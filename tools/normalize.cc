/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/Common.h"
#include "mirtk/Options.h"
#include "mirtk/IOConfig.h"

#include "mirtk/Image.h"
#include "mirtk/ImageHistogram1D.h"
//#include "mirtk/SegmentationFunction.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/NormalizeNyul.h"

using namespace mirtk;

void PrintHelp(const char *name)
{
  std::cout << "usage: " << name << " <target> <source> <output_source> [options]" << std::endl;
  std::cout << "       " << name << " <input> <output> -equalize <padding>" << std::endl;
  std::cout << std::endl;
  std::cout << "Normalizes the intensity distribution of an image to be similar to" << std::endl;
  std::cout << "the intensity distribution of a given reference image. Moreover," << std::endl;
  std::cout << "this tool can be used to equalize the histograms of either a single" << std::endl;
  std::cout << "given image or two images using the same transfer function." << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -Tp <value>   Target padding value" << std::endl;
  std::cout << "  -Sp <value>   Source padding value" << std::endl;
  std::cout << "  -piecewise    Use a piecewise linear function as suggested by Nyul et al." << std::endl;
  std::cout << "  -equalize <padding> [<target_output>]   Equalize histograms before normalization." << std::endl;
  PrintCommonOptions(std::cout);
}

int main(int argc, char **argv)
{
  // Positional arguments
  REQUIRES_POSARGS(2);
  InitializeIOLibrary();

  const char *target_name   = NULL;
  const char *target_output = NULL;
  const char *source_name   = NULL;
  const char *source_output = NULL;

  if (NUM_POSARGS == 2) {
    source_name   = POSARG(1);
    source_output = POSARG(2);
  } else if (NUM_POSARGS == 3) {
    target_name   = POSARG(1);
    source_name   = POSARG(2);
    source_output = POSARG(3);
  } else {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Optional arguments
  double source_padding   = MIN_GREY;
  double target_padding   = MIN_GREY;
  double equalize_padding = .0;
  bool   equalize         = false;
  bool   piecewise        = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-Tp")) target_padding = atof(ARGUMENT);
    else if (OPTION("-Sp")) source_padding = atof(ARGUMENT);
    else if (OPTION("-piecewise")) piecewise = true;
    else if (OPTION("-equalize")) {
      equalize = true;
      equalize_padding = atof(ARGUMENT);
      if (NUM_POSARGS > 2) {
        target_output = (HAS_ARGUMENT ? ARGUMENT : NULL);
      }
    }
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (target_name == NULL) {
    if (target_output) {
      PrintHelp(EXECNAME);
      exit(1);
    }
    equalize = true;
  }

  // Read input image(s)
  RealImage target;
  if (target_name) {
    if (verbose) std::cout << "Reading target image ... ", std::cout.flush();
    target.Read(target_name);
    if (verbose) std::cout << "done" << std::endl;
  }

  if (verbose) std::cout << "Reading source image ... ", std::cout.flush();
  RealImage source(source_name);
  if (verbose) std::cout << "done" << std::endl;

  // Equalize histograms
  if (equalize) {
    double min, max;
    RealImage *reference = (target.IsEmpty() ? &source : &target);

    if (verbose) {
      std::cout << "Equalize histogram" << ((reference == &target) ? "s" : "") << " ... ";
      std::cout.flush();
    }

    reference->GetMinMaxAsDouble(min, max);
    if (min < equalize_padding) min = equalize_padding;

		ImageHistogram1D<RealPixel> histogram;
		histogram.Evaluate(reference, equalize_padding);
		//NormalizeNyul::histogramImage(&histogram, reference, equalize_padding);
		histogram.Equalize(min, max);
		histogram.BackProject(reference);

    if (reference != &source) {
      source.GetMinMaxAsDouble(min, max);
      if (min < equalize_padding) min = equalize_padding;
      histogram.Evaluate(&source, equalize_padding);
      //NormalizeNyul::histogramImage(&histogram, &source, equalize_padding);
      histogram.Equalize(min, max);
      histogram.BackProject(&source);
    }

    if (verbose) std::cout << "done" << std::endl;
  }

  // Stop if source image equalization is done only
  if (target.IsEmpty()) {
    source.Write(source_output);
    return 0;
  }

  if (verbose) std::cout << "Normalize histogram ... ";

  // Normalize histogram
  if (piecewise) {
    NormalizeNyul nn(source, target);
    nn.SetPadding(source_padding, target_padding);
    nn.Run();
    source = nn.GetOutput();
  } else {
    double a, b, cov, var, x_avg, y_avg, x, y, z;

	  int n = 0;
	  x_avg = 0;
	  y_avg = 0;
	  for (int k = 0; k < source.GetZ(); k++)
		for (int j = 0; j < source.GetY(); j++)
		for (int i = 0; i < source.GetX(); i++) {
      x = i; y = j; z = k;
      source.ImageToWorld(x,y,z);
      target.WorldToImage(x,y,z);
      x = round(x); y = round(y); z = round(z);
      if (x >= 0 && x < target.GetX() &&
          y >= 0 && y < target.GetY() &&
          z >= 0 && z < target.GetZ()) {
        if ((source(i, j, k) > source_padding) && (target(x,y,z) > target_padding)) {
          n++;
          x_avg += source(i, j, k);
          y_avg += target(x, y, z);
        }
      }
	  }
	  if (n == 0) {
      std::cout << "failed" << std::endl;
      std::cerr << EXECNAME << ": Number of samples should be larger than zero" << std::endl;
      exit(1);
	  }

	  cov = 0;
	  var = 0;
	  for (int k = 0; k < source.GetZ(); k++)
    for (int j = 0; j < source.GetY(); j++)
    for (int i = 0; i < source.GetX(); i++) {
      x = i; y = j; z = k;
      source.ImageToWorld(x,y,z);
      target.WorldToImage(x,y,z);
      x = round(x); y = round(y); z = round(z);
      if (x >= 0 && x < target.GetX() &&
          y >= 0 && y < target.GetY() &&
          z >= 0 && z < target.GetZ()) {
        if ((source(i, j, k) > source_padding) && (target(x, y, z) > target_padding)) {
          cov += (source(i, j, k) - x_avg) * (target(x, y, z) - y_avg);
          var += (source(i, j, k) - x_avg) * (source(i, j, k) - x_avg);
        }
      }
	  }
	  cov /= n;
	  var /= n;
	  b = cov / var;
	  a = y_avg - b * x_avg;

    if (verbose) {
      std::cout << "scaling = " << b << ", offset = " << b;
    }

	  for (int k = 0; k < source.GetZ(); k++)
		for (int j = 0; j < source.GetY(); j++)
		for (int i = 0; i < source.GetX(); i++) {
      if (source(i, j, k) > source_padding) {
        source(i, j, k) = round(a + b * source(i, j, k));
      }
	  }

    if (verbose) std::cout << "done" << std::endl;
  }

  // Write normalized image
  source.Write(source_output);
  if (target_output) target.Write(target_output);
}

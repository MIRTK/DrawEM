/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Antonios Makropoulos
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
#include "mirtk/GenericImage.h"

#include <algorithm>   
#include <vector>      

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input> <options>" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
 	std::cout << "  Calculates statistics by filtering with a kernel." << std::endl;
	std::cout << std::endl;
	std::cout << "Options: " << std::endl;
	std::cout << "  -kernel <number>       kernel size: number^3 (number must be even, and >=3 !), default: 3 " << std::endl;
	std::cout << std::endl;
	std::cout << "Operations: " << std::endl;
	std::cout << "  -min <output>          calculate min" << std::endl;
	std::cout << "  -max <output>          calculate max" << std::endl;
	std::cout << "  -mean <output>         calculate mean" << std::endl;
	std::cout << "  -median <output>       calculate median" << std::endl;
	std::cout << "  -std <output>          calculate std" << std::endl;
	std::cout << "  -std-median <output>   calculate std based on median" << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------


int main(int argc, char **argv)
{
  EXPECTS_POSARGS(1);

  InitializeIOLibrary();
  RealImage img(POSARG(1));

  const char *mask_name        = nullptr;
  const char *min_name         = nullptr;
  const char *max_name         = nullptr;
  const char *mean_name        = nullptr;
  const char *median_name      = nullptr;
  const char *std_name         = nullptr;
  const char *std_median_name  = nullptr;
  bool        have_output_name = false;
  int         kernel           = 1;

  for (ALL_OPTIONS) {
    if (OPTION("-kernel")){
      PARSE_ARGUMENT(kernel);
      if (kernel % 2 == 0 || kernel < 3) {
        FatalError("Invalid -kernel width");
      }
      kernel = (kernel - 1) / 2;
    }
    else if (OPTION("-min")) {
      min_name = ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-max")){
      max_name=ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-mean")){
      mean_name=ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-median")){
      median_name=ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-std")){
      std_name=ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-std_median")){
      std_median_name=ARGUMENT;
      have_output_name = true;
    }
    else if (OPTION("-mask")) {
      mask_name = ARGUMENT;
    }
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }
  if (!have_output_name) {
    FatalError("At least one output file name option must be given!");
  }

  RealImage min_img, max_img, mean_img, median_img, std_img, std_median_img;
  if (min_name) min_img.Initialize(img.Attributes());
  if (max_name) max_img.Initialize(img.Attributes());
  if (mean_name) mean_img.Initialize(img.Attributes());
  if (median_name) median_img.Initialize(img.Attributes());
  if (std_name) std_img.Initialize(img.Attributes());
  if (std_median_name) std_median_img.Initialize(img.Attributes());

  BinaryImage mask;
  if (mask_name) mask.Read(mask_name);

  Array<RealPixel> values;
  values.reserve((2 * kernel + 1) * (2 * kernel + 1) * (2 * kernel + 1));

  int x,x1,x2,xn, y,y1,y2,yn, z,z1,z2,zn;
  double minv, maxv, meanv, medianv, sum2;

  for (z = 0; z < img.Z(); ++z)
  for (y = 0; y < img.Y(); ++y)
  for (x = 0; x < img.X(); ++x) {

    x1 = max(0,x-kernel);
    y1 = max(0,y-kernel);
    z1 = max(0,z-kernel);
    x2 = min(x+kernel, img.X()-1);
    y2 = min(y+kernel, img.Y()-1);
    z2 = min(z+kernel, img.Z()-1);

    values.clear();
    for (xn = x1; xn <= x2; ++xn)
    for (yn = y1; yn <= y2; ++yn)
    for (zn = z1; zn <= z2; ++zn) {
      if (!mask_name || mask(xn, yn, zn) != 0) {
        values.push_back(img.Get(xn, yn, zn));
      }
    }

    if (values.empty()) {
      if (min_name) min_img(x, y, z) = 0.;
      if (max_name) max_img(x, y, z) = 0.;
      if (median_name) median_img(x, y, z) = 0.;
      if (mean_name) mean_img(x, y, z) = 0.;
      if (std_name) std_img(x, y, z) = 0.;
      if (std_median_name) std_median_img(x, y, z) = 0;
      continue;
    }

    minv = maxv = values.front(), meanv = 0.;
    for (auto v : values) {
      minv = min(minv, v);
      maxv = max(maxv, v);
      meanv += v;
    }
    meanv /= values.size();

    if (min_name) min_img(x, y, z) = minv;
    if (max_name) max_img(x, y, z) = maxv;
    if (mean_name) mean_img(x, y, z) = meanv;
    if (median_name || std_median_name) {
      medianv = NthElement(values, static_cast<int>(values.size()) / 2);
      if (median_name) median_img(x,y,z) = medianv;
    }
    if (std_name) {
      sum2 = 0.;
      for (auto v : values) {
        v -= meanv;
        sum2 += v * v;
      }
      std_img(x, y, z) = sqrt(sum2 / static_cast<double>(values.size()));
    }
    if (std_median_name) {
      sum2 = 0.;
      for (auto v : values) {
        v -= medianv;
        sum2 += v * v;
      }
      std_median_img(x, y, z) = sqrt(sum2 / static_cast<double>(values.size()));
    }
  }

  if (min_name) min_img.Write(min_name);
  if (max_name) max_img.Write(max_name);
  if (mean_name) mean_img.Write(mean_name);
  if (median_name) median_img.Write(median_name);
  if (std_name) std_img.Write(std_name);
  if (std_median_name) std_median_img.Write(std_median_name);

  return 0;
}

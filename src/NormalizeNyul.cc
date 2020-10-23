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

#include "mirtk/NormalizeNyul.h"

#include "mirtk/ImageHistogram1D.h"
#include "mirtk/CommonExport.h"


namespace mirtk {


// Verbose flag declared in Options.h of Common module
MIRTK_Common_EXPORT int verbose;


NormalizeNyul::NormalizeNyul(RealImage source, RealImage target)
{
	_source = source;
	_target = target;
	_source_padding = 0.;
	_target_padding = 0.;
}


void NormalizeNyul::SetMask(RealImage source_mask, RealImage target_mask)
{
	RealPixel * iPtr, *mPtr;
	iPtr = _source.GetPointerToVoxels();
	mPtr = source_mask.GetPointerToVoxels();
	for (int i = 0; i < _source.GetNumberOfVoxels(); ++i) {
		if (*mPtr < 1) {
			*iPtr = voxel_cast<RealPixel>(_source_padding);
		}
		iPtr++;
		mPtr++;
	}
	iPtr = _target.GetPointerToVoxels();
	mPtr = source_mask.GetPointerToVoxels();
	for (int i = 0; i < _target.GetNumberOfVoxels(); ++i) {
		if (*mPtr < 1) {
			*iPtr = voxel_cast<RealPixel>(_target_padding);
		}
		iPtr++;
		mPtr++;
	}
}

RealImage NormalizeNyul::GetOutput()
{
	return _source;
}

RealImage NormalizeNyul::GetTarget()
{
	return _target;
}

void NormalizeNyul::SetPadding(double source_padding, double target_padding)
{
	_source_padding = source_padding;
	_target_padding = target_padding;
}

void NormalizeNyul::Run()
{
  const int    hist_bins = 4000;
  const int    steps     = 10;
	const double cut_off   = .01;
  const double step      = (1.0 - 2.0*cut_off / 100.0) / steps;

	ImageHistogram1D<RealPixel> src_hist_s, trg_hist_s;
	src_hist_s.PutNumberOfBins(hist_bins);
	trg_hist_s.PutNumberOfBins(hist_bins);
	src_hist_s.Evaluate(&_source, static_cast<RealPixel>(_source_padding));
	trg_hist_s.Evaluate(&_target, static_cast<RealPixel>(_target_padding));

	const double src_min = src_hist_s.GetMin();
  const double src_max = src_hist_s.GetMax();
  const double trg_min = trg_hist_s.GetMin();
  const double trg_max = trg_hist_s.GetMax();

	//provide mapping for background
  Array<double> src_levels_s;
  Array<double> trg_levels_s;
	src_levels_s.push_back(src_min);
	trg_levels_s.push_back(trg_min);
  if (verbose) {
    cout << "[ min ] " << src_levels_s[0] << " => " << trg_levels_s[0] << endl;
    cout << "nr tgt samples: " << trg_hist_s.NumberOfSamples() << ". nr src samples: " << src_hist_s.NumberOfSamples() << endl;
    cout << "nr tgt bins: " << trg_hist_s.NumberOfBins() << ". nr src bins: " << src_hist_s.NumberOfBins() << endl;
  }
	for (int i = 0; i <= steps; ++i) {
		double pct =  i * step + cut_off/100.;
		double src_lev_s = src_hist_s.CDFToVal(pct);
		double trg_lev_s = trg_hist_s.CDFToVal(pct);
    if (verbose) {
      cout << "step " << i << " perc: " << pct << " src perc: " << src_lev_s << " trg perc: " << trg_lev_s << endl;
    }
		if (trg_lev_s-trg_levels_s[trg_levels_s.size()-1]<(trg_max-trg_min)/100000.0) {
		  cerr << "Warning: " << pct*100 <<" percentile collapses in target, skipping" << endl;
		} else {
		  src_levels_s.push_back(src_lev_s);
		  trg_levels_s.push_back(trg_lev_s);
      if (verbose) {
        cout << "[ " << pct*100.0 << " ] " << src_levels_s[src_levels_s.size() - 1] << " => " << trg_levels_s[trg_levels_s.size() - 1] << endl;
      }
		}
	}

	//provide mapping for upper range
	src_levels_s.push_back(src_max);
	trg_levels_s.push_back(trg_max);
  if (verbose) {
    cout << "[ max ] " << src_levels_s[src_levels_s.size() - 1] << " => " << trg_levels_s[trg_levels_s.size() - 1] << endl;
  }
  if (verbose) {
    cout << "Recalculating intensities..." << flush;
  }
	RealPixel *ptr = _source.GetPointerToVoxels();
	for (int i = 0; i < _source.GetNumberOfVoxels(); ++i, ++ptr) {
	  //use LUT to map the intensities
	  double input  = *ptr;
	  double output = input;
    size_t bin;
    for (bin = 0; bin < src_levels_s.size(); bin++) {
      if (input <= src_levels_s[bin]) break;
    }
    if (input <= this->_source_padding) {
      output = this->_source_padding;
    } else if (bin == 0) { // first bin ?
      output = trg_levels_s[0];
    } else if (bin >= (src_levels_s.size() - 1)) {
		  output = trg_levels_s[trg_levels_s.size()-1];
    } else {
      output = (input - src_levels_s[bin - 1])
             / (src_levels_s[bin] - src_levels_s[bin - 1])
             *(trg_levels_s[bin] - trg_levels_s[bin - 1]) + trg_levels_s[bin - 1];
    }
	  *ptr = output;
	}
  if (verbose) {
    cout << " done" << endl;
  }
}


} // namespace mirtk

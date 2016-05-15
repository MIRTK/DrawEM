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


#ifndef MIRTKNORMALIZENYUL_H_
#define MIRTKNORMALIZENYUL_H_

#include "mirtk/Image.h"
//#include "mirtk/Registration.h"
#include "mirtk/ImageHistogram1D.h"


namespace mirtk {

class NormalizeNyul{

private:
	RealImage _target;
	RealImage _source;
	int _source_padding;
	int _target_padding;

public:
	NormalizeNyul(RealImage source, RealImage target);
	void SetMask(RealImage source_mask, RealImage target_mask);
	void SetPadding(int source_padding, int target_padding);
	static void histogramImage(Histogram1D<RealPixel>* histogram, RealImage* image, double padding);
	void Run();
	RealImage GetOutput();
	RealImage GetTarget();
};

}
#endif

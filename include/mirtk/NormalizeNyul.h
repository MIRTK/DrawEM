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


#ifndef MIRTK_NORMALIZENYUL_H_
#define MIRTK_NORMALIZENYUL_H_

#include "mirtk/GenericImage.h"


namespace mirtk {


/**
 * This filter implements the intensity normalization algorithm published in
 *
 * Nyul, L.G.; Udupa, J.K.; Xuan Zhang, "New variants of a method of MRI scale standardization",
 * Medical Imaging, IEEE Transactions on , vol.19, no.2, pp.143-150, Feb 2000
 * http://dx.doi.org/10.1109/42.836373 "
 *
 * Source code adapted from an implementation by Vladimir Fonov in EZminc (https://github.com/vfonov/EZminc)
 */
class NormalizeNyul
{
private:
	RealImage _target;
	RealImage _source;
  double _source_padding;
  double _target_padding;

public:
	NormalizeNyul(RealImage source, RealImage target);
	void SetMask(RealImage source_mask, RealImage target_mask);
	void SetPadding(double source_padding, double target_padding);
	void Run();
	RealImage GetOutput();
	RealImage GetTarget();
};


} // namespace mirtk

#endif // MIRTK_NORMALIZENYUL_H_

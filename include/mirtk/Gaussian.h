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


#ifndef _MIRTKGAUSSIAN_H

#define _MIRTKGAUSSIAN_H

#include "mirtk/Image.h"

/**

multivariate gaussian probability distribution

 */

namespace mirtk {

class Gaussian : public Object
{
	mirtkObjectMacro(Gaussian);

protected:

	double _mi;
	double _sigma;
	double _norm;

public:

	void Initialise(const double &mi, const double &sigma);
	double Evaluate(const double &x);
	double GetNorm();
};

inline double Gaussian::Evaluate(const double &x)
{
	return _norm * exp(-((x - _mi) * (x - _mi)) / (2.0 * _sigma));
}

inline double Gaussian::GetNorm()
{
	return _norm;
}

}
#endif

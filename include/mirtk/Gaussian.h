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

#ifndef _MIRTK_GAUSSIAN_H
#define _MIRTK_GAUSSIAN_H

#include "mirtk/Object.h"
#include "mirtk/Math.h"


namespace mirtk {

/**
 * multivariate gaussian probability distribution
 */
class Gaussian : public Object
{
	mirtkObjectMacro(Gaussian);

protected:

	double _mean;
	double _var;
	double _norm;

public:

	void Initialise(double mean, double var);
	double Evaluate(double x) const;
	double GetNorm() const;
};


inline double Gaussian::Evaluate(double x) const
{
  x -= _mean;
	return _norm * exp(-.5 * x*x / _var);
}

inline double Gaussian::GetNorm() const
{
	return _norm;
}


} // namespace mirtk

#endif // _MIRTK_GAUSSIAN_H

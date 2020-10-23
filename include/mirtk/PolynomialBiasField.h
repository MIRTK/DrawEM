/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Christian Ledig
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


#ifndef MIRTKPOLYNOMIALBIASFIELD_H_
#define MIRTKPOLYNOMIALBIASFIELD_H_

#include "mirtk/BiasField.h"

namespace mirtk {

class PolynomialBiasField : public BiasField
{
	mirtkObjectMacro(PolynomialBiasField);

private:
	int _dop;
	double* _coeff;
	int _numOfCoefficients;

public:
	PolynomialBiasField();

	/**
	 * @param dop max degree of polynomial
	 */
	PolynomialBiasField(const GreyImage &image, int dop);
	~PolynomialBiasField();

	/// Calculate weighted least square fit of polynomial to data
	virtual void WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no);

	double Bias(double, double, double);

	double Approximate(double *, double *, double *, double *, int);
	void Interpolate(double* dbias);

	/// Subdivide FFD
	void Subdivide();

	/// Reads FFD from file
	void Read (char *);

	/// Writes FFD to file
	virtual void Write(char *);

	/// Print info
	virtual void Print();

private:
	double evaluatePolynomial(double x, double y, double z);
	int getNumberOfCoefficients(int dop);
};

}

#endif /* MIRTKPOLYNOMIALBIASFIELD_H_ */

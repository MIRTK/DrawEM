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


#include "mirtk/Matrix.h"
#include "mirtk/Vector.h"
#include "mirtk/Cifstream.h"

#include "mirtk/PolynomialBiasField.h"
#include "mirtk/BiasField.h"


namespace mirtk {

PolynomialBiasField::PolynomialBiasField()
{
}

PolynomialBiasField::PolynomialBiasField(const GreyImage &image, int dop)
{
	_dop = dop;
	_numOfCoefficients = getNumberOfCoefficients(dop);
	_coeff = new double[_numOfCoefficients];
	memset( _coeff, 0, sizeof(double) * _numOfCoefficients );
}

PolynomialBiasField::~PolynomialBiasField()
{
	delete _coeff;
}

//void PolynomialBiasField::WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)
//{
//	// just consider each eachs voxel...
//	int each = 1;
//	no /= each;
//
//	Matrix A(no, _numOfCoefficients);
//	Vector vecB(no);
//
//	std::cout << "numOfCoefficients: " << _numOfCoefficients << std::endl;
//	std::cout << "num of voxels: " << no << std::endl;
//
//	for( int rr = 0; rr < no; ++rr )
//	{
//		int r = rr * each;
//		vecB.Put(rr,bias[r]);
//
//		double x = x1[r];
//		double y = y1[r];
//		double z = z1[r];
//
//		int c = 0;
//
//		double cur_x = pow(x, _dop);
//		for( int xd = _dop; xd >= 0; --xd)
//		{
//			double cur_y = pow(y,_dop-xd);
//			for( int yd = _dop-xd; yd >= 0; --yd)
//			{
//				double cur_z = cur_x*cur_y * pow(z,_dop-xd-yd);
//				for( int zd = _dop-xd-yd; zd >= 0; --zd)
//				{
//					A.Put(rr,c, cur_z);
//					c++;
//					if( z != 0) cur_z /= z;
//				}
//				if( y != 0) cur_y /= y;
//			}
//			if( x != 0) cur_x /= x;
//		}
//	}
//
//	Matrix AtW = A;
//	AtW.Transpose();
//
//	for( int r = 0; r < _numOfCoefficients; ++r)
//	{
//		for( int c= 0; c < no; ++c )
//		{
//			AtW.Put(r,c, AtW(r,c) * weights[c*each] );
//		}
//	}
//
// 	Vector rightSide = AtW * vecB;
// 	Matrix leftSide = AtW * A;
//
// 	leftSide.Invert();
//
//	Vector vecC = leftSide * rightSide;
//
//	for( int r = 0; r < _numOfCoefficients; ++r )
//	{
//		_coeff[r] = vecC.Get(r);
//	}
//}

// ALTERNATIVE IMPLEMENTATION USING symmetry of A'WA, should be by a factor of 2 faster than above implementation
void PolynomialBiasField::WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)
{
	// just consider each eachs voxel...
	int each = 3;
	no /= each;

	Matrix A(no, _numOfCoefficients);
	Vector vecB(no);

	std::cout << "numOfCoefficients: " << _numOfCoefficients << std::endl;
	std::cout << "num of voxels: " << no << std::endl;

	double* AtWA = new double[_numOfCoefficients * _numOfCoefficients];
	memset(AtWA, 0, sizeof(double) * _numOfCoefficients * _numOfCoefficients );

	double* Basis = new double[_numOfCoefficients];

	for( int rr = 0; rr < no; ++rr )
	{
		int r = rr * each;
		double weight = weights[rr*each];
		vecB.Put(rr,bias[r]);

		double x = x1[r];
		double y = y1[r];
		double z = z1[r];

		int c = 0;

		double cur_x = 1.0;
		double cur_y = 1.0;

		for( int xd = 0; xd <= _dop; ++xd )
		{
			cur_y = 1.0;
			for( int yd = 0; yd <= _dop-xd; ++yd )
			{
				double tmp = cur_x*cur_y;
				for( int zd = 0; zd <= _dop-xd-yd; ++zd )
				{
					A.Put(rr,c, tmp);
					Basis[c] = tmp;
					c++;
					tmp *= z;
				}
				cur_y *= y;
			}
			cur_x *= x;
		}

		double* Aptr = (double*) AtWA;
		double* Basisptr1 = (double*) Basis;
		double* Basisptr2;

		for( int j2 = 0; j2 < _numOfCoefficients; j2++, Basisptr1++)
		{
			Basisptr2= &Basis[j2];
			Aptr= &AtWA[j2+j2*_numOfCoefficients];
			for(int i2=j2; i2<_numOfCoefficients; i2++, Aptr++, Basisptr2++){
				(*Aptr)+=(*Basisptr2)*(weight)*(*Basisptr1);
			}
		}
	}

	Matrix AtW = A;
	AtW.Transpose();

	for( int r = 0; r < _numOfCoefficients; ++r)
	{
		for( int c= 0; c < no; ++c )
		{
			AtW.Put(r,c, AtW(r,c) * weights[c*each] );
		}
	}

	Matrix leftSide(_numOfCoefficients, _numOfCoefficients);
	for(int j2=0; j2<_numOfCoefficients; j2++)
	{
		for(int i2=j2; i2<_numOfCoefficients; i2++)
		{
			leftSide.Put(i2,j2,(double)(AtWA[i2+j2*_numOfCoefficients]));
			leftSide.Put(j2,i2,(double)(AtWA[i2+j2*_numOfCoefficients]));
		}
	}

	Vector rightSide = AtW * vecB;
	leftSide.Invert();

	Vector vecC = leftSide * rightSide;

	for( int r = 0; r < _numOfCoefficients; ++r )
	{
        _coeff[r] = vecC.Get(r);
    }
}

double PolynomialBiasField::evaluatePolynomial(double x, double y, double z)
{
	double res = 0;
	int n = 0;

	double cur_x = 1.0;
	double cur_y = 1.0;

	for( int xd = 0; xd <= _dop; ++xd )
	{
		cur_y = 1.0;
		for( int yd = 0; yd <= _dop-xd; ++yd )
		{
			double tmp = cur_x * cur_y;
			for( int zd = 0; zd <= _dop-xd-yd; ++zd )
			{
				res +=   tmp * _coeff[n];
				n++;
				tmp *= z;
			}
			cur_y *= y;
		}
		cur_x *= x;
	}
	return res;
}


double PolynomialBiasField::Bias(double x, double y, double z)
{
	return evaluatePolynomial(x, y, z);
}

void PolynomialBiasField::Interpolate(double* dbias)
{

}

void PolynomialBiasField::Subdivide()
{

}

void PolynomialBiasField::Read(char *name)
{
	unsigned int magic_no;
	unsigned int trans_type;

	// Open file
	Cifstream from;
	from.Open(name);

	// Read magic no. for transformations
	from.ReadAsUInt(&magic_no, 1, 0);

	// Read transformation type
	from.ReadAsUInt(&trans_type, 1);

	if ((magic_no != MIRTKBIASFIELD_MAGIC) && (trans_type != MIRTKBIASFIELD_POLYNOMIAL)) {
		std::cerr << "PolynomialBiasField::Read: File format not recognized" << std::endl;
		exit(1);
	}

	if( _coeff ) delete _coeff;

	// Write data
	from.ReadAsInt(&_dop, 1);
	from.ReadAsInt(&_numOfCoefficients, 1);
	_coeff = new double[_numOfCoefficients];

	from.ReadAsDouble(_coeff, _numOfCoefficients);

	// Close file stream
	from.Close();
}


void PolynomialBiasField::Write(char *name)
{
	// Open file
	Cofstream to;
	to.Open(name);

	// Write magic no. for transformations
	unsigned int magic_no = MIRTKBIASFIELD_MAGIC;
	to.WriteAsUInt(&magic_no, 1, 0);

	// Write transformation type
	unsigned int trans_type = MIRTKBIASFIELD_POLYNOMIAL;
	to.WriteAsUInt(&trans_type, 1);

	// Write data
	to.WriteAsInt(&_dop, 1);
	to.WriteAsInt(&_numOfCoefficients, 1);
	to.WriteAsDouble(_coeff, _numOfCoefficients);

	// Close file stream
	to.Close();
}

void PolynomialBiasField::Print()
{
	std::cerr<<std::endl<<"Polynomial Bias Field info:" << std::endl;
	// Write no. of control points
	std::cout << "Control points: " << _x << " x " << _y << " x " << _z << std::endl;
	std::cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << std::endl;
	std::cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << " " << std::endl;
	std::cout << "Orientation: " << _xaxis[0] << " " << _xaxis[1] << " "
			<< _xaxis[2] << " " << _yaxis[0] << " " << _yaxis[1] << " "
			<< _yaxis[2] << " " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << std::endl;
	std::cout << _numOfCoefficients << " Coefficients:" << std::endl;

	int n;
	for (n=0; n<_numOfCoefficients; n++) {
		std::cout << _coeff[n] << " ";
	}
	std::cout<<std::endl;
}


double PolynomialBiasField::Approximate(double *x1, double *y1, double *z1, double *bias, int no)
{
	return .0;
}

int PolynomialBiasField::getNumberOfCoefficients(int dop)
{
	int n = 0;
	for( int x = dop; x >=0; --x )
	{
		for( int y = dop-x; y >=0; --y )
		{
			for( int z = dop-x-y; z >=0; --z )
			{
				n++;
			}
		}
	}
	return n;
}

}

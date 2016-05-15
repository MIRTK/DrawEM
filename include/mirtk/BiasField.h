/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Maria Murgasova
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

#ifndef _MIRTKBIASFIELD_H

#define _MIRTKBIASFIELD_H

//#include "mirtk/Geometry.h"
#include "mirtk/Point.h"
#include "mirtk/Matrix.h"

#include "mirtk/Transformation.h"

#define MIRTKBIASFIELD_MAGIC            815008

#define MIRTKBIASFIELD_BSPLINE          1
#define MIRTKBIASFIELD_POLYNOMIAL       2

namespace mirtk {

class BiasField : public Object
{
	mirtkAbstractMacro(BiasField);

protected:

	/// Number of control points in x
	int _x;

	/// Number of control points in y
	int _y;

	/// Number of control points in z
	int _z;

	/// Spacing of control points in x (in mm)
	double _dx;

	/// Spacing of control points in y (in mm)
	double _dy;

	/// Spacing of control points in z (in mm)
	double _dz;

	/// Direction of x-axis
	double _xaxis[3];

	/// Direction of y-axis
	double _yaxis[3];

	/// Direction of z-axis
	double _zaxis[3];

	/// Origin
	Point _origin;

	/// Transformation matrix from lattice coordinates to world coordinates
	Matrix _matL2W;

	/// Transformation matrix from world coordinates to lattice coordinates
	Matrix _matW2L;

	/// Displacement at the control points
	double ***_data;

	/// Allocate memory for control points
	static double ***Allocate  (double ***, int, int, int);

	/// Deallocate memory for control points
	static double ***Deallocate(double ***, int, int, int);

	/// Update transformation matrix
	virtual void UpdateMatrix();

public:

	/** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
	virtual double Approximate(double *, double *, double *, double *, int) = 0;

	/// Calculate weighted least square fit to data
	virtual void WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)=0;


	/** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param bias The bias at each control point. */
	virtual void Interpolate(double* bias) = 0;

	/// Subdivide FFD
	virtual void Subdivide() = 0;

	/// Returns the of control points in x
	virtual int GetX() const;

	/// Returns the of control points in y
	virtual int GetY() const;

	/// Returns the of control points in z
	virtual int GetZ() const;

	/// Returns the number of parameters of the transformation
	virtual int NumberOfDOFs() const;

	/// Get the control point spacing (in mm)
	virtual void GetSpacing(double &, double &, double &) const;

	/// Put orientation of free-form deformation
	virtual void  PutOrientation(double *, double *, double *);

	/// Get orientation of free-form deformation
	virtual void  GetOrientation(double *, double *, double *) const;

	/// Puts a control point value
	virtual void   Put(int, double);

	/// Puts a control point value
	virtual void   Put(int, int, int, double);

	/// Gets a control point value
	virtual double Get(int) const;

	/// Gets a control point value
	virtual double Get(int, int, int) const;

	/// Transforms world coordinates (in mm) to FFD coordinates
	virtual void WorldToLattice(double &, double &, double &) const;

	/// Transforms world coordinates (in mm) to FFD coordinates
	virtual void WorldToLattice(Point &) const;

	/// Transforms FFD coordinates to world coordinates (in mm)
	virtual void LatticeToWorld(double &, double &, double &) const;

	/// Transforms FFD coordinates to world coordinates (in mm)
	virtual void LatticeToWorld(Point &) const;

	/// Transforms index of control points to FFD coordinates
	virtual void IndexToLattice(int index, int& i, int& j, int& k) const;

	/// Transforms  FFD coordinates to index of control point
	virtual int  LatticeToIndex(int i, int j, int k) const;

	/// Returns the control point location (in mm)
	virtual void  ControlPointLocation(int, double &, double &, double &) const;

	/// Returns the control point location (in mm)
	virtual Point ControlPointLocation(int) const;

	/// Put the bounding box for FFD (in mm)
	virtual void PutBoundingBox(Point, Point);

	/// Put the bounding box for FFD (in mm)
	virtual void PutBoundingBox(double, double, double,
			double, double, double);

	/// Returns the bounding box for FFD (in mm)
	virtual void BoundingBox(Point &, Point &) const;

	/// Returns the bounding box for FFD (in mm)
	virtual void BoundingBox(double &, double &, double &,
			double &, double &, double &) const;

	/** Returns the bounding box for a control point (in mm). The last
	 *  parameter specifies what fraction of the bounding box to return. The
	 *  default is 1 which equals 100% of the bounding box.
	 */
	virtual void BoundingBox(int, Point &, Point &, double = 1) const;

	/** Returns the bounding box for a control point (in mm). The last
	 *  parameter specifies what fraction of the bounding box to return. The
	 *  default is 1 which equals 100% of the bounding box.
	 */
	virtual void BoundingBox(int, double &, double &, double &,
			double &, double &, double &, double = 1) const;

	/** Returns the bounding box for a control point (in pixels). The last
	 *  parameter specifies what fraction of the bounding box to return. The
	 *  default is 1 which equals 100% of the bounding box.
	 */
	virtual void BoundingBox(GreyImage *, int, int &, int &, int &,
			int &, int &, int &, double = 1) const;

	/// Calculate value of bias field at a point
	virtual double Bias(double, double, double) = 0;

	/// Reads a transformation from a file (abstract)
	virtual void Read (char *) = 0;

	/// Writes a transformation to a file (abstract)
	virtual void Write(char *) = 0;

	/// Prints info (abstract)
	virtual void Print() = 0;

};

inline int BiasField::GetX() const
{
	return _x;
}

inline int BiasField::GetY() const
{
	return _y;
}

inline int BiasField::GetZ() const
{
	return _z;
}

inline int BiasField::NumberOfDOFs() const
{
	return _x*_y*_z;
}

inline double BiasField::Get(int index) const
{
	int i, j, k;

	if (index >= _x*_y*_z) {
		std::cerr << "BiasField::Get: No such dof" << std::endl;
		exit(1);
	}
	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
	return _data[k][j][i];
}

inline double BiasField::Get(int i, int j, int k) const
{
	if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z)) {
		std::cerr << "BiasField::Get: No such dof" << std::endl;
		exit(1);
	}
	return _data[k][j][i];
}

inline void BiasField::Put(int index, double x)
{
	Point p;
	int i, j, k;

	if (index >= _x*_y*_z) {
		std::cerr << "BiasField::Put: No such dof" << std::endl;
		exit(1);
	}
	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
	_data[k][j][i] = x;
}

inline void BiasField::Put(int i, int j, int k, double x)
{
	if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z)) {
		std::cerr << "BiasField::Put: No such dof" << std::endl;
		exit(1);
	}
	_data[k][j][i] = x;
}

inline void BiasField::GetSpacing(double &dx, double &dy, double &dz) const
{
	dx = _dx;
	dy = _dy;
	dz = _dz;
}

inline void BiasField::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
{
	_xaxis[0] = xaxis[0];
	_xaxis[1] = xaxis[1];
	_xaxis[2] = xaxis[2];
	_yaxis[0] = yaxis[0];
	_yaxis[1] = yaxis[1];
	_yaxis[2] = yaxis[2];
	_zaxis[0] = zaxis[0];
	_zaxis[1] = zaxis[1];
	_zaxis[2] = zaxis[2];

	this->UpdateMatrix();
}

inline void BiasField::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
{
	xaxis[0] = _xaxis[0];
	xaxis[1] = _xaxis[1];
	xaxis[2] = _xaxis[2];
	yaxis[0] = _yaxis[0];
	yaxis[1] = _yaxis[1];
	yaxis[2] = _yaxis[2];
	zaxis[0] = _zaxis[0];
	zaxis[1] = _zaxis[1];
	zaxis[2] = _zaxis[2];
}

inline void BiasField::PutBoundingBox(double x1, double y1, double z1, double x2, double y2, double z2)
{
	// Initialize control point domain
	_origin._x = (x2 + x1) / 2.0;
	_origin._y = (y2 + y1) / 2.0;
	_origin._z = (z2 + z1) / 2.0;

	double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
	double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
	double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
	x1 = a;
	y1 = b;
	z1 = c;
	a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
	b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
	c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
	x2 = a;
	y2 = b;
	z2 = c;

	// Initialize control point spacing

	_x = round((x2-x1)/_dx) +1;
	_y = round((y2-y1)/_dy) +1;
	_z = round((z2-z1)/_dz) +1;
	_dx = (x2 - x1) / (_x - 1);
	_dy = (y2 - y1) / (_y - 1);
	if (z2 > z1) {
		_dz = (z2 - z1) / (_z - 1);
	} else {
		_dz = 1;
	}

	this->UpdateMatrix();
}

inline void BiasField::PutBoundingBox(Point p1, Point p2)
{
	this->PutBoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

inline void BiasField::WorldToLattice(double &x, double &y, double &z) const
{
	double a, b, c;

	// Pre-multiply point with transformation matrix
	a = _matW2L(0, 0)*x+_matW2L(0, 1)*y+_matW2L(0, 2)*z+_matW2L(0, 3);
	b = _matW2L(1, 0)*x+_matW2L(1, 1)*y+_matW2L(1, 2)*z+_matW2L(1, 3);
	c = _matW2L(2, 0)*x+_matW2L(2, 1)*y+_matW2L(2, 2)*z+_matW2L(2, 3);

	// Copy result back
	x = a;
	y = b;
	z = c;
}

inline void BiasField::WorldToLattice(Point &p) const
{
	this->WorldToLattice(p._x, p._y, p._z);
}

inline void BiasField::LatticeToWorld(double &x, double &y, double &z) const
{
	double a, b, c;

	// Pre-multiply point with transformation matrix
	a = _matL2W(0, 0)*x+_matL2W(0, 1)*y+_matL2W(0, 2)*z+_matL2W(0, 3);
	b = _matL2W(1, 0)*x+_matL2W(1, 1)*y+_matL2W(1, 2)*z+_matL2W(1, 3);
	c = _matL2W(2, 0)*x+_matL2W(2, 1)*y+_matL2W(2, 2)*z+_matL2W(2, 3);

	// Copy result back
	x = a;
	y = b;
	z = c;
}

inline void BiasField::LatticeToWorld(Point &p) const
{
	this->LatticeToWorld(p._x, p._y, p._z);
}

inline void BiasField::IndexToLattice(int index, int& i, int& j, int& k) const
{

	if (index >= _x*_y*_z) {
		index -= _x*_y*_z;
		if (index >= _x*_y*_z) {
			index -= _x*_y*_z;
		}
	}
	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
}

inline int BiasField::LatticeToIndex(int i, int j, int k) const
{
	return i * _y * _z + j * _z + k;
}

//#include "BSplineBiasField.h"

}

#endif

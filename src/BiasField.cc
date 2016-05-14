/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Daniel Rueckert
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


#include "mirtk/BiasField.h"

namespace mirtk {

void BiasField::UpdateMatrix()
{
	// Update image to world coordinate system matrix
	_matL2W.Ident();
	_matL2W(0, 0) = _xaxis[0];
	_matL2W(1, 0) = _xaxis[1];
	_matL2W(2, 0) = _xaxis[2];
	_matL2W(0, 1) = _yaxis[0];
	_matL2W(1, 1) = _yaxis[1];
	_matL2W(2, 1) = _yaxis[2];
	_matL2W(0, 2) = _zaxis[0];
	_matL2W(1, 2) = _zaxis[1];
	_matL2W(2, 2) = _zaxis[2];

	Matrix tmp1(4, 4);
	tmp1.Ident();
	tmp1(0, 3) = - (_x - 1) / 2.0;
	tmp1(1, 3) = - (_y - 1) / 2.0;
	tmp1(2, 3) = - (_z - 1) / 2.0;

	Matrix tmp2(4, 4);
	tmp2.Ident();
	tmp2(0, 0) = _dx;
	tmp2(1, 1) = _dy;
	tmp2(2, 2) = _dz;

	Matrix tmp3(4, 4);
	tmp3.Ident();
	tmp3(0, 3) = _origin._x;
	tmp3(1, 3) = _origin._y;
	tmp3(2, 3) = _origin._z;

	_matL2W = tmp3 * (_matL2W * (tmp2 * tmp1));

	// Update world to image coordinate system matrix
	_matW2L.Ident();
	_matW2L(0, 0) = _xaxis[0];
	_matW2L(0, 1) = _xaxis[1];
	_matW2L(0, 2) = _xaxis[2];
	_matW2L(1, 0) = _yaxis[0];
	_matW2L(1, 1) = _yaxis[1];
	_matW2L(1, 2) = _yaxis[2];
	_matW2L(2, 0) = _zaxis[0];
	_matW2L(2, 1) = _zaxis[1];
	_matW2L(2, 2) = _zaxis[2];

	tmp1.Ident();
	tmp1(0, 3) = (_x - 1) / 2.0;
	tmp1(1, 3) = (_y - 1) / 2.0;
	tmp1(2, 3) = (_z - 1) / 2.0;

	tmp2.Ident();
	tmp2(0, 0) = 1.0 / _dx;
	tmp2(1, 1) = 1.0 / _dy;
	tmp2(2, 2) = 1.0 / _dz;

	tmp3.Ident();
	tmp3(0, 3) = - _origin._x;
	tmp3(1, 3) = - _origin._y;
	tmp3(2, 3) = - _origin._z;

	_matW2L = tmp1 * (tmp2 * (_matW2L * tmp3));
}

double ***BiasField::Allocate(double ***data, int x, int y, int z)
{
	int i, j, k;

	if ((x == 0) || (y == 0) || (z == 0)) {
		return NULL;
	}

	if ((data = new double **[z+8]) == NULL) {
		std::cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
		std::cerr << z << "\n";
		exit(1);
	}
	data += 4;

	if ((data[-4] = new double *[(z+8)*(y+8)]) == NULL) {
		std::cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
		std::cerr << z << "\n";
		exit(1);
	}
	data[-4] += 4;

	for (i = -3; i < z+4; i++) {
		data[i] = data[i-1] + (y+8);
	}

	if ((data[-4][-4] = new double[(z+8)*(y+8)*(x+8)]) == NULL) {
		std::cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
		std::cerr << z << "\n";
		exit(1);
	}
	data[-4][-4] += 4;

	for (i = -4; i < z+4; i++) {
		for (j = -4; j < y+4; j++) {
			data[i][j] = data[-4][-4] + ((i+4)*(y+8)+(j+4))*(x+8);
		}
	}

	for (i = -4; i < z+4; i++) {
		for (j = -4; j < y+4; j++) {
			for (k = -4; k < x+4; k++) {
				data[i][j][k] = 0;
			}
		}
	}

	return data;
}

double ***BiasField::Deallocate(double ***data, int, int, int)
{
	if (data != NULL) {
		delete [](data[-4][-4]-4);
		delete [](data[-4]-4);
		delete [](data-4);
	}
	return NULL;
}

Point BiasField::ControlPointLocation(int index) const
{
	int i, j, k;

	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
	Point  p(i, j, k);
	this->LatticeToWorld(p);
	return p;
}

void BiasField::ControlPointLocation(int index, double &x, double &y, double &z) const
{
	x = index/(_y*_z);
	y = index%(_y*_z)/_z;
	z = index%(_y*_z)%_z;
	this->LatticeToWorld(x, y, z);
}

void BiasField::BoundingBox(Point &p1, Point &p2) const
{
	p1._x = 0;
	p1._y = 0;
	p1._z = 0;
	p2._x = _x - 1;
	p2._y = _y - 1;
	p2._z = _z - 1;
	this->LatticeToWorld(p1);
	this->LatticeToWorld(p2);
}

void BiasField::BoundingBox(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2) const
{
	x1 = 0;
	y1 = 0;
	z1 = 0;
	x2 = _x - 1;
	y2 = _y - 1;
	z2 = _z - 1;
	this->LatticeToWorld(x1, y1, z1);
	this->LatticeToWorld(x2, y2, z2);
}

void BiasField::BoundingBox(int index, Point &p1, Point &p2, double fraction) const
{
	int i, j, k;

	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
	p1 = Point(i-2*fraction, j-2*fraction, k-2*fraction);
	p2 = Point(i+2*fraction, j+2*fraction, k+2*fraction);
	this->LatticeToWorld(p1);
	this->LatticeToWorld(p2);
}

void BiasField::BoundingBox(int index, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, double fraction) const
{
	x1 = index/(_y*_z)-2*fraction;
	y1 = index%(_y*_z)/_z-2*fraction;
	z1 = index%(_y*_z)%_z-2*fraction;
	x2 = index/(_y*_z)+2*fraction;
	y2 = index%(_y*_z)/_z+2*fraction;
	z2 = index%(_y*_z)%_z+2*fraction;
	this->LatticeToWorld(x1, y1, z1);
	this->LatticeToWorld(x2, y2, z2);
}

void BiasField::BoundingBox(GreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	double x1, y1, z1, x2, y2, z2;

	// Calculate bounding box in world coordinates
	this->BoundingBox(index, x1, y1, z1, x2, y2, z2, fraction);

	// Transform world coordinates to image coordinates
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);

	// Calculate bounding box in image coordinates
	i1 = (x1 < 0) ? 0 : int(x1)+1;
	j1 = (y1 < 0) ? 0 : int(y1)+1;
	k1 = (z1 < 0) ? 0 : int(z1)+1;
	i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
	j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
	k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

}

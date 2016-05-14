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


#include "mirtk/BiasCorrection.h"
namespace mirtk {

BiasCorrection::BiasCorrection()
{
	// Set parameters
	_Padding   = MIN_GREY;

	// Set inputs
	_target    = NULL;
	_reference = NULL;

	// Set output
	_biasfield = NULL;

	// mask
	_mask = NULL;
}

BiasCorrection::~BiasCorrection()
{
}

void BiasCorrection::Initialize()
{
}

void BiasCorrection::SetMask( ByteImage* imagePtr)
{
	_mask = imagePtr;
}

void BiasCorrection::Finalize()
{
}

void BiasCorrection::Run()
{
	int i, j, k, n;
	RealPixel *ptr2target, *ptr2ref, *ptr2w;

	if (_reference == NULL) {
		std::cerr << "BiasCorrection::Run: Filter has no reference input" << std::endl;
		exit(1);
	}

	if (_target == NULL) {
		std::cerr << "BiasCorrection::Run: Filter has no target input" << std::endl;
		exit(1);
	}

	if (_biasfield == NULL) {
		std::cerr << "BiasCorrection::Run: Filter has no transformation output" << std::endl;
		exit(1);
	}

	// Do the initial set up for all levels
	this->Initialize();

	// Compute no of unpadded voxels
	n = 0;
	ptr2target = _target->GetPointerToVoxels();
	BytePixel *pm = _mask->GetPointerToVoxels();
	for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
		if (*ptr2target != _Padding && *pm == 1 ) {
			n++;
		}
		pm++;
		ptr2target++;
	}

	double *x = new double[n];
	double *y = new double[n];
	double *z = new double[n];
	double *b = new double[n];
	double *w = new double[n];
	n = 0;
	ptr2target = _target->GetPointerToVoxels();
	ptr2ref    = _reference->GetPointerToVoxels();
	ptr2w      = _weights->GetPointerToVoxels();
	pm = _mask->GetPointerToVoxels();
	for (k = 0; k < _target->GetZ(); k++) {
		for (j = 0; j < _target->GetY(); j++) {
			for (i = 0; i < _target->GetX(); i++) {
				if (*ptr2target != _Padding && *pm == 1) {
					x[n] = i;
					y[n] = j;
					z[n] = k;
					_target->ImageToWorld(x[n], y[n], z[n]);
					b[n] = *ptr2target - (double) *ptr2ref;
                    w[n] = *ptr2w;
					n++;
				}
				pm++;
				ptr2target++;
				ptr2ref++;
				ptr2w++;
			}
		}
	}

	std::cout << "Computing bias field ... ";
	std::cout.flush();
	// _biasfield->Approximate(x, y, z, b, n);
	_biasfield->WeightedLeastSquares(x, y, z, b, w, n);
	std::cout << "done" << std::endl;

	delete []x;
	delete []y;
	delete []z;
	delete []b;

	// Do the final cleaning up for all levels
	this->Finalize();
}

void BiasCorrection::Apply(RealImage &image)
{
	int i, j, k;
	double x, y, z, bias;

	// Default is the target image
	image = *_target;

	for (k = 0; k < image.GetZ(); k++) {
		for (j = 0; j < image.GetY(); j++) {
			for (i = 0; i < image.GetX(); i++) {
				x = i;
				y = j;
				z = k;
				image.ImageToWorld(x, y, z);
				bias = _biasfield->Bias(x, y, z);
                if (_target->Get(i, j, k) != _Padding ) {
                //if (_mask->Get(i, j, k)  && _target->Get(i, j, k) != _Padding) {
					//image(i, j, k) = round(image(i, j, k) - bias);
					image(i, j, k) = image(i, j, k) - bias;
					//if( image(i, j, k) < 0.0 ) image(i, j, k) = 0.0;
				}
			}
		}
	}
}

void BiasCorrection::ApplyToImage(RealImage &image)
{
	int i, j, k;
	double x, y, z, bias;
	std::cerr<<"Applying bias ...";
	//std::cerr<<_biasfield;

	for (k = 0; k < image.GetZ(); k++) {
		for (j = 0; j < image.GetY(); j++) {
			for (i = 0; i < image.GetX(); i++) {
				x = i;
				y = j;
				z = k;
				image.ImageToWorld(x, y, z);
				bias = _biasfield->Bias(x, y, z);
				if (image.Get(i, j, k) != _Padding) {
					image(i, j, k) = round(image(i, j, k) - bias);
				}
			}

		}
	}

	std::cerr<<"done."<<std::endl;
}

void BiasCorrection::ApplyToImage(GreyImage &image)
{
	int i, j, k;
	double x, y, z, bias;

	for (k = 0; k < image.GetZ(); k++) {
		for (j = 0; j < image.GetY(); j++) {
			for (i = 0; i < image.GetX(); i++) {
				x = i;
				y = j;
				z = k;
				image.ImageToWorld(x, y, z);
				bias = _biasfield->Bias(x, y, z);
				if (image.Get(i, j, k) != _Padding) {
					image(i, j, k) = round(image(i, j, k) / exp(bias/1000));
				}
			}
		}
	}
}

}

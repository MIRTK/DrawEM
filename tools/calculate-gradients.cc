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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/GaussianBlurring.h"

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input> <output> <sigma> [<options>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
 	std::cout << "  A blur with S.D. sigma is appled before the gradient is estimated" << std::endl;
	std::cout << std::endl;
	std::cout << "Optional: " << std::endl;
	std::cout << std::endl;
	std::cout << "  -sep basename    Write separate components to basename-x.nii.gz, basename-y.nii.gz, basename-z.nii.gz" << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------


int main(int argc, char **argv)
{
        REQUIRES_POSARGS(3);

        InitializeIOLibrary();

	RealImage input;
	RealImage output;
	RealImage gradX;
	RealImage gradY;
	RealImage gradZ;
	double voxelx, voxely, voxelz;
	double dx, dy, dz;
	int i, j, k, xdim, ydim, zdim;
	RealPixel tmp = 0;
	double ssq, sigma;
	Matrix w2i;

	char *input_name  = POSARG(1);
	char *output_name = POSARG(2);
	sigma = atof(POSARG(3));
	char *sepBasename = NULL;

	    for (ALL_OPTIONS) {
		if (OPTION("-sep")){
		    sepBasename = ARGUMENT;
		}
		else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
	    }


	// Read input
	input.Read(input_name);
	output.Read(input_name);
	gradX.Read(input_name);
	gradY.Read(input_name);
	gradZ.Read(input_name);

	// Blur input image
	if (sigma > 0){
		GaussianBlurring<RealPixel> gaussianBlurring( sigma );
		gaussianBlurring.Input (&input);
		gaussianBlurring.Output(&input);
		gaussianBlurring.Run();
	}

	w2i = input.GetWorldToImageMatrix();

	xdim = input.GetX();
	ydim = input.GetY();
	zdim = input.GetZ();

	input.GetPixelSize(&voxelx, &voxely, &voxelz);
	voxelx = abs(voxelx);
	voxely = abs(voxely);
	voxelz = abs(voxelz);

	for (k = 0; k < zdim; ++k) {
		for (j = 0; j < ydim; ++j){
			for (i = 0; i < xdim; ++i){

				if (i == 0 || i == xdim - 1 ||
						j == 0 || j == ydim - 1 ||
						k == 0 || k == zdim - 1) {

					gradX.Put(i, j, k, 0);
					gradY.Put(i, j, k, 0);
					gradZ.Put(i, j, k, 0);
					output.Put(i, j, k, 0);
					continue;

				}

				ssq = 0.0;

				// Gradient in image coordinates.
				dx = (input.Get(i+1, j  , k  ) - input.Get(i-1, j  , k  )) / 2.0;
				dy = (input.Get(i  , j+1, k  ) - input.Get(i  , j-1, k  )) / 2.0;
				dz = (input.Get(i  , j  , k+1) - input.Get(i  , j  , k-1)) / 2.0;

				// Convert to world coordinates.
				// Assume scalar function is f
				// World coords are (w1, w2, w3) and image coordinates are (im1, im2, im3).
				// e.g for w1,
				//
				// d f    d f   d im1   d f   d im2   d f   d im3
				// ---- = ----- ----- + ----- ----- + ----- -----
				// d w1   d im1 d w1    d im2 d w1    d im3 d w1
				//
				// I.e. can obtain world gradient by multiplying image gradient by world
				// to image matrix which expresses image coordinates as a function of
				// world coordinates.
				tmp = w2i(0, 0) * dx + w2i(0, 1) * dy + w2i(0, 2) * dz;
				gradX.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				tmp = w2i(1, 0) * dx + w2i(1, 1) * dy + w2i(1, 2) * dz;
				gradY.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				tmp = w2i(2, 0) * dx + w2i(2, 1) * dy + w2i(2, 2) * dz;
				gradZ.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				output.Put(i, j, k, sqrt(ssq));

			}
		}
	}


	output.Write(output_name);

	if (sepBasename != NULL){
		char buffer[250];

		sprintf(buffer, "%s-x.nii.gz", sepBasename);
		gradX.Write(buffer);

		sprintf(buffer, "%s-y.nii.gz", sepBasename);
		gradY.Write(buffer);

		sprintf(buffer, "%s-z.nii.gz", sepBasename);
		gradZ.Write(buffer);

	}
}


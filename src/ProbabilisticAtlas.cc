/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Wenzhe Shi
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

#include "mirtk/ProbabilisticAtlas.h"
#include "mirtk/GenericImage.h"


namespace mirtk {


ProbabilisticAtlas::ProbabilisticAtlas()
{
	_number_of_voxels=0;
	_number_of_tissues=0;
}

void ProbabilisticAtlas::SwapImages(int a, int b)
{
	if( a >= _number_of_tissues || b >= _number_of_tissues )
	{
		cerr << "cannot swap images, index out of bounds!" << endl;
		return;
	}
	RealImage tmpimage = _images[a];
	_images[a] = _images[b];
	_images[b] = tmpimage;
	RealPixel* tmpptr = _pointers[a];
	_pointers[a] = _pointers[b];
	_pointers[a] = tmpptr;
}

void ProbabilisticAtlas::AddImage(RealImage image)
{
	if (_images.size() == 0) {
		_number_of_voxels = image.GetNumberOfVoxels();
	} else {
		if (_number_of_voxels != image.GetNumberOfVoxels()) {
			cerr << "Image sizes mismatch" << endl;
			exit(1);
		}
	}
	_images.push_back(image);
	_pointers.push_back(image.GetPointerToVoxels());
	_number_of_tissues = static_cast<int>(_images.size());
}

void ProbabilisticAtlas::NormalizeAtlas(RealImage background)
{
	int i, j;
	RealPixel norm;

	// Add extra image
	if (_number_of_tissues == 0) {
		cerr << "ProbabilisticAtlas::NormalizeAtlas: No probability maps found" << endl;
		exit(1);
	} else {
		this->AddImage(background);
		_number_of_tissues = static_cast<int>(_images.size());

	}

	// normalize atlas to 0 to 1
	this->First();
	for (i = 0; i < _number_of_voxels; i++) {
		norm = 0;
		for (j = 0; j < _number_of_tissues; j++) {
			norm += *(_pointers[j]);
		}
		if (norm>0) {
			for (j = 0; j < _number_of_tissues; j++)
				*(_pointers[j]) = *(_pointers[j])/norm;
		} else {
			for (j = 0; j < _number_of_tissues-1; j++)
				*(_pointers[j]) = 0;
			*(_pointers[_number_of_tissues-1]) = 1;
		}
		this->Next();
	}
}

void ProbabilisticAtlas::NormalizeAtlas()
{
	int i, j;
	RealPixel norm;

	// Add extra image
	if (_number_of_tissues == 0) {
		cerr << "ProbabilisticAtlas::NormalizeAtlas: No probability maps found" << endl;
		exit(1);
	}
	_number_of_tissues = static_cast<int>(_images.size());


	// normalize atlas to 0 to 1
	this->First();
	for (i = 0; i < _number_of_voxels; i++) {
		norm = 0;

		for (j = 0; j < _number_of_tissues-1; j++) {
			if (*(_pointers[j]) < .0) {
				//cerr << " atlas prob < 0: " << *(_pointers[j]) << "  tissue " << j << " voxel" << i << endl;
				*(_pointers[j]) = .0;
			}
		}

		for (j = 0; j < _number_of_tissues; j++) {
			norm += *(_pointers[j]);
		}
		if (norm>0) {
      for (j = 0; j < _number_of_tissues; j++) {
        *(_pointers[j]) = *(_pointers[j]) / norm;
      }
		} else {
			int vx,vy,vz;
			int index=i;
			vz = index / (_images[0].GetX() * _images[0].GetY());
			index -= vz * (_images[0].GetX() * _images[0].GetY());
			vx = index % _images[0].GetX();
			vy = index / _images[0].GetX();
      for (j = 0; j < _number_of_tissues; j++) {
        *(_pointers[j]) = 0;
      }
		}
		this->Next();
	}
}


void ProbabilisticAtlas::AddBackground()
{
	int i, j;
	RealPixel norm, min, max;

	// Add extra image
	if (_number_of_tissues == 0) {
		cerr << "ProbabilisticAtlas::AddBackground: No probability maps found" << endl;
		exit(1);
	} else {
		this->AddImage(_images[0]);
		_number_of_tissues = static_cast<int>(_images.size());
	}

	// normalize atlas to 0 to 1
	RealImage other;
	other = _images[0];
	for (i = 1; i < _number_of_tissues-1; i++) {
		cerr<<i<<endl;
		other += _images[i];
	}
	other.GetMinMax(&min, &max);

	this->First();
	for (i = 0; i < _number_of_voxels; i++) {
		norm = 0;
		for (j = 0; j < _number_of_tissues-1; j++) {
			*(_pointers[j]) = (*(_pointers[j]) - min) / (max - min);
			norm += *(_pointers[j]);
		}
		double value = 1 - norm;
		if (value < 0) value = 0;
		if (value > 1) value = 1;
		*(_pointers[_number_of_tissues-1]) = value;
		this->Next();
	}
}

void ProbabilisticAtlas::AddProbabilityMaps(int n, RealImage **atlas)
{
	for (int i = 0; i < n; i++) {
		this->AddImage(*(atlas[i]));
	}
};

void ProbabilisticAtlas::Write(int i, const char *filename)
{
	if  (i < _number_of_tissues) {
		_images[i].Write(filename);
	} else {
		cerr << "ProbabilisticAtlas::Write: No such probability map" << endl;
		exit(1);
	}
}

RealImage ProbabilisticAtlas::ComputeHardSegmentation()
{
	int tissue;
	double max;
	_segmentation = _images[0];
	First();
	RealPixel *ptr = _segmentation.GetPointerToVoxels();
	for (int i = 0; i < _number_of_voxels; ++i) {
		max = 0.;
		tissue = -1;
		for (int j = 0; j < GetNumberOfTissues(); ++j) {
			if (GetValue(j) > max) {
				tissue = j;
				max = GetValue(j);
			}
		}
		*ptr = tissue;
		Next();
		ptr++;
	}
	return _segmentation;
}

RealImage ProbabilisticAtlas::GetImage(int i)
{
	return _images[i];
}

void ProbabilisticAtlas::ExtractLabel(int label, RealImage& image)
{
	RealPixel *ptr = _segmentation.GetPointerToVoxels();
	RealPixel *ptr_im = image.GetPointerToVoxels();
	for (int i = 0; i < _number_of_voxels; ++i) {
		if (*ptr == label) *ptr_im = 1;
		else *ptr_im = 0;
		ptr++;
		ptr_im++;
	}
}

void ProbabilisticAtlas::WriteHardSegmentation(const char *filename)
{
	_segmentation.Write(filename);
}


} // namespace mirtk

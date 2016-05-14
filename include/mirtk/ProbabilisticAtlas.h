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


#ifndef _MIRTKPROBABILISTICATLAS_H

#define _MIRTKPROBABILISTICATLAS_H

#include "mirtk/Object.h"

#include "mirtk/Image.h"

#include <vector>


/**

Atlas probability map class

 */
using namespace std;
namespace mirtk {


class ProbabilisticAtlas : public Object
{
	mirtkObjectMacro(ProbabilisticAtlas);

	// Vector of probability maps
	vector<RealImage> _images;

	// Vector of pointers to probability maps
	vector<RealPixel *> _pointers;

	///Number of voxels
	int _number_of_voxels;

	///Number of voxels
	int _number_of_tissues;

	/// Hard segmentation
	RealImage _segmentation;


public:

	/// Constructor
	ProbabilisticAtlas();

	/// swap images within atlas
	void SwapImages(int, int);

	/// Adds new image (channel)
	void AddImage(RealImage image);

	/// Moves pointers in all images to the first voxel
	void First();

	/// Moves pointers in all images to the next voxel
	void Next();

	///Returns intensity value at pointer
	RealPixel GetValue(unsigned int channel);

	///Returns intensity value at position
	RealPixel GetValue(int x, int y, int z, unsigned int tissue);

	///Sets intensity value at pointer
	void SetValue(unsigned int tissue, RealPixel value);

	///Sets intensity value at position
	void SetValue(int x, int y, int z, unsigned int tissue, RealPixel value);

	/// Returns number of voxels
	int GetNumberOfVoxels();

	/// Add probability maps
	void AddProbabilityMaps(int, RealImage **atlas);

	/// Computes probability map for the last tissue type
	void AddBackground();

	/// normalize atlas without adding background
	void NormalizeAtlas();

	/// normalize atlas while adding background
	void NormalizeAtlas(RealImage background);

	///Returns number of tissues
	int GetNumberOfTissues();

	/// Write
	void Write(int, const char *);

	/// Computes hard segmentation
	RealImage ComputeHardSegmentation();

	/// Extract a label from hard segmentation
	void ExtractLabel(int label, RealImage& image);

	/// Writes hard segmentation into a file
	void WriteHardSegmentation(const char *filename);

	RealImage GetImage(int);

};

inline void ProbabilisticAtlas::First()
{
	unsigned int i;
	for (i=0; i<_pointers.size(); i++) _pointers[i] = _images[i].GetPointerToVoxels();
}

inline void ProbabilisticAtlas::Next()
{
	unsigned int i;
	for (i=0; i<_pointers.size(); i++) _pointers[i]++;
}

inline RealPixel ProbabilisticAtlas::GetValue(unsigned int channel)
{
	if (channel < _images.size()) return *_pointers[channel];
	else {
		std::cerr << "Channel identificator " << channel <<" out of range." <<std::endl;
		exit(1);
	}
}

inline RealPixel ProbabilisticAtlas::GetValue(int x, int y, int z, unsigned int channel)
{
	if (channel < _images.size()) return _images[channel].Get(x,y,z);
	else {
		std::cerr << "Channel identificator " << channel <<" out of range." <<std::endl;
		exit(1);
	}
}

inline void ProbabilisticAtlas::SetValue(unsigned int channel, RealPixel value)
{
	if (channel < _images.size()) *_pointers[channel] = value;
	else {
		std::cerr << "Channel identificator " << channel << " out of range." <<std::endl;
		exit(1);
	}
}

inline void ProbabilisticAtlas::SetValue(int x, int y, int z, unsigned int channel, RealPixel value)
{
	if (channel < _images.size()) _images[channel].Put( x, y, z, value);
	else {
		std::cerr << "Channel identificator " << channel <<" out of range." <<std::endl;
		exit(1);
	}
}

inline int ProbabilisticAtlas::GetNumberOfVoxels()
{
	return _number_of_voxels;
}

inline int ProbabilisticAtlas::GetNumberOfTissues()
{
	return _number_of_tissues;
}

}

#endif

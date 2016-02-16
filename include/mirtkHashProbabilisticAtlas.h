/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Antonios Makropoulos
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


#ifndef _mirtkHashProbabilisticAtlas_H

#define _mirtkHashProbabilisticAtlas_H

#include <mirtkObject.h>

#include <mirtkImage.h>
#include <mirtkHashImage.h>

#include <vector>


/**

Atlas probability mapnr class

 */
using namespace std;
namespace mirtk{


class mirtkHashProbabilisticAtlas : public Object
{
	typedef HashImage<int> HashIntegerImage;

	mirtkObjectMacro(mirtkHashProbabilisticAtlas);

	// Vector of probability maps
	vector<HashRealImage*> _images;

	// Number of voxels
	int _number_of_voxels;

	// Number of maps
	int _number_of_maps;

	// Hard segmentation
	HashIntegerImage *_segmentation;

	// _position of iterator
	int _position;

	// whether we have background
 	bool _has_background;

public:

	// Constructor
	mirtkHashProbabilisticAtlas();

	// Destructor
	~mirtkHashProbabilisticAtlas();

	// Copy operator
	mirtkHashProbabilisticAtlas& operator=(const mirtkHashProbabilisticAtlas &atlas);

	// swap images within atlas
	void SwapImages(int, int);

	// Adds new image
	template <class ImageType>
    void AddImage(ImageType image);

	// Moves pointers in all images to the first voxel
	void First();

	// Moves pointers in all images to the next voxel
	void Next();

	// Returns intensity value at pointer
	RealPixel GetValue(unsigned int mapnr);

	// Returns intensity value at _position
	RealPixel GetValue(int x, int y, int z, unsigned int mapnr);

	// Sets intensity value at pointer
	void SetValue(unsigned int mapnr, RealPixel value);

	// Sets intensity value at _position
	void SetValue(int x, int y, int z, unsigned int mapnr, RealPixel value);

	// Returns number of voxels
    int GetNumberOfVoxels() const;

	// Add probability maps
	template <class ImageType>
    void AddProbabilityMaps(int, ImageType **atlas);

	// Add background map as (1-rest maps)
	void AddBackground();

	// Add background map 
	template <class ImageType>
    void AddBackground(ImageType image);

	// Whether it has background
	bool HasBackground() const;

	// normalize atlas
	void NormalizeAtlas();

	//Returns number of maps
	int GetNumberOfMaps() const;

	// Write
	void Write(int, const char *);

	// Computes hard segmentation
	HashIntegerImage ComputeHardSegmentation();

	// Extract a label from hard segmentation
	void ExtractLabel(int label, HashByteImage &image);

	// Writes hard segmentation into a file
    void WriteHardSegmentation(const char *filename);

    // Get image
    HashRealImage GetImage(unsigned int mapnr) const;

    // Get data
    HashRealImage::DataIterator Begin(unsigned int mapnr) const;
    HashRealImage::DataIterator End(unsigned int mapnr) const;

};

inline void mirtkHashProbabilisticAtlas::First(){
	_position=0;
}

inline void mirtkHashProbabilisticAtlas::Next(){
	_position++;
}

inline RealPixel mirtkHashProbabilisticAtlas::GetValue(unsigned int mapnr){
	if (mapnr < _images.size()) return _images[mapnr]->Get(_position);
	else {
		cerr << "map identificator " << mapnr <<" out of range." <<endl;
		exit(1);
	}
}

inline RealPixel mirtkHashProbabilisticAtlas::GetValue(int x, int y, int z, unsigned int mapnr){
	if (mapnr < _images.size()) return _images[mapnr]->Get(x,y,z);
	else {
		cerr << "map identificator " << mapnr <<" out of range." <<endl;
		exit(1);
	}
}

inline void mirtkHashProbabilisticAtlas::SetValue(unsigned int mapnr, RealPixel value){
	if (mapnr < _images.size()) _images[mapnr]->Put(_position, value);
	else {
		cerr << "map identificator " << mapnr << " out of range." <<endl;
		exit(1);
	}
}

inline void mirtkHashProbabilisticAtlas::SetValue(int x, int y, int z, unsigned int mapnr, RealPixel value){
	if (mapnr < _images.size()) _images[mapnr]->Put(x,y,z, value);
	else {
		cerr << "map identificator " << mapnr <<" out of range." <<endl;
		exit(1);
	}
}

inline int mirtkHashProbabilisticAtlas::GetNumberOfVoxels() const{
	return _number_of_voxels;
}

inline int mirtkHashProbabilisticAtlas::GetNumberOfMaps() const{
	return _number_of_maps;
}

inline bool mirtkHashProbabilisticAtlas::HasBackground() const{
	return _has_background;
}

template <class ImageType>
inline void mirtkHashProbabilisticAtlas::AddBackground(ImageType image){
	AddImage(image);
	_has_background=true;
}

inline HashRealImage mirtkHashProbabilisticAtlas::GetImage(unsigned int mapnr) const{
    if (mapnr < _images.size()) return *_images[mapnr];
    else {
        cerr << "map identificator " << mapnr <<" out of range." <<endl;
        exit(1);
    }
}

inline HashRealImage::DataIterator mirtkHashProbabilisticAtlas::Begin(unsigned int mapnr) const
{
    if (mapnr < _images.size()) return _images[mapnr]->Begin();
    else {
        cerr << "map identificator " << mapnr <<" out of range." <<endl;
        exit(1);
    }
}
inline HashRealImage::DataIterator mirtkHashProbabilisticAtlas::End(unsigned int mapnr) const
{
    if (mapnr < _images.size()) return _images[mapnr]->End();
    else {
        cerr << "map identificator " << mapnr <<" out of range." <<endl;
        exit(1);
    }
}

}

#endif

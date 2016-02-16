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

#include <mirtkHashProbabilisticAtlas.h>

namespace mirtk{

mirtkHashProbabilisticAtlas::mirtkHashProbabilisticAtlas(){
	_number_of_voxels = 0;
	_number_of_maps = 0;
	_position = 0;
	_has_background = false;
	_segmentation = NULL;
}

mirtkHashProbabilisticAtlas::~mirtkHashProbabilisticAtlas(){
	if (_segmentation) delete _segmentation;
	for(int i=0; i<_images.size(); i++) delete _images[i];
}

mirtkHashProbabilisticAtlas& mirtkHashProbabilisticAtlas::operator=(const mirtkHashProbabilisticAtlas &atlas)
{
  if (this != &atlas) {
	if (_segmentation) delete _segmentation;
	for(int i=0; i<_images.size(); i++) delete _images[i];
	int N=atlas.GetNumberOfMaps();
	for(int i=0; i<N; i++) AddImage(atlas.GetImage(i));
	_has_background = atlas.HasBackground();
  }
  return *this;
}

void mirtkHashProbabilisticAtlas::SwapImages(int a, int b){
	if( a >= _number_of_maps || b >= _number_of_maps ){
		cerr << "cannot swap images, index out of bounds!" << endl;
		return;
	}
	HashRealImage *tmpimage = _images[a];
	_images[a] = _images[b];
	_images[b] = tmpimage;
}

template <class ImageType>
void mirtkHashProbabilisticAtlas::AddImage(ImageType image){
	if (_images.size() == 0) {
		_number_of_voxels = image.GetNumberOfVoxels();
	} else {
		if (_number_of_voxels != image.GetNumberOfVoxels()) {
			cerr << "Image sizes mismatch" << endl;
			exit(1);
		}
	}
	_images.push_back(new HashRealImage(image));
	if(_has_background) SwapImages(_images.size()-2, _images.size()-1);
	_number_of_maps = _images.size();

}

void mirtkHashProbabilisticAtlas::NormalizeAtlas(){
	int i, j;

	// Add extra image
	if (_number_of_maps == 0) {
		cerr << "mirtkHashProbabilisticAtlas::NormalizeAtlas: No probability maps found" << endl;
		exit(1);
	} 

	// normalize atlas to 0 to 1
	this->First();
	RealPixel norm;
	RealPixel values[_number_of_maps];
	for (i = 0; i < _number_of_voxels; i++) {
		norm = 0;
		for (j = 0; j < _number_of_maps; j++) {
			values[j]= _images[j]->Get(i);
			if(values[j]<0){
				_images[j]->Put(i, 0);
				values[j]=0;
			}else norm += values[j];
		}
		if (norm>0) {
			for (j = 0; j < _number_of_maps; j++)
				_images[j]->Put(i, values[j]/norm);
		} else {
			for (j = 0; j < _number_of_maps-1; j++)
				_images[j]->Put(i, 0);
			if(_has_background) _images[_number_of_maps-1]->Put(i, 1);
		}
		this->Next();
	}
}


void mirtkHashProbabilisticAtlas::AddBackground(){
	int i, j;
	// Add extra image
	if (_number_of_maps == 0) {
		cerr << "mirtkHashProbabilisticAtlas::AddBackground: No probability maps found" << endl;
		exit(1);
	} 
	this->AddImage(*(_images[0]));

	// normalize atlas to 0 to 1
	RealPixel norm, min, max, div, value;
	RealPixel values[_number_of_maps];

	HashRealImage *other;
	other = _images[0];
	for (i = 1; i < _number_of_maps-1; i++) {
		(*other) += (*_images[i]);
	}
	other->GetMinMax(&min, &max);
	div=(max - min);

	this->First();
	for (i = 0; i < _number_of_voxels; i++) {
		norm = 0;
		for (j = 0; j < _number_of_maps-1; j++) {
			values[j]=(_images[j]->Get(i)- min) / div;
			_images[j]->Put(i, values[j]);
			norm += values[j];
		}
		value = 1 - norm;
		if (value < 0) value = 0;
		if (value > 1) value = 1;
		 _images[_number_of_maps-1]->Put(i, value);
		this->Next();
	}

	_has_background=true;
}

template <class ImageType>
void mirtkHashProbabilisticAtlas::AddProbabilityMaps(int n, ImageType **atlas){
	for (int i = 0; i < n; i++) {
		this->AddImage(*(atlas[i]));
	}
}

void mirtkHashProbabilisticAtlas::Write(int i, const char *filename){
	if  (i < _number_of_maps) {
		_images[i]->Write(filename);
	} else {
		cerr << "mirtkHashProbabilisticAtlas::Write: No such probability map" << endl;
		exit(1);
	}
}

HashImage<int> mirtkHashProbabilisticAtlas::ComputeHardSegmentation(){
	int i, mapnr, j = 0;
    RealPixel max = 0;

	if(!_segmentation) _segmentation = new HashImage<int>(_images[0]->Attributes());
    (*_segmentation) = -1;

	First();
	for (i = 0; i < _number_of_voxels; i++) {
		max = 0;
		mapnr = -1;
		for (j=0; j< GetNumberOfMaps(); j++) {
			if (_images[j]->Get(i) > max) {
				mapnr = j;
				max = _images[j]->Get(i);
			}
        }
		_segmentation->Put(i, mapnr);
		Next();
	}
	return *_segmentation;
}

void mirtkHashProbabilisticAtlas::ExtractLabel(int label, HashByteImage& image){
	if(!_segmentation) ComputeHardSegmentation();
	for (int i = 0; i < _number_of_voxels; i++) {
		if (_segmentation->Get(i) == label) image.Put(i, 1);
		else image.Put(i, 0);
	}
}

void mirtkHashProbabilisticAtlas::WriteHardSegmentation(const char *filename){
	_segmentation->Write(filename);
}


template void mirtkHashProbabilisticAtlas::AddImage(RealImage image);
template void mirtkHashProbabilisticAtlas::AddImage(HashRealImage image);
template void mirtkHashProbabilisticAtlas::AddProbabilityMaps(int, RealImage **atlas);
template void mirtkHashProbabilisticAtlas::AddProbabilityMaps(int, HashRealImage **atlas);
template void mirtkHashProbabilisticAtlas::AddBackground(RealImage image);
template void mirtkHashProbabilisticAtlas::AddBackground(HashRealImage image);

}

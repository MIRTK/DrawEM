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


#include "mirtk/ImageHistogram1D.h"

namespace mirtk {

template <class VoxelType> void ImageHistogram1D<VoxelType>::Evaluate(GenericImage<VoxelType> *image, VoxelType padding)
{  
	double value,min,max;
	int i,j,k,l;
	image->GetMinMaxAsDouble(&min,&max);
	this->PutMin(min);
	this->PutMax(max);
	this->PutNumberOfBins(512);
	for (l = 0; l < image->GetT(); l++){
		for (k = 0; k < image->GetZ(); k++){
			for (j = 0; j < image->GetY(); j++){
				for (i = 0; i < image->GetX(); i++){
					value = image->GetAsDouble(i, j, k, l);
					if(value > padding)
						this->AddSample(value);
				}
			}
		}
	}
}

template <class VoxelType> void ImageHistogram1D<VoxelType>::BackProject(GenericImage<VoxelType> *image)
{  
	VoxelType value;
	int i,j,k,l;
	for (l = 0; l < image->GetT(); l++){
		for (k = 0; k < image->GetZ(); k++){
			for (j = 0; j < image->GetY(); j++){
				for (i = 0; i < image->GetX(); i++){
					value = image->GetAsDouble(i, j, k, l);
					value = this->ValToBin(value);
					value = this->_bins[(int)round(value)];
					image->PutAsDouble(i,j,k,l,value);
				}
			}
		}
	}
}

template <class VoxelType> void ImageHistogram1D<VoxelType>::Equalize(VoxelType min,VoxelType max)
{
	int i;
	double count = 0,current;
	for(i=0;i<this->_nbins;i++){
        current = this->BinToPDF(i);
		this->_bins[i] = (count+current/2.0)*(max - min) + min;
        count += current;
	}
    this->_emin = min;
    this->_emax = max;
}

template class ImageHistogram1D<unsigned char>;
template class ImageHistogram1D<short>;
template class ImageHistogram1D<unsigned short>;
template class ImageHistogram1D<float>;
template class ImageHistogram1D<double>;

}

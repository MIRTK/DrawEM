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


#ifndef _ImageHistogram1D_H

#define _ImageHistogram1D_H

#include "mirtk/Histogram1D.h"
#include "mirtk/GenericImage.h"

namespace mirtk {

template <class VoxelType> 
class ImageHistogram1D : public Histogram1D<VoxelType>
{
protected:
    /// min for equalize
    VoxelType _emin;
    /// max for equalize
    VoxelType _emax;

public:
	/// Evaluate the histogram from a given image with padding value
	virtual void Evaluate(GenericImage<VoxelType> *, VoxelType padding = -10000);
	/// Histogram Equalization
	virtual void Equalize(VoxelType min,VoxelType max);
	/// Back project the equalized histogram to image
	virtual void BackProject(GenericImage<VoxelType> *); 
};

}
#endif

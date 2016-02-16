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


#ifndef _MIRTKBIASCORRECTION_H

#define _MIRTKBIASCORRECTIO_H

#include <mirtkImage.h>

#include <mirtkResampling.h>

#include <mirtkTransformation.h>

#include <mirtkBiasField.h>

namespace mirtk{

class mirtkBiasCorrection : public Object
{
	mirtkObjectMacro(mirtkBiasCorrection);

protected:

	RealImage *_target;

	RealImage *_reference;

	RealImage *_weights;
	ByteImage *_mask;

	/// Output
	mirtkBiasField *_biasfield;

	/// Padding value
	GreyPixel _Padding;

	/// Initial set up for the registration
	virtual void Initialize();

	/// Final set up for the registration
	virtual void Finalize();

public:

	/// Constructor
	mirtkBiasCorrection();

	/// Destructor
	virtual ~mirtkBiasCorrection();

	/// Sets input for the bias correction filter
	virtual void SetInput (RealImage *, RealImage *);

	/// Sets weights for the bias correction filter
	virtual void SetWeights (RealImage *);

	/// Sets output for the bias correction filter
	virtual void SetOutput(mirtkBiasField *);

	virtual void SetMask( ByteImage *);

	/// Runs the bias correction filter
	virtual void Run();

	/// Apply bias correction to _input
	virtual void Apply(RealImage &);

	/// Apply bias correction to any image
	virtual void ApplyToImage(RealImage &);

	/// Apply bias correction to any image including logarithmic transform
	virtual void ApplyToImage(GreyImage &);

	// Access parameters
	virtual void SetPadding(short Padding);
	virtual short GetPadding();

};

inline void mirtkBiasCorrection::SetInput(RealImage *target, RealImage *reference)
{
	_target    = target;
	_reference = reference;
}

inline void mirtkBiasCorrection::SetWeights(RealImage *weights)
{
	_weights    = weights;
}

inline void mirtkBiasCorrection::SetOutput(mirtkBiasField *biasfield)
{
	_biasfield = biasfield;
}

inline void mirtkBiasCorrection::SetPadding(short Padding)
{
	_Padding = Padding;
}
inline short mirtkBiasCorrection::GetPadding(){
	return _Padding;
}

}
#endif

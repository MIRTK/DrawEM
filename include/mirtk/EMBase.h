/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Christian Ledig
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

#ifndef _MIRTKEMCLASSIFICATIONNEO_H
#define _MIRTKEMCLASSIFICATIONNEO_H

#include "mirtk/Image.h"
#include "mirtk/Object.h"
#include "mirtk/HashProbabilisticAtlas.h"
#include "mirtk/Gaussian.h"
#include "mirtk/Histogram1D.h"
#include "mirtk/MeanShift.h"

/*

Expectation Maximisation Algorithm

 */
namespace mirtk {

class EMBase : public Object
{

protected:
	typedef GenericImage<int> IntegerImage;

	/// Input image
	RealImage _input;

	/// Brain tissue
	RealImage _brain;

	/// weights image
    RealImage _weights;

	/// image estimate
	RealImage _estimate;

	/// Posterior probability maps -  segmentation
	HashProbabilisticAtlas _output;

	/// Partial volume soft segmentation
	HashProbabilisticAtlas _pv_output;

	/// Probability maps (atlas)
	HashProbabilisticAtlas _atlas;

	/// image segmentation
    IntegerImage _segmentation;

    /// posterior penalty
    RealImage _postpenalty;
    bool _postpen;

	/// Number of tissues
	int _number_of_tissues;

	/// Number of voxels
	int _number_of_voxels;

	/// Gaussian distribution parameters for each tissue type
    Gaussian *_G;
	double *_mi;
	double *_sigma;
	/// mixing coefficients for GMM
	double *_c;

	/// Likelihood
	double _f;

	/// Padding value
	RealPixel _padding;

    /// whether we have additional background tissue
	bool _has_background;

    /// superlabels
	int *_super;
	bool _superlabels;

    /// whether a mask is set
    bool _mask_set;

    /// whether initial posteriors is set
    bool _posteriors_set;

public:
	/// Input mask
	ByteImage _mask;

	/// Empty constructor
	EMBase();

	/// Constructor when adding background
	template <class ImageType>
	EMBase(int noTissues, ImageType **atlas, ImageType *background);

	/// Constructor without adding background
	template <class ImageType>
	EMBase(int noTissues, ImageType **atlas);

	/// Constructor without adding background + initialise posteriors
	template <class ImageType>
	EMBase(int noTissues, ImageType **atlas, ImageType **initposteriors ) ;

	/// Destructor
	virtual ~EMBase();

	// add a probability map
	template <class ImageType>
	void addProbabilityMap(ImageType image);

	// add background
	void addBackground();
	template <class ImageType>
	void addBackground(ImageType image);

	// normalise the atlas
	void NormalizeAtlas();

	/// Initialize filter
	virtual void Initialise();

	/// Initialize parameters
	void InitialiseParameters();

	/// Initialize filter
    virtual void InitialiseGMM();

	/// Initialize Gaussian Mixture Parameters and calculate initial posteriors
	void InitialiseGMMParameters(int n, double *m, double *s, double *c);

	/// Returns the name of the class
	virtual const char* NameOfClass() const;

	/// Estimates posterior probabilities
	virtual void EStep();

	/// Estimates parameters
	virtual void MStep();

	/// Estimates posterior probabilities
	virtual void EStepGMM(bool uniform_prior = false);

	/// Estimates parameters in GMM
	virtual void MStepGMM(bool uniform_prior = false);

	/// Estimates parameters in GMM with equal variance
	virtual void MStepVarGMM(bool uniform_prior = false);

	/// Computes log likelihood for current parameters
	virtual double LogLikelihood();

	/// Computes log likelihood for GMM for current parameters
	virtual double LogLikelihoodGMM();

	/// Get ProbMap
	virtual void GetProbMap(int i,RealImage& image);

	/// Print parameters
	virtual void Print();

	/// Print parameters
	virtual void PrintGMM();

	/// Calculate weights
	virtual void WStep();

	/// Set image
	virtual void SetInput(const RealImage &);

    /// set mask
    void SetMask(ByteImage &mask);

	/// Set padding value
	virtual void SetPadding(RealPixel);

    /// Set posterior penalty
    void setPostPenalty(RealImage &postpenalty);
    /// Set superlabels
    void setSuperlabels(int *superlabels);

	/// Execute one iteration and return log likelihood
	virtual double Iterate(int iteration);

	/// Execute one iteration and return log likelihood
	virtual double IterateGMM(int iteration, bool equal_var = false, bool uniform_prior = false);

    /// Get GMM memberships
    virtual void GetProportions(double *);

	/// Construct segmentation based on current posterior probabilities
	virtual void ConstructSegmentation(IntegerImage &);
	/// Construct segmentation based on current posterior probabilities
	virtual void ConstructSegmentation();

    /// Write probability map into a file
	void WriteProbMap(int i, const char *filename);
    /// Write Gaussian parameters into a file
	void WriteGaussianParameters(const char *file_name, int flag = 0);
    /// Write image estimate
	void WriteEstimate(const char *filename);
    /// Write weights
	void WriteWeights(const char *filename);
    /// Write input
	void WriteInput(const char *filename);
    /// Write segmentation
	void WriteSegmentation(const char *filename);


    /// Returns log likelihood for given intensity value
    double PointLogLikelihoodGMM(double x);
    /// Initialize gaussians with current GMM parameters
	void GInit();
    /// create mask taking into account priors and padding
    void CreateMask();
    /// return means
	virtual void GetMean(double *);
    /// return variances
    virtual void GetVariance(double *);

    /// initialise GMM parameters
	void InitialiseGMMParameters(int n);
    /// set uniform prior
    void UniformPrior();
};


inline void EMBase::addBackground(){
	_atlas.AddBackground();
	_has_background = true;
	_number_of_tissues = _atlas.GetNumberOfMaps();
}

template <class ImageType>
inline void EMBase::addBackground(ImageType background){
	_atlas.AddBackground(background);
	_has_background = true;
	_number_of_tissues = _atlas.GetNumberOfMaps();
}


template <class ImageType>
inline void EMBase::addProbabilityMap(ImageType image){
	_atlas.AddImage(image);
	_number_of_tissues = _atlas.GetNumberOfMaps();
}

inline void EMBase::NormalizeAtlas(){
	_atlas.NormalizeAtlas();
}

inline void EMBase::SetPadding(RealPixel padding)
{
	_padding = padding;
}

inline void EMBase::WriteEstimate(const char *filename)
{
	_estimate.Write(filename);
}

inline void EMBase::WriteInput(const char *filename)
{
	_input.Write(filename);
}

inline void EMBase::WriteSegmentation(const char *filename)
{
	_segmentation.Write(filename);
}

inline const char* EMBase::NameOfClass() const { return "EMBase"; } 

}

#endif

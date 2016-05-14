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

#ifndef DrawEM_H_
#define DrawEM_H_

#include "mirtk/EMBase.h"

#include "mirtk/PolynomialBiasField.h"
#include "mirtk/Image.h"
#include "mirtk/HashProbabilisticAtlas.h"
#include "mirtk/BiasField.h"
#include "mirtk/BiasCorrection.h"
#include "mirtk/Image.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/EuclideanDistanceTransform.h"
#include "mirtk/ConnectedComponents.h"

#include <set>
#include <map>
#include <vector>
#include <utility>

namespace mirtk {

class DrawEM : public EMBase
{
    mirtkObjectMacro(DrawEM);

protected:
    //

    /// local weights for MRF field
    RealImage _MRF_weights;

    /// Uncorrected image
    RealImage _uncorrected;

    /// Bias field correction filter
    BiasCorrection _biascorrection;

    /// Bias field
    BiasField *_biasfield;

    /// MRF connectivity
    Matrix _connectivity;

    /// PV classes
    map<int,int> pv_classes;
    vector< pair<int, int> > pv_connections;
    vector<double> pv_fc;

    /// use 26-nn MRF
    bool bignn;
    /// do PV correction similar to Hui et al.
    bool huipvcorr;
    /// weight of the MRF
    double mrfweight;
    /// intra MRF
    double beta,betainter;

    /// inter MRF
    bool intermrf;
    RealImage **_MRF_inter;

    /// the tissue class of each label
    int *tissuelabels;
    int csflabel,wmlabel,gmlabel,outlabel;

private:
    bool isPVclass(int pvclass);
    double getMRFenergy(int index, int tissue);
    double getMRFInterEnergy(int index, int tissue);

public:
    /// Constructor
    DrawEM();
    template <class ImageType>
    DrawEM(int noTissues, ImageType **atlas, ImageType *background);
    template <class ImageType>
    DrawEM(int noTissues, ImageType **atlas);
    template <class ImageType>
    DrawEM(int noTissues, ImageType **atlas, ImageType **initposteriors);

    /// Initialize parameters
    void InitialiseParameters();

    /// Get the bias field
    void GetBiasField(RealImage &image);

    /// add partial volume between classes classA and classB
    int AddPartialVolumeClass(int classA, int classB, int huiclass=0);

    /// estimate probabilities
    void EStepMRF(void);

    /// relax priors
    void RStep(void);
    void RStep(double rf);

    /// Estimates bias field
    virtual void BStep();

protected:

    using EMBase::SetInput;

public:

    /// Set image
    virtual void SetInput(const RealImage &, const Matrix &);

    /// Set bias field
    virtual void SetBiasField(BiasField *);

    /// Execute one iteration and return log likelihood
    virtual double Iterate(int iteration);

    /// Compute the bias corrected image
    virtual void GetBiasCorrectedImage(RealImage &);


    /// set the MRF strength
    virtual void setMRFstrength(double mrfw);
    /// set a 26-neighborhood in the MRF
    virtual void setbignn(bool bnn);
    /// computes the MRF with the 26-neighborhood
    double getMRFenergy_diag(int index, int tissue);

    /// removes the PV classes!
    //void removePVclasses(double threshold = 0.1);

    /// sets the tissue class of each label
    void setTissueLabels(int num,int *atisslabels);

    /// set hui-style PV correction
    void setHui(bool hui);
    /// hui-style PV correction - HACKY!!
    void huiPVCorrection(bool changePosterior=false);
    /// construct seg with the tissue class of each label instead of the label itself
    void ConstructSegmentationHui(IntegerImage &segmentation);
    /// finds the overall tissue probabilities at the voxel (x,y,z)
    /// by adding the probability of the different labels belonging to the tissue class
    void getHuiValues(double &outval,double &csfval,double &gmval,double &wmval,int x,int y,int z,bool atlas);
    /// modifies the probability of the labels according to overall tissue probabilities at the voxel (x,y,z)
    /// the overall tissue probability is divided to its labels according to their "contribution" to the tissue
    void setHuiValues(double &outval,double &csfval,double &gmval,double &wmval,int x,int y,int z,bool atlas);

    /// set beta_intra
    virtual void setBeta(double beta);
    /// set beta_inter
    virtual void setBetaInter(double betainter);
    /// set mrf_inter term
    virtual void setMRFInterAtlas(RealImage **&atlas);

};

inline void DrawEM::setHui(bool hui){huipvcorr=hui;}
inline void DrawEM::setbignn(bool bnn){bignn=bnn;}
inline void DrawEM::setMRFstrength(double mrfw){mrfweight=mrfw;}
inline void DrawEM::setMRFInterAtlas(RealImage **&atlas){	_MRF_inter=atlas; intermrf=true;}
inline void DrawEM::setBeta(double b){beta=b;}
inline void DrawEM::setBetaInter(double b){betainter=b;}
inline void DrawEM::setTissueLabels(int num,int *atisslabels){
    tissuelabels=new int[num];
    for(int i=0;i<num;i++)  tissuelabels[i]=atisslabels[i];
}

}




#endif /* DrawEM_H_ */

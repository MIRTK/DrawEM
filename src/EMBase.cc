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


#include "mirtk/EMBase.h"

namespace mirtk {

EMBase::EMBase(){
	InitialiseParameters();
}

template <class ImageType>
EMBase::EMBase(int noTissues, ImageType **atlas, ImageType *background)
{
	InitialiseParameters();
	_atlas.AddProbabilityMaps(noTissues, atlas);
	addBackground(*background);
	_number_of_tissues = _atlas.GetNumberOfMaps();
}

template <class ImageType>
EMBase::EMBase(int noTissues, ImageType **atlas)
{
	InitialiseParameters();
	_atlas.AddProbabilityMaps(noTissues, atlas);
	_number_of_tissues = _atlas.GetNumberOfMaps();
}

template <class ImageType>
EMBase::EMBase(int noTissues, ImageType **atlas, ImageType **initposteriors)
{
	InitialiseParameters();
	_atlas.AddProbabilityMaps(noTissues, atlas);
	_output.AddProbabilityMaps(noTissues, initposteriors);
	_posteriors_set=true;
	_number_of_tissues = _atlas.GetNumberOfMaps();
}


EMBase::~EMBase()
{
	if (_G!=NULL) delete []_G;
	delete []_mi;
	delete []_sigma;
	delete []_c;
}

void EMBase::InitialiseParameters(){
	_padding = MIN_GREY;
	_number_of_voxels = 0;
	_number_of_tissues = 0;
	_mi = NULL;
	_sigma = NULL;
	_c = NULL;
	_f = 0;
    _G = NULL;
	_superlabels=false;
    _postpen=false;
	_posteriors_set=false;
	_has_background=false;
	_mask_set=false;
}

void EMBase::SetInput(const RealImage &image)
{
	_input = image;
	_estimate = image;
    _weights = image;
    _number_of_voxels=_input.GetNumberOfVoxels();
}

void EMBase::CreateMask()
{
	_atlas.First();
    if(!_mask_set){
            _mask.Initialize(_input.Attributes());
    }
	RealPixel *p=_input.GetPointerToVoxels();
	BytePixel *m=_mask.GetPointerToVoxels();
	for (int i=0; i<_number_of_voxels; i++) {
        if(!_mask_set || *m==1){
            if (*p!=_padding){
                bool flag = true;
                for (int j = 0; j < _number_of_tissues; j++) {
                    if (_atlas.GetValue(j) > 0) {
                        flag=false;break;
                    }
                }
                if(flag) *m=0;
                else *m=1;
            }
            else *m=0;
        }

		p++;
		m++;
		_atlas.Next();
	}
    _mask_set = true;
}

void EMBase::SetMask(ByteImage &mask)
{
	_mask=mask;
	_mask_set = true;
}

void EMBase::Initialise()
{
	_number_of_tissues = _atlas.GetNumberOfMaps();

	_atlas.NormalizeAtlas();
	if(!_posteriors_set) _output = _atlas;
	else _output.NormalizeAtlas();

    CreateMask();

	_mi = new double[_number_of_tissues];
	_sigma = new double[_number_of_tissues];
	_c = new double[_number_of_tissues];

	this->MStep();
	Print();
}

void EMBase::InitialiseGMM()
{
	_atlas.NormalizeAtlas();
	if(!_posteriors_set) _output = _atlas;
	else _output.NormalizeAtlas();

	this->MStepGMM();
	Print();
}


void EMBase::setPostPenalty(RealImage &pp){
    _postpenalty=pp;
    _postpen=true;
}


void EMBase::UniformPrior()
{
	int i, k;

	_atlas.First();
	BytePixel *pm = _mask.GetPointerToVoxels();
    for (i=0; i< _number_of_voxels; i++) {
		if (*pm == 1) {
			for (k = 0; k < _number_of_tissues; k++) {
				_atlas.SetValue(k,1.0/_number_of_tissues);
			}
		}
		pm++;
		_atlas.Next();
	}
}

void EMBase::InitialiseGMMParameters(int n)
{
	int i;
    std::cout <<"Estimating GMM parameters ... ";
	RealPixel imin, imax;
    _input.GetMinMaxPad(&imin, &imax,_padding);
	_number_of_tissues=n;
	_mi=new double[n];
	_sigma=new double[n];
	_c=new double[n];
	for (i=0; i<n; i++) {
		_mi[i]=imin+i*(imax-imin)/(n-1);
		_sigma[i]=((imax-imin)/(n-1))*((imax-imin)/(n-1));
		_c[i]=1.0/n;
	}
	PrintGMM();
}


void EMBase::InitialiseGMMParameters(int n, double *m, double *s, double *c)
{
	int i;

	_number_of_tissues = n;
	_mi = new double[_number_of_tissues];
	_sigma = new double[_number_of_tissues];
	_c = new double[_number_of_tissues];

	for ( i=0; i<n; i++) {
		_mi[i]=m[i];
		_sigma[i]=s[i];
		_c[i]=c[i];
	}

    PrintGMM();
	_output = _atlas;
	EStepGMM();
}

void EMBase::MStep()
{
    std::cout << "M-step" << std::endl;
    int k;
    double mi_num[_number_of_tissues];
    double sigma_num[_number_of_tissues];
    double denom[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
        mi_num[k] = 0;
        sigma_num[k] = 0;
		denom[k] = 0;
	}

	//superlabels
    double mi_num_super[_number_of_tissues];
    double sigma_num_super[_number_of_tissues];
    double denom_super[_number_of_tissues];
    if(_superlabels){
        for (k = 0; k < _number_of_tissues; k++) {
            mi_num_super[k] = 0;
            sigma_num_super[k] = 0;
            denom_super[k]=0;
        }
    }

    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                mi_num[k] += it->second * _input.Get(it->first);
                denom[k]  += it->second;
            }
        }
    }


    /*
	_output.First();
	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
            for (k = 0; k < _number_of_tissues; k++) {
				mi_num[k] += _output.GetValue(k) * *ptr;
				denom[k]  += _output.GetValue(k);
			}
		}
		ptr++;
		pm++;

		_output.Next();
    }*/


	//superlabels
    if(_superlabels){
		for (k = 0; k < _number_of_tissues; k++) {
            mi_num_super[_super[k]] +=  mi_num[k];
			denom_super[_super[k]] += denom[k];
		}

		for (k = 0; k < _number_of_tissues; k++) {
			mi_num[k] =  mi_num_super[_super[k]];
			denom[k] = denom_super[_super[k]];
		}
	}

	for (k = 0; k < _number_of_tissues; k++) {
		if (denom[k] != 0) {
			_mi[k] = mi_num[k] / denom[k];
		} else {
            std::cerr << "Division by zero while computing tissue mean!" << std::endl;
			exit(1);
		}
	}


    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                sigma_num[k] += it->second * pow(_input.Get(it->first) - _mi[k],2);
            }
        }
    }


    /*
	_output.First();
	ptr = _input.GetPointerToVoxels();
	pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
			for (k = 0; k <_number_of_tissues; k++) {
				sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
			}
		}
		ptr++;
		pm++;

		_output.Next();
    }*/

	//superlabels
    if(_superlabels){
		for (k = 0; k < _number_of_tissues; k++) {
			sigma_num_super[_super[k]] +=  sigma_num[k];
		}

		for (k = 0; k < _number_of_tissues; k++) {
			sigma_num[k] =  sigma_num_super[_super[k]];
		}
	}

    for (k = 0; k <_number_of_tissues; k++) {
		_sigma[k] = sigma_num[k] / denom[k];
		_sigma[k] = max( _sigma[k], 0.005 );

	}
}

void EMBase::EStep()
{
    std::cout << "E-step" << std::endl;
	int i, k;
	double x;

	RealImage segmentation;

	Gaussian* G = new Gaussian[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
		G[k].Initialise( _mi[k], _sigma[k]);
	}

	RealPixel *pptr;
    if(_postpen) pptr= _postpenalty.GetPointerToVoxels();


	double likelihood[_number_of_tissues];double sumlike;

	_atlas.First();
	_output.First();
	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();
    double* numerator = new double[_number_of_tissues];
    double denominator, temp;
    for (i=0; i< _number_of_voxels; i++) {
		if (*pm == 1) {

			x = *ptr;
            sumlike=0;
			for (k = 0; k < _number_of_tissues; k++) {
				likelihood[k] = G[k].Evaluate(x);
                sumlike+=likelihood[k];
			}
			for (k = 0; k < _number_of_tissues; k++) {
				likelihood[k] /= sumlike;
			}

            denominator=0;
			for (k = 0; k < _number_of_tissues; k++) {
				temp = likelihood[k] * _atlas.GetValue(k);
				numerator[k] = temp;
				denominator += temp;
			}

			//model averaging
            if (_postpen && denominator != 0) {
				double olddenom=denominator;
				denominator=0;
				for (k = 0; k < _number_of_tissues; k++) {
					double value = numerator[k]/olddenom;
					double priorvalue=_atlas.GetValue(k);
					value=(1-*pptr)*value +*pptr *  priorvalue;
					numerator[k] = value;
					denominator += value;
				}
			}



			for (k = 0; k < _number_of_tissues; k++) {
				if (denominator != 0) {
					double value = numerator[k]/denominator;
					if ((value < 0) || (value > 1)) {
						int x,y,z;
						_input.IndexToVoxel(i, x, y, z);
						std::cerr << "Probability value = " << value <<" @ Estep at voxel "<< x<<" "<<y<<" "<<z<< ", structure " << k << std::endl;
						if (value < 0)value=0;
						if (value > 1)value=1;
					}_output.SetValue(k, value);
                		} else {
                   			 _output.SetValue(k,_atlas.GetValue(k));
				}
			}

            if (denominator == 0) {
                int x,y,z;
                _input.IndexToVoxel(i, x, y, z);
                std::cerr<<"Division by 0 while computing probabilities at voxel "<<x<<","<<y<<","<<z<<std::endl;
            }
		} else {
			for (k = 0; k < _number_of_tissues ; k++) {
				_output.SetValue(k, 0);
			}
		}
		ptr++;
		pm++;
        if(_postpen)pptr++;
		_atlas.Next();
        _output.Next();
	}
    delete[] numerator;
	delete[] G;
}


void EMBase::WStep()
{
    std::cout << "W-step" << std::endl;
	int i,k;
	double num, den;
    std::cout<<"Calculating weights ...";
	RealPixel *pi=_input.GetPointerToVoxels();
    RealPixel *pw=_weights.GetPointerToVoxels();
	RealPixel *pe=_estimate.GetPointerToVoxels();
	BytePixel *pm=_mask.GetPointerToVoxels();
	_output.First();
	_atlas.First();

    for (i=0; i< _number_of_voxels; i++) {
        if (*pm == 1){
			num=0;
			den=0;
			for (k=0; k<_number_of_tissues; k++) {
				num += _output.GetValue(k)*_mi[k]/_sigma[k];
				den += _output.GetValue(k)/_sigma[k];
			}
            if (den!=0){
                *pw=den;
                *pe=num/den;
            }else{
                *pw=_padding;
                *pe=_padding;
            }
		} else {
			*pw=_padding;
			*pe=_padding;
		}

		pi++;
		pm++;
        pw++;
		pe++;
		_output.Next();
		_atlas.Next();
	}
    std::cout<<"done."<<std::endl;
}

void EMBase::GetMean(double *mean){
	int i;
	for(i=0;i<_number_of_tissues;i++){
		mean[i] = _mi[i];
	}
}

void EMBase::GetVariance(double *variance){
	int i;
	for(i=0;i<_number_of_tissues;i++){
		variance[i] = sqrt(_sigma[i]);
	}
}

void EMBase::MStepGMM(bool uniform_prior)
{
    std::cout << "M-step GMM" << std::endl;
    int k;
    double mi_num[_number_of_tissues];
    double sigma_num[_number_of_tissues];
    double denom[_number_of_tissues];
    double num_vox[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
        mi_num[k] = 0;
        sigma_num[k] = 0;
        denom[k] = 0;
		num_vox[k] = 0;
	}


    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                mi_num[k] += it->second * _input.Get(it->first);
                denom[k]  += it->second;
                num_vox[k]+= 1;
            }
        }
    }

    /*
	_output.First();
	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		// Check for backgound
		if (*pm == 1) {
			for (k = 0; k < _number_of_tissues; k++) {
				mi_num[k] += _output.GetValue(k) * *ptr;
				denom[k]  += _output.GetValue(k);
				num_vox[k]+= 1;
			}
		}
		ptr++;
		pm++;
		_output.Next();
    }*/

	for (k = 0; k < _number_of_tissues; k++) {
		if (denom[k] != 0) {
			_mi[k] = mi_num[k] / denom[k];
		} else {
			std::cerr <<"Tissue "<< k <<": Division by zero while computing tissue mean!" << std::endl;
			exit(1);
		}
		if (uniform_prior) _c[k]=1.0/_number_of_tissues;
		else _c[k]=denom[k]/num_vox[k];
	}


    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                sigma_num[k] += it->second * pow(_input.Get(it->first) - _mi[k],2);
            }
        }
    }

    /*_output.First();
	ptr = _input.GetPointerToVoxels();
	pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
			for (k = 0; k <_number_of_tissues; k++) {
				sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
			}
		}
		ptr++;
		pm++;
		_output.Next();
    }*/

	for (k = 0; k <_number_of_tissues; k++) {
		_sigma[k] = sigma_num[k] / denom[k];
		if(_sigma[k]<1)
			_sigma[k] = 1;
	}
}

void EMBase::MStepVarGMM(bool uniform_prior)
{
    std::cout << "M-step VarGMM" << std::endl;
    int k;
    double mi_num[_number_of_tissues];
    double denom[_number_of_tissues];
    double num_vox[_number_of_tissues];
    double sigma_num = 0;

	for (k = 0; k < _number_of_tissues; k++) {
		mi_num[k] = 0;
        denom[k] = 0;
        num_vox[k] = 0;
    }

    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                mi_num[k] += it->second * _input.Get(it->first);
                denom[k]  += it->second;
                num_vox[k]+= 1;
            }
        }
    }

    /*
	_output.First();
	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
			for (k = 0; k < _number_of_tissues; k++) {
				mi_num[k] += _output.GetValue(k) * *ptr;
				denom[k]  += _output.GetValue(k);
				num_vox[k]+= 1;
			}
		}
		ptr++;
		pm++;
		_output.Next();
    }*/


	for (k = 0; k < _number_of_tissues; k++) {
		if (denom[k] != 0) {
			_mi[k] = mi_num[k] / denom[k];
		} else {
           		std::cerr << "Division by zero while computing tissue mean!" << std::endl;
			exit(1);
		}
		if (uniform_prior) _c[k]=1.0/_number_of_tissues;
		else _c[k]=denom[k]/num_vox[k];
	}



    for (k = 0; k < _number_of_tissues; k++) {
        HashRealImage::DataIterator begin=_output.Begin(k), end=_output.End(k);
        for ( auto it = begin; it != end; ++it ){
            if (_mask.Get(it->first)==1){
                sigma_num += it->second * pow(_input.Get(it->first) - _mi[k],2);
            }
        }
    }

    /*
	_output.First();
	ptr = _input.GetPointerToVoxels();
	pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
			for (k = 0; k <_number_of_tissues; k++) {
				sigma_num += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
			}
		}
		ptr++;
		pm++;
		_output.Next();
    }*/

	double sum =0;
	for (k = 0; k <_number_of_tissues; k++) sum += denom[k];
	for (k = 0; k <_number_of_tissues; k++) {
		if (sum>0) _sigma[k] = sigma_num / sum;
	}
}



void EMBase::EStepGMM(bool uniform_prior)
{
    std::cout << "E-step GMM" << std::endl;
	int i, k;
	double x;
	double *gv = new double[_number_of_tissues];
	double *numerator = new double[_number_of_tissues];
	Gaussian *G = new Gaussian[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
		G[k].Initialise( _mi[k], _sigma[k]);
	}

	_output.First();
	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();

    for (i=0; i< _number_of_voxels; i++) {
		double denominator=0, temp=0;
		if (*pm == 1) {
			x = *ptr;
			for (k = 0; k < _number_of_tissues; k++) {
				temp = G[k].Evaluate(x);
				gv[k] = temp;
				if (!uniform_prior) temp = temp * _c[k];
				numerator[k] = temp;
				denominator += temp;
			}
			for (k = 0; k < _number_of_tissues; k++) {
				if (denominator != 0) {
					double value = numerator[k]/denominator;

					if ((value < 0) || (value > 1)) {
						int x,y,z;
						_input.IndexToVoxel(i, x, y, z);
						std::cerr << "Probability value = " << value <<" @ Estep gmm at voxel "<< x<<" "<<y<<" "<<z<< ", structure " << k << std::endl;
						if (value < 0)value=0;
						if (value > 1)value=1;
					}_output.SetValue(k, value);
                } else {
                    _output.SetValue(k,_atlas.GetValue(k));
				}
            }

            if (denominator <= 0) {
                int x,y,z;
                _input.IndexToVoxel(i, x, y, z);
                 std::cerr<<"Division by 0 while computing probabilities at voxel "<<x<<","<<y<<","<<z<<std::endl;
            }

		} else {
			for (k = 0; k < _number_of_tissues ; k++) {
				_output.SetValue(k, 0);
            }
		}
		ptr++;
		pm++;
		_output.Next();
	}

	delete[] G;
	delete[] gv;
	delete[] numerator;
}

void EMBase::Print()
{
	int k;

    std::cout << "mean:";
	for (k = 0;  k <_number_of_tissues; k++) {
        std::cout << " " << k << ": " << _mi[k];
    }
    std::cout << std::endl;

    std::cout << "sigma:";
	for (k = 0; k < _number_of_tissues; k++) {
        std::cout << " " << k << ": " << sqrt(_sigma[k]);
	}
    std::cout << std::endl;
}

void EMBase::PrintGMM()
{
	int k;
	Print();
    std::cout << "c:";
	for (k = 0; k < _number_of_tissues; k++) {
        std::cout << " " << k << ": " << _c[k];
    }
    std::cout << std::endl;

}

double EMBase::Iterate(int)
{
	this->EStep();
	this->MStep();
	Print();
	return LogLikelihood();
}

double EMBase::IterateGMM(int iteration, bool equal_var, bool uniform_prior)
{
	if (iteration > 1) this->EStepGMM();
	if (equal_var) this->MStepVarGMM(uniform_prior);
	else this->MStepGMM(uniform_prior);
	PrintGMM();

	return LogLikelihoodGMM();
}

double EMBase::LogLikelihood()
{
	int i, k;
	double temp, f;
    std::cout<< "Log likelihood: ";
	Gaussian* G = new Gaussian[_number_of_tissues];
	double* gv = new double[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
		G[k].Initialise( _mi[k], _sigma[k]);
	}

	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();
	_output.First();
	f = 0;
    for (i = 0; i < _number_of_voxels; i++) {
        if (*pm == 1) {
			temp = 0;
			double max = 0;
			int max_k = 0;

			for (k=0; k < _number_of_tissues; k++) {
				// Estimation of gaussian probability of intensity (*ptr) for tissue k
				gv[k] = G[k].Evaluate(*ptr);
				if( gv[k] > 1 ) gv[k] = 1.0;

				if( max < gv[k] )
				{
					max_k = k;
					max = gv[k];
				}

				// Probability that current voxel is from tissue k
				temp += gv[k] * _output.GetValue(k);
			}

			if ((temp > 0) && (temp <= 1)) {
				f += log(temp);
			}
		}
		ptr++;
		pm++;
		_output.Next();
	}

	f = -f;
	double diff, rel_diff;
	diff = _f-f;

	if (_f == 0) rel_diff = 1;
	else rel_diff = diff/_f;

	_f=f;

    std::cout << "f= "<< f << " diff = " << diff << " rel_diff = " << rel_diff <<std::endl;
	delete[] G;
	delete[] gv;

	return rel_diff;
}

double EMBase::LogLikelihoodGMM()
{
	int i, k;
	double temp, f;
    std::cout<< "Log likelihood GMM: ";
	Gaussian* G = new Gaussian[_number_of_tissues];
	double* gv = new double[_number_of_tissues];

	for (k = 0; k < _number_of_tissues; k++) {
		G[k].Initialise( _mi[k], _sigma[k]);
	}

	RealPixel *ptr = _input.GetPointerToVoxels();
	BytePixel *pm = _mask.GetPointerToVoxels();
	_output.First();
	f = 0;
    for (i = 0; i < _number_of_voxels; i++) {
		if (*pm == 1) {
			temp = 0;
			for (k=0; k < _number_of_tissues; k++) {
				// Estimation of gaussian probability of intensity (*ptr) for tissue k
				gv[k] = G[k].Evaluate(*ptr);
				// Probability that current voxel is from tissue k
				temp += gv[k] * _c[k];
			}

			if ((temp > 0) && (temp <= 1)) {
				f += log(temp);
			}
		}
		ptr++;
		pm++;
		_output.Next();
	}

	f = -f;
	double diff, rel_diff;
	diff = _f-f;

	if (_f == 0) rel_diff = 1;
	else rel_diff = diff/_f;

	_f=f;

    std::cout << "f= "<< f << " diff = " << diff << " rel_diff = " << rel_diff <<std::endl;
	delete[] G;
	delete[] gv;

	return rel_diff;
}

void EMBase::ConstructSegmentation(IntegerImage &segmentation)
{
    int i, j, m;
    RealPixel max;

    std::cout<<"Constructing segmentation"<<std::endl;

    // Initialize pointers of probability maps
    _output.First();

    // Initialize segmentation to same size as input
    segmentation = IntegerImage(_input.Attributes());
    RealPixel *ptr = _input.GetPointerToVoxels();
    int *sptr = segmentation.GetPointerToVoxels();
    BytePixel *pm = _mask.GetPointerToVoxels();

    for (i = 0; i < _number_of_voxels; i++) {
        m = 0;
        if (*pm == 1) {
            max  = 0;
            for (j = 0; j < _number_of_tissues; j++) {
                if (_output.GetValue(j) > max) {
                    max  = _output.GetValue(j);
                    m = j+1;
                    if ( _has_background && (j+1) == _number_of_tissues) m=0;
                }
            }
            if(max==0){
                int x,y,z;
                int index=i;
                z = index / (_input.GetX() * _input.GetY());
                index -= z * (_input.GetX() * _input.GetY());
                x = index % _input.GetX();
                y = index / _input.GetX();
                std::cerr<<"voxel at "<<x<<","<<y<<","<<z<<" has 0 prob"<<std::endl;
            }
        }
        *sptr = m;
        sptr++;
        ptr++;
        pm++;
        _output.Next();
    }
}

void EMBase::ConstructSegmentation()
{
	ConstructSegmentation(_segmentation);
}


void EMBase::GetProbMap(int i,RealImage& image){
	if  (i < _number_of_tissues) {
		image= _output.GetImage(i);
	} else {
		std::cerr << "HashProbabilisticAtlas::Write: No such probability map" << std::endl;
		exit(1);
	}
}

void EMBase::WriteProbMap(int i, const char *filename)
{
	if  (i < _number_of_tissues) {
		HashRealImage image= _output.GetImage(i);
		image.Write(filename);
	} else {
		std::cerr << "HashProbabilisticAtlas::Write: No such probability map" << std::endl;
		exit(1);
	}
}

void EMBase::WriteGaussianParameters(const char *file_name, int flag)
{
    std::cout << "Writing GaussianDistributionParameters: " << file_name << std::endl;

	ofstream fileOut(file_name);

	if (!fileOut) {
		std::cerr << "Can't open file " << file_name << std::endl;
		exit(1);
	}

	int k,l,m;

	if(flag){
		// out put without names
		for (k=0; k<_number_of_tissues; k++) {
			fileOut << _mi[k] << " " << _sigma[k] << std::endl;
		}
	}else{
		// out put with names
		fileOut << "mi: " <<std::endl;
		for (k=0; k<_number_of_tissues; k++) {
			fileOut << "Tissue " << k << ": (";
			for (l=0; l < 1/*_input.GetNumberOfChannels()*/; l++) {
				fileOut << _mi[k];//.Get(l);
				if (l == 0/*_input.GetNumberOfChannels() - 1*/) fileOut << ")" << std::endl;
				else fileOut << ", ";
			}
		}

		fileOut << "sigma: " << std::endl;
		for (k=0; k<_number_of_tissues; k++) {
			fileOut << "Tissue " << k << ": (";
			//<<std::endl << "(";

			for (l=0; l < 1/*_input.GetNumberOfChannels()*/; l++) {
				//fileOut << "(";
				for (m = 0; m < 1/*_input.GetNumberOfChannels()*/; m++) {
					double s = _sigma[k];//.Get(m,l);
					if ( s >= 0) fileOut << sqrt(s);
                    else fileOut << -sqrt(-s);
				}
			}
			//if (l == 0/*_input.GetNumberOfChannels() - 1*/) fileOut << ")" << std::endl;
			fileOut <<  ")" << std::endl;
		}
	}
}

void EMBase::WriteWeights(const char *filename)
{
	RealImage w(_weights);
	RealPixel *pw = w.GetPointerToVoxels();
	int i;
	double sigma_min=_sigma[0];

	for (i=1; i<_number_of_tissues; i++) if (sigma_min > _sigma[i]) sigma_min = _sigma[i];
	std::cerr<<"sigma min = "<<sigma_min<<std::endl;

	for (i=0; i<w.GetNumberOfVoxels(); i++) {
		if (*pw != _padding) *pw=(*pw) * sigma_min * 100;
		pw++;
	}

	w.Write(filename);
}

double EMBase::PointLogLikelihoodGMM(double x)
{
	int k;
	double temp=0;

	for (k=0; k < _number_of_tissues; k++) temp+= _G[k].Evaluate(x) * _c[k];
	if (-log(temp)> 1000000) exit(1);

	if ((temp > 1) || (temp < 0)) {
		std::cerr << "Could not compute likelihood, probability out of range = " << temp << std::endl;
		exit(1);
	}
	return -log(temp);
}


void EMBase::GInit()
{
	if (_G!=NULL) delete[] _G;
	_G = new Gaussian[_number_of_tissues];

	for (int k = 0; k < _number_of_tissues; k++) {
		_G[k].Initialise( _mi[k], _sigma[k]);
	}
}



void EMBase::setSuperlabels(int *superlabels){
	_super=superlabels;
	_superlabels=true;
}


void EMBase::GetProportions(double *proportions){
	int i;
	for(i=0;i<_number_of_tissues;i++){
		proportions[i] =  _c[i];
	}
}




template EMBase::EMBase(int, RealImage **, RealImage *);
template EMBase::EMBase(int, HashRealImage **, HashRealImage *);
template EMBase::EMBase(int, RealImage **);
template EMBase::EMBase(int, HashRealImage **);
template EMBase::EMBase(int, RealImage **, RealImage **);
template EMBase::EMBase(int, HashRealImage **, HashRealImage **);

template void EMBase::addProbabilityMap(RealImage image);
template void EMBase::addProbabilityMap(HashRealImage image);
template void EMBase::addBackground(RealImage image);
template void EMBase::addBackground(HashRealImage image);

}

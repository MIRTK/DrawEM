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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

#include "mirtk/PolynomialBiasField.h"
#include "mirtk/DrawEM.h"
#include "mirtk/Matrix.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input> <N> <prob1> .. <probN> <output> [options]" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
    std::cout << "  Runs the DrawEM segmentation at the input image with the provided N probability maps of structures. " << std::endl;
	std::cout << "  The main algorithm is outlined in [1]. " << std::endl;
	std::cout << std::endl;

	std::cout << "Input options:" << std::endl;
	std::cout << "GENERAL EM PARAMETERS:" << std::endl;
	std::cout << " -biasfielddegree <number>       polynomial degree (of one dimension) of biasfield (default = 4)" << std::endl;
	std::cout << " -mask <mask>                    mask image" << std::endl;
    std::cout << " -padding <number>               padding value (default is min intensity)" << std::endl;
    std::cout << " -iterations <number>            max number of iterations (default: 20)" << std::endl;
	std::cout << " -reldiff <double>               min relative difference that assumes convergence" << std::endl;
    std::cout << " -pv <class1> <class2>           add partial volume class between class class1 and class2" << std::endl;
	std::cout << std::endl;

	std::cout << "MRF PARAMETERS:" << std::endl;
	std::cout << " -mrf <file>                     use MRF, file contains a nxn connection matrix defining in at entry (i,j) if class i is adjacent, distant or identic to class j. " << std::endl;
	std::cout << " -mrfstrength <double>           MRF strength (default 1) " << std::endl;
    std::cout << " -bigmrf                         use 26-connectivity in the MRF neighborhood (default no)" << std::endl;
    std::cout << " -mrftimes <number>              max number of times mrf will be performed (default == number of max iterations)" << std::endl;
	std::cout << std::endl;

	std::cout << "RELAXATION PARAMETERS:" << std::endl;
    std::cout << " -relax                          perform prior relaxation (default:no) " << std::endl;
    std::cout << " -relaxfactor <double>           relaxation factor (weighting of smooth posteriors) (default: 0.5)" << std::endl;
    std::cout << " -relaxtimes <number>            number of times relaxation will be performed (default: 1)" << std::endl;
	std::cout << std::endl;

	std::cout << "SEGMENTATION PARAMETERS SPECIFIC IN [1]:" << std::endl;
	std::cout << " -postpenalty <file>              image with the GMM weighting. Note that the weight is inverse to the AM IEEE TMI paper (0:use posterior, 1:use prior)" << std::endl;
	std::cout << " -superlabel <N> <class1> .. <classN>    "
			"                                  classes defining a superclass" << std::endl;
	std::cout << " -tissues <Nout> <out1> .. <outN> <Ncsf> <csf1> .. <csfN> <Ngm> <gm1> .. <gmN> <Nwm> <wm1> .. <wmN>   "
			"                                  tissue classes" << std::endl;
	std::cout << " -hui                             hui-style PV correction (must have set tissue classes)" << std::endl;
	std::cout << " -pvh <class1> <class2> <tissue>  add partial volume class between class number1 and number2, set tissue class for the pv" << std::endl;
	std::cout << "                                  tissue is 1 for outlier, 2 for csf, 3 for gm, 4 for wm" << std::endl;
    std::cout << std::endl;

	std::cout << "OUTPUT PARAMETERS:" << std::endl;
	std::cout << " -corrected <file>                save corrected image" << std::endl;
	std::cout << " -savepv <file>                   save segmentation with pvs" << std::endl;
	std::cout << " -saveprob <number> <file>        save posterior probability of tissue to file"<<std::endl;
	std::cout << " -biasfield <file>                save final bias field (for log transformed intensities) to file. " << std::endl;
    std::cout << std::endl;
    PrintStandardOptions(std::cout);
    std::cout << std::endl;

	std::cout << "References:" << std::endl;
    std::cout << "[1] A. Makropoulos et al. Automatic whole brain MRI segmentation of the developing neonatal brain, IEEE TMI, 2014 "<<std::endl;

}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

void readConnectivityMatrix( Matrix* mat, char* filename)
{
	// Open file
	ifstream from;
	from.open(filename);

	int value;
    std::cout<<"Connectivity matrix: "<<filename<<std::endl;
	for( int i = 0; i < mat->Rows(); ++i ){
		for( int j = 0; j < mat->Cols(); ++j ){
			from >> value;
			mat->Put(i,j,value);
		}
    }
    std::cout<<std::endl;
    from.close();
}



int main(int argc, char **argv)
{
	REQUIRES_POSARGS(4);
	InitializeIOLibrary();
	int a=1;

	clock_t begin = clock();
	char *output_segmentation, *output_biascorrection, *output_biasfield, *connections, *output_pvsegmentation, *mask;
	int no_mrf_correction = 1;
        bool relax = false;
	vector< pair<int,int> > pv_classes;
	vector<int> hpv;
	int i, n, padding, maxIterations;
	int biasfield_degree = 4;
	output_biasfield = NULL;
	connections = NULL;
	output_pvsegmentation = NULL;
	mask = NULL;
	output_biascorrection=NULL;

	// Input image
	std::cout<<"reading "<<argv[1]<<std::endl;
	RealImage image;
	image.Read(POSARG(a++));
	image.Print();

	// Number of tissues
	n = atoi(POSARG(a++));
	std::cout<<n<<" atlases"<<std::endl;

	// Probabilistic atlas
	char **atlas_names=new char*[n];
	for (i = 0; i < n; i++) {
		atlas_names[i] = POSARG(a++);
        }

	// File name for segmentation
	output_segmentation = POSARG(a);

	// Default parameters
	maxIterations = 20;
	padding    = MIN_GREY;
	int number_of_classes = n;
	double rfactor=0.5;
	int relaxtimes=1;
	int mrftimes=maxIterations;
	bool mrfdecl=false;
	double mrfstrength=1;
    double reldiff=0.005;
	bool bignn=false,hui=false,superlbls=false;
	bool postpen=false,settissues=false;
	RealImage postpenalty;
	int *tissuelabels, *superlabels;
	int ss=0;
	vector<string> savesegs;
	vector<int> savesegsnr;
    char *output_pv=NULL;



	for (ALL_OPTIONS) {
		if (OPTION("-mrf")){
			connections = ARGUMENT;
		}
		else if (OPTION("-iterations")){
			maxIterations=atoi(ARGUMENT);
			if(!mrfdecl) mrftimes=maxIterations;
		}
		else if (OPTION("-reldiff")){
			reldiff = atof(ARGUMENT);
        }
		else if (OPTION("-corrected")){
			output_biascorrection = ARGUMENT;
		}
		else if (OPTION("-savepv")){
			output_pv = ARGUMENT;
		}
		else if (OPTION("-mask")){
			mask=ARGUMENT;
		}
		else if (OPTION("-pv")){
			int a, b;
			a = atoi(ARGUMENT);
			b = atoi(ARGUMENT);
			pv_classes.push_back(make_pair(a,b));
			hpv.push_back(0);
		}
		else if (OPTION("-pvh")){
			int a, b, c;
			a = atoi(ARGUMENT);
			b = atoi(ARGUMENT);
			c = atoi(ARGUMENT);
			pv_classes.push_back(make_pair(a,b));
			hpv.push_back(c);
        }
		else if (OPTION("-padding")){
			padding=atoi(ARGUMENT);
		}
		else if (OPTION("-biasfielddegree")){
			biasfield_degree=atoi(ARGUMENT);
			std::cout << "Degree of biasfield polynomial: " << biasfield_degree << std::endl;
		}
		else if (OPTION("-biasfield")){
			output_biasfield=ARGUMENT;
			std::cout << "Output biasfield to: " << output_biasfield << std::endl;
		}
		else if (OPTION("-postpenalty")){
			char *ppf=ARGUMENT;
			postpenalty.Read(ppf);
			std::cout<<"will use postpenalty "<<ppf<<std::endl;
			postpen=true;
		}
        else if (OPTION("-relax")){
            relax = true;
        }
		else if (OPTION("-relaxfactor")){
			rfactor=atof(ARGUMENT);
            relax = true;
		}
        else if (OPTION("-relaxtimes")){
            relaxtimes=atoi(ARGUMENT);
            relax = true;
        }
		else if (OPTION("-saveprobs")) {
			char* probsBase = ARGUMENT;
			savesegsnr.clear();
			savesegs.clear();
			for (i = 0; i < n; i++){
				ostringstream sstr;
				sstr<<probsBase<<i<<".nii.gz";
				savesegsnr.push_back(i);
				savesegs.push_back(sstr.str());
			}
			ss=n;
		}
		else if (OPTION("-saveprob")) {
			savesegsnr.push_back(atoi(ARGUMENT));
			savesegs.push_back(ARGUMENT);
			ss++;
		}
		else if (OPTION("-mrfstrength")){
			mrfstrength=atof(ARGUMENT);
		}
        else if (OPTION("-bigmrf")){
			bignn=true;
		}
		else if (OPTION("-mrftimes")){
			mrftimes=atoi(ARGUMENT);
			mrfdecl=true;
        }
		else if (OPTION("-hui")){
			hui=true;
		}
		else if (OPTION("-tissues")){
			tissuelabels=new int[n];
			for(int i=0;i<n;i++)tissuelabels[i]=0;

			for(int j=1;j<=4;j++){
				int numtiss=atoi(ARGUMENT);
				for(int i=0;i<numtiss;i++){
					int settiss=atoi(ARGUMENT);
                    tissuelabels[settiss]=j;
				}
			}
			settissues=true;
		}
		else if (OPTION("-superlabel")){
			int a = atoi(ARGUMENT);
			if(!superlbls){
				superlbls=true;
				superlabels=new int[number_of_classes];
				for(int i=0;i<number_of_classes;i++)superlabels[i]=i;
			}
			int superlbl=atoi(ARGUMENT);
            superlabels[superlbl]=superlbl;
            for(int i=1;i<a;i++){
				superlabels[atoi(ARGUMENT)]=superlbl;
			}
		}
		else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
	}



	// logtransform image
	// be careful log transformation might transform intensities to the actual padding! --> change padding value
	int min_log_voxel = 1000;
	RealImage padding_map = image;

    for( int i = 0; i < image.GetX(); ++i){
        for( int j = 0; j < image.GetY(); ++j){
            for( int k = 0; k < image.GetZ(); ++k){
                if(image.Get(i,j,k) == padding ){
					padding_map.Put(i,j,k, 1);
                }else{
					padding_map.Put(i,j,k, 0);
				}
                // no log transformation possible for 0 intensity values, treat like padding voxels
                if( image.Get(i,j,k) == 0 ){
					image.Put(i,j,k, padding);
					padding_map.Put(i,j,k,1);
				}
                if( image.Get(i,j,k) > 0 && image.Get(i,j,k) != padding ){
					image.Put(i,j,k, log(image.Get(i,j,k)));
					min_log_voxel = min(min_log_voxel,  static_cast<int>( image.Get(i,j,k)) );
				}
			}
		}
	}
	int original_padding = padding;
	padding = min_log_voxel-1;

    for( int i = 0; i < image.GetX(); ++i){
        for( int j = 0; j < image.GetY(); ++j){
            for( int k = 0; k < image.GetZ(); ++k){
                if( padding_map.Get(i,j,k) == 1 ){
					image.Put(i,j,k, padding);
				}
			}
        }
    }


    // Create classification object
    std::cout<<"initialize segmentation"<<std::endl;
    DrawEM *classification = new DrawEM();
	double atlasmin, atlasmax;
	for (i = 0; i < n; i++) {
		std::cout << "Image " << i <<" = " << atlas_names[i];
		RealImage atlas(atlas_names[i]);
		classification->addProbabilityMap(atlas);
		atlas.GetMinMaxAsDouble(&atlasmin, &atlasmax);
		std::cout << " with range: "<<  atlasmin <<" - "<<atlasmax<<std::endl;
	}

	Matrix* G;
	if( connections != NULL ){
		G = new Matrix(n,n);
		readConnectivityMatrix(G, connections);
		no_mrf_correction = 0;
	}
	else{
		// no MRF correction if theres no MRF
		no_mrf_correction = 1;
		G = new Matrix(1,1);
	}
	classification->SetInput(image, *G);

	if(bignn)classification->setbignn(bignn);
	if(superlbls)classification->setSuperlabels(superlabels);
	if(settissues)	classification->setTissueLabels(n,tissuelabels);
    else if(hui){ std::cerr<<"need to set tissues for pv correction"<<std::endl; PrintHelp(EXECNAME); exit(1);}
	if(hui)	classification->setHui(hui);
	if(mrfstrength!=1)classification->setMRFstrength(mrfstrength);

    classification->SetPadding(padding);
	if ( mask != NULL ){
		RealImage maskImage(mask);
		ByteImage maskByteImage(maskImage.Attributes());
		RealPixel* ptr = maskImage.GetPointerToVoxels();
		BytePixel* bptr = maskByteImage.GetPointerToVoxels();
		for( int i = 0; i < maskImage.GetNumberOfVoxels(); ++i ){
			if( *ptr > 0 ) *bptr = 1;
			else *bptr = 0;
			ptr++; bptr++;
		}
		classification->SetMask(maskByteImage);
    }
    classification->Initialise();
    std::cout<<"initialization: success"<<std::endl;
	std::cout.setf(std::ios::unitbuf);


	double rel_diff = 1.0;
	bool stop = false;
	int curr_biasfield_degree = 1;
	int iter = 0;
    int improvePhase = 0;

	// Create bias field
	PolynomialBiasField *biasfield = new PolynomialBiasField(image, curr_biasfield_degree);
	classification->SetBiasField(biasfield);

	bool BFupdate = false;
	bool MRFupdate = false;
	vector<int> pv_positions;
	int number_current_iterations = 0;
	bool PVon = false;
	bool relaxed = false;
	int modlabel=1;


    while( !stop && iter < maxIterations ){

		if( !MRFupdate || mrftimes<=0)
			classification->EStep();
		else{
			classification->EStepMRF();	
			mrftimes--;
		}


		if(hui && iter%2==modlabel)
			classification->huiPVCorrection();

		if( BFupdate ) {
			classification->WStep();
            classification->BStep();
		}
		classification->MStep();

		classification->Print();
		rel_diff = classification->LogLikelihood();


        if( rel_diff < reldiff && iter < maxIterations ){
			//once we have small rel_diff
            switch( improvePhase ){
			case 0:
				//gradually increase the biasfield_degree until desired
                if( curr_biasfield_degree < biasfield_degree && number_current_iterations){
					curr_biasfield_degree++;
					delete biasfield;
					biasfield = new PolynomialBiasField(image, curr_biasfield_degree);
					classification->SetBiasField(biasfield);

					BFupdate = true;
                    if( curr_biasfield_degree == biasfield_degree ){
						improvePhase++;
					}
					break;
                }else{
					improvePhase++;
				}
			case 1:
				improvePhase++;
				if(postpen){
					classification->setPostPenalty(postpenalty);
					break;
				}
			case 2:
				//MRF
                if( !MRFupdate && !no_mrf_correction ){
					MRFupdate = true;
					improvePhase++;
					break;
                }else{
					improvePhase++;
				}
			case 3:
                if(( relax && !relaxed )){
					// relax priors
					std::cout << "relaxing priors now..." << std::endl;
					classification->RStep(rfactor);
                    std::cout << "done!" << std::endl;
					relaxtimes--;
					if(relaxtimes==0){
						improvePhase++;
						relaxed = true;
					}
					break;
                }else{
					improvePhase++;
				}
			case 4:
				//add pv class
                if( pv_classes.size() && !PVon ){
                    for( unsigned int i = 0; i < pv_classes.size(); ++i ){
						int a, b;
						pair<int, int> tmp = pv_classes[i];
						a = tmp.first;
						b = tmp.second;

						int h = hpv[i];

						// add partial volume classes
						std::cout << "adding partial volume classes between " << a << " " << b<< " with hui class"<< h<< std::endl;
						pv_positions.push_back( classification->AddPartialVolumeClass( a, b, h) );
						std::cout << "New PV Class at position: " << pv_positions[i] << std::endl;
						std::cout << "done!" << std::endl;
					}
					PVon = true;
					improvePhase++;
					break;
                }else{
					improvePhase++;
				}
            case 5:
				stop = true;
				break;
			}
			number_current_iterations = 0;

        }else{
			number_current_iterations++;
		}

		iter++;
	}

	if(hui )classification->huiPVCorrection(true);



	GenericImage<int> output_image = image;
    classification->ConstructSegmentation(output_image);

	// Save segmentation
	std::cout<<"saving segmentation to "<<output_segmentation<<std::endl;
	output_image.Write(output_segmentation);

	RealImage bias(image.Attributes());
	if (output_biascorrection != NULL) {
		// Bias corrected image
		std::cout<<"preparing bias corrected image"<<std::endl;
		classification->GetBiasCorrectedImage(bias);

        for( int k = 0; k < bias.GetZ(); ++k){
            for( int j = 0; j < bias.GetY(); ++j){
                for( int i = 0; i < bias.GetX(); ++i){
                    if( bias.Get(i,j,k) == padding ){
						bias.Put(i,j,k, original_padding);
                    }else{
						bias.Put(i,j,k, exp(bias.Get(i,j,k) ));
					}
				}
			}
		}
		std::cout<<"saving bias corrected image to "<<output_biascorrection<<std::endl;
		bias.Write(output_biascorrection);
	}

	if (output_biasfield != NULL) {
		classification->GetBiasField( bias );
		std::cout<<"saving bias field to "<<output_biasfield<<std::endl;
		bias.Write(output_biasfield);
	}

	for( unsigned int i = 0; i < ss; ++i ){
		std::cout<<"saving probability map of structure "<<savesegsnr[i]<<" to "<<savesegs[i]<<std::endl;
		classification->WriteProbMap(savesegsnr[i],savesegs[i].c_str());
	}

	delete G;
	delete classification;


	clock_t end = clock();
	int elapsed_secs = round( double(end - begin) / CLOCKS_PER_SEC);
	int elapsed_mins = elapsed_secs/60;
	elapsed_secs%=elapsed_mins;
	std::cout<<"elapsed time: "<<elapsed_mins<<" min "<<elapsed_secs<<" sec"<<std::endl;

	return 0;
}


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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/GaussianBlurring.h"

#include "mirtk/Gaussian.h"
#include "mirtk/KMeans.h"


using namespace mirtk;
using namespace std;



// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input> <K> <output> [options]" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
	std::cout << "  Runs k-means clustering at the input image with the provided K number of classes. " << std::endl;
	std::cout << "  e.g. " << name << " input.nii.gz 10 output.nii.gz" << std::endl;
	std::cout << std::endl;
	std::cout << "Input options:" << std::endl;
	std::cout << "  -mask <mask>               run k-means inside the provided mask" << std::endl;
	std::cout <<	"  -replicates <number>       number of random initializations, the best is retained (default: 10)." << std::endl;
	std::cout <<	"  -iterations <number>       maximum number of iterations per replicate  (default: 100)." << std::endl;
	std::cout << "  -saveprobs <basename>      output probabilities. These are computed by initialising a Gaussian distribution for each cluster with the mean and variance of the points."<< std::endl;
	std::cout <<	"  -blurprobs <sigma>         The output probabilities are blurred with a kernel of size sigma." << std::endl;
	std::cout << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

int main(int argc, char **argv){
	REQUIRES_POSARGS(3);
	InitializeIOLibrary();
	int a=1;

	RealImage input;
	input.Read(POSARG(a++));
	int k = atoi(POSARG(a++));
	const char* outputname = POSARG(a);

	bool usemask = false;
	BinaryImage mask;
	const char* probsBase = NULL;
	double sigma = 0.0;
	int iterations=100;
	int replicates=10;

	for (ALL_OPTIONS) {
		if (OPTION("-mask")){
			usemask = true;
			mask.Read(ARGUMENT);
		}
		else if (OPTION("-saveprobs")) {
			probsBase = ARGUMENT;
		}
		else if (OPTION("-blurprobs")) {
			sigma = atof(ARGUMENT);
		}
		else if (OPTION("-replicates")) {
			replicates = atoi(ARGUMENT);
		}
		else if (OPTION("-iterations")) {
			iterations = atoi(ARGUMENT);
		}
		else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
	}


	if (!usemask){
		mask.Initialize (input.GetImageAttributes());
	}



	int numintensities=0;
	for(int x = 0; x < input.GetX(); x++){
		for(int y = 0; y < input.GetY(); y++){
			for(int z = 0; z < input.GetZ(); z++){
				if(!usemask){
					if(input.Get(x, y, z, 0)>0){
						mask.Put(x, y, z, 1);
						numintensities++;
					}
				}else{
					if(mask.Get(x, y, z, 0))
						numintensities++;
				}
			}
		}
	}

	int i=0;
	double * intensities=new double[numintensities];
	double sumintensities=0;
	for(int x = 0; x < input.GetX(); x++){
		for(int y = 0; y < input.GetY(); y++){
			for(int z = 0; z < input.GetZ(); z++){
				if(mask.Get(x, y, z, 0)){
					double val=input.Get(x, y, z, 0);
					intensities[i++]=val;
					sumintensities+=val;
				}
			}
		}
	}

	kmeans km(intensities,numintensities,k,iterations,replicates);
	int *intensitiesCentroids=km.getPointClusters();


	GreyImage output;
	output.Initialize (input.GetImageAttributes());
	int c=0;
	for(int x = 0; x < input.GetX(); x++){
		for(int y = 0; y < input.GetY(); y++){
			for(int z = 0; z < input.GetZ(); z++){
				if(mask.Get(x, y, z, 0)){
					output.Put(x, y, z, 0, intensitiesCentroids[c++] + 1);
				}
			}
		}
	}
	output.Write(outputname);
	output.Clear();




	if (probsBase != NULL){
		double *centroids=km.getCentroids();
		double *vars=new double[k];
		double *denom=new double[k];
		Gaussian *G=new Gaussian[k];


		for (i=0;i<k;i++){ 
			denom[i]=0;
			vars[i]=0;
		}
		for(int x=0;x<numintensities;x++){
			i=intensitiesCentroids[x];
			vars[i] += (intensities[x] - centroids[i]) * (intensities[x] - centroids[i]);
			denom[i]++;
		}
		for (i=0;i<k;i++){ 
			vars[i]/=denom[i];
		}


		RealImage *probs = new RealImage[k];
		for (int i = 0; i < k; i++){
			probs[i].Initialize (input.GetImageAttributes());
			G[i].Initialise (centroids[i], vars[i]);
		}

		double val, probval[k];
		for(int x = 0; x < input.GetX(); x++){
			for(int y = 0; y < input.GetY(); y++){
				for(int z = 0; z < input.GetZ(); z++){
					if(mask.Get(x, y, z, 0)){
						double sumprob=0;
						for (int i = 0; i < k; i++){
							val=input.Get(x, y, z, 0);
							probval[i]=G[i].Evaluate(input.Get(x, y, z, 0));
							sumprob+=probval[i];
						}
						if(sumprob==0)continue;
						for (int i = 0; i < k; i++){
							probval[i]/=sumprob;
							probs[i].Put(x,y,z,0,probval[i]);
						}

					}

				}
			}
		}

		for (int i=0;i<k;i++){
			if(sigma>0){
				GaussianBlurring<RealPixel> filter(sigma);
				filter.Input(&probs[i]);
				filter.Output(&probs[i]);
				filter.Run();
			}

			ostringstream ss;
			ss<<probsBase<<i<<".nii.gz";
			probs[i].Write((ss.str()).c_str());
			std::cout<<"c["<<i<<"]="<<centroids[i]<<" , "<<vars[i]<<std::endl;
		}

		delete[] probs;
		delete[] vars;
		delete[] G;
		delete[] intensities;
		delete[] centroids;
	}

	exit(0);
}




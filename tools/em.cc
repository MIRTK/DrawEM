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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

#include "mirtk/EMBase.h"
#include <vector>
#include <sstream>

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
	std::cout << "  Runs EM segmentation at the input image with the provided N probability maps of structures. " << std::endl;
	std::cout << "  e.g. " << name << " input.nii.gz 5 bg.nii.gz csf.nii.gz gm.nii.gz wm.nii.gz dgm.nii.gz segmentation.nii.gz" << std::endl;
	std::cout << std::endl;
	std::cout << "Input options:" << std::endl;
	std::cout << "  -mask <mask>               run EM inside the provided mask" << std::endl;
	std::cout << "  -padding <number>          run EM where input > padding value" << std::endl;
	std::cout << "  -iterations <number>       maximum number of iterations (default: 50)" << std::endl;
	std::cout << "  -saveprob <number> <file>  save posterior probability of structure with number <number> to file "<<std::endl;
	std::cout <<	"                             (0-indexed i.e. structure 1 has number 0)"<<std::endl;
	std::cout << "  -saveprobs <basename>      save posterior probability of structures to files with basename <basename>"<<std::endl;
	std::cout << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------




int main(int argc, char **argv)
{
	REQUIRES_POSARGS(4);
	InitializeIOLibrary();
	int a=1;

	int i, n, padding, iterations;

	// Input image
	RealImage image;
	image.Read(POSARG(a++));
	char *output_name;

	// Number of tissues
	n = atoi(POSARG(a++));

	// Probabilistic atlas
	char **atlas_names=new char*[n];
	for (i = 0; i < n; i++) {
		atlas_names[i] = POSARG(a++);
	}
	// Probabilistic atlas
	// File name for output
	output_name = POSARG(a);

	// Default parameters
	iterations = 50;
	padding    = -1;//MIN_GREY;
	bool usemask = false;
	ByteImage mask;

	int ss=0;
	vector<string> savesegs;
	vector<int> savesegsnr;

	for (ALL_OPTIONS) {
		if (OPTION("-mask")){
			usemask = true;
			mask.Read(ARGUMENT);
		}
		else if (OPTION("-padding")){
			padding=atoi(ARGUMENT);
		}
		else if (OPTION("-iterations")){
			iterations=atoi(ARGUMENT);
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
		else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
	}

	EMBase *classification = new EMBase();
	double atlasmin, atlasmax;
	for (i = 0; i < n; i++) {
		std::cout << "Image " << i <<" = " << atlas_names[i];
		RealImage atlas(atlas_names[i]);
		classification->addProbabilityMap(atlas);
		atlas.GetMinMaxAsDouble(&atlasmin, &atlasmax);
		std::cout << " with range: "<<  atlasmin <<" - "<<atlasmax<<std::endl;
	}
	if (usemask) classification->SetMask(mask);
	classification->SetPadding(padding);
	classification->SetInput(image);
	classification->Initialise();

	double rel_diff;
	i=0;
	do {
		std::cout << "Iteration = " << i+1 << " / " << iterations << std::endl;
		rel_diff = classification->Iterate(i);
		i++;
	} while ((rel_diff>0.001)&&(i<iterations));

	GenericImage<int> segmentation;
	classification->ConstructSegmentation(segmentation);
	segmentation.Write(output_name);


	for( unsigned int i = 0; i < ss; ++i ){
		std::cout<<"saving probability map of structure "<<savesegsnr[i]<<" to "<<savesegs[i]<<std::endl;
		classification->WriteProbMap(savesegsnr[i],savesegs[i].c_str());
	}

	delete classification;


}


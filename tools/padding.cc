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

#include <vector>
#include <fstream>
#include <sstream>
#include <string>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage 1: " << name << " <inputA> <inputB> <output> <Value in inputB> <Padding Value in output> [<-invert>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
	std::cout << "  Changes the inputA according to inputB.  Sets the padding value where the value is found in inputB. " << std::endl;
	std::cout << "  e.g. " << name << " inputA.nii.gz inputB.nii.gz output.nii.gz 1 0" << std::endl;
	std::cout << "  If the <-invert> flag is supplied, the padding value is set where the value is NOT found. " << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------"<<std::endl;
	std::cout << std::endl;
	std::cout << "Usage 2: " << name << " <inputA> <inputB> <output> <N> <Value 1 in inputB> .. <Value N in inputB> <Padding Value in output> [<-invert>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
	std::cout << "  Changes the inputA according to inputB.  Sets the padding value where the N values are found in inputB. " << std::endl;
	std::cout << "  e.g. " << name << " inputA.nii.gz inputB.nii.gz output.nii.gz 2 1 2 100" << std::endl;
	std::cout << "  Different paddings can be entered consecutively to change different values" << std::endl;
	std::cout << "  e.g. " << name << " inputA.nii.gz inputB.nii.gz output.nii.gz 1 1 100  1 2 200  1 3 300" << std::endl;
	std::cout << "  If the <-invert> flag is supplied, the padding value is set where the value is NOT found. " << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------"<<std::endl;
	std::cout << std::endl;
	std::cout << "Usage 3: " << name << " <inputA> <inputB> <output> <csv>" << std::endl;
	std::cout << std::endl;
	std::cout << "  The csv file contains the values to change with rows:" << std::endl;
	std::cout << "    value padding_value" << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc < 5) {
		PrintHelp(EXECNAME);
	}


	const char *inputA_name = argv[1];
	const char *inputB_name = argv[2];
	const char *output_name = argv[3];

	InitializeIOLibrary();
	unique_ptr<BaseImage> inputA(BaseImage::New(inputA_name));
	unique_ptr<BaseImage> inputB(BaseImage::New(inputB_name));

	int x, y, z, t;
	double i, j, k;
	vector<int> numvalues;
	vector<double*> values;
	vector<double> padding;
	vector<bool> invert;
	int setsvalues=0;

	int restargs = 4;
	if (argc-1 == restargs){
		//read csv file
		char* csv_name=argv[restargs];
		ifstream csv(csv_name);
		if (!csv.is_open()){
			std::cerr << "Could not open file " << csv_name << std::endl;
			exit(1);
		}

		string line;
		double val1, val2;
		while ( getline(csv, line) ){
			istringstream ssline(line);
			ssline>>val1; ssline>>val2; 
			if ( val1==val2 ) continue;

			numvalues.push_back(1);
			values.push_back(new double[1]);
			values[setsvalues][0]=val1;
			padding.push_back(val2);
			invert.push_back(false);

			setsvalues++;
		}
		csv.close();
	}else{
		//process parameters
		bool onepair = false;
		if ( (argc-1 == restargs + 1) || 
				((argc-1 == restargs + 2) && (strcmp(argv[argc-1], "-invert") == 0)) )
			onepair = true;

		int nv;
		int a = restargs;
		while (a < argc){
			if(onepair){ 
				nv=1;
			}
			else{
				nv=atoi(argv[a]);
				a++;
			}
			numvalues.push_back(nv);

			values.push_back(new double[nv]);
			for(int i=0;i<nv;i++){
				values[setsvalues][i] = atof(argv[a]);
				a++;
			}

			padding.push_back(atof(argv[a]));
			a++;

			invert.push_back(false);
			if ((a <argc) && (strcmp(argv[a], "-invert") == 0)) {
				invert[setsvalues] = true;
				a++;
			}
			setsvalues++;
		}
	}

	if(setsvalues==0){
		PrintHelp(EXECNAME);
	}


	for (t = 0; t < inputA->GetT(); t++) {
		for (z = 0; z < inputA->GetZ(); z++) {
			for (y = 0; y < inputA->GetY(); y++) {
				for (x = 0; x < inputA->GetX(); x++) {
					i = x; j = y; k = z;
					inputA->ImageToWorld(i,j,k);
					inputB->WorldToImage(i,j,k);
					i = round(i); j = round(j); k = round(k);

					if (! (i >= 0 && i < inputB->GetX()
							&& j >= 0 && j < inputB->GetY()
							&& k >= 0 && k < inputB->GetZ()
							&& t >= 0 && t < inputB->GetT()) ) continue;


					double val=inputB->GetAsDouble(i, j, k, t);
					for(int s=0; s<setsvalues; s++){
						bool paddit=false;
						for(int th=0;th<numvalues[s];th++){
							if (val == values[s][th]){
								paddit=true;break;
							}
						}
						if(invert[s]) paddit=!paddit;
						if(paddit) inputA->PutAsDouble(x, y, z, t, padding[s]);
					}

				}
			}
		}
	}

	inputA->Write(output_name);

	return 0;
}


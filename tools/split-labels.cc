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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkImageIOConfig.h>
#include <mirtkGenericImage.h>
#include <mirtkHashProbabilisticAtlas.h>
#include <string>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	cout << endl;
	cout << "Usage: " << name << " <N> <labelmap_1> .. <labelmap_N>  <weight_1> .. <weight_N>    <R> <label_1> .. <label_R> <probmap_1> .. <probmap_R>" << endl;
	cout << endl;
	cout << "Description:" << endl;
	cout << "  Measures the probability of the different labels in the N label maps <labelmap_1> .. <labelmap_N> " <<endl;
	cout << "  according to the weights (maps) <weight_1> .. <weight_N> (based on occurence)." << endl;
	cout << "  It then outputs the probability of the specified R labels <label_1> .. <label_R> to <probmap_1> .. <probmap_R>  "<<endl;
	cout << endl;
	PrintStandardOptions(cout);
	cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

int main(int argc, char **argv){

	REQUIRES_POSARGS(6);
	InitializeImageIOLibrary();

	int a = 1;
	int numAtlases = atoi(POSARG(a)); a++;
	string atlasNames[numAtlases], weightNames[numAtlases];

	// atlases
	for(int j = 0; j < numAtlases; j++){
		atlasNames[j]=POSARG(a++);
	}
	// weights
	for(int j = 0; j < numAtlases; j++){
		weightNames[j]=POSARG(a++);
	}

	// structures
	int numStructuresToDo= atoi(POSARG(a)); a++;
	string names[numStructuresToDo];
	int values[numStructuresToDo];
	for(int j = 0; j < numStructuresToDo; j++){
		values[j] = atoi(POSARG(a)); a++;
	}
	for(int j = 0; j < numStructuresToDo; j++){
		names[j] = POSARG(a); a++;
	}

	// output
	GreyImage atlas(atlasNames[0].c_str());
	RealImage weight(weightNames[0].c_str());
	mirtkHashProbabilisticAtlas probs;
	for(int j = 0; j < numStructuresToDo; j++){
		probs.AddImage( HashRealImage(atlas.Attributes()) );
	}

	// process atlases
	GreyPixel *ptr;
	RealPixel *wptr, *swptr;
	RealImage sumweight(weight.Attributes());
	for(int a = 0; a < numAtlases; a++){
		if(a>0){
			atlas.Read(atlasNames[a].c_str());
			weight.Read(weightNames[a].c_str());
		}

		ptr=atlas.GetPointerToVoxels();
		wptr=weight.GetPointerToVoxels();
		swptr=sumweight.GetPointerToVoxels();
		probs.First();

		for(int i = 0; i < atlas.GetNumberOfVoxels(); i++){
			if(*wptr>0){
				for(int j = 0; j < numStructuresToDo; j++){
					if(*ptr == values[j]){
						//if(a==1) cout<<*wptr + probs.GetValue(j)<<"/"<<*swptr+*wptr<<endl;
						probs.SetValue(j, *wptr + probs.GetValue(j));
						break;
					}
				}
				*swptr+=*wptr;
			}
			ptr++; wptr++; swptr++;
			probs.Next();
		}
	}

	// normalise output
	probs.First();
	swptr=sumweight.GetPointerToVoxels();
	for(int i = 0; i < atlas.GetNumberOfVoxels(); i++){
		if(*swptr>0){
			for(int j = 0; j < numStructuresToDo; j++){
				probs.SetValue(j, probs.GetValue(j) * 100.0 / *swptr );
			}
		}
		swptr++;
		probs.Next();
	}

	// write output
	for(int j = 0; j < numStructuresToDo; j++){
		probs.Write(j, names[j].c_str());
	}

	return 0;
}


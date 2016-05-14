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

#include <map>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input>" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
 	std::cout << "  Measures the volume of each label in the input image" << std::endl;
	std::cout << std::endl;
	std::cout << "Optional: " << std::endl;
	std::cout << std::endl;
	std::cout << "  -voxels    count the number of voxels instead" << std::endl;
	PrintStandardOptions(std::cout);
	std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------


int main(int argc, char **argv)
{
  REQUIRES_POSARGS(1);

  InitializeIOLibrary();
  RealImage img(POSARG(1));

  bool involume=true;
    for (ALL_OPTIONS) {
        if (OPTION("-voxels")){
            involume=false;
        }
        else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }

  map<int, int> voxels;

    for(int x = 0; x < img.GetX(); x++){
		for(int y = 0; y < img.GetY(); y++){
			for(int z = 0; z < img.GetZ(); z++){
			    int label=img.Get(x, y, z, 0);
			    if(label==0)continue;
			    voxels[label]++;
			}
		}
    }

    double xvsize, yvsize, zvsize;
    img.GetPixelSize(&xvsize,&yvsize,&zvsize); 
    double vsize=xvsize* yvsize* zvsize;	

    for( map<int,int>::iterator ii=voxels.begin(); ii!=voxels.end(); ++ii){
		int label=(*ii).first;
		if(label==0)continue;
		double volume=(*ii).second ;
		if(involume) volume=volume*vsize;
	       std::cout << label << " " <<std::setprecision(8)<< volume << std::endl;
    }
    
    return 0;
}

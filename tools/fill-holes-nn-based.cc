/*
 * Medical Image Registration ToolKit (MIRTK)
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

// Source code adapted from FSL fslmaths -fillh


#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

#include "mirtk/ConnectedComponents.h"
#include <set>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
    std::cout << "Usage: " << name << " <input> <suspected-holes> <output>" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
    std::cout << "  Fills holes in the input."<<std::endl;
    std::cout << "  The surrounding voxels of the suspected-holes are measured and if the majority belongs to the input they are filled." <<std::endl;
    std::cout << std::endl;
    std::cout << "Input options:" << std::endl;
    std::cout << "  -connectivity <number>  connectivity (6/26), default: 6" <<std::endl;
    std::cout << "  -majority <double>      majority factor (0-1), default: 0.9" <<std::endl;
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

    int a = 1;
    GreyImage image(POSARG(a++));
    GreyImage suspected(POSARG(a++));
    char *ouput_name=POSARG(a);
    ConnectivityType connectivity= (ConnectivityType) 6;
    double majority=0.9;

    for (ALL_OPTIONS) {
        if (OPTION("-connectivity")){
            connectivity = (ConnectivityType) atoi(ARGUMENT);
	    if(connectivity!=6 && connectivity!=26){ std::cerr<<"connectivity has to be either 6 or 26!"<<std::endl; exit(1);}
        }
        else if (OPTION("-majority")){
            majority=atof(ARGUMENT);
        }
        else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }

    GreyImage suspectedcc(image.Attributes());
    ConnectedComponents<GreyPixel> cc(CC_LargestFirst, connectivity);
    cc.Input(&suspected);
    cc.Output(&suspectedcc);
    cc.Run();

    int numccs=cc.NumberOfComponents();
    int numccNNs[numccs];
    int numccMaskNNs[numccs];
    for(int i=0; i<numccs; i++){
	numccNNs[i]=0;
	numccMaskNNs[i]=0;
    }


    int conn[27];
    for(int i=0; i<27; i++) conn[i]=0;

    if(connectivity==6){
      double dist;
      int i=0;
      for(int z=-1; z<=1; z++)
      for(int y=-1; y<=1; y++)
      for(int x=-1; x<=1; x++){
	dist=x*x+y*y+z*z;
	if(dist==1) conn[i]=1;
	i++;
      }
    }else {
    	for(int i=0; i<27; i++) conn[i]=1;
    	conn[13]=0;
    }


    int xmax=suspectedcc.GetX()-1;
    int ymax=suspectedcc.GetY()-1;
    int zmax=suspectedcc.GetZ()-1;
    int vlabel, vnnlabel, wlabel;
	for( int x = 0; x < xmax; x++){
	  for( int y = 0; y < ymax; y++){
	    for( int z = 0; z < zmax; z++){
                vlabel=suspectedcc.Get( x ,y, z)-1;
		if(vlabel==-1) continue;

	      int nn=0;
	      for(int zd=-1; zd<=1; zd++)
	      for(int yd=-1; yd<=1; yd++)
	      for(int xd=-1; xd<=1; xd++){
		if(conn[nn] && x+xd>=0 && y+yd>=0 && z+zd>=0 && x-xd<=xmax && y+yd<=ymax && z+zd<=zmax ){
                    vnnlabel=suspectedcc.Get( x+xd , y+yd, z+zd)-1;
		    if(vnnlabel != vlabel){
			numccNNs[vlabel]++;
                    	wlabel=image.Get( x+xd , y+yd, z+zd);
		    	if(wlabel>0) numccMaskNNs[vlabel]++;
		    }
		}
		nn++;
	      }
	  }
	}
    }

    bool fill[numccs];
    for(int i=0; i<numccs; i++)
	fill[i]=( ((double)numccMaskNNs[i]/(double)numccNNs[i]) >= majority);
	

    GreyPixel *p=image.GetPointerToVoxels();
    GreyPixel *c=suspectedcc.GetPointerToVoxels();
    int number_of_voxels=image.GetNumberOfVoxels();
    for (int i=0; i<number_of_voxels; i++) {
	*p=(*p>0);
	if(*c>0 && fill[*c-1]) *p=1;
        p++; c++;
    }

    image.Write(ouput_name);

	return 0;
}


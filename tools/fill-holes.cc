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

// Source code adapted from FSL fslmaths -fillh


#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkImageIOConfig.h>
#include <mirtkGenericImage.h>

#include <mirtkConnectedComponents.h>
#include <set>

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	cout << endl;
    cout << "Usage: " << name << " <input> <output>" << endl;
	cout << endl;
	cout << "Description:" << endl;
    cout << "  Fills holes in the input."<<endl;
    cout << "  Note: The code is adapted from fslmaths -fillh" <<endl;
    cout << endl;
    cout << "Input options:" << endl;
    cout << "  -connectivity <number>  voxel connectivity for finding holes - 6 or 26 (default: 6)" <<endl;
    cout << endl;
	PrintStandardOptions(cout);
	cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

int main(int argc, char **argv){

    REQUIRES_POSARGS(2);
	InitializeImageIOLibrary();

    int a = 1;

    RealImage image(POSARG(a++));
    char *ouput_name=POSARG(a);
    ConnectivityType connectivity = (ConnectivityType) 6;

    for (ALL_OPTIONS) {
        if (OPTION("-connectivity")){
            connectivity = (ConnectivityType) atoi(ARGUMENT);
        }
        else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }


    GreyImage mask(image.Attributes());
    RealPixel *p=image.GetPointerToVoxels();
    GreyPixel *m=mask.GetPointerToVoxels();
    int number_of_voxels=image.GetNumberOfVoxels();
    for (int i=0; i<number_of_voxels; i++) {
        *m=(*p==0);
        p++; m++;
    }

    ConnectedComponents<GreyPixel> cc(CC_LargestFirst, connectivity);
    cc.Input(&mask);
    cc.Output(&mask);
    cc.Run();

    int maxx=mask.GetX()-1;
    int maxy=mask.GetY()-1;
    int maxz=mask.GetZ()-1;

    set<int> edgelabs;
      for (int z=0; z<=maxz; z++) {
        for (int y=0; y<=maxy; y++) {
          if (mask.Get(0,y,z)!=0) edgelabs.insert(mask.Get(0,y,z));
          if (mask.Get(maxx,y,z)!=0) edgelabs.insert(mask.Get(maxx,y,z));
        }
      }
      for (int z=0; z<=maxz; z++) {
        for (int x=0; x<=maxx; x++) {
          if (mask.Get(x,0,z)!=0) edgelabs.insert(mask.Get(x,0,z));
          if (mask.Get(x,maxy,z)!=0) edgelabs.insert(mask.Get(x,maxy,z));
        }
      }
      for (int y=0; y<=maxy; y++) {
        for (int x=0; x<=maxx; x++) {
          if (mask.Get(x,y,0)!=0) edgelabs.insert(mask.Get(x,y,0));
          if (mask.Get(x,y,maxz)!=0) edgelabs.insert(mask.Get(x,y,maxz));
        }
      }

      while (!edgelabs.empty()) {
        int val=*(edgelabs.begin());
        for (int z=0; z<=maxz; z++) {
          for (int y=0; y<=maxy; y++) {
        for (int x=0; x<=maxx; x++) {
          if (mask.Get(x,y,z)==val) mask.Put(x,y,z,0);
        }
          }
        }
        edgelabs.erase(val);
      }

      for (int z=0; z<=maxz; z++) {
        for (int y=0; y<=maxy; y++) {
          for (int x=0; x<=maxx; x++) {
             if (mask.Get(x,y,z)>0.5) image.Put(x,y,z,1);
          }
        }
      }

    image.Write(ouput_name);

	return 0;
}


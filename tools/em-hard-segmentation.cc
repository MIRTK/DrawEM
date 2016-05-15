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

#include "mirtk/HashProbabilisticAtlas.h"


using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
    std::cout << std::endl;
    std::cout << "Usage: " << name << " <N> <atlas1> .. <atlasN> <output> [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Description:" << std::endl;
    std::cout << "  Computes the hard segmentation of the N atlases. " << std::endl;
    std::cout << "  Optionally MRF smoothing can be additionally applied." << std::endl;
    std::cout << std::endl;
    std::cout << "Input options:" << std::endl;
    std::cout << "  -mrftimes <number>              number of times the mrf term will be applied (default 0)" << std::endl;
    std::cout <<	"  -mrfweight <double>             weight of the mrf term (default 1/3)." << std::endl;
    std::cout <<	"  -posteriors <post1> .. <postN>  write posteriors (useful when MRF is used)." << std::endl;
    std::cout << std::endl;
    PrintStandardOptions(std::cout);
    std::cout << std::endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------

double dx,dy,dz;
double sx, sy,sz;
int maxx,maxy,maxz;
double mrfweight=1.0/3.0;
int n;
HashProbabilisticAtlas atlas;





double getMRFenergyTerm(int x, int y, int z, int k)
{
	int lx = max(x-1,0);
    int rx = min(x+1, maxx-1 );
	int ly = max(y-1,0);
    int ry = min(y+1, maxy-1 );
	int lz = max(z-1,0);
    int rz = min(z+1, maxz-1 );
    double term= sx * ( atlas.GetValue(rx,y,z,k) + atlas.GetValue(lx,y,z,k) ) + sy * ( atlas.GetValue(x,ry,z,k) + atlas.GetValue(x,ly,z,k) ) + sz * ( atlas.GetValue(x,y,rz,k) + atlas.GetValue(x,y,lz,k) );
	return term;
}


void MRFStep(){
	double mrfterms[n];
	double numerator[n];
	double denominator;
        for(int x = 0; x < maxx; x++){
          for(int y = 0; y < maxy; y++){
            for(int z = 0; z < maxz; z++){

				for (int k = 0; k < n; k++)
					mrfterms[k] = getMRFenergyTerm(x,y,z,k);
				
				denominator=0;
				for (int k = 0; k < n; k++){
					double energy = .0;
					for (int t = 0; t < n; t++){
						if(t==k)continue;
						energy += mrfterms[t];
					}
                    double expo= exp(-mrfweight * energy);
                    numerator[k] = atlas.GetValue(x,y,z,k) *expo;
					denominator+=numerator[k];
				}

				if(denominator>0)
                for (int k = 0; k < n; k++){
                    atlas.SetValue(x,y,z,k,numerator[k]/denominator);
				}
			}
		}
	}
	      
}



int main(int argc, char **argv){
    REQUIRES_POSARGS(4);
    InitializeIOLibrary();
    int a=1;

  // Number of tissues
  n = atoi(POSARG(a++));

  // Read atlas for each tissue
  double min, max;
  for (int i = 0; i < n; i++) {
    RealImage image(POSARG(a));
    std::cerr << "Image " << i <<" = " << POSARG(a);
    image.GetMinMaxAsDouble(&min, &max);
    std::cout << " with range: "<<  min <<" - "<<max<<std::endl;
    atlas.AddImage(image);
    a++;
  }

  // File name for output
  char *output_name = POSARG(a++);

  string *post = new string[n];
  bool saveposts=false;
  int mrftimes=0;
  for (ALL_OPTIONS) {
      if (OPTION("-mrfweight")){
          mrfweight = atof(ARGUMENT);
      }
      else if (OPTION("-mrftimes")) {
          mrftimes = atoi(ARGUMENT);
      }
      else if (OPTION("-posteriors")) {
          saveposts = true;
          for (int i = 0; i < n; i++) {
            post[i] = ARGUMENT;
          }
      }
      else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  RealImage sum;
  if(mrftimes>0){
      sum=atlas.GetImage(0).ToGenericImage();
      for (int i = 1; i < n; i++){
          HashRealImage::DataIterator begin=atlas.Begin(i), end=atlas.End(i);
          for ( auto it = begin; it != end; ++it ){
              sum.Put(it->first, sum.Get(it->first)+ it->second);
          }
      }

      maxx=sum.GetX();
      maxy=sum.GetY();
      maxz=sum.GetZ();

      sum.GetPixelSize(&dx,&dy,&dz);
      sx = 1.0/dx;
      sy = 1.0/dy;
      sz = 1.0/dz;

    for(int i=1;i<=mrftimes;i++){
      std::cout<<"MRF iteration "<<i<<"...";
      MRFStep();
      std::cout<<"done"<<std::endl;
    }
  }

  HashImage<int> seg = atlas.ComputeHardSegmentation();
  seg+=1;
  seg.Write(output_name);



  if(saveposts){
      for (int i = 0; i < n; i++){
          RealImage image=atlas.GetImage(i).ToGenericImage();
          if(mrftimes>0){
              RealPixel *ptr = image.GetPointerToVoxels();
              RealPixel *sptr = sum.GetPointerToVoxels();
              for (int v = 0; v < image.GetNumberOfVoxels(); v++) {
                  *ptr *= *sptr;
                  ptr++; sptr++;
              }
          }
          image.Write(post[i].c_str());
       }
  }
	  
  delete[] post;


}


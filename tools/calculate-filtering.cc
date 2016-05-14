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

#include <algorithm>   
#include <vector>      

using namespace mirtk;
using namespace std;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	std::cout << std::endl;
	std::cout << "Usage: " << name << " <input> <options>" << std::endl;
	std::cout << std::endl;
	std::cout << "Description:" << std::endl;
 	std::cout << "  Calculates statistics by filtering with a kernel." << std::endl;
	std::cout << std::endl;
	std::cout << "Options: " << std::endl;
	std::cout << "  -kernel <number>       kernel size: number^3 (number must be even, and >=3 !), default: 3 " << std::endl;
	std::cout << std::endl;
	std::cout << "Operations: " << std::endl;
	std::cout << "  -min <output>          calculate min" << std::endl;
	std::cout << "  -max <output>          calculate max" << std::endl;
	std::cout << "  -mean <output>         calculate mean" << std::endl;
	std::cout << "  -median <output>       calculate median" << std::endl;
	std::cout << "  -std <output>          calculate std" << std::endl;
	std::cout << "  -std-median <output>   calculate std based on median" << std::endl;
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

  char *min_name, *max_name, *mean_name, *median_name, *std_name, *std_median_name;
  RealImage *min=NULL, *max=NULL, *mean=NULL, *median=NULL, *std=NULL, *std_median=NULL;
  bool ok=false;
  int kernel=1;

    for (ALL_OPTIONS) {
        if (OPTION("-kernel")){
           kernel=atoi(ARGUMENT); 
	   if(kernel%2==0 || kernel<3) HANDLE_STANDARD_OR_UNKNOWN_OPTION();
	   kernel=float((kernel-1)/2);
        }
        else if (OPTION("-min")){
            min_name=ARGUMENT; min=new RealImage(img.Attributes());
	    ok=true;
        }
        else if (OPTION("-max")){
            max_name=ARGUMENT; max=new RealImage(img.Attributes());
	    ok=true;
        }
        else if (OPTION("-mean")){
            mean_name=ARGUMENT; mean=new RealImage(img.Attributes());
	    ok=true;
        }
        else if (OPTION("-median")){
            median_name=ARGUMENT; median=new RealImage(img.Attributes());
	    ok=true;
        }
        else if (OPTION("-std")){
            std_name=ARGUMENT; std=new RealImage(img.Attributes());
	    ok=true;
        }
        else if (OPTION("-std_median")){
            std_median_name=ARGUMENT; std_median=new RealImage(img.Attributes());
	    ok=true;
        }
        else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }
    if(!ok){ PrintHelp(EXECNAME); exit(1); }


    int x,x1,x2,xn, y,y1,y2,yn, z,z1,z2,zn, i;
    double minv,maxv,meanv,medianv,stdv,std_medianv,v;
    for(x = 0; x < img.GetX(); x++){
      for(y = 0; y < img.GetY(); y++){
	for(z = 0; z < img.GetZ(); z++){
	  x1=std::max(0,x-kernel);
	  y1=std::max(0,y-kernel);
	  z1=std::max(0,z-kernel);
	  x2=std::min(x+kernel,img.GetX()-1);
	  y2=std::min(y+kernel,img.GetY()-1);
	  z2=std::min(z+kernel,img.GetZ()-1);

	  vector<RealPixel> values;
	  for(xn = x1; xn <=x2; xn++){
      	    for(yn = y1; yn <=y2; yn++){
	      for(zn = z1; zn <=z2; zn++){
		values.push_back(img.Get(xn, yn, zn));
	      }
	    }
	  }

	  minv=DBL_MAX; maxv=DBL_MIN; meanv=0;
	  for(i=0; i<values.size(); i++){
	     v=values[i];
	     if(min) minv=std::min(minv, v);
	     if(max) maxv=std::max(maxv, v);
	     if(mean || std) meanv+=v;
	  }
	  if(mean || std) meanv/=(double)values.size();
	  if(median || std_median){ 
	     sort(values.begin(), values.end());
	     medianv=values[floor(values.size()/2)];
	  }

	  if(std){
	     stdv=0;
	     for(i=0; i<values.size(); i++){
    		stdv += pow(values[i]-meanv, 2);
	     }
	     stdv/=values.size();
	  }
	  if(std_median){
	     std_medianv=0;
	     for(i=0; i<values.size(); i++){
    		std_medianv += pow(values[i]-medianv, 2);
	     }
	     std_medianv/=values.size();
	  }

	  if(min) min->Put(x,y,z, minv);
	  if(max) max->Put(x,y,z, maxv);
	  if(mean) mean->Put(x,y,z, meanv);
	  if(median) median->Put(x,y,z, medianv);
	  if(std) std->Put(x,y,z, stdv);
	  if(std_median) std_median->Put(x,y,z, std_medianv);

	}
      }
    }


    if(min){ min->Write(min_name); delete min;}
    if(max){ max->Write(max_name); delete max;}
    if(mean){ mean->Write(mean_name); delete mean;}
    if(median){ median->Write(median_name); delete median;}
    if(std){ std->Write(std_name); delete std;}
    if(std_median){ std_median->Write(std_median_name); delete std_median;}
    
    return 0;
}

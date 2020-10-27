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

#include <iostream>
#include <fstream>
#include <map>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

using namespace mirtk;


void usage()
{
  std::cerr << "Usage: label-connectivity [num_labels] [N] [atlas1_labels.nii.gz atlas1_labels.nii.gz .. atlasN_labels.nii.gz] [output]" << endl;
  exit(1);
}



int main(int argc, char **argv){

  if (argc < 5) {
    usage();
  }
  InitializeIOLibrary();

    int x,y,z,k,j,n;

    int numlabels =atoi(argv[1]);
    argc--;
    argv++;

    int conn[numlabels+1][numlabels+1];
    for (int i=0;i<numlabels+1;i++)
        for (int j=0;j<numlabels+1;j++){
            conn[i][j]=0;
        }


    int numatlases =atoi(argv[1]);
    argc--;
    argv++;

    int vlabel,rxlabel,lxlabel,rylabel,lylabel,rzlabel,lzlabel,wlabel;
    int lx,rx,ly,ry,rz,lz;


    for (int a=0;a<numatlases;a++){
        GreyImage atlas(argv[1]);
        argc--;
        argv++;


    int xmax=atlas.GetX()-1;
    int ymax=atlas.GetY()-1;
    int zmax=atlas.GetZ()-1;

	for( int x = 0; x < atlas.GetX(); x++){
	  for( int y = 0; y < atlas.GetY(); y++){
		  for( int z = 0; z < atlas.GetZ(); z++){
                vlabel=atlas.Get( x ,y, z);

                if(x-1>=0){
                    wlabel=atlas.Get( x-1 ,y,z);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }
                if(x+1<=xmax){
                    wlabel=atlas.Get( x+1 ,y,z);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }
                if(y-1>=0){
                    wlabel=atlas.Get( x,y-1,z);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }
                if(y+1<=ymax){
                    wlabel=atlas.Get( x,y+1,z);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }if(z-1>=0){
                    wlabel=atlas.Get( x,y,z-1);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }
                if(z+1<=zmax){
                    wlabel=atlas.Get( x,y,z+1);
                    if(vlabel!=wlabel)conn[vlabel][wlabel]++;
                }

		  }
	  }
	}
    }


    char *filename=argv[1];
    ofstream outfile(filename);
    for (int i=1;i<numlabels+1;i++){
        for (int j=1;j<numlabels+1;j++){
            outfile<<conn[i][j];
            if(j!=numlabels)outfile<<" ";
        }outfile<<endl;
    }
    outfile.close();

  return 0;
}

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

#ifndef _MIRTKKMEANS_H
#define _MIRTKKMEANS_H

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <algorithm>

using namespace std;

class kmeans
{

protected:
	int k,iter;
	double * centroid;
	double * point;
	int * pointcluster, * oldpointcluster;
	int n;
	bool converged;
	double * bestcentroid;
	int * bestpointcluster;
	double bestsumdistances;
	int *nb;

	void KMeansAssign();
	void KMeansCluster();

public:
	//method to run
	kmeans();

	kmeans(double *thepoints,int numpoints,int numclusters,int numiters=100,int replicates=10);
	kmeans(double *thepoints,int numpoints,int numclusters,double* centroids_init,int numiters=100,int replicates=10);
	double *getCentroids();
	int *getPointClusters();

};



inline double * kmeans::getCentroids(){return bestcentroid;}
inline int * kmeans::getPointClusters(){return bestpointcluster;}




#endif







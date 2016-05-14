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


#include "mirtk/KMeans.h"

// Hard membership
void kmeans::KMeansAssign()
{
	int winner;
	double dist,distmin;

	for(int i=0;i<n;i++)
	{
		distmin=DBL_MAX;
		for(int j=0;j<k;j++)
		{
			dist=fabs(point[i]-centroid[j]);

			if (dist<distmin) {distmin=dist; winner=j;}
		}
		pointcluster[i]=winner;
		if(pointcluster[i] != oldpointcluster[i]){
			converged=false;
			oldpointcluster[i]=pointcluster[i];
		}
	}
}


void kmeans::KMeansCluster()
{
	int i,j,cl;

	for(j=0;j<k;j++)
	{
		nb[j]=0;
		centroid[j]=0;
	}

	for(i=0;i<n;i++)
	{
		cl=pointcluster[i];
		centroid[cl]+=point[i];
		nb[cl]++;
	}

	for(j=0;j<k;j++){
		if (nb[j]>0) centroid[j]/=nb[j];
	}
}



kmeans::kmeans(double *thepoints,int numpoints,int numclusters,int numiters,int replicates){
	n=numpoints;
	k=numclusters;
	point=thepoints;


	nb=new int[k];
	centroid=new double[k];
	pointcluster=new int[n];
	oldpointcluster=new int[n];


	bestsumdistances=DBL_MAX;
	bestcentroid=new double[k];
	bestpointcluster=new int[n];

	double bestsumdistances0=DBL_MAX;
	double *bestcentroid0=new double[k];
	int *bestpointcluster0=new int[n];


	for(int rep=1;rep<=replicates;rep++){
		for(int i=0;i<n;i++)
			oldpointcluster[i]=-1;

		for(int i=0;i<k;i++){
			centroid[i]=point[rand()%n];
		}

		sort(centroid, centroid + k);



		for(int iter=1;iter<=numiters;iter++){
			converged=true;
			KMeansAssign();
			KMeansCluster();

			if(converged)break;
		}

		std::sort(centroid, centroid + k);
		KMeansAssign();


		double sumdistances=0;
		for(int i=0;i<n;i++)
			sumdistances+=fabs(point[i]-centroid[pointcluster[i]]);


		if(sumdistances<bestsumdistances0){
			bestsumdistances0=sumdistances;
			for(int i=0;i<n;i++) bestpointcluster0[i]=pointcluster[i];
			for(int i=0;i<k;i++) bestcentroid0[i]=centroid[i];
		}

		bool zerofail=0;
		for(int i=0;i<k;i++)if(nb[i]==0){
			zerofail=true;
			break;}
		if(sumdistances<bestsumdistances && !zerofail){
			bestsumdistances=sumdistances;
			for(int i=0;i<n;i++) bestpointcluster[i]=pointcluster[i];
			for(int i=0;i<k;i++) bestcentroid[i]=centroid[i];
		}
	}

	if(bestsumdistances==DBL_MAX){
		bestsumdistances=bestsumdistances0;
		for(int i=0;i<n;i++) bestpointcluster[i]=bestpointcluster0[i];
		for(int i=0;i<k;i++) bestcentroid[i]=bestcentroid0[i];
	}
}











kmeans::kmeans(double *thepoints,int numpoints,int numclusters,double* centroids_init,int numiters,int replicates){
	n=numpoints;
	k=numclusters;
	point=thepoints;




	nb=new int[k];
	centroid=new double[k];
	pointcluster=new int[n];
	oldpointcluster=new int[n];


	bestsumdistances=DBL_MAX;
	bestcentroid=new double[k];
	bestpointcluster=new int[n];

	double bestsumdistances0=DBL_MAX;
	double *bestcentroid0=new double[k];
	int *bestpointcluster0=new int[n];

	double mindiff=1e10;
	for(int i=0;i<k;i++){
		for(int j=0;j<k;j++){
			if(i==j)continue;
			double diff=fabs(centroids_init[i]-centroids_init[j]);
			if(diff<mindiff) mindiff=diff;
		}
	}


	int change=round(mindiff);
	for(int rep=1;rep<=replicates;rep++){
		for(int i=0;i<n;i++)
			oldpointcluster[i]=-1;

		for(int i=0;i<k;i++){
			centroid[i]=centroids_init[i] + rand() % change;
		}

		sort(centroid, centroid + k);



		for(int iter=1;iter<=numiters;iter++){
			converged=true;
			KMeansAssign();
			KMeansCluster();

			if(converged)break;
		}

		std::sort(centroid, centroid + k);
		KMeansAssign();


		double sumdistances=0;
		for(int i=0;i<n;i++)
			sumdistances+=fabs(point[i]-centroid[pointcluster[i]]);


		if(sumdistances<bestsumdistances0){
			bestsumdistances0=sumdistances;
			for(int i=0;i<n;i++) bestpointcluster0[i]=pointcluster[i];
			for(int i=0;i<k;i++) bestcentroid0[i]=centroid[i];
		}

		bool zerofail=0;
		for(int i=0;i<k;i++)if(nb[i]==0){
			zerofail=true;
			break;}
		if(sumdistances<bestsumdistances && !zerofail){
			bestsumdistances=sumdistances;
			for(int i=0;i<n;i++) bestpointcluster[i]=pointcluster[i];
			for(int i=0;i<k;i++) bestcentroid[i]=centroid[i];
		}

	}

	if(bestsumdistances==DBL_MAX){
		bestsumdistances=bestsumdistances0;
		for(int i=0;i<n;i++) bestpointcluster[i]=bestpointcluster0[i];
		for(int i=0;i<k;i++) bestcentroid[i]=bestcentroid0[i];
	}
}

kmeans::kmeans(){}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


//find in which grid box the point at the given coordinates is at. returns an array with 3 int elements.
void findinwhichbox(float x,float y,float z,int n,int m, int k, int *ret){ 
	ret[0] = (int)(x*n);
	ret[1] = (int)(y*m);
	ret[2] = (int)(z*k);
	return;
}

//free allocated 3d memory - not used for now, might be useful in the future
// this is copy-paste code, and it just works
void free_3d(float ***data, size_t xlen, size_t ylen)
{
	size_t i, j;

	for (i=0; i < xlen; ++i) {
		if (data[i] != NULL) {
			for (j=0; j < ylen; ++j)
				free(data[i][j]);
			free(data[i]);
		}
	}
	free(data);
}

// allocate a 3d float table, with fixed dimensions (ie no different line size for each line)
// this is copy-paste code, and it just works
float ***alloc_3d(size_t xlen, size_t ylen, size_t zlen)
{
	float ***p;
	size_t i, j;

	if ((p = malloc(xlen * sizeof *p)) == NULL) {
		perror("malloc 1");
		return NULL;
	}

	for (i=0; i < xlen; ++i)
		p[i] = NULL;

	for (i=0; i < xlen; ++i)
		if ((p[i] = malloc(ylen * sizeof *p[i])) == NULL) {
			perror("malloc 2");
			free_3d(p, xlen, ylen);
			return NULL;
		}

	for (i=0; i < xlen; ++i)
		for (j=0; j < ylen; ++j)
			p[i][j] = NULL;

	for (i=0; i < xlen; ++i)
		for (j=0; j < ylen; ++j)
			if ((p[i][j] = malloc(zlen * sizeof *p[i][j])) == NULL) {
				perror("malloc 3");
				free_3d(p, xlen, ylen);
				return NULL;
			}

	return p;
}


int main(){

	srand (time(NULL));  //such randomness wow
	
	int Pprocesses=4; //from 0 to 7
	int processes=1<<Pprocesses;
	int Pnumberofpoints=15; // from 0(?) to 25
	int numberofpoints=1<<Pnumberofpoints;
	int Pnumberboxes=7; //from 12 to 16
	int numboxes=1<<Pnumberboxes;
	int Pnumboxesperprocess=Pnumberboxes-Pprocesses; //2^x number of grid boxes per process
	int numboxesperprocess=1<<Pnumboxesperprocess;
	
	int processid=2; //for testing purposes - this should be dynamically allocated
	
	
	// coordinates will be in here. data example: {x1, y1, z1, x2, y2, z2, ...}
	float *q = malloc(sizeof(float) * (3 * numberofpoints / processes));
	float *c = malloc(sizeof(float) * (3 * numberofpoints / processes));
	
	
	for (int i=0;i<numberofpoints*3/processes;i++){
		q[i]=(float)rand() / RAND_MAX; //random value from 0 to 1
		c[i]=(float)rand() / RAND_MAX;
	}
	
	/*for (int i = 0; i < numberofpoints/processes; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f ", q[3*i+j]);
		}
	printf("\n");
	}*/
		
	//This will simulate our box grid, creating 2^pnumberboxes boxes. 
	//Each boxdimension shows how many boxes are on each dimension of our cube.
	//n<=m<=k is kept.
	//Trust me, self of the future, I did the math.
	int boxdimensions[3];
	boxdimensions[0]=1<<(Pnumberboxes/3);
	if (Pnumberboxes%3 > 0)	boxdimensions[0]*=2; //double the boxes in this row if mod3 is 1 or 2
	boxdimensions[1]=1<<(Pnumberboxes/3);
	if (Pnumberboxes%3 == 2)	boxdimensions[1]*=2; //double the boxes in this row if mod3 is 2
	boxdimensions[2]=1<<(Pnumberboxes/3);
	printf("grid dimensions: %d, %d, %d\n",boxdimensions[0],boxdimensions[1],boxdimensions[2]);
	/*
	if (boxdimensions[0]*boxdimensions[1]*boxdimensions[2]!=numboxes) {//fuckup in the math
		printf("mathfuckup\n"); 
		return 1;
	} 
	*/
	
	//same idea, but with splits. each splitted part of the grid will go to another process
	int splitdimensions[3];
	int gridsplitsize[3];
	splitdimensions[0]=1<<(Pprocesses/3);
	if (Pprocesses%3 > 0)	splitdimensions[0]*=2; //double the boxes in this row if mod3 is 1 or 2
	splitdimensions[1]=1<<(Pprocesses/3);
	if (Pprocesses%3 == 2)	splitdimensions[1]*=2; //double the boxes in this row if mod3 is 2
	splitdimensions[2]=1<<(Pprocesses/3);	
	for(int i=0;i<3;i++) gridsplitsize[i]=boxdimensions[i]/splitdimensions[i];
	printf("grid split size: %d, %d, %d\n",gridsplitsize[0],gridsplitsize[1],gridsplitsize[2]);
	
		
	//these could be int, but then i'd need another alloc_3d for int type - maybe ill do it later	
	//these will keep tabs on how many points are in each grid box
	float ***qgridcount; //how many points in each grid box
	float ***cgridcount; //how many points in each grid box
	if ((qgridcount = alloc_3d((size_t)boxdimensions[0], (size_t)boxdimensions[1], (size_t)boxdimensions[2])) == NULL) 
		return 1; //malloc fail	
	if ((cgridcount = alloc_3d((size_t)boxdimensions[0], (size_t)boxdimensions[1], (size_t)boxdimensions[2])) == NULL) 
		return 1; //malloc fail
	for (int i=0; i < boxdimensions[0]; ++i)
		for (int j=0; j < boxdimensions[1]; ++j)
			for (int k=0; k < boxdimensions[2]; ++k){
				qgridcount[i][j][k]=cgridcount[i][j][k]=0;
			}
			
	//find where each point belongs
	int *qbox = malloc(sizeof(int) * (3 * numberofpoints / processes));
	int *cbox = malloc(sizeof(int) * (3 * numberofpoints / processes));
	
	for (int i=0;i<numberofpoints/processes;i++){
		findinwhichbox(q[3*i],q[3*i+1],q[3*i+2],boxdimensions[0],boxdimensions[1],boxdimensions[2], &qbox[3*i]);
		findinwhichbox(c[3*i],c[3*i+1],c[3*i+2],boxdimensions[0],boxdimensions[1],boxdimensions[2], &cbox[3*i]);
		//printf("Point %d is at grid box: %d, %d, %d\n",i+1, qbox[3*i],qbox[3*i+1],qbox[3*i+2]);
		//count number of points in each box
		qgridcount[qbox[3*i]][qbox[3*i+1]][qbox[3*i+2]]+=1;
		cgridcount[cbox[3*i]][cbox[3*i+1]][cbox[3*i+2]]+=1;
	}
	
	/*
	for (int i=0; i < boxdimensions[0]; ++i)
		for (int j=0; j < boxdimensions[1]; ++j)
			for (int k=0; k < boxdimensions[2]; ++k){
				printf("Box %d,%d,%d\nQ: %f\nC: %f\n\n",i+1,j+1,k+1,qgridcount[i][j][k],cgridcount[i][j][k]);
			}*/
	return 0;
}

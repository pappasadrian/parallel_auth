#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

//global variables
int boxdimensions[3];

//get box coordinates from box id
void getboxcoords(int id, int *coords){
	coords[0]=id%boxdimensions[0];
	coords[1]=(id/boxdimensions[0])%boxdimensions[1];
	coords[2]=((id/boxdimensions[0])/boxdimensions[1])%boxdimensions[2];
}
//get box id from coordinates
int getboxid(int x, int y, int z){
	int i= x+y*(boxdimensions[0])+z*boxdimensions[0]*boxdimensions[1];
	return i;
}

//find in which grid box the point at the given coordinates is at. returns box id
int findinwhichbox(float x,float y,float z){ 
	int ret[3];
	ret[0] = (int)(x*boxdimensions[0]); //rounds down to the nearest int
	ret[1] = (int)(y*boxdimensions[1]);
	ret[2] = (int)(z*boxdimensions[2]);
	return getboxid(ret[0],ret[1],ret[2]);
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


int main(int argc, char **argv){

	if (argc != 4) {
		printf("Usage: %s Q P B\n",argv[0]);
		printf("Q: number of points (power of two - from 20 to 25)\n");
		printf("P: number of processes (power of two - from 0 to 7)\n");
		printf("B: number of boxes (power of two - from 12 to 16)\n");
		return 1;
	}
  
	srand (time(NULL));  //such randomness wow
	
	int Pprocesses=atoi(argv[2]); //from 0 to 7
	int processes=1<<Pprocesses;
	int Pnumberofpoints=atoi(argv[1]); // from 0(?) to 25 -in all processes. this process has numberofpoints/processes points.
	int numberofpoints=1<<Pnumberofpoints;
	int Pnumberboxes=atoi(argv[3]); //from 12 to 16
	int numboxes=1<<Pnumberboxes;
	int Pnumboxesperprocess=Pnumberboxes-Pprocesses; //2^x number of grid boxes per process, aka splits
	int numboxesperprocess=1<<Pnumboxesperprocess;
	
	int processid=2; //for testing purposes - this should be dynamically allocated
	
	
	// coordinates will be in here. data example: {x1, y1, z1, x2, y2, z2, ...}
	float *q = malloc(sizeof(float) * (3 * numberofpoints / processes)); //numberofpoints/processes, because the points are spread through the active processes
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
	//generally, we can refer to these boxes either by ID or by coordinates. there are functions that do this. 
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
	//i am not sure why i did this. but it might get handy. 
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
		int idq = findinwhichbox(q[3*i],q[3*i+1],q[3*i+2]);
		int idc = findinwhichbox(c[3*i],c[3*i+1],c[3*i+2]);
		getboxcoords(idq, &qbox[3*i]);
		getboxcoords(idc, &cbox[3*i]);
		//printf("Point %d is at grid box: %d, %d, %d\n",i+1, qbox[3*i],qbox[3*i+1],qbox[3*i+2]);
		//count number of points in each box
		qgridcount[qbox[3*i]][qbox[3*i+1]][qbox[3*i+2]]++;
		cgridcount[cbox[3*i]][cbox[3*i+1]][cbox[3*i+2]]++;
	}
	
	
	//these will have the coordinates of each point in a specific box
	int *qpointsinbox[numboxes][3];
	int *cpointsinbox[numboxes][3];
	
	for (int i=0; i < numboxes; ++i){
		//get box number, coordinates from boxid. boxid=i
		int coord[3];
		getboxcoords(i, &coord[0]);
		//printf("\ndef id = %d \n%d,%d,%d\n",i,coord[0],coord[1],coord[2]);
		//if( getboxid(coord[0],coord[1],coord[2]) != i) printf("MATHERRORHERE");
		for (int j=0; j < 3; ++j){
			if ((qpointsinbox[i][j] = malloc(qgridcount[coord[0]][coord[1]][coord[2]] * sizeof *qpointsinbox[i][j])) == NULL) {
				perror("malloc 3");
				return 1;
			}
			if ((cpointsinbox[i][j] = malloc(cgridcount[coord[0]][coord[1]][coord[2]] * sizeof *cpointsinbox[i][j])) == NULL) {
				perror("malloc 3");
				return 1;
			}
		}
	}
	
	//put points in boxes, and keep the data in different tables - this will be useful for passing boxes around
	int *ccountforboxes = malloc(numboxes * sizeof(int));//keep tabs
	int *qcountforboxes = malloc(numboxes * sizeof(int));
	for (int i=0;i<numboxes;i++) {qcountforboxes[i]=0; ccountforboxes[i]=0; }
	
	for (int i=0;i<numberofpoints/processes;i++){
		int qtempid=getboxid(qbox[3*i],qbox[3*i+1],qbox[3*i+2]);
		int ctempid=getboxid(cbox[3*i],cbox[3*i+1],cbox[3*i+2]);
		for (int j=0;j<3;j++){
			qpointsinbox[qtempid][j][qcountforboxes[qtempid]]=q[3*i+j];
			cpointsinbox[ctempid][j][ccountforboxes[ctempid]]=q[3*i+j];	
		}
		qcountforboxes[qtempid]++;
		ccountforboxes[ctempid]++;
		//if (qcountforboxes[qtempid]>qgridcount[qbox[3*i]][qbox[3*i+1]][qbox[3*i+2]]) printf("mathfuckup\n");
		//if (ccountforboxes[ctempid]>cgridcount[cbox[3*i]][cbox[3*i+1]][cbox[3*i+2]]) printf("mathfuckup2\n");
	}
	
	
	
	//data tests - can be ignored
	/*for (int i=0;i<numberofpoints/processes;i++){
		
		qcountforboxes[0]++;
		ccountforboxes[0]++;

	}*/
	
	/*
	int qtemp=0;
	int qtemp2=0;
	int ctemp=0;
	int ctemp2=0;
	int ct[3];
	for (int i=0; i < numboxes; ++i){
		getboxcoords(i, &ct[0]);
		qtemp+=(int)qgridcount[ct[0]][ct[1]][ct[2]];
		ctemp+=(int)cgridcount[ct[0]][ct[1]][ct[2]];
		int tempid=getboxid(ct[0],ct[1],ct[2]);
		qtemp2+=qcountforboxes[tempid];
		ctemp2+=ccountforboxes[tempid];
		printf("Box %d,%d,%d \nID: %d == %d\nQ: %d %d\nC: %d %d\n\n",ct[0]+1,ct[1]+1,ct[2]+1,tempid, i,(int)qgridcount[ct[0]][ct[1]][ct[2]],qcountforboxes[tempid],(int)cgridcount[ct[0]][ct[1]][ct[2]],ccountforboxes[tempid]);
	}
	printf("Qtot: %d\nCtot: %d\n",qtemp,ctemp);
	printf("Qtot: %d\nCtot: %d\n",qtemp2,ctemp2);
	printf("%d %d\n",numberofpoints/processes,numboxes);
*/
	return 0;
}

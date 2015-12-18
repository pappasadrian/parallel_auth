#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

//global variables
int boxdimensions[3];
int splitdimensions[3];
int boxespersplit[3];
int Pprocesses;
int processes;
int Pnumberofpoints;
int numberofpoints;
int Pnumberboxes;
int numboxes;
int Pnumboxesperprocess;
int numboxesperprocess;

int processid=3; //for testing purposes - this should be dynamically allocated

//get split coordinates from process id
void getsplitcoords(int id, int *coords){
	coords[0]=id%splitdimensions[0];
	coords[1]=(id/splitdimensions[0])%splitdimensions[1];
	coords[2]=((id/splitdimensions[0])/splitdimensions[1])%splitdimensions[2];
}

//get process id from split coordinates
int getprocessid(int x, int y, int z){
	int i= x+y*(splitdimensions[0])+z*splitdimensions[0]*splitdimensions[1];
	return i;
}

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

//find which process's is the given box (by id)
int whosebox(int boxid){
	int boxcoords[3];
	getboxcoords(boxid,&boxcoords[0]);
	int splitcoords[3];
	for (int i=0;i<3;i++) splitcoords[i]=boxcoords[i]/boxespersplit[i];
	int pro;
	pro = getprocessid(splitcoords[0],splitcoords[1],splitcoords[2]);
	return pro;
}	

//check if given boxid belongs to current process
int ismybox(int id){
	if (whosebox(id)==processid) return 1;
	else return 0;
}

//get the id of a neigboring box. direction can be from 1 to 6, think the boxes
//this disregards diagonally adjacent boxes. this could be changed to 14 or 26 adjacent boxes
//must be discussed
int getneigborid(int id, int dir){
	int boxcoords[3];
	getboxcoords(id,&boxcoords[0]);
	switch (dir){
		case 1:
			if (boxcoords[0]<boxdimensions[0]) boxcoords[0]++;
			else return -1;
			break;
		case 2:
			if (boxcoords[1]<boxdimensions[1]) boxcoords[1]++;
			else return -1;
			break;
		case 3:
			if (boxcoords[2]<boxdimensions[2]) boxcoords[2]++;
			else return -1;
			break;
		case 4:
			if (boxcoords[0]>0)boxcoords[0]--;
			else return -1;
			break;
		case 5:
			if (boxcoords[1]>0)boxcoords[1]--;
			else return -1;
			break;
		case 6: 
			if (boxcoords[2]>0)boxcoords[2]--;
			else return -1;
			break;
		default:
			return -1;
		} 
	int temp;
	temp=getboxid(boxcoords[0],boxcoords[1],boxcoords[2]);
	return temp;
}

//check if given box is adjacent to a specified process's split. will be used for C. 
//there must be some smarter, more math-based way to perform this rather than brute forcing it
//but fuck it
void getadjacentboxesofaprocess(int pid, int *nb){	
	int count=0;
	for (int i=0;i<numboxes;i++){	
		if(whosebox(i)==pid) {
			for (int j=1;j<7;j++){
				int testneigbor = getneigborid(i,j);
				if(testneigbor!=-1&&ismybox(testneigbor)==0) 
				{
					nb[count]=testneigbor;
					count++;
					//printf("%d\n",testneigbor);
				}
			}
		}
	}
	nb[count]=-1; //shows the end of the neigbors.
}

//find if a given box is adjacent to another split, and return that split's processid
//else (if it has no neigbor split, return -1
int whoseneigboristhis(int thisid){
	for (int j=1;j<7;j++){
		int testneigbor = getneigborid(thisid,j);
		if(testneigbor!=-1&&ismybox(testneigbor)==0) 
		{
			int temp = whosebox(testneigbor);
			return temp;
		}
	}
	return -1;
}

//evaluate the euclidean distance between two points
float euclidean(float x1,float y1,float z1,float x2,float y2,float z2){
	float d = sqrtf((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	return d;
}

void checkqpointsinbox(float **c, float **q, int clength, int qlength){
	for (int i=0;i<qlength;i++) i++;
}

//check if given box is closer than the candidate distance of q from candidate c
//this is done by using this generic code, that sees if a sphere intersects a cuboid.
//returns true if it intersect, false if it doesnt.
float squared(float v) { return v * v; }
//bool doescubeintersectsphere(vec3 C1, vec3 C2, vec3 S, float R)
int doescubeintersectsphere(int boxid, float qx, float qy, float qz, float R)
{
	float distancesquared = R * R;
	/* assume C1 and C2 are element-wise sorted, if not, do that now */
	int boxcoords[3];
	getboxcoords(boxid, &boxcoords[0]);
	float xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=(float)boxcoords[0]/(float)boxdimensions[0];
	xmax=(float)(boxcoords[0]+1)/(float)boxdimensions[0];
	ymin=(float)boxcoords[1]/(float)boxdimensions[1];
	ymax=(float)(boxcoords[1]+1)/(float)boxdimensions[1];
	zmin=(float)boxcoords[2]/(float)boxdimensions[2];
	zmax=(float)(boxcoords[2]+1)/(float)boxdimensions[2];
	
	//copy paste from here on
	//lets trust the copy paste
	if (qx < xmin) distancesquared -= squared(qx - xmin);
	else if (qx > xmax) distancesquared -= squared(qx - xmax);
	if (qy < ymin) distancesquared -= squared(qy - ymin);
	else if (qy > ymax) distancesquared -= squared(qy - ymax);
	if (qz < zmin) distancesquared -= squared(qz - zmin);
	else if (qz > zmax) distancesquared -= squared(qz - zmax);
	return distancesquared > 0;
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

//main function
int main(int argc, char **argv){

	if (argc != 4) {
		printf("Usage: %s Q P B\n",argv[0]);
		printf("Q: number of points (power of two - from 20 to 25)\n");
		printf("P: number of processes (power of two - from 0 to 7)\n");
		printf("B: number of boxes (power of two - from 12 to 16)\n");
		return 1;
	}
  
	srand (time(NULL));  //such randomness wow
	
	Pprocesses=atoi(argv[2]); //from 0 to 7
	processes=1<<Pprocesses;
	Pnumberofpoints=atoi(argv[1]); // from 0(?) to 25 -in all processes. this process has numberofpoints/processes points.
	numberofpoints=1<<Pnumberofpoints;
	Pnumberboxes=atoi(argv[3]); //from 12 to 16
	numboxes=1<<Pnumberboxes;
	Pnumboxesperprocess=Pnumberboxes-Pprocesses; //2^x number of grid boxes per process, aka splits
	numboxesperprocess=1<<Pnumboxesperprocess;
	
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
	//printf("grid dimensions: %d, %d, %d\n",boxdimensions[0],boxdimensions[1],boxdimensions[2]);
	/*
	if (boxdimensions[0]*boxdimensions[1]*boxdimensions[2]!=numboxes) {//fuckup in the math
		printf("mathfuckup\n"); 
		return 1;
	} 
	*/
	
	//same idea, but with splits. each split contains many boxes and belongs to a process
	splitdimensions[0]=1<<(Pprocesses/3);
	if (Pprocesses%3 > 0)	splitdimensions[0]*=2; //double the boxes in this row if mod3 is 1 or 2
	splitdimensions[1]=1<<(Pprocesses/3);
	if (Pprocesses%3 == 2)	splitdimensions[1]*=2; //double the boxes in this row if mod3 is 2
	splitdimensions[2]=1<<(Pprocesses/3);
		
	for(int i=0;i<3;i++) boxespersplit[i]=boxdimensions[i]/splitdimensions[i];
	
	//printf("grid split size: %d, %d, %d\n",boxespersplit[0],boxespersplit[1],boxespersplit[2]);
	
	//these could be int, but then i'd need another alloc_3d for int type - maybe ill do it later	
	//these will keep tabs on how many points are in each grid box
	//this could be done with box IDs and not coordinates, and it wouldnt require the 3d allocation, it was written before the convention
	//im too bored to fix it
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
	float *qpointsinbox[numboxes][3];
	float *cpointsinbox[numboxes][3];
	
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
	for (int i=0;i<numboxes;i++) {qcountforboxes[i]=0; ccountforboxes[i]=0; }//init with 0s
	
	for (int i=0;i<numberofpoints/processes;i++){//for each point in this process
		int qtempid=getboxid(qbox[3*i],qbox[3*i+1],qbox[3*i+2]);
		int ctempid=getboxid(cbox[3*i],cbox[3*i+1],cbox[3*i+2]);
		for (int j=0;j<3;j++){
			qpointsinbox[qtempid][j][qcountforboxes[qtempid]]=q[3*i+j];
			cpointsinbox[ctempid][j][ccountforboxes[ctempid]]=c[3*i+j];	
		}
		qcountforboxes[qtempid]++;
		ccountforboxes[ctempid]++;
		//if (qcountforboxes[qtempid]>qgridcount[qbox[3*i]][qbox[3*i+1]][qbox[3*i+2]]) printf("mathfuckup\n");
		//if (ccountforboxes[ctempid]>cgridcount[cbox[3*i]][cbox[3*i+1]][cbox[3*i+2]]) printf("mathfuckup2\n");
	}
	
	//find neinbors of boxes of this split, in order to keep them for C. 
	int maxneigbors;
	maxneigbors=(boxespersplit[0]*boxespersplit[1]+boxespersplit[0]*boxespersplit[2]+boxespersplit[1]*boxespersplit[2])*2;
	int *neigbors = malloc(maxneigbors * sizeof(int));
	getadjacentboxesofaprocess(processid, &neigbors[0]);
	
	//from here on, the data passing must commence.
	//each process will keep the boxes that it owns, and pass the other ones to the respective processes.
	//all processes must know which box goes to which process.
	for (int i=0;i<numboxes;i++){
		int tempid=whosebox(i);
		if (tempid!=processid){
			//passbox(**qpointsinbox[i],toprocess(processid))
			//same for C
		}
	}
	
	//see whose neigbor is each box of this process
	//so, if we're searching in this box, we have to request from the neigbor process to pass the box to us.
	// i failed to pass my awkward 3d matrices to a function, so i did this monstrocity...
	for (int i=0;i<numboxes;i++){
		if (ismybox(i)){
			//printf("%d q points\n",ccountforboxes[i]);
			for (int j=0;j<qcountforboxes[i];j++){
				float qcoordtemp[3];
				float ccandidate[3];
				float cfinal[3];
				qcoordtemp[0]=qpointsinbox[i][0][j];
				qcoordtemp[1]=qpointsinbox[i][1][j];
				qcoordtemp[2]=qpointsinbox[i][2][j];
				
				float bestdistance;
				int tempcandidate;
				float tempdistance;
				//printf("%d c points\n",ccountforboxes[j]);
				for (int cp=0;cp<ccountforboxes[j];cp++){
					tempdistance = euclidean(qcoordtemp[0],qcoordtemp[1],qcoordtemp[2],cpointsinbox[i][0][cp],cpointsinbox[i][1][cp],cpointsinbox[i][2][cp]);
					if (tempdistance<bestdistance||cp==0) {
						//printf("-|-");
						bestdistance=tempdistance;
						tempcandidate=cp;
					}
					
				}
				cfinal[0]=cpointsinbox[i][0][tempcandidate];
				cfinal[1]=cpointsinbox[i][1][tempcandidate];
				cfinal[2]=cpointsinbox[i][2][tempcandidate];
				
				for (int dir=1;dir<7;dir++){
					int tempid=getneigborid(i, dir);//new c box to search at
					if (tempid!=-1){//if there is a neigbor
						if (ismybox(tempid)){
							//printf("LOOKINGATNEIGBOR\n");
							for (int cp=0;cp<ccountforboxes[tempid];cp++){
								tempdistance = euclidean(qcoordtemp[0],qcoordtemp[1],qcoordtemp[2],cpointsinbox[tempid][0][cp],cpointsinbox[tempid][1][cp],cpointsinbox[tempid][2][cp]);
								if (tempdistance<bestdistance) {
									//printf("BETTERATNEIGBOR\n");
									bestdistance=tempdistance;
									tempcandidate=cp;
									cfinal[0]=cpointsinbox[i][0][tempcandidate];
									cfinal[1]=cpointsinbox[i][1][tempcandidate];
									cfinal[2]=cpointsinbox[i][2][tempcandidate];
								}
							}
						}
						else{//must look in neigbor process
							//lets ignore this for now
							//printf("could have been at neigbor process\n");
						}
					}
					
				}	
			//printf("Point Q at coords %f,%f,%f is nearest to point C at coords %f,%f,%f\n",qcoordtemp[0],qcoordtemp[1],qcoordtemp[2],cfinal[0],cfinal[1],cfinal[2]);
			}
		//printf("\nNEXTBOX\n");
		}
	}
	
	//NOTE THAT THE FOLLOWING ARE NOT IN ACCORDANCE TO THE EKFONISI
	//suggestion/hack/slacking around: each process should process its boxes, meaning that it will check the Q points in its boxes
	//however, it is a good idea that it keeps C points that are not only in its boxes, but also to the boxes adjacent to it. 
	//like that, when it is searching at the neighbour boxes of a box at an edge, it wont have to mpi call another process to get the data back
	//this will only be done once at the beginning.
	//ie each process will keep and receive the Q points in its boxes, and the C points in its boxes and adjacent boxes.
	//this is a hack but it will save us from a lot of mpi calls during search, which might lead to many locks that will drop performance 
	
	//another hack would be to completely disregard searches of nearest neigbor in different boxes outside our split
	//but it is bad, it will lead to mistakes at the borders.
	//however, it will make message passing much much easier and ill use it if all else fails.
	
	//note that there is a shitload of messages to be passed around during point search, so this might actually be a good idea...
	
	//data tests - can be ignored
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
	/*see points in a certain box
	int testbox=2000;
	for (int i=0;i<qcountforboxes[testbox];i++){
		printf("\nPoint %d: ",i);
		for (int j=0;j<3;j++)
			printf("%f,",qpointsinbox[testbox][j][i]);
	}
	*/
	//which box is in each split
	//for (int i=0;i<numboxes;i++)	printf("Box %d is in split %d\n",i,whosebox(i)+1);
	//find my boxes
	//for (int i=0;i<numboxes;i++)	if(ismybox(i)==1)	printf("Box %d is mine :D\n",i);
	
	//test for box with sphere intersections
	int count=0;
	for (int i=0;i<numboxes;i++){
		float testx=0.1;
		float testy=0.5;
		float testz=0.4;
		float testdist=0.04;
		if (findinwhichbox(testx,testy,testz)!=i){//if this is not Q's box
			if (doescubeintersectsphere(i, testx, testy, testz, testdist)){
				//needed just for print purposes
				int boxcoords[3];
				getboxcoords(i, &boxcoords[0]);
				float xmin, xmax, ymin, ymax, zmin, zmax;
				xmin=(float)boxcoords[0]/(float)boxdimensions[0];
				xmax=(float)(boxcoords[0]+1)/(float)boxdimensions[0];
				ymin=(float)boxcoords[1]/(float)boxdimensions[1];
				ymax=(float)(boxcoords[1]+1)/(float)boxdimensions[1];
				zmin=(float)boxcoords[2]/(float)boxdimensions[2];
				zmax=(float)(boxcoords[2]+1)/(float)boxdimensions[2];
				printf("Cube %d is within search range\n",i);
				printf("X from %f to %f\n",xmin,xmax);
				printf("Y from %f to %f\n",ymin,ymax);
				printf("Z from %f to %f\n\n",zmin,zmax);
				count++;
				//needed just for print purposes
			}
			
		}
		else{//else, we are in the same box as Q - probably do nothing
			//needed just for print purposes
			int boxcoords[3];
			getboxcoords(i, &boxcoords[0]);
			float xmin, xmax, ymin, ymax, zmin, zmax;
			xmin=(float)boxcoords[0]/(float)boxdimensions[0];
			xmax=(float)(boxcoords[0]+1)/(float)boxdimensions[0];
			ymin=(float)boxcoords[1]/(float)boxdimensions[1];
			ymax=(float)(boxcoords[1]+1)/(float)boxdimensions[1];
			zmin=(float)boxcoords[2]/(float)boxdimensions[2];
			zmax=(float)(boxcoords[2]+1)/(float)boxdimensions[2];
			printf("Cube %d is container of Q\n",i);
			printf("X from %f to %f\n",xmin,xmax);
			printf("Y from %f to %f\n",ymin,ymax);
			printf("Z from %f to %f\n\n",zmin,zmax);
			//needed just for print purposes
			
		}
	}
	printf("found %d neighbor boxes\n",count);
	
}


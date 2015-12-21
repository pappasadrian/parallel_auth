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

//struct for point coordinates
struct pointcoords {
   float x;
   float y;
   float z;
};

//struct for box coordinates
//can be used for splits too
struct boxcoords {
   int x;
   int y;
   int z;
};

//structs to be passed around via mpi
struct messagebox {
	int numpoints;
	int boxid;
	struct pointcoords *point;
};

int processid=6; //for testing purposes - this should be dynamically allocated

//get split coordinates from process id
struct boxcoords procID_2_split(int id){
	struct boxcoords coords;
	coords.x=id%splitdimensions[0];
	coords.y=(id/splitdimensions[0])%splitdimensions[1];
	coords.z=((id/splitdimensions[0])/splitdimensions[1])%splitdimensions[2];
	return coords;
}

//get process id from split coordinates
int split_2_procID(struct boxcoords a){
	int i= a.x+a.y*(splitdimensions[0])+a.z*splitdimensions[0]*splitdimensions[1];
	return i;
}

//get box coordinates from box id
struct boxcoords get_box_coords(int id){
	struct boxcoords out;
	out.x=id%boxdimensions[0];
	out.y=(id/boxdimensions[0])%boxdimensions[1];
	out.z=((id/boxdimensions[0])/boxdimensions[1])%boxdimensions[2];
	return out;
}

//get box id from coordinates
int get_box_id(struct boxcoords box){
	int i= box.x+box.y*(boxdimensions[0])+box.z*boxdimensions[0]*boxdimensions[1];
	return i;
}

//find in which grid box the point at the given coordinates is at. returns box id
int find_in_which_box(struct pointcoords a){ 
	struct boxcoords box;
	box.x = (int)(a.x*boxdimensions[0]); //rounds down to the nearest int
	box.y = (int)(a.y*boxdimensions[1]);
	box.z = (int)(a.z*boxdimensions[2]);
	return get_box_id(box);
}

//find which process's is the given box (by id)
int get_box_owner(int boxid){
	struct boxcoords boxc;
	boxc=get_box_coords(boxid);
	struct boxcoords splitcoords;
	splitcoords.x=boxc.x/boxespersplit[0];
	splitcoords.y=boxc.y/boxespersplit[1];
	splitcoords.z=boxc.z/boxespersplit[2];
	int pro;
	pro = split_2_procID(splitcoords);
	return pro;
}	

//check if given boxid belongs to current process
int is_my_box(int id){
	if (get_box_owner(id)==processid) return 1;
	else return 0;
}


//new getneighborid, which can return any of the 26 adjacent boxes.
int get_neighbor_boxID_coords(int id, int dirx, int diry, int dirz){
	struct boxcoords boxc;
	boxc=get_box_coords(id);
	switch (dirx){
		case 1:
			if (boxc.x<boxdimensions[0]) boxc.x++;
			else return -1;//out of bounds
			break;
		case 0:
			break;
		case -1:
			if (boxc.x>0) boxc.x--;
			else return -1;//out of bounds
			break;
		default://input not correct
			return -1;
		} 
	switch (diry){
		case 1:
			if (boxc.y<boxdimensions[1]) boxc.y++;
			else return -1;//out of bounds
			break;
		case 0:
			break;
		case -1:
			if (boxc.y>0) boxc.y--;
			else return -1;//out of bounds
			break;
		default://input not correct
			return -1;
		}
	switch (dirz){
		case 1:
			if (boxc.z<boxdimensions[2])boxc.z++;
			else return -1;//out of bounds
			break;
		case 0:
			break;
		case -1:
			if (boxc.z>0) boxc.z--;
			else return -1;//out of bounds
			break;
		default://input not correct
			return -1;
		}
	int temp;
	temp=get_box_id(boxc);
	//printf("%d ",temp);
	return temp;
}

//helper function for using getneigborid(int id, int dirx, int diry, int dirz)
//in for loops from 0 to 26, for checking all neigbors
//one of the iterations of the loop will return Q's box. in that case, we return -2 for reference
int get_neighbor_id(int loopnum, int boxid){
	int x, y, z;
	int neighbor;
	//loopnum++;
	x=(loopnum%3) - 1;
	z=(loopnum/9)%3-1;
	y=(loopnum/3)%3-1;
	//printf("%d: %d,%d,%d\n",loopnum,x,y,z);
	//printf("%d,%d,%d\n",x,y,z);
	if (x==0&&y==0&&z==0){
		return -2;
	}
	else{
		neighbor=get_neighbor_boxID_coords(boxid,x,y,z);
		return neighbor;
	}
}


//check if given box is adjacent to a specified process's split. will be used for C. 
//there must be some smarter, more math-based way to perform this rather than brute forcing it
//but fuck it
//this is not working and need to be checked
/*
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
*/

//evaluate the euclidean distance between two points
float euclidean(struct pointcoords a, struct pointcoords b){
	float d = sqrtf((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)+(b.z-a.z)*(b.z-a.z));
	return d;
}


//check if given box is closer than the candidate distance of q from candidate c
//this is done by using this generic code, that sees if a sphere intersects a cuboid.
//returns true if it intersect, false if it doesnt.
float squared(float v) { return v * v; }
//bool doescubeintersectsphere(vec3 C1, vec3 C2, vec3 S, float R)
int does_box_intersect_sphere(int boxid, struct pointcoords q, float R)
{
	if (boxid==-1) return 0;
	float distancesquared = R * R;
	/* assume C1 and C2 are element-wise sorted, if not, do that now */
	struct boxcoords boxc;
	boxc=get_box_coords(boxid);
	float xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=(float)boxc.x/(float)boxdimensions[0];
	xmax=(float)(boxc.x+1)/(float)boxdimensions[0];
	ymin=(float)boxc.y/(float)boxdimensions[1];
	ymax=(float)(boxc.y+1)/(float)boxdimensions[1];
	zmin=(float)boxc.z/(float)boxdimensions[2];
	zmax=(float)(boxc.z+1)/(float)boxdimensions[2];
	
	//copy paste from here on
	//lets trust the copy paste
	if (q.x < xmin) distancesquared -= squared(q.x - xmin);
	else if (q.x > xmax) distancesquared -= squared(q.x - xmax);
	if (q.y < ymin) distancesquared -= squared(q.y - ymin);
	else if (q.y > ymax) distancesquared -= squared(q.y - ymax);
	if (q.z < zmin) distancesquared -= squared(q.z - zmin);
	else if (q.z > zmax) distancesquared -= squared(q.z - zmax);
	return distancesquared > 0;
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
	struct pointcoords *q = malloc(sizeof(struct pointcoords) * (numberofpoints / processes)); //numberofpoints/processes, because the points are spread through the active processes
	struct pointcoords *c = malloc(sizeof(struct pointcoords) * (numberofpoints / processes));
	
	for (int i=0;i<numberofpoints/processes;i++){
		q[i].x=(float)rand() / RAND_MAX; //random value from 0 to 1
		q[i].y=(float)rand() / RAND_MAX;
		q[i].z=(float)rand() / RAND_MAX;
		c[i].x=(float)rand() / RAND_MAX;
		c[i].y=(float)rand() / RAND_MAX;
		c[i].z=(float)rand() / RAND_MAX;
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
	
	//these will keep tabs on how many points are in each grid box
	int *qgridcount=malloc(numboxes*sizeof(int)); //how many points in each grid box
	int *cgridcount=malloc(numboxes*sizeof(int)); //how many points in each grid box
	for (int i=0;i<numboxes;i++){
		qgridcount[i]=cgridcount[i]=0;
	}
			
	//find where each point belongs
	struct boxcoords *qbox = malloc(sizeof(struct boxcoords) * (numberofpoints / processes));
	struct boxcoords *cbox = malloc(sizeof(struct boxcoords) * (numberofpoints / processes));
	//and count how many q and c points are in this split
	int qcountinthissplit=0;
	int ccountinthissplit=0;
	for (int i=0;i<numberofpoints/processes;i++){
		int idq = find_in_which_box(q[i]);
		int idc = find_in_which_box(c[i]);
		qbox[i]=get_box_coords(idq);
		cbox[i]=get_box_coords(idc);
		//printf("Point %d is at grid box: %d, %d, %d\n",i+1, qbox[3*i],qbox[3*i+1],qbox[3*i+2]);
		//count number of points in each box
		qgridcount[idq]++;
		cgridcount[idc]++;
		if (is_my_box(idq)) qcountinthissplit++;
		if (is_my_box(idc)) ccountinthissplit++;
	}
	
	//these will have the coordinates of each point in a specific box
	struct pointcoords *qpointsinbox[numboxes];
	struct pointcoords *cpointsinbox[numboxes];
	
	for (int i=0; i < numboxes; ++i){
		//get box number, coordinates from boxid. boxid=i
		struct boxcoords coord;
		coord=get_box_coords(i);
		//printf("\ndef id = %d \n%d,%d,%d\n",i,coord[0],coord[1],coord[2]);
		//if( getboxid(coord[0],coord[1],coord[2]) != i) printf("MATHERRORHERE");
		if ((qpointsinbox[i] = malloc(qgridcount[i] * sizeof *qpointsinbox[i])) == NULL) {
			perror("malloc 3");
			return 1;
		}
		if ((cpointsinbox[i] = malloc(cgridcount[i] * sizeof *cpointsinbox[i])) == NULL) {
			perror("malloc 3");
			return 1;
		}
	}
	
	//put points in boxes, and keep the data in different tables - this will be useful for passing boxes around
	int *ccountforboxes = malloc(numboxes * sizeof(int));//keep tabs
	int *qcountforboxes = malloc(numboxes * sizeof(int));
	for (int i=0;i<numboxes;i++) {qcountforboxes[i]=0; ccountforboxes[i]=0; }//init with 0s
	
	for (int i=0;i<numberofpoints/processes;i++){//for each point in this process
		int qtempid=get_box_id(qbox[i]);
		int ctempid=get_box_id(cbox[i]);
		qpointsinbox[qtempid][qcountforboxes[qtempid]]=q[i];
		cpointsinbox[ctempid][ccountforboxes[ctempid]]=c[i];
		qcountforboxes[qtempid]++;
		ccountforboxes[ctempid]++;
		//if (qcountforboxes[qtempid]>qgridcount[qbox[3*i]][qbox[3*i+1]][qbox[3*i+2]]) printf("mathfuckup\n");
		//if (ccountforboxes[ctempid]>cgridcount[cbox[3*i]][cbox[3*i+1]][cbox[3*i+2]]) printf("mathfuckup2\n");
	}
	
	
	//from here on, the data passing must commence.
	//each process will keep the boxes that it owns, and pass the other ones to the respective processes.
	//all processes must know which box goes to which process.
	
	//here, we create the message structs, each containing a box to be sent to another process.
	int numboxesnotmine=numboxes-numboxesperprocess;
	struct messagebox *qmessage = malloc (sizeof(struct messagebox) * (numboxesnotmine));
	struct messagebox *cmessage = malloc (sizeof(struct messagebox) * (numboxesnotmine));
	int boxesnotminecounter=0;
	for (int i=0;i<numboxes;i++){
		int tempid=get_box_owner(i);
		if (!(is_my_box(tempid))){
			qmessage[boxesnotminecounter].numpoints=qgridcount[i];
			qmessage[boxesnotminecounter].boxid=i;
			qmessage[boxesnotminecounter].point = malloc(qmessage[boxesnotminecounter].numpoints * sizeof *qmessage[boxesnotminecounter].point);  
			for (int j=0;j<qmessage[boxesnotminecounter].numpoints;j++){
				qmessage[boxesnotminecounter].point[j]=qpointsinbox[i][j];//fix this
			}
			cmessage[boxesnotminecounter].numpoints=cgridcount[i];
			cmessage[boxesnotminecounter].boxid=i;
			cmessage[boxesnotminecounter].point = malloc(cmessage[boxesnotminecounter].numpoints * sizeof *cmessage[boxesnotminecounter].point);  
			for (int j=0;j<cmessage[boxesnotminecounter].numpoints;j++){
				cmessage[boxesnotminecounter].point[j]=cpointsinbox[i][j];//fix this
			}
			boxesnotminecounter++;	
		}
	}
	
	//prepare an array of box structs to be sent to a specific process
	//mpi should send the two 1d arrays qmessagetoprocess and cmessagetoprocess
	//to the processidtosendto process.
	int processidtosendto=4;//this should be chosen dynamically somehow
	struct messagebox *qmessagetoprocess = malloc (sizeof(struct messagebox) * (numboxesperprocess));
	struct messagebox *cmessagetoprocess = malloc (sizeof(struct messagebox) * (numboxesperprocess));
	int ccounterformessagetoprocess=0;
	int qcounterformessagetoprocess=0;
	for (int i=0;i<numboxesnotmine;i++){
		if(get_box_owner(qmessage[i].boxid)==processidtosendto){
			qmessagetoprocess[qcounterformessagetoprocess]=qmessage[i];
			qcounterformessagetoprocess++;
		}
		if(get_box_owner(cmessage[i].boxid)==processidtosendto){
			cmessagetoprocess[ccounterformessagetoprocess]=cmessage[i];
			ccounterformessagetoprocess++;
		}
	}
	
	//receive an array of box structs
	//temporarily use the qmessagetoprocess table from above
	receiveddata=1;
	if(receiveddata){//do this when we receive an array of data
		for (int i=0;i<numboxesperprocess;i++){
			qmessagetoprocess[i];
			cmessagetoprocess[i];
			qcountforboxes[qmessagetoprocess[i].boxid]+=qmessagetoprocess[i].numboxes;
			//realloc qpointsinbox[qmessagetoprocess[i].boxid];
			//concatenate qpointsinbox with qcountforboxes[i].pointcoords
			ccountforboxes[cmessagetoprocess[i].boxid]+=cmessagetoprocess[i].numboxes;
			//realloc cpointsinbox[cmessagetoprocess[i].boxid];
			//concatenate cpointsinbox with ccountforboxes[i].pointcoords
		}
	}
	
	/* USE THIS FOR CONCATENATION OF ARRAYS - TEST CODE
	float x[4] = { 1, 1, 1, 1 };
	float y[4] = { 2, 2, 2, 2 };

	float* total = malloc(8 * sizeof(float)); // array to hold the result

	memcpy(total,     x, 4 * sizeof(float)); // copy 4 floats from x to total[0]...total[3]
	memcpy(total + 4, y, 4 * sizeof(float)); // copy 4 floats from y to total[4]...total[7]
	
	*/
	
	
	//SEARCH PART OF THE PROGRAM
	//array to gather the results
	struct pointcoords *results = malloc(sizeof(struct pointcoords) * (2 * qcountinthissplit));
	int pointertoresults=0;
	//see whose neigbor is each box of this process
	//so, if we're searching in this box, we have to request from the neigbor process to pass the box to us.
	// i failed to pass my awkward 3d matrices to a function, so i did this monstrocity...
	for (int i=0;i<numboxes;i++){
		if (is_my_box(i)){
			//printf("%d q points\n",ccountforboxes[i]);
			for (int j=0;j<qcountforboxes[i];j++){
				struct pointcoords qcoordtemp;
				struct pointcoords cfinal;
				qcoordtemp=qpointsinbox[i][j];
				
				float bestdistance=-1;
				int tempcandidate;
				float tempdistance;
				//printf("%d c points\n",ccountforboxes[j]);
				for (int cp=0;cp<ccountforboxes[j];cp++){
					tempdistance = euclidean(qcoordtemp,cpointsinbox[i][cp]);
					if (tempdistance<bestdistance||cp==0) {
						//printf("-|-");
						bestdistance=tempdistance;
						tempcandidate=cp;
					}
				}
				//just in case there are no C points in this box, use a very large bestdistance, so that the following searches will work
				if (bestdistance==-1){
					bestdistance=1;
				}
				else{
					cfinal=cpointsinbox[i][tempcandidate];
				}
				int checkedneighborcounter=0;
				int inneighborprocessescounter=0;
				for (int dir=0;dir<27;dir++){
					int tempid=get_neighbor_id(dir, i);
					int shouldwecheckneighbor=0;
					if (bestdistance<1){
						shouldwecheckneighbor=does_box_intersect_sphere(tempid, cfinal, bestdistance);
						//printf("%d",shouldwecheckneighbor);
					}
					else shouldwecheckneighbor=1;
					if (tempid>-1&&shouldwecheckneighbor){//if there is a neigbor and the box in question is not Q's box and the box in question is within range of the sphere
						checkedneighborcounter++;
						if (is_my_box(tempid)){//this is happening at a box of the same process
							//printf("LOOKING IN THIS PROCESS");
							for (int cp=0;cp<ccountforboxes[tempid];cp++){
								tempdistance = euclidean(qcoordtemp,cpointsinbox[tempid][cp]);
								if (tempdistance<bestdistance) {
									//printf("BETTER AT NEIGBOR %d\n",dir);
									bestdistance=tempdistance;
									tempcandidate=cp;
									cfinal=cpointsinbox[i][tempcandidate];
								}
							}
						}
						else{//must look in neigbor process
							//lets ignore this for now
							inneighborprocessescounter++;
							//printf("could have been at neigbor process\n");
						}
					}
				}
			printf("Point Q at coords %f,%f,%f is nearest to point C at coords %f,%f,%f\n%d neighbor box%s of the same process checked for this result.\n%d candidate box%s of neighbor processes should have been checked.\n\n",qcoordtemp.x,qcoordtemp.y,qcoordtemp.z,cfinal.x,cfinal.y,cfinal.z,checkedneighborcounter,(checkedneighborcounter!=1)?"es":"",inneighborprocessescounter,(inneighborprocessescounter!=1)?"es":"");
			results[pointertoresults]=qcoordtemp;
			results[pointertoresults+1]=cfinal;
			pointertoresults+=2;
			}
		//printf("\nNEXTBOX\n");
		}
	}
	
	//independent print of results
	//they are rounded to 3 decimal places, just for facilitating the view
	for (int i=0;i<qcountinthissplit;i++){
		//printf("#%d Q point at coords (%.3f,%.3f,%.3f) is nearest to point C at coords (%.3f,%.3f,%.3f)\n\n",i+1,results[2*i].x,results[2*i].y,results[2*i].z,results[2*i+1].x,results[2*i+1].y,results[2*i+1].z);
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

}


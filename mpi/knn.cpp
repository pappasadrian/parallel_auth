#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <string.h>
#include <iostream>

  using std::vector;
  using std::cout;
  using std::endl;
  
  
int processid=0; //for testing purposes - this should be dynamically allocated
//global variables
int boxdimensions[3];
int splitdimensions[3];
int boxespersplit[3];
int Pprocesses;
int processes;
//different amount of qpoints and cpoints
int PQnumberofpoints;
int qnumberofpoints;
int PCnumberofpoints;
int cnumberofpoints;
int Pnumberboxes;
int numboxes;
int Pnumboxesperprocess;
int numboxesperprocess;

//struct for point coordinates
struct Point{
  float x;
  float y;
  float z;
  Point(){};
  Point(float x, float y, float z);
};

Point::Point(float x, float y, float z){
  this->x=x;
  this->y=y;
  this->z=z;
}

struct QPoint : public Point{
  Point nn;
  bool found_nn;
  QPoint();
  ~QPoint(){};
  QPoint(float x, float y, float z);
};

QPoint::QPoint()
{
  found_nn=false;
}

QPoint::QPoint(float x, float y, float z) : Point( x,  y,  z)
{
  found_nn=false;
}

std::vector<Point> naive_search(QPoint q, std::vector<Point> search_space);

//struct for box coordinates
//can be used for splits too
struct boxcoords {
  boxcoords(){}
  int x;
  int y;
  int z;
};

template<class T> class Box{
public:
  boxcoords coords;
  int id;
  int owner;
  std::vector<T> point_cloud;
};

//structs to be passed around via mpi
struct messagebox {
  int numpoints;
  int boxid;
  std::vector<Point> point_cloud;
  messagebox(){};
};

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
int find_in_which_box(struct Point a){
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

//evaluate the euclidean distance between two points
float euclidean(struct Point a, struct Point b){
  float d = sqrtf((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)+(b.z-a.z)*(b.z-a.z));
  return d;
}

//simple printout of a point. just for testing reasons
void printcoords (struct Point point){
  cout<<"("<<point.x<<", "<<point.y<<", "<<point.z<<")";
}


//check if given box is closer than the candidate distance of q from candidate c
//this is done by using this generic code, that sees if a sphere intersects a cuboid.
//returns true if it intersect, false if it doesnt.
float squared(float v) { return v * v; }
//bool doescubeintersectsphere(vec3 C1, vec3 C2, vec3 S, float R)
int does_box_intersect_sphere(int boxid, struct Point q, float R)
{
  if (boxid==-1) return 0;
  float distancesquared = R * R;
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

  if (argc != 5) {
      printf("Usage: %s Q C P B\n",argv[0]);
      printf("Q: number of Q points (power of two - from 20 to 25)\n");
      printf("C: number of C points (power of two - from 20 to 25)\n");
      printf("P: number of processes (power of two - from 0 to 7)\n");
      printf("B: number of boxes (power of two - from 12 to 16)\n");
      return 1;
    }

  srand (time(NULL));  //such randomness wow
  Pprocesses=atoi(argv[3]); //from 0 to 7
  processes=1<<Pprocesses;
  PQnumberofpoints=atoi(argv[1]); // from 0(?) to 25 -in all processes. this process has numberofpoints/processes Q points.
  qnumberofpoints=1<<PQnumberofpoints;
  PCnumberofpoints=atoi(argv[1]); // from 0(?) to 25 -in all processes. this process has numberofpoints/processes C points.
  cnumberofpoints=1<<PCnumberofpoints;
  Pnumberboxes=atoi(argv[4]); //from 12 to 16
  numboxes=1<<Pnumberboxes;
  Pnumboxesperprocess=Pnumberboxes-Pprocesses; //2^x number of grid boxes per process, aka splits
  numboxesperprocess=1<<Pnumboxesperprocess;

  // coordinates will be in here. data example: {x1, y1, z1, x2, y2, z2, ...}
  //numberofpoints/processes, because the points are spread through the active processes
  int q_pts_per_process=(qnumberofpoints / processes);
  int c_pts_per_process=(cnumberofpoints / processes);
  //This will simulate our box grid, creating 2^pnumberboxes boxes.
  //Each boxdimension shows how many boxes are on each dimension of our cube.
  //n<=m<=k is kept.
  //generally, we can refer to these boxes either by ID or by coordinates. there are functions that do this.
  boxdimensions[0]=1<<(Pnumberboxes/3);
  if (Pnumberboxes%3 > 0)	boxdimensions[0]*=2; //double the boxes in this row if mod3 is 1 or 2
  boxdimensions[1]=1<<(Pnumberboxes/3);
  if (Pnumberboxes%3 == 2)	boxdimensions[1]*=2; //double the boxes in this row if mod3 is 2
  boxdimensions[2]=1<<(Pnumberboxes/3);

  //same idea, but with splits. each split contains many boxes and belongs to a process
  splitdimensions[0]=1<<(Pprocesses/3);
  if (Pprocesses%3 > 0)	splitdimensions[0]*=2; //double the boxes in this row if mod3 is 1 or 2
  splitdimensions[1]=1<<(Pprocesses/3);
  if (Pprocesses%3 == 2)	splitdimensions[1]*=2; //double the boxes in this row if mod3 is 2
  splitdimensions[2]=1<<(Pprocesses/3);
  for(int i=0;i<3;i++) boxespersplit[i]=boxdimensions[i]/splitdimensions[i];
  
  
  vector<Box<Point> > c_boxes;
  vector<Box<QPoint> > q_boxes;
  q_boxes.resize(numboxes);
  c_boxes.resize(numboxes);
  
  //initialize q_boxes and c_boxes
  for (uint i=0;i<q_boxes.size();i++){
      q_boxes[i].id=i;
      q_boxes[i].coords=get_box_coords(i);
      q_boxes[i].owner=get_box_owner(i);
    };
  for (uint i=0;i<c_boxes.size();i++){
      c_boxes[i].id=i;
      c_boxes[i].coords=get_box_coords(i);
      c_boxes[i].owner=get_box_owner(i);
    };

//generate q points
  for (int i=0;i<q_pts_per_process;i++){
      QPoint qtemp((float)rand() / RAND_MAX,(float)rand() / RAND_MAX,(float)rand() / RAND_MAX);
      q_boxes[find_in_which_box(qtemp)].point_cloud.push_back(qtemp);
    }
//generate c points    
  for (int i=0;i<c_pts_per_process;i++){
      Point ctemp((float)rand() / RAND_MAX,(float)rand() / RAND_MAX,(float)rand() / RAND_MAX);
      c_boxes[find_in_which_box(ctemp)].point_cloud.push_back(ctemp);
    }
    
//data passing should happen here

//make Q box for passing around
//same exactly for C (copy paste or something
for (int pro=0;pro<processes;pro++){
	if (pro!=processid){
		vector<Box<QPoint> > sendout;
		int sendcount=0;
		for (uint i=0;i<q_boxes.size();i++){
			if(get_box_owner(q_boxes[i].id==pro)){
				//ADD THIS BOX TO THE VECTOR TO BE SENT
				//IM GUESSING SOMETHING LIKE
				sendout[sendcount]=q_boxes[i];
				sendcount++;
			}
		}
		//here, we must send the sendout to the neighbor
	}
}

//POINT SEARCH
  //for each Q box
  for (uint i=0; i<q_boxes.size();i++){
      // if this box is not empty of Q points, and it belongs to me
      //cout<<"Searching box "<<i<<"\n";
      if (!q_boxes[i].point_cloud.empty() && is_my_box(q_boxes[i].id)){
          //for each Q point in this box
          for (uint j=0;j<q_boxes[i].point_cloud.size();j++){
              //cout<<"=========\nNext Point:\n";
              vector<Point> tentative_nn;
              tentative_nn=naive_search(q_boxes[i].point_cloud[j],c_boxes[i].point_cloud);
              
              //cout<<"Dist:"<<euclidean(tentative_nn[0],q_boxes[i].point_cloud[j])<<endl;
              
              //ANTONI!! IS THIS NECESSARY? WHAT IF WE DONT FIND A POSSIBLE C POINT IN THIS BOX BUT THERE IS ONE IN THE NEIGHBORS?
              if (!tentative_nn.empty()){
              
                  float cur_dist=euclidean(tentative_nn[0],q_boxes[i].point_cloud[j]);
                  
                  //for each possible direction on all sides of the box (26)
                  for (int dir=0;dir<27;dir++){
                      int temp_id=get_neighbor_id(dir, i);
                      bool should_we_check_neighbors=does_box_intersect_sphere(temp_id,q_boxes[i].point_cloud[j],cur_dist);
                      //check neighbor if a box exists in that direction
                      //and if there is some point in searching there
                      //else, dont check it/ignore it
                      if (temp_id>-1 && should_we_check_neighbors){
                      	  //if that box belongs to this process
                          if (is_my_box(temp_id)){
                              //we only have this info for C boxes in this process
                              if (!c_boxes[temp_id].point_cloud.empty()){
		                      //cout<<"Looking in neighbor "<<temp_id<<endl;
		                      vector<Point> temp=naive_search(q_boxes[i].point_cloud[j],c_boxes[temp_id].point_cloud);
		                      float d_temp=euclidean(q_boxes[i].point_cloud[j],temp[0]);
		                      if (d_temp<euclidean(q_boxes[i].point_cloud[j],tentative_nn[0]) ){
		                          tentative_nn=temp;
		                          cur_dist=d_temp;
		                          //cout<<"Found better C in neighbor "<<temp_id<<" at direction "<<dir<<endl;
		                        }
                                }
                            }
                          //else, if C is not my box, but belongs to a neighbor
                          else{
                              //make contact with temp_id neighbor box of neighbor process
                              //we must search in that box
                              //cout<<temp_id<<" belongs to another process. IMPLEMENT THIS!"<<endl;
                              //neighbor_proc_count++;
                              //neighbor receives the comm and must do the search. this has not been implemented either
                            }
                        }
                    }
                }
                /* testing of evaluation
                cout<<"Point Q ";
                printcoords(q_boxes[i].point_cloud[j]);
                cout<<" closest to \nPoint C ";
                printcoords(tentative_nn[0]);
                cout<<"\n=========\n\n";
                */
            }
        }
    }
}

std::vector<Point> naive_search(QPoint q, std::vector<Point> search_space){
  using std::vector;
  using std::cout;
  using std::endl;
  vector<Point> nn;
  uint siz=search_space.size();
  if (siz==0)
    {
      cout<<"Search space is empty"<<endl; }
  else {
      int min_index=0;
      float min_dist=2.0;
      float temp_dist;
      for (uint i=0;i<siz;i++)
        {
          temp_dist= euclidean(q,search_space[i]);
          if (temp_dist<min_dist)
            {
              min_dist=temp_dist;
              min_index=i;
            }
        }
      nn.push_back(search_space[min_index]);
    }
  return nn;
}


////independent print of results
////they are rounded to 3 decimal places, just for facilitating the view
//for (int i=0;i<qcountinthissplit;i++){
//  //printf("#%d Q point at coords (%.3f,%.3f,%.3f) is nearest to point C at coords (%.3f,%.3f,%.3f)\n\n",i+1,results[2*i].x,results[2*i].y,results[2*i].z,results[2*i+1].x,results[2*i+1].y,results[2*i+1].z);
//}

////NOTE THAT THE FOLLOWING ARE NOT IN ACCORDANCE TO THE EKFONISI
////suggestion/hack/slacking around: each process should process its boxes, meaning that it will check the Q points in its boxes
////however, it is a good idea that it keeps C points that are not only in its boxes, but also to the boxes adjacent to it.
////like that, when it is searching at the neighbour boxes of a box at an edge, it wont have to mpi call another process to get the data back
////this will only be done once at the beginning.
////ie each process will keep and receive the Q points in its boxes, and the C points in its boxes and adjacent boxes.
////this is a hack but it will save us from a lot of mpi calls during search, which might lead to many locks that will drop performance

////another hack would be to completely disregard searches of nearest neigbor in different boxes outside our split
////but it is bad, it will lead to mistakes at the borders.
////however, it will make message passing much much easier and ill use it if all else fails.

////note that there is a shitload of messages to be passed around during point search, so this might actually be a good idea...




#ifndef KNN_HPP
#define KNN_HPP

#endif // KNN_HPP
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <iomanip>
int rank; //for testing purposes - this should be dynamically allocated
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
  Point();
  Point(float x, float y, float z);
  std::vector<float> to_vector();
};
Point::Point(){
  //  x=(float)rand() / RAND_MAX;
  //  y=(float)rand()/ RAND_MAX;
  //  z=(float)rand()/ RAND_MAX;
}

Point::Point(float x, float y, float z){
  this->x=x;
  this->y=y;
  this->z=z;
}

std::vector<float> Point::to_vector(){
  std::vector<float> vec;
  vec.push_back(x);
  vec.push_back(y);
  vec.push_back(z);
  return vec;
}
struct PointMsg{
  float x;
  float y;
  float z;
};

struct QPoint : public Point{
  Point nn;
  bool found_nn;
  QPoint();
  ~QPoint(){}
  QPoint(float x, float y, float z);
};

QPoint::QPoint() : Point()
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
int get_owning_process(Point p){
  return get_box_owner(find_in_which_box(p));
}

//check if given boxid belongs to current process
int is_my_box(int id){
  if (get_box_owner(id)==rank) return 1;
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
  //printf(s"%d ",temp);
  return temp;
}

//helper function for using getneigborid(int id, int dirx, int diry, int dirz)
//in for loops from 0 to 26, for checking all neigbors
//one of the iterations of the loop will return Q's box. in that case,
//we return -2 for reference
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


//check if given box is closer than the candidate distance of q from candidate c
//this is done by using this generic code, that sees if a sphere intersects
//a cuboid. returns true if it intersect, false if it doesnt.
float squared(float v) { return v * v; }
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
  if (q.x < xmin) distancesquared -= squared(q.x - xmin);
  else if (q.x > xmax) distancesquared -= squared(q.x - xmax);
  if (q.y < ymin) distancesquared -= squared(q.y - ymin);
  else if (q.y > ymax) distancesquared -= squared(q.y - ymax);
  if (q.z < zmin) distancesquared -= squared(q.z - zmin);
  else if (q.z > zmax) distancesquared -= squared(q.z - zmax);
  return distancesquared > 0;
}


std::vector<Point> naive_search(QPoint q, std::vector<Point> search_space){
  using std::vector;
  using std::cout;
  using std::endl;
  vector<Point> nn;
  uint siz=search_space.size();
  if (siz==0)
    {
      cout<<"Search space is empty"<<endl;
    }
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

void generate_random_points(std::vector<Point>& pts,const int num ){
  pts.reserve(num);
  for (int i=0;i<num;i++){
      pts.push_back(
            Point((float)rand() / RAND_MAX,(float)rand() / RAND_MAX,(float)rand() / RAND_MAX));
    }
}

void assign_points_to_proccesses(const std::vector<Point>& pts,
                                 const int proc_count,
                                 std::vector<std::vector<Point> >& sorted_points){
  sorted_points.resize(proc_count);
  for (uint i=0;i<pts.size();i++){
      sorted_points[get_owning_process(pts[i])].push_back(pts[i]);
    }
}

void prepare_scatterv_msg(const std::vector<std::vector<Point> >& proc_array,
                          std::vector<float>& sendbuf,
                          std::vector<int>& count,
                          std::vector<int>& displ){
  using namespace std;
  count.reserve(proc_array.size());
  displ.reserve(proc_array.size());
  displ[0]=0;
  for (uint i=0;i<proc_array.size();i++)  {
      cout<<endl;
      for (uint j=0;j<proc_array[i].size();j++){
          sendbuf.push_back(proc_array[i][j].x);
          sendbuf.push_back(proc_array[i][j].y);
          sendbuf.push_back(proc_array[i][j].z);
//          if (i==1)cout<<std::setprecision(4)<<proc_array[0][j].x<<" "<<proc_array[0][j].y<<" "<<proc_array[0][j].z<<" "<<proc_array[i][j].x<<" "<<proc_array[i][j].y<<" "<<proc_array[i][j].z<<endl;
          cout<<setprecision(3);
//          cout<<i<<","<<j<<" "<<proc_array[i][j].x<<" "<<proc_array[i][j].y<<" "<<proc_array[i][j].z<<endl;
        }
      count[i]=3*proc_array[i].size();
      if (i) displ[i]=displ[i-1]+count[i-1];
    }
}

struct BoundaryMsg{
  float pt[3];
  int box;
  void set(std::vector<float> vec, int box_id);
};

void BoundaryMsg::set(std::vector<float> vec, int box_id){
  for (int i=0;i<3;i++) pt[i]=vec[i];
  box=box_id;
}

template <class T> std::vector<T> flatten(const std::vector<std::vector<T> >& vec_2d) {
  std::vector<T> flat_vec;
  for (uint i=0;i<vec_2d.size();i++){
      for (uint j=0;j<vec_2d[i].size();j++){
          flat_vec.push_back(vec_2d[i][j]);
        }
    }
  return flat_vec;
}

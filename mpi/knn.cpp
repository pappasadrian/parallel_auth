#include "knn.hpp"

int main(int argc, char **argv){
  using std::vector;
  using std::cout;
  using std::endl;
  if (argc != 5) {
      cout<<"Usage:"<<argv[0]<<" Q C P B\n";
      cout<<"Q: number of Q points (power of two - from 20 to 25)\n";
      cout<<"C: number of C points (power of two - from 20 to 25)\n";
      cout<<"P: number of processes (power of two - from 0 to 7)\n";
      cout<<"B: number of boxes (power of two - from 12 to 16)\n";
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

  for (uint i=0; i<q_boxes.size();i++){
      if (!q_boxes[i].point_cloud.empty() ){
          for (uint j=0;j<q_boxes[i].point_cloud.size();j++){
              vector<Point> tentative_nn;
              tentative_nn=naive_search(q_boxes[i].point_cloud[j],c_boxes[i].point_cloud);
              //cout<<"Dist:"<<euclidean(tentative_nn[0],q_boxes[i].point_cloud[j])<<endl;
              if (!tentative_nn.empty()){
                  float cur_dist=euclidean(tentative_nn[0],q_boxes[i].point_cloud[j]);
                  for (int dir=0;dir<27;dir++){
                      int temp_id=get_neighbor_id(dir, i);
                      bool should_we_check_neighbors=does_box_intersect_sphere(temp_id,q_boxes[i].point_cloud[j],cur_dist);
                      if (temp_id>-1 && should_we_check_neighbors ){
                          if (is_my_box(temp_id) && !c_boxes[temp_id].point_cloud.empty()){
                              cout<<"LOOKING IN THE PROCESS"<<endl;
                              vector<Point> temp=naive_search(q_boxes[i].point_cloud[j],c_boxes[temp_id].point_cloud);
                              float d_temp=euclidean(q_boxes[i].point_cloud[j],temp[0]);
                              if (d_temp<euclidean(q_boxes[i].point_cloud[j],tentative_nn[0]) ){
                                  tentative_nn=temp;
                                  cur_dist=d_temp;
                                  cout<<"Found better in neighbor "<<dir<<endl;
                                }
                            }
                          else{
                              //cout<<is_my_box(temp_id)<<":In other process. IMPLEMENT THIS!"<<endl;
                              //neighbor_proc_count++;
                            }
                        }
                    }
                }
            }
        }
    }
}



////independent print of results
////they are rounded to 3 decimal places, just for facilitating the view
//for (int i=0;i<qcountinthissplit;i++){
//  //printf("#%d Q point at coords (%.3f,%.3f,%.3f) is nearest to point C at coords (%.3f,%.3f,%.3f)\n\n",i+1,results[2*i].x,results[2*i].y,results[2*i].z,results[2*i+1].x,results[2*i+1].y,results[2*i+1].z);
//}

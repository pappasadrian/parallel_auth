#include "knn.hpp"
#include <mpi/mpi.h>

int main(int argc, char **argv){
  MPI_Init(&argc,&argv);
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
  MPI_Comm_size(MPI_COMM_WORLD,&processes);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, rank, processes);

  Pprocesses=atoi(argv[3]);
  int my_proc_num=1<<Pprocesses;
  if (my_proc_num!=processes){
      cout<<"ERROR: Inconsistent process numbers. Exiting";
      cout<<my_proc_num<<"!="<<processes;
      MPI_Abort(MPI_COMM_WORLD,1);
      exit(123);
    }

  // from 0(?) to 25 -in all processes. this process has numberofpoints/processes Q points.
  PQnumberofpoints=atoi(argv[1]);
  qnumberofpoints=1<<PQnumberofpoints;
  cout<<"Qnum="<<qnumberofpoints<<endl;
  // from 0(?) to 25 -in all processes. this process has numberofpoints/processes C points.
  PCnumberofpoints=atoi(argv[2]);
  cnumberofpoints=1<<PCnumberofpoints;
  cout<<"cnum="<<cnumberofpoints<<endl;
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
  //generally, we can refer to these boxes either by ID or by coordinates. there
  //are functions that do this.
  boxdimensions[0]=1<<(Pnumberboxes/3);
  if (Pnumberboxes%3 > 0)  boxdimensions[0]*=2;
  boxdimensions[1]=1<<(Pnumberboxes/3);
  if (Pnumberboxes%3 == 2)	boxdimensions[1]*=2;
  boxdimensions[2]=1<<(Pnumberboxes/3);

  splitdimensions[0]=1<<(Pprocesses/3);
  if (Pprocesses%3 > 0)	splitdimensions[0]*=2;
  splitdimensions[1]=1<<(Pprocesses/3);
  if (Pprocesses%3 == 2)	splitdimensions[1]*=2;
  splitdimensions[2]=1<<(Pprocesses/3);
  for(int i=0;i<3;i++) boxespersplit[i]=boxdimensions[i]/splitdimensions[i];

  vector<Point> c_points;
  vector<vector<Point> >  c_pts_in_proc;
  vector<float> c_sendbuff;
  vector<int> c_count,c_displ;
  vector<float> c_recvbuf;
  c_recvbuf.reserve(4*c_pts_per_process);
  c_count.reserve(processes);
  if (rank==0){
      generate_random_points(c_points,cnumberofpoints);
      assign_points_to_proccesses(c_points,processes,c_pts_in_proc);
      prepare_scatterv_msg(c_pts_in_proc,c_sendbuff,c_count,c_displ);
      cout<<"Sending "<<c_pts_in_proc[0][0].x<<endl;
    }
  //Sending counts of each process' search space, and then scattering the points
  //of said space to each process. Should be done with MPI_Struct or MPI_Pack but
  //AINT NOBODY GOT TIME FOR THAT
  MPI_Bcast(&c_count[0],processes,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatterv(&c_sendbuff[0],&c_count[0],&c_displ[0],MPI_FLOAT,&c_recvbuf[0],4*c_pts_per_process,MPI_FLOAT,0,MPI_COMM_WORLD);
  //  if (rank){
  cout<<std::setprecision(3)<<rank<< " received " <<c_recvbuf[0]<<endl;
  cout<<std::setprecision(3)<<rank<< " count " <<c_count[rank]<<endl;
  //  for (int i=0;i<c_count[rank];i+=3)
  //    cout<<rank<<" "<<c_recvbuf[i]<<" "<<c_recvbuf[i+1]<<" "<<c_recvbuf[i+2]<<endl;
  MPI_Barrier(MPI_COMM_WORLD);


  cout<<"================ Q POINTS ======================="<<endl;
  vector<Point> q_points;
  vector<vector<Point> >  q_pts_in_proc;
  vector<float> q_sendbuff;
  vector<int> q_count,q_displ;
  vector<float> q_recvbuf;
  q_recvbuf.reserve(6*q_pts_per_process);
  q_count.reserve(processes);
  if (rank==0){
      generate_random_points(q_points,qnumberofpoints);
      assign_points_to_proccesses(q_points,processes,q_pts_in_proc);
      prepare_scatterv_msg(q_pts_in_proc,q_sendbuff,q_count,q_displ);
      cout<<"Sending "<<q_pts_in_proc[0][0].x<<endl;
    }
  //Sending counts of each process' search space, and then scattering the points
  //of said space to each process. Should be done with MPI_Struct or MPI_Pack but
  //AINT NOBODY GOT TIME FOR THAT
  MPI_Bcast(&q_count[0],processes,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatterv(&q_sendbuff[0],&q_count[0],&q_displ[0],MPI_FLOAT,&q_recvbuf[0],6*q_pts_per_process,MPI_FLOAT,0,MPI_COMM_WORLD);
  //  if (rank){
  cout<<std::setprecision(3)<<rank<< " received " <<q_recvbuf[0]<<endl;
  cout<<std::setprecision(3)<<rank<< " count " <<q_count[rank]<<endl;
    for (int i=0;i<q_count[rank];i+=3)
      cout<<rank<<" "<<q_recvbuf[i]<<" "<<q_recvbuf[i+1]<<" "<<q_recvbuf[i+2]<<endl;



  vector<Box<Point> > c_boxes;
  c_boxes.resize(numboxes);
  for (uint i=0;i<c_boxes.size();i++){
      c_boxes[i].id=i;
      c_boxes[i].coords=get_box_coords(i);
      c_boxes[i].owner=get_box_owner(i);
    };


  vector<Box<QPoint> > q_boxes;
  q_boxes.resize(numboxes);
  for (uint i=0;i<q_boxes.size();i++){
      if (is_my_box(i) ){
          q_boxes[i].id=i;
          q_boxes[i].coords=get_box_coords(i);
          q_boxes[i].owner=get_box_owner(i);
        }
    };

  for (int i=0; i<c_count[rank];i+=3){
      Point temp=Point(c_recvbuf[i],c_recvbuf[i+1],c_recvbuf[i+2]);
      c_boxes[find_in_which_box(temp)].point_cloud.push_back(temp);
    }

  for (int i=0; i<q_count[rank];i+=3){
      QPoint temp=QPoint(q_recvbuf[i],q_recvbuf[i+1],q_recvbuf[i+2]);
      q_boxes[find_in_which_box(temp)].point_cloud.push_back(temp);
    }

  for (uint i=0; i<q_boxes.size();i++){
      if (!q_boxes[i].point_cloud.empty() && rank==q_boxes[i].id){
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


  MPI_Finalize();
}



////independent print of results
////they are rounded to 3 decimal places, just for facilitating the view
//for (int i=0;i<qcountinthissplit;i++){
//  //printf("#%d Q point at coords (%.3f,%.3f,%.3f) is nearest to point C at
//coords (%.3f,%.3f,%.3f)\n\n",i+1,results[2*i].x,results[2*i].y,results[2*i].z,results[2*i+1].x,results[2*i+1].y,results[2*i+1].z);
//}





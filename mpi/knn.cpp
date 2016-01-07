#include "knn.hpp"
#include <mpi.h>
#include <ctime>
#include <sys/time.h>
#include <time.h>

int main(int argc, char **argv){
  MPI_Init(&argc,&argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  using std::vector;
  using std::cout;
  using std::endl;
  struct timeval startwtime, endwtime;
  double execute_time;
  if (argc != 5) {
      cout<<"Usage:"<<argv[0]<<" Q C P B\n";
      cout<<"Q: number of Q points (power of two - from 20 to 25)\n";
      cout<<"C: number of C points (power of two - from 20 to 25)\n";
      cout<<"P: number of processes (power of two - from 0 to 7)\n";
      cout<<"B: number of boxes (power of two - from 12 to 16)\n";
      return 1;
    }
  MPI_Comm_size(comm,&processes);
  MPI_Comm_rank(comm,&rank);
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
      MPI_Abort(comm,1);
      exit(123);
    }

  // from 0(?) to 25 -in all processes. this process has numberofpoints/processes Q points.
  PQnumberofpoints=atoi(argv[1]);
  qnumberofpoints=1<<PQnumberofpoints;
//  cout<<"Qnum="<<qnumberofpoints<<endl;
  // from 0(?) to 25 -in all processes. this process has numberofpoints/processes C points.
  PCnumberofpoints=atoi(argv[2]);
  cnumberofpoints=1<<PCnumberofpoints;
//  cout<<"cnum="<<cnumberofpoints<<endl;
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
      gettimeofday(&startwtime,NULL); //Implicit MPI_Barrier below
    }
  //Sending counts of each process' search space, and then scattering the points
  //of said space to each process. Should be done with MPI_Struct or MPI_Pack but
  //AINT NOBODY GOT TIME FOR THAT
  MPI_Bcast(&c_count[0],processes,MPI_INT,0,comm);
  MPI_Scatterv(&c_sendbuff[0],&c_count[0],&c_displ[0],MPI_FLOAT,&c_recvbuf[0],4*c_pts_per_process,MPI_FLOAT,0,comm);



  vector<Point> q_points;
  vector<vector<Point> >  q_pts_in_proc;
  vector<float> q_sendbuff;
  vector<int> q_count,q_displ;

  q_count.reserve(processes);
  if (rank==0){
      generate_random_points(q_points,qnumberofpoints);
      assign_points_to_proccesses(q_points,processes,q_pts_in_proc);
      prepare_scatterv_msg(q_pts_in_proc,q_sendbuff,q_count,q_displ);
      //      cout<<"Sending "<<q_pts_in_proc[0][0].x<<endl;
    }
  //Sending counts of each process' search space, and then scattering the points
  //of said space to each process. Should be done with MPI_Struct or MPI_Pack but
  //AINT NOBODY GOT TIME FOR THAT
  MPI_Bcast(&q_count[0],processes,MPI_INT,0,comm);
  vector<float> q_recvbuf;
  q_recvbuf.resize(3*q_count[rank]);
  MPI_Scatterv(&q_sendbuff[0],&q_count[0],&q_displ[0],MPI_FLOAT,&q_recvbuf[0],6*q_pts_per_process,MPI_FLOAT,0,comm);

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


  vector<float> rslts_sendbuf;
  vector<vector<BoundaryMsg> > boundary_pts;
  boundary_pts.resize(processes);
  vector<vector<int> > boundary_pts_index;
  boundary_pts_index.resize(2);

  ;
  for (uint i=0; i<q_boxes.size();i++){
      if (!q_boxes[i].point_cloud.empty() && rank==q_boxes[i].owner){
          for (uint j=0;j<q_boxes[i].point_cloud.size();j++){
              vector<Point> tentative_nn;
              tentative_nn=naive_search(q_boxes[i].point_cloud[j],c_boxes[i].point_cloud);
              float cur_dist=2;
              if (!tentative_nn.empty()){
                  cur_dist=euclidean(tentative_nn[0],q_boxes[i].point_cloud[j]);
                }
              for (int dir=0;dir<27;dir++){
                  int temp_id=get_neighbor_id(dir, i);
                  if (temp_id>numboxes-1) temp_id=-2;
                  bool should_we_check_neighbors=does_box_intersect_sphere(temp_id,q_boxes[i].point_cloud[j],cur_dist);
                  if (temp_id>-1 && should_we_check_neighbors ){
                      if (is_my_box(temp_id) && !c_boxes[temp_id].point_cloud.empty()){
                          vector<Point> temp=naive_search(q_boxes[i].point_cloud[j],c_boxes[temp_id].point_cloud);
                          float d_temp=euclidean(q_boxes[i].point_cloud[j],temp[0]);
                          if (d_temp<euclidean(q_boxes[i].point_cloud[j],tentative_nn[0]) ){
                              tentative_nn=temp;
                              cur_dist=d_temp;
                            }
                        }
                      else{
                          int split_id=get_box_owner(temp_id);
                          BoundaryMsg temp;
                          temp.set(q_boxes[i].point_cloud[j].to_vector(),temp_id);
                          boundary_pts[split_id].push_back(temp);
                          boundary_pts_index[0].push_back(i);
                          boundary_pts_index[1].push_back(j);
                        }
                    }
                }

              if (!tentative_nn.empty()){
                  q_boxes[i].point_cloud[j].nn=tentative_nn[0];
                }
              else {
                  q_boxes[i].point_cloud[j].nn.x=2;
                  q_boxes[i].point_cloud[j].nn.y=2;
                  q_boxes[i].point_cloud[j].nn.z=2;
                }
            }
        }
    }

  vector<int> pts_sendcount,pts_sdispl;
  pts_sendcount.resize(processes);
  pts_sdispl.resize(processes);
  pts_sdispl[0]=0;
  for (int i=0;i<processes;i++){
      if (i) pts_sdispl[i]=pts_sdispl[i-1]+pts_sendcount[i-1];
      pts_sendcount[i]=(int)boundary_pts[i].size();
    }
  vector<int> pts_recvcount;
  pts_recvcount.resize(processes);
  MPI_Alltoall(&pts_sendcount[0],1,MPI_INT,&pts_recvcount[0],1,MPI_INT,comm);
  vector<int> pts_rdispl;
  pts_rdispl.resize(processes);
  int pts_recvbuf_len=0;
  for (int i=0;i<processes;i++){
      pts_recvbuf_len+=pts_recvcount[i];
      if(i) pts_rdispl[i]=pts_rdispl[i-1]+pts_recvcount[i-1];
    }

  vector<BoundaryMsg> pts_sendbuf=flatten(boundary_pts);
  vector<BoundaryMsg> pts_recvbuf;
  pts_recvbuf.resize(pts_recvbuf_len);


  MPI_Datatype BoundaryMsg_MPI;
  MPI_Datatype type[2]; //I miss c++11
  type[0]=MPI_FLOAT;
  type[1]=MPI_INT;
  int blocklen[2];
  blocklen[0]=3;
  blocklen[1]=1;
  MPI_Aint disp[2];
  disp[0]=0;
  disp[1]=3*sizeof(float);
  MPI_Type_struct(2,blocklen,disp,type,&BoundaryMsg_MPI);
  MPI_Type_commit(&BoundaryMsg_MPI);

  MPI_Alltoallv(&pts_sendbuf[0],&pts_sendcount[0],&pts_sdispl[0],BoundaryMsg_MPI,
      &pts_recvbuf[0],&pts_recvcount[0],&pts_rdispl[0],BoundaryMsg_MPI,comm);

  vector<PointMsg> nn_sendbuff;
  nn_sendbuff.resize(pts_recvbuf.size());
  for (uint i=0;i<pts_recvbuf.size();i++){
      QPoint temp(pts_recvbuf[i].pt[0],pts_recvbuf[i].pt[1],pts_recvbuf[i].pt[2]);
      vector<Point> nn=naive_search(temp,c_boxes[pts_recvbuf[i].box].point_cloud);
      if (nn.empty() ) nn.push_back(Point(2,2,2));
      nn_sendbuff[i]=nn[0].to_msg();
    }
  vector<PointMsg> nn_recvbuff;
  nn_recvbuff.resize(pts_sendbuf.size());


  MPI_Datatype PointMsg_MPI;
  MPI_Type_contiguous(3,MPI_FLOAT,&PointMsg_MPI);
  MPI_Type_commit(&PointMsg_MPI);

  MPI_Alltoallv(&nn_sendbuff[0],&pts_recvcount[0],&pts_rdispl[0],PointMsg_MPI,
      &nn_recvbuff[0],&pts_sendcount[0],&pts_sdispl[0],PointMsg_MPI,comm);

  for (uint i=0;i<boundary_pts_index[0].size();i++){
      int ii=boundary_pts_index[0][i];
      int jj=boundary_pts_index[1][i];
      Point foreign_nn=nn_recvbuff[i].to_point();
      float d_old=euclidean(q_boxes[ii].point_cloud[jj],q_boxes[ii].point_cloud[jj].nn);
      float d_new=euclidean(q_boxes[ii].point_cloud[jj], foreign_nn);
      if (d_new<d_old) q_boxes[ii].point_cloud[jj].nn=foreign_nn;
    }



  for (uint i=0;i<q_boxes.size();i++){
      for (uint j=0;j<q_boxes[i].point_cloud.size();j++){
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].x);
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].y);
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].z);
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].nn.x);
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].nn.y);
          rslts_sendbuf.push_back(q_boxes[i].point_cloud[j].nn.z);
        }
    }
  vector<float> rslts_recvbuf;
  vector<int> rslts_cnt;
  vector<int> rslts_displs;

  if (rank==0){
      rslts_recvbuf.reserve(6*qnumberofpoints);
      rslts_cnt.reserve(qnumberofpoints);
      for (int i=0;i<processes;i++) rslts_cnt[i]=2*q_count[i];
      rslts_displs.reserve(qnumberofpoints);
      rslts_displs[0]=0;
      for (int i=0;i<processes;i++){
          rslts_displs[i]=rslts_cnt[i-1]+rslts_displs[i-1];
        }
    }

  MPI_Gatherv(&rslts_sendbuf[0],(int)rslts_sendbuf.size(),MPI_FLOAT,
      &rslts_recvbuf[0],&rslts_cnt[0],&rslts_displs[0],MPI_FLOAT,0,
      comm);
  if (rank==0){
 gettimeofday(&endwtime,NULL);
  execute_time = (double)( ( endwtime.tv_usec - startwtime.tv_usec ) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec );
 cout<<"MPI duration"<<execute_time<<endl;
      //      MPI_t
      //      double single_time= MPi
 gettimeofday(&startwtime,NULL);
 vector<Box<QPoint> > single_results=single_threaded_search(c_points,q_points,numboxes);
 gettimeofday(&endwtime,NULL);
 execute_time = (double)( ( endwtime.tv_usec - startwtime.tv_usec ) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec );
 cout<<"Single threaded duration "<<execute_time<<endl;
      cout<<single_results.size();

    /**  Validation code
       * Simple naive search across the whole space
       * Because it is too slow uncomment only when making changes*/
//      int wrong_results=0;
//      vector<QPoint> results;
//      results.resize(qnumberofpoints);
//      int j=0;
//      for (int i=0;i<6*qnumberofpoints;i+=6){
//          results[j]=QPoint(rslts_recvbuf[i],rslts_recvbuf[i+1],rslts_recvbuf[i+2]);
//          results[j].nn=Point(rslts_recvbuf[i+3],rslts_recvbuf[i+4],rslts_recvbuf[i+5]);
//          j++;
//        }
//      for (int i=0;i<qnumberofpoints;i++){
//          vector<Point> nn_real=naive_search(results[i],c_points);
//          float d_found=euclidean(results[i],results[i].nn);
//          float d_real=euclidean(results[i],nn_real[0]);
//          if (abs(d_real-d_found)>0.00001) wrong_results++;
//        }
//      cout<<"Misclassification rate: "<< (float)wrong_results/qnumberofpoints<<endl;
//      cout<<"wrong"<<wrong_results<<endl;

    }
  MPI_Finalize();
}






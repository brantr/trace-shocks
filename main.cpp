#include <stdio.h>
#include <math.h>
#include <vector>
#include "shock_data_types.hpp"       /*tracer and shock data structers*/
#include "write_shock_catalogues.hpp" //correct catalogues
#include "timer.h"


tracer tin;                       /* buffer for adding tracers to tracer vectors */

struct interaction
{
  int snap_A;
  int snap_B;
  long idx_A;
  long id_A;
  long l_A;
  long o_A;
  float d_A;
  long idx_B;
  long id_B;
  long l_B;
  long o_B;
  float d_B;
  float frac_A;
  float frac_B;
  float frac_A_dense;
  float frac_B_dense;
  float frac_dense;
  float x_A[3];
  float x_B[3];
  float x_C[3];
  long n;
};
bool interaction_sort(interaction ia, interaction ib)
{
  //sort by decreasing frac_A
  return ia.frac_A>ib.frac_A;
}
bool interaction_sort_dense(interaction ia, interaction ib)
{
  //sort by decreasing frac_A
  return ia.frac_A_dense>ib.frac_A_dense;
}
void show_shocks(vector<shock> s);
void show_interactions(vector<interaction> ia);
void write_interactions(char fname[], vector<interaction>  ia, float *times);
void read_interactions( char fname[], vector<interaction> *ia);
void keep_duplicates(vector<long> iunion, vector<long> *ioverlap);
void write_interaction_counts(int iA, int iB, long n_total, vector<long> n_ints, vector<long> o_ints);
void read_interaction_counts(char fname[], vector<long> *n_ints, vector<long> *o_ints);
void write_branches(char fname[], vector<long> ia, vector<long> ib);
int main(int argc, char **argv)
{
  char fdir[200];
  char fdint[200];
  char fbint[200];
  char fint[200];
  char fbase[200];
  char fcount[200];
  char flist_A[200];
  char flist_B[200];
  char fdata_A[200];
  char fdata_B[200];
  char foutput[200];
  char fnints[200];
  char fdtrace[200];
  char fbranch[200];


  int  ishock = 0;
  int imin = 700;
  int imax = 901;

  int  iA;
  int  iB;

  long nA;
  long nB;

  long tt;
  long ss;

  float dr = sqrt(3)/512.;

  int nsearch = 10;

  float frac_A;
  float frac_B;
  float d_A;
  float d_B;
  float x_A[3];
  float x_B[3];

  vector<shock>  sA;
  vector<shock>  sB;
  vector<tracer> tA;
  vector<tracer> tB;

  vector<long> n_ints;
  vector<long> o_ints;
  long n_ints_total;
  long n_ints_snap;
  long n_ints_in;

  int fflag = 0; //forward mode?

  //vector<long> l_ia;
  //vector<long> o_ia;
  vector<interaction> ia;
  vector<interaction> ia_tmp;
  vector<interaction> ia_shock;
  interaction ia_new;

  vector<long> nbranches;
  vector<long> nbinteractions;


  kdtree2 *bp_tree;
  kdtree2_result_vector res;
  array2dfloat bp_tree_data;
  vector<float> xc(3);

  long ioff = 0;

  float fdr = 3.; //how many cell distances?

  double time_A;
  double time_B;

  const int max_red = 5;

  array<vector<long>, max_red> islist;
  vector<long>::iterator il;

  int *snap_list; //list of snapshots to use in tracking
  int n_snaps;    //total number of snapshots to use

  FILE *fp_times;
  char fname_times[200];
  int n_times;
  float *times;

  int n_redundancy = 1;

  time_A = timer();

  if(argc<4)
  {
    printf("./trace_shocks imin imax ishock [fname] [n_redundancy] [forward_flag]\n");
    exit(-1);
  }

  //if we've supplied a range of snapshots
  //to search, use that
  if(argc>=4)
  {
  	imin   = atoi(argv[1]);
    imax   = atoi(argv[2]);
  	ishock = atoi(argv[3]);
  }
  if(argc>=6)
  {
    n_redundancy = atoi(argv[5]);
  }
  printf("n_redundancy = %d\n",n_redundancy);

  if(argc>=5)
  {
  	char fname_snap_list[200];



  	sprintf(fname_snap_list,"%s",argv[4]);

   	FILE *fp_snap_list;

   	if(!(fp_snap_list = fopen(fname_snap_list,"r")))
   	{
   		printf("Error opening %s\n",fname_snap_list);
   		exit(-1);
   	}

   	fscanf(fp_snap_list,"%d\n",&n_snaps);
    printf("n_snaps = %d\n",n_snaps);

    snap_list = (int *)malloc(n_snaps*sizeof(int));

   	for(int i=0;i<n_snaps;i++)
   	{
   		fscanf(fp_snap_list,"%d\n",&snap_list[i]);
   		//printf("snap_list[%d] = %d\n",i,snap_list[i]);
   	}
    fclose(fp_snap_list);

    if(argc>=6)
      fflag = atoi(argv[5]);

    if(fflag)
    {
      //forward
      imax = snap_list[n_snaps-1];
      imin = snap_list[0];
    }else{
      //backward
      imin = snap_list[n_snaps-1];
      imax = snap_list[0];
    }
  }else{

    n_snaps = imax-imin+1;
    snap_list = (int *)malloc(n_snaps*sizeof(int));

    for(int i=0;i<n_snaps;i++)
      snap_list[i] = imax - i;
  }

  //load snapshot times
  sprintf(fname_times,"times.txt");
  if(!(fp_times = fopen(fname_times,"r")))
  {
    printf("Error opening %s.\n",fname_times);
    exit(-1);
  }
  fscanf(fp_times,"%d\n",&n_times);
  times = (float *) malloc(n_times*sizeof(float));
  for(int i=0;i<n_times;i++)
  	fscanf(fp_times,"%f\n",&times[i]);
  fclose(fp_times);



  printf("ishock  = %d\n",ishock);
  printf("imin    = %d %d\n",imin,snap_list[n_snaps-1]);
  printf("imax    = %d %d\n",imax,snap_list[0]);

  sprintf(fdir,"data/");
  sprintf(fdint,"interactions/");
  sprintf(fbint,"interactions");
  sprintf(fbase,"peak.blended");
  sprintf(fdtrace,"traces/");

  
  sprintf(fint,"%s%s.%04d.%04d.%d.txt",fdint,fbint,imin,imax,n_redundancy);
  sprintf(foutput,"%strace.%04d.%04d.%d.%08d.txt",fdtrace,imin,imax,n_redundancy,ishock);
  sprintf(fbranch,"%sbranches.%04d.%04d.%d.%08d.txt",fdtrace,imin,imax,n_redundancy,ishock);

  //Read interactions
  printf("Reading %s\n",fint);
  read_interactions(fint,&ia);

  printf("ia.size() %ld\n",ia.size());

  long m;
//  for(iA = imax; iA>imin; iA--)

  //need to adjust islist to allow for 
  //tracking of shocks in each of the redundant 
  //snapshot times

  //the first list of shocks is just the
  //ishock of interest
  islist[0].push_back(ishock);

  for(int i_snap = 0;i_snap<n_snaps-n_redundancy; i_snap++)
  {
  	iA = snap_list[i_snap]; 

    for(int i_red=1;i_red <= n_redundancy; i_red++)
    {
      iB = snap_list[i_snap+i_red];

      printf("i_red %d iA %d iB %d\n",i_red,iA,iB);
      if(fflag)
      {
        sprintf(fcount,"%sinteraction_count.%04d.%04d.txt",fdint,iA,iB);
      }else{
        sprintf(fcount,"%sinteraction_count.%04d.%04d.txt",fdint,iB,iA);
      }
      printf("reading %s\n",fcount);
      read_interaction_counts(fcount,&n_ints,&o_ints);
      printf("******\n");
      for(int ss=0;ss<islist[0].size();ss++)
      {
        ishock = islist[0][ss];

        if(ishock>=n_ints.size())
        {
          printf("ERROR iA %04d iB %04d ishock %d n_ints.size() %ld\n",iA,iB,ishock,n_ints.size());
          exit(-1);
        }

        //print information about interaction counts
        printf("iA %04d iB %04d ishock %d n_ins %ld o_ints %ld\n",iA,iB,ishock,n_ints[ishock],o_ints[ishock]);

        for(int i=0;i<n_ints[ishock];i++)
        {
          ia_tmp.push_back(ia[ioff + o_ints[ishock] + i]);
          m = ia_tmp.size()-1;
          printf("iA %04d idxA %6ld idxB %6ld nints %4ld idA %10ld idB %10ld nex %8ld fA %5.4e fB %5.4e fdA %5.4e fdB %5.4e f %5.4e\n",ia_tmp[m].snap_A,ia_tmp[m].idx_A,ia_tmp[m].idx_B,n_ints[ishock],ia_tmp[m].id_A,ia_tmp[m].id_B,ia_tmp[m].n,ia_tmp[m].frac_A,ia_tmp[m].frac_B,ia_tmp[m].frac_A_dense,ia_tmp[m].frac_B_dense,ia_tmp[m].frac_dense);
        } //end loop over interactions
      }// end loop over list of traced shocks

      //reset the shock list
      //vector<long>().swap(islist[i_red]);
      for(int i=0;i<ia_tmp.size();i++)
      {
        ia_shock.push_back(ia_tmp[i]);
        islist[i_red].push_back(ia_tmp[i].idx_B);
      }

      //remember the number of interactions
      nbinteractions.push_back(ia_tmp.size());

      //keep unique
      std::sort(islist[i_red].begin(), islist[i_red].end());
      il = std::unique(islist[i_red].begin(), islist[i_red].end());
      islist[i_red].resize( std::distance(islist[i_red].begin(), il) );

      for(int i=0;i<islist[i_red].size();i++)
        printf("i_red %d i %d islist[%d] %ld\n",i_red,i,i,islist[i_red][i]);

      //if there are no interactions, break the loop
      if(ia_tmp.size()==0)
        break;

      printf("HERE\n");

      //remember the current number of branches
      //nbranches.push_back(ia_tmp.size());
      nbranches.push_back(islist[i_red].size());
      printf("nbranches %ld\n",nbranches[nbranches.size()-1]);

      //reset temporary interaction list
      vector<interaction>().swap(ia_tmp);

      //advance ioff
      long n_int_sum = 0;
      for(int i=0;i<n_ints.size();i++)
        n_int_sum += n_ints[i];

      //ioff += n_ints.size();
      ioff += n_int_sum;
      printf("ioff = %ld n_int_sum %ld n_ints.size() %ld\n",ioff,n_int_sum,n_ints.size());

      vector<long>().swap(n_ints);
      vector<long>().swap(o_ints);
    }//end loop over second snapshots

    for(int i=0;i<n_redundancy-1;i++)
      islist[i].swap(islist[i+1]);
    vector<long>().swap(islist[n_redundancy-1]);

  }//end loop over first snapshots



  //for(int i=0;i<ia.size();i++)
	//printf("snap_A %04d\tia %10ld\tib %10ld\tn %10ld\tfrac A %5.4e\tfrac B %5.4e\txA %e %e %e\txB %e %e %e\n",ia[i].snap_A,ia[i].id_A,ia[i].id_B,ia[i].n,ia[i].frac_A,ia[i].frac_B,ia[i].x_A[0],ia[i].x_A[1],ia[i].x_A[2],ia[i].x_B[0],ia[i].x_B[1],ia[i].x_B[2]);

  //save the interactions to a file
  write_interactions(foutput,ia_shock,times);

  write_branches(fbranch,nbranches,nbinteractions);

  time_B = timer();

  printf("Total time = %es.\n",time_B-time_A);

  return 0;
}
void show_shocks(vector<shock> s)
{
  int nlim = 10;
  printf("***************\n");
  for(int i=0;i<nlim;i++)
	  printf("i %6d\tl %10ld\to %10ld\td %15.14e\tid %10ld\tmin %5.4e %5.4e %5.4e\tmax %5.4e %5.4e %5.4e\n",i,s[i].l,s[i].o,s[i].d,s[i].id,s[i].min[0],s[i].min[1],s[i].min[2],s[i].max[0],s[i].max[1],s[i].max[2]);
  printf("***************\n");
  for(int i=s.size()-nlim;i<s.size();i++)
	  printf("i %6d\tl %10ld\to %10ld\td %15.14e\tid %10ld\tmin %5.4e %5.4e %5.4e\tmax %5.4e %5.4e %5.4e\n",i,s[i].l,s[i].o,s[i].d,s[i].id,s[i].min[0],s[i].min[1],s[i].min[2],s[i].max[0],s[i].max[1],s[i].max[2]);

}
void show_interactions(vector<interaction> ia, vector<long> n_ints, vector<long> o_ints)
{
	int nlim = ia.size();
	printf("**********\n");
	for(int i=0;i<nlim;i++)
	{
		printf("snap_A %04d\tsnap_B %04d\tidA %10ld\tidB %10lddA %5.4e\tdB %5.4e\n",ia[i].snap_A,ia[i].snap_B,ia[i].id_A,ia[i].id_B,ia[i].d_A,ia[i].d_B);
	}
}
void keep_duplicates(vector<long> iunion, vector<long> *ioverlap)
{
  vector<long>::iterator ia;

  ia = std::adjacent_find(iunion.begin(), iunion.end());
  if(ia!=iunion.end())
  {
    ioverlap->push_back(*ia);
    while(ia!=iunion.end())
    {
      ia = std::adjacent_find(++ia, iunion.end());
      if(ia!=iunion.end())
        ioverlap->push_back(*ia);
    }
  }
}
void write_branches(char fname[], vector<long> ia, vector<long> ib)
{
  FILE *fp;
  if(!(fp=fopen(fname,"w")))
  {
    printf("Error opening %s.\n",fname);
    exit(-1);
  }
  //fprintf(fp,"%ld\n",ia.size());
  for(size_t i=0;i<ia.size();i++)
  {
    fprintf(fp,"%ld\t%ld\n",ia[i],ib[i]);
  }
  fclose(fp);
}
void write_interactions(char fname[], vector<interaction> ia, float *times)
{
  FILE *fp;
  if(!(fp=fopen(fname,"w")))
  {
  	printf("Error opening %s.\n",fname);
  	exit(-1);
  }
  //fprintf(fp,"%ld\n",ia.size());
  for(size_t i=0;i<ia.size();i++)
  {
  	//fprintf(fp,"%04d %04d %9.8f %9.8f %8ld %8ld %8ld %5.4e %5.4e %10ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e %10ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e\n",ia[i].snap_A,ia[i].snap_B,times[ia[i].snap_A],times[ia[i].snap_B],ia[i].idx_A,ia[i].idx_B,ia[i].n,ia[i].frac_A,ia[i].frac_B,ia[i].id_A,ia[i].l_A,ia[i].o_A,ia[i].d_A,ia[i].x_A[0],ia[i].x_A[1],ia[i].x_A[2],ia[i].id_B,ia[i].l_B,ia[i].o_B,ia[i].d_B,ia[i].x_B[0],ia[i].x_B[1],ia[i].x_B[2]);
    fprintf(fp,"%04d %04d %9.8f %9.8f %8ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e %5.4e %10ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e %10ld %8ld %8ld %5.4e % 5.4e % 5.4e % 5.4e % 5.4e % 5.4e % 5.4e\n",ia[i].snap_A,ia[i].snap_B,times[ia[i].snap_A],times[ia[i].snap_B],ia[i].idx_A,ia[i].idx_B,ia[i].n,ia[i].frac_A,ia[i].frac_B,ia[i].frac_A_dense,ia[i].frac_B_dense,ia[i].frac_dense,ia[i].id_A,ia[i].l_A,ia[i].o_A,ia[i].d_A,ia[i].x_A[0],ia[i].x_A[1],ia[i].x_A[2],ia[i].id_B,ia[i].l_B,ia[i].o_B,ia[i].d_B,ia[i].x_B[0],ia[i].x_B[1],ia[i].x_B[2],ia[i].x_C[0],ia[i].x_C[1],ia[i].x_C[2]);

  }
  fclose(fp);
}
void read_interactions(char fname[], vector<interaction> *ia)
{
  FILE *fp;
  if(!(fp=fopen(fname,"r")))
  {
  	printf("Error opening %s.\n",fname);
  	exit(-1);
  }
  int snap_A;
  int snap_B;
  long idx_A;
  long id_A;
  long l_A;
  long o_A;
  float d_A;
  long idx_B;
  long id_B;
  long l_B;
  long o_B;
  float d_B;
  float frac_A;
  float frac_B;
  float frac_A_dense;
  float frac_B_dense;
  float frac_dense;

  float x_A[3];
  float x_B[3];
  float x_C[3];

  long n;
  fscanf(fp,"%ld\n",&n);
  ia->resize(n);
  for(size_t i=0;i<ia->size();i++)
  {
  	fscanf(fp,"%04d %04d %8ld %8ld %8ld %f %f %f %f %f %10ld %8ld %8ld %f %f %f %f %10ld %8ld %8ld %f %f %f %f %f %f %f\n",&snap_A,&snap_B,&idx_A,&idx_B,&n,&frac_A,&frac_B,&frac_A_dense,&frac_B_dense,&frac_dense,&id_A,&l_A,&o_A,&d_A,&x_A[0],&x_A[1],&x_A[2],&id_B,&l_B,&o_B,&d_B,&x_B[0],&x_B[1],&x_B[2],&x_C[0],&x_C[1],&x_C[2]);
  	(*ia)[i].n = n;
  	(*ia)[i].frac_A = frac_A;
  	(*ia)[i].frac_B = frac_B;
    (*ia)[i].frac_A_dense = frac_A_dense;
    (*ia)[i].frac_B_dense = frac_B_dense;
    (*ia)[i].frac_dense = frac_dense;
  	(*ia)[i].snap_A = snap_A;
  	(*ia)[i].snap_B = snap_B;
  	(*ia)[i].idx_A = idx_A;
  	(*ia)[i].id_A = id_A;
  	(*ia)[i].l_A = l_A;
  	(*ia)[i].o_A = o_A;
  	(*ia)[i].d_A = d_A;
  	(*ia)[i].idx_B = idx_B;
  	(*ia)[i].id_B = id_B;
  	(*ia)[i].l_B = l_B;
  	(*ia)[i].o_B = o_B;
  	(*ia)[i].d_B = d_B;
  	for(int k=0;k<3;k++)
  	{
  		(*ia)[i].x_A[k] = x_A[k];
  		(*ia)[i].x_B[k] = x_B[k];
      (*ia)[i].x_C[k] = x_C[k];
  	}
  }
  fclose(fp);
}
void read_interaction_counts(char fname[], vector<long> *n_ints, vector<long> *o_ints)
{
  long ntot;
  long nin, oin;
  FILE *fp;
  if(!(fp=fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(-1);
  }
  fscanf(fp,"%ld\n",&ntot);
  n_ints->resize(ntot);
  o_ints->resize(ntot);
  for(int i=0;i<ntot;i++)
  {
    fscanf(fp,"%ld\t%ld\n",&nin,&oin);
    (*n_ints)[i] = nin;
    (*o_ints)[i] = oin;
  }
  fclose(fp);
}

void write_interaction_counts(int iA, int iB, long n_total, vector<long> n_ints, vector<long> o_ints)
{
  FILE *fp;
  char fname[200];

  sprintf(fname,"interactions/interaction_count.%04d.%04d.txt",iA,iB);
  if(!(fp=fopen(fname,"w")))
  {
  	printf("Error opening %s.\n",fname);
  	exit(-1);
  }
  fprintf(fp,"%ld\n",n_total);
  for(int i=0;i<n_total;i++)
  	fprintf(fp,"%ld\t%ld\n",n_ints[i],o_ints[i]);
  fclose(fp);
}

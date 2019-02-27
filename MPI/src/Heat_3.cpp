#include <stdio.h>
#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

int main(int argc, char **argv)
{
	int i,j,t,ierr,num_grid,tsteps,proc_id,num_proc,start,end;
	double *u,*u_new,*final;
	double r=0.25,t_beg,t_end; //constant defined by Brian
	//MPI_Status stat;
	//argv[0] is program name
	t_beg=atoi(argv[1]);
	t_end=atoi(argv[2]);
	num_grid=atoi(argv[3]);
	tsteps=atoi(argv[4]);

	ierr=MPI_Init(&argc,&argv);
	ierr=MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);
	ierr=MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	int each=num_grid/num_proc;
	int extra=num_grid%num_proc;
	
if(num_proc==1)
{
	u=(double*)malloc((num_grid+2)*sizeof(double));
  	u_new=(double*)malloc((num_grid+2)*sizeof(double));
  	final=(double*)malloc((num_grid+2)*sizeof(double));

  	start=1;
	end=num_grid+1;

	u[end]=t_end;			//set_const
	u_new[end]=t_end;
		
	u[start-1]=t_beg;			//set_const
  	u_new[start-1]=t_beg;

  	for(i=start;i<end;i++)			//initialize everything else
	{
		u[i]=0.0;
		u_new[i]=0.0;
	}

	for(t=0;t<tsteps;t++)			//time loop
	{
		for(i=start;i<end;i++)			//compute
		{
			u_new[i]=(1-2*r)*u[i]+r*u[i-1]+r*u[i+1];
		}

		for(i=start;i<end;i++)			//transfer computed values
		{
			u[i]=u_new[i];
		}

		u[end]=t_end;
		u[start-1]=t_beg;
	}

	for(i=start;i<end;i++)
	{
		final[i]=u[i];
	}

	ofstream outf;
	outf.open("output.csv");
	for(i=1;i<num_grid;i++)
	{
		outf<<final[i]<<", ";
		cout<<final[i]<<" ";
	}
	outf<<final[num_grid];
	cout<<final[num_grid]<<"\n";
	outf.close();

	free(final);
}
else
{
	if(proc_id==num_proc-1)
	{
		u=(double*)malloc((each+extra+2)*sizeof(double));
  		u_new=(double*)malloc((each+extra+2)*sizeof(double));
  		// start=proc_id*each+1;
  		// end=(proc_id+1)*each+extra+1;
  		start=1;
  		end=each+extra+1;
  		u[end]=t_end;			//set_const
  		u_new[end]=t_end;
  		for(i=start-1;i<end;i++)			//initialize last set
		{
			u[i]=0.0;
			u_new[i]=0.0;
		}
	}
	else
	{
		u=(double*)malloc((each+2)*sizeof(double));
  		u_new=(double*)malloc((each+2)*sizeof(double));
  		// start=proc_id*each+1;
  		// end=(proc_id+1)*each+1;
  		start=1;
  		end=each+1;
  		if(proc_id==0)				
  		{
  			u[start-1]=t_beg;			//set_const
  			u_new[start-1]=t_beg;
  			for(i=start;i<end+1;i++)			//initialize first set
			{
				u[i]=0.0;
				u_new[i]=0.0;
			}

			final=(double*)malloc((num_grid+2)*sizeof(double));	//create final array
  		}
  		else
  		{
  			for(i=start-1;i<end+1;i++)			//initialize everything else
			{
				u[i]=0.0;
				u_new[i]=0.0;
			}
  		}
	}


	for(t=0;t<tsteps;t++)			//time loop
	{
		if(proc_id != num_proc-1)
		{
			//send rigth_val(at end-1) to R, tag=1
			MPI_Send(&u[end-1],1,MPI_DOUBLE,proc_id+1,1,MPI_COMM_WORLD);
		}

		if(proc_id != 0)
		{
			//rec right_val(into start-1) from L, tag=1
			MPI_Recv(&u[start-1],1,MPI_DOUBLE,proc_id-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			//send left_val(at start) to L, tag=2
			MPI_Send(&u[start],1,MPI_DOUBLE,proc_id-1,2,MPI_COMM_WORLD);
		}

		if(proc_id != num_proc-1)
		{
			//rec left_val(into end) from R, tag=2
			MPI_Recv(&u[end],1,MPI_DOUBLE,proc_id+1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}

		for(i=start;i<end;i++)			//compute
		{
			u_new[i]=(1-2*r)*u[i]+r*u[i-1]+r*u[i+1];
		}

		for(i=start;i<end;i++)			//transfer computed values
		{
			u[i]=u_new[i];
		}

		if(proc_id==0)					//set_const again
		{
			u[start-1]=t_beg;
		}

		if(proc_id==num_proc-1)
		{
			u[end]=t_end;
		}

	}

	if(proc_id != 0)
	{
		//send from start to end, to proc 0, tag=3
		for(i=start;i<end;i++)	
		{
			MPI_Send(&u[i],1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
	}
	else
	{	
		for(i=start;i<end;i++)
		{
			final[i]=u[i];
		}

		for(j=1;j<num_proc;j++)		//rec from all proc
		{
			if(j==num_proc-1)
			{
				start=j*each+1;
				end=(j+1)*each+extra+1;
			}
			else
			{
				start=j*each+1;
				end=(j+1)*each+1;
			}

			for(i=start;i<end;i++)
			{
				//rec from all grids, tag=3
				MPI_Recv(&final[i],1,MPI_DOUBLE,j,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}

		//write to csv
		// cout<<"\nOUTPUT:\n";
		// for(i=1;i<num_grid+1;i++)
		// {
		// 	cout<<" "<<final[i]<<" ";
		// }

		ofstream outf;
		outf.open("output.csv");
		for(i=1;i<num_grid;i++)
		{
			outf<<final[i]<<", ";	
		}
		outf<<final[num_grid];
		outf.close();

		free(final);
	}
}

	free(u);
	free(u_new);

	MPI_Finalize();
	return 0;
}
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<cstdio>
#include<vector>
#include<stdlib.h>
#include<string.h>
#include<cstring>
#include<sys/time.h>

using namespace std;

const int TPB=1024;
int NX;
int NY;
int NZ;
int tsteps;
float l;
int dim;
float init_t;

void init(float *T)
{
	if(dim==2)
	{
		for(int j=0;j<NY;j++)
			for(int i=0;i<NX;i++)
			{
				int index=i+j*NX;
				T[index]=init_t;
			}
	}

	else if(dim==3)
	{
		for(int k=0;k<NZ;k++)
		{
			for(int j=0;j<NY;j++)
			{
				for(int i=0;i<NX;i++)
				{
					int index=i+j*NX+k*(NX*NY);
					T[index]=init_t;
				}
			}
		}
	}
}

 void set_const(float *T, float *T_const)
 {
 	if(dim==2)
 	{
 	for(int j=0;j<NY;j++)
		for(int i=0;i<NX;i++)
		{
			int index=i+j*NX;
			if(T_const[index]!=(-1))
				T[index]=T_const[index];
		}
 	}
 	else if(dim==3)
 	{
 		for(int k=0;k<NZ;k++)
		{
			for(int j=0;j<NY;j++)
			{
				for(int i=0;i<NX;i++)
				{
					int index=i+j*NX+k*(NX*NY);
					if(T_const[index]!=(-1))
						T[index]=T_const[index];
				}
			}
		}
 	}	
 }

// __global__ void print_kernel(float *T_old, float *T_new,float *T_const_GPU,int *dim_GPU,float *l_GPU,int *NX_GPU,int *NY_GPU,int *NZ_GPU)
// {

// 		//printf("\ndim_GPU: %d, l_GPU: %f, NX_GPU: %d, NY_GPU: %d, NZ_GPU: %d\n",*dim_GPU,*l_GPU,*NX_GPU,*NY_GPU,*NZ_GPU);
// 		printf("\nT1 - %d")
	
// }

__global__ void swap_kernel(float *T_old, float *T_new,int *dim_GPU,int *NX_GPU,int *NY_GPU,int *NZ_GPU)
{

	//T1=T2
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(*dim_GPU==2)
	{
		if(index<(*NX_GPU)*(*NY_GPU))
		{
			// printf("\nT_old before swap:\n");
			// printf("index: %d, T_old[index]: %d",index,T_old[index]);

			T_old[index]=T_new[index];

			// printf("\nT_old after swap:\n");
			// printf("index: %d, T_old[index]: %d",index,T_old[index]);
		}
	}

	else if(*dim_GPU==3)
	{
		if(index<(*NX_GPU)*(*NY_GPU)*(*NZ_GPU))
		{
			// printf("\nINSIDE 3D SWAP\n");
			// printf("\nT_old before swap:\n");
			// printf("index: %d, T_old[index]: %d",index,T_old[index]);

			T_old[index]=T_new[index];

			// printf("\nT_old after swap:\n");
			// printf("index: %d, T_old[index]: %d",index,T_old[index]);
		}
	}
}

 __global__ void update_const(float *T, float *T_const_GPU, int *dim_GPU, int *NX_GPU,int *NY_GPU,int *NZ_GPU)
 {
 	if(*dim_GPU==2)
 	{
	 	int index = threadIdx.x + blockIdx.x * blockDim.x;
		 if(index<(*NX_GPU)*(*NY_GPU))
		 {
		 	//printf("index: %d, %f", index, T[index]);
		 	if(T_const_GPU[index]!=(-1))
		 	{
				T[index]=T_const_GPU[index];
				//printf("index: %d, %f ---", index, T_const_GPU[index]);
		 	}
	 	}
	}
	else if(*dim_GPU==3)
	{
		int index = threadIdx.x + blockIdx.x * blockDim.x;	
		if(index<(*NX_GPU)*(*NY_GPU)*(*NZ_GPU))
		{
		if(T_const_GPU[index]!=(-1))
			T[index]=T_const_GPU[index];
		}
	}
 }

__global__ void update(int *t_GPU,float *T_old, float *T_new,float *T_const_GPU,int *dim_GPU,float *l_GPU,int *NX_GPU,int *NY_GPU,int *NZ_GPU)
{
	int index=0,left,right,top,bott,front,back;
	// printf("\ntimestep: %d\n",(*t_GPU));
	index = threadIdx.x + blockIdx.x * blockDim.x;

	if(*dim_GPU==2)
	{
		if(index<((*NX_GPU)*(*NY_GPU)))
		{
			if(index % (*NX_GPU) == 0)
				left=index;
			else
				left=index-1;

			if((index + 1) % (*NX_GPU) == 0)
				right=index;
			else
				right=index+1;

			if(index / (*NX_GPU) == 0)
				top=index;
			else
				top=index-(*NX_GPU);

			if(((((*NX_GPU)*(*NY_GPU))-index-1) / (*NX_GPU)) == 0)
				bott=index;
			else
				bott=index+(*NX_GPU);

			//printf("index: %d, right: %f:, left: %f, top: %f, bottom: %f\n", index, T_old[right], T_old[left], T_old[top], T_old[bott]);

			T_new[index]=T_old[index] + (*l_GPU) *(T_old[left] + T_old[right] + T_old[top] + T_old[bott] - (4*T_old[index]));
		}
	}
	else if(*dim_GPU==3)
	{				
		if(index<(*NX_GPU)*(*NY_GPU)*(*NZ_GPU))
		{
			if(index % *NX_GPU == 0)
				left=index;
			else
				left=index-1;

			if((index + 1) % *NX_GPU == 0)
				right=index;
			else
				right=index+1;

			if(index / ((*NX_GPU) * (*NY_GPU)) == 0)
				front=index;
			else
				front=index-((*NX_GPU)*(*NY_GPU));

			if(( ((*NX_GPU)*(*NY_GPU)*(*NZ_GPU)) - index - 1) / ((*NX_GPU)*(*NY_GPU)) == 0)
				back=index;
			else
				back=index+((*NX_GPU)*(*NY_GPU));

			if((index%((*NX_GPU)*(*NY_GPU)))/(*NX_GPU) == 0)
				top=index;
			else
				top=index-(*NX_GPU);

			if((( ((*NX_GPU)*(*NY_GPU)*(*NZ_GPU))-index-1 ) % ((*NX_GPU)*(*NY_GPU))) / *NX_GPU == 0)
				bott=index;
			else
				bott=index+(*NX_GPU);

			T_new[index]=T_old[index] + (*l_GPU) *(T_old[left] + T_old[right] + T_old[top] + T_old[bott] + T_old[front] + T_old[back] - 6*T_old[index]);
		}
	}
}


int main(int argc,char **argv)
{
	int i,j,k,index,t,num_block,ctr=0;
	float *T_const;
	char line[500];
	char *token;
	char *filename=argv[1];
	FILE * outf=fopen("heatOutput.csv","w");

	FILE * file;
	file=fopen(filename,"r");
	while(fgets(line,sizeof(line),file))
	{
		if(line[0]!='#')
		{
			ctr++;
			if(ctr==1)
			{
				dim=line[0]-'0';
				//cout<<"\ndim: "<<dim;
			}

			if(ctr==2)
			{
				token=strtok(line,"\n");
				l=atof(token);
				//cout<<"\nl: "<<l;
			}

			if(ctr==3)
			{
				token=strtok(line,"\n");
				tsteps=atoi(token);
				//cout<<"\ntsteps: "<<tsteps;
			}

			if(ctr==4)
			{
				if(dim==2)
				{
					token=strtok(line,",");
					if(token!=NULL)
					{
						NX=atoi(token);
						//cout<<"\n2D NX: "<<NX;
						token=strtok(NULL,"\n");
					}

					if(token!=NULL)
					{
						NY=atoi(token);
						//cout<<"\n2D NY: "<<NY;
					}

						NZ=0;

					T_const=new float[NX*NY];
					for(j=0;j<NY;j++)	//Initializing T_const to -1
						for(i=0;i<NX;i++)
						{
							index=i+j*NX;
							T_const[index]=-1;
						}
				}
				else if(dim==3)
				{
					token=strtok(line,",");
					if(token!=NULL)
					{
						NX=atoi(token);
						//cout<<"\n3D NX: "<<NX;
						token=strtok(NULL,",");
					}
					if(token!=NULL)
					{
						NY=atoi(token);
						//cout<<"\n3D NY: "<<NY;
						token=strtok(NULL,"\n");
					}
					if(token!=NULL)
					{
						NZ=atoi(token);
						//cout<<"\n3D NZ: "<<NZ;
					}
					T_const=new float[NX*NY*NZ];
					for(k=0;k<NZ;k++)	//Initializing T_const to -1
					{
						for(j=0;j<NY;j++)
						{
							for(i=0;i<NX;i++)
							{
								index=i+j*NX+k*(NX*NY);
								T_const[index]=-1;
							}
						}
					}
				}
			}

			if(ctr==5)
			{
				token=strtok(line,"\n");
				init_t=atof(token);
				//cout<<"\ninit_t: "<<init_t;
			}

			if(ctr>5)
			{
				if(dim==2)
				{
					int x2,y2,w2,h2;
					float t2;

					token=strtok(line,",");
					if(token!=NULL)
						{ x2=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ y2=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ w2=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ h2=atoi(token); token=strtok(NULL,"\n"); }
					if(token!=NULL)
						{ t2=atof(token); }

					for(j=y2;j<y2+h2;j++)	//creating T_const matrix 
						for(i=x2;i<x2+w2;i++)
						{
							index=i+j*NX;
							T_const[index]=t2;
						}	

				}
				else if(dim==3)
				{
					int x3,y3,z3,w3,h3,d3;
					float t3;

					token=strtok(line,",");
					if(token!=NULL)
						{ x3=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ y3=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ z3=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ w3=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ h3=atoi(token); token=strtok(NULL,","); }
					if(token!=NULL)
						{ d3=atoi(token); token=strtok(NULL,"\n"); }
					if(token!=NULL)
						{ t3=atof(token); }

					for(k=z3;k<z3+d3;k++)	//Create const matrix
					{
						for(j=y3;j<y3+h3;j++)
						{
							for(i=x3;i<x3+w3;i++)
							{
								index=i+j*NX+k*(NX*NY);
								T_const[index]=t3;
							}
						}
					}
				}
			}
		}
	}
				

	if(dim==2)
	{
		float *T=new float[NX*NY];	 		//CPU Memory
		float *T1, *T2, *T_const_GPU;		//GPU Memory
		//float *T_temp=new float[NX*NY];
		int *dim_CPU,*NX_CPU,*NY_CPU,*NZ_CPU,*t_CPU;	//t_CPU and t_GPU created only to print timesteps from kernel
		float *l_CPU;

		dim_CPU=&dim; l_CPU=&l; NX_CPU=&NX; NY_CPU=&NY; NZ_CPU=&NZ; 

		int *dim_GPU,*NX_GPU,*NY_GPU,*NZ_GPU,*t_GPU;
		float *l_GPU;

		init(T);

		set_const(T,T_const);	//copying constants from T_const to T

		cudaMalloc((void**)&T1,NX*NY*sizeof(float));	//T_old
		cudaMalloc((void**)&T2,NX*NY*sizeof(float));	//T_new
		cudaMalloc((void**)&T_const_GPU,NX*NY*sizeof(float)); //T_const

		cudaMalloc((void**)&dim_GPU,sizeof(int));
		cudaMalloc((void**)&t_GPU,sizeof(int));
		cudaMalloc((void**)&l_GPU,sizeof(float));
		cudaMalloc((void**)&NX_GPU,sizeof(int));
		cudaMalloc((void**)&NY_GPU,sizeof(int));
		cudaMalloc((void**)&NZ_GPU,sizeof(int));

		cudaMemcpy(dim_GPU,dim_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(l_GPU,l_CPU,sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(NX_GPU,NX_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(NY_GPU,NY_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(NZ_GPU,NZ_CPU,sizeof(int),cudaMemcpyHostToDevice);

		cudaMemcpy(T1,T,NX*NY*sizeof(float),cudaMemcpyHostToDevice);
		//cudaMemcpy(T2,T,NX*NY*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(T_const_GPU,T_const,NX*NY*sizeof(float),cudaMemcpyHostToDevice);

		num_block=((NX*NY)+TPB-1)/TPB;	//TPB = Threads Per Block = 1024

		//print_kernel<<<1,1>>>(T1,T2,T_const_GPU,dim_GPU,l_GPU,NX_GPU,NY_GPU,NZ_GPU);

		//struct timeval start,end;
		//double diff;
		//gettimeofday(&start,NULL);

		for(t=0;t<tsteps;t++)
		{
			t_CPU=&t;
			cudaMemcpy(t_GPU,t_CPU,sizeof(int),cudaMemcpyHostToDevice);
	
			update<<<num_block,TPB>>>(t_GPU,T1,T2,T_const_GPU,dim_GPU,l_GPU,NX_GPU,NY_GPU,NZ_GPU);



			update_const<<<num_block,TPB>>>(T2,T_const_GPU,dim_GPU,NX_GPU,NY_GPU,NZ_GPU);
			
			// cudaMemcpy(T_temp, T2, NX*NY*sizeof(float), cudaMemcpyDeviceToHost);	//SWAP IN CPU IMPLEMENTATION
			// cout<<"\nT_TEMP:\n";
			// for(int j=0;j<NY;j++)
			// {
			// 	for(int i=0;i<NX;i++)
			// 	{
			// 		int index=i+j*NX;
			// 		cout<<" "<<T_temp[index]<<" ";
			// 	}
			// 	cout<<"\n";
			// }
			//cudaMemcpy(T1, T_temp, NX*NY*sizeof(float), cudaMemcpyHostToDevice);
			//update_const<<<num_block,TPB>>>(T1,T_const_GPU,dim_GPU);
			

			swap_kernel<<<num_block,TPB>>>(T1,T2,dim_GPU,NX_GPU,NY_GPU,NZ_GPU);	//KERNEL SWAP IMPLEMENTATION


			//T1=T2;	//SIMPLE(WRONG?) SWAP IMPLEMENTATION

			//DOUBLE UPDATE CALL IMPLEMENTATION
		}

	//	gettimeofday(&end,NULL);
	//	diff=(end.tv_sec - start.tv_sec)*1000.0;	//sec to ms conversion

	//	cout<<"\n2D TIME: "<<diff<<"\n";

		// if(tsteps%2==0)	//FOR DOUBLE UPDATE CALL IMPLEMENTATION
		// {
			cudaMemcpy(T,T1,NX*NY*sizeof(float),cudaMemcpyDeviceToHost);
		// }
		// else
		// {
		// 	cudaMemcpy(T,T1,NX*NY*sizeof(float),cudaMemcpyDeviceToHost);
		// }
		
		//cout<<"\nFINAL:\n";

		for(int j=0;j<NY;j++)
		{
			for(int i=0;i<NX;i++)
			{
				index=i+j*NX;
				if((index+1) % NX == 0)
				{
					fprintf(outf,"%f",T[index]);
					printf("%f",T[index]);
				}
				else
				{
					fprintf(outf,"%f, ",T[index]);
					printf("%f",T[index]);
				}
			}
			fprintf(outf,"\n");
			printf("\n");
		}

		delete T;
		delete T_const;

		cudaFree(T1);
		cudaFree(T2);
		cudaFree(T_const_GPU);
		cudaFree(l_GPU);
		cudaFree(dim_GPU);
		cudaFree(NX_GPU);
		cudaFree(NY_GPU);
		cudaFree(NZ_GPU);
		cudaFree(t_GPU);
	}

	else if(dim==3)
	{
		float *T=new float[NX*NY*NZ];	 		//CPU Memory
		float *T1, *T2, *T_const_GPU;		//GPU Memory
		//float *T_temp=new float[NX*NY];
		int *dim_CPU,*NX_CPU,*NY_CPU,*NZ_CPU,*t_CPU;	//t_CPU and t_GPU created only to print timesteps from kernel
		float *l_CPU;

		dim_CPU=&dim; l_CPU=&l; NX_CPU=&NX; NY_CPU=&NY; NZ_CPU=&NZ; 

		int *dim_GPU,*NX_GPU,*NY_GPU,*NZ_GPU,*t_GPU;
		float *l_GPU;

		init(T);

		set_const(T,T_const);	//copying constants from T_const to T

		cudaMalloc((void**)&T1,NX*NY*NZ*sizeof(float));	//T_old
		cudaMalloc((void**)&T2,NX*NY*NZ*sizeof(float));	//T_new
		cudaMalloc((void**)&T_const_GPU,NX*NY*NZ*sizeof(float));

		cudaMalloc((void**)&dim_GPU,sizeof(int));
		cudaMalloc((void**)&l_GPU,sizeof(float));
		cudaMalloc((void**)&NX_GPU,sizeof(int));
		cudaMalloc((void**)&NY_GPU,sizeof(int));
		cudaMalloc((void**)&NZ_GPU,sizeof(int));
		cudaMalloc((void**)&t_GPU,sizeof(int));

		cudaMemcpy(dim_GPU,dim_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(l_GPU,l_CPU,sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(NX_GPU,NX_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(NY_GPU,NY_CPU,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(NZ_GPU,NZ_CPU,sizeof(int),cudaMemcpyHostToDevice);

		cudaMemcpy(T1,T,NX*NY*NZ*sizeof(float),cudaMemcpyHostToDevice);
		//cudaMemcpy(T2,T,NX*NY*NZ*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(T_const_GPU,T_const,NX*NY*NZ*sizeof(float),cudaMemcpyHostToDevice);

		num_block=((NX*NY*NZ)+TPB-1)/TPB;	//TPB = Threads Per Block = 1024

	//	struct timeval start,end;
	//	double diff;
	//	gettimeofday(&start,NULL);

		for(t=0;t<tsteps;t++)
		{
			t_CPU=&t;
			cudaMemcpy(t_GPU,t_CPU,sizeof(int),cudaMemcpyHostToDevice);

			update<<<num_block,TPB>>>(t_GPU,T1,T2,T_const_GPU,dim_GPU,l_GPU,NX_GPU,NY_GPU,NZ_GPU);
			update_const<<<num_block,TPB>>>(T2,T_const_GPU,dim_GPU,NX_GPU,NY_GPU,NZ_GPU);

			//update_const<<<num_block,TPB>>>(T1,T_const_GPU,dim_GPU);
			//swap
			//T1=T2;

			//swap_kernel<<<num_block,TPB>>>(T1,T2,dim_GPU,NX_GPU,NY_GPU,NZ_GPU);

			// update<<<num_block,TPB>>>(t_GPU,T2,T1,T_const_GPU,dim_GPU,l_GPU,NX_GPU,NY_GPU,NZ_GPU);
			// update_const<<<num_block,TPB>>>(T1,T_const_GPU,dim_GPU);
			// update_const<<<num_block,TPB>>>(T2,T_const_GPU,dim_GPU);
		}

	//	gettimeofday(&end,NULL);
	//	diff=(end.tv_sec - start.tv_sec)*1000.0;	//sec to ms conversion

	//	cout<<"\n3D TIME: "<<diff<<"\n";

		// if(tsteps%2==0)
		// {
			cudaMemcpy(T,T2,NX*NY*NZ*sizeof(float),cudaMemcpyDeviceToHost);
		// }
		// else
		// {
		// 	cudaMemcpy(T,T1,NX*NY*NZ*sizeof(float),cudaMemcpyDeviceToHost);
		// }
		
		/*for(int i=0;i<NX*NY*NZ;i++)
		{
		if(i%NX==0)
			fprintf(outf,"\n");
		if(i%(NX*NY)==0)
			fprintf(outf,"\n");
		if(i%NX==NX-1)
			fprintf(outf,"%f ",T[i]);
		else
			fprintf(outf,"%f, ",T[i]);
		}*/

		for(int k=0;k<NZ;k++)
		{
			for(int j=0;j<NY;j++)
			{
				for(int i=0;i<NX;i++)
				{
					index=i+j*NX+k*(NX*NY);
					if((index+1) % NX == 0)
					{
						fprintf(outf,"%f",T[index]);
						printf("%f",T[index]);
					}
					else
					{
						fprintf(outf,"%f, ",T[index]);
						printf("%f",T[index]);
					}
				}
				fprintf(outf,"\n");
				printf("\n");
			}
			fprintf(outf,"\n");
			printf("\n");
		}

		delete T;
		delete T_const;

		cudaFree(T1);
		cudaFree(T2);
		cudaFree(T_const_GPU);
		cudaFree(l_GPU);
		cudaFree(dim_GPU);
		cudaFree(NX_GPU);
		cudaFree(NY_GPU);
		cudaFree(NZ_GPU);
		cudaFree(t_GPU);
	}

		return 0;
}

	//GPU JOB: qsub -I -q coc-ice -l nodes=1:ppn=2:gpus=1,walltime=2:00:00,pmem=2gb
	//CPU JOB: qsub -I -q coc-ice -l nodes=1:ppn=12,walltime=2:00:00,pmem=2gb
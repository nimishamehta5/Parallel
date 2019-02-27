#include<iostream>
#include<thread>
#include<fstream>
#include<stdio.h>
#include<cstdio>
#include<stdlib.h>
#include<string.h>
#include<cstring>
#include<ctime>
#include<chrono>
#include<sys/time.h>
using namespace std;

int NX;
int NY;
int NZ;
int tsteps;
float l;
int dim;
float init_t;

void set_const(float *T, float *T_const);
void thread_compute(int t_num, float *T, float *T_temp);
void init(float *T);

int main(int argc,char **argv)
{
	int i,j,k,index,t,num_block,ctr=0;
	float *T_const;
	char line[500];
	char *token;
	char *filename=argv[1];
	float *T;
	float *T_temp;
	float *T_swap;
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
					for(j=0;j<NY;j++)	//Initializing T_const to 0.0
						for(i=0;i<NX;i++)
						{
							index=i+j*NX;
							T_const[index]=(-1);
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
					for(k=0;k<NZ;k++)	//Initializing T_const to 0.0
					{
						for(j=0;j<NY;j++)
						{
							for(i=0;i<NX;i++)
							{
								index=i+j*NX+k*(NX*NY);
								T_const[index]=(-1);
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

	std::thread threads[11];

	if(dim==2)
	{
		T=new float[NX*NY];
		T_temp=new float[NX*NY];	
		T_swap=new float[NX*NY];
	}
	else if(dim==3)
	{
		T=new float[NX*NY*NZ];
		T_temp=new float[NX*NY*NZ];
		T_swap=new float[NX*NY*NZ];
	}

	int m_start,m_end,m_interval,m,left,right,top,bott;

	init(T);

	set_const(T,T_const);

	if(dim==2)
		{
			m_interval=(NX*NY)/12;
		}
	else if(dim==3)
		{
			m_interval=(NX*NY*NZ)/12;
		}

	//m_start=0; m_end=m_interval;

	//  time_t start,end;
	//  double sec;
	
	// start=time(NULL);

	//std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

	// struct timeval start,end;
	// double diff;
	// gettimeofday(&start,NULL);

	for(t=0;t<tsteps;t++)
	{
		// printf("\n******** TIMESTEP ********: %d\n",t);
		thread_compute(1,T,T_temp);

		for(i=1;i<12;i++)
		{
			threads[i]=std::thread(thread_compute,i+1,T,T_temp);
		}

		for(i=1;i<12;i++)
		{
			threads[i].join();
		}

		T_swap=T;
		T=T_temp;
		T_temp=T_swap;

		set_const(T,T_const);
		set_const(T_temp,T_const);
	}

	// gettimeofday(&end,NULL);
	// diff=(end.tv_sec - start.tv_sec)*1000.0;	//sec to ms conversion

	// cout<<"\n2D TIME: "<<diff<<"\n";

	 //end=time(NULL);

	 //sec=difftime(end,start);
	//cout<<"TIME: "<<sec;

	// std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	// std::chrono::duration<double> difference = (end - start);
	// cout<<"\nTIME: "<<difference.count();


	if(dim==2)
	{
		// cout<<"\nFINAL:\n";

		// for(int j=0;j<NY;j++)
		// {
		// 	for(int i=0;i<NX;i++)
		// 	{
		// 		int index=i+j*NX;
		// 		cout<<" "<<T[index]<<" ";
		// 	}
		// 	cout<<"\n";
		// }

		for(int j=0;j<NY;j++)
		{
			for(int i=0;i<NX;i++)
			{
				index=i+j*NX;
				if((index+1) % NX == 0)
				{
					fprintf(outf,"%f",T[index]);
				}
				else
				{
					fprintf(outf,"%f, ",T[index]);
				}
			}
			fprintf(outf,"\n");
		}
	}

	else if(dim==3)
	{
	// 	cout<<"\nFINAL:\n";

	// 	for(int k=0;k<NZ;k++)
	// 	{
	// 		for(int j=0;j<NY;j++)
	// 		{
	// 			for(int i=0;i<NX;i++)
	// 			{
	// 				int index=i+j*NX;
	// 				cout<<" "<<T[index]<<" ";
	// 			}
	// 			cout<<"\n";
	// 		}
	// 		cout<<"\n";
	// 	}
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
					}
					else
					{
						fprintf(outf,"%f, ",T[index]);
					}
				}
				fprintf(outf,"\n");
			}
			fprintf(outf,"\n");
		}
	 }
	
}

	
void thread_compute(int t_num, float *T, float *T_temp)
{
	int start,end,interval=0,extra=0,i;
	int left,right,top,bott,front,back;
	// std::thread::id t_id=std::this_thread::get_id();
	// std::cout << "\nThis is the %d thread.\n",t_id;
	if(dim==2)
	{
		interval=(NX*NY)/12;
		extra=(NX*NY)%12;	
	}
	else if(dim==3)
	{
		interval=(NX*NY*NZ)/12;
		extra=(NX*NY*NZ)%12;
	}

	//compute T_temp from T

	//printf("\nThread number: %d\n",t_num);

	start=(t_num-1)*interval;
	end=(t_num)*interval;

	if(t_num==12)	{end+=extra;}

	//printf("\nStart: %d, End: %d\n",start,end);
	
	if(dim==2)
	{
		for(i=start;i<end;i++)
		{
			if(i % (NX) == 0)
				left=i;
			else
				left=i-1;

			if((i + 1) % (NX) == 0)
				right=i;
			else
				right=i+1;

			if(i / (NX) == 0)
				top=i;
			else
				top=i-(NX);

			if(((((NX)*(NY))-i-1) / (NX)) == 0)
				bott=i;
			else
				bott=i+(NX);

			//printf("\nleft,right,top,bottom: %f, %f, %f, %f\n",T[left],T[right],T[top],T[bott]);

			T_temp[i]=T[i] + (l) *(T[left] + T[right] + T[top] + T[bott] - (4*T[i]));

			//printf("\nT_temp[%d]: %f\n",i,T_temp[i]);
		}
	}
	else if(dim==3)
	{
		for(i=start;i<end;i++)
		{
			if(i % NX == 0)
				left=i;
			else
				left=i-1;

			if((i + 1) % NX == 0)
				right=i;
			else
				right=i+1;

			if(i / ((NX) * (NY)) == 0)
				front=i;
			else
				front=i-((NX)*(NY));

			if(( ((NX)*(NY)*(NZ)) - i - 1) / ((NX)*(NY)) == 0)
				back=i;
			else
				back=i+((NX)*(NY));

			if((i%((NX)*(NY)))/(NX) == 0)
				top=i;
			else
				top=i-(NX);

			if((( ((NX)*(NY)*(NZ))-i-1 ) % ((NX)*(NY))) / NX == 0)
				bott=i;
			else
				bott=i+(NX);

			T_temp[i]=T[i] + (l) *(T[left] + T[right] + T[top] + T[bott] + T[front] + T[back] - 6*T[i]);
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
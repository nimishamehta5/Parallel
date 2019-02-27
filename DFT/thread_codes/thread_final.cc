#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include"complex.cc"
#include"input_image.cc"
//#include"input_image.h"
#include<thread>
#include<mutex>
#include<cstring>
#include<time.h>
#include<sys/time.h>

void compute_dft(Complex *D,Complex *D_calc,int t_id,int n,int nt);
void transpose(Complex *D_in,Complex *D_out,int n);
void compute_idft(Complex *ID,Complex *ID_calc,int t_id,int n,int nt);

void compute_dft(Complex *D,Complex *D_calc,int t_id,int n,int nt)
{
  int start,end,k,x,row_st,row_end;
  float s,c;

  x=n/nt;
Complex sum;
  row_st=t_id*x;
  row_end=(t_id+1)*x;

  for(k=row_st;k<row_end;k++)
  {
  	start=k*n;
  	end=start+n-1;

	  for(int i=start, i1 = 0;i<=end, i1 < n;i++, i1++)
	  {
        sum.real = 0;
        sum.imag= 0;
	     for(int j=start, j1 =0;j<=end, j1<n;j++, j1++)
	     {
	        c=cos(-2.*PI*i1*j1/n);
	        s=sin(-2.*PI*i1*j1/n);
	        sum.real=sum.real+(c*D[j].real-s*D[j].imag);
	        sum.imag=sum.imag+(s*D[j].real+D[j].imag*c);
	     }
	     D_calc[i] = sum;
	  }
  }

}

void compute_dft_row(Complex *D,Complex *D_calc,int t_id,int n,int nt)
{
  int start,end,i,j,k,x,row_st,row_end;
  Complex sum;
  float s,c;

  x=n/nt;

  row_st=t_id*x;
  row_end=(t_id+1)*x;

  for(k=row_st;k<row_end;k++)
  {
  	start=k*n;
  	end=start+n;

	  for(int i=start, i1=0;i<= (start+end)/2, i1<=n/2;i++, i1++)
	  {
	     sum.real=0.0;
	     sum.imag=0.0;
	     for(int j=start, j1 =0;j<=end, j1<n;j++, j1++)
	     {
	        c=cos(-2.*PI*i1*j1/n);
	        s=sin(-2.*PI*i1*j1/n);
	        sum.real=sum.real+(c*D[j].real-s*D[j].imag);
	        sum.imag=sum.imag+(s*D[j].real+D[j].imag*c);
	     }
	     D_calc[i] = sum;
	  }

	  for(i = ((start+end)/2) + 1; i < end; i++){
        D_calc[i] = D_calc[(end)-(i-start)].conj();
	  }
  }

}

void compute_idft(Complex *ID,Complex *ID_calc,int t_id,int n,int nt)
{
  int start,end,i,j,k,x,row_st,row_end;
  Complex sum;

  x=n/nt;

  row_st=t_id*x;
  row_end=(t_id+1)*x;

  float c,s,n_f;
  n_f=(float)n;

  for(k=row_st;k<row_end;k++)
  {
  	start=(t_id*x*n)+((k%x)*n);
  	end=start+n-1;

	  for(i=start;i<=end;i++)
	  {
	     sum.real=0.0;
	     sum.imag=0.0;
	     for(j=start;j<=end;j++)
	     {
	        c=cos(2*PI*(i%n)*(j%n)/n);
	        s=sin(2*PI*(i%n)*(j%n)/n);
	        //Complex factor(c,s);
	        sum.real=sum.real+(c*ID[j].real-s*ID[j].imag);
	        sum.imag=sum.imag+(c*ID[j].imag+s*ID[j].real);
	        //sum=sum+(factor*ID[j]);
	     }
	     ID_calc[i].real=sum.real/n_f;
	     ID_calc[i].imag=sum.imag/n_f;
	  }
  }
}

void transpose(Complex *D_in,Complex *D_out,int n)
{
	int i,j,rhs,lhs;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			rhs=i+j*n;
			lhs=j+i*n;
			D_out[lhs]=D_in[rhs];
		}
	}
}

int main(int argc,char **argv)
{


	int n,i,j,nt;
	char* type;
	type=argv[1];
	//cout<<type;

	nt=16;
	//nt=2;

	InputImage input_file(argv[2]);
	n = input_file.get_height();

	if(n < nt) nt = n;


	if(strcmp(type,"forward")==0)
	{
		//n=2;

		Complex *D_in=(Complex*)malloc(sizeof(Complex)*n*n);
		Complex *D_out=(Complex*)malloc(sizeof(Complex)*n*n);

		D_in=input_file.get_image_data();

		// D_in[0]=1;
		// D_in[1]=2;
		// D_in[2]=3;
		// D_in[3]=4;


		std::thread threads[nt];
		for(i=0;i<nt;i++)
		{
			threads[i]=std::thread(compute_dft_row,D_in,D_out,i,n,nt);
		}

		for(i=0;i<nt;i++)
		{
			threads[i].join();
		}



		transpose(D_out,D_in,n);



		for(int i=0;i<nt;i++)
		{
			threads[i]=std::thread(compute_dft,D_in,D_out,i,n,nt);
		}

		for(i=0;i<nt;i++)
		{
			threads[i].join();
		}



		// cout<<"\nFINAL:\n";
		// for(int i = 0; i < n; i++)
		// {
		// 	for(int j = 0; j < n;j++)
		// 	{
		// 		cout<<D_out[i*n+j].real;
		// 	}
		// 	cout<<"\n";
		// }

		transpose(D_out,D_in,n);



		input_file.save_image_data(argv[3],D_in,n,n);



		// for(int i = 0; i < n; i++)
		// {
		// 	for(int j = 0; j < n;j++)
		// 	{
		// 		fprintf(fp,"(%f,%f) ",D_out[i*n+j].real,D_out[i*n+j].imag);
		// 	}
		// 	fprintf(fp, "\n");
		// }

	}

	else if(strcmp(type,"reverse")==0)
	{
		//n=2;
		Complex *ID_in=(Complex*)malloc(sizeof(Complex)*n*n);
		Complex *ID_out=(Complex*)malloc(sizeof(Complex)*n*n);

		ID_in=input_file.get_image_data();

		// ID_in[0].real=10;
		// ID_in[1].real=-4;
		// ID_in[2].real=-2;
		// ID_in[3].real=0;
		// ID_in[0].imag=0;
		// ID_in[1].imag=0;
		// ID_in[2].imag=0;
		// ID_in[3].imag=0;

		std::thread threads[nt];

		for(i=0;i<nt;i++)
		{
			threads[i]=std::thread(compute_idft,ID_in,ID_out,i,n,nt);
		}

		for(i=0;i<nt;i++)
		{
			threads[i].join();
		}

		transpose(ID_out,ID_in,n);

		for(int i=0;i<nt;i++)
		{
			threads[i]=std::thread(compute_idft,ID_in,ID_out,i,n,nt);
		}

		for(i=0;i<nt;i++)
		{
			threads[i].join();
		}

		// cout<<"\nFINAL:\n";
		// for(i=0;i<n*n;i++)
		// {
		// 	cout<<ID_out[i].real<<"+"<<ID_out[i].imag<<"i"<<"\t";
		// }

		transpose(ID_out,ID_in,n);

		input_file.save_image_data(argv[3],ID_in,n,n);

		// for(int i = 0; i < n; i++)
		// {
		// 	for(int j = 0; j < n;j++)
		// 	{
		// 		fprintf(fp,"(%f,%f) ",ID_out[i*n+j].real,ID_out[i*n+j].imag);
		// 	}
		// 	fprintf(fp, "\n");
		// }
	}
}

//qsub -I -q coc-ice -l nodes=1:ppn=8,walltime=2:00:00,pmem=2gb

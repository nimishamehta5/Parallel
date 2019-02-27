

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cuda.h>
#include "complex.h"
#include "input_image.h"
#include <bits/stdc++.h>
#include <string>
#include "complex.cc"
#include "input_image.cc"
#include <chrono>

//const float PI = 3.14159265358979f;
//threadsperblock 1024
__global__ void compute_gpu_dft(Complex *arr_old, Complex *arr_new, int *d_num_blocks, float *d_PI_val)
{
    int num_blocks = *d_num_blocks;
    float PI_val = *d_PI_val;
	int row_id = blockIdx.x/num_blocks;
	int row_elements = blockDim.x*num_blocks;
	int x_value = (blockIdx.x%num_blocks)*blockDim.x+threadIdx.x;
	int idx = x_value*row_elements + row_id;

	float rl = 0;
	float im = 0;

	for( int j = 0; j < row_elements; j++)
	{
	    float sinval, cosval;
	    //coef[i*numx+j] = Complex(cos(-2.*PI*j*x_value/row_elements), sin(-2.*PI*j*x_value/row_elements));
	    __sincosf ( -2.*PI_val*j*x_value/row_elements, &sinval, &cosval );
	    float old_real = arr_old[j+row_id*row_elements].real;
	    float old_imag = arr_old[j+row_id*row_elements].imag;

//		rl += arr_old[j+row_id*row_elements].real*coef[x_value+j*row_elements].real - arr_old[j+row_id*row_elements].imag*coef[ x_value+j*row_elements].imag;
//		im += arr_old[j+row_id*row_elements].real*coef[ x_value+j*row_elements].imag + arr_old[j+row_id*row_elements].imag*coef[ x_value+j*row_elements].real;

		rl += old_real*cosval - old_imag*sinval;
		im += old_real*sinval + old_imag*cosval;

	}
	arr_new[idx].real = rl;
	arr_new[idx].imag = im;
}


int main( int argc, char ** argv)
{

int numberofblocks=2;
//auto start = std::chrono::system_clock::now();
//FILE *fp, *fp1;
//fp1 = std::fopen(argv[3], "w");
//fp = std::fopen("heatOutput.txt", "w");
Complex *arr,*arr_old;
//Complex *W;
InputImage input_file(argv[2]);
int numx = input_file.get_width();
int T_P_B = input_file.get_height()/numberofblocks;
Complex *freq;
int sizearr = numx*numx*sizeof(Complex);
freq=(Complex*)malloc(sizearr);
Complex *final_freq=(Complex*)malloc(sizearr);
Complex *coef = (Complex*)malloc(sizearr);
freq=input_file.get_image_data();

//auto end = std::chrono::system_clock::now();
  //  auto elapsed1 = end - start;

int numberofrows = input_file.get_height();





//cpu allocation



//char * cstr = new char [argv[1].length()+1];
//std::strcpy (cstr, argv[1].c_str());

//std::string str ="Tower1024.txt";
//char c[20] = "Tower1024.txt";
//InputImage input_file("Tower1024.txt");


// freq[0].real = 0;
// freq[1].real = 0;
// freq[2].real = 1;
// freq[3].real = 0;

// freq[0].imag = -1;
// freq[1].imag = 2;
// freq[2].imag = 0;
// freq[3].imag = 0;


//for(int i = 0; i < numx; i++)
//{
//	for( int j = 0; j < numx; j++)
//	{
//		coef[i*numx+j] = Complex(cos(-2.*PI*i*j/numx), sin(-2.*PI*i*j/numx));
//	}
//}

// std::cout << "printing the coeficient array" << std::endl;
// for(int i = 0; i < numx; i++)
// {
// 	for( int j = 0; j < numx; j++)
// 	{
// 		std::cout << coef[i*numx+j].real << " " << coef[i*numx+j].imag << std::endl;
// 	}
// }

// //gpu allocation
int *d_numblocks;float *d_pi;
int sizeint = sizeof(int);
int sizefloat = sizeof(float);
cudaMalloc((void**)&arr,sizearr);
cudaMalloc((void**)&arr_old,sizearr);
cudaMalloc((void**)&d_numblocks,sizeint);
cudaMalloc((void**)&d_pi,sizefloat);

//cudaMalloc((void**)&W,sizearr);



// for(int i = 0; i < numx; i++)
// {
// 	for( int j = 0; j < numx; j++)
// 	{
// 		std::cout << freq[i*numx+j].real << " " << freq[i*numx+j].imag << std::endl;
// 	}
// }

//host to device
//_host__ ​ __device__ ​cudaError_t cudaMemcpyAsync ( arr_old, freq, sizearr, cudaMemcpyHostToDevice, cudaStream_t stream = 0 )
cudaMemcpy(arr_old,freq,sizearr,cudaMemcpyHostToDevice);
cudaMemcpy(d_pi,&PI,sizefloat,cudaMemcpyHostToDevice);
cudaMemcpy(d_numblocks, &numberofblocks,sizeint,cudaMemcpyHostToDevice);
//cudaMemcpy(arr,final_freq,sizearr,cudaMemcpyHostToDevice);
//cudaMemcpy(W,coef,sizearr,cudaMemcpyHostToDevice);


  //std::cout << elapsed.count() << '\n';
  //auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);


	//printf("Printing I %d", i);
	compute_gpu_dft<<<numx*numberofrows / T_P_B , T_P_B>>>(arr_old,arr, d_numblocks,d_pi);



// cudaMemcpy(freq,arr,sizearr,cudaMemcpyDeviceToHost);
// for(int i = 0; i < numx; i++)
// {
// 	for( int j = 0; j < numx; j++)
// 	{
// 		std::cout << freq[i*numx+j].real << " " << freq[i*numx+j].imag << std::endl;
// 	}
// }


	compute_gpu_dft<<<numx*numberofrows / T_P_B, T_P_B>>>(arr,arr_old, d_numblocks,d_pi);





//cudaMemcpy(freq,arr,sizearr,cudaMemcpyDeviceToHost);
cudaMemcpy(final_freq,arr_old,sizearr,cudaMemcpyDeviceToHost);

	// for(int i = 0; i < numx; i++)
	// {
	// 	for(int j = 0; j < numx;j++)
	// 	{
	// 		std::cout << final_freq[i*numx+j].real << std::endl;
	// 	}
	// }

	//end = std::chrono::system_clock::now();
	//auto elapsed2 = end - start;
  //std::cout << elapsed.count() << '\n';
  //auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  //std::cout << "computation taken in nano seconds: " << (elapsed2-elapsed1).count() << '\n';

    //	for(int i = 0; i < numx; i++)
    //	{
    //		for(int j = 0; j < numx;j++)
    //		{
    //			//std::cout << final_freq[i*numx+j].real/(numx*numx) << std::endl;
    //			fprintf(fp1,"(%f,",final_freq[i*numx+j].real );
    //			fprintf(fp1,"%f) ",final_freq[i*numx+j].imag );
    //		}
    //		fprintf(fp1, "\n");
    //	}

    input_file.save_image_data(argv[3], final_freq, numx, numx);
	// int i,j,rhs,lhs;
	// for(i=0;i<numx;i++)
	// {
	// 	for(j=0;j<numx;j++)
	// 	{
	// 		rhs=i+j*numx;
	// 		lhs=j+i*numx;
	// 		freq[lhs]=final_freq[rhs];
	// 	}
	// }

	// for(int i = 0; i < numx; i++)
	// {
	// 	for( int j = 0; j < numx; j++)
	// 	{
	// 		coef[i*numx+j] = Complex(cos(2.*PI*i*j/numx), sin(2.*PI*i*j/numx));
	// 	}
	// }

	// // std::cout << "printing the coeficient array" << std::endl;
	// // for(int i = 0; i < numx; i++)
	// // {
	// // 	for( int j = 0; j < numx; j++)
	// // 	{
	// // 		std::cout << coef[i*numx+j].real << " " << coef[i*numx+j].imag << std::endl;
	// // 	}
	// // }

	// cudaMemcpy(arr_old,freq,sizearr,cudaMemcpyHostToDevice);
	// cudaMemcpy(arr,final_freq,sizearr,cudaMemcpyHostToDevice);
	// cudaMemcpy(W,coef,sizearr,cudaMemcpyHostToDevice);
	// for( int i = 0; i < numx; i++)
	// {

	// 	//printf("Printing I %d", i);
	// 	compute_gpu_dft<<<numx / T_P_B, T_P_B>>>(arr_old,arr, W , i);
	// }

	// cudaDeviceSynchronize();

	// cudaMemcpy(freq,arr,sizearr,cudaMemcpyDeviceToHost);
	// // for(int i = 0; i < numx; i++)
	// // {
	// // 	for( int j = 0; j < numx; j++)
	// // 	{
	// // 		std::cout << freq[i*numx+j].real << " " << freq[i*numx+j].imag << std::endl;
	// // 	}
	// // }

	// for( int i = 0; i < numx; i++)
	// {
	// 	compute_gpu_dft<<<numx / T_P_B, T_P_B>>>(arr,arr_old, W,i);
	// }

	// cudaMemcpy(freq,arr,sizearr,cudaMemcpyDeviceToHost);
	// cudaMemcpy(final_freq,arr_old,sizearr,cudaMemcpyDeviceToHost);

	// for(i=0;i<numx;i++)
	// {
	// 	for(j=0;j<numx;j++)
	// 	{
	// 		rhs=i+j*numx;
	// 		lhs=j+i*numx;
	// 		freq[lhs]=final_freq[rhs];
	// 	}
	// }

	// for(int i = 0; i < numx; i++)
	// {
	// 	for(int j = 0; j < numx;j++)
	// 	{
	// 		//std::cout << final_freq[i*numx+j].real/(numx*numx) << std::endl;
	// 		fprintf(fp,"%f ",freq[i*numx+j].real/(numx*numx) );
	// 		//fprintf(fp,"%f ",final_freq[i*numx+j].imag/(numx*numx) );
	// 	}
	// 	fprintf(fp, "\n");
	// }

  cudaFree(final_freq);
  //cudaFree(W);
  cudaFree(freq);
  cudaFree(d_numblocks);
  cudaFree(d_pi);

  //end = std::chrono::system_clock::now();
  //auto elapsed = end - start;
  //std::cout << elapsed.count() << '\n';
  //auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  //std::cout << "Time taken in nano seconds: " << elapsed.count() << '\n';
  //fclose(fp);
  //fclose(fp1);
  return 0;
}

#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "complex.cc"
#include <sstream>
#include <sys/time.h>
#include <cstddef>

//dft11 with rows also combined into single send

void separate (Complex* a, int n) {
    Complex* b = new Complex[n/2];  // get temp heap storage
    for(int i=0; i<n/2; i++)    // copy all odd elements to heap storage
        b[i] = a[i*2+1];
    for(int i=0; i<n/2; i++)    // copy all even elements to lower-half of a[]
        a[i] = a[i*2];
    for(int i=0; i<n/2; i++)    // copy all odd (from heap) to upper-half of a[]
        a[i+n/2] = b[i];
    delete[] b;                 // delete heap storage
}

void fft2 (Complex* X, int N) {
    if(N < 2) {
        // bottom of recursion.
        // Do nothing here, because already X[0] = x[0]
    } else {
        separate(X,N);      // all evens to lower half, all odds to upper half
        fft2(X,     N/2);   // recurse even items
        fft2(X+N/2, N/2);   // recurse odd  items
        // combine results of two half recursions
        for(int k=0; k<N/2; k++) {
            Complex e = X[k    ];   // even
            Complex o = X[k+N/2];   // odd
                         // w is the "twiddle-factor"
            //Complex w = exp( Complex(0,-2.*PI*k/N) );
            Complex w = Complex (cos(-2.*PI*k/N), sin(-2.*PI*k/N));
            X[k    ] = e + w * o;
            X[k+N/2] = e - w * o;
            //cout << std::setprecision(3)<< N << " k: " << k << " w " << w << " e: " << e << " o " << o << " X[k]  "  << X[k] << "  X[k+n/2]   " << X[k+N/2] << endl;
        }
    }
}

void fft2D(Complex **X, Complex **result, int W, int H){
    for(int i = 0; i < H; i++){
        fft2(X[i], W);
    }
    cerr << "rows done\n";
    Complex *column = new Complex [H];
    cout << "after half fft\n";
//    for(int i = 0; i < height; i++){
//        for (int j = 0; j < W; j++){
//            cout << X[i][j] << "   ";
//        }
//        cout << endl;
//    }
    for(int j = 0; j < W; j++){
        cerr << "column: " << j << endl;
        for(int i = 0; i < H; i++){
            column[i] = X[i][j];
        }
        fft2(column, H);
        for(int i = 0; i < H; i++){
            result[i][j] = column[i];
        }
    }
    delete column;
}

void ifft2 (Complex* X, int N) {
    if(N < 2) {
        // bottom of recursion.
        // Do nothing here, because already X[0] = x[0]
    } else {
        separate(X,N);      // all evens to lower half, all odds to upper half
        ifft2(X,     N/2);   // recurse even items
        ifft2(X+N/2, N/2);   // recurse odd  items
        // combine results of two half recursions
        for(int k=0; k<N/2; k++) {
            Complex e = X[k    ];   // even
            Complex o = X[k+N/2];   // odd
                         // w is the "twiddle-factor"
            //Complex w = exp( Complex(0,-2.*PI*k/N) );
            Complex w = Complex (cos(2.*PI*k/N) , sin(2.*PI*k/N));
            X[k    ] = e + w * o;
            X[k+N/2] = e - w * o;
            //cout << std::setprecision(3)<< N << " k: " << k << " w " << w << " e: " << e << " o " << o << " X[k]  "  << X[k] << "  X[k+n/2]   " << X[k+N/2] << endl;
        }
    }
}

void ifft2D(Complex *X, int W, int H){
    Complex *column = new Complex [H];
    for(int j = 0; j < W; j++){
       // cerr << "column: " << j << endl;
        for(int i = 0; i < H; i++){
            column[i] = X[i*W+j];
        }
        ifft2(column, H);
        for(int i = 0; i < H; i++){
            X[i*W + j] = Complex(column[i].real/(W*H), column[i].imag/(W*H));
        }
    }
    delete column;

//    cout << "after half ift\n";
//    for(int i = 0; i < height; i++){
//        for (int j = 0; j < width; j++){
//            cout << X[i][j] << "   ";
//        }
//        cout << endl;
//    }

    for(int i = 0; i < H; i++){
        ifft2(&X[i*W], W);
    }


//    cerr << "rows done\n";
}

int main(int argc, char** argv){
    int numtasks, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype comp_datatype;
    int comp_count = 1;
    MPI_Datatype types[1] = {MPI_INT};
    int          block_length[1] = {2};
    MPI_Aint     displacement[1] = {offsetof(Complex, real)};

    MPI_Type_create_struct(comp_count, block_length, displacement, types, &comp_datatype);
    MPI_Type_commit(&comp_datatype);

//    struct timeval  curtime;
//    gettimeofday(&curtime, NULL);
//    double time_in_mill = (curtime.tv_sec) * 1000 + (curtime.tv_usec) / 1000 ;

    Complex *X, *Y;

    ifstream myFile(argv[2]);

    string line;
    getline(myFile, line);

    stringstream ss( line );


   // while( ss.good() )
    //{
    int height, width;
    ss >> height >> width;

    //cout << "rank : " << my_rank << "height: " << height << "width : " << width<< endl;
    if(numtasks > width){
        numtasks = width;
    }
    if(my_rank < width){
        int row_start[numtasks], row_stop[numtasks], column_start[numtasks], column_stop[numtasks];
        int rows_per_node = height/numtasks;
        int columns_per_node = width/numtasks;


        for(int i = 0; i < (numtasks-1); i++){
            row_start[i] = rows_per_node*i;
            row_stop[i] = rows_per_node*(i+1);
        }
        row_start[(numtasks-1)] = rows_per_node*(numtasks-1);
        row_stop[(numtasks-1)] = height;



        for(int i = 0; i < (numtasks-1); i++){
            column_start[i] = columns_per_node*i;
            column_stop[i] = columns_per_node*(i+1);
        }
        column_start[(numtasks-1)] = columns_per_node*(numtasks-1);
        column_stop[(numtasks-1)] = width;

        X = new Complex [rows_per_node*width];
        Y = new Complex [columns_per_node*height];

        int row = 0;
        while(getline(myFile, line)){
            if(row<row_start[my_rank]) {
                row += 1;
                continue;
            }
            if(row>=row_stop[my_rank]) break;
            stringstream ss1( line );
            float rl;
            for(int i = 0 ; i < width; i++){
                ss1 >> rl;
                X[(row-row_start[my_rank])*width+i].real = rl;
            }
            row += 1;
        }

//if(my_rank == 0){
//        ofstream ofs("Tower20481.txt");
//        ofs << width << " " << height << std::endl;
//
//    for(int r = 0; r < height; ++r) {
//        for(int c = 0; c < width; ++c) {
//            ofs << X[r * width + c] << " ";
//        }
//        ofs << std::endl;
//    }
//    }
//        gettimeofday(&curtime, NULL);
//        double bc_in_mill = (curtime.tv_sec) * 1000 + (curtime.tv_usec) / 1000 ;
//
//        double duration_1 = (bc_in_mill - time_in_mill)/1000;
//        if(my_rank == 0) std::cout << "before comp: " << duration_1 << std::endl;
        for(int i = 0; i < rows_per_node; i++){//compute fft of rows corresponding to the node
            //fft2realimag(X[i], width, rl[i], im[i]);
            fft2(&X[i*width], width);
        }
    //    MPI_Status send_status[numtasks], rec_status[numtasks];
    //    MPI_Request send_request[numtasks], rec_request[numtasks];

//        gettimeofday(&curtime, NULL);
//        double ac_in_mill = (curtime.tv_sec) * 1000 + (curtime.tv_usec) / 1000 ;
//
//         duration_1 = (ac_in_mill - time_in_mill)/1000;
//        if(my_rank == 0) std::cout << "after row comp: " << duration_1 << std::endl;

        int comp_columns = column_stop[my_rank] - column_start[my_rank];
        Complex *columns = new Complex [comp_columns*height];



        MPI_Status *send_status = new MPI_Status [numtasks];
        MPI_Request *send_request = new MPI_Request [numtasks];

        MPI_Status *rec_status = new MPI_Status [numtasks];
        MPI_Request *rec_request = new MPI_Request [numtasks];

        MPI_Datatype row_comp_datatype;
        int vec_count = row_stop[my_rank] - row_start[my_rank];
        int vec_blocklength = column_stop[my_rank] - column_start[my_rank];
        int vec_stride = width;
        MPI_Type_vector(vec_count, vec_blocklength, vec_stride, comp_datatype, &row_comp_datatype);
        MPI_Type_commit(&row_comp_datatype);

        for(int i = 0; i < numtasks; i++){//create array of things to be sent and send it
            if(i == my_rank) continue;
            MPI_Isend(&X[column_start[i]], 1, row_comp_datatype, i, 0, MPI_COMM_WORLD, &send_request[i]);
            MPI_Irecv(&Y[row_start[i]*columns_per_node], columns_per_node*rows_per_node, comp_datatype, i, 0, MPI_COMM_WORLD, &rec_request[i]);
        }

        for(int m = 0; m < numtasks; m ++){
            if(m == my_rank) continue;
            MPI_Wait(&send_request[m], &send_status[m]);
            MPI_Wait(&rec_request[m], &rec_status[m]);
        }

        delete send_request; delete send_status;
        delete rec_request; delete rec_status;

        for(int i = row_start[my_rank]; i < row_stop[my_rank]; i++){
            for(int j = 0; j < columns_per_node; j++){
                Y[i*columns_per_node + j] = X[(i-row_start[my_rank])*width + (j+column_start[my_rank])];
            }
        }

//        gettimeofday(&curtime, NULL);
//        double after_row = (curtime.tv_sec) * 1000 + (curtime.tv_usec) / 1000 ;
//
//        duration_1 = (after_row - time_in_mill)/1000;
//        if(my_rank == 0) std::cout << "after row: " << duration_1 << std::endl;

        MPI_Request lsend_req; MPI_Status lsend_status;

        for(int i = 0; i < comp_columns; i++){
            for(int j = 0; j < height; j++){
                //columns[i*height + j] = X[j*width + i+column_start[my_rank]];
                columns[i*height + j] = Y[j*comp_columns + i];
            }
            fft2(&columns[i*height], height);
        }
        if(my_rank != 0){
                MPI_Isend(columns, height*comp_columns, comp_datatype, 0, 0, MPI_COMM_WORLD, &lsend_req);
                MPI_Wait(&lsend_req, &lsend_status);
        }

        if(my_rank == 0){

            int col_work0 = column_stop[0] - column_start[0];
            Complex *rec_columns = new Complex[(width-col_work0)*height];
            MPI_Request *lrec_req = new MPI_Request [numtasks];
            MPI_Status *lrec_status = new MPI_Status [numtasks];

            //rec_columns[0][0] = 65;
            for(int i = 1; i < numtasks; i++){
                //for(int k = column_start[i]; k < column_stop[i]; k++){
                    //cout << "receiving from node " << i << " k -col_work0= " << k-col_work0 << endl;
                    MPI_Irecv(&rec_columns[(column_start[i]-col_work0)*height], height*comp_columns, comp_datatype, i, 0, MPI_COMM_WORLD, &lrec_req[i]);
                //}
            }
            for(int i = 1; i < numtasks; i++){
                MPI_Wait(&lrec_req[i], &lrec_status[i]);
                //cout << "wait rec rank " << my_rank << "  " << i-col_work0 << endl;
            }
            delete lrec_req; delete lrec_status;

//            gettimeofday(&curtime, NULL);
//            double after_col = (curtime.tv_sec) * 1000 + (curtime.tv_usec) / 1000 ;
//
//            duration_1 = (after_col- time_in_mill)/1000;
//            if(my_rank == 0) std::cout << "after column: " << duration_1 << std::endl;

//            Complex *Z = new Complex [width*height];

            ofstream outFile(argv[3]);
            outFile << width << " " << height << endl;
            for(int i = 0; i < height; i++){
                for(int j = 0; j < comp_columns; j++){
                    //cout << columns[j*height + i] << "  ";
                    outFile << columns[j*height + i] << " ";
//                    Z[i*width + j] = columns[j*height + i];
                }

                for(int j = column_start[1]; j < width; j++){
                    //cout << rec_columns[(j-column_start[1])*height + i] << "  ";
                    outFile << rec_columns[(j-column_start[1])*height + i] << " ";
//                    Z[i*width + j] = columns[j*height + i];
                }
        //cout << endl;
        outFile << endl;
            }
//            ifft2D(Z, width, height);
//            ofstream orr("idft_mpi.txt");
//            for(int i = 0; i < height; i++){
//                for(int j = 0; j < width; j++){
//                    orr << Z[i*width + j] << " ";
//                }
//                orr << endl;
//            }

        }
        delete columns;

    //    if(my_rank == 0){
    //        for(int i = 0; i < comp_columns; i++){
    //            for(int j = 0; j < height; j++){
    //                cout << columns[i][j] << "   ";
    //            }
    //            cout << endl;
    //        }
    //    }
    //
    //    if(my_rank == 0){
    //        for(int i = 0; i < height; i++){
    //            for(int j = 0; j < width; j++){
    //                cout << X[i][j] << "  ";
    //            }
    //	   cout << endl;
    //        }
    //    }

        delete X;
        delete Y;
        //if(my_rank == 0) std::cout << "runtime: " << duration_1 << std::endl;
    }
//cout << "ending rank " << my_rank << endl;

    MPI_Finalize();
}

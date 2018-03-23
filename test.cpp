#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

#include "C2Decomp.hpp"

int main(int argc, char *argv[]){
   int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if(!mpiRank){
        cout << endl;
        cout << "-------------------" << endl;
    	cout << " C2Decomp Testing " << endl;
    	cout << "-------------------" << endl;
    	cout << endl;

	cout << "attempting C2Decomp" << endl;
    }
 
    int nx = 50, ny = 50, nz = 50;
    int pRow = 2, pCol = 2;
    bool periodicBC[3] = {false, false, false};

    C2Decomp *c2d;
    c2d = new C2Decomp(nx, ny, nz, pRow, pCol, periodicBC);

    if(!mpiRank){
	cout << "Probably we blew up didn't we?" << endl;
    }


    int i = 3;
	if(mpiRank == i){
		cout << "Rank " << i << endl;		
	        cout << c2d->xStart[0] << " " << c2d->xStart[1] << " " << c2d->xStart[2] << endl;
	        cout << c2d->yStart[0] << " " << c2d->yStart[1] << " " << c2d->yStart[2] << endl;
	        cout << c2d->zStart[0] << " " << c2d->zStart[1] << " " << c2d->zStart[2] << endl;
	        cout << c2d->xEnd[0] << " " << c2d->xEnd[1] << " " << c2d->xEnd[2] << endl;
	        cout << c2d->yEnd[0] << " " << c2d->yEnd[1] << " " << c2d->yEnd[2] << endl;
	        cout << c2d->zEnd[0] << " " << c2d->zEnd[1] << " " << c2d->zEnd[2] << endl;
	        cout << c2d->xSize[0] << " " << c2d->xSize[1] << " " << c2d->xSize[2] << endl;
	        cout << c2d->ySize[0] << " " << c2d->ySize[1] << " " << c2d->ySize[2] << endl;
	        cout << c2d->zSize[0] << " " << c2d->zSize[1] << " " << c2d->zSize[2] << endl;

	}



    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}










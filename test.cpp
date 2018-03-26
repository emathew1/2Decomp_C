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

    }
 
    int nx = 31, ny = 17, nz = 11;
    int pRow = 2, pCol = 3;
    bool periodicBC[3] = {false, false, false};

if(!mpiRank) cout << "initializing " << endl;
    C2Decomp *c2d;
    c2d = new C2Decomp(nx, ny, nz, pRow, pCol, periodicBC);
    MPI_Barrier(MPI_COMM_WORLD);
if(!mpiRank) cout << "done initializing " << endl;
    
    int m = 0;
    double data1[nz][ny][nx];
    for(int kp = 0; kp < nz; kp++){
	for(int jp = 0; jp < ny; jp++){
	    for(int ip = 0; ip < nx; ip++){
		data1[kp][jp][ip] = (double)m;
		m++;
	    }
	}
    }

    double xSize[3], ySize[3], zSize[3];
    xSize[0] = c2d->xSize[0];
    xSize[1] = c2d->xSize[1];
    xSize[2] = c2d->xSize[2];
    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];
    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];

    double *u1, *u2;


    c2d->allocX(u1);
    c2d->allocY(u2);

    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 0; ip < xSize[0]; ip++){
		int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		u1[ii] = data1[c2d->xStart[2]+kp][c2d->xStart[1]+jp][c2d->xStart[0]+ip];
	    }
	}
    }

    for(int ip = 0; ip < totRank; ip++){
	if(mpiRank == ip){
	   cout << "---"<< mpiRank << "---" << endl;
	   for(int jp = 0; jp < xSize[0]*xSize[1]*xSize[2]; jp++){
	       cout << u1[jp] << endl;
	   } 
	}
	MPI_Barrier(MPI_COMM_WORLD);
    } 

    c2d->transposeX2Y(u1, u2);

  if(!mpiRank){
    //Testing transposition
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 0; jp < ySize[1]; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
		int ii = kp*ySize[1]*ySize[0] + jp*ySize[0] + ip;
		cout << u2[ii] << " " << data1[c2d->yStart[2]+kp][c2d->yStart[1]+jp][c2d->yStart[0]+ip] << endl;;
	    }
	}
    }
  }

    int i = 0;
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
	if(mpiRank == i){
		cout << endl;
	}

	if(mpiRank == i){
	    for(int j = 0; j < c2d->dims[0]; j++){
		cout << c2d->decompMain.x1cnts[j] << " " << c2d->decompMain.x1disp[j] << endl;
	    }
	    cout << endl;
	    for(int j = 0; j < c2d->dims[0]; j++){
		cout << c2d->decompMain.y1cnts[j] << " " << c2d->decompMain.y1disp[j] << endl;
	    }
	    cout << endl;

	    for(int j = 0; j < c2d->dims[1]; j++){
		cout << c2d->decompMain.y2dist[j] << " " << c2d->decompMain.z2dist[j] << endl;
	    }
	}


    c2d->deallocXYZ(u1);
    c2d->deallocXYZ(u2);

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}










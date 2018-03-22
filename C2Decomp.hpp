#ifndef _2DECOMPCH_
#define _2DECOMPCH_

#include "mpi.h"
#include <cstdlib>
#include <cstddef>
#include <iostream>
#include "math.h"
#include <string>

using namespace::std;

class C2Decomp{

    public:	
	//Just assume that we're using double precision all of the time
	typedef double myType; 
	MPI_Datatype realType;
	int myTypeBytes;
	
	//Global Size
	int nxGlobal, nyGlobal, nzGlobal;

	//MPI rank info
	int nRank, nProc;


    private:
	//parameters for 2D Cartesian Topology
	int dims[2], coord[2];
	int periodic[2];


    public:
	MPI_Comm DECOMP_2D_COMM_CART_X, DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z;
	MPI_Comm DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL;


    private:
	//Defining neighboring blocks 
	int neighbor[3][6];
	//Flags for periodic condition in 3D
	bool periodicX, periodicY, periodicZ;


    public:
	//Struct used to store decomposition info for a given global data size
	typedef struct decompinfo{
	    int xst[3], xen[3], xsz[3];
	    int yst[3], yen[3], ysz[3];
	    int zst[3], zen[3], zsz[3];

	    int *x1dist, *y1dist, *y2dist, *z2dist;
	    int *x1cnts, *y1cnts, *y2cnts, *z2cnts;
	    int *x1disp, *y1disp, *y2disp, *z2disp;

	    int x1count, y1count, y2count, z2count;

	    bool even;
	} DecompInfo;


    public:
	//main default decomposition information for global size nx*ny*nz
	DecompInfo decompMain;

	
    public:
	//Starting/ending index and size of data held by the current processor	
	//duplicate 'decompMain', needed by apps to define data structure
	int xStart[3], xEnd[3], xSize[3]; //x-pencil
	int yStart[3], yEnd[3], ySize[3]; //y-pencil
	int zStart[3], zEnd[3], zSize[3]; //z-pencil


    private:
	//These are the buffers used by MPI_ALLTOALL(V) calls
	int decompBufSize;
	myType *work1_r, *work2_r; //Only implementing real for now... 
        

    public:

	C2Decomp(int nx, int ny, int nz, int pRow, int pCol, bool periodicBC[3]){

	    nxGlobal = nx;
	    nyGlobal = ny;
	    nzGlobal = nz;

	    periodicX = periodicBC[0];	
	    periodicY = periodicBC[1];	
	    periodicZ = periodicBC[2];	
	
	    decompBufSize = 0;
	    work1_r = NULL;
	    work2_r = NULL;
	
	    realType = MPI_DOUBLE_PRECISION;

	    decomp2DInit(pRow, pCol);
	} 

	void decomp2DInit(int pRow, int pCol);

	void decomp2DFinalize();

	//Just get it running without the optional decomp for now...
	void transposeX2Y(myType ***src, myType ***dst);
	void transposeY2Z(myType ***src, myType ***dst);
	void transposeZ2Y(myType ***src, myType ***dst);
	void transposeY2X(myType ***src, myType ***dst);
	
	//calls for overlapping communication and computation...
	void transposeX2Y_Start(int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	void transposeX2Y_Wait (int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);

	void transposeY2Z_Start(int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	void transposeY2Z_Wait (int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);

	void transposeZ2Y_Start(int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	void transposeZ2Y_Wait (int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	
	void transposeY2X_Start(int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	void transposeY2X_Wait (int handle, myType ***src, myType ***dst, myType ***sbuf, myType ***rbuf);
	
	void decompInfoInit();
	void decompInfoFinalize(DecompInfo decompinfo_in);


	//only doing real 
	void allocX(myType ***var); 
	void allocY(myType ***var); 
	void allocZ(myType ***var); 

	void updateHalo(myType ***in, myType ***out, int level);

	void decomp2DAbort(int errorCode, string msg);
	void initNeighbor();	
	void getDist();
	void distribute(int data1, int proc, int *st, int *en, int *sz);
	void partition(int nx, int ny, int nz, int pdim[3], int lstart[3], int lend[3], int lsize[3]);
	void prepareBuffer(DecompInfo *dii);

	void getDecompInfo(DecompInfo dcompinfo_in);
	

};


#endif
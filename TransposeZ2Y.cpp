#include "C2Decomp.hpp"

void C2Decomp::transposeZ2Y(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Alltoallv(src,     decompMain.z2cnts, decompMain.z2disp, realType, 
		  work2_r, decompMain.y2cnts, decompMain.y2disp, realType, 		   
		  DECOMP_2D_COMM_ROW);

    memMergeZY(work2_r, d1, d2, d3, dst, dims[1], decompMain.y2dist);

}



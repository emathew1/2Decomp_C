#include "C2Decomp.hpp"

void C2Decomp::transposeY2Z(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    s1 = decompMain.ysz[0];
    s2 = decompMain.ysz[1];
    s3 = decompMain.ysz[2];

    d1 = decompMain.zsz[0];
    d2 = decompMain.zsz[1];
    d3 = decompMain.zsz[2];

    memSplitYZ(src, s1, s2, s3, work1_r, dims[1], decompMain.y2dist);

    MPI_Alltoallv(work1_r, decompMain.y2cnts, decompMain.y2disp, realType, 
		  dst,     decompMain.z2cnts, decompMain.z2disp, realType, 		   
		  DECOMP_2D_COMM_ROW);


}



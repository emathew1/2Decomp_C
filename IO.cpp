#include "C2Decomp.hpp"

void C2Decomp::writeOne(int ipencil, double *var, string filename){

    MPI_Offset disp, filesize;
    MPI_File fh;
    MPI_Datatype data_type, new_type;

    int sizes[3], subsizes[3], starts[3];
    
    sizes[0] = decompMain.xsz[0];
    sizes[1] = decompMain.ysz[1];
    sizes[2] = decompMain.zsz[2];

    if(ipencil == 0){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.xsz[ip];
	    starts[ip] = decompMain.xst[ip]-1;
	}
    }else if(ipencil == 1){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.ysz[ip];
	    starts[ip] = decompMain.yst[ip]-1;
	}
    }else if(ipencil == 2){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.zsz[ip];
	    starts[ip] = decompMain.zst[ip]-1;
	}
    }

    data_type = MPI_DOUBLE;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);
    MPI_Type_commit(&new_type);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    filesize = 0;
    MPI_File_set_size(fh, filesize);

    disp = 0;
    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(fh, var, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    MPI_Type_free(&new_type);

}

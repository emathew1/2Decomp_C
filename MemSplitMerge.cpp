#include "C2Decomp.hpp"

void C2Decomp::memSplitXY(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;

    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2= i1+dist[m]-1;
	}

	pos = decompMain.x1disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = 0; j < n2; j++){
		for(int i = (i1-1); i < i2; i++){ 
		    int ii = k*n2*n1 + j*n1 + i;
		    out[pos] = in[ii];
		    pos++;
		}
	    }
	}

    }

};


void C2Decomp::memMergeXY(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;

    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.y1disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = (i1-1); j < i2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[ii] = in[pos];
		    pos++;

		}
	    }
	}

    }

}


void C2Decomp::memSplitYZ(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.y2disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = (i1-1); j < i2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[pos] = in[ii];
		    pos++;
		}
	    }
	}

    }

};

void C2Decomp::memMergeYZ(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.z2disp[m];

	for(int k = (i1-1); k < i2; k++){
	    for(int j = 0; j < n2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[ii] = in[pos];
		    pos++;
		}
	    }
	}

    }

};

void C2Decomp::memSplitZY(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.z2disp[m];

	for(int k = (i1-1); k < i2; k++){
	    for(int j = 0; j < n2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[pos] = in[ii];
		    pos++;
		}
	    }
	}

    }

};

void C2Decomp::memMergeZY(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.y2disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = (i1-1); j < i2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[ii] = in[pos];
		    pos++;
		}
	    }
	}

    }

}

void C2Decomp::memSplitYX(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.y1disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = (i1-1); j < i2; j++){
		for(int i = 0; i < n1; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[pos] = in[ii];
		    pos++;
		}
	    }
	}

    }

};

void C2Decomp::memMergeYX(double *in, int n1, int n2, int n3, double *out, int iproc, int *dist){

    int i1, i2, pos;
    
    for(int m = 0; m < iproc; m++){
	if(m == 0){
	    i1 = 1;
	    i2 = dist[0];
	}else{
	    i1 = i2+1;
	    i2 = i1+dist[m]-1;
	}

	pos = decompMain.x1disp[m];

	for(int k = 0; k < n3; k++){
	    for(int j = 0; j < n2; j++){
		for(int i = (i1-1); i < i2; i++){
		    int ii = k*n2*n1 + j*n1 + i;
		    out[ii] = in[pos];
		    pos++;
		}
	    }
	}

    }

}



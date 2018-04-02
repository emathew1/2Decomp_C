#include "C2Decomp.hpp"

void C2Decomp::FindFactor(int num, int *factors, int &nfact){

    int m; 

    //Finding factors <= sqrt(num)
    m = (int)sqrt((double)num);
    nfact = 1;
    for(int ip = 1; ip < m+1; ip++){
        if(num/ip*ip == num){
	    factors[nfact-1] = ip;
	    nfact++;
        }

    }
    nfact--;

    //Finding factors > sqrt(num)
    if(factors[nfact-1]*factors[nfact-1] != num){
	for(int ip = nfact+1; ip < 2*nfact+1; ip++){
	    factors[ip-1] = num/factors[2*nfact-ip];
	}
	nfact *= 2;
    }else{
	for(int ip = nfact+1; ip < 2*nfact; ip++){
	    factors[ip-1] = num/factors[2*nfact-ip];
	}	
	nfact = nfact*2 - 1;
    }

   
}

void C2Decomp::best2DGrid(int nProc, int &row, int &col){

    if(!nRank){
	int *factors = new int[20];
 	int nfact = 0;
	cout << "HERE" << endl;

	FindFactor(96, factors, nfact);

	for(int ip = 0; ip < nfact; ip++){
	    cout << factors[ip] << endl;
	}
    }
}



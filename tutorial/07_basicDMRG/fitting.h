//fitting.h


#ifndef ___FITTING
#define ___FITTING


//#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "sweep_utility.h"
#include "parameters.h"


using namespace std;
namespace itensor {
    //------begin: fitting sweeps--------------------------------------------------------
    std::vector<ITensor> fitting_sweep(std::vector<ITensor> N, std::vector<ITensor> M, double chi){
        int sizeN=N.size(), sizeM=M.size();
        int size=0;
        if (sizeN==sizeM) {
            size=sizeN;
        }
        else {
            cout<<"Error: Sizes of M and N in the fitting process do not Match!!!"<<endl;
            exit(0);
        }
        
        /*
        double diff;
        ITensor product=(N[0]-M[0])*dag( prime(N[0]-M[0],(N[0]-M[0]).indices()[0]) );
        
        for (int x=1; x<size-1; x++) {//int x=1; x<size; x++
            product=product*((N[x]-M[x])*dag(prime(prime(N[x]-M[x], (N[x]-M[x]).indices()[0]),(N[x]-M[x]).indices()[1])));
        }
        product=product*(N[size-1]-M[size-1])*dag( prime(N[size-1]-M[size-1],(N[size-1]-M[size-1]).indices()[0]) );
        
        PrintDat(product);
        //Print( prime(prime(N[1]-M[1], (N[1]-M[1]).indices()[0]),(N[1]-M[1]).indices()[1]) );
        */
        //PrintDat(N[])
        return M;
    }
//------end: fitting sweeps---------------------------------------------------------
}
#endif





















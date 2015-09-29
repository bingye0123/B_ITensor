//this is the Heisenberg model on square lattice with open boundary condition

#ifndef ___HEISENBERG_SL
#define ___HEISENBERG_SL


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "tensor_utility.h"
#include "parameters.h"
//#include "fitting.h"


//#define NX 4
//#define NY 3
//#define NS 3
//#define NB 3

//static int N=NX*NY;

using namespace std;
namespace itensor {
    
    //mapping bond index number to the two site indices numbers that the bond has, all sites and bonds are labeled by  a indexed number
    std::vector<int> mapping(int i){
        std::vector<int> temp(2);
        
        temp[0]=(N-1)/2;
        temp[1]=i;
        
        return temp;
    }
    
    
    //--------begin: setting up the tensor network for hamilnonian-------------------------------------------------------
    void setup_correlator_tensors(int hamx, int nx, int ny){
        
        for (int i=0; i<N; i++) {
            //U[i]=Index(nameint("U",i),NB*NB);
            //D[i]=Index(nameint("D",i),NB*NB);
            L[i]=Index(nameint("L",i),NB*NB);
            R[i]=Index(nameint("R",i),NB*NB);
        }
        for (int i=0; i<N; i++) {
            if ((i+1)%nx!=0) {
                R[i]=L[i+1];
            }
            //if (i/nx!=ny-1) {
            //    D[i]=U[i+nx];
            //}
        }
        
        
        if (NS==3) {
            //---------begin: Define tensors Hi_ss'=<s|Hi|s'>------------------------------
            if (hamx<(terms/3)) { //the first quater: S^z*S^z
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamx)[0]) && i!=(mapping(hamx)[1])) {
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                if (mapping(hamx)[0]!=mapping(hamx)[1]) {
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](1), primed(s[mapping(hamx)[0]])(1))=1;
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](3), primed(s[mapping(hamx)[0]])(3))=-1;
                    O[N*hamx+mapping(hamx)[1]](s[mapping(hamx)[1]](1), primed(s[mapping(hamx)[1]])(1))=1;
                    O[N*hamx+mapping(hamx)[1]](s[mapping(hamx)[1]](3), primed(s[mapping(hamx)[1]])(3))=-1;
                }
                else {
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](1), primed(s[mapping(hamx)[0]])(1))=1;
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](2), primed(s[mapping(hamx)[0]])(2))=1;
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](3), primed(s[mapping(hamx)[0]])(3))=1;
                }
                
            }
            else if (hamx<(2*terms/3)){ //the second quater: 1/2*S^+*S^-
                int hamxtemp=hamx-terms/3;
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamxtemp)[0]) && i!=(mapping(hamxtemp)[1])) {//sites that are not in the specific bond: identity
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                
                if (mapping(hamxtemp)[0]!=mapping(hamxtemp)[1]) {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](1), primed(s[mapping(hamxtemp)[0]])(2))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(3))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](2), primed(s[mapping(hamxtemp)[1]])(1))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](3), primed(s[mapping(hamxtemp)[1]])(2))=_1oversq2;
                }
                else {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](1), primed(s[mapping(hamxtemp)[0]])(1))=0.5;
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(2))=0.5;
                }
            }
            else {
                int hamxtemp=hamx-2*terms/3;
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamxtemp)[0]) && i!=(mapping(hamxtemp)[1])) {//sites that are not in the specific bond: identity
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                
                if (mapping(hamxtemp)[0]!=mapping(hamxtemp)[1]) {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(1))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](3), primed(s[mapping(hamxtemp)[0]])(2))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](1), primed(s[mapping(hamxtemp)[1]])(2))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](2), primed(s[mapping(hamxtemp)[1]])(3))=_1oversq2;
                }
                else {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(2))=0.5;
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](3), primed(s[mapping(hamxtemp)[0]])(3))=0.5;

                }
                
            }
        }
        if (NS==2) {
            //---------begin: Define tensors Hi_ss'=<s|Hi|s'>------------------------------
            if (hamx<(terms/3)) { //the first quater: S^z*S^z
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamx)[0]) && i!=(mapping(hamx)[1])) {
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                if (mapping(hamx)[0]!=mapping(hamx)[1]) {
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](1), primed(s[mapping(hamx)[0]])(1))=1;
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](2), primed(s[mapping(hamx)[0]])(2))=-1;
                    O[N*hamx+mapping(hamx)[1]](s[mapping(hamx)[1]](1), primed(s[mapping(hamx)[1]])(1))=1;
                    O[N*hamx+mapping(hamx)[1]](s[mapping(hamx)[1]](2), primed(s[mapping(hamx)[1]])(2))=-1;
                }
                else {
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](1), primed(s[mapping(hamx)[0]])(1))=1;
                    O[N*hamx+mapping(hamx)[0]](s[mapping(hamx)[0]](2), primed(s[mapping(hamx)[0]])(2))=1;
                }
            }
            else if (hamx<(2*terms/3)){ //the second quater: 1/2*S^+*S^-
                int hamxtemp=hamx-terms/3;
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamxtemp)[0]) && i!=(mapping(hamxtemp)[1])) {//sites that are not in the specific bond: identity
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                
                if (mapping(hamxtemp)[0]!=mapping(hamxtemp)[1]) {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](1), primed(s[mapping(hamxtemp)[0]])(2))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](2), primed(s[mapping(hamxtemp)[1]])(1))=_1oversq2;
                }
                else {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](1), primed(s[mapping(hamxtemp)[0]])(1))=0.5;
                }
            }
            else {
                int hamxtemp=hamx-2*terms/3;
                for (int i=0; i<N; i++) {
                    O[N*hamx+i]=ITensor(s[i],primed(s[i]));
                    if (i!=(mapping(hamxtemp)[0]) && i!=(mapping(hamxtemp)[1])) {//sites that are not in the specific bond: identity
                        for (int k=1; k<=NS; k++) {
                            O[N*hamx+i](s[i](k),primed(s[i])(k))=1;
                        }
                    }
                }
                
                if (mapping(hamxtemp)[0]!=mapping(hamxtemp)[1]) {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(1))=_1oversq2;
                    O[N*hamx+mapping(hamxtemp)[1]](s[mapping(hamxtemp)[1]](1), primed(s[mapping(hamxtemp)[1]])(2))=_1oversq2;
                }
                else {
                    O[N*hamx+mapping(hamxtemp)[0]](s[mapping(hamxtemp)[0]](2), primed(s[mapping(hamxtemp)[0]])(2))=0.5;
                }
            }
        }
        
//        for (int i=0; i<N; i++) {
//            cout<<"hamx="<<hamx<<","<<i;
//            PrintDat(O[N*hamx+i]);
//        }
        
        
        //---------begin: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=y*nx+x;
                //                if (backflag==1) {
                //                    cout<<"ind="<<ind;
                //                    Print(a[ind]);
                //                    Print(O[ind+N*hamx]);
                //                }
                
                //                if (hamx==1) {
                //                    cout<<"x="<<x<<", y="<<y;
                //                    Print(O[ind+N*hamx]);
                //                    Print(a[ind]);
                //                }
                
                if (y==0) {
                    if (x==0) {
                        //cout<<"193~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice4((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], R[ind]);
                        //cout<<"195~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"198~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice4((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind]);
                        //cout<<"200~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"203~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice6((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"205~~~~~~~~"<<endl;
                    }
                }
                else if (y==ny-1) {
                    if (x==0) {
                        //cout<<"210~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice4((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], R[ind]);
                        //cout<<"212~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"215~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice4((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind]);
                        //cout<<"217~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"220~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice6((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"222~~~~~~~~"<<endl;
                    }
                }
                else {
                    if (x==0) {
                        //cout<<"227~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice6((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"229~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"232~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice6((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], L[ind]);
                        //cout<<"234~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"237~~~~~~~~"<<endl;
                        W[ind+nx*ny*hamx]=group_indice8((O[ind+N*hamx]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], L[ind], a[ind].indices()[4], dag(prime(a[ind])).indices()[4], R[ind]);
                        //cout<<"239~~~~~~~~"<<endl;
                    }
                }
            }
        }
        
    }
    //--------end: setting up the tensor network for hamiltonian-------------------------------------------------------
    
    
    
    
    //--------begin: setting up the tensor network for identity operator-------------------------------------------------------
    void setup_identity_tensors(int nx, int ny){
        
        for (int i=0; i<N; i++) {
            //U[i]=Index(nameint("U",i),NB*NB);
            //D[i]=Index(nameint("D",i),NB*NB);
            L[i]=Index(nameint("L",i),NB*NB);
            R[i]=Index(nameint("R",i),NB*NB);
        }
        for (int i=0; i<N; i++) {
            if ((i+1)%nx!=0) {
                R[i]=L[i+1];
            }
            //if (i/nx!=ny-1) {
            //D[i]=U[i+nx];
            //}
        }
        
        //---------begin: Define tensors N_ss'=<s|N|s'>------------------------------
        for (int i=0; i<N; i++) {
            Identity[i]=ITensor(s[i],primed(s[i]));
            for (int k=1; k<=NS; k++) {
                Identity[i](s[i](k),primed(s[i])(k))=1;
            }
        }
        //---------end: Define tensors O_ss'=<s|O|s'>--------------------------------
        /*
         for (int i=0; i<N; i++) {
         U[i]=Index(nameint("U",i),NB*NB);
         D[i]=Index(nameint("D",i),NB*NB);
         L[i]=Index(nameint("L",i),NB*NB);
         R[i]=Index(nameint("R",i),NB*NB);
         }
         for (int i=0; i<N; i++) {
         if ((i+1)%nx!=0) {
         R[i]=L[i+1];
         }
         if (i/nx!=ny-1) {
         D[i]=U[i+nx];
         }
         }
         */
        //---------begin: group indices to get IA, IB and IW: IA in MPS bra, IB in MPS ket, and IW in MPO-------------
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=y*nx+x;
                
                if (y==0) {
                    if (x==0) {
                        //cout<<"299~~~~~~~~"<<endl;
                        IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], R[ind]);
                        //cout<<"301~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"304~~~~~~~~"<<endl;
                        IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind]);
                        //cout<<"306~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"309~~~~~~~~"<<endl;
                        IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], D[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"311~~~~~~~~"<<endl;
                    }
                }
                else if (y==ny-1) {
                    if (x==0) {
                        //cout<<"316~~~~~~~~"<<endl;
                        IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], R[ind]);
                        //cout<<"318~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"321~~~~~~~~"<<endl;
                        IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind]);
                        //cout<<"323~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"326~~~~~~~~"<<endl;
                        IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], L[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"328~~~~~~~~"<<endl;
                    }
                }
                else {
                    if (x==0) {
                        //cout<<"333~~~~~~~~"<<endl;
                        
                        //                        cout<<"x="<<x<<", y="<<y;
                        //                        Print(Identity[ind]);
                        //                        Print(a[ind]);
                        //                        Print(u[ind]);
                        //                        Print(U[ind]);
                        //                        Print(d[ind]);
                        //                        Print(D[ind]);
                        //                        Print(r[ind]);
                        //                        Print(R[ind]);
                        
                        IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], R[ind]);
                        //cout<<"335~~~~~~~~"<<endl;
                    }
                    else if (x==nx-1) {
                        //cout<<"338~~~~~~~~"<<endl;
                        IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], L[ind]);
                        //cout<<"340~~~~~~~~"<<endl;
                    }
                    else {
                        //cout<<"343~~~~~~~~"<<endl;
                        IW[ind]=group_indice8((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], dag(prime(a[ind])).indices()[1], U[ind], a[ind].indices()[2], dag(prime(a[ind])).indices()[2], D[ind], a[ind].indices()[3], dag(prime(a[ind])).indices()[3], L[ind], a[ind].indices()[4], dag(prime(a[ind])).indices()[4], R[ind]);
                        //cout<<"345~~~~~~~~"<<endl;
                    }
                }
            }
        }
        
        //start: need to be comment out!!!
        
        //PrintDat(a[1]);
        //for (int x=1; x<2; x++) {
        //    cout<<"x="<<x;
        //    PrintDat(IA[x]);
        //Print(IB[x]);
        //Print(IW[x]);
        //}
        
        //end: need to be comment out!!!
        //---------end: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
        
        //PrintDat(A[2]);
    }
    //--------end: setting up the tensor network for identity operator-------------------------------------------------------
    
}
#endif

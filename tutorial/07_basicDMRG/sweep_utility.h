//sweep_utility.h
//include: setup_tensors(), compressor2


#ifndef ___SWEEP_UTILITY
#define ___SWEEP_UTILITY


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "tensor_utility.h"
#include "parameters.h"
//#include "fitting.h"


//#define nx 4
//#define ny 3
//#define NS 3
//#define NB 3

//static int N=nx*ny;

using namespace std;
namespace itensor {
    
    //static std::vector<Index> s(N),u(N),d(N),l(N),r(N);
    //static std::vector<ITensor> O(N);
    //static std::vector<ITensor> a(N); //site tensors
    //static std::vector<ITensor> A(nx), B(nx),W(nx*(ny-2));
    //static std::vector<Index> U(N),D(N),L(N),R(N);
    /*
    //--------begin: setting up the tensor network-------------------------------------------------------
    void setup_tensors(){
        for (int i=0; i<N; i++) {
            s[i]=Index(nameint("s",i),NS,Site);
            u[i]=Index(nameint("u",i),NB);
            d[i]=Index(nameint("d",i),NB);
            l[i]=Index(nameint("l",i),NB);
            r[i]=Index(nameint("r",i),NB);
        }
        for (int i=0; i<N; i++) {
            if ((i+1)%nx!=0) {
                r[i]=l[i+1];
            }
            if (i/nx!=ny-1) {
                d[i]=u[i+nx];
            }
        }
        //---------begin: Define tensors O_ss'=<s|O|s'>------------------------------
        ////trial: only one term JS1zS2z in Heisenberg model///
        
        for (int i=0; i<N; i++) {
            O[i]=ITensor(s[i],primed(s[i]));
            for (int k=1; k<=NS; k++) {
                O[i](s[i](k),primed(s[i])(k))=1;
            }
        }
        O[0](s[0](3),primed(s[0])(3))=-1;
        O[1](s[1](3),primed(s[1])(3))=-1;
        //---------end: Define tensors O_ss'=<s|O|s'>--------------------------------
        
        //---------begin: Define and randomize site tensors A^s_udlr-----------------
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=y*nx+x;
                if (x==0) {
                    if (y==0) {
                        a[ind]=ITensor(s[ind],d[ind],r[ind]);
                    }
                    else if (y==ny-1) {
                        a[ind]=ITensor(s[ind],u[ind],r[ind]);
                    }
                    else{
                        a[ind]=ITensor(s[ind],u[ind],d[ind],r[ind]);
                    }
                }
                else if (x==nx-1) {
                    if (y==0) {
                        a[ind]=ITensor(s[ind],d[ind],l[ind]);
                    }
                    else if (y==ny-1) {
                        a[ind]=ITensor(s[ind],u[ind],l[ind]);
                    }
                    else{
                        a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind]);
                    }
                }
                else if (x>0&&x<nx-1&&y==0) {
                    a[ind]=ITensor(s[ind],d[ind],l[ind],r[ind]);
                }
                else if (x>0&&x<nx-1&&y==ny-1) {
                    a[ind]=ITensor(s[ind],u[ind],l[ind],r[ind]);
                }
                else{
                    a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind],r[ind]);
                }
                a[ind].randomize();
                
                //????????????????????????????????????????????????????????????????????????????????????
                //a[ind] /= a[ind].norm();            //?????????normalized????????????????????????????
                //????????????????????????????????????????????????????????????????????????????????????
                ////E[ind]=O[ind]*a[ind]*dag(prime(a[ind]));
                //PrintDat(A[ind]);
            }
        }
        //---------end: Define and randomize site tensors A^s_udlr-----------------
        
        //---------begin: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
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
        
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=y*nx+x;
                int nind=ind-nx;
                if (y==0) {
                    if (x==0) {
                        A[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                    else if (x==nx-1) {
                        A[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind]);
                    }
                    else {
                        A[x]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                }
                else if (y==ny-1) {
                    if (x==0) {
                        B[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                    else if (x==nx-1) {
                        B[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], l[ind], primed(l[ind]), L[ind]);
                    }
                    else {
                        B[x]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], l[ind], primed(l[ind]), L[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                }
                else {
                    if (x==0) {
                        W[nind]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], d[ind], primed(d[ind]), D[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                    else if (x==nx-1) {
                        W[nind]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind]);
                    }
                    else {
                        W[nind]=group_indice8((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind], r[ind], primed(r[ind]), R[ind]);
                    }
                }
            }
        }
        //---------end: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
        //PrintDat(A[2]);
    }
    //--------end: setting up the tensor network-------------------------------------------------------
    */
    
    

    
    //---------------parameters:   1.input bra: temp1  2.MPO: gate----------------------------------------------------------//
    //----------------------3.indix: i, with y=ny-i where y is the vertical coordinate in 2D system-------------------------//
    //----------------------(the coordinate system: x axis pointing right, y axis pointing down, origin at up-left corner)--//
    //----------------------4.SVD cutoff for zip-up algorithm---------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------//
    //---------------The algorithm:-----------------------------------------------------------------------------------------//
    //---------------(1)set the othogonality center of temp1 and gate to the first site: sweep from right to left-----------//
    //---------------(2)perform zip-up algorithm: sweep from left to right--------------------------------------------------//
    //---------------(3)perform fitting algorithm to fit the output bra with K=||output_bra-gate*temp1||^2 minimum----------//
    
    //--------begin: compressor of bra, starting from the bottom of the system----------------------------------------------
    std::vector<ITensor> compressor2(std::vector<ITensor> temp1, std::vector<ITensor> gate, int i, double cutoff, int m, int dcut, int nx, int ny){
        
        int y=ny-1-i;
        
        std::vector<ITensor> temp(nx);    //this is the original mps
        std::vector<ITensor> gate_original(nx);    //this is the original mps
        for (int x=0; x<nx; x++) {
            temp[x]=temp1[x];
            gate_original[x]=gate[x];
        }
        
        gate=mpo_otho_center_1(gate, y, nx, ny);    //move the othogonality center of MPO to the first site
        //temp1=bra_otho_center_1(temp1, y);
        //-----begin: zip-up algorithm-----------------------
        
        
        int bind[nx],wind[nx];  //B index and W index

        for (int x=0; x<nx; x++) {
            bind[x]=x+nx*(ny-i);
            wind[x]=x+nx*(ny-1-i);
        }
        
        std::vector<ITensor> c(nx);//, S(nx);
        
        c[0]=gate[0]*temp1[0];  ////zip-up first step: starting from the first site in this row
        
        
        
        OptSet opts;
        opts.add("Cutoff", cutoff);
        ITensor V0,D0;
        
//        PrintDat(temp1[0]);
//        PrintDat(gate[0])
//        PrintDat(c[0]);
        
        
        temp1[0]=ITensor(gate_original[0].indices()[0]);
        svd(c[0],temp1[0],V0,D0,opts);        ////svd: with truncation error cutoff
        R[wind[0]]=Index(nameint("R",wind[0]),temp1[0].indices()[1].m());    ///rename indice from svd to beta
        L[wind[1]]=Index(nameint("L",wind[1]),temp1[0].indices()[1].m());    ///rename indice from svd to beta
        R[wind[0]]=L[wind[1]];
        temp1[0]=change_index_name2(temp1[0], gate_original[0].indices()[0], R[wind[0]]);          ///temp[0] is the U0 in the written notes
        
        
        //--------
        for (int x=1; x<nx-1; x++) {          ////int x=1; x<nx-1; x++
            c[x]=dag(prime(temp1[x-1],gate_original[x-1].indices()[0]))*prime(c[x-1],gate_original[x-1].indices()[0])*gate[x]*temp1[x];

            temp1[x]=ITensor(L[wind[x]],gate_original[x].indices()[0]);
            svd(c[x], temp1[x], V0, D0,opts);

            R[wind[x]]=Index(nameint("R",wind[x]),temp1[x].indices()[2].m());
            L[wind[x+1]]=Index(nameint("L",wind[x+1]),temp1[x].indices()[2].m());
            R[wind[x]]=L[wind[x+1]];
            temp1[x]=change_index_name3(temp1[x], L[wind[x]], gate_original[x].indices()[0], R[wind[x]]);
            
        }
        
        c[nx-1]=dag(prime(temp1[nx-2],gate_original[nx-2].indices()[0]))*prime(c[nx-2],gate_original[nx-2].indices()[0])*gate[nx-1]*temp1[nx-1];
        temp1[nx-1]=c[nx-1];
        
        //for (int x=0; x<nx; x++) {
        //    Print(temp1[x]);
        //    if (x==1) {
        //        PrintDat(temp1[x]);
        //    }
        //}
        //cout<<"~~~~~~~"<<endl;
        
        /////temp1=bra_otho_center_1m(temp1, y-1, m);//move the othogonality center of input bra to the first site with size cutoff m
        //for (int x=0; x<nx; x++) {
        //    Print(temp1[x]);
        //}
        
        temp1[nx-1]=sort_indices_braket(temp1[nx-1]);
        for (int x=1; x<nx-1; x++) {
            temp1[x]=sort_indices_compressor(temp1[x]);
        }
        
        temp1=bra_otho_center_1(temp1, y, nx, ny);//move the othogonality center of input bra to the first site without cutoff
        //-----end: zip-up algorithm-------------------------
        
        //for (int x=0; x<nx; x++) {
        //    Print(temp1[x]);
            //if (x==1) {
            //    PrintDat(temp1[x]);
            //}
        //}
        //cout<<">>>>>>>>>>>>>>>>"<<endl;
        
//        for (int x=0; x<nx; x++) {
//            PrintDat(temp1[x]);
//        }
//        for (int x=0; x<nx; x++) {
//            PrintDat(gate[x]);
//        }
        
        if (nx>=3) {
            //-----begin: fitting algorithm----------------------
            temp1=fitting_algorithm(temp1, gate, temp, dcut, i, nx, ny);
            //-----end: fitting algorithm------------------------
        }

        //temp1=bra_otho_center_1(temp1,y,nx,ny);//~~~~~~~~~~~~~~
        
         
         
        return temp1;   //return to the output bra=MPO*bra
    }
    //--------end: compressor of bra, starting from the bottom of the system----------------------------------------------
    
    //--------begin: compressor of ket, starting from the top of the system----------------------------------------------
    std::vector<ITensor> compressor1(std::vector<ITensor> temp1, std::vector<ITensor> gate, int y, double cutoff, int m, int dcut, int nx, int ny){
        std::vector<ITensor> temp(nx);    //this is the original mps
        std::vector<ITensor> gate_original(nx);    //this is the original mps
        for (int x=0; x<nx; x++) {
            temp[x]=temp1[x];
            gate_original[x]=gate[x];
        }
        //cout<<"287~~~~"<<endl;
        gate=mpo_otho_center_1(gate, y, nx, ny);    //move the othogonality center of MPO to the first site
        //cout<<"289~~~~"<<endl;
        //temp1=bra_otho_center_1(temp1, y);
        //-----begin: zip-up algorithm-----------------------
        int bind[nx],wind[nx];  //B index and W index
        
        for (int x=0; x<nx; x++) {
            bind[x]=x+nx*(y-1);
            wind[x]=x+nx*y;
        }
        
        std::vector<ITensor> c(nx);//, S(nx);
        //cout<<"300~~~~"<<endl;
        c[0]=gate[0]*temp1[0];  ////zip-up first step: starting from the first site in this row
        //Print(c[0]);
        OptSet opts;
        opts.add("Cutoff", cutoff);
        ITensor V0,D0;
        //cout<<"305~~~~"<<endl;
        temp1[0]=ITensor(gate_original[0].indices()[1]);
        svd(c[0],temp1[0],V0,D0,opts);        ////svd: with truncation error cutoff
        
//        cout<<"~~~~~~~";
//        Print(c[0]);
//        Print(temp1[0]);
//        Print(V0);
//        Print(D0);
        
        
        R[wind[0]]=Index(nameint("R",wind[0]),temp1[0].indices()[1].m());    ///rename indice from svd to beta
        L[wind[1]]=Index(nameint("L",wind[1]),temp1[0].indices()[1].m());    ///rename indice from svd to beta
        R[wind[0]]=L[wind[1]];
        temp1[0]=change_index_name2(temp1[0], gate_original[0].indices()[1], R[wind[0]]);          ///temp[0] is the U0 in the written notes
        
        //Print(temp1[0]);
        
        //cout<<"313~~~~"<<endl;
        for (int x=1; x<nx-1; x++) {          ////int x=1; x<nx-1; x++
            c[x]=dag(temp1[x-1])*c[x-1]*gate[x]*temp1[x];
            
            temp1[x]=ITensor(L[wind[x]],gate_original[x].indices()[1]);
            svd(c[x], temp1[x], V0, D0,opts);
            
            R[wind[x]]=Index(nameint("R",wind[x]),temp1[x].indices()[2].m());
            L[wind[x+1]]=Index(nameint("L",wind[x+1]),temp1[x].indices()[2].m());
            R[wind[x]]=L[wind[x+1]];
            temp1[x]=change_index_name3(temp1[x], L[wind[x]], gate_original[x].indices()[1], R[wind[x]]);
            
        }
        
        //cout<<"327~~~~"<<endl;
        c[nx-1]=dag(temp1[nx-2])*c[nx-2]*gate[nx-1]*temp1[nx-1];
        temp1[nx-1]=c[nx-1];
        
        //for (int x=0; x<nx; x++) {
        //    Print(temp1[x]);
        //    if (x==1) {
        //        PrintDat(temp1[x]);
        //    }
        //}
        //cout<<"~~~~~~~"<<endl;
        
        /////temp1=bra_otho_center_1m(temp1, y-1, m);//move the othogonality center of input bra to the first site with size cutoff m
        //for (int x=0; x<nx; x++) {
        //    Print(temp1[x]);
        //}
        //cout<<"343~~~~"<<endl;
        temp1[nx-1]=sort_indices_braket(temp1[nx-1]);
        for (int x=1; x<nx-1; x++) {
            temp1[x]=sort_indices_compressor(temp1[x]);
        }
        //cout<<"348~~~~"<<endl;
        temp1=ket_otho_center_1(temp1, y, nx, ny);//move the othogonality center of input bra to the first site without cutoff
        //cout<<"350~~~~"<<endl;
        //-----end: zip-up algorithm-------------------------
        
        
//        for (int x=0; x<nx; x++) {
//            cout<<"x="<<x;
//            Print(temp1[x]);
//            Print(gate[x]);
//            Print(temp[x]);
//        }
        //cout<<">>>>>>>>>>>>>>>>"<<endl;
        
        if (nx>=3) {
            //-----begin: fitting algorithm----------------------
            temp1=fitting_algorithm(temp1, gate, temp, dcut, ny-1-y, nx, ny);
            //cout<<"365~~~~"<<endl;
            //-----end: fitting algorithm------------------------
        }
        
        
        return temp1;   //return to the output bra=MPO*bra
    }
    //--------end: compressor of ket, starting from the top of the system----------------------------------------------
    
    
    

    

}
#endif

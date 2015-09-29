#ifndef ___PARAMETERS
#define ___PARAMETERS


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
//#include "sweep_utility.h"

#define SQRTJ 0.3

#define COMPLEX 0   //COMPLEX=1: make all tensors' elements complex number
#define NX 7
#define NY 7
#define NS 3
#define NB 3
#define CORRELATIONS 1
#define SEARCH 0
//NEIGHBOR=1: nearst neighbor
#define NEIGHBOR 1
//LATTICE=1: square lattice
#define LATTICE 1
//IGG=0: trivial; IGG=1: Z2
#define IGG 0

#define NORM 0


#define CHIC4 1

using namespace std;
namespace itensor {
    int number_of_terms(){
        if (NEIGHBOR==1&&LATTICE==1&&SEARCH==1) {
            return 2*NX*NY-NX-NY;
        }
        else if(LATTICE==1&&CORRELATIONS==1){
            //return NX+(NX-1)*(NX-3)/8;//number of correlators for 1/8 of sites
            return NX*NY;
        }
    }
    /*
    static int terms=number_of_terms(); //number of terms in hamiltonian
    static int N=NX*NY;
    static int m=9;     //tensor size cutoff in zip-up process: may not be used
    static int dcut=18;  //tensor size cutoff in fitting process
    static double cutoff=1E-10; //truncation error cutoff in the zip-up process, preferred than m
    static double chi=1E-6;   //difference between the fitted MPS and original MPS chi=||fitted_MPS-MPO*MPS||^2
    
    static std::vector<Index> s(N),u(N),d(N),l(N),r(N);
    static std::vector<ITensor> O(N);
    static std::vector<ITensor> a(N); //site tensors
    static std::vector<ITensor> A(NX), B(NX),W(NX*(NY-2));
    static std::vector<Index> U(N),D(N),L(N),R(N);
    static std::vector<ITensor> bra(NX*(NY-1)), ket(NX*(NY-1));
    */
    static int nsweep=100;
    static int nxsweep=1;
    static int nfittingsweep=1;
    static int fitting_process_print=0;
    static double renormalizor=1000000;
    static int terms=3*number_of_terms(); //number of terms in hamiltonian, S*S=S^z*S^z+S^+*S^- two terms in hamiltonian
    static int N=NX*NY;
    //static int cn=(NX-1)+(NX-1)*(NX-3)/8; //number of correlators
    static int m=9;     //tensor size cutoff in zip-up process: may not be used
    static int dcut=18;  //tensor size cutoff in fitting process
    static double cutoff=1E-10; //truncation error cutoff in the zip-up process, preferred than m
    static double chi=1E-6;   //difference between the fitted MPS and original MPS chi=||fitted_MPS-MPO*MPS||^2
    static double _1oversq2=0.7071067812;
    static int printflag=0;//printflag=1 print energy only in scan process, printflag=0 print more details including energy with optimazation position indicators etc
    static int endsweepflag=0;//endsweepflag=1 only print energy at the end of the sweep, otherwise endsweepflag=0
    static int backflag=0;//backflag=0 gives forward sweep, backflag=1 gives backward sweep
    static double scaleup=10;
    
    static std::vector<Index> s(N),u(N),d(N),l(N),r(N);
    static std::vector<ITensor> O(N*terms);
    static std::vector<ITensor> Identity(N);
    static std::vector<ITensor> a(N),areal(N),aimag(N); //site tensors
    static std::vector<ITensor> bh((NX+1)*NY), bv(NX*(NY+1));
    static std::vector<ITensor> aold(N);
    static std::vector<ITensor> W(NX*NY*terms), IW(NX*NY);
    static std::vector<Index> U(N),D(N),L(N),R(N);
    static std::vector<ITensor> bra(NX*(NY-1)*terms), ket(NX*(NY-1)*terms), Ibra(NX*(NY-1)), Iket(NX*(NY-1));
    static std::vector<ITensor> bbra(NY*(NX-1)*terms), bket(NY*(NX-1)*terms), bIbra(NY*(NX-1)), bIket(NY*(NX-1));
    static Index gamma,gammap,delta,deltap;
    //static std::vector<ITensor> AA(NX),BB(NX),WW(NX*(NY-2));
    static std::vector<ITensor> WW(NX*NY),WWW(NX*NY);
    static int mystart, myend;
    static std::vector<ITensor> _Iket(NX), _Ibra(NX);
    static std::vector<ITensor> b_Iket(NY), b_Ibra(NY);
    static std::vector<ITensor> _ket(NX*terms),_bra(NX*terms);
    static std::vector<ITensor> b_ket(NY*terms),b_bra(NY*terms);
    static double startwtime, endwtime;
    static int root=0;
    
    static ITensor IF;
    static std::vector<ITensor> bra_temp(NX), ket_temp(NX);
    static std::vector<ITensor> bbra_temp(NY), bket_temp(NY);
    static ITensor FF;
    static Index delta2=Index("delta2",NB*NB*NS);
    static Index delta2p=Index("delta2p",NB*NB*NS);
    static Index delta3=Index("delta3",NB*NB*NB*NS);
    static Index delta3p=Index("delta3p",NB*NB*NB*NS);
    static Index delta4=Index("delta4",NB*NB*NB*NB*NS);
    static Index delta4p=Index("delta4p",NB*NB*NB*NB*NS);
    
    static Matrix AF,BF,eigenvectors;
    static Vector eigenvalues,evector_real,evector_imag,evector_temp;
    static std::vector<ITensor> IA(NX),AA(NX);
    static std::vector<ITensor> bIA(NY),bAA(NY);
    static int cholesky_flag;
    static Vector dr,di;
    static Matrix vl,vr,diag;
    
    static ITensor temp_real,temp_imag;
    static int nx,ny;
    static int realflag,imagflag;
    static int xmin,xmax;
    static int _switch;
    static int min_index;
    static int skip_flag;
    
    static Vector rrr;
    static Vector iii;
    static Vector beta;
    static Matrix zzl;
    static Matrix zzr;
    static double energy=1000000;
    static int xtemp,ytemp;
}
#endif





















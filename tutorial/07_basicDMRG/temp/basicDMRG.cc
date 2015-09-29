#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "sites/spinone.h"
#include "hams/Heisenberg.h"
//#include <mpi.h>

#define NX 4
#define NY 3
#define NS 3
#define NB 3

using namespace itensor;
using namespace std;

ITensor group_indice2(ITensor Z, Index a1, Index a2, Index a12){
    ITensor temp=ITensor(a12);
    int ma1=a1.m(),ma2=a2.m(),ma12=a12.m();
    if (ma1*ma2!=ma12) {
        cout<<"Error: Tensor sizes do not match!!!"<<endl;
        exit(0);
    }
    for (int i1=1; i1<=ma1; i1++) {
        for (int i2=1; i2<=ma2; i2++) {
            int i12=(i1-1)*ma2+i2;
            temp(a12(i12))=Z(a1(i1),a2(i2));
        }
    }
    return temp;
}

ITensor group_indice4(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34){
    ITensor temp=ITensor(a12,a34);
    int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m();
    if (ma1*ma2!=ma12||ma3*ma4!=ma34) {
        cout<<"Error: Tensor sizes do not match!!!"<<endl;
        exit(0);
    }
    for (int i1=1; i1<=ma1; i1++) {
        for (int i2=1; i2<=ma2; i2++) {
            int i12=(i1-1)*ma2+i2;
            for (int i3=1; i3<=ma3; i3++) {
                for (int i4=1; i4<=ma4; i4++) {
                    int i34=(i3-1)*ma4+i4;
                    temp(a12(i12),a34(i34))=Z(a1(i1),a2(i2),a3(i3),a4(i4));
                }
            }
        }
    }
    return temp;
}

ITensor group_indice6(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5, Index a6, Index a56){
    ITensor temp=ITensor(a12,a34,a56);
    
    int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m(), ma5=a5.m(),ma6=a6.m(),ma56=a56.m();
    if (ma1*ma2!=ma12||ma3*ma4!=ma34||ma5*ma6!=ma56) {
        cout<<"Error: Tensor sizes do not match!!!"<<endl;
        exit(0);
    }
    
    for (int i1=1; i1<=ma1; i1++) {
        for (int i2=1; i2<=ma2; i2++) {
            int i12=(i1-1)*ma2+i2;
            for (int i3=1; i3<=ma3; i3++) {
                for (int i4=1; i4<=ma4; i4++) {
                    int i34=(i3-1)*ma4+i4;
                    for (int i5=1; i5<=ma5; i5++) {
                        for (int i6=1; i6<=ma6; i6++) {
                            int i56=(i5-1)*ma6+i6;
                            temp(a12(i12),a34(i34),a56(i56))=Z(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6));
                        }
                    }
                }
            }
        }
    }
    
    return temp;
}

ITensor group_indice8(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5, Index a6, Index a56, Index a7, Index a8, Index a78){
    ITensor temp=ITensor(a12,a34,a56,a78);
    int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m(), ma5=a5.m(),ma6=a6.m(),ma56=a56.m(), ma7=a7.m(),ma8=a8.m(),ma78=a78.m();
    if (ma1*ma2!=ma12||ma3*ma4!=ma34||ma5*ma6!=ma56||ma7*ma8!=ma78) {
        cout<<"Error: Tensor sizes do not match!!!"<<endl;
        exit(0);
    }
    for (int i1=1; i1<=ma1; i1++) {
        for (int i2=1; i2<=ma2; i2++) {
            int i12=(i1-1)*ma2+i2;
            for (int i3=1; i3<=ma3; i3++) {
                for (int i4=1; i4<=ma4; i4++) {
                    int i34=(i3-1)*ma4+i4;
                    for (int i5=1; i5<=ma5; i5++) {
                        for (int i6=1; i6<=ma6; i6++) {
                            int i56=(i5-1)*ma6+i6;
                            for (int i7=1; i7<=ma7; i7++) {
                                for (int i8=1; i8<=ma8; i8++) {
                                    int i78=(i7-1)*ma8+i8;
                                    temp(a12(i12),a34(i34),a56(i56),a78(i78))=Z(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6),a7(i7),a8(i8));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return temp;
}

int
main(int argc, char* argv[])
{
//    MPI::Init(argc, argv); // initialize MPI environment
//    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
//    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    int N=NX*NY;
    std::vector<Index> s(N),u(N),d(N),l(N),r(N);
    for (int i=0; i<N; i++) {
        s[i]=Index(nameint("s",i),NS,Site);
        u[i]=Index(nameint("u",i),NB);
        d[i]=Index(nameint("d",i),NB);
        l[i]=Index(nameint("l",i),NB);
        r[i]=Index(nameint("r",i),NB);
    }
    for (int i=0; i<N; i++) {
        if ((i+1)%NX!=0) {
            r[i]=l[i+1];
        }
        if (i/NX!=NY-1) {
            d[i]=u[i+NX];
        }
    }
    
    //---------begin: Define tensors O_ss'=<s|O|s'>------------------------------
    ////trial: only one term JS1zS2z in Heisenberg model///
    std::vector<ITensor> O(N);
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
    std::vector<ITensor> a(N); //site tensors
    //std::vector<ITensor> Ap(N); //site tensors
    //std::vector<ITensor> E(N); //expectation value tensors E_uu'dd'll'rr'
    //for (int i=0; i<N; i++) {
    //    E[i].randomize("Complex");/////////////////////////
    //}
    for (int x=0; x<NX; x++) {
        for (int y=0; y<NY; y++) {
            int ind=y*NX+x;
            if (x==0) {
                if (y==0) {
                    a[ind]=ITensor(s[ind],d[ind],r[ind]);
                }
                else if (y==NY-1) {
                    a[ind]=ITensor(s[ind],u[ind],r[ind]);
                }
                else{
                    a[ind]=ITensor(s[ind],u[ind],d[ind],r[ind]);
                }
            }
            else if (x==NX-1) {
                if (y==0) {
                    a[ind]=ITensor(s[ind],d[ind],l[ind]);
                }
                else if (y==NY-1) {
                    a[ind]=ITensor(s[ind],u[ind],l[ind]);
                }
                else{
                    a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind]);
                }
            }
            else if (x>0&&x<NX-1&&y==0) {
                a[ind]=ITensor(s[ind],d[ind],l[ind],r[ind]);
            }
            else if (x>0&&x<NX-1&&y==NY-1) {
                a[ind]=ITensor(s[ind],u[ind],l[ind],r[ind]);
            }
            else{
                a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind],r[ind]);
            }
            a[ind].randomize();
            ////E[ind]=O[ind]*a[ind]*dag(prime(a[ind]));
            //PrintDat(A[ind]);
        }
    }
    //---------end: Define and randomize site tensors A^s_udlr-----------------
    
    //---------begin: group indices to get A and W: A in MPS bra, B in MPS ket, and W in MPO-------------
    std::vector<ITensor> A(NX), B(NX),W(NX*(NY-2));
    std::vector<Index> U(N),D(N),L(N),R(N);
    for (int i=0; i<N; i++) {
        U[i]=Index(nameint("U",i),NB*NB);
        D[i]=Index(nameint("D",i),NB*NB);
        L[i]=Index(nameint("L",i),NB*NB);
        R[i]=Index(nameint("R",i),NB*NB);
    }
    for (int i=0; i<N; i++) {
        if ((i+1)%NX!=0) {
            R[i]=L[i+1];
        }
        if (i/NX!=NY-1) {
            D[i]=U[i+NX];
        }
    }
    
    for (int x=0; x<NX; x++) {
        for (int y=0; y<NY; y++) {
            int ind=y*NX+x;
            int nind=ind-NX;
            if (y==0) {
                if (x==0) {
                    A[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], r[ind], primed(r[ind]), R[ind]);
                }
                else if (x==NX-1) {
                    A[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind]);
                }
                else {
                    A[x]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind], r[ind], primed(r[ind]), R[ind]);
                }
            }
            else if (y==NY-1) {
                if (x==0) {
                    B[x]=group_indice4((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], r[ind], primed(r[ind]), R[ind]);
                }
                else if (x==NX-1) {
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
                else if (x==NX-1) {
                    W[nind]=group_indice6((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind]);
                }
                else {
                    W[nind]=group_indice8((O[ind]*a[ind]*dag(prime(a[ind]))), u[ind], primed(u[ind]), U[ind], d[ind], primed(d[ind]), D[ind], l[ind], primed(l[ind]), L[ind], r[ind], primed(r[ind]), R[ind]);
                }
            }
        }
    }
    //---------end: group indices to get A and W: A in MPS bra, B in MPS ket, and W in MPO-------------
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    //---------test
    for (int x=0; x<NX; x++) {
        int ind=0*NX+x;
        PrintDat(O[ind]*a[ind]*dag(prime(a[ind])));
        PrintDat(A[x]);
    }
    cout<<"<<<<<<<<<<<<<<<<<<<<"<<endl;
    for (int x=0; x<NX; x++) {
        int ind=(NY-1)*NX+x;
        PrintDat(O[ind]*a[ind]*dag(prime(a[ind])));
        PrintDat(B[x]);
    }
    cout<<"<<<<<<<<<<<<<<<<<<<<"<<endl;
    for (int x=0; x<NX; x++) {
        for (int y=1; y<NY-1; y++) {
            int ind=y*NX+x;
            PrintDat(O[ind]*a[ind]*dag(prime(a[ind])));
            PrintDat(W[ind-NX]);
        }
    }
    //---------
    */
    
    
    
    
    
    
    
    
    
    
    //////PrintDat(A[1]);
    //for (int i=0;i<N; i++) {
        //A[i].randomize("Complex");
        //PrintDat(O[i]);
    //}
    
    //PrintDat(A[1]);
    
    
    //PrintDat(A[0]*A[2]);
    //int N=3;
    
    //SpinOne model(N);
    //MPS psi(model);
    
    
    ///PrintDat(psi);
    /*
    std::vector<ITensor> aa(6);
    Index s1("site 1", 2, Site), b1("bond 1",3);
    aa[1]=ITensor(s1,b1);
    aa[1].randomize("Complex");
    PrintDat(aa[1]);

    PrintDat(psi.A(1));
    PrintDat(psi.A(2));
     */
    
    //Index s1("site 1", 2, Site), b1("bond 1",3);
    //Index s2("site 2", 2, Site);
    //ITensor a1(s1,b1),a2(b1,s2);
    //a1.randomize("Complex");
    //a2.randomize("Complex");
    
    //PrintDat(a1);
    
    //PrintDat(a1.takeRealPart());
    //PrintDat(psi.A(1));
    //PrintDat(psi.A(2));
    
    //ITensor bondWF = psi.A(1)*psi.A(2);
    //PrintDat(bondWF);
    //cout<<psi.A(1)<<endl;
    //cout<<"Hello!"<<endl;
    
    
    return 0;
}




///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
/*
std::vector<int> siteindex(int i){
    std::vector<int> temp(2);
    temp[0]=i%NX;
    temp[1]=i/NX;
    return temp;
}
std::vector<int> bondhindex(int i){
    std::vector<int> temp(2);
    temp[0]=i%(NX-1);
    temp[1]=i/(NX-1);
    return temp;
}
std::vector<int> bondvindex(int i){
    std::vector<int> temp(2);
    temp[0]=i%NX;
    temp[1]=i/NX;
    return temp;
}
*/

/*int N=NX*NY;
int Bnh=(NX-1)*NY,Bnv=NX*(NY-1);
std::vector<Index> s(N);
////Index ss("sites", NS, Site);
for (int i=0; i<N; i++) {
    ////s[i]=ss;
    s[i]=Index(nameint("s",i),NS,Site);
}

std::vector<Index> bh(Bnh);
std::vector<Index> bv(Bnv);
//Index bb("bonds", NB);
for (int i=0; i<Bnh; i++) {
    ///bh[i]=bb;
    bh[i]=Index(nameint("bondh",i),NB);
}
for (int i=0; i<Bnv; i++) {
    bv[i]=Index(nameint("bondv",i),NB);
}

//---------begin: Define and randomize all the site tensors A^s_udlr-----------
std::vector<ITensor> A(N); //site tensors
for (int i=0; i<N; i++) {
    int xx=siteindex(i)[0];
    int yy=siteindex(i)[1];
    
    if (xx==0&&yy==0) {
        A[i]=ITensor(s[i],bv[i],bh[i-yy]);
    }
    if (xx==NX-1&&yy==0) {
        A[i]=ITensor(s[i],bv[i],bh[i-1-yy]);
    }
    if (xx==0&&yy==NY-1) {
        A[i]=ITensor(s[i],bv[i-NX],bh[i-yy]);
    }
    if (xx==NX-1&&yy==NY-1) {
        A[i]=ITensor(s[i],bv[i-NX],bh[i-1-yy]);
    }
    if (xx>0&&xx<NX-1) {
        if (yy==0) {
            A[i]=ITensor(s[i],bv[i],bh[i-1-yy],bh[i-yy]);
        }
        if (yy==NY-1) {
            A[i]=ITensor(s[i],bv[i-NX],bh[i-1-yy],bh[i-yy]);
        }
    }
    if (yy>0&&yy<NY-1) {
        if (xx==0) {
            A[i]=ITensor(s[i],bv[i-NX],bv[i],bh[i-yy]);
        }
        if (xx==NX-1) {
            A[i]=ITensor(s[i],bv[i-NX],bv[i],bh[i-1-yy]);
        }
    }
    if (xx>0&&xx<NX-1&&yy>0&&yy<NY-1) {
        A[i]=ITensor(s[i],bv[i-NX],bv[i],bh[i-1-yy],bh[i-yy]);
    }
    A[i].randomize("Complex");
}
*/
//---------end: Define and randomize all the site tensors A^s_udlr-----------
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


/*
#include "symMonteCarlo.h"

#include "core.h"

#include "model/tjmod.h"
#include "hams/tJham.h"
//#include "model/hubbardmod.h"
//#include "hams/ExtendedHubbardham.h"

//LATTICE=0 indicates Honeycomb lattice
//LATTICE=1 indicates Triangular lattice
//#define LATTICE 1
//SYSSHAPE=0 indicates system has rhombus shape
//SYSSHAPE=1 indicates system has hexagonal shape
//#define SYSSHAPE 0
#define RANDOM 1
#define NUMBERING 1


using boost::format;
using namespace std;

typedef tJ
//typedef Hubbard
states;


//number of electrons
const extern int Nferm=8;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm=500;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure=1000;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure=40;

int
main(int argc, char* argv[])
{
    MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    
    int N=16;
    
	//
    // Initialize the site degrees of freedom.
    //
    states model(N);    // make a chain of N spins
    
    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    const IQMPO H = tJChain(model,Opt("J",0.5)&Opt("t",-1.0));
    //const IQMPO H = ExtendedHubbard(model,Opt("U",10.0)& Opt("t1",1)& Opt("t2",0.0)& Opt("V1",0));
	
    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(model);
	//up and down spins
    ITiVec init_state_vec=initialize_random_state_vec_tj(N,Nferm);
    //ITiVec init_state_vec=initialize_random_state_vec_hubbard(N,Nferm);
    if (RANDOM==0) {

        init_state_vec(0)=1;
        init_state_vec(1)=2;
        init_state_vec(2)=1;
        init_state_vec(3)=3;
        init_state_vec(4)=2;
        init_state_vec(5)=1;
        init_state_vec(6)=2;
        init_state_vec(7)=3;
        init_state_vec(8)=3;
        init_state_vec(9)=1;
        init_state_vec(10)=2;
        init_state_vec(11)=1;
        init_state_vec(12)=1;
        init_state_vec(13)=1;
        init_state_vec(14)=3;
        init_state_vec(15)=1;
    }
	
	if (rank==0) {
		cout<<"initial state: "<<init_state_vec<<endl;
	}
    
    for(int i=1;i<=N;i++)
    {
		if(init_state_vec(i-1)==1)
		{
			initState.set(i,"Emp");
		}
		else if(init_state_vec(i-1)==2)
		{
			initState.set(i,"Up");
		}
		else if(init_state_vec(i-1)==3)
		{
			initState.set(i,"Dn");
		}
		else if(init_state_vec(i-1)==4)
		{
			initState.set(i,"UpDn");
		}
		else
		{
			cout<<"InitState Error"<<endl;
			exit(0);
		}
    }
    
    IQMPS psi(initState);
	
    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
	if (rank==0) {
		cout << totalQN(psi) << endl;
		cout << format("Initial energy = %.10f") % psiHphi(psi,H,psi) << endl;
	}
	
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    
	int Nsweeps=75;
	Sweeps sweeps(Nsweeps);
	sweeps.maxm() = 10,20,20,20,20,20,50,50,50,100,100,100,100,200,200,200,200,200,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,4000,4000,4000,4000,4000,4000,4000,4000,4000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000;
	sweeps.cutoff() = 1E-10,1E-10,0.;
	sweeps.niter() = 2;
	sweeps.noise() = 1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.;
	cout << sweeps;
    
    //
    // Begin the DMRG calculation
    //
	
	
    Real En = dmrg(psi,H,sweeps,Quiet());
    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;
    cout  <<"Energy: " << psiHphi(psi,H,psi)<<endl;
    
    
    
	psi.position(1);
	if (rank==0) { cout << "Norm2: " << Dot(conj(psi.A(1)), psi.A(1))<<endl; }
	cout.flush();
    
	ITdVec nup_vec(N),ndn_vec(N);
	
	for(int j=1; j<=N; j++) {
		IQTensor nup_op = model.op("Nup",j);
		IQTensor ndn_op = model.op("Ndn",j);
		psi.position(j);
		IQTensor ket = psi.A(j);
		IQTensor bra = conj(primed(ket,Site));
		nup_vec(j-1) = Dot(bra, nup_op*ket);
		ndn_vec(j-1) = Dot(bra, ndn_op*ket);
	}
	
	if (rank==0) {
        cout<<"Nup: "<<setprecision(10)<<nup_vec<<endl;
        cout<<"Ndn: "<<setprecision(10)<<ndn_vec<<endl;
        cout<<"Total N: "<<setprecision(10)<<nup_vec+ndn_vec<<endl;
        cout<<"Sz: "<<setprecision(10)<<0.5*(nup_vec-ndn_vec)<<endl;
	}
	
	psi.position(1);
    
    
    cout<<"Writing to disk:"<<endl;
  	writeToFile("psi_dat",psi);
	writeToFile("model_dat",model);
	writeToFile("Ham_dat",H);
    cout<<"Writing done."<<endl;
    
	ITiVec inv_sym_mapping_table, T1_sym_mapping_table_16_site, T2_sym_mapping_table_16_site, C6_sym_mapping_table_16_site;
    if (NUMBERING==0) {
        inv_sym_mapping_table="8 11 10 9 4 7 6 5 0 3 2 1 12 15 14 13";
        T1_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
        T2_sym_mapping_table_16_site="12 13 14 15 0 1 2 3 4 5 6 7 8 9 10 11";
        C6_sym_mapping_table_16_site="9 5 1 13 14 10 6 2 3 15 11 7 4 0 12 8";
    }
    if (NUMBERING==1) {
        inv_sym_mapping_table="0 3 2 1 12 15 14 13 8 11 10 9 4 7 6 5";
        T1_sym_mapping_table_16_site="4 5 6 7 8 9 10 11 12 13 14 15 0 1 2 3";
        T2_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
        C6_sym_mapping_table_16_site="2 15 8 5 3 12 9 6 0 13 10 7 1 14 11 4";
    }
	ITiVecArray C6_generator_mapping_table(1);
	C6_generator_mapping_table(0)=C6_sym_mapping_table_16_site;
	ITcVec Ang_mom(1);
	Ang_mom(0)=std::exp(Complex_i*2.*Pi/3.);
	//Ang_mom(0)=1.;
	
    
	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_16_site;
	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_16_site;
	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_16_site;
    //	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_32_site;
    //	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_32_site;
    //	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_32_site;
    
	ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table;
	translation_generator_mapping_table(1)=T2_sym_mapping_table;
	ITcVec COM_mom="1. 1.";
	
	ITiVecArray transl_C6_generator_mapping_table(3);
	transl_C6_generator_mapping_table(0)=T1_sym_mapping_table;
	transl_C6_generator_mapping_table(1)=T2_sym_mapping_table;
	transl_C6_generator_mapping_table(2)=C6_sym_mapping_table;
	
	ITcVec COM_Ang_mom="1. 1. 0.";
	COM_Ang_mom(2)=std::exp(-Complex_i*2.*Pi/3.);
	
	std::string desc("16_J075_m2P3");
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table,translation_generator_mapping_table,COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
    //cout<<"T1<<<<<<<<<<<<<<<<<<<"<<endl;
    //cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(desc,psi,T1_sym_mapping_table)<<endl;
    //cout<<"T2<<<<<<<<<<<<<<<<<<<"<<endl;
    //cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(desc,psi,T2_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_tj_MonteCarlo(psi,C6_generator_mapping_table,Ang_mom)<<endl;
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_hubbard_MonteCarlo(psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_hubbard_MonteCarlo(psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi, C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, C6_sym_mapping_table_8_site, C6_generator_mapping_table, Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,T1_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, T1_sym_mapping_table_8_site, translation_generator_mapping_table, COM_mom)<<endl;
	
	ITiVec triv_sym_mapping_table="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15";
	ITiVecArray triv_generator_mapping_table(1);
	triv_generator_mapping_table(0)=triv_sym_mapping_table;
	ITcVec triv_mom="1.";
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,triv_generator_mapping_table,triv_mom)<<endl;
    
	
    MPI::Finalize();
    
    return 0;
}

*/



/*
int 
main(int argc, char* argv[])
    {
    int N = 6;

    SpinOne sites(N);

    MPO H = Heisenberg(sites);

    Sweeps sweeps(5); //number of sweeps is 5
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
        
    MPS psi(sites);

    Real energy = dmrg(psi,H,sweeps);

    cout << "Ground State Energy = " << energy << endl;

    return 0;
}
*/
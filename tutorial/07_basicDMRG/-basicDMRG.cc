#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "iqcombiner.h"
#include "sweeps.h"
#include "sites/spinone.h"
#include "hams/Heisenberg.h"
#include <mpi.h>
#include "sweep_utility.h"
#include "parameters.h"
#include "tensor_utility.h"
#include "tensor_sweep.h"
#include "heisenberg_sl.h"
#include "main_loop.h"

#include "generalizedev.h"

using namespace itensor;
using namespace std;


int
main(int argc, char* argv[])
{
    
    MPI_Status status;
    MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process

    //int root=0;
    if (rank==root) {
        MPI_Barrier(MPI_COMM_WORLD);
        startwtime = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        startwtime = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    
    if (rank==root) {
        cout<<"System size: "<<NX<<"*"<<NY<<endl;
        cout<<"Number of terms (including Identity operator)="<<terms+1<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if (rank==root) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else {   //use root processor to do Identity operator related calculations, other processors to do hamiltonian related calculations
        mystart = (terms / (size-1)) * (rank-1);
        if ((terms) % (size-1) > (rank-1)){
            mystart += (rank-1);
            myend = mystart + ((terms) / (size-1)) + 1;
        }else{
            mystart += (terms) % (size-1);
            myend = mystart + ((terms) / (size-1));
        }
        printf("CPU%d %d ~ %d\n",rank,mystart,myend);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    //-----begin: setup random site tensors---------------------------------------
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
    
    //---------begin: Define, randomize and broadcast site tensors a^s_udlr-----------------
    
    //---the hidden part of code is for randomize REAL a tensors
    if (COMPLEX==0) {
    for (int x=0; x<NX; x++) {
        for (int y=0; y<NY; y++) {
            int ind=y*NX+x;
            if (x==0) {
                if (y==0) {
                    a[ind]=ITensor(s[ind],d[ind],r[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else if (y==NY-1) {
                    a[ind]=ITensor(s[ind],u[ind],r[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else{
                    a[ind]=ITensor(s[ind],u[ind],d[ind],r[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                    MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            else if (x==NX-1) {
                if (y==0) {
                    a[ind]=ITensor(s[ind],d[ind],l[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else if (y==NY-1) {
                    a[ind]=ITensor(s[ind],u[ind],l[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else{
                    a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind]);
                    if (rank==0) {
                        a[ind].randomize();
                        if (NORM==1) {
                            a[ind] /= a[ind].norm();
                        }
                    }
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                    MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            else if (x>0&&x<NX-1&&y==0) {
                a[ind]=ITensor(s[ind],d[ind],l[ind],r[ind]);
                if (rank==0) {
                    a[ind].randomize();
                    if (NORM==1) {
                        a[ind] /= a[ind].norm();
                    }
                }
                for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
            else if (x>0&&x<NX-1&&y==NY-1) {
                a[ind]=ITensor(s[ind],u[ind],l[ind],r[ind]);
                if (rank==0) {
                    a[ind].randomize();
                    if (NORM==1) {
                        a[ind] /= a[ind].norm();
                    }
                }
                for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
            else{
                a[ind]=ITensor(s[ind],u[ind],d[ind],l[ind],r[ind]);
                if (rank==0) {
                    a[ind].randomize();
                    if (NORM==1) {
                        a[ind] /= a[ind].norm();
                    }
                }
                for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                for (int i4=1; i4<=a[ind].indices()[4].m(); i4++) {
                                    MPI_Bcast(&(a[ind](a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3),a[ind].indices()[4](i4))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            //PrintDat(a[ind]);
        }
    }
    }
    //---the hidden part of code is for randomize COMPLEX a tensors
    if (COMPLEX==1) {
    for (int x=0; x<NX; x++) {
        for (int y=0; y<NY; y++) {
            int ind=y*NX+x;
            if (x==0) {
                if (y==0) {
                    areal[ind]=ITensor(s[ind],d[ind],r[ind]);
                    aimag[ind]=ITensor(s[ind],d[ind],r[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else if (y==NY-1) {
                    areal[ind]=ITensor(s[ind],u[ind],r[ind]);
                    aimag[ind]=ITensor(s[ind],u[ind],r[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else{
                    areal[ind]=ITensor(s[ind],u[ind],d[ind],r[ind]);
                    aimag[ind]=ITensor(s[ind],u[ind],d[ind],r[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=areal[ind].indices()[3].m(); i3++) {
                                    MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2),areal[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                    MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2),aimag[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            else if (x==NX-1) {
                if (y==0) {
                    areal[ind]=ITensor(s[ind],d[ind],l[ind]);
                    aimag[ind]=ITensor(s[ind],d[ind],l[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else if (y==NY-1) {
                    areal[ind]=ITensor(s[ind],u[ind],l[ind]);
                    aimag[ind]=ITensor(s[ind],u[ind],l[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
                else{
                    areal[ind]=ITensor(s[ind],u[ind],d[ind],l[ind]);
                    aimag[ind]=ITensor(s[ind],u[ind],d[ind],l[ind]);
                    if (rank==0) {
                        areal[ind].randomize();
                        aimag[ind].randomize();
                    }
                    for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=areal[ind].indices()[3].m(); i3++) {
                                    MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2),areal[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                    MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2),aimag[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            else if (x>0&&x<NX-1&&y==0) {
                areal[ind]=ITensor(s[ind],d[ind],l[ind],r[ind]);
                aimag[ind]=ITensor(s[ind],d[ind],l[ind],r[ind]);
                if (rank==0) {
                    areal[ind].randomize();
                    aimag[ind].randomize();
                }
                for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=areal[ind].indices()[3].m(); i3++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2),areal[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2),aimag[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
            else if (x>0&&x<NX-1&&y==NY-1) {
                areal[ind]=ITensor(s[ind],u[ind],l[ind],r[ind]);
                aimag[ind]=ITensor(s[ind],u[ind],l[ind],r[ind]);
                if (rank==0) {
                    areal[ind].randomize();
                    aimag[ind].randomize();
                }
                for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=areal[ind].indices()[3].m(); i3++) {
                                MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2),areal[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2),aimag[ind].indices()[3](i3))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
            else{
                areal[ind]=ITensor(s[ind],u[ind],d[ind],l[ind],r[ind]);
                aimag[ind]=ITensor(s[ind],u[ind],d[ind],l[ind],r[ind]);
                if (rank==0) {
                    areal[ind].randomize();
                    aimag[ind].randomize();
                }
                for (int i0=1; i0<=areal[ind].indices()[0].m(); i0++) {
                    for (int i1=1; i1<=areal[ind].indices()[1].m(); i1++) {
                        for (int i2=1; i2<=areal[ind].indices()[2].m(); i2++) {
                            for (int i3=1; i3<=areal[ind].indices()[3].m(); i3++) {
                                for (int i4=1; i4<=areal[ind].indices()[4].m(); i4++) {
                                    MPI_Bcast(&(areal[ind](areal[ind].indices()[0](i0),areal[ind].indices()[1](i1),areal[ind].indices()[2](i2),areal[ind].indices()[3](i3),areal[ind].indices()[4](i4))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                    MPI_Bcast(&(aimag[ind](aimag[ind].indices()[0](i0),aimag[ind].indices()[1](i1),aimag[ind].indices()[2](i2),aimag[ind].indices()[3](i3),aimag[ind].indices()[4](i4))), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                                }
                            }
                        }
                    }
                }
            }
            //PrintDat(a[ind]);
            a[ind]=areal[ind]+Complex_i*aimag[ind];
            if (NORM==1) {
                a[ind] /= a[ind].norm();
            }
        }
    }
    }
    //---------end: Define, randomize and broadcast site tensors a^s_udlr-----------------
    

    
    backflag=0;
    if (rank==root) {
        setup_identity_tensors(NX,NY);   //find all IA, IB, IW
        identity_find_bras(IW,NX,NY);  //find all Ibra's
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
        for (int i=mystart; i<myend; i++) {
            setup_ham_tensors(i, NX,NY);   //find all A, B, W
            for (int x=0; x<NX*NY; x++) {
                WW[x]=W[x+NX*NY*i];
            }
            ham_find_bras(WW, i, NX,NY);   //find all bra's
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    //----begin: main loop, sweep---------------------------------------------------------------------------------
    for (int sweep_ind=1; sweep_ind<=nsweep; sweep_ind++) {//int sweep_ind=1; sweep_ind<=nsweep; sweep_ind++
        main_loop(NX, NY, backflag, root, sweep_ind, MPI_COMM_WORLD);
    }
    
    //----end: main loop, sweep---------------------------------------------------------------------------------
    if (rank==root) {
        for (int i=0; i<N; i++) {
            cout<<"i="<<i;
            PrintDat(a[i]);
        }
        if (NX==2&&NY==2) {
            PrintDat(a[0]*a[1]*a[2]*a[3]);
        }
    }
    
	
    if (rank==root) {
        endwtime=MPI_Wtime();
        printf("wall clock time = %f\n", endwtime-startwtime);
    }

    MPI::Finalize();
   
    
    return 0;
}










































///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
/*//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



/*//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//---------begin: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
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
//---------end: group indices to get A, B and W: A in MPS bra, B in MPS ket, and W in MPO-------------
*///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
*///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
























































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

/*
if (rank==root) {
    setup_identity_tensors();   //find all IA, IB, IW
    //cout<<"this is end signal from rank "<<rank<<"!"<<endl;
}
else {
    for (int i=mystart; i<myend; i++) {
        setup_ham_tensors(i);   //find all A, B, W
        //cout<<"this is end signal from rank "<<rank<<" and i="<<i<<"!"<<endl;
    }
}


std::vector<ITensor> _Iket(NX), _Ibra(NX),_ket(NX*terms),_bra(NX*terms); //this is the temp ket/bra used in the loop, obtained from Iket/Ibra
//---begin: find bra's and ket(0)-----------------------
if (rank==root) {
    identity_find_bras(IA, IB, IW);  //find all Ibra's
    IA=ket_otho_center_1(IA, 0);//change the othogonality center of A to site 1
    for (int x=0; x<NX; x++) {  //find Iket(0);
        Iket[x]=IA[x];
    }
    cout<<"this is end signal from rank "<<rank<<"!"<<endl;
    }
    else {
        std::vector<ITensor> AA(NX),BB(NX),WW(NX*(NY-2));
        for (int i=mystart; i<myend; i++) {
            //if (i==61) {/////////////
            
            for (int x=0; x<NX; x++) {
                AA[x]=A[x+NX*i];
                BB[x]=B[x+NX*i];
            }
            for (int x=0; x<NX*(NY-2); x++) {
                WW[x]=W[x+NX*(NY-2)*i];
            }
            
            ham_find_bras(AA, BB, WW, i);   //find all bra's
            
            AA=ket_otho_center_1(AA, 0);
            for (int x=0; x<NX; x++) {  //find ket(0)
                ket[x+NX*(NY-1)*i]=AA[x];
            }
            
            cout<<"this is end signal from rank "<<rank<<" and i="<<i<<"!"<<endl;
            //}/////////////
        }
    }
    //---end: find bra's and ket(0)-----------------------

*/

/*
if (rank==root) {
    for (int x=0; x<NX; x++) {
        if (Y<NY-1) {
            
            if (Ibra[x+NX*(NY-2-Y)].r()==2) {
                _Ibra[x]=ITensor(Ibra[x+NX*(NY-2-Y)].indices()[0], Ibra[x+NX*(NY-2-Y)].indices()[1]);
            }
            else if (Ibra[x+NX*(NY-2-Y)].r()==3){
                _Ibra[x]=ITensor(Ibra[x+NX*(NY-2-Y)].indices()[0], Ibra[x+NX*(NY-2-Y)].indices()[1], Ibra[x+NX*(NY-2-Y)].indices()[2]);
            }
            
            _Ibra[x]=Ibra[x+NX*(NY-2-Y)];
        }
        if (Y>0) {
            
            if (Iket[x+NX*(Y-1)].r()==2) {
                _Iket[x]=ITensor(Iket[x+NX*(Y-1)].indices()[0], Iket[x+NX*(Y-1)].indices()[1]);
            }
            else if (Iket[x+NX*(Y-1)].r()==3){
                _Iket[x]=ITensor(Iket[x+NX*(Y-1)].indices()[0], Iket[x+NX*(Y-1)].indices()[1], Iket[x+NX*(Y-1)].indices()[2]);
            }
            
            _Iket[x]=Iket[x+NX*(Y-1)];
        }
    }
    
    //                for (int x=0; x<NX; x++) {
    //                    cout<<"rank="<<rank<<", Y="<<Y<<", x="<<x;
    //                    Print(_Iket[x]);
    //                }
    
    MPI_Barrier(MPI_COMM_WORLD);
}
else {
    for (int i=mystart; i<myend; i++) {
        for (int x=0; x<NX; x++) {
            if (Y<NY-1) {
                
                if (bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].r()==2) {
                    _bra[x+NX*i]=ITensor(bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].indices()[0], bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].indices()[1]);
                }
                else if (bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].r()==3){
                    _bra[x+NX*i]=ITensor(bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].indices()[0], bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].indices()[1], bra[x+NX*(NY-2-Y)+NX*(NY-1)*i].indices()[2]);
                }
                
                _bra[x+NX*i]=bra[x+NX*(NY-2-Y)+NX*(NY-1)*i];
            }
            if (Y>0) {
                
                if (ket[x+NX*(Y-1)+NX*(NY-1)*i].r()==2) {
                    _ket[x+NX*i]=ITensor(ket[x+NX*(Y-1)+NX*(NY-1)*i].indices()[0], ket[x+NX*(Y-1)+NX*(NY-1)*i].indices()[1]);
                }
                else if (ket[x+NX*(Y-1)+NX*(NY-1)*i].r()==3){
                    _ket[x+NX*i]=ITensor(ket[x+NX*(Y-1)+NX*(NY-1)*i].indices()[0], ket[x+NX*(Y-1)+NX*(NY-1)*i].indices()[1], ket[x+NX*(Y-1)+NX*(NY-1)*i].indices()[2]);
                }
                
                _ket[x+NX*i]=ket[x+NX*(Y-1)+NX*(NY-1)*i];
            }
            
            //if (i==mystart) {
            //    cout<<"rank="<<rank<<", Y="<<Y<<", x="<<x;
            //    Print(_bra[x+NX*i]);
            //    Print(_ket[x+NX*i]);
            //}
            
        }
        
        //                    for (int x=0; x<NX; x++) {
        //                        cout<<"rank="<<rank<<", i="<<i<<", Y="<<Y<<", x="<<x;
        //                        Print(_ket[x+NX*i]);
        //                    }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}
*/

/*
if (rank==root) {
    if (Y==0) {
        //std::vector<ITensor> IA(NX);
        for (int x=0; x<NX; x++) {
            IA[x]=IW[x];
        }
        
        IA=ket_otho_center_1(IA, 0);//change the othogonality center of A to site 1
        for (int x=0; x<NX; x++) {  //find Iket(0);
            Iket[x+NX*Y]=IA[x];
        }
    }
    else if (Y<NY-1&&Y>0){
        std::vector<ITensor> ket_temp(NX),gate_temp(NX),temp2(NX);
        for (int x=0; x<NX; x++) {
            ket_temp[x]=Iket[x+NX*(Y-1)];
            gate_temp[x]=IW[x+NX*Y];
        }
        temp2=compressor1(ket_temp,gate_temp,Y,cutoff, m, dcut);
        for (int x=0; x<NX; x++) {
            Iket[x+NX*Y]=temp2[x];
        }
    }
    
    //                for (int x=0; x<NX; x++) {
    //                    cout<<"rank="<<rank<<", Y="<<Y<<", x="<<x;
    //                    Print(Iket[x+NX*Y]);
    //                }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    //endwtime=MPI_Wtime();
    //printf("wall clock time = %f of processor %d\n", endwtime-startwtime, rank);
}
else {
    for (int i=mystart; i<myend; i++) {
        
        //cout<<"rank="<<rank<<", i="<<i<<", Y="<<Y<<endl;
        if (i==mystart) {
            if (Y==0) {
                //std::vector<ITensor> AA(NX);
                for (int x=0; x<NX; x++) {
                    AA[x]=W[x+NX*NY*i];
                }
                
                AA=ket_otho_center_1(AA, 0);
                for (int x=0; x<NX; x++) {  //find ket(0)
                    ket[x+NX*Y+NX*(NY-1)*i]=AA[x];
                }
                
                cout<<"rank="<<rank<<", i="<<i<<", Y="<<Y<<endl;
            }
            else if (Y<NY-1&&Y>0){
                std::vector<ITensor> ket_temp(NX),gate_temp(NX),temp2(NX);
                for (int x=0; x<NX; x++) {
                    ket_temp[x]=ket[x+NX*(Y-1)+NX*(NY-1)*i];
                    gate_temp[x]=W[x+NX*Y+N*i];
                }
                temp2=compressor1(ket_temp,gate_temp,Y,cutoff, m, dcut);
                for (int x=0; x<NX; x++) {
                    ket[x+NX*Y+NX*(NY-1)*i]=temp2[x];
                }
            }
        }
                
                
                //cout<<"rank="<<rank<<", i="<<i<<", Y="<<Y<<endl;
                
                //                    for (int x=0; x<NX; x++) {
                //                        cout<<"rank="<<rank<<", i="<<i<<", Y="<<Y<<", x="<<x;
                //                        Print(ket[x+NX*Y+NX*(NY-1)*i]);
                //                    }
                
                
                }
                MPI_Barrier(MPI_COMM_WORLD);
                //endwtime=MPI_Wtime();
                //printf("wall clock time = %f of processor %d\n", endwtime-startwtime, rank);
            }
*/
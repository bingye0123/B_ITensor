//fitting.h


#ifndef ___MAIN_LOOP
#define ___MAIN_LOOP


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
#include <math.h>

#include "generalizedev.h"


using namespace std;
namespace itensor {
    void main_loop(int nx, int ny, int _backflag, int root, int sweep_ind, MPI_Comm communicator){
        
        int rank;
        MPI_Comm_rank(communicator, &rank);
        int size;
        MPI_Comm_size(communicator, &size);
        MPI_Status status;
        
        //----begin: row index looping, from Y=0 to Y=ny-1------
        for (int Y=0; Y<ny; Y++) {    //int Y=0; Y<ny; Y++ //int Y=1; Y<ny-1; Y++
            if (rank==root && printflag==0) {
                cout<<"Sweep "<<sweep_ind<<", Y="<<Y<<endl;
            }
            //----begin: build up _bra and _ket for the loop------------------------------------
            if (rank==root) {
                if (Y<ny-1) {
                    for (int x=0; x<nx; x++) {
                        if (_backflag==0) {
                            _Ibra[x]=Ibra[x+nx*(ny-2-Y)];
                        }
                        else {
                            b_Ibra[x]=bIbra[x+nx*(ny-2-Y)];
                        }
                    }
                }
                if (Y>0) {
                    if (Y==1) {
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                IA[x]=IW[x];
                            }
                            else {
                                bIA[x]=IW[x];
                            }
                        }
                        if (_backflag==0) {
                            IA=ket_otho_center_1(IA, Y-1, nx, ny);//change the othogonality center of A to site 1
                        }
                        else {
                            bIA=ket_otho_center_1(bIA, Y-1, nx, ny);//change the othogonality center of A to site 1
                        }
                        for (int x=0; x<nx; x++) {  //find Iket(0);
                            if (_backflag==0) {
                                Iket[x+nx*(Y-1)]=IA[x];
                                _Iket[x]=IA[x];
                            }
                            else {
                                bIket[x+nx*(Y-1)]=bIA[x];
                                b_Iket[x]=bIA[x];
                            }
                        }
                    }
                    else {
                        std::vector<ITensor> ket_temp(nx),gate_temp(nx),temp2(nx);
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                ket_temp[x]=Iket[x+nx*(Y-2)];
                            }
                            else {
                                ket_temp[x]=bIket[x+nx*(Y-2)];
                            }
                            gate_temp[x]=IW[x+nx*(Y-1)];
                        }
                        temp2=compressor1(ket_temp,gate_temp,Y-1,cutoff, m, dcut,nx,ny);
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                Iket[x+nx*(Y-1)]=temp2[x];
                                _Iket[x]=temp2[x];
                            }
                            else {
                                bIket[x+nx*(Y-1)]=temp2[x];
                                b_Iket[x]=temp2[x];
                            }
                            
                        }
                    }
                }
                
                MPI_Barrier(communicator);
            }
            else {
                //if (Y!=2) {//~~~~~
                
                for (int i=mystart; i<myend; i++) {
                    //cout<<"<><><><><><><>i="<<i<<endl;
                    
                    if (Y<ny-1) {
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                _bra[x+nx*i]=bra[x+nx*(ny-2-Y)+nx*(ny-1)*i];
                            }
                            else {
                                b_bra[x+nx*i]=bbra[x+nx*(ny-2-Y)+nx*(ny-1)*i];
                            }
                        }
                    }
                    if (Y>0) {
                        if (Y==1) {
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    AA[x]=W[x+nx*ny*i];
                                }
                                else {
                                    bAA[x]=W[x+nx*ny*i];
                                }
                            }
                            if (_backflag==0) {
                                AA=ket_otho_center_1(AA, Y-1,nx,ny);
                            }
                            else {
                                bAA=ket_otho_center_1(bAA, Y-1,nx,ny);
                            }
                            for (int x=0; x<nx; x++) {  //find ket(0)
                                if (_backflag==0) {
                                    ket[x+nx*(Y-1)+nx*(ny-1)*i]=AA[x];
                                    _ket[x+nx*i]=AA[x];
                                }
                                else {
                                    bket[x+nx*(Y-1)+nx*(ny-1)*i]=bAA[x];
                                    b_ket[x+nx*i]=bAA[x];
                                }
                            }
                        }
                        else {
                            std::vector<ITensor> ket_temp(nx),gate_temp(nx),temp2(nx);
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    ket_temp[x]=ket[x+nx*(Y-2)+nx*(ny-1)*i];
                                }
                                else {
                                    ket_temp[x]=bket[x+nx*(Y-2)+nx*(ny-1)*i];
                                }
                                gate_temp[x]=W[x+nx*(Y-1)+N*i];
                            }
                            temp2=compressor1(ket_temp,gate_temp,Y-1,cutoff, m, dcut,nx,ny);
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    ket[x+nx*(Y-1)+nx*(ny-1)*i]=temp2[x];
                                    _ket[x+nx*i]=temp2[x];
                                }
                                else {
                                    bket[x+nx*(Y-1)+nx*(ny-1)*i]=temp2[x];
                                    b_ket[x+nx*i]=temp2[x];
                                }
                            }
                        }
                        
                    }
                    
                    
                }
                
                //}//~~~~~~~
                MPI_Barrier(communicator);
            }
            //----end: build up _bra and _ket for the loop------------------------------------
            
            //if (Y==2) {//~~~~~
            
            
            //----begin: sweep of a row------------------------
            for (int xsweep=1; xsweep<=nxsweep; xsweep++) {
                for (int X=0; X<nx; X++) {  //int X=0; X<nx; X++
                    
                    //skip_flag=0;
                    
                    std::vector<ITensor> F(size);
                    if (rank==root) {
                        if (printflag==0) {
                            cout<<"--------------------------------------"<<endl;
                            cout<<"Working on ("<<X<<","<<Y<<"):"<<endl;
                            //PrintDat(a[0]);//~~~~~~~
                        }
                        if (_backflag==0) {
                            IF=sweepx(_Ibra,_Iket,X,Y,rank,0,nx,ny);
                        }
                        else {
                            IF=sweepx(b_Ibra,b_Iket,X,Y,rank,0,nx,ny);
                        }
                        for (int xxx=0; xxx<size; xxx++) {
                            F[xxx]=IF;
                        }
                        MPI_Barrier(communicator);
                    }
                    else {
                        for (int i=mystart; i<myend; i++) {
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    bra_temp[x]=_bra[x+nx*i];
                                    ket_temp[x]=_ket[x+nx*i];
                                }
                                else {
                                    bbra_temp[x]=b_bra[x+nx*i];
                                    bket_temp[x]=b_ket[x+nx*i];
                                }
                            }
                            if (i==mystart) {
                                if (_backflag==0) {
                                    FF=sweepx(bra_temp,ket_temp,X,Y,rank,i,nx,ny);
                                }
                                else {
                                    FF=sweepx(bbra_temp,bket_temp,X,Y,rank,i,nx,ny);
                                }
                            }
                            else {
                                if (_backflag==0) {
                                    FF+=sweepx(bra_temp,ket_temp,X,Y,rank,i,nx,ny);
                                }
                                else {
                                    FF+=sweepx(bbra_temp,bket_temp,X,Y,rank,i,nx,ny);
                                }
                            }
                        }
                        temp_real=realPart(FF);
                        temp_imag=imagPart(FF);
                        //------
                        for (int dd0=1; dd0<=FF.indices()[0].m(); dd0++) { //int dd0=1; dd0<=F[rank].indices()[0].m(); dd0++
                            for (int dd1=1; dd1<=FF.indices()[1].m(); dd1++) { //int dd1=1; dd1<=F[rank].indices()[1].m(); dd1++
                                MPI_Send(&(temp_real(FF.indices()[0](dd0),FF.indices()[1](dd1))), 1, MPI_DOUBLE, root, (X+Y*nx)*1000+rank*100+dd0*10+dd1, communicator);
                                MPI_Send(&(temp_imag(FF.indices()[0](dd0),FF.indices()[1](dd1))), 1, MPI_DOUBLE, root, 10000+(X+Y*nx)*1000+rank*100+dd0*10+dd1, communicator);
                            }
                        }
                        MPI_Barrier(communicator);
                    }
                    
                    
                    if (rank==root) {
                        
                        for (int rr=1; rr<size; rr++) {
                            temp_real=realPart(F[rr]);
                            temp_imag=imagPart(F[rr]);
                            for (int dd0=1; dd0<=F[rr].indices()[0].m(); dd0++) { //int dd0=1; dd0<=F[rank].indices()[0].m(); dd0++
                                for (int dd1=1; dd1<=F[rr].indices()[1].m(); dd1++) { //int dd1=1; dd1<=F[rank].indices()[1].m(); dd1++
                                    MPI_Recv(&(temp_real(F[rr].indices()[0](dd0),F[rr].indices()[1](dd1))), 1, MPI_DOUBLE, rr, (X+Y*nx)*1000+rr*100+dd0*10+dd1, communicator, &status);
                                    MPI_Recv(&(temp_imag(F[rr].indices()[0](dd0),F[rr].indices()[1](dd1))), 1, MPI_DOUBLE, rr, 10000+(X+Y*nx)*1000+rr*100+dd0*10+dd1, communicator, &status);
                                }
                            }
                            F[rr]=temp_real+Complex_i*temp_imag;
                        }
                        AF=Matrix(2*IF.indices()[0].m(), 2*IF.indices()[1].m());
                        BF=Matrix(2*IF.indices()[0].m(), 2*IF.indices()[1].m());
                        
                        evector_real=Vector(IF.indices()[0].m());
                        evector_imag=Vector(IF.indices()[0].m());
                        
                        for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,1)
                            for (int y=1; y<=IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=realPart(IF)(IF.indices()[0](x),IF.indices()[1](y));
                            }
                        }
                        for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,2)
                            for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=realPart(IF)(IF.indices()[0](x-IF.indices()[0].m()),IF.indices()[1](y-IF.indices()[1].m()));
                            }
                        }
                        for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,2)
                            for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=-imagPart(IF)(IF.indices()[0](x),IF.indices()[1](y-IF.indices()[1].m()));
                            }
                        }
                        for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,1)
                            for (int y=1; y<=IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=imagPart(IF)(IF.indices()[0](x-IF.indices()[0].m()),IF.indices()[1](y));
                            }
                        }
                        for (int rr=1; rr<size; rr++) {
                            for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,1)
                                for (int y=1; y<=IF.indices()[1].m(); y++) {
                                    AF(x,y)+=realPart(F[rr])(F[rr].indices()[0](x),F[rr].indices()[1](y));
                                }
                            }
                            for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,2)
                                for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                    AF(x,y)+=realPart(F[rr])(F[rr].indices()[0](x-IF.indices()[0].m()),F[rr].indices()[1](y-IF.indices()[1].m()));
                                }
                            }
                            for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,2)
                                for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                    AF(x,y)+=-imagPart(F[rr])(F[rr].indices()[0](x),F[rr].indices()[1](y-IF.indices()[1].m()));
                                }
                            }
                            for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,1)
                                for (int y=1; y<=IF.indices()[1].m(); y++) {
                                    AF(x,y)+=imagPart(F[rr])(F[rr].indices()[0](x-IF.indices()[0].m()),F[rr].indices()[1](y));
                                }
                            }
                        }
                        
                        
                        
                        
                        //                        const int Nrow = 3;
                        //                        const int Ncol = 3;
                        //                        const int maxm = min(Nrow,Ncol);
                        //
                        //                        Matrix M(Nrow,Ncol);
                        //                        M(1,1) = 0.1; M(1,2) = 0.223707; M(1,3) = 0.10;
                        //                        M(2,1) = 0.223707; M(2,2) = 0.213707; M(2,3) = -0.10;
                        //                        M(3,1) = 0.10; M(3,2) = -0.10; M(3,3) = 0.9;
                        //                        Print(M);
                        //
                        //                        Matrix N(Nrow,Ncol);
                        //                        N(1,1) = 2; N(1,2) = -1; N(1,3) = 0;
                        //                        N(2,1) = -1; N(2,2) = 2; N(2,3) = -1;
                        //                        N(3,1) = 0; N(3,2) = -1; N(3,3) = 2;
                        //                        Print(N);
                        //
                        //                        Vector rr;
                        //                        Vector ii;
                        //                        Vector beta;
                        //                        Matrix zzl;
                        //                        Matrix zzr;
                        //                        GeneralizedNonSymmetricEigenProblem(M, N, rr, ii, beta, zzl, zzr);
                        //                        Print(rr);
                        //                        Print(ii);
                        //                        Print(beta);
                        //                        Print(zzl);
                        //                        Print(zzr);
                        //
                        //                        double en;
                        //                        Vector evector(3);
                        //                        for (int x=1; x<=3; x++) {
                        //                            if (ii(x)==0) {
                        //                                en=rr(x)/beta(x);
                        //                                for (int y=1; y<=3; y++) {
                        //                                    evector(y)=zzr(y,x);
                        //                                }
                        //                                cout<<en<<endl;
                        //                                Print(evector);
                        //                            }
                        //                        }
                        
                        if (BF(1,1)!=0) {
                            for (int x=2*IF.indices()[0].m(); x>=1; x--) {
                                for (int y=2*IF.indices()[0].m(); y>=1; y--) {
                                    AF(x,y)=AF(x,y)/BF(1,1);
                                    BF(x,y)=BF(x,y)/BF(1,1);
                                }
                            }
                        }
                        else {
                            cout<<"BF(1,1)==0, AF and BF matrices are not normalized!!!"<<endl;
                        }
                        
                        Matrix BFF=BF;
                        CholeskyFactorization(BFF);
                        
                        if (cholesky_flag==0) {
                            EigenValues(BF, eigenvalues, eigenvectors);
                            Print(eigenvalues);
                            
                            for (int i=0; i<N; i++) {
                                cout<<"i="<<i;
                                PrintDat(a[i]);
                            }
                        }
                        
                        GeneralizedEigenProblem(AF, BF, eigenvalues, eigenvectors);
                        
                        min_index=1;
                        
                        //------begin: print format, controled by several parameters--
                        energy=eigenvalues(min_index);
                        if (endsweepflag==1) {
                            if (X==nx-1&&Y==ny-1&&xsweep==nxsweep) {
                                printf("Energy= %.10f\n",energy);
                            }
                        }
                        else {
                            printf("Energy= %.10f\n",energy);
                        }
                        
                        if (printflag==0) {
                            cout<<"--------------------------------------"<<endl;
                            cout<<endl;
                        }
                        //------end: print format, controled by several parameters--
                        
                        //------begin: update a tensors-----------------------------
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            evector_real(x)=eigenvectors(x,min_index);
                            evector_imag(x)=eigenvectors(x+IF.indices()[0].m(),min_index);
                        }
                        
                        realflag=0;
                        imagflag=0;
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            if (evector_real(x)!=0) {
                                realflag=1;
                                break;
                            }
                            //                            realflag+=evector_real(x);
                            //                            imagflag+=evector_imag(x);
                        }
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            if (evector_imag(x)!=0) {
                                imagflag=1;
                                break;
                            }
                        }
                        if (realflag==0&&imagflag==1) {
                            evector_temp=Vector(IF.indices()[0].m());
                            for (int x=1; x<=IF.indices()[0].m(); x++) {
                                evector_temp(x)=evector_real(x);
                            }
                            for (int x=1; x<=IF.indices()[0].m(); x++) {
                                evector_real(x)=evector_imag(x);
                                evector_imag(x)=evector_temp(x);
                            }
                        }
                        
                        a[Y*nx+X]=mapping_complex_vector_to_tensor(evector_real,evector_imag,a[Y*nx+X],X,Y,nx,ny);
                        
                        //                        cout<<"X="<<X<<", Y="<<Y;
                        //                        PrintDat(a[Y*nx+X]);
                        
                        //~~~~~~~~~~~~if (NORM==1) {
                        //a[Y*nx+X] /= a[Y*nx+X].norm();
                        //~~~~~~~~~~~~}
                        areal[Y*nx+X]=realPart(a[Y*nx+X]);
                        aimag[Y*nx+X]=imagPart(a[Y*nx+X]);
                        //------end: update a tensors-----------------------------
                        //---end: GE problem when AF,BF are symmetric and BF positive definite--
                        //}
                        
                        MPI_Barrier(communicator);
                    }
                    else {
                        
                        //MPI_Recv(&skip_flag, 1, MPI_INT, root, 0, communicator, &status);
                        
                        //if (skip_flag==0) {
                        areal[Y*nx+X]=realPart(a[Y*nx+X]);
                        aimag[Y*nx+X]=imagPart(a[Y*nx+X]);
                        //}
                        MPI_Barrier(communicator);
                        
                    }
                    
                    
                    //if (skip_flag==0) {
                    //broading the updated a to all processors
                    if (a[Y*nx+X].r()==3) {
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2))), 1, MPI_DOUBLE, 0, communicator);
                                    MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2))), 1, MPI_DOUBLE, 0, communicator);
                                }
                            }
                        }
                    }
                    else if (a[Y*nx+X].r()==4){
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    for (int i3=1; i3<=a[Y*nx+X].indices()[3].m(); i3++) {
                                        MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2),areal[Y*nx+X].indices()[3](i3))), 1, MPI_DOUBLE, 0, communicator);
                                        MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2),aimag[Y*nx+X].indices()[3](i3))), 1, MPI_DOUBLE, 0, communicator);
                                    }
                                }
                            }
                        }
                    }
                    else if (a[Y*nx+X].r()==5){
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    for (int i3=1; i3<=a[Y*nx+X].indices()[3].m(); i3++) {
                                        for (int i4=1; i4<=a[Y*nx+X].indices()[4].m(); i4++) {
                                            MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2),areal[Y*nx+X].indices()[3](i3),areal[Y*nx+X].indices()[4](i4))), 1, MPI_DOUBLE, 0, communicator);
                                            MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2),aimag[Y*nx+X].indices()[3](i3),aimag[Y*nx+X].indices()[4](i4))), 1, MPI_DOUBLE, 0, communicator);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    a[Y*nx+X]=areal[Y*nx+X]+Complex_i*aimag[Y*nx+X];
                    
                    int ind=Y*nx+X;
                    
                    if (rank==root) {
                        //                        cout<<"IW before updating at ("<<X<<","<<Y<<")";
                        //                        Print(IW[ind]);
                        if (Y==0) {
                            if (X==0) { //X==0, Y==0
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else if (X==nx-1){ //X==nx-1, Y==0
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else { //0<X<nx-1, Y==0
                                //cout<<"409~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"411~~~~~~~~~~"<<endl;
                            }
                        }
                        else if (Y==ny-1){
                            if (X==0) { //X==0, Y==ny-1
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else if (X==nx-1){ //X==nx-1, Y==ny-1
                                //cout<<"419~~~~~~~"<<endl;
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                                //cout<<"421~~~~~~~"<<endl;
                            }
                            else { //0<X<nx-1, Y==ny-1
                                //cout<<"422~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"424~~~~~~~~~~"<<endl;
                            }
                        }
                        else {
                            if (X==0) { //X==0, 0<Y<ny-1
                                //cout<<"429~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"431~~~~~~~~~~"<<endl;
                            }
                            else if (X==nx-1){ //X==nx-1, 0<Y<ny-1
                                //cout<<"434~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"436~~~~~~~~~~"<<endl;
                            }
                            else { //0<X<nx-1, 0<Y<ny-1
                                IW[ind]=group_indice8((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2], a[ind].indices()[4], primed(a[ind].indices()[4]), IW[ind].indices()[3]);
                            }
                        }
                        MPI_Barrier(communicator);
                    }
                    
                    else {  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        for (int i=mystart; i<myend; i++) {
                            
                            if (Y==0) {
                                if (X==0) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else if (X==nx-1) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else {
                                    //cout<<"456~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"458~~~~~~~~~~"<<endl;
                                }
                            }
                            else if (Y==ny-1) {
                                if (X==0) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else if (X==nx-1) {
                                    //cout<<"468~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                    //cout<<"470~~~~~~~"<<endl;
                                }
                                else {
                                    //cout<<"469~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"471~~~~~~~~~~"<<endl;
                                }
                            }
                            else {
                                if (X==0) {
                                    //cout<<"476~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"478~~~~~~~~~~"<<endl;
                                }
                                else if (X==nx-1) {
                                    //cout<<"481~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"483~~~~~~~~~~"<<endl;
                                }
                                else {
                                    W[ind+N*i]=group_indice8((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2], a[ind].indices()[4], primed(a[ind].indices()[4]), W[ind+N*i].indices()[3]);
                                }
                            }
                            
                        }
                        MPI_Barrier(communicator);
                    }
                    //}
                    
                }
            }
            //----end: sweep of a row------------------------
            
            //}//~~~~~~
        }
        //----end: row index looping, from Y=0 to Y=ny-1------
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        //----begin: row index looping backward, from Y=ny-1 to Y=0------
        for (int Y=ny-1; Y>=0; Y--) {    //int Y=0; Y<ny; Y++ //int Y=1; Y<ny-1; Y++
            
            
            if (Y==ny-1) {
                xmin=0;
                xmax=nx-2;
            }
            else if(Y==0){
                xmin=1;
                xmax=nx-1;
            }
            else {
                xmin=0;
                xmax=nx-1;
            }
            
            
            if (rank==root && printflag==0) {
                cout<<"Sweep "<<sweep_ind<<", Y="<<Y<<" (backward)"<<endl;
            }
            //----begin: build up _bra and _ket for the loop------------------------------------
            if (rank==root) {
                if (Y>0) {
                    for (int x=0; x<nx; x++) {
                        if (_backflag==0) {
                            _Iket[x]=Iket[x+nx*(Y-1)];
                        }
                        else {
                            b_Iket[x]=bIket[x+nx*(Y-1)];
                        }
                    }
                }
                if (Y<ny-1) {
                    if (Y==ny-2) {
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                IA[x]=IW[x+nx*(Y+1)];
                            }
                            else {
                                bIA[x]=IW[x+nx*(Y+1)];
                            }
                        }
                        if (_backflag==0) {
                            IA=bra_otho_center_1(IA,Y+1,nx,ny);
                        }
                        else {
                            bIA=bra_otho_center_1(bIA,Y+1,nx,ny);
                        }
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                Ibra[x]=IA[x];
                                _Ibra[x]=IA[x];
                            }
                            else {
                                bIbra[x]=bIA[x];
                                b_Ibra[x]=bIA[x];
                            }
                        }
                    }
                    else {
                        std::vector<ITensor> bra_temp(nx),gate_temp(nx),temp2(nx);
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                bra_temp[x]=Ibra[x+nx*(ny-3-Y)];
                            }
                            else {
                                bra_temp[x]=bIbra[x+nx*(ny-3-Y)];
                            }
                            gate_temp[x]=IW[x+nx*(Y+1)];
                        }
                        temp2=compressor2(bra_temp,gate_temp,ny-2-Y,cutoff,m,dcut,nx,ny);
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                Ibra[x+nx*(ny-2-Y)]=temp2[x];
                                _Ibra[x]=temp2[x];
                            }
                            else {
                                bIbra[x+nx*(ny-2-Y)]=temp2[x];
                                b_Ibra[x]=temp2[x];
                            }
                        }
                        
                    }
                }
                
                MPI_Barrier(communicator);
            }
            else {
                for (int i=mystart; i<myend; i++) {
                    if (Y>0) {
                        for (int x=0; x<nx; x++) {
                            if (_backflag==0) {
                                _ket[x+nx*i]=ket[x+nx*(Y-1)+nx*(ny-1)*i];
                            }
                            else {
                                b_ket[x+nx*i]=bket[x+nx*(Y-1)+nx*(ny-1)*i];
                            }
                        }
                    }
                    if (Y<ny-1) {
                        if (Y==ny-2) {
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    AA[x]=W[x+nx*(Y+1)+N*i];
                                }
                                else {
                                    bAA[x]=W[x+nx*(Y+1)+N*i];
                                }
                            }
                            if (_backflag==0) {
                                AA=bra_otho_center_1(AA,Y+1,nx,ny);
                            }
                            else {
                                bAA=bra_otho_center_1(bAA,Y+1,nx,ny);
                            }
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    bra[x+nx*(ny-1)*i]=AA[x];
                                    _bra[x+nx*i]=AA[x];
                                }
                                else {
                                    bbra[x+nx*(ny-1)*i]=bAA[x];
                                    b_bra[x+nx*i]=bAA[x];
                                }
                            }
                        }
                        else {
                            std::vector<ITensor> bra_temp(nx),gate_temp(nx),temp2(nx);
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    bra_temp[x]=bra[x+nx*(ny-3-Y)+nx*(ny-1)*i];
                                }
                                else {
                                    bra_temp[x]=bbra[x+nx*(ny-3-Y)+nx*(ny-1)*i];
                                }
                                gate_temp[x]=W[x+nx*(Y+1)+N*i];
                            }
                            temp2=compressor2(bra_temp,gate_temp,ny-2-Y,cutoff,m,dcut,nx,ny);
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    bra[x+nx*(ny-2-Y)+nx*(ny-1)*i]=temp2[x];
                                    _bra[x+nx*i]=temp2[x];
                                }
                                else {
                                    bbra[x+nx*(ny-2-Y)+nx*(ny-1)*i]=temp2[x];
                                    b_bra[x+nx*i]=temp2[x];
                                }
                            }
                            
                        }
                    }
                }
                
                MPI_Barrier(communicator);
            }
            //----end: build up _bra and _ket for the loop------------------------------------
            
            //if (Y==2) {//~~~~~
            
            //----begin: sweep of a row------------------------
            for (int xsweep=1; xsweep<=nxsweep; xsweep++) {
                for (int X=xmax; X>=xmin; X--) {  //int X=0; X<nx; X++
                    
                    //skip_flag=0;
                    
                    std::vector<ITensor> F(size);
                    if (rank==root) {
                        if (printflag==0) {
                            cout<<"--------------------------------------"<<endl;
                            cout<<"Working on ("<<X<<","<<Y<<") (backward):"<<endl;
                            //PrintDat(a[0]);//~~~~~~~
                        }
                        if (_backflag==0) {
                            IF=sweepx(_Ibra,_Iket,X,Y,rank,0,nx,ny);
                        }
                        else {
                            IF=sweepx(b_Ibra,b_Iket,X,Y,rank,0,nx,ny);
                        }
                        for (int xxx=0; xxx<size; xxx++) {
                            F[xxx]=IF;
                        }
                        MPI_Barrier(communicator);
                    }
                    else {
                        for (int i=mystart; i<myend; i++) {
                            for (int x=0; x<nx; x++) {
                                if (_backflag==0) {
                                    bra_temp[x]=_bra[x+nx*i];
                                    ket_temp[x]=_ket[x+nx*i];
                                }
                                else {
                                    bbra_temp[x]=b_bra[x+nx*i];
                                    bket_temp[x]=b_ket[x+nx*i];
                                }
                            }
                            if (i==mystart) {
                                if (_backflag==0) {
                                    FF=sweepx(bra_temp,ket_temp,X,Y,rank,i,nx,ny);
                                }
                                else {
                                    FF=sweepx(bbra_temp,bket_temp,X,Y,rank,i,nx,ny);
                                }
                            }
                            else {
                                if (_backflag==0) {
                                    FF+=sweepx(bra_temp,ket_temp,X,Y,rank,i,nx,ny);
                                }
                                else {
                                    FF+=sweepx(bbra_temp,bket_temp,X,Y,rank,i,nx,ny);
                                }
                            }
                        }
                        temp_real=realPart(FF);
                        temp_imag=imagPart(FF);
                        //------
                        for (int dd0=1; dd0<=FF.indices()[0].m(); dd0++) { //int dd0=1; dd0<=F[rank].indices()[0].m(); dd0++
                            for (int dd1=1; dd1<=FF.indices()[1].m(); dd1++) { //int dd1=1; dd1<=F[rank].indices()[1].m(); dd1++
                                MPI_Send(&(temp_real(FF.indices()[0](dd0),FF.indices()[1](dd1))), 1, MPI_DOUBLE, root, (X+Y*nx)*1000+rank*100+dd0*10+dd1, communicator);
                                MPI_Send(&(temp_imag(FF.indices()[0](dd0),FF.indices()[1](dd1))), 1, MPI_DOUBLE, root, 10000+(X+Y*nx)*1000+rank*100+dd0*10+dd1, communicator);
                            }
                        }
                        MPI_Barrier(communicator);
                    }
                    
                    
                    if (rank==root) {
                        
                        for (int rr=1; rr<size; rr++) {
                            temp_real=realPart(F[rr]);
                            temp_imag=imagPart(F[rr]);
                            for (int dd0=1; dd0<=F[rr].indices()[0].m(); dd0++) { //int dd0=1; dd0<=F[rank].indices()[0].m(); dd0++
                                for (int dd1=1; dd1<=F[rr].indices()[1].m(); dd1++) { //int dd1=1; dd1<=F[rank].indices()[1].m(); dd1++
                                    MPI_Recv(&(temp_real(F[rr].indices()[0](dd0),F[rr].indices()[1](dd1))), 1, MPI_DOUBLE, rr, (X+Y*nx)*1000+rr*100+dd0*10+dd1, communicator, &status);
                                    MPI_Recv(&(temp_imag(F[rr].indices()[0](dd0),F[rr].indices()[1](dd1))), 1, MPI_DOUBLE, rr, 10000+(X+Y*nx)*1000+rr*100+dd0*10+dd1, communicator, &status);
                                }
                            }
                            F[rr]=temp_real+Complex_i*temp_imag;
                        }
                        AF=Matrix(2*IF.indices()[0].m(), 2*IF.indices()[1].m());
                        BF=Matrix(2*IF.indices()[0].m(), 2*IF.indices()[1].m());
                        
                        evector_real=Vector(IF.indices()[0].m());
                        evector_imag=Vector(IF.indices()[0].m());
                        for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,1)
                            for (int y=1; y<=IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=realPart(IF)(IF.indices()[0](x),IF.indices()[1](y));
                            }
                        }
                        for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,2)
                            for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=realPart(IF)(IF.indices()[0](x-IF.indices()[0].m()),IF.indices()[1](y-IF.indices()[1].m()));
                            }
                        }
                        for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,2)
                            for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=-imagPart(IF)(IF.indices()[0](x),IF.indices()[1](y-IF.indices()[1].m()));
                            }
                        }
                        for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,1)
                            for (int y=1; y<=IF.indices()[1].m(); y++) {
                                AF(x,y)=0;
                                BF(x,y)=imagPart(IF)(IF.indices()[0](x-IF.indices()[0].m()),IF.indices()[1](y));
                            }
                        }
                        for (int rr=1; rr<size; rr++) {
                            for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,1)
                                for (int y=1; y<=IF.indices()[1].m(); y++) {
                                    AF(x,y)+=realPart(F[rr])(F[rr].indices()[0](x),F[rr].indices()[1](y));
                                }
                            }
                            for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,2)
                                for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                    AF(x,y)+=realPart(F[rr])(F[rr].indices()[0](x-IF.indices()[0].m()),F[rr].indices()[1](y-IF.indices()[1].m()));
                                }
                            }
                            for (int x=1; x<=IF.indices()[0].m(); x++) {    //block: (1,2)
                                for (int y=1+IF.indices()[1].m(); y<=2*IF.indices()[1].m(); y++) {
                                    AF(x,y)+=-imagPart(F[rr])(F[rr].indices()[0](x),F[rr].indices()[1](y-IF.indices()[1].m()));
                                }
                            }
                            for (int x=1+IF.indices()[0].m(); x<=2*IF.indices()[0].m(); x++) {    //block: (2,1)
                                for (int y=1; y<=IF.indices()[1].m(); y++) {
                                    AF(x,y)+=imagPart(F[rr])(F[rr].indices()[0](x-IF.indices()[0].m()),F[rr].indices()[1](y));
                                }
                            }
                        }
                        
                        if (BF(1,1)!=0) {
                            for (int x=2*IF.indices()[0].m(); x>=1; x--) {
                                for (int y=2*IF.indices()[0].m(); y>=1; y--) {
                                    AF(x,y)=AF(x,y)/BF(1,1);
                                    BF(x,y)=BF(x,y)/BF(1,1);
                                }
                            }
                        }
                        else {
                            cout<<"BF(1,1)==0, AF and BF matrices are not normalized!!!"<<endl;
                        }
                        
                        Matrix BFF=BF;
                        CholeskyFactorization(BFF);
                        
                        if (cholesky_flag==0) {
                            EigenValues(BF, eigenvalues, eigenvectors);
                            Print(eigenvalues);
                            
                            for (int i=0; i<N; i++) {
                                cout<<"i="<<i;
                                PrintDat(a[i]);
                            }
                        }
                        
                        GeneralizedEigenProblem(AF, BF, eigenvalues, eigenvectors);
                        
                        min_index=1;
                        
                        //------begin: print format, controled by several parameters--
                        energy=eigenvalues(min_index);
                        if (endsweepflag==1) {
                            if (X==nx-1&&Y==ny-1&&xsweep==nxsweep) {
                                printf("Energy= %.10f\n",energy);
                            }
                        }
                        else {
                            printf("Energy= %.10f\n",energy);
                        }
                        
                        if (printflag==0) {
                            cout<<"--------------------------------------"<<endl;
                            cout<<endl;
                        }
                        //------end: print format, controled by several parameters--
                        
                        //------begin: update a tensors-----------------------------
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            evector_real(x)=eigenvectors(x,min_index);
                            evector_imag(x)=eigenvectors(x+IF.indices()[0].m(),min_index);
                        }
                        
                        realflag=0;
                        imagflag=0;
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            if (evector_real(x)!=0) {
                                realflag=1;
                                break;
                            }
                            //                            realflag+=evector_real(x);
                            //                            imagflag+=evector_imag(x);
                        }
                        for (int x=1; x<=IF.indices()[0].m(); x++) {
                            if (evector_imag(x)!=0) {
                                imagflag=1;
                                break;
                            }
                        }
                        if (realflag==0&&imagflag==1) {
                            evector_temp=Vector(IF.indices()[0].m());
                            for (int x=1; x<=IF.indices()[0].m(); x++) {
                                evector_temp(x)=evector_real(x);
                            }
                            for (int x=1; x<=IF.indices()[0].m(); x++) {
                                evector_real(x)=evector_imag(x);
                                evector_imag(x)=evector_temp(x);
                            }
                        }
                
                        
                        
                        a[Y*nx+X]=mapping_complex_vector_to_tensor(evector_real,evector_imag,a[Y*nx+X],X,Y,nx,ny);
                        
                        //                        cout<<"X="<<X<<", Y="<<Y;
                        //                        PrintDat(a[Y*nx+X]);
                        
                        //~~~~~~~~~~~~if (NORM==1) {
                        //a[Y*nx+X] /= a[Y*nx+X].norm();
                        //~~~~~~~~~~~~}
                        areal[Y*nx+X]=realPart(a[Y*nx+X]);
                        aimag[Y*nx+X]=imagPart(a[Y*nx+X]);
                        //------end: update a tensors-----------------------------
                        //---end: GE problem when AF,BF are symmetric and BF positive definite--
                        //}
                        
                        MPI_Barrier(communicator);
                        
                    }
                    else {
                        //MPI_Recv(&skip_flag, 1, MPI_INT, root, 0, communicator, &status);
                        
                        //if (skip_flag==0) {//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        areal[Y*nx+X]=realPart(a[Y*nx+X]);
                        aimag[Y*nx+X]=imagPart(a[Y*nx+X]);
                        //}
                        MPI_Barrier(communicator);
                    }
                    
                    
                    //if (skip_flag==0) {///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    //broading the updated a to all processors
                    if (a[Y*nx+X].r()==3) {
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2))), 1, MPI_DOUBLE, 0, communicator);
                                    MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2))), 1, MPI_DOUBLE, 0, communicator);
                                }
                            }
                        }
                    }
                    else if (a[Y*nx+X].r()==4){
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    for (int i3=1; i3<=a[Y*nx+X].indices()[3].m(); i3++) {
                                        MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2),areal[Y*nx+X].indices()[3](i3))), 1, MPI_DOUBLE, 0, communicator);
                                        MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2),aimag[Y*nx+X].indices()[3](i3))), 1, MPI_DOUBLE, 0, communicator);
                                    }
                                }
                            }
                        }
                    }
                    else if (a[Y*nx+X].r()==5){
                        for (int i0=1; i0<=a[Y*nx+X].indices()[0].m(); i0++) {
                            for (int i1=1; i1<=a[Y*nx+X].indices()[1].m(); i1++) {
                                for (int i2=1; i2<=a[Y*nx+X].indices()[2].m(); i2++) {
                                    for (int i3=1; i3<=a[Y*nx+X].indices()[3].m(); i3++) {
                                        for (int i4=1; i4<=a[Y*nx+X].indices()[4].m(); i4++) {
                                            MPI_Bcast(&(areal[Y*nx+X](areal[Y*nx+X].indices()[0](i0),areal[Y*nx+X].indices()[1](i1),areal[Y*nx+X].indices()[2](i2),areal[Y*nx+X].indices()[3](i3),areal[Y*nx+X].indices()[4](i4))), 1, MPI_DOUBLE, 0, communicator);
                                            MPI_Bcast(&(aimag[Y*nx+X](aimag[Y*nx+X].indices()[0](i0),aimag[Y*nx+X].indices()[1](i1),aimag[Y*nx+X].indices()[2](i2),aimag[Y*nx+X].indices()[3](i3),aimag[Y*nx+X].indices()[4](i4))), 1, MPI_DOUBLE, 0, communicator);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    a[Y*nx+X]=areal[Y*nx+X]+Complex_i*aimag[Y*nx+X];
                    
                    int ind=Y*nx+X;
                    
                    if (rank==root) {
                        //                        cout<<"IW before updating at ("<<X<<","<<Y<<")";
                        //                        Print(IW[ind]);
                        if (Y==0) {
                            if (X==0) { //X==0, Y==0
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else if (X==nx-1){ //X==nx-1, Y==0
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else { //0<X<nx-1, Y==0
                                //cout<<"409~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"411~~~~~~~~~~"<<endl;
                            }
                        }
                        else if (Y==ny-1){
                            if (X==0) { //X==0, Y==ny-1
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                            }
                            else if (X==nx-1){ //X==nx-1, Y==ny-1
                                //cout<<"419~~~~~~~"<<endl;
                                IW[ind]=group_indice4((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1]);
                                //cout<<"421~~~~~~~"<<endl;
                            }
                            else { //0<X<nx-1, Y==ny-1
                                //cout<<"422~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"424~~~~~~~~~~"<<endl;
                            }
                        }
                        else {
                            if (X==0) { //X==0, 0<Y<ny-1
                                //cout<<"429~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"431~~~~~~~~~~"<<endl;
                            }
                            else if (X==nx-1){ //X==nx-1, 0<Y<ny-1
                                //cout<<"434~~~~~~~~~~"<<endl;
                                IW[ind]=group_indice6((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2]);
                                //cout<<"436~~~~~~~~~~"<<endl;
                            }
                            else { //0<X<nx-1, 0<Y<ny-1
                                IW[ind]=group_indice8((Identity[ind]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), IW[ind].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), IW[ind].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), IW[ind].indices()[2], a[ind].indices()[4], primed(a[ind].indices()[4]), IW[ind].indices()[3]);
                            }
                        }
                        MPI_Barrier(communicator);
                    }
                    
                    else {  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        for (int i=mystart; i<myend; i++) {
                            
                            if (Y==0) {
                                if (X==0) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else if (X==nx-1) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else {
                                    //cout<<"456~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"458~~~~~~~~~~"<<endl;
                                }
                            }
                            else if (Y==ny-1) {
                                if (X==0) {
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                }
                                else if (X==nx-1) {
                                    //cout<<"468~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice4((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1]);
                                    //cout<<"470~~~~~~~"<<endl;
                                }
                                else {
                                    //cout<<"469~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"471~~~~~~~~~~"<<endl;
                                }
                            }
                            else {
                                if (X==0) {
                                    //cout<<"476~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"478~~~~~~~~~~"<<endl;
                                }
                                else if (X==nx-1) {
                                    //cout<<"481~~~~~~~~~~"<<endl;
                                    W[ind+N*i]=group_indice6((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2]);
                                    //cout<<"483~~~~~~~~~~"<<endl;
                                }
                                else {
                                    W[ind+N*i]=group_indice8((O[ind+N*i]*a[ind]*dag(prime(a[ind]))), a[ind].indices()[1], primed(a[ind].indices()[1]), W[ind+N*i].indices()[0], a[ind].indices()[2], primed(a[ind].indices()[2]), W[ind+N*i].indices()[1], a[ind].indices()[3], primed(a[ind].indices()[3]), W[ind+N*i].indices()[2], a[ind].indices()[4], primed(a[ind].indices()[4]), W[ind+N*i].indices()[3]);
                                }
                            }
                            
                        }
                        MPI_Barrier(communicator);
                    }
                    //}
                }
            }
            //----end: sweep of a row------------------------
            
            //}//~~~~~~
        }
        //----end: row index looping backward, from Y=ny-1 to Y=0------
        
        
        
        //                if (rank==root) {
        //                    identity_find_bras(IW,nx,ny);  //find all Ibra's
        //                    MPI_Barrier(communicator);
        //                }
        //                else {
        //                    std::vector<ITensor> WW(nx*ny);
        //                    for (int i=mystart; i<myend; i++) {
        //                        for (int x=0; x<nx*ny; x++) {
        //                            WW[x]=W[x+nx*ny*i];
        //                        }
        //                        ham_find_bras(WW, i, nx,ny);   //find all bra's
        //                    }
        //
        //                    MPI_Barrier(communicator);
        //                }
    }
    
    
    void define_indices(int nx, int ny, MPI_Comm communicator){
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
        
        MPI_Barrier(communicator);
    }
    
    void backward_transform(int root, int nx, int ny, MPI_Comm communicator){
        int rank;
        MPI_Comm_rank(communicator, &rank);
        int size;
        MPI_Comm_size(communicator, &size);
        MPI_Status status;
        
        int x_new,y_new;
        
        //        for (int i=0; i<N; i++) {
        //            s[i]=Index(nameint("s",i),NS,Site);
        //            u[i]=Index(nameint("u",i),NB);
        //            d[i]=Index(nameint("d",i),NB);
        //            l[i]=Index(nameint("l",i),NB);
        //            r[i]=Index(nameint("r",i),NB);
        //        }
        //        for (int i=0; i<N; i++) {
        //            if ((i+1)%ny!=0) {
        //                r[i]=l[i+1];
        //            }
        //            if (i/ny!=nx-1) {
        //                d[i]=u[i+ny];
        //            }
        //        }
        
        
        //if (rank==0) {
        std::vector<ITensor> a_real(N), a_imag(N);
        //            for (int k=0; k<N; k++) {
        //                a_real[k]=realPart(a[k]);
        //                a_imag[k]=imagPart(a[k]);
        //            }
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=x+nx*y;
                x_new=ny-1-y;
                y_new=nx-1-x;
                int newind=x_new+ny*y_new;
                
                if (y_new==0) {
                    if (x_new==0) {
                        a_real[newind]=ITensor(s[newind],d[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],d[newind],r[newind]);
                    }
                    else if (x_new==ny-1){
                        a_real[newind]=ITensor(s[newind],d[newind],l[newind]);
                        a_imag[newind]=ITensor(s[newind],d[newind],l[newind]);
                    }
                    else {
                        a_real[newind]=ITensor(s[newind],d[newind],l[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],d[newind],l[newind],r[newind]);
                    }
                }
                else if (y_new==nx-1){
                    if (x_new==0) {
                        a_real[newind]=ITensor(s[newind],u[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],r[newind]);
                    }
                    else if (x_new==ny-1){
                        a_real[newind]=ITensor(s[newind],u[newind],l[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],l[newind]);
                    }
                    else {
                        a_real[newind]=ITensor(s[newind],u[newind],l[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],l[newind],r[newind]);
                    }
                }
                else {
                    if (x_new==0) {
                        a_real[newind]=ITensor(s[newind],u[newind],d[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],d[newind],r[newind]);
                    }
                    else if (x_new==ny-1){
                        a_real[newind]=ITensor(s[newind],u[newind],d[newind],l[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],d[newind],l[newind]);
                    }
                    else {
                        a_real[newind]=ITensor(s[newind],u[newind],d[newind],l[newind],r[newind]);
                        a_imag[newind]=ITensor(s[newind],u[newind],d[newind],l[newind],r[newind]);
                    }
                }
                if (a[ind].r()==3) {
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                a_real[newind](a_real[newind].indices()[0](i0),a_real[newind].indices()[1](i2),a_real[newind].indices()[2](i1))=realPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2));
                                a_imag[newind](a_imag[newind].indices()[0](i0),a_imag[newind].indices()[1](i2),a_imag[newind].indices()[2](i1))=imagPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2));
                            }
                        }
                    }
                }
                else if (a[ind].r()==4){
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                    a_real[newind](a_real[newind].indices()[0](i0),a_real[newind].indices()[1](i3),a_real[newind].indices()[2](i2),a_real[newind].indices()[3](i1))=realPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3));
                                    a_imag[newind](a_imag[newind].indices()[0](i0),a_imag[newind].indices()[1](i3),a_imag[newind].indices()[2](i2),a_imag[newind].indices()[3](i1))=imagPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3));
                                }
                            }
                        }
                    }
                }
                else if (a[ind].r()==5){
                    for (int i0=1; i0<=a[ind].indices()[0].m(); i0++) {
                        for (int i1=1; i1<=a[ind].indices()[1].m(); i1++) {
                            for (int i2=1; i2<=a[ind].indices()[2].m(); i2++) {
                                for (int i3=1; i3<=a[ind].indices()[3].m(); i3++) {
                                    for (int i4=1; i4<=a[ind].indices()[4].m(); i4++) {
                                        a_real[newind](a_real[newind].indices()[0](i0),a_real[newind].indices()[1](i4),a_real[newind].indices()[2](i3),a_real[newind].indices()[3](i2),a_real[newind].indices()[4](i1))=realPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3),a[ind].indices()[4](i4));
                                        a_imag[newind](a_imag[newind].indices()[0](i0),a_imag[newind].indices()[1](i4),a_imag[newind].indices()[2](i3),a_imag[newind].indices()[3](i2),a_imag[newind].indices()[4](i1))=imagPart(a[ind])(a[ind].indices()[0](i0),a[ind].indices()[1](i1),a[ind].indices()[2](i2),a[ind].indices()[3](i3),a[ind].indices()[4](i4));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int ind=x+nx*y;
                x_new=ny-1-y;
                y_new=nx-1-x;
                int newind=x_new+ny*y_new;
                
                //cout<<"ind="<<ind<<", newind="<<newind;
                //PrintDat(a[ind]);
                a[newind]=a_real[newind]+Complex_i*a_imag[newind];
                //PrintDat(a[newind]);
            }
        }
        
        MPI_Barrier(communicator);
    }
}
#endif





















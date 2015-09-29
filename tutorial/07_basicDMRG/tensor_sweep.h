//tensor_sweep.h
//include tensor_sweep

#ifndef ___TENSOR_SWEEP
#define ___TENSOR_SWEEP


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "tensor_utility.h"
#include "sweep_utility.h"
#include "parameters.h"
//#include "fitting.h"


using namespace std;
namespace itensor {
    //tensor_sweep(A,B,W,bra(ny)):
    void identity_find_bras(std::vector<ITensor> _W, int nx, int ny){
        

        //-----------begin: find all the fitted bra's and store them----------------------
        std::vector<ITensor> _B(nx);
        for (int x=0; x<nx; x++) {
            _B[x]=_W[x+nx*(ny-1)];
        }
        _B=bra_otho_center_1(_B, ny-1, nx, ny);//change the othogonality center of B to site 1
        //PrintDat(B[0]);
        for (int x=0; x<nx; x++) {
            if (backflag==0) {
                Ibra[x]=_B[x];
            }
            else {
                bIbra[x]=_B[x];
            }
            
        }
        
        
        std::vector<ITensor> temp1(nx),temp2(nx),gate(nx);
        for (int i=1; i<=ny-2; i++) {   //////////int i=1; i<=ny-2; i++
            for (int x=0; x<nx; x++) {
                if (backflag==0) {
                    temp1[x]=Ibra[x+nx*(i-1)];
                }
                else {
                    temp1[x]=bIbra[x+nx*(i-1)];
                }
                gate[x]=_W[x+nx*(ny-1-i)];
                
                //Print(W[x+nx*(ny-2-i)+nx*(ny-2)*hamx]);
                //cout<<x+nx*(ny-2-i)+nx*(ny-2)*hamx<<endl;
            }
            
            
            
            temp2=compressor2(temp1,gate,i,cutoff, m, dcut, nx, ny);  //to make variable m effective, make change in compressor in sweep_utility.h
            
            for (int x=0; x<nx; x++) {
                if (backflag==0) {
                    Ibra[nx*i+x]=temp2[x];
                }
                else {
                    bIbra[nx*i+x]=temp2[x];
                }
                
            }
            
            
            //---for test only------
            //if (i==1) {
            //    for (int xx=0; xx<nx; xx++) {
            //        cout<<xx<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            //        PrintDat(temp2[xx]);
            //    }
            //}
            //---for test only------
        }
        
        
        //-----------end: find all the fitted bra's and store them----------------------
        
        
        
        //bra_otho_center_1(temp2, ny-2);
        
        //for (int x=0; x<nx; x++) {
        //    cout<<x<<"temp1~~~~~~~~~"<<endl;
        //    PrintDat(temp2[x]);
        //cout<<x<<"gate~~~~~~~~~"<<endl;
        //PrintDat(gate[x]);
        //}
        
        
        
//         for (int i=0; i<=ny-2; i++) {
//         cout<<"UU"<<i<<"~~~~~~~~~"<<endl;
//         for (int x=0; x<nx; x++) {
//         cout<<x<<"~~~~~~~~~~~"<<endl;
//         PrintDat(Ibra[nx*i+x]);
//         }
//         }
        
    }
    
    
    void ham_find_bras(std::vector<ITensor> _W, int hamx, int nx, int ny){
        
            
        //-----------begin: find all the fitted bra's and store them----------------------
        std::vector<ITensor> _B(nx);
        for (int x=0; x<nx; x++) {
            _B[x]=_W[x+nx*(ny-1)];
            
//            if (backflag==1) {
//                cout<<"x="<<x;
//                PrintDat(_B[x]);
//            }
        }
        
        //if (backflag==0) {//~~~~~~
            
        _B=bra_otho_center_1(_B, ny-1, nx, ny);//change the othogonality center of B to site 1
        //if (backflag==0) {//~~~~~
        for (int x=0; x<nx; x++) {
            if (backflag==0) {
                bra[x+nx*(ny-1)*hamx]=_B[x];
            }
            else {
                bbra[x+nx*(ny-1)*hamx]=_B[x];
            }
        }
        
        std::vector<ITensor> temp1(nx),temp2(nx),gate(nx);
        for (int i=1; i<=ny-2; i++) {   //////////int i=1; i<=ny-2; i++
            for (int x=0; x<nx; x++) {
                if (backflag==0) {
                    temp1[x]=bra[x+nx*(i-1)+nx*(ny-1)*hamx];
                }
                else {
                    temp1[x]=bbra[x+nx*(i-1)+nx*(ny-1)*hamx];
                }
                
                gate[x]=_W[x+nx*(ny-1-i)];
                //Print(W[x+nx*(ny-2-i)+nx*(ny-2)*hamx]);
                //cout<<x+nx*(ny-2-i)+nx*(ny-2)*hamx<<endl;
            }
            
            temp2=compressor2(temp1,gate,i,cutoff, m, dcut, nx, ny);  //to make variable m effective, make change in compressor in sweep_utility.h
            
            for (int x=0; x<nx; x++) {
                if (backflag==0) {
                    bra[nx*i+x+nx*(ny-1)*hamx]=temp2[x];
                }
                else {
                    bbra[nx*i+x+nx*(ny-1)*hamx]=temp2[x];
                }
                
            }
            
            //---for test only------
            //if (i==1) {
            //    for (int xx=0; xx<nx; xx++) {
            //        cout<<xx<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            //        PrintDat(temp2[xx]);
            //    }
            //}
            //---for test only------
        }
            
        //}//~~~~~~~
        //-----------end: find all the fitted bra's and store them----------------------
        
        
    }
    
    
    
/*
    void bham_find_bras(std::vector<ITensor> _W, int hamx, int nx, int ny){
        
        
        //-----------begin: find all the fitted bra's and store them----------------------
        std::vector<ITensor> _B(nx);
        for (int x=0; x<nx; x++) {
            _B[x]=_W[x+nx*(ny-1)];
        }
        
        
        _B=bbra_otho_center_1(_B, ny-1, nx, ny);//change the othogonality center of B to site 1
        
//        for (int x=0; x<nx; x++) {
//            if (backflag==0) {
//                bra[x+nx*(ny-1)*hamx]=_B[x];
//            }
//            else {
//                bbra[x+nx*(ny-1)*hamx]=_B[x];
//            }
//        }
//        
//        std::vector<ITensor> temp1(nx),temp2(nx),gate(nx);
//        for (int i=1; i<=ny-2; i++) {   //////////int i=1; i<=ny-2; i++
//            for (int x=0; x<nx; x++) {
//                if (backflag==0) {
//                    temp1[x]=bra[x+nx*(i-1)+nx*(ny-1)*hamx];
//                }
//                else {
//                    temp1[x]=bbra[x+nx*(i-1)+nx*(ny-1)*hamx];
//                }
//                
//                gate[x]=_W[x+nx*(ny-1-i)];
//            }
//            
//            temp2=compressor2(temp1,gate,i,cutoff, m, dcut, nx, ny);  //to make variable m effective, make change in compressor in sweep_utility.h
//            
//            for (int x=0; x<nx; x++) {
//                if (backflag==0) {
//                    bra[nx*i+x+nx*(ny-1)*hamx]=temp2[x];
//                }
//                else {
//                    bbra[nx*i+x+nx*(ny-1)*hamx]=temp2[x];
//                }
//                
//            }
//            
//            //---for test only------
//            //if (i==1) {
//            //    for (int xx=0; xx<nx; xx++) {
//            //        cout<<xx<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//            //        PrintDat(temp2[xx]);
//            //    }
//            //}
//            //---for test only------
//        }
        
        //-----------end: find all the fitted bra's and store them----------------------
        
        
    }
*/
    
    void identity_find_kets(std::vector<ITensor> _W, int nx, int ny){
        
        for (int Y=1; Y<ny; Y++) {
                if (Y==1) {
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            IA[x]=IW[x];
//                        }
//                        else {
//                            bIA[x]=IW[x];
//                        }
                    }
//                    if (_backflag==0) {
                        IA=ket_otho_center_1(IA, Y-1, nx, ny);//change the othogonality center of A to site 1
//                    }
//                    else {
//                        bIA=ket_otho_center_1(bIA, Y-1, nx, ny);//change the othogonality center of A to site 1
//                    }
                    for (int x=0; x<nx; x++) {  //find Iket(0);
//                        if (_backflag==0) {
                            Iket[x+nx*(Y-1)]=IA[x];
                            _Iket[x]=IA[x];
//                        }
//                        else {
//                            bIket[x+nx*(Y-1)]=bIA[x];
//                            b_Iket[x]=bIA[x];
//                        }
                    }
                }
                else {
                    std::vector<ITensor> ket_temp(nx),gate_temp(nx),temp2(nx);
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            ket_temp[x]=Iket[x+nx*(Y-2)];
//                        }
//                        else {
//                            ket_temp[x]=bIket[x+nx*(Y-2)];
//                        }
                        gate_temp[x]=IW[x+nx*(Y-1)];
                    }
                    temp2=compressor1(ket_temp,gate_temp,Y-1,cutoff, m, dcut,nx,ny);
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            Iket[x+nx*(Y-1)]=temp2[x];
                            _Iket[x]=temp2[x];
//                        }
//                        else {
//                            bIket[x+nx*(Y-1)]=temp2[x];
//                            b_Iket[x]=temp2[x];
//                        }
                        
                    }
                }
        }
    }
    
    
    void ham_find_kets(std::vector<ITensor> _W, int i, int nx, int ny){
        
        for (int Y=1; Y<ny; Y++) {
                if (Y==1) {
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            AA[x]=W[x+nx*ny*i];
//                        }
//                        else {
//                            bAA[x]=W[x+nx*ny*i];
//                        }
                    }
//                    if (_backflag==0) {
                        AA=ket_otho_center_1(AA, Y-1,nx,ny);
//                    }
//                    else {
//                        bAA=ket_otho_center_1(bAA, Y-1,nx,ny);
//                    }
                    for (int x=0; x<nx; x++) {  //find ket(0)
//                        if (_backflag==0) {
                            ket[x+nx*(Y-1)+nx*(ny-1)*i]=AA[x];
                            _ket[x+nx*i]=AA[x];
//                        }
//                        else {
//                            bket[x+nx*(Y-1)+nx*(ny-1)*i]=bAA[x];
//                            b_ket[x+nx*i]=bAA[x];
//                        }
                    }
                }
                else {
                    std::vector<ITensor> ket_temp(nx),gate_temp(nx),temp2(nx);
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            ket_temp[x]=ket[x+nx*(Y-2)+nx*(ny-1)*i];
//                        }
//                        else {
//                            ket_temp[x]=bket[x+nx*(Y-2)+nx*(ny-1)*i];
//                        }
                        gate_temp[x]=W[x+nx*(Y-1)+N*i];
                    }
                    temp2=compressor1(ket_temp,gate_temp,Y-1,cutoff, m, dcut,nx,ny);
                    for (int x=0; x<nx; x++) {
//                        if (_backflag==0) {
                            ket[x+nx*(Y-1)+nx*(ny-1)*i]=temp2[x];
                            _ket[x+nx*i]=temp2[x];
//                        }
//                        else {
//                            bket[x+nx*(Y-1)+nx*(ny-1)*i]=temp2[x];
//                            b_ket[x+nx*i]=temp2[x];
//                        }
                    }
                }
        }
    }
    
    
}
#endif





/*
#ifndef ___TENSOR_SWEEP
#define ___TENSOR_SWEEP


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "tensor_utility.h"
#include "sweep_utility.h"
#include "parameters.h"
//#include "fitting.h"


using namespace std;
namespace itensor {

    //tensor_sweep(A,B,W,bra(ny)):
    void tensor_sweep(std::vector<ITensor> A, std::vector<ITensor> B, std::vector<ITensor> W){
        //cout<<"Hello!"<<endl;
        //---------begin: realizing main algorithm: to be wrapped into .h file--------------------------------
        std::vector<ITensor> bra(nx*(ny-1)), ket(nx*(ny-1));
        //for (int x=0; x<nx; x++) {
        //    PrintDat(B[x]);
        //}
        //Print(B[0].norm());
        //PrintDat(B[0]);
        
        B=bra_otho_center_1(B, ny-1);//change the othogonality center of B to site 1
        
        //PrintDat(B[0]);
        for (int x=0; x<nx; x++) {
            bra[x]=B[x];
        }

        
        std::vector<ITensor> temp1(nx),temp2(nx),gate(nx);
        for (int i=1; i<=ny-2; i++) {   //////////int i=1; i<=ny-2; i++
            for (int x=0; x<nx; x++) {
                temp1[x]=bra[x+nx*(i-1)];
                gate[x]=W[x+nx*(ny-2-i)];
            }
            
            temp2=compressor2(temp1,gate,i,cutoff, m, dcut);  //to make variable m effective, make change in compressor in sweep_utility.h
            
            for (int x=0; x<nx; x++) {
                bra[nx*i+x]=temp2[x];
            }
            
            //---for test only------
            //if (i==1) {
            //    for (int xx=0; xx<nx; xx++) {
            //        cout<<xx<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            //        PrintDat(temp2[xx]);
            //    }
            //}
            //---for test only------
        }
        
        //bra_otho_center_1(temp2, ny-2);
        
        //for (int x=0; x<nx; x++) {
        //    cout<<x<<"temp1~~~~~~~~~"<<endl;
        //    PrintDat(temp2[x]);
        //cout<<x<<"gate~~~~~~~~~"<<endl;
        //PrintDat(gate[x]);
        //}
        
        
        
        //for (int i=0; i<=ny-2; i++) {
        //    cout<<"UU"<<i<<"~~~~~~~~~"<<endl;
        //    for (int x=0; x<nx; x++) {
        //        cout<<x<<"~~~~~~~~~~~"<<endl;
        //        PrintDat(bra[nx*i+x]);
        //    }
        //}
         
        //---------end: realizing main algorithm--------------------------------------------------------------
    }

}
#endif
*/

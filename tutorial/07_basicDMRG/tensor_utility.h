//tensor_utility.h
//include change_index_name2/3/4, group_indice2/4/6/8, bra_otho_center_1, mpo_otho_center_1

#ifndef ___TENSOR_UTILITY
#define ___TENSOR_UTILITY


#include "core.h"
#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "sweep_utility.h"
#include "parameters.h"
#include "fitting.h"
#include <mpi.h>

#include "generalizedev.h"

using namespace std;
namespace itensor {
    //    if (COMPLEX==1){
    //
    //        //-------begin: resize tensors----------------------------------------------------
    //        ITensor resize2(ITensor N,Index b,int dcut){  //change the first index size to be dcut and rename it into b
    //            return N;
    //        }
    //        //-------begin: resize tensors----------------------------------------------------
    //
    //        //------
    //        ITensor sort_indices_braket(ITensor Z){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            if (Z.r()==2) {
    //                tempreal=ITensor(Z.indices()[1],Z.indices()[0]);
    //                tempimag=ITensor(Z.indices()[1],Z.indices()[0]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        tempreal(tempreal.indices()[0](i1),tempreal.indices()[1](i0))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1));
    //                        tempimag(tempimag.indices()[0](i1),tempimag.indices()[1](i0))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1));
    //                    }
    //                }
    //            }
    //            else if (Z.r()==3){
    //                tempreal=ITensor(Z.indices()[2],Z.indices()[0],Z.indices()[1]);
    //                tempimag=ITensor(Z.indices()[2],Z.indices()[0],Z.indices()[1]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
    //                            tempreal(tempreal.indices()[0](i2),tempreal.indices()[1](i0),tempreal.indices()[2](i1))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                            tempimag(tempimag.indices()[0](i2),tempimag.indices()[1](i0),tempimag.indices()[2](i1))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor sort_indices_compressor(ITensor Z){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            if (Z.r()==2) {
    //                tempreal=ITensor(Z.indices()[1],Z.indices()[0]);
    //                tempimag=ITensor(Z.indices()[1],Z.indices()[0]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        tempreal(tempreal.indices()[0](i1),tempreal.indices()[1](i0))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1));
    //                        tempimag(tempimag.indices()[0](i1),tempimag.indices()[1](i0))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1));
    //                    }
    //                }
    //            }
    //            else if (Z.r()==3){
    //                tempreal=ITensor(Z.indices()[1],Z.indices()[0],Z.indices()[2]);
    //                tempimag=ITensor(Z.indices()[1],Z.indices()[0],Z.indices()[2]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
    //                            tempreal(tempreal.indices()[0](i1),tempreal.indices()[1](i0),tempreal.indices()[2](i2))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                            tempimag(tempimag.indices()[0](i1),tempimag.indices()[1](i0),tempimag.indices()[2](i2))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor sort_indices_mpo(ITensor Z){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            if (Z.r()==3) {
    //                tempreal=ITensor(Z.indices()[1],Z.indices()[2],Z.indices()[0]);
    //                tempimag=ITensor(Z.indices()[1],Z.indices()[2],Z.indices()[0]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
    //                            tempreal(tempreal.indices()[0](i1),tempreal.indices()[1](i2),tempreal.indices()[2](i0))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                            tempimag(tempimag.indices()[0](i1),tempimag.indices()[1](i2),tempimag.indices()[2](i0))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
    //                        }
    //                    }
    //                }
    //            }
    //            else if (Z.r()==4){
    //                tempreal=ITensor(Z.indices()[2],Z.indices()[3],Z.indices()[0],Z.indices()[1]);
    //                tempimag=ITensor(Z.indices()[2],Z.indices()[3],Z.indices()[0],Z.indices()[1]);
    //                for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
    //                    for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
    //                        for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
    //                            for (int i3=1; i3<=Z.indices()[3].m(); i3++) {
    //                                tempreal(tempreal.indices()[0](i2),tempreal.indices()[1](i3),tempreal.indices()[2](i0),tempreal.indices()[3](i1))=realPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2),Z.indices()[3](i3));
    //                                tempimag(tempimag.indices()[0](i2),tempimag.indices()[1](i3),tempimag.indices()[2](i0),tempimag.indices()[3](i1))=imagPart(Z)(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2),Z.indices()[3](i3));
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor sort_indices_otho_2(ITensor Z, int a0, int a1){
    //            if (Z.r()!=2) {
    //                cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_2)"<<endl;
    //                exit(0);
    //            }
    //            ITensor tempreal(Z.indices()[a0],Z.indices()[a1]);
    //            ITensor tempimag(Z.indices()[a0],Z.indices()[a1]);
    //            for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
    //                for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
    //                    tempreal(tempreal.indices()[0](i0),tempreal.indices()[1](i1))=realPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1));
    //                    tempimag(tempimag.indices()[0](i0),tempimag.indices()[1](i1))=imagPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1));
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        ITensor sort_indices_otho_3(ITensor Z, int a0, int a1, int a2){
    //            if (Z.r()!=3) {
    //                cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_3)"<<endl;
    //                exit(0);
    //            }
    //            ITensor tempreal(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2]);
    //            ITensor tempimag(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2]);
    //            for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
    //                for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
    //                    for (int i2=1; i2<=Z.indices()[a2].m(); i2++) {
    //                        tempreal(tempreal.indices()[0](i0),tempreal.indices()[1](i1),tempreal.indices()[2](i2))=realPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2));
    //                        tempimag(tempimag.indices()[0](i0),tempimag.indices()[1](i1),tempimag.indices()[2](i2))=imagPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2));
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        ITensor sort_indices_otho_4(ITensor Z, int a0, int a1, int a2, int a3){
    //            if (Z.r()!=4) {
    //                cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_3)"<<endl;
    //                exit(0);
    //            }
    //            ITensor tempreal(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2], Z.indices()[a3]);
    //            ITensor tempimag(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2], Z.indices()[a3]);
    //            for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
    //                for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
    //                    for (int i2=1; i2<=Z.indices()[a2].m(); i2++) {
    //                        for (int i3=1; i3<=Z.indices()[a3].m(); i3++) {
    //                            tempreal(tempreal.indices()[0](i0),tempreal.indices()[1](i1),tempreal.indices()[2](i2),tempreal.indices()[3](i3))=realPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2),Z.indices()[a3](i3));
    //                            tempimag(tempimag.indices()[0](i0),tempimag.indices()[1](i1),tempimag.indices()[2](i2),tempimag.indices()[3](i3))=imagPart(Z)(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2),Z.indices()[a3](i3));
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        //------
    //
    //        //-----
    //        ITensor change_index_sequence2(ITensor N){
    //            if (N.r()!=2) {
    //                cout<<"Error: Number of indices in target tensor in not 2!!! (In func: change_index_sequence2) "<<endl;
    //                exit(0);
    //            }
    //            ITensor tempreal(N.indices()[1],N.indices()[0]);
    //            ITensor tempimag(N.indices()[1],N.indices()[0]);
    //            for (int i1=1; i1<=N.indices()[1].m(); i1++) {
    //                for (int i0=1; i0<=N.indices()[0].m(); i0++) {
    //                    tempreal(tempreal.indices()[0](i1),tempreal.indices()[1](i0))=realPart(N)(N.indices()[0](i0),N.indices()[1](i1));
    //                    tempimag(tempimag.indices()[0](i1),tempimag.indices()[1](i0))=imagPart(N)(N.indices()[0](i0),N.indices()[1](i1));
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        //-----
    //
    //        //-------begin: change names of indice--------------------------------------------
    //        ITensor change_index_name2(ITensor A, Index a1, Index a2){
    //            ITensor tempreal=ITensor(a1, a2);
    //            ITensor tempimag=ITensor(a1, a2);
    //            if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m()) {
    //                cout<<"Error: New indice do not match old indice!!! (in func: change_index_name2)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=a1.m(); i1++) {
    //                for (int i2=1; i2<=a2.m(); i2++) {
    //                    tempreal(a1(i1),a2(i2))=realPart(A)(A.indices()[0](i1),A.indices()[1](i2));
    //                    tempimag(a1(i1),a2(i2))=imagPart(A)(A.indices()[0](i1),A.indices()[1](i2));
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor change_index_name3(ITensor A, Index a1, Index a2, Index a3){
    //            ITensor tempreal=ITensor(a1, a2, a3);
    //            ITensor tempimag=ITensor(a1, a2, a3);
    //            if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m() || A.indices()[2].m()!=a3.m()) {
    //                cout<<"Error: New indice do not match old indice!!! (in func: change_index_name3)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=a1.m(); i1++) {
    //                for (int i2=1; i2<=a2.m(); i2++) {
    //                    for (int i3=1; i3<=a3.m(); i3++) {
    //                        tempreal(a1(i1),a2(i2),a3(i3))=realPart(A)(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3));
    //                        tempimag(a1(i1),a2(i2),a3(i3))=imagPart(A)(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3));
    //                    }
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor change_index_name4(ITensor A, Index a1, Index a2, Index a3, Index a4){
    //            ITensor tempreal=ITensor(a1, a2, a3, a4);
    //            ITensor tempimag=ITensor(a1, a2, a3, a4);
    //            if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m() || A.indices()[2].m()!=a3.m() || A.indices()[3].m()!=a4.m()) {
    //                cout<<"Error: New indice do not match old indice!!! (in func: change_index_name4)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=a1.m(); i1++) {
    //                for (int i2=1; i2<=a2.m(); i2++) {
    //                    for (int i3=1; i3<=a3.m(); i3++) {
    //                        for (int i4=1; i4<=a4.m(); i4++) {
    //                            tempreal(a1(i1),a2(i2),a3(i3),a4(i4))=realPart(A)(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3),A.indices()[3](i4));
    //                            tempimag(a1(i1),a2(i2),a3(i3),a4(i4))=imagPart(A)(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3),A.indices()[3](i4));
    //                        }
    //                    }
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        //-------end: change names of indice--------------------------------------------
    //
    //        //-------begin: tensor indice combiners-------------------------------------------
    //        ITensor group_indice2(ITensor Z, Index a1, Index a2, Index a12){
    //            ITensor tempreal=ITensor(a12);
    //            ITensor tempimag=ITensor(a12);
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m();
    //            if (ma1*ma2!=ma12) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice2)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    tempreal(a12(i12))=realPart(Z)(a1(i1),a2(i2));
    //                    tempimag(a12(i12))=imagPart(Z)(a1(i1),a2(i2));
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor group_indice3(ITensor Z, Index a1, Index a2, Index a12, Index a3){  //group index a1 and a2, leave a3 inact: Z_a1,a2,a3=M_(a1,a2),a3
    //            ITensor tempreal=ITensor(a12,a3);
    //            ITensor tempimag=ITensor(a12,a3);
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(),ma3=a3.m();
    //            if (ma1*ma2!=ma12) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice3)"<<endl;
    //                exit(0);
    //            }
    //
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    for (int i3=1; i3<=ma3; i3++) {
    //                        tempreal(a12(i12),a3(i3))=realPart(Z)(a1(i1),a2(i2),a3(i3));
    //                        tempimag(a12(i12),a3(i3))=imagPart(Z)(a1(i1),a2(i2),a3(i3));
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor group_indice4(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34){
    //            ITensor tempreal=ITensor(a12,a34);
    //            ITensor tempimag=ITensor(a12,a34);
    //
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m();
    //            if (ma1*ma2!=ma12||ma3*ma4!=ma34) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice4)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    for (int i3=1; i3<=ma3; i3++) {
    //                        for (int i4=1; i4<=ma4; i4++) {
    //                            int i34=(i3-1)*ma4+i4;
    //                            tempreal(a12(i12),a34(i34))=realPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4));
    //                            tempimag(a12(i12),a34(i34))=imagPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4));
    //                            //temp(a12(i12),a34(i34))=Z(a1(i1),a2(i2),a3(i3),a4(i4));
    //                        }
    //                    }
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor group_indice5(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5){  //group index a1 and a2, a3 and a4, leave a5 inact: Z_a1,a2,a3,a4,a5=M_(a1,a2),(a3,a4),a5
    //            ITensor tempreal=ITensor(a12,a34,a5);
    //            ITensor tempimag=ITensor(a12,a34,a5);
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(),ma3=a3.m(),ma4=a4.m(),ma34=a34.m(),ma5=a5.m();
    //            if (ma1*ma2!=ma12 || ma3*ma4!=ma34) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice5)"<<endl;
    //                exit(0);
    //            }
    //
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    for (int i3=1; i3<=ma3; i3++) {
    //                        for (int i4=1; i4<=ma4; i4++) {
    //                            int i34=(i3-1)*ma4+i4;
    //
    //                            for (int i5=1; i5<=ma5; i5++) {
    //                                tempreal(a12(i12),a34(i34),a5(i5))=realPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5));
    //                                tempimag(a12(i12),a34(i34),a5(i5))=imagPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5));
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor group_indice6(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5, Index a6, Index a56){
    //            ITensor tempreal=ITensor(a12,a34,a56);
    //            ITensor tempimag=ITensor(a12,a34,a56);
    //
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m(), ma5=a5.m(),ma6=a6.m(),ma56=a56.m();
    //            if (ma1*ma2!=ma12||ma3*ma4!=ma34||ma5*ma6!=ma56) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice6)"<<endl;
    //                exit(0);
    //            }
    //
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    for (int i3=1; i3<=ma3; i3++) {
    //                        for (int i4=1; i4<=ma4; i4++) {
    //                            int i34=(i3-1)*ma4+i4;
    //                            for (int i5=1; i5<=ma5; i5++) {
    //                                for (int i6=1; i6<=ma6; i6++) {
    //                                    int i56=(i5-1)*ma6+i6;
    //                                    tempreal(a12(i12),a34(i34),a56(i56))=realPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6));
    //                                    tempimag(a12(i12),a34(i34),a56(i56))=imagPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6));
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor group_indice8(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5, Index a6, Index a56, Index a7, Index a8, Index a78){
    //            ITensor tempreal=ITensor(a12,a34,a56,a78);
    //            ITensor tempimag=ITensor(a12,a34,a56,a78);
    //            int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m(), ma5=a5.m(),ma6=a6.m(),ma56=a56.m(), ma7=a7.m(),ma8=a8.m(),ma78=a78.m();
    //            if (ma1*ma2!=ma12||ma3*ma4!=ma34||ma5*ma6!=ma56||ma7*ma8!=ma78) {
    //                cout<<"Error: Tensor sizes do not match!!! (in func: group_indice8)"<<endl;
    //                exit(0);
    //            }
    //            for (int i1=1; i1<=ma1; i1++) {
    //                for (int i2=1; i2<=ma2; i2++) {
    //                    int i12=(i1-1)*ma2+i2;
    //                    for (int i3=1; i3<=ma3; i3++) {
    //                        for (int i4=1; i4<=ma4; i4++) {
    //                            int i34=(i3-1)*ma4+i4;
    //                            for (int i5=1; i5<=ma5; i5++) {
    //                                for (int i6=1; i6<=ma6; i6++) {
    //                                    int i56=(i5-1)*ma6+i6;
    //                                    for (int i7=1; i7<=ma7; i7++) {
    //                                        for (int i8=1; i8<=ma8; i8++) {
    //                                            int i78=(i7-1)*ma8+i8;
    //                                            tempreal(a12(i12),a34(i34),a56(i56),a78(i78))=realPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6),a7(i7),a8(i8));
    //                                            tempimag(a12(i12),a34(i34),a56(i56),a78(i78))=imagPart(Z)(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5),a6(i6),a7(i7),a8(i8));
    //                                        }
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //            return tempreal+Complex_i*tempimag;
    //        }
    //        //-------end: tensor indice combiners-------------------------------------------
    //
    //        //-------begin: move othogonality center to the first site for bra, without cutoff----------------------
    //        std::vector<ITensor> bra_otho_center_1(std::vector<ITensor> Z, int y){
    //
    //            int size=Z.size();
    //            int u_ind=-1,l_ind=-1,r_ind=-1;
    //            ITensor S, V, D;
    //            int aaa;
    //
    //            std::vector<ITensor> original(size);
    //            for (int x=0; x<size; x++) {
    //                original[x]=Z[x];
    //            }
    //
    //            D=ITensor(original[size-1].indices()[0]);
    //            svd(Z[size-1], S, V, D);    //svd
    //
    //            //find out indice names and sequence
    //            for (int i=0; i<2; i++) {
    //                if (D.indices()[i]==original[size-1].indices()[0]) {
    //                    u_ind=i;
    //                }
    //            }
    //            l_ind=1-u_ind;
    //
    //            L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
    //
    //            //Z=D, and rename indices as U, L, R indice
    //            if (u_ind==0) {
    //                Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[size-1]=change_index_name2(D, L[size-1+nx*y], original[size-1].indices()[0]);
    //            }
    //
    //            Z[size-1]=sort_indices_otho_2(Z[size-1], u_ind, l_ind);
    //
    //            //        D=sort_indices_otho_2(D, u_ind, l_ind);
    //            //        Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
    //
    //            for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[x]=Z[x]*S*V;
    //                //if (x==size-3) {
    //                //    PrintDat(Z[x]);
    //                //}
    //
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (Z[x].indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (Z[x].indices()[i]==original[x].indices()[1]) {
    //                        l_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                r_ind=3-aaa;
    //                R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
    //
    //
    //                //PrintDat(Z[x]);
    //                //rename indices of Z as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
    //                    }
    //                }
    //
    //                //PrintDat(Z[x]);
    //                S=ITensor(original[x].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
    //                svd(Z[x], S, V, D); //svd
    //                //PrintDat(S);
    //                //PrintDat(V);
    //                //PrintDat(D);
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (D.indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (D.indices()[i]==R[x+nx*y]) {
    //                        r_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                l_ind=3-aaa;
    //
    //                L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
    //
    //                //Z=D, and rename indices as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
    //                    }
    //                }
    //
    //                Z[x]=sort_indices_otho_3(Z[x], u_ind, l_ind, r_ind);
    //
    //                //            D=sort_indices_otho_3(D, u_ind, l_ind, r_ind);
    //                //            Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
    //                //PrintDat(D);
    //                //PrintDat(Z[x]);
    //
    //            }
    //
    //            Z[0]=Z[0]*S*V;
    //            for (int i=0; i<2; i++) {
    //                if (Z[0].indices()[i]==original[0].indices()[0]) {
    //                    u_ind=i;
    //                }
    //            }
    //            r_ind=1-u_ind;
    //            R[0+nx*y]=L[1+nx*y];
    //            if (u_ind==0) {
    //                Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[0]=change_index_name2(Z[0], R[0+nx*y], original[0].indices()[0]);
    //            }
    //
    //            Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
    //
    //            //        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
    //            //        Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<x<<"~~~~~~~~~"<<endl;
    //            //    PrintDat(Z[x]);
    //            //}
    //            //        for (int x=0; x<size; x++) {
    //            //            Z[x]=sort_indices_braket(Z[x]);
    //            //        }
    //
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //
    //        //-------begin: move othogonality center to the first site for bra, without cutoff----------------------
    //        std::vector<ITensor> bra_otho_center_1_b(std::vector<ITensor> Z, int x){
    //
    //            int size=Z.size();
    //            int l_ind=-1,u_ind=-1,d_ind=-1;
    //            ITensor _S, _V, _D;
    //            int aaa;
    //
    //            std::vector<ITensor> original(size);
    //            for (int y=0; y<size; y++) {
    //                original[y]=Z[y];
    //            }
    //
    //            _D=ITensor(original[size-1].indices()[1]);
    //            svd(Z[size-1], _S, _V, _D);    //svd
    //
    //            U[x+nx*(size-1)]=Index(nameint("U",x+nx*(size-1)),_D.indices()[0].m());  //resize index, to be used for the updated Z
    //
    //            //Z=D, and rename indices as U, L, R indice
    //            Z[size-1]=change_index_name2(_D, U[x+nx*(size-1)], Z[size-1].indices()[1]);
    //
    //            for (int y=size-2; y>0; y--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[y]=Z[y]*_S*_V;
    //
    //                D[x+nx*y]=U[x+nx*(y+1)];
    //                Z[y]=change_index_name3(Z[y], D[x+nx*y], Z[y].indices()[1], Z[y].indices()[2]);
    //                _S=ITensor(Z[y].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
    //                svd(Z[y], _S, _V, _D); //svd
    //
    //                U[x+nx*y]=Index(nameint("U",x+nx*y),_D.indices()[0].m());  //resize index, to be used for the updated Z
    //                Z[y]=change_index_name3(_D, U[x+nx*y], _D.indices()[1], _D.indices()[2]);
    //
    //            }
    //
    //            Z[0]=Z[0]*_S*_V;
    //            D[x]=U[x+nx*1];
    //            Z[0]=change_index_name2(Z[0], D[x], Z[0].indices()[1]);
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //
    //        //-------begin: move othogonality center to the first site for bra, with dimension cutoff m----------------------
    //        std::vector<ITensor> bra_otho_center_1m(std::vector<ITensor> Z, int y, int m){
    //            int size=Z.size();
    //            int u_ind=-1,l_ind=-1,r_ind=-1;
    //            ITensor S, V, D;
    //            int aaa;
    //
    //            D=ITensor(U[size-1+nx*y]);
    //
    //            OptSet opts;
    //            opts.add("Maxm", m);
    //            svd(Z[size-1], S, V, D, opts);    //svd
    //
    //            //find out indice names and sequence
    //            for (int i=0; i<2; i++) {
    //                if (D.indices()[i]==U[size-1+nx*y]) {
    //                    u_ind=i;
    //                }
    //            }
    //            l_ind=1-u_ind;
    //
    //            L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
    //
    //            //Z=D, and rename indices as U, L, R indice
    //            if (u_ind==0) {
    //                Z[size-1]=change_index_name2(D, U[size-1+nx*y], L[size-1+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[size-1]=change_index_name2(D, L[size-1+nx*y], U[size-1+nx*y]);
    //            }
    //
    //            for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[x]=Z[x]*S*V;
    //
    //                //cout<<"x="<<x<<endl;
    //                //if (x==size-3) {
    //                //    Print(Z[x]);
    //                //}
    //
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (Z[x].indices()[i]==U[x+nx*y]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (Z[x].indices()[i]==L[x+nx*y]) {
    //                        l_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                r_ind=3-aaa;
    //                R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
    //
    //
    //                //PrintDat(Z[x]);
    //                //rename indices of Z as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(Z[x], U[x+nx*y], L[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], U[x+nx*y], R[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], L[x+nx*y], U[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], U[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], L[x+nx*y], R[x+nx*y], U[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], L[x+nx*y], U[x+nx*y]);
    //                    }
    //                }
    //
    //                //PrintDat(Z[x]);
    //                S=ITensor(L[x+nx*y]);   //let L index be in S tensor, U and R be in D tensor
    //                svd(Z[x], S, V, D, opts); //svd
    //                //if (x==2) {
    //                //    Print(Z[x]);
    //                //    Print(S);
    //                //    Print(V);
    //                //    Print(D);
    //                //}
    //
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (D.indices()[i]==U[x+nx*y]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (D.indices()[i]==R[x+nx*y]) {
    //                        r_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                l_ind=3-aaa;
    //
    //                L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
    //
    //                //Z=D, and rename indices as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(D, U[x+nx*y], L[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, U[x+nx*y], R[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], U[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], U[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], U[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], U[x+nx*y]);
    //                    }
    //                }
    //                //PrintDat(D);
    //                //PrintDat(Z[x]);
    //
    //            }
    //
    //            Z[0]=Z[0]*S*V;
    //            for (int i=0; i<2; i++) {
    //                if (Z[0].indices()[i]==U[0+nx*y]) {
    //                    u_ind=i;
    //                }
    //            }
    //            r_ind=1-u_ind;
    //            R[0+nx*y]=L[1+nx*y];
    //            if (u_ind==0) {
    //                Z[0]=change_index_name2(Z[0], U[0+nx*y], R[0+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[0]=change_index_name2(Z[0], R[0+nx*y], U[0+nx*y]);
    //            }
    //
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<x<<"~~~~~~~~~"<<endl;
    //            //    PrintDat(Z[x]);
    //            //}
    //
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //        //-------begin: move othogonality center to the first site for ket, without cutoff----------------------
    //        std::vector<ITensor> ket_otho_center_1(std::vector<ITensor> Z, int y){
    //            int size=Z.size();
    //            int u_ind=-1,l_ind=-1,r_ind=-1;
    //            ITensor S, V, D;
    //            int aaa;
    //
    //
    //            std::vector<ITensor> original(size);
    //            for (int x=0; x<size; x++) {
    //                original[x]=Z[x];
    //            }
    //
    //            D=ITensor(original[size-1].indices()[0]);   //for ket: U[size-1+nx+nx*y]=D[size-1+nx*y]
    //            svd(Z[size-1], S, V, D);    //svd
    //
    //            //find out indice names and sequence
    //            for (int i=0; i<2; i++) {
    //                if (D.indices()[i]==original[size-1].indices()[0]) {
    //                    u_ind=i;
    //                }
    //            }
    //            l_ind=1-u_ind;
    //
    //            L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
    //
    //            if (u_ind==0) {
    //                Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[size-1]=change_index_name2(D, L[size-1+nx*y], original[size-1].indices()[0]);
    //            }
    //
    //            Z[size-1]=sort_indices_otho_2(Z[size-1], u_ind, l_ind);
    //
    //
    //            //        D=sort_indices_otho_2(D, u_ind, l_ind);
    //            //        Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
    //
    //
    //            for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[x]=Z[x]*S*V;
    //
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (Z[x].indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (Z[x].indices()[i]==original[x].indices()[1]) {
    //                        l_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                r_ind=3-aaa;
    //                R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
    //
    //
    //                //rename indices of Z as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
    //                    }
    //                }
    //
    //                S=ITensor(original[x].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
    //                svd(Z[x], S, V, D); //svd
    //
    //                aaa=0;
    //                for (int i=0; i<3; i++) {
    //                    if (D.indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (D.indices()[i]==R[x+nx*y]) {
    //                        r_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                l_ind=3-aaa;
    //
    //                L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
    //
    //
    //                //Z=D, and rename indices as U,L,R
    //                if (u_ind==0) {
    //                    if (l_ind==1) {
    //                        Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
    //                    }
    //                }
    //                else {
    //                    if (l_ind==0) {
    //                        Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
    //                    }
    //                    else {
    //                        Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
    //                    }
    //                }
    //
    //                Z[x]=sort_indices_otho_3(Z[x], u_ind, l_ind, r_ind);
    //
    //
    //                //            D=sort_indices_otho_3(D, u_ind, l_ind, r_ind);
    //                //            Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
    //
    //
    //
    //            }
    //
    //
    //            Z[0]=Z[0]*S*V;
    //
    //            for (int i=0; i<2; i++) {
    //                if (Z[0].indices()[i]==original[0].indices()[0]) {
    //                    u_ind=i;
    //                }
    //            }
    //            r_ind=1-u_ind;
    //            R[0+nx*y]=L[1+nx*y];
    //
    //
    //            if (u_ind==0) {
    //                Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
    //            }
    //            if (u_ind==1) {
    //                Z[0]=change_index_name2(Z[0], R[0+nx*y], original[0].indices()[0]);
    //            }
    //
    //            Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
    //
    //
    //            //        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
    //            //        Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
    //
    //
    //            /*
    //             for (int x=0; x<nx; x++) {
    //             cout<<"x="<<x;
    //             Print(Z[x]);
    //             }
    //             */
    //
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //
    //
    //
    //        //-------begin: move othogonality center to the first site for mpo----------------------
    //        std::vector<ITensor> mpo_otho_center_1(std::vector<ITensor> Z, int y){
    //            int size=Z.size();
    //            int u_ind=-1,d_ind=-1,l_ind=-1,r_ind=-1;
    //            ITensor S, V, DD;
    //
    //            std::vector<ITensor> original(size);
    //            for (int x=0; x<size; x++) {
    //                original[x]=Z[x];
    //            }
    //
    //            S=ITensor(original[size-1].indices()[2]);
    //            svd(Z[size-1], S, V, DD);
    //            int aaa=0;
    //            for (int i=0; i<3; i++) {
    //                if (DD.indices()[i]==original[size-1].indices()[0]) {
    //                    u_ind=i;
    //                    aaa+=i;
    //                }
    //                if (DD.indices()[i]==original[size-1].indices()[1]) {
    //                    d_ind=i;
    //                    aaa+=i;
    //                }
    //            }
    //            l_ind=3-aaa;
    //
    //            L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),DD.indices()[l_ind].m());  //resize index, to be used for the updated Z
    //
    //            if (u_ind==0) {
    //                if (d_ind==1) {
    //                    Z[size-1]=change_index_name3(DD, original[size-1].indices()[0], original[size-1].indices()[1], L[size-1+nx*y]);
    //                }
    //                else {
    //                    Z[size-1]=change_index_name3(DD, original[size-1].indices()[0], L[size-1+nx*y], original[size-1].indices()[1]);
    //                }
    //            }
    //            else if (u_ind==1) {
    //                if (d_ind==0) {
    //                    Z[size-1]=change_index_name3(DD, original[size-1].indices()[1], original[size-1].indices()[0], L[size-1+nx*y]);
    //                }
    //                else {
    //                    Z[size-1]=change_index_name3(DD, L[size-1+nx*y], original[size-1].indices()[0], original[size-1].indices()[1]);
    //                }
    //            }
    //            else {
    //                if (d_ind==0) {
    //                    Z[size-1]=change_index_name3(DD, original[size-1].indices()[1], L[size-1+nx*y], original[size-1].indices()[0]);
    //                }
    //                else {
    //                    Z[size-1]=change_index_name3(DD, L[size-1+nx*y], original[size-1].indices()[1], original[size-1].indices()[0]);
    //                }
    //            }
    //            //PrintDat(S);
    //            //PrintDat(V);
    //
    //
    //            for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[x]=Z[x]*S*V;
    //                //PrintDat(Z[x]);
    //
    //                //if (x==size-3) {
    //                //    PrintDat(Z[x]);
    //                //}
    //
    //
    //                int aaa=0;
    //                for (int i=0; i<4; i++) {
    //                    if (Z[x].indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (Z[x].indices()[i]==original[x].indices()[1]) {
    //                        d_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (Z[x].indices()[i]==original[x].indices()[2]) {
    //                        l_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                r_ind=6-aaa;
    //
    //                R[x+nx*y]=L[x+1+nx*y];
    //                //PrintDat(Z[x]);
    //
    //                ///rename indices of Z[x] to the wanted format
    //                if (u_ind==0) {
    //                    if (d_ind==1) {
    //                        if (l_ind==2) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[1], original[x].indices()[2], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else if (d_ind==2) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[2], original[x].indices()[1], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[2], R[x+nx*y], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[2], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (d_ind==0) {
    //                        if (l_ind==2) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[0], original[x].indices()[2], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else if (d_ind==2) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[2], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else if (u_ind==2) {
    //                    if (d_ind==0) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[2], original[x].indices()[0], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else if (d_ind==1) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0], original[x].indices()[2]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[2], original[x].indices()[0], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else {
    //                    if (d_ind==0) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[2], R[x+nx*y], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[2], original[x].indices()[0]);
    //                        }
    //                    }
    //                    else if (d_ind==1) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[2], original[x].indices()[0]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(Z[x], original[x].indices()[2], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[2], original[x].indices()[1], original[x].indices()[0]);
    //                        }
    //                    }
    //                }
    //
    //                //PrintDat(Z[x]);
    //
    //                S=ITensor(original[x].indices()[2]);
    //                svd(Z[x], S, V, DD);
    //                //PrintDat(S);
    //                //PrintDat(V);
    //                //PrintDat(DD);
    //
    //                aaa=0;
    //                for (int i=0; i<4; i++) {
    //                    if (DD.indices()[i]==original[x].indices()[0]) {
    //                        u_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (DD.indices()[i]==original[x].indices()[1]) {
    //                        d_ind=i;
    //                        aaa+=i;
    //                    }
    //                    if (DD.indices()[i]==R[x+nx*y]) {
    //                        r_ind=i;
    //                        aaa+=i;
    //                    }
    //                }
    //                l_ind=6-aaa;
    //
    //                L[x+nx*y]=Index(nameint("L",x+nx*y),DD.indices()[l_ind].m());    //resize index, to be used for the updated Z
    //
    //                if (u_ind==0) {
    //                    if (d_ind==1) {
    //                        if (l_ind==2) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], original[x].indices()[1], L[x+nx*y], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], original[x].indices()[1], R[x+nx*y], L[x+nx*y]);
    //                        }
    //                    }
    //                    else if (d_ind==2) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], L[x+nx*y], original[x].indices()[1], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], R[x+nx*y], original[x].indices()[1], L[x+nx*y]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], L[x+nx*y], R[x+nx*y], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[0], R[x+nx*y], L[x+nx*y], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else if (u_ind==1) {
    //                    if (d_ind==0) {
    //                        if (l_ind==2) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
    //                        }
    //                    }
    //                    else if (d_ind==2) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[0], original[x].indices()[1], L[x+nx*y]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[0], L[x+nx*y], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else if (u_ind==2) {
    //                    if (d_ind==0) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
    //                        }
    //                    }
    //                    else if (d_ind==1) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[1], original[x].indices()[0], L[x+nx*y]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], L[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
    //                        }
    //                    }
    //                }
    //                else {
    //                    if (d_ind==0) {
    //                        if (l_ind==1) {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, original[x].indices()[1], R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
    //                        }
    //                    }
    //                    else if (d_ind==1) {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[1], L[x+nx*y], original[x].indices()[0]);
    //                        }
    //                    }
    //                    else {
    //                        if (l_ind==0) {
    //                            Z[x]=change_index_name4(DD, L[x+nx*y], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
    //                        }
    //                        else {
    //                            Z[x]=change_index_name4(DD, R[x+nx*y], L[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
    //                        }
    //                    }
    //                }
    //
    //                //PrintDat(DD);
    //                //PrintDat(Z[x]);
    //
    //
    //            }
    //
    //            Z[0]=Z[0]*S*V;
    //            aaa=0;
    //            for (int i=0; i<3; i++) {
    //                if (Z[0].indices()[i]==original[0].indices()[0]) {
    //                    u_ind=i;
    //                    aaa+=i;
    //                }
    //                if (Z[0].indices()[i]==original[0].indices()[1]) {
    //                    d_ind=i;
    //                    aaa+=i;
    //                }
    //            }
    //            r_ind=3-aaa;
    //
    //            R[0+nx*y]=L[1+nx*y];
    //
    //            if (u_ind==0) {
    //                if (d_ind==1) {
    //                    Z[0]=change_index_name3(Z[0], original[0].indices()[0], original[0].indices()[1], R[0+nx*y]);
    //                }
    //                else {
    //                    Z[0]=change_index_name3(Z[0], original[0].indices()[0], R[0+nx*y], original[0].indices()[1]);
    //                }
    //            }
    //            else if (u_ind==1) {
    //                if (d_ind==0) {
    //                    Z[0]=change_index_name3(Z[0], original[0].indices()[1], original[0].indices()[0], R[0+nx*y]);
    //                }
    //                else {
    //                    Z[0]=change_index_name3(Z[0], R[0+nx*y], original[0].indices()[0], original[0].indices()[1]);
    //                }
    //            }
    //            else {
    //                if (d_ind==0) {
    //                    Z[0]=change_index_name3(Z[0], original[0].indices()[1], R[0+nx*y], original[0].indices()[0]);
    //                }
    //                else {
    //                    Z[0]=change_index_name3(Z[0], R[0+nx*y], original[0].indices()[1], original[0].indices()[0]);
    //                }
    //            }
    //
    //            for (int x=0; x<size; x++) {
    //                Z[x]=sort_indices_mpo(Z[x]);
    //            }
    //
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //
    //        //-------begin: move othogonality center to the first site for mpo----------------------
    //        std::vector<ITensor> mpo_otho_center_1b(std::vector<ITensor> Z, int x){
    //            int size=Z.size();
    //            ITensor S, V, DD;
    //
    //            std::vector<ITensor> original(size);
    //            for (int y=0; y<size; y++) {
    //                original[y]=Z[y];
    //            }
    //
    //            S=ITensor(original[size-1].indices()[0]);
    //            svd(Z[size-1], S, V, DD);
    //
    //            U[x+nx*(size-1)]=Index(nameint("U",x+nx*(size-1)),DD.indices()[0].m());  //resize index, to be used for the updated Z
    //            //PrintDat(S);
    //            //PrintDat(V);
    //            Z[size-1]=change_index_name3(DD, U[x+nx*(size-1)], DD.indices()[1], DD.indices()[2]);
    //
    //            for (int y=size-2; y>0; y--) {  //int x=size-2; x>0; x--
    //                ///change index name for Z[x]//////
    //                Z[y]=Z[y]*S*V;
    //                //PrintDat(Z[x]);
    //                //if (x==size-3) {
    //                //    PrintDat(Z[x]);
    //                //}
    //
    //                D[x+nx*y]=U[x+nx*(y+1)];
    //                //PrintDat(Z[x]);
    //
    //                Z[y]=change_index_name4(Z[y], D[x+nx*y], Z[y].indices()[1], Z[y].indices()[2], Z[y].indices()[3]);
    //
    //                S=ITensor(Z[y].indices()[1]);
    //                svd(Z[y], S, V, DD);
    //
    //                U[x+nx*y]=Index(nameint("U",x+nx*y),DD.indices()[0].m());
    //
    //                Z[y]=change_index_name4(DD, U[x+nx*y], DD.indices()[1], DD.indices()[2], DD.indices()[3]);
    //
    //                //PrintDat(DD);
    //            }
    //
    //
    //            Z[0]=Z[0]*S*V;
    //            D[x]=U[x+nx];
    //            Z[0]=change_index_name3(Z[0], D[x], Z[0].indices()[1], Z[0].indices()[2]);
    //
    //            return Z;
    //        }
    //        //-------end: move othogonality center to the first site----------------------
    //
    //
    //        //-------begin: fitting algorithm---------------------------------------------
    //        std::vector<ITensor> fitting_algorithm(std::vector<ITensor> zipped_mps, std::vector<ITensor> mpo, std::vector<ITensor> original_mps, int dcut, int i){
    //            int zsize=zipped_mps.size();
    //            int msize=mpo.size();
    //            int osize=original_mps.size();
    //            int size;
    //            if (zsize==msize&&msize==osize) {
    //                size=zsize;
    //            }
    //            else {
    //                cout<<"ERROR: Sizes of zipped mps, mpo and the original mps does not match!!! (in func: fitting_algorithm)"<<endl;
    //                exit(0);
    //            }
    //
    //            int wind[nx];  //B index and W index
    //            for (int x=0; x<nx; x++) {
    //                wind[x]=x+nx*(ny-1-i);
    //            }
    //
    //            std::vector<ITensor> RRR(size-2), LLL(size-2);
    //            ITensor total;
    //
    //            ITensor SSS,VVV,DDD;
    //
    //
    //            ITensor psi_a,psi_b,k_aa,k_bb,k_ab;
    //            /*
    //             //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //             if (size<8) {
    //             cout<<"Sweep step 0:    ";
    //             psi_a=mpo[0]*original_mps[0];
    //             psi_b=zipped_mps[0];
    //             //Print(psi_a);
    //             //Print(psi_b);
    //             for (int x=1; x<size; x++) {
    //             psi_a=psi_a*(mpo[x]*original_mps[x]);
    //             psi_b=psi_b*zipped_mps[x];
    //             }
    //             //Print(psi_a.norm());
    //             //Print(psi_a);
    //             //Print(psi_b);
    //             //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //             printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //             cout<<endl;
    //             //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //             }
    //             //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //             */
    //            //cout<<"sweep0"<<endl;
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<"x="<<x<<endl;
    //            //    Print(zipped_mps[x]);
    //            //}
    //
    //
    //            for (int sw=1; sw<=1; sw++) {
    //
    //                //~~~~~~~~cout<<"Forward scanning......"<<endl;
    //                //forward sweep
    //                RRR[0]=zipped_mps[size-1]*mpo[size-1]*original_mps[size-1];
    //                for (int k=1; k<size-2; k++) {
    //                    RRR[k]=RRR[k-1]*(zipped_mps[size-1-k]*mpo[size-1-k]*original_mps[size-1-k]);
    //                }
    //
    //                total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);
    //
    //                SSS=ITensor(total.indices()[2]);
    //                svd(total, SSS,VVV,DDD);
    //                zipped_mps[0]=SSS;
    //                zipped_mps[1]=VVV*DDD;
    //
    //                LLL[0]=zipped_mps[0]*mpo[0]*original_mps[0];
    //                /*
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 if (size<8) {
    //                 cout<<"Sweep step "<<1+(2*size-3)*(sw-1)<<":    ";
    //                 psi_a=mpo[0]*original_mps[0];
    //                 psi_b=zipped_mps[0];
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 for (int x=1; x<size; x++) {
    //                 psi_a=psi_a*(mpo[x]*original_mps[x]);
    //                 psi_b=psi_b*zipped_mps[x];
    //                 }
    //                 //Print(psi_a.norm());
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                 printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                 cout<<endl;
    //                 //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                 }
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 */
    //                for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
    //                    total=RRR[size-3-k]*(mpo[1+k]*original_mps[1+k])*(mpo[0+k]*original_mps[0+k])*LLL[k-1];
    //                    SSS=ITensor(total.indices()[2], total.indices()[3]);
    //                    svd(total,SSS,VVV,DDD);
    //                    zipped_mps[0+k]=SSS;
    //                    zipped_mps[1+k]=VVV*DDD;
    //                    LLL[0+k]=LLL[0+k-1]*(zipped_mps[0+k]*mpo[0+k]*original_mps[0+k]);
    //                    /*
    //                     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                     if (size<8) {
    //                     cout<<"Sweep step "<<k+1+(2*size-3)*(sw-1)<<":    ";
    //                     psi_a=mpo[0]*original_mps[0];
    //                     psi_b=zipped_mps[0];
    //                     //Print(psi_a);
    //                     //Print(psi_b);
    //                     for (int x=1; x<size; x++) {
    //                     psi_a=psi_a*(mpo[x]*original_mps[x]);
    //                     psi_b=psi_b*zipped_mps[x];
    //                     }
    //                     //Print(psi_a.norm());
    //                     //Print(psi_a);
    //                     //Print(psi_b);
    //                     //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                     printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                     cout<<endl;
    //                     //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                     }
    //                     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                     */
    //                }
    //
    //                total=LLL[size-3]*(mpo[size-2]*original_mps[size-2])*(mpo[size-1]*original_mps[size-1]);
    //                SSS=ITensor(total.indices()[0],total.indices()[1]);
    //                svd(total,SSS,VVV,DDD);
    //                zipped_mps[size-2]=SSS*VVV;
    //                zipped_mps[size-1]=DDD;
    //                /*
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 if (size<8) {
    //                 cout<<"Sweep step "<<(size-1)+(2*size-3)*(sw-1)<<":    ";
    //                 psi_a=mpo[0]*original_mps[0];
    //                 psi_b=zipped_mps[0];
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 for (int x=1; x<size; x++) {
    //                 psi_a=psi_a*(mpo[x]*original_mps[x]);
    //                 psi_b=psi_b*zipped_mps[x];
    //                 }
    //                 //Print(psi_a.norm());
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                 printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                 cout<<endl;
    //                 //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                 }
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 */
    //                //cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //                //for (int x=0; x<size; x++) {
    //                //    cout<<"x="<<x<<endl;
    //                //    Print(zipped_mps[x]);
    //                //}
    //                //~~~~~~~~~~~~cout<<"Backward scanning......"<<endl;
    //                //backward sweep
    //                for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
    //                    if (k==1) {
    //                        RRR[k-1]=zipped_mps[size-k]*mpo[size-k]*original_mps[size-k];
    //                    }
    //                    else {
    //                        RRR[k-1]=RRR[k-2]*(zipped_mps[size-k]*mpo[size-k]*original_mps[size-k]);
    //                    }
    //                    total=RRR[k-1]*(mpo[size-1-k]*original_mps[size-1-k])*(mpo[size-2-k]*original_mps[size-2-k])*LLL[size-3-k];
    //                    SSS=ITensor(total.indices()[2],total.indices()[3]);
    //                    svd(total,SSS,VVV,DDD);
    //                    zipped_mps[size-2-k]=SSS*VVV;
    //                    zipped_mps[size-1-k]=DDD;
    //                    /*
    //                     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                     if (size<8) {
    //                     cout<<"Sweep step "<<(size-1+k)+(2*size-3)*(sw-1)<<":    ";
    //                     psi_a=mpo[0]*original_mps[0];
    //                     psi_b=zipped_mps[0];
    //                     //Print(psi_a);
    //                     //Print(psi_b);
    //                     for (int x=1; x<size; x++) {
    //                     psi_a=psi_a*(mpo[x]*original_mps[x]);
    //                     psi_b=psi_b*zipped_mps[x];
    //                     }
    //                     //Print(psi_a.norm());
    //                     //Print(psi_a);
    //                     //Print(psi_b);
    //                     //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                     printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                     cout<<endl;
    //                     //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                     }
    //                     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                     */
    //                }
    //
    //                if (size>3) {
    //                    RRR[size-3]=RRR[size-4]*(zipped_mps[2]*mpo[2]*original_mps[2]);
    //                }
    //                else if (size==3){
    //                    RRR[size-3]=(zipped_mps[2]*mpo[2]*original_mps[2]);
    //                }
    //
    //                total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);
    //                SSS=ITensor(total.indices()[2]);
    //                svd(total,SSS,VVV,DDD);
    //                zipped_mps[0]=SSS*VVV;
    //                zipped_mps[1]=DDD;
    //                /*
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 if (size<8) {
    //                 cout<<"Sweep step "<<(2*size-3)+(2*size-3)*(sw-1)<<":    ";
    //                 psi_a=mpo[0]*original_mps[0];
    //                 psi_b=zipped_mps[0];
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 for (int x=1; x<size; x++) {
    //                 psi_a=psi_a*(mpo[x]*original_mps[x]);
    //                 psi_b=psi_b*zipped_mps[x];
    //                 }
    //                 //Print(psi_a.norm());
    //                 //Print(psi_a);
    //                 //Print(psi_b);
    //                 //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                 printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                 cout<<endl;
    //                 //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                 }
    //                 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                 */
    //                //zipped_mps[0]=change_index_sequence2(zipped_mps[0]);
    //
    //                //            cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //                //            for (int x=0; x<size; x++) {
    //                //                cout<<"x="<<x<<endl;
    //                //                Print(zipped_mps[x]);
    //                //            }
    //
    //                //rename index names makes it easy to debug
    //                for (int x=0; x<size-1; x++) {
    //                    R[wind[x]]=Index(nameint("R",wind[x]),zipped_mps[x+1].indices()[0].m());
    //                    L[wind[x+1]]=Index(nameint("L",wind[x+1]),zipped_mps[x+1].indices()[0].m());
    //                    R[wind[x]]=L[wind[x+1]];
    //                }
    //                zipped_mps[0]=change_index_name2(zipped_mps[0], R[wind[0]], zipped_mps[0].indices()[1]);
    //                for (int x=1; x<size-1; x++) {
    //                    zipped_mps[x]=change_index_name3(zipped_mps[x], L[wind[x]], R[wind[x]], zipped_mps[x].indices()[2]);
    //                }
    //                zipped_mps[size-1]=change_index_name2(zipped_mps[size-1], L[wind[size-1]], zipped_mps[size-1].indices()[1]);
    //
    //                for (int x=0; x<size; x++) {
    //                    zipped_mps[x]=sort_indices_braket(zipped_mps[x]);
    //                }
    //
    //                //cout<<"<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
    //                //cout<<"sweep"<<sw<<endl;
    //                //for (int x=0; x<size; x++) {
    //                //    cout<<"x="<<x<<endl;
    //                //    Print(zipped_mps[x]);
    //                //}
    //
    //
    //            }
    //
    //
    //
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<"x="<<x<<endl;
    //            //    Print(zipped_mps[x]);
    //            //PrintDat(zipped_mps[x]);
    //            //Print(mpo[x]);
    //            //Print(original_mps[x]);
    //            //}
    //
    //
    //            return zipped_mps;
    //        }
    //        //-------end: fitting algorithm-----------------------------------------------
    //
    //
    //        //-------begin: fitting algorithm---------------------------------------------
    //        std::vector<ITensor> fitting_algorithmb(std::vector<ITensor> zipped_mps, std::vector<ITensor> mpo, std::vector<ITensor> original_mps, int dcut, int i){
    //
    //            int zsize=zipped_mps.size();
    //            int msize=mpo.size();
    //            int osize=original_mps.size();
    //            int size;
    //            if (zsize==msize&&msize==osize) {
    //                size=zsize;
    //            }
    //            else {
    //                cout<<"ERROR: Sizes of zipped mps, mpo and the original mps does not match!!! (in func: fitting_algorithm)"<<endl;
    //                exit(0);
    //            }
    //
    //            int wind[ny];  //B index and W index
    //            for (int y=0; y<ny; y++) {
    //                wind[y]=i+nx*y;
    //            }
    //
    //            std::vector<ITensor> RRR(size-2), LLL(size-2);
    //            ITensor total;
    //
    //            ITensor SSS,VVV,DDD;
    //
    //
    //            ITensor psi_a,psi_b,k_aa,k_bb,k_ab;
    //
    //            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //            if (size<8) {
    //                cout<<"Sweep step 0:    ";
    //                psi_a=mpo[0]*original_mps[0];
    //                psi_b=zipped_mps[0];
    //                //Print(psi_a);
    //                //Print(psi_b);
    //                for (int y=1; y<size; y++) {
    //                    psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                    psi_b=psi_b*zipped_mps[y];
    //                }
    //                //Print(psi_a.norm());
    //                //Print(psi_a);
    //                //Print(psi_b);
    //                //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                cout<<endl;
    //                //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //            }
    //            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //            //cout<<"sweep0"<<endl;
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<"x="<<x<<endl;
    //            //    Print(zipped_mps[x]);
    //            //}
    //
    //
    //            for (int sw=1; sw<=5; sw++) {
    //
    //                //~~~~~~~~cout<<"Forward scanning......"<<endl;
    //                //forward sweep
    //                RRR[0]=zipped_mps[size-1]*mpo[size-1]*original_mps[size-1];
    //
    //                for (int k=1; k<size-2; k++) {
    //                    RRR[k]=RRR[k-1]*(zipped_mps[size-1-k]*mpo[size-1-k]*original_mps[size-1-k]);
    //                }
    //
    //                total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);
    //
    //                SSS=ITensor(total.indices()[2]);
    //                svd(total, SSS,VVV,DDD);
    //                zipped_mps[0]=SSS;
    //                zipped_mps[1]=VVV*DDD;
    //
    //                LLL[0]=zipped_mps[0]*mpo[0]*original_mps[0];
    //
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                if (size<8) {
    //                    cout<<"Sweep step "<<1+(2*size-3)*(sw-1)<<":    ";
    //                    psi_a=mpo[0]*original_mps[0];
    //                    psi_b=zipped_mps[0];
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    for (int y=1; y<size; y++) {
    //                        psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                        psi_b=psi_b*zipped_mps[y];
    //                    }
    //                    //Print(psi_a.norm());
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                    cout<<endl;
    //                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                }
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //                for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
    //                    total=RRR[size-3-k]*(mpo[1+k]*original_mps[1+k])*(mpo[0+k]*original_mps[0+k])*LLL[k-1];
    //                    SSS=ITensor(total.indices()[2], total.indices()[3]);
    //                    svd(total,SSS,VVV,DDD);
    //                    zipped_mps[0+k]=SSS;
    //                    zipped_mps[1+k]=VVV*DDD;
    //                    LLL[0+k]=LLL[0+k-1]*(zipped_mps[0+k]*mpo[0+k]*original_mps[0+k]);
    //
    //                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                    if (size<8) {
    //                        cout<<"Sweep step "<<k+1+(2*size-3)*(sw-1)<<":    ";
    //                        psi_a=mpo[0]*original_mps[0];
    //                        psi_b=zipped_mps[0];
    //                        //Print(psi_a);
    //                        //Print(psi_b);
    //                        for (int y=1; y<size; y++) {
    //                            psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                            psi_b=psi_b*zipped_mps[y];
    //                        }
    //                        //Print(psi_a.norm());
    //                        //Print(psi_a);
    //                        //Print(psi_b);
    //                        //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                        printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                        cout<<endl;
    //                        //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                    }
    //                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //                }
    //
    //                total=LLL[size-3]*(mpo[size-2]*original_mps[size-2])*(mpo[size-1]*original_mps[size-1]);
    //                SSS=ITensor(total.indices()[0],total.indices()[1]);
    //                svd(total,SSS,VVV,DDD);
    //                zipped_mps[size-2]=SSS*VVV;
    //                zipped_mps[size-1]=DDD;
    //
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                if (size<8) {
    //                    cout<<"Sweep step "<<(size-1)+(2*size-3)*(sw-1)<<":    ";
    //                    psi_a=mpo[0]*original_mps[0];
    //                    psi_b=zipped_mps[0];
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    for (int y=1; y<size; y++) {
    //                        psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                        psi_b=psi_b*zipped_mps[y];
    //                    }
    //                    //Print(psi_a.norm());
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                    cout<<endl;
    //                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                }
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //                //cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //                //for (int x=0; x<size; x++) {
    //                //    cout<<"x="<<x<<endl;
    //                //    Print(zipped_mps[x]);
    //                //}
    //                //~~~~~~~~~~~~cout<<"Backward scanning......"<<endl;
    //                //backward sweep
    //                for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
    //                    if (k==1) {
    //                        RRR[k-1]=zipped_mps[size-k]*mpo[size-k]*original_mps[size-k];
    //                    }
    //                    else {
    //                        RRR[k-1]=RRR[k-2]*(zipped_mps[size-k]*mpo[size-k]*original_mps[size-k]);
    //                    }
    //                    total=RRR[k-1]*(mpo[size-1-k]*original_mps[size-1-k])*(mpo[size-2-k]*original_mps[size-2-k])*LLL[size-3-k];
    //                    SSS=ITensor(total.indices()[2],total.indices()[3]);
    //                    svd(total,SSS,VVV,DDD);
    //                    zipped_mps[size-2-k]=SSS*VVV;
    //                    zipped_mps[size-1-k]=DDD;
    //
    //                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                    if (size<8) {
    //                        cout<<"Sweep step "<<(size-1+k)+(2*size-3)*(sw-1)<<":    ";
    //                        psi_a=mpo[0]*original_mps[0];
    //                        psi_b=zipped_mps[0];
    //                        //Print(psi_a);
    //                        //Print(psi_b);
    //                        for (int y=1; y<size; y++) {
    //                            psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                            psi_b=psi_b*zipped_mps[y];
    //                        }
    //                        //Print(psi_a.norm());
    //                        //Print(psi_a);
    //                        //Print(psi_b);
    //                        //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                        printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                        cout<<endl;
    //                        //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                    }
    //                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //                }
    //
    //                if (size>3) {
    //                    RRR[size-3]=RRR[size-4]*(zipped_mps[2]*mpo[2]*original_mps[2]);
    //                }
    //                else if (size==3){
    //                    RRR[size-3]=(zipped_mps[2]*mpo[2]*original_mps[2]);
    //                }
    //
    //                total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);
    //                SSS=ITensor(total.indices()[2]);
    //                svd(total,SSS,VVV,DDD);
    //                zipped_mps[0]=SSS*VVV;
    //                zipped_mps[1]=DDD;
    //
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                if (size<8) {
    //                    cout<<"Sweep step "<<(2*size-3)+(2*size-3)*(sw-1)<<":    ";
    //                    psi_a=mpo[0]*original_mps[0];
    //                    psi_b=zipped_mps[0];
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    for (int y=1; y<size; y++) {
    //                        psi_a=psi_a*(mpo[y]*original_mps[y]);
    //                        psi_b=psi_b*zipped_mps[y];
    //                    }
    //                    //Print(psi_a.norm());
    //                    //Print(psi_a);
    //                    //Print(psi_b);
    //                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
    //                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
    //                    cout<<endl;
    //                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
    //                }
    //                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //                //zipped_mps[0]=change_index_sequence2(zipped_mps[0]);
    //
    //                //                        cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //                //                        for (int x=0; x<size; x++) {
    //                //                            cout<<"y="<<x<<endl;
    //                //                            Print(zipped_mps[x]);
    //                //                        }
    //
    //
    //
    //
    //                //rename index names makes it easy to debug
    //                for (int y=0; y<size-1; y++) {
    //                    D[wind[y]]=Index(nameint("D",wind[y]),zipped_mps[y+1].indices()[0].m());
    //                    U[wind[y+1]]=Index(nameint("U",wind[y+1]),zipped_mps[y+1].indices()[0].m());
    //                    D[wind[y]]=U[wind[y+1]];
    //                }
    //                zipped_mps[0]=change_index_name2(zipped_mps[0], D[wind[0]], zipped_mps[0].indices()[1]);
    //                for (int y=1; y<size-1; y++) {
    //                    zipped_mps[y]=change_index_name3(zipped_mps[y], U[wind[y]], D[wind[y]], zipped_mps[y].indices()[2]);
    //                }
    //                zipped_mps[size-1]=change_index_name2(zipped_mps[size-1], U[wind[size-1]], zipped_mps[size-1].indices()[1]);
    //
    //                //            for (int x=0; x<size; x++) {
    //                //                zipped_mps[x]=sort_indices_braket(zipped_mps[x]);
    //                //            }
    //
    //                //cout<<"<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
    //                //cout<<"sweep"<<sw<<endl;
    //                //            for (int y=0; y<size; y++) {
    //                //                cout<<"y="<<y<<endl;
    //                //                Print(zipped_mps[y]);
    //                //            }
    //
    //
    //            }
    //
    //
    //
    //            //for (int x=0; x<size; x++) {
    //            //    cout<<"x="<<x<<endl;
    //            //    Print(zipped_mps[x]);
    //            //PrintDat(zipped_mps[x]);
    //            //Print(mpo[x]);
    //            //Print(original_mps[x]);
    //            //}
    //
    //
    //            return zipped_mps;
    //        }
    //        //-------end: fitting algorithm-----------------------------------------------
    //
    //
    //
    //        ITensor reform2(ITensor R){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            std::vector<Index> i(2);
    //            for (int k=0; k<R.r(); k++) {
    //                i[k]=R.indices()[k];
    //            }
    //            gamma=Index("gamma",NB*NB);
    //            gammap=Index("gammap",NB*NB);
    //            tempreal=ITensor(gamma,gammap);
    //            tempimag=ITensor(gamma,gammap);
    //
    //            for (int i0=1; i0<=NB; i0++) {
    //                for (int i1=1; i1<=NB; i1++) {
    //                    for (int i2=1; i2<=NB; i2++) {
    //                        for (int i3=1; i3<=NB; i3++) {
    //                            tempreal(gamma((i0-1)*NB+i2),gammap((i1-1)*NB+i3))=realPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3));
    //                            tempimag(gamma((i0-1)*NB+i2),gammap((i1-1)*NB+i3))=imagPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3));
    //                        }
    //                    }
    //                }
    //            }
    //
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor reform3(ITensor R){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            std::vector<Index> i(3);
    //            for (int k=0; k<R.r(); k++) {
    //                i[k]=R.indices()[k];
    //            }
    //            gamma=Index("gamma",NB*NB*NB);
    //            gammap=Index("gammap",NB*NB*NB);
    //            tempreal=ITensor(gamma,gammap);
    //            tempimag=ITensor(gamma,gammap);
    //            for (int i0=1; i0<=NB; i0++) {
    //                for (int i1=1; i1<=NB; i1++) {
    //                    for (int i2=1; i2<=NB; i2++) {
    //                        for (int i3=1; i3<=NB; i3++) {
    //                            for (int i4=1; i4<=NB; i4++) {
    //                                for (int i5=1; i5<=NB; i5++) {
    //                                    tempreal(gamma((i0-1)*NB*NB+(i2-1)*NB+i4),gammap((i1-1)*NB*NB+(i3-1)*NB+i5))=realPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5));
    //                                    tempimag(gamma((i0-1)*NB*NB+(i2-1)*NB+i4),gammap((i1-1)*NB*NB+(i3-1)*NB+i5))=imagPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5));
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        ITensor reform4(ITensor R){
    //            ITensor tempreal;
    //            ITensor tempimag;
    //            std::vector<Index> i(4);
    //            for (int k=0; k<R.r(); k++) {
    //                i[k]=R.indices()[k];
    //            }
    //            gamma=Index("gamma",NB*NB*NB*NB);
    //            gammap=Index("gammap",NB*NB*NB*NB);
    //            tempreal=ITensor(gamma,gammap);
    //            tempimag=ITensor(gamma,gammap);
    //            for (int i0=1; i0<=NB; i0++) {
    //                for (int i1=1; i1<=NB; i1++) {
    //                    for (int i2=1; i2<=NB; i2++) {
    //                        for (int i3=1; i3<=NB; i3++) {
    //                            for (int i4=1; i4<=NB; i4++) {
    //                                for (int i5=1; i5<=NB; i5++) {
    //                                    for (int i6=1; i6<=NB; i6++) {
    //                                        for (int i7=1; i7<=NB; i7++) {
    //                                            tempreal(gamma((i0-1)*NB*NB*NB+(i2-1)*NB*NB+(i4-1)*NB+i6),gammap((i1-1)*NB*NB*NB+(i3-1)*NB*NB+(i5-1)*NB+i7))=realPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5),i[3]((i6-1)*NB+i7));
    //                                            tempimag(gamma((i0-1)*NB*NB*NB+(i2-1)*NB*NB+(i4-1)*NB+i6),gammap((i1-1)*NB*NB*NB+(i3-1)*NB*NB+(i5-1)*NB+i7))=imagPart(R)(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5),i[3]((i6-1)*NB+i7));
    //                                        }
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //
    //            return tempreal+Complex_i*tempimag;
    //        }
    //
    //        //------begin: finding out the F tensors, to be used to optimize a(X,Y)-------
    //        ITensor sweepx(std::vector<ITensor> bra, std::vector<ITensor> ket, int X, int Y, int rank, int i){
    //            if (bra.size()!=ket.size()) {
    //                cout<<"ERROR: Cannot find F tensors, sizes of bra and ket do not match!!! (in func: sweepx)"<<endl;
    //            }
    //            ITensor L,R,temp;
    //            ITensor abuddy; //abuddy indicates the indice-grouped matrice abuddy_gamma2=a^s_udlr
    //            if (rank==0) {
    //                if (Y==0) {
    //                    if (X==0) { //rank==root, X==0, Y==0
    //                        R=bra[nx-1]*IW[nx-1+nx*Y];
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*bra[x]*IW[x+nx*Y];
    //                        }
    //                        R=R*bra[0];
    //                        R=reform2(R);
    //                        R=R*Identity[X+nx*Y];
    //                        temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta2, R.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
    //                    }
    //                    else if (X==nx-1){  //rank==root, X==nx-1, Y==0
    //                        L=bra[0]*IW[0+nx*Y];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*bra[x]*IW[x+nx*Y];
    //                        }
    //                        L=L*bra[nx-1];
    //                        //Print(L);
    //                        L=reform2(L);
    //                        L=L*Identity[X+nx*Y];
    //                        temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta2, L.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
    //                    }
    //                    else {  //rank==root, 0<X<nx-1, Y==0
    //                        L=bra[0]*IW[0+nx*Y];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*bra[x]*IW[x+nx*Y];
    //                        }
    //                        R=bra[nx-1]*IW[nx-1+nx*Y];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*bra[x]*IW[x+nx*Y];
    //                        }
    //                        temp=L*bra[X]*R;
    //                        //Print(temp);
    //                        temp=reform3(temp);
    //                        temp=temp*Identity[X+nx*Y];
    //                        //Print(temp);
    //                        temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta3, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
    //                        //Print(temp);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //                else if (Y==ny-1){
    //                    if (X==0) { //rank==root, X==0, Y==ny-1
    //                        R=ket[nx-1]*IW[nx-1+nx*Y];
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*ket[x]*IW[x+nx*Y];
    //                        }
    //                        R=R*ket[0];
    //                        R=reform2(R);
    //                        R=R*Identity[X+nx*Y];
    //                        temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta2, R.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
    //                    }
    //                    else if (X==nx-1){  //rank==root, X==nx-1, Y==ny-1
    //                        L=ket[0]*IW[0+nx*Y];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*ket[x]*IW[x+nx*Y];
    //                        }
    //                        L=L*ket[nx-1];
    //                        L=reform2(L);
    //                        L=L*Identity[X+nx*Y];
    //                        temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta2, L.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
    //                    }
    //                    else {  //rank==root, 0<X<nx-1, Y==ny-1
    //                        L=ket[0]*IW[0+nx*Y];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*ket[x]*IW[x+nx*Y];
    //                        }
    //                        R=ket[nx-1]*IW[nx-1+nx*Y];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*ket[x]*IW[x+nx*Y];
    //                        }
    //                        temp=L*ket[X]*R;
    //                        //Print(temp);
    //                        temp=reform3(temp);
    //                        temp=temp*Identity[X+nx*Y];
    //                        //Print(temp);
    //                        temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta3, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
    //                        //Print(temp);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //                else {
    //                    if (X==0) { //rank==root, X==0, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        R=bra[nx-1]*IW[nx-1+nx*Y]*ket[nx-1];
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*bra[x]*IW[x+nx*Y]*ket[x];
    //                        }
    //                        R=R*bra[0]*ket[0];
    //                        R=reform3(R);
    //                        R=R*Identity[X+nx*Y];
    //                        temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta3, R.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
    //                        //Print(R);
    //                        //PrintDat(R);
    //                    }
    //                    else if (X==nx-1){  //rank==root, X==nx-1, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        L=bra[0]*IW[0+nx*Y]*ket[0];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*bra[x]*IW[x+nx*Y]*ket[x];
    //                        }
    //                        L=L*bra[nx-1]*ket[nx-1];
    //                        L=reform3(L);
    //                        L=L*Identity[X+nx*Y];
    //                        temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta3, L.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
    //                        //PrintDat(L);
    //                        //PrintDat(L);
    //                    }
    //                    else {  //rank==root, 0<X<nx-1, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NB*NS);
    //                        L=bra[0]*IW[0+nx*Y]*ket[0];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*bra[x]*IW[x+nx*Y]*ket[x];
    //                        }
    //                        R=bra[nx-1]*IW[nx-1+nx*Y]*ket[nx-1];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*bra[x]*IW[x+nx*Y]*ket[x];
    //                        }
    //                        temp=L*bra[X]*ket[X]*R;
    //                        //Print(temp);
    //                        temp=reform4(temp);
    //                        temp=temp*Identity[X+nx*Y];
    //                        temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta4, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta4p);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //            }
    //
    //            else {  //rank!=0
    //
    //                if (Y==0) {
    //                    if (X==0) { //rank!=root, X==0, Y==0
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        R=bra[nx-1]*W[nx-1+nx*Y+N*i];
    //                        //Print(R);
    //
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*bra[x]*W[x+nx*Y+N*i];
    //                        }
    //                        R=R*bra[0];
    //                        //Print(R);
    //                        R=reform2(R);
    //                        //                    cout<<"rank="<<rank<<", i="<<i;//~~~~~~~
    //                        //                    PrintDat(R);//~~~~~~~~~
    //                        R=R*O[X+nx*Y+N*i];
    //                        //cout<<"i="<<i;
    //                        //PrintDat(O[X+nx*Y+N*i]);
    //
    //                        temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
    //
    //                        //Print(R);
    //                        //PrintDat(R);
    //
    //                    }
    //                    else if (X==nx-1){  //rank!=root, X==nx-1, Y==0
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        L=bra[0]*W[0+nx*Y+N*i];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*bra[x]*W[x+nx*Y+N*i];
    //                        }
    //                        L=L*bra[nx-1];
    //                        L=reform2(L);
    //                        L=L*O[X+nx*Y+N*i];
    //                        temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
    //                        //PrintDat(L);
    //                        //PrintDat(L);
    //                    }
    //                    else {  //rank!=root, 0<X<nx-1, Y==0
    //                        //delta=Index("delta",NB*NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NB*NS);
    //                        L=bra[0]*W[0+nx*Y+N*i];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*bra[x]*W[x+nx*Y+N*i];
    //                        }
    //                        R=bra[nx-1]*W[nx-1+nx*Y+N*i];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*bra[x]*W[x+nx*Y+N*i];
    //                        }
    //                        temp=L*bra[X]*R;
    //                        //Print(temp);
    //                        temp=reform3(temp);
    //                        temp=temp*O[X+nx*Y+N*i];
    //                        temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //
    //                else if (Y==ny-1){
    //                    if (X==0) { //rank!=root, X==0, Y==ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //
    //                        R=ket[nx-1]*W[nx-1+nx*Y+N*i];
    //                        //Print(R);
    //
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*ket[x]*W[x+nx*Y+N*i];
    //                        }
    //
    //                        R=R*ket[0];
    //
    //                        //Print(R);
    //                        R=reform2(R);
    //                        R=R*O[X+nx*Y+N*i];
    //                        temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
    //                        //Print(R);
    //                        //PrintDat(R);
    //
    //
    //                    }
    //
    //                    else if (X==nx-1){  //rank!=root, X==nx-1, Y==ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        L=ket[0]*W[0+nx*Y+N*i];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*ket[x]*W[x+nx*Y+N*i];
    //                        }
    //                        L=L*ket[nx-1];
    //                        L=reform2(L);
    //                        L=L*O[X+nx*Y+N*i];
    //                        temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
    //                        //PrintDat(L);
    //                        //PrintDat(L);
    //                    }
    //                    else {  //rank!=root, 0<X<nx-1, Y==ny-1
    //                        //delta=Index("delta",NB*NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NB*NS);
    //                        L=ket[0]*W[0+nx*Y+N*i];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*ket[x]*W[x+nx*Y+N*i];
    //                        }
    //                        R=ket[nx-1]*W[nx-1+nx*Y+N*i];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*ket[x]*W[x+nx*Y+N*i];
    //                        }
    //                        temp=L*ket[X]*R;
    //                        //Print(temp);
    //                        temp=reform3(temp);
    //                        temp=temp*O[X+nx*Y+N*i];
    //                        temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //                else {
    //                    if (X==0) { //rank!=root, X==0, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        R=bra[nx-1]*W[nx-1+nx*Y+N*i]*ket[nx-1];
    //                        //Print(R);
    //
    //                        for (int x=nx-2; x>0; x--) {
    //                            R=R*bra[x]*W[x+nx*Y+N*i]*ket[x];
    //                        }
    //                        R=R*bra[0]*ket[0];
    //                        //Print(R);
    //                        R=reform3(R);
    //                        R=R*O[X+nx*Y+N*i];
    //                        temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
    //                        //Print(R);
    //                        //PrintDat(R);
    //
    //                    }
    //                    else if (X==nx-1){  //rank!=root, X==nx-1, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NS);
    //                        L=bra[0]*W[0+nx*Y+N*i]*ket[0];
    //                        for (int x=1; x<nx-1; x++) {
    //                            L=L*bra[x]*W[x+nx*Y+N*i]*ket[x];
    //                        }
    //                        L=L*bra[nx-1]*ket[nx-1];
    //                        L=reform3(L);
    //                        L=L*O[X+nx*Y+N*i];
    //                        temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
    //                        //PrintDat(L);
    //                        //PrintDat(L);
    //                    }
    //                    else {  //rank!=root, 0<X<nx-1, 0<Y<ny-1
    //                        //delta=Index("delta",NB*NB*NB*NB*NS);
    //                        //deltap=Index("deltap",NB*NB*NB*NB*NS);
    //                        L=bra[0]*W[0+nx*Y+N*i]*ket[0];
    //                        for (int x=1; x<X; x++) {
    //                            L=L*bra[x]*W[x+nx*Y+N*i]*ket[x];
    //                        }
    //                        R=bra[nx-1]*W[nx-1+nx*Y+N*i]*ket[nx-1];
    //                        for (int x=nx-2; x>X; x--) {
    //                            R=R*bra[x]*W[x+nx*Y+N*i]*ket[x];
    //                        }
    //                        temp=L*bra[X]*ket[X]*R;
    //                        //Print(temp);
    //                        temp=reform4(temp);
    //                        temp=temp*O[X+nx*Y+N*i];
    //                        temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta4, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta4p);
    //                        //PrintDat(temp);
    //                    }
    //                }
    //                //------
    //
    //            }
    //
    //            return temp;
    //        }
    //        //------end: finding out the F tensors, to be used to optimize a(X,Y)-------
    //
    //
    //        //------begin:
    //        ITensor mapping_vector_to_tensor(Vector vector, ITensor tensor, int x, int y){
    //            if (y==0||y==ny-1) {
    //                if (x==0||x==nx-1) {    //x==0,y==0 or x==nx-1,y==0 or x==0,y==ny-1 or x==nx-1,y==ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int is=1; is<=NS; is++) {
    //                                tensor(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector((((i1-1)*NB+i2)-1)*NS+is);
    //                            }
    //                        }
    //                    }
    //                }
    //                else {  //0<x<nx-1,y==0 or 0<x<nx-1,y==ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //            else {
    //                if (x==0) { //x==0, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //                else if (x==nx-1){  //x==nx-1, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //                else {  //0<x<nx-1, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int i4=1; i4<=NB; i4++) {
    //                                    for (int is=1; is<=NS; is++) {
    //                                        tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //
    //            return tensor;
    //        }
    //        //------end:
    //
    //        //------begin:
    //        ITensor mapping_complex_vector_to_tensor(Vector vector_real, Vector vector_imag, ITensor tensor, int x, int y){
    //
    //            ITensor tensor_real,tensor_imag;
    //            tensor_real=realPart(tensor);
    //            tensor_imag=imagPart(tensor);
    //
    //            if (y==0||y==ny-1) {
    //                if (x==0||x==nx-1) {    //x==0,y==0 or x==nx-1,y==0 or x==0,y==ny-1 or x==nx-1,y==ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int is=1; is<=NS; is++) {
    //                                tensor_real(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector_real((((i1-1)*NB+i2)-1)*NS+is);
    //                                tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector_imag((((i1-1)*NB+i2)-1)*NS+is);
    //                            }
    //                        }
    //                    }
    //                }
    //                else {  //0<x<nx-1,y==0 or 0<x<nx-1,y==ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor_real(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                    tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //            else {
    //                if (x==0) { //x==0, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                    tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //                else if (x==nx-1){  //x==nx-1, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int is=1; is<=NS; is++) {
    //                                    tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                    tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //                else {  //0<x<nx-1, 0<y<ny-1
    //                    for (int i1=1; i1<=NB; i1++) {
    //                        for (int i2=1; i2<=NB; i2++) {
    //                            for (int i3=1; i3<=NB; i3++) {
    //                                for (int i4=1; i4<=NB; i4++) {
    //                                    for (int is=1; is<=NS; is++) {
    //                                        tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector_real((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
    //                                        tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector_imag((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //            tensor=tensor_real+Complex_i*tensor_imag;
    //
    //            return tensor;
    //        }
    //        //------end:
    //
    //
    //        std::vector<ITensor> move_otho_next_bra(std::vector<ITensor> tensors, int otho_center, int y){  //y<ny-1
    //            int size=tensors.size();
    //            ITensor S,V,D;
    //
    //            if (otho_center==0) {
    //                D=ITensor(tensors[otho_center].indices()[1]);
    //                svd(tensors[otho_center], S, V, D);    //svd
    //
    //                R[otho_center+nx*(y+1)]=Index(nameint("L",otho_center+1+nx*(y+1)),S.indices()[1].m());
    //                tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*(y+1)]);
    //                L[otho_center+1+nx*(y+1)]=R[otho_center+nx*(y+1)];
    //
    //                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y+1)], tensors[otho_center+1].indices()[2]);
    //            }
    //            else if(otho_center<size-1){
    //                D=ITensor(tensors[otho_center].indices()[2]);
    //                svd(tensors[otho_center], S, V, D);    //svd
    //
    //                R[otho_center+nx*(y+1)]=Index(nameint("L",otho_center+1+nx*(y+1)),S.indices()[2].m());
    //                tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*(y+1)]);
    //                L[otho_center+1+nx*(y+1)]=R[otho_center+nx*(y+1)];
    //
    //                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y+1)], tensors[otho_center+1].indices()[2]);
    //            }
    //
    //            return tensors;
    //        }
    //
    //        std::vector<ITensor> move_otho_next_ket(std::vector<ITensor> tensors, int otho_center, int y){  //y<0
    //            int size=tensors.size();
    //            ITensor S,V,D;
    //
    //            if (otho_center==0) {
    //                D=ITensor(tensors[otho_center].indices()[1]);
    //                svd(tensors[otho_center], S, V, D);    //svd
    //
    //                R[otho_center+nx*(y-1)]=Index(nameint("L",otho_center+1+nx*(y-1)),S.indices()[1].m());
    //                tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*(y-1)]);
    //                L[otho_center+1+nx*(y-1)]=R[otho_center+nx*(y-1)];
    //
    //                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y-1)], tensors[otho_center+1].indices()[2]);
    //            }
    //            else if(otho_center<size-1){
    //                D=ITensor(tensors[otho_center].indices()[2]);
    //                svd(tensors[otho_center], S, V, D);    //svd
    //
    //                R[otho_center+nx*(y-1)]=Index(nameint("L",otho_center+1+nx*(y-1)),S.indices()[2].m());
    //                tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*(y-1)]);
    //                L[otho_center+1+nx*(y-1)]=R[otho_center+nx*(y-1)];
    //
    //                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y-1)], tensors[otho_center+1].indices()[2]);
    //            }
    //
    //            return tensors;
    //        }
    //
    //        std::vector<ITensor> move_otho_next_mpo(std::vector<ITensor> tensors, int otho_center, int y){
    //            int size=tensors.size();
    //            ITensor S,V,D;
    //
    //            if (y==0||y==ny-1) {
    //                if (otho_center==0) {
    //                    D=ITensor(tensors[otho_center].indices()[1]);
    //                    svd(tensors[otho_center], S, V, D);    //svd
    //
    //                    R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[1].m());
    //                    tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*y]);
    //                    L[otho_center+1+nx*y]=R[otho_center+nx*y];
    //
    //                    tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                    tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                    tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[2]);
    //                }
    //                else if(otho_center<size-1){
    //                    D=ITensor(tensors[otho_center].indices()[2]);
    //                    svd(tensors[otho_center], S, V, D);    //svd
    //
    //                    R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[2].m());
    //                    tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*y]);
    //                    L[otho_center+1+nx*y]=R[otho_center+nx*y];
    //
    //                    tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                    tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
    //                    tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[2]);
    //                }
    //            }
    //            else {
    //                if (otho_center==0) {
    //                    D=ITensor(tensors[otho_center].indices()[2]);
    //                    svd(tensors[otho_center], S, V, D);    //svd
    //
    //                    R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[2].m());
    //                    tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*y]);
    //                    L[otho_center+1+nx*y]=R[otho_center+nx*y];
    //
    //                    tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                    tensors[otho_center+1]=sort_indices_otho_4(tensors[otho_center+1], 0, 1, 3, 2);
    //                    tensors[otho_center+1]=change_index_name4(tensors[otho_center+1], tensors[otho_center+1].indices()[0], tensors[otho_center+1].indices()[1], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[3]);
    //                }
    //                else if(otho_center<size-1){
    //                    D=ITensor(tensors[otho_center].indices()[3]);
    //                    svd(tensors[otho_center], S, V, D);    //svd
    //
    //                    R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[3].m());
    //                    tensors[otho_center]=change_index_name4(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], tensors[otho_center].indices()[2], R[otho_center+nx*y]);
    //                    L[otho_center+1+nx*y]=R[otho_center+nx*y];
    //
    //                    tensors[otho_center+1]=tensors[otho_center+1]*V*D;
    //                    tensors[otho_center+1]=sort_indices_otho_4(tensors[otho_center+1], 0, 1, 3, 2);
    //                    tensors[otho_center+1]=change_index_name4(tensors[otho_center+1], tensors[otho_center+1].indices()[0], tensors[otho_center+1].indices()[1], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[3]);
    //                }
    //            }
    //
    //            return tensors;
    //        }
    //
    //    }
    
    
    
    
    //    if (COMPLEX==0){
    //-------begin: resize tensors----------------------------------------------------
    ITensor resize2(ITensor N,Index b,int dcut){  //change the first index size to be dcut and rename it into b
        return N;
    }
    //-------begin: resize tensors----------------------------------------------------
    
    //------
    ITensor sort_indices_braket(ITensor Z){
        ITensor temp;
        if (Z.r()==2) {
            temp=ITensor(Z.indices()[1],Z.indices()[0]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    temp(temp.indices()[0](i1),temp.indices()[1](i0))=Z(Z.indices()[0](i0),Z.indices()[1](i1));
                }
            }
        }
        else if (Z.r()==3){
            temp=ITensor(Z.indices()[2],Z.indices()[0],Z.indices()[1]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
                        temp(temp.indices()[0](i2),temp.indices()[1](i0),temp.indices()[2](i1))=Z(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
                    }
                }
            }
        }
        
        return temp;
    }
    
    ITensor sort_indices_compressor(ITensor Z){
        ITensor temp;
        if (Z.r()==2) {
            temp=ITensor(Z.indices()[1],Z.indices()[0]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    temp(temp.indices()[0](i1),temp.indices()[1](i0))=Z(Z.indices()[0](i0),Z.indices()[1](i1));
                }
            }
        }
        else if (Z.r()==3){
            temp=ITensor(Z.indices()[1],Z.indices()[0],Z.indices()[2]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
                        temp(temp.indices()[0](i1),temp.indices()[1](i0),temp.indices()[2](i2))=Z(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
                    }
                }
            }
        }
        
        return temp;
    }
    
    ITensor sort_indices_mpo(ITensor Z){
        ITensor temp;
        if (Z.r()==3) {
            temp=ITensor(Z.indices()[1],Z.indices()[2],Z.indices()[0]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
                        temp(temp.indices()[0](i1),temp.indices()[1](i2),temp.indices()[2](i0))=Z(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2));
                    }
                }
            }
        }
        else if (Z.r()==4){
            temp=ITensor(Z.indices()[2],Z.indices()[3],Z.indices()[0],Z.indices()[1]);
            for (int i0=1; i0<=Z.indices()[0].m(); i0++) {
                for (int i1=1; i1<=Z.indices()[1].m(); i1++) {
                    for (int i2=1; i2<=Z.indices()[2].m(); i2++) {
                        for (int i3=1; i3<=Z.indices()[3].m(); i3++) {
                            temp(temp.indices()[0](i2),temp.indices()[1](i3),temp.indices()[2](i0),temp.indices()[3](i1))=Z(Z.indices()[0](i0),Z.indices()[1](i1),Z.indices()[2](i2),Z.indices()[3](i3));
                        }
                    }
                }
            }
        }
        
        return temp;
    }
    
    ITensor sort_indices_otho_2(ITensor Z, int a0, int a1){
        if (Z.r()!=2) {
            cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_2)"<<endl;
            exit(0);
        }
        ITensor temp(Z.indices()[a0],Z.indices()[a1]);
        for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
            for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
                temp(temp.indices()[0](i0),temp.indices()[1](i1))=Z(Z.indices()[a0](i0),Z.indices()[a1](i1));
            }
        }
        
        return temp;
    }
    ITensor sort_indices_otho_3(ITensor Z, int a0, int a1, int a2){
        if (Z.r()!=3) {
            cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_3)"<<endl;
            exit(0);
        }
        ITensor temp(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2]);
        for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
            for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
                for (int i2=1; i2<=Z.indices()[a2].m(); i2++) {
                    temp(temp.indices()[0](i0),temp.indices()[1](i1),temp.indices()[2](i2))=Z(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2));
                }
            }
        }
        
        return temp;
    }
    ITensor sort_indices_otho_4(ITensor Z, int a0, int a1, int a2, int a3){
        if (Z.r()!=4) {
            cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_3)"<<endl;
            exit(0);
        }
        ITensor temp(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2], Z.indices()[a3]);
        for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
            for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
                for (int i2=1; i2<=Z.indices()[a2].m(); i2++) {
                    for (int i3=1; i3<=Z.indices()[a3].m(); i3++) {
                        temp(temp.indices()[0](i0),temp.indices()[1](i1),temp.indices()[2](i2),temp.indices()[3](i3))=Z(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2),Z.indices()[a3](i3));
                    }
                }
            }
        }
        
        return temp;
    }
    ITensor sort_indices_otho_5(ITensor Z, int a0, int a1, int a2, int a3, int a4){
        if (Z.r()!=5) {
            cout<<"Error: Number of indices does not match!!! (In func: sort_indices_otho_3)"<<endl;
            exit(0);
        }
        ITensor temp(Z.indices()[a0], Z.indices()[a1], Z.indices()[a2], Z.indices()[a3], Z.indices()[a4]);
        for (int i0=1; i0<=Z.indices()[a0].m(); i0++) {
            for (int i1=1; i1<=Z.indices()[a1].m(); i1++) {
                for (int i2=1; i2<=Z.indices()[a2].m(); i2++) {
                    for (int i3=1; i3<=Z.indices()[a3].m(); i3++) {
                        for (int i4=1; i4<=Z.indices()[a4].m(); i4++) {
                            temp(temp.indices()[0](i0),temp.indices()[1](i1),temp.indices()[2](i2),temp.indices()[3](i3),temp.indices()[4](i4))=Z(Z.indices()[a0](i0),Z.indices()[a1](i1),Z.indices()[a2](i2),Z.indices()[a3](i3),Z.indices()[a4](i4));
                        }
                    }
                }
            }
        }
        
        return temp;
    }
    //------
    
    //-----
    ITensor change_index_sequence2(ITensor N){
        if (N.r()!=2) {
            cout<<"Error: Number of indices in target tensor in not 2!!! (In func: change_index_sequence2) "<<endl;
            exit(0);
        }
        ITensor temp(N.indices()[1],N.indices()[0]);
        for (int i1=1; i1<=N.indices()[1].m(); i1++) {
            for (int i0=1; i0<=N.indices()[0].m(); i0++) {
                temp(temp.indices()[0](i1),temp.indices()[1](i0))=N(N.indices()[0](i0),N.indices()[1](i1));
            }
        }
        
        return temp;
    }
    //-----
    
    //-------begin: change names of indice--------------------------------------------
    ITensor change_index_name2(ITensor A, Index a1, Index a2){
        ITensor temp=ITensor(a1, a2);
        
        //cout<<"indicater~~~"<<endl;
        
        if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m()) {
            cout<<"Error: New indice do not match old indice!!! (in func: change_index_name2)"<<endl;
            //Print(A);
            //Print(a1);
            //Print(a2);
            //            cout<<"~~~~~~~~~~"<<endl;
            exit(0);
        }
        for (int i1=1; i1<=a1.m(); i1++) {
            for (int i2=1; i2<=a2.m(); i2++) {
                temp(a1(i1),a2(i2))=A(A.indices()[0](i1),A.indices()[1](i2));
            }
        }
        return temp;
    }
    
    ITensor change_index_name3(ITensor A, Index a1, Index a2, Index a3){
        ITensor temp=ITensor(a1, a2, a3);
        if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m() || A.indices()[2].m()!=a3.m()) {
            cout<<"Error: New indice do not match old indice!!! (in func: change_index_name3)"<<endl;
            exit(0);
        }
        for (int i1=1; i1<=a1.m(); i1++) {
            for (int i2=1; i2<=a2.m(); i2++) {
                for (int i3=1; i3<=a3.m(); i3++) {
                    temp(a1(i1),a2(i2),a3(i3))=A(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3));
                }
            }
        }
        return temp;
    }
    
    ITensor change_index_name4(ITensor A, Index a1, Index a2, Index a3, Index a4){
        ITensor temp=ITensor(a1, a2, a3, a4);
        if (A.indices()[0].m()!=a1.m() || A.indices()[1].m()!=a2.m() || A.indices()[2].m()!=a3.m() || A.indices()[3].m()!=a4.m()) {
            cout<<"Error: New indice do not match old indice!!! (in func: change_index_name4)"<<endl;
            exit(0);
        }
        for (int i1=1; i1<=a1.m(); i1++) {
            for (int i2=1; i2<=a2.m(); i2++) {
                for (int i3=1; i3<=a3.m(); i3++) {
                    for (int i4=1; i4<=a4.m(); i4++) {
                        temp(a1(i1),a2(i2),a3(i3),a4(i4))=A(A.indices()[0](i1),A.indices()[1](i2),A.indices()[2](i3),A.indices()[3](i4));
                    }
                }
            }
        }
        return temp;
    }
    //-------end: change names of indice--------------------------------------------
    
    //-------begin: tensor indice combiners-------------------------------------------
    ITensor group_indice2(ITensor Z, Index a1, Index a2, Index a12){
        ITensor temp=ITensor(a12);
        int ma1=a1.m(),ma2=a2.m(),ma12=a12.m();
        if (ma1*ma2!=ma12) {
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice2)"<<endl;
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
    
    ITensor group_indice3(ITensor Z, Index a1, Index a2, Index a12, Index a3){  //group index a1 and a2, leave a3 inact: Z_a1,a2,a3=M_(a1,a2),a3
        ITensor temp=ITensor(a12,a3);
        int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(),ma3=a3.m();
        if (ma1*ma2!=ma12) {
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice3)"<<endl;
            exit(0);
        }
        
        for (int i1=1; i1<=ma1; i1++) {
            for (int i2=1; i2<=ma2; i2++) {
                int i12=(i1-1)*ma2+i2;
                for (int i3=1; i3<=ma3; i3++) {
                    temp(a12(i12),a3(i3))=Z(a1(i1),a2(i2),a3(i3));
                }
            }
        }
        
        return temp;
    }
    
    ITensor group_indice4(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34){
        ITensor temp=ITensor(a12,a34);
        int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(), ma3=a3.m(),ma4=a4.m(),ma34=a34.m();
        if (ma1*ma2!=ma12||ma3*ma4!=ma34) {
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice4)"<<endl;
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
    
    ITensor group_indice5(ITensor Z, Index a1, Index a2, Index a12, Index a3, Index a4, Index a34, Index a5){  //group index a1 and a2, a3 and a4, leave a5 inact: Z_a1,a2,a3,a4,a5=M_(a1,a2),(a3,a4),a5
        ITensor temp=ITensor(a12,a34,a5);
        int ma1=a1.m(),ma2=a2.m(),ma12=a12.m(),ma3=a3.m(),ma4=a4.m(),ma34=a34.m(),ma5=a5.m();
        if (ma1*ma2!=ma12 || ma3*ma4!=ma34) {
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice5)"<<endl;
            exit(0);
        }
        
        for (int i1=1; i1<=ma1; i1++) {
            for (int i2=1; i2<=ma2; i2++) {
                int i12=(i1-1)*ma2+i2;
                for (int i3=1; i3<=ma3; i3++) {
                    for (int i4=1; i4<=ma4; i4++) {
                        int i34=(i3-1)*ma4+i4;
                        
                        for (int i5=1; i5<=ma5; i5++) {
                            temp(a12(i12),a34(i34),a5(i5))=Z(a1(i1),a2(i2),a3(i3),a4(i4),a5(i5));
                        }
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
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice6)"<<endl;
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
            cout<<"Error: Tensor sizes do not match!!! (in func: group_indice8)"<<endl;
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
    //-------end: tensor indice combiners-------------------------------------------
    
    //-------begin: move othogonality center to the first site for bra, without cutoff----------------------
    std::vector<ITensor> bra_otho_center_1(std::vector<ITensor> Z, int y, int nx, int ny){
        
        //        if (backflag==1) {
        //            for (int x=0; x<NY; x++) {
        //                cout<<"x="<<x;
        //                PrintDat(Z[x]);
        //            }
        //        }
        
        int size=Z.size();
        int u_ind=-1,l_ind=-1,r_ind=-1;
        ITensor S, V, D;
        int aaa;
        
        std::vector<ITensor> original(size);
        for (int x=0; x<size; x++) {
            original[x]=Z[x];
        }
        
        D=ITensor(original[size-1].indices()[0]);
        svd(Z[size-1], S, V, D);    //svd
        
        //find out indice names and sequence
        for (int i=0; i<2; i++) {
            if (D.indices()[i]==original[size-1].indices()[0]) {
                u_ind=i;
            }
        }
        l_ind=1-u_ind;
        //        if (backflag==1) {
        //            Print(Z[size-1]);
        //        }
        
        //if (backflag==0) {//~~~~~
        L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
        //if (backflag==0) {//~~~~~~~~~
        //Z=D, and rename indices as U, L, R indice
        if (u_ind==0) {
            Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        }
        if (u_ind==1) {
            Z[size-1]=change_index_name2(D, L[size-1+nx*y], original[size-1].indices()[0]);
        }
        
        Z[size-1]=sort_indices_otho_2(Z[size-1], u_ind, l_ind);
        
        //        D=sort_indices_otho_2(D, u_ind, l_ind);
        //        Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        
        for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
            ///change index name for Z[x]//////
            Z[x]=Z[x]*S*V;
            //if (x==size-3) {
            //    PrintDat(Z[x]);
            //}
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (Z[x].indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==original[x].indices()[1]) {
                    l_ind=i;
                    aaa+=i;
                }
            }
            r_ind=3-aaa;
            R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
            
            
            //PrintDat(Z[x]);
            //rename indices of Z as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                }
            }
            
            //PrintDat(Z[x]);
            S=ITensor(original[x].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
            svd(Z[x], S, V, D); //svd
            //PrintDat(S);
            //PrintDat(V);
            //PrintDat(D);
            aaa=0;
            for (int i=0; i<3; i++) {
                if (D.indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (D.indices()[i]==R[x+nx*y]) {
                    r_ind=i;
                    aaa+=i;
                }
            }
            l_ind=3-aaa;
            
            L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
            
            //Z=D, and rename indices as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
                }
            }
            
            Z[x]=sort_indices_otho_3(Z[x], u_ind, l_ind, r_ind);
            
            //            D=sort_indices_otho_3(D, u_ind, l_ind, r_ind);
            //            Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
            //PrintDat(D);
            //PrintDat(Z[x]);
            
        }
        
        Z[0]=Z[0]*S*V;
        for (int i=0; i<2; i++) {
            if (Z[0].indices()[i]==original[0].indices()[0]) {
                u_ind=i;
            }
        }
        r_ind=1-u_ind;
        R[0+nx*y]=L[1+nx*y];
        if (u_ind==0) {
            Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        }
        if (u_ind==1) {
            Z[0]=change_index_name2(Z[0], R[0+nx*y], original[0].indices()[0]);
        }
        
        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        
        //        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        //        Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        //for (int x=0; x<size; x++) {
        //    cout<<x<<"~~~~~~~~~"<<endl;
        //    PrintDat(Z[x]);
        //}
        //        for (int x=0; x<size; x++) {
        //            Z[x]=sort_indices_braket(Z[x]);
        //        }
        //}//~~~~~~~~
        
        return Z;
    }
    //-------end: move othogonality center to the first site----------------------
    
    
    
    
    //-------begin: move othogonality center to the first site for bra, without cutoff----------------------
    std::vector<ITensor> bbra_otho_center_1(std::vector<ITensor> Z, int y, int nx, int ny){
        
        //        if (backflag==1) {
        //            for (int x=0; x<NY; x++) {
        //                cout<<"x="<<x;
        //                PrintDat(Z[x]);
        //            }
        //        }
        
        int size=Z.size();
        int u_ind=-1,l_ind=-1,r_ind=-1;
        ITensor S, V, D;
        int aaa;
        
        std::vector<ITensor> original(size);
        for (int x=0; x<size; x++) {
            original[x]=Z[x];
            //            cout<<"x="<<x;
            //            Print(original[x]);
        }
        
        D=ITensor(original[size-1].indices()[0]);
        svd(Z[size-1], S, V, D);    //svd
        
        //find out indice names and sequence
        for (int i=0; i<2; i++) {
            if (D.indices()[i]==original[size-1].indices()[0]) {
                u_ind=i;
            }
        }
        l_ind=1-u_ind;
        //        if (backflag==1) {
        //            Print(Z[size-1]);
        //        }
        for (int x=0; x<size; x++) {
            cout<<x<<"~~~~~~~~~"<<endl;
            PrintDat(Z[x]);
        }
        //if (backflag==0) {//~~~~~
        L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
        //if (backflag==0) {//~~~~~~~~~
        //Z=D, and rename indices as U, L, R indice
        if (u_ind==0) {
            Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        }
        if (u_ind==1) {
            Z[size-1]=change_index_name2(D, L[size-1+nx*y], original[size-1].indices()[0]);
        }
        
        Z[size-1]=sort_indices_otho_2(Z[size-1], u_ind, l_ind);
        
        //        D=sort_indices_otho_2(D, u_ind, l_ind);
        //        Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        
        for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
            ///change index name for Z[x]//////
            Z[x]=Z[x]*S*V;
            //if (x==size-3) {
            //    PrintDat(Z[x]);
            //}
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (Z[x].indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==original[x].indices()[1]) {
                    l_ind=i;
                    aaa+=i;
                }
            }
            r_ind=3-aaa;
            R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
            
            
            //PrintDat(Z[x]);
            //rename indices of Z as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                }
            }
            
            //PrintDat(Z[x]);
            S=ITensor(original[x].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
            svd(Z[x], S, V, D); //svd
            //PrintDat(S);
            //PrintDat(V);
            //PrintDat(D);
            aaa=0;
            for (int i=0; i<3; i++) {
                if (D.indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (D.indices()[i]==R[x+nx*y]) {
                    r_ind=i;
                    aaa+=i;
                }
            }
            l_ind=3-aaa;
            
            L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
            
            //Z=D, and rename indices as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
                }
            }
            
            Z[x]=sort_indices_otho_3(Z[x], u_ind, l_ind, r_ind);
            
            //            D=sort_indices_otho_3(D, u_ind, l_ind, r_ind);
            //            Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
            //PrintDat(D);
            //PrintDat(Z[x]);
            
        }
        
        Z[0]=Z[0]*S*V;
        for (int i=0; i<2; i++) {
            if (Z[0].indices()[i]==original[0].indices()[0]) {
                u_ind=i;
            }
        }
        r_ind=1-u_ind;
        R[0+nx*y]=L[1+nx*y];
        if (u_ind==0) {
            Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        }
        if (u_ind==1) {
            Z[0]=change_index_name2(Z[0], R[0+nx*y], original[0].indices()[0]);
        }
        
        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        
        //        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        //        Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        //        for (int x=0; x<size; x++) {
        //            cout<<x<<"~~~~~~~~~"<<endl;
        //            PrintDat(Z[x]);
        //        }
        //        for (int x=0; x<size; x++) {
        //            Z[x]=sort_indices_braket(Z[x]);
        //        }
        //}//~~~~~~~~
        
        return Z;
    }
    //-------end: move othogonality center to the first site----------------------
    
    
    
    
    
    
    //-------begin: move othogonality center to the first site for bra, with dimension cutoff m----------------------
    std::vector<ITensor> bra_otho_center_1m(std::vector<ITensor> Z, int y, int m, int nx, int ny){
        int size=Z.size();
        int u_ind=-1,l_ind=-1,r_ind=-1;
        ITensor S, V, D;
        int aaa;
        
        D=ITensor(U[size-1+nx*y]);
        
        OptSet opts;
        opts.add("Maxm", m);
        svd(Z[size-1], S, V, D, opts);    //svd
        
        //find out indice names and sequence
        for (int i=0; i<2; i++) {
            if (D.indices()[i]==U[size-1+nx*y]) {
                u_ind=i;
            }
        }
        l_ind=1-u_ind;
        
        L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
        
        //Z=D, and rename indices as U, L, R indice
        if (u_ind==0) {
            Z[size-1]=change_index_name2(D, U[size-1+nx*y], L[size-1+nx*y]);
        }
        if (u_ind==1) {
            Z[size-1]=change_index_name2(D, L[size-1+nx*y], U[size-1+nx*y]);
        }
        
        for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
            ///change index name for Z[x]//////
            Z[x]=Z[x]*S*V;
            
            //cout<<"x="<<x<<endl;
            //if (x==size-3) {
            //    Print(Z[x]);
            //}
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (Z[x].indices()[i]==U[x+nx*y]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==L[x+nx*y]) {
                    l_ind=i;
                    aaa+=i;
                }
            }
            r_ind=3-aaa;
            R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
            
            
            //PrintDat(Z[x]);
            //rename indices of Z as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(Z[x], U[x+nx*y], L[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], U[x+nx*y], R[x+nx*y], L[x+nx*y]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], L[x+nx*y], U[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], U[x+nx*y], L[x+nx*y]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], L[x+nx*y], R[x+nx*y], U[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], L[x+nx*y], U[x+nx*y]);
                }
            }
            
            //PrintDat(Z[x]);
            S=ITensor(L[x+nx*y]);   //let L index be in S tensor, U and R be in D tensor
            svd(Z[x], S, V, D, opts); //svd
            //if (x==2) {
            //    Print(Z[x]);
            //    Print(S);
            //    Print(V);
            //    Print(D);
            //}
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (D.indices()[i]==U[x+nx*y]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (D.indices()[i]==R[x+nx*y]) {
                    r_ind=i;
                    aaa+=i;
                }
            }
            l_ind=3-aaa;
            
            L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
            
            //Z=D, and rename indices as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(D, U[x+nx*y], L[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, U[x+nx*y], R[x+nx*y], L[x+nx*y]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], U[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], U[x+nx*y], L[x+nx*y]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], U[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], U[x+nx*y]);
                }
            }
            //PrintDat(D);
            //PrintDat(Z[x]);
            
        }
        
        Z[0]=Z[0]*S*V;
        for (int i=0; i<2; i++) {
            if (Z[0].indices()[i]==U[0+nx*y]) {
                u_ind=i;
            }
        }
        r_ind=1-u_ind;
        R[0+nx*y]=L[1+nx*y];
        if (u_ind==0) {
            Z[0]=change_index_name2(Z[0], U[0+nx*y], R[0+nx*y]);
        }
        if (u_ind==1) {
            Z[0]=change_index_name2(Z[0], R[0+nx*y], U[0+nx*y]);
        }
        
        //for (int x=0; x<size; x++) {
        //    cout<<x<<"~~~~~~~~~"<<endl;
        //    PrintDat(Z[x]);
        //}
        
        return Z;
    }
    //-------end: move othogonality center to the first site----------------------
    
    //-------begin: move othogonality center to the first site for ket, without cutoff----------------------
    std::vector<ITensor> ket_otho_center_1(std::vector<ITensor> Z, int y, int nx, int ny){
        int size=Z.size();
        int u_ind=-1,l_ind=-1,r_ind=-1;
        ITensor S, V, D;
        int aaa;
        
        
        std::vector<ITensor> original(size);
        for (int x=0; x<size; x++) {
            original[x]=Z[x];
        }
        
        D=ITensor(original[size-1].indices()[0]);   //for ket: U[size-1+nx+nx*y]=D[size-1+nx*y]
        svd(Z[size-1], S, V, D);    //svd
        
        //find out indice names and sequence
        for (int i=0; i<2; i++) {
            if (D.indices()[i]==original[size-1].indices()[0]) {
                u_ind=i;
            }
        }
        l_ind=1-u_ind;
        
        L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),D.indices()[l_ind].m());  //resize index, to be used for the updated Z
        
        if (u_ind==0) {
            //cout<<"3568~~~~~~~~~~"<<endl;
            Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        }
        if (u_ind==1) {
            //cout<<"3572~~~~~~~~~~"<<endl;
            Z[size-1]=change_index_name2(D, L[size-1+nx*y], original[size-1].indices()[0]);
        }
        
        Z[size-1]=sort_indices_otho_2(Z[size-1], u_ind, l_ind);
        
        
        //        D=sort_indices_otho_2(D, u_ind, l_ind);
        //        Z[size-1]=change_index_name2(D, original[size-1].indices()[0], L[size-1+nx*y]);
        
        
        for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
            ///change index name for Z[x]//////
            Z[x]=Z[x]*S*V;
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (Z[x].indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==original[x].indices()[1]) {
                    l_ind=i;
                    aaa+=i;
                }
            }
            r_ind=3-aaa;
            R[x+nx*y]=L[x+1+nx*y];  //let R indices match neighbour L indices
            
            
            //rename indices of Z as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                }
            }
            
            S=ITensor(original[x].indices()[1]);   //let L index be in S tensor, U and R be in D tensor
            svd(Z[x], S, V, D); //svd
            
            aaa=0;
            for (int i=0; i<3; i++) {
                if (D.indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (D.indices()[i]==R[x+nx*y]) {
                    r_ind=i;
                    aaa+=i;
                }
            }
            l_ind=3-aaa;
            
            L[x+nx*y]=Index(nameint("L",x+nx*y),D.indices()[l_ind].m());    //resize index, to be used for the updated Z
            
            
            //Z=D, and rename indices as U,L,R
            if (u_ind==0) {
                if (l_ind==1) {
                    Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
                }
            }
            else if (u_ind==1) {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
                }
            }
            else {
                if (l_ind==0) {
                    Z[x]=change_index_name3(D, L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
                }
                else {
                    Z[x]=change_index_name3(D, R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
                }
            }
            
            Z[x]=sort_indices_otho_3(Z[x], u_ind, l_ind, r_ind);
            
            
            //            D=sort_indices_otho_3(D, u_ind, l_ind, r_ind);
            //            Z[x]=change_index_name3(D, original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
            
            
            
        }
        
        
        Z[0]=Z[0]*S*V;
        
        for (int i=0; i<2; i++) {
            if (Z[0].indices()[i]==original[0].indices()[0]) {
                u_ind=i;
            }
        }
        r_ind=1-u_ind;
        R[0+nx*y]=L[1+nx*y];
        
        
        if (u_ind==0) {
            //cout<<"3696~~~~~~~~~~"<<endl;
            Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        }
        if (u_ind==1) {
            //cout<<"3700~~~~~~~~~~"<<endl;
            Z[0]=change_index_name2(Z[0], R[0+nx*y], original[0].indices()[0]);
        }
        
        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        
        
        //        Z[0]=sort_indices_otho_2(Z[0], u_ind, r_ind);
        //        Z[0]=change_index_name2(Z[0], original[0].indices()[0], R[0+nx*y]);
        
        
        /*
         for (int x=0; x<nx; x++) {
         cout<<"x="<<x;
         Print(Z[x]);
         }
         */
        
        return Z;
    }
    //-------end: move othogonality center to the first site----------------------
    
    
    
    
    //-------begin: move othogonality center to the first site for mpo----------------------
    std::vector<ITensor> mpo_otho_center_1(std::vector<ITensor> Z, int y, int nx, int ny){
        int size=Z.size();
        int u_ind=-1,d_ind=-1,l_ind=-1,r_ind=-1;
        ITensor S, V, DD;
        
        std::vector<ITensor> original(size);
        for (int x=0; x<size; x++) {
            original[x]=Z[x];
        }
        
        S=ITensor(original[size-1].indices()[2]);
        svd(Z[size-1], S, V, DD);
        int aaa=0;
        for (int i=0; i<3; i++) {
            if (DD.indices()[i]==original[size-1].indices()[0]) {
                u_ind=i;
                aaa+=i;
            }
            if (DD.indices()[i]==original[size-1].indices()[1]) {
                d_ind=i;
                aaa+=i;
            }
        }
        l_ind=3-aaa;
        
        L[size-1+nx*y]=Index(nameint("L",size-1+nx*y),DD.indices()[l_ind].m());  //resize index, to be used for the updated Z
        
        if (u_ind==0) {
            if (d_ind==1) {
                Z[size-1]=change_index_name3(DD, original[size-1].indices()[0], original[size-1].indices()[1], L[size-1+nx*y]);
            }
            else {
                Z[size-1]=change_index_name3(DD, original[size-1].indices()[0], L[size-1+nx*y], original[size-1].indices()[1]);
            }
        }
        else if (u_ind==1) {
            if (d_ind==0) {
                Z[size-1]=change_index_name3(DD, original[size-1].indices()[1], original[size-1].indices()[0], L[size-1+nx*y]);
            }
            else {
                Z[size-1]=change_index_name3(DD, L[size-1+nx*y], original[size-1].indices()[0], original[size-1].indices()[1]);
            }
        }
        else {
            if (d_ind==0) {
                Z[size-1]=change_index_name3(DD, original[size-1].indices()[1], L[size-1+nx*y], original[size-1].indices()[0]);
            }
            else {
                Z[size-1]=change_index_name3(DD, L[size-1+nx*y], original[size-1].indices()[1], original[size-1].indices()[0]);
            }
        }
        //PrintDat(S);
        //PrintDat(V);
        
        
        for (int x=size-2; x>0; x--) {  //int x=size-2; x>0; x--
            ///change index name for Z[x]//////
            Z[x]=Z[x]*S*V;
            //PrintDat(Z[x]);
            
            //if (x==size-3) {
            //    PrintDat(Z[x]);
            //}
            
            
            int aaa=0;
            for (int i=0; i<4; i++) {
                if (Z[x].indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==original[x].indices()[1]) {
                    d_ind=i;
                    aaa+=i;
                }
                if (Z[x].indices()[i]==original[x].indices()[2]) {
                    l_ind=i;
                    aaa+=i;
                }
            }
            r_ind=6-aaa;
            
            R[x+nx*y]=L[x+1+nx*y];
            //PrintDat(Z[x]);
            
            ///rename indices of Z[x] to the wanted format
            if (u_ind==0) {
                if (d_ind==1) {
                    if (l_ind==2) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[1], original[x].indices()[2], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[1], R[x+nx*y], original[x].indices()[2]);
                    }
                }
                else if (d_ind==2) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[2], original[x].indices()[1], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[1], original[x].indices()[2]);
                    }
                }
                else {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], original[x].indices()[2], R[x+nx*y], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[0], R[x+nx*y], original[x].indices()[2], original[x].indices()[1]);
                    }
                }
            }
            else if (u_ind==1) {
                if (d_ind==0) {
                    if (l_ind==2) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[0], original[x].indices()[2], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[0], R[x+nx*y], original[x].indices()[2]);
                    }
                }
                else if (d_ind==2) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[1], original[x].indices()[2]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[0], original[x].indices()[2], original[x].indices()[1]);
                    }
                }
            }
            else if (u_ind==2) {
                if (d_ind==0) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[2], original[x].indices()[0], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[0], original[x].indices()[2]);
                    }
                }
                else if (d_ind==1) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[0], original[x].indices()[2]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[2], original[x].indices()[0], original[x].indices()[1]);
                    }
                }
            }
            else {
                if (d_ind==0) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], original[x].indices()[2], R[x+nx*y], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[1], R[x+nx*y], original[x].indices()[2], original[x].indices()[0]);
                    }
                }
                else if (d_ind==1) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[1], original[x].indices()[2], original[x].indices()[0]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(Z[x], original[x].indices()[2], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(Z[x], R[x+nx*y], original[x].indices()[2], original[x].indices()[1], original[x].indices()[0]);
                    }
                }
            }
            
            //PrintDat(Z[x]);
            
            S=ITensor(original[x].indices()[2]);
            svd(Z[x], S, V, DD);
            //PrintDat(S);
            //PrintDat(V);
            //PrintDat(DD);
            
            aaa=0;
            for (int i=0; i<4; i++) {
                if (DD.indices()[i]==original[x].indices()[0]) {
                    u_ind=i;
                    aaa+=i;
                }
                if (DD.indices()[i]==original[x].indices()[1]) {
                    d_ind=i;
                    aaa+=i;
                }
                if (DD.indices()[i]==R[x+nx*y]) {
                    r_ind=i;
                    aaa+=i;
                }
            }
            l_ind=6-aaa;
            
            L[x+nx*y]=Index(nameint("L",x+nx*y),DD.indices()[l_ind].m());    //resize index, to be used for the updated Z
            
            if (u_ind==0) {
                if (d_ind==1) {
                    if (l_ind==2) {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], original[x].indices()[1], L[x+nx*y], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], original[x].indices()[1], R[x+nx*y], L[x+nx*y]);
                    }
                }
                else if (d_ind==2) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], L[x+nx*y], original[x].indices()[1], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], R[x+nx*y], original[x].indices()[1], L[x+nx*y]);
                    }
                }
                else {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], L[x+nx*y], R[x+nx*y], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[0], R[x+nx*y], L[x+nx*y], original[x].indices()[1]);
                    }
                }
            }
            else if (u_ind==1) {
                if (d_ind==0) {
                    if (l_ind==2) {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], original[x].indices()[0], L[x+nx*y], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], original[x].indices()[0], R[x+nx*y], L[x+nx*y]);
                    }
                }
                else if (d_ind==2) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[0], original[x].indices()[1], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[0], original[x].indices()[1], L[x+nx*y]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[0], R[x+nx*y], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[0], L[x+nx*y], original[x].indices()[1]);
                    }
                }
            }
            else if (u_ind==2) {
                if (d_ind==0) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], L[x+nx*y], original[x].indices()[0], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], R[x+nx*y], original[x].indices()[0], L[x+nx*y]);
                    }
                }
                else if (d_ind==1) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[1], original[x].indices()[0], R[x+nx*y]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[1], original[x].indices()[0], L[x+nx*y]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], R[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], L[x+nx*y], original[x].indices()[0], original[x].indices()[1]);
                    }
                }
            }
            else {
                if (d_ind==0) {
                    if (l_ind==1) {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], L[x+nx*y], R[x+nx*y], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, original[x].indices()[1], R[x+nx*y], L[x+nx*y], original[x].indices()[0]);
                    }
                }
                else if (d_ind==1) {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], original[x].indices()[1], R[x+nx*y], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], original[x].indices()[1], L[x+nx*y], original[x].indices()[0]);
                    }
                }
                else {
                    if (l_ind==0) {
                        Z[x]=change_index_name4(DD, L[x+nx*y], R[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                    }
                    else {
                        Z[x]=change_index_name4(DD, R[x+nx*y], L[x+nx*y], original[x].indices()[1], original[x].indices()[0]);
                    }
                }
            }
            
            //PrintDat(DD);
            //PrintDat(Z[x]);
            
            
        }
        
        Z[0]=Z[0]*S*V;
        aaa=0;
        for (int i=0; i<3; i++) {
            if (Z[0].indices()[i]==original[0].indices()[0]) {
                u_ind=i;
                aaa+=i;
            }
            if (Z[0].indices()[i]==original[0].indices()[1]) {
                d_ind=i;
                aaa+=i;
            }
        }
        r_ind=3-aaa;
        
        R[0+nx*y]=L[1+nx*y];
        
        if (u_ind==0) {
            if (d_ind==1) {
                Z[0]=change_index_name3(Z[0], original[0].indices()[0], original[0].indices()[1], R[0+nx*y]);
            }
            else {
                Z[0]=change_index_name3(Z[0], original[0].indices()[0], R[0+nx*y], original[0].indices()[1]);
            }
        }
        else if (u_ind==1) {
            if (d_ind==0) {
                Z[0]=change_index_name3(Z[0], original[0].indices()[1], original[0].indices()[0], R[0+nx*y]);
            }
            else {
                Z[0]=change_index_name3(Z[0], R[0+nx*y], original[0].indices()[0], original[0].indices()[1]);
            }
        }
        else {
            if (d_ind==0) {
                Z[0]=change_index_name3(Z[0], original[0].indices()[1], R[0+nx*y], original[0].indices()[0]);
            }
            else {
                Z[0]=change_index_name3(Z[0], R[0+nx*y], original[0].indices()[1], original[0].indices()[0]);
            }
        }
        
        for (int x=0; x<size; x++) {
            Z[x]=sort_indices_mpo(Z[x]);
        }
        
        return Z;
    }
    //-------end: move othogonality center to the first site----------------------
    
    
    //-------begin: fitting algorithm---------------------------------------------
    std::vector<ITensor> fitting_algorithm(std::vector<ITensor> zipped_mps, std::vector<ITensor> mpo, std::vector<ITensor> original_mps, int dcut, int i, int nx, int ny){
        
        int zsize=zipped_mps.size();
        int msize=mpo.size();
        int osize=original_mps.size();
        int size;
        if (zsize==msize&&msize==osize) {
            size=zsize;
        }
        else {
            cout<<"ERROR: Sizes of zipped mps, mpo and the original mps does not match!!! (in func: fitting_algorithm)"<<endl;
            exit(0);
        }
        
//        if (i==2) {
//            for (int x=0; x<size; x++) {
//                PrintDat(zipped_mps[x]);
//            }
//        }
        
        
        int index[size],oindex[size];
        for (int x=0; x<size; x++) {
            index[x]=0;
            oindex[x]=0;
        }
        

        int lflag=0,rflag=0;
        if (original_mps[0].indices()[0].m()==1) {
            oindex[0]=1;
        }
        

        for (int x=1; x<size-1; x++) {
            if (oindex[x-1]==1) {
                if (original_mps[x].indices()[1].m()==1&&original_mps[x].indices()[2].m()==1) {
                    oindex[x]=1;
                }
            }
            else {
                if (original_mps[x].indices()[2].m()==1) {
                    oindex[x]=1;
                }
            }
        }
        
        
        int wind[nx];  //B index and W index
        for (int x=0; x<nx; x++) {
            wind[x]=x+nx*(ny-1-i);
        }
        

        std::vector<ITensor> RRR(size-2), LLL(size-2);
        ITensor total;
        
        ITensor SSS,VVV,DDD;
        
        ITensor psi_a,psi_b,k_aa,k_bb,k_ab;
        
        if (fitting_process_print==1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (size<8) {
                cout<<"Sweep step 0:    ";
                psi_a=mpo[0]*original_mps[0];
                psi_b=zipped_mps[0];
                //Print(psi_a);
                //Print(psi_b);
                for (int x=1; x<size; x++) {
                    psi_a=psi_a*(mpo[x]*original_mps[x]);
                    psi_b=psi_b*zipped_mps[x];
                }
                psi_a=psi_a/psi_a.norm();
                psi_b=psi_b/psi_b.norm();
                //Print(psi_a.norm());
                //Print(psi_a);
                //Print(psi_b);
                //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                cout<<endl;
                //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
            }
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        }
        

        for (int sw=1; sw<=nfittingsweep; sw++) {
            
            //~~~~~~~~cout<<"Forward scanning......"<<endl;
            //forward sweep
            RRR[0]=zipped_mps[size-1]*mpo[size-1]*original_mps[size-1];
            
            //cout<<"4173~~~~~~~~~"<<endl;
            
//            cout<<"4231~~~~~~"<<endl;
            
//            cout<<"4239<><><><><><><><"<<endl;
            
            //Print(zipped_mps[2]*mpo[2]*original_mps[2]);
//            PrintDat(zipped_mps[2]);
//            PrintDat(mpo[2]);
//            PrintDat(original_mps[2]);

//            cout<<"4242<><><><><><><><"<<endl;
            for (int k=1; k<size-2; k++) {
//                RRR[k]=RRR[k-1]*(zipped_mps[size-1-k]*mpo[size-1-k]*original_mps[size-1-k]);
                RRR[k]=RRR[k-1]*zipped_mps[size-1-k]*mpo[size-1-k]*original_mps[size-1-k];
            }
//            cout<<"4237~~~~~~"<<endl;

            //cout<<"4177~~~~~~~~~"<<endl;
            total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);
            

            if (total.indices()[2].m()==1) {
                SSS=ITensor(total.indices()[1]);
            }
            else {
                SSS=ITensor(total.indices()[2]);
            }
            

            //Print(SSS);//~~~~
            
            //cout<<"4179~~~~~~~~~"<<endl;
            //SSS=ITensor(total.indices()[2]);
            //SSS=ITensor(total.indices()[2]);
            svd(total, SSS,VVV,DDD);
            

            //            Print(SSS);
            //            Print(VVV*DDD);
            //cout<<"4182~~~~~~~~~"<<endl;
            zipped_mps[0]=SSS;
            zipped_mps[1]=VVV*DDD;
            //cout<<"4185~~~~~~~~~"<<endl;
            LLL[0]=zipped_mps[0]*mpo[0]*original_mps[0];
            

            //cout<<"4187~~~~~~~~~"<<endl;
            if (fitting_process_print==1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (size<8) {
                    cout<<"Sweep step "<<1+(2*size-3)*(sw-1)<<":    ";
                    psi_a=mpo[0]*original_mps[0];
                    psi_b=zipped_mps[0];
                    //Print(psi_a);
                    //Print(psi_b);
                    for (int x=1; x<size; x++) {
                        psi_a=psi_a*(mpo[x]*original_mps[x]);
                        psi_b=psi_b*zipped_mps[x];
                    }
                    psi_a=psi_a/psi_a.norm();
                    psi_b=psi_b/psi_b.norm();
                    //Print(psi_a.norm());
                    //Print(psi_a);
                    //Print(psi_b);
                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                    //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                    cout<<endl;
                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
                }
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            }
            

            for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
                //cout<<"4211~~~~~~~~~"<<endl;
                total=RRR[size-3-k]*(mpo[1+k]*original_mps[1+k])*(mpo[0+k]*original_mps[0+k])*LLL[k-1];
                //cout<<"4213~~~~~~~~~"<<endl;
                
                lflag=(oindex[k-1]==1&&LLL[k-1].indices()[1].m()==1&&LLL[k-1].indices()[2].m()==1)||(oindex[k-1]==0&&LLL[k-1].indices()[2].m()==1);
                rflag=(oindex[k+1]==1&&RRR[size-3-k].indices()[1].m()==1&&RRR[size-3-k].indices()[2].m()==1)||(oindex[k+1]==0&&RRR[size-3-k].indices()[2].m()==1);
                
                if (lflag!=1&&rflag!=1) {
                    SSS=ITensor(total.indices()[2], total.indices()[3]);
                }
                else if(lflag!=1&&rflag==1){
                    SSS=ITensor(total.indices()[1], total.indices()[2]);
                }
                else if(lflag==1&&rflag!=1){
                    SSS=ITensor(total.indices()[2], total.indices()[3]);
                }
                else if(lflag==1&&rflag==1){
                    SSS=ITensor(total.indices()[1], total.indices()[3]);
                }
                //Print(SSS);//~~~~~~
                //SSS=ITensor(total.indices()[2], total.indices()[3]);
                svd(total,SSS,VVV,DDD);
                //                if (SSS.indices()[2].m()==1&&SSS.indices()[1].m()!=1) {
                //                    index[k]=1;
                //                }
                zipped_mps[0+k]=SSS;
                zipped_mps[1+k]=VVV*DDD;
                //cout<<"4218~~~~~~~~~"<<endl;
                LLL[0+k]=LLL[0+k-1]*zipped_mps[0+k]*mpo[0+k]*original_mps[0+k];
                //cout<<"4221~~~~~~~~~"<<endl;
                if (fitting_process_print==1) {
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if (size<8) {
                        cout<<"Sweep step "<<k+1+(2*size-3)*(sw-1)<<":    ";
                        psi_a=mpo[0]*original_mps[0];
                        psi_b=zipped_mps[0];
                        //Print(psi_a);
                        //Print(psi_b);
                        for (int x=1; x<size; x++) {
                            psi_a=psi_a*(mpo[x]*original_mps[x]);
                            psi_b=psi_b*zipped_mps[x];
                        }
                        psi_a=psi_a/psi_a.norm();
                        psi_b=psi_b/psi_b.norm();
                        //Print(psi_a.norm());
                        //Print(psi_a);
                        //Print(psi_b);
                        //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                        printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                        //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                        cout<<endl;
                        //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
                    }
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                }
                
                
            }

            //cout<<"4244~~~~~~~~~"<<endl;
            total=LLL[size-3]*(mpo[size-2]*original_mps[size-2])*(mpo[size-1]*original_mps[size-1]);
            if (total.indices()[2].m()==1) {
                SSS=ITensor(total.indices()[0],total.indices()[2]);
            }
            else {
                SSS=ITensor(total.indices()[0],total.indices()[1]);
            }
            

            //Print(SSS);//~~~~~~
            //SSS=ITensor(total.indices()[0],total.indices()[1]);
            svd(total,SSS,VVV,DDD);
            
            zipped_mps[size-2]=SSS*VVV;
            zipped_mps[size-1]=DDD;
            if (zipped_mps[size-1].indices()[1].m()==1) {
                index[size-1]=1;
            }

            //cout<<"4250~~~~~~~~~"<<endl;
            if (fitting_process_print==1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (size<8) {
                    cout<<"Sweep step "<<(size-1)+(2*size-3)*(sw-1)<<":    ";
                    psi_a=mpo[0]*original_mps[0];
                    psi_b=zipped_mps[0];
                    //Print(psi_a);
                    //Print(psi_b);
                    for (int x=1; x<size; x++) {
                        psi_a=psi_a*(mpo[x]*original_mps[x]);
                        psi_b=psi_b*zipped_mps[x];
                    }
                    psi_a=psi_a/psi_a.norm();
                    psi_b=psi_b/psi_b.norm();
                    //Print(psi_a.norm());
                    //Print(psi_a);
                    //Print(psi_b);
                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                    //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                    cout<<endl;
                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
                }
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            }
            
            

            //cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            //for (int x=0; x<size; x++) {
            //    cout<<"x="<<x<<endl;
            //    Print(zipped_mps[x]);
            //}
            //~~~~~~~~~~~~cout<<"Backward scanning......"<<endl;
            //backward sweep
            for (int k=1; k<size-2; k++) {  //int k=1; k<size-2; k++
                //cout<<"4281~~~~~~~~~"<<endl;
                if (k==1) {
                    RRR[k-1]=zipped_mps[size-k]*mpo[size-k]*original_mps[size-k];
                }
                else {
                    RRR[k-1]=RRR[k-2]*zipped_mps[size-k]*mpo[size-k]*original_mps[size-k];
                }
                total=RRR[k-1]*(mpo[size-1-k]*original_mps[size-1-k])*(mpo[size-2-k]*original_mps[size-2-k])*LLL[size-3-k];
                
                lflag=(oindex[size-3-k]==1&&LLL[size-3-k].indices()[1].m()==1&&LLL[size-3-k].indices()[2].m()==1)||(oindex[size-3-k]==0&&LLL[size-3-k].indices()[2].m()==1);
                rflag=(oindex[size-1-k]==1&&RRR[k-1].indices()[1].m()==1&&RRR[k-1].indices()[2].m()==1)||(oindex[size-1-k]==0&&RRR[k-1].indices()[2].m()==1);
                
                if (lflag!=1&&rflag!=1) {
                    SSS=ITensor(total.indices()[2], total.indices()[3]);
                }
                else if(lflag!=1&&rflag==1){
                    SSS=ITensor(total.indices()[1], total.indices()[2]);
                }
                else if(lflag==1&&rflag!=1){
                    SSS=ITensor(total.indices()[2], total.indices()[3]);
                }
                else if(lflag==1&&rflag==1){
                    SSS=ITensor(total.indices()[1], total.indices()[3]);
                }
                
                //SSS=ITensor(total.indices()[2],total.indices()[3]);
                svd(total,SSS,VVV,DDD);
                zipped_mps[size-2-k]=SSS*VVV;
                zipped_mps[size-1-k]=DDD;
                
                if (lflag!=1&&rflag!=1&&zipped_mps[size-1-k].indices()[2].m()==1) {
                    index[size-1-k]=1;
                }
                else if(lflag!=1&&rflag==1){
                    if (zipped_mps[size-2-k].indices()[2].m()==1) {
                        index[size-1-k]=3;
                    }
                    else {
                        index[size-1-k]=2;
                    }
                }
                else if(lflag==1&&rflag!=1&&zipped_mps[size-1-k].indices()[2].m()==1){
                    index[size-1-k]=1;
                }
                else if(lflag==1&&rflag==1){
                    if (zipped_mps[size-1-k].indices()[2].m()==1&&zipped_mps[size-1-k].indices()[1].m()==1) {
                        index[size-1-k]=3;
                    }
                    else {
                        index[size-1-k]=2;
                    }
                }
                //cout<<"4293~~~~~~~~~"<<endl;
                if (fitting_process_print==1) {
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if (size<8) {
                        cout<<"Sweep step "<<(size-1+k)+(2*size-3)*(sw-1)<<":    ";
                        psi_a=mpo[0]*original_mps[0];
                        psi_b=zipped_mps[0];
                        //Print(psi_a);
                        //Print(psi_b);
                        for (int x=1; x<size; x++) {
                            psi_a=psi_a*(mpo[x]*original_mps[x]);
                            psi_b=psi_b*zipped_mps[x];
                        }
                        psi_a=psi_a/psi_a.norm();
                        psi_b=psi_b/psi_b.norm();
                        //Print(psi_a.norm());
                        //Print(psi_a);
                        //Print(psi_b);
                        //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                        printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                        //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                        cout<<endl;
                        //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
                    }
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                }
                
                
            }

            //cout<<"4317~~~~~~~~~"<<endl;
            if (size>3) {
                RRR[size-3]=RRR[size-4]*zipped_mps[2]*mpo[2]*original_mps[2];
            }
            else if (size==3){
                RRR[size-3]=zipped_mps[2]*mpo[2]*original_mps[2];
            }
            
            total=RRR[size-3]*(mpo[1]*original_mps[1])*(mpo[0]*original_mps[0]);

            
            if (total.indices()[2].m()==1) {
                SSS=ITensor(total.indices()[1]);
            }
            else {
                SSS=ITensor(total.indices()[2]);
            }
            

            svd(total,SSS,VVV,DDD);
            zipped_mps[0]=SSS*VVV;
            zipped_mps[1]=DDD;
            

            if (total.indices()[2].m()==1) {
                if (zipped_mps[1].indices()[2].m()==1&&zipped_mps[1].indices()[1].m()==1) {
                    index[1]=3;
                }
                else {
                    index[1]=2;
                }
            }
            else {
                if (zipped_mps[1].indices()[2].m()==1) {
                    index[1]=1;
                }
            }
            

            //cout<<"4330~~~~~~~~~"<<endl;
            if (fitting_process_print==1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (size<8) {
                    cout<<"Sweep step "<<(2*size-3)+(2*size-3)*(sw-1)<<":    ";
                    psi_a=mpo[0]*original_mps[0];
                    psi_b=zipped_mps[0];
                    //Print(psi_a);
                    //Print(psi_b);
                    for (int x=1; x<size; x++) {
                        psi_a=psi_a*(mpo[x]*original_mps[x]);
                        psi_b=psi_b*zipped_mps[x];
                    }
                    psi_a=psi_a/psi_a.norm();
                    psi_b=psi_b/psi_b.norm();
                    //Print(psi_a.norm());
                    //Print(psi_a);
                    //Print(psi_b);
                    //PrintDat(dag(psi_a-psi_b)*(psi_a-psi_b));
                    printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1)) );
                    //printf("%.10f", (dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))/psi_a.norm() );
                    cout<<endl;
                    //cout<<(dag(psi_a-psi_b)*(psi_a-psi_b))((dag(psi_a-psi_b)*(psi_a-psi_b)).indices()[0](1))<<endl;
                }
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            }
            
            
            //zipped_mps[0]=change_index_sequence2(zipped_mps[0]);
            
            //            cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            //            for (int x=0; x<size; x++) {
            //                cout<<"x="<<x<<endl;
            //                Print(zipped_mps[x]);
            //            }
            
            /*
             if (index[size-1]==0) {
             R[wind[size-2]]=Index(nameint("R",wind[size-2]),zipped_mps[size-1].indices()[0].m());
             L[wind[size-1]]=Index(nameint("L",wind[size-1]),zipped_mps[size-1].indices()[0].m());
             R[wind[size-2]]=L[wind[size-1]];
             zipped_mps[size-1]=change_index_name2(zipped_mps[size-1], L[wind[size-1]], zipped_mps[size-1].indices()[1]);
             }
             else if (index[size-1]==1) {
             R[wind[size-2]]=Index(nameint("R",wind[size-2]),zipped_mps[size-1].indices()[1].m());
             L[wind[size-1]]=Index(nameint("L",wind[size-1]),zipped_mps[size-1].indices()[1].m());
             R[wind[size-2]]=L[wind[size-1]];
             zipped_mps[size-1]=change_index_name2(zipped_mps[size-1], zipped_mps[size-1].indices()[0], L[wind[size-1]]);
             }
             
             
             for (int x=size-2; x>1; x--) {
             if (index[x]==0) {
             R[wind[x-1]]=Index(nameint("R",wind[x-1]),zipped_mps[x].indices()[0].m());
             L[wind[x]]=Index(nameint("L",wind[x]),zipped_mps[x].indices()[0].m());
             R[wind[x-1]]=L[wind[x]];
             zipped_mps[x]=change_index_name3(zipped_mps[x], L[wind[x]], R[wind[x]], zipped_mps[x].indices()[2]);
             zipped_mps[x]=sort_indices_otho_3(zipped_mps[x],1,2,0);
             }
             else if (index[x]==1) {
             R[wind[x-1]]=Index(nameint("R",wind[x-1]),zipped_mps[x].indices()[2].m());
             L[wind[x]]=Index(nameint("L",wind[x]),zipped_mps[x].indices()[2].m());
             R[wind[x-1]]=L[wind[x]];
             zipped_mps[x]=change_index_name3(zipped_mps[x], R[wind[x]], zipped_mps[x].indices()[1], L[wind[x]]);
             zipped_mps[x]=sort_indices_otho_3(zipped_mps[x],1,0,2);
             }
             else if (index[x]==2) {
             R[wind[x-1]]=Index(nameint("R",wind[x-1]),zipped_mps[x].indices()[0].m());
             L[wind[x]]=Index(nameint("L",wind[x]),zipped_mps[x].indices()[0].m());
             R[wind[x-1]]=L[wind[x]];
             zipped_mps[x]=change_index_name3(zipped_mps[x], L[wind[x]], zipped_mps[x].indices()[1], R[wind[x]]);
             zipped_mps[x]=sort_indices_otho_3(zipped_mps[x],1,0,2);
             }
             else if (index[x]==3) {
             R[wind[x-1]]=Index(nameint("R",wind[x-1]),zipped_mps[x].indices()[2].m());
             L[wind[x]]=Index(nameint("L",wind[x]),zipped_mps[x].indices()[2].m());
             R[wind[x-1]]=L[wind[x]];
             zipped_mps[x]=change_index_name3(zipped_mps[x], zipped_mps[x].indices()[0], R[wind[x]], L[wind[x]]);
             }
             }
             
             
             if (index[1]==0) {
             R[wind[0]]=Index(nameint("R",wind[0]),zipped_mps[1].indices()[0].m());
             L[wind[1]]=Index(nameint("L",wind[1]),zipped_mps[1].indices()[0].m());
             R[wind[0]]=L[wind[1]];
             zipped_mps[1]=change_index_name3(zipped_mps[1], L[wind[1]], R[wind[1]], zipped_mps[1].indices()[2]);
             zipped_mps[1]=sort_indices_otho_3(zipped_mps[1],1,2,0);
             }
             else if (index[1]==1) {
             R[wind[0]]=Index(nameint("R",wind[0]),zipped_mps[1].indices()[2].m());
             L[wind[1]]=Index(nameint("L",wind[1]),zipped_mps[1].indices()[2].m());
             R[wind[0]]=L[wind[1]];
             zipped_mps[1]=change_index_name3(zipped_mps[1], R[wind[1]], zipped_mps[1].indices()[1], L[wind[1]]);
             zipped_mps[1]=sort_indices_otho_3(zipped_mps[1],1,0,2);
             }
             else if (index[1]==2) {
             R[wind[0]]=Index(nameint("R",wind[0]),zipped_mps[1].indices()[0].m());
             L[wind[1]]=Index(nameint("L",wind[1]),zipped_mps[1].indices()[0].m());
             R[wind[0]]=L[wind[1]];
             zipped_mps[1]=change_index_name3(zipped_mps[1], L[wind[1]], zipped_mps[1].indices()[1], R[wind[1]]);
             zipped_mps[1]=sort_indices_otho_3(zipped_mps[1],1,0,2);
             }
             else if (index[1]==3) {
             R[wind[0]]=Index(nameint("R",wind[0]),zipped_mps[1].indices()[2].m());
             L[wind[1]]=Index(nameint("L",wind[1]),zipped_mps[1].indices()[2].m());
             R[wind[0]]=L[wind[1]];
             zipped_mps[1]=change_index_name3(zipped_mps[1], zipped_mps[1].indices()[0], R[wind[1]], L[wind[1]]);
             }
             
             zipped_mps[0]=change_index_name2(zipped_mps[0], zipped_mps[0].indices()[0], R[wind[0]]);
             */
            
            
            
            //            for (int x=0; x<size; x++) {
            //                cout<<"x="<<x<<endl;
            //                Print(zipped_mps[x]);
            //            }
            /*
             //rename index names makes it easy to debug
             for (int x=0; x<size-1; x++) {
             R[wind[x]]=Index(nameint("R",wind[x]),zipped_mps[x+1].indices()[0].m());
             L[wind[x+1]]=Index(nameint("L",wind[x+1]),zipped_mps[x+1].indices()[0].m());
             R[wind[x]]=L[wind[x+1]];
             }
             //cout<<"4367~~~~~~~~~"<<endl;////////////////////////////////////////
             zipped_mps[0]=change_index_name2(zipped_mps[0], R[wind[0]], zipped_mps[0].indices()[1]);
             //cout<<"4369~~~~~~~~~"<<endl;////////////////////////////////////////
             for (int x=1; x<size-1; x++) {
             zipped_mps[x]=change_index_name3(zipped_mps[x], L[wind[x]], R[wind[x]], zipped_mps[x].indices()[2]);
             }
             //cout<<"4373~~~~~~~~~"<<endl;
             zipped_mps[size-1]=change_index_name2(zipped_mps[size-1], L[wind[size-1]], zipped_mps[size-1].indices()[1]);
             //cout<<"4375~~~~~~~~~"<<endl;
             for (int x=0; x<size; x++) {
             zipped_mps[x]=sort_indices_braket(zipped_mps[x]);
             }
             //cout<<"4379~~~~~~~~~"<<endl;
             //cout<<"<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
             //cout<<"sweep"<<sw<<endl;
             */
            
        }
        

        
        //for (int x=0; x<size; x++) {
        //    cout<<"x="<<x<<endl;
        //    Print(zipped_mps[x]);
        //PrintDat(zipped_mps[x]);
        //Print(mpo[x]);
        //Print(original_mps[x]);
        //}
        
        
        return zipped_mps;
    }
    //-------end: fitting algorithm-----------------------------------------------
    
    ITensor reform2(ITensor R){
        ITensor temp;
        std::vector<Index> i(2);
        for (int k=0; k<R.r(); k++) {
            i[k]=R.indices()[k];
        }
        gamma=Index("gamma",NB*NB);
        gammap=Index("gammap",NB*NB);
        temp=ITensor(gamma,gammap);
        
        for (int i0=1; i0<=NB; i0++) {
            for (int i1=1; i1<=NB; i1++) {
                for (int i2=1; i2<=NB; i2++) {
                    for (int i3=1; i3<=NB; i3++) {
                        temp(gamma((i0-1)*NB+i2),gammap((i1-1)*NB+i3))=R(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3));
                    }
                }
            }
        }
        
        
        return temp;
    }
    
    ITensor reform3(ITensor R){
        ITensor temp;
        std::vector<Index> i(3);
        for (int k=0; k<R.r(); k++) {
            i[k]=R.indices()[k];
        }
        gamma=Index("gamma",NB*NB*NB);
        gammap=Index("gammap",NB*NB*NB);
        temp=ITensor(gamma,gammap);
        for (int i0=1; i0<=NB; i0++) {
            for (int i1=1; i1<=NB; i1++) {
                for (int i2=1; i2<=NB; i2++) {
                    for (int i3=1; i3<=NB; i3++) {
                        for (int i4=1; i4<=NB; i4++) {
                            for (int i5=1; i5<=NB; i5++) {
                                temp(gamma((i0-1)*NB*NB+(i2-1)*NB+i4),gammap((i1-1)*NB*NB+(i3-1)*NB+i5))=R(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5));
                            }
                        }
                    }
                }
            }
        }
        
        
        return temp;
    }
    
    ITensor reform4(ITensor R){
        ITensor temp;
        std::vector<Index> i(4);
        for (int k=0; k<R.r(); k++) {
            i[k]=R.indices()[k];
        }
        gamma=Index("gamma",NB*NB*NB*NB);
        gammap=Index("gammap",NB*NB*NB*NB);
        temp=ITensor(gamma,gammap);
        for (int i0=1; i0<=NB; i0++) {
            for (int i1=1; i1<=NB; i1++) {
                for (int i2=1; i2<=NB; i2++) {
                    for (int i3=1; i3<=NB; i3++) {
                        for (int i4=1; i4<=NB; i4++) {
                            for (int i5=1; i5<=NB; i5++) {
                                for (int i6=1; i6<=NB; i6++) {
                                    for (int i7=1; i7<=NB; i7++) {
                                        temp(gamma((i0-1)*NB*NB*NB+(i2-1)*NB*NB+(i4-1)*NB+i6),gammap((i1-1)*NB*NB*NB+(i3-1)*NB*NB+(i5-1)*NB+i7))=R(i[0]((i0-1)*NB+i1),i[1]((i2-1)*NB+i3),i[2]((i4-1)*NB+i5),i[3]((i6-1)*NB+i7));
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
    
    //------begin: finding out the F tensors, to be used to optimize a(X,Y)-------
    ITensor sweepx(std::vector<ITensor> bra, std::vector<ITensor> ket, int X, int Y, int rank, int i, int nx, int ny){
        if (bra.size()!=ket.size()) {
            cout<<"ERROR: Cannot find F tensors, sizes of bra and ket do not match!!! (in func: sweepx)"<<endl;
        }
        ITensor L,R,temp;
        ITensor abuddy; //abuddy indicates the indice-grouped matrice abuddy_gamma2=a^s_udlr
        if (rank==0) {
            if (Y==0) {
                if (X==0) { //rank==root, X==0, Y==0
                    R=bra[nx-1]*IW[nx-1+nx*Y];
                    for (int x=nx-2; x>0; x--) {
                        R=R*bra[x]*IW[x+nx*Y];
                    }
                    R=R*bra[0];
                    R=reform2(R);
                    R=R*Identity[X+nx*Y];
                    temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta2, R.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
                }
                else if (X==nx-1){  //rank==root, X==nx-1, Y==0
                    L=bra[0]*IW[0+nx*Y];
                    for (int x=1; x<nx-1; x++) {
                        L=L*bra[x]*IW[x+nx*Y];
                    }
                    L=L*bra[nx-1];
                    //Print(L);
                    L=reform2(L);
                    L=L*Identity[X+nx*Y];
                    temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta2, L.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
                }
                else {  //rank==root, 0<X<nx-1, Y==0
                    L=bra[0]*IW[0+nx*Y];
                    for (int x=1; x<X; x++) {
                        L=L*bra[x]*IW[x+nx*Y];
                    }
                    R=bra[nx-1]*IW[nx-1+nx*Y];
                    for (int x=nx-2; x>X; x--) {
                        R=R*bra[x]*IW[x+nx*Y];
                    }
                    temp=L*bra[X]*R;
                    //Print(temp);
                    temp=reform3(temp);
                    temp=temp*Identity[X+nx*Y];
                    //Print(temp);
                    temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta3, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
                    //Print(temp);
                    //PrintDat(temp);
                }
            }
            else if (Y==ny-1){
                if (X==0) { //rank==root, X==0, Y==ny-1
                    R=ket[nx-1]*IW[nx-1+nx*Y];
                    for (int x=nx-2; x>0; x--) {
                        R=R*ket[x]*IW[x+nx*Y];
                    }
                    R=R*ket[0];
                    R=reform2(R);
                    R=R*Identity[X+nx*Y];
                    temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta2, R.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
                }
                else if (X==nx-1){  //rank==root, X==nx-1, Y==ny-1
                    L=ket[0]*IW[0+nx*Y];
                    for (int x=1; x<nx-1; x++) {
                        L=L*ket[x]*IW[x+nx*Y];
                    }
                    L=L*ket[nx-1];
                    L=reform2(L);
                    L=L*Identity[X+nx*Y];
                    temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta2, L.indices()[1], Identity[X+nx*Y].indices()[1], delta2p);
                }
                else {  //rank==root, 0<X<nx-1, Y==ny-1
                    L=ket[0]*IW[0+nx*Y];
                    for (int x=1; x<X; x++) {
                        L=L*ket[x]*IW[x+nx*Y];
                    }
                    R=ket[nx-1]*IW[nx-1+nx*Y];
                    for (int x=nx-2; x>X; x--) {
                        R=R*ket[x]*IW[x+nx*Y];
                    }
                    temp=L*ket[X]*R;
                    //Print(temp);
                    temp=reform3(temp);
                    temp=temp*Identity[X+nx*Y];
                    //Print(temp);
                    temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta3, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
                    //Print(temp);
                    //PrintDat(temp);
                }
            }
            else {
                if (X==0) { //rank==root, X==0, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    R=bra[nx-1]*IW[nx-1+nx*Y]*ket[nx-1];
                    for (int x=nx-2; x>0; x--) {
                        R=R*bra[x]*IW[x+nx*Y]*ket[x];
                    }
                    R=R*bra[0]*ket[0];
                    R=reform3(R);
                    R=R*Identity[X+nx*Y];
                    temp=group_indice4(R, R.indices()[0], Identity[X+nx*Y].indices()[0], delta3, R.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
                    //Print(R);
                    //PrintDat(R);
                }
                else if (X==nx-1){  //rank==root, X==nx-1, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    L=bra[0]*IW[0+nx*Y]*ket[0];
                    for (int x=1; x<nx-1; x++) {
                        L=L*bra[x]*IW[x+nx*Y]*ket[x];
                    }
                    L=L*bra[nx-1]*ket[nx-1];
                    L=reform3(L);
                    L=L*Identity[X+nx*Y];
                    temp=group_indice4(L, L.indices()[0], Identity[X+nx*Y].indices()[0], delta3, L.indices()[1], Identity[X+nx*Y].indices()[1], delta3p);
                    //PrintDat(L);
                    //PrintDat(L);
                }
                else {  //rank==root, 0<X<nx-1, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NB*NS);
                    L=bra[0]*IW[0+nx*Y]*ket[0];
                    for (int x=1; x<X; x++) {
                        L=L*bra[x]*IW[x+nx*Y]*ket[x];
                    }
                    R=bra[nx-1]*IW[nx-1+nx*Y]*ket[nx-1];
                    for (int x=nx-2; x>X; x--) {
                        R=R*bra[x]*IW[x+nx*Y]*ket[x];
                    }
                    temp=L*bra[X]*ket[X]*R;
                    //Print(temp);
                    temp=reform4(temp);
                    temp=temp*Identity[X+nx*Y];
                    temp=group_indice4(temp, temp.indices()[0], Identity[X+nx*Y].indices()[0], delta4, temp.indices()[1], Identity[X+nx*Y].indices()[1], delta4p);
                    //PrintDat(temp);
                }
            }
        }
        
        else {  //rank!=0
            
            if (Y==0) {
                if (X==0) { //rank!=root, X==0, Y==0
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    R=bra[nx-1]*W[nx-1+nx*Y+N*i];
                    //Print(R);
                    
                    for (int x=nx-2; x>0; x--) {
                        R=R*bra[x]*W[x+nx*Y+N*i];
                    }
                    R=R*bra[0];
                    //Print(R);
                    R=reform2(R);
                    //                    cout<<"rank="<<rank<<", i="<<i;//~~~~~~~
                    //                    PrintDat(R);//~~~~~~~~~
                    R=R*O[X+nx*Y+N*i];
                    //cout<<"i="<<i;
                    //PrintDat(O[X+nx*Y+N*i]);
                    
                    temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
                    
                    //Print(R);
                    //PrintDat(R);
                    
                }
                else if (X==nx-1){  //rank!=root, X==nx-1, Y==0
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    L=bra[0]*W[0+nx*Y+N*i];
                    for (int x=1; x<nx-1; x++) {
                        L=L*bra[x]*W[x+nx*Y+N*i];
                    }
                    L=L*bra[nx-1];
                    L=reform2(L);
                    L=L*O[X+nx*Y+N*i];
                    temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
                    //PrintDat(L);
                    //PrintDat(L);
                }
                else {  //rank!=root, 0<X<nx-1, Y==0
                    //delta=Index("delta",NB*NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NB*NS);
                    L=bra[0]*W[0+nx*Y+N*i];
                    for (int x=1; x<X; x++) {
                        L=L*bra[x]*W[x+nx*Y+N*i];
                    }
                    R=bra[nx-1]*W[nx-1+nx*Y+N*i];
                    for (int x=nx-2; x>X; x--) {
                        R=R*bra[x]*W[x+nx*Y+N*i];
                    }
                    temp=L*bra[X]*R;
                    //Print(temp);
                    temp=reform3(temp);
                    temp=temp*O[X+nx*Y+N*i];
                    temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
                    //PrintDat(temp);
                }
            }
            
            else if (Y==ny-1){
                if (X==0) { //rank!=root, X==0, Y==ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    
                    R=ket[nx-1]*W[nx-1+nx*Y+N*i];
                    //Print(R);
                    
                    for (int x=nx-2; x>0; x--) {
                        R=R*ket[x]*W[x+nx*Y+N*i];
                    }
                    
                    R=R*ket[0];
                    
                    //Print(R);
                    R=reform2(R);
                    R=R*O[X+nx*Y+N*i];
                    temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
                    //Print(R);
                    //PrintDat(R);
                    
                    
                }
                
                else if (X==nx-1){  //rank!=root, X==nx-1, Y==ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    L=ket[0]*W[0+nx*Y+N*i];
                    for (int x=1; x<nx-1; x++) {
                        L=L*ket[x]*W[x+nx*Y+N*i];
                    }
                    L=L*ket[nx-1];
                    L=reform2(L);
                    L=L*O[X+nx*Y+N*i];
                    temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta2, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta2p);
                    //PrintDat(L);
                    //PrintDat(L);
                }
                else {  //rank!=root, 0<X<nx-1, Y==ny-1
                    //delta=Index("delta",NB*NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NB*NS);
                    L=ket[0]*W[0+nx*Y+N*i];
                    for (int x=1; x<X; x++) {
                        L=L*ket[x]*W[x+nx*Y+N*i];
                    }
                    R=ket[nx-1]*W[nx-1+nx*Y+N*i];
                    for (int x=nx-2; x>X; x--) {
                        R=R*ket[x]*W[x+nx*Y+N*i];
                    }
                    temp=L*ket[X]*R;
                    //Print(temp);
                    temp=reform3(temp);
                    temp=temp*O[X+nx*Y+N*i];
                    temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
                    //PrintDat(temp);
                }
            }
            else {
                if (X==0) { //rank!=root, X==0, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    R=bra[nx-1]*W[nx-1+nx*Y+N*i]*ket[nx-1];
                    //Print(R);
                    
                    for (int x=nx-2; x>0; x--) {
                        R=R*bra[x]*W[x+nx*Y+N*i]*ket[x];
                    }
                    R=R*bra[0]*ket[0];
                    //Print(R);
                    R=reform3(R);
                    R=R*O[X+nx*Y+N*i];
                    temp=group_indice4(R, R.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, R.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
                    //Print(R);
                    //PrintDat(R);
                    
                }
                else if (X==nx-1){  //rank!=root, X==nx-1, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NS);
                    L=bra[0]*W[0+nx*Y+N*i]*ket[0];
                    for (int x=1; x<nx-1; x++) {
                        L=L*bra[x]*W[x+nx*Y+N*i]*ket[x];
                    }
                    L=L*bra[nx-1]*ket[nx-1];
                    L=reform3(L);
                    L=L*O[X+nx*Y+N*i];
                    temp=group_indice4(L, L.indices()[0], O[X+nx*Y+N*i].indices()[0], delta3, L.indices()[1], O[X+nx*Y+N*i].indices()[1], delta3p);
                    //PrintDat(L);
                    //PrintDat(L);
                }
                else {  //rank!=root, 0<X<nx-1, 0<Y<ny-1
                    //delta=Index("delta",NB*NB*NB*NB*NS);
                    //deltap=Index("deltap",NB*NB*NB*NB*NS);
                    L=bra[0]*W[0+nx*Y+N*i]*ket[0];
                    for (int x=1; x<X; x++) {
                        L=L*bra[x]*W[x+nx*Y+N*i]*ket[x];
                    }
                    R=bra[nx-1]*W[nx-1+nx*Y+N*i]*ket[nx-1];
                    for (int x=nx-2; x>X; x--) {
                        R=R*bra[x]*W[x+nx*Y+N*i]*ket[x];
                    }
                    temp=L*bra[X]*ket[X]*R;
                    //Print(temp);
                    temp=reform4(temp);
                    temp=temp*O[X+nx*Y+N*i];
                    temp=group_indice4(temp, temp.indices()[0], O[X+nx*Y+N*i].indices()[0], delta4, temp.indices()[1], O[X+nx*Y+N*i].indices()[1], delta4p);
                    //PrintDat(temp);
                }
            }
            //------
            
        }
        
        return temp;
    }
    //------end: finding out the F tensors, to be used to optimize a(X,Y)-------
    
    
    //------begin:
    ITensor mapping_vector_to_tensor(Vector vector, ITensor tensor, int x, int y, int nx, int ny){
        if (y==0||y==ny-1) {
            if (x==0||x==nx-1) {    //x==0,y==0 or x==nx-1,y==0 or x==0,y==ny-1 or x==nx-1,y==ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int is=1; is<=NS; is++) {
                            tensor(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector((((i1-1)*NB+i2)-1)*NS+is);
                        }
                    }
                }
            }
            else {  //0<x<nx-1,y==0 or 0<x<nx-1,y==ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
        }
        else {
            if (x==0) { //x==0, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
            else if (x==nx-1){  //x==nx-1, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
            else {  //0<x<nx-1, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int i4=1; i4<=NB; i4++) {
                                for (int is=1; is<=NS; is++) {
                                    tensor(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return tensor;
    }
    //------end:
    
    //------begin:
    ITensor mapping_complex_vector_to_tensor(Vector vector_real, Vector vector_imag, ITensor tensor, int x, int y, int nx, int ny){
        
        ITensor tensor_real,tensor_imag;
        tensor_real=realPart(tensor);
        tensor_imag=imagPart(tensor);
        
        if (y==0||y==ny-1) {
            if (x==0||x==nx-1) {    //x==0,y==0 or x==nx-1,y==0 or x==0,y==ny-1 or x==nx-1,y==ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int is=1; is<=NS; is++) {
                            tensor_real(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector_real((((i1-1)*NB+i2)-1)*NS+is);
                            tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1))=vector_imag((((i1-1)*NB+i2)-1)*NS+is);
                        }
                    }
                }
            }
            else {  //0<x<nx-1,y==0 or 0<x<nx-1,y==ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor_real(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                                tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i2),tensor.indices()[2](i1),tensor.indices()[3](i3))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
        }
        else {
            if (x==0) { //x==0, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                                tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
            else if (x==nx-1){  //x==nx-1, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int is=1; is<=NS; is++) {
                                tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_real((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                                tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1))=vector_imag((((i1-1)*NB*NB+(i2-1)*NB+i3)-1)*NS+is);
                            }
                        }
                    }
                }
            }
            else {  //0<x<nx-1, 0<y<ny-1
                for (int i1=1; i1<=NB; i1++) {
                    for (int i2=1; i2<=NB; i2++) {
                        for (int i3=1; i3<=NB; i3++) {
                            for (int i4=1; i4<=NB; i4++) {
                                for (int is=1; is<=NS; is++) {
                                    tensor_real(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector_real((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
                                    tensor_imag(tensor.indices()[0](is),tensor.indices()[1](i3),tensor.indices()[2](i2),tensor.indices()[3](i1),tensor.indices()[4](i4))=vector_imag((((i1-1)*NB*NB*NB+(i2-1)*NB*NB+(i3-1)*NB+i4)-1)*NS+is);
                                }
                            }
                        }
                    }
                }
            }
        }
        tensor=tensor_real+Complex_i*tensor_imag;
        
        return tensor;
    }
    //------end:
    
    
    std::vector<ITensor> move_otho_next_bra(std::vector<ITensor> tensors, int otho_center, int y, int nx, int ny){  //y<ny-1
        int size=tensors.size();
        ITensor S,V,D;
        
        if (otho_center==0) {
            D=ITensor(tensors[otho_center].indices()[1]);
            svd(tensors[otho_center], S, V, D);    //svd
            
            R[otho_center+nx*(y+1)]=Index(nameint("L",otho_center+1+nx*(y+1)),S.indices()[1].m());
            tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*(y+1)]);
            L[otho_center+1+nx*(y+1)]=R[otho_center+nx*(y+1)];
            
            tensors[otho_center+1]=tensors[otho_center+1]*V*D;
            tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
            tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y+1)], tensors[otho_center+1].indices()[2]);
        }
        else if(otho_center<size-1){
            D=ITensor(tensors[otho_center].indices()[2]);
            svd(tensors[otho_center], S, V, D);    //svd
            
            R[otho_center+nx*(y+1)]=Index(nameint("L",otho_center+1+nx*(y+1)),S.indices()[2].m());
            tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*(y+1)]);
            L[otho_center+1+nx*(y+1)]=R[otho_center+nx*(y+1)];
            
            tensors[otho_center+1]=tensors[otho_center+1]*V*D;
            tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
            tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y+1)], tensors[otho_center+1].indices()[2]);
        }
        
        return tensors;
    }
    
    std::vector<ITensor> move_otho_next_ket(std::vector<ITensor> tensors, int otho_center, int y, int nx, int ny){  //y<0
        int size=tensors.size();
        ITensor S,V,D;
        
        if (otho_center==0) {
            D=ITensor(tensors[otho_center].indices()[1]);
            svd(tensors[otho_center], S, V, D);    //svd
            
            R[otho_center+nx*(y-1)]=Index(nameint("L",otho_center+1+nx*(y-1)),S.indices()[1].m());
            tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*(y-1)]);
            L[otho_center+1+nx*(y-1)]=R[otho_center+nx*(y-1)];
            
            tensors[otho_center+1]=tensors[otho_center+1]*V*D;
            tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
            tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y-1)], tensors[otho_center+1].indices()[2]);
        }
        else if(otho_center<size-1){
            D=ITensor(tensors[otho_center].indices()[2]);
            svd(tensors[otho_center], S, V, D);    //svd
            
            R[otho_center+nx*(y-1)]=Index(nameint("L",otho_center+1+nx*(y-1)),S.indices()[2].m());
            tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*(y-1)]);
            L[otho_center+1+nx*(y-1)]=R[otho_center+nx*(y-1)];
            
            tensors[otho_center+1]=tensors[otho_center+1]*V*D;
            tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
            tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*(y-1)], tensors[otho_center+1].indices()[2]);
        }
        
        return tensors;
    }
    
    std::vector<ITensor> move_otho_next_mpo(std::vector<ITensor> tensors, int otho_center, int y, int nx, int ny){
        int size=tensors.size();
        ITensor S,V,D;
        
        if (y==0||y==ny-1) {
            if (otho_center==0) {
                D=ITensor(tensors[otho_center].indices()[1]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[1].m());
                tensors[otho_center]=change_index_name2(S, tensors[otho_center].indices()[0], R[otho_center+nx*y]);
                L[otho_center+1+nx*y]=R[otho_center+nx*y];
                
                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[2]);
            }
            else if(otho_center<size-1){
                D=ITensor(tensors[otho_center].indices()[2]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[2].m());
                tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*y]);
                L[otho_center+1+nx*y]=R[otho_center+nx*y];
                
                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
                tensors[otho_center+1]=sort_indices_otho_3(tensors[otho_center+1], 0, 2, 1);
                tensors[otho_center+1]=change_index_name3(tensors[otho_center+1], tensors[otho_center+1].indices()[0], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[2]);
            }
        }
        else {
            if (otho_center==0) {
                D=ITensor(tensors[otho_center].indices()[2]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[2].m());
                tensors[otho_center]=change_index_name3(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], R[otho_center+nx*y]);
                L[otho_center+1+nx*y]=R[otho_center+nx*y];
                
                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
                tensors[otho_center+1]=sort_indices_otho_4(tensors[otho_center+1], 0, 1, 3, 2);
                tensors[otho_center+1]=change_index_name4(tensors[otho_center+1], tensors[otho_center+1].indices()[0], tensors[otho_center+1].indices()[1], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[3]);
            }
            else if(otho_center<size-1){
                D=ITensor(tensors[otho_center].indices()[3]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                R[otho_center+nx*y]=Index(nameint("L",otho_center+1+nx*y),S.indices()[3].m());
                tensors[otho_center]=change_index_name4(S, tensors[otho_center].indices()[0], tensors[otho_center].indices()[1], tensors[otho_center].indices()[2], R[otho_center+nx*y]);
                L[otho_center+1+nx*y]=R[otho_center+nx*y];
                
                tensors[otho_center+1]=tensors[otho_center+1]*V*D;
                tensors[otho_center+1]=sort_indices_otho_4(tensors[otho_center+1], 0, 1, 3, 2);
                tensors[otho_center+1]=change_index_name4(tensors[otho_center+1], tensors[otho_center+1].indices()[0], tensors[otho_center+1].indices()[1], L[otho_center+1+nx*y], tensors[otho_center+1].indices()[3]);
            }
        }
        
        return tensors;
    }
    
    
    std::vector<ITensor> move_otho_previous_mpo(std::vector<ITensor> tensors, int otho_center, int y, int nx, int ny){
        int size=tensors.size();
        ITensor S,V,D;
        
        if (y==0||y==ny-1) {
            if (otho_center==size-1) {
                D=ITensor(tensors[otho_center].indices()[0]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                D=sort_indices_otho_2(D,1,0);
                L[otho_center+nx*y]=Index(nameint("L",otho_center+nx*y),D.indices()[1].m());
                tensors[otho_center]=change_index_name2(D, D.indices()[0],L[otho_center+nx*y]);
                
                
                tensors[otho_center-1]=tensors[otho_center]*S*V;
                tensors[otho_center-1]=change_index_name3(tensors[otho_center-1], tensors[otho_center-1].indices()[0], tensors[otho_center-1].indices()[1], L[otho_center+nx*y]);
            }
            else if(otho_center>0){
                S=ITensor(tensors[otho_center].indices()[1]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                D=sort_indices_otho_3(D,1,0,2);
                L[otho_center+nx*y]=Index(nameint("L",otho_center+nx*y),D.indices()[1].m());
                tensors[otho_center]=change_index_name3(D, D.indices()[0], L[otho_center+nx*y], D.indices()[2]);
                
                tensors[otho_center-1]=tensors[otho_center-1]*S*V;
                tensors[otho_center-1]=change_index_name3(tensors[otho_center-1], tensors[otho_center-1].indices()[0], tensors[otho_center-1].indices()[1], L[otho_center+nx*y]);
            }
        }
        else {
            if (otho_center==size-1) {
                S=ITensor(tensors[otho_center].indices()[2]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                D=sort_indices_otho_3(D,1,2,0);
                L[otho_center+nx*y]=Index(nameint("L",otho_center+nx*y),D.indices()[2].m());
                tensors[otho_center]=change_index_name3(D, D.indices()[0], D.indices()[1], L[otho_center+nx*y]);
                
                tensors[otho_center-1]=tensors[otho_center-1]*S*V;
                tensors[otho_center-1]=change_index_name4(tensors[otho_center-1], tensors[otho_center-1].indices()[0], tensors[otho_center-1].indices()[1], tensors[otho_center-1].indices()[2],L[otho_center+nx*y]);
            }
            else if(otho_center>0){
                S=ITensor(tensors[otho_center].indices()[2]);
                svd(tensors[otho_center], S, V, D);    //svd
                
                D=sort_indices_otho_4(D,1,2,0,3);
                L[otho_center+nx*y]=Index(nameint("L",otho_center+nx*y),D.indices()[2].m());
                tensors[otho_center]=change_index_name4(D, D.indices()[0], D.indices()[1], L[otho_center+nx*y], D.indices()[3]);
                
                tensors[otho_center-1]=tensors[otho_center-1]*S*V;
                tensors[otho_center-1]=change_index_name4(tensors[otho_center-1], tensors[otho_center-1].indices()[0], tensors[otho_center-1].indices()[1], tensors[otho_center-1].indices()[2], L[otho_center+nx*y]);
            }
        }
        
        return tensors;
    }
    
    
//    Copyright (c) 2013, John D'Errico
//    All rights reserved.
//    
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are met:
    
//    * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in
//    the documentation and/or other materials provided with the distribution
//    
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
//    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//                           SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//                           INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//    POSSIBILITY OF SUCH DAMAGE.
    Matrix shift_matrix_k(Matrix A, int nn, int k){
        Vector d;
        Matrix U;
        EigenValues(A, d, U);
        
        //adding a tiny multiple (the lowest eigenvalue and an eps) of an identity matrix
        Matrix dd(nn,nn), eps(nn,nn);
        dd=0;
        eps=0;
        for(int j = 1; j <= nn; ++j){
            dd(j,j) = -d(0)*k*k;
            eps(j,j)=0.000001;
            //eps(j,j)=std::numeric_limits<double>::epsilon();
        }
        
        return A+dd+eps;
        //        U*(dd+eps)*Inverse(U);
    }
    
    Matrix make_positive_definite(Matrix A, int nn){
        
        Matrix U,V;
        Vector d;
        SVD(A,U,d,V);
        
        //Compute the symmetric polar factor
        Matrix dd(nn,nn);
        dd=0;
        for(int j = 1; j <= nn; ++j){
            dd(j,j) = d(j);
        }
        Matrix VV=V;
        for (int x=1; x<=nn; x++) {
            for (int y=1; y<=nn; y++) {
                VV(y,x)=V(x,y);
            }
        }
        Matrix H=VV*dd*V;
        
        //get H
        H=0.5*(A+H);
        
        //ensure the symmetry of H
        Matrix HT=H;
        for (int x=1; x<=nn; x++) {
            for (int y=1; y<=nn; y++) {
                HT(y,x)=H(x,y);
            }
        }
        H=0.5*(H+HT);
        
        int k=0;
        cholesky_flag=0;
        while (cholesky_flag==0) {
            //check whether H is Positive Definite
            Matrix HH=H;
            CholeskyFactorization(HH);
            
            k++;
            
            //if not, tweak by adding a tiny multiple of an identity matrix
            if (cholesky_flag==0) {
                H=shift_matrix_k(H, nn, k);
            }
        }
        
        return H;
    }
    
    
    
    Matrix shift_matrix(Matrix A, int nn){
        Vector d;
        Matrix U;
        EigenValues(A, d, U);
        
        //adding a tiny multiple (the lowest eigenvalue and an eps) of an identity matrix
        Matrix dd(nn,nn), eps(nn,nn);
        dd=0;
        eps=0;
        for(int j = 1; j <= nn; ++j){
            dd(j,j) = -d(0);
            eps(j,j)=0.000001;
        }
        
        return A+dd+eps;
        //        U*(dd+eps)*Inverse(U);
    }
    //    }
    /*
     std::vector<ITensor> fitting_algorithm(std::vector<ITensor> zipped_mps, std::vector<ITensor> mpo, std::vector<ITensor> original_mps, int dcut, int i){
     int zsize=zipped_mps.size();
     int msize=mpo.size();
     int osize=original_mps.size();
     int size;
     if (zsize==msize&&msize==osize) {
     size=zsize;
     }
     else {
     cout<<"ERROR: Sizes of zipped mps, mpo and the original mps does not match!!!"<<endl;
     exit(0);
     }
     
     int wind[nx];  //W index
     for (int x=0; x<nx; x++) {
     wind[x]=x+nx*(ny-1-i);
     }
     
     std::vector<Index> alpha(size),Alpha(size);
     
     for (int x=0; x<size-1; x++) {
     alpha[x]=Index(nameint("alpha",x),R[wind[x]].m()*R[wind[x]+size].m());
     }
     //PrintDat(original_mps[0]);
     //PrintDat(mpo[0]);
     
     //PrintDat(alpha[0]);
     //PrintDat(R[wind[0]+nx]);
     //PrintDat(U[wind[0]]);
     
     //mpo*original_mps
     std::vector<ITensor> N(size);
     
     for (int x=0; x<size; x++) {
     N[x]=mpo[x]*original_mps[x];
     if (x==0) {
     alpha[x]=Index(nameint("alpha",x),N[x].indices()[0].m()*N[x].indices()[2].m());
     N[x]=group_indice3(N[x], N[x].indices()[0], N[x].indices()[2], alpha[x], N[x].indices()[1]);
     Alpha[x]=Index(nameint("Alpha",x),dcut);
     N[x]=resize2(N[x], Alpha[x], dcut);
     }
     else if(x!=size-1) {
     alpha[x]=Index(nameint("alpha",x),N[x].indices()[1].m()*N[x].indices()[4].m());
     N[x]=group_indice5(N[x], N[x].indices()[0], N[x].indices()[3], alpha[x-1], N[x].indices()[1], N[x].indices()[4], alpha[x], N[x].indices()[2]);
     }
     else {
     N[x]=group_indice3(N[x], N[x].indices()[0], N[x].indices()[2], alpha[x-1], N[x].indices()[1]);
     }
     }
     
     //PrintDat(N[0]);
     
     //define trial tensors M and initialize M using zipped_mps
     std::vector<ITensor> M(size);
     int aaa;
     for (int x=0; x<size; x++) {
     if (x==0 || x==size-1) {
     M[x]=ITensor(N[x].indices()[0],N[x].indices()[1]);
     if (dcut<=M[x].indices()[0].m()) {
     aaa=dcut;
     }
     else {
     aaa=M[x].indices()[0].m();
     }
     for (int i0=1; i0<=aaa; i0++) {
     for (int i1=1; i1<=M[x].indices()[1].m(); i1++) {
     M[x](M[x].indices()[0](i0),M[x].indices()[1](i1))=zipped_mps[x](zipped_mps[x].indices()[0](i0),zipped_mps[x].indices()[1](i1));
     }
     }
     }
     else {
     M[x]=ITensor(N[x].indices()[0],N[x].indices()[1],N[x].indices()[2]);
     if (dcut<=M[x].indices()[0].m() && dcut<=M[x].indices()[1].m()) {
     aaa=dcut;
     }
     else if (M[x].indices()[0].m()<=dcut && M[x].indices()[0].m()<=M[x].indices()[1].m()) {
     aaa=M[x].indices()[0].m();
     }
     else {
     aaa=M[x].indices()[1].m();
     }
     for (int i0=1; i0<=aaa; i0++) {
     for (int i1=1; i1<=aaa; i1++) {
     for (int i2=1; i2<=M[x].indices()[2].m(); i2++) {
     M[x](M[x].indices()[0](i0),M[x].indices()[1](i1),M[x].indices()[2](i2))=zipped_mps[x](zipped_mps[x].indices()[0](i0),zipped_mps[x].indices()[1](i1),zipped_mps[x].indices()[2](i2));
     }
     }
     }
     }
     //PrintDat(M[x]);
     //Print(zipped_mps[x]);
     }
     
     //M=fitting_sweep(N, M, chi);
     
     
     
     return zipped_mps;
     }
     */
    
    
    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }
    
    ITensor define_site_tensor(ITensor a){
        if (LATTICE==1 && IGG==0 && NS==3 && NB==3) { //square lattice, trivial IGG, spin 1, bond dimension D=3
            srand (time(NULL));
            int np=6;
            double para[np];
            for (int i=0; i<np; i++) {
                para[i]=fRand(-1,1);
                //cout<<para[i]<<endl;
            }
            
            //1 -> -    2 -> 0    3 -> +;
            
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=para[0]-para[2]-3*para[3]+para[4]+para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=-para[0]-para[1]+2*para[3]+2*para[4];
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=-para[0]+2*para[4]+2*para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=para[0]+4*para[4];
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=para[0]+para[2]+3*para[3]+para[4]+para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=-para[0]+para[1]-2*para[3]+2*para[4];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=para[1]-para[2]+para[3]-3*para[4]+para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=-para[1]+para[2]-para[3]-3*para[4]+para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=para[1]+para[2]+para[3]-3*para[4]-para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 1 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=-para[1]-para[2]-para[3]-3*para[4]-para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=para[2]-3*para[3]-3*para[4]-para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=-para[2]+3*para[3]-3*para[4]-para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=6*para[3]-2*para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=-6*para[3]-2*para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 1 ))=6*para[4]-2*para[5];
            a(a.indices()[0]( 3 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 1 ))=4*para[5];
            
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=para[0]+para[2]-3*para[3]+para[4]-para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=-para[0]+para[2]-3*para[3]-para[4]+para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=para[0]-2*para[4]+2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=-para[0]+2*para[4]-2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=para[0]-para[2]+3*para[3]+para[4]-para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=-para[0]-para[2]+3*para[3]-para[4]+para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=para[1]+4*para[3];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=para[1]-2*para[3]+2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=para[1]-2*para[3]-2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=-para[1]+2*para[3]+2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=-para[1]+2*para[3]-2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=-para[1]-4*para[3];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=-para[2]-3*para[3]+3*para[4]-para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=para[2]+3*para[3]+3*para[4]-para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=para[2]+3*para[3]-3*para[4]+para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=-para[2]-3*para[3]-3*para[4]+para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 3 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 1 ))=-6*para[4]-2*para[5];
            a(a.indices()[0]( 2 ),a.indices()[1]( 1 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 3 ))=6*para[4]+2*para[5];
            
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=para[0]-para[1]+2*para[3]-2*para[4];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 3 ),a.indices()[4]( 1 ))=-para[0]-para[2]-3*para[3]-para[4]-para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=para[0]-2*para[4]-2*para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 2 ))=-para[0]-4*para[4];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 2 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=para[0]+para[1]-2*para[3]-2*para[4];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 3 ))=-para[0]+para[2]+3*para[3]-para[4]-para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 1 ))=para[1]+para[2]+para[3]+3*para[4]+para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=para[1]-para[2]+para[3]+3*para[4]-para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 3 ),a.indices()[3]( 1 ),a.indices()[4]( 2 ))=-para[1]-para[2]-para[3]+3*para[4]+para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 3 ),a.indices()[2]( 1 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=-para[1]+para[2]-para[3]+3*para[4]-para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=para[2]-3*para[3]+3*para[4]+para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 2 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=-para[2]+3*para[3]+3*para[4]+para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 3 ),a.indices()[4]( 2 ))=6*para[3]+2*para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 1 ),a.indices()[2]( 3 ),a.indices()[3]( 2 ),a.indices()[4]( 3 ))=-6*para[3]+2*para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 2 ),a.indices()[2]( 1 ),a.indices()[3]( 3 ),a.indices()[4]( 3 ))=-6*para[4]+2*para[5];
            a(a.indices()[0]( 1 ),a.indices()[1]( 1 ),a.indices()[2]( 2 ),a.indices()[3]( 3 ),a.indices()[4]( 3 ))=-4*para[5];
            //a(a.indices()[0](  ),a.indices()[1](  ),a.indices()[2](  ),a.indices()[3](  ),a.indices()[4](  ))=
        }
        
        a/=a.norm();
        //PrintDat(a);
        return a;
    }
    
    ITensor define_bond_h_tensor(ITensor b, int x, int y){
        if (LATTICE==1 && IGG==0 && NS==3 && NB==3 && CHIC4==1) { //square lattice, trivial IGG, spin 1, bond dimension D=3
            if (x==0 || x==NX) {
                b(b.indices()[0]( 2 ))=-1;
            }
            else {
                b(b.indices()[0]( 1 ),b.indices()[1]( 3 ))=1;
                b(b.indices()[0]( 2 ),b.indices()[1]( 2 ))=-1;
                b(b.indices()[0]( 3 ),b.indices()[1]( 1 ))=1;
            }
        }
        
        b/=b.norm();
        return b;
    }
    
    ITensor define_bond_v_tensor(ITensor b, int x, int y){
        if (LATTICE==1 && IGG==0 && NS==3 && NB==3 && CHIC4==1) { //square lattice, trivial IGG, spin 1, bond dimension D=3
            if (y==0 || y==NY) {
                b(b.indices()[0]( 2 ))=-1;
            }
            else {
                b(b.indices()[0]( 1 ),b.indices()[1]( 3 ))=1;
                b(b.indices()[0]( 2 ),b.indices()[1]( 2 ))=-1;
                b(b.indices()[0]( 3 ),b.indices()[1]( 1 ))=1;
            }
        }
        
        b/=b.norm();
        return b;
    }
}

#endif





















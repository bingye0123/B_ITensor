//tensor_utility.h
//include change_index_name2/3/4, group_indice2/4/6/8, bra_otho_center_1, mpo_otho_center_1

#ifndef ___TENSOR_UTILITY
#define ___TENSOR_UTILITY


#include "core.h"
#include "parameters.h"
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

using namespace std;
namespace itensor {
    
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
        
        //a/=a.norm();
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
        
        ///////////b/=(b.norm()*1000);
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
        
        //////////b/=(b.norm()*1000);
        return b;
    }
    
    
}

#endif





















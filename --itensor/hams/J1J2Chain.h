//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_J1J2CHAIN_H
#define __ITENSOR_HAMS_J1J2CHAIN_H
#include "../mpo.h"

class J1J2Chain
    {
    public:

    J1J2Chain(const Model& model, 
              const OptSet& opts = Global::opts());

    operator const IQMPO&() { init(); return H; }

    operator MPO() { init(); return H.toMPO(); }

    private:

    /////////////
    //
    // Data Members

    const Model& model_;

    Real J1_, 
         J2_;

    IQMPO H;

    bool initted_;

    //
    /////////////

    void
    init();

    };

inline J1J2Chain::
J1J2Chain(const Model& model,
          const OptSet& opts)
    :
    model_(model),
    initted_(false)
    {
    J1_ = opts.getReal("J1",1.);
    J2_ = opts.getReal("J2",0.);
    }
    

void inline J1J2Chain::
init()
    {
    if(initted_) return;

    H = IQMPO(model_);
    const int N = model_.N();

    const int kk = 2,
              ds = 5,    //start of diagonal part
              kd = kk+2, //bond dimension of diagonal part
              k  = (ds-1)+kd;    //total bond dimension

    std::vector<IQIndex> iqlinks(N+1);

    //The names of these indices refer to their Sz quantum numbers (plus or minus), 
    std::vector<Index> q0(N+1),
                       qP(N+1), 
                       qM(N+1);

    for(int i = 0; i <= N; ++i)
        {
        qP.at(i) = Index(nameint("qP_",i),kk);
        qM.at(i) = Index(nameint("qM_",i),kk);
        q0.at(i) = Index(nameint("q0_",i),kd);

        iqlinks.at(i) = IQIndex(nameint("hl",i),
                                qP[i],QN(-2),
                                qM[i],QN(+2),
                                q0[i],QN( 0));
        }

    for(int j = 1; j <= N; ++j)
        {
        //Create j^th A (an IQTensor)
        IQTensor &W = H.Anc(j);
        IQIndex row = conj(iqlinks.at(j-1)),
                col = iqlinks.at(j);

        W = IQTensor(conj(model_.si(j)),model_.siP(j),row,col);

        //Identity string operators
        W += model_.op("Id",j) * row(ds) * col(ds);
        W += model_.op("Id",j) * row(k) * col(k);

        //S+ S- terms
        W += model_.op("Sp",j) * row(k) * col(1) * (J1_/2.);
        W += model_.op("Sp",j) * row(k) * col(2) * (J2_/2.);
        W += model_.op("Id",j) * row(2) * col(1);
        W += model_.op("Sm",j) * row(1) * col(ds);

        //S- S+ terms
        W += model_.op("Sm",j) * row(k) * col(3) * (J1_/2.);
        W += model_.op("Sm",j) * row(k) * col(4) * (J2_/2.);
        W += model_.op("Id",j) * row(4) * col(3);
        W += model_.op("Sp",j) * row(3) * col(ds);

        //Sz Sz terms
        W += model_.op("Sz",j) * row(k) * col(6) * J1_;
        W += model_.op("Sz",j) * row(k) * col(7) * J2_;
        W += model_.op("Id",j) * row(7) * col(6);
        W += model_.op("Sz",j) * row(6) * col(ds);
        }

    H.Anc(1) *= IQTensor(iqlinks.at(0)(k));
    H.Anc(N) *= IQTensor(conj(iqlinks.at(N)(ds)));

    initted_ = true;
    }


#endif

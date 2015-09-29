//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMBUILDER_H
#define __ITENSOR_HAMBUILDER_H
#include "mpo.h"

//
//
// HamBuilder
//
// Class for creating product-operator MPOs,
// usually to be combined into a more complex
// MPO such as a Hamiltonian.
//

class HamBuilder
    {
    public:

    HamBuilder(const Model& mod);

    int
    N() const { return N_; }

    //Sets res to an identity MPO
    template <class Tensor>
    void
    getMPO(MPOt<Tensor>& res, Real fac=1) const;

    //Sets the j1'th operator of res to op1
    //and the rest to the identity
    template <class Tensor>
    void
    getMPO(int j1, const Tensor& op1, 
           MPOt<Tensor>& res, Real fac=1) const;

    //Sets the j1'th operator of res to op1,
    //the j2'th operator of res to op2,
    //and the rest to the identity
    template <class Tensor>
    void
    getMPO(int j1, const Tensor& op1,
           int j2, const Tensor& op2,
           MPOt<Tensor>& res, Real fac=1) const;

    // etc.

    template <class Tensor>
    void
    getMPO(int j1, const Tensor& op1,
           int j2, const Tensor& op2,
           int j3, const Tensor& op3,
           MPOt<Tensor>& res, Real fac=1) const;

    template <class Tensor>
    void
    getMPO(int j1, const Tensor& op1,
           int j2, const Tensor& op2,
           int j3, const Tensor& op3,
           int j4, const Tensor& op4,
           MPOt<Tensor>& res, Real fac=1) const;

    private:

    /////////////////
    //
    // Data Members

    const Model& mod_;
    const int N_;

    //
    /////////////////

    void
    putLinks(MPOt<ITensor>& res) const;
    void
    putLinks(MPOt<IQTensor>& res) const;

    template <class Tensor>
    void
    initialize(MPOt<Tensor>& res) const;

    static int
    hamNumber()
        {
        static int num_ = 0;
        ++num_;
        return num_;
        }


    };

inline HamBuilder::
HamBuilder(const Model& mod)
    :
    mod_(mod),
    N_(mod.N())
    { }

template <class Tensor>
void HamBuilder::
getMPO(MPOt<Tensor>& res, Real fac) const
    {
    initialize(res);
    putLinks(res);
    res.Anc(1) *= fac;
    }

template <class Tensor>
void HamBuilder::
getMPO(int j1, const Tensor& op1, 
       MPOt<Tensor>& res, Real fac) const
    {
    initialize(res);
#ifdef DEBUG
    if(!hasindex(op1,mod_.si(j1)))
        {
        Print(j1);
        Print(op1.indices());
        Error("Tensor does not have correct Site index");
        }
#endif
    res.Anc(j1) = op1;
    res.Anc(j1) *= fac;
    putLinks(res);
    }

template <class Tensor>
void HamBuilder::
getMPO(int j1, const Tensor& op1,
       int j2, const Tensor& op2,
       MPOt<Tensor>& res, Real fac) const
    {
    initialize(res);
#ifdef DEBUG
    if(!hasindex(op1,mod_.si(j1)))
        {
        Print(j1);
        Print(op1.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op2,mod_.si(j2)))
        {
        Print(j2);
        Print(op2.indices());
        Error("Tensor does not have correct Site index");
        }
#endif
    res.Anc(j1) = op1;
    res.Anc(j2) = op2;
    res.Anc(j1) *= fac;
    putLinks(res);
    }

template <class Tensor>
void HamBuilder::
getMPO(int j1, const Tensor& op1,
       int j2, const Tensor& op2,
       int j3, const Tensor& op3,
       MPOt<Tensor>& res, Real fac) const
    {
    initialize(res);
#ifdef DEBUG
    if(!hasindex(op1,mod_.si(j1)))
        {
        Print(j1);
        Print(op1.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op2,mod_.si(j2)))
        {
        Print(j2);
        Print(op2.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op3,mod_.si(j3)))
        {
        Print(j3);
        Print(op3.indices());
        Error("Tensor does not have correct Site index");
        }
#endif
    res.Anc(j1) = op1;
    res.Anc(j2) = op2;
    res.Anc(j3) = op3;
    res.Anc(j1) *= fac;
    putLinks(res);
    }

template <class Tensor>
void HamBuilder::
getMPO(int j1, const Tensor& op1,
       int j2, const Tensor& op2,
       int j3, const Tensor& op3,
       int j4, const Tensor& op4,
       MPOt<Tensor>& res, Real fac) const
    {
    initialize(res);
#ifdef DEBUG
    if(!hasindex(op1,mod_.si(j1)))
        {
        Print(j1);
        Print(op1.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op2,mod_.si(j2)))
        {
        Print(j2);
        Print(op2.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op3,mod_.si(j3)))
        {
        Print(j3);
        Print(op3.indices());
        Error("Tensor does not have correct Site index");
        }
    if(!hasindex(op4,mod_.si(j4)))
        {
        Print(j4);
        Print(op4.indices());
        Error("Tensor does not have correct Site index");
        }
#endif
    res.Anc(j1) = op1;
    res.Anc(j2) = op2;
    res.Anc(j3) = op3;
    res.Anc(j4) = op4;
    res.Anc(j1) *= fac;
    putLinks(res);
    }


template <class Tensor>
void HamBuilder::
initialize(MPOt<Tensor>& res) const
    {
    res = MPOt<Tensor>(mod_);
    for(int j = 1; j <= N_; ++j)
        res.Anc(j) = mod_.op("Id",j);
    }


void inline HamBuilder::
putLinks(MPOt<ITensor>& res) const
    {
    int ver = hamNumber();
    std::vector<Index> links(N_);
    for(int i = 1; i < N_; ++i)
        {
        boost::format nm = boost::format("h%d-%d") % ver % i;
        links.at(i) = Index(nm.str());
        }
    res.Anc(1) *= links.at(1)(1);
    for(int i = 1; i < N_; ++i)
        {
        res.Anc(i) *= links.at(i-1)(1);
        res.Anc(i) *= links.at(i)(1);
        }
    res.Anc(N_) *= links.at(N_-1)(1);
    }

void inline HamBuilder::
putLinks(MPOt<IQTensor>& res) const
    {
    QN q;

    int ver = hamNumber();
    std::vector<IQIndex> links(N_);
    for(int i = 1; i < N_; ++i)
        {
        boost::format nm = boost::format("h%d-%d") % ver % i,
                      Nm = boost::format("H%d-%d") % ver % i;
        q += div(res.A(i));
        links.at(i) = IQIndex(Nm.str(),
                             Index(nm.str()),q);
        }

    res.Anc(1) *= links.at(1)(1);
    for(int i = 2; i < N_; ++i)
        {
        res.Anc(i) *= conj(links.at(i-1)(1));
        res.Anc(i) *= links.at(i)(1);
        }
    res.Anc(N_) *= conj(links.at(N_-1)(1));
    }

#endif

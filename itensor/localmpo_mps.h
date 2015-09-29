//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO_MPS
#define __ITENSOR_LOCALMPO_MPS
#include "localmpo.h"

template <class Tensor>
class LocalMPO_MPS
    {
    public:

    //
    // Constructors
    //

    LocalMPO_MPS();

    LocalMPO_MPS(const MPOt<Tensor>& Op, 
                 const std::vector<MPSt<Tensor> >& psis,
                 const OptSet& opts = Global::opts());

    //
    // Typedefs
    //

    typedef typename Tensor::IndexT
    IndexT;

    typedef typename Tensor::CombinerT
    CombinerT;

    typedef LocalMPO<Tensor>
    LocalMPOType;


    //
    // Sparse matrix methods
    //

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const { return lmpo_.expect(phi); }

    Tensor
    deltaRho(const Tensor& AA, 
             const CombinerT& comb, Direction dir) const
        { return lmpo_.deltaRho(AA,comb,dir); }

    Tensor
    diag() const { return lmpo_.diag(); }

    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    int
    size() const { return lmpo_.size(); }

    bool
    isNull() const { return Op_ == 0; }

    Real
    weight() const { return weight_; }
    void
    weight(Real val) { weight_ = val; }

    bool
    doWrite() const { return lmpo_.doWrite(); }
    void
    doWrite(bool val) { lmpo_.doWrite(val); }

    private:

    /////////////////
    //
    // Data Members
    //

    const MPOt<Tensor>* Op_;
    const std::vector<MPSt<Tensor> >* psis_;
    //LocalMPO object representing projected version
    //of the MPO Op_
    LocalMPO<Tensor> lmpo_; 
    //LocalMPO objects representing projected version
    //of each MPS in psis_
    std::vector<LocalMPO<Tensor> > lmps_;

    Real weight_;

    //
    /////////////////

    };

template <class Tensor>
inline LocalMPO_MPS<Tensor>::
LocalMPO_MPS()
    : 
    Op_(0),
    psis_(0),
    weight_(1)
    { }

template <class Tensor>
inline LocalMPO_MPS<Tensor>::
LocalMPO_MPS(const MPOt<Tensor>& Op,
             const std::vector<MPSt<Tensor> >& psis,
             const OptSet& opts)
    : 
    Op_(&Op),
    psis_(&psis),
    lmps_(psis.size()),
    weight_(1)
    { 
    lmpo_ = LocalMPOType(Op);

    for(size_t j = 0; j < lmps_.size(); ++j)
        lmps_[j] = LocalMPOType(psis[j]);

    if(opts.defined("Weight"))
        weight(opts.getReal("Weight"));
    }

template <class Tensor>
void inline LocalMPO_MPS<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    lmpo_.product(phi,phip);

    Tensor outer;
    for(size_t j = 0; j < lmps_.size(); ++j)
        {
        lmps_[j].product(phi,outer);
        outer *= weight_;
        phip += outer;
        }
    }

template <class Tensor>
template <class MPSType> 
void inline LocalMPO_MPS<Tensor>::
position(int b, const MPSType& psi)
    {
    lmpo_.position(b,psi);
    for(size_t j = 0; j < lmps_.size(); ++j)
        lmps_[j].position(b,psi);
    }

#endif

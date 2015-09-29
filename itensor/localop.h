//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCAL_OP
#define __ITENSOR_LOCAL_OP
#include "iqtensor.h"

#define Cout std::cout
#define Endl std::endl

//
// The LocalOp class represents
// an MPO or other operator that
// has been projected into the
// reduced Hilbert space of 
// two sites of an MPS.
//
//   .-              -.
//   |    |      |    |
//   L - Op1 -- Op2 - R
//   |    |      |    |
//   '-              -'
//
// (Note that L, Op1, Op2 and R
//  are not required to have this
//  precise structure. L and R
//  can even be null in which case
//  they will not be used.)
//


template <class Tensor>
class LocalOp
    {
    public:

    typedef typename Tensor::IndexT
    IndexT;

    typedef typename Tensor::CombinerT
    CombinerT;

    //
    // Constructors
    //

    LocalOp(const OptSet& opts = Global::opts());

    LocalOp(const Tensor& Op1, const Tensor& Op2,
            const OptSet& opts = Global::opts());

    LocalOp(const Tensor& Op1, const Tensor& Op2, 
            const Tensor& L, const Tensor& R,
            const OptSet& opts = Global::opts());

    //
    // Sparse Matrix Methods
    //

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const;

    Tensor
    deltaRho(const Tensor& rho, 
             const CombinerT& comb, Direction dir) const;

    Tensor
    deltaPhi(const Tensor& phi) const;

    Tensor
    diag() const;

    int
    size() const;

    //
    // Accessor Methods
    //

    void
    update(const Tensor& Op1, const Tensor& Op2);

    void
    update(const Tensor& Op1, const Tensor& Op2, 
           const Tensor& L, const Tensor& R);

    const Tensor&
    Op1() const 
        { 
        if(isNull()) Error("LocalOp is null");
        return *Op1_;
        }

    const Tensor&
    Op2() const 
        { 
        if(isNull()) Error("LocalOp is null");
        return *Op2_;
        }

    const Tensor&
    L() const 
        { 
        if(isNull()) Error("LocalOp is null");
        return *L_;
        }

    const Tensor&
    R() const 
        { 
        if(isNull()) Error("LocalOp is null");
        return *R_;
        }

    const Tensor&
    bondTensor() const 
        { 
        makeBond();
        return bond_;
        }

    bool
    combineMPO() const { return combine_mpo_; }
    void
    combineMPO(bool val) { combine_mpo_ = val; }

    bool
    isNull() const { return Op1_ == 0; }

    bool
    LIsNull() const;

    bool
    RIsNull() const;

    static const LocalOp& Null()
        {
        static LocalOp Null_;
        return Null_;
        }

    void
    operator=(const LocalOp& other)
        {
        Op1_ = other.Op1_;
        Op2_ = other.Op2_;
        L_ = other.L_;
        R_ = other.R_;
        combine_mpo_ = other.combine_mpo_;
        bond_ = other.bond_;
        }

    private:

    /////////////////
    //
    // Data Members
    //

    const Tensor *Op1_, *Op2_; 
    const Tensor *L_, *R_; 
    bool combine_mpo_;
    mutable int size_;
    mutable Tensor bond_;

    //
    /////////////////

    void
    makeBond() const;

    void
    processOpts(const OptSet& opts)
        {
        combine_mpo_ = opts.getBool("CombineMPO",true);
        }

    };

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const OptSet& opts)
    :
    Op1_(0),
    Op2_(0),
    L_(0),
    R_(0),
    combine_mpo_(true),
    size_(-1)
    { 
    processOpts(opts);
    }

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const Tensor& Op1, const Tensor& Op2,
        const OptSet& opts)
    : 
    Op1_(0),
    Op2_(0),
    L_(0),
    R_(0),
    combine_mpo_(true),
    size_(-1)
    {
    processOpts(opts);
    update(Op1,Op2);
    }

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const Tensor& Op1, const Tensor& Op2, 
        const Tensor& L, const Tensor& R,
        const OptSet& opts)
    : 
    Op1_(0),
    Op2_(0),
    L_(0),
    R_(0),
    combine_mpo_(true),
    size_(-1)
    {
    processOpts(opts);
    update(Op1,Op2,L,R);
    }

template <class Tensor>
void inline LocalOp<Tensor>::
update(const Tensor& Op1, const Tensor& Op2)
    {
    Op1_ = &Op1;
    Op2_ = &Op2;
    L_ = 0;
    R_ = 0;
    size_ = -1;
    bond_ = Tensor();
    }

template <class Tensor>
void inline LocalOp<Tensor>::
update(const Tensor& Op1, const Tensor& Op2, 
       const Tensor& L, const Tensor& R)
    {
    update(Op1,Op2);
    L_ = &L;
    R_ = &R;
    }

template <class Tensor>
bool inline LocalOp<Tensor>::
LIsNull() const
    {
    if(L_ == 0) return true;
    return L_->isNull();
    }

template <class Tensor>
bool inline LocalOp<Tensor>::
RIsNull() const
    {
    if(R_ == 0) return true;
    return R_->isNull();
    }

template <class Tensor>
void inline LocalOp<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    if(this->isNull()) Error("LocalOp is null");

    const Tensor& Op1 = *Op1_;
    const Tensor& Op2 = *Op2_;

    if(LIsNull())
        {
        phip = phi;

        if(!RIsNull()) 
            phip *= R(); //m^3 k d

        if(combine_mpo_)
            {
            phip *= bondTensor();
            }
        else
            {
            phip *= Op2; //m^2 k^2
            phip *= Op1; //m^2 k^2
            }
        }
    else
        {
        phip = phi * L(); //m^3 k d

        if(combine_mpo_)
            {
            phip *= bondTensor();
            }
        else
            {
            phip *= Op1; //m^2 k^2
            phip *= Op2; //m^2 k^2
            }

        if(!RIsNull()) 
            phip *= R();
        }

    phip.mapprime(1,0);
    }

template <class Tensor>
Real inline LocalOp<Tensor>::
expect(const Tensor& phi) const
    {
    Tensor phip;
    product(phi,phip);
    return Dot(phip,phi);
    }

template <class Tensor>
Tensor inline LocalOp<Tensor>::
deltaRho(const Tensor& AA, const CombinerT& comb, Direction dir) const
    {
    Tensor delta(AA);
    if(dir == Fromleft)
        {
        if(!LIsNull()) delta *= L();
        delta *= (*Op1_);
        }
    else //dir == Fromright
        {
        if(!RIsNull()) delta *= R();
        delta *= (*Op2_);
        }

    delta.noprime();
    delta = comb * delta;
    
    delta *= conj(primed(delta,comb.right()));

    return delta;
    }

template <class Tensor>
Tensor inline LocalOp<Tensor>::
deltaPhi(const Tensor& phi) const
    {
    Tensor deltaL(phi),
           deltaR(phi);

    if(!LIsNull()) 
        {
        deltaL *= L();
        }

    if(!RIsNull()) 
        {
        deltaR *= R();
        }

    const Tensor& Op1 = *Op1_;
    const Tensor& Op2 = *Op2_;

    deltaL *= Op1;
    deltaR *= Op2;

    IndexT hl = commonIndex(Op1,Op2,Link);

    deltaL.trace(hl);
    deltaL.mapprime(1,0);

    deltaR.trace(hl);
    deltaR.mapprime(1,0);

    deltaL += deltaR;

    return deltaL;
    }

template <>
IQTensor inline LocalOp<IQTensor>::
deltaPhi(const IQTensor& phi) const
    {
    IQTensor deltaL(phi),
           deltaR(phi);

    if(!LIsNull()) 
        {
        deltaL *= L();
        }

    if(!RIsNull()) 
        {
        deltaR *= R();
        }

    const IQTensor& Op1 = *Op1_;
    const IQTensor& Op2 = *Op2_;

    deltaL *= Op1;
    deltaR *= Op2;

    IndexT hl = commonIndex(Op1,Op2,Link);

    deltaL.trace(hl);
    deltaL.mapprime(1,0);

    deltaR.trace(hl);
    deltaR.mapprime(1,0);

    deltaL += deltaR;

    IQTensor delta(deltaL);
    delta *= 0;

    const QN targetQn = div(phi);

    Foreach(const ITensor& block, deltaL.blocks())
        {
        QN div;
        Foreach(const Index& I, block.indices())
            div += qn(deltaL,I);

        if(div == targetQn)
            delta += block;
        }

    return delta;
    }

template <class Tensor>
Tensor inline LocalOp<Tensor>::
diag() const
    {
    if(this->isNull()) Error("LocalOp is null");

    const Tensor& Op1 = *Op1_;
    const Tensor& Op2 = *Op2_;

    IndexT toTie;
    bool found = false;

    Foreach(const IndexT& s, Op1.indices())
        {
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) 
        {
        Print(Op1);
        Error("Couldn't find Index");
        }

    Tensor Diag = tieIndices(Op1,toTie,primed(toTie),toTie);

    found = false;
    Foreach(const IndexT& s, Op2.indices())
        {
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag *= tieIndices(Op2,toTie,primed(toTie),toTie);

    if(!LIsNull())
        {
        found = false;
        Foreach(const IndexT& ll, L().indices())
            {
            if(ll.primeLevel() == 0 && hasindex(L(),primed(ll)))
                {
                toTie = ll;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= tieIndices(L(),toTie,primed(toTie),toTie);
        else
            Diag *= L();
        }

    if(!RIsNull())
        {
        found = false;
        Foreach(const IndexT& rr, R().indices())
            {
            if(rr.primeLevel() == 0 && hasindex(R(),primed(rr)))
                {
                toTie = rr;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= tieIndices(R(),toTie,primed(toTie),toTie);
        else
            Diag *= R();
        }

    Diag.conj();
    return Diag;
    }

template <class Tensor>
int inline LocalOp<Tensor>::
size() const
    {
    if(this->isNull()) Error("LocalOp is null");
    if(size_ == -1)
        {
        //Calculate linear size of this 
        //op as a square matrix
        size_ = 1;
        if(!LIsNull()) 
            {
            Foreach(const IndexT& I, L().indices())
                {
                if(I.primeLevel() > 0)
                    {
                    size_ *= I.m();
                    break;
                    }
                }
            }
        if(!RIsNull()) 
            {
            Foreach(const IndexT& I, R().indices())
                {
                if(I.primeLevel() > 0)
                    {
                    size_ *= I.m();
                    break;
                    }
                }
            }

        size_ *= findtype(*Op1_,Site).m();
        size_ *= findtype(*Op2_,Site).m();
        }
    return size_;
    }

template <class Tensor>
void inline LocalOp<Tensor>::
makeBond() const
    {
    if(bond_.isNull()) 
        {
        if(!combine_mpo_) Error("combineMPO is false");
        bond_ = (*Op1_) * (*Op2_);
        }
    }

#undef Cout
#undef Endl

#endif

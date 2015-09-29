//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINLESS_H
#define __ITENSOR_SPINLESS_H
#include "../model.h"

class Spinless : public Model
    {
    public:

    Spinless();

    Spinless(int N, const OptSet& opts = Global::opts());

    IQIndexVal
    Emp(int i) const;

    IQIndexVal
    Occ(int i) const;

    IQIndexVal
    EmpP(int i) const;

    IQIndexVal
    OccP(int i) const;

    IQTensor
    projEmp(int i) const { return getOp(i,"ProjEmp"); }

    IQTensor
    projOcc(int i) const { return getOp(i,"ProjOcc"); }

    private:

    int
    getN() const;

    const IQIndex&
    getSi(int i) const;
    
    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts = Global::opts()) const;

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

    void
    constructSites();

    IQTensor
    makeN(int i) const;
    IQTensor
    makeC(int i) const;
    IQTensor
    makeCdag(int i) const;
    IQTensor
    makeFermiPhase(int i) const;
    IQTensor
    makeProjEmp(int i) const;
    IQTensor
    makeProjOcc(int i) const;

        
    //Data members -----------------

    int N_;

    bool odd_even_up_down_;

    bool conserve_Nf_;

    std::vector<IQIndex> site_;

    };

inline Spinless::
Spinless()
    : N_(-1),
      odd_even_up_down_(false),
      conserve_Nf_(true)
    { 
    }

inline Spinless::
Spinless(int N, const OptSet& opts)
    : N_(N),
      site_(N_+1)
    { 
    odd_even_up_down_ = opts.getBool("OddEvenUpDown",false);
    conserve_Nf_ = opts.getBool("ConserveNf",true);
    constructSites();
    }

void inline Spinless::
constructSites()
    {
    const int occ = (conserve_Nf_ ? 1 : 0);
    if(odd_even_up_down_)
        {
        for(int i = 1; i <= N_; ++i)
            {
            if(i%2==1)
                {
                site_.at(i) = IQIndex(nameint("Spinless Up site=",i),
                Index(nameint("Emp for Up site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Up site",i),1,Site),QN(+1,occ,1));
                }
            else
                {
                site_.at(i) = IQIndex(nameint("Spinless Dn site=",i),
                Index(nameint("Emp for Dn site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Dn site",i),1,Site),QN(-1,occ,1));
                }
            }
        }
    else
        {
        for(int i = 1; i <= N_; ++i)
            {
            site_.at(i) = IQIndex(nameint("Spinless site=",i),
            Index(nameint("Emp for site",i),1,Site),QN(0,0,0),
            Index(nameint("Occ for site",i),1,Site),QN(0,occ,1));
            }
        }
    }

void inline Spinless::
doRead(std::istream& s)
    {
    s.read((char*) &odd_even_up_down_,sizeof(odd_even_up_down_));
    s.read((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

void inline Spinless::
doWrite(std::ostream& s) const
    {
    s.write((char*) &odd_even_up_down_,sizeof(odd_even_up_down_));
    s.write((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline Spinless::
getN() const
    { return N_; }

inline 
const IQIndex& Spinless::
getSi(int i) const
    { return site_.at(i); }


IQIndexVal inline Spinless::
Emp(int i) const
    {
    return getSi(i)(1);
    }

IQIndexVal inline Spinless::
Occ(int i) const
    {
    return getSi(i)(2);
    }

IQIndexVal inline Spinless::
EmpP(int i) const
    {
    return primed(getSi(i))(1);
    }

IQIndexVal inline Spinless::
OccP(int i) const
    {
    return primed(getSi(i))(2);
    }

inline IQIndexVal Spinless::
getState(int i, const String& state) const
    {
    if(state == "Emp" || state == "0") 
        {
        return getSi(i)(1);
        }
    else 
    if(state == "Occ" || state == "1") 
        {
        return getSi(i)(2);
        }
    else
        {
        Error("State " + state + " not recognized");
        return getSi(i)(1);
        }
    }

inline IQTensor Spinless::
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Emp(s(1)),
               EmpP(sP(1)),
               Occ(s(2)),
               OccP(sP(2));

    IQTensor Op(conj(s),sP);

    if(opname == "N" || opname == "n")
        {
        Op(Occ,OccP) = 1;
        }
    else
    if(opname == "C")
        {
        Op(Occ,EmpP) = 1;
        }
    else
    if(opname == "Cdag")
        {
        Op(Emp,OccP) = 1;
        }
    else
    if(opname == "A")
        {
        Op(Occ,EmpP) = 1;
        }
    else
    if(opname == "Adag")
        {
        Op(Emp,OccP) = 1;
        }
    else
    if(opname == "F" || opname == "FermiPhase")
        {
        Op(Emp,EmpP) = 1;
        Op(Occ,OccP) = -1;
        }
    else
    if(opname == "projEmp")
        {
        Op(Emp,EmpP) = 1;
        }
    else
    if(opname == "projOcc")
        {
        Op(Occ,OccP) = 1; 
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

IQTensor inline Spinless::
makeN(int i) const
    {
    IQTensor N(conj(si(i)),siP(i));
    N(Occ(i),OccP(i)) = 1;
    return N;
    }

IQTensor inline Spinless::
makeC(int i) const
    {
    IQTensor C(conj(si(i)),siP(i));
    C(Occ(i),EmpP(i)) = 1;
    return C;
    }

IQTensor inline Spinless::
makeCdag(int i) const
    {
    IQTensor Cdag(conj(si(i)),siP(i));
    Cdag(Emp(i),OccP(i)) = 1;
    return Cdag;
    }

IQTensor inline Spinless::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Occ(i),OccP(i)) = -1;
    return fermiPhase;
    }

IQTensor inline Spinless::
makeProjEmp(int i) const
    {
    IQTensor projEmp(conj(si(i)),siP(i));
    projEmp(Emp(i),EmpP(i)) = 1;
    return projEmp;
    }

IQTensor inline Spinless::
makeProjOcc(int i) const
    {
    IQTensor projOcc(conj(si(i)),siP(i));
    projOcc(Occ(i),OccP(i)) = 1;
    return projOcc;
    }

#endif

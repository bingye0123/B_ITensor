//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TJ_H
#define __ITENSOR_TJ_H
#include "model.h"

class tJ : public Model
    {
    public:

    tJ();

    tJ(int N);

    tJ(std::ifstream& s) { doRead(s); }

    IQIndexVal
    Emp(int i) const;

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Dn(int i) const;

    IQIndexVal
    EmpP(int i) const;

    IQIndexVal
    UpP(int i) const;

    IQIndexVal
    DnP(int i) const;
        
    IQTensor
    getOp(int i, const String& opname, const OptSet& opts) const;
        
    IQIndexVal
    getState(int i, const String& state) const;


    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();
        
    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

inline tJ::
tJ()
    : N_(-1)
    { }

inline tJ::
tJ(int N)
    : N_(N),
      site_(N_+1)
    { 
    constructSites();
    }

inline void tJ::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        site_.at(j) = IQIndex(nameint("tJ site=",j),
            Index(nameint("Emp for site ",j),1,Site),  QN( 0,0,0),
            Index(nameint("Up for site ",j),1,Site),   QN(+1,1,1),
            Index(nameint("Dn for site ",j),1,Site),   QN(-1,1,1));
    }

inline void tJ::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void tJ::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int tJ::
getN() const
    { return N_; }

inline const IQIndex& tJ::
getSi(int i) const
    { return site_.at(i); }

inline IQIndex tJ::
getSiP(int i) const
{ return primed(site_.at(i)); }

IQIndexVal inline tJ::
getState(int i, const String& state) const
    {
    if(state == "0" || state == "Emp") 
        {
        return getSi(i)(1);
        }
    else 
    if(state == "+" || state == "Up") 
        {
        return getSi(i)(2);
        }
    else 
    if(state == "-" || state == "Dn") 
        {
        return getSi(i)(3);
        }
    else
        {
        Error("State " + state + " not recognized");
        return IQIndexVal();
        }
    }

inline IQIndexVal tJ::
Emp(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal tJ::
Up(int i) const
    {
    return getSi(i)(2);
    }

inline IQIndexVal tJ::
Dn(int i) const
    {
    return getSi(i)(3);
    }

inline IQIndexVal tJ::
EmpP(int i) const
    {
    return primed(getSi(i))(1);
    }

inline IQIndexVal tJ::
UpP(int i) const
    {
    return primed(getSi(i))(2);
    }

inline IQIndexVal tJ::
DnP(int i) const
    {
    return primed(getSiP(i))(3);
    }

inline IQTensor tJ::
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Em(s(1)),
               EmP(sP(1)),
               Up(s(2)),
               UpP(sP(2)),
               Dn(s(3)),
               DnP(sP(3));

    IQTensor Op(conj(s),sP);

    if(opname == "TReverse")
        {
        Op(Em,EmP) = +1;
        Op(Dn,UpP) = -1; //correct?
        Op(Up,DnP) = +1;
        }
    else
    if(opname == "Nup")
        {
        Op(Up,UpP) = 1;
        }
    else
    if(opname == "Ndn")
        {
        Op(Dn,DnP) = 1;
        }
    else
    if(opname == "Ntot")
        {
        Op(Up,UpP) = 1;
        Op(Dn,DnP) = 1;
        }
    else
    if(opname == "Cup")
        {
        Op(Up,EmP) = 1; 
        }
    else
    if(opname == "Cdagup")
        {
        Op(Em,UpP) = 1; 
        }
    else
    if(opname == "Cdn")
        {
        Op(Dn,EmP) = 1; 
        }
    else
    if(opname == "Cdagdn")
        {
        Op(Em,DnP) = 1; 
        }
    else
    if(opname == "Aup")
        {
        Op(Up,EmP) = 1; 
        }
    else
    if(opname == "Adagup")
        {
        Op(Em,UpP) = 1; 
        }
    else
    if(opname == "Adn")
        {
        Op(Dn,EmP) = 1; 
        }
    else
    if(opname == "Adagdn")
        {
        Op(Em,DnP) = 1; 
        }
    else
    if(opname == "FermiPhase" || opname == "F")
        {
        Op(Em,EmP) = +1; 
        Op(Up,UpP) = -1;
        Op(Dn,DnP) = -1;
        }
    else
    if(opname == "Sz")
        {
        Op(Up,UpP) = +0.5; 
        Op(Dn,DnP) = -0.5;
        }
    else
    if(opname == "Sx")
        {
        Op(Up,DnP) = 1; 
        Op(Dn,UpP) = 1;
        }
    else
    if(opname == "Sp")
        {
        Op(Dn,UpP) = 1;
        }
    else
    if(opname == "Sm")
        {
        Op(Up,DnP) = 1;
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }

#endif

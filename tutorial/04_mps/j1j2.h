#ifndef __J1J2_H
#define __J1J2_H

#include "core.h"
#include "sites/spinhalf.h"
#include "hams/J1J2Chain.h"

namespace itensor {

MPS inline
computeGroundState(const SpinHalf& sites, Real J2)
    {
    MPO H = J1J2Chain(sites,Opt("J2",J2));

    MPS psi(sites);

    Sweeps sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-9;

    println("Starting ground state calculation for J2 = ",J2);

    dmrg(psi,H,sweeps,"Quiet");

    println("Done with ground state calculation.");

    return psi;
    }

}; //namespace itensor

#endif

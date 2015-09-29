//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INTEGRATORS_H
#define __ITENSOR_INTEGRATORS_H

#include "global.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// 4th Order Runge-Kutta.
// Assumes a time-independent Force.
// The vector v should be 1 indexed.
//
template <class Tensor, typename Deriv>
void
rungeKutta4(const Deriv& D, Real tstep, std::vector<Tensor>& v, 
            const OptSet& opts = Global::opts())
    {
    int N = int(v.size())-1;
    while(v.at(N).isNull() && N > 1) --N;

    if(N <= 0) Error("Empty vector v (v should be 1-indexed)");

    std::vector<Tensor> k1,k2,k3,k4;

    k1 = D(v);

    //d = v + (tstep/2)*k1
    std::vector<Tensor> d(v);
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k1.at(j);
        }
    k2 = D(d);

    //d = v + (tstep/2)*k2
    d = v;
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k2.at(j);
        }
    k3 = D(d);

    //d = v + (tstep)*k3
    d = v;
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep)*k3.at(j);
        }
    k4 = D(d);


    for(int j = 1; j <= N; ++j)
        {
        v.at(j) += (tstep/6.)*(k1.at(j) + 2*k2.at(j) + 2*k3.at(j) + k4.at(j));
        }
    }

template <class Tensor, typename Deriv>
void
midpointMethod(const Deriv& D, Real tstep, std::vector<Tensor>& v, 
               const OptSet& opts = Global::opts())
    {
    int N = int(v.size())-1;
    while(v.at(N).isNull() && N > 1) --N;

    if(N <= 0) Error("Empty vector v (v should be 1-indexed)");

    std::vector<Tensor> k1,k2,k3,k4;

    k1 = D(v);

    std::vector<Tensor> d(v);
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k1.at(j);
        }
    k2 = D(d);

    for(int j = 1; j <= N; ++j)
        {
        v.at(j) += tstep*k2.at(j);
        }
    }

#undef Cout
#undef Endl
#undef Format

#endif

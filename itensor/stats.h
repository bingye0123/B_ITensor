//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STATS_H
#define __ITENSOR_STATS_H

#include "real.h"

class Stats
    {
    public:

    std::vector<Real> dat;
    Stats() { dat.reserve(100); }

    void putin(Real x) { dat.push_back(x); }

    Real 
    avg() const
        {
        if(dat.empty()) return 0;
        Real s = 0;
        for(size_t j = 0; j < dat.size(); ++j) s += dat.at(j);
        return s/dat.size();
        }

    Real 
    sigma() const
        {
        if(dat.size() < 2) return 0;
        Real av = avg(), s = 0;
        for(size_t j = 0; j < dat.size(); ++j) s += sqr(dat.at(j)-av);
        return sqrt(s/dat.size());
        }

    Real 
    err() const
        {
        Real n = dat.size();
        if(n > 1) return sigma()/sqrt(n-1);
        else return sigma()/sqrt(n);
        }

    Real 
    binerr(int bs) const
        {
        int n = dat.size();
        if(n < bs) return 0;
        int nbin = n / bs;
        Vector bin(nbin);
        bin = 0.0;
        for(int i = 0; i < n; i++)
            {
            if(i/bs > (nbin-1)) break;
            bin.el(i/bs) += dat.at(i);
            }
        bin *= 1.0/bs;
        Stats bmeas;
        for(int i = 1; i <= nbin; ++i) bmeas.putin(bin(i));
        return bmeas.err();
        }

    Real 
    bootstrap(int nresamples = 1000)
        { //Uses the bootstrap procedure to estimate the std deviation
        Stats boot;
        int n = dat.size();
        if(n < 2) return 0.0;
        for(int sample = 1; sample <= nresamples; ++sample)
            {
            Real avg = 0.0;
            for(int i = 1; i <= n; ++i)
                {
                int which = int(ran1()*n);
                avg += dat.at(which);
                }
            avg /= n;
            boot.putin(avg);
            }
        return boot.sigma();
        }

    Real 
    bootstrap_2c(Real prefactor,Stats firstc, bool correlated = true, int nresamples = 1000)
        { //Uses the bootstrap procedure to estimate the second cumulant error
        Stats boot2c;
        Real n = dat.size();
        if(n < 2) return 0.0;
        if(firstc.dat.size() != n) error("Stats instances don't have the same size data sets.");
        for(int sample = 1; sample <= nresamples; ++sample)
            {
            Real avg2 = 0.0;
            Real avg = 0.0;
            int which;
            for(int i = 1; i <= n; ++i)
                { //Sample with replacement (replacement means 'which' can take the same value more than once)
                which = int(ran1()*n);
                avg2 += dat.at(which);
                if(!correlated) which = int(ran1()*n); //use <H^2> and <H> from the same METTS
                avg += firstc.dat.at(which);
                }
            avg2 /= n; avg /= n;
            boot2c.putin(prefactor*(avg2-avg*avg));
            }
        return boot2c.sigma();
        }

    };

#endif

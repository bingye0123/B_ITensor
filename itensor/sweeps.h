//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SWEEPS_HEADER_H
#define __ITENSOR_SWEEPS_HEADER_H
#include "global.h"
#include "input.h"
#include "boost/function.hpp"

template <typename T>
class SweepSetter;

//
// To use the InputGroup / table constructor,
// the format required in the input file is:
//
// sweep_table_name
//      {
//      maxm   minm  cutoff  niter  noise
//      20     20    1E-8    4      1E-8
//      40     20    1E-8    3      1E-9
//      80     20    1E-10   2      1E-10
//      160    20    1E-12   2      0
//      240    20    1E-12   2      0
//      }
//
class Sweeps
    {
    public:

    //Constructors --------------

    Sweeps();

    Sweeps(int nsweeps, 
           int minm = 1, 
           int maxm = 500, 
           Real cutoff = 1E-8);

    Sweeps(int nsweeps, InputGroup& sweep_table);
    
    //Accessor methods ----------

    int 
    nsweep() const { return nsweep_; }
    void 
    nsweep(int val);

    int 
    minm(int sw) const { return minm_.at(sw); }
    void 
    setminm(int sw, int val) { minm_.at(sw) = val; }

    //Use as sweeps.minm() = 20,20,10; (all remaining set to 10)
    SweepSetter<int> 
    minm();

    int 
    maxm(int sw) const { return maxm_.at(sw); }
    void 
    setmaxm(int sw, int val) { maxm_.at(sw) = val; }

    //Use as sweeps.maxm() = 50,50,100,100,500; (all remaining set to 500)
    SweepSetter<int> 
    maxm();

    Real 
    cutoff(int sw) const { return cutoff_.at(sw); }
    void 
    setcutoff(int sw, Real val) { cutoff_.at(sw) = val; }

    //Use as sweeps.cutoff() = 1E-8; (cutoff set to 1E-8 for all sweeps)
    //or as sweeps.cutoff() = 1E-4,1E-4,1E-5,1E-5,1E-10; (all remaining set to 1E-10)
    SweepSetter<Real> 
    cutoff();

    Real 
    noise(int sw) const { return noise_.at(sw); }
    void 
    setnoise(int sw, Real val) { noise_.at(sw) = val; }
    void 
    setnoise(Real val) { noise_.assign(nsweep_+1,val); }

    //Use as sweeps.noise() = 1E-10; (noise set to 1E-10 for all sweeps)
    //or as sweeps.noise() = 1E-8,1E-9,1E-10,0.0; (all remaining set to 0)
    SweepSetter<Real> 
    noise();

    int 
    niter(int sw) const { return niter_.at(sw); }
    void 
    setniter(int sw, int val) { niter_.at(sw) = val; }
    void 
    setniter(int val) { niter_.assign(nsweep_+1,val); }

    //Use as sweeps.niter() = 5,4,3,2; (all remaining set to 2)
    SweepSetter<int> 
    niter();

    private:

    void 
    init(int min_m, int max_m, Real cut);

    void 
    tableInit(InputGroup& table);

    std::vector<int> maxm_,
                     minm_,
                     niter_;
    std::vector<Real> cutoff_,
                      noise_;
    int nsweep_;
    };

//
// Helper class for Sweeps accessor methods.
// Accumulates a comma separated list of 
// values of type T, storing them in the
// vector v passed to its constructor.
//
template <typename T>
class SweepSetter
    {
    public:

    SweepSetter(std::vector<T>& v)
        :
        nsweep_(int(v.size())-1),
        started_(false),
        v_(v),
        j_(1)
        { 
        last_val_ = v_[j_];
        }
    
    ~SweepSetter()
        {
        for(; j_ <= nsweep_; ++j_)
            {
            v_[j_] = last_val_;
            }
        }

    template <typename Arg>
    SweepSetter& 
    operator=(Arg val)
        {
        started_ = true;
        return operator,(val);
        }

    SweepSetter& 
    operator,(T val)
        {
        checkStarted();
        if(j_ > nsweep_) return *this;
        v_[j_] = val;
        ++j_;
        last_val_ = val;
        return *this;
        }

    SweepSetter& 
    operator,(const Opt& opt)
        { 
        checkStarted();
        if(opt.name() == "Repeat")
            {
            if(j_ == 1) Error("No value to repeat");
            for(int n = 1; n < opt.intVal(); ++n, ++j_)
                {
                if(j_ > nsweep_) return *this;
                v_[j_] = last_val_;
                }
            }
        return *this;
        }

    typedef boost::function2<T,int,int> 
    Func;

    //
    //Sets the remaining values of v_ by 
    //accepting (the address of a) function 
    //with signature:
    //
    //  T f(int sw, int nsweep) { ... }
    //
    //or an instance of a function object
    //which defines a method:
    //
    //  T operator()(int sw, int nsweep) const { ... } 
    //
    SweepSetter&
    operator,(Func f)
        {
        checkStarted();
        for(; j_ <= nsweep_ ; ++j_)
            {
            v_[j_] = f(j_,nsweep_);
            }
        last_val_ = v_.back();
        return *this;
        }

    private:

    const int nsweep_;
    bool started_;
    std::vector<T>& v_;
    int j_;
    T last_val_;

    void
    checkStarted() const
        {
        if(!started_) Error("SweepSetter notation is setter() = #,#,...;");
        }
    };

struct RampM
    {
    RampM(int start_m, int end_m,
          const OptSet& opts = Global::opts())
        :
        start_m_(start_m),
        end_m_(end_m),
        nwarm_(opts.getInt("Warmup",-1))
        { }

    int
    operator()(int sw, int nsweep) const 
        { 
        const int actual_nwarm = (nwarm_ < 0 ? nsweep : min(nwarm_+1,nsweep));
        if(sw <= actual_nwarm)
            return (int) (start_m_ + (sw-1.)/(actual_nwarm-1.)*(end_m_-start_m_));
        else
            return end_m_;
        }

    private:
    const
    int start_m_,
        end_m_,
        nwarm_;
    };

struct ExpM
    {
    ExpM(int start_m, int end_m,
         const OptSet& opts = Global::opts())
        :
        start_m_(start_m),
        end_m_(end_m),
        exp_base_(opts.getReal("ExpBase",2.))
        { }

    int
    operator()(int sw, int nsweep) const 
        { 
        int expm = start_m_*int(pow(exp_base_,sw-1));
        if(expm <= 0) return end_m_; //catch overflow
        return min(expm,end_m_);
        }

    private:
    const
    int start_m_,
        end_m_;
    const
    Real exp_base_;
    };

inline Sweeps::
Sweeps()
    :
    nsweep_(0)
    { }

inline Sweeps::
Sweeps(int nsw, int min_m, int max_m, Real cut)
    :
    nsweep_(nsw)
    {
    init(min_m,max_m,cut);
    }

inline Sweeps::
Sweeps(int nsw, InputGroup& sweep_table)
    : 
    nsweep_(nsw)
    {
    tableInit(sweep_table);
    }

SweepSetter<int> inline Sweeps::
minm() { return SweepSetter<int>(minm_); }

SweepSetter<int> inline Sweeps::
maxm() { return SweepSetter<int>(maxm_); }

SweepSetter<Real> inline Sweeps::
cutoff() { return SweepSetter<Real>(cutoff_); }

SweepSetter<Real> inline Sweeps::
noise() { return SweepSetter<Real>(noise_); }

SweepSetter<int> inline Sweeps::
niter() { return SweepSetter<int>(niter_); }

void inline Sweeps::
nsweep(int val)
    { 
    if(val > nsweep_) 
        Error("Can't use nsweep accessor to increase number of sweeps.");
    nsweep_ = val; 
    }


void inline Sweeps::
init(int min_m, int max_m, Real cut)
    {
    minm_ = std::vector<int>(nsweep_+1,min_m);
    maxm_ = std::vector<int>(nsweep_+1,max_m);
    cutoff_ = std::vector<Real>(nsweep_+1,cut);
    niter_ = std::vector<int>(nsweep_+1,2);
    noise_ = std::vector<Real>(nsweep_+1,0);

    //Set number of Davidson iterations
    const int Max_niter = 9;
    for(int s = 1; s <= min(4,nsweep_); ++s)
        {
        niter_.at(s) = max(Max_niter-s+1,2);
        }

    } //Sweeps::init

void inline Sweeps::
tableInit(InputGroup& table)
    {
    if(!table.GotoGroup()) 
        Error("Couldn't find table " + table.name);

    minm_ = std::vector<int>(nsweep_+1);
    maxm_ = std::vector<int>(nsweep_+1);
    cutoff_ = std::vector<Real>(nsweep_+1);
    niter_ = std::vector<int>(nsweep_+1);
    noise_ = std::vector<Real>(nsweep_+1);

    table.SkipLine(); //SkipLine so we can have a table key
    for(int i = 1; i <= nsweep_; i++)
        {
        table.infile.file >> maxm_[i] >> minm_[i] >> cutoff_[i] >> niter_[i] >> noise_[i];
        }

    } //Sweeps::tableInit

inline std::ostream&
operator<<(std::ostream& s, const Sweeps& swps)
    {
    s << "Sweeps:\n";
    for(int sw = 1; sw <= swps.nsweep(); ++sw)
        {
        s << boost::format("%d  Maxm=%d, Minm=%d, Cutoff=%.1E, Niter=%d, Noise=%.1E\n")
             % sw % swps.maxm(sw) % swps.minm(sw) % swps.cutoff(sw) %swps.niter(sw) % swps.noise(sw);
        }
    return s;
    }

void inline
sweepnext(int &b, int &ha, int N, int min_b = 1)
    {
    const int inc = (ha==1 ? +1 : -1);
    b += inc;
    if(b == (ha==1 ? N : min_b-1))
        {
        b -= inc;
        ++ha;
        }
    }


#endif //__ITENSOR_SWEEPS_HEADER_H

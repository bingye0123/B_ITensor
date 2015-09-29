#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/Heisenberg.h"
using boost::format;
using namespace std;

//typedef SpinHalf
//Spin;           //use S=1/2 degrees of freedom

//Un-comment above typedef and comment this one to switch spin type
typedef SpinOne
Spin;             //use S=1 degrees of freedom

int 
main(int argc, char* argv[])
    {
    int N = 100;

    //
    // Initialize the site degrees of freedom.
    //
    Spin model(N);    // make a chain of N spins

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    MPO H = Heisenberg(model);

    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(model);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    MPS psi(initState);

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    cout << format("Initial energy = %.5f") % psiHphi(psi,H,psi) << endl;

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    cout << sweeps;

    //
    // Begin the DMRG calculation
    //

    Real En = dmrg(psi,H,sweeps,Quiet());

    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;

    return 0;
    }

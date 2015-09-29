
//#include "symMonteCarlo.h"
#include <mpi.h>
#include "core.h"
#include "itensor/model/spinone.h"
#include "itensor/hams/Heisenberg.h"
//#include "itensor/itensor.h"

using boost::format;
using namespace std;
//using namespace itensor;

//typedef tJ
//typedef Hubbard
//states;


int
main(int argc, char* argv[])
{
    MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    
    int N=8;
    
    SpinOne model(N);
    MPS psi(model);
    
    cout<<psi.A(1)<<endl;
    //cout<<"Hello!"<<endl;
    
    
    return 0;
}



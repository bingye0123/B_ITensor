
#include "symMonteCarlo.h"

#include "core.h"

#include "model/tjmod.h"
#include "hams/tJham.h"
//#include "model/hubbardmod.h"
//#include "hams/ExtendedHubbardham.h"

using boost::format;
using namespace std;

typedef tJ
//typedef Hubbard
states;



int
main(int argc, char* argv[])
{
    int N = 32;

	states model(N);
	readFromFile("model_dat",model);
	IQMPS psi(model);
	//readFromFile("psi_dat",psi);
	psi.read("psi_d9qvXn");
	//read("PH_p4u8mi");
    //const IQMPO H = tJChain(model,Opt("J",0.5)& Opt("t",1.0));
	IQMPO H=tJChain(model);
	readFromFile("Ham_dat",H);
	
	int Nsweeps=4;
	Sweeps sweeps(Nsweeps);
	sweeps.maxm() = 10000,10000,10000,10000;
	sweeps.cutoff() =0.;
	sweeps.niter() = 2;
	sweeps.noise() = 1E-7,1E-7,0.;
	    cout << sweeps;

    psi.doWrite(true);
    H.doWrite(true);
    Real En = dmrg(psi,H,sweeps,Quiet());
    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;

    psi.doWrite(false);
    H.doWrite(false);
      cout<<"Writing to disk:"<<endl;
    	writeToFile("psiJ05_dat",psi);
	writeToFile("modelJ05_dat",model);
	writeToFile("HamJ05_dat",H);
      cout<<"Writing done."<<endl;
    
    
    //cout  <<"Energy: " << psiHphi(psi,H,psi)<<endl;

	psi.position(1);
	cout << "Norm2: " << Dot(conj(psi.A(1)), psi.A(1))<<endl;
	
	ITdVec nup_vec(N),ndn_vec(N);
	
	for(int j=1; j<=N; j++) {
		IQTensor nup_op = model.op("Nup",j);
	        IQTensor ndn_op = model.op("Ndn",j);
		psi.position(j);
		IQTensor ket = psi.A(j);
		IQTensor bra = conj(primed(ket,Site));
		nup_vec(j-1) = Dot(bra, nup_op*ket);
		ndn_vec(j-1) = Dot(bra, ndn_op*ket);
	}
	
	cout<<"Nup: "<<nup_vec<<endl;
	cout<<"Ndn: "<<ndn_vec<<endl;
	cout<<"Total N: "<<nup_vec+ndn_vec<<endl;
	cout<<"Sz: "<<0.5*(nup_vec-ndn_vec)<<endl;
	
	psi.position(1);

    
	
	
    return 0;
}



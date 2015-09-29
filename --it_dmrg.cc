
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

//number of electrons
const extern int Nferm=6;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm=500;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure=1000;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure=40;


int
main(int argc, char* argv[])
{
  MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process

    int N=8;
    
	//
    // Initialize the site degrees of freedom.
    //
    states model(N);    // make a chain of N spins
    
    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    const IQMPO H = tJChain(model,Opt("J",0.2)&Opt("t",1.0));
    //const IQMPO H = ExtendedHubbard(model,Opt("U",10.0)& Opt("t1",1)& Opt("t2",0.0)& Opt("V1",0));
	
    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(model);
	//up and down spins
    ITiVec init_state_vec=initialize_random_state_vec_tj(N,Nferm);
    //ITiVec init_state_vec=initialize_random_state_vec_hubbard(N,Nferm);
	
	if (rank==0) {
		cout<<"initial state: "<<init_state_vec<<endl;
	}

    for(int i=1;i<=N;i++)
    {
		if(init_state_vec(i-1)==1)
		{
			initState.set(i,"Emp");
		}
		else if(init_state_vec(i-1)==2)
		{
			initState.set(i,"Up");
		}
		else if(init_state_vec(i-1)==3)
		{
			initState.set(i,"Dn");
		}
		else if(init_state_vec(i-1)==4)
		{
			initState.set(i,"UpDn");
		}
		else
		{
			cout<<"InitState Error"<<endl;
			exit(0);
		}
    }
    
    IQMPS psi(initState);
	
    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
	if (rank==0) {
		cout << totalQN(psi) << endl;
		cout << format("Initial energy = %.5f") % psiHphi(psi,H,psi) << endl;
	}
	
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
	   
	int Nsweeps=45;
	Sweeps sweeps(Nsweeps);
	sweeps.maxm() = 10,20,20,20,20,20,50,50,50,100,100,100,100,200,200,200,200,200,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,4000,4000,4000,4000,4000,4000,4000,4000,4000;
	sweeps.cutoff() = 1E-10,1E-10,0.;
	sweeps.niter() = 2;
	sweeps.noise() = 1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,0.;
	cout << sweeps;
	
    //
    // Begin the DMRG calculation
    //
	
	
    Real En = dmrg(psi,H,sweeps,Quiet());
    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;
    cout  <<"Energy: " << psiHphi(psi,H,psi)<<endl;

    

	psi.position(1);
	if (rank==0) { cout << "Norm2: " << Dot(conj(psi.A(1)), psi.A(1))<<endl; }
	cout.flush();

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
	
	if (rank==0) {
	cout<<"Nup: "<<nup_vec<<endl;
	cout<<"Ndn: "<<ndn_vec<<endl;
	cout<<"Total N: "<<nup_vec+ndn_vec<<endl;
	cout<<"Sz: "<<0.5*(nup_vec-ndn_vec)<<endl;
	}
	
	psi.position(1);


	
	ITiVec inv_sym_mapping_table(N);
	//inversion: site -> N-1-site
	for(int i=0;i<inv_sym_mapping_table.size();i++)
	{
		inv_sym_mapping_table(i)=N-1-i;
	}
	
	ITiVec T1_sym_mapping_table_8_site="2 3 0 1 6 7 4 5";
	ITiVec T2_sym_mapping_table_8_site="4 5 6 7 0 1 2 3";
	ITiVec C6_sym_mapping_table_8_site="1 2 7 4 3 0 5 6";
	ITiVecArray C6_generator_mapping_table(1);
	C6_generator_mapping_table(0)=C6_sym_mapping_table_8_site;
	ITcVec Ang_mom(1);
	Ang_mom(0)=std::exp(Complex_i*2.*Pi/3.);
	//Ang_mom(0)=1.;
	
	ITiVec T1_sym_mapping_table_32_site="3 4 5 6 7 8 1 2 11 12 13 14 15 16 9 10 19 20 21 22 23 24 17 18 27 28 29 30 31 32 25 26";
	T1_sym_mapping_table_32_site=T1_sym_mapping_table_32_site-1;
	ITiVec T2_sym_mapping_table_32_site="25 26 27 28 29 30 31 32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24";
	T2_sym_mapping_table_32_site=T2_sym_mapping_table_32_site-1;
	ITiVec C6_sym_mapping_table_32_site="4 5 14 15 24 17 26 27 2 3 12 13 22 23 32 25 8 1 10 11 20 21 30 31 6 7 16 9 18 19 28 29";
	C6_sym_mapping_table_32_site=C6_sym_mapping_table_32_site-1;
	
	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_8_site;
	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_8_site;
	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_8_site;
//	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_32_site;
//	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_32_site;
//	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_32_site;

	ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table;
	translation_generator_mapping_table(1)=T2_sym_mapping_table;
	ITcVec COM_mom="1. 1.";
	
	ITiVecArray transl_C6_generator_mapping_table(3);
	transl_C6_generator_mapping_table(0)=T1_sym_mapping_table;
	transl_C6_generator_mapping_table(1)=T2_sym_mapping_table;
	transl_C6_generator_mapping_table(2)=C6_sym_mapping_table;
	
	ITcVec COM_Ang_mom="1. 1. 0.";
	COM_Ang_mom(2)=std::exp(-Complex_i*2.*Pi/3.);
	
	std::string desc("8_J02_m2P3");
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table,translation_generator_mapping_table,COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_tj_MonteCarlo(psi,C6_generator_mapping_table,Ang_mom)<<endl;
	
	cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_hubbard_MonteCarlo(psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_hubbard_MonteCarlo(psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi, C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, C6_sym_mapping_table_8_site, C6_generator_mapping_table, Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,T1_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, T1_sym_mapping_table_8_site, translation_generator_mapping_table, COM_mom)<<endl;
	
	ITiVec triv_sym_mapping_table="0 1 2 3 4 5 6 7";
	ITiVecArray triv_generator_mapping_table(1);
	triv_generator_mapping_table(0)=triv_sym_mapping_table;
	ITcVec triv_mom="1.";
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,triv_generator_mapping_table,triv_mom)<<endl;

	
    MPI::Finalize();
    
    return 0;
}




#include "symMonteCarlo.h"

#include "core.h"

#include "model/tjmod.h"
#include "hams/tJham.h"
//#include "model/hubbardmod.h"
//#include "hams/ExtendedHubbardham.h"

//LATTICE=0 indicates Honeycomb lattice
//LATTICE=1 indicates Triangular lattice
//#define LATTICE 1
//SYSSHAPE=0 indicates system has rhombus shape
//SYSSHAPE=1 indicates system has hexagonal shape
//#define SYSSHAPE 0
#define RANDOM 1
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the central row to edge rows
#define NUMBERING 0


using boost::format;
using namespace std;

typedef tJ
//typedef Hubbard
states;


//number of electrons
const extern int Nferm=14;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm=500;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure=1000;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure=40;


int
main(int argc, char* argv[])
{
    double startwtime = 0.0,endwtime;
    MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process

    
    if(rank == 0)
        startwtime = MPI_Wtime();
    int N=28;
    
	//
    // Initialize the site degrees of freedom.
    //
    states model(N);    // make a chain of N spins
    
    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    IQMPO H = tJChain(model,Opt("J",2.0)&Opt("t",-1.0));
    //const IQMPO H = ExtendedHubbard(model,Opt("U",10.0)& Opt("t1",1)& Opt("t2",0.0)& Opt("V1",0));
	
    
    
    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(model);
	//up and down spins
    ITiVec init_state_vec=initialize_random_state_vec_tj(N,Nferm);
    //ITiVec init_state_vec=initialize_random_state_vec_hubbard(N,Nferm);
    if (RANDOM==0) {
        for (int i=0; i<=17; i++) {
            init_state_vec(i)=1;
        }
        for (int i=18; i<=35; i++) {
            init_state_vec(i)=1;
        }
        /*
        for (int i=0; i<=16; i+=2) {
            init_state_vec(i+18)=2;
            init_state_vec(i+19)=3;
        };
         */
        
        for (int i=0; i<=10; i+=2) {
            init_state_vec(i)=2;
            init_state_vec(i+1)=3;
        };
        for (int i=30; i<=34; i+=2) {
            init_state_vec(i)=2;
            init_state_vec(i+1)=3;
        }
         
    }
	
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
		cout << format("Initial energy = %.10f") % psiHphi(psi,H,psi) << endl;
	}
	
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
	   
	int Nsweeps=70;
	Sweeps sweeps(Nsweeps);
	sweeps.maxm() = 10,20,20,20,20,20,50,50,50,100,100,100,100,200,200,200,200,200,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,4000,4000,4000,4000,4000,4000,4000,4000,4000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000;
	sweeps.cutoff() = 1E-10,1E-10,0.;
	sweeps.niter() = 2;
	sweeps.noise() = 1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.;
	cout << sweeps;
        
    //
    // Begin the DMRG calculation
    //
	psi.doWrite(true);
    H.doWrite(true);
	
    Real En = dmrg(psi,H,sweeps,Quiet());
    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;
    cout  <<"Energy: " << psiHphi(psi,H,psi)<<endl;

    psi.doWrite(false);
    H.doWrite(false);

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
        cout<<"Nup: "<<setprecision(10)<<nup_vec<<endl;
        cout<<"Ndn: "<<setprecision(10)<<ndn_vec<<endl;
        cout<<"Total N: "<<setprecision(10)<<nup_vec+ndn_vec<<endl;
        cout<<"Sz: "<<setprecision(10)<<0.5*(nup_vec-ndn_vec)<<endl;
	}
	
	psi.position(1);


    cout<<"Writing to disk:"<<endl;
  	writeToFile("psi_dat",psi);
	writeToFile("model_dat",model);
	writeToFile("Ham_dat",H);
    cout<<"Writing done."<<endl;
    
    ITiVec inv_sym_mapping_table, T1_sym_mapping_table_28_site, T2_sym_mapping_table_28_site, C6_sym_mapping_table_28_site;
	if (NUMBERING==0) {
        inv_sym_mapping_table="0 27 26 25 24 23 6 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 22 5 4 3 2 1";
        T1_sym_mapping_table_28_site="1 17 3 4 5 6 23 8 9 10 11 27 13 14 15 16 2 18 19 20 21 22 7 24 25 26 12 0";
        T2_sym_mapping_table_28_site="24 25 11 27 0 1 17 2 3 4 5 6 22 7 8 9 10 26 12 13 14 15 16 18 19 20 21 23";
        C6_sym_mapping_table_28_site="22 2 17 12 7 3 0 18 13 8 4 1 23 19 14 9 5 27 24 20 15 10 6 25 21 16 11 26";
    }
    if (NUMBERING==1) {
        inv_sym_mapping_table="0 27 26 25 24 23 6 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 22 5 4 3 2 1";
        T1_sym_mapping_table_28_site="1 17 3 4 5 6 23 8 9 10 11 27 13 14 15 16 2 18 19 20 21 22 7 24 25 26 12 0";
        T2_sym_mapping_table_28_site="24 25 11 27 0 1 17 2 3 4 5 6 22 7 8 9 10 26 12 13 14 15 16 18 19 20 21 23";
        C6_sym_mapping_table_28_site="22 2 17 12 7 3 0 18 13 8 4 1 23 19 14 9 5 27 24 20 15 10 6 25 21 16 11 26";
    }
	ITiVecArray C6_generator_mapping_table(1);
	C6_generator_mapping_table(0)=C6_sym_mapping_table_28_site;
	ITcVec Ang_mom(1);
	Ang_mom(0)=std::exp(Complex_i*2.*Pi/3.);
	//Ang_mom(0)=1.;
	
    
	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_28_site;
	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_28_site;
	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_28_site;
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
	
	std::string desc("28_J2_m2P3");
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table,translation_generator_mapping_table,COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_tj_MonteCarlo(psi,C6_generator_mapping_table,Ang_mom)<<endl;
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
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
	
	ITiVec triv_sym_mapping_table="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27";
	ITiVecArray triv_generator_mapping_table(1);
	triv_generator_mapping_table(0)=triv_sym_mapping_table;
	ITcVec triv_mom="1.";
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,triv_generator_mapping_table,triv_mom)<<endl;

    if(rank==0){
		endwtime=MPI_Wtime();
		printf("wall clock time = %f \n", endwtime-startwtime);
	}
	
    MPI::Finalize();
    
    return 0;
}



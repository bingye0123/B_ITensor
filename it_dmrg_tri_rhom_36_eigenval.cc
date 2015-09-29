
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
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the central row to edge rows
#define NUMBERING 1
#define DO_NORM 1


using boost::format;
using namespace std;

typedef tJ
//typedef Hubbard
states;


//number of electrons
const extern int Nferm=18;
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
    int N=36;
    //cout<<"loading info"<<endl;
    if (rank==0) { cout<<"Reading from disk:"<<endl; }
    states model(N);
	readFromFile("model_dat",model);
	IQMPS psi(model);
	readFromFile("psi_dat",psi);
	//psi.read("psi_d9qvXn");
	//read("PH_p4u8mi");
    //const IQMPO H = tJChain(model,Opt("J",0.5)& Opt("t",1.0));
	IQMPO H=tJChain(model);
	readFromFile("Ham_dat",H);
    if (rank==0) { cout<<"Reading done."<<endl; }
	
    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
	if (rank==0) {
		cout << totalQN(psi) << endl;
		cout << format("Initial energy from loaded psi = %.10f") % psiHphi(psi,H,psi) << endl;
	}
	
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
	   
	//int Nsweeps=45;
	//Sweeps sweeps(Nsweeps);
	//sweeps.maxm() = 10,20,20,20,20,20,50,50,50,100,100,100,100,200,200,200,200,200,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,4000,4000,4000,4000,4000,4000,4000,4000,4000;
    //int Nsweeps=10;
	//Sweeps sweeps(Nsweeps);
    //sweeps.maxm() = 4000,4000,4000,4000,4000,8000,8000,8000,8000,8000;//,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000;
	//sweeps.cutoff() = 1E-10,1E-10,0.;
	//sweeps.niter() = 2;
	//sweeps.noise() = 1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,0.;
	//cout << sweeps;
    	//int Nsweeps=50;
	//Sweeps sweeps(Nsweeps);
	//sweeps.maxm() = 10,20,20,20,20,20,50,50,50,100,100,100,100,200,200,200,200,200,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000;
	//sweeps.cutoff() = 1E-10,1E-10,0.;
	//sweeps.niter() = 2;
	//sweeps.noise() = 1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,0.,0.,0.,0.,0.,0.;
	//cout << sweeps;
	
    //
    // Begin the DMRG calculation
    //
	
	
    //Real En = dmrg(psi,H,sweeps,Quiet());
    //
    // Print the final energy reported by DMRG
    //
    //cout << format("\nGround State Energy = %.10f")%En << endl;
    //cout  <<"Energy: " << psiHphi(psi,H,psi)<<endl;

    

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

    //cout<<"Writing to disk:"<<endl;
  	//writeToFile("psi_dat",psi);
	//writeToFile("model_dat",model);
	//writeToFile("Ham_dat",H);
    //cout<<"Writing done."<<endl;

	//cout<<"measuring norm now"<<endl;
    ITiVec inv_sym_mapping_table, T1_sym_mapping_table_36_site, T2_sym_mapping_table_36_site, C6_sym_mapping_table_36_site;
	if (NUMBERING==0) {
        inv_sym_mapping_table="24 29 28 27 26 25 18 23 22 21 20 19 12 17 16 15 14 13 6 11 10 9 8 7 0 5 4 3 2 1 30 35 34 33 32 31";
        T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
        T2_sym_mapping_table_36_site="30 31 32 33 34 35 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29";
        C6_sym_mapping_table_36_site="19 13 7 1 31 25 26 20 14 8 2 32 33 27 21 15 9 3 4 34 28 22 16 10 11 5 35 29 23 17 12 6 0 30 24 18";
    }
    if (NUMBERING==1) {
        inv_sym_mapping_table="0 5 4 3 2 1 12 17 16 15 14 13 6 11 10 9 8 7 24 29 28 27 26 25 18 23 22 21 20 19 30 35 34 33 32 31";
        T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
        T2_sym_mapping_table_36_site="12 13 14 15 16 17 0 1 2 3 4 5 24 25 26 27 28 29 6 7 8 9 10 11 30 31 32 33 34 35 18 19 20 21 22 23";
        C6_sym_mapping_table_36_site="33 21 9 3 15 27 28 34 22 10 4 16 20 8 2 14 26 32 17 29 35 23 11 5 7 1 13 25 31 19 0 12 24 30 18 6";
    }
	ITiVecArray C6_generator_mapping_table(1);
	C6_generator_mapping_table(0)=C6_sym_mapping_table_36_site;
	ITcVec Ang_mom(1);
	Ang_mom(0)=std::exp(Complex_i*2.*Pi/3.);
	//Ang_mom(0)=1.;
	
   
	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_36_site;
	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_36_site;
	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_36_site;
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
	COM_Ang_mom(2)=std::exp(Complex_i*2.*Pi/3.);
    
    ITiVec triv_sym_mapping_table="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35";
	ITiVecArray triv_generator_mapping_table(1);
	triv_generator_mapping_table(0)=triv_sym_mapping_table;
	ITcVec triv_mom="1.";
	
	std::string desc("36");
	
    if (DO_NORM==1) {
        ITcVec mom="1. 1.";
        int count=0;
        for (int p=0; p<6; p++) {
            for (int q=0; q<6; q++) {
                mom(0)=std::exp(Complex_i*(double) p*Pi/3.);
                mom(1)=std::exp(Complex_i*(double) q*Pi/3.);
                //cout<<mom<<endl;
                cout<<count<<"k "<<mom<<"   "<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_norm_measurement_tj_MonteCarlo(desc, psi, translation_generator_mapping_table, mom)<<endl;
                count++;
            }
        }
    }
     
    
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(desc,psi,C6_sym_mapping_table,translation_generator_mapping_table,COM_mom)<<endl;
    //cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_norm_measurement_tj_MonteCarlo(desc, psi, translation_generator_mapping_table, COM_mom)<<endl;
    //cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_norm_measurement_tj_MonteCarlo(desc, psi, triv_generator_mapping_table, triv_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
    //////cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(desc,psi,T2_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_tj_MonteCarlo(psi,C6_generator_mapping_table,Ang_mom)<<endl;
    //cout<<"From rank="<<rank<<" Monte Carlo measure:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	
	///////cout<<"From rank="<<rank<<" Monte Carlo measur e:"<<symmetry_sector_measurement_tj_MonteCarlo(desc,psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_hubbard_MonteCarlo(psi,transl_C6_generator_mapping_table,COM_Ang_mom)<<endl;
	
	
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_hubbard_MonteCarlo(psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	
	//////cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(desc,psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(desc,psi, inv_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_tj_MonteCarlo(desc,psi, C6_sym_mapping_table_16_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(desc,psi, C6_sym_mapping_table, C6_generator_mapping_table, Ang_mom)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,C6_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_eigval_hubbard_MonteCarlo(psi,T1_sym_mapping_table_8_site)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval: "<<symmetry_sector_eigval_tj_MonteCarlo(psi, T1_sym_mapping_table_8_site, translation_generator_mapping_table, COM_mom)<<endl;
    
    if(rank==0){
		endwtime=MPI_Wtime();
		printf("wall clock time = %f \n", endwtime-startwtime);
	}
	
    MPI::Finalize();
    
    return 0;
}



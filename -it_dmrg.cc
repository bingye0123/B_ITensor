
#include "symMonteCarlo.h"

#include "core.h"

//#include "model/tjmod.h"
//#include "hams/tJham.h"
#include "model/hubbardmod.h"
#include "hams/ExtendedHubbardham.h"

using boost::format;
using namespace std;

//typedef tJ
typedef Hubbard
states;

//number of electrons
const extern int Nferm=24;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm=500;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure=500;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure=40;


int
main(int argc, char* argv[])
{
  MPI::Init(argc, argv); // initialize MPI environment
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process

    int N=24;
    
	if (rank==0) { cout<<"Reading from disk:"<<endl; }
	states model(N);
	readFromFile("model_dat",model);
	IQMPS psi(model);
	psi.read("psi_folder");
	// readFromFile("psi_dat",psi);
	if (rank==0) { cout<<"Reading done."<<endl; }


	psi.position(1);
	if (rank==0) { cout << "Norm2: " << Dot(conj(psi.A(1)), psi.A(1))<<endl; }
	cout.flush();
/*
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
	
	if (rank==1) {
	cout<<"Nup: "<<nup_vec<<endl;
	cout<<"Ndn: "<<ndn_vec<<endl;
	cout<<"Total N: "<<nup_vec+ndn_vec<<endl;
	cout<<"Sz: "<<0.5*(nup_vec-ndn_vec)<<endl;
	}
	
	psi.position(1);
*/

	
	ITiVec inv_sym_mapping_table(N);
	//inversion: site -> N-1-site
	for(int i=0;i<inv_sym_mapping_table.size();i++)
	{
		inv_sym_mapping_table(i)=N-1-i;
	}
	
	ITiVec T1_sym_mapping_table_8_site="2 3 0 1 6 7 4 5";
	ITiVec T2_sym_mapping_table_8_site="4 5 6 7 0 1 2 3";
	ITiVec C6_sym_mapping_table_8_site="1 2 7 4 3 0 5 6";
	
	ITiVec T1_sym_mapping_table_32_site="3 4 5 6 7 8 1 2 11 12 13 14 15 16 9 10 19 20 21 22 23 24 17 18 27 28 29 30 31 32 25 26";
	T1_sym_mapping_table_32_site=T1_sym_mapping_table_32_site-1;
	ITiVec T2_sym_mapping_table_32_site="25 26 27 28 29 30 31 32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24";
	T2_sym_mapping_table_32_site=T2_sym_mapping_table_32_site-1;
	ITiVec C6_sym_mapping_table_32_site="4 5 14 15 24 17 26 27 2 3 12 13 22 23 32 25 8 1 10 11 20 21 30 31 6 7 16 9 18 19 28 29";
	C6_sym_mapping_table_32_site=C6_sym_mapping_table_32_site-1;
	
	ITiVec T1_sym_mapping_table_24_site="3 4 5 13 14 8 9 10 11 12 20 21 15 16 17 18 19 1 2 22 23 24 6 7";
	T1_sym_mapping_table_24_site=T1_sym_mapping_table_24_site-1;
	ITiVec T2_sym_mapping_table_24_site="21 22 23 24 6 1 2 3 4 5 13 14 7 8 9 10 11 12 20 15 16 17 18 19";
	T2_sym_mapping_table_24_site=T2_sym_mapping_table_24_site-1;
	ITiVec C6_sym_mapping_table_24_site="4 5 11 12 19 2 3 9 10 17 18 24 1 7 8 15 16 22 23 6 13 14 20 21";
	C6_sym_mapping_table_24_site=C6_sym_mapping_table_24_site-1;
	
//	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_8_site;
//	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_8_site;
//	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_8_site;
	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_32_site;
	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_32_site;
	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_32_site;
//	ITiVec T1_sym_mapping_table=T1_sym_mapping_table_24_site;
//	ITiVec T2_sym_mapping_table=T2_sym_mapping_table_24_site;
//	ITiVec C6_sym_mapping_table=C6_sym_mapping_table_24_site;
	
	ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table;
	translation_generator_mapping_table(1)=T2_sym_mapping_table;
	ITcVec COM_mom="1. 1.";
	
	ITiVecArray transl_C6_generator_mapping_table(3);
	transl_C6_generator_mapping_table(0)=T1_sym_mapping_table;
	transl_C6_generator_mapping_table(1)=T2_sym_mapping_table;
	transl_C6_generator_mapping_table(2)=C6_sym_mapping_table;
	ITcVec COM_Ang_mom="1. 1. 0.";
	COM_Ang_mom(2)=std::exp(Complex_i*Pi/3.);

	ITiVecArray C6_generator_mapping_table(1);
	C6_generator_mapping_table(0)=C6_sym_mapping_table;
	ITcVec Ang_mom(1);
	Ang_mom(0)=std::exp(-Complex_i*2.*Pi/3.);

	//std::string descsector("noGm2Pi3_10000");
	std::string descsectoreig("G_c6_500");
	//std::string descnorm("G_10000");
	//std::string desceig("inv_10000");
	//std::string desceig("C6_10000");
	
	//cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_measurement_tj_MonteCarlo(descsector,psi,C6_generator_mapping_table,Ang_mom)<<endl;

	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval:"<<symmetry_norm_measurement_tj_MonteCarlo(descnorm,psi,translation_generator_mapping_table, COM_mom)<<endl;

	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval:"<<symmetry_eigval_tj_MonteCarlo(desceig,psi,inv_sym_mapping_table)<<endl;
	//cout<<"From rank="<<rank<<" Monte Carlo sym eigval:"<<symmetry_eigval_tj_MonteCarlo(desceig,psi,C6_sym_mapping_table)<<endl;
	cout<<"From rank="<<rank<<" Monte Carlo measure: "<<symmetry_sector_eigval_hubbard_MonteCarlo(descsectoreig,psi, C6_sym_mapping_table, translation_generator_mapping_table, COM_mom)<<endl;


    MPI::Finalize();
    
    return 0;
}



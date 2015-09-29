#include "symMonteCarlo.h"

#define H24SAMPLE 0
#define DO_SZSZ 1
#define DO_CHIR 0
#define DO_PAIRS 0
//if 1 do all bonds on lattice, THIS NOW INCLUDES Sz(i)*Sz(i)
//if 0 cover maximal distance from fixed unit-cell (sites 0 and 1)
//if 2 all distances from fixed unit-cell (sites 0 and 1)
//#define ALL_BONDS 1
//if 1, gives 6 small and 2 big triangles per hexagon (128 triangs)
//else, gives triangles with all NN pairs ij and arbitrary k (1440 triangs)
#define CHIR_8TRIANGS_PERHEXAGON 1
//write to files the lattice states at each measurement step
#define WRITE_STATES 1
//norm <sympsi|psi>/<sympsi|sympsi>
#define DO_SECTORNORM 0

//LATTICE=0 indicates Honeycomb lattice
//LATTICE=1 indicates Triangular lattice
#define LATTICE 1
//SYSSHAPE=0 indicates system has rhombus shape
//SYSSHAPE=1 indicates system has hexagonal shape
#define SYSSHAPE 1

#define SYSSIZE 28

//For 36 sites
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the central row to edge rows
//For 16 sites
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the left col to right col
#define NUMBERING 0

//number of electrons
const extern int Nferm;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure;

//PBCs for the snake numbering of honeycomb
int pbca1r(int i, int j,int r){
	int y;
	if (div(j,r).quot>div(i,r).quot) { y=j-r; } else y=j;
	return y;
}

int pbca1l(int i,int r){
	int y;
	if ((div(i,r).quot<div(i+1,r).quot) || (i<0)) { y=i+r; } else y=i;
	return y;
}

int pbca2(int i, int Ns){
	int y;
	if (i>Ns) { y=i-Ns; } else if (i<0) { y=i+Ns; } else { y=i; }
	return y;
}


void make_bonds(ITiMatrix & bond_site_list, int Nsite){
    //SzSz measurement, for all bonds
	bond_site_list.set_size(Nsite*(Nsite+1)/2,2);
	int count=0;
	
	for(int i=0;i<Nsite;i++)
	{
		for(int j=i;j<Nsite;j++)
		{
			bond_site_list(count,0)=i;
			bond_site_list(count,1)=j;
			count++;
		}
	}
    /*
#if ALL_BONDS==2
	//SzSz measurement, bonds from first unit-cell
	bond_site_list.set_size(Nsite-1+Nsite-2,2);
	int count=0;
	for(int i=1;i<Nsite;i++)
	{
		bond_site_list(count,0)=0;
		bond_site_list(count,1)=i;
		count++;
		if (i!=1) {//site 1 does not bond to 0 nor to itself
			bond_site_list(count,0)=1;
			bond_site_list(count,1)=i;
			count++;
		}
	}
#endif
#if ALL_BONDS==1
	//SzSz measurement, for all bonds
	bond_site_list.set_size(Nsite*(Nsite+1)/2,2);
	int count=0;
	
	for(int i=0;i<Nsite;i++)
	{
		for(int j=i;j<Nsite;j++)
		{
			bond_site_list(count,0)=i;
			bond_site_list(count,1)=j;
			count++;
		}
	}
#elif ALL_BONDS==0
	//SzSz measurement, for up to two unit-cells in a1 and -a2 directions, starting from fixed unit-cell (sites 0 and 1)
	int d,L1;
	if (Nsite==32) { d=4/2;L1=8; } else if (Nsite==8) { d=2/2;L1=4; } else { cout << "Check bond choice for SzSz!!!"<<endl; d=0;L1=4; }
	bond_site_list.set_size(2*(d+1)*(d+1)-1+2*(d+1)*(d+1)-2,2);
	int count=0;
	
	for(int i=0;i<=2*d+1;i++)
	{//along a1 (horizontally)
		for(int j=0;j<=d;j++)
		{//along -a2
			if (!(i==0 && j==0)) {//site 0 does not bond to itself
				bond_site_list(count,0)=0;
				bond_site_list(count,1)=i+j*L1;
				count++;
			}
			if (!(i==0 && j==0) && !(i==1 && j==0)) {//site 1 does not bond to 0 nor to itself
				bond_site_list(count,0)=1;
				bond_site_list(count,1)=i+j*L1;
				count++;
			}
		}
	}
#endif
    */
}

#if (LATTICE==0 && SYSSHAPE==0 && SYSSIZE==8)
//for Si.(SjxSk) measurement; 8 triangles per hexagon, each clockwise;
//starting from middle site in small triangles; from upper-left in big-down triangles, bottom-left in big-up
void make_triangles_8perhexagon(ITiMatrix & tri_site_list,int Nsite){
	tri_site_list.set_size(4*Nsite,3);
	int rrow;//number of sites per row
	if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
	int count=0,j,k;
	for(int i=0;i<Nsite;i++)
	{
		if (!(i%2)) {//A site
			j=i+1;k=pbca2(i+rrow+1, Nsite);//triang 1
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=pbca2(i+rrow+1,Nsite);k=pbca1l(i-1,rrow);//triang 2
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=pbca1l(i-1,rrow);k=i+1;//triang 3
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=pbca1r(i,i+2,rrow);k=pbca1r(pbca2(i+rrow+1, Nsite),pbca2(i+rrow+2, Nsite),rrow);//big DOWN triangle
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
		} else {//B site
			j=pbca1r(i,i+1,rrow);k=i-1;//triang 4
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=i-1;k=pbca2(i-(rrow+1),Nsite);//triang 5
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=pbca2(i-(rrow+1),Nsite);k=pbca1r(i,i+1,rrow);//triang 6
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
			j=pbca2(i-rrow,Nsite);k=pbca1r(i,i+2,rrow);//big UP triangle
			tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
			count++;
		}
	}
}
#endif


#if (LATTICE==0 && SYSSHAPE==0 && SYSSIZE==8)
//for Si.(SjxSk) measurement; NN bonds ij and arbitrary k
void make_triangles_allfromNNbonds(ITiMatrix & tri_site_list, int Nsite){
	tri_site_list.set_size(3*Nsite*(Nsite-2)/2,3);
	int rrow;//number of sites per row
	if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
	int count=0,j;
	
	for(int i=0;i<Nsite;i++)
	{
		if (!(i%2)) {//A site has two NN
			j=i+1;//NN is B in same UC, no pbc
			for(int k=0;k<Nsite;k++)
			{
				if (k!=i && k!=j) {
					tri_site_list(count,0)=i;
					tri_site_list(count,1)=j;
					tri_site_list(count,2)=k;
					count++;
				}
			}
			j=(i+rrow+1)%Nsite;//NN is B in row below, pbc along a2
			for(int k=0;k<Nsite;k++)
			{
				if (k!=i && k!=j) {
					tri_site_list(count,0)=i;
					tri_site_list(count,1)=j;
					tri_site_list(count,2)=k;
					count++;
				}
			}
		} else {//B site has one NN in next UC along row
			if (div(i+1,rrow).quot>div(i,rrow).quot) {//pbc along a1
				j=i+1-rrow;
			} else {
				j=i+1;
			}
			for(int k=0;k<Nsite;k++)
			{
				if (k!=i && k!=j) {
					tri_site_list(count,0)=i;
					tri_site_list(count,1)=j;
					tri_site_list(count,2)=k;
					count++;
				}
			}
		}
	}
}
#endif

#if (LATTICE==0 && SYSSHAPE==0 && SYSSIZE==8)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(3*Nsite/2,2);
    ITiVec T1_sym_mapping_table_8_site="2 3 0 1 6 7 4 5";
	ITiVec T2_sym_mapping_table_8_site="4 5 6 7 0 1 2 3";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_8_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_8_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<4; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(5);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(1);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(2);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<3*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<3*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(3*(Nsite/2)*(3*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==16 && NUMBERING==0)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
	ITiVec T2_sym_mapping_table_16_site="12 13 14 15 0 1 2 3 4 5 6 7 8 9 10 11";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_16_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_16_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<16; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(4);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(15);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==16 && NUMBERING==1)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_16_site="4 5 6 7 8 9 10 11 12 13 14 15 0 1 2 3";
    ITiVec T2_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_16_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_16_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<16; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(7);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(12);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif


#if (LATTICE==1 && SYSSHAPE==1 && SYSSIZE==12)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_12_site="1 5 3 4 9 6 7 8 0 10 11 2";
	ITiVec T2_sym_mapping_table_12_site="10 11 0 1 5 2 3 4 9 6 7 8";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_12_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_12_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<12; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(2);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(9);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif


#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==36 && NUMBERING==0)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
	ITiVec T2_sym_mapping_table_36_site="30 31 32 33 34 35 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_36_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_36_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<36; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(6);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(35);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==36 && NUMBERING==1)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
	ITiVec T2_sym_mapping_table_36_site="12 13 14 15 16 17 0 1 2 3 4 5 24 25 26 27 28 29 6 7 8 9 10 11 30 31 32 33 34 35 18 19 20 21 22 23";
    
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_36_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_36_site;
	ITcVec COM_mom="1. 1.";
    COM_mom(0)=std::exp(Complex_i*2.0*Pi/3.);
    COM_mom(1)=std::exp(Complex_i*4.0*Pi/3.);
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<36; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(6);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(17);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif

#if (LATTICE==1 && SYSSHAPE==1 && SYSSIZE==28 && NUMBERING==0)
//NN bonds; unordered/ordered pairs of bonds, because <B\dag_kl B_ij> is just complex conjugate of <B\dag_ij B_kl>
//create shared list to identify pairs that share a site; order sites so that site j=k is shared; FAILS if both sites shared, so create bond_pairs carefully
void make_bond_pairs_FULL(ITiMatrix & bond_pairs_list, ITiVec & shared, int Nsite){
    
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_28_site="1 17 3 4 5 6 23 8 9 10 11 27 13 14 15 16 2 18 19 20 21 22 7 24 25 26 12 0";
	ITiVec T2_sym_mapping_table_28_site="24 25 11 27 0 1 17 2 3 4 5 6 22 7 8 9 10 26 12 13 14 15 16 18 19 20 21 23";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_28_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_28_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
    int count=0,j;
    for (int i=0; i<28; i++) {
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(1);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(3);
        count++;
        NNbonds(count,0)=translation_group_site_mapping_table(i)(0);
        NNbonds(count,1)=translation_group_site_mapping_table(i)(24);
        count++;
    }
    
    //cout<<"NNbonds: "<<endl<<NNbonds<<endl; cout.flush();
    
	//create (un)ordered pairs of NN bonds without pairing bond to itself
    //	bond_pairs_list.set_size((3*(Nsite/2)*(3*(Nsite/2)-1))/2,4);
	bond_pairs_list.set_size((6*(Nsite/2)*(6*(Nsite/2)-1)),4);
	count=0;
	for (int b=0; b<6*Nsite/2; b++) {
        //		for (int b2=b+1; b2<3*Nsite/2; b2++) {
        //		for (int b2=0; b2<b; b2++) {
		for (int b2=0; b2<6*Nsite/2; b2++) {
			if (b2!=b) {
				bond_pairs_list(count,0)=NNbonds(b,0);
				bond_pairs_list(count,1)=NNbonds(b,1);
				bond_pairs_list(count,2)=NNbonds(b2,0);
				bond_pairs_list(count,3)=NNbonds(b2,1);
				count++;
			}
		}
	}
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
    //	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2)*(6*(Nsite/2)-1));
	shared.zeros();
	int i,k,l;
	for (int t=0; t<bond_pairs_list.rows(); t++) {
		//check if site shared; if yes, make it the j=k site; FAILS if both sites shared, so create bond_pairs carefully
		if (bond_pairs_list(t,0)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,0)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,1); j=bond_pairs_list(t,0); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,2)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,2); l=bond_pairs_list(t,3); }
		if (bond_pairs_list(t,1)==bond_pairs_list(t,3)) { shared(t)=1; i=bond_pairs_list(t,0); j=bond_pairs_list(t,1); k=bond_pairs_list(t,3); l=bond_pairs_list(t,2); }
		if (shared(t)==1) {	bond_pairs_list(t,0)=i; bond_pairs_list(t,1)=j; bond_pairs_list(t,2)=k; bond_pairs_list(t,3)=l;  }
	}
    
}
#endif


#if (LATTICE==0 && SYSSHAPE==1 && SYSSIZE==24)
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(3*Nsite/2,4);
	int count=0;
	
	bps(count,0)=0;	bps(count,1)=6;	bps(count,2)=4;	bps(count,3)=10;	count++;
	bps(count,0)=2;	bps(count,1)=8;	bps(count,2)=13;bps(count,3)=19;	count++;
	bps(count,0)=4;	bps(count,1)=10;bps(count,2)=15;bps(count,3)=21;	count++;
	bps(count,0)=13;bps(count,1)=19;bps(count,2)=17;bps(count,3)=23;	count++;
	bps(count,0)=15;bps(count,1)=21;bps(count,2)=0;bps(count,3)=6;	count++;
	bps(count,0)=17;bps(count,1)=23;bps(count,2)=2;bps(count,3)=8;	count++;

	bps(count,0)=5;bps(count,1)=12;bps(count,2)=9;bps(count,3)=16;	count++;
	bps(count,0)=7;bps(count,1)=14;bps(count,2)=11;bps(count,3)=18;	count++;
	bps(count,0)=9;bps(count,1)=16;bps(count,2)=20;bps(count,3)=1;	count++;
	bps(count,0)=11;bps(count,1)=18;bps(count,2)=22;bps(count,3)=3;	count++;
	bps(count,0)=1;bps(count,1)=20;bps(count,2)=5;bps(count,3)=12;	count++;
	bps(count,0)=3;bps(count,1)=22;bps(count,2)=7;bps(count,3)=14;	count++;

	bps(count,0)=0;	bps(count,1)=6;bps(count,2)=3;bps(count,3)=4;	count++;
	bps(count,0)=2;	bps(count,1)=8;bps(count,2)=12;bps(count,3)=13;	count++;
	bps(count,0)=4;	bps(count,1)=10;bps(count,2)=14;bps(count,3)=15;	count++;
	bps(count,0)=13;bps(count,1)=19;bps(count,2)=16;bps(count,3)=17;	count++;
	bps(count,0)=15;bps(count,1)=21;bps(count,2)=18;bps(count,3)=0;	count++;
	bps(count,0)=17;bps(count,1)=23;bps(count,2)=1;bps(count,3)=2;	count++;
	
	bps(count,0)=5;bps(count,1)=12;bps(count,2)=8;bps(count,3)=9;	count++;
	bps(count,0)=7;bps(count,1)=14;bps(count,2)=10;bps(count,3)=11;	count++;
	bps(count,0)=9;bps(count,1)=16;bps(count,2)=19;bps(count,3)=20;	count++;
	bps(count,0)=11;bps(count,1)=18;bps(count,2)=21;bps(count,3)=22;	count++;
	bps(count,0)=1;bps(count,1)=20;bps(count,2)=23;bps(count,3)=5;	count++;
	bps(count,0)=3;bps(count,1)=22;bps(count,2)=6;bps(count,3)=7;	count++;
	
	bps(count,0)=0;	bps(count,1)=6;	bps(count,2)=9;bps(count,3)=10;	count++;
	bps(count,0)=2;	bps(count,1)=8;	bps(count,2)=11;bps(count,3)=19;	count++;
	bps(count,0)=4;	bps(count,1)=10;bps(count,2)=20;bps(count,3)=21;	count++;
	bps(count,0)=13;bps(count,1)=19;bps(count,2)=22;bps(count,3)=23;	count++;
	bps(count,0)=15;bps(count,1)=21;bps(count,2)=5;bps(count,3)=6;	count++;
	bps(count,0)=17;bps(count,1)=23;bps(count,2)=7;bps(count,3)=8;	count++;
	
	bps(count,0)=5;bps(count,1)=12;bps(count,2)=15;bps(count,3)=16;	count++;
	bps(count,0)=7;bps(count,1)=14;bps(count,2)=17;bps(count,3)=18;	count++;
	bps(count,0)=9;bps(count,1)=16;bps(count,2)=0;bps(count,3)=1;	count++;
	bps(count,0)=11;bps(count,1)=18;bps(count,2)=2;bps(count,3)=3;	count++;
	bps(count,0)=1;bps(count,1)=20;bps(count,2)=4;bps(count,3)=12;	count++;
	bps(count,0)=3;bps(count,1)=22;bps(count,2)=13;bps(count,3)=14;	count++;

	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
	//	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(3*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==0 && SYSSHAPE==0 && SYSSIZE==8)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(3*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiVec T1_sym_mapping_table_8_site="2 3 0 1 6 7 4 5";
	ITiVec T2_sym_mapping_table_8_site="4 5 6 7 0 1 2 3";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_8_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_8_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);

	int count=0;//,j,k,l;
    for (int i=0; i<4; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(7);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(5);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(7);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(3);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(7);
        count++;
    }
    
    
    
    
    /*
	//create all NN bonds
	for(int i=0;i<Nsite;i++)
	{
		if (!(i%2)) {//A site
			//NN is B in same UC, no pbc
			k=pbca1r(i,i+rrow/2,rrow);//furthest vertical bond
			l=pbca2(k+rrow+1,Nsite);
			j=i+1;
			bps(count,0)=i;
			bps(count,1)=j;
			bps(count,2)=k;
			bps(count,3)=l;
			count++;
			//NN is B in row below, pbc along a2
			j=(i+rrow+1)%Nsite;
			bps(count,0)=i;
			bps(count,1)=j;
			bps(count,2)=k;
			bps(count,3)=l;
			count++;
			//B site from previous UC along row
			j=pbca1l(i-1,rrow);
			bps(count,0)=i;
			bps(count,1)=j;
			bps(count,2)=k;
			bps(count,3)=l;
			count++;
		}
	}
    */
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
	//	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(3*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==16 && NUMBERING==0)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
	ITiVec T2_sym_mapping_table_16_site="12 13 14 15 0 1 2 3 4 5 6 7 8 9 10 11";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_16_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_16_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<16; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(8);
        bps(count,1)=translation_group_site_mapping_table(i)(12);
        bps(count,2)=translation_group_site_mapping_table(i)(5);
        bps(count,3)=translation_group_site_mapping_table(i)(6); //////10,11
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(12);
        bps(count,1)=translation_group_site_mapping_table(i)(13);
        bps(count,2)=translation_group_site_mapping_table(i)(5);
        bps(count,3)=translation_group_site_mapping_table(i)(6);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(13);
        bps(count,1)=translation_group_site_mapping_table(i)(8);
        bps(count,2)=translation_group_site_mapping_table(i)(5);
        bps(count,3)=translation_group_site_mapping_table(i)(6);
        count++;
    }
    
    
    
    
    /*
     //create all NN bonds
     for(int i=0;i<Nsite;i++)
     {
     if (!(i%2)) {//A site
     //NN is B in same UC, no pbc
     k=pbca1r(i,i+rrow/2,rrow);//furthest vertical bond
     l=pbca2(k+rrow+1,Nsite);
     j=i+1;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //NN is B in row below, pbc along a2
     j=(i+rrow+1)%Nsite;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //B site from previous UC along row
     j=pbca1l(i-1,rrow);
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     }
     }
     */
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
	//	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==16 && NUMBERING==1)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_16_site="4 5 6 7 8 9 10 11 12 13 14 15 0 1 2 3";
    ITiVec T2_sym_mapping_table_16_site="1 2 3 0 5 6 7 4 9 10 11 8 13 14 15 12";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_16_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_16_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<16; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(1);
        bps(count,1)=translation_group_site_mapping_table(i)(0);
        bps(count,2)=translation_group_site_mapping_table(i)(6);
        bps(count,3)=translation_group_site_mapping_table(i)(10); //////10,11
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(4);
        bps(count,2)=translation_group_site_mapping_table(i)(6);
        bps(count,3)=translation_group_site_mapping_table(i)(10);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(4);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(6);
        bps(count,3)=translation_group_site_mapping_table(i)(10);
        count++;
    }
    
    
    
    
    /*
     //create all NN bonds
     for(int i=0;i<Nsite;i++)
     {
     if (!(i%2)) {//A site
     //NN is B in same UC, no pbc
     k=pbca1r(i,i+rrow/2,rrow);//furthest vertical bond
     l=pbca2(k+rrow+1,Nsite);
     j=i+1;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //NN is B in row below, pbc along a2
     j=(i+rrow+1)%Nsite;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //B site from previous UC along row
     j=pbca1l(i-1,rrow);
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     }
     }
     */
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
	//	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==1 && SYSSIZE==12)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_12_site="1 5 3 4 9 6 7 8 0 10 11 2";
	ITiVec T2_sym_mapping_table_12_site="10 11 0 1 5 2 3 4 9 6 7 8";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_12_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_12_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<12; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(7);
        bps(count,3)=translation_group_site_mapping_table(i)(8);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(2);
        bps(count,2)=translation_group_site_mapping_table(i)(7);
        bps(count,3)=translation_group_site_mapping_table(i)(8);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(9);
        bps(count,2)=translation_group_site_mapping_table(i)(7);
        bps(count,3)=translation_group_site_mapping_table(i)(8);
        count++;
    }
    
    
    
    
    /*
     //create all NN bonds
     for(int i=0;i<Nsite;i++)
     {
     if (!(i%2)) {//A site
     //NN is B in same UC, no pbc
     k=pbca1r(i,i+rrow/2,rrow);//furthest vertical bond
     l=pbca2(k+rrow+1,Nsite);
     j=i+1;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //NN is B in row below, pbc along a2
     j=(i+rrow+1)%Nsite;
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     //B site from previous UC along row
     j=pbca1l(i-1,rrow);
     bps(count,0)=i;
     bps(count,1)=j;
     bps(count,2)=k;
     bps(count,3)=l;
     count++;
     }
     }
     */
	//create list identifying pairs that share a site; in that case, order sites in bond_pairs so that j=k
	//	shared.set_size(3*Nsite/2*(3*Nsite/2-1)/2);
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif


#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==4)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_4_site="1 0 3 2";
	ITiVec T2_sym_mapping_table_4_site="2 3 0 1";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_4_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_4_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<4; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(3);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(2);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(3);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(3);
        bps(count,2)=translation_group_site_mapping_table(i)(2);
        bps(count,3)=translation_group_site_mapping_table(i)(3);
        count++;
    }

	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==36 && NUMBERING==0)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
	ITiVec T2_sym_mapping_table_36_site="30 31 32 33 34 35 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_36_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_36_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<4; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(21);
        bps(count,3)=translation_group_site_mapping_table(i)(22);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(6);
        bps(count,2)=translation_group_site_mapping_table(i)(21);
        bps(count,3)=translation_group_site_mapping_table(i)(22);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(35);
        bps(count,2)=translation_group_site_mapping_table(i)(21);
        bps(count,3)=translation_group_site_mapping_table(i)(22);
        count++;
    }
    
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==0 && SYSSIZE==36 && NUMBERING==1)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_36_site="1 2 3 4 5 0 7 8 9 10 11 6 13 14 15 16 17 12 19 20 21 22 23 18 25 26 27 28 29 24 31 32 33 34 35 30";
	ITiVec T2_sym_mapping_table_36_site="12 13 14 15 16 17 0 1 2 3 4 5 24 25 26 27 28 29 6 7 8 9 10 11 30 31 32 33 34 35 18 19 20 21 22 23";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_36_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_36_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<4; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(1);
        bps(count,2)=translation_group_site_mapping_table(i)(33);
        bps(count,3)=translation_group_site_mapping_table(i)(34);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(6);
        bps(count,2)=translation_group_site_mapping_table(i)(33);
        bps(count,3)=translation_group_site_mapping_table(i)(34);
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(0);
        bps(count,1)=translation_group_site_mapping_table(i)(17);
        bps(count,2)=translation_group_site_mapping_table(i)(33);
        bps(count,3)=translation_group_site_mapping_table(i)(34);
        count++;
    }
    
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif

#if (LATTICE==1 && SYSSHAPE==1 && SYSSIZE==28 && NUMBERING==0)
//NN bonds; three inequivalent bonds furthest from fixed bond, one per unit-cell
//create shared list, although no pairs share a site
void make_bond_pairs(ITiMatrix & bps, ITiVec & shared, int Nsite){
	bps.set_size(6*(Nsite/2),4);
	//int rrow;//number of sites per row
	//if (Nsite==32) { rrow=8; } else if (Nsite==8) { rrow=4;	}
    ITiMatrix NNbonds(6*Nsite/2,2);
    ITiVec T1_sym_mapping_table_28_site="1 17 3 4 5 6 23 8 9 10 11 27 13 14 15 16 2 18 19 20 21 22 7 24 25 26 12 0";
	ITiVec T2_sym_mapping_table_28_site="24 25 11 27 0 1 17 2 3 4 5 6 22 7 8 9 10 26 12 13 14 15 16 18 19 20 21 23";
    ITiVecArray translation_generator_mapping_table(2);
	translation_generator_mapping_table(0)=T1_sym_mapping_table_28_site;
	translation_generator_mapping_table(1)=T2_sym_mapping_table_28_site;
	ITcVec COM_mom="1. 1.";
    
    ITiVecArray translation_group_site_mapping_table;
    ITcVec translation_group_eigval_table;
    generate_symmetry_group(translation_generator_mapping_table, COM_mom, translation_group_site_mapping_table, translation_group_eigval_table);
    
	int count=0;//,j,k,l;
    for (int i=0; i<4; i++) {
        bps(count,0)=translation_group_site_mapping_table(i)(14);
        bps(count,1)=translation_group_site_mapping_table(i)(15);
        bps(count,2)=translation_group_site_mapping_table(i)(6);
        bps(count,3)=translation_group_site_mapping_table(i)(23); ///the bond -
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(14);
        bps(count,1)=translation_group_site_mapping_table(i)(15);
        bps(count,2)=translation_group_site_mapping_table(i)(18);
        bps(count,3)=translation_group_site_mapping_table(i)(23); ///the bond \
        count++;
        bps(count,0)=translation_group_site_mapping_table(i)(14);
        bps(count,1)=translation_group_site_mapping_table(i)(15);
        bps(count,2)=translation_group_site_mapping_table(i)(6);
        bps(count,3)=translation_group_site_mapping_table(i)(18); ///the bond /
        count++;
    }
    
	shared.set_size(6*(Nsite/2));
	shared.zeros();
}
#endif


ITcVecArray symmetry_sector_measurement_tj_MonteCarlo(std::string desc,const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table)
{
  ITcVecArray result(1);
  
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
    int Nsite=psi.N();
    
    
  ITiVecArray sym_sector_group_site_mapping_table;
  ITcVec sym_sector_group_eigval_table;
  generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
/////////////////////////////////////
  state_storage_tj cur_state=initialize_random_state_vec_tj_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);

//either do all bonds (fine for 2x2x2) or
//fix single unit-cell (OK if projected to CoM momentum Gamma) and consider its bonds up to largest possible length
	//ITiMatrix bond_site_list;
#if DO_SZSZ==1
	ITiMatrix bond_site_list;
	make_bonds(bond_site_list, Nsite);
/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_site_list for S.S: "<<endl<<bond_site_list<<endl;
		cout.flush();
	}
/////////////////////////////////////
	ITdMatrix SzSz_list(Monte_Carlo_nmeasure,bond_site_list.rows());
#endif
#if DO_CHIR==1
#if CHIR_8TRIANGS_PERHEXAGON==1
	ITiMatrix tri_site_list;
	//for Si.(SjxSk) measurement; 8 triangles per hexagon, each clockwise;
	//starting from middle site in small triangles; from upper-left in big-down triangles, bottom-left in big-up
	make_triangles_8perhexagon(tri_site_list,Nsite);
#else
	//for Si.(SjxSk) measurement; NN bonds ij and arbitrary k
	ITiMatrix tri_site_list;
	make_triangles_allfromNNbonds(tri_site_list,Nsite);
#endif
	/////////////////////////////////////
	if (rank==0) {
		cout<<"tri_site_list for S.(SxS): "<<endl<<tri_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITcMatrix Chir_list(Monte_Carlo_nmeasure,tri_site_list.rows());
#endif
#if DO_PAIRS==1
	ITiMatrix bond_pairs_list;
	ITiVec shared_list;
	make_bond_pairs_FULL(bond_pairs_list, shared_list,Nsite);
	/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_pairs_list for pair-pair correlation: "<<endl<<bond_pairs_list<<endl; cout.flush();
	}
	/////////////////////////////////////
	ITcMatrix pairpair_list(Monte_Carlo_nmeasure,bond_pairs_list.rows());
#endif

	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			if(i%10==0)cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////

		try_one_step_tj_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	}
#if WRITE_STATES==1
/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "states"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str ();
	std::ofstream stout;
	stout.open (ime1.c_str ());
/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				if(j%10==0)cout<<"j:"<<j<<endl;
				cout.flush();
			}
			/////////////////////////////////////
			try_one_step_tj_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
		}
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
#if DO_SZSZ==1
			ITdVec cur_SzSz=SzSzcorr_tj( psi, state_vec, bond_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			//dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			SzSz_list.set_row(measure_count,cur_SzSz);
#endif
#if DO_CHIR==1
			ITcVec cur_Chir=SSS_tj( psi, state_vec, tri_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			Chir_list.set_row(measure_count,cur_Chir);
/////////////////////////////////////
//			if (rank==0) {
//				cout<<endl<<state_vec<<endl;
//				cout<<endl<<cur_Chir<<endl;
//			}
/////////////////////////////////////
#endif
#if DO_PAIRS==1
			ITcVec cur_pairpair=pairpair_tj( psi, state_vec, bond_pairs_list, shared_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			pairpair_list.set_row(measure_count,cur_pairpair);
#endif
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif
#if DO_SZSZ==1
	ITdVec SzSz_result(bond_site_list.rows()),SzSz_dev(bond_site_list.rows());
	for(int i=0;i<SzSz_result.size();i++)
	{
		SzSz_result(i)=itpp::mean(SzSz_list.get_col(i));
		SzSz_dev(i)=std::sqrt(itpp::variance(SzSz_list.get_col(i))/SzSz_list.rows());
	}
	cout<<"From rank="<<rank<<": SzSz_mean="<<SzSz_result<<endl<<" stderr="<<SzSz_dev<<endl;

/////////////////////////////////////
	std::stringstream fim;
	fim << "SzSz"<<desc<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime = fim.str ();
	std::ofstream finout;
	finout.open (ime.c_str ());
	finout<<bond_site_list.rows()<<endl;
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<bond_site_list(pair,0)<<endl;
		finout<<bond_site_list(pair,1)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_dev(pair)<<endl;
	}
	finout.close();
/////////////////////////////////////
	result(0)=itpp::to_cvec(SzSz_result);
#endif
#if DO_CHIR==1
	ITcVec Chir_result(tri_site_list.rows()),Chir_dev(tri_site_list.rows());
	for(int i=0;i<Chir_result.size();i++)
	{
		Chir_result(i)=itpp::mean(Chir_list.get_col(i));
		Chir_dev(i)=std::sqrt(itpp::variance(Chir_list.get_col(i))/Chir_list.rows());
	}
	cout<<"From rank="<<rank<<": S.(SxS)_mean="<<Chir_result<<endl<<" stderr="<<Chir_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream cfim;
	cfim << "Chir"<<desc<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string cime = cfim.str();
	std::ofstream cfinout;
	cfinout.open(cime.c_str());
	cfinout<<tri_site_list.rows()<<endl;
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<tri_site_list(pair,0)<<endl;
		cfinout<<tri_site_list(pair,1)<<endl;
		cfinout<<tri_site_list(pair,2)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_result(pair)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_dev(pair)<<endl;
	}
	cfinout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(Chir_result);
#endif
#if DO_PAIRS==1
	ITcVec pairpair_result(bond_pairs_list.rows()),pairpair_dev(bond_pairs_list.rows());
	for(int i=0;i<pairpair_result.size();i++)
	{
		pairpair_result(i)=itpp::mean(pairpair_list.get_col(i));
		pairpair_dev(i)=std::sqrt(itpp::variance(pairpair_list.get_col(i))/pairpair_list.rows());
	}
	cout<<"From rank="<<rank<<": pair-pair_mean="<<pairpair_result<<endl<<" stderr="<<pairpair_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream pcfim;
	pcfim << "PP"<<desc<<"rank"<<rank<<".txt";
	std::string pcime = pcfim.str();
	std::ofstream pcfinout;
	pcfinout.open(pcime.c_str());
	pcfinout<<bond_pairs_list.rows()<<endl;
	for (int pair=0; pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<bond_pairs_list(pair,0)<<endl;
		pcfinout<<bond_pairs_list(pair,1)<<endl;
		pcfinout<<bond_pairs_list(pair,2)<<endl;
		pcfinout<<bond_pairs_list(pair,3)<<endl;
	}
	for (int pair=0;pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<pairpair_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<pairpair_dev(pair)<<endl;
	}
	pcfinout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(pairpair_result);
#endif

	return result;
}


ITcVecArray symmetry_sector_measurement_hubbard_MonteCarlo(std::string desc,const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table)
{
	ITcVecArray result(1);
	
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
	//cout<<size<<endl<<rank<<endl;
    int Nsite=psi.N();
    
    
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	state_storage_hubbard cur_state=initialize_random_state_vec_hubbard_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);

//either do all bonds (fine for 2x2x2) or
//fix single unit-cell (OK if projected to CoM momentum Gamma) and consider its bonds up to largest possible length
	//ITiMatrix bond_site_list;
#if DO_SZSZ==1
	ITiMatrix bond_site_list;
	make_bonds(bond_site_list, Nsite);
/*	if (rank==0) {
		std::stringstream fim;
		fim << "bonds.txt";
		std::string ime = fim.str ();
		std::ofstream finout;
		finout.open (ime.c_str ());
		finout<<endl<<bond_site_list<<endl;
		finout.close();
	}
*/
	/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_site_list for S.S: "<<endl<<bond_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITdMatrix SzSz_list(Monte_Carlo_nmeasure,bond_site_list.rows());
#endif
	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			if(i%10==0)cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		try_one_step_hubbard_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	}
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "states"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str();
	//cout<<ime1<<endl;
	std::ofstream stout;
	stout.open(ime1.c_str());
	/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				if(j%10==0)cout<<"j:"<<j<<endl;
				cout.flush();
			}
			/////////////////////////////////////
			try_one_step_hubbard_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
		}
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
#if DO_SZSZ==1
			ITdVec cur_SzSz=SzSzcorr_hubbard( psi, state_vec, bond_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			//dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			SzSz_list.set_row(measure_count,cur_SzSz);
#endif
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif
#if DO_SZSZ==1
	ITdVec SzSz_result(bond_site_list.rows()),SzSz_dev(bond_site_list.rows());
	for(int i=0;i<SzSz_result.size();i++)
	{
		SzSz_result(i)=itpp::mean(SzSz_list.get_col(i));
		SzSz_dev(i)=std::sqrt(itpp::variance(SzSz_list.get_col(i))/SzSz_list.rows());
	}
	cout<<"From rank="<<rank<<": SzSz_mean="<<SzSz_result<<endl<<" stderr="<<SzSz_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream fim;
	fim << "SzSz"<<desc<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	finout<<bond_site_list.rows()<<endl;
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<bond_site_list(pair,0)<<endl;
		finout<<bond_site_list(pair,1)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_dev(pair)<<endl;
	}
	finout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(SzSz_result);
#endif
	return result;
}



ITdVec SzSzcorr_tj(const IQMPS & psi, const ITiVec & state_vec, const ITiMatrix & bond_site_list,const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  
  ITdVec Sz_vec="0. 1. -1.";
  Sz_vec*=0.5;
  ITdVec result(bond_site_list.rows());
  for(int i=0;i<bond_site_list.rows();i++)
  { 
    result(i)=Sz_vec(state_vec(bond_site_list(i,0))-1)*Sz_vec(state_vec(bond_site_list(i,1))-1);
  }
//     cout<<"cur state vec:"<<state_vec<<endl;
//     cout<<"S.S:"<<result<<endl;
  return result;
}


ITdVec SzSzcorr_hubbard(const IQMPS & psi, const ITiVec & state_vec, const ITiMatrix & bond_site_list,const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
	
	ITdVec Sz_vec="0. 1. -1. 0.";
	Sz_vec*=0.5;
	ITdVec result(bond_site_list.rows());
	for(int i=0;i<bond_site_list.rows();i++)
	{
		result(i)=Sz_vec(state_vec(bond_site_list(i,0))-1)*Sz_vec(state_vec(bond_site_list(i,1))-1);
	}
	//     cout<<"cur state vec:"<<state_vec<<endl;
	//     cout<<"S.S:"<<result<<endl;
	return result;
}


ITcVec SSS_tj(const IQMPS & psi, const ITiVec & st, const ITiMatrix & tri, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	ITcVec result(tri.rows());
	result.zeros();
	dComplex old_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	//dComplex old_amp=quantum_amplitude(psi,st);
	/////////////////////////////////////
	int rank=MPI::COMM_WORLD.Get_rank();
	//	if (rank==0) { cout<<old_amp<<endl<<tri.rows()<<endl; }
	/////////////////////////////////////
	for (int t=0; t<tri.rows(); t++) {
		if (st[tri(t,0)]!=1 && st[tri(t,1)]!=1 && st[tri(t,2)]!=1) {
			if ((st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=6 && (st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=9) {
				ITiVec nst1(st);
				nst1[tri(t,0)]=st[tri(t,1)];
				nst1[tri(t,1)]=st[tri(t,2)];
				nst1[tri(t,2)]=st[tri(t,0)];
				dComplex new_amp1=symmetrized_quantum_amplitude_tj(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
				//dComplex new_amp1=quantum_amplitude(psi,nst1);
				result[t]+=-0.25*Complex_i*new_amp1/old_amp;
				
				ITiVec nst2(st);
				nst2[tri(t,0)]=st[tri(t,2)];
				nst2[tri(t,1)]=st[tri(t,0)];
				nst2[tri(t,2)]=st[tri(t,1)];
				dComplex new_amp2=symmetrized_quantum_amplitude_tj(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
				//dComplex new_amp2=quantum_amplitude(psi,nst2);
				result[t]+=0.25*Complex_i*new_amp2/old_amp;
/////////////////////////////////////
//				if (rank==0) { cout<<"T"<<tri.get_row(t)<<"  :  "<<st<<" "<<old_amp<<" ;  "<<nst1<<" "<<sgn1*new_amp1<<" ;  "<<nst2<<" "<<sgn2*new_amp2<<" ,"<<result[t]<<endl; cout.flush();}
/////////////////////////////////////
			}
		}
		/////////////////////////////////////
		if (rank==0) { cout<<t<<","; cout.flush();}
		/////////////////////////////////////
	}
	/////////////////////////////////////
	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return result;
}


ITcVec SSS_hubbard(const IQMPS & psi, const ITiVec & st, const ITiMatrix & tri, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	ITcVec result(tri.rows());
	result.zeros();
	dComplex old_amp=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	/////////////////////////////////////
	int rank=MPI::COMM_WORLD.Get_rank();
//	if (rank==0) { cout<<old_amp<<endl<<tri.rows()<<endl; }
	/////////////////////////////////////
	for (int t=0; t<tri.rows(); t++) {
		if (st[tri(t,0)]!=1 && st[tri(t,1)]!=1 && st[tri(t,2)]!=1 && st[tri(t,0)]!=4 && st[tri(t,1)]!=4 && st[tri(t,2)]!=4) {
			if ((st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=6 && (st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=9) {
				ITiVec nst1(st);
				nst1[tri(t,0)]=st[tri(t,1)];
				nst1[tri(t,1)]=st[tri(t,2)];
				nst1[tri(t,2)]=st[tri(t,0)];
				dComplex new_amp1=symmetrized_quantum_amplitude_hubbard(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
				result[t]+=-0.25*Complex_i*std::conj(new_amp1)/old_amp;
				
				ITiVec nst2(st);
				nst2[tri(t,0)]=st[tri(t,2)];
				nst2[tri(t,1)]=st[tri(t,0)];
				nst2[tri(t,2)]=st[tri(t,1)];
				dComplex new_amp2=symmetrized_quantum_amplitude_hubbard(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
				result[t]+=0.25*Complex_i*std::conj(new_amp2)/old_amp;
			}
		}
		/////////////////////////////////////
		if (rank==0) { cout<<t<<","; cout.flush();}
		/////////////////////////////////////
	}
	/////////////////////////////////////
	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return result;
}


dComplex sectornorm_tj(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	dComplex sym_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex amp=quantum_amplitude(psi,st);
	/////////////////////////////////////
//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return amp/sym_amp;
}


dComplex sectornorm_hubbard(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	dComplex sym_amp=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex amp=quantum_amplitude(psi,st);
	/////////////////////////////////////
	//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return amp/sym_amp;
}


int hopsign(const ITiVec & state, int site, int neigh){
	int x=0,sign;
	if (site<neigh) { for(int m=site+1;m<neigh;m++){ if(state(m)>1 && state(m)<4) { x++; } }
	} else { for(int m=neigh+1;m<site;m++){ if(state(m)>1 && state(m)<4) { x++; } } }
	if(x%2) { sign=-1; } else { sign=1; }
	return sign;
}


ITcVec pairpair_tj(const IQMPS & psi, const ITiVec & st, const ITiMatrix & bps,  const ITiVec & shared, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	ITcVec result(bps.rows());
	result.zeros();
	dComplex old_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	/////////////////////////////////////
	int rank=MPI::COMM_WORLD.Get_rank();
		//if (rank==0) { cout<<old_amp<<endl<<bps.rows()<<endl; }
	/////////////////////////////////////
	int i,j,k,l,sign,temp;
	for (int t=0; t<bps.rows(); t++) {
		if (shared(t)==1) {//bonds share site j=k
			if (st[bps(t,0)]==1) {//non-zero matrix element only if site-i is hole
				if (st[bps(t,1)]*st[bps(t,3)]==6) {//j=k is up, l is down, OR vice versa
					//hop (i,l) times n_j
					i=bps(t,0); j=bps(t,1); k=bps(t,2); l=bps(t,3);
					ITiVec nst1(st);
					nst1[i]=st[l];
					nst1[l]=st[i];
					sign=hopsign(st,i,l);
					dComplex new_amp1=symmetrized_quantum_amplitude_tj(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
					result[t]+=(1.*sign)*new_amp1/old_amp;

					//hop (i,l) times flip j
					ITiVec nst2(st);
					nst2[i]=st[j];
					nst2[j]=st[l];
					nst2[l]=st[i];//this is just the hole
					dComplex new_amp2=symmetrized_quantum_amplitude_tj(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
					//the fermion sign here is always (for any order of i,j,l) opposite to hopping sign between (i,l) calculated above
					result[t]+=-(1.*sign)*new_amp2/old_amp;
				}
			}
		} else {//bonds do not overlap
			if (st[bps(t,0)]==1 && st[bps(t,1)]==1) {//non-zero only if both sites i,j are holes
				if (st[bps(t,2)]*st[bps(t,3)]==6) {//k is up, l is down, OR vice versa
					i=bps(t,0); j=bps(t,1); k=bps(t,2); l=bps(t,3);
					ITiVec nst1(st);
					nst1[i]=st[l];
					nst1[l]=st[i];//this is just a hole
					temp=hopsign(st,i,l);
					sign=hopsign(nst1,j,k);
					nst1[j]=nst1[k];
					nst1[k]=1;//this is just the hole nst1[j]=st[j]
					dComplex new_amp1=symmetrized_quantum_amplitude_tj(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
					result[t]+=(1.*sign*temp)*new_amp1/old_amp;
					//if (rank==0) { cout<<st<<endl<<bps.get_row(t)<<st[i]<<st[j]<<st[k]<<st[l]<<sign*temp<<endl; }
					
					ITiVec nst2(st);
					nst2[i]=st[k];
					nst2[k]=st[i];//this is just a hole
					temp=hopsign(st,i,k);
					sign=hopsign(nst2,j,l);
					nst2[j]=nst2[l];
					nst2[l]=1;//this is just the hole nst2[j]=st[j]
					dComplex new_amp2=symmetrized_quantum_amplitude_tj(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
					result[t]+=(1.*sign*temp)*new_amp2/old_amp;
					//if (rank==0) { cout<<sign*temp<<endl; }
				}
			}
		}
		/////////////////////////////////////
		if (rank==0) { cout<<t<<","; cout.flush();}
		/////////////////////////////////////
	}
	/////////////////////////////////////
	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return result;
}




ITcVecArray symmetry_sector_measurement_fromfiles_tj_MonteCarlo(std::ifstream& inputstates, std::string desc,int nmeasure,const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table){
	ITcVecArray result(1);
	
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
	//cout<<size<<endl<<rank<<endl;
	int Nsite=psi.N();
	
	//Load states from file
	if (!inputstates) { cout<<"!No input file on "<<rank<<endl; exit(1); }
	std::vector<int> loadstates;
	loadstates.reserve(10000*Nsite);
	std::string str;
	int tempi,co=1;
	while (!std::getline(inputstates,str).eof()) {
		std::istringstream ss(str);
		if(!(ss>>tempi)) { cout<<"Bad input!"; }
		else loadstates.push_back(tempi);
	}
	cout<<"End of loading states"<<endl;
	inputstates.close();
	int nmeas,ndata=loadstates.size();
	cout<<ndata<<endl;
	if (ndata%Nsite) { cout<<"Partial data!!!"<<endl; }
	nmeas=ndata/Nsite;
	if (nmeasure>0) { nmeas=nmeasure; }

	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
#if DO_SZSZ==1
	ITiMatrix bond_site_list;
	make_bonds(bond_site_list, Nsite);
	/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_site_list for S.S: "<<endl<<bond_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITdMatrix SzSz_list(nmeas,bond_site_list.rows());
#endif
#if DO_CHIR==1
#if CHIR_8TRIANGS_PERHEXAGON==1
	ITiMatrix tri_site_list;
	//for Si.(SjxSk) measurement; 8 triangles per hexagon, each clockwise;
	//starting from middle site in small triangles; from upper-left in big-down triangles, bottom-left in big-up
	make_triangles_8perhexagon(tri_site_list,Nsite);
#else
	//for Si.(SjxSk) measurement; NN bonds ij and arbitrary k
	ITiMatrix tri_site_list;
	make_triangles_allfromNNbonds(tri_site_list,Nsite);
#endif
	/////////////////////////////////////
	if (rank==0) {
		cout<<"tri_site_list for S.(SxS): "<<endl<<tri_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITcMatrix Chir_list(nmeas,tri_site_list.rows());
#endif
#if DO_SECTORNORM==1
	ITcVec SecNorm_list(nmeas);
#endif
#if DO_PAIRS==1
	ITiMatrix bond_pairs_list;
	ITiVec shared_list;
	make_bond_pairs(bond_pairs_list, shared_list,Nsite);
	/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_pairs_list for pair-pair correlation: "<<endl<<bond_pairs_list<<endl; cout.flush();
	}
	/////////////////////////////////////
	ITcMatrix pairpair_list(nmeas,bond_pairs_list.rows());
#endif
	
	int measure_count=0;
	
	for(int i=0;i<nmeas;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec(Nsite);
			for (int ss=0; ss<Nsite; ss++) {
				state_vec(ss)=loadstates[Nsite*measure_count+ss];
			}
#if DO_SZSZ==1
			ITdVec cur_SzSz=SzSzcorr_tj( psi, state_vec, bond_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			SzSz_list.set_row(measure_count,cur_SzSz);
#endif
#if DO_CHIR==1
			ITcVec cur_Chir=SSS_tj( psi, state_vec, tri_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			Chir_list.set_row(measure_count,cur_Chir);
#endif
#if DO_SECTORNORM==1
			SecNorm_list(measure_count)=sectornorm_tj( psi, state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
#endif
#if DO_PAIRS==1
			ITcVec cur_pairpair=pairpair_tj( psi, state_vec, bond_pairs_list, shared_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			pairpair_list.set_row(measure_count,cur_pairpair);
#endif
			//dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			measure_count++;
		}
	}
	
#if DO_SZSZ==1
	ITdVec SzSz_result(bond_site_list.rows()),SzSz_dev(bond_site_list.rows());
	for(int i=0;i<SzSz_result.size();i++)
	{
		SzSz_result(i)=itpp::mean(SzSz_list.get_col(i));
		SzSz_dev(i)=std::sqrt(itpp::variance(SzSz_list.get_col(i))/SzSz_list.rows());
	}
	cout<<"From rank="<<rank<<": SzSz_mean="<<SzSz_result<<endl<<" stderr="<<SzSz_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream fim;
	fim << "SzSzload"<<desc<<"rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	finout<<bond_site_list.rows()<<endl;
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<bond_site_list(pair,0)<<endl;
		finout<<bond_site_list(pair,1)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_dev(pair)<<endl;
	}
	finout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(SzSz_result);
#endif
#if DO_CHIR==1
	ITcVec Chir_result(tri_site_list.rows()),Chir_dev(tri_site_list.rows());
	for(int i=0;i<Chir_result.size();i++)
	{
		Chir_result(i)=itpp::mean(Chir_list.get_col(i));
		Chir_dev(i)=std::sqrt(itpp::variance(Chir_list.get_col(i))/Chir_list.rows());
	}
	cout<<"From rank="<<rank<<": S.(SxS)_mean="<<Chir_result<<endl<<" stderr="<<Chir_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream cfim;
	cfim << "Chirload"<<desc<<"rank"<<rank<<".txt";
	std::string cime = cfim.str();
	std::ofstream cfinout;
	cfinout.open(cime.c_str());
	cfinout<<tri_site_list.rows()<<endl;
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<tri_site_list(pair,0)<<endl;
		cfinout<<tri_site_list(pair,1)<<endl;
		cfinout<<tri_site_list(pair,2)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_result(pair)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_dev(pair)<<endl;
	}
	cfinout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(Chir_result);
#endif
#if DO_SECTORNORM==1
	dComplex SecNorm_result=itpp::mean(SecNorm_list);
	dComplex SecNorm_dev=std::sqrt(itpp::variance(SecNorm_list));
	
	cout<<"From rank="<<rank<<": SecNorm_mean="<<SecNorm_result<<endl<<" stderr="<<SecNorm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "SecNormload"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<SecNorm_result<<endl<<SecNorm_dev<<endl;
	for (int pair=0; pair<SecNorm_list.size(); pair++) {
		nfinout<<SecNorm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(SecNorm_result);
	result(0)=temp;
#endif
#if DO_PAIRS==1
	ITcVec pairpair_result(bond_pairs_list.rows()),pairpair_dev(bond_pairs_list.rows());
	for(int i=0;i<pairpair_result.size();i++)
	{
		pairpair_result(i)=itpp::mean(pairpair_list.get_col(i));
		pairpair_dev(i)=std::sqrt(itpp::variance(pairpair_list.get_col(i))/pairpair_list.rows());
	}
	cout<<"From rank="<<rank<<": pair-pair_mean="<<pairpair_result<<endl<<" stderr="<<pairpair_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream pcfim;
	pcfim << "PPload"<<desc<<"rank"<<rank<<".txt";
	std::string pcime = pcfim.str();
	std::ofstream pcfinout;
	pcfinout.open(pcime.c_str());
	pcfinout<<bond_pairs_list.rows()<<endl;
	for (int pair=0; pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<bond_pairs_list(pair,0)<<endl;
		pcfinout<<bond_pairs_list(pair,1)<<endl;
		pcfinout<<bond_pairs_list(pair,2)<<endl;
		pcfinout<<bond_pairs_list(pair,3)<<endl;
	}
	for (int pair=0;pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<pairpair_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_pairs_list.rows(); pair++) {
		pcfinout<<pairpair_dev(pair)<<endl;
	}
	pcfinout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(pairpair_result);
#endif

	return result;
}


ITcVecArray symmetry_sector_measurement_fromfiles_hubbard_MonteCarlo(std::ifstream& inputstates, std::string desc,int nmeasure, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table){
	ITcVecArray result(1);

    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
	//cout<<size<<endl<<rank<<endl;
	int Nsite=psi.N();

//Load states from file
	if (!inputstates) { cout<<"!No input file on "<<rank<<endl; exit(1); }
	std::vector<int> loadstates;
	loadstates.reserve(10000*Nsite);
	std::string str;
	int tempi,co=1;
	while (!std::getline(inputstates,str).eof()) {
		std::istringstream ss(str);
		if(!(ss>>tempi)) { cout<<"Bad input!"; }
		else loadstates.push_back(tempi);
	}
	cout<<"End of loading states"<<endl;
	inputstates.close();
	int nmeas,ndata=loadstates.size();
	cout<<ndata<<endl;
	if (ndata%Nsite) { cout<<"Partial data!!!"<<endl; }
	nmeas=ndata/Nsite;
	if (nmeasure>0) { nmeas=nmeasure; }

	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
#if DO_SZSZ==1
	ITiMatrix bond_site_list;
	make_bonds(bond_site_list, Nsite);
	/////////////////////////////////////
	if (rank==0) {
		cout<<"bond_site_list for S.S: "<<endl<<bond_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITdMatrix SzSz_list(nmeas,bond_site_list.rows());
#endif
#if DO_CHIR==1
#if CHIR_8TRIANGS_PERHEXAGON==1
	ITiMatrix tri_site_list;
	//for Si.(SjxSk) measurement; 8 triangles per hexagon, each clockwise;
	//starting from middle site in small triangles; from upper-left in big-down triangles, bottom-left in big-up
	make_triangles_8perhexagon(tri_site_list,Nsite);
#else
	//for Si.(SjxSk) measurement; NN bonds ij and arbitrary k
	ITiMatrix tri_site_list;
	make_triangles_allfromNNbonds(tri_site_list,Nsite);
#endif
	/////////////////////////////////////
	if (rank==0) {
		cout<<"tri_site_list for S.(SxS): "<<endl<<tri_site_list<<endl;
		cout.flush();
	}
	/////////////////////////////////////
	ITcMatrix Chir_list(nmeas,tri_site_list.rows());
#endif
#if DO_SECTORNORM==1
	ITcVec SecNorm_list(nmeas);
#endif

	int measure_count=0;

	for(int i=0;i<nmeas;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec(Nsite);
			for (int ss=0; ss<Nsite; ss++) {
				state_vec(ss)=loadstates[Nsite*measure_count+ss];
			}
#if DO_SZSZ==1
			ITdVec cur_SzSz=SzSzcorr_hubbard( psi, state_vec, bond_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			SzSz_list.set_row(measure_count,cur_SzSz);
#endif
#if DO_CHIR==1
			ITcVec cur_Chir=SSS_hubbard( psi, state_vec, tri_site_list, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			Chir_list.set_row(measure_count,cur_Chir);
#endif
#if DO_SECTORNORM==1
			SecNorm_list(measure_count)=sectornorm_hubbard( psi, state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
#endif
			//dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			measure_count++;
		}
	}

#if DO_SZSZ==1
	ITdVec SzSz_result(bond_site_list.rows()),SzSz_dev(bond_site_list.rows());
	for(int i=0;i<SzSz_result.size();i++)
	{
		SzSz_result(i)=itpp::mean(SzSz_list.get_col(i));
		SzSz_dev(i)=std::sqrt(itpp::variance(SzSz_list.get_col(i))/SzSz_list.rows());
	}
	cout<<"From rank="<<rank<<": SzSz_mean="<<SzSz_result<<endl<<" stderr="<<SzSz_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream fim;
	fim << "SzSzload"<<desc<<"rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	finout<<bond_site_list.rows()<<endl;
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<bond_site_list(pair,0)<<endl;
		finout<<bond_site_list(pair,1)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_result(pair)<<endl;
	}
	for (int pair=0; pair<bond_site_list.rows(); pair++) {
		finout<<SzSz_dev(pair)<<endl;
	}
	finout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(SzSz_result);
#endif
#if DO_CHIR==1
	ITcVec Chir_result(tri_site_list.rows()),Chir_dev(tri_site_list.rows());
	for(int i=0;i<Chir_result.size();i++)
	{
		Chir_result(i)=itpp::mean(Chir_list.get_col(i));
		Chir_dev(i)=std::sqrt(itpp::variance(Chir_list.get_col(i))/Chir_list.rows());
	}
	cout<<"From rank="<<rank<<": S.(SxS)_mean="<<Chir_result<<endl<<" stderr="<<Chir_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream cfim;
	cfim << "Chirload"<<desc<<"rank"<<rank<<".txt";
	std::string cime = cfim.str();
	std::ofstream cfinout;
	cfinout.open(cime.c_str());
	cfinout<<tri_site_list.rows()<<endl;
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<tri_site_list(pair,0)<<endl;
		cfinout<<tri_site_list(pair,1)<<endl;
		cfinout<<tri_site_list(pair,2)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_result(pair)<<endl;
	}
	for (int pair=0; pair<tri_site_list.rows(); pair++) {
		cfinout<<Chir_dev(pair)<<endl;
	}
	cfinout.close();
	/////////////////////////////////////
	result(0)=itpp::to_cvec(Chir_result);
#endif
#if DO_SECTORNORM==1
	dComplex SecNorm_result=itpp::mean(SecNorm_list);
	dComplex SecNorm_dev=std::sqrt(itpp::variance(SecNorm_list));

	cout<<"From rank="<<rank<<": SecNorm_mean="<<SecNorm_result<<endl<<" stderr="<<SecNorm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "SecNormload"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<SecNorm_result<<endl<<SecNorm_dev<<endl;
	for (int pair=0; pair<SecNorm_list.size(); pair++) {
		nfinout<<SecNorm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(SecNorm_result);
	result(0)=temp;
#endif

	return result;
}



dComplex norm_tj(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	dComplex sym_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex amp=quantum_amplitude(psi,st);
	/////////////////////////////////////
	//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return sym_amp/amp;
}


dComplex norm_hubbard(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
	dComplex sym_amp=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex amp=quantum_amplitude(psi,st);
	/////////////////////////////////////
	//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return sym_amp/amp;
}


//This is <sympsi|psi> / <psi|psi>, different from above "sectornorm" <psi|sympsi> / <sympsi|sympsi>
ITcVecArray symmetry_norm_measurement_tj_MonteCarlo(std::string desc, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table){
	ITcVecArray result(1);
	
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
    int Nsite=psi.N();
    
    
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	state_storage_tj cur_state=initialize_random_state_vec_tj_rengine(psi,Nferm,_rengine);

	ITcVec Norm_list(Monte_Carlo_nmeasure);

	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			//if(i%10==0)cout<<"i:"<<i<<endl;             //uncomment if you want to see i printed
			//cout.flush();
		}
		/////////////////////////////////////
		
		try_one_step_tj(psi,cur_state, _rengine);
	}
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "statesUNSYM"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str ();
	std::ofstream stout;
	stout.open (ime1.c_str ());
	/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				//if(j%10==0)cout<<"j:"<<j<<endl;           //uncomment if you want to see j printed
				//cout.flush();
			}
			/////////////////////////////////////
			try_one_step_tj(psi,cur_state, _rengine);
		}
		/////////////////////////////////////
		if (rank==0) {
			//cout<<"i:"<<i<<endl;                         //uncomment if you want to see i printed
			//cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
			Norm_list(measure_count)=norm_tj( psi, state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif

	dComplex Norm_result=itpp::mean(Norm_list);
	dComplex Norm_dev=std::sqrt(itpp::variance(Norm_list));
	
	cout<<"From rank="<<rank<<": Norm_mean="<<Norm_result<<endl<<" stderr="<<Norm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "NormUNSYM"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<Norm_result<<endl<<Norm_dev<<endl;
	for (int pair=0; pair<Norm_list.size(); pair++) {
		nfinout<<Norm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(Norm_result);
	result(0)=temp;

	return result;
}


//This is <sympsi|psi> / <psi|psi>, different from above "sectornorm" <psi|sympsi> / <sympsi|sympsi>
ITcVecArray symmetry_norm_measurement_hubbard_MonteCarlo(std::string desc, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table){
	ITcVecArray result(1);
	
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
    int Nsite=psi.N();
    
    
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	state_storage_hubbard cur_state=initialize_random_state_vec_hubbard_rengine(psi,Nferm,_rengine);
	
	ITcVec Norm_list(Monte_Carlo_nmeasure);
	
	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			if(i%10==0)cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		
		try_one_step_hubbard(psi,cur_state, _rengine);
	}
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "statesUNSYM"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str ();
	std::ofstream stout;
	stout.open (ime1.c_str ());
	/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				if(j%10==0)cout<<"j:"<<j<<endl;
				cout.flush();
			}
			/////////////////////////////////////
			try_one_step_hubbard(psi,cur_state, _rengine);
		}
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
			Norm_list(measure_count)=norm_hubbard( psi, state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif
	
	dComplex Norm_result=itpp::mean(Norm_list);
	dComplex Norm_dev=std::sqrt(itpp::variance(Norm_list));
	
	cout<<"From rank="<<rank<<": Norm_mean="<<Norm_result<<endl<<" stderr="<<Norm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "NormUNSYM"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<Norm_result<<endl<<Norm_dev<<endl;
	for (int pair=0; pair<Norm_list.size(); pair++) {
		nfinout<<Norm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(Norm_result);
	result(0)=temp;
	
	return result;
}


dComplex norm2sym_tj(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table, const ITiVecArray & sym_sector_group_site_mapping_table2, const ITcVec & sym_sector_group_eigval_table2){
	dComplex sym_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex sym_amp2=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table2, sym_sector_group_eigval_table2);
	/////////////////////////////////////
	//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return sym_amp2/sym_amp;
}


dComplex norm2sym_hubbard(const IQMPS & psi, const ITiVec & st, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table, const ITiVecArray & sym_sector_group_site_mapping_table2, const ITcVec & sym_sector_group_eigval_table2){
	dComplex sym_amp=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	dComplex sym_amp2=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table2, sym_sector_group_eigval_table2);
	/////////////////////////////////////
	//	if (rank==0) { cout<<endl;}
	/////////////////////////////////////
	return sym_amp2/sym_amp;
}


//This is <sym2psi|sympsi> / <sympsi|sympsi>
ITcVecArray symmetry2_norm_measurement_tj_MonteCarlo(std::string desc, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table, const ITiVecArray & sym_sector_generator_site_mapping_table2, const ITcVec & sym_sector_generator_eigval_table2){
	ITcVecArray result(1);
	
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
    int Nsite=psi.N();
    
    
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	ITiVecArray sym_sector_group_site_mapping_table2;
	ITcVec sym_sector_group_eigval_table2;
	generate_symmetry_group(sym_sector_generator_site_mapping_table2, sym_sector_generator_eigval_table2, sym_sector_group_site_mapping_table2, sym_sector_group_eigval_table2);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	state_storage_tj cur_state=initialize_random_state_vec_tj_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	ITcVec Norm_list(Monte_Carlo_nmeasure);
	
	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			if(i%10==0)cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		
		try_one_step_tj_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	}
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "statesSYM1_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str ();
	std::ofstream stout;
	stout.open (ime1.c_str ());
	/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				if(j%10==0)cout<<"j:"<<j<<endl;
				cout.flush();
			}
			/////////////////////////////////////
		try_one_step_tj_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
		}
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
			Norm_list(measure_count)=norm2sym_tj(psi,state_vec,sym_sector_group_site_mapping_table,sym_sector_group_eigval_table,sym_sector_group_site_mapping_table2,sym_sector_group_eigval_table2);
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif
	
	dComplex Norm_result=itpp::mean(Norm_list);
	dComplex Norm_dev=std::sqrt(itpp::variance(Norm_list));
	
	cout<<"From rank="<<rank<<": Norm_mean="<<Norm_result<<endl<<" stderr="<<Norm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "NormSYM1SYM2"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<Norm_result<<endl<<Norm_dev<<endl;
	for (int pair=0; pair<Norm_list.size(); pair++) {
		nfinout<<Norm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(Norm_result);
	result(0)=temp;
	
	return result;
}


//This is <sym2psi|sympsi> / <sympsi|sympsi>
ITcVecArray symmetry2_norm_measurement_hubbard_MonteCarlo(std::string desc, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table, const ITiVecArray & sym_sector_generator_site_mapping_table2, const ITcVec & sym_sector_generator_eigval_table2){
	ITcVecArray result(1);
	
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);
    
    int Nsite=psi.N();
    
    
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	ITiVecArray sym_sector_group_site_mapping_table2;
	ITcVec sym_sector_group_eigval_table2;
	generate_symmetry_group(sym_sector_generator_site_mapping_table2, sym_sector_generator_eigval_table2, sym_sector_group_site_mapping_table2, sym_sector_group_eigval_table2);
	
	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	state_storage_hubbard cur_state=initialize_random_state_vec_hubbard_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	
	ITcVec Norm_list(Monte_Carlo_nmeasure);
	
	int measure_count=0;
	
	for(int i=0;i<Monte_Carlo_ntherm;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			if(i%10==0)cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		
		try_one_step_hubbard_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	}
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "statesSYM1_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1 = fim1.str ();
	std::ofstream stout;
	stout.open (ime1.c_str ());
	/////////////////////////////////////
#endif
	for(int i=0;i<Monte_Carlo_nmeasure;i++)
	{
		for(int j=0;j<Monte_Carlo_n_between_measure;j++)
		{
			/////////////////////////////////////
			if (rank==0) {
				if(j%10==0)cout<<"j:"<<j<<endl;
				cout.flush();
			}
			/////////////////////////////////////
			try_one_step_hubbard_symmetrized(psi,cur_state, _rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
		}
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec=cur_state._state_vec;
#if WRITE_STATES==1
			for (int u=0; u<Nsite; u++) {
				stout<<state_vec[u]<<endl;
			}
#endif
			Norm_list(measure_count)=norm2sym_hubbard(psi,state_vec,sym_sector_group_site_mapping_table,sym_sector_group_eigval_table,sym_sector_group_site_mapping_table2,sym_sector_group_eigval_table2);
			measure_count++;
		}
	}
#if WRITE_STATES==1
	stout.close();
#endif
	
	dComplex Norm_result=itpp::mean(Norm_list);
	dComplex Norm_dev=std::sqrt(itpp::variance(Norm_list));
	
	cout<<"From rank="<<rank<<": Norm_mean="<<Norm_result<<endl<<" stderr="<<Norm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "NormSYM1SYM2"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<Norm_result<<endl<<Norm_dev<<endl;
	for (int pair=0; pair<Norm_list.size(); pair++) {
		nfinout<<Norm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(Norm_result);
	result(0)=temp;
	
	return result;
}


//This is <sym2psi|sympsi> / <sympsi|sympsi>, loads sympsi MC states from files; "SYM1"=sym, "SYM2"=sym2
ITcVecArray symmetry2_norm_measurement_fromfiles_tj_MonteCarlo(std::ifstream& inputstates, std::string desc,int nmeasure, const IQMPS & psi, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table, const ITiVecArray & sym_sector_generator_site_mapping_table2, const ITcVec & sym_sector_generator_eigval_table2){
	ITcVecArray result(1);
	
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
	//cout<<size<<endl<<rank<<endl;
	int Nsite=psi.N();
	
	//Load states from file
	if (!inputstates) { cout<<"!No input file on "<<rank<<endl; exit(1); }
	std::vector<int> loadstates;
	loadstates.reserve(10000*Nsite);
	std::string str;
	int tempi,co=1;
	while (!std::getline(inputstates,str).eof()) {
		std::istringstream ss(str);
		if(!(ss>>tempi)) { cout<<"Bad input!"; }
		else loadstates.push_back(tempi);
	}
	cout<<"End of loading states"<<endl;
	inputstates.close();
	int nmeas,ndata=loadstates.size();
	cout<<ndata<<endl;
	if (ndata%Nsite) { cout<<"Partial data!!!"<<endl; }
	nmeas=ndata/Nsite;
	if (nmeasure>0) { nmeas=nmeasure; }
	
	ITiVecArray sym_sector_group_site_mapping_table;
	ITcVec sym_sector_group_eigval_table;
	generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
	ITiVecArray sym_sector_group_site_mapping_table2;
	ITcVec sym_sector_group_eigval_table2;
	generate_symmetry_group(sym_sector_generator_site_mapping_table2, sym_sector_generator_eigval_table2, sym_sector_group_site_mapping_table2, sym_sector_group_eigval_table2);

	/////////////////////////////////////
	if (rank==0) {
		cout<<sym_sector_group_site_mapping_table<<endl;
		cout<<sym_sector_group_eigval_table<<endl;
	}
	/////////////////////////////////////
	ITcVec Norm_list(Monte_Carlo_nmeasure);
	
	int measure_count=0;
	
	for(int i=0;i<nmeas;i++)
	{
		/////////////////////////////////////
		if (rank==0) {
			cout<<"i:"<<i<<endl;
			cout.flush();
		}
		/////////////////////////////////////
		//do one measurement
		{
			ITiVec state_vec(Nsite);
			for (int ss=0; ss<Nsite; ss++) {
				state_vec(ss)=loadstates[Nsite*measure_count+ss];
			}
			Norm_list(measure_count)=norm2sym_tj(psi,state_vec,sym_sector_group_site_mapping_table,sym_sector_group_eigval_table,sym_sector_group_site_mapping_table2,sym_sector_group_eigval_table2);
			//dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
			measure_count++;
		}
	}
	
	dComplex Norm_result=itpp::mean(Norm_list);
	dComplex Norm_dev=std::sqrt(itpp::variance(Norm_list));
	
	cout<<"From rank="<<rank<<": Norm_mean="<<Norm_result<<endl<<" stderr="<<Norm_dev<<endl;
	
	/////////////////////////////////////
	std::stringstream nfim;
	nfim << "NormloadSYM1SYM2"<<desc<<"rank"<<rank<<".txt";
	std::string nime = nfim.str();
	std::ofstream nfinout;
	nfinout.open(nime.c_str());
	nfinout<<Norm_result<<endl<<Norm_dev<<endl;
	for (int pair=0; pair<Norm_list.size(); pair++) {
		nfinout<<Norm_list(pair)<<endl;
	}
	nfinout.close();
	/////////////////////////////////////
	ITcVec temp(1);
	temp(0)=std::abs(Norm_result);
	result(0)=temp;
	
	return result;
}



/* OLD CODE for making chirality triangles, only for 32 site
 #if DO_CHIR==1
 #if CHIR_8TRIANGS_PERHEXAGON==1
 //Si.(SjxSk) measurement; 8 triangles per hexagon, each clockwise;
 //starting from middle site in small triangles; from upper-left in big-down triangles, bottom-left in big-up
 ITiMatrix tri_site_list(4*Nsite,3);
 int count=0,j,k;
 for(int i=0;i<Nsite;i++)
 {
 if (!(i%2)) {//A site
 j=i+1;k=pbca2(i+9, Nsite);//triang 1
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=pbca2(i+9,Nsite);k=pbca1l(i-1);//triang 2
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=pbca1l(i-1);k=i+1;//triang 3
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=pbca1r(i,i+2);k=pbca1r(pbca2(i+9, Nsite),pbca2(i+10, Nsite));//big DOWN triangle
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 } else {//B site
 j=pbca1r(i,i+1);k=i-1;//triang 4
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=i-1;k=pbca2(i-9,Nsite);//triang 5
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=pbca2(i-9,Nsite);k=pbca1r(i,i+1);//triang 6
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 j=pbca2(i-8,Nsite);k=pbca1r(i,i+2);//big UP triangle
 tri_site_list(count,0)=i; tri_site_list(count,1)=j; tri_site_list(count,2)=k;
 count++;
 }
 }
 #else
 //Si.(SjxSk) measurement; NN bonds ij and arbitrary k
 ITiMatrix tri_site_list(3*Nsite*(Nsite-2)/2,3);
 int count=0,j;
 
 for(int i=0;i<Nsite;i++)
 {
 if (!(i%2)) {//A site has two NN
 j=i+1;//NN is B in same UC, no pbc
 for(int k=0;k<Nsite;k++)
 {
 if (k!=i && k!=j) {
 tri_site_list(count,0)=i;
 tri_site_list(count,1)=j;
 tri_site_list(count,2)=k;
 count++;
 }
 }
 j=(i+9)%32;//NN is B in row below, pbc along a2
 for(int k=0;k<Nsite;k++)
 {
 if (k!=i && k!=j) {
 tri_site_list(count,0)=i;
 tri_site_list(count,1)=j;
 tri_site_list(count,2)=k;
 count++;
 }
 }
 } else {//B site has one NN in next UC along row
 if (div(i+1,8).quot>div(i,8).quot) {//pbc along a1
 j=i+1-8;
 } else {
 j=i+1;
 }
 for(int k=0;k<Nsite;k++)
 {
 if (k!=i && k!=j) {
 tri_site_list(count,0)=i;
 tri_site_list(count,1)=j;
 tri_site_list(count,2)=k;
 count++;
 }
 }
 }
 }
 #endif
 /////////////////////////////////////
 if (rank==0) {
 cout<<"tri_site_list for S.(SxS): "<<endl<<tri_site_list<<endl;
 cout.flush();
 }
 /////////////////////////////////////
 ITcMatrix Chir_list(nmeas,tri_site_list.rows());
 #endif
*/


/*CHIRALITY with fermion sign, NOT HERMITIAN!
 ITcVec SSS_tj(const IQMPS & psi, const ITiVec & st, const ITiMatrix & tri, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
 ITcVec result(tri.rows());
 result.zeros();
 ITiVec ord(st.size());
 for (int y=0; y<st.size(); y++) { ord(y)=y; }
 dComplex old_amp=symmetrized_quantum_amplitude_tj(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 //dComplex old_amp=quantum_amplitude(psi,st);
 /////////////////////////////////////
 int rank=MPI::COMM_WORLD.Get_rank();
 //	if (rank==0) { cout<<old_amp<<endl<<tri.rows()<<endl; }
 /////////////////////////////////////
 for (int t=0; t<tri.rows(); t++) {
 if (st[tri(t,0)]!=1 && st[tri(t,1)]!=1 && st[tri(t,2)]!=1) {
 if ((st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=6 && (st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=9) {
 ITiVec perm1(ord),nst1(st);
 nst1[tri(t,0)]=st[tri(t,1)];
 nst1[tri(t,1)]=st[tri(t,2)];
 nst1[tri(t,2)]=st[tri(t,0)];
 perm1[tri(t,0)]=ord[tri(t,1)];
 perm1[tri(t,1)]=ord[tri(t,2)];
 perm1[tri(t,2)]=ord[tri(t,0)];
 dComplex new_amp1=symmetrized_quantum_amplitude_tj(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 //dComplex new_amp1=quantum_amplitude(psi,nst1);
 //dP sgn1=sym_mapping_fermion_sign_tj(st,perm1);
 dP sgn1=1.;
 result[t]+=-0.25*Complex_i*sgn1*new_amp1/old_amp;
 
 ITiVec nst2(st),perm2(ord);
 nst2[tri(t,0)]=st[tri(t,2)];
 nst2[tri(t,1)]=st[tri(t,0)];
 nst2[tri(t,2)]=st[tri(t,1)];
 perm2[tri(t,0)]=ord[tri(t,2)];
 perm2[tri(t,1)]=ord[tri(t,0)];
 perm2[tri(t,2)]=ord[tri(t,1)];
 dComplex new_amp2=symmetrized_quantum_amplitude_tj(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 //dComplex new_amp2=quantum_amplitude(psi,nst2);
 //dP sgn2=sym_mapping_fermion_sign_tj(st,perm2);
 dP sgn2=1.;
 result[t]+=0.25*Complex_i*sgn2*new_amp2/old_amp;
 /////////////////////////////////////
 //				if (rank==0) { cout<<"T"<<tri.get_row(t)<<"  :  "<<st<<" "<<old_amp<<" ;  "<<nst1<<" "<<sgn1*new_amp1<<" ;  "<<nst2<<" "<<sgn2*new_amp2<<" ,"<<result[t]<<endl; cout.flush();}
 /////////////////////////////////////
 }
 }
 /////////////////////////////////////
 //		if (rank==0) { cout<<t<<","; cout.flush();}
 /////////////////////////////////////
 }
 /////////////////////////////////////
 //	if (rank==0) { cout<<endl;}
 /////////////////////////////////////
 return result;
 }
 
 
 ITcVec SSS_hubbard(const IQMPS & psi, const ITiVec & st, const ITiMatrix & tri, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table){
 ITcVec result(tri.rows());
 result.zeros();
 ITiVec ord(st.size());
 for (int y=0; y<st.size(); y++) { ord(y)=y; }
 dComplex old_amp=symmetrized_quantum_amplitude_hubbard(psi,st,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 /////////////////////////////////////
 int rank=MPI::COMM_WORLD.Get_rank();
 //	if (rank==0) { cout<<old_amp<<endl<<tri.rows()<<endl; }
 /////////////////////////////////////
 for (int t=0; t<tri.rows(); t++) {
 if (st[tri(t,0)]!=1 && st[tri(t,1)]!=1 && st[tri(t,2)]!=1 && st[tri(t,0)]!=4 && st[tri(t,1)]!=4 && st[tri(t,2)]!=4) {
 if ((st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=6 && (st[tri(t,0)]+st[tri(t,1)]+st[tri(t,2)])!=9) {
 ITiVec perm1(ord),nst1(st);
 nst1[tri(t,0)]=st[tri(t,1)];
 nst1[tri(t,1)]=st[tri(t,2)];
 nst1[tri(t,2)]=st[tri(t,0)];
 perm1[tri(t,0)]=ord[tri(t,1)];
 perm1[tri(t,1)]=ord[tri(t,2)];
 perm1[tri(t,2)]=ord[tri(t,0)];
 dComplex new_amp1=symmetrized_quantum_amplitude_hubbard(psi,nst1,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 dP sgn1=sym_mapping_fermion_sign_hubbard(st,perm1);
 result[t]+=-0.25*Complex_i*sgn1*std::conj(new_amp1)/old_amp;
 
 ITiVec nst2(st),perm2(ord);
 nst2[tri(t,0)]=st[tri(t,2)];
 nst2[tri(t,1)]=st[tri(t,0)];
 nst2[tri(t,2)]=st[tri(t,1)];
 perm2[tri(t,0)]=ord[tri(t,2)];
 perm2[tri(t,1)]=ord[tri(t,0)];
 perm2[tri(t,2)]=ord[tri(t,1)];
 dComplex new_amp2=symmetrized_quantum_amplitude_hubbard(psi,nst2,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
 dP sgn2=sym_mapping_fermion_sign_hubbard(st,perm2);
 result[t]+=0.25*Complex_i*sgn2*std::conj(new_amp2)/old_amp;
 }
 }
 /////////////////////////////////////
 if (rank==0) { cout<<t<<","; cout.flush();}
 /////////////////////////////////////
 }
 /////////////////////////////////////
 if (rank==0) { cout<<endl;}
 /////////////////////////////////////
 return result;
 }
*/


#include "symMonteCarlo.h"


const extern int Nferm;

dComplex quantum_amplitude(const IQMPS & psi, const ITiVec & state_vec)
{
  int N=psi.N();

  IQTensor iqtt;
  for(int i=1;i<=N;i++)
  {
    IQTensor cur_site_iqt=psi.A(i);
    //cout<<"site: "<<i<<" is complex? "<<cur_site_iqt.isComplex()<<endl;
    IQIndex cur_site_ind=findtype(cur_site_iqt, Site);
    IQTensor cur_site_state(conj(cur_site_ind));
    cur_site_state(cur_site_ind(state_vec(i-1)))=1.;
    cur_site_iqt=cur_site_iqt*cur_site_state;
    if(i==1)
    {
      iqtt=cur_site_iqt;
    }
    else
    {
      iqtt=iqtt*cur_site_iqt;
    }
  }
  return iqtt.toComplex();
}

// dP sym_mapping_fermion_sign_hubbard(const ITiVec & old_state_vec, const ITiVec & sym_site_mapping_table)
// {
//   int Nsite=old_state_vec.size();
//   //state-1: empty;state-2/3: up/down; state-4: double occupy
//   ITiVec Nf_vec="0,1,1,2";
//   int Nf=0;
//   for(int i=0;i<Nsite;i++)
//   {
//     Nf+=Nf_vec(old_state_vec(i)-1);
//   }
//   ITiVec mapped_new_state(Nf);
//   int count=0;
//   for(int i=0;i<Nsite;i++)
//   {
//     if(old_state_vec(i)==2)//if old state is up
//     {
//       mapped_new_state(count)=2*sym_site_mapping_table(i);
//       count++;
//     }
//     else if(old_state_vec(i)==3)//if old state is down
//     {
//       mapped_new_state(count)=2*sym_site_mapping_table(i)+1;
//       count++;
//     }
//     else if(old_state_vec(i)==4)//if old state is downup
//     {
//       mapped_new_state(count)=2*sym_site_mapping_table(i);
//       count++;
//       mapped_new_state(count)=2*sym_site_mapping_table(i)+1;
//       count++;
//     }
//   }
//   
//   //cout<< mapped_new_state<<endl;
//   
//   //compute sign of permutation
//   int N = mapped_new_state.size();
//   bool neg = true;
//   for(int i = 0; i<N-1; i++) {
//     for(int j = i+1; j<N; j++) {
//         if(mapped_new_state(i) > mapped_new_state(j)) neg = !neg;
//     }
//   }
//   
//   if(neg)
//   {
//     return 1.;
//   }
//   else
//   {
//     return -1.;
//   }
// }

// dP sym_mapping_fermion_sign_tj(const ITiVec & old_state_vec, const ITiVec & sym_site_mapping_table)
// {
//   return sym_mapping_fermion_sign_hubbard(old_state_vec,sym_site_mapping_table);
// }



dP sym_mapping_fermion_sign_hubbard(const ITiVec & old_state_vec, const ITiVec & sym_site_mapping_table)
{
        int N=old_state_vec.size();
	std::vector<int> sss(N),perm(N);
	for(int i=0;i<N;i++)
	{
	  sss[i]=old_state_vec(i);
	  perm[i]=sym_site_mapping_table(i);
	}
	
  	int sig=1,i,temp;
	int j,inc,k;
	std::vector<int> st(N),order(N),porder(N),ord(3*N/4);
	//need state vector, having 0 for hole or doubleocc, and +/-1 for up/down, for fermion sign; sign varies, but for tJ-like states should be always -1 for 2x2x2 and +1 for 4x4x2 due to reversed permutation
	for (int y=0; y<N; y++) {
		if (sss[y]==3) { st[y]=-1; } else if (sss[y]==4) { st[y]=0; } else { st[y]=sss[y]-1; }
	}

	int c=1;
	for (k=0; k<N; k++) {
		switch (st[k]) {
			case 0:
				order[k]=0;
				break;
			default:
				order[k]=c;
				c++;
				break;
		}
	}
	
	for (k=0; k<N; k++) {
		porder[k]=order[perm[k]];
	}
	
	int d=0,v;
	for (k=0; k<N; k++) {
		v=porder[k];
		if (v!=0) { ord[d]=v; d++; }
	}
	
	for(inc=Nferm/2; inc>0; inc /= 2){
        for(i=inc; i<Nferm; i++){
            temp=ord[i];
            for(j=i; j>=inc; j -= inc){
                if(temp < ord[j-inc]) {
                    ord[j] = ord[j-inc];
					sig*=-1;}
                else break;
            }
            ord[j] = temp;
        }
    }
	return(1.*sig);

}


dP sym_mapping_fermion_sign_tj(const ITiVec & old_state_vec, const ITiVec & sym_site_mapping_table){
/////////////////////////////////////
//	if (MPI::COMM_WORLD.Get_rank()==14) {
//		cout<<old_state_vec<<endl;
//		cout<<sym_site_mapping_table<<endl;
//	}
/////////////////////////////////////
	int N=old_state_vec.size();
	std::vector<int> sss(N),perm(N);
	for(int i=0;i<N;i++)
	{
	  sss[i]=old_state_vec(i);
	  perm[i]=sym_site_mapping_table(i);
	}
	int sig=1,i,temp;
	int j,inc,k;
	std::vector<int> st(N),order(N),porder(N),ord(Nferm);
	//need state vector, having 0 for hole, and +/-1 for up/down, for fermion sign; sign should be always -1 for 2x2x2 and +1 for 4x4x2 due to reversed permutation
	for (int y=0; y<N; y++) {
		if (sss[y]==3) { st[y]=-1; } else { st[y]=sss[y]-1; }
	}
/////////////////////////////////////
//	if (MPI::COMM_WORLD.Get_rank()==14) {
//		for (int tr=0;tr<N;tr++){ cout<<st[tr]<<",";}; cout<<endl;
//	}
/////////////////////////////////////
	int c=1;
	for (k=0; k<N; k++) {
		switch (st[k]) {
			case 0:
				order[k]=0;
				break;
			default:
				order[k]=c;
				c++;
				break;
		}
	}
/////////////////////////////////////
//	if (MPI::COMM_WORLD.Get_rank()==14) {
//		for (int tr=0;tr<N;tr++){ cout<<order[tr]<<",,";}; cout<<endl;
//	}
/////////////////////////////////////
	
	for (k=0; k<N; k++) {
		porder[k]=order[perm[k]];
	}
	
	int d=0,v;
	for (k=0; k<N; k++) {
		v=porder[k];
		if (v!=0) { ord[d]=v; d++; }
	}
/////////////////////////////////////
//	if (MPI::COMM_WORLD.Get_rank()==14) {
//		for (int tr=0;tr<ord.size();tr++){ cout<<ord[tr]<<";";}; cout<<endl;
//	}
/////////////////////////////////////

	for(inc=Nferm/2; inc>0; inc /= 2){
        for(i=inc; i<Nferm; i++){
            temp=ord[i];
            for(j=i; j>=inc; j -= inc){
                if(temp < ord[j-inc]) {
                    ord[j] = ord[j-inc];
					sig*=-1;}
                else break;
            }
            ord[j] = temp;
        }
    }
	return(1.*sig);
}



//total Sz is set to 0.
ITiVec initialize_random_state_vec_hubbard(int Nsite, int Nf)
{
  //state-1: empty;state-2/3: up/down; state-4: double occupy
  ITiVec Nf_vec="0,1,1,2";
  
  itpp::RNG_randomize();
  //itpp::RNG_reset(0);
  
  
  
  //up site initializing
  ITiVec up_state_vec(Nsite);
  up_state_vec.zeros();
  int count=0;
  do{
    //choose a site
    int cur_site=itpp::randi(0,Nsite-1);
    //set its state
    if(up_state_vec(cur_site)==0)
    {
      up_state_vec(cur_site)=1;
      count++;
    }
  }
  while(count<Nf/2);
  //down site initializing
  ITiVec down_state_vec(Nsite);
  down_state_vec.zeros();
  count=0;
  do{
    //choose a site
    int cur_site=itpp::randi(0,Nsite-1);
    //set its state
    if(down_state_vec(cur_site)==0)
    {
      down_state_vec(cur_site)=1;
      count++;
    }
  }
  while(count<Nf/2);
  
  ITiVec state_vec(Nsite);
  for(int i=0;i<Nsite;i++)
  {
    if(up_state_vec(i)==0)
    {
      if(down_state_vec(i)==0)
      {
	state_vec(i)=1;
      }
      else
      {
	state_vec(i)=3;
      }
    }
    else if(up_state_vec(i)==1)
    {
      if(down_state_vec(i)==0)
      {
	state_vec(i)=2;
      }
      else
      {
	state_vec(i)=4;
      }
    }
  }
  
  return state_vec;
}

ITiVec initialize_random_state_vec_tj(int Nsite, int Nf)
{
  
    itpp::RNG_randomize();
    //itpp::RNG_reset(0);

  
  //up site initializing
  ITiVec state_vec(Nsite);
  state_vec.ones();
  int count=0;
  do{
    //choose a site
    int cur_site=itpp::randi(0,Nsite-1);
    //set its state
    if(state_vec(cur_site)==1)
    {
      state_vec(cur_site)=2;
      count++;
    }
  }
  while(count<Nf/2);
  //down site initializing
  count=0;
  do{
    //choose a site
    int cur_site=itpp::randi(0,Nsite-1);
    //set its state
    if(state_vec(cur_site)==1)
    {
      state_vec(cur_site)=3;
      count++;
    }
  }
  while(count<Nf/2);
  
  return state_vec;
}

dComplex symmetry_eigval_hubbard(const IQMPS & psi, const ITiVec & sym_site_mapping_table)
{
  
  ITiVec state_vec=initialize_random_state_vec_hubbard(psi.N(),Nferm);
  
  dComplex old_amp=quantum_amplitude(psi, state_vec);
  cout<<"old amp:"<<old_amp<<endl;
  
  ITiVec new_state_vec(state_vec.size());
  for(int i=0;i<new_state_vec.size();i++)
  {
    new_state_vec(i)=state_vec(sym_site_mapping_table(i));
  }
  
  dComplex new_amp=quantum_amplitude(psi,new_state_vec);
  cout<<"new amp:"<<new_amp<<endl;
  
  //cout<<"sign of permutation: "<< sym_mapping_fermion_sign_hubbard (state_vec, sym_site_mapping_table)<<endl;
  return sym_mapping_fermion_sign_hubbard (state_vec, sym_site_mapping_table)*new_amp/old_amp;
}

//return the mapping table by firstly applying table1,then table2.
//result(i)=table2(table1(i));
ITiVec sym_multiply_site_mapping_table(const ITiVec & table1, const ITiVec & table2)
{
  ITiVec result(table1.size());
  for(int i=0;i<result.size();i++)
  {
    result(i)=table2(table1(i));
  }
  return result;
}

void generate_symmetry_group(const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table, ITiVecArray & sym_sector_group_site_mapping_table, ITcVec & sym_sector_group_eigval_table)
{
	
	ITiVecArrayArray sym_sector_each_generator_site_mapping_table(sym_sector_generator_site_mapping_table.size());
	ITcVecArray sym_sector_each_generator_eigval_table(sym_sector_generator_site_mapping_table.size());
	for(int i=0;i<sym_sector_generator_site_mapping_table.size();i++)
	{
		ITiVecArray cur_sym_sector_generator_site_mapping_table(0);
		ITcVec cur_sym_sector_generator_eigval_table(0);
		//initial sym element=identity:
		ITiVec id_site_mapping_list(sym_sector_generator_site_mapping_table(i).size());
		for(int j=0;j<id_site_mapping_list.size();j++)
		{
			id_site_mapping_list(j)=j;
		}
		ITiVec cur_site_mapping_list=id_site_mapping_list;
		dComplex cur_eigval=1.;
		do
		{
			cur_sym_sector_generator_site_mapping_table=itpp::concat(cur_sym_sector_generator_site_mapping_table,cur_site_mapping_list);
			cur_sym_sector_generator_eigval_table=itpp::concat(cur_sym_sector_generator_eigval_table,cur_eigval);
			//apply the current generator once
			cur_site_mapping_list=sym_multiply_site_mapping_table(cur_site_mapping_list,sym_sector_generator_site_mapping_table(i));
			
			cur_eigval*=sym_sector_generator_eigval_table(i);
		}
		while(cur_site_mapping_list!=id_site_mapping_list);
		sym_sector_each_generator_site_mapping_table(i)=cur_sym_sector_generator_site_mapping_table;
		sym_sector_each_generator_eigval_table(i)=cur_sym_sector_generator_eigval_table;
	}
	
	//cout<<sym_sector_each_generator_site_mapping_table<<endl;
	//cout<<sym_sector_each_generator_eigval_table<<endl;
	
	if(sym_sector_generator_site_mapping_table.size()==1)//if one generator
	{
		sym_sector_group_site_mapping_table=sym_sector_each_generator_site_mapping_table(0);
		sym_sector_group_eigval_table=sym_sector_each_generator_eigval_table(0);
	}
	else if(sym_sector_generator_site_mapping_table.size()==2)//if two generators, do product
	{
		sym_sector_group_site_mapping_table.set_size(sym_sector_each_generator_site_mapping_table(0).size()*sym_sector_each_generator_site_mapping_table(1).size());
		sym_sector_group_eigval_table.set_size(sym_sector_group_site_mapping_table.size());
        cout<<"sym_sector_each_generator_group:"<<endl<<sym_sector_each_generator_site_mapping_table<<endl;
		int count=0;
		for(int i=0;i<sym_sector_each_generator_site_mapping_table(0).size();i++)
		{
			for(int j=0;j<sym_sector_each_generator_site_mapping_table(1).size();j++)
			{
				sym_sector_group_site_mapping_table(count)=sym_multiply_site_mapping_table(sym_sector_each_generator_site_mapping_table(0)(i),sym_sector_each_generator_site_mapping_table(1)(j));
				sym_sector_group_eigval_table(count)=sym_sector_each_generator_eigval_table(0)(i)*sym_sector_each_generator_eigval_table(1)(j);
				count++;
			}
		}
	}
	else if(sym_sector_generator_site_mapping_table.size()==3)//if three generators, do product
	{
		sym_sector_group_site_mapping_table.set_size(sym_sector_each_generator_site_mapping_table(0).size()*sym_sector_each_generator_site_mapping_table(1).size()*sym_sector_each_generator_site_mapping_table(2).size());
		sym_sector_group_eigval_table.set_size(sym_sector_group_site_mapping_table.size());
		cout<<"sym_sector_each_generator_group:"<<endl<<sym_sector_each_generator_site_mapping_table<<endl;
		int count=0;
		for(int i=0;i<sym_sector_each_generator_site_mapping_table(0).size();i++)
		{
			for(int j=0;j<sym_sector_each_generator_site_mapping_table(1).size();j++)
			{
				for(int k=0;k<sym_sector_each_generator_site_mapping_table(2).size();k++)
				{
					sym_sector_group_site_mapping_table(count)=sym_multiply_site_mapping_table(sym_sector_each_generator_site_mapping_table(0)(i),sym_sector_each_generator_site_mapping_table(1)(j));
					sym_sector_group_site_mapping_table(count)=sym_multiply_site_mapping_table(sym_sector_group_site_mapping_table(count),sym_sector_each_generator_site_mapping_table(2)(k));
					
					sym_sector_group_eigval_table(count)=sym_sector_each_generator_eigval_table(0)(i)*sym_sector_each_generator_eigval_table(1)(j)*sym_sector_each_generator_eigval_table(2)(k);
					
					count++;
				}
			}
		}
	}
	else
	{
		cout<<"ERROR: too many generators in generate_symmetry_group function."<<endl;
		exit(0);
	}
	
	//check symmetry group elements are independent.
	//remove duplicates
	ITiVec duplicate_q(sym_sector_group_site_mapping_table.size());
	duplicate_q.zeros();
	for(int i=0;i<sym_sector_group_site_mapping_table.size();i++)
	{
		for(int j=i+1;j<sym_sector_group_site_mapping_table.size();j++)
		{
			if(sym_sector_group_site_mapping_table(i)==sym_sector_group_site_mapping_table(j))
			{
				cout<<"Warning: different group elements same in generate_symmetry_group function."<<endl;
				duplicate_q(j)=1;
			}
		}
	}
	int n_removed=itpp::sum(duplicate_q);
	ITiVecArray new_sym_sector_group_site_mapping_table(sym_sector_group_site_mapping_table.size()-n_removed);
	ITcVec new_sym_sector_group_eigval_table(sym_sector_group_site_mapping_table.size()-n_removed);
	int count=0;
	for(int i=0;i<sym_sector_group_site_mapping_table.size();i++)
	{
		if(duplicate_q(i)==0)
		{
			new_sym_sector_group_site_mapping_table(count)=sym_sector_group_site_mapping_table(i);
			new_sym_sector_group_eigval_table(count)=sym_sector_group_eigval_table(i);
			count++;
		}
	}
	
	sym_sector_group_site_mapping_table=new_sym_sector_group_site_mapping_table;
	sym_sector_group_eigval_table=new_sym_sector_group_eigval_table;
	
	cout<<"Group elements generated."<<endl;
	return;
}


void  find_sym_states_with_coeff_hubbard(const ITiVec & state_vec, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec &sym_sector_group_eigval_table, ITiVecArray & sym_state_vec_list, ITiVec & sym_state_vec_union_ind, ITcVec & sym_state_coeff)
{
  sym_state_vec_list.set_size(sym_sector_group_site_mapping_table.size());
  for(int i=0;i<sym_state_vec_list.size();i++)
  {
    sym_state_vec_list(i)=sym_multiply_site_mapping_table(sym_sector_group_site_mapping_table(i),state_vec);
  }

  
  
  sym_state_vec_union_ind.set_size(0);
  sym_state_coeff.set_size(0);
  
  ITiVec sym_state_ind_list(sym_state_vec_list.size());
  for(int i=0;i<sym_state_ind_list.size();i++)
  {
    sym_state_ind_list(i)=i;
  }
  do{
    ITiVec cur_state_vec=sym_state_vec_list(sym_state_ind_list(0));
    ITiVec find_q(sym_state_ind_list.size());
    find_q.zeros();
    for(int i=0;i<sym_state_ind_list.size();i++)
    {
		if(sym_state_vec_list(sym_state_ind_list(i))==cur_state_vec)
		{
			find_q(i)=1;
		}
    }
    int n_found=itpp::sum(find_q);
    dComplex coeff=0.;
    for(int i=0;i<find_q.size();i++)
    {
      if(find_q(i)==1)
	
      coeff+=sym_mapping_fermion_sign_hubbard(state_vec, sym_sector_group_site_mapping_table(sym_state_ind_list(i)))*sym_sector_group_eigval_table(sym_state_ind_list(i));
    }
    sym_state_vec_union_ind=itpp::concat(sym_state_vec_union_ind,sym_state_ind_list(0));
    sym_state_coeff=itpp::concat(sym_state_coeff,coeff);
    
    ITiVec new_sym_state_ind_list(sym_state_ind_list.size()-n_found);
    int count=0;
    for(int i=0;i<find_q.size();i++)
    {
		if(find_q(i)==0)
		{
			new_sym_state_ind_list(count)=sym_state_ind_list(i);
			count++;
		}
    }
    
    sym_state_ind_list=new_sym_state_ind_list;
  }
  while(sym_state_ind_list.size()>0);
  
  return;
}

void  find_sym_states_with_coeff_tj(const ITiVec & state_vec, const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec &sym_sector_group_eigval_table, ITiVecArray & sym_state_vec_list, ITiVec & sym_state_vec_union_ind, ITcVec & sym_state_coeff)
{
  sym_state_vec_list.set_size(sym_sector_group_site_mapping_table.size());
  for(int i=0;i<sym_state_vec_list.size();i++)
  {
    sym_state_vec_list(i)=sym_multiply_site_mapping_table(sym_sector_group_site_mapping_table(i),state_vec);
  }

  
  
  sym_state_vec_union_ind.set_size(0);
  sym_state_coeff.set_size(0);
  
  ITiVec sym_state_ind_list(sym_state_vec_list.size());
  for(int i=0;i<sym_state_ind_list.size();i++)
  {
    sym_state_ind_list(i)=i;
  }
  do{
    ITiVec cur_state_vec=sym_state_vec_list(sym_state_ind_list(0));
    ITiVec find_q(sym_state_ind_list.size());
    find_q.zeros();
    for(int i=0;i<sym_state_ind_list.size();i++)
    {
      if(sym_state_vec_list(sym_state_ind_list(i))==cur_state_vec)
      {
	find_q(i)=1;
      }
    }
    int n_found=itpp::sum(find_q);
    dComplex coeff=0.;
    for(int i=0;i<find_q.size();i++)
    {
      if(find_q(i)==1)
	
      coeff+=sym_mapping_fermion_sign_tj(state_vec, sym_sector_group_site_mapping_table(sym_state_ind_list(i)))*sym_sector_group_eigval_table(sym_state_ind_list(i));
    }
    sym_state_vec_union_ind=itpp::concat(sym_state_vec_union_ind,sym_state_ind_list(0));
    sym_state_coeff=itpp::concat(sym_state_coeff,coeff);
    
    ITiVec new_sym_state_ind_list(sym_state_ind_list.size()-n_found);
    int count=0;
    for(int i=0;i<find_q.size();i++)
    {
      if(find_q(i)==0)
      {
	new_sym_state_ind_list(count)=sym_state_ind_list(i);
	count++;
      }
    }
    
    sym_state_ind_list=new_sym_state_ind_list;
  }
  while(sym_state_ind_list.size()>0);
  
  return;
}


dComplex symmetrized_quantum_amplitude_hubbard(const IQMPS & psi, const ITiVec & state_vec, const ITiVecArray &sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  
  ITiVecArray sym_state_vec_list;  
  ITiVec sym_state_vec_union_ind;
  ITcVec sym_state_coeff;
  
  find_sym_states_with_coeff_hubbard(state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table, sym_state_vec_list,sym_state_vec_union_ind,sym_state_coeff);
  
  dComplex result=0.;
  
  for(int i=0;i<sym_state_coeff.size();i++)
  {
    result+=sym_state_coeff(i)*quantum_amplitude(psi,sym_state_vec_list(sym_state_vec_union_ind(i)));
  }
  if(itpp::norm(sym_state_coeff)<1.E-14)
  {
    return 0.;
  }
  else
  {
    return sym_state_coeff(0)/std::pow(itpp::norm(sym_state_coeff),2)*result;
  }
}

dComplex symmetrized_quantum_amplitude_tj(const IQMPS & psi, const ITiVec & state_vec, const ITiVecArray &sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  
  ITiVecArray sym_state_vec_list;  
  ITiVec sym_state_vec_union_ind;
  ITcVec sym_state_coeff;
  
  find_sym_states_with_coeff_tj(state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table, sym_state_vec_list,sym_state_vec_union_ind,sym_state_coeff);
  
  dComplex result=0.;
  
  for(int i=0;i<sym_state_coeff.size();i++)
  {
    result+=sym_state_coeff(i)*quantum_amplitude(psi,sym_state_vec_list(sym_state_vec_union_ind(i)));
  }
  if(itpp::norm(sym_state_coeff)<1.E-14)
  {
    return 0.;
  }
  else
  {
    return sym_state_coeff(0)/std::pow(itpp::norm(sym_state_coeff),2)*result;
  }
}

//assuming total symmetry group is generated by direct product of those abelian groups generated by each generator.
dComplex symmetry_sector_eigval_hubbard(const IQMPS & psi, const ITiVec & sym_site_mapping_table, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table)
{
  ITiVecArray sym_sector_group_site_mapping_table;
  ITcVec sym_sector_group_eigval_table;
  generate_symmetry_group(sym_sector_generator_site_mapping_table, sym_sector_generator_eigval_table, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
  
  ITiVec state_vec=initialize_random_state_vec_hubbard(psi.N(),Nferm);
  
  dComplex old_amp=symmetrized_quantum_amplitude_hubbard( psi,  state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
  ITiVec new_state_vec=sym_multiply_site_mapping_table(sym_site_mapping_table,state_vec);
  dComplex new_amp=symmetrized_quantum_amplitude_hubbard( psi,  new_state_vec, sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
  cout<<"old_amp:"<<old_amp<<endl;
  return sym_mapping_fermion_sign_hubbard (state_vec, sym_site_mapping_table)*new_amp/old_amp;
}

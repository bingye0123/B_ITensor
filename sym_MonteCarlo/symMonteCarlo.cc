
#include "symMonteCarlo.h"


//number of electrons
const extern int Nferm;
//number of Monte Carlo thermalization steps
const extern int Monte_Carlo_ntherm;
//number of Monte Carlo measurement
const extern int Monte_Carlo_nmeasure;
//number of steps between measurement
const extern int Monte_Carlo_n_between_measure;

//write into file the MC states used for measurement
#define WRITE_STATES 0

void generate_state_lists_using_state_vec_hubbard(state_storage_hubbard & state)
{
  ITiVec & state_vec=state._state_vec;
  ITiVec & up_fill_list=state._up_fill_list;
  ITiVec & down_fill_list=state._down_fill_list;
  ITiVec & up_emp_list=state._up_emp_list;
  ITiVec & down_emp_list=state._down_emp_list;
  
  int Nsite=state_vec.size();
  int Nf=0;
  ITiVec Nf_vec="0,1,1,2";
  for(int i=0;i<state_vec.size();i++)
  {
    Nf+=Nf_vec(state_vec(i)-1);
  }
  up_fill_list.set_size(Nf/2);
  down_fill_list.set_size(Nf/2);
  up_emp_list.set_size(Nsite-Nf/2);
  down_emp_list.set_size(Nsite-Nf/2);
  
  int count_up_fill=0,count_down_fill=0,count_up_emp=0,count_down_emp=0;
  
  for(int i=0;i<state_vec.size();i++)
  {
    if(state_vec(i)==1)
    {
      up_emp_list(count_up_emp)=i;
      count_up_emp++;
      down_emp_list(count_down_emp)=i;
      count_down_emp++;
    }
    else if(state_vec(i)==2)
    {
      up_fill_list(count_up_fill)=i;
      count_up_fill++;
      down_emp_list(count_down_emp)=i;
      count_down_emp++;
    }
    else if(state_vec(i)==3)
    {
      up_emp_list(count_up_emp)=i;
      count_up_emp++;
      down_fill_list(count_down_fill)=i;
      count_down_fill++;
    }
    else if(state_vec(i)==4)
    {
      up_fill_list(count_up_fill)=i;
      count_up_fill++;
      down_fill_list(count_down_fill)=i;
      count_down_fill++;
    }
    else
    {
      cout<<"ERROR: wrong state in generate_state_lists_using_state_vec_hubbard."<<endl;
      exit(0);
    }
  }
  
  return;
}

void generate_state_lists_using_state_vec_tj(state_storage_tj & state)
{
  ITiVec & state_vec=state._state_vec;
  ITiVec & up_fill_list=state._up_fill_list;
  ITiVec & down_fill_list=state._down_fill_list;
  ITiVec & emp_list=state._emp_list;
  
  int Nsite=state_vec.size();
  int Nf=0;
  ITiVec Nf_vec="0,1,1";
  for(int i=0;i<state_vec.size();i++)
  {
    Nf+=Nf_vec(state_vec(i)-1);
  }
  up_fill_list.set_size(Nf/2);
  down_fill_list.set_size(Nf/2);
  emp_list.set_size(Nsite-Nf);
  
  int count_up_fill=0,count_down_fill=0,count_emp=0;
  
  for(int i=0;i<state_vec.size();i++)
  {
    if(state_vec(i)==1)
    {
      emp_list(count_emp)=i;
      count_emp++;
    }
    else if(state_vec(i)==2)
    {
      up_fill_list(count_up_fill)=i;
      count_up_fill++;
    }
    else if(state_vec(i)==3)
    {
      down_fill_list(count_down_fill)=i;
      count_down_fill++;
    }
    else
    {
      cout<<"ERROR: wrong state in generate_state_lists_using_state_vec_tj."<<endl;
      exit(0);
    }
  }
  
  return;
}


state_storage_hubbard initialize_random_state_vec_hubbard_rengine(const IQMPS & psi, int Nf, trng::yarn3 & _rengine)
{
  int Nsite=psi.N();
  //state-1: empty;state-2/3: up/down; state-4: double occupy
  ITiVec Nf_vec="0,1,1,2";
  
  trng::uniform_int_dist U(0, Nsite);
  
  //up site initializing
  ITiVec up_state_vec(Nsite);
  up_state_vec.zeros();
  int count=0;
  do{
    //choose a site
    int cur_site=U(_rengine);
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
    int cur_site=U(_rengine);;
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
  
  state_storage_hubbard state;
  state._state_vec=state_vec;
  generate_state_lists_using_state_vec_hubbard(state);
  state._amp=quantum_amplitude(psi,state_vec);
  return state;
}


state_storage_tj initialize_random_state_vec_tj_rengine(const IQMPS & psi, int Nf, trng::yarn3 & _rengine)
{
  int Nsite=psi.N();
  
  trng::uniform_int_dist U(0, Nsite);
  
  //up site initializing
  ITiVec state_vec(Nsite);
  state_vec.ones();
  int count=0;
  do{
    //choose a site
    int cur_site=U(_rengine);
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
    int cur_site=U(_rengine);
    //set its state
    if(state_vec(cur_site)==1)
    {
      state_vec(cur_site)=3;
      count++;
    }
  }
  while(count<Nf/2);
  
  
  
  state_storage_tj state;
  state._state_vec=state_vec;
  generate_state_lists_using_state_vec_tj(state);
  state._amp=quantum_amplitude(psi,state_vec);
  return state;
}


state_storage_hubbard initialize_random_state_vec_hubbard_rengine_symmetrized(const IQMPS & psi, int Nf, trng::yarn3 & _rengine, const ITiVecArray &sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  state_storage_hubbard state=initialize_random_state_vec_hubbard_rengine(psi, Nf, _rengine);
  state._amp=symmetrized_quantum_amplitude_hubbard(psi,state._state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  return state;
}

state_storage_tj initialize_random_state_vec_tj_rengine_symmetrized(const IQMPS & psi, int Nf, trng::yarn3 & _rengine, const ITiVecArray &sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  state_storage_tj state=initialize_random_state_vec_tj_rengine(psi, Nf, _rengine);
  state._amp=symmetrized_quantum_amplitude_tj(psi,state._state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  return state;
}

bool acceptq(const dP probratio, trng::yarn3 & _rengine)
{
  if(probratio>1.)
  {
    return true;
  }
  else
  {
    trng::uniform01_dist<> u;
    dP r=u(_rengine);
    if(r<probratio)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}



SwapTypeHubbard try_one_step_hubbard(const IQMPS & psi, state_storage_hubbard & state, trng::yarn3 & _rengine)
{
  SwapTypeHubbard swt=NONMOVE;
  
  int wayofnewupswaps=state._up_fill_list.size()*state._up_emp_list.size();
  int wayofnewdownswaps=state._down_fill_list.size()*state._down_emp_list.size();
  
  int totwayofswaps=wayofnewupswaps+wayofnewdownswaps;
  
  trng::uniform_int_dist U(0, totwayofswaps);
  int a=U(_rengine);

  int new_up_fill_ind, new_down_fill_ind, new_up_emp_ind, new_down_emp_ind, old_up_fill_ind,old_down_fill_ind,old_up_emp_ind,old_down_emp_ind;
  
  const ITiVec & old_state_vec=state._state_vec;
  ITiVec new_state_vec=old_state_vec;
  
  dComplex new_amp=0.;
  
  if(a<wayofnewupswaps)
  {
    new_up_fill_ind=a/state._up_fill_list.size();
    new_up_emp_ind=a%state._up_fill_list.size();
    old_up_emp_ind=new_up_fill_ind;
    old_up_fill_ind=new_up_emp_ind;
    
    if(old_state_vec(state._up_fill_list(old_up_fill_ind))==2&&old_state_vec(state._up_emp_list(old_up_emp_ind))==1)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=2;
      swt=UP_EMP_to_EMP_UP;
    }
    else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==2&&old_state_vec(state._up_emp_list(old_up_emp_ind))==3)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=4;
      swt=UP_DOWN_to_EMP_DL;
    }
     else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==4&&old_state_vec(state._up_emp_list(old_up_emp_ind))==1)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=2;
      swt=DL_EMP_to_DOWN_UP;
    }
     else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==4&&old_state_vec(state._up_emp_list(old_up_emp_ind))==3)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=4;
      swt=DL_DOWN_to_DOWN_DL;
    }
    else
    {
      cout<<"Error: WRONG "<<swt<< " Swap in try_one_step_hubbard."<<endl;
      exit(0);
    }
    new_amp=quantum_amplitude(psi,new_state_vec);
  }
  
  else
  {
    a-=(wayofnewupswaps);
    new_down_fill_ind=a/state._down_fill_list.size();
    new_down_emp_ind=a%state._down_fill_list.size();
    old_down_emp_ind=new_down_fill_ind;
    old_down_fill_ind=new_down_emp_ind;
    
     if(old_state_vec(state._down_fill_list(old_down_fill_ind))==3&&old_state_vec(state._down_emp_list(old_down_emp_ind))==1)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=3;
      swt=DOWN_EMP_to_EMP_DOWN;
    }
    else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==3&&old_state_vec(state._down_emp_list(old_down_emp_ind))==2)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=4;
      swt=DOWN_UP_to_EMP_DL;
    }
     else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==4&&old_state_vec(state._down_emp_list(old_down_emp_ind))==1)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=3;
      swt=DL_EMP_to_UP_DOWN;
    }
     else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==4&&old_state_vec(state._down_emp_list(old_down_emp_ind))==2)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=4;
      swt=DL_UP_to_UP_DL;
    }
    else
    {
      cout<<"Error: WRONG "<<swt<< " Swap in try_one_step_hubbard."<<endl;
      exit(0);
    }
    new_amp=quantum_amplitude(psi,new_state_vec);
  }
  
  dP probratio=pow(std::abs(new_amp/state._amp),2.);
  if(acceptq(probratio,_rengine))
  {
    state._state_vec=new_state_vec;
    generate_state_lists_using_state_vec_hubbard(state);
    state._amp=new_amp;
    return swt;
  }
  else
  {
    return NONMOVE;
  }
}


SwapTypetJ try_one_step_tj(const IQMPS & psi, state_storage_tj & state, trng::yarn3 & _rengine)
{
  SwapTypetJ swt=NONMOVEtJ;
  
  int wayofnewupswaps=state._up_fill_list.size()*state._emp_list.size();
  int wayofnewdownswaps=state._down_fill_list.size()*state._emp_list.size();
  int waysofupdownswaps=state._up_fill_list.size()*state._down_fill_list.size();
  int totwayofswaps=wayofnewupswaps+wayofnewdownswaps+waysofupdownswaps;
  
  trng::uniform_int_dist U(0, totwayofswaps);
  int a=U(_rengine);

  int new_up_fill_ind, new_down_fill_ind, new_emp_ind, old_up_fill_ind,old_down_fill_ind,old_emp_ind;
  
  const ITiVec & old_state_vec=state._state_vec;
  ITiVec new_state_vec=old_state_vec;
  
  dComplex new_amp=0.;
  
  if(a<wayofnewupswaps)
  {
    new_up_fill_ind=a/state._up_fill_list.size();
    new_emp_ind=a%state._up_fill_list.size();
    old_emp_ind=new_up_fill_ind;
    old_up_fill_ind=new_emp_ind;
    new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
    new_state_vec(state._emp_list(old_emp_ind))=2;
    swt=UP_EMP;
    
    new_amp=quantum_amplitude(psi,new_state_vec);
  }
  
  else if(a<wayofnewupswaps+wayofnewdownswaps)
  {
    a-=(wayofnewupswaps);
    new_down_fill_ind=a/state._down_fill_list.size();
    new_emp_ind=a%state._down_fill_list.size();
    old_emp_ind=new_down_fill_ind;
    old_down_fill_ind=new_emp_ind;
    new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
    new_state_vec(state._emp_list(old_emp_ind))=3;
    swt=DOWN_EMP;
    
    new_amp=quantum_amplitude(psi,new_state_vec);
  }
  else
  {
    a-=wayofnewupswaps+wayofnewdownswaps;
    new_down_fill_ind=a/state._down_fill_list.size();
    new_up_fill_ind=a%state._down_fill_list.size();
    old_up_fill_ind=new_down_fill_ind;
    old_down_fill_ind=new_up_fill_ind;
    new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
    new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
    swt=UP_DOWN;

    new_amp=quantum_amplitude(psi,new_state_vec);
  }
  
  dP probratio=pow(std::abs(new_amp/state._amp),2.);
  if(acceptq(probratio,_rengine))
  {
    state._state_vec=new_state_vec;
    generate_state_lists_using_state_vec_tj(state);
    state._amp=new_amp;
    return swt;
  }
  else
  {
    return NONMOVEtJ;
  }
}

SwapTypeHubbard try_one_step_hubbard_symmetrized(const IQMPS & psi, state_storage_hubbard & state, trng::yarn3 & _rengine,const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  SwapTypeHubbard swt=NONMOVE;
  
  int wayofnewupswaps=state._up_fill_list.size()*state._up_emp_list.size();
  int wayofnewdownswaps=state._down_fill_list.size()*state._down_emp_list.size();
  
  int totwayofswaps=wayofnewupswaps+wayofnewdownswaps;
  
  trng::uniform_int_dist U(0, totwayofswaps);
  int a=U(_rengine);

  int new_up_fill_ind, new_down_fill_ind, new_up_emp_ind, new_down_emp_ind, old_up_fill_ind,old_down_fill_ind,old_up_emp_ind,old_down_emp_ind;
  
  const ITiVec & old_state_vec=state._state_vec;
  ITiVec new_state_vec=old_state_vec;
  
  dComplex new_amp=0.;
  
  if(a<wayofnewupswaps)
  {
    new_up_fill_ind=a/state._up_fill_list.size();
    new_up_emp_ind=a%state._up_fill_list.size();
    old_up_emp_ind=new_up_fill_ind;
    old_up_fill_ind=new_up_emp_ind;
    
    if(old_state_vec(state._up_fill_list(old_up_fill_ind))==2&&old_state_vec(state._up_emp_list(old_up_emp_ind))==1)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=2;
      swt=UP_EMP_to_EMP_UP;
    }
    else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==2&&old_state_vec(state._up_emp_list(old_up_emp_ind))==3)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=4;
      swt=UP_DOWN_to_EMP_DL;
    }
     else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==4&&old_state_vec(state._up_emp_list(old_up_emp_ind))==1)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=2;
      swt=DL_EMP_to_DOWN_UP;
    }
     else if(old_state_vec(state._up_fill_list(old_up_fill_ind))==4&&old_state_vec(state._up_emp_list(old_up_emp_ind))==3)
    {
      new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
      new_state_vec(state._up_emp_list(old_up_emp_ind))=4;
      swt=DL_DOWN_to_DOWN_DL;
    }
    else
    {
      cout<<"Error: WRONG "<<swt<< " Swap in try_one_step_hubbard."<<endl;
      exit(0);
    }
    new_amp=symmetrized_quantum_amplitude_hubbard(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  }
  
  else
  {
    a-=(wayofnewupswaps);
    new_down_fill_ind=a/state._down_fill_list.size();
    new_down_emp_ind=a%state._down_fill_list.size();
    old_down_emp_ind=new_down_fill_ind;
    old_down_fill_ind=new_down_emp_ind;
    
     if(old_state_vec(state._down_fill_list(old_down_fill_ind))==3&&old_state_vec(state._down_emp_list(old_down_emp_ind))==1)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=3;
      swt=DOWN_EMP_to_EMP_DOWN;
    }
    else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==3&&old_state_vec(state._down_emp_list(old_down_emp_ind))==2)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=4;
      swt=DOWN_UP_to_EMP_DL;
    }
     else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==4&&old_state_vec(state._down_emp_list(old_down_emp_ind))==1)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=3;
      swt=DL_EMP_to_UP_DOWN;
    }
     else if(old_state_vec(state._down_fill_list(old_down_fill_ind))==4&&old_state_vec(state._down_emp_list(old_down_emp_ind))==2)
    {
      new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
      new_state_vec(state._down_emp_list(old_down_emp_ind))=4;
      swt=DL_UP_to_UP_DL;
    }
    else
    {
      cout<<"Error: WRONG "<<swt<< " Swap in try_one_step_hubbard."<<endl;
      exit(0);
    }
    new_amp=symmetrized_quantum_amplitude_hubbard(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  }
  
  dP probratio=pow(std::abs(new_amp/state._amp),2.);
  if(acceptq(probratio,_rengine))
  {
    state._state_vec=new_state_vec;
    generate_state_lists_using_state_vec_hubbard(state);
    state._amp=new_amp;
    return swt;
  }
  else
  {
    return NONMOVE;
  }
}

SwapTypetJ try_one_step_tj_symmetrized(const IQMPS & psi, state_storage_tj & state, trng::yarn3 & _rengine,const ITiVecArray & sym_sector_group_site_mapping_table, const ITcVec & sym_sector_group_eigval_table)
{
  SwapTypetJ swt=NONMOVEtJ;
  
  int wayofnewupswaps=state._up_fill_list.size()*state._emp_list.size();
  int wayofnewdownswaps=state._down_fill_list.size()*state._emp_list.size();
  int waysofupdownswaps=state._up_fill_list.size()*state._down_fill_list.size();
  int totwayofswaps=wayofnewupswaps+wayofnewdownswaps+waysofupdownswaps;
  
  trng::uniform_int_dist U(0, totwayofswaps);
  int a=U(_rengine);

  int new_up_fill_ind, new_down_fill_ind, new_emp_ind, old_up_fill_ind,old_down_fill_ind,old_emp_ind;
  
  const ITiVec & old_state_vec=state._state_vec;
  ITiVec new_state_vec=old_state_vec;
  
  dComplex new_amp=0.;
  
  if(a<wayofnewupswaps)
  {
    new_up_fill_ind=a/state._up_fill_list.size();
    new_emp_ind=a%state._up_fill_list.size();
    old_emp_ind=new_up_fill_ind;
    old_up_fill_ind=new_emp_ind;
    new_state_vec(state._up_fill_list(old_up_fill_ind))=1;
    new_state_vec(state._emp_list(old_emp_ind))=2;
    swt=UP_EMP;
    
    new_amp=symmetrized_quantum_amplitude_tj(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  }
  
  else if(a<wayofnewupswaps+wayofnewdownswaps)
  {
    a-=(wayofnewupswaps);
    new_down_fill_ind=a/state._down_fill_list.size();
    new_emp_ind=a%state._down_fill_list.size();
    old_emp_ind=new_down_fill_ind;
    old_down_fill_ind=new_emp_ind;
    new_state_vec(state._down_fill_list(old_down_fill_ind))=1;
    new_state_vec(state._emp_list(old_emp_ind))=3;
    swt=DOWN_EMP;
    
    new_amp=symmetrized_quantum_amplitude_tj(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  }
  else
  {
    a-=wayofnewupswaps+wayofnewdownswaps;
    new_down_fill_ind=a/state._down_fill_list.size();
    new_up_fill_ind=a%state._down_fill_list.size();
    old_up_fill_ind=new_down_fill_ind;
    old_down_fill_ind=new_up_fill_ind;
    new_state_vec(state._down_fill_list(old_down_fill_ind))=2;
    new_state_vec(state._up_fill_list(old_up_fill_ind))=3;
    swt=UP_DOWN;

    new_amp=symmetrized_quantum_amplitude_tj(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  }
  
  dP probratio=pow(std::abs(new_amp/state._amp),2.);
  if(acceptq(probratio,_rengine))
  {
    state._state_vec=new_state_vec;
    generate_state_lists_using_state_vec_tj(state);
    state._amp=new_amp;
    return swt;
  }
  else
  {
    return NONMOVEtJ;
  }
}

dComplex symmetry_eigval_hubbard_MonteCarlo(std::string desc,const IQMPS & psi, const ITiVec & sym_site_mapping_table)
{
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);

	int Nsite=psi.N();

  state_storage_hubbard cur_state=initialize_random_state_vec_hubbard_rengine(psi,Nferm,_rengine);
  
  ITcVec eigval_list(Monte_Carlo_nmeasure);
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
	/////////////////////////////////////
	std::stringstream fim;
	fim << "dataEigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	/////////////////////////////////////
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1w;
	fim1w << "statesUNSYM_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1w = fim1w.str ();
	std::ofstream stout;
	stout.open(ime1w.c_str());
	/////////////////////////////////////
#endif

  for(int i=0;i<Monte_Carlo_nmeasure;i++)
  {
	  /////////////////////////////////////
	  if (rank==0) {
		  if(i%10==0)cout<<"i:"<<i<<endl;
		  cout.flush();
	  }
	  /////////////////////////////////////
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
      
     //do one measurement
     {
     ITiVec state_vec=cur_state._state_vec;
     dComplex old_amp=quantum_amplitude(psi, state_vec);
     if(std::abs((cur_state._amp-old_amp)/old_amp)>1.E-10)
     {
       cout<<"ERROR: amp wrong in symmetry_eigval_hubbard_MonteCarlo"<<endl;
       exit(0);
     }
     finout<<"From rank="<<rank<<"old amp:"<<old_amp<<endl;
  
     ITiVec new_state_vec(state_vec.size());
     for(int i=0;i<new_state_vec.size();i++)
     {
        new_state_vec(i)=state_vec(sym_site_mapping_table(i));
     }
  
     dComplex new_amp=quantum_amplitude(psi,new_state_vec);
     finout<<"From rank="<<rank<<"new amp:"<<new_amp<<endl;
#if WRITE_STATES==1
		 for (int u=0; u<Nsite; u++) { stout<<state_vec[u]<<endl; }
#endif
     dComplex cur_eigval=sym_mapping_fermion_sign_hubbard (state_vec, sym_site_mapping_table)*new_amp/old_amp;
     eigval_list(measure_count)=cur_eigval;
     measure_count++;
         
     }
  }
#if WRITE_STATES==1
	stout.close();
#endif
	finout.close();
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "Eigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime1 = fim1.str();
	std::ofstream finout1;
	finout1.open(ime1.c_str());
	/////////////////////////////////////
	
	finout1<<itpp::mean(eigval_list)<<endl<<std::sqrt(itpp::variance(eigval_list))<<endl;
	finout1<<eigval_list<<endl;
	//  finout1<<"From rank="<<rank<<": eigval_mean="<<itpp::mean(eigval_list)<<" stderr="<<std::sqrt(itpp::variance(eigval_list)/eigval_list.size())<<endl;
	finout1.close();
  return itpp::mean(eigval_list);
}

dComplex symmetry_eigval_tj_MonteCarlo(std::string desc,const IQMPS & psi, const ITiVec & sym_site_mapping_table)
{
    trng::yarn3 _rengine;
    _rengine.seed(0);
    itpp::RNG_randomize();
    _rengine.seed(itpp::randi(0 , 3200));
    //initialize parallel random engine using leapfrog method
    int size=MPI::COMM_WORLD.Get_size(); // get total number of processes
    int rank=MPI::COMM_WORLD.Get_rank(); // get rank of current process
    _rengine.split(size, rank);

	int Nsite=psi.N();

  state_storage_tj cur_state=initialize_random_state_vec_tj_rengine(psi,Nferm,_rengine);
  
  ITcVec eigval_list(Monte_Carlo_nmeasure);
  int measure_count=0;
  
  for(int i=0;i<Monte_Carlo_ntherm;i++)
  {
	  /////////////////////////////////////
	  if (rank==0) {
		  if(i%10==0)cout<<"i:"<<i<<endl;
		  cout.flush();
	  }
	  /////////////////////////////////////
    try_one_step_tj(psi,cur_state, _rengine);
  }
  	/////////////////////////////////////
	std::stringstream fim;
	fim << "dataEigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	/////////////////////////////////////
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1w;
	fim1w << "statesUNSYM_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1w = fim1w.str ();
	std::ofstream stout;
	stout.open(ime1w.c_str());
	/////////////////////////////////////
#endif

  for(int i=0;i<Monte_Carlo_nmeasure;i++)
  {
	  /////////////////////////////////////
	  if (rank==0) {
		  if(i%10==0)cout<<"i:"<<i<<endl;
		  cout.flush();
	  }
	  /////////////////////////////////////
     for(int j=0;j<Monte_Carlo_n_between_measure;j++)
     {
		 /////////////////////////////////////
		 if (rank==0) {
			 if(j%10==0)cout<<"j:"<<j<<endl;
			 cout.flush();
		 }
		 /////////////////////////////////////
        try_one_step_tj(psi,cur_state, _rengine);
     }
     
     //do one measurement
     {
     ITiVec state_vec=cur_state._state_vec;
     dComplex old_amp=quantum_amplitude(psi, state_vec);
     if(std::abs((cur_state._amp-old_amp)/old_amp)>1.E-10)
     {
       cout<<"ERROR: amp wrong in symmetry_eigval_tj_MonteCarlo"<<endl;
       exit(0);
     }
     finout<<"From rank="<<rank<<"old amp:"<<old_amp<<endl;
  
     ITiVec new_state_vec(state_vec.size());
     for(int i=0;i<new_state_vec.size();i++)
     {
        new_state_vec(i)=state_vec(sym_site_mapping_table(i));
     }
  
     dComplex new_amp=quantum_amplitude(psi,new_state_vec);
     finout<<"From rank="<<rank<<"new amp:"<<new_amp<<endl;
#if WRITE_STATES==1
		 for (int u=0; u<Nsite; u++) { stout<<state_vec[u]<<endl; }
#endif

     dComplex cur_eigval=sym_mapping_fermion_sign_tj (state_vec, sym_site_mapping_table)*new_amp/old_amp;

     eigval_list(measure_count)=cur_eigval;
     measure_count++;
     }
  }
#if WRITE_STATES==1
	stout.close();
#endif
	finout.close();
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "Eigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime1 = fim1.str();
	std::ofstream finout1;
	finout1.open(ime1.c_str());
	/////////////////////////////////////
	
	finout1<<itpp::mean(eigval_list)<<endl<<std::sqrt(itpp::variance(eigval_list))<<endl;
	finout1<<eigval_list<<endl;
	//  finout1<<"From rank="<<rank<<": eigval_mean="<<itpp::mean(eigval_list)<<" stderr="<<std::sqrt(itpp::variance(eigval_list)/eigval_list.size())<<endl;
	finout1.close();

	return itpp::mean(eigval_list);
}

dComplex symmetry_sector_eigval_hubbard_MonteCarlo(std::string desc, const IQMPS & psi, const ITiVec & sym_site_mapping_table, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table)
{
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
    
  state_storage_hubbard cur_state=initialize_random_state_vec_hubbard_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
  ITcVec eigval_list(Monte_Carlo_nmeasure);
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
  
	/////////////////////////////////////
	std::stringstream fim;
	fim << "dataEigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	/////////////////////////////////////
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1w;
	fim1w << "statesSYM1_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1w = fim1w.str ();
	std::ofstream stout;
	stout.open(ime1w.c_str());
	/////////////////////////////////////
#endif

  for(int i=0;i<Monte_Carlo_nmeasure;i++)
  {
	  /////////////////////////////////////
	  if (rank==0) {
		  if(i%10==0)cout<<"i:"<<i<<endl;
		  cout.flush();
	  }
	  /////////////////////////////////////
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
     
     //do one measurement
     {
     ITiVec state_vec=cur_state._state_vec;
     dComplex old_amp=symmetrized_quantum_amplitude_hubbard(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
     if(std::abs((cur_state._amp-old_amp)/old_amp)>1.E-10)
     {
       cout<<"ERROR: amp wrong in symmetry_eigval_hubbard_MonteCarlo"<<endl;
       exit(0);
     }
     finout<<"From rank="<<rank<<"old amp:"<<old_amp<<endl;
  
     ITiVec new_state_vec(state_vec.size());
     for(int i=0;i<new_state_vec.size();i++)
     {
        new_state_vec(i)=state_vec(sym_site_mapping_table(i));
     }
  
     dComplex new_amp=symmetrized_quantum_amplitude_hubbard(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
     finout<<"From rank="<<rank<<"new amp:"<<new_amp<<endl;
#if WRITE_STATES==1
		 for (int u=0; u<Nsite; u++) { stout<<state_vec[u]<<endl; }
#endif

     dComplex cur_eigval=sym_mapping_fermion_sign_hubbard (state_vec, sym_site_mapping_table)*new_amp/old_amp;

     eigval_list(measure_count)=cur_eigval;
     measure_count++;
     }
  }
#if WRITE_STATES==1
	stout.close();
#endif
	finout.close();
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "Eigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime1 = fim1.str();
	std::ofstream finout1;
	finout1.open(ime1.c_str());
	/////////////////////////////////////

	finout1<<itpp::mean(eigval_list)<<endl<<std::sqrt(itpp::variance(eigval_list))<<endl;
	finout1<<eigval_list<<endl;
//  finout1<<"From rank="<<rank<<": eigval_mean="<<itpp::mean(eigval_list)<<" stderr="<<std::sqrt(itpp::variance(eigval_list)/eigval_list.size())<<endl;
	finout1.close();
  return itpp::mean(eigval_list);
  
  
}



dComplex symmetry_sector_eigval_tj_MonteCarlo(std::string desc,const IQMPS & psi, const ITiVec & sym_site_mapping_table, const ITiVecArray & sym_sector_generator_site_mapping_table, const ITcVec & sym_sector_generator_eigval_table)
{
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
  
  cout<<sym_sector_group_site_mapping_table<<endl;
  cout<<sym_sector_group_eigval_table<<endl;
  
  state_storage_tj cur_state=initialize_random_state_vec_tj_rengine_symmetrized(psi,Nferm,_rengine,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
  
  ITcVec eigval_list(Monte_Carlo_nmeasure);
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
  	/////////////////////////////////////
	std::stringstream fim;
	fim << "dataEigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime = fim.str();
	std::ofstream finout;
	finout.open(ime.c_str());
	/////////////////////////////////////
#if WRITE_STATES==1
	/////////////////////////////////////
	std::stringstream fim1w;
	fim1w << "statesSYM1_"<<Nsite<<desc<<"_"<<Monte_Carlo_ntherm<<"_"<<Monte_Carlo_nmeasure<<"x"<<Monte_Carlo_n_between_measure<<"rank"<<rank<<".txt";
	std::string ime1w = fim1w.str ();
	std::ofstream stout;
	stout.open(ime1w.c_str());
	/////////////////////////////////////
#endif

  for(int i=0;i<Monte_Carlo_nmeasure;i++)
  {
	  /////////////////////////////////////
	  if (rank==0) {
		  if(i%10==0)cout<<"i:"<<i<<endl;
		  cout.flush();
	  }
	  /////////////////////////////////////
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
     
     //do one measurement
     {
     ITiVec state_vec=cur_state._state_vec;
     dComplex old_amp=symmetrized_quantum_amplitude_tj(psi, state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
     if(std::abs((cur_state._amp-old_amp)/old_amp)>1.E-10)
     {
       finout<<"ERROR: amp wrong in symmetry_eigval_tj_MonteCarlo"<<endl;
       exit(0);
     }
     finout<<"From rank="<<rank<<"old amp:"<<old_amp<<endl;
  
     ITiVec new_state_vec(state_vec.size());
     for(int i=0;i<new_state_vec.size();i++)
     {
        new_state_vec(i)=state_vec(sym_site_mapping_table(i));
     }
  
     dComplex new_amp=symmetrized_quantum_amplitude_tj(psi,new_state_vec,sym_sector_group_site_mapping_table, sym_sector_group_eigval_table);
     finout<<"From rank="<<rank<<"new amp:"<<new_amp<<endl;
#if WRITE_STATES==1
		 for (int u=0; u<Nsite; u++) { stout<<state_vec[u]<<endl; }
#endif
     
     dComplex cur_eigval=sym_mapping_fermion_sign_tj (state_vec, sym_site_mapping_table)*new_amp/old_amp;

     eigval_list(measure_count)=cur_eigval;
     measure_count++;
     }
  }
#if WRITE_STATES==1
	stout.close();
#endif
  //cout<<eigval_list<<endl;
  //cout<<"From rank="<<rank<<": eigval_mean="<<itpp::mean(eigval_list)<<" stderr="<<std::sqrt(itpp::variance(eigval_list)/eigval_list.size())<<endl;
	finout.close();
	/////////////////////////////////////
	std::stringstream fim1;
	fim1 << "Eigval"<<desc<<"_rank"<<rank<<".txt";
	std::string ime1 = fim1.str();
	std::ofstream finout1;
	finout1.open(ime1.c_str());
	/////////////////////////////////////
	
	finout1<<itpp::mean(eigval_list)<<endl<<std::sqrt(itpp::variance(eigval_list))<<endl;
	finout1<<eigval_list<<endl;
	//  finout1<<"From rank="<<rank<<": eigval_mean="<<itpp::mean(eigval_list)<<" stderr="<<std::sqrt(itpp::variance(eigval_list)/eigval_list.size())<<endl;
	finout1.close();

	return itpp::mean(eigval_list);
  
  
}


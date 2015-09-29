//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#define __ITENSOR_HAMS_EXTENDEDHUBBARD_H
#include "../mpo.h"
#include "../model/hubbard.h"

//LATTICE=0 indicates Honeycomb lattice
//LATTICE=1 indicates Triangular lattice
#define LATTICE 1
//SYSSHAPE=0 indicates system has rhombus shape
//SYSSHAPE=1 indicates system has hexagonal shape
#define SYSSHAPE 1
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the central row to edge rows
#define NUMBERING 1
#define BREAKING 0

class ExtendedHubbard
    {
    public:

    ExtendedHubbard(const Hubbard& model, 
                    const OptSet& opts = Global::opts());

    Real
    t1() const { return t1_; }
    void
    t1(Real val) { initted_ = false; t1_ = val; }

    Real
    t2() const { return t2_; }
    void
    t2(Real val) { initted_ = false; t2_ = val; }

    Real
    U() const { return U_; }
    void
    U(Real val) { initted_ = false; U_ = val; }

    Real
    V1() const { return V1_; }
    void
    V1(Real val) { initted_ = false; V1_ = val; }
        
    Real
    Udif() const { return Udif_; }
    void
    Udif(Real val) { initted_ = false; Udif_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    const Hubbard& model_;
    Real U_,
         t1_,
         t2_,
         V1_,
         Udif_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void init_();

    }; //class HubbardChain

inline ExtendedHubbard::
ExtendedHubbard(const Hubbard& model,
                const OptSet& opts)
    :
    model_(model), 
    initted_(false)
    { 
    U_ = opts.getReal("U",19.5);
    t1_ = opts.getReal("t1",-1.0);
    t2_ = opts.getReal("t2",0);
    V1_ = opts.getReal("V1",0);
    Udif_ = opts.getReal("Udif",0.0);
    }

void inline ExtendedHubbard::
init_()
    {
    if(initted_) return;
        
    cout << "U is " << U_ << endl;
    cout << "t1 is " << t1_ << endl;
    cout << "t2 is " << t2_ << endl;
    cout << "V1 is " << V1_ << endl;
    cout << "Udif is " << Udif_ << endl;

    H = MPO(model_);

    const int Ns = model_.N();
		//dimension of MPO matrices; equals 2+4*nn, with nn longest range hopping
//			const int k = 78;//honeycomb hexagonal 24 sites
			//const int k = 94;//4x4x2
			//const int k = 22;//2x2x2
//			const int k = 7 + (t2_ == 0 ? 0 : 4);//other lattices NN and NNN
        //else if(LATTICE==1 && SYSSHAPE==1 && Ns==12){
        const int k = 42;//triangular lattice hexagon shape with 12 sites
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==16){
        //const int k = 62;//triangular lattice rhombus shape with 16 sites
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==24){
        //const int k = 46;//triangular lattice rhombus shape with 24 sites
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==4){
        //const int k = 14;//triangular lattice rhombus shape with 4 sites
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==0){
        //const int k = 142;//triangular lattice rhombus shape with 36 sites, numbering 0
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==1){
        //const int k = 70;//triangular lattice rhombus shape with 36 sites, numbering 1
        //}
        //else if(LATTICE==1 && SYSSHAPE==0 && Ns==36){
        //const int k = 70;//triangular lattice rhombus shape with 36 sites
        //}

		std::vector<Index> links(Ns+1);
		for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);
		
		ITensor W;
		
        
		if (LATTICE==0 && SYSSHAPE==0 && Ns==8) {//honeycomb 2x2x2
			Real tt1,tt3,tt5;
			for(int n = 1; n <= Ns; ++n)
			{
				if (n==1 || n==8) {
					tt1=t1_;tt3=t1_;tt5=t1_;
				} else if (n==2 || n==5){
					tt1=t1_;tt3=t1_;tt5=0;
				} else if (n==3){
					tt1=t1_;tt3=0;tt5=t1_;
				} else if (n==4){
					tt1=0;tt3=t1_;tt5=0;
				} else if (n==6 || n==7){
					tt1=t1_;tt3=0;tt5=0;
				}
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;

				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
				W += model_.Cdagdn(n)*row(7)*col(1)*(-1.0);
				W += model_.Cup(n)*row(12)*col(1)*(-1.0);
				W += model_.Cdn(n)*row(17)*col(1)*(-1.0);

				for (int j=0; j<4; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+5+j)*col(2+5+j);
					W += model_.fermiPhase(n)*row(3+10+j)*col(2+10+j);
					W += model_.fermiPhase(n)*row(3+15+j)*col(2+15+j);
				}

				if (tt1>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * tt1;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(7) * tt1;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(12) * tt1;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(17) * tt1;
				}
				if (tt3>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(4) * tt3;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(9) * tt3;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(14) * tt3;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(19) * tt3;
				}
				if (tt5>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(6) * tt5;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(11) * tt5;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(16) * tt5;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(21) * tt5;
				}
			}
		} else if (LATTICE==0 && SYSSHAPE==0 && Ns==32) {//honeycomb 4x4x2
				for(int n = 1; n <= Ns; ++n) {
					std::vector<Real> ts(23,0);
					int nn=23;//longest range hopping
					
					if (n==8) {
						ts[22]=t1_;
					} else if (n==25){
						ts[0]=t1_;ts[6]=t1_;
					} else if (n==2 ||n==4 ||n==6){
						ts[0]=t1_;ts[22]=t1_;
					} else if (n==1||n==9||n==17){
						ts[0]=t1_;ts[6]=t1_;ts[8]=t1_;
					} else if (n==3 || n==5|| n==7|| n==11|| n==13|| n==15|| n==19|| n==21|| n==23){
						ts[0]=t1_;ts[8]=t1_;
					} else if (n!=16 && n!=24 && n!=32) {
						ts[0]=t1_;
					}
					//for(int ert=0;ert<nn;ert++){ std::cout<<ts[ert]; }
					//std::cout<<std::endl;
					
					ITensor& W = H.Anc(n);
					Index &row = links[n-1], &col = links[n];
					
					W = ITensor(model_.si(n),model_.siP(n),row,col);
					
					//Identity strings
					W += model_.id(n) * row(1) * col(1);
					W += model_.id(n) * row(k) * col(k);
					
					//Hubbard U term
					W += model_.Nupdn(n) * row(k) * col(1) * U_;
					
					W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
					W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);
					W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);
					W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);
					
					for (int j=0; j<nn-1; j++) {
						W += model_.fermiPhase(n)*row(3+j)*col(2+j);
						W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
						W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
						W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					}
					for (int j=1; j<=nn; j++) {
						if (ts[j-1]>0.00001) {
							W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
							W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
							W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
							W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
						}
					}
			}
		} else if(LATTICE==0 && SYSSHAPE==1 && Ns==24){//honeycomb hexagonal 24 sites
			for(int n = 1; n <= Ns; ++n) {
				std::vector<Real> ts(19,0);
				int nn=19;//longest range hopping
				
				if (n==1) {
					ts[0]=t1_;ts[5]=t1_;ts[17]=t1_;
				} else if (n==2 || n==4){
					ts[0]=t1_;ts[18]=t1_;
				} else if (n==3 ||n==14 ||n==16 ||n==18){
					ts[0]=t1_;ts[5]=t1_;
				} else if (n==5){
					ts[5]=t1_;ts[7]=t1_;
				} else if (n==6){
					ts[0]=t1_;ts[6]=t1_;ts[17]=t1_;
				} else if (n==7|| n==9||n==11|| n==13|| n==15|| n==17|| n==20||n==21||n==22|| n==23){
					ts[0]=t1_;
				} else if (n==8 || n==10) {
					ts[0]=t1_;ts[6]=t1_;
				} else if (n==12) {
					ts[6]=t1_;ts[7]=t1_;
				}
				//for(int ert=0;ert<nn;ert++){ std::cout<<ts[ert]; }
				//std::cout<<std::endl;
				
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				
				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
				}
				for (int j=1; j<=nn; j++) {
					if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					}
				}
			}
		}else if(LATTICE==1 && SYSSHAPE==0 && Ns==16){//triangular lattice rhombus shape with 16 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(15,0);
				int nn=15;
                
                if (BREAKING==0) {
                    if (n==1) {
                        ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[4]=t1_;ts[11]=t1_;ts[14]=t1_;   //0,2,3,4,11,14
                    } else if(n==2 || n==3){
                        ts[0]=t1_;ts[3]=t1_;ts[4]=t1_;ts[10]=t1_;ts[11]=t1_;  //0,3,4,10,11
                    } else if(n==4){
                        ts[0]=t1_;ts[3]=t1_;ts[10]=t1_;ts[11]=t1_;  //0,3,10,11
                    } else if(n==5 || n==9){
                        ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[4]=t1_;  //0,2,3,4
                    } else if(n==6 || n==7 || n==10 || n==11){
                        ts[0]=t1_;ts[3]=t1_;ts[4]=t1_;
                    } else if(n==8 || n==12){
                        ts[0]=t1_;ts[3]=t1_;
                    } else if(n==13){
                        ts[0]=t1_;ts[2]=t1_;
                    } else if(n==14 || n==15){
                        ts[0]=t1_;
                    }
                }else if(BREAKING==1){
                    
                    if (n==1) {
                        ts[0]=2.0*t1_;//ts[2]=2.0*t1_;ts[11]=2.0*t1_;ts[14]=2.0*t1_;//2,11,14
                    }else if (n==4) {
                        ts[8]=2.0*t1_;ts[11]=2.0*t1_;//8,11
                    }else if (n==13) {
                        ts[2]=2.0*t1_;//2
                    }else if (n==2) {
                        ts[0]=2.0*t1_;ts[11]=2.0*t1_;ts[12]=2.0*t1_;//0,11,12
                    }else if (n==3) {
                        ts[10]=2.0*t1_;ts[11]=2.0*t1_;//10,11
                    }else if (n==14) {
                        ts[0]=2.0*t1_;//0
                    }else if (n==5) {
                        ts[2]=2.0*t1_;ts[3]=2.0*t1_;ts[6]=2.0*t1_;//2,3,6
                    }else if (n==8) {
                        ts[0]=2.0*t1_;ts[3]=2.0*t1_;//0,3
                    }else if (n==9) {
                        ts[2]=2.0*t1_;//2
                    }else if (n==6) {
                        ts[0]=2.0*t1_;ts[3]=2.0*t1_;ts[4]=2.0*t1_;//0,3,4
                    }else if (n==7) {
                        ts[2]=2.0*t1_;ts[3]=2.0*t1_;//2,3
                    }else if (n==10) {
                        ts[0]=2.0*t1_;//0
                    }
                     
                }
                
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				//Hubbard U term
				//W += model_.Nupdn(n) * row(k) * col(1) * U_;
				/////////////////////////////////////////////////W += model_.Nupdn(n) * row(2+4*nn) * col(1) * U_;//////////
				W += (model_.Nupdn(n)) * row(2+4*nn) * col(1) * U_;
                W += (model_.Ndif(n)) * row(2+4*nn) * col(1) * Udif_;
                W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);//
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);//
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);//
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);//
                
                //Hubbard V1 term
                W += model_.op("Ntot",n) * row(k-1) * col(1);
                W += model_.op("Ntot",n) * row(k) * col(k-1) * V1_;
                
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
				}
			}
        
        }else if(LATTICE==1 && SYSSHAPE==0 && Ns==24){//triangular lattice rhombus shape with 16 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(11,0);
				int nn=11;
                
                
                    if (n==1) {
                        ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[4]=t1_;ts[7]=t1_;ts[10]=t1_;   //0,2,3,4,7,10
                    } else if(n==2 || n==3){
                        ts[0]=t1_;ts[3]=t1_;ts[4]=t1_;ts[6]=t1_;ts[7]=t1_;  //0,3,4,6,7
                    } else if(n==4){
                        ts[0]=t1_;ts[3]=t1_;ts[6]=t1_;ts[7]=t1_;  //0,3,6,7
                    } else if(n==5 || n==13){
                        ts[0]=t1_;ts[2]=t1_;ts[7]=t1_;ts[8]=t1_;  //0,2,7,8
                    } else if(n==6 || n==7 || n==14 || n==15){
                        ts[0]=t1_;ts[7]=t1_;ts[8]=t1_;//0,7,8
                    } else if(n==8 || n==16){
                        ts[4]=t1_;ts[7]=t1_;//4,7
                    } else if(n==9){
                        ts[0]=t1_;ts[2]=t1_;ts[7]=t1_;ts[10]=t1_;//0,2,7,10
                    } else if(n==10 || n==11){
                        ts[0]=t1_;ts[6]=t1_;ts[7]=t1_;//0,6,7
                    } else if(n==12){
                        ts[6]=t1_;ts[7]=t1_;//6,7
                    } else if(n==17){
                        ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[6]=t1_;//0,2,3,6
                    } else if(n==18 || n==19){
                        ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;//0,2,3
                    } else if(n==20){
                        ts[2]=t1_;ts[3]=t1_;//2,3
                    } else if(n==21){
                        ts[0]=t1_;ts[2]=t1_;//0,2
                    } else if(n==22 || n==23){
                        ts[0]=t1_;//0
                    }
                
                
                
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				//Hubbard U term
				//W += model_.Nupdn(n) * row(k) * col(1) * U_;
				/////////////////////////////////////////////////W += model_.Nupdn(n) * row(2+4*nn) * col(1) * U_;//////////
				W += (model_.Nupdn(n)) * row(2+4*nn) * col(1) * U_;
                W += (model_.Ndif(n)) * row(2+4*nn) * col(1) * Udif_;
                W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);//
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);//
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);//
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);//
                
                //Hubbard V1 term
                W += model_.op("Ntot",n) * row(k-1) * col(1);
                W += model_.op("Ntot",n) * row(k) * col(k-1) * V1_;
                
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
				}
			}
            
        }
        /*
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==36){//triangular lattice rhombus shape with 36 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(17,0);
				int nn=17;
                
                
                if (n==1) {
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[6]=t1_;ts[11]=t1_;ts[16]=t1_;   //0,4,5,6,11,16
                } else if(n>=2 && n<=5){
                    ts[0]=t1_;ts[5]=t1_;ts[6]=t1_;ts[10]=t1_;ts[11]=t1_;  //0,5,6,10,11
                } else if(n==6){
                    ts[0]=t1_;ts[5]=t1_;ts[10]=t1_;ts[11]=t1_;  //0,5,10,11
                } else if(n==7 || n==19){
                    ts[0]=t1_;ts[4]=t1_;ts[11]=t1_;ts[12]=t1_;  //0,4,11,12
                } else if((n>=8 && n<=11) || (n>=20 && n<=23)){
                    ts[0]=t1_;ts[11]=t1_;ts[12]=t1_;//0,11,12
                } else if(n==12){
                    ts[6]=t1_;ts[11]=t1_;//6,11
                } else if(n==13){
                    ts[0]=t1_;ts[4]=t1_;ts[11]=t1_;ts[16]=t1_;//0,4,11,16
                } else if(n>=14 && n<=17){
                    ts[0]=t1_;ts[10]=t1_;ts[11]=t1_;//0,10,11
                }else if(n==18){
                    ts[10]=t1_;ts[11]=t1_;//10,11
                }else if(n==24){
                    ts[6]=t1_;ts[11]=t1_;//6,11
                }else if(n==25){
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[10]=t1_;//0,4,5,10
                }else if(n>=26 && n<=29){
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;//0,4,5
                }else if(n==30){
                    ts[4]=t1_;ts[5]=t1_;//4,5
                }else if(n==31){
                    ts[0]=t1_;ts[4]=t1_;//0,4
                }else if(n>=32 && n<=35){
                    ts[0]=t1_;//0
                }
                
                
                
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				//Hubbard U term
				//W += model_.Nupdn(n) * row(k) * col(1) * U_;
				/////////////////////////////////////////////////W += model_.Nupdn(n) * row(2+4*nn) * col(1) * U_;//////////
				W += (model_.Nupdn(n)) * row(2+4*nn) * col(1) * U_;
                W += (model_.Ndif(n)) * row(2+4*nn) * col(1) * Udif_;
                W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);//
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);//
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);//
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);//
                
                //Hubbard V1 term
                W += model_.op("Ntot",n) * row(k-1) * col(1);
                W += model_.op("Ntot",n) * row(k) * col(k-1) * V1_;
                
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
				}
			}
            
        }
        */
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==4){//triangular lattice rhombus shape with 4 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(3,0);
				int nn=3;
                
                if (n==1) {
                    ts[0]=2.0*t1_;ts[1]=2.0*t1_;ts[2]=2.0*t1_;   //0,2,3,4,11,14
                } else if(n==2){
                    ts[0]=2.0*t1_;ts[1]=2.0*t1_; //0,3,4,10,11
                } else if(n==3){
                    ts[0]=2.0*t1_;  //0,3,10,11
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				//W += model_.Nupdn(n) * row(2+4*nn) * col(1) * U_;//////////
				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);//
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);//
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);//
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);//
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
				}
			}
        }else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==0){//triangular lattice rhombus shape with 36 sites, numbering 0/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(35,0);
				int nn=35;
                
                if (n==1) {
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[6]=t1_;ts[29]=t1_;ts[34]=t1_;  //0,4,5,6,29,34
                } else if(n>=2 && n<=5){
                    ts[0]=t1_;ts[5]=t1_;ts[6]=t1_;ts[28]=t1_;ts[29]=t1_;  //0,5,6,28,29
                } else if(n==6){
                    ts[0]=t1_;ts[5]=t1_;ts[28]=t1_;ts[29]=t1_;  //0,5,28,29
                } else if(n==7 || n==13 || n==19 || n==25){
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[6]=t1_;  //0,4,5,6
                } else if((n>=8 && n<=11) || (n>=14 && n<=17) || (n>=20 && n<=23) || (n>=26 && n<=29)){
                    ts[0]=t1_;ts[5]=t1_;ts[6]=t1_;//0,5,6
                } else if(n==12 || n==18 || n==24 || n==30){
                    ts[0]=t1_;ts[5]=t1_;//0,5
                } else if(n==31){
                    ts[0]=t1_;ts[4]=t1_;//0,4
                } else if(n>=32 && n<=35){
                    ts[0]=t1_;
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				
				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
				}
			}
        }else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==1){//triangular lattice rhombus shape with 36 sites, numbering 1/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(17,0);
				int nn=17;
                
                if (n==1) {
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[6]=t1_;ts[11]=t1_;ts[16]=t1_;   //0,4,5,6,11,16
                } else if(n>=2 && n<=5){
                    ts[0]=t1_;ts[5]=t1_;ts[6]=t1_;ts[10]=t1_;ts[11]=t1_; //0,5,6,10,11
                } else if(n==6){
                    ts[0]=t1_;ts[5]=t1_;ts[10]=t1_;ts[11]=t1_;  //0,5,10,11
                } else if(n==7 || n==19){
                    ts[0]=t1_;ts[4]=t1_;ts[11]=t1_;ts[12]=t1_;  //0,4,11,12
                } else if((n>=8 && n<=11) || (n>=20 && n<=23)){
                    ts[0]=t1_;ts[11]=t1_;ts[12]=t1_;//0,11,12
                } else if(n==12){
                    ts[6]=t1_;ts[11]=t1_;//6,11
                } else if(n==13){
                    ts[0]=t1_;ts[4]=t1_;ts[11]=t1_;ts[16]=t1_;//0,4,11,16
                } else if(n>=14 && n<=17){
                    ts[0]=t1_;ts[10]=t1_;ts[11]=t1_;//0,10,11
                }else if(n==18){
                    ts[10]=t1_;ts[11]=t1_;//10,11
                }else if(n==24){
                    ts[6]=t1_;ts[11]=t1_;//6,11
                }else if(n==25){
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;ts[10]=t1_;//0,4,5,10
                }else if(n>=26 && n<=29){
                    ts[0]=t1_;ts[4]=t1_;ts[5]=t1_;//0,4,5
                }else if(n==30){
                    ts[4]=t1_;ts[5]=t1_;//4,5
                }else if(n==31){
                    ts[0]=t1_;ts[4]=t1_;//0,4
                }else if(n>=32 && n<=35){
                    ts[0]=t1_;//0
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				
				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
				}
			}
        }else if(LATTICE==1 && SYSSHAPE==1 && Ns==12){//triangular lattice hexagon shape with 12 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(10,0);
				int nn=10;
                
                if (n==1) {
                    ts[0]=t1_;ts[1]=t1_;ts[2]=t1_;ts[7]=t1_;ts[8]=t1_;ts[9]=t1_;   //0,1,2,7,8,9
                } else if(n==2){
                    ts[1]=t1_;ts[2]=t1_;ts[3]=t1_;ts[8]=t1_;ts[9]=t1_; //1,2,3,8,9
                } else if(n==3){
                    ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[5]=t1_;ts[8]=t1_;  //0,2,3,5,8
                } else if(n==4 || n==7 || n==8){
                    ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;  //0,2,3
                } else if(n==5){
                    ts[0]=t1_;ts[2]=t1_;ts[3]=t1_;ts[4]=t1_;  //0,2,3,4
                } else if(n==6){
                    ts[0]=t1_;ts[3]=t1_;ts[5]=t1_;  //0,3,5
                } else if(n==9){
                    ts[0]=t1_;ts[2]=t1_;  //0,2
                } else if(n==10 || n==11){
                    ts[0]=t1_; //0
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				//W += model_.Nupdn(n) * row(2+4*nn) * col(1) * U_;//////////
				W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);//
				W += model_.Cdagdn(n)*row(2+nn)*col(1)*(-1.0);//
				W += model_.Cup(n)*row(2+2*nn)*col(1)*(-1.0);//
				W += model_.Cdn(n)*row(2+3*nn)*col(1)*(-1.0);//
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
				}
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
				}
			}
        }else {
			for(int n = 1; n <= Ns; ++n)
			{
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				//Identity strings
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				//Hubbard U term
				W += model_.Nupdn(n) * row(k) * col(1) * U_;
				
				//Hubbard V1 term
				W += model_.Ntot(n) * row(k-1) * col(1);
				W += model_.Ntot(n) * row(k) * col(k-1) * V1_;
				
				if(t2_ == 0)
				{
					//Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t1_;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(3) * t1_;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(4) * t1_;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(5) * t1_;
					
					W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
					W += model_.Cdagdn(n) * row(3) * col(1) * (-1.0);
					W += model_.Cup(n) * row(4) * col(1) * (-1.0);
					W += model_.Cdn(n) * row(5) * col(1) * (-1.0);
				}
				else // t2_ != 0
				{
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t1_;
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(3) * t2_;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(4) * t1_;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(5) * t2_;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(6) * t1_;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(7) * t2_;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(8) * t1_;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(9) * t2_;
					
					W += model_.Cdagup(n)*row(2)*col(1)*(-1.0);
					W += model_.fermiPhase(n)*row(3)*col(2);
					W += model_.Cdagdn(n)*row(4)*col(1)*(-1.0);
					W += model_.fermiPhase(n)*row(5)*col(4);
					W += model_.Cup(n)*row(6)*col(1)*(-1.0);
					W += model_.fermiPhase(n)*row(7)*col(6);
					W += model_.Cdn(n)*row(8)*col(1)*(-1.0);
					W += model_.fermiPhase(n)*row(9)*col(8);
				}
			}
		}

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#endif

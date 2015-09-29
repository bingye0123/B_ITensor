//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_TJCHAIN_H
#define __ITENSOR_HAMS_TJCHAIN_H
#include "../mpo.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//LATTICE=0 indicates Honeycomb lattice
//LATTICE=1 indicates Triangular lattice
#define LATTICE 1
//SYSSHAPE=0 indicates system has rhombus shape
//SYSSHAPE=1 indicates system has hexagonal shape
#define SYSSHAPE 1
//For 36 sites
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the central row to edge rows
//For 16 sites
//NUMBERING=0 indicates numbering from the top row to the bottom row
//NUMBERING=1 indicates numbering from the left col to right col
#define NUMBERING 0

class tJChain 
    {
    public:

    tJChain(const Model& model,
            const OptSet& opts = Global::opts());

    Real
    t() const { return t_; }

    void
    t(Real val) { initted = false; t_ = val; }

    Real
    J() const { return J_; }

    void
    J(Real val) { initted = false; J_ = val; }

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    const Model& model_;
    Real t_,J_;
    bool initted;
    MPO H;

    void init_();

    }; //class tJChain

inline tJChain::
tJChain(const Model& model, const OptSet& opts)
    : 
    model_(model), 
    initted(false)
    { 
    J_ = opts.getReal("J",2.0);
    t_ = opts.getReal("t",-1.0);
    }


void inline tJChain::
init_()
    {
    if(initted) return;
    Cout << "J is " << J_ << Endl;
    Cout << "t is " << t_ << Endl;

    H = MPO(model_);

    const int Ns = model_.N();
    //const int k = 10;//other lattices NN
    //if (LATTICE==0 && SYSSHAPE==0 && Ns==8) {
    //    const int k = 42;//honeycomb lattice rhombus shape with 8 sites 2x2x2
    //}
    //else if(LATTICE==0 && SYSSHAPE==0 && Ns==32){
    //    const int k = 186;//honeycomb lattice rhombus shape with 32 sites 4x4x2
    //}
    //else if(LATTICE==0 && SYSSHAPE==1 && Ns==24){
    //    const int k = 154;//honeycomb lattice hexagonal shape with 24 sites
    //}
	//else if(LATTICE==1 && SYSSHAPE==0 && Ns==16 && NUMBERING==0){
        //const int k = 122;//triangular lattice rhombus shape with 16 sites
    //}
    //else if(LATTICE==1 && SYSSHAPE==0 && Ns==16 && NUMBERING==1){
        //const int k = 106;//triangular lattice rhombus shape with 16 sites
    //}
    //else if(LATTICE==1 && SYSSHAPE==0 && Ns==4){
        //const int k = 26;//triangular lattice rhombus shape with 4 sites
    //}
    //else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==0){
        //const int k = 282;//triangular lattice rhombus shape with 36 sites, numbering 0
    //}
    //else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==1){
        //const int k = 138;//triangular lattice rhombus shape with 36 sites, numbering 1
    //}
    //else if(LATTICE==1 && SYSSHAPE==1 && Ns==28 && NUMBERING==0){
    const int k = 218;//triangular lattice hexagon shape with 28 sites, numbering 0
    //}
		
    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("tjl",l),k);

    ITensor W;
		if (LATTICE==0 && SYSSHAPE==0 && Ns==8) {//honeycomb lattice rhombus shape with 8 sites 2x2x2
			Real tt1,tt3,tt5,J1,J3,J5;
			for(int n = 1; n <= Ns; ++n)
			{
				if (n==1 || n==8) {
					tt1=t_;tt3=t_;tt5=t_;J1=J_;J3=J_;J5=J_;
				} else if (n==2 || n==5){
					tt1=t_;tt3=t_;tt5=0;J1=J_;J3=J_;J5=0;
				} else if (n==3){
					tt1=t_;tt3=0;tt5=t_;J1=J_;J3=0;J5=J_;
				} else if (n==4){
					tt1=0;tt3=t_;tt5=0;J1=0;J3=J_;J5=0;
				} else if (n==6 || n==7){
					tt1=t_;tt3=0;tt5=0;J1=J_;J3=0;J5=0;
				}
				
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);

				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(7) * col(1) * (-1.0);
				W += model_.Cup(n) * row(12) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(17) * col(1) * (-1.0);
				W += model_.sz(n) * row(22) * col(1);
				W += model_.sm(n) * row(27) * col(1);
				W += model_.sp(n) * row(32) * col(1);
				W += model_.Ntot(n) * row(37) * col(1);

				for (int j=0; j<4; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+5+j)*col(2+5+j);
					W += model_.fermiPhase(n)*row(3+10+j)*col(2+10+j);
					W += model_.fermiPhase(n)*row(3+15+j)*col(2+15+j);

					W += model_.id(n)*row(3+20+j)*col(2+20+j);
					W += model_.id(n)*row(3+25+j)*col(2+25+j);
					W += model_.id(n)*row(3+30+j)*col(2+30+j);
					W += model_.id(n)*row(3+35+j)*col(2+35+j);
				}

				//if (tt1>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * tt1;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(7) * tt1;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(12) * tt1;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(17) * tt1;
					W += model_.sz(n) * row(k) * col(22) * J1;
					W += model_.sp(n) * row(k) * col(27) * J1/2;
					W += model_.sm(n) * row(k) * col(32) * J1/2;
					W += model_.Ntot(n) * row(k) * col(37) * (-J1/4);
				//}
				//if (tt3>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(4) * tt3;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(9) * tt3;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(14) * tt3;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(19) * tt3;
					W += model_.sz(n) * row(k) * col(24) * J3;
					W += model_.sp(n) * row(k) * col(29) * J3/2;
					W += model_.sm(n) * row(k) * col(34) * J3/2;
					W += model_.Ntot(n) * row(k) * col(39) * (-J3/4);
				//}
				//if (tt5>0.00001) {
					W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(6) * tt5;
					W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(11) * tt5;
					W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(16) * tt5;
					W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(21) * tt5;
					W += model_.sz(n) * row(k) * col(26) * J5;
					W += model_.sp(n) * row(k) * col(31) * J5/2;
					W += model_.sm(n) * row(k) * col(36) * J5/2;
					W += model_.Ntot(n) * row(k) * col(41) * (-J5/4);
				//}
			}
		} else if(LATTICE==0 && SYSSHAPE==0 && Ns==32) {//honeycomb lattice rhombus shape with 32 sites 4x4x2
			for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(23,0);
				std::vector<Real> Js(23,0);
				int nn=23;

				if (n==8) {
					ts[22]=t_;Js[22]=J_;
				} else if (n==25){
					ts[0]=t_;Js[0]=J_;ts[6]=t_;Js[6]=J_;
				} else if (n==2 ||n==4 ||n==6){
					ts[0]=t_;Js[0]=J_;ts[22]=t_;Js[22]=J_;
				} else if (n==1||n==9||n==17){
					ts[0]=t_;Js[0]=J_;ts[6]=t_;Js[6]=J_;ts[8]=t_;Js[8]=J_;
				} else if (n==3 || n==5|| n==7|| n==11|| n==13|| n==15|| n==19|| n==21|| n==23){
					ts[0]=t_;Js[0]=J_;ts[8]=t_;Js[8]=J_;
				} else if (n!=16 && n!=24 && n!=32) {
					ts[0]=t_;Js[0]=J_;
				}

				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
						W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
						W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
						W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
						W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
		} else if(LATTICE==0 && SYSSHAPE==1 && Ns==24) {//honeycomb lattice hexagonal shape with 24 sites
			for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(19,0);
				std::vector<Real> Js(19,0);
				int nn=19;
				
				if (n==1) {
					ts[0]=t_;ts[5]=t_;ts[17]=t_;Js[0]=J_;Js[5]=J_;Js[17]=J_;
				} else if (n==2 || n==4){
					ts[0]=t_;ts[18]=t_;Js[0]=J_;Js[18]=J_;
				} else if (n==3 ||n==14 ||n==16 ||n==18){
					ts[0]=t_;ts[5]=t_;Js[0]=J_;Js[5]=J_;
				} else if (n==5){
					ts[5]=t_;ts[7]=t_;Js[5]=J_;Js[7]=J_;
				} else if (n==6){
					ts[0]=t_;ts[6]=t_;ts[17]=t_;Js[0]=J_;Js[6]=J_;Js[17]=J_;
				} else if (n==7|| n==9||n==11|| n==13|| n==15|| n==17|| n==20||n==21||n==22|| n==23){
					ts[0]=t_;Js[0]=J_;
				} else if (n==8 || n==10) {
					ts[0]=t_;ts[6]=t_;Js[0]=J_;Js[6]=J_;
				} else if (n==12) {
					ts[6]=t_;ts[7]=t_;Js[6]=J_;Js[7]=J_;
				}
				
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
						W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
						W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
						W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
						W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
		}
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==16 && NUMBERING==0){//triangular lattice rhombus shape with 16 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(15,0);
				std::vector<Real> Js(15,0);
				int nn=15;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;ts[11]=t_;Js[11]=J_;ts[14]=t_;Js[14]=J_;   //0,2,3,4,11,14
                } else if(n==2 || n==3){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;  //0,3,4,10,11
                } else if(n==4){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;  //0,3,10,11
                } else if(n==5 || n==9){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;  //0,2,3,4
                } else if(n==6 || n==7 || n==10 || n==11){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;
                } else if(n==8 || n==12){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;
                } else if(n==13){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;
                } else if(n==14 || n==15){
                    ts[0]=t_;Js[0]=J_;
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);//
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);//
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);//
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);//
				W += model_.sz(n) * row(2+4*nn) * col(1);//
				W += model_.sm(n) * row(2+5*nn) * col(1);//
				W += model_.sp(n) * row(2+6*nn) * col(1);//
				W += model_.Ntot(n) * row(2+7*nn) * col(1);//
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);//
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);//
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);//
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);//
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
					//if (Js[j-1]>0.00001) {
						W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];//
						W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;//
						W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;//
						W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);//
					//}
				}
			}
        }
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==16 && NUMBERING==1){//triangular lattice rhombus shape with 16 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(13,0);
				std::vector<Real> Js(13,0);
				int nn=13;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[6]=t_;Js[6]=J_;ts[11]=t_;Js[11]=J_;ts[12]=t_;Js[12]=J_;   //0,2,3,6,11,12
                } else if(n==2 || n==3){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[11]=t_;Js[11]=J_;ts[12]=t_;Js[12]=J_;  //0,2,3,11,12
                } else if(n==4){
                    ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[8]=t_;Js[8]=J_;ts[11]=t_;Js[11]=J_;  //2,3,8,11
                } else if(n==5 || n==9){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[6]=t_;Js[6]=J_;  //0,2,3,6
                } else if(n==6 || n==7 || n==10 || n==11){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;
                } else if(n==8 || n==12){
                    ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;
                } else if(n==13){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;
                } else if(n==14 || n==15){
                    ts[0]=t_;Js[0]=J_;
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);//
				W += model_.id(n) * row(k) * col(k);//
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);//
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);//
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);//
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);//
				W += model_.sz(n) * row(2+4*nn) * col(1);//
				W += model_.sm(n) * row(2+5*nn) * col(1);//
				W += model_.sp(n) * row(2+6*nn) * col(1);//
				W += model_.Ntot(n) * row(2+7*nn) * col(1);//
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);//
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);//
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);//
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);//
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);//
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);//
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);//
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);//
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];//
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];//
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];//
					//}
					//if (Js[j-1]>0.00001) {
                    W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];//
                    W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;//
                    W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;//
                    W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);//
					//}
				}
			}
        }
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==4){//triangular lattice rhombus shape with 4 sites/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(15,0);
				std::vector<Real> Js(15,0);
				int nn=3;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[1]=t_;Js[1]=J_;ts[2]=t_;Js[2]=J_;   //0,2,3,4,11,14
                } else if(n==2){
                    ts[0]=t_;Js[0]=J_;ts[1]=t_;Js[1]=J_; //0,3,4,10,11
                } else if(n==3){
                    ts[0]=t_;Js[0]=J_;  //0,3,10,11
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
						W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
						W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
						W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
						W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
        }
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==0){//triangular lattice rhombus shape with 36 sites, numbering 0/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(35,0);
				std::vector<Real> Js(35,0);
				int nn=35;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;ts[29]=t_;Js[29]=J_;ts[34]=t_;Js[34]=J_;   //0,4,5,6,29,34
                } else if(n>=2 && n<=5){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;ts[28]=t_;Js[28]=J_;ts[29]=t_;Js[29]=J_;  //0,5,6,28,29
                } else if(n==6){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;ts[28]=t_;Js[28]=J_;ts[29]=t_;Js[29]=J_;  //0,5,28,29
                } else if(n==7 || n==13 || n==19 || n==25){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;  //0,4,5,6
                } else if((n>=8 && n<=11) || (n>=14 && n<=17) || (n>=20 && n<=23) || (n>=26 && n<=29)){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;//0,5,6
                } else if(n==12 || n==18 || n==24 || n==30){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;//0,5
                } else if(n==31){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;//0,4
                } else if(n>=32 && n<=35){
                    ts[0]=t_;Js[0]=J_;
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
						W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
						W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
						W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
						W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
						W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
						W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
						W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
						W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
        }
        else if(LATTICE==1 && SYSSHAPE==1 && Ns==28 && NUMBERING==0){//triangular lattice hexagon shape with 28 sites, numbering 0/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(27,0);
				std::vector<Real> Js(27,0);
				int nn=27;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[23]=t_;Js[23]=J_;ts[24]=t_;Js[24]=J_;ts[26]=t_;Js[26]=J_;   //0,2,3,23,24,26
                } else if(n==2){
                    ts[2]=t_;Js[2]=J_;ts[3]=t_;Js[3]=J_;ts[15]=t_;Js[15]=J_;ts[23]=t_;Js[23]=J_;ts[24]=t_;Js[24]=J_;  //2,3,15,23,24
                } else if(n==3){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[8]=t_;Js[8]=J_;ts[13]=t_;Js[13]=J_;ts[19]=t_;Js[19]=J_;ts[24]=t_;Js[24]=J_;  //0,4,8,13,19,24
                } else if(n==4){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;ts[23]=t_;Js[23]=J_;  //0,3,4,23
                } else if(n==5||n==20||n==21||n==22){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;//0,3,4
                } else if(n==6){
                    ts[0]=t_;Js[0]=J_;ts[1]=t_;Js[1]=J_;ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;//0,1,3,4
                } else if(n==7){
                    ts[3]=t_;Js[3]=J_;ts[4]=t_;Js[4]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;ts[16]=t_;Js[16]=J_;//3,4,10,11,16
                } else if(n==8){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[14]=t_;Js[14]=J_;//0,4,5,14
                } else if((n>=9&&n<=11)||(n>=14&&n<=16)){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;//0,4,5
                } else if(n==12){
                    ts[4]=t_;Js[4]=J_;ts[11]=t_;Js[11]=J_;ts[15]=t_;Js[15]=J_;//4,11,15
                } else if(n==13){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[9]=t_;Js[9]=J_;ts[13]=t_;Js[13]=J_;  //0,4,5,9,13
                } else if(n==17){
                    ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;  //4,5
                } else if(n==18){
                    ts[0]=t_;Js[0]=J_;ts[8]=t_;Js[8]=J_;  //0,8
                } else if(n==19){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;  //0,4
                } else if(n==23){
                    ts[3]=t_;Js[3]=J_;  //3
                } else if(n==24){
                    ts[0]=t_;Js[0]=J_;ts[3]=t_;Js[3]=J_;  //0,3
                } else if(n==25){
                    ts[0]=t_;Js[0]=J_;ts[2]=t_;Js[2]=J_;  //0,2
                } else if(n==26){
                    ts[0]=t_;Js[0]=J_;  //0
                }
                
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
                    W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
                    W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
                    W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
                    W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
        }
        else if(LATTICE==1 && SYSSHAPE==0 && Ns==36 && NUMBERING==1){//triangular lattice rhombus shape with 36 sites, numbering 1/////////////////////
            for(int n = 1; n <= Ns; ++n)
			{
				std::vector<Real> ts(17,0);
				std::vector<Real> Js(17,0);
				int nn=17;
                
                if (n==1) {
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;ts[11]=t_;Js[11]=J_;ts[16]=t_;Js[16]=J_;   //0,4,5,6,11,16
                } else if(n>=2 && n<=5){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;ts[6]=t_;Js[6]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;  //0,5,6,10,11
                } else if(n==6){
                    ts[0]=t_;Js[0]=J_;ts[5]=t_;Js[5]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;  //0,5,10,11
                } else if(n==7 || n==19){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[11]=t_;Js[11]=J_;ts[12]=t_;Js[12]=J_;  //0,4,11,12
                } else if((n>=8 && n<=11) || (n>=20 && n<=23)){
                    ts[0]=t_;Js[0]=J_;ts[11]=t_;Js[11]=J_;ts[12]=t_;Js[12]=J_;//0,11,12
                } else if(n==12){
                    ts[6]=t_;Js[6]=J_;ts[11]=t_;Js[11]=J_;//6,11
                } else if(n==13){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[11]=t_;Js[11]=J_;ts[16]=t_;Js[16]=J_;//0,4,11,16
                } else if(n>=14 && n<=17){
                    ts[0]=t_;Js[0]=J_;ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;//0,10,11
                }else if(n==18){
                    ts[10]=t_;Js[10]=J_;ts[11]=t_;Js[11]=J_;//10,11
                }else if(n==24){
                    ts[6]=t_;Js[6]=J_;ts[11]=t_;Js[11]=J_;//6,11
                }else if(n==25){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;ts[10]=t_;Js[10]=J_;//0,4,5,10
                }else if(n>=26 && n<=29){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;//0,4,5
                }else if(n==30){
                    ts[4]=t_;Js[4]=J_;ts[5]=t_;Js[5]=J_;//4,5
                }else if(n==31){
                    ts[0]=t_;Js[0]=J_;ts[4]=t_;Js[4]=J_;//0,4
                }else if(n>=32 && n<=35){
                    ts[0]=t_;Js[0]=J_;//0
                }
                
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(2+nn) * col(1) * (-1.0);
				W += model_.Cup(n) * row(2+2*nn) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(2+3*nn) * col(1) * (-1.0);
				W += model_.sz(n) * row(2+4*nn) * col(1);
				W += model_.sm(n) * row(2+5*nn) * col(1);
				W += model_.sp(n) * row(2+6*nn) * col(1);
				W += model_.Ntot(n) * row(2+7*nn) * col(1);
				
				for (int j=0; j<nn-1; j++) {
					W += model_.fermiPhase(n)*row(3+j)*col(2+j);
					W += model_.fermiPhase(n)*row(3+nn+j)*col(2+nn+j);
					W += model_.fermiPhase(n)*row(3+2*nn+j)*col(2+2*nn+j);
					W += model_.fermiPhase(n)*row(3+3*nn+j)*col(2+3*nn+j);
					
					W += model_.id(n)*row(3+4*nn+j)*col(2+4*nn+j);
					W += model_.id(n)*row(3+5*nn+j)*col(2+5*nn+j);
					W += model_.id(n)*row(3+6*nn+j)*col(2+6*nn+j);
					W += model_.id(n)*row(3+7*nn+j)*col(2+7*nn+j);
				}
				
				for (int j=1; j<=nn; j++) {
					//if (ts[j-1]>0.00001) {
                    W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(1+j) * ts[j-1];
                    W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(1+j+nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(1+j+2*nn) * ts[j-1];
                    W += multSiteOps(model_.Cdagdn(n),model_.fermiPhase(n)) * row(k) * col(1+j+3*nn) * ts[j-1];
					//}
					//if (Js[j-1]>0.00001) {
                    W += model_.sz(n) * row(k) * col(1+j+4*nn) * Js[j-1];
                    W += model_.sp(n) * row(k) * col(1+j+5*nn) * Js[j-1]/2;
                    W += model_.sm(n) * row(k) * col(1+j+6*nn) * Js[j-1]/2;
                    W += model_.Ntot(n) * row(k) * col(1+j+7*nn) * (-Js[j-1]/4);
					//}
				}
			}
        }
        else {//other lattices
			for(int n = 1; n <= Ns; ++n)
			{
				ITensor& W = H.Anc(n);
				Index &row = links[n-1], &col = links[n];
				
				W = ITensor(model_.si(n),model_.siP(n),row,col);
				
				// fermiPhase will be needed for longer range hopping, but not for this nn chain
				
				//W += model_.Nupdn(n) * row(k) * col(1) * U_;
				//W += multSiteOps(model_.fermiPhase(n),model_.Cup(n)) * row(k) * col(2) * t_;
				//W += multSiteOps(model_.fermiPhase(n),model_.Cdn(n)) * row(k) * col(3) * t_;
				//W += multSiteOps(model_.Cdagup(n),model_.fermiPhase(n)) * row(k) * col(4) * t_;
				//W += multSiteOps(model_.Cdagdn(n),model.fermiPhase(n)) * row(k) * col(5) * t_;
				
				W += model_.Cup(n) * row(k) * col(2) * t_;
				W += model_.Cdn(n) * row(k) * col(3) * t_;
				W += model_.Cdagup(n) * row(k) * col(4) * t_;
				W += model_.Cdagdn(n) * row(k) * col(5) * t_;
				
				W += model_.sz(n) * row(k) * col(6) * J_;
				W += model_.sp(n) * row(k) * col(7) * J_/2;
				W += model_.sm(n) * row(k) * col(8) * J_/2;
				W += model_.Ntot(n) * row(k) * col(9) * (-J_/4);
				W += model_.id(n) * row(k) * col(k);
				
				W += model_.id(n) * row(1) * col(1);
				W += model_.Cdagup(n) * row(2) * col(1) * (-1.0);
				W += model_.Cdagdn(n) * row(3) * col(1) * (-1.0);
				W += model_.Cup(n) * row(4) * col(1) * (-1.0);
				W += model_.Cdn(n) * row(5) * col(1) * (-1.0);
				W += model_.sz(n) * row(6) * col(1);
				W += model_.sm(n) * row(7) * col(1);
				W += model_.sp(n) * row(8) * col(1);
				W += model_.Ntot(n) * row(9) * col(1);
			}
		}
			

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted = true;
    }

#undef Cout
#undef Endl
#undef Format

#endif

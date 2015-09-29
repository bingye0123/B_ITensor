//tensor_utility.h
//include change_index_name2/3/4, group_indice2/4/6/8, bra_otho_center_1, mpo_otho_center_1

#ifndef ___GENERALIZEDEV
#define ___GENERALIZEDEV


#include "matrix.h"
#include "tarray1.h"
#include "minmax.h"
#include <math.h>
#include <fstream>

#include "/home/yeba/B_ITensor/matrix/ev_lapack_wrap.h"

using namespace std;
namespace itensor {
    void
    GeneralizedEigenProblem(const MatrixRef& A, Matrix B, Vector& D, Matrix& Z)
    {
        LAPACK_INT N = A.Ncols();
        if(A.Nrows() < 1)
            _merror("GeneralizedEigenProblem: 0 dimensions matrix");
        if (N != A.Nrows() || A.Nrows() < 1)
            _merror("GeneralizedEigenProblem: Input Matrix must be square");
        
        char jobz = 'V';
        char uplo = 'U';
        LAPACK_INT info;
        
        D.ReDimension(N);
        Z = A;
        
        dsygv_wrapper(&jobz,&uplo,&N,Z.Store(),B.Store(),D.Store(),&info);
        
        if(info != 0)
        {
            cerr << "info is " << info << endl;
            cout << "info is " << info << endl;
            Error("Error in call to dsygv");
            return;
        }
        
        //Transpose Z before return
        Z = Z.t();
    }
    
    void
    GeneralizedNonSymmetricEigenProblem(const Matrix A, Matrix B, Vector& R, Vector& I, Vector& BETA, Matrix& ZL, Matrix& ZR)
    {
        LAPACK_INT N = A.Ncols();
        if(A.Nrows() < 1)
            _merror("GeneralizedNonSymmetricEigenProblem: 0 dimensions matrix");
        if (N != A.Nrows() || A.Nrows() < 1)
            _merror("GeneralizedNonSymmetricEigenProblem: Input Matrix must be square");
        
        char jobvl = 'V';
        char jobvr = 'V';
        LAPACK_INT info;
        
        R.ReDimension(N);
        I.ReDimension(N);
        BETA.ReDimension(N);
        //BETA=R;
        ZL = A;
        ZR = A;
        dggev_wrapper(&jobvl, &jobvr, &N, A.Store(), B.Store(), R.Store(), I.Store(), BETA.Store(), ZL.Store(), ZR.Store(), &info);
        //dsygv_wrapper(&jobz,&uplo,&N,Z.Store(),B.Store(),D.Store(),&info);
        if(info != 0)
        {
            cerr << "info is " << info << endl;
            cout << "info is " << info << endl;
            Error("Error in call to dggev");
            return;
        }
        
        //Transpose Z before return
        ZL = ZL.t();
        ZR = ZR.t();
    }
    
//    void
//    GeneralizedNonSymmetricEigenProblem(const MatrixRef& A, Matrix B, Vector& R, Vector& I, Vector BETA, Matrix& Z)
//    {
//        LAPACK_INT N = A.Ncols();
//        if(A.Nrows() < 1)
//            _merror("GeneralizedNonSymmetricEigenProblem: 0 dimensions matrix");
//        if (N != A.Nrows() || A.Nrows() < 1)
//            _merror("GeneralizedNonSymmetricEigenProblem: Input Matrix must be square");
//        
//        char jobvl = 'N';
//        char jobvr = 'V';
//        LAPACK_INT info;
//        
//        R.ReDimension(N);
//        I.ReDimension(N);
//        BETA.ReDimension(N);
//        Z = A;
//        
//        dggev_wrapper(&jobvl, &jobvr, &N, A.Store(), B.Store(), R.Store(), I.Store(), BETA.Store(), Z.Store(), &info);
//        //dsygv_wrapper(&jobz,&uplo,&N,Z.Store(),B.Store(),D.Store(),&info);
//        
//        if(info != 0)
//        {
//            cerr << "info is " << info << endl;
//            cout << "info is " << info << endl;
//            Error("Error in call to dggev");
//            return;
//        }
//        
//        //Transpose Z before return
//        Z = Z.t();
//    }
    
    void
    CholeskyFactorization(const MatrixRef& A)
    {
        LAPACK_INT N = A.Ncols();
        if(A.Nrows() < 1)
            _merror("CholeskyFactorization: 0 dimensions matrix");
        if (N != A.Nrows() || A.Nrows() < 1)
            _merror("CholeskyFactorization: Input Matrix must be square");
        char uplo = 'U';
        LAPACK_INT info;
        
        dpotrf_wrapper(&uplo,&N,A.Store(),&info);
        
        if(info < 0)
        {
            //cerr << "info is " << info << endl;
            cout << "Cholesky Factorization: the "<<-info<<"th argument had an illegal value!" << endl;
            cholesky_flag=0;
            //Error("dpotrf");
            //return;
        }
        else if (info>0){
            cout<< "Cholesky Factorization: the leading minor of order "<<info<<" is not positive definite!"<<endl;
            cholesky_flag=0;
        }
        else {
            if (printflag==0) {
                cout<<"Cholesky Factorization: successful!"<<endl;
            }
            
            cholesky_flag=1;
        }
    }
     
    void
    DiagonalizeRealMatrix(const MatrixRef& A, Vector& DR, Vector& DI, Matrix& VL, Matrix& VR ){
        LAPACK_INT N = A.Ncols();
        if(A.Nrows() < 1)
            _merror("make_positive_definite: 0 dimensions matrix");
        if (N != A.Nrows() || A.Nrows() < 1)
            _merror("make_positive_definite: Input Matrix must be square");
        char jobvl = 'V';
        char jobvr = 'V';
        LAPACK_INT info;
        
        dgeev_wrapper(&jobvl,&jobvr,&N,A.Store(),DR.Store(),DI.Store(),VL.Store(),VR.Store(),&info);
        
        if(info != 0)
        {
            cerr << "info is " << info << endl;
            cout << "info is " << info << endl;
            Error("Error in call to dgeev");
            return;
        }
        
        VL=VL.t();
        VR=VR.t();
        
    }
}
#endif
























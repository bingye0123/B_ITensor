#include "core.h"
#include "/home/yeba/B_ITensor/matrix/ev_lapack_wrap.h"
//#include "/home/yeba/B_ITensor/tutorial/07_basicDMRG/generalizedev.h"

using namespace itensor;
using namespace std;


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

int
main(int argc, char* argv[])
    {

    //
    // SVD of matrix M
    //

        const int Nrow = 4;
        const int Ncol = 4;
        const int maxm = min(Nrow,Ncol);
        
        Matrix M(Nrow,Ncol);
        M(1,1) = 3.9; M(1,2) = 12.5; M(1,3) = -34.5; M(1,4) = -0.5;
        M(2,1) = 4.3; M(2,2) = 21.5; M(2,3) = -47.5; M(2,4) = 7.5;
        M(3,1) = 4.3; M(3,2) = 21.5; M(3,3) = -43.5; M(3,4) = 3.5;
        M(4,1) = 4.4; M(4,2) = 26.0; M(4,3) = -46.0; M(4,4) = 6.0;
        Print(M);
        
        Matrix N(Nrow,Ncol);
        N(1,1) = 1.0; N(1,2) = 2.0; N(1,3) = -3.0; N(1,4) = 1.0;
        N(2,1) = 1.0; N(2,2) = 3.0; N(2,3) = -5.0; N(2,4) = 4.0;
        N(3,1) = 1.0; N(3,2) = 3.0; N(3,3) = -4.0; N(3,4) = 3.0;
        N(4,1) = 1.0; N(4,2) = 3.0; N(4,3) = -4.0; N(4,4) = 4.0;
        Print(N);
        
        Vector rr;
        Vector ii;
        Vector beta;
        Matrix zzl;
        Matrix zzr;
        GeneralizedNonSymmetricEigenProblem(M, N, rr, ii, beta, zzl, zzr);
        Print(rr);
        Print(ii);
        Print(beta);
        Print(zzl);
        Print(zzr);

        
        

//    Matrix U,V;
//    Vector d;
//    SVD(M,U,d,V);
//
//    Print(U);
//    Print(d);
//    Print(V);
//
//    Matrix Dtrunc(maxm,maxm);
//    Dtrunc = 0;
//
//    const int nkeep = 2;
//    for(int j = 1; j <= nkeep; ++j)
//        Dtrunc(j,j) = d(j);
//
//cout<<"D truncation"<<endl;
//Print(Dtrunc);
//
//    Matrix MM = U*Dtrunc*V;
//
//    Matrix Diff = MM-M;
//    Matrix D2 = Diff.t()*Diff;
//    Real n2 = D2(1,1) + D2(2,2) + D2(3,3);
//
//    Print(n2);
//
//    println();
//    
//
//    //
//    // SVD of two-site wavefunction
//    //
//    
//    Index s1("s1",2,Site),
//          s2("s2",2,Site);
//
//    ITensor sing(s1,s2),
//            prod(s1,s2);
//
//    //Make sing a singlet
//    sing(s1(1),s2(2)) =  1./sqrt(2);
//    sing(s1(2),s2(1)) = -1./sqrt(2);
//
//    //Make prod a product state
//    prod(s1(1),s2(2)) =  1.;
//
//    for(Real mix = 0; mix <= 1.; mix += 0.1)
//        {
//        //
//        // Create a new wavefunction that is 
//        // (1-mix) times a product state plus (mix)
//        // times a singlet (i.e. maximally entangled state).
//        //
//        // SVD this wavefunction and analyze the results.
//        // Try computing and plotting the entanglement entropy.
//        //
//        }


    return 0;
    }

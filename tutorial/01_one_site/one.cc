#include "core.h"
using namespace std;
using namespace itensor;

int
main(int argc, char* argv[]){
    Index a("L7",1000000000);
    Index b("U10",2);
    Index c("U9",2);
    Index d("L5",1);
    
    ITensor A(a);
    //A(a(1),b(1),c(1),d(1))=0.435839;
    A.randomize();
    
    Print(A);
    
    //ITensor S,V,D;
    //S=ITensor(A.indices()[1],A.indices()[3]);
    
    //svd(A, S, V, D);
    
    //Print(S);
    
    
    
//    for (int i=0; i<1; i++) {
//        if (A.indices()[i].m()==1) {
//            <#statements#>
//        }
//    }
    //Print(A.indices()[0]);
    //Print(A.indices()[1]);
    
//    a.randomize();
//    b.randomize();
    
//    ITensor x(s,primed(s)),xx(s,primed(s));
//    x.randomize("Complex");
//    xx.randomize("Complex");
//    PrintDat(x);
//    PrintDat(xx);
////    ITensor x=a+Complex_i*b;
//    Matrix rr(2,2),ii(2,2);
//    for (int i=1; i<=2; i++) {
//        for (int j=1; j<=2; j++) {
//            rr(i,j)=realPart(x)(s(i),primed(s)(i))+realPart(xx)(s(i),primed(s)(i));
//            ii(i,j)=imagPart(x)(s(i),primed(s)(i))+imagPart(xx)(s(i),primed(s)(i));
//        }
//    }
//    
//    Print(rr);

    //xx(z(1))=realPart(x)(s(1),primed(s)(1))+Complex_i*imagPart(x)(s(1),primed(s)(1));
//    Print(realPart(x)(s(1),primed(s)(1)));
//    realPart(x)(s(1),primed(s)(1))=1;
//    imagPart(x)(s(1),primed(s)(1))=2;
    
//    PrintDat(xx);
//    PrintDat(realPart(x));
//    PrintDat(imagPart(x));
//    Print(realPart(x)(s(1),primed(s)(1)));
//    Print(imagPart(x)(s(1),primed(s)(1)));
//    ITensor S,V,D;
//    
//    S=ITensor(s);
//    
//    svd(x, S, V, D);
//    PrintDat(S);
//    PrintDat(V);
//    PrintDat(V*D);
    
}


//#include "core.h"
//
//using namespace itensor;
//
//int
//main(int argc, char* argv[])
//    {
//    //
//    // Single-site wavefunction
//    //
//    
//    //Make a dimension 2 Index
//    Index s("s",2);
//
//    //Construct an ITensor
//    ITensor psi(s); //default initialized to zero
//
//
//    //
//    // Initialize up spin
//    //
//
//    //Set first element to 1.
//    psi(s(1)) = 2;
//
//    PrintData(psi);
//    
//    //exit(0); //uncomment to exit here
//
//    //
//    // Operators 
//    //
//
//    ITensor Sz(s,prime(s)),
//            Sx(s,prime(s));
//
//    commaInit(Sz,s,prime(s)) = 0.5, 0.0,
//                               0.0,-0.5;
//
//    commaInit(Sx,s,prime(s)) = 0.0, 0.5,
//                               0.5, 0.0;
//
//    PrintData(Sz);
//    PrintData(Sx);
//
//    //exit(0); //uncomment to exit here
//
//    //
//    // Product Sx * phi 
//    //
//
//    ITensor phi = Sx * psi;
//
//    phi.noprime();
//
//    PrintData(phi);
//
//    //exit(0); //uncomment to exit here
//
//    //
//    // 45* angle spin
//    //
//
//    const Real theta = Pi/4;
//
//    //Extra factors of two come from S=1/2 representation
//    psi(s(1)) = cos(theta/2.);
//    psi(s(2)) = sin(theta/2.);
//
//    PrintData(psi);
//
//    //exit(0); //uncomment to exit here
//
//    //
//    // Expectation values
//    //
//
//    ITensor cpsi = dag(prime(psi));
//
//    Real zz = (cpsi * Sz * psi).toReal();
//    Real xx = (cpsi * Sx * psi).toReal();
//
//    println("<Sz> = ", zz);
//    println("<Sx> = ", xx);
//
//    return 0;
//    }

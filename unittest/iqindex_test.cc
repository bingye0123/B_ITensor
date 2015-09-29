#include "test.h"
#include "iqindex.h"

using namespace itensor;
using namespace std;

TEST_CASE("IQIndexTest")
    {

SECTION("Null")
    {
    IQIndex i1;
    CHECK(!i1);
    CHECK_EQUAL(1,i1.m());

    IQIndex I("I",Index("i"),QN());
    CHECK(I);
    }

SECTION("Arrows")
    {
    CHECK_EQUAL(-In,Out);
    CHECK_EQUAL(-Out,In);
    }

SECTION("Primes")
    {
    IQIndex I("I",Index("i"),QN());

    I = prime(I);
    CHECK_EQUAL(I.primeLevel(),1);

    I = prime(I);
    CHECK_EQUAL(I.primeLevel(),2);

    I = prime(I,7);
    CHECK_EQUAL(I.primeLevel(),9);

    I = prime(I,-7);
    CHECK_EQUAL(I.primeLevel(),2);
    }


}

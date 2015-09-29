//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//

// input.h -- classes for reading from input files

#ifndef _input_h
#define _input_h
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include "String.h"
#include <cctype>
#define SP << " " <<
typedef double Real;
void error(const std::string& s);

class InputFile
    {
public:
    std::string filename;
    std::ifstream file;
    int opened;
    InputFile(std::string fname) : filename(fname), opened(0) {}
    void open();
    void close();
    };
std::ostream & operator << (std::ostream &s, InputFile &a);

typedef long int lint;

class InputGroup
    {
public:
    InputFile & infile;
    InputGroup * parent;
    std::string name;
    bool quiet;
    InputGroup(InputFile& inf, std::string nam,const char* c = 0)
		: infile(inf), parent(0), name(nam), quiet(false)
	{
	std::cout << "Making input group " << name;
	if(c) std::cout << ": " << c;
	std::cout << std::endl;
	}
    InputGroup(InputGroup& par, std::string nam,const char* c = 0)
		: infile(par.infile), parent(&par), name(nam), quiet(false)
	{
	std::cout << "Making input group " << parent->name << "." << name;
	if(c) std::cout << ": " << c;
	std::cout << std::endl;
	}

    int GotoGroup();		// Goes to group, then eats "{" + whitespace
    int GotoToken(std::string s);	// Goes to Token, then eats "=" + whitespace

// The following go to s, and read into i,r,t, or yes, printing c.

    int GetInt(std::string s, int& i,const char* c = 0);
    int GetLong(std::string s,lint& i,const char* c = 0);
    int GetReal(std::string s, Real& r,const char* c = 0);	
    int GetString(std::string s, std::string& t,const char* c = 0);
    int GetYesNo(std::string s, int& yes,const char* c = 0);	 // understands yes/no

// The following are mandatory versions; if they doesn't get it, we quit
    void GetIntM(std::string s, int& i,const char* c = 0);	
    void GetLongM(std::string s, lint& i,const char* c = 0);	
    void GetRealM(std::string s, Real& r,const char* c = 0);
    void GetStringM(std::string s, std::string& t,const char* c = 0);
    void GetYesNoM(std::string s, int& yes,const char* c = 0);

    void SkipLine();
    };

    /*
InputGroup(InputFile& inf, String nam,const char* c)
	    : infile(inf), name(nam), parent(0) 
    {
    std::cout << "Making input group " << name;
    if(c) std::cout << ", " << c;
    std::cout << std::endl;
    }
    */

/* 
To read in a table:

in input file:

tablename
    {
    a	x	y	z
    1	4.0	5.0	7.0
    2	4.0	2.0	7.0
    4	3.0	5.0	7.0
    5	4.0	5.0	1.0
    }

Then in program:
    InputGroup table(parent,"tablename");
    if(table.GotoGroup())
	{
	table.SkipLine();
	for(int i = 1; i <= n; i++)
	    table.infile.file >> a[i] >> x[i] >> y[i] >> z[i];
	}
*/

int gettoken(std::istream& is, std::string& s);

#endif

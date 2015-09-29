//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "input.h"
#include <algorithm>
using std::ostream;
using std::istream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;


void InputFile::open()
    {
    close();
    file.open(filename.c_str());
    if(!file) 
	error("can't open input file");
    opened = 1;
    }

void InputFile::close()
    {
    if(opened)
	{ file.clear(); file.close(); }
    opened = 0;
    }

ostream & operator << (ostream &s, InputFile &a)
    {
    a.open();
    char c;
    s << "Input filename is " << a.filename << endl;
    while(a.file.get(c))
	s << c;
    a.close();
    return s;
    }

void eatwhite(istream& is)
    {
    char c;
    while(is.get(c))
	{
	if(!isspace(c))
	    {
	    is.putback(c);
	    break;
	    }
	}
    }

void matchbracket(istream& is)
    {
    int nest = 1;
    char c;
    while(is >> c)
	{
	if(c == '{')
	    nest++;
	else if(c == '}')
	    nest--;
	// cout << "match: " << c << endl;
	if(nest == 0) break;
	}
    if(nest != 0)
	error("unterminated bracket");
    }

int gettoken(istream& is, string& s)
    {
    char c[512], ch;
    eatwhite(is);
    for(int i = 0; i < 512; i++)
	{
	if(!is.get(ch))
	    return 0;
	// cout << "gettoken: " << ch << endl;
	if(isalpha(ch) || ch == '_' || (i != 0 && isdigit(ch)))
	    c[i] = ch;
	else
	    {
	    if(i == 0)
		{
		c[0] = ch;
		c[1] = '\0';
		s = c;
		return 1;
		}
	    is.putback(ch);
	    c[i] = '\0';
	    break;
	    }
	}
    s = c;
    return 1;
    }

int InputGroup::GotoGroup()
    {
    if(parent != 0)
	parent->GotoGroup();
    else
	infile.open();
    eatwhite(infile.file);
    while(1)
	{
	string s;
	if(!gettoken(infile.file,s))
	    return 0;
	// cout << "GotoGroup: got string " << s << endl;
	if(s == name)
	    {
	    eatwhite(infile.file);
	    char c;
	    infile.file.get(c);
	    if(c != '{')
		error("bracket does not follow group name");
	    eatwhite(infile.file);
	    return 1;
	    }
	else
	    {
	    if(s[0] == '{')
		{
		matchbracket(infile.file);
		}

	    while(1)
		{
		char c;
		infile.file.get(c);
		if(c == '{')
		    matchbracket(infile.file);
		if(c == '\n')
		    {
		    eatwhite(infile.file);
		    break;
		    }
		}
	    }
	}
    }

int InputGroup::GotoToken(string s)
    {
    if(!GotoGroup()) 
	return 0;
    while(1)
	{
	string t;
	if(!gettoken(infile.file,t))
	    {
	    // cout << "no token, returning from GotoToken" << endl;
	    return 0;
	    }
	// cout << "GotoToken: got string " << t << endl;
	if(t[0] == '}') return 0;
	if(t == s) 
	    {
	    eatwhite(infile.file);
	    char c;
	    if(!(infile.file.get(c))) return 0;
	    if(c != '=')
		return 0;
	    eatwhite(infile.file);
	    return 1;
	    }
	if(t[0] == '{')
	    {
	    matchbracket(infile.file);
	    }
	while(1)
	    {
	    char c;
	    if(!(infile.file.get(c))) return 0;
	    if(c == '\n' || c == ',' || c == ';')
		{
		eatwhite(infile.file);
		break;
		}
	    if(c == '{')
		matchbracket(infile.file);
	    }
	}
    }

int InputGroup::GetInt(string s, int& res,const char* c)
    {
    if(!GotoToken(s) || !(infile.file >> res)) 
	{
	if(!quiet) 
    {
    cout << "Couldnt get " << name << "." << s;
	if(c) cout  << ": " << c;
	cout << endl;
    }
	return 0;
	}
    if(!quiet) 
    {
    cout << "Got " << name << "." << s << " = " << res;
    if(c) cout  << ": " << c;
    cout << endl;
    }
    return 1;
    }

int InputGroup::GetLong(string s, lint& res,const char* c)
    {
    if(!GotoToken(s) || !(infile.file >> res)) 
	{
	if(!quiet) 
    {
    cout << "Couldnt get " << name << "." << s;
	if(c) cout  << ": " << c;
	cout << endl;
    }
	return 0;
	}
    if(!quiet) 
    {
    cout << "Got " << name << "." << s << " = " << res;
    if(c) cout  << ": " << c;
    cout << endl;
    }
    return 1;
    }

int InputGroup::GetReal(string s, Real& res,const char* c)
    {
    if(!GotoToken(s) || !(infile.file >> res)) 
	{
	if(!quiet) 
    {
    cout << "Couldnt get " << name << "." << s;
	if(c) cout  << ": " << c;
	cout << endl;
    }
	return 0;
	}
    if(!quiet) 
    {
    cout << "Got " << name << "." << s << " = " << res;
    if(c) cout  << ": " << c;
    cout << endl;
    }
    return 1;
    }

int InputGroup::GetString(string s, string& res,const char* c)
    {
    if(!GotoToken(s) || !(infile.file >> res)) 
	{
	if(!quiet) 
    {
    cout << "Couldnt get " << name << "." << s;
	if(c) cout  << ": " << c;
	cout << endl;
    }
	return 0;
	}
    if(!quiet) 
    {
    cout << "Got " << name << "." << s << " = " << res;
    if(c) cout  << ": " << c;
    cout << endl;
    }
    return 1;
    }

char mydolower(char c) { return tolower(c); }

int InputGroup::GetYesNo(string s, int& yes,const char* c)
    {
    string res;
    if(!GotoToken(s) || !(infile.file >> res)) 
	{
	if(!quiet) 
    {
    cout << "Couldnt get " << name << "." << s;
	if(c) cout  << ": " << c;
	cout << endl;
    }
	return 0;
	}
    if(!quiet) 
    {
    cout << "Got " << name << "." << s << " = " << res;
    if(c) cout  << ": " << c;
    cout << endl;
    }
    transform(res.begin(),res.end(),res.begin(),mydolower);
    if(res == "yes" || res == "y")
	{
	yes = 1;
	return 1;
	}
    if(res == "no" || res == "n")
	{
	yes = 0;
	return 1;
	}
    return 0;
    }

void InputGroup::SkipLine()
    {
    char c = '\0';
    while(c != '\n')
	infile.file.get(c);
    eatwhite(infile.file);
    }

void InputGroup::GetIntM(string s, int& res,const char* c)
    {
    if(!GetInt(s,res,c))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::GetLongM(string s, lint& res,const char* c)
    {
    if(!GetLong(s,res,c))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::GetRealM(string s, Real& res,const char* c)
    {
    if(!GetReal(s,res,c))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::GetStringM(string s, string& res,const char* c)
    {
    if(!GetString(s,res,c))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::GetYesNoM(string s, int& yes,const char* c)
    {
    if(!GetYesNo(s,yes,c))
        error("mandatory item: " + s + ", exiting");
    }


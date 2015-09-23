#ifndef Paramizer_IS_INCLUDED
#define Paramizer_IS_INCLUDED


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <list>

using namespace std;


class Parameter {

 protected:
  string str;
  string name;
  string value;
  int    _parsed;

 public:
  Parameter(string _str);
  int     parsed()   {return _parsed;}
  string  getStr()   {return str;}
  string  getName()  {return name;}
  string& getValue() {return value;}
  void    summary();

 protected:
  int    parse();
};



class Paramizer {

 protected:
  list<Parameter> inputParamList;
  list<Parameter> outputParamList;

 public:
  Paramizer(int argn, char* args[]);
  Paramizer();
  ~Paramizer();
  void init(int argn, char* args[]);
  void getInt(int& var, char* paramName, int def = 0);
  void getDouble(double& var, char* paramName, double def = 0.);
  void getString(string& var, char* paramName, string def = "");
  void summary();

 protected:
  void parse();
  void parseFile(string filename);
  void setParam(Parameter param);

};


#endif

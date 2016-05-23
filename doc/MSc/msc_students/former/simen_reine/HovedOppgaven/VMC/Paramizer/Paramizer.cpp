#include "Paramizer.h"


Parameter::Parameter(string _str) : str(_str) {
  _parsed = parse();
}


int Parameter::parse() {
  int n;
  if ((n = str.find_first_of('=')) >= 0) {
    name  = str.substr(0,n);
    value = str.substr(n+1);
    return 1;
  }
  return 0;
}


void Parameter::summary() {
  cerr << endl;
  if (_parsed == 1)
    cerr << endl << "Name:  " << "\"" << name  << "\""
	 << endl << "Value: " << "\"" << value << "\""
	 << endl;
  else
    cerr << str << endl;    
}


/************************************************************/




Paramizer::Paramizer(int argn, char* args[]) {
  for (int i = 1; i < argn; i++)
    inputParamList.push_back(Parameter(string(args[i])));
  parse();
}


Paramizer::Paramizer() {
}

Paramizer::~Paramizer() {
}


void Paramizer::init(int argn, char* args[]) {
  for (int i = 1; i < argn; i++)
    inputParamList.push_back(Parameter(string(args[i])));
  parse();
}


void Paramizer::parse() {

  for ( list<Parameter>::iterator iter = inputParamList.begin(); iter != inputParamList.end(); iter++ )
    if ( iter->parsed() == 0 )
      if ( iter->getStr() == string("-f") ) {
	iter++;
	if (iter == inputParamList.end()) {
	  cerr << "\nFilename missing.\n";
	  exit(1);
	}
	parseFile(iter->getStr());
      }
      else {
	cerr << "\nCannot resolve input string \"" << iter->getStr() << "\"\n";
	exit(1);
      }
    else
      setParam(*iter);

}


void Paramizer::parseFile(string filename) {
  ifstream file( filename.c_str() );

  if(!file) {
    cerr << "\nFile \"" << filename << "\" not found.\n";
    exit(1);
  }

  string s;
  
  while( (file >> s).good() ) {
    Parameter param(s);
    if (param.parsed() == 0) {
      cerr << "\nError while parsing file \"" << filename << "\".\n";
      exit(1);
    }
    setParam(param);
  }
}


void Paramizer::setParam(Parameter param) {
  int found = 0;
  for ( list<Parameter>::iterator iter = outputParamList.begin(); iter != outputParamList.end(); iter++ )
    if (iter->getName() == param.getName()) {
      iter->getValue() = param.getValue();      
      found = 1;
      break;
    }
  if (found == 0)
    outputParamList.push_back(param);
}


void Paramizer::getInt(int& var, char* paramName, int def) {
  int found = 0;
  for ( list<Parameter>::iterator iter = outputParamList.begin(); iter != outputParamList.end(); iter++ )
    if (iter->getName() == paramName) {
      var = atoi(iter->getValue().c_str());
      found = 1;
      break;
    }
  if (found == 0)
    var = def;
}


void Paramizer::getDouble(double& var, char* paramName, double def) {
  int found = 0;
  for ( list<Parameter>::iterator iter = outputParamList.begin(); iter != outputParamList.end(); iter++ )
    if (iter->getName() == paramName) {
      var = atof(iter->getValue().c_str());
      found = 1;
      break;
    }
  if (found == 0)
    var = def;
}


void Paramizer::getString(string& var, char* paramName, string def) {
  int found = 0;
  for ( list<Parameter>::iterator iter = outputParamList.begin(); iter != outputParamList.end(); iter++ )
    if (iter->getName() == paramName) {
      var = iter->getValue();
      found = 1;
      break;
    }
  if (found == 0)
    var = def;
}


void Paramizer::summary() {
  cerr << endl << "Number of parameters: " << outputParamList.size() << endl;
  for( list<Parameter>::iterator iter = outputParamList.begin(); iter != outputParamList.end(); iter++ )
    iter->summary();
}


// Spørsmål: Blir s allokert hver gang/ noen gang i det hele tatt?
//           Hva med var i getString?

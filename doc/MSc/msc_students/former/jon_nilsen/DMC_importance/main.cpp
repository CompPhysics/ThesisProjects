#include "DMC_importance.hpp"

using namespace std;

int main(int argc, const char* argv[]){

  if(argc!=3){
    cerr << "usage: gpp.app infile outfile" << endl;
    return 1;
  }

  string infile=argv[1];
  
  string outfile=argv[2];
  
  DMC obj(argv[1], argv[2]);
  
  obj.diffMC();
  
  return 0;

}

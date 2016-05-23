#include "functions.hpp"

double Helium::valuePt(QickArray &pos, QickArray &params){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  if(dim!=2 || cols!=2)
    std::cerr << "pos must be a 2d 2 particles x pos array" << std::endl;

  static QickArray dummy(rows);
  double ans=1;
  for(int i=0; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      dummy(j)=pos(i,j);
    ans*=wave1s(dummy, params);
  }
  
  for(int j=0; j!=rows; j++)
    dummy(j)=pos(1,j)-pos(0,j);

  return ans*corr(dummy,params);

}


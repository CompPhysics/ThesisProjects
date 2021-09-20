// package of functors

#ifndef functions_hpp_IS_INCLUDED
#define functions_hpp_IS_INCLUDED

#include <cmath>
#include "QickArray/QickArray.hpp"

using namespace std;

class Func{
  
public:

  virtual double valuePt(QickArray *pos){ return 42; }

  virtual double valuePt(QickArray &pos, QickArray &params){ return 42; }

  virtual double valuePt(Func &func, QickArray &pos, 
			 QickArray &params){ return 42; }

  virtual double valuePt(Func &func1, Func &func2, 
			 QickArray &pos, QickArray &params){ return 42; }

  virtual double valuePt(Func &func1, Func &func2, QickArray &pos1,
			 QickArray &pos2, QickArray &params,
			 double D, double tau){ return 42; }

  virtual double valuePt(Func &func1, Func &func2, Func &Func3,
			 QickArray &pos, QickArray &params){ return 42; }

  virtual QickArray& valuePt(Func &func, QickArray &pos, 
			     QickArray &dummy, QickArray &params){ 
    return dummy; 
  }

  virtual inline double operator()(QickArray &pos){ return 42; }

  virtual inline double operator()(QickArray &pos, QickArray &params){ 
    return 42; }

  virtual inline double operator()(Func &func, QickArray &pos
				   , QickArray &params){ return 42; }

  virtual inline double operator()(Func &func1, Func &func2, 
				   QickArray &pos, QickArray &params)
  { return 42; }

  virtual inline double operator()(Func &func1, Func &func2, QickArray &pos1,
				   QickArray &pos2, QickArray &params,
				   double D, double tau)
  { return 42; }

  virtual inline double operator()(Func &func1, Func &func2, Func &func3,
				   QickArray &pos, QickArray &params)
  { return 42; }

  virtual inline QickArray& operator()(Func &func, QickArray &pos, 
				      QickArray &dummy, QickArray &params){ 
    return dummy; 
  }

  virtual inline double sqr(double x){ return x*x; }

};


class Wave : public virtual Func{

public:

  virtual double valuePt(QickArray &pos, QickArray &params){return 42.;}

  virtual inline double operator()(QickArray &pos, QickArray &params){ 
    return valuePt(pos, params);
  }

};

class Wave1s : public virtual Wave{

public:

  virtual double valuePt(QickArray &pos, QickArray &params);

  virtual inline double operator()(QickArray &pos, QickArray &params){ 
    return valuePt(pos, params);
  }
};

// exp(-alpha*(r))
inline double Wave1s:: valuePt(QickArray &pos, QickArray &params){

  int dim=pos.no_of_elements();
  double ans=0;
  for(int i=0; i!=dim; i++)
    ans+=sqr(pos(i));

  return exp(-params(0)*sqrt(ans));

}

class Correlation : public virtual Wave{

public:

  virtual double valuePt(QickArray &pos, QickArray &params);

  virtual inline double operator()(QickArray &pos, QickArray &params){ 
    return valuePt(pos, params);
  }
};

// e^(r_ij/(1+params(1)*r_ij)) nb! pos=r_ij
inline double Correlation:: valuePt(QickArray &pos, QickArray &params){

  int dim=pos.no_of_elements();
  double r_ij=0;
  for(int i=0; i!=dim; i++)
    r_ij+=sqr(pos(i));
  r_ij = sqrt(r_ij);

  return exp(.5*r_ij/(1+params(1)*r_ij));

}

class Helium : public virtual Wave{

private:
  Wave1s wave1s;
  Correlation corr;
public:

  virtual double valuePt(QickArray &pos, QickArray &params);

  virtual inline double operator()(QickArray &pos, QickArray &params){ 
    return valuePt(pos, params);
  }
};

class Potential : public virtual Func{

public:

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }

};

// V = -sum(number_particles/r_i) + sum(1/r_ij)
inline double Potential:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e=0;
  /* 
  ** adding potential energy; 
  ** V = -sum(number_particles/r_i) + sum(1/r_ij)
  */
  double r_single_particle, r_temp;
  /* contribution from electron-proton potential  */
  for (int i = 0; i != cols; i++){
    r_single_particle = 0;
    for (int j = 0; j != rows; j++) {
      r_temp = pos(i, j);
      r_single_particle += r_temp*r_temp;
    }
    potential_e -= cols/sqrt(r_single_particle);
  }
  
  /* contribution from electron-electron potential  */
  double r_ij, r_old_ik, r_old_jk;
  for (int i = 0; i != cols-1; i++) { 
    for (int j = i+1; j != cols; j++) {
      r_ij = 0;  
      for (int k = 0; k != rows; k++) { 
	r_old_ik = pos(i,k);
	r_old_jk = pos(j,k);
	r_ij += sqr(r_old_ik-r_old_jk);
      }

      potential_e += 1/sqrt(r_ij);
    }
  }

  return potential_e;

}

// functor evaluating first derivate returnin positions x dim array
class Nabla : public virtual Func{

public:

  virtual QickArray& valuePt(Func &func, QickArray &pos, 
			     QickArray &dummy, QickArray &params);
  
  virtual inline QickArray& operator()(Func &func, QickArray &pos,
				      QickArray &dummy, QickArray &params){ 
    
    return valuePt(func, pos, dummy, params);
  }

};

inline QickArray& Nabla:: valuePt(Func &func, QickArray &pos, 
				  QickArray &dummy, QickArray &params){

  const double h=1e-7;

  int cols,rows,dim;

  dim = pos.get_dimension_info(&cols,&rows);

  double ans=0;

  dummy *= 0;
  // using centered difference
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      // f(x-h)
      pos(i,j) -= h;
      dummy(i,j) -= func(pos, params);
      // f(x+h)
      pos(i,j) += 2*h;
      dummy(i,j) += func(pos, params);
      // f(x)
      pos(i,j) -= h;
    }
  
  dummy /= 2*h;

  return dummy;

}

// functor evaluating second derivate summing over all positions
class Nabla2 : public virtual Func{

public:

  virtual double valuePt(Func &func, QickArray &pos, QickArray &params);

  virtual inline double operator()(Func &func, QickArray &pos, 
				   QickArray &params){ 

    return valuePt(func, pos, params);
  }

};

inline double Nabla2:: valuePt(Func &func, QickArray &pos, QickArray &params){

  const double h=5e-5;

  int cols,rows,dim;

  dim = pos.get_dimension_info(&cols,&rows);

  double ans=0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      // f(x-h)
      pos(i,j)-=h;
      ans += func(pos, params);
      // f(x+h)
      pos(i,j)+=2*h;
      ans += func(pos, params);
      // f(x)
      pos(i,j)-=h;
      ans -= 2*func(pos, params);
    }
  
  return ans/h/h;

}


// functor evaluating local energy
class LocalEnergy : public virtual Func{

public:

  virtual double valuePt(Func &wf, Func &ham, 
			 Func &pot, QickArray &pos, QickArray &params);

  virtual inline double operator()(Func &wf, Func &ham,
				   Func &pot, QickArray &pos, 
				   QickArray &params){ 

    return valuePt(wf, ham, pot, pos, params);
  }

};

inline double LocalEnergy:: valuePt(Func &wf, Func &ham, 
				    Func &pot, QickArray &pos, 
				    QickArray &params){

  return ham(wf, pot, pos, params)/wf(pos, params);
}

// functor evaluating quantum force
class QForce : public virtual Func{

public:

  virtual QickArray& valuePt(Func &wf, QickArray &pos, 
			     QickArray &dummy, QickArray &params);
  
  virtual inline QickArray& operator()(Func &wf, QickArray &pos,
				      QickArray &dummy, QickArray &params){ 
    
    return valuePt(wf, pos, dummy, params);
  }

};

inline QickArray& QForce:: valuePt(Func &wf, QickArray &pos, 
				   QickArray &dummy, QickArray &params){

  Nabla nabla;
  
  dummy = nabla(wf, pos, dummy, params);

  dummy *= 2/wf(pos, params);

  return dummy;

}

// functor evaluating hamiltonian
class Hamilton : public virtual Func{

public:

  virtual double valuePt(Func &wf, Func &pot, QickArray &pos, 
			 QickArray &params);

  virtual inline double operator()(Func &wf, Func &pot, QickArray &pos,
				   QickArray &params){ 

    return valuePt(wf, pot, pos, params);
  }

};

inline double Hamilton:: valuePt(Func &wf, Func &pot, 
				 QickArray &pos, QickArray &params){

  Nabla2 kin;

  // H*psi=-.5*nabla2*psi+V*psi
  return -.5*kin(wf, pos, params)+pot(pos)*wf(pos, params);

}

// functor evaluating green's function
class Greens : public virtual Func{

public:

  virtual double valuePt(Func &wf, Func &q_force, QickArray &pos, 
			 QickArray &new_pos, QickArray &params,
			 double D, double tau);

  virtual inline double operator()(Func &wf, Func &q_force, QickArray &pos,
				   QickArray &new_pos, QickArray &params,
				   double D, double tau){ 

    return valuePt(wf, q_force, pos, new_pos, params, D, tau);
  }

};

inline double Greens:: valuePt(Func &wf, Func &q_force, 
			       QickArray &pos, QickArray &new_pos, 
			       QickArray &params, double D, double tau){

  int cols,rows,dim;

  dim = pos.get_dimension_info(&cols,&rows);

  static QickArray qf(cols, rows);

  qf = q_force(wf, pos, qf, params);

  double G = 0;
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++)
      G += sqr(pos(i,j) - new_pos(i,j) - D*tau*qf(i,j));
  // H*psi=-.5*nabla2*psi+V*psi
  return exp(-G*D/tau);

}

#endif

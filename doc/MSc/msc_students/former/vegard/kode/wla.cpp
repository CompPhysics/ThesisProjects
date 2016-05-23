#include "wla.h"
#include <iostream>//Just for testing purposes...
#include <cmath>

dwt::dwt(int width, int filter){
  l = width;
  type = filter;
  g = new double[l];
  h = new double[l];
  init();
}

dwt::~dwt(){
  delete[] g;
  delete[] h;
}

void dwt::init(){
  if(type == 1){             
    for(int i = 0; i < l; i++){ 
      g[i] = get_daub(l, i);
    }
  }else if(type == 2){
    for(int i = 0; i < l; i++){
      g[i] = get_la(l, i);
    }
    nu = get_nu(l, type);
  }else if(type == 3){
    for(int i = 0; i < l; i++){
      g[i] = get_bl(l, i);
    }
    nu = get_nu(l, type);
  }else if(type == 4){
    for(int i = 0; i < l; i++){
      g[i] = get_coiflet(l, i);
    }
    nu = get_nu(l, type);
  }

  for(int i = 0; i < l; i++){
    h[i] = pow(-1,i)*g[l-1-i];
  }
}

void dwt::transform(int n, double *v_in, double *v_out, double *w_out){
  int k;
  int nh = n/2;

  for(int i = 0; i < nh; i++){
    k = 2*i+1;
    v_out[i] = g[0]*v_in[k];
    w_out[i] = h[0]*v_in[k];
    for(int j = 1; j < l; j++){
      k--;
      if(k < 0){
	k = n-1;
      }
      v_out[i] += g[j]*v_in[k];
      w_out[i] += h[j]*v_in[k];
    }
  }
}

void dwt::inverse_transform(int n, double *v_in, double *w_in, double *x_out){
  int k, m, o;
  int p = -2; 
  int q = -1;  
  int lh = l/2;
  
  for(int i = 0; i < n; i++){
    k = i;
    m = 0;
    o = 1;
    p += 2;
    q += 2;
    x_out[p] = h[o]*w_in[k]+g[o]*v_in[k];
    x_out[q] = h[m]*w_in[k]+g[m]*v_in[k];
    if(l > 2){
      for(int j = 1; j < lh; j++){
	k++;
	if(k >= n){
	  k = 0;
	}
	m += 2;
	o += 2;
	x_out[p] += h[o]*w_in[k]+g[o]*v_in[k];
	x_out[q] += h[m]*w_in[k]+g[m]*v_in[k];
      }
    }
  }
}

void dwt::partial_transform(int n, int j_max, double *x_in, double **wv){
  int nh = n/2;
  double *v_in = new double[n];
  double *v_out = new double[nh];
  double *w = new double[nh];

  for(int i = 0; i < n; i++){
    v_in[i] = x_in[i];
  }

  for(int i = 0; i < j_max; i++){
    dwt::transform(n/pow(2,i), v_in, v_out, w);
    for(int j = 0; j < n; j++){
      v_in[j] = v_out[j];
      wv[i][j] = w[j];
    }
    n /= 2;
  }
  n *= 2;
    for(int j = 0; j < n; j++){
    wv[j_max][j] = v_out[j];
  }

  delete[] v_in;
  delete[] v_out;
  delete[] w;
}

void dwt::mra(int n, int j_max, double **wv, double **ds){//Take w_matrix as input!!!
  double *x_out = new double[n];
  double *zero = new double[n];

  for(int i = 0; i < n; i++){
    zero[i] = 0;
  }

  double *tmp = new double[n];
  //Compute details D_j:
  for(int i = 0; i < j_max; i++){      
    dwt::inverse_transform(n/pow(2,i), zero, wv[i], x_out);//Check n/pow(2,x) throughout...
    for(int j = i-1; j >= 0; j--){
      for(int k = 0; k < n; k++){
	tmp[k] = x_out[k];
      }
      dwt::inverse_transform(n/pow(2,j), tmp, zero, x_out);
    }
    for(int j = 0; j < n; j++){
      ds[i][j] = x_out[j];
    }
  }

  //Compute smooth S_j0:
  dwt::inverse_transform(n/pow(2,j_max-1), wv[j_max], zero, x_out);
  for(int k = j_max-2; k >= 0; k--){
    for(int i = 0; i < n; i++){
      tmp[i] = x_out[i];
    }
    dwt::inverse_transform(n/pow(2,k), tmp, zero, x_out);
  }
  for(int j = 0; j < n; j++){
    ds[j_max][j] = x_out[j];
  }

  delete[] tmp;
  delete[] x_out;
  delete[] zero;
}

void dwt::wavevar(int bias, int n, int j_max, int p, double *edof, double *wvar, double *ci, double **w){
  int l_j, m_j, n_j, level, ii;                                   
  double nusq, phi, eta, q_lower, q_upper;                              

  if(bias == 0){
    for(int i = 0; i < j_max; i++){
      nusq = .0;
      level = i+1;
      l_j = ceil((l-2)*(1-(1./pow(2,level))));
      n_j = floor(double(n)/pow(2,level)-1);
      m_j = n_j-l_j+1;
      for(int j = l_j; j < n_j; j++){
	nusq += w[i][j]*w[i][j];
      }
      wvar[i] = nusq/(m_j*pow(2,level));
    }
  }else if(bias == 1){
    for(int i = 0; i < j_max; i++){
      nusq = .0;
      level = i+1;
      n_j = floor(double(n)/pow(2,level)-1)+1;
      for(int j = 0; j < n_j; j++){
	nusq += w[i][j]*w[i][j];
      }
      wvar[i] = nusq/(n_j*pow(2,level));
    }
  }

  //Calculate confidence interval
  if(p == 90){
    phi = 1.6449;
  }else if(p == 95){
    phi = 1.96;
  }else if(p == 99){
    phi = 2.5758;
  }
  
  for(int i = 0; i < j_max; i++){
    eta = edof[i];
    q_lower = eta*pow((1-2./(9*eta)+phi*sqrt(2./(9*eta))),3);
    phi *= -1;
    q_upper = fabs(eta*pow((1-2./(9*eta)+phi*sqrt(2./(9*eta))),3));
    phi *= -1;
    ii = j_max+i;
    ci[i] = log(eta*wvar[i]/q_lower);
    ci[ii] = log(eta*wvar[i]/q_upper);
  }
}

void dwt::edof(int bias, int n, int j_max, double *eta){//Sjekk likning øverst side 345. Ikke dele på 2^j???
  int j, l_j;
  double m_j, n_j;
  
  if(bias == 0){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      l_j = ceil((l-2)*(1-(1./pow(2,j))));
      n_j = floor(double(n)/pow(2,j)-1);
      m_j = n_j-l_j+1;
      eta[i] = m_j/pow(2,j);
      if(eta[i] < 1){
	eta[i] = 1;
      }
    }
  }else if(bias == 1){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      n_j = floor(double(n)/pow(2,j)-1);
      eta[i] = n_j/pow(2,j);
      if(eta[i] < 1){
	eta[i] = 1;
      }
    }
  }
}

int dwt::get_lj(int j){ //All three classes needs this one...
  return ((pow(2,j)-1)*(l-1)+1);
}

int dwt::g_shift(int j){        //Make these private
  j++;
  int l_j = ((pow(2,j)-1)*(l-1)+1);//dwt::get_lj(j);
  int nu_j = nu*(l_j-1)/(l-1);//Check with formulaes p147
  int gamma_j = ceil(((fabs(nu_j)+1)/pow(2,j))-1); 

  return gamma_j;
}

int dwt::h_shift(int j){        //Make these private
  j++;
  int l_j = ((pow(2,j)-1)*(l-1)+1);
  int nu_j = -(l_j/2+l/2+nu-1);
  int gamma_j = ceil(((fabs(nu_j)+1)/pow(2,j))-1);
  
  return gamma_j;
}

void dwt::shift(int j_max, int n, double **wv){
  int shift, n_j;                                          
  double *tmp = new double[n/2];

  for(int i = 0; i < j_max; i++){
    n_j = n/pow(2,i+1);
    shift = dwt::h_shift(i); //wavelet
    for(int j = shift; j < n_j; j++){
      tmp[j-shift] = wv[i][j];
    }
    for(int j = 0; j < shift; j++){
      tmp[n_j-shift+j] = wv[i][j];
    }
    for(int j = 0; j < n_j; j++){
      wv[i][j] = tmp[j];
    }
  }
  shift = dwt::g_shift(j_max-1);//scaling
  for(int j = shift; j < n_j; j++){//Turn i and j around to make it more clear...
    tmp[j-shift] = wv[j_max][j];
  }
  for(int j = 0; j < shift; j++){
    tmp[n_j-shift+j] = wv[j_max][j];
  }
  for(int j = 0; j < n_j; j++){
    wv[j_max][j] = tmp[j];
  }
  

  delete[] tmp;
}

void dwt::get_boundaries(int param, int j_max, int n, double *boundaries){
  int j, l_j, n_j, gamma, gamma_hat;
  
  for(int i = 0; i <= 2*j_max+1; i++){//initialize boundaries
    boundaries[i] = -1;
  }

  if(param == 0){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      l_j = ceil((l-2)*(1-(1./pow(2,j))));                                //require l_j <= n_j
      boundaries[i] = l_j-1;
    }
    boundaries[2*j_max] = l_j-1;
  }else if(param == 1){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      l_j = ceil((l-2)*(1-(1./pow(2,j))));
      n_j = n/pow(2,j);
      gamma = dwt::h_shift(i);
      gamma_hat = l_j-gamma;
      boundaries[i] = gamma_hat-1;//Koeffisient til og med gamma-1 er påvirket av grenser...
      boundaries[i+j_max] = n_j-gamma;//Fra og med
    }
    gamma = dwt::g_shift(j_max-1);
    gamma_hat = l_j-gamma;
    boundaries[2*j_max] = gamma_hat-1;
    boundaries[2*j_max+1] = n_j-gamma;
  }else if(param == 2){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      l_j = ceil((l-2)*(1-(1./pow(2,j))));
      boundaries[i] = pow(2,j)*l_j-1;
      l_j = ((pow(2,j)-1)*(l-1)+1);
      boundaries[i+j_max] = n-(l_j-pow(2,j));
    }
    l_j = ceil((l-2)*(1-(1./pow(2,j))));
    boundaries[2*j_max] = pow(2,j)*l_j-1;
    l_j = ((pow(2,j)-1)*(l-1)+1);
    boundaries[2*j_max+1] = n-(l_j-pow(2,j));
  }
}

modwt::modwt(int width, int filter){
  l = width;
  type = filter;
  g = new double[l];
  h = new double[l];
  init();
}

modwt::~modwt(){
  delete[] g;
  delete[] h;
}

void modwt::init(){
  if(type == 1){              //Make switch case instead of if else?
    for(int i = 0; i < l; i++){
      g[i] = get_daub(l, i)/sqrt(2);
    }
  }else if(type == 2){
    for(int i = 0; i < l; i++){
      g[i] = get_la(l, i)/sqrt(2);
    }
    nu = get_nu(l, type);
  }else if(type == 3){
    for(int i = 0; i < l; i++){
      g[i] = get_bl(l, i)/sqrt(2);
    }
    nu = get_nu(l, type);
  }else if(type == 4){
    for(int i = 0; i < l; i++){
      g[i] = get_coiflet(l, i)/sqrt(2);
    }
    nu = get_nu(l, type);
  }

  for(int i = 0; i < l; i++){
    h[i] = pow(-1,i)*g[l-1-i];
  }
}

void modwt::transform(int n, int it, double *v_in, double *v_out, double *w_out){
  int k;
  
  for(int i = 0; i < n; i++){
    k = i;
    v_out[i] = g[0]*v_in[k];
    w_out[i] = h[0]*v_in[k];
    for(int j = 1; j < l; j++){
      k -= pow(2,it);
      if(k < 0){
	k += n;
      }
      v_out[i] += g[j]*v_in[k];
      w_out[i] += h[j]*v_in[k];
    }
  }
}

void modwt::inverse_transform(int n, int it, double *v_in, double *w_in, double *x_out){
  int k;
  
  for(int i = 0; i < n; i++){
    k = i;
    x_out[i] = h[0]*w_in[k]+g[0]*v_in[k];
    for(int j = 1; j < l; j++){
      k += pow(2, it);
      if(k >= n){
	k -= n;
      }
      x_out[i] += h[j]*w_in[k]+g[j]*v_in[k];
    }
  }
}

void modwt::mra(int n, int j_max, double **wv, double **ds){
  double *x_out = new double[n];
  double *zero = new double[n];
  double *tmp = new double[n];
  
  for(int i = 0; i < n; i++){
    zero[i] = 0;
  }

  //Compute details D_j:
  for(int i = 0; i < j_max; i++){
    modwt::inverse_transform(n, i, zero, wv[i], x_out);
    for(int j = i-1; j >= 0; j--){
      for(int k = 0; k < n; k++){
	tmp[k] = x_out[k];
      }
      modwt::inverse_transform(n, j, tmp, zero, x_out);
    }
    for(int j = 0; j < n; j++){
      ds[i][j] = x_out[j];
    }
  }
  
  //Compute smooth S_j0:
  modwt::inverse_transform(n, j_max-1, wv[j_max], zero, x_out);
  for(int k = j_max-2; k >= 0; k--){
    for(int i = 0; i < n; i++){
      tmp[i] = x_out[i];
    }
    modwt::inverse_transform(n, k, tmp, zero, x_out);
  }
  for(int j = 0; j < n; j++){
    ds[j_max][j] = x_out[j];
  }

  delete[] x_out;
  delete[] zero;
  delete[] tmp;
}

void modwt::partial_transform(int n, int j_max, double *x_in, double **wv){
  double *v_in = new double[n];
  double *v_out = new double[n];
  double *w = new double[n];

  for(int i = 0; i < n; i++){
    v_in[i] = x_in[i];
  }
  for(int i = 0; i < j_max; i++){ 
    modwt::transform(n, i, v_in, v_out, w);
    for(int j = 0; j < n; j++){
      v_in[j] = v_out[j];
      wv[i][j] = w[j];
    }
  }
  for(int j = 0; j < n; j++){
    wv[j_max][j] = v_out[j];
  }
  
  delete[] v_in;
  delete[] v_out;
  delete[] w;
}

void modwt::wavevar(int bias, int n, int j_max,  int p, double *edof, double *wvar, double *ci, double **w){
  int l_j, m_j, ii;
  double nusq, phi, eta, q_lower, q_upper;

  if(bias == 0){//unbiased estimator requires N - L_j >= 0
    for(int i = 0; i < j_max; i++){
      nusq = .0;
      l_j = get_lj(i);
      m_j = n-l_j+1;
      l_j--;
      for(int j = l_j; j < n; j++){
	nusq += w[i][j]*w[i][j];
      }
      wvar[i] = nusq/m_j;
    }
  }else if(bias == 1){
    for(int i = 0; i < j_max; i++){
      nusq = .0;
      for(int j = 0; j < n; j++){
	nusq += w[i][j]*w[i][j];
      }
      wvar[i] = nusq/n;
    }
  }

  //Calculate confidence interval
  if(p == 90){
    phi = 1.6449;
  }else if(p == 95){
    phi = 1.96;
  }else if(p == 99){
    phi = 2.5758;
  }
  for(int i = 0; i < j_max; i++){
    eta = edof[i];
    q_lower = eta*pow((1-2./(9*eta)+phi*sqrt(2./(9*eta))),3);
    phi *= -1;
    q_upper = fabs(eta*pow((1-2./(9*eta)+phi*sqrt(2./(9*eta))),3));
    phi *= -1;
    ii = j_max+i;
    ci[i] = log(eta*wvar[i]/q_lower);
    ci[ii] = log(eta*wvar[i]/q_upper);
  }
}

void modwt::edof(int bias, int n, int j_max, double *eta){
  int j, l_j, m_j;

  if(bias == 0){
    for(int i = 0; i < j_max; i++){
      j = i;
      l_j = get_lj(j);
      m_j = n-l_j+1;
      j++;
      eta[i] = m_j/pow(2,j);
      if(eta[i] < 1){
	eta[i] = 1;
      }
    }
  }else if(bias == 1){
    for(int i = 0; i < j_max; i++){
      j = i+1;
      eta[i] = n/pow(2,j);
      if(eta[i] < 1){
	eta[i] = 1;
      }
    }
  }
}

int modwt::get_lj(int j){//Fjerne disse
  j++; //I will call when j = 0,1,2... But that means j = 1,2,3 etc
  return ((pow(2,j)-1)*(l-1)+1);
}

int modwt::g_shift(int j){
  int l_j = modwt::get_lj(j);
  int nu_j = nu*(l_j-1)/(l-1);

  return nu_j;
}

int modwt::h_shift(int j){
  int l_j = modwt::get_lj(j);
  int nu_j = -(l_j/2+l/2+nu-1);

  return nu_j;
}

void modwt::shift(int j_max, int n, double **wv){
  int shift;
  double *tmp = new double[n];
  
  for(int j = 0; j < j_max; j++){
    shift = fabs(modwt::h_shift(j)); //wavelet
    for(int i = shift; i < n; i++){//Turn i and j around to make it more clear...
      tmp[i-shift] = wv[j][i];
    }
    for(int i = 0; i < shift; i++){
      tmp[n-shift+i] = wv[j][i];
    }
    for(int i = 0; i < n; i++){
      wv[j][i] = tmp[i];
    }
  }  
  shift = fabs(modwt::g_shift(j_max-1));//scaling
  for(int i = shift; i < n; i++){//Turn i and j around to make it more clear...
    tmp[i-shift] = wv[j_max][i];
  }
  for(int i = 0; i < shift; i++){
    tmp[n-shift+i] = wv[j_max][i];
  }
  for(int i = 0; i < n; i++){
    wv[j_max][i] = tmp[i];
  }

  delete[] tmp;
}

void modwt::get_boundaries(int param, int j_max, int n, double *boundaries){
  int l_j, nu;

  for(int i = 0; i <= 2*j_max+1; i++){//initialize boundaries
    boundaries[i] = -1;
  }
  
  if(param == 0){
    for(int i = 0; i < j_max; i++){
      l_j = modwt::get_lj(i);
      boundaries[i] = l_j-2;
    }
    boundaries[2*j_max] = l_j-2;
  }else if(param == 1){
    for(int i = 0; i < j_max; i++){
      l_j = modwt::get_lj(i);
      nu = modwt::h_shift(i);
      boundaries[i] = l_j-2-fabs(nu);
      boundaries[i+j_max] = n-fabs(nu);
    }
    nu = modwt::g_shift(j_max-1);
    boundaries[2*j_max] = l_j-2-fabs(nu);
    boundaries[2*j_max+1] = n-fabs(nu);
  }else if(param == 2){
    for(int i = 0; i < j_max; i++){
      l_j = modwt::get_lj(i);
      boundaries[i] = l_j-2;
      boundaries[i+j_max] = n-l_j+1;
    }
    boundaries[2*j_max] = l_j-2;
    boundaries[2*j_max+1] = n-l_j+1;
  }
}

void toolbox::reflection(int n, double *&x){
  double *tmp = new double[2*n];

  for(int i = 0; i < n; i++){
    tmp[i] = x[i];
    tmp[i+n] = x[n-i-1];
  }

  delete[] x;
  x = tmp;
}

void toolbox::padding(int padding,int n, int &m, double *&x){
  int j = 0;
  double value;
  double *tmp;

  if(padding == 0){//pad with zeros
    value = 0;
  }else if(padding == 1){//pad with sample mean
    value = .0;
    for(int i = 0; i < n; i++){
      value += x[i];
    }
    value /= n;
  }

  while(pow(2,j) < n){
    j++;
  }
  m = pow(2,j);
  if(m > n){
    tmp = new double[m];
    for(int i = 0; i < n; i++){
      tmp[i] = x[i];
    }
    for(int i = n; i < m; i++){
      tmp[i] = value;
    }
    delete[] x;
    x = tmp;
  }
}

void toolbox::truncate(int n, double *&x){
  double *tmp = new double[n];
  
  for(int i = 0; i < n; i++){
    tmp[i] = x[i];
  }
  delete[] x;
  x = tmp;
}

void toolbox::wlse(int j_0, int j_min, int j_max, double *edof, double *wvar, double *y, double *y_hat, double *output){
  double arg, beta, eta, zeta, var;
  double *tau = new double[j_0];
  double *w = new double[j_0];

  for(int i = 0; i < j_0; i++){
    arg = edof[i]/2;
    w[i] = 1./trigamma(arg);
    y[i] = log(wvar[i])-digamma(arg)+log(arg);
    tau[i] = pow(2,i);
  }

  double a = .0;
  double b = .0;
  double c = .0;
  double d = .0;
  double e = .0;
  double f = .0;
  
  j_min--;
  for(int i = j_min; i < j_max; i++){
    a += w[i];
    b += w[i]*log(tau[i])*y[i];
    c += w[i]*log(tau[i]);
    d += w[i]*y[i];
    e += w[i]*pow(log(tau[i]),2);
    f += w[i]*log(tau[i]);
  }
  
  beta = (a*b-c*d)/(a*e-f*f);
  zeta = (e*d-c*b)/(a*e-f*f);
  var = a/(a*e-f*f);
  output[0] = beta/2;
  output[1] = var;
  output[2] = sqrt(.25*var);
  
  for(int i = j_min; i < j_max; i++){
    y_hat[i] = zeta+beta*log(tau[i]);
    //std::cout<<"i: "<<i+1<<"; "<<y_hat[i]<<std::endl;
  }
  
  delete[] tau;
  delete[] w;
}

void toolbox::iwlse(int j_min, int j_max, int n, double *hurst, double *output, double **w){
  int jj = j_max-j_min+1;                                      
  double eul = 0.5772156649015328606; //Euler-Mascheroni constant
  double pi = 3.14159265358979323846; //The well known PI
  double a, b, c, d, var;
  double avgh = .0;
  double *tau = new double[j_max];
  double *y = new double[j_max]; 
  double test;
  j_min--;
  
  for(int i = 0; i < n; i++){//time //Should not allow to sum over boundary coefficients! p109 solution guide  
    for(int j = j_min; j < j_max; j++){//scale
      tau[j] = pow(2,j);
      y[j] = log(pow(w[j][i],2))+log(2)+eul; 
    }
    a = .0;
    b = .0;
    c = .0;
    d = .0;
    for(int j = j_min; j < j_max; j++){
      a += log(tau[j])*y[j];
      b += log(tau[j]);
      c += y[j];
      d += pow(log(tau[j]),2); 
    }
	   
    hurst[i] = .5*(jj*a-b*c)/(jj*d-b*b);
    avgh += hurst[i];
  }
  var = (jj*pow(pi,2))/(8*(jj*d-b*b));
  output[0] = avgh/n;
  output[1] = var;
  output[2] = sqrt(.25*var);

  delete[] tau;
  delete[] y;
}

double toolbox::digamma(double x){ //http://people.sc.fsu.edu/~burkardt/cpp_src/prob/prob.C
  double s3 = 0.08333333333; //Must be checked for accuracy and speed
  double s4 = 0.0083333333333;
  double s5 = 0.003968253968;
  double value = .0;
  double y, z;

  z = x;
  while(z < 10){
    value -= 1./z;
    z += 1.;
  }
  y = 1./pow(z,2);
  value += log(z)-.5/z - y*(s3-y*(s4-y*s5));

  return value;
}

double toolbox::trigamma(double x){ //Must be checked for accuracy
  double b2 = 1/6;
  double b4 = -1/30;
  double b6 = 1/42;
  double b8 = -1/30;
  double b10 = 5/66;
  double value = .0;
  double y, z;

  z = x;
  while(z < 10000){//Works okay with z < 10
    value += 1./pow(z,2);
    z += 1.;
  }
  y = 1./pow(z,2);
  value += .5*y+(1.+y*(b2+y*(b4+y*(b6+y*(b8+y*b10)))))/z;
  
  return value;
}

void toolbox::ahurst(int transtype, int l, int type, int bias, int refl, int n, int j_0, int j_min, int j_max, int p,
		     double *x_in, double *y, double *y_hat, double *ci, double *output){
  int pad = 0;
  double *edof = new double[j_0];
  double *wvar = new double[j_0];
  double **wv = new double *[j_0+1];

  if(transtype == 0){
    int m;//m is length of padded array (m = n if array not padded)
    dwt object(l, type);
    if(refl == 0){
      toolbox::padding(pad, n, m, x_in);
      for(int j = 0; j <= j_0; j++){ //Har nå satt nh i ptrans til å være n 
	wv[j] = new double[m/2];//Hei!!!!!! Kan jo gjøre hver av disse kortere for hver iterasjon!!
      }
      object.partial_transform(m, j_0, x_in, wv);//Sjekk kommentar 3 side 367!!!
      object.edof(bias, n, j_0, edof);
      object.wavevar(bias, n, j_0, p, edof, wvar, ci, wv);
      toolbox::wlse(j_0, j_min, j_max, edof, wvar, y, y_hat, output);
      toolbox::truncate(n, x_in);
    }else if(refl == 1){
      toolbox::padding(pad, n, m, x_in);
      toolbox::reflection(m, x_in);                   
      for(int j = 0; j <= j_0; j++){
	wv[j] = new double[2*m];
      }
      object.partial_transform(2*m, j_0, x_in, wv);
      object.edof(bias, n, j_0, edof);
      object.wavevar(bias, 2*n, j_0, p, edof, wvar, ci, wv);
      toolbox::wlse(j_0, j_min, j_max, edof, wvar, y, y_hat, output);
      toolbox::truncate(n, x_in);
    }
  }else if(transtype == 1){
    modwt object(l, type);
    if(refl == 0){
      for(int j = 0; j <= j_0; j++){
	wv[j] = new double[n];
      }
      object.partial_transform(n, j_0, x_in, wv);
      //object.shift(j_0, n, wv);   //REMOVE
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 55; j++){
	std::cout<<wv[i][j]<<std::endl;
      }
      std::cout<<"----------"<<std::endl;
    }
      object.edof(bias, n, j_0, edof);
      object.wavevar(bias, n, j_0, p, edof, wvar, ci, wv);
      toolbox::wlse(j_0, j_min, j_max, edof, wvar, y, y_hat, output);
    }else if(refl == 1){
      for(int j = 0; j <= j_0; j++){
	wv[j] = new double[2*n];
      }
      toolbox::reflection(n, x_in);
      object.partial_transform(2*n, j_0, x_in, wv);
      object.edof(bias, n, j_0, edof);
      object.wavevar(bias, 2*n, j_0, p, edof, wvar, ci, wv);
      toolbox::wlse(j_0, j_min, j_max, edof, wvar, y, y_hat, output);
      toolbox::truncate(n, x_in);
    }
  }

  for(int j = 0; j <= j_max; j++){
    delete[] wv[j];
  }
  delete[] wv;
  delete[] wvar;
  delete[] edof;
}

void toolbox::ihurst(int l, int type, int refl, int j_0, int j_min, int j_max, int n, 
		     double *x_in, double *hurst, double *output){
  modwt object(l, type);
  double **ds = new double *[j_0+1];
  double **wv = new double *[j_0+1];

  if(refl == 0){
    for(int j = 0; j <= j_0; j++){
      wv[j] = new double[n];
    }
    object.partial_transform(n, j_0, x_in, wv);
    object.shift(j_0, n, wv);                      
    toolbox::iwlse(j_min, j_max, n, hurst, output, wv);//Add j_0 ass parameter, error with j_min in this function...
  }else if(refl == 1){
    for(int j = 0; j <= j_max; j++){
      wv[j] = new double[2*n];
    }
    toolbox::reflection(n, x_in);
    object.partial_transform(2*n, j_0, x_in, wv);
    for(int i = 0; i <= j_max; i++){
      toolbox::truncate(n, wv[i]);
    }
    object.shift(j_0, n, wv);
    toolbox::iwlse(j_min, j_max, n, hurst, output, wv);
  }
  
  //Smooooooth mr Hurst
  std::cout<<"J0: "<<j_0<<std::endl;
  std::cout<<"JMAX: "<<j_max<<std::endl;
  for(int j = 0; j <= j_max; j++){
    ds[j] = new double[n];
  }
  object.partial_transform(n, j_max, hurst, wv);
  object.mra(n, j_max, wv, ds);
  for(int i = 0; i < n; i++){
    hurst[i] = ds[j_max][i]; 
  }
  
  for(int j = 0; j <= j_max; j++){
    delete[] wv[j];
    delete[] ds[j];
  }
  delete[] ds;
  delete[] wv;
}

void toolbox::mra(int transtype, int l, int type, int refl, int j_max, int n, 
		  double *boundaries, double *x_in, double **ds){//Check boundaries
  int pad = 1;
  int mra = 2;
  double **wv = new double *[j_max+1];

  if(transtype == 0){
    int m;
    dwt object(l, type);

    if(refl == 0){
      toolbox::padding(pad, n, m, x_in);
      for(int j = 0; j <= j_max; j++){
	wv[j] = new double[m]; //Can adjust this to get correct length of each level!
      }
      object.partial_transform(m, j_max, x_in, wv);
      object.mra(m, j_max, wv, ds);
      object.get_boundaries(mra, j_max, m, boundaries);
      toolbox::truncate(n, x_in);
    }else if(refl == 1){
      toolbox::padding(pad, n, m, x_in);
      toolbox::reflection(m, x_in);
      for(int j = 0; j <= j_max; j++){
	wv[j] = new double[2*m];
      }
      object.partial_transform(2*m, j_max, x_in, wv);
      object.mra(2*m, j_max, wv, ds);
      object.get_boundaries(mra, j_max, m, boundaries);
      for(int j = 0; j <= j_max; j++){
	toolbox::truncate(m, ds[j]);
      }
      toolbox::truncate(n, x_in);
    }
  }else if(transtype == 1){
    modwt object(l, type);

    if(refl == 0){
      for(int j = 0; j <= j_max; j++){
	wv[j] = new double[n];
      }
      object.partial_transform(n, j_max, x_in, wv);
      object.mra(n, j_max, wv, ds);
      object.get_boundaries(mra, j_max, n, boundaries);//Check that boundaries are ok!
    }else if(refl == 1){
      for(int j = 0; j <= j_max; j++){
	wv[j] = new double[2*n];//Not sure?
      }
      toolbox::reflection(n, x_in);
      object.partial_transform(2*n, j_max, x_in, wv);
      object.mra(2*n, j_max, wv, ds);
      object.get_boundaries(mra, j_max, n, boundaries);
      for(int j = 0; j <= j_max; j++){
	toolbox::truncate(n, ds[j]);
      }
      toolbox::truncate(n, x_in);
    }
  }
  /*
  for(int i = 0; i < j_max; i++){
    for(int j = 0; j < 6; j++){
      std::cout<<ds[i][j]<<std::endl;
    }
    std::cout<<"----------"<<std::endl;
  }
  for(int i = 0; i < 6; i++){
    std::cout<<ds[j_max][i]<<std::endl;
  }
  */
  for(int j = 0; j <= j_max; j++){
    delete[] wv[j];
  }
  delete[] wv;
}

void toolbox::wa(int transtype, int l, int type, int refl, int j_max, int n, 
		 double *boundaries, double *x_in, double **wv){//Needs testing for reflecting boundaries!!!
  int shifted = 0;
  int pad = 1;
  
  if(transtype == 0){
    int m;
    dwt object(l, type);
    
    if(refl == 0){
    toolbox::padding(pad, n, m, x_in);
    object.partial_transform(m, j_max, x_in, wv);
    if(type != 1){
      shifted = 1;
      object.shift(j_max, m, wv);
    }
    object.get_boundaries(shifted, j_max, n, boundaries);
    toolbox::truncate(n, x_in);
    }else if(refl == 1){
      toolbox::padding(pad, n, m, x_in);
      toolbox::reflection(m, x_in);
      object.partial_transform(2*m, j_max, x_in, wv);
      for(int i = 0; i < j_max; i++){ //Ikke truncate på V??????????????
	toolbox::truncate(m/pow(2,i), wv[i]);
      }
      toolbox::truncate(n, wv[j_max]);
      if(type != 1){
	object.shift(j_max, m, wv);
      }
      object.get_boundaries(shifted, j_max, n, boundaries);//Should use m (and change to n_j in the function
      toolbox::truncate(n, x_in);                          //DWT is allways of power 2**n at this stage
    }
  }else if(transtype == 1){
    modwt object(l, type);

    if(refl == 0){
      object.partial_transform(n, j_max, x_in, wv);
      if(type != 1){
	shifted = 1;
	object.shift(j_max, n, wv);
      }
      object.get_boundaries(shifted, j_max, n, boundaries);
    }else if(refl == 1){
      toolbox::reflection(n, x_in);
      object.partial_transform(2*n, j_max, x_in, wv); //Make sure that input length of matrix is 2n when refl == 1
      for(int i = 0; i <= j_max; i++){
	toolbox::truncate(n, wv[i]);
      }
      if(type != 1){
	object.shift(j_max, n, wv);
      }
      object.get_boundaries(shifted, j_max, n, boundaries);
      toolbox::truncate(n, x_in);
    }
  }
}

void toolbox::shurst(int l, int type, int refl, int j_0, int j_min, int j_max, int n, int n_b, 
		     double *x_in){//, double *output){
  modwt object(l, type);
  int n_s = floor(n/n_b);
  int start = 0;
  int end = n_s;
  int bias = 1;
  int p = 95;
  double avg = .0;
  double nusq = .0;
  double *ci = new double[2*j_0+2];
  double *edof = new double[j_0];
  double *wvar = new double[j_0];
  double **wv = new double *[j_0+1];
  double *tmp = new double[3];
  double *useless = new double[j_0];
  double *useless2 = new double[j_0];

  for(int i = 0; i <= j_0; i++){
    wv[i] = new double[n];
  }

  object.partial_transform(n, j_0, x_in, wv);
  object.edof(bias, n_s, j_0, edof);
  object.shift(j_0, n, wv);
  for(int i = 0; i < n_b; i++){
    for(int j = 0; j < j_0; j++){
      for(int k = start; k < end; k++){
      	nusq += wv[j][k]*wv[j][k];
      }
      nusq /= n_s;
      wvar[j] = nusq;
    }
    toolbox::wlse(j_0, j_min, j_max, edof, wvar, useless, useless2, tmp);
    start += n_s;
    end += n_s;
    nusq = .0;
    avg += tmp[0];
    std::cout<<tmp[0]<<std::endl;
  }
  avg /= n_b;
  std::cout<<"Avg: "<<avg<<std::endl;
  
  delete[] useless2;
  delete[] useless;
  delete[] ci;
  delete[] edof;
  delete[] wvar;
  delete[] wv;
  delete[] tmp;
}

int get_nu(int l, int type){
  if(type == 2){         //LA
    if(l == 14){
      return -l/2+2;
    }else if(l == 10 || l == 18){
      return -l/2;
    }else{
      return -l/2+1;
    }
  }else if(type == 3){   //BL
    if(l == 14){
      return -5;
    }else if(l == 18){
      return -11;
    }else if(l == 20){
      return -9;
    }
  }else if(type == 4){   //Coiflet
    return -2*l/3+1;
  }
}

double get_daub(int l, int i){ //Er dette dumt?
  const double d2[2] = {            
    0.7071067811865475,            
    0.7071067811865475
  };
  const double d4[4] = {
    0.4829629131445341,
    0.8365163037378077,
    0.2241438680420134,
    -0.1294095225512603
  };
  const double d6[6] = {
    0.3326705529500827,
    0.8068915093110928,
    0.4598775021184915,
    -0.1350110200102546,
    -0.0854412738820267,
    0.0352262918857096,
  };
  const double d8[8] = {
    0.2303778133074431,
    0.7148465705484058,
    0.6308807679358788,
    -0.0279837694166834,
    -0.1870348117179132,
    0.0308413818353661,
    0.0328830116666778,
    -0.0105974017850021
  };
  const double d10[10] = {
    0.1601023979741930,
    0.6038292697971898,
    0.7243085284377729,
    0.1384281459013204,
    -0.2422948870663824,
    -0.0322448695846381,
    0.0775714938400459,
    -0.0062414902127983,
    -0.0125807519990820,
    0.0033357252854738
  };
  const double d12[12] = {
    0.1115407433501094,
    0.4946238903984530,
    0.7511339080210954,
    0.3152503517091980,
    -0.2262646939654399,
    -0.1297668675672624,
    0.0975016055873224,
    0.0275228655303053,
    -0.0315820393174862,
    0.0005538422011614,
    0.0047772575109455,
    -0.0010773010853085
  };
  const double d14[14] = {
    0.0778520540850081,
    0.3965393194819136,
    0.7291320908462368,
    0.4697822874052154,
    -0.1439060039285293,
    -0.2240361849938538,
    0.0713092192668312,
    0.0806126091510820,
    -0.0380299369350125,
    -0.0165745416306664,
    0.0125509985560993,
    0.0004295779729214,
    -0.0018016407040474,
    0.0003537137999745
  };
  const double d16[16] = {
    0.0544158422431049,
    0.3128715909143031,
    0.6756307362972904,
    0.5853546836541907,
    -0.0158291052563816,
    -0.2840155429615702,
    0.0004724845739124,
    0.1287474266204837,
    -0.0173693010018083,
    -0.0440882539307952,
    0.0139810279173995,
    0.0087460940474061,
    -0.0048703529934518,
    -0.0003917403733770,
    0.0006754494064506,
    -0.0001174767841248
  };
  const double d18[18] = {
    0.0380779473638791,
    0.2438346746125939,
    0.6048231236901156,
    0.6572880780512955,
    0.1331973858249927,
    -0.2932737832791761,
    -0.0968407832229524,
    0.1485407493381306,
    0.0307256814793395,
    -0.0676328290613302,
    0.0002509471148340,
    0.0223616621236805,
    -0.0047232047577520,
    -0.0042815036824636,
    0.0018476468830564,
    0.0002303857635232,
    -0.0002519631889427,
    0.0000393473203163
  };
  const double d20[20] = {
    0.0266700579005546,
    0.1881768000776863,
    0.5272011889317202,
    0.6884590394536250,
    0.2811723436606485,
    -0.2498464243272283,
    -0.1959462743773399,
    0.1273693403357890,
    0.0930573646035802,
    -0.0713941471663697,
    -0.0294575368218480,
    0.0332126740593703,
    0.0036065535669880,
    -0.0107331754833036,
    0.0013953517470692,
    0.0019924052951930,
    -0.0006858566949566,
    -0.0001164668551285,
    0.0000935886703202,
    -0.0000132642028945
    };
  
  if(l == 2){
    return d2[i];
  }else if(l == 4){
    return d4[i];
  }else if(l == 6){
    return d6[i];
  }else if(l == 8){
    return d8[i];
  }else if(l == 10){
    return d10[i];
  }else if(l == 12){
    return d12[i];
  }else if(l == 14){
    return d14[i];
  }else if(l == 16){
    return d16[i];
  }else if(l == 18){
    return d18[i];
  }else if(l == 20){
    return d20[i];
  }
}

double get_la(int l, int i){
  const double la8[8] = {
    -0.0757657147893407,
    -0.0296355276459541,
    0.4976186676324578,
    0.8037387518052163,
    0.2978577956055422,
    -0.0992195435769354,
    -0.0126039672622612,
    0.0322231006040713
  };
  const double la10[10] = {
    0.0195388827353869,
    -0.0211018340249298,
    -0.1753280899081075,
    0.0166021057644243,
    0.6339789634569490,
    0.7234076904038076,
    0.1993975339769955,
    -0.0391342493025834,
    0.0295194909260734,
    0.0273330683451645
  };
  const double la12[12] = {
    0.0154041093273377,
    0.0034907120843304,
    -0.1179901111484105,
    -0.0483117425859981,
    0.4910559419276396,
    0.7876411410287941,
    0.3379294217282401,
    -0.0726375227866000,
    -0.0210602925126954,
    0.0447249017707482,
    0.0017677118643983,
    -0.0078007083247650
  };
  const double la14[14] = {
    0.0102681767084968,
    0.0040102448717033,
    -0.1078082377036168,
    -0.1400472404427030,
    0.2886296317509833,
    0.7677643170045710,
    0.5361019170907720,
    0.0174412550871099,
    -0.0495528349370410,
    0.0678926935015971,
    0.0305155131659062,
    -0.0126363034031526,
    -0.0010473848889657,
    0.0026818145681164
  };
  const double la16[16] = {
    -0.0033824159513594,
    -0.0005421323316355,
    0.0316950878103452,
    0.0076074873252848,
    -0.1432942383510542,
    -0.0612733590679088,
    0.4813596512592012,
    0.7771857516997478,
    0.3644418948359564,
    -0.0519458381078751,
    -0.0272190299168137,
    0.0491371796734768,
    0.0038087520140601,
    -0.0149522583367926,
    -0.0003029205145516,
    0.0018899503329007
  };
  const double la18[18] = {
    0.0010694900326538,
    -0.0004731544985879,
    -0.0102640640276849,
    0.0088592674935117,
    0.0620777893027638,
    -0.0182337707798257,
    -0.1915508312964873,
    0.0352724880359345,
    0.6173384491413523,
    0.7178970827642257,
    0.2387609146074182,
    -0.0545689584305765,
    0.0005834627463312,
    0.0302248788579895,
    -0.0115282102079848,
    -0.0132719677815332,
    0.0006197808890549,
    0.0014009155255716
  };
  const double la20[20] = {
    0.0007701598091030,
    0.0000956326707837,
    -0.0086412992759401,
    -0.0014653825833465,
    0.0459272392237649,
    0.0116098939129724,
    -0.1594942788575307,
    -0.0708805358108615,
    0.4716906668426588,
    0.7695100370143388,
    0.3838267612253823,
    -0.0355367403054689,
    -0.0319900568281631,
    0.0499949720791560,
    0.0057649120455518,
    -0.0203549398039460,
    -0.0008043589345370,
    0.0045931735836703,
    0.0000570360843390,
    -0.0004593294205481
  };

  if(l == 8){
    return la8[i];
  }else if(l == 10){
    return la10[i];
  }else if(l == 12){
    return la12[i];
  }else if(l == 14){
    return la14[i];
  }else if(l == 16){
    return la16[i];
  }else if(l == 18){
    return la18[i];
  }else if(l == 20){
    return la20[i];
  }
}

double get_bl(int l, int i){
    const double bl8[8] = {
    -0.0757657147893407,
    -0.0296355276459541,
    0.4976186676324578,
    0.8037387518052163,
    0.2978577956055422,
    -0.0992195435769354,
    -0.0126039672622612,
    0.0322231006040713
  };
  const double bl10[10] = {
    0.0195388827353869,
    -0.0211018340249298,
    -0.1753280899081075,
    0.0166021057644243,
    0.6339789634569490,
    0.7234076904038076,
    0.1993975339769955,
    -0.0391342493025834,
    0.0295194909260734,
    0.0273330683451645
  };
  const double bl12[12] = {
    0.0154041093273377,
    0.0034907120843304,
    -0.1179901111484105,
    -0.0483117425859981,
    0.4910559419276396,
    0.7876411410287941,
    0.3379294217282401,
    -0.0726375227866000,
    -0.0210602925126954,
    0.0447249017707482,
    0.0017677118643983,
    -0.0078007083247650
  };
  const double bl14[14] = {
    0.0120154192834842,
    0.0172133762994439,
    -0.0649080035533744,
    -0.0641312898189170,
    0.3602184608985549,
    0.7819215932965554,
    0.4836109156937821,
    -0.0568044768822707,
    -0.1010109208664125,
    0.0447423494687405,
    0.0204642075778225,
    -0.0181266051311065,
    -0.0032832978473081,
    0.0022918339541009
  };
  const double bl16[16] = {
    -0.0033824159513594,
    -0.0005421323316355,
    0.0316950878103452,
    0.0076074873252848,
    -0.1432942383510542,
    -0.0612733590679088,
    0.4813596512592012,
    0.7771857516997478,
    0.3644418948359564,
    -0.0519458381078751,
    -0.0272190299168137,
    0.0491371796734768,
    0.0038087520140601,
    -0.0149522583367926,
    -0.0003029205145516,
    0.0018899503329007
  };
  const double bl18[18] = {
    0.0002594576266544,
    -0.0006273974067728,
    -0.0019161070047557,
    0.0059845525181721,
    0.0040676562965785,
    -0.0295361433733604,
    -0.0002189514157348,
    0.0856124017265279,
    -0.0211480310688774,
    -0.1432929759396520,
    0.2337782900224977,
    0.7374707619933686,
    0.5926551374433956,
    0.0805670008868546,
    -0.1143343069619310,
    -0.0348460237698368,
    0.0139636362487191,
    0.0057746045512475
  };
  const double bl20[20] = {
    0.0008625782242896,
    0.0007154205305517,
    -0.0070567640909701,
    0.0005956827305406,
    0.0496861265075979,
    0.0262403647054251,
    -0.1215521061578162,
    -0.0150192395413644,
    0.5137098728334054,
    0.7669548365010849,
    0.3402160135110789,
    -0.0878787107378667,
    -0.0670899071680668,
    0.0338423550064691,
    -0.0008687519578684,
    -0.0230054612862905,
    -0.0011404297773324,
    0.0050716491945793,
    0.0003401492622332,
    -0.0004101159165852
  };
  if(l == 8){
    return bl8[i];
  }else if(l == 10){
    return bl10[i];
  }else if(l == 12){
    return bl12[i];
  }else if(l == 14){
    return bl14[i];
  }else if(l == 16){
    return bl16[i];
  }else if(l == 18){
    return bl18[i];
  }else if(l == 20){
    return bl20[i];
  }
}

double get_coiflet(int l, int i){
  const double c6[6] = {
    -0.0156557285289848,
    -0.0727326213410511,
    0.3848648565381134,
    0.8525720416423900,
    0.3378976709511590,
    -0.0727322757411889
  };
  const double c12[12] = {
    -0.0007205494453679,
    -0.0018232088707116,
    0.0056114348194211,
    0.0236801719464464,
    -0.0594344186467388,
    -0.0764885990786692,
    0.4170051844236707,
    0.8127236354493977,
    0.3861100668229939,
    -0.0673725547222826,
    -0.0414649367819558,
    0.0163873364635998
  };
  const double c18[18] = {
    -0.0000345997728362,
    -0.0000709833031381,
    0.0004662169601129,
    0.0011175187708906,
    -0.0025745176887502,
    -0.0090079761366615,
    0.0158805448636158,
    0.0345550275730615,
    -0.0823019271068856,
    -0.0717998216193117,
    0.4284834763776168,
    0.7937772226256169,
    0.4051769024096150,
    -0.0611233900026726,
    -0.0657719112818552,
    0.0234526961418362,
    0.0077825964273254,
    -0.0037935128644910
  };
  const double c24[24] = {
    -0.0000017849850031,
    -0.0000032596802369,
    0.0000312298758654,
    0.0000623390344610,
    -0.0002599745524878,
    -0.0005890207562444,
    0.0012665619292991,
    0.0037514361572790,
    -0.0056582866866115,
    -0.0152117315279485,
    0.0250822618448678,
    0.0393344271233433,
    -0.0962204420340021,
    -0.0666274742634348,
    0.4343860564915321,
    0.7822389309206135,
    0.4153084070304910,
    -0.0560773133167630,
    -0.0812666996808907,
    0.0266823001560570,
    0.0160689439647787,
    -0.0073461663276432,
    -0.0016294920126020,
    0.0008923136685824
  };
  const double c30[30] = {
    -0.0000000951765727,
    -0.0000001674428858,
    0.0000020637618516,
    0.0000037346551755,
    -0.0000213150268122,
    -0.0000413404322768,
    0.0001405411497166,
    0.0003022595818445,
    -0.0006381313431115,
    -0.0016628637021860,
    0.0024333732129107,
    0.0067641854487565,
    -0.0091642311634348,
    -0.0197617789446276,
    0.0326835742705106,
    0.0412892087544753,
    -0.1055742087143175,
    -0.0620359639693546,
    0.4379916262173834,
    0.7742896037334738,
    0.4215662067346898,
    -0.0520431631816557,
    -0.0919200105692549,
    0.0281680289738655,
    0.0234081567882734,
    -0.0101311175209033,
    -0.0041593587818186,
    0.0021782363583355,
    0.0003585896879330,
    -0.0002120808398259
  };

  if(l == 6){
    return c6[i];
  }else if(l == 12){
    return c12[i];
  }else if(l == 18){
    return c18[i];
  }else if(l == 24){
    return c24[i];
  }else if(l == 30){
    return c30[i];
  }
}


//Boundaries: l_j' = gamma_j + gamma_j*              Tables 137
//D/S : x = 2**j*l_j' (0,...,x V N-l_j'+2**j,..,N-1) (Eq 139) 

/*  
//WA printout with boundaries

  for(int i = 0; i < j_max; i++){
    //std::cout<<"W_"<<i+1<<":"<<std::endl;
    std::cout<<"Boundaries: "<<boundaries[i]<<" , "<<boundaries[i+j_max]<<std::endl;
    for(int j = 0; j < 6; j++){
      std::cout<<w[i][j]<<std::endl;
    }
    std::cout<<"------------"<<std::endl;
  }
  std::cout<<"V_"<<j_max<<":"<<std::endl;
  std::cout<<"Boundaries: "<<boundaries[2*j_max]<<" , "<<boundaries[2*j_max+1]<<std::endl;
  for(int i = 0; i < 6; i++){
    std::cout<<x_in[i]<<std::endl;
  }

//AHURST options printout
  if(transtype == 0){
    std::cout<<"DWT"<<std::endl;
  }else if(transtype == 1){
    std::cout<<"MODWT"<<std::endl;
  }
  if(refl == 0){
    std::cout<<"Circular"<<std::endl;
  }else if(refl == 1){
    std::cout<<"Reflection"<<std::endl;
  }
  if(bias == 0){
    std::cout<<"Unbiased"<<std::endl;
  }else if(bias == 1){
    std::cout<<"Biased"<<std::endl;
  }
  std::cout<<"j_min: "<<j_min<<" j_max: "<<j_max<<std::endl;

*/

  /*
  i = 0;
  for(int j = j_min; j < j_max; i++){ //Ok if y_hat is included. Should
    a += w[i];
    b += w[i]*log(tau[i])*y[i];
    c += w[i]*log(tau[i]);
    d += w[i]*y[i];
    e += w[i]*pow(log(tau[i]),2);
    f += w[i]*log(tau[i]);
    i++;
  }
  */
  /*
  i = 0;
  for(int j = j_min; j < j_max; j++){ //Rewrite to get y in all points...
    arg = edof[j]/2;
    w[i] = 1./trigamma(arg);
    y[i] = log(wvar[j])-digamma(arg)+log(arg);
    tau[i] = pow(2,j);
    std::cout<<y[i]<<std::endl;
    i++;
  }
  */

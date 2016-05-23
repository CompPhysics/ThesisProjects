class dwt{                  //Class names should start with capital letter
 private:
  int l, type, nu;          
  double *g, *h;

  void transform(int, double *, double *, double *);
  void inverse_transform(int, double *, double *, double *);
  void init();
  int get_lj(int); //Not needed?
  int g_shift(int); 
  int h_shift(int);	
  
 public:
  dwt(int, int);
  ~dwt();
  void partial_transform(int, int, double *, double **); //rename trans, itrans and ptrans
  void mra(int, int, double **, double **);
  void wavevar(int, int, int, int, double *, double *, double *, double **);
  void edof(int, int, int, double *);
  void shift(int, int, double **);
  void get_boundaries(int, int, int, double *);
  friend int get_nu(int, int);
  friend double get_daub(int, int);
  friend double get_la(int, int);
  friend double get_bl(int, int);
  friend double get_coiflet(int, int);
};

class modwt{
 private:
  int l, type, nu;
  double *g, *h;

  int get_lj(int); //Not needed?
  int g_shift(int);
  int h_shift(int);
  void transform(int, int, double *, double *, double *);
  void inverse_transform(int, int, double *, double *, double *);
  void init();
  
 public:
  modwt(int, int);
  ~modwt();
  void partial_transform(int, int, double *, double **);
  void mra(int, int, double **, double **);
  void wavevar(int, int, int, int, double *, double *, double *, double **);
  void edof(int, int, int, double *);
  void shift(int, int, double **);
  void get_boundaries(int, int, int, double *);
  friend int get_nu(int, int);
  friend double get_daub(int, int);
  friend double get_la(int, int);
  friend double get_bl(int, int);
  friend double get_coiflet(int, int);
};

class toolbox{
 private:
  void reflection(int, double *&);
  void padding(int, int, int&, double *&);
  void truncate(int, double *&);
  void wlse(int, int, int, double *, double *, double *, double *, double *);
  void iwlse(int, int, int, double *, double *, double **);
  double digamma(double);
  double trigamma(double);
 public:
  void ahurst(int, int, int, int, int, int, int, int, int, int, double *, double *, double *, double *, double *);
  void ihurst(int, int, int, int, int, int, int, double *, double *, double *);
  void shurst(int, int, int, int, int, int, int, int, double *);
  void mra(int, int, int, int, int, int, double *, double *, double **);
  void wa(int, int, int, int, int, int, double *, double *, double **);
};

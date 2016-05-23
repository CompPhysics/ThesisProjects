class dwt{                  //Class names should start with capital letter
 private:
  int w, family;
  double *g, *h;
  
  void init_coeff(int, int, double *, double *);
  friend double get_daub(int, int);
  friend double get_la(int, int);
  friend double get_bl(int, int);
  friend double get_coiflet(int, int);
  
 public:
  dwt(int, int);
  ~dwt();
  void transform(int, double *, double *, double *);
  void inverse_transform(int, double *, double *, double *);
};

class modwt{
 private:
  int w, family, it;
  double *g, *h;

  void init_coeff(int, int, double *, double *);
  friend double get_daub(int, int);
  friend double get_la(int, int);
  friend double get_bl(int, int);
  friend double get_coiflet(int, int);

 public:
  modwt(int, int);
  ~modwt();
  void transform(int, double *, double *, double *);
  void inverse_transform(int, double *, double *, double *);
};

class lifting{
 public:
  lifting();
  //split
  //predict
  //update
  //merge
};

class cwt{
 public:
  cwt();
};

void linear_regression(int, double *, double *);
double awc(int, double *);

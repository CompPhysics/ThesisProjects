#include "QickArray.hpp"

using namespace std;

QickArray::QickArray(void) {
#ifdef _DEBUG_
  cout << "entering QickArray()" << endl;
#endif
  data_ = NULL;
  for (int i=0; i<5; i++)
    dimensions_[i] = 0;
  n_dims_ = 0;
  is_initialized_ = _FALSE_;
#ifdef _DEBUG_
  cout << "exiting QickArray()" << endl;
#endif

}


QickArray::QickArray(int i1, int i2, int i3, int i4, int i5, int i6) {
#ifdef _DEBUG_
  cout << "entering QickArray(int i1, etc.)" << endl;
#endif

  data_ = NULL;
  is_initialized_ = 0; // we _know_ memory is not allocated right now.
  n_dims_ = 0;

  for (int i=0; i<6; i++) 
    dimensions_[0] = 0;

  int error = redim(i1, i2, i3, i4, i5, i6);

  if (!error)
    ;


#ifdef _DEBUG_
  if (error)
    cout << "Initialization did not succeed!" << endl;
  cout << "exiting QickArray(int i1, etc.)" << endl;
#endif

}


//
// this function allocates memory and shapes array.
// uesd in e.g. constructor.
//
// it is the ONLY function that makes a matrix ready for use.
//
int QickArray::redim(int i1, int i2, int i3, int i4, int i5, int i6) {

#ifdef _DEBUG_
  cout << "entering redim(int i1, etc.)" << endl;
#endif

  //delete old array if any
  delete_memory();

  //if an index is zero, then the dimension is ignored.

  n_dims_ = 6;

  if (i6==0) n_dims_--;
  if (i5==0) n_dims_--;
  if (i4==0) n_dims_--;
  if (i3==0) n_dims_--;
  if (i2==0) n_dims_--;
  if (i1==0) n_dims_--;

  //n_dims_ now contains correct number of dimensions.
  
  if (n_dims_ == 0) {
#ifdef _DEBUG_
    cout << "You just tried to create a zero-dimensional array!" << endl;
    cout << "exiting(int i1, etc.)" << endl;
#endif
    return 1;
    // here, is_initialized_ == 0, and memory is not allocated. ok.
  } 
  
  //we have a dimension not zero.
  
  //place dimensions.
  dimensions_[0] = i1; dimensions_[1] = i2; dimensions_[2] = i3;
  dimensions_[3] = i4; dimensions_[4] = i5; dimensions_[5] = i6;
  
  //check for negative dimensions. if found, bail out!
  
  for (int i=0; i<6; i++)
    if (dimensions_[i]<0) {
#ifdef _DEBUG_
      cout << "Found negative dimension! You bastard!" << endl;
      cout << "entering redim(int i1, etc.)" << endl;
#endif
      return 1;
    }
  
  //remove possible zeros from dimensions_[]
  //this got somewhat complicated, but it works...
  for (int i=0; i<6; i++) {
    int temp = 0;
    for (int j=i; j<6; j++)
      temp += dimensions_[j];
    if (temp)
      while (dimensions_[i]==0) {
	for (int j=i; j<5; j++)
	  dimensions_[j] = dimensions_[j+1];
	dimensions_[5] = 0;
      }
  }
  
#ifdef _DEBUG_
  
  cout << "dimensions_ after zero-removal:" << endl;
  for (int i=0; i<6; i++)
    cout << dimensions_[i] <<" ";
  cout << endl;
  
#endif
  
  //dimensions_ array now okay, even if the user was incredibly stupid.
  //n_dims_ also okay.


  //calculate total number of elements.

  no_of_elements_ = 1;

  for (int i=0; i<6; i++)
    if (dimensions_[i]>0)
      no_of_elements_ *= dimensions_[i];


  //store indexing-steps.
  step1_ = dimensions_[0];
  step2_ = step1_ * dimensions_[1];
  step3_ = step2_ * dimensions_[2];
  step4_ = step3_ * dimensions_[3];
  step5_ = step4_ * dimensions_[4];



  //allocate memory, then set is_initialized_. array ready for use!

  allocate_memory();  //allocates memory and sets is_initialized_ flag.

#ifdef _DEBUG_
  cout << "exiting redim(int i1, etc.)" << endl;
#endif

  return 0;
}




QickArray::~QickArray(void) {

#ifdef _DEBUG_
  cout << "entering ~QickArray()." << endl;
#endif

  //just erase memory.

  delete_memory();

#ifdef _DEBUG_
  cout << "exiting ~QickArray()." << endl;
#endif

}


//this little functions calculates the product of all dimensions, thus
//returning total number of elements in array.
int QickArray::no_of_elements(void) const {

  return no_of_elements_;

  // this was a real bottleneck! better to store in variable.

  /*
  int temp=1;

  for (int i=0; i<6; i++)
    if (dimensions_[i]>0)
      temp *= dimensions_[i];

  return temp;
  */
}


//simply allocate needed memory for data_.
int QickArray::allocate_memory(void) {
#ifdef _DEBUG_
  cout << "entering allocate_memory()." << endl;
#endif
  
  if (!is_initialized_) {
#ifdef _DEBUG_
    cout << "allocating..." << endl;
#endif
    data_ = new double[no_of_elements_]; //allocate
    is_initialized_ = _TRUE_;             //say that mem is allocated.
    fill_data(0);                         //fill
  }
  else {
#ifdef _DEBUG_
    cout << "exiting allocate_memory() with error-status." << endl;
#endif
    return 1;
  }
#ifdef _DEBUG_
    cout << "exiting allocate_memory()" << endl;
#endif

  return 0;
}


//simply allocate needed memory for data_.
int QickArray::delete_memory(void) {
#ifdef _DEBUG_
  cout << "entering delete_memory()." << endl;
#endif


  //delete data if present. relies om is_initialized_!
  if (is_initialized_) {
    delete[] data_;
    is_initialized_ = _FALSE_;
  }
  else {
#ifdef _DEBUG_
    cout << "exiting delete_memory() prematurely." << endl;
#endif
    return 1;
  }

#ifdef _DEBUG_
  cout << "exiting delete_memory()." << endl;
#endif

  return 0;
}


// this function fills data_ with the value given.
int QickArray::fill_data(double value) {
#ifdef _DEBUG_
  cout << "entering fill_data(double value)." << endl;
#endif
  
  if (is_initialized_) {
    for (int i=0; i<no_of_elements_; i++)
      data_[i] = value;
    return 0;
  } else {
    return 1;
  }
#ifdef _DEBUG_
  cout << "exiting fill_data(double value)." << endl;
#endif

}





//just some quick info-functions for the outside world...
int* QickArray::get_dimensions_ptr(void) { return dimensions_; }

int QickArray::get_dimension_info(int *i1, int *i2, int *i3, int *i4, int *i5, int *i6) {
  if (i1) *i1 = dimensions_[0];
  if (i2) *i2 = dimensions_[1];
  if (i3) *i3 = dimensions_[2];
  if (i4) *i4 = dimensions_[3];
  if (i5) *i5 = dimensions_[4];
  if (i6) *i6 = dimensions_[5];
  return n_dims_;
}


int QickArray::is_initialized(void) { return is_initialized_; }

double* QickArray::get_data_ptr(void) { return data_; }


//
// returns logical index from absolute index.
//
int QickArray::logical_index(int idx, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6){

#ifdef _DEBUG_
  cout << "entering logical_index." << endl;
#endif

  int temp = idx;
  int i = 0;

  if (n_dims_ >= 1) {
    i = temp % dimensions_[0];
    if (i1) *i1 = i;

    temp -= i;
    temp /= dimensions_[0];
  }
  if (n_dims_ >= 2) {
    i = temp % dimensions_[1];
    if (i2) *i2 = i;

    temp -= i;
    temp /= dimensions_[1];
  }
  if (n_dims_ >= 3) {
    i = temp % dimensions_[2];
    if (i3) *i3 = i;

    temp -= i;
    temp /= dimensions_[2];
  }
  if (n_dims_ >= 4) {
    i = temp % dimensions_[3];
    if (i4) *i4 = i;

    temp -= i;
    temp /= dimensions_[3];
  }
  if (n_dims_ >= 5) {
    i = temp % dimensions_[4];
    if (i5) *i5 = i;

    temp -= i;
    temp /= dimensions_[4];
  }
  if (n_dims_ >= 6) {
    i = temp % dimensions_[5];
    if (i6) *i6 = i;

    temp -= i;
    temp /= dimensions_[5];
  }
  
  return 0;

}



//get element indexed by i1 ... i6
const double& QickArray::get(int i1, int i2, int i3, int i4, int i5, int i6) const{
  return (*this)(i1,i2,i3,i4,i5,i6);
}

//put value into element indexed by i1 ... i6
int QickArray::put(double value, int i1, int i2, int i3, int i4, int i5, int i6) {
  data_[index(i1,i2,i3,i4,i5,i6)] = value;
  return 0;
}


//
// overloaded assignment operator. clones object.
//
QickArray& QickArray::operator=(const QickArray &clone_me) {
#ifdef _DEBUG_
  cout << "entering operator=(const QickArray&)." << endl;
#endif

  // no need to copy if this == &clone_me.
  if (this == &clone_me) {
#ifdef _DEBUG_
    cout << "exiting operator=(const QickArray&); it was a trivial task." << endl;
#endif
    return *this;

  }

  if ((is_initialized_) && (no_of_elements_ != clone_me.no_of_elements_))
    delete_memory();

  for (int i=0; i<6; i++)
    dimensions_[i] = clone_me.dimensions_[i];
  
  n_dims_ = clone_me.n_dims_;
  no_of_elements_ = clone_me.no_of_elements_;

  step1_ = clone_me.step1_;
  step2_ = clone_me.step2_;
  step3_ = clone_me.step3_;
  step4_ = clone_me.step4_;
  step5_ = clone_me.step5_;

  //copy data in clone_me if clone_me is initialized.
  if (clone_me.is_initialized_) {

    allocate_memory();

    for (int i=0; i<no_of_elements_; i++)
      data_[i] = clone_me.data_[i];
  } 


#ifdef _DEBUG_
  cout << "exiting operator=(const QickArray&)." << endl;
#endif

  return *this; 
}

//
// overloaded +=
//
QickArray& QickArray::operator+=(const QickArray &A) {  //memberwise addition


  int A_no_of_elements = A.no_of_elements_;
  int n = this->no_of_elements_;

  //add as many elements as possible.

  if (A_no_of_elements < n)
    n = A_no_of_elements;

  
  for (int i=0; i<n; i++)
    this->data_[i] += A.data_[i];

  return *this;
}


//
// overloaded -=
//
QickArray& QickArray::operator-=(const QickArray &A) {  //memberwise subtr


  int A_no_of_elements = A.no_of_elements_;
  int n = this->no_of_elements_;

  //add as many elements as possible.

  if (A_no_of_elements < n)
    n = A_no_of_elements;

  
  for (int i=0; i<n; i++)
    this->data_[i] -= A.data_[i];

  return *this;
}

//
// overloaded *=
//
QickArray& QickArray::operator*=(const QickArray &A) {  //memberwise mult


  int A_no_of_elements = A.no_of_elements_;
  int n = this->no_of_elements_;

  //add as many elements as possible.

  if (A_no_of_elements < n)
    n = A_no_of_elements;

  
  for (int i=0; i<n; i++)
    this->data_[i] *= A.data_[i];

  return *this;
}

//
// overloaded /=
//
QickArray& QickArray::operator/=(const QickArray &A) {  //memberwise division


  int A_no_of_elements = A.no_of_elements_;
  int n = this->no_of_elements_;

  //add as many elements as possible.

  if (A_no_of_elements < n)
    n = A_no_of_elements;

  
  for (int i=0; i<n; i++)
    this->data_[i] /= A.data_[i];

  return *this;
}




//
// overloaded += for scalars
//
QickArray& QickArray::operator+=(const double x) {

  int n = this->no_of_elements_;

  for (int i=0; i<n; i++)
    this->data_[i] += x;

  return *this;
}



//
// overloaded -= for scalars
//
QickArray& QickArray::operator-=(const double x) {

  int n = this->no_of_elements_;

  for (int i=0; i<n; i++)
    this->data_[i] -= x;

  return *this;
}



//
// overloaded *= for scalars
//
QickArray& QickArray::operator*=(const double x) {

  int n = this->no_of_elements_;

  for (int i=0; i<n; i++)
    this->data_[i] *= x;

  return *this;
}



//
// overloaded /= for scalars
//
QickArray& QickArray::operator/=(const double x) {

  int n = this->no_of_elements_;

  for (int i=0; i<n; i++)
    this->data_[i] /= x;

  return *this;
}





//
// constructor that copies object.
//
QickArray::QickArray(const QickArray& clone_me) { *this = clone_me; }


//
// function that returns a callable python string.
// if compact==_FALSE_ it is with newlines.
//
string QickArray::python_string(int compact) const {

#ifdef _DEBUG_
  cout << "entering python_string(int=0)." << endl;
#endif

  int indices[] = {0, 0, 0, 0, 0, 0};
  
  string return_me;
  if (is_initialized_)
    return_me = "[" + sub_python_string(n_dims_, indices, compact) + "]";
  else
    return_me = "[ array not initialized ]";

#ifdef _DEBUG_
  cout << "exiting python_string(int=0)." << endl;
#endif


  return return_me;

}


string QickArray::sub_python_string(int dim_no, int indices[], int compact ) const {

#ifdef _DEBUG_
  cout << "entering sub_python_string()." << endl;
#endif

  int i;
  char number[50];

  string newl = ((compact) ? "" : "\n");

  string return_me = "";

  if (dim_no==1) {
    for (i=0; i<dimensions_[0]; i++) {
      sprintf(number, "%10g", get(i, indices[1], indices[2], indices[3], indices[4], indices[5]));
      return_me += number;
      if (i<dimensions_[0]-1)
        return_me += ", ";
    }
  } 
  else {
 
    
    for (i=0; i<dimensions_[dim_no-1]; i++) {

      indices[dim_no-1] = i;
      return_me += newl + "[" + sub_python_string(dim_no-1, indices, compact) + "]" + newl;
      if (i<dimensions_[dim_no-1]-1)
        return_me += ", ";
    }
    
  }

#ifdef _DEBUG_
  cout << "exiting sub_python_string()." << endl;
#endif

  return return_me;
}

//
// skikkelig skitten rapport... litt ubrukelig.
//
string QickArray::summary(void) {

#ifdef _DEBUG_
  cout << "entering summary()"<<endl;
#endif

  string ret;
  char temp[2000];

  if (is_initialized_)
    ret  = "Summary of initialized array:\n";
  else
    ret = "Summary of uninitialized array:\n";

  sprintf(temp, "   min == %g, max == %g, average == %g, \n", min(), max(), average());
  ret += temp;

  sprintf(temp, "   variance == %g, stddev == %g. \n", variance(), stddev());

  ret = ret + temp;

  sprintf(temp, "{ %d, %d, %d, %d, %d, %d }, no_of_elements == %d", dimensions_[0], 
	  dimensions_[1], dimensions_[2],
	  dimensions_[3], dimensions_[4], dimensions_[5], no_of_elements_);
  
  ret += "   dimensions = ";
  ret += temp;
  ret += "\n";
  ret += python_string(_FALSE_);
  ret += "\n";

#ifdef _DEBUG_
  cout << "entering summary()"<<endl;
#endif


  return ret;
  

}

/*

THESE WERE DANGEROUS! MEMORY LEAK! :(

QickArray& QickArray::operator*(const QickArray &B) {

#ifdef _DEBUG_
  cout << "Entering operator*(const QickArray &B)" << endl;
#endif
  QickArray *Result;
  int A_rows = dimensions_[1],   A_cols = dimensions_[0];  //actually "this"...
  int B_rows = B.dimensions_[1], B_cols = B.dimensions_[0]; 

  //both matrices must be 2D...

  if ((n_dims_!=2) || (B.n_dims_!=2)) {
#ifdef _DEBUG_
    cout << "matmult(): matrices are not of correct dimension! returning self.\n";
    cout << "exiting operator*(const QickArray &)" << endl;
#endif
    return (*this);
  }
    
  // no. of rows in B must equal no. of cols in this.

  if (A_cols!=B_rows) {
#ifdef _DEBUG_
    cout << "matmult(): no. of rows in parameter must match no. of cols in self. returning self.\n";
    cout << "exiting operator*(const QickArray &B)" << endl;
#endif

    return (*this);
  }

  //create matrix
  Result = new QickArray(B_cols, A_rows);

  for (int j=0; j<A_rows; j++)
    for (int i=0; i<B_cols; i++) 
      for (int k=0; k<B_rows; k++) {
        (*Result)(i, j) += (*this)(k,j)*B(i,k);
      }

#ifdef _DEBUG_
  cout << "exiting operator*(const QickArray &B)" << endl;
#endif

  return *Result;
}

//
// matrix addition
//
QickArray& QickArray::operator+(const QickArray &B) {
#ifdef _DEBUG_
  cout << "entering operator+(const QickArray &B)" << endl;
#endif

  // A and B must be of same size.

  if (no_of_elements_ != B.no_of_elements_) {
#ifdef _DEBUG_
    cout << "Arrays must match in number og elements during addition. returning self.\n";
    cout << "exiting operator+(const QickArray &B)" << endl;
#endif

    return (*this);

  }

  // add!

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[0], dimensions_[1], dimensions_[2]);


  for (int i=0; i<no_of_elements_; i++)
    Result->data_[i] = data_[i] + B.data_[i];

#ifdef _DEBUG_
    cout << "exiting operator+(const QickArray &B)" << endl;
#endif

  return *Result;
  
}

//
// matrix subtraction
//
QickArray& QickArray::operator-(const QickArray &B) {
#ifdef _DEBUG_
  cout << "entering operator-(const QickArray &B)" << endl;
#endif


  if (no_of_elements_ != B.no_of_elements_) {
#ifdef _DEBUG_
    cout << "Arrays must match in number of elements during subtraction. returning self.\n";
    cout << "exiting operator-(const QickArray &B)" << endl;
#endif

    return (*this);

  }

  // subtract!

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[0], dimensions_[1], dimensions_[2]);


  for (int i=0; i<no_of_elements_; i++)
    Result->data_[i] = data_[i] - B.data_[i];

#ifdef _DEBUG_
    cout << "exiting operator-(const QickArray &B)" << endl;
#endif

  return *Result;
  
}





//
// four simple overloaded operators for manip. of
// array with a double scalar.
//

QickArray& QickArray::operator*(const double x) {
#ifdef _DEBUG_
  cout << "entering operator*(const double)" << endl;
#endif

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[3], dimensions_[4], dimensions_[5]);

  if (is_initialized_)
    for (int i=0; i<no_of_elements_; i++) {
      Result->data_[i] = data_[i] * x;
    }
  
#ifdef _DEBUG_
  cout << "exiting operator*(const double)" << endl;
#endif

  return *Result;
}

QickArray& QickArray::operator/(const double x) {
#ifdef _DEBUG_
  cout << "entering operator/(const double)" << endl;
#endif

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[3], dimensions_[4], dimensions_[5]);

  if (is_initialized_)
    for (int i=0; i<no_of_elements_; i++) {
      Result->data_[i] = data_[i] / x;
    }
  
#ifdef _DEBUG_
  cout << "exiting operator/(const double)" << endl;
#endif

  return *Result;
}

QickArray& QickArray::operator+(const double x) {
#ifdef _DEBUG_
  cout << "entering operator+(const double)" << endl;
#endif

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[3], dimensions_[4], dimensions_[5]);

  if (is_initialized_)
    for (int i=0; i<no_of_elements_; i++) {
      Result->data_[i] = data_[i] + x;
    }
  
#ifdef _DEBUG_
  cout << "exiting operator+(const double)" << endl;
#endif

  return *Result;
}


QickArray& QickArray::operator-(const double x) {
#ifdef _DEBUG_
  cout << "entering operator-(const double)" << endl;
#endif

  QickArray *Result = new QickArray(dimensions_[0], dimensions_[1], dimensions_[2],
				    dimensions_[3], dimensions_[4], dimensions_[5]);

  if (is_initialized_)
    for (int i=0; i<no_of_elements_; i++) {
      Result->data_[i] = data_[i] - x;
    }
  
#ifdef _DEBUG_
  cout << "exiting operator-(const double)" << endl;
#endif

  return *Result;
}


END OF DANGEROUS METHODS!


*/

//
// return minimum of all elements.
//
double QickArray::min(int *i1, int *i2, int *i3, int *i4, int *i5, int *i6) {
  double temp;
  int ind=0;

#ifdef _DEBUG_
  cout << "entering min()"<<endl;
#endif
  
  temp = data_[0];
  for (int i=1; i<no_of_elements_; i++)
    temp = (data_[i]<temp) ? data_[ind=i] : temp ;  //Hey! :P

  logical_index(ind, i1, i2, i3, i4, i5, i6);

  return temp;
}

//
// return (behold!) maximum of all elements.
//
double QickArray::max(int *i1, int *i2, int *i3, int *i4, int *i5, int *i6) {
#ifdef _DEBUG_
  cout << "entering max()"<<endl;
#endif

  double temp;
  int ind=0;
  
  temp = data_[0];
  for (int i=1; i<no_of_elements_; i++)
    temp = (data_[i]>temp) ? data_[ind=i] : temp ;   //Watch this! :P

  logical_index(ind, i1, i2, i3, i4, i5, i6);

  return temp;
}

//
// return average of all elements.
//
double QickArray::average(void) {
#ifdef _DEBUG_
  cout << "entering average()"<<endl;
#endif

  double sum = 0;

  for (int i=0; i<no_of_elements_; i++)
    sum += data_[i];

  return sum/(double)no_of_elements_;

}

//
// return variance.
//
double QickArray::variance(void) {
#ifdef _DEBUG_
  cout << "entering variance()"<<endl;
#endif
  double avg = average();
  double sumsq = 0;

  for (int i=0; i<no_of_elements_; i++)
    sumsq += (data_[i] - avg)*(data_[i] - avg);

  return sumsq/(double)no_of_elements_;

}

//
// return standard deviation.
//
double QickArray::stddev(void) {
#ifdef _DEBUG_
  cout << "entering stddev()"<<endl;
#endif
  return sqrt(variance());
}


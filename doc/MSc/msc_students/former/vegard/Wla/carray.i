%module carray 
%inline %{ 
double *new_carray(int size) { 
       return (double *) malloc(size*sizeof(double)); 
} 
double carray_get(double *a, int i) { 
       return a[i]; 
} 
void carray_set(double *a, int i, double value) { 
     a[i] = value; 
} 
%} 
%rename(delete_darray) free(void *); 
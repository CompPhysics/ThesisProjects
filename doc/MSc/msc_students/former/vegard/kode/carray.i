%module carray 
%inline %{ 
double **new_cmatrix(int row, int col){
       double **matrix = new double *[row];
       for(int j = 0; j < row; j++){
       	       matrix[j] = new double[col];
       }
       return matrix;
}
double cmatrix_get(double **matrix, int i, int j){
       return matrix[i][j];
}
double cmatrix_set(double **matrix, int i, int j, double value){
       matrix[i][j] = value;
}
double *new_carray(int n){//Just for testing...
       double *array = new double[n];
       return array;
}
double carray_get(double *a, int i) { 
       return a[i]; 
} 
void carray_set(double *a, int i, double value) { 
     a[i] = value; 
} 
%} 
//Free arrays???????????????
//Most of these should be void???

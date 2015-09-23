#ifndef InlineJastrow_IS_INCLUDED
#define InlineJastrow_IS_INCLUDED

#ifndef Jastrow_IS_INCLUDED
#error Must be included through Jastrow.h
#endif

/*
void Jastrow::updateProposedMove() {
  _f_matrix = f_matrix-1+current_particle;
  _rij_new_column = rij_new_column-1;
  _f_new_column = f_new_column-1;
  int k=N-1;
  difference = 0;
  for ( int i=0; i < current_particle; i++ ) {
    (*++_f_new_column)=(*f)( (*++_rij_new_column) );
    difference += (*_f_new_column) - (*_f_matrix);
    k--;
    _f_matrix+=k;
  }
  for ( int i=current_particle+1; i<N; i++) {
    (*++_f_new_column)=(*f)( (*++_rij_new_column) );
    difference += (*_f_new_column) - (*++_f_matrix);
  }
}
*/

#endif

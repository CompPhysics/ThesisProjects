      ! dett er en test

      program test
      
      implicit none
      double precision :: c,k,f
      

      write(*,*) 'skrin inn c'
      read(*,*) c
      k = c+273.15
      f = c*1.8+32.
      
      write(*,*) 'k og f'
      write(*,*) k,f
      
      end

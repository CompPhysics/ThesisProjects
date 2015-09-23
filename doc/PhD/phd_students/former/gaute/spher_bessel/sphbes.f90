!     double precision allows for 16 digit accuracy
!     real allows for 8 digits accuracy



      SUBROUTINE sphbes(n,x,sj)
      INTEGER :: n
      DOUBLE PRECISION ::  sj,x      
!     USES bessjy
      DOUBLE PRECISION :: factor,order,rj,RTPIO2    
      PARAMETER (RTPIO2= 1.2533141d0)
!     PARAMETER (RTPIO2= 1.2533141374588012)  
      if(n.lt.0.or.x.le.0.)pause 'bad arguments in sphbes'
      order=n+0.5d0
      CALL bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      return
      END SUBROUTINE sphbes

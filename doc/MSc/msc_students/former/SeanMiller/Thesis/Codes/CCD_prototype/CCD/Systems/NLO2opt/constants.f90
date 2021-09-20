module mconstants
implicit none
integer,parameter :: wp=selected_real_kind(15)
public :: pi,pi2,twopi
public :: hbarc
private

!mathematical constants
real(wp), parameter :: pi    = 3.14159265358979_wp
real(wp), parameter :: pi2   = 9.86960440108936_wp
real(wp), parameter :: twopi = 6.28318530717958_wp

!physics constants
real(wp), parameter :: hbarc = 197.3269718_wp

end module mconstants

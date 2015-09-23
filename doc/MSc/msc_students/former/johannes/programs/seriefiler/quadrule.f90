! documentation is found on 
! http://www.csit.fsu.edu/~burkardt/f_src/quadrule/quadrule.html


subroutine bashforth_set ( norder, xtab, weight )

!*******************************************************************************
!
!! BASHFORTH_SET sets abscissas and weights for Adams-Bashforth quadrature.
!
!  Definition:
!
!    Adams-Bashforth quadrature formulas are normally used in solving
!    ordinary differential equations, and are not really suitable for
!    general quadrature computations.  However, an Adams-Bashforth formula
!    is equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an explicit formula that relies only on known values
!    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
!    have been included here.
!
!    Suppose the unknown function is denoted by Y(X), with derivative
!    F(Y(X)), and that approximate values of the function are known at a
!    series of X values, which we write as X(1), X(2), ..., X(M).  We write
!    the value Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
!             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+1-I)) approximately.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!  Integration interval:
!
!    [ 0, 1 ].
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( 1 - I ),
!
!  Note:
!
!    The Adams-Bashforth formulas require equally spaced data.
!
!    Here is how the formula is applied in the case with non-unit spacing:
!
!      Integral ( A <= X <= A+H ) F(X) dX =
!      H * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( A - (I-1)*H ),
!      approximately.
!
!    The reference lists the second coefficient of the order 8 Adams-Bashforth
!    formula as
!      weight(2) =  -1162169.0D+00 / 120960.0D+00
!    but this should be
!      weight(2) =  -1152169.0D+00 / 120960.0D+00
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Jean Lapidus and John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.  NORDER should be
!    between 1 and 8.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    WEIGHT(1) is the weight at X = 0, WEIGHT(2) the weight at X = -1,
!    and so on.  The weights are rational, and should sum to 1.  Some
!    weights may be negative.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 1 ) then

    weight(1) = 1.0D+00

  else if ( norder == 2 ) then

    d = 2.0D+00

    weight(1) =   3.0D+00 / d
    weight(2) = - 1.0D+00 / d

  else if ( norder == 3 ) then

    d = 12.0D+00

    weight(1) =   23.0D+00 / d
    weight(2) = - 16.0D+00 / d
    weight(3) =    5.0D+00 / d

  else if ( norder == 4 ) then

    d = 24.0D+00

    weight(1) =   55.0D+00 / d
    weight(2) = - 59.0D+00 / d
    weight(3) =   37.0D+00 / d
    weight(4) =  - 9.0D+00 / d

  else if ( norder == 5 ) then

    d = 720.0D+00

    weight(1) =   1901.0D+00 / d
    weight(2) = - 2774.0D+00 / d
    weight(3) =   2616.0D+00 / d
    weight(4) = - 1274.0D+00 / d
    weight(5) =    251.0D+00 / d

  else if ( norder == 6 ) then

    d = 1440.0D+00

    weight(1) =   4277.0D+00 / d
    weight(2) = - 7923.0D+00 / d
    weight(3) =   9982.0D+00 / d
    weight(4) = - 7298.0D+00 / d
    weight(5) =   2877.0D+00 / d
    weight(6) =  - 475.0D+00 / d

  else if ( norder == 7 ) then

    d = 60480.0D+00

    weight(1) =    198721.0D+00 / d
    weight(2) =  - 447288.0D+00 / d
    weight(3) =    705549.0D+00 / d
    weight(4) =  - 688256.0D+00 / d
    weight(5) =    407139.0D+00 / d
    weight(6) =  - 134472.0D+00 / d
    weight(7) =     19087.0D+00 / d

  else if ( norder == 8 ) then

    d = 120960.0D+00

    weight(1) =     434241.0D+00 / d
    weight(2) =  - 1152169.0D+00 / d
    weight(3) =    2183877.0D+00 / d
    weight(4) =  - 2664477.0D+00 / d
    weight(5) =    2102243.0D+00 / d
    weight(6) =  - 1041723.0D+00 / d
    weight(7) =     295767.0D+00 / d
    weight(8) =    - 36799.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASHFORTH_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 through 8.'
    stop

  end if

  do i = 1, norder
    xtab(i) = real ( 1 - i, kind = 8 )
  end do

  return
end
subroutine bdf_set ( norder, alpha, beta, gamma )

!*******************************************************************************
!
!! BDF_SET sets weights for backward differentiation ODE weights.
!
!  Discussion:
!
!    GAMMA * Y(N+1) = Sum ( 1 <= I <= NORDER ) ALPHA(I) * Y(N+1-I)
!                     + dX * BETA * Y'(X(N+1),Y(N+1))
!
!    This is equivalent to the backward differentiation corrector formulas.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, between 1 and 6.
!
!    Output, real ( kind = 8 ) ALPHA(NORDER), BETA, GAMMA, the weights.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha(norder)
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
!
  if ( norder == 1 ) then
    beta =     1.0D+00
    gamma =    1.0D+00
    alpha(1) = 1.0D+00
  else if ( norder == 2 ) then
    beta =       2.0D+00
    gamma =      3.0D+00
    alpha(1) =   4.0D+00
    alpha(2) = - 1.0D+00
  else if ( norder == 3 ) then
    beta =       6.0D+00
    gamma =     11.0D+00
    alpha(1) =  18.0D+00
    alpha(2) = - 9.0D+00
    alpha(3) =   2.0D+00
  else if ( norder == 4 ) then
    beta =       12.0D+00
    gamma =      25.0D+00
    alpha(1) =   48.0D+00
    alpha(2) = - 36.0D+00
    alpha(3) =   16.0D+00
    alpha(4) =  - 3.0D+00
  else if ( norder == 5 ) then
    beta =        60.0D+00
    gamma =      137.0D+00
    alpha(1) =   300.0D+00
    alpha(2) = - 300.0D+00
    alpha(3) =   200.0D+00
    alpha(4) =  - 75.0D+00
    alpha(5) =    12.0D+00
  else if ( norder == 6 ) then
    beta =        60.0D+00
    gamma =      147.0D+00
    alpha(1) =   360.0D+00
    alpha(2) = - 450.0D+00
    alpha(3) =   400.0D+00
    alpha(4) = - 225.0D+00
    alpha(5) =    72.0D+00
    alpha(6) =  - 10.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BDF_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal order requested = ', norder
    stop
  end if

  return
end
subroutine bdfc_set ( norder, weight, xtab )
!
!*******************************************************************************
!
!! BDFC_SET sets weights for backward differentiation corrector quadrature.
!
!
!  Definition:
!
!    A backward differentiation corrector formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(1), and the backward differences at X(1) that approximate the
!    derivatives there.
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BD**(I-1) F ( 1 ).
!
!    Here, "BD**(I-1) F ( 1 )" denotes the (I-1)st backward difference
!    of F at X = 1, using a spacing of 1.  In particular,
!
!    BD**0 F(1) = F(1)
!    BD**1 F(1) = F(1) - F(0)
!    BD**2 F(1) = F(1) - 2 * F(0) + F(-1 )
!
!  Note:
!
!    The relationship between a backward difference corrector and the
!    corresponding Adams-Moulton formula may be illustrated for the
!    BDF corrector of order 4:
!
!      BD**0 F(1) - 1/2 * BD**1 F(1) - 1/12 * BD**2 F(1) - 1/24 * BDF**3 F(1)
!      =            F(1)
!        -  1/2 * ( F(1) -         F(0) )
!        - 1/12 * ( F(1) - 2     * F(0) +        F(-1) )
!        - 1/24 * ( F(1) - 3     * F(0) + 3    * F(-1) -        F(-2) )
!      =   9/24 *   F(1) + 19/24 * F(0) - 5/24 * F(-1) + 1/24 * F(-2)
!
!    which is the Adams-Moulton formula of order 4.
! 
!  Modified:
!
!    28 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary Differential
!      Equations,
!    Academic Press, 1988.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which can be
!    any value from 1 to 19.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
  implicit none
!
  integer, parameter :: maxord = 19
!
  integer norder
!
  integer i
  real ( kind = 8 ) w(maxord)
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  w(1) =                 1.0D+00
  w(2) =               - 1.0D+00 /               2.0D+00
  w(3) =               - 1.0D+00 /              12.0D+00
  w(4) =               - 1.0D+00 /              24.0D+00
  w(5) =              - 19.0D+00 /             720.0D+00
  w(6) =               - 3.0D+00 /             160.0D+00
  w(7) =             - 863.0D+00 /           60480.0D+00
  w(8) =             - 275.0D+00 /           24792.0D+00
  w(9) =           - 33953.0D+00 /         3628800.0D+00
  w(10) =           - 8183.0D+00 /         1036800.0D+00
  w(11) =        - 3250433.0D+00 /       479001600.0D+00
  w(12) =           - 4671.0D+00 /          788480.0D+00
  w(13) =    - 13695779093.0D+00 /   2615348736000.0D+00
  w(14) =     - 2224234463.0D+00 /    475517952000.0D+00
  w(15) =   - 132282840127.0D+00 /  31384184832000.0D+00
  w(16) =     - 2639651053.0D+00 /    689762304000.0D+00
  w(17) =  111956703448001.0D+00 /   3201186852864.0D+00
  w(18) =         50188465.0D+00 /     15613165568.0D+00
  w(19) = 2334028946344463.0D+00 / 786014494949376.0D+00

  do i = 1, min ( norder, maxord )
    weight(i) = w(i)
  end do

  do i = 1, norder
    xtab(i) = real ( 2 - i, kind = 8 )
  end do

  return
end
subroutine bdfp_set ( norder, weight, xtab )
!
!*******************************************************************************
!
!! BDFP_SET sets weights for backward differentiation predictor quadrature.
!
!
!  Definition:
!
!    A backward differentiation predictor formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(2), and the backward differences at X(2) that approximate the
!    derivatives there.  A backward differentiation predictor formula
!    is equivalent to an Adams-Bashforth formula of the same order.
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BD**(I-1) F ( 0 ),
!
!    Here, "BD**(I-1) F ( 0 )" denotes the (I-1)st backward difference
!    of F at X = 0, using a spacing of 1.  In particular,
!
!    BD**0 F(0) = F(0)
!    BD**1 F(0) = F(0) - F(-1)
!    BD**2 F(0) = F(0) - 2 * F(-1) + F(-2 )
!
!  Note:
!
!    The relationship between a backward difference predictor and the
!    corresponding Adams-Bashforth formula may be illustrated for the
!    BDF predictor of order 3:
!
!      BD**0 F(0) + 0.5 * BD**1 F(0) + 5/12 * BD**2 F(0)
!      =            F(0)
!        + 1/2  * ( F(0) -         F(1) )
!        + 5/12 * ( F(0) - 2     * F(-1) +      F(-2) )
!      =  23/12 *   F(0) - 16/12 * F(-1) + 5/12 F(-2)
!
!    which is the Adams-Bashforth formula of order 3.
! 
!  Modified:
!
!    29 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary Differential
!      Equations,
!    Academic Press, 1988.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which can be
!    any value from 1 to 19.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weight of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
  implicit none
!
  integer, parameter :: maxord = 19
!
  integer norder
!
  integer i
  real ( kind = 8 ) w(maxord)
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  w(1) =                       1.0D+00
  w(2) =                       1.0D+00 /                2.0D+00
  w(3) =                       5.0D+00 /               12.0D+00
  w(4) =                       3.0D+00 /                8.0D+00
  w(5) =                     251.0D+00 /              720.0D+00
  w(6) =                      95.0D+00 /              288.0D+00
  w(7) =                   19087.0D+00 /            60480.0D+00
  w(8) =                    5257.0D+00 /            17280.0D+00
  w(9) =                 1070017.0D+00 /          3628800.0D+00
  w(10) =                  25713.0D+00 /            89600.0D+00
  w(11) =               26842253.0D+00 /         95800320.0D+00
  w(12) =                4777223.0D+00 /         17418240.0D+00
  w(13) =           703604254357.0D+00 /    2615348736000.0D+00
  w(14) =           106364763817.0D+00 /     402361344000.0D+00
  w(15) =          1166309819657.0D+00 /    4483454976000.0D+00
  w(16) =               25221445.0D+00 /         98402304.0D+00
  w(17) =       8092989203533249.0D+00 /    3201186852864.0D+00
  w(18) =         85455477715379.0D+00 /      34237292544.0D+00
  w(19) =   12600467236042756559.0D+00 / 5109094217170944.0D+00

  do i = 1, min ( norder, maxord )
    weight(i) = w(i)
  end do

  do i = 1, norder
    xtab(i) = real ( 1 - i, kind = 8 )
  end do

  return
end
subroutine bdf_sum ( func, norder, weight, xtab, diftab, result )
!
!*******************************************************************************
!
!! BDF_SUM carries out an explicit backward difference quadrature rule for [0,1].
!
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX
!
!  Formula:
!
!    RESULT = Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BDF**(I-1) FUNC ( 0 )
!
!  Note:
!
!    The integral from 0 to 1 is approximated using data at X = 0,
!    -1, -2, ..., -NORDER+1.  This is a form of extrapolation, and
!    the approximation can become poor as NORDER increases.
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which evaluates
!    the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
!    Input, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Workspace, real ( kind = 8 ) DIFTAB(NORDER).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) diftab(norder)
  real ( kind = 8 ), external :: func
  integer i
  integer j
  real ( kind = 8 ) result
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  do i = 1, norder
    diftab(i) = func ( xtab(i) )
  end do

  do i = 2, norder
    do j = i, norder
      diftab(norder+i-j) = ( diftab(norder+i-j-1) - diftab(norder+i-j) )
    end do
  end do

  result = dot_product ( weight(1:norder), diftab(1:norder) )

  return
end
subroutine cheb_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! CHEB_SET sets abscissas and weights for Chebyshev quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Note:
!
!    The Chebyshev rule is distinguished by using equal weights.
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    H Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
!    There are NO other Chebyshev rules with real abscissas.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule,
!    which are symmetric in [-1,1].
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule,
!    which should each equal 2 / NORDER.
!
  implicit none
!
  integer norder
!
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 1 ) then

    xtab(1) = 0.0D+00

  else if ( norder == 2 ) then

    xtab(1) = - 1.0D+00 / sqrt ( 3.0D+00 )
    xtab(2) =   1.0D+00 / sqrt ( 3.0D+00 )

  else if ( norder == 3 ) then

    xtab(1) = - 1.0D+00 / sqrt ( 2.0D+00 )
    xtab(2) =   0.0D+00
    xtab(3) =   1.0D+00 / sqrt ( 2.0D+00 )

  else if ( norder == 4 ) then

    xtab(1) =   - sqrt ( ( 1.0D+00 + 2.0D+00/ sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    xtab(2) =   - sqrt ( ( 1.0D+00 - 2.0D+00/ sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    xtab(3) =     sqrt ( ( 1.0D+00 - 2.0D+00/ sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    xtab(4) =     sqrt ( ( 1.0D+00 + 2.0D+00/ sqrt ( 5.0D+00 ) ) / 3.0D+00 )

  else if ( norder == 5 ) then

    xtab(1) = - sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00) ) / 12.0D+00 )
    xtab(2) = - sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00) ) / 12.0D+00 )
    xtab(3) =   0.0D+00
    xtab(4) =   sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00) ) / 12.0D+00 )
    xtab(5) =   sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00) ) / 12.0D+00 )

  else if ( norder == 6 ) then

    xtab(1) = - 0.866246818107820591383598D+00
    xtab(2) = - 0.422518653761111529118546D+00
    xtab(3) = - 0.266635401516704720331534D+00
    xtab(4) =   0.266635401516704720331534D+00
    xtab(5) =   0.422518653761111529118546D+00
    xtab(6) =   0.866246818107820591383598D+00

  else if ( norder == 7 ) then

    xtab(1) = - 0.883861700758049035704224D+00
    xtab(2) = - 0.529656775285156811385048D+00
    xtab(3) = - 0.323911810519907637519673D+00
    xtab(4) =   0.0D+00
    xtab(5) =   0.323911810519907637519673D+00
    xtab(6) =   0.529656775285156811385048D+00
    xtab(7) =   0.883861700758049035704224D+00

  else if ( norder == 9 ) then

    xtab(1) = - 0.911589307728434473664949D+00
    xtab(2) = - 0.601018655380238071428128D+00
    xtab(3) = - 0.528761783057879993260181D+00
    xtab(4) = - 0.167906184214803943068031D+00
    xtab(5) =   0.0D+00
    xtab(6) =   0.167906184214803943068031D+00
    xtab(7) =   0.528761783057879993260181D+00
    xtab(8) =   0.601018655380238071428128D+00
    xtab(9) =   0.911589307728434473664949D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEB_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 through 7, and 9.'
    stop

  end if

  weight(1:norder) = 2.0D+00 / real ( norder, kind = 8 )

  return
end
subroutine cheb_to_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! CHEB_TO_SET sets up open Gauss-Chebyshev (first kind) quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 / SQRT ( 1 - X**2 )
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) / SQRT ( 1 - X**2 ) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-1 or less.
!
!  Note:
!
!    The abscissas of the rule are zeroes of the Chebyshev polynomials
!    of the first kind, T(NORDER)(X).
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule,
!    which are all equal to PI / NORDER.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) angle
  real ( kind = 8 ) d_pi
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  do i = 1, norder
    angle = real ( 2 * i - 1, kind = 8 ) * d_pi ( ) &
          / real ( 2 * norder, kind = 8 )
    xtab(i) = cos ( angle )
  end do

  weight(1:norder) = d_pi ( ) / real ( norder, kind = 8 )

  return
end
subroutine cheb_tc_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! CHEB_TC_SET sets up closed Gauss-Chebyshev (first kind) quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 / SQRT ( 1 - X**2 )
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) / SQRT ( 1 - X**2 ) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-3 or less.
!
!  Note:
!
!    The abscissas include -1 and 1.
!
!    If the order is doubled, the abscissas of the new rule include
!    all the points of the old rule.  This fact can be used to
!    efficiently implement error estimation.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which must be at least 2.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The first and last weights are 0.5 * PI / ( NORDER - 1),
!    and all other weights are PI / ( NORDER - 1 ).
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) angle
  real ( kind = 8 ) d_pi
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEB_TC_SET - Fatal error!'
    write ( *, '(a)' ) '  NORDER must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value was NORDER = ', norder
    stop
  end if

  do i = 1, norder

    angle = real (      i - 1, kind = 8 ) * d_pi ( ) &
          / real ( norder - 1, kind = 8 )
    xtab(i) = cos ( angle )

  end do

  weight(1) = d_pi ( ) / real ( 2 * ( norder - 1 ), kind = 8 )
  weight(2:norder-1) = d_pi ( ) / real ( norder - 1, kind = 8 )
  weight(norder) = d_pi ( ) / real ( 2 * ( norder - 1 ), kind = 8 )

  return
end
subroutine cheb_u_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! CHEB_U_SET sets abscissas and weights for Gauss-Chebyshev quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    SQRT ( 1 - X**2 )
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) SQRT ( 1 - X**2 ) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-1 or less.
!
!  Note:
!
!    The abscissas are zeroes of the Chebyshev polynomials
!    of the second kind, U(NORDER)(X).
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule,
!    which are all equal to PI / NORDER.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) angle
  real ( kind = 8 ) d_pi
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  do i = 1, norder
    angle = real ( i, kind = 8 ) * d_pi ( ) / real ( norder + 1, kind = 8 )
    xtab(i) = cos ( angle )
    weight(i) = d_pi ( ) * ( sin ( angle ) )**2 / real ( norder + 1, kind = 8 )
  end do

  return
end
subroutine d_swap ( x, y )

!*******************************************************************************
!
!! D_SWAP switches two double precision values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
  z = x
  x = y
  y = z

  return
end
function d_pi ( )

!*******************************************************************************
!
!! D_PI returns the value of pi as a double precision quantity.
!
!  Modified:
!
!    28 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) D_PI, the value of pi.
!
  implicit none
!
  real ( kind = 8 ) d_pi
!
  d_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine dvec_reverse ( n, a )

!*******************************************************************************
!
!! DVEC_REVERSE reverses the elements of a double precision vector.
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0D+00 )
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0D+00 )
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none
!
  integer n
!
  real ( kind = 8 ) a(n)
  integer i
!
  do i = 1, n/2
    call d_swap ( a(i), a(n+1-i) )
  end do

  return
end
function gamma ( x )
!
!*******************************************************************************
!
!! GAMMA computes the gamma function using Hastings's approximation.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument at which the gamma function
!    is to be evaluated.  X must be greater than 0, and less than 70.
!
!    Output, real ( kind = 8 ) GAMMA, the gamma function at X.
!
  implicit none
!
  real ( kind = 8 ) gam
  real ( kind = 8 ) gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) za
!
  gam ( y ) = ((((((( &
          0.035868343D+00   * y &
        - 0.193527818D+00 ) * y &
        + 0.482199394D+00 ) * y &
        - 0.756704078D+00 ) * y &
        + 0.918206857D+00 ) * y &
        - 0.897056937D+00 ) * y &
        + 0.988205891D+00 ) * y &
        - 0.577191652D+00 ) * y + 1.0D+00

  if ( x <= 0.0D+00 ) then
    gamma = 0.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA - Fatal error!'
    write ( *, '(a)' ) '  Input argument X <= 0.'
    stop
  end if

  if ( 70.0D+00 <= x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA - Fatal error!'
    write ( *, '(a)' ) '  Input argument 70 <= X.'
    stop
  end if

  if ( x == 1.0D+00 ) then
    gamma = 1.0D+00
    return
  end if

  if ( x <= 1.0D+00 ) then
    gamma = gam ( x ) / x
    return
  end if

  z = x

  za = 1.0D+00

  do

    z = z - 1.0D+00

    if ( z < 1.0D+00 ) then
      gamma = za * gam ( z )
      exit
    else if ( z == 1.0D+00 ) then
      gamma = za
      exit
    end if

    za = za * z

  end do

  return
end
subroutine hermite_com ( norder, xtab, weight )
!
!*******************************************************************************
!
!! HERMITE_COM computes the abscissa and weights for Gauss-Hermite quadrature.
!
!
!  Discussion:
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!  Integration interval:
!
!    ( -Infinity, +Infinity )
!
!  Weight function:
!
!    exp ( - X**2 )
!
!  Integral to approximate:
!
!    Integral ( -INFINITY < X < +INFINITY ) exp ( - X**2 ) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the formula to be computed.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the Gauss-Hermite abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the Gauss-Hermite weights.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  real ( kind = 8 ) gamma
  integer i
  real ( kind = 8 ) p1
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(norder)
!
  cc = 1.7724538509D+00 * gamma ( real ( norder, kind = 8 ) ) &
    / ( 2.0D+00**( norder - 1 ) )

  s = ( 2.0D+00 * real ( norder, kind = 8 ) + 1.0D+00 )**( 1.0D+00 / 6.0D+00 )

  do i = 1, ( norder + 1 ) / 2

    if ( i == 1 ) then

      x = s**3 - 1.85575D+00 / s

    else if ( i == 2 ) then

      x = x - 1.14D+00 * ( ( real ( norder, kind = 8 ) )**0.426D+00 ) / x

    else if ( i == 3 ) then

      x = 1.86D+00 * x - 0.86D+00 * xtab(1)

    else if ( i == 4 ) then

      x = 1.91D+00 * x - 0.91D+00 * xtab(2)

    else

      x = 2.0D+00 * x - xtab(i-2)

    end if

    call hermite_root ( x, norder, dp2, p1 )

    xtab(i) = x
    weight(i) = ( cc / dp2 ) / p1

    xtab(norder-i+1) = - x
    weight(norder-i+1) = weight(i)

  end do
!
!  Reverse the order of the XTAB values.
!
  call dvec_reverse ( norder, xtab )

  return
end
subroutine hermite_recur ( p2, dp2, p1, x, norder )
!
!*******************************************************************************
!
!! HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of H(NORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of H'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
  implicit none
!
  integer i
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer norder
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, norder

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2  = x * p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * p0
    dp2 = x * dp1 + p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * dp0

  end do

  return
end
subroutine hermite_root ( x, norder, dp2, p1 )
!
!*******************************************************************************
!
!! HERMITE_ROOT improves an approximate root of a Hermite polynomial.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the Hermite polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of H'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(NORDER-1)(X).
!
  implicit none
!
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ), parameter :: eps = 1.0D-12
  integer i
  integer, parameter :: maxstep = 10
  integer norder
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  do i = 1, maxstep

    call hermite_recur ( p2, dp2, p1, x, norder )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine hermite_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! HERMITE_SET sets abscissas and weights for Hermite quadrature.
!
!
!  Integration interval:
!
!    ( -Infinity, +Infinity )
!
!  Weight function:
!
!    exp ( - X**2 )
!
!  Integral to approximate:
!
!    Integral ( -INFINITY < X < +INFINITY ) exp ( - X**2 ) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20, or one of the values
!    30, 32, 40, 50, 60 or 64.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule,
!    which are symmetrically placed around 0.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive and symmetric, and should sum
!    to SQRT(PI).
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d_pi
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) = 0.0D+00

    weight(1) = sqrt ( d_pi ( ) )

  else if ( norder == 2 )then

    xtab(1) = - 0.707106781186547524400844362105D+00
    xtab(2) =   0.707106781186547524400844362105D+00

    weight(1) = 0.886226925452758013649083741671D+00
    weight(2) = 0.886226925452758013649083741671D+00

  else if ( norder == 3 ) then

    xtab(1) = - 0.122474487139158904909864203735D+01
    xtab(2) =   0.0D+00
    xtab(3) =   0.122474487139158904909864203735D+01

    weight(1) = 0.295408975150919337883027913890D+00
    weight(2) = 0.118163590060367735153211165556D+01
    weight(3) = 0.295408975150919337883027913890D+00

  else if ( norder == 4 ) then

    xtab(1) = - 0.165068012388578455588334111112D+01
    xtab(2) = - 0.524647623275290317884060253835D+00
    xtab(3) =   0.524647623275290317884060253835D+00
    xtab(4) =   0.165068012388578455588334111112D+01

    weight(1) = 0.813128354472451771430345571899D-01
    weight(2) = 0.804914090005512836506049184481D+00
    weight(3) = 0.804914090005512836506049184481D+00
    weight(4) = 0.813128354472451771430345571899D-01

  else if ( norder == 5 ) then

    xtab(1) = - 0.202018287045608563292872408814D+01
    xtab(2) = - 0.958572464613818507112770593893D+00
    xtab(3) =   0.0D+00
    xtab(4) =   0.958572464613818507112770593893D+00
    xtab(5) =   0.202018287045608563292872408814D+01

    weight(1) = 0.199532420590459132077434585942D-01
    weight(2) = 0.393619323152241159828495620852D+00
    weight(3) = 0.945308720482941881225689324449D+00
    weight(4) = 0.393619323152241159828495620852D+00
    weight(5) = 0.199532420590459132077434585942D-01

  else if ( norder == 6 ) then

    xtab(1) = - 0.235060497367449222283392198706D+01
    xtab(2) = - 0.133584907401369694971489528297D+01
    xtab(3) = - 0.436077411927616508679215948251D+00
    xtab(4) =   0.436077411927616508679215948251D+00
    xtab(5) =   0.133584907401369694971489528297D+01
    xtab(6) =   0.235060497367449222283392198706D+01

    weight(1) = 0.453000990550884564085747256463D-02
    weight(2) = 0.157067320322856643916311563508D+00
    weight(3) = 0.724629595224392524091914705598D+00
    weight(4) = 0.724629595224392524091914705598D+00
    weight(5) = 0.157067320322856643916311563508D+00
    weight(6) = 0.453000990550884564085747256463D-02

  else if ( norder == 7 ) then

    xtab(1) = - 0.265196135683523349244708200652D+01
    xtab(2) = - 0.167355162876747144503180139830D+01
    xtab(3) = - 0.816287882858964663038710959027D+00
    xtab(4) =   0.0D+00
    xtab(5) =   0.816287882858964663038710959027D+00
    xtab(6) =   0.167355162876747144503180139830D+01
    xtab(7) =   0.265196135683523349244708200652D+01

    weight(1) = 0.971781245099519154149424255939D-03
    weight(2) = 0.545155828191270305921785688417D-01
    weight(3) = 0.425607252610127800520317466666D+00
    weight(4) = 0.810264617556807326764876563813D+00
    weight(5) = 0.425607252610127800520317466666D+00
    weight(6) = 0.545155828191270305921785688417D-01
    weight(7) = 0.971781245099519154149424255939D-03

  else if ( norder == 8 ) then

    xtab(1) = - 0.293063742025724401922350270524D+01
    xtab(2) = - 0.198165675669584292585463063977D+01
    xtab(3) = - 0.115719371244678019472076577906D+01
    xtab(4) = - 0.381186990207322116854718885584D+00
    xtab(5) =   0.381186990207322116854718885584D+00
    xtab(6) =   0.115719371244678019472076577906D+01
    xtab(7) =   0.198165675669584292585463063977D+01
    xtab(8) =   0.293063742025724401922350270524D+01

    weight(1) = 0.199604072211367619206090452544D-03
    weight(2) = 0.170779830074134754562030564364D-01
    weight(3) = 0.207802325814891879543258620286D+00
    weight(4) = 0.661147012558241291030415974496D+00
    weight(5) = 0.661147012558241291030415974496D+00
    weight(6) = 0.207802325814891879543258620286D+00
    weight(7) = 0.170779830074134754562030564364D-01
    weight(8) = 0.199604072211367619206090452544D-03

  else if ( norder == 9 ) then

    xtab(1) = - 0.319099320178152760723004779538D+01
    xtab(2) = - 0.226658058453184311180209693284D+01
    xtab(3) = - 0.146855328921666793166701573925D+01
    xtab(4) = - 0.723551018752837573322639864579D+00
    xtab(5) =   0.0D+00
    xtab(6) =   0.723551018752837573322639864579D+00
    xtab(7) =   0.146855328921666793166701573925D+01
    xtab(8) =   0.226658058453184311180209693284D+01
    xtab(9) =   0.319099320178152760723004779538D+01

    weight(1) = 0.396069772632643819045862946425D-04
    weight(2) = 0.494362427553694721722456597763D-02
    weight(3) = 0.884745273943765732879751147476D-01
    weight(4) = 0.432651559002555750199812112956D+00
    weight(5) = 0.720235215606050957124334723389D+00
    weight(6) = 0.432651559002555750199812112956D+00
    weight(7) = 0.884745273943765732879751147476D-01
    weight(8) = 0.494362427553694721722456597763D-02
    weight(9) = 0.396069772632643819045862946425D-04

  else if ( norder == 10 ) then

    xtab(1) =  - 0.343615911883773760332672549432D+01
    xtab(2) =  - 0.253273167423278979640896079775D+01
    xtab(3) =  - 0.175668364929988177345140122011D+01
    xtab(4) =  - 0.103661082978951365417749191676D+01
    xtab(5) =  - 0.342901327223704608789165025557D+00
    xtab(6) =    0.342901327223704608789165025557D+00
    xtab(7) =    0.103661082978951365417749191676D+01
    xtab(8) =    0.175668364929988177345140122011D+01
    xtab(9) =    0.253273167423278979640896079775D+01
    xtab(10) =   0.343615911883773760332672549432D+01

    weight(1) =  0.764043285523262062915936785960D-05
    weight(2) =  0.134364574678123269220156558585D-02
    weight(3) =  0.338743944554810631361647312776D-01
    weight(4) =  0.240138611082314686416523295006D+00
    weight(5) =  0.610862633735325798783564990433D+00
    weight(6) =  0.610862633735325798783564990433D+00
    weight(7) =  0.240138611082314686416523295006D+00
    weight(8) =  0.338743944554810631361647312776D-01
    weight(9) =  0.134364574678123269220156558585D-02
    weight(10) = 0.764043285523262062915936785960D-05

  else if ( norder == 11 ) then

    xtab(1) =  - 0.366847084655958251845837146485D+01
    xtab(2) =  - 0.278329009978165177083671870152D+01
    xtab(3) =  - 0.202594801582575533516591283121D+01
    xtab(4) =  - 0.132655708449493285594973473558D+01
    xtab(5) =  - 0.656809566882099765024611575383D+00
    xtab(6) =    0.0D+00
    xtab(7) =    0.656809566882099765024611575383D+00
    xtab(8) =    0.132655708449493285594973473558D+01
    xtab(9) =    0.202594801582575533516591283121D+01
    xtab(10) =   0.278329009978165177083671870152D+01
    xtab(11) =   0.366847084655958251845837146485D+01

    weight(1) =  0.143956039371425822033088366032D-05
    weight(2) =  0.346819466323345510643413772940D-03
    weight(3) =  0.119113954449115324503874202916D-01
    weight(4) =  0.117227875167708503381788649308D+00
    weight(5) =  0.429359752356125028446073598601D+00
    weight(6) =  0.654759286914591779203940657627D+00
    weight(7) =  0.429359752356125028446073598601D+00
    weight(8) =  0.117227875167708503381788649308D+00
    weight(9) =  0.119113954449115324503874202916D-01
    weight(10) = 0.346819466323345510643413772940D-03
    weight(11) = 0.143956039371425822033088366032D-05

  else if ( norder == 12 ) then

    xtab(1) =  - 0.388972489786978191927164274724D+01
    xtab(2) =  - 0.302063702512088977171067937518D+01
    xtab(3) =  - 0.227950708050105990018772856942D+01
    xtab(4) =  - 0.159768263515260479670966277090D+01
    xtab(5) =  - 0.947788391240163743704578131060D+00
    xtab(6) =  - 0.314240376254359111276611634095D+00
    xtab(7) =    0.314240376254359111276611634095D+00
    xtab(8) =    0.947788391240163743704578131060D+00
    xtab(9) =    0.159768263515260479670966277090D+01
    xtab(10) =   0.227950708050105990018772856942D+01
    xtab(11) =   0.302063702512088977171067937518D+01
    xtab(12) =   0.388972489786978191927164274724D+01

    weight(1) =  0.265855168435630160602311400877D-06
    weight(2) =  0.857368704358785865456906323153D-04
    weight(3) =  0.390539058462906185999438432620D-02
    weight(4) =  0.516079856158839299918734423606D-01
    weight(5) =  0.260492310264161129233396139765D+00
    weight(6) =  0.570135236262479578347113482275D+00
    weight(7) =  0.570135236262479578347113482275D+00
    weight(8) =  0.260492310264161129233396139765D+00
    weight(9) =  0.516079856158839299918734423606D-01
    weight(10) = 0.390539058462906185999438432620D-02
    weight(11) = 0.857368704358785865456906323153D-04
    weight(12) = 0.265855168435630160602311400877D-06

  else if ( norder == 13 ) then

    xtab(1) =  - 0.410133759617863964117891508007D+01
    xtab(2) =  - 0.324660897837240998812205115236D+01
    xtab(3) =  - 0.251973568567823788343040913628D+01
    xtab(4) =  - 0.185310765160151214200350644316D+01
    xtab(5) =  - 0.122005503659074842622205526637D+01
    xtab(6) =  - 0.605763879171060113080537108602D+00
    xtab(7) =    0.0D+00
    xtab(8) =    0.605763879171060113080537108602D+00
    xtab(9) =    0.122005503659074842622205526637D+01
    xtab(10) =   0.185310765160151214200350644316D+01
    xtab(11) =   0.251973568567823788343040913628D+01
    xtab(12) =   0.324660897837240998812205115236D+01
    xtab(13) =   0.410133759617863964117891508007D+01

    weight(1) =  0.482573185007313108834997332342D-07
    weight(2) =  0.204303604027070731248669432937D-04
    weight(3) =  0.120745999271938594730924899224D-02
    weight(4) =  0.208627752961699392166033805050D-01
    weight(5) =  0.140323320687023437762792268873D+00
    weight(6) =  0.421616296898543221746893558568D+00
    weight(7) =  0.604393187921161642342099068579D+00
    weight(8) =  0.421616296898543221746893558568D+00
    weight(9) =  0.140323320687023437762792268873D+00
    weight(10) = 0.208627752961699392166033805050D-01
    weight(11) = 0.120745999271938594730924899224D-02
    weight(12) = 0.204303604027070731248669432937D-04
    weight(13) = 0.482573185007313108834997332342D-07

  else if ( norder == 14 ) then

    xtab(1) =  - 0.430444857047363181262129810037D+01
    xtab(2) =  - 0.346265693360227055020891736115D+01
    xtab(3) =  - 0.274847072498540256862499852415D+01
    xtab(4) =  - 0.209518325850771681573497272630D+01
    xtab(5) =  - 0.147668273114114087058350654421D+01
    xtab(6) =  - 0.878713787329399416114679311861D+00
    xtab(7) =  - 0.291745510672562078446113075799D+00
    xtab(8) =    0.291745510672562078446113075799D+00
    xtab(9) =    0.878713787329399416114679311861D+00
    xtab(10) =   0.147668273114114087058350654421D+01
    xtab(11) =   0.209518325850771681573497272630D+01
    xtab(12) =   0.274847072498540256862499852415D+01
    xtab(13) =   0.346265693360227055020891736115D+01
    xtab(14) =   0.430444857047363181262129810037D+01

    weight(1) =  0.862859116812515794532041783429D-08
    weight(2) =  0.471648435501891674887688950105D-05
    weight(3) =  0.355092613551923610483661076691D-03
    weight(4) =  0.785005472645794431048644334608D-02
    weight(5) =  0.685055342234652055387163312367D-01
    weight(6) =  0.273105609064246603352569187026D+00
    weight(7) =  0.536405909712090149794921296776D+00
    weight(8) =  0.536405909712090149794921296776D+00
    weight(9) =  0.273105609064246603352569187026D+00
    weight(10) = 0.685055342234652055387163312367D-01
    weight(11) = 0.785005472645794431048644334608D-02
    weight(12) = 0.355092613551923610483661076691D-03
    weight(13) = 0.471648435501891674887688950105D-05
    weight(14) = 0.862859116812515794532041783429D-08

  else if ( norder == 15 ) then

    xtab(1) =  - 0.449999070730939155366438053053D+01
    xtab(2) =  - 0.366995037340445253472922383312D+01
    xtab(3) =  - 0.296716692790560324848896036355D+01
    xtab(4) =  - 0.232573248617385774545404479449D+01
    xtab(5) =  - 0.171999257518648893241583152515D+01
    xtab(6) =  - 0.113611558521092066631913490556D+01
    xtab(7) =  - 0.565069583255575748526020337198D+00
    xtab(8) =    0.0D+00
    xtab(9) =    0.565069583255575748526020337198D+00
    xtab(10) =   0.113611558521092066631913490556D+01
    xtab(11) =   0.171999257518648893241583152515D+01
    xtab(12) =   0.232573248617385774545404479449D+01
    xtab(13) =   0.296716692790560324848896036355D+01
    xtab(14) =   0.366995037340445253472922383312D+01
    xtab(15) =   0.449999070730939155366438053053D+01

    weight(1) =  0.152247580425351702016062666965D-08
    weight(2) =  0.105911554771106663577520791055D-05
    weight(3) =  0.100004441232499868127296736177D-03
    weight(4) =  0.277806884291277589607887049229D-02
    weight(5) =  0.307800338725460822286814158758D-01
    weight(6) =  0.158488915795935746883839384960D+00
    weight(7) =  0.412028687498898627025891079568D+00
    weight(8) =  0.564100308726417532852625797340D+00
    weight(9) =  0.412028687498898627025891079568D+00
    weight(10) = 0.158488915795935746883839384960D+00
    weight(11) = 0.307800338725460822286814158758D-01
    weight(12) = 0.277806884291277589607887049229D-02
    weight(13) = 0.100004441232499868127296736177D-03
    weight(14) = 0.105911554771106663577520791055D-05
    weight(15) = 0.152247580425351702016062666965D-08

  else if ( norder == 16 ) then

    xtab(1) =  - 0.468873893930581836468849864875D+01
    xtab(2) =  - 0.386944790486012269871942409801D+01
    xtab(3) =  - 0.317699916197995602681399455926D+01
    xtab(4) =  - 0.254620215784748136215932870545D+01
    xtab(5) =  - 0.195178799091625397743465541496D+01
    xtab(6) =  - 0.138025853919888079637208966969D+01
    xtab(7) =  - 0.822951449144655892582454496734D+00
    xtab(8) =  - 0.273481046138152452158280401965D+00
    xtab(9) =    0.273481046138152452158280401965D+00
    xtab(10) =   0.822951449144655892582454496734D+00
    xtab(11) =   0.138025853919888079637208966969D+01
    xtab(12) =   0.195178799091625397743465541496D+01
    xtab(13) =   0.254620215784748136215932870545D+01
    xtab(14) =   0.317699916197995602681399455926D+01
    xtab(15) =   0.386944790486012269871942409801D+01
    xtab(16) =   0.468873893930581836468849864875D+01

    weight(1) =  0.265480747401118224470926366050D-09
    weight(2) =  0.232098084486521065338749423185D-06
    weight(3) =  0.271186009253788151201891432244D-04
    weight(4) =  0.932284008624180529914277305537D-03
    weight(5) =  0.128803115355099736834642999312D-01
    weight(6) =  0.838100413989858294154207349001D-01
    weight(7) =  0.280647458528533675369463335380D+00
    weight(8) =  0.507929479016613741913517341791D+00
    weight(9) =  0.507929479016613741913517341791D+00
    weight(10) = 0.280647458528533675369463335380D+00
    weight(11) = 0.838100413989858294154207349001D-01
    weight(12) = 0.128803115355099736834642999312D-01
    weight(13) = 0.932284008624180529914277305537D-03
    weight(14) = 0.271186009253788151201891432244D-04
    weight(15) = 0.232098084486521065338749423185D-06
    weight(16) = 0.265480747401118224470926366050D-09

  else if ( norder == 17 ) then

    xtab(1) =  - 0.487134519367440308834927655662D+01
    xtab(2) =  - 0.406194667587547430689245559698D+01
    xtab(3) =  - 0.337893209114149408338327069289D+01
    xtab(4) =  - 0.275776291570388873092640349574D+01
    xtab(5) =  - 0.217350282666662081927537907149D+01
    xtab(6) =  - 0.161292431422123133311288254454D+01
    xtab(7) =  - 0.106764872574345055363045773799D+01
    xtab(8) =  - 0.531633001342654731349086553718D+00
    xtab(9) =    0.0D+00
    xtab(10) =   0.531633001342654731349086553718D+00
    xtab(11) =   0.106764872574345055363045773799D+01
    xtab(12) =   0.161292431422123133311288254454D+01
    xtab(13) =   0.217350282666662081927537907149D+01
    xtab(14) =   0.275776291570388873092640349574D+01
    xtab(15) =   0.337893209114149408338327069289D+01
    xtab(16) =   0.406194667587547430689245559698D+01
    xtab(17) =   0.487134519367440308834927655662D+01

    weight(1) =  0.458057893079863330580889281222D-10
    weight(2) =  0.497707898163079405227863353715D-07
    weight(3) =  0.711228914002130958353327376218D-05
    weight(4) =  0.298643286697753041151336643059D-03
    weight(5) =  0.506734995762753791170069495879D-02
    weight(6) =  0.409200341495762798094994877854D-01
    weight(7) =  0.172648297670097079217645196219D+00
    weight(8) =  0.401826469470411956577635085257D+00
    weight(9) =  0.530917937624863560331883103379D+00
    weight(10) = 0.401826469470411956577635085257D+00
    weight(11) = 0.172648297670097079217645196219D+00
    weight(12) = 0.409200341495762798094994877854D-01
    weight(13) = 0.506734995762753791170069495879D-02
    weight(14) = 0.298643286697753041151336643059D-03
    weight(15) = 0.711228914002130958353327376218D-05
    weight(16) = 0.497707898163079405227863353715D-07
    weight(17) = 0.458057893079863330580889281222D-10

  else if ( norder == 18 ) then

    xtab(1) =  - 0.504836400887446676837203757885D+01
    xtab(2) =  - 0.424811787356812646302342016090D+01
    xtab(3) =  - 0.357376906848626607950067599377D+01
    xtab(4) =  - 0.296137750553160684477863254906D+01
    xtab(5) =  - 0.238629908916668600026459301424D+01
    xtab(6) =  - 0.183553160426162889225383944409D+01
    xtab(7) =  - 0.130092085838961736566626555439D+01
    xtab(8) =  - 0.776682919267411661316659462284D+00
    xtab(9) =  - 0.258267750519096759258116098711D+00
    xtab(10) =   0.258267750519096759258116098711D+00
    xtab(11) =   0.776682919267411661316659462284D+00
    xtab(12) =   0.130092085838961736566626555439D+01
    xtab(13) =   0.183553160426162889225383944409D+01
    xtab(14) =   0.238629908916668600026459301424D+01
    xtab(15) =   0.296137750553160684477863254906D+01
    xtab(16) =   0.357376906848626607950067599377D+01
    xtab(17) =   0.424811787356812646302342016090D+01
    xtab(18) =   0.504836400887446676837203757885D+01

    weight(1) =  0.782819977211589102925147471012D-11
    weight(2) =  0.104672057957920824443559608435D-07
    weight(3) =  0.181065448109343040959702385911D-05
    weight(4) =  0.918112686792940352914675407371D-04
    weight(5) =  0.188852263026841789438175325426D-02
    weight(6) =  0.186400423875446519219315221973D-01
    weight(7) =  0.973017476413154293308537234155D-01
    weight(8) =  0.284807285669979578595606820713D+00
    weight(9) =  0.483495694725455552876410522141D+00
    weight(10) = 0.483495694725455552876410522141D+00
    weight(11) = 0.284807285669979578595606820713D+00
    weight(12) = 0.973017476413154293308537234155D-01
    weight(13) = 0.186400423875446519219315221973D-01
    weight(14) = 0.188852263026841789438175325426D-02
    weight(15) = 0.918112686792940352914675407371D-04
    weight(16) = 0.181065448109343040959702385911D-05
    weight(17) = 0.104672057957920824443559608435D-07
    weight(18) = 0.782819977211589102925147471012D-11

  else if ( norder == 19 ) then

    xtab(1) =  - 0.522027169053748216460967142500D+01
    xtab(2) =  - 0.442853280660377943723498532226D+01
    xtab(3) =  - 0.376218735196402009751489394104D+01
    xtab(4) =  - 0.315784881834760228184318034120D+01
    xtab(5) =  - 0.259113378979454256492128084112D+01
    xtab(6) =  - 0.204923170985061937575050838669D+01
    xtab(7) =  - 0.152417061939353303183354859367D+01
    xtab(8) =  - 0.101036838713431135136859873726D+01
    xtab(9) =  - 0.503520163423888209373811765050D+00
    xtab(10) =   0.0D+00
    xtab(11) =   0.503520163423888209373811765050D+00
    xtab(12) =   0.101036838713431135136859873726D+01
    xtab(13) =   0.152417061939353303183354859367D+01
    xtab(14) =   0.204923170985061937575050838669D+01
    xtab(15) =   0.259113378979454256492128084112D+01
    xtab(16) =   0.315784881834760228184318034120D+01
    xtab(17) =   0.376218735196402009751489394104D+01
    xtab(18) =   0.442853280660377943723498532226D+01
    xtab(19) =   0.522027169053748216460967142500D+01

    weight(1) =  0.132629709449851575185289154385D-11
    weight(2) =  0.216305100986355475019693077221D-08
    weight(3) =  0.448824314722312295179447915594D-06
    weight(4) =  0.272091977631616257711941025214D-04
    weight(5) =  0.670877521407181106194696282100D-03
    weight(6) =  0.798886677772299020922211491861D-02
    weight(7) =  0.508103869090520673569908110358D-01
    weight(8) =  0.183632701306997074156148485766D+00
    weight(9) =  0.391608988613030244504042313621D+00
    weight(10) = 0.502974888276186530840731361096D+00
    weight(11) = 0.391608988613030244504042313621D+00
    weight(12) = 0.183632701306997074156148485766D+00
    weight(13) = 0.508103869090520673569908110358D-01
    weight(14) = 0.798886677772299020922211491861D-02
    weight(15) = 0.670877521407181106194696282100D-03
    weight(16) = 0.272091977631616257711941025214D-04
    weight(17) = 0.448824314722312295179447915594D-06
    weight(18) = 0.216305100986355475019693077221D-08
    weight(19) = 0.132629709449851575185289154385D-11

  else if ( norder == 20 ) then

    xtab(1) =  - 0.538748089001123286201690041068D+01
    xtab(2) =  - 0.460368244955074427307767524898D+01
    xtab(3) =  - 0.394476404011562521037562880052D+01
    xtab(4) =  - 0.334785456738321632691492452300D+01
    xtab(5) =  - 0.278880605842813048052503375640D+01
    xtab(6) =  - 0.225497400208927552308233334473D+01
    xtab(7) =  - 0.173853771211658620678086566214D+01
    xtab(8) =  - 0.123407621539532300788581834696D+01
    xtab(9) =  - 0.737473728545394358705605144252D+00
    xtab(10) = - 0.245340708300901249903836530634D+00
    xtab(11) =   0.245340708300901249903836530634D+00
    xtab(12) =   0.737473728545394358705605144252D+00
    xtab(13) =   0.123407621539532300788581834696D+01
    xtab(14) =   0.173853771211658620678086566214D+01
    xtab(15) =   0.225497400208927552308233334473D+01
    xtab(16) =   0.278880605842813048052503375640D+01
    xtab(17) =   0.334785456738321632691492452300D+01
    xtab(18) =   0.394476404011562521037562880052D+01
    xtab(19) =   0.460368244955074427307767524898D+01
    xtab(20) =   0.538748089001123286201690041068D+01

    weight(1) =  0.222939364553415129252250061603D-12
    weight(2) =  0.439934099227318055362885145547D-09
    weight(3) =  0.108606937076928169399952456345D-06
    weight(4) =  0.780255647853206369414599199965D-05
    weight(5) =  0.228338636016353967257145917963D-03
    weight(6) =  0.324377334223786183218324713235D-02
    weight(7) =  0.248105208874636108821649525589D-01
    weight(8) =  0.109017206020023320013755033535D+00
    weight(9) =  0.286675505362834129719659706228D+00
    weight(10) = 0.462243669600610089650328639861D+00
    weight(11) = 0.462243669600610089650328639861D+00
    weight(12) = 0.286675505362834129719659706228D+00
    weight(13) = 0.109017206020023320013755033535D+00
    weight(14) = 0.248105208874636108821649525589D-01
    weight(15) = 0.324377334223786183218324713235D-02
    weight(16) = 0.228338636016353967257145917963D-03
    weight(17) = 0.780255647853206369414599199965D-05
    weight(18) = 0.108606937076928169399952456345D-06
    weight(19) = 0.439934099227318055362885145547D-09
    weight(20) = 0.222939364553415129252250061603D-12

  else if ( norder == 30 ) then

    xtab( 1) =   -6.86334529352989158106110835756D+00
    xtab( 2) =   -6.13827922012393462039499237854D+00
    xtab( 3) =   -5.53314715156749572511833355558D+00
    xtab( 4) =   -4.98891896858994394448649710633D+00
    xtab( 5) =   -4.48305535709251834188703761971D+00
    xtab( 6) =   -4.00390860386122881522787601332D+00
    xtab( 7) =   -3.54444387315534988692540090217D+00
    xtab( 8) =   -3.09997052958644174868873332237D+00
    xtab( 9) =   -2.66713212453561720057110646422D+00
    xtab(10) =   -2.24339146776150407247297999483D+00
    xtab(11) =   -1.82674114360368803883588048351D+00
    xtab(12) =   -1.41552780019818851194072510555D+00
    xtab(13) =   -1.00833827104672346180498960870D+00
    xtab(14) =   -0.603921058625552307778155678757D+00
    xtab(15) =   -0.201128576548871485545763013244D+00
    xtab(16) =    0.201128576548871485545763013244D+00
    xtab(17) =    0.603921058625552307778155678757D+00
    xtab(18) =    1.00833827104672346180498960870D+00
    xtab(19) =    1.41552780019818851194072510555D+00
    xtab(20) =    1.82674114360368803883588048351D+00
    xtab(21) =    2.24339146776150407247297999483D+00
    xtab(22) =    2.66713212453561720057110646422D+00
    xtab(23) =    3.09997052958644174868873332237D+00
    xtab(24) =    3.54444387315534988692540090217D+00
    xtab(25) =    4.00390860386122881522787601332D+00
    xtab(26) =    4.48305535709251834188703761971D+00
    xtab(27) =    4.98891896858994394448649710633D+00
    xtab(28) =    5.53314715156749572511833355558D+00
    xtab(29) =    6.13827922012393462039499237854D+00
    xtab(30) =    6.86334529352989158106110835756D+00

    weight( 1) =   0.290825470013122622941102747365D-20
    weight( 2) =   0.281033360275090370876277491534D-16
    weight( 3) =   0.287860708054870606219239791142D-13
    weight( 4) =   0.810618629746304420399344796173D-11
    weight( 5) =   0.917858042437852820850075742492D-09
    weight( 6) =   0.510852245077594627738963204403D-07
    weight( 7) =   0.157909488732471028834638794022D-05
    weight( 8) =   0.293872522892298764150118423412D-04
    weight( 9) =   0.348310124318685523420995323183D-03
    weight(10) =   0.273792247306765846298942568953D-02
    weight(11) =   0.147038297048266835152773557787D-01
    weight(12) =   0.551441768702342511680754948183D-01
    weight(13) =   0.146735847540890099751693643152D+00
    weight(14) =   0.280130930839212667413493211293D+00
    weight(15) =   0.386394889541813862555601849165D+00
    weight(16) =   0.386394889541813862555601849165D+00
    weight(17) =   0.280130930839212667413493211293D+00
    weight(18) =   0.146735847540890099751693643152D+00
    weight(19) =   0.551441768702342511680754948183D-01
    weight(20) =   0.147038297048266835152773557787D-01
    weight(21) =   0.273792247306765846298942568953D-02
    weight(22) =   0.348310124318685523420995323183D-03
    weight(23) =   0.293872522892298764150118423412D-04
    weight(24) =   0.157909488732471028834638794022D-05
    weight(25) =   0.510852245077594627738963204403D-07
    weight(26) =   0.917858042437852820850075742492D-09
    weight(27) =   0.810618629746304420399344796173D-11
    weight(28) =   0.287860708054870606219239791142D-13
    weight(29) =   0.281033360275090370876277491534D-16
    weight(30) =   0.290825470013122622941102747365D-20

  else if ( norder == 32 ) then

    xtab( 1) =   -7.12581390983D+00
    xtab( 2) =   -6.40949814927D+00
    xtab( 3) =   -5.81222594952D+00
    xtab( 4) =   -5.27555098652D+00
    xtab( 5) =   -4.77716450350D+00
    xtab( 6) =   -4.30554795335D+00
    xtab( 7) =   -3.85375548547D+00
    xtab( 8) =   -3.41716749282D+00
    xtab( 9) =   -2.99249082500D+00
    xtab(10) =   -2.57724953773D+00
    xtab(11) =   -2.16949918361D+00
    xtab(12) =   -1.76765410946D+00
    xtab(13) =   -1.37037641095D+00
    xtab(14) =  -0.976500463590D+00
    xtab(15) =  -0.584978765436D+00
    xtab(16) =  -0.194840741569D+00
    xtab(17) =   0.194840741569D+00
    xtab(18) =   0.584978765436D+00
    xtab(19) =   0.976500463590D+00
    xtab(20) =    1.37037641095D+00
    xtab(21) =    1.76765410946D+00
    xtab(22) =    2.16949918361D+00
    xtab(23) =    2.57724953773D+00
    xtab(24) =    2.99249082500D+00
    xtab(25) =    3.41716749282D+00
    xtab(26) =    3.85375548547D+00
    xtab(27) =    4.30554795335D+00
    xtab(28) =    4.77716450350D+00
    xtab(29) =    5.27555098652D+00
    xtab(30) =    5.81222594952D+00
    xtab(31) =    6.40949814927D+00
    xtab(32) =    7.12581390983D+00

    weight( 1) =   0.731067642736D-22
    weight( 2) =   0.923173653649D-18
    weight( 3) =   0.119734401709D-14
    weight( 4) =   0.421501021125D-12
    weight( 5) =   0.593329146300D-10
    weight( 6) =   0.409883216476D-08
    weight( 7) =   0.157416779254D-06
    weight( 8) =   0.365058512955D-05
    weight( 9) =   0.541658406172D-04
    weight(10) =   0.536268365526D-03
    weight(11) =   0.365489032664D-02
    weight(12) =   0.175534288315D-01
    weight(13) =   0.604581309557D-01
    weight(14) =   0.151269734076D+00
    weight(15) =   0.277458142302D+00
    weight(16) =   0.375238352592D+00
    weight(17) =   0.375238352592D+00
    weight(18) =   0.277458142302D+00
    weight(19) =   0.151269734076D+00
    weight(20) =   0.604581309557D-01
    weight(21) =   0.175534288315D-01
    weight(22) =   0.365489032664D-02
    weight(23) =   0.536268365526D-03
    weight(24) =   0.541658406172D-04
    weight(25) =   0.365058512955D-05
    weight(26) =   0.157416779254D-06
    weight(27) =   0.409883216476D-08
    weight(28) =   0.593329146300D-10
    weight(29) =   0.421501021125D-12
    weight(30) =   0.119734401709D-14
    weight(31) =   0.923173653649D-18
    weight(32) =   0.731067642736D-22

  else if ( norder == 40 ) then

    xtab( 1) =   -8.09876113925D+00
    xtab( 2) =   -7.41158253149D+00
    xtab( 3) =   -6.84023730525D+00
    xtab( 4) =   -6.32825535122D+00
    xtab( 5) =   -5.85409505603D+00
    xtab( 6) =   -5.40665424797D+00
    xtab( 7) =   -4.97926097855D+00
    xtab( 8) =   -4.56750207284D+00
    xtab( 9) =   -4.16825706683D+00
    xtab(10) =   -3.77920675344D+00
    xtab(11) =   -3.39855826586D+00
    xtab(12) =   -3.02487988390D+00
    xtab(13) =   -2.65699599844D+00
    xtab(14) =   -2.29391714188D+00
    xtab(15) =   -1.93479147228D+00
    xtab(16) =   -1.57886989493D+00
    xtab(17) =   -1.22548010905D+00
    xtab(18) =  -0.874006612357D+00
    xtab(19) =  -0.523874713832D+00
    xtab(20) =  -0.174537214598D+00
    xtab(21) =   0.174537214598D+00
    xtab(22) =   0.523874713832D+00
    xtab(23) =   0.874006612357D+00
    xtab(24) =    1.22548010905D+00
    xtab(25) =    1.57886989493D+00
    xtab(26) =    1.93479147228D+00
    xtab(27) =    2.29391714188D+00
    xtab(28) =    2.65699599844D+00
    xtab(29) =    3.02487988390D+00
    xtab(30) =    3.39855826586D+00
    xtab(31) =    3.77920675344D+00
    xtab(32) =    4.16825706683D+00
    xtab(33) =    4.56750207284D+00
    xtab(34) =    4.97926097855D+00
    xtab(35) =    5.40665424797D+00
    xtab(36) =    5.85409505603D+00
    xtab(37) =    6.32825535122D+00
    xtab(38) =    6.84023730525D+00
    xtab(39) =    7.41158253149D+00
    xtab(40) =    8.09876113925D+00

    weight( 1) =   0.259104371384D-28
    weight( 2) =   0.854405696375D-24
    weight( 3) =   0.256759336540D-20
    weight( 4) =   0.198918101211D-17
    weight( 5) =   0.600835878947D-15
    weight( 6) =   0.880570764518D-13
    weight( 7) =   0.715652805267D-11
    weight( 8) =   0.352562079135D-09
    weight( 9) =   0.112123608322D-07
    weight(10) =   0.241114416359D-06
    weight(11) =   0.363157615067D-05
    weight(12) =   0.393693398108D-04
    weight(13) =   0.313853594540D-03
    weight(14) =   0.187149682959D-02
    weight(15) =   0.846088800823D-02
    weight(16) =   0.293125655361D-01
    weight(17) =   0.784746058652D-01
    weight(18) =   0.163378732713D+00
    weight(19) =   0.265728251876D+00
    weight(20) =   0.338643277425D+00
    weight(21) =   0.338643277425D+00
    weight(22) =   0.265728251876D+00
    weight(23) =   0.163378732713D+00
    weight(24) =   0.784746058652D-01
    weight(25) =   0.293125655361D-01
    weight(26) =   0.846088800823D-02
    weight(27) =   0.187149682959D-02
    weight(28) =   0.313853594540D-03
    weight(29) =   0.393693398108D-04
    weight(30) =   0.363157615067D-05
    weight(31) =   0.241114416359D-06
    weight(32) =   0.112123608322D-07
    weight(33) =   0.352562079135D-09
    weight(34) =   0.715652805267D-11
    weight(35) =   0.880570764518D-13
    weight(36) =   0.600835878947D-15
    weight(37) =   0.198918101211D-17
    weight(38) =   0.256759336540D-20
    weight(39) =   0.854405696375D-24
    weight(40) =   0.259104371384D-28

  else if ( norder == 50 ) then

    xtab( 1) =   -9.18240695813D+00
    xtab( 2) =   -8.52277103092D+00
    xtab( 3) =   -7.97562236821D+00
    xtab( 4) =   -7.48640942986D+00
    xtab( 5) =   -7.03432350977D+00
    xtab( 6) =   -6.60864797386D+00
    xtab( 7) =   -6.20295251927D+00
    xtab( 8) =   -5.81299467542D+00
    xtab( 9) =   -5.43578608722D+00
    xtab(10) =   -5.06911758492D+00
    xtab(11) =   -4.71129366617D+00
    xtab(12) =   -4.36097316045D+00
    xtab(13) =   -4.01706817286D+00
    xtab(14) =   -3.67867706252D+00
    xtab(15) =   -3.34503831394D+00
    xtab(16) =   -3.01549776957D+00
    xtab(17) =   -2.68948470227D+00
    xtab(18) =   -2.36649390430D+00
    xtab(19) =   -2.04607196869D+00
    xtab(20) =   -1.72780654752D+00
    xtab(21) =   -1.41131775490D+00
    xtab(22) =   -1.09625112896D+00
    xtab(23) =  -0.782271729555D+00
    xtab(24) =  -0.469059056678D+00
    xtab(25) =  -0.156302546889D+00
    xtab(26) =   0.156302546889D+00
    xtab(27) =   0.469059056678D+00
    xtab(28) =   0.782271729555D+00
    xtab(29) =    1.09625112896D+00
    xtab(30) =    1.41131775490D+00
    xtab(31) =    1.72780654752D+00
    xtab(32) =    2.04607196869D+00
    xtab(33) =    2.36649390430D+00
    xtab(34) =    2.68948470227D+00
    xtab(35) =    3.01549776957D+00
    xtab(36) =    3.34503831394D+00
    xtab(37) =    3.67867706252D+00
    xtab(38) =    4.01706817286D+00
    xtab(39) =    4.36097316045D+00
    xtab(40) =    4.71129366617D+00
    xtab(41) =    5.06911758492D+00
    xtab(42) =    5.43578608722D+00
    xtab(43) =    5.81299467542D+00
    xtab(44) =    6.20295251927D+00
    xtab(45) =    6.60864797386D+00
    xtab(46) =    7.03432350977D+00
    xtab(47) =    7.48640942986D+00
    xtab(48) =    7.97562236821D+00
    xtab(49) =    8.52277103092D+00
    xtab(50) =    9.18240695813D+00

    weight( 1) =   0.183379404857D-36
    weight( 2) =   0.167380166790D-31
    weight( 3) =   0.121524412340D-27
    weight( 4) =   0.213765830835D-24
    weight( 5) =   0.141709359957D-21
    weight( 6) =   0.447098436530D-19
    weight( 7) =   0.774238295702D-17
    weight( 8) =   0.809426189344D-15
    weight( 9) =   0.546594403180D-13
    weight(10) =   0.250665552389D-11
    weight(11) =   0.811187736448D-10
    weight(12) =   0.190904054379D-08
    weight(13) =   0.334679340401D-07
    weight(14) =   0.445702996680D-06
    weight(15) =   0.458168270794D-05
    weight(16) =   0.368401905377D-04
    weight(17) =   0.234269892109D-03
    weight(18) =   0.118901178175D-02
    weight(19) =   0.485326382616D-02
    weight(20) =   0.160319410684D-01
    weight(21) =   0.430791591566D-01
    weight(22) =   0.945489354768D-01
    weight(23) =   0.170032455676D+00
    weight(24) =   0.251130856331D+00
    weight(25) =   0.305085129203D+00
    weight(26) =   0.305085129203D+00
    weight(27) =   0.251130856331D+00
    weight(28) =   0.170032455676D+00
    weight(29) =   0.945489354768D-01
    weight(30) =   0.430791591566D-01
    weight(31) =   0.160319410684D-01
    weight(32) =   0.485326382616D-02
    weight(33) =   0.118901178175D-02
    weight(34) =   0.234269892109D-03
    weight(35) =   0.368401905377D-04
    weight(36) =   0.458168270794D-05
    weight(37) =   0.445702996680D-06
    weight(38) =   0.334679340401D-07
    weight(39) =   0.190904054379D-08
    weight(40) =   0.811187736448D-10
    weight(41) =   0.250665552389D-11
    weight(42) =   0.546594403180D-13
    weight(43) =   0.809426189344D-15
    weight(44) =   0.774238295702D-17
    weight(45) =   0.447098436530D-19
    weight(46) =   0.141709359957D-21
    weight(47) =   0.213765830835D-24
    weight(48) =   0.121524412340D-27
    weight(49) =   0.167380166790D-31
    weight(50) =   0.183379404857D-36

  else if ( norder == 60 ) then

    xtab( 1) =   -10.1591092462D+00
    xtab( 2) =   -9.52090367701D+00
    xtab( 3) =   -8.99239800140D+00
    xtab( 4) =   -8.52056928412D+00
    xtab( 5) =   -8.08518865425D+00
    xtab( 6) =   -7.67583993750D+00
    xtab( 7) =   -7.28627659440D+00
    xtab( 8) =   -6.91238153219D+00
    xtab( 9) =   -6.55125916706D+00
    xtab(10) =   -6.20077355799D+00
    xtab(11) =   -5.85929019639D+00
    xtab(12) =   -5.52552108614D+00
    xtab(13) =   -5.19842653458D+00
    xtab(14) =   -4.87715007747D+00
    xtab(15) =   -4.56097375794D+00
    xtab(16) =   -4.24928643596D+00
    xtab(17) =   -3.94156073393D+00
    xtab(18) =   -3.63733587617D+00
    xtab(19) =   -3.33620465355D+00
    xtab(20) =   -3.03780333823D+00
    xtab(21) =   -2.74180374807D+00
    xtab(22) =   -2.44790690231D+00
    xtab(23) =   -2.15583787123D+00
    xtab(24) =   -1.86534153123D+00
    xtab(25) =   -1.57617901198D+00
    xtab(26) =   -1.28812467487D+00
    xtab(27) =   -1.00096349956D+00
    xtab(28) =  -0.714488781673D+00
    xtab(29) =  -0.428500064221D+00
    xtab(30) =  -0.142801238703D+00
    xtab(31) =   0.142801238703D+00
    xtab(32) =   0.428500064221D+00
    xtab(33) =   0.714488781673D+00
    xtab(34) =    1.00096349956D+00
    xtab(35) =    1.28812467487D+00
    xtab(36) =    1.57617901198D+00
    xtab(37) =    1.86534153123D+00
    xtab(38) =    2.15583787123D+00
    xtab(39) =    2.44790690231D+00
    xtab(40) =    2.74180374807D+00
    xtab(41) =    3.03780333823D+00
    xtab(42) =    3.33620465355D+00
    xtab(43) =    3.63733587617D+00
    xtab(44) =    3.94156073393D+00
    xtab(45) =    4.24928643596D+00
    xtab(46) =    4.56097375794D+00
    xtab(47) =    4.87715007747D+00
    xtab(48) =    5.19842653458D+00
    xtab(49) =    5.52552108614D+00
    xtab(50) =    5.85929019639D+00
    xtab(51) =    6.20077355799D+00
    xtab(52) =    6.55125916706D+00
    xtab(53) =    6.91238153219D+00
    xtab(54) =    7.28627659440D+00
    xtab(55) =    7.67583993750D+00
    xtab(56) =    8.08518865425D+00
    xtab(57) =    8.52056928412D+00
    xtab(58) =    8.99239800140D+00
    xtab(59) =    9.52090367701D+00
    xtab(60) =    10.1591092462D+00

    weight( 1) =   0.110958724796D-44
    weight( 2) =   0.243974758810D-39
    weight( 3) =   0.377162672698D-35
    weight( 4) =   0.133255961176D-31
    weight( 5) =   0.171557314767D-28
    weight( 6) =   0.102940599693D-25
    weight( 7) =   0.334575695574D-23
    weight( 8) =   0.651256725748D-21
    weight( 9) =   0.815364047300D-19
    weight(10) =   0.692324790956D-17
    weight(11) =   0.415244410968D-15
    weight(12) =   0.181662457614D-13
    weight(13) =   0.594843051597D-12
    weight(14) =   0.148895734905D-10
    weight(15) =   0.289935901280D-09
    weight(16) =   0.445682277521D-08
    weight(17) =   0.547555461926D-07
    weight(18) =   0.543351613419D-06
    weight(19) =   0.439428693625D-05
    weight(20) =   0.291874190415D-04
    weight(21) =   0.160277334681D-03
    weight(22) =   0.731773556963D-03
    weight(23) =   0.279132482894D-02
    weight(24) =   0.893217836028D-02
    weight(25) =   0.240612727660D-01
    weight(26) =   0.547189709320D-01
    weight(27) =   0.105298763697D+00
    weight(28) =   0.171776156918D+00
    weight(29) =   0.237868904958D+00
    weight(30) =   0.279853117522D+00
    weight(31) =   0.279853117522D+00
    weight(32) =   0.237868904958D+00
    weight(33) =   0.171776156918D+00
    weight(34) =   0.105298763697D+00
    weight(35) =   0.547189709320D-01
    weight(36) =   0.240612727660D-01
    weight(37) =   0.893217836028D-02
    weight(38) =   0.279132482894D-02
    weight(39) =   0.731773556963D-03
    weight(40) =   0.160277334681D-03
    weight(41) =   0.291874190415D-04
    weight(42) =   0.439428693625D-05
    weight(43) =   0.543351613419D-06
    weight(44) =   0.547555461926D-07
    weight(45) =   0.445682277521D-08
    weight(46) =   0.289935901280D-09
    weight(47) =   0.148895734905D-10
    weight(48) =   0.594843051597D-12
    weight(49) =   0.181662457614D-13
    weight(50) =   0.415244410968D-15
    weight(51) =   0.692324790956D-17
    weight(52) =   0.815364047300D-19
    weight(53) =   0.651256725748D-21
    weight(54) =   0.334575695574D-23
    weight(55) =   0.102940599693D-25
    weight(56) =   0.171557314767D-28
    weight(57) =   0.133255961176D-31
    weight(58) =   0.377162672698D-35
    weight(59) =   0.243974758810D-39
    weight(60) =   0.110958724796D-44

  else if ( norder == 64 ) then

    xtab( 1) =   -10.5261231680D+00
    xtab( 2) =   -9.89528758683D+00
    xtab( 3) =   -9.37315954965D+00
    xtab( 4) =   -8.90724909996D+00
    xtab( 5) =   -8.47752908338D+00
    xtab( 6) =   -8.07368728501D+00
    xtab( 7) =   -7.68954016404D+00
    xtab( 8) =   -7.32101303278D+00
    xtab( 9) =   -6.96524112055D+00
    xtab(10) =   -6.62011226264D+00
    xtab(11) =   -6.28401122877D+00
    xtab(12) =   -5.95566632680D+00
    xtab(13) =   -5.63405216435D+00
    xtab(14) =   -5.31832522463D+00
    xtab(15) =   -5.00777960220D+00
    xtab(16) =   -4.70181564741D+00
    xtab(17) =   -4.39991716823D+00
    xtab(18) =   -4.10163447457D+00
    xtab(19) =   -3.80657151395D+00
    xtab(20) =   -3.51437593574D+00
    xtab(21) =   -3.22473129199D+00
    xtab(22) =   -2.93735082300D+00
    xtab(23) =   -2.65197243543D+00
    xtab(24) =   -2.36835458863D+00
    xtab(25) =   -2.08627287988D+00
    xtab(26) =   -1.80551717147D+00
    xtab(27) =   -1.52588914021D+00
    xtab(28) =   -1.24720015694D+00
    xtab(29) =  -0.969269423071D+00
    xtab(30) =  -0.691922305810D+00
    xtab(31) =  -0.414988824121D+00
    xtab(32) =  -0.138302244987D+00
    xtab(33) =   0.138302244987D+00
    xtab(34) =   0.414988824121D+00
    xtab(35) =   0.691922305810D+00
    xtab(36) =   0.969269423071D+00
    xtab(37) =    1.24720015694D+00
    xtab(38) =    1.52588914021D+00
    xtab(39) =    1.80551717147D+00
    xtab(40) =    2.08627287988D+00
    xtab(41) =    2.36835458863D+00
    xtab(42) =    2.65197243543D+00
    xtab(43) =    2.93735082300D+00
    xtab(44) =    3.22473129199D+00
    xtab(45) =    3.51437593574D+00
    xtab(46) =    3.80657151395D+00
    xtab(47) =    4.10163447457D+00
    xtab(48) =    4.39991716823D+00
    xtab(49) =    4.70181564741D+00
    xtab(50) =    5.00777960220D+00
    xtab(51) =    5.31832522463D+00
    xtab(52) =    5.63405216435D+00
    xtab(53) =    5.95566632680D+00
    xtab(54) =    6.28401122877D+00
    xtab(55) =    6.62011226264D+00
    xtab(56) =    6.96524112055D+00
    xtab(57) =    7.32101303278D+00
    xtab(58) =    7.68954016404D+00
    xtab(59) =    8.07368728501D+00
    xtab(60) =    8.47752908338D+00
    xtab(61) =    8.90724909996D+00
    xtab(62) =    9.37315954965D+00
    xtab(63) =    9.89528758683D+00
    xtab(64) =    10.5261231680D+00

    weight( 1) =   0.553570653584D-48
    weight( 2) =   0.167974799010D-42
    weight( 3) =   0.342113801099D-38
    weight( 4) =   0.155739062462D-34
    weight( 5) =   0.254966089910D-31
    weight( 6) =   0.192910359546D-28
    weight( 7) =   0.786179778889D-26
    weight( 8) =   0.191170688329D-23
    weight( 9) =   0.298286278427D-21
    weight(10) =   0.315225456649D-19
    weight(11) =   0.235188471067D-17
    weight(12) =   0.128009339117D-15
    weight(13) =   0.521862372645D-14
    weight(14) =   0.162834073070D-12
    weight(15) =   0.395917776693D-11
    weight(16) =   0.761521725012D-10
    weight(17) =   0.117361674232D-08
    weight(18) =   0.146512531647D-07
    weight(19) =   0.149553293672D-06
    weight(20) =   0.125834025103D-05
    weight(21) =   0.878849923082D-05
    weight(22) =   0.512592913577D-04
    weight(23) =   0.250983698512D-03
    weight(24) =   0.103632909950D-02
    weight(25) =   0.362258697852D-02
    weight(26) =   0.107560405098D-01
    weight(27) =   0.272031289536D-01
    weight(28) =   0.587399819634D-01
    weight(29) =   0.108498349306D+00
    weight(30) =   0.171685842349D+00
    weight(31) =   0.232994786062D+00
    weight(32) =   0.271377424940D+00
    weight(33) =   0.271377424940D+00
    weight(34) =   0.232994786062D+00
    weight(35) =   0.171685842349D+00
    weight(36) =   0.108498349306D+00
    weight(37) =   0.587399819634D-01
    weight(38) =   0.272031289536D-01
    weight(39) =   0.107560405098D-01
    weight(40) =   0.362258697852D-02
    weight(41) =   0.103632909950D-02
    weight(42) =   0.250983698512D-03
    weight(43) =   0.512592913577D-04
    weight(44) =   0.878849923082D-05
    weight(45) =   0.125834025103D-05
    weight(46) =   0.149553293672D-06
    weight(47) =   0.146512531647D-07
    weight(48) =   0.117361674232D-08
    weight(49) =   0.761521725012D-10
    weight(50) =   0.395917776693D-11
    weight(51) =   0.162834073070D-12
    weight(52) =   0.521862372645D-14
    weight(53) =   0.128009339117D-15
    weight(54) =   0.235188471067D-17
    weight(55) =   0.315225456649D-19
    weight(56) =   0.298286278427D-21
    weight(57) =   0.191170688329D-23
    weight(58) =   0.786179778889D-26
    weight(59) =   0.192910359546D-28
    weight(60) =   0.254966089910D-31
    weight(61) =   0.155739062462D-34
    weight(62) =   0.342113801099D-38
    weight(63) =   0.167974799010D-42
    weight(64) =   0.553570653584D-48

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 to 20,'
    write ( *, '(a)' ) '  30, 32, 40, 50, 60 and 64.'
    stop

  end if

  return
end
subroutine jacobi_com ( norder, xtab, weight, alpha, beta )
!
!*******************************************************************************
!
!! JACOBI_COM computes the abscissa and weights for Gauss-Jacobi quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) (1+X)**ALPHA * (1-X)**BETA * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    09 September 2000
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the quadrature rule to be computed.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the Gauss-Jacobi abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the Gauss-Jacobi weights.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) an
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) beta
  real ( kind = 8 ) bn
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) cc
  real ( kind = 8 ) delta
  real ( kind = 8 ) dp2
  real ( kind = 8 ) log_gamma
  integer i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) temp
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(norder)
!
!  Set the recursion coefficients.
!
  do i = 1, norder

    if ( alpha + beta == 0.0D+00 .or. beta - alpha == 0.0D+00) then

      b(i) = 0.0D+00

    else

      b(i) = ( alpha + beta ) * ( beta - alpha ) / &
            ( ( alpha + beta + real ( 2 * i, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) ) )

    end if

    if ( i == 1 ) then

      c(i) = 0.0D+00

    else

      c(i) = 4.0D+00 * real ( i - 1, kind = 8 ) &
            * ( alpha + real ( i - 1, kind = 8 ) ) &
            * ( beta + real ( i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( i - 1, kind = 8 ) ) / &
            ( ( alpha + beta + real ( 2 * i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) )**2 &
            * ( alpha + beta + real ( 2 * i - 3, kind = 8 ) ) )

    end if

  end do

  delta = exp ( log_gamma ( alpha + 1.0D+00 ) + log_gamma ( beta + 1.0D+00) &
    + log_gamma ( alpha + beta + 2.0D+00 ) )

  cc = delta * 2.0D+00**( alpha + beta + 1.0D+00 ) * product ( c(2:norder) )

  do i = 1, norder

    if ( i == 1 ) then

      an = alpha / real ( norder, kind = 8 )
      bn = beta / real ( norder, kind = 8 )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00+ real ( norder**2, kind = 8 ) ) &
        + 0.768D+00 * an / real ( norder, kind = 8 ) )

      r2 = 1.0D+00 + 1.48D+00 * an + 0.96D+00 * bn &
        + 0.452D+00 * an**2 + 0.83D+00 * an * bn

      x = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( norder, kind = 8 ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( norder, kind = 8 )

      r3 = 1.0D+00 + 0.012D+00 * beta * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( norder, kind = 8 )

      x = x - r1 * r2 * r3 * ( 1.0D+00 - x )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( norder, kind = 8 ) - 8.0D+00 ) &
        / real ( norder, kind = 8 )

      r3 = 1.0D+00 + 8.0D+00 * beta / &
        ( ( 6.28D+00 + beta ) * real ( norder**2, kind = 8 ) )

      x = x - r1 * r2 * r3 * ( xtab(1) - x )

    else if ( i < norder - 1 ) then

      x = 3.0D+00 * xtab(i-1) - 3.0D+00 * xtab(i-2) + xtab(i-3)

    else if ( i == norder - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * beta ) / ( 0.766D+00 + 0.119D+00 * beta )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( norder, kind = 8 ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( norder, kind = 8 ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( norder**2, kind = 8 ) ) )

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

    else if ( i == norder ) then

      r1 = ( 1.0D+00 + 0.37D+00 * beta ) / ( 1.67D+00 + 0.28D+00 * beta )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( norder, kind = 8 ) - 8.0D+00 ) &
        / real ( norder, kind = 8 ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( norder**2, kind = 8 ) ) )

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

    end if

    call jacobi_root ( x, norder, alpha, beta, dp2, p1, b, c )

    xtab(i) = x
    weight(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the order of the XTAB values.
!
  call dvec_reverse ( norder, xtab )

  return
end
subroutine jacobi_recur ( p2, dp2, p1, x, norder, alpha, beta, b, c )
!
!*******************************************************************************
!
!! JACOBI_RECUR finds the value and derivative of a Jacobi polynomial.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of J(NORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.
!
!    Input, real ( kind = 8 ) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0D+00 )
  dp2 = 1.0D+00

  do i = 2, norder

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine jacobi_root ( x, norder, alpha, beta, dp2, p1, b, c )
!
!*******************************************************************************
!
!! JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
!
!
!  Modified:
!
!    09 December 2000
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.
!
!    Output, real ( kind = 8 ) DP2, the value of J'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(NORDER), C(NORDER), the recursion coefficients.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  integer i
  integer, parameter :: maxstep = 10
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  eps = epsilon ( x )

  do i = 1, maxstep

    call jacobi_recur ( p2, dp2, p1, x, norder, alpha, beta, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine kronrod_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! KRONROD_SET sets abscissas and weights for Gauss-Kronrod quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Note:
!
!    A Kronrod rule is used in conjunction with a lower order
!    Gauss rule, and provides an efficient error estimation.
!
!    The error may be estimated as the difference in the two integral
!    approximations.
!
!    The efficiency comes about because the Kronrod uses the abscissas
!    of the Gauss rule, thus saving on the number of function evaluations
!    necessary.  If the Kronrod rule were replaced by a Gauss rule of
!    the same order, a higher precision integral estimate would be
!    made, but the function would have to be evaluated at many more
!    points.
!
!    The Gauss Kronrod pair of rules involves an ( NORDER + 1 ) / 2
!    point Gauss-Legendre rule and an NORDER point Kronrod rule.
!    Thus, the 15 point Kronrod rule should be paired with the
!    Gauss-Legendre 7 point rule.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner,
!    QUADPACK, A Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which may be
!    15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
!    order 7, 10, 15 or 20.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule, which
!    are symmetrically places in [-1,1].
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 15 ) then

    xtab(1) =  - 0.9914553711208126D+00
    xtab(2) =  - 0.9491079123427585D+00
    xtab(3) =  - 0.8648644233597691D+00
    xtab(4) =  - 0.7415311855993944D+00
    xtab(5) =  - 0.5860872354676911D+00
    xtab(6) =  - 0.4058451513773972D+00
    xtab(7) =  - 0.2077849550789850D+00
    xtab(8) =    0.0D+00
    xtab(9) =    0.2077849550789850D+00
    xtab(10) =   0.4058451513773972D+00
    xtab(11) =   0.5860872354676911D+00
    xtab(12) =   0.7415311855993944D+00
    xtab(13) =   0.8648644233597691D+00
    xtab(14) =   0.9491079123427585D+00
    xtab(15) =   0.9914553711208126D+00

    weight(1) =  0.2293532201052922D-01
    weight(2) =  0.6309209262997855D-01
    weight(3) =  0.1047900103222502D+00
    weight(4) =  0.1406532597155259D+00
    weight(5) =  0.1690047266392679D+00
    weight(6) =  0.1903505780647854D+00
    weight(7) =  0.2044329400752989D+00
    weight(8) =  0.2094821410847278D+00
    weight(9) =  0.2044329400752989D+00
    weight(10) = 0.1903505780647854D+00
    weight(11) = 0.1690047266392679D+00
    weight(12) = 0.1406532597155259D+00
    weight(13) = 0.1047900103222502D+00
    weight(14) = 0.6309209262997855D-01
    weight(15) = 0.2293532201052922D-01

  else if ( norder == 21 ) then

    xtab(1) =  - 0.9956571630258081D+00
    xtab(2) =  - 0.9739065285171717D+00
    xtab(3) =  - 0.9301574913557082D+00
    xtab(4) =  - 0.8650633666889845D+00
    xtab(5) =  - 0.7808177265864169D+00
    xtab(6) =  - 0.6794095682990244D+00
    xtab(7) =  - 0.5627571346686047D+00
    xtab(8) =  - 0.4333953941292472D+00
    xtab(9) =  - 0.2943928627014602D+00
    xtab(10) = - 0.1488743389816312D+00
    xtab(11) =   0.0D+00
    xtab(12) =   0.1488743389816312D+00
    xtab(13) =   0.2943928627014602D+00
    xtab(14) =   0.4333953941292472D+00
    xtab(15) =   0.5627571346686047D+00
    xtab(16) =   0.6794095682990244D+00
    xtab(17) =   0.7808177265864169D+00
    xtab(18) =   0.8650633666889845D+00
    xtab(19) =   0.9301574913557082D+00
    xtab(20) =   0.9739065285171717D+00
    xtab(21) =   0.9956571630258081D+00

    weight(1) =  0.1169463886737187D-01
    weight(2) =  0.3255816230796473D-01
    weight(3) =  0.5475589657435200D-01
    weight(4) =  0.7503967481091995D-01
    weight(5) =  0.9312545458369761D-01
    weight(6) =  0.1093871588022976D+00
    weight(7) =  0.1234919762620659D+00
    weight(8) =  0.1347092173114733D+00
    weight(9) =  0.1427759385770601D+00
    weight(10) = 0.1477391049013385D+00
    weight(11) = 0.1494455540029169D+00
    weight(12) = 0.1477391049013385D+00
    weight(13) = 0.1427759385770601D+00
    weight(14) = 0.1347092173114733D+00
    weight(15) = 0.1234919762620659D+00
    weight(16) = 0.1093871588022976D+00
    weight(17) = 0.9312545458369761D-01
    weight(18) = 0.7503967481091995D-01
    weight(19) = 0.5475589657435200D-01
    weight(20) = 0.3255816230796473D-01
    weight(21) = 0.1169463886737187D-01

  else if ( norder == 31 ) then

    xtab(1) =  - 0.9980022986933971D+00
    xtab(2) =  - 0.9879925180204854D+00
    xtab(3) =  - 0.9677390756791391D+00
    xtab(4) =  - 0.9372733924007059D+00
    xtab(5) =  - 0.8972645323440819D+00
    xtab(6) =  - 0.8482065834104272D+00
    xtab(7) =  - 0.7904185014424659D+00
    xtab(8) =  - 0.7244177313601700D+00
    xtab(9) =  - 0.6509967412974170D+00
    xtab(10) = - 0.5709721726085388D+00
    xtab(11) = - 0.4850818636402397D+00
    xtab(12) = - 0.3941513470775634D+00
    xtab(13) = - 0.2991800071531688D+00
    xtab(14) = - 0.2011940939974345D+00
    xtab(15) = - 0.1011420669187175D+00
    xtab(16) =   0.0D+00
    xtab(17) =   0.1011420669187175D+00
    xtab(18) =   0.2011940939974345D+00
    xtab(19) =   0.2991800071531688D+00
    xtab(20) =   0.3941513470775634D+00
    xtab(21) =   0.4850818636402397D+00
    xtab(22) =   0.5709721726085388D+00
    xtab(23) =   0.6509967412974170D+00
    xtab(24) =   0.7244177313601700D+00
    xtab(25) =   0.7904185014424659D+00
    xtab(26) =   0.8482065834104272D+00
    xtab(27) =   0.8972645323440819D+00
    xtab(28) =   0.9372733924007059D+00
    xtab(29) =   0.9677390756791391D+00
    xtab(30) =   0.9879925180204854D+00
    xtab(31) =   0.9980022986933971D+00

    weight(1) =  0.5377479872923349D-02
    weight(2) =  0.1500794732931612D-01
    weight(3) =  0.2546084732671532D-01
    weight(4) =  0.3534636079137585D-01
    weight(5) =  0.4458975132476488D-01
    weight(6) =  0.5348152469092809D-01
    weight(7) =  0.6200956780067064D-01
    weight(8) =  0.6985412131872826D-01
    weight(9) =  0.7684968075772038D-01
    weight(10) = 0.8308050282313302D-01
    weight(11) = 0.8856444305621177D-01
    weight(12) = 0.9312659817082532D-01
    weight(13) = 0.9664272698362368D-01
    weight(14) = 0.9917359872179196D-01
    weight(15) = 0.1007698455238756D+00
    weight(16) = 0.1013300070147915D+00
    weight(17) = 0.1007698455238756D+00
    weight(18) = 0.9917359872179196D-01
    weight(19) = 0.9664272698362368D-01
    weight(20) = 0.9312659817082532D-01
    weight(21) = 0.8856444305621177D-01
    weight(22) = 0.8308050282313302D-01
    weight(23) = 0.7684968075772038D-01
    weight(24) = 0.6985412131872826D-01
    weight(25) = 0.6200956780067064D-01
    weight(26) = 0.5348152469092809D-01
    weight(27) = 0.4458975132476488D-01
    weight(28) = 0.3534636079137585D-01
    weight(29) = 0.2546084732671532D-01
    weight(30) = 0.1500794732931612D-01
    weight(31) = 0.5377479872923349D-02

  else if ( norder == 41 ) then

    xtab(1) =  - 0.9988590315882777D+00
    xtab(2) =  - 0.9931285991850949D+00
    xtab(3) =  - 0.9815078774502503D+00
    xtab(4) =  - 0.9639719272779138D+00
    xtab(5) =  - 0.9408226338317548D+00
    xtab(6) =  - 0.9122344282513259D+00
    xtab(7) =  - 0.8782768112522820D+00
    xtab(8) =  - 0.8391169718222188D+00
    xtab(9) =  - 0.7950414288375512D+00
    xtab(10) = - 0.7463319064601508D+00
    xtab(11) = - 0.6932376563347514D+00
    xtab(12) = - 0.6360536807265150D+00
    xtab(13) = - 0.5751404468197103D+00
    xtab(14) = - 0.5108670019508271D+00
    xtab(15) = - 0.4435931752387251D+00
    xtab(16) = - 0.3737060887154196D+00
    xtab(17) = - 0.3016278681149130D+00
    xtab(18) = - 0.2277858511416451D+00
    xtab(19) = - 0.1526054652409227D+00
    xtab(20) = - 0.7652652113349733D-01
    xtab(21) =   0.0D+00
    xtab(22) =   0.7652652113349733D-01
    xtab(23) =   0.1526054652409227D+00
    xtab(24) =   0.2277858511416451D+00
    xtab(25) =   0.3016278681149130D+00
    xtab(26) =   0.3737060887154196D+00
    xtab(27) =   0.4435931752387251D+00
    xtab(28) =   0.5108670019508271D+00
    xtab(29) =   0.5751404468197103D+00
    xtab(30) =   0.6360536807265150D+00
    xtab(31) =   0.6932376563347514D+00
    xtab(32) =   0.7463319064601508D+00
    xtab(33) =   0.7950414288375512D+00
    xtab(34) =   0.8391169718222188D+00
    xtab(35) =   0.8782768112522820D+00
    xtab(36) =   0.9122344282513259D+00
    xtab(37) =   0.9408226338317548D+00
    xtab(38) =   0.9639719272779138D+00
    xtab(39) =   0.9815078774502503D+00
    xtab(40) =   0.9931285991850949D+00
    xtab(41) =   0.9988590315882777D+00

    weight(1) =  0.3073583718520532D-02
    weight(2) =  0.8600269855642942D-02
    weight(3) =  0.1462616925697125D-01
    weight(4) =  0.2038837346126652D-01
    weight(5) =  0.2588213360495116D-01
    weight(6) =  0.3128730677703280D-01
    weight(7) =  0.3660016975820080D-01
    weight(8) =  0.4166887332797369D-01
    weight(9) =  0.4643482186749767D-01
    weight(10) = 0.5094457392372869D-01
    weight(11) = 0.5519510534828599D-01
    weight(12) = 0.5911140088063957D-01
    weight(13) = 0.6265323755478117D-01
    weight(14) = 0.6583459713361842D-01
    weight(15) = 0.6864867292852162D-01
    weight(16) = 0.7105442355344407D-01
    weight(17) = 0.7303069033278667D-01
    weight(18) = 0.7458287540049919D-01
    weight(19) = 0.7570449768455667D-01
    weight(20) = 0.7637786767208074D-01
    weight(21) = 0.7660071191799966D-01
    weight(22) = 0.7637786767208074D-01
    weight(23) = 0.7570449768455667D-01
    weight(24) = 0.7458287540049919D-01
    weight(25) = 0.7303069033278667D-01
    weight(26) = 0.7105442355344407D-01
    weight(27) = 0.6864867292852162D-01
    weight(28) = 0.6583459713361842D-01
    weight(29) = 0.6265323755478117D-01
    weight(30) = 0.5911140088063957D-01
    weight(31) = 0.5519510534828599D-01
    weight(32) = 0.5094457392372869D-01
    weight(33) = 0.4643482186749767D-01
    weight(34) = 0.4166887332797369D-01
    weight(35) = 0.3660016975820080D-01
    weight(36) = 0.3128730677703280D-01
    weight(37) = 0.2588213360495116D-01
    weight(38) = 0.2038837346126652D-01
    weight(39) = 0.1462616925697125D-01
    weight(40) = 0.8600269855642942D-02
    weight(41) = 0.3073583718520532D-02

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KRONROD_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 15, 21, 31 or 41.'
    stop

  end if

  return
end
subroutine laguerre_com ( norder, xtab, weight, alpha )
!
!*******************************************************************************
!
!! LAGUERRE_COM computes the abscissa and weights for Gauss-Laguerre quadrature.
!
!
!  Discussion:
!
!    In the simplest case, ALPHA is 0, and we are approximating the
!    integral from 0 to INFINITY of EXP(-X) * F(X).  When this is so,
!    it is easy to modify the rule to approximate the integral from
!    A to INFINITY as well.
!
!    If ALPHA is nonzero, then there is no simple way to extend the
!    rule to approximate the integral from A to INFINITY.  The simplest
!    procedures would be to approximate the integral from 0 to A.
!
!  Integration interval:
!
!      [ A, +Infinity ) or [ 0, +Infinity )
!
!  Weight function:
!
!      EXP ( - X ) or EXP ( - X ) * X**ALPHA
!
!  Integral to approximate:
!
!      Integral ( A <= X < +INFINITY ) EXP ( - X ) * F(X) dX
!    or
!      Integral ( 0 <= X < +INFINITY ) EXP ( - X ) * X**ALPHA * F(X) dX
!
!  Approximate integral:
!
!      EXP ( - A ) * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( A+XTAB(I) )
!    or
!      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the quadrature rule to be computed.
!    NORDER must be at least 1.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the Gauss-Laguerre abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the Gauss-Laguerre weights.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    Set ALPHA = 0.0D+00 for the simplest rule.
!    ALPHA must be nonnegative.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  real ( kind = 8 ) gamma
  integer i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) ratio
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(norder)
!
!  Set the recursion coefficients.
!
  do i = 1, norder
    b(i) = ( alpha + real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, norder
    c(i) = real ( i - 1, kind = 8 ) * ( alpha + real ( i - 1, kind = 8 ) )
  end do

  cc = gamma ( alpha + 1.0D+00 ) * product ( c(2:norder) )

  do i = 1, norder
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      x = ( 1.0D+00 + alpha ) * ( 3.0D+00+ 0.92 * alpha ) / &
        ( 1.0D+00 + 2.4D+00 * real ( norder, kind = 8 ) + 1.8D+00 * alpha )

    else if ( i == 2 ) then

      x = x + ( 15.0D+00 + 6.25D+00 * alpha ) / &
        ( 1.0D+00 + 0.9D+00 * alpha + 2.5D+00 * real ( norder, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      r2 = 1.26D+00 * real ( i - 2, kind = 8 ) * alpha / &
        ( 1.0D+00 + 3.5D+00 * real ( i - 2, kind = 8 ) )

      ratio = ( r1 + r2 ) / ( 1.0D+00 + 0.3D+00 * alpha )

      x = x + ratio * ( x - xtab(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call laguerre_root ( x, norder, alpha, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    xtab(i) = x
    weight(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine laguerre_recur ( p2, dp2, p1, x, norder, alpha, b, c )
!
!*******************************************************************************
!
!! LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(NORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Input, real ( kind = 8 ) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - alpha - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, norder

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine laguerre_root ( x, norder, alpha, dp2, p1, b, c )
!
!*******************************************************************************
!
!! LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
!
!
!  Modified:
!
!    09 December 2000
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(NORDER), C(NORDER), the recursion coefficients.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) c(norder)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  integer i
  integer, parameter :: maxstep = 10
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  eps = epsilon ( x )

  do i = 1, maxstep

    call laguerre_recur ( p2, dp2, p1, x, norder, alpha, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine laguerre_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LAGUERRE_SET sets abscissas and weights for Laguerre quadrature.
!
!
!  Integration interval:
!
!    [ 0, +Infinity )
!
!  Weight function:
!
!    EXP ( - X )
!
!  Integral to approximate:
!
!    Integral ( 0 <= X < +INFINITY ) EXP ( - X ) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Note:
!
!    The abscissas are the zeroes of the Laguerre polynomial L(NORDER)(X).
!
!  Modified:
!
!    17 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, and should add to 1.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 1 ) then

    xtab(1) =    1.0D+00
    weight(1) =  1.0D+00

  else if ( norder == 2 ) then

    xtab(1) =    0.585786437626904951198311275790D+00
    xtab(2) =    0.341421356237309504880168872421D+01

    weight(1) =  0.853553390593273762200422181052D+00
    weight(2) =  0.146446609406726237799577818948D+00

  else if ( norder == 3 ) then

    xtab(1) =    0.415774556783479083311533873128D+00
    xtab(2) =    0.229428036027904171982205036136D+01
    xtab(3) =    0.628994508293747919686641576551D+01

    weight(1) =  0.711093009929173015449590191143D+00
    weight(2) =  0.278517733569240848801444888457D+00
    weight(3) =  0.103892565015861357489649204007D-01

  else if ( norder == 4 ) then

    xtab(1) =    0.322547689619392311800361943361D+00
    xtab(2) =    0.174576110115834657568681671252D+01
    xtab(3) =    0.453662029692112798327928538496D+01
    xtab(4) =    0.939507091230113312923353644342D+01

    weight(1) =  0.603154104341633601635966023818D+00
    weight(2) =  0.357418692437799686641492017458D+00
    weight(3) =  0.388879085150053842724381681562D-01
    weight(4) =  0.539294705561327450103790567621D-03

  else if ( norder == 5 ) then

    xtab(1) =    0.263560319718140910203061943361D+00
    xtab(2) =    0.141340305910651679221840798019D+01
    xtab(3) =    0.359642577104072208122318658878D+01
    xtab(4) =    0.708581000585883755692212418111D+01
    xtab(5) =    0.126408008442757826594332193066D+02

    weight(1) =  0.521755610582808652475860928792D+00
    weight(2) =  0.398666811083175927454133348144D+00
    weight(3) =  0.759424496817075953876533114055D-01
    weight(4) =  0.361175867992204845446126257304D-02
    weight(5) =  0.233699723857762278911490845516D-04

  else if ( norder == 6 ) then

    xtab(1) =    0.222846604179260689464354826787D+00
    xtab(2) =    0.118893210167262303074315092194D+01
    xtab(3) =    0.299273632605931407769132528451D+01
    xtab(4) =    0.577514356910451050183983036943D+01
    xtab(5) =    0.983746741838258991771554702994D+01
    xtab(6) =    0.159828739806017017825457915674D+02

    weight(1) =  0.458964673949963593568284877709D+00
    weight(2) =  0.417000830772120994113377566193D+00
    weight(3) =  0.113373382074044975738706185098D+00
    weight(4) =  0.103991974531490748989133028469D-01
    weight(5) =  0.261017202814932059479242860001D-03
    weight(6) =  0.898547906429621238825292052825D-06

  else if ( norder == 7 ) then

    xtab(1) =    0.193043676560362413838247885004D+00
    xtab(2) =    0.102666489533919195034519944317D+01
    xtab(3) =    0.256787674495074620690778622666D+01
    xtab(4) =    0.490035308452648456810171437810D+01
    xtab(5) =    0.818215344456286079108182755123D+01
    xtab(6) =    0.127341802917978137580126424582D+02
    xtab(7) =    0.193957278622625403117125820576D+02

    weight(1) =  0.409318951701273902130432880018D+00
    weight(2) =  0.421831277861719779929281005417D+00
    weight(3) =  0.147126348657505278395374184637D+00
    weight(4) =  0.206335144687169398657056149642D-01
    weight(5) =  0.107401014328074552213195962843D-02
    weight(6) =  0.158654643485642012687326223234D-04
    weight(7) =  0.317031547899558056227132215385D-07

  else if ( norder == 8 ) then

    xtab(1) =    0.170279632305100999788861856608D+00
    xtab(2) =    0.903701776799379912186020223555D+00
    xtab(3) =    0.225108662986613068930711836697D+01
    xtab(4) =    0.426670017028765879364942182690D+01
    xtab(5) =    0.704590540239346569727932548212D+01
    xtab(6) =    0.107585160101809952240599567880D+02
    xtab(7) =    0.157406786412780045780287611584D+02
    xtab(8) =    0.228631317368892641057005342974D+02

    weight(1) =  0.369188589341637529920582839376D+00
    weight(2) =  0.418786780814342956076978581333D+00
    weight(3) =  0.175794986637171805699659866777D+00
    weight(4) =  0.333434922612156515221325349344D-01
    weight(5) =  0.279453623522567252493892414793D-02
    weight(6) =  0.907650877335821310423850149336D-04
    weight(7) =  0.848574671627253154486801830893D-06
    weight(8) =  0.104800117487151038161508853552D-08

  else if ( norder == 9 ) then

    xtab(1) =    0.152322227731808247428107073127D+00
    xtab(2) =    0.807220022742255847741419210952D+00
    xtab(3) =    0.200513515561934712298303324701D+01
    xtab(4) =    0.378347397333123299167540609364D+01
    xtab(5) =    0.620495677787661260697353521006D+01
    xtab(6) =    0.937298525168757620180971073215D+01
    xtab(7) =    0.134662369110920935710978818397D+02
    xtab(8) =    0.188335977889916966141498992996D+02
    xtab(9) =    0.263740718909273767961410072937D+02

    weight(1) =  0.336126421797962519673467717606D+00
    weight(2) =  0.411213980423984387309146942793D+00
    weight(3) =  0.199287525370885580860575607212D+00
    weight(4) =  0.474605627656515992621163600479D-01
    weight(5) =  0.559962661079458317700419900556D-02
    weight(6) =  0.305249767093210566305412824291D-03
    weight(7) =  0.659212302607535239225572284875D-05
    weight(8) =  0.411076933034954844290241040330D-07
    weight(9) =  0.329087403035070757646681380323D-10

  else if ( norder == 10 ) then

    xtab(1) =    0.137793470540492430830772505653D+00
    xtab(2) =    0.729454549503170498160373121676D+00
    xtab(3) =    0.180834290174031604823292007575D+01
    xtab(4) =    0.340143369785489951448253222141D+01
    xtab(5) =    0.555249614006380363241755848687D+01
    xtab(6) =    0.833015274676449670023876719727D+01
    xtab(7) =    0.118437858379000655649185389191D+02
    xtab(8) =    0.162792578313781020995326539358D+02
    xtab(9) =    0.219965858119807619512770901956D+02
    xtab(10) =   0.299206970122738915599087933408D+02

    weight(1) =  0.308441115765020141547470834678D+00
    weight(2) =  0.401119929155273551515780309913D+00
    weight(3) =  0.218068287611809421588648523475D+00
    weight(4) =  0.620874560986777473929021293135D-01
    weight(5) =  0.950151697518110055383907219417D-02
    weight(6) =  0.753008388587538775455964353676D-03
    weight(7) =  0.282592334959956556742256382685D-04
    weight(8) =  0.424931398496268637258657665975D-06
    weight(9) =  0.183956482397963078092153522436D-08
    weight(10) = 0.991182721960900855837754728324D-12

  else if ( norder == 11 ) then

    xtab(1) =    0.125796442187967522675794577516D+00
    xtab(2) =    0.665418255839227841678127839420D+00
    xtab(3) =    0.164715054587216930958700321365D+01
    xtab(4) =    0.309113814303525495330195934259D+01
    xtab(5) =    0.502928440157983321236999508366D+01
    xtab(6) =    0.750988786380661681941099714450D+01
    xtab(7) =    0.106059509995469677805559216457D+02
    xtab(8) =    0.144316137580641855353200450349D+02
    xtab(9) =    0.191788574032146786478174853989D+02
    xtab(10) =   0.252177093396775611040909447797D+02
    xtab(11) =   0.334971928471755372731917259395D+02

    weight(1) =  0.284933212894200605056051024724D+00
    weight(2) =  0.389720889527849377937553508048D+00
    weight(3) =  0.232781831848991333940223795543D+00
    weight(4) =  0.765644535461966864008541790132D-01
    weight(5) =  0.143932827673506950918639187409D-01
    weight(6) =  0.151888084648487306984777640042D-02
    weight(7) =  0.851312243547192259720424170600D-04
    weight(8) =  0.229240387957450407857683270709D-05
    weight(9) =  0.248635370276779587373391491114D-07
    weight(10) = 0.771262693369132047028152590222D-10
    weight(11) = 0.288377586832362386159777761217D-13

  else if ( norder == 12 ) then

    xtab(1) =    0.115722117358020675267196428240D+00
    xtab(2) =    0.611757484515130665391630053042D+00
    xtab(3) =    0.151261026977641878678173792687D+01
    xtab(4) =    0.283375133774350722862747177657D+01
    xtab(5) =    0.459922763941834848460572922485D+01
    xtab(6) =    0.684452545311517734775433041849D+01
    xtab(7) =    0.962131684245686704391238234923D+01
    xtab(8) =    0.130060549933063477203460524294D+02
    xtab(9) =    0.171168551874622557281840528008D+02
    xtab(10) =   0.221510903793970056699218950837D+02
    xtab(11) =   0.284879672509840003125686072325D+02
    xtab(12) =   0.370991210444669203366389142764D+02

    weight(1) =  0.264731371055443190349738892056D+00
    weight(2) =  0.377759275873137982024490556707D+00
    weight(3) =  0.244082011319877564254870818274D+00
    weight(4) =  0.904492222116809307275054934667D-01
    weight(5) =  0.201023811546340965226612867827D-01
    weight(6) =  0.266397354186531588105415760678D-02
    weight(7) =  0.203231592662999392121432860438D-03
    weight(8) =  0.836505585681979874533632766396D-05
    weight(9) =  0.166849387654091026116989532619D-06
    weight(10) = 0.134239103051500414552392025055D-08
    weight(11) = 0.306160163503502078142407718971D-11
    weight(12) = 0.814807746742624168247311868103D-15

  else if ( norder == 13 ) then

    xtab(1) =    0.107142388472252310648493376977D+00
    xtab(2) =    0.566131899040401853406036347177D+00
    xtab(3) =    0.139856433645101971792750259921D+01
    xtab(4) =    0.261659710840641129808364008472D+01
    xtab(5) =    0.423884592901703327937303389926D+01
    xtab(6) =    0.629225627114007378039376523025D+01
    xtab(7) =    0.881500194118697804733348868036D+01
    xtab(8) =    0.118614035888112425762212021880D+02
    xtab(9) =    0.155107620377037527818478532958D+02
    xtab(10) =   0.198846356638802283332036594634D+02
    xtab(11) =   0.251852638646777580842970297823D+02
    xtab(12) =   0.318003863019472683713663283526D+02
    xtab(13) =   0.407230086692655795658979667001D+02

    weight(1) =  0.247188708429962621346249185964D+00
    weight(2) =  0.365688822900521945306717530893D+00
    weight(3) =  0.252562420057658502356824288815D+00
    weight(4) =  0.103470758024183705114218631672D+00
    weight(5) =  0.264327544155616157781587735702D-01
    weight(6) =  0.422039604025475276555209292644D-02
    weight(7) =  0.411881770472734774892472527082D-03
    weight(8) =  0.235154739815532386882897300772D-04
    weight(9) =  0.731731162024909910401047197761D-06
    weight(10) = 0.110884162570398067979150974759D-07
    weight(11) = 0.677082669220589884064621459082D-10
    weight(12) = 0.115997995990507606094507145382D-12
    weight(13) = 0.224509320389275841599187226865D-16

  else if ( norder == 14 ) then

    xtab(1) =    0.997475070325975745736829452514D-01
    xtab(2) =    0.526857648851902896404583451502D+00
    xtab(3) =    0.130062912125149648170842022116D+01
    xtab(4) =    0.243080107873084463616999751038D+01
    xtab(5) =    0.393210282229321888213134366778D+01
    xtab(6) =    0.582553621830170841933899983898D+01
    xtab(7) =    0.814024014156514503005978046052D+01
    xtab(8) =    0.109164995073660188408130510904D+02
    xtab(9) =    0.142108050111612886831059780825D+02
    xtab(10) =   0.181048922202180984125546272083D+02
    xtab(11) =   0.227233816282696248232280886985D+02
    xtab(12) =   0.282729817232482056954158923218D+02
    xtab(13) =   0.351494436605924265828643121364D+02
    xtab(14) =   0.443660817111174230416312423666D+02

    weight(1) =  0.231815577144864977840774861104D+00
    weight(2) =  0.353784691597543151802331301273D+00
    weight(3) =  0.258734610245428085987320561144D+00
    weight(4) =  0.115482893556923210087304988673D+00
    weight(5) =  0.331920921593373600387499587137D-01
    weight(6) =  0.619286943700661021678785967675D-02
    weight(7) =  0.739890377867385942425890907080D-03
    weight(8) =  0.549071946684169837857331777667D-04
    weight(9) =  0.240958576408537749675775256553D-05
    weight(10) = 0.580154398167649518088619303904D-07
    weight(11) = 0.681931469248497411961562387084D-09
    weight(12) = 0.322120775189484793980885399656D-11
    weight(13) = 0.422135244051658735159797335643D-14
    weight(14) = 0.605237502228918880839870806281D-18

  else if ( norder == 15 ) then

    xtab(1) =    0.933078120172818047629030383672D-01
    xtab(2) =    0.492691740301883908960101791412D+00
    xtab(3) =    0.121559541207094946372992716488D+01
    xtab(4) =    0.226994952620374320247421741375D+01
    xtab(5) =    0.366762272175143727724905959436D+01
    xtab(6) =    0.542533662741355316534358132596D+01
    xtab(7) =    0.756591622661306786049739555812D+01
    xtab(8) =    0.101202285680191127347927394568D+02
    xtab(9) =    0.131302824821757235640991204176D+02
    xtab(10) =   0.166544077083299578225202408430D+02
    xtab(11) =   0.207764788994487667729157175676D+02
    xtab(12) =   0.256238942267287801445868285977D+02
    xtab(13) =   0.314075191697539385152432196202D+02
    xtab(14) =   0.385306833064860094162515167595D+02
    xtab(15) =   0.480260855726857943465734308508D+02

    weight(1) =  0.218234885940086889856413236448D+00
    weight(2) =  0.342210177922883329638948956807D+00
    weight(3) =  0.263027577941680097414812275022D+00
    weight(4) =  0.126425818105930535843030549378D+00
    weight(5) =  0.402068649210009148415854789871D-01
    weight(6) =  0.856387780361183836391575987649D-02
    weight(7) =  0.121243614721425207621920522467D-02
    weight(8) =  0.111674392344251941992578595518D-03
    weight(9) =  0.645992676202290092465319025312D-05
    weight(10) = 0.222631690709627263033182809179D-06
    weight(11) = 0.422743038497936500735127949331D-08
    weight(12) = 0.392189726704108929038460981949D-10
    weight(13) = 0.145651526407312640633273963455D-12
    weight(14) = 0.148302705111330133546164737187D-15
    weight(15) = 0.160059490621113323104997812370D-19

  else if ( norder == 16 ) then

    xtab(1) =    0.876494104789278403601980973401D-01
    xtab(2) =    0.462696328915080831880838260664D+00
    xtab(3) =    0.114105777483122685687794501811D+01
    xtab(4) =    0.212928364509838061632615907066D+01
    xtab(5) =    0.343708663389320664523510701675D+01
    xtab(6) =    0.507801861454976791292305830814D+01
    xtab(7) =    0.707033853504823413039598947080D+01
    xtab(8) =    0.943831433639193878394724672911D+01
    xtab(9) =    0.122142233688661587369391246088D+02
    xtab(10) =   0.154415273687816170767647741622D+02
    xtab(11) =   0.191801568567531348546631409497D+02
    xtab(12) =   0.235159056939919085318231872752D+02
    xtab(13) =   0.285787297428821403675206137099D+02
    xtab(14) =   0.345833987022866258145276871778D+02
    xtab(15) =   0.419404526476883326354722330252D+02
    xtab(16) =   0.517011603395433183643426971197D+02

    weight(1) =  0.206151714957800994334273636741D+00
    weight(2) =  0.331057854950884165992983098710D+00
    weight(3) =  0.265795777644214152599502020650D+00
    weight(4) =  0.136296934296377539975547513526D+00
    weight(5) =  0.473289286941252189780623392781D-01
    weight(6) =  0.112999000803394532312490459701D-01
    weight(7) =  0.184907094352631086429176783252D-02
    weight(8) =  0.204271915308278460126018338421D-03
    weight(9) =  0.148445868739812987713515067551D-04
    weight(10) = 0.682831933087119956439559590327D-06
    weight(11) = 0.188102484107967321388159920418D-07
    weight(12) = 0.286235024297388161963062629156D-09
    weight(13) = 0.212707903322410296739033610978D-11
    weight(14) = 0.629796700251786778717446214552D-14
    weight(15) = 0.505047370003551282040213233303D-17
    weight(16) = 0.416146237037285519042648356116D-21

  else if ( norder == 17 ) then

    xtab(1) =    0.826382147089476690543986151980D-01
    xtab(2) =    0.436150323558710436375959029847D+00
    xtab(3) =    0.107517657751142857732980316755D+01
    xtab(4) =    0.200519353164923224070293371933D+01
    xtab(5) =    0.323425612404744376157380120696D+01
    xtab(6) =    0.477351351370019726480932076262D+01
    xtab(7) =    0.663782920536495266541643929703D+01
    xtab(8) =    0.884668551116980005369470571184D+01
    xtab(9) =    0.114255293193733525869726151469D+02
    xtab(10) =   0.144078230374813180021982874959D+02
    xtab(11) =   0.178382847307011409290658752412D+02
    xtab(12) =   0.217782682577222653261749080522D+02
    xtab(13) =   0.263153178112487997766149598369D+02
    xtab(14) =   0.315817716804567331343908517497D+02
    xtab(15) =   0.377960938374771007286092846663D+02
    xtab(16) =   0.453757165339889661829258363215D+02
    xtab(17) =   0.553897517898396106640900199790D+02

    weight(1) =  0.195332205251770832145927297697D+00
    weight(2) =  0.320375357274540281336625631970D+00
    weight(3) =  0.267329726357171097238809604160D+00
    weight(4) =  0.145129854358758625407426447473D+00
    weight(5) =  0.544369432453384577793805803066D-01
    weight(6) =  0.143572977660618672917767247431D-01
    weight(7) =  0.266282473557277256843236250006D-02
    weight(8) =  0.343679727156299920611775097985D-03
    weight(9) =  0.302755178378287010943703641131D-04
    weight(10) = 0.176851505323167689538081156159D-05
    weight(11) = 0.657627288681043332199222748162D-07
    weight(12) = 0.146973093215954679034375821888D-08
    weight(13) = 0.181691036255544979555476861323D-10
    weight(14) = 0.109540138892868740297645078918D-12
    weight(15) = 0.261737388222337042155132062413D-15
    weight(16) = 0.167293569314615469085022374652D-18
    weight(17) = 0.106562631627404278815253271162D-22

  else if ( norder == 18 ) then

    xtab(1) =    0.781691666697054712986747615334D-01
    xtab(2) =    0.412490085259129291039101536536D+00
    xtab(3) =    0.101652017962353968919093686187D+01
    xtab(4) =    0.189488850996976091426727831954D+01
    xtab(5) =    0.305435311320265975115241130719D+01
    xtab(6) =    0.450420553888989282633795571455D+01
    xtab(7) =    0.625672507394911145274209116326D+01
    xtab(8) =    0.832782515660563002170470261564D+01
    xtab(9) =    0.107379900477576093352179033397D+02
    xtab(10) =   0.135136562075550898190863812108D+02
    xtab(11) =   0.166893062819301059378183984163D+02
    xtab(12) =   0.203107676262677428561313764553D+02
    xtab(13) =   0.244406813592837027656442257980D+02
    xtab(14) =   0.291682086625796161312980677805D+02
    xtab(15) =   0.346279270656601721454012429438D+02
    xtab(16) =   0.410418167728087581392948614284D+02
    xtab(17) =   0.488339227160865227486586093290D+02
    xtab(18) =   0.590905464359012507037157810181D+02

    weight(1) =  0.185588603146918805623337752284D+00
    weight(2) =  0.310181766370225293649597595713D+00
    weight(3) =  0.267866567148536354820854394783D+00
    weight(4) =  0.152979747468074906553843082053D+00
    weight(5) =  0.614349178609616527076780103487D-01
    weight(6) =  0.176872130807729312772600233761D-01
    weight(7) =  0.366017976775991779802657207890D-02
    weight(8) =  0.540622787007735323128416319257D-03
    weight(9) =  0.561696505121423113817929049294D-04
    weight(10) = 0.401530788370115755858883625279D-05
    weight(11) = 0.191466985667567497969210011321D-06
    weight(12) = 0.583609526863159412918086289717D-08
    weight(13) = 0.107171126695539012772851317562D-09
    weight(14) = 0.108909871388883385562011298291D-11
    weight(15) = 0.538666474837830887608094323164D-14
    weight(16) = 0.104986597803570340877859934846D-16
    weight(17) = 0.540539845163105364356554467358D-20
    weight(18) = 0.269165326920102862708377715980D-24

  else if ( norder == 19 ) then

    xtab(1) =    0.741587837572050877131369916024D-01
    xtab(2) =    0.391268613319994607337648350299D+00
    xtab(3) =    0.963957343997958058624879377130D+00
    xtab(4) =    0.179617558206832812557725825252D+01
    xtab(5) =    0.289365138187378399116494713237D+01
    xtab(6) =    0.426421553962776647436040018167D+01
    xtab(7) =    0.591814156164404855815360191408D+01
    xtab(8) =    0.786861891533473373105668358176D+01
    xtab(9) =    0.101324237168152659251627415800D+02
    xtab(10) =   0.127308814638423980045092979656D+02
    xtab(11) =   0.156912783398358885454136069861D+02
    xtab(12) =   0.190489932098235501532136429732D+02
    xtab(13) =   0.228508497608294829323930586693D+02
    xtab(14) =   0.271606693274114488789963947149D+02
    xtab(15) =   0.320691222518622423224362865906D+02
    xtab(16) =   0.377129058012196494770647508283D+02
    xtab(17) =   0.443173627958314961196067736013D+02
    xtab(18) =   0.523129024574043831658644222420D+02
    xtab(19) =   0.628024231535003758413504690673D+02

    weight(1) =  0.176768474915912502251035479815D+00
    weight(2) =  0.300478143607254379482156807712D+00
    weight(3) =  0.267599547038175030772695440648D+00
    weight(4) =  0.159913372135580216785512147895D+00
    weight(5) =  0.682493799761491134552355368344D-01
    weight(6) =  0.212393076065443249244062193091D-01
    weight(7) =  0.484162735114839596725013121019D-02
    weight(8) =  0.804912747381366766594647138204D-03
    weight(9) =  0.965247209315350170843161738801D-04
    weight(10) = 0.820730525805103054408982992869D-05
    weight(11) = 0.483056672473077253944806671560D-06
    weight(12) = 0.190499136112328569993615674552D-07
    weight(13) = 0.481668463092806155766936380273D-09
    weight(14) = 0.734825883955114437684376840171D-11
    weight(15) = 0.620227538757261639893719012423D-13
    weight(16) = 0.254143084301542272371866857954D-15
    weight(17) = 0.407886129682571235007187465134D-18
    weight(18) = 0.170775018759383706100412325084D-21
    weight(19) = 0.671506464990818995998969111749D-26

  else if ( norder == 20 ) then

    xtab(1) =    0.705398896919887533666890045842D-01
    xtab(2) =    0.372126818001611443794241388761D+00
    xtab(3) =    0.916582102483273564667716277074D+00
    xtab(4) =    0.170730653102834388068768966741D+01
    xtab(5) =    0.274919925530943212964503046049D+01
    xtab(6) =    0.404892531385088692237495336913D+01
    xtab(7) =    0.561517497086161651410453988565D+01
    xtab(8) =    0.745901745367106330976886021837D+01
    xtab(9) =    0.959439286958109677247367273428D+01
    xtab(10) =   0.120388025469643163096234092989D+02
    xtab(11) =   0.148142934426307399785126797100D+02
    xtab(12) =   0.179488955205193760173657909926D+02
    xtab(13) =   0.214787882402850109757351703696D+02
    xtab(14) =   0.254517027931869055035186774846D+02
    xtab(15) =   0.299325546317006120067136561352D+02
    xtab(16) =   0.350134342404790000062849359067D+02
    xtab(17) =   0.408330570567285710620295677078D+02
    xtab(18) =   0.476199940473465021399416271529D+02
    xtab(19) =   0.558107957500638988907507734445D+02
    xtab(20) =   0.665244165256157538186403187915D+02

    weight(1) =  0.168746801851113862149223899689D+00
    weight(2) =  0.291254362006068281716795323812D+00
    weight(3) =  0.266686102867001288549520868998D+00
    weight(4) =  0.166002453269506840031469127816D+00
    weight(5) =  0.748260646687923705400624639615D-01
    weight(6) =  0.249644173092832210728227383234D-01
    weight(7) =  0.620255084457223684744754785395D-02
    weight(8) =  0.114496238647690824203955356969D-02
    weight(9) =  0.155741773027811974779809513214D-03
    weight(10) = 0.154014408652249156893806714048D-04
    weight(11) = 0.108648636651798235147970004439D-05
    weight(12) = 0.533012090955671475092780244305D-07
    weight(13) = 0.175798117905058200357787637840D-08
    weight(14) = 0.372550240251232087262924585338D-10
    weight(15) = 0.476752925157819052449488071613D-12
    weight(16) = 0.337284424336243841236506064991D-14
    weight(17) = 0.115501433950039883096396247181D-16
    weight(18) = 0.153952214058234355346383319667D-19
    weight(19) = 0.528644272556915782880273587683D-23
    weight(20) = 0.165645661249902329590781908529D-27

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 to 20.'
    stop

  end if

  return
end
subroutine laguerre_sum ( func, a, norder, xtab, weight, result )
!
!*******************************************************************************
!
!! LAGUERRE_SUM carries out Laguerre quadrature over [ A, +Infinity ).
!
!
!  Discussion:
!
!    The simplest Laguerre integral to approximate is the
!    integral from 0 to INFINITY of EXP(-X) * F(X).  When this is so,
!    it is easy to modify the rule to approximate the integral from
!    A to INFINITY as well.
!
!    Another common Laguerre integral to approximate is the
!    integral from 0 to Infinity of EXP(-X) * X**ALPHA * F(X).
!    This routine may be used to sum up the terms of the Laguerre
!    rule for such an integral as well.  However, if ALPHA is nonzero,
!    then there is no simple way to extend the rule to approximate the
!    integral from A to INFINITY.  The simplest procedures would be
!    to approximate the integral from 0 to A.
!
!  Integration interval:
!
!    [ A, +Infinity ) or [ 0, +Infinity )
!
!  Weight function:
!
!    EXP ( - X ) or EXP ( - X ) * X**ALPHA
!
!  Integral to approximate:
!
!    Integral ( A <= X <= +Infinity ) EXP ( -X ) * F(X) dX or
!    Integral ( 0 <= X <= +Infinity ) EXP ( -X ) * X**ALPHA * F(X) dX
!
!  Approximate integral:
!
!    EXP ( - A ) * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) + A ) 
!
!    or
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, the beginning of the integration interval.
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ), external :: func
  integer i
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_SUM - Fatal error!'
    write ( *, '(a,i6)' ) '  Nonpositive NORDER = ', norder
    stop
  end if

  result = 0.0D+00
  do i = 1, norder
    result = result + weight(i) * func ( xtab(i) + a )
  end do
  result = exp ( - a ) * result

  return
end
subroutine legendre_com ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    16 September 1998
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) d_pi
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer i
  integer iback
  integer k
  integer m
  integer mp1mi
  integer ncopy
  integer nmove
  real ( kind = 8 ) p
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) weight(norder)
!
  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COM - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop
  end if

  e1 = real ( norder * ( norder + 1 ), kind = 8 )

  m = ( norder + 1 ) / 2

  do i = 1, ( norder + 1 ) / 2

    mp1mi = m + 1 - i
    t = real ( 4 * i - 1, kind = 8 ) * d_pi ( ) &
      / real ( 4 * norder + 2, kind = 8 )
    x0 = cos(t) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( norder, kind = 8 ) ) &
      / real ( 8 * norder**2, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, norder
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( norder, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0**2 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0**2 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0**2 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0**2 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v**2 - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    xtab(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp**2 ) / fx**2

  end do

  if ( mod ( norder, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( norder + 1 ) / 2
  ncopy = norder - nmove

  do i = 1, nmove
    iback = norder + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, norder - nmove
    xtab(i) = - xtab(norder+1-i)
    weight(i) = weight(norder+1-i)
  end do

  return
end
subroutine legendre_recur ( p2, dp2, p1, x, norder )
!
!*******************************************************************************
!
!! LEGENDRE_RECUR finds the value and derivative of a Legendre polynomial.
!
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of P(NORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of P'(NORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of P(NORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
!
  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, norder

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( real ( 2 * i - 1, kind = 8 ) * x * p1 &
         + real (   - i + 1, kind = 8 )     * p0 ) &
         / real (     i,     kind = 8 )

    dp2 = ( real ( 2 * i - 1, kind = 8 ) * ( p1 + x * dp1 ) &
          - real (     i - 1, kind = 8 ) * dp0 ) &
          / real (     i,     kind = 8 )

  end do

  return
end
subroutine legendre_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-1).
!
!  Note:
!
!    The abscissas of the rule are the zeroes of the Legendre polynomial
!    P(NORDER)(X).
!
!    The integral produced by a Gauss-Legendre rule is equal to the
!    integral of the unique polynomial of degree NORDER-1 which
!    agrees with the function at the NORDER abscissas of the rule.
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20, 32 or 64.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =   0.0D+00

    weight(1) = 2.0D+00

  else if ( norder == 2 ) then

    xtab(1) = - 0.577350269189625764509148780502D+00
    xtab(2) =   0.577350269189625764509148780502D+00

    weight(1) = 1.0D+00
    weight(2) = 1.0D+00

  else if ( norder == 3 ) then

    xtab(1) = - 0.774596669241483377035853079956D+00
    xtab(2) =   0.0D+00
    xtab(3) =   0.774596669241483377035853079956D+00

    weight(1) = 5.0D+00 / 9.0D+00
    weight(2) = 8.0D+00 / 9.0D+00
    weight(3) = 5.0D+00 / 9.0D+00

  else if ( norder == 4 ) then

    xtab(1) = - 0.861136311594052575223946488893D+00
    xtab(2) = - 0.339981043584856264802665759103D+00
    xtab(3) =   0.339981043584856264802665759103D+00
    xtab(4) =   0.861136311594052575223946488893D+00

    weight(1) = 0.347854845137453857373063949222D+00
    weight(2) = 0.652145154862546142626936050778D+00
    weight(3) = 0.652145154862546142626936050778D+00
    weight(4) = 0.347854845137453857373063949222D+00

  else if ( norder == 5 ) then

    xtab(1) = - 0.906179845938663992797626878299D+00
    xtab(2) = - 0.538469310105683091036314420700D+00
    xtab(3) =   0.0D+00
    xtab(4) =   0.538469310105683091036314420700D+00
    xtab(5) =   0.906179845938663992797626878299D+00

    weight(1) = 0.236926885056189087514264040720D+00
    weight(2) = 0.478628670499366468041291514836D+00
    weight(3) = 0.568888888888888888888888888889D+00
    weight(4) = 0.478628670499366468041291514836D+00
    weight(5) = 0.236926885056189087514264040720D+00

  else if ( norder == 6 ) then

    xtab(1) = - 0.932469514203152027812301554494D+00
    xtab(2) = - 0.661209386466264513661399595020D+00
    xtab(3) = - 0.238619186083196908630501721681D+00
    xtab(4) =   0.238619186083196908630501721681D+00
    xtab(5) =   0.661209386466264513661399595020D+00
    xtab(6) =   0.932469514203152027812301554494D+00

    weight(1) = 0.171324492379170345040296142173D+00
    weight(2) = 0.360761573048138607569833513838D+00
    weight(3) = 0.467913934572691047389870343990D+00
    weight(4) = 0.467913934572691047389870343990D+00
    weight(5) = 0.360761573048138607569833513838D+00
    weight(6) = 0.171324492379170345040296142173D+00

  else if ( norder == 7 ) then

    xtab(1) = - 0.949107912342758524526189684048D+00
    xtab(2) = - 0.741531185599394439863864773281D+00
    xtab(3) = - 0.405845151377397166906606412077D+00
    xtab(4) =   0.0D+00
    xtab(5) =   0.405845151377397166906606412077D+00
    xtab(6) =   0.741531185599394439863864773281D+00
    xtab(7) =   0.949107912342758524526189684048D+00

    weight(1) = 0.129484966168869693270611432679D+00
    weight(2) = 0.279705391489276667901467771424D+00
    weight(3) = 0.381830050505118944950369775489D+00
    weight(4) = 0.417959183673469387755102040816D+00
    weight(5) = 0.381830050505118944950369775489D+00
    weight(6) = 0.279705391489276667901467771424D+00
    weight(7) = 0.129484966168869693270611432679D+00

  else if ( norder == 8 ) then

    xtab(1) = - 0.960289856497536231683560868569D+00
    xtab(2) = - 0.796666477413626739591553936476D+00
    xtab(3) = - 0.525532409916328985817739049189D+00
    xtab(4) = - 0.183434642495649804939476142360D+00
    xtab(5) =   0.183434642495649804939476142360D+00
    xtab(6) =   0.525532409916328985817739049189D+00
    xtab(7) =   0.796666477413626739591553936476D+00
    xtab(8) =   0.960289856497536231683560868569D+00

    weight(1) = 0.101228536290376259152531354310D+00
    weight(2) = 0.222381034453374470544355994426D+00
    weight(3) = 0.313706645877887287337962201987D+00
    weight(4) = 0.362683783378361982965150449277D+00
    weight(5) = 0.362683783378361982965150449277D+00
    weight(6) = 0.313706645877887287337962201987D+00
    weight(7) = 0.222381034453374470544355994426D+00
    weight(8) = 0.101228536290376259152531354310D+00

  else if ( norder == 9 ) then

    xtab(1) = - 0.968160239507626089835576202904D+00
    xtab(2) = - 0.836031107326635794299429788070D+00
    xtab(3) = - 0.613371432700590397308702039341D+00
    xtab(4) = - 0.324253423403808929038538014643D+00
    xtab(5) =   0.0D+00
    xtab(6) =   0.324253423403808929038538014643D+00
    xtab(7) =   0.613371432700590397308702039341D+00
    xtab(8) =   0.836031107326635794299429788070D+00
    xtab(9) =   0.968160239507626089835576202904D+00

    weight(1) = 0.812743883615744119718921581105D-01
    weight(2) = 0.180648160694857404058472031243D+00
    weight(3) = 0.260610696402935462318742869419D+00
    weight(4) = 0.312347077040002840068630406584D+00
    weight(5) = 0.330239355001259763164525069287D+00
    weight(6) = 0.312347077040002840068630406584D+00
    weight(7) = 0.260610696402935462318742869419D+00
    weight(8) = 0.180648160694857404058472031243D+00
    weight(9) = 0.812743883615744119718921581105D-01

  else if ( norder == 10 ) then

    xtab(1) =  - 0.973906528517171720077964012084D+00
    xtab(2) =  - 0.865063366688984510732096688423D+00
    xtab(3) =  - 0.679409568299024406234327365115D+00
    xtab(4) =  - 0.433395394129247190799265943166D+00
    xtab(5) =  - 0.148874338981631210884826001130D+00
    xtab(6) =    0.148874338981631210884826001130D+00
    xtab(7) =    0.433395394129247190799265943166D+00
    xtab(8) =    0.679409568299024406234327365115D+00
    xtab(9) =    0.865063366688984510732096688423D+00
    xtab(10) =   0.973906528517171720077964012084D+00

    weight(1) =  0.666713443086881375935688098933D-01
    weight(2) =  0.149451349150580593145776339658D+00
    weight(3) =  0.219086362515982043995534934228D+00
    weight(4) =  0.269266719309996355091226921569D+00
    weight(5) =  0.295524224714752870173892994651D+00
    weight(6) =  0.295524224714752870173892994651D+00
    weight(7) =  0.269266719309996355091226921569D+00
    weight(8) =  0.219086362515982043995534934228D+00
    weight(9) =  0.149451349150580593145776339658D+00
    weight(10) = 0.666713443086881375935688098933D-01

  else if ( norder == 11 ) then

    xtab(1) =  - 0.978228658146056992803938001123D+00
    xtab(2) =  - 0.887062599768095299075157769304D+00
    xtab(3) =  - 0.730152005574049324093416252031D+00
    xtab(4) =  - 0.519096129206811815925725669459D+00
    xtab(5) =  - 0.269543155952344972331531985401D+00
    xtab(6) =    0.0D+00
    xtab(7) =    0.269543155952344972331531985401D+00
    xtab(8) =    0.519096129206811815925725669459D+00
    xtab(9) =    0.730152005574049324093416252031D+00
    xtab(10) =   0.887062599768095299075157769304D+00
    xtab(11) =   0.978228658146056992803938001123D+00

    weight(1) =  0.556685671161736664827537204425D-01
    weight(2) =  0.125580369464904624634694299224D+00
    weight(3) =  0.186290210927734251426097641432D+00
    weight(4) =  0.233193764591990479918523704843D+00
    weight(5) =  0.262804544510246662180688869891D+00
    weight(6) =  0.272925086777900630714483528336D+00
    weight(7) =  0.262804544510246662180688869891D+00
    weight(8) =  0.233193764591990479918523704843D+00
    weight(9) =  0.186290210927734251426097641432D+00
    weight(10) = 0.125580369464904624634694299224D+00
    weight(11) = 0.556685671161736664827537204425D-01

  else if ( norder == 12 ) then

    xtab(1) =  - 0.981560634246719250690549090149D+00
    xtab(2) =  - 0.904117256370474856678465866119D+00
    xtab(3) =  - 0.769902674194304687036893833213D+00
    xtab(4) =  - 0.587317954286617447296702418941D+00
    xtab(5) =  - 0.367831498998180193752691536644D+00
    xtab(6) =  - 0.125233408511468915472441369464D+00
    xtab(7) =    0.125233408511468915472441369464D+00
    xtab(8) =    0.367831498998180193752691536644D+00
    xtab(9) =    0.587317954286617447296702418941D+00
    xtab(10) =   0.769902674194304687036893833213D+00
    xtab(11) =   0.904117256370474856678465866119D+00
    xtab(12) =   0.981560634246719250690549090149D+00

    weight(1) =  0.471753363865118271946159614850D-01
    weight(2) =  0.106939325995318430960254718194D+00
    weight(3) =  0.160078328543346226334652529543D+00
    weight(4) =  0.203167426723065921749064455810D+00
    weight(5) =  0.233492536538354808760849898925D+00
    weight(6) =  0.249147045813402785000562436043D+00
    weight(7) =  0.249147045813402785000562436043D+00
    weight(8) =  0.233492536538354808760849898925D+00
    weight(9) =  0.203167426723065921749064455810D+00
    weight(10) = 0.160078328543346226334652529543D+00
    weight(11) = 0.106939325995318430960254718194D+00
    weight(12) = 0.471753363865118271946159614850D-01

  else if ( norder == 13 ) then

    xtab(1) =  - 0.984183054718588149472829448807D+00
    xtab(2) =  - 0.917598399222977965206547836501D+00
    xtab(3) =  - 0.801578090733309912794206489583D+00
    xtab(4) =  - 0.642349339440340220643984606996D+00
    xtab(5) =  - 0.448492751036446852877912852128D+00
    xtab(6) =  - 0.230458315955134794065528121098D+00
    xtab(7) =    0.0D+00
    xtab(8) =    0.230458315955134794065528121098D+00
    xtab(9) =    0.448492751036446852877912852128D+00
    xtab(10) =   0.642349339440340220643984606996D+00
    xtab(11) =   0.801578090733309912794206489583D+00
    xtab(12) =   0.917598399222977965206547836501D+00
    xtab(13) =   0.984183054718588149472829448807D+00

    weight(1) =  0.404840047653158795200215922010D-01
    weight(2) =  0.921214998377284479144217759538D-01
    weight(3) =  0.138873510219787238463601776869D+00
    weight(4) =  0.178145980761945738280046691996D+00
    weight(5) =  0.207816047536888502312523219306D+00
    weight(6) =  0.226283180262897238412090186040D+00
    weight(7) =  0.232551553230873910194589515269D+00
    weight(8) =  0.226283180262897238412090186040D+00
    weight(9) =  0.207816047536888502312523219306D+00
    weight(10) = 0.178145980761945738280046691996D+00
    weight(11) = 0.138873510219787238463601776869D+00
    weight(12) = 0.921214998377284479144217759538D-01
    weight(13) = 0.404840047653158795200215922010D-01

  else if ( norder == 14 ) then

    xtab(1) =  - 0.986283808696812338841597266704D+00
    xtab(2) =  - 0.928434883663573517336391139378D+00
    xtab(3) =  - 0.827201315069764993189794742650D+00
    xtab(4) =  - 0.687292904811685470148019803019D+00
    xtab(5) =  - 0.515248636358154091965290718551D+00
    xtab(6) =  - 0.319112368927889760435671824168D+00
    xtab(7) =  - 0.108054948707343662066244650220D+00
    xtab(8) =    0.108054948707343662066244650220D+00
    xtab(9) =    0.319112368927889760435671824168D+00
    xtab(10) =   0.515248636358154091965290718551D+00
    xtab(11) =   0.687292904811685470148019803019D+00
    xtab(12) =   0.827201315069764993189794742650D+00
    xtab(13) =   0.928434883663573517336391139378D+00
    xtab(14) =   0.986283808696812338841597266704D+00

    weight(1) =  0.351194603317518630318328761382D-01
    weight(2) =  0.801580871597602098056332770629D-01
    weight(3) =  0.121518570687903184689414809072D+00
    weight(4) =  0.157203167158193534569601938624D+00
    weight(5) =  0.185538397477937813741716590125D+00
    weight(6) =  0.205198463721295603965924065661D+00
    weight(7) =  0.215263853463157790195876443316D+00
    weight(8) =  0.215263853463157790195876443316D+00
    weight(9) =  0.205198463721295603965924065661D+00
    weight(10) = 0.185538397477937813741716590125D+00
    weight(11) = 0.157203167158193534569601938624D+00
    weight(12) = 0.121518570687903184689414809072D+00
    weight(13) = 0.801580871597602098056332770629D-01
    weight(14) = 0.351194603317518630318328761382D-01

  else if ( norder == 15 ) then

    xtab(1) =  - 0.987992518020485428489565718587D+00
    xtab(2) =  - 0.937273392400705904307758947710D+00
    xtab(3) =  - 0.848206583410427216200648320774D+00
    xtab(4) =  - 0.724417731360170047416186054614D+00
    xtab(5) =  - 0.570972172608538847537226737254D+00
    xtab(6) =  - 0.394151347077563369897207370981D+00
    xtab(7) =  - 0.201194093997434522300628303395D+00
    xtab(8) =    0.0D+00
    xtab(9) =    0.201194093997434522300628303395D+00
    xtab(10) =   0.394151347077563369897207370981D+00
    xtab(11) =   0.570972172608538847537226737254D+00
    xtab(12) =   0.724417731360170047416186054614D+00
    xtab(13) =   0.848206583410427216200648320774D+00
    xtab(14) =   0.937273392400705904307758947710D+00
    xtab(15) =   0.987992518020485428489565718587D+00

    weight(1) =  0.307532419961172683546283935772D-01
    weight(2) =  0.703660474881081247092674164507D-01
    weight(3) =  0.107159220467171935011869546686D+00
    weight(4) =  0.139570677926154314447804794511D+00
    weight(5) =  0.166269205816993933553200860481D+00
    weight(6) =  0.186161000015562211026800561866D+00
    weight(7) =  0.198431485327111576456118326444D+00
    weight(8) =  0.202578241925561272880620199968D+00
    weight(9) =  0.198431485327111576456118326444D+00
    weight(10) = 0.186161000015562211026800561866D+00
    weight(11) = 0.166269205816993933553200860481D+00
    weight(12) = 0.139570677926154314447804794511D+00
    weight(13) = 0.107159220467171935011869546686D+00
    weight(14) = 0.703660474881081247092674164507D-01
    weight(15) = 0.307532419961172683546283935772D-01

  else if ( norder == 16 ) then

    xtab(1) =  - 0.989400934991649932596154173450D+00
    xtab(2) =  - 0.944575023073232576077988415535D+00
    xtab(3) =  - 0.865631202387831743880467897712D+00
    xtab(4) =  - 0.755404408355003033895101194847D+00
    xtab(5) =  - 0.617876244402643748446671764049D+00
    xtab(6) =  - 0.458016777657227386342419442984D+00
    xtab(7) =  - 0.281603550779258913230460501460D+00
    xtab(8) =  - 0.950125098376374401853193354250D-01
    xtab(9) =    0.950125098376374401853193354250D-01
    xtab(10) =   0.281603550779258913230460501460D+00
    xtab(11) =   0.458016777657227386342419442984D+00
    xtab(12) =   0.617876244402643748446671764049D+00
    xtab(13) =   0.755404408355003033895101194847D+00
    xtab(14) =   0.865631202387831743880467897712D+00
    xtab(15) =   0.944575023073232576077988415535D+00
    xtab(16) =   0.989400934991649932596154173450D+00

    weight(1) =  0.271524594117540948517805724560D-01
    weight(2) =  0.622535239386478928628438369944D-01
    weight(3) =  0.951585116824927848099251076022D-01
    weight(4) =  0.124628971255533872052476282192D+00
    weight(5) =  0.149595988816576732081501730547D+00
    weight(6) =  0.169156519395002538189312079030D+00
    weight(7) =  0.182603415044923588866763667969D+00
    weight(8) =  0.189450610455068496285396723208D+00
    weight(9) =  0.189450610455068496285396723208D+00
    weight(10) = 0.182603415044923588866763667969D+00
    weight(11) = 0.169156519395002538189312079030D+00
    weight(12) = 0.149595988816576732081501730547D+00
    weight(13) = 0.124628971255533872052476282192D+00
    weight(14) = 0.951585116824927848099251076022D-01
    weight(15) = 0.622535239386478928628438369944D-01
    weight(16) = 0.271524594117540948517805724560D-01

  else if ( norder == 17 ) then

    xtab(1) =  - 0.990575475314417335675434019941D+00
    xtab(2) =  - 0.950675521768767761222716957896D+00
    xtab(3) =  - 0.880239153726985902122955694488D+00
    xtab(4) =  - 0.781514003896801406925230055520D+00
    xtab(5) =  - 0.657671159216690765850302216643D+00
    xtab(6) =  - 0.512690537086476967886246568630D+00
    xtab(7) =  - 0.351231763453876315297185517095D+00
    xtab(8) =  - 0.178484181495847855850677493654D+00
    xtab(9) =    0.0D+00
    xtab(10) =   0.178484181495847855850677493654D+00
    xtab(11) =   0.351231763453876315297185517095D+00
    xtab(12) =   0.512690537086476967886246568630D+00
    xtab(13) =   0.657671159216690765850302216643D+00
    xtab(14) =   0.781514003896801406925230055520D+00
    xtab(15) =   0.880239153726985902122955694488D+00
    xtab(16) =   0.950675521768767761222716957896D+00
    xtab(17) =   0.990575475314417335675434019941D+00

    weight(1) =  0.241483028685479319601100262876D-01
    weight(2) =  0.554595293739872011294401653582D-01
    weight(3) =  0.850361483171791808835353701911D-01
    weight(4) =  0.111883847193403971094788385626D+00
    weight(5) =  0.135136368468525473286319981702D+00
    weight(6) =  0.154045761076810288081431594802D+00
    weight(7) =  0.168004102156450044509970663788D+00
    weight(8) =  0.176562705366992646325270990113D+00
    weight(9) =  0.179446470356206525458265644262D+00
    weight(10) = 0.176562705366992646325270990113D+00
    weight(11) = 0.168004102156450044509970663788D+00
    weight(12) = 0.154045761076810288081431594802D+00
    weight(13) = 0.135136368468525473286319981702D+00
    weight(14) = 0.111883847193403971094788385626D+00
    weight(15) = 0.850361483171791808835353701911D-01
    weight(16) = 0.554595293739872011294401653582D-01
    weight(17) = 0.241483028685479319601100262876D-01

  else if ( norder == 18 ) then

    xtab(1) =  - 0.991565168420930946730016004706D+00
    xtab(2) =  - 0.955823949571397755181195892930D+00
    xtab(3) =  - 0.892602466497555739206060591127D+00
    xtab(4) =  - 0.803704958972523115682417455015D+00
    xtab(5) =  - 0.691687043060353207874891081289D+00
    xtab(6) =  - 0.559770831073947534607871548525D+00
    xtab(7) =  - 0.411751161462842646035931793833D+00
    xtab(8) =  - 0.251886225691505509588972854878D+00
    xtab(9) =  - 0.847750130417353012422618529358D-01
    xtab(10) =   0.847750130417353012422618529358D-01
    xtab(11) =   0.251886225691505509588972854878D+00
    xtab(12) =   0.411751161462842646035931793833D+00
    xtab(13) =   0.559770831073947534607871548525D+00
    xtab(14) =   0.691687043060353207874891081289D+00
    xtab(15) =   0.803704958972523115682417455015D+00
    xtab(16) =   0.892602466497555739206060591127D+00
    xtab(17) =   0.955823949571397755181195892930D+00
    xtab(18) =   0.991565168420930946730016004706D+00

    weight(1) =  0.216160135264833103133427102665D-01
    weight(2) =  0.497145488949697964533349462026D-01
    weight(3) =  0.764257302548890565291296776166D-01
    weight(4) =  0.100942044106287165562813984925D+00
    weight(5) =  0.122555206711478460184519126800D+00
    weight(6) =  0.140642914670650651204731303752D+00
    weight(7) =  0.154684675126265244925418003836D+00
    weight(8) =  0.164276483745832722986053776466D+00
    weight(9) =  0.169142382963143591840656470135D+00
    weight(10) = 0.169142382963143591840656470135D+00
    weight(11) = 0.164276483745832722986053776466D+00
    weight(12) = 0.154684675126265244925418003836D+00
    weight(13) = 0.140642914670650651204731303752D+00
    weight(14) = 0.122555206711478460184519126800D+00
    weight(15) = 0.100942044106287165562813984925D+00
    weight(16) = 0.764257302548890565291296776166D-01
    weight(17) = 0.497145488949697964533349462026D-01
    weight(18) = 0.216160135264833103133427102665D-01

  else if ( norder == 19 ) then

    xtab(1) =  - 0.992406843843584403189017670253D+00
    xtab(2) =  - 0.960208152134830030852778840688D+00
    xtab(3) =  - 0.903155903614817901642660928532D+00
    xtab(4) =  - 0.822714656537142824978922486713D+00
    xtab(5) =  - 0.720966177335229378617095860824D+00
    xtab(6) =  - 0.600545304661681023469638164946D+00
    xtab(7) =  - 0.464570741375960945717267148104D+00
    xtab(8) =  - 0.316564099963629831990117328850D+00
    xtab(9) =  - 0.160358645640225375868096115741D+00
    xtab(10) =   0.0D+00
    xtab(11) =   0.160358645640225375868096115741D+00
    xtab(12) =   0.316564099963629831990117328850D+00
    xtab(13) =   0.464570741375960945717267148104D+00
    xtab(14) =   0.600545304661681023469638164946D+00
    xtab(15) =   0.720966177335229378617095860824D+00
    xtab(16) =   0.822714656537142824978922486713D+00
    xtab(17) =   0.903155903614817901642660928532D+00
    xtab(18) =   0.960208152134830030852778840688D+00
    xtab(19) =   0.992406843843584403189017670253D+00

    weight(1) =  0.194617882297264770363120414644D-01
    weight(2) =  0.448142267656996003328381574020D-01
    weight(3) =  0.690445427376412265807082580060D-01
    weight(4) =  0.914900216224499994644620941238D-01
    weight(5) =  0.111566645547333994716023901682D+00
    weight(6) =  0.128753962539336227675515784857D+00
    weight(7) =  0.142606702173606611775746109442D+00
    weight(8) =  0.152766042065859666778855400898D+00
    weight(9) =  0.158968843393954347649956439465D+00
    weight(10) = 0.161054449848783695979163625321D+00
    weight(11) = 0.158968843393954347649956439465D+00
    weight(12) = 0.152766042065859666778855400898D+00
    weight(13) = 0.142606702173606611775746109442D+00
    weight(14) = 0.128753962539336227675515784857D+00
    weight(15) = 0.111566645547333994716023901682D+00
    weight(16) = 0.914900216224499994644620941238D-01
    weight(17) = 0.690445427376412265807082580060D-01
    weight(18) = 0.448142267656996003328381574020D-01
    weight(19) = 0.194617882297264770363120414644D-01

  else if ( norder == 20 ) then

    xtab(1) =  - 0.993128599185094924786122388471D+00
    xtab(2) =  - 0.963971927277913791267666131197D+00
    xtab(3) =  - 0.912234428251325905867752441203D+00
    xtab(4) =  - 0.839116971822218823394529061702D+00
    xtab(5) =  - 0.746331906460150792614305070356D+00
    xtab(6) =  - 0.636053680726515025452836696226D+00
    xtab(7) =  - 0.510867001950827098004364050955D+00
    xtab(8) =  - 0.373706088715419560672548177025D+00
    xtab(9) =  - 0.227785851141645078080496195369D+00
    xtab(10) = - 0.765265211334973337546404093988D-01
    xtab(11) =   0.765265211334973337546404093988D-01
    xtab(12) =   0.227785851141645078080496195369D+00
    xtab(13) =   0.373706088715419560672548177025D+00
    xtab(14) =   0.510867001950827098004364050955D+00
    xtab(15) =   0.636053680726515025452836696226D+00
    xtab(16) =   0.746331906460150792614305070356D+00
    xtab(17) =   0.839116971822218823394529061702D+00
    xtab(18) =   0.912234428251325905867752441203D+00
    xtab(19) =   0.963971927277913791267666131197D+00
    xtab(20) =   0.993128599185094924786122388471D+00

    weight(1) =  0.176140071391521183118619623519D-01
    weight(2) =  0.406014298003869413310399522749D-01
    weight(3) =  0.626720483341090635695065351870D-01
    weight(4) =  0.832767415767047487247581432220D-01
    weight(5) =  0.101930119817240435036750135480D+00
    weight(6) =  0.118194531961518417312377377711D+00
    weight(7) =  0.131688638449176626898494499748D+00
    weight(8) =  0.142096109318382051329298325067D+00
    weight(9) =  0.149172986472603746787828737002D+00
    weight(10) = 0.152753387130725850698084331955D+00
    weight(11) = 0.152753387130725850698084331955D+00
    weight(12) = 0.149172986472603746787828737002D+00
    weight(13) = 0.142096109318382051329298325067D+00
    weight(14) = 0.131688638449176626898494499748D+00
    weight(15) = 0.118194531961518417312377377711D+00
    weight(16) = 0.101930119817240435036750135480D+00
    weight(17) = 0.832767415767047487247581432220D-01
    weight(18) = 0.626720483341090635695065351870D-01
    weight(19) = 0.406014298003869413310399522749D-01
    weight(20) = 0.176140071391521183118619623519D-01

  else if ( norder == 32 ) then

    xtab(1) =  - 0.997263861849481563544981128665D+00
    xtab(2) =  - 0.985611511545268335400175044631D+00
    xtab(3) =  - 0.964762255587506430773811928118D+00
    xtab(4) =  - 0.934906075937739689170919134835D+00
    xtab(5) =  - 0.896321155766052123965307243719D+00
    xtab(6) =  - 0.849367613732569970133693004968D+00
    xtab(7) =  - 0.794483795967942406963097298970D+00
    xtab(8) =  - 0.732182118740289680387426665091D+00
    xtab(9) =  - 0.663044266930215200975115168663D+00
    xtab(10) = - 0.587715757240762329040745476402D+00
    xtab(11) = - 0.506899908932229390023747474378D+00
    xtab(12) = - 0.421351276130635345364119436172D+00
    xtab(13) = - 0.331868602282127649779916805730D+00
    xtab(14) = - 0.239287362252137074544603209166D+00
    xtab(15) = - 0.144471961582796493485186373599D+00
    xtab(16) = - 0.483076656877383162348125704405D-01
    xtab(17) =   0.483076656877383162348125704405D-01
    xtab(18) =   0.144471961582796493485186373599D+00
    xtab(19) =   0.239287362252137074544603209166D+00
    xtab(20) =   0.331868602282127649779916805730D+00
    xtab(21) =   0.421351276130635345364119436172D+00
    xtab(22) =   0.506899908932229390023747474378D+00
    xtab(23) =   0.587715757240762329040745476402D+00
    xtab(24) =   0.663044266930215200975115168663D+00
    xtab(25) =   0.732182118740289680387426665091D+00
    xtab(26) =   0.794483795967942406963097298970D+00
    xtab(27) =   0.849367613732569970133693004968D+00
    xtab(28) =   0.896321155766052123965307243719D+00
    xtab(29) =   0.934906075937739689170919134835D+00
    xtab(30) =   0.964762255587506430773811928118D+00
    xtab(31) =   0.985611511545268335400175044631D+00
    xtab(32) =   0.997263861849481563544981128665D+00

    weight(1) =  0.701861000947009660040706373885D-02
    weight(2) =  0.162743947309056706051705622064D-01
    weight(3) =  0.253920653092620594557525897892D-01
    weight(4) =  0.342738629130214331026877322524D-01
    weight(5) =  0.428358980222266806568786466061D-01
    weight(6) =  0.509980592623761761961632446895D-01
    weight(7) =  0.586840934785355471452836373002D-01
    weight(8) =  0.658222227763618468376500637069D-01
    weight(9) =  0.723457941088485062253993564785D-01
    weight(10) = 0.781938957870703064717409188283D-01
    weight(11) = 0.833119242269467552221990746043D-01
    weight(12) = 0.876520930044038111427714627518D-01
    weight(13) = 0.911738786957638847128685771116D-01
    weight(14) = 0.938443990808045656391802376681D-01
    weight(15) = 0.956387200792748594190820022041D-01
    weight(16) = 0.965400885147278005667648300636D-01
    weight(17) = 0.965400885147278005667648300636D-01
    weight(18) = 0.956387200792748594190820022041D-01
    weight(19) = 0.938443990808045656391802376681D-01
    weight(20) = 0.911738786957638847128685771116D-01
    weight(21) = 0.876520930044038111427714627518D-01
    weight(22) = 0.833119242269467552221990746043D-01
    weight(23) = 0.781938957870703064717409188283D-01
    weight(24) = 0.723457941088485062253993564785D-01
    weight(25) = 0.658222227763618468376500637069D-01
    weight(26) = 0.586840934785355471452836373002D-01
    weight(27) = 0.509980592623761761961632446895D-01
    weight(28) = 0.428358980222266806568786466061D-01
    weight(29) = 0.342738629130214331026877322524D-01
    weight(30) = 0.253920653092620594557525897892D-01
    weight(31) = 0.162743947309056706051705622064D-01
    weight(32) = 0.701861000947009660040706373885D-02

  else if ( norder == 64 ) then

    xtab(1) =  - 0.999305041735772139456905624346D+00
    xtab(2) =  - 0.996340116771955279346924500676D+00
    xtab(3) =  - 0.991013371476744320739382383443D+00
    xtab(4) =  - 0.983336253884625956931299302157D+00
    xtab(5) =  - 0.973326827789910963741853507352D+00
    xtab(6) =  - 0.961008799652053718918614121897D+00
    xtab(7) =  - 0.946411374858402816062481491347D+00
    xtab(8) =  - 0.929569172131939575821490154559D+00
    xtab(9) =  - 0.910522137078502805756380668008D+00
    xtab(10) = - 0.889315445995114105853404038273D+00
    xtab(11) = - 0.865999398154092819760783385070D+00
    xtab(12) = - 0.840629296252580362751691544696D+00
    xtab(13) = - 0.813265315122797559741923338086D+00
    xtab(14) = - 0.783972358943341407610220525214D+00
    xtab(15) = - 0.752819907260531896611863774886D+00
    xtab(16) = - 0.719881850171610826848940217832D+00
    xtab(17) = - 0.685236313054233242563558371031D+00
    xtab(18) = - 0.648965471254657339857761231993D+00
    xtab(19) = - 0.611155355172393250248852971019D+00
    xtab(20) = - 0.571895646202634034283878116659D+00
    xtab(21) = - 0.531279464019894545658013903544D+00
    xtab(22) = - 0.489403145707052957478526307022D+00
    xtab(23) = - 0.446366017253464087984947714759D+00
    xtab(24) = - 0.402270157963991603695766771260D+00
    xtab(25) = - 0.357220158337668115950442615046D+00
    xtab(26) = - 0.311322871990210956157512698560D+00
    xtab(27) = - 0.264687162208767416373964172510D+00
    xtab(28) = - 0.217423643740007084149648748989D+00
    xtab(29) = - 0.169644420423992818037313629748D+00
    xtab(30) = - 0.121462819296120554470376463492D+00
    xtab(31) = - 0.729931217877990394495429419403D-01
    xtab(32) = - 0.243502926634244325089558428537D-01
    xtab(33) =   0.243502926634244325089558428537D-01
    xtab(34) =   0.729931217877990394495429419403D-01
    xtab(35) =   0.121462819296120554470376463492D+00
    xtab(36) =   0.169644420423992818037313629748D+00
    xtab(37) =   0.217423643740007084149648748989D+00
    xtab(38) =   0.264687162208767416373964172510D+00
    xtab(39) =   0.311322871990210956157512698560D+00
    xtab(40) =   0.357220158337668115950442615046D+00
    xtab(41) =   0.402270157963991603695766771260D+00
    xtab(42) =   0.446366017253464087984947714759D+00
    xtab(43) =   0.489403145707052957478526307022D+00
    xtab(44) =   0.531279464019894545658013903544D+00
    xtab(45) =   0.571895646202634034283878116659D+00
    xtab(46) =   0.611155355172393250248852971019D+00
    xtab(47) =   0.648965471254657339857761231993D+00
    xtab(48) =   0.685236313054233242563558371031D+00
    xtab(49) =   0.719881850171610826848940217832D+00
    xtab(50) =   0.752819907260531896611863774886D+00
    xtab(51) =   0.783972358943341407610220525214D+00
    xtab(52) =   0.813265315122797559741923338086D+00
    xtab(53) =   0.840629296252580362751691544696D+00
    xtab(54) =   0.865999398154092819760783385070D+00
    xtab(55) =   0.889315445995114105853404038273D+00
    xtab(56) =   0.910522137078502805756380668008D+00
    xtab(57) =   0.929569172131939575821490154559D+00
    xtab(58) =   0.946411374858402816062481491347D+00
    xtab(59) =   0.961008799652053718918614121897D+00
    xtab(60) =   0.973326827789910963741853507352D+00
    xtab(61) =   0.983336253884625956931299302157D+00
    xtab(62) =   0.991013371476744320739382383443D+00
    xtab(63) =   0.996340116771955279346924500676D+00
    xtab(64) =   0.999305041735772139456905624346D+00

    weight(1) =  0.178328072169643294729607914497D-02
    weight(2) =  0.414703326056246763528753572855D-02
    weight(3) =  0.650445796897836285611736039998D-02
    weight(4) =  0.884675982636394772303091465973D-02
    weight(5) =  0.111681394601311288185904930192D-01
    weight(6) =  0.134630478967186425980607666860D-01
    weight(7) =  0.157260304760247193219659952975D-01
    weight(8) =  0.179517157756973430850453020011D-01
    weight(9) =  0.201348231535302093723403167285D-01
    weight(10) = 0.222701738083832541592983303842D-01
    weight(11) = 0.243527025687108733381775504091D-01
    weight(12) = 0.263774697150546586716917926252D-01
    weight(13) = 0.283396726142594832275113052002D-01
    weight(14) = 0.302346570724024788679740598195D-01
    weight(15) = 0.320579283548515535854675043479D-01
    weight(16) = 0.338051618371416093915654821107D-01
    weight(17) = 0.354722132568823838106931467152D-01
    weight(18) = 0.370551285402400460404151018096D-01
    weight(19) = 0.385501531786156291289624969468D-01
    weight(20) = 0.399537411327203413866569261283D-01
    weight(21) = 0.412625632426235286101562974736D-01
    weight(22) = 0.424735151236535890073397679088D-01
    weight(23) = 0.435837245293234533768278609737D-01
    weight(24) = 0.445905581637565630601347100309D-01
    weight(25) = 0.454916279274181444797709969713D-01
    weight(26) = 0.462847965813144172959532492323D-01
    weight(27) = 0.469681828162100173253262857546D-01
    weight(28) = 0.475401657148303086622822069442D-01
    weight(29) = 0.479993885964583077281261798713D-01
    weight(30) = 0.483447622348029571697695271580D-01
    weight(31) = 0.485754674415034269347990667840D-01
    weight(32) = 0.486909570091397203833653907347D-01
    weight(33) = 0.486909570091397203833653907347D-01
    weight(34) = 0.485754674415034269347990667840D-01
    weight(35) = 0.483447622348029571697695271580D-01
    weight(36) = 0.479993885964583077281261798713D-01
    weight(37) = 0.475401657148303086622822069442D-01
    weight(38) = 0.469681828162100173253262857546D-01
    weight(39) = 0.462847965813144172959532492323D-01
    weight(40) = 0.454916279274181444797709969713D-01
    weight(41) = 0.445905581637565630601347100309D-01
    weight(42) = 0.435837245293234533768278609737D-01
    weight(43) = 0.424735151236535890073397679088D-01
    weight(44) = 0.412625632426235286101562974736D-01
    weight(45) = 0.399537411327203413866569261283D-01
    weight(46) = 0.385501531786156291289624969468D-01
    weight(47) = 0.370551285402400460404151018096D-01
    weight(48) = 0.354722132568823838106931467152D-01
    weight(49) = 0.338051618371416093915654821107D-01
    weight(50) = 0.320579283548515535854675043479D-01
    weight(51) = 0.302346570724024788679740598195D-01
    weight(52) = 0.283396726142594832275113052002D-01
    weight(53) = 0.263774697150546586716917926252D-01
    weight(54) = 0.243527025687108733381775504091D-01
    weight(55) = 0.222701738083832541592983303842D-01
    weight(56) = 0.201348231535302093723403167285D-01
    weight(57) = 0.179517157756973430850453020011D-01
    weight(58) = 0.157260304760247193219659952975D-01
    weight(59) = 0.134630478967186425980607666860D-01
    weight(60) = 0.111681394601311288185904930192D-01
    weight(61) = 0.884675982636394772303091465973D-02
    weight(62) = 0.650445796897836285611736039998D-02
    weight(63) = 0.414703326056246763528753572855D-02
    weight(64) = 0.178328072169643294729607914497D-02

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 to 20, 32 or 64.'
    stop

  end if

  return
end
subroutine legendre_set_cos ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_COS sets a Gauss-Legendre rule for COS(X) * F(X) on [-PI/2,PI/2].
!
!
!  Integration interval:
!
!    [ -PI/2, PI/2 ]
!
!  Weight function:
!
!    COS(X) * F(X)
!
!  Integral to approximate:
!
!    Integral ( -PI/2 <= X <= PI/2 ) COS(X) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Discussion:
!
!    The same rule can be used to approximate
!
!      Integral ( 0 <= X <= PI ) SIN(X) * F(X) dX
!
!    as
!
!      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) + PI/2 )
!
!  Modified:
!
!    23 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993, QA299.3E93, page 310.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1, 2, 4, 8 or 16.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) = 0.0D+00

    weight(1) = 2.0D+00

  else if ( norder == 2 ) then

    xtab(1) = - 0.68366739008990304094D+00
    xtab(2) =   0.68366739008990304094D+00

    weight(1) = 1.0D+00
    weight(2) = 1.0D+00

  else if ( norder == 4 ) then

    xtab(1) = - 1.1906765638948557415D+00
    xtab(2) = - 0.43928746686001514756D+00
    xtab(3) =   0.43928746686001514756D+00
    xtab(4) =   1.1906765638948557415D+00

    weight(1) = 0.22407061812762016065D+00
    weight(2) = 0.77592938187237983935D+00
    weight(3) = 0.77592938187237983935D+00
    weight(4) = 0.22407061812762016065D+00

  else if ( norder == 8 ) then

    xtab(1) = - 1.4414905401823575701D+00
    xtab(2) = - 1.1537256454567275850D+00
    xtab(3) = - 0.74346864787549244989D+00
    xtab(4) = - 0.25649650741623123020D+00
    xtab(5) =   0.25649650741623123020D+00
    xtab(6) =   0.74346864787549244989D+00
    xtab(7) =   1.1537256454567275850D+00
    xtab(8) =   1.4414905401823575701D+00

    weight(1) = 0.027535633513767011149D+00
    weight(2) = 0.14420409203022750950D+00
    weight(3) = 0.33626447785280459621D+00
    weight(4) = 0.49199579660320088314D+00
    weight(5) = 0.49199579660320088314D+00
    weight(6) = 0.33626447785280459621D+00
    weight(7) = 0.14420409203022750950D+00
    weight(8) = 0.027535633513767011149D+00

  else if ( norder == 16 ) then

    xtab( 1) = - 1.5327507132362304779D+00
    xtab( 2) = - 1.4446014873666514608D+00
    xtab( 3) = - 1.3097818904452936698D+00
    xtab( 4) = - 1.1330068786005003695D+00
    xtab( 5) = - 0.92027786206637096497D+00
    xtab( 6) = - 0.67861108097560545347D+00
    xtab( 7) = - 0.41577197673418943962D+00
    xtab( 8) = - 0.14003444424696773778D+00
    xtab( 9) =   0.14003444424696773778D+00
    xtab(10) =   0.41577197673418943962D+00
    xtab(11) =   0.67861108097560545347D+00
    xtab(12) =   0.92027786206637096497D+00
    xtab(13) =   1.1330068786005003695D+00
    xtab(14) =   1.3097818904452936698D+00
    xtab(15) =   1.4446014873666514608D+00
    xtab(16) =   1.5327507132362304779D+00

    weight( 1) = 0.0024194677567615628193D+00
    weight( 2) = 0.014115268156854008264D+00
    weight( 3) = 0.040437893946503669410D+00
    weight( 4) = 0.083026647573217742131D+00
    weight( 5) = 0.13834195526951273359D+00
    weight( 6) = 0.19741148870253455567D+00
    weight( 7) = 0.24763632094635522403D+00
    weight( 8) = 0.27661095764826050408D+00
    weight( 9) = 0.27661095764826050408D+00
    weight(10) = 0.24763632094635522403D+00
    weight(11) = 0.19741148870253455567D+00
    weight(12) = 0.13834195526951273359D+00
    weight(13) = 0.083026647573217742131D+00
    weight(14) = 0.040437893946503669410D+00
    weight(15) = 0.014115268156854008264D+00
    weight(16) = 0.0024194677567615628193D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_COS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_cos2 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_COS2 sets a Gauss-Legendre rule for COS(X) * F(X) on [0,PI/2].
!
!
!  Integration interval:
!
!    [ 0, PI/2 ]
!
!  Weight function:
!
!    COS(X) * F(X)
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= PI/2 ) COS(X) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Discussion:
!
!    The same rule can be used to approximate
!
!      Integral ( 0 <= X <= PI/2 ) SIN(X) * F(X) dX
!
!    as
!
!      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( PI/2 - XTAB(I) )
!
!  Modified:
!
!    24 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993, QA299.3E93, page 311.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 2, 4, 8 or 16.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 2 ) then

    xtab(1) = 0.26587388056307823382D+00
    xtab(2) = 1.0351526093171315182D+00

    weight(1) = 0.60362553280827113087D+00
    weight(2) = 0.39637446719172886913D+00

  else if ( norder == 4 ) then

    xtab(1) = 0.095669389196858636773D+00
    xtab(2) = 0.45240902327067096554D+00
    xtab(3) = 0.93185057672024082424D+00
    xtab(4) = 1.3564439599666466230D+00

    weight( 1) = 0.23783071419515504517D+00
    weight( 2) = 0.40265695523581253512D+00
    weight( 3) = 0.28681737948564715225D+00
    weight( 4) = 0.072694951083385267446D+00

  else if ( norder == 8 ) then

    xtab(1) = 0.029023729768913933432D+00
    xtab(2) = 0.14828524404581819442D+00
    xtab(3) = 0.34531111151664787488D+00
    xtab(4) = 0.59447696797658360178D+00
    xtab(5) = 0.86538380686123504827D+00
    xtab(6) = 1.1263076093187456632D+00
    xtab(7) = 1.3470150460281258016D+00
    xtab(8) = 1.5015603622059195568D+00

    weight( 1) = 0.073908998095117384985D+00
    weight( 2) = 0.16002993702338006099D+00
    weight( 3) = 0.21444434341803549108D+00
    weight( 4) = 0.21979581268851903339D+00
    weight( 5) = 0.17581164478209568886D+00
    weight( 6) = 0.10560448025308322171D+00
    weight( 7) = 0.042485497299217201089D+00
    weight( 8) = 0.0079192864405519178899D+00

  else if ( norder == 16 ) then

    xtab( 1) = 0.0080145034906295973494D+00
    xtab( 2) = 0.041893031354246254797D+00
    xtab( 3) = 0.10149954486757579459D+00
    xtab( 4) = 0.18463185923836617507D+00
    xtab( 5) = 0.28826388487760574589D+00
    xtab( 6) = 0.40870579076464794191D+00
    xtab( 7) = 0.54176054986913847463D+00
    xtab( 8) = 0.68287636658719416893D+00
    xtab( 9) = 0.82729287620416833520D+00
    xtab(10) = 0.97018212594829367065D+00
    xtab(11) = 1.1067865150286247873D+00
    xtab(12) = 1.2325555697227748824D+00
    xtab(13) = 1.3432821921580721861D+00
    xtab(14) = 1.4352370549295032923D+00
    xtab(15) = 1.5052970876794669248D+00
    xtab(16) = 1.5510586944086135769D+00

    weight( 1) = 0.020528714977215248902D+00
    weight( 2) = 0.046990919853597958123D+00
    weight( 3) = 0.071441021312218541698D+00
    weight( 4) = 0.092350338329243052271D+00
    weight( 5) = 0.10804928026816236935D+00
    weight( 6) = 0.11698241243306261791D+00
    weight( 7) = 0.11812395361762037649D+00
    weight( 8) = 0.11137584940420091049D+00
    weight( 9) = 0.097778236145946543110D+00
    weight(10) = 0.079418758985944482077D+00
    weight(11) = 0.059039620053768691402D+00
    weight(12) = 0.039458876783728165671D+00
    weight(13) = 0.022987785677206847531D+00
    weight(14) = 0.011010405600421536861D+00
    weight(15) = 0.0038123928030499915653D+00
    weight(16) = 0.00065143375461266656171D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_COS2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_log ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_LOG sets a Gauss-Legendre rule for - LOG(X) * F(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    - LOG(X) * F(X)
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) - LOG(X) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    05 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964, page 920.
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993, QA299.3E93, page 309.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 through 8, or 16.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) = 0.25D+00

    weight(1) = 1.0D+00

  else if ( norder == 2 ) then

    xtab(1) = 0.112008806166976182957205488948D+00
    xtab(2) = 0.602276908118738102757080225338D+00

    weight(1) = 0.718539319030384440665510200891D+00
    weight(2) = 0.281460680969615559334489799109D+00

  else if ( norder == 3 ) then

    xtab(1) = 0.0638907930873254049961166031363D+00
    xtab(2) = 0.368997063715618765546197645857D+00
    xtab(3) = 0.766880303938941455423682659817D+00

    weight(1) = 0.513404552232363325129300497567D+00
    weight(2) = 0.391980041201487554806287180966D+00
    weight(3) = 0.0946154065661491200644123214672D+00

  else if ( norder == 4 ) then

    xtab(1) = 0.0414484801993832208033213101564D+00
    xtab(2) = 0.245274914320602251939675759523D+00
    xtab(3) = 0.556165453560275837180184354376D+00
    xtab(4) = 0.848982394532985174647849188085D+00

    weight(1) = 0.383464068145135124850046522343D+00
    weight(2) = 0.386875317774762627336008234554D+00
    weight(3) = 0.190435126950142415361360014547D+00
    weight(4) = 0.0392254871299598324525852285552D+00

  else if ( norder == 5 ) then

    xtab(1) = 0.0291344721519720533037267621154D+00
    xtab(2) = 0.173977213320897628701139710829D+00
    xtab(3) = 0.411702520284902043174931924646D+00
    xtab(4) = 0.677314174582820380701802667998D+00
    xtab(5) = 0.894771361031008283638886204455D+00

    weight(1) = 0.297893471782894457272257877888D+00
    weight(2) = 0.349776226513224180375071870307D+00
    weight(3) = 0.234488290044052418886906857943D+00
    weight(4) = 0.0989304595166331469761807114404D+00
    weight(5) = 0.0189115521431957964895826824218D+00

  else if ( norder == 6 ) then

    xtab(1) = 0.0216340058441169489956958558537D+00
    xtab(2) = 0.129583391154950796131158505009D+00
    xtab(3) = 0.314020449914765508798248188420D+00
    xtab(4) = 0.538657217351802144548941893993D+00
    xtab(5) = 0.756915337377402852164544156139D+00
    xtab(6) = 0.922668851372120237333873231507D+00

    weight(1) = 0.238763662578547569722268303330D+00
    weight(2) = 0.308286573273946792969383109211D+00
    weight(3) = 0.245317426563210385984932540188D+00
    weight(4) = 0.142008756566476685421345576030D+00
    weight(5) = 0.0554546223248862900151353549662D+00
    weight(6) = 0.0101689586929322758869351162755D+00

  else if ( norder == 7 ) then

    xtab(1) = 0.0167193554082585159416673609320D+00
    xtab(2) = 0.100185677915675121586885031757D+00
    xtab(3) = 0.246294246207930599046668547239D+00
    xtab(4) = 0.433463493257033105832882482601D+00
    xtab(5) = 0.632350988047766088461805812245D+00
    xtab(6) = 0.811118626740105576526226796782D+00
    xtab(7) = 0.940848166743347721760134113379D+00

    weight(1) = 0.196169389425248207525427377585D+00
    weight(2) = 0.270302644247272982145271719533D+00
    weight(3) = 0.239681873007690948308072785252D+00
    weight(4) = 0.165775774810432906560869687736D+00
    weight(5) = 0.0889432271376579644357238403458D+00
    weight(6) = 0.0331943043565710670254494111034D+00
    weight(7) = 0.00593278701512592399918517844468D+00

  else if ( norder == 8 ) then

    xtab(1) = 0.0133202441608924650122526725243D+00
    xtab(2) = 0.0797504290138949384098277291424D+00
    xtab(3) = 0.197871029326188053794476159516D+00
    xtab(4) = 0.354153994351909419671463603538D+00
    xtab(5) = 0.529458575234917277706149699996D+00
    xtab(6) = 0.701814529939099963837152670310D+00
    xtab(7) = 0.849379320441106676048309202301D+00
    xtab(8) = 0.953326450056359788767379678514D+00

    weight(1) = 0.164416604728002886831472568326D+00
    weight(2) = 0.237525610023306020501348561960D+00
    weight(3) = 0.226841984431919126368780402936D+00
    weight(4) = 0.175754079006070244988056212006D+00
    weight(5) = 0.112924030246759051855000442086D+00
    weight(6) = 0.0578722107177820723985279672940D+00
    weight(7) = 0.0209790737421329780434615241150D+00
    weight(8) = 0.00368640710402761901335232127647D+00

  else if ( norder == 16 ) then

    xtab( 1) = 0.00389783448711591592405360527037D+00
    xtab( 2) = 0.0230289456168732398203176309848D+00
    xtab( 3) = 0.0582803983062404123483532298394D+00
    xtab( 4) = 0.108678365091054036487713613051D+00
    xtab( 5) = 0.172609454909843937760843776232D+00
    xtab( 6) = 0.247937054470578495147671753047D+00
    xtab( 7) = 0.332094549129917155984755859320D+00
    xtab( 8) = 0.422183910581948600115088366745D+00
    xtab( 9) = 0.515082473381462603476277704052D+00
    xtab(10) = 0.607556120447728724086384921709D+00
    xtab(11) = 0.696375653228214061156318166581D+00
    xtab(12) = 0.778432565873265405203868167732D+00
    xtab(13) = 0.850850269715391083233822761319D+00
    xtab(14) = 0.911086857222271905418818994060D+00
    xtab(15) = 0.957025571703542157591520509383D+00
    xtab(16) = 0.987047800247984476758697436516D+00

    weight( 1) = 0.0607917100435912328511733871235D+00
    weight( 2) = 0.102915677517582144387691736210D+00
    weight( 3) = 0.122355662046009193557547513197D+00
    weight( 4) = 0.127569246937015988717042209239D+00
    weight( 5) = 0.123013574600070915423123365137D+00
    weight( 6) = 0.111847244855485722621848903429D+00
    weight( 7) = 0.0965963851521243412529681650802D+00
    weight( 8) = 0.0793566643514731387824416770520D+00
    weight( 9) = 0.0618504945819652070951360853113D+00
    weight(10) = 0.0454352465077266686288299526509D+00
    weight(11) = 0.0310989747515818064092528627927D+00
    weight(12) = 0.0194597659273608420780860268669D+00
    weight(13) = 0.0107762549632055256455393162159D+00
    weight(14) = 0.00497254289008764171250524951646D+00
    weight(15) = 0.00167820111005119451503546419059D+00
    weight(16) = 0.000282353764668436321778085987413D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_LOG - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_sqrtx_01 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_SQRTX_01 sets a Gauss-Legendre rule for SQRT(X) * F(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    SQRT ( X )
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) SQRT ( X ) * F(X) dX =
!    Integral ( 0 <= Y <= 1 ) 2 * Y**2 * F(Y**2) dY.
!    (using Y = SQRT(X) )
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    23 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 30th Edition, 2000, page 696.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  integer norder2
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) xtab2(2*norder+1)
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) weight2(2*norder+1)
!
  norder2 = 2 * norder + 1

  call legendre_set ( norder2, xtab2, weight2 )

  xtab(1:norder) = xtab2(norder+2:2*norder+1)**2
  weight(1:norder) = 2.0D+00 * weight2(norder+2:2*norder+1) * xtab(1:norder)

  return
end
subroutine legendre_set_sqrtx2_01 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_SQRTX2_01 sets a Gauss-Legendre rule for F(X) / SQRT(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    1 / SQRT ( X )
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) / SQRT ( X ) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 30th Edition, 2000, page 696.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  integer norder2
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) xtab2(2*norder+1)
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) weight2(2*norder+1)
!
  norder2 = 2 * norder + 1

  call legendre_set ( norder2, xtab2, weight2 )

  xtab(1:norder) = xtab2(norder+2:2*norder+1)**2
  weight(1:norder) = 2.0D+00 * weight2(norder+2:2*norder+1)

  return
end
subroutine legendre_set_x0_01 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_X0_01 sets a Gauss-Legendre rule for F(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964, page 921.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 8.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =   0.5D+00

    weight(1) = 1.0D+00

  else if ( norder == 2 ) then

    xtab(1) = 0.2113248654D+00
    xtab(2) = 0.7886751346D+00

    weight(1) = 0.5D+00
    weight(2) = 0.5D+00

  else if ( norder == 3 ) then

    xtab(1) = 0.1127016654D+00
    xtab(2) = 0.5000000000D+00
    xtab(3) = 0.8872983346D+00

    weight(1) = 5.0D+00 / 18.0D+00
    weight(2) = 8.0D+00 / 18.0D+00
    weight(3) = 5.0D+00 / 18.0D+00

  else if ( norder == 4 ) then

    xtab(1) = 0.0694318442D+00
    xtab(2) = 0.3300094782D+00
    xtab(3) = 0.6699905218D+00
    xtab(4) = 0.9305681558D+00

    weight(1) = 0.1739274226D+00
    weight(2) = 0.3260725774D+00
    weight(3) = 0.3260725774D+00
    weight(4) = 0.1739274226D+00

  else if ( norder == 5 ) then

    xtab(1) = 0.0469100770D+00
    xtab(2) = 0.2307653449D+00
    xtab(3) = 0.5000000000D+00
    xtab(4) = 0.7692346551D+00
    xtab(5) = 0.9530899230D+00

    weight(1) = 0.1184634425D+00
    weight(2) = 0.2393143352D+00
    weight(3) = 0.2844444444D+00
    weight(4) = 0.2393143352D+00
    weight(5) = 0.1184634425D+00

  else if ( norder == 6 ) then

    xtab(1) = 0.0337652429D+00
    xtab(2) = 0.1693953068D+00
    xtab(3) = 0.3806904070D+00
    xtab(4) = 0.6193095930D+00
    xtab(5) = 0.8306046932D+00
    xtab(6) = 0.9662347571D+00

    weight(1) = 0.0856622462D+00
    weight(2) = 0.1803807865D+00
    weight(3) = 0.2339569673D+00
    weight(4) = 0.2339569673D+00
    weight(5) = 0.1803807865D+00
    weight(6) = 0.0856622462D+00

  else if ( norder == 7 ) then

    xtab(1) = 0.0254460438D+00
    xtab(2) = 0.1292344072D+00
    xtab(3) = 0.2970774243D+00
    xtab(4) = 0.5000000000D+00
    xtab(5) = 0.7029225757D+00
    xtab(6) = 0.8707655928D+00
    xtab(7) = 0.9745539562D+00

    weight(1) = 0.0647424831D+00
    weight(2) = 0.1398526957D+00
    weight(3) = 0.1909150253D+00
    weight(4) = 0.2089795918D+00
    weight(5) = 0.1909150253D+00
    weight(6) = 0.1398526957D+00
    weight(7) = 0.0647424831D+00

  else if ( norder == 8 ) then

    xtab(1) = 0.0198550718D+00
    xtab(2) = 0.1016667613D+00
    xtab(3) = 0.2372337950D+00
    xtab(4) = 0.4082826788D+00
    xtab(5) = 0.5917173212D+00
    xtab(6) = 0.7627662050D+00
    xtab(7) = 0.8983332387D+00
    xtab(8) = 0.9801449282D+00

    weight(1) = 0.0506142681D+00
    weight(2) = 0.1111905172D+00
    weight(3) = 0.1568533229D+00
    weight(4) = 0.1813418917D+00
    weight(5) = 0.1813418917D+00
    weight(6) = 0.1568533229D+00
    weight(7) = 0.1111905172D+00
    weight(8) = 0.0506142681D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X0_01 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_x1 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 + X
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966, Table #3.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 9.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =  0.333333333333333333333333333333D+00

    weight(1) = 2.0D+00

  else if ( norder == 2 ) then

    xtab(1) = -0.289897948556635619639456814941D+00
    xtab(2) =  0.689897948556635619639456814941D+00

    weight(1) =  0.727834473024091322422523991699D+00
    weight(2) =  1.27216552697590867757747600830D+00

  else if ( norder == 3 ) then

    xtab(1) = -0.575318923521694112050483779752D+00
    xtab(2) =  0.181066271118530578270147495862D+00
    xtab(3) =  0.822824080974592105208907712461D+00

    weight(1) =  0.279307919605816490135525088716D+00
    weight(2) =  0.916964425438344986775682378225D+00
    weight(3) =  0.803727654955838523088792533058D+00

  else if ( norder == 4 ) then

    xtab(1) = -0.720480271312438895695825837750D+00
    xtab(2) = -0.167180864737833640113395337326D+00
    xtab(3) =  0.446313972723752344639908004629D+00
    xtab(4) =  0.885791607770964635613757614892D+00

    weight(1) =  0.124723883800032328695500588386D+00
    weight(2) =  0.519390190432929763305824811559D+00
    weight(3) =  0.813858272041085443165617903743D+00
    weight(4) =  0.542027653725952464833056696312D+00

  else if ( norder == 5 ) then

    xtab(1) = -0.802929828402347147753002204224D+00
    xtab(2) = -0.390928546707272189029229647442D+00
    xtab(3) =  0.124050379505227711989974959990D+00
    xtab(4) =  0.603973164252783654928415726409D+00
    xtab(5) =  0.920380285897062515318386619813D+00

    weight(1) =  0.0629916580867691047411692662740D+00
    weight(2) =  0.295635480290466681402532877367D+00
    weight(3) =  0.585547948338679234792151477424D+00
    weight(4) =  0.668698552377478261966702492391D+00
    weight(5) =  0.387126360906606717097443886545D+00

  else if ( norder == 6 ) then

    xtab(1) = -0.853891342639482229703747931639D+00
    xtab(2) = -0.538467724060109001833766720231D+00
    xtab(3) = -0.117343037543100264162786683611D+00
    xtab(4) =  0.326030619437691401805894055838D+00
    xtab(5) =  0.703842800663031416300046295008D+00
    xtab(6) =  0.941367145680430216055899446174D+00

    weight(1) =  0.0349532072544381270240692132496D+00
    weight(2) =  0.175820662202035902032706497222D+00
    weight(3) =  0.394644603562621056482338042193D+00
    weight(4) =  0.563170215152795712476307356284D+00
    weight(5) =  0.542169988926074467362761586552D+00
    weight(6) =  0.289241322902034734621817304499D+00

  else if ( norder == 7 ) then

    xtab(1) = -0.887474878926155707068695617935D+00
    xtab(2) = -0.639518616526215270024840114382D+00
    xtab(3) = -0.294750565773660725252184459658D+00
    xtab(4) =  0.0943072526611107660028971153047D+00
    xtab(5) =  0.468420354430821063046421216613D+00
    xtab(6) =  0.770641893678191536180719525865D+00
    xtab(7) =  0.955041227122575003782349000858D+00

    weight(1) =  0.0208574488112296163587654972151D+00
    weight(2) =  0.109633426887493901777324193433D+00
    weight(3) =  0.265538785861965879934591955055D+00
    weight(4) =  0.428500262783494679963649011999D+00
    weight(5) =  0.509563589198353307674937943100D+00
    weight(6) =  0.442037032763498409684482945478D+00
    weight(7) =  0.223869453693964204606248453720D+00

  else if ( norder == 8 ) then

    xtab(1) = -0.910732089420060298533757956283D+00
    xtab(2) = -0.711267485915708857029562959544D+00
    xtab(3) = -0.426350485711138962102627520502D+00
    xtab(4) = -0.0903733696068532980645444599064D+00
    xtab(5) =  0.256135670833455395138292079035D+00
    xtab(6) =  0.571383041208738483284917464837D+00
    xtab(7) =  0.817352784200412087992517083851D+00
    xtab(8) =  0.964440169705273096373589797925D+00

    weight(1) =  0.0131807657689951954189692640444D+00
    weight(2) =  0.0713716106239448335742111888042D+00
    weight(3) =  0.181757278018795592332221684383D+00
    weight(4) =  0.316798397969276640481632757440D+00
    weight(5) =  0.424189437743720042818124385645D+00
    weight(6) =  0.450023197883549464687088394417D+00
    weight(7) =  0.364476094545494505382889847132D+00
    weight(8) =  0.178203217446223725304862478136D+00

  else if ( norder == 9 ) then

    xtab(1) = -0.927484374233581078117671398464D+00
    xtab(2) = -0.763842042420002599615429776011D+00
    xtab(3) = -0.525646030370079229365386614293D+00
    xtab(4) = -0.236234469390588049278459503207D+00
    xtab(5) =  0.0760591978379781302337137826389D+00
    xtab(6) =  0.380664840144724365880759065541D+00
    xtab(7) =  0.647766687674009436273648507855D+00
    xtab(8) =  0.851225220581607910728163628088D+00
    xtab(9) =  0.971175180702246902734346518378D+00

    weight(1) =  0.00872338834309252349019620448007D+00
    weight(2) =  0.0482400171391415162069086091476D+00
    weight(3) =  0.127219285964216005046760427743D+00
    weight(4) =  0.233604781180660442262926091607D+00
    weight(5) =  0.337433287379681397577000079834D+00
    weight(6) =  0.401235236773473158616600898930D+00
    weight(7) =  0.394134968689382820640692081477D+00
    weight(8) =  0.304297020437232650320317215016D+00
    weight(9) =  0.145112014093119485838598391765D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X1 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal input value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_x1_01 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_X1_01 sets a Gauss-Legendre rule for X * F(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    X
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) X * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964, page 921.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 8.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =   0.6666666667D+00

    weight(1) = 0.5000000000D+00

  else if ( norder == 2 ) then

    xtab(1) = 0.3550510257D+00
    xtab(2) = 0.8449489743D+00

    weight(1) = 0.1819586183D+00
    weight(2) = 0.3180413817D+00

  else if ( norder == 3 ) then

    xtab(1) = 0.2123405382D+00
    xtab(2) = 0.5905331356D+00
    xtab(3) = 0.9114120405D+00

    weight(1) = 0.0698269799D+00
    weight(2) = 0.2292411064D+00
    weight(3) = 0.2009319137D+00

  else if ( norder == 4 ) then

    xtab(1) = 0.1397598643D+00
    xtab(2) = 0.4164095676D+00
    xtab(3) = 0.7231569864D+00
    xtab(4) = 0.9428958039D+00

    weight(1) = 0.0311809710D+00
    weight(2) = 0.1298475476D+00
    weight(3) = 0.2034645680D+00
    weight(4) = 0.1355069134D+00

  else if ( norder == 5 ) then

    xtab(1) = 0.0985350858D+00
    xtab(2) = 0.3045357266D+00
    xtab(3) = 0.5620251898D+00
    xtab(4) = 0.8019865821D+00
    xtab(5) = 0.9601901429D+00

    weight(1) = 0.0157479145D+00
    weight(2) = 0.0739088701D+00
    weight(3) = 0.1463888701D+00
    weight(4) = 0.1671746381D+00
    weight(5) = 0.0967815902D+00

  else if ( norder == 6 ) then

    xtab(1) = 0.0730543287D+00
    xtab(2) = 0.2307661380D+00
    xtab(3) = 0.4413284812D+00
    xtab(4) = 0.6630153097D+00
    xtab(5) = 0.8519214003D+00
    xtab(6) = 0.9706835728D+00

    weight(1) = 0.0087383108D+00
    weight(2) = 0.0439551656D+00
    weight(3) = 0.0986611509D+00
    weight(4) = 0.1407925538D+00
    weight(5) = 0.1355424972D+00
    weight(6) = 0.0723103307D+00

  else if ( norder == 7 ) then

    xtab(1) = 0.0562625605D+00
    xtab(2) = 0.1802406917D+00
    xtab(3) = 0.3526247171D+00
    xtab(4) = 0.5471536263D+00
    xtab(5) = 0.7342101772D+00
    xtab(6) = 0.8853209468D+00
    xtab(7) = 0.9775206136D+00

    weight(1) = 0.0052143622D+00
    weight(2) = 0.0274083567D+00
    weight(3) = 0.0663846965D+00
    weight(4) = 0.1071250657D+00
    weight(5) = 0.1273908973D+00
    weight(6) = 0.1105092582D+00
    weight(7) = 0.0559673634D+00

  else if ( norder == 8 ) then

    xtab(1) = 0.0446339553D+00
    xtab(2) = 0.1443662570D+00
    xtab(3) = 0.2868247571D+00
    xtab(4) = 0.4548133152D+00
    xtab(5) = 0.6280678354D+00
    xtab(6) = 0.7856915206D+00
    xtab(7) = 0.9086763921D+00
    xtab(8) = 0.9822200849D+00

    weight(1) = 0.0032951914D+00
    weight(2) = 0.0178429027D+00
    weight(3) = 0.0454393195D+00
    weight(4) = 0.0791995995D+00
    weight(5) = 0.1060473594D+00
    weight(6) = 0.1125057995D+00
    weight(7) = 0.0911190236D+00
    weight(8) = 0.0445508044D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X1_01 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_x2 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_X2 sets a Gauss-Legendre rule for ( 1 + X )**2 * F(X) on [-1,1].
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    ( 1 + X )**2
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) ( 1 + X )**2 * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966, Table #3.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 9.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =  0.5D+00

    weight(1) =  2.66666666666666666666666666666D+00

  else if ( norder == 2 ) then

    xtab(1) = -0.0883036880224505775998524725910D+00
    xtab(2) =  0.754970354689117244266519139258D+00

    weight(1) =  0.806287056638603444666851075928D+00
    weight(2) =  1.86037961002806322199981559074D+00

  else if ( norder == 3 ) then

    xtab(1) = -0.410004419776996766244796955168D+00
    xtab(2) =  0.305992467923296230556472913192D+00
    xtab(3) =  0.854011951853700535688324041976D+00

    weight(1) =  0.239605624068645584091811926047D+00
    weight(2) =  1.16997015407892817602809616291D+00
    weight(3) =  1.25709088851909290654675857771D+00

  else if ( norder == 4 ) then

    xtab(1) = -0.591702835793545726606755921586D+00
    xtab(2) = -0.0340945902087350046811467387661D+00
    xtab(3) =  0.522798524896275389882037174551D+00
    xtab(4) =  0.902998901106005341405865485802D+00

    weight(1) =  0.0828179259993445222751812523731D+00
    weight(2) =  0.549071097383384602539010760334D+00
    weight(3) =  1.14767031839371367238662411421D+00
    weight(4) =  0.887107324890223869465850539752D+00

  else if ( norder == 5 ) then

    xtab(1) = -0.702108425894032836232448374820D+00
    xtab(2) = -0.268666945261773544694327777841D+00
    xtab(3) =  0.220227225868961343518209179230D+00
    xtab(4) =  0.653039358456608553790815164028D+00
    xtab(5) =  0.930842120163569816951085142737D+00

    weight(1) =  0.0329106016247920636689299329544D+00
    weight(2) =  0.256444805783695354037991444453D+00
    weight(3) =  0.713601289772720001490035944563D+00
    weight(4) =  1.00959169519929190423066348132D+00
    weight(5) =  0.654118274286167343239045863379D+00

  else if ( norder == 6 ) then

    xtab(1) = -0.773611232355123732602532012021D+00
    xtab(2) = -0.431362254623427837535325249187D+00
    xtab(3) = -0.0180728263295041680220798103354D+00
    xtab(4) =  0.395126163954217534500188844163D+00
    xtab(5) =  0.736872116684029732026178298518D+00
    xtab(6) =  0.948190889812665614490712786006D+00

    weight(1) =  0.0146486064549543818622276447204D+00
    weight(2) =  0.125762377479560410622810097040D+00
    weight(3) =  0.410316569036929681761034600615D+00
    weight(4) =  0.756617493988329628546336413760D+00
    weight(5) =  0.859011997894245060846045458784D+00
    weight(6) =  0.500309621812647503028212451747D+00

  else if ( norder == 7 ) then

    xtab(1) = -0.822366333126005527278634734418D+00
    xtab(2) = -0.547034493182875002223997992852D+00
    xtab(3) = -0.200043026557985860387937545780D+00
    xtab(4) =  0.171995710805880507163425502299D+00
    xtab(5) =  0.518891747903884926692601716998D+00
    xtab(6) =  0.793821941703901970495546427988D+00
    xtab(7) =  0.959734452453198985538996625765D+00

    weight(1) =  0.00714150426951365443207221475404D+00
    weight(2) =  0.0653034050584375560578544725498D+00
    weight(3) =  0.235377690316228918725962815880D+00
    weight(4) =  0.505171029671130381676271523850D+00
    weight(5) =  0.733870426238362032891332767175D+00
    weight(6) =  0.725590596901489156295739839779D+00
    weight(7) =  0.394212014211504966587433032679D+00

  else if ( norder == 8 ) then

    xtab(1) = -0.857017929919813794402037235698D+00
    xtab(2) = -0.631543407166567521509503573952D+00
    xtab(3) = -0.339104543648722903660229021109D+00
    xtab(4) = -0.0111941563689783438801237300122D+00
    xtab(5) =  0.316696017045595559454075475675D+00
    xtab(6) =  0.609049663022520165351466780939D+00
    xtab(7) =  0.834198765028697794599267293239D+00
    xtab(8) =  0.967804480896157932935972899807D+00

    weight(1) =  0.00374814227227757804631954025851D+00
    weight(2) =  0.0357961737041152639660521680263D+00
    weight(3) =  0.137974910241879862433949246199D+00
    weight(4) =  0.326515411108352185491692769217D+00
    weight(5) =  0.547577467373226177976217604887D+00
    weight(6) =  0.682278153375510121675529810121D+00
    weight(7) =  0.614544746137780998436053880546D+00
    weight(8) =  0.318231662453524478640851647411D+00

  else if ( norder == 9 ) then

    xtab(1) = -0.882491728426548422828684254270D+00
    xtab(2) = -0.694873684026474640346360850039D+00
    xtab(3) = -0.446537143480670863635920316400D+00
    xtab(4) = -0.159388112702326252531544826624D+00
    xtab(5) =  0.141092709224374414981503995427D+00
    xtab(6) =  0.428217823321559204544020866175D+00
    xtab(7) =  0.676480966471850715860378175342D+00
    xtab(8) =  0.863830940812464825046988286026D+00
    xtab(9) =  0.973668228805771018909618924364D+00

    weight(1) =  0.00209009877215570354392734918986D+00
    weight(2) =  0.0205951891648697848186537272448D+00
    weight(3) =  0.0832489326348178964194106978875D+00
    weight(4) =  0.210746247220398685903797568021D+00
    weight(5) =  0.388325022916052063676224499399D+00
    weight(6) =  0.554275165518437673725822282791D+00
    weight(7) =  0.621388553284444032628761363828D+00
    weight(8) =  0.523916296267173054255512857631D+00
    weight(9) =  0.262081160888317771694556320674D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal input value of NORDER = ', norder
    stop

  end if

  return
end
subroutine legendre_set_x2_01 ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LEGENDRE_SET_X2_01 sets a Gauss-Legendre rule for X**2 * F(X) on [0,1].
!
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    X**2
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) X*X * F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964, page 921.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 8.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =   0.75D+00

    weight(1) = 1.0D+00 / 3.0D+00

  else if ( norder == 2 ) then

    xtab(1) = 0.4558481560D+00
    xtab(2) = 0.8774851773D+00

    weight(1) = 0.1007858821D+00
    weight(2) = 0.2325474513D+00

  else if ( norder == 3 ) then

    xtab(1) = 0.2949977901D+00
    xtab(2) = 0.6529962340D+00
    xtab(3) = 0.9270059759D+00

    weight(1) = 0.0299507030D+00
    weight(2) = 0.1462462693D+00
    weight(3) = 0.1571363611D+00

  else if ( norder == 4 ) then

    xtab(1) = 0.2041485821D+00
    xtab(2) = 0.4829527049D+00
    xtab(3) = 0.7613992624D+00
    xtab(4) = 0.9514994506D+00

    weight(1) = 0.0103522408D+00
    weight(2) = 0.0686338872D+00
    weight(3) = 0.1434587898D+00
    weight(4) = 0.1108884156D+00

  else if ( norder == 5 ) then

    xtab(1) = 0.1489457871D+00
    xtab(2) = 0.3656665274D+00
    xtab(3) = 0.6101136129D+00
    xtab(4) = 0.8265196792D+00
    xtab(5) = 0.9654210601D+00

    weight(1) = 0.0041138252D+00
    weight(2) = 0.0320556007D+00
    weight(3) = 0.0892001612D+00
    weight(4) = 0.1261989619D+00
    weight(5) = 0.0817647843D+00

  else if ( norder == 6 ) then

    xtab(1) = 0.1131943838D+00
    xtab(2) = 0.2843188727D+00
    xtab(3) = 0.4909635868D+00
    xtab(4) = 0.6975630820D+00
    xtab(5) = 0.8684360583D+00
    xtab(6) = 0.9740954449D+00

    weight(1) = 0.0018310758D+00
    weight(2) = 0.0157202972D+00
    weight(3) = 0.0512895711D+00
    weight(4) = 0.0945771867D+00
    weight(5) = 0.1073764997D+00
    weight(6) = 0.0625387027D+00

  else if ( norder == 7 ) then

    xtab(1) = 0.0888168334D+00
    xtab(2) = 0.2264827534D+00
    xtab(3) = 0.3999784867D+00
    xtab(4) = 0.5859978554D+00
    xtab(5) = 0.7594458740D+00
    xtab(6) = 0.8969109709D+00
    xtab(7) = 0.9798672262D+00

    weight(1) = 0.0008926880D+00
    weight(2) = 0.0081629256D+00
    weight(3) = 0.0294222113D+00
    weight(4) = 0.0631463787D+00
    weight(5) = 0.0917338033D+00
    weight(6) = 0.0906988246D+00
    weight(7) = 0.0492765018D+00

  else if ( norder == 8 ) then

    xtab(1) = 0.0714910350D+00
    xtab(2) = 0.1842282964D+00
    xtab(3) = 0.3304477282D+00
    xtab(4) = 0.4944029218D+00
    xtab(5) = 0.6583480085D+00
    xtab(6) = 0.8045248315D+00
    xtab(7) = 0.9170993825D+00
    xtab(8) = 0.9839022404D+00

    weight(1) = 0.0004685178D+00
    weight(2) = 0.0044745217D+00
    weight(3) = 0.0172468638D+00
    weight(4) = 0.0408144264D+00
    weight(5) = 0.0684471834D+00
    weight(6) = 0.0852847692D+00
    weight(7) = 0.0768180933D+00
    weight(8) = 0.0397789578D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X2_01 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    stop

  end if

  return
end
subroutine lobatto_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! LOBATTO_SET sets abscissas and weights for Lobatto quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-3).
!
!  Note:
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas of the rule.
!
!  Modified:
!
!    20 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 2 and 20.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas for the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 2 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =    1.0D+00

    weight(1) =  1.0D+00
    weight(2) =  1.0D+00

  else if ( norder == 3 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =    0.0D+00
    xtab(3) =    1.0D+00

    weight(1) =  1.0D+00 / 3.0D+00
    weight(2) =  4.0D+00 / 3.0D+00
    weight(3) =  1.0D+00 / 3.0D+00

  else if ( norder == 4 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.447213595499957939281834733746D+00
    xtab(3) =    0.447213595499957939281834733746D+00
    xtab(4) =    1.0D+00

    weight(1) =  1.0D+00 / 6.0D+00
    weight(2) =  5.0D+00 / 6.0D+00
    weight(3) =  5.0D+00 / 6.0D+00
    weight(4) =  1.0D+00 / 6.0D+00

  else if ( norder == 5 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.654653670707977143798292456247D+00
    xtab(3) =    0.0D+00
    xtab(4) =    0.654653670707977143798292456247D+00
    xtab(5) =    1.0D+00

    weight(1) =  9.0D+00 / 90.0D+00
    weight(2) = 49.0D+00 / 90.0D+00
    weight(3) = 64.0D+00 / 90.0D+00
    weight(4) = 49.0D+00 / 90.0D+00
    weight(5) =  9.0D+00 / 90.0D+00

  else if ( norder == 6 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.765055323929464692851002973959D+00
    xtab(3) =  - 0.285231516480645096314150994041D+00
    xtab(4) =    0.285231516480645096314150994041D+00
    xtab(5) =    0.765055323929464692851002973959D+00
    xtab(6) =    1.0D+00

    weight(1) =  0.066666666666666666666666666667D+00
    weight(2) =  0.378474956297846980316612808212D+00
    weight(3) =  0.554858377035486353016720525121D+00
    weight(4) =  0.554858377035486353016720525121D+00
    weight(5) =  0.378474956297846980316612808212D+00
    weight(6) =  0.066666666666666666666666666667D+00

  else if ( norder == 7 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.830223896278566929872032213967D+00
    xtab(3) =  - 0.468848793470714213803771881909D+00
    xtab(4) =    0.0D+00
    xtab(5) =    0.468848793470714213803771881909D+00
    xtab(6) =    0.830223896278566929872032213967D+00
    xtab(7) =    1.0D+00

    weight(1) =  0.476190476190476190476190476190D-01
    weight(2) =  0.276826047361565948010700406290D+00
    weight(3) =  0.431745381209862623417871022281D+00
    weight(4) =  0.487619047619047619047619047619D+00
    weight(5) =  0.431745381209862623417871022281D+00
    weight(6) =  0.276826047361565948010700406290D+00
    weight(7) =  0.476190476190476190476190476190D-01

  else if ( norder == 8 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.871740148509606615337445761221D+00
    xtab(3) =  - 0.591700181433142302144510731398D+00
    xtab(4) =  - 0.209299217902478868768657260345D+00
    xtab(5) =    0.209299217902478868768657260345D+00
    xtab(6) =    0.591700181433142302144510731398D+00
    xtab(7) =    0.871740148509606615337445761221D+00
    xtab(8) =    1.0D+00

    weight(1) =  0.357142857142857142857142857143D-01
    weight(2) =  0.210704227143506039382991065776D+00
    weight(3) =  0.341122692483504364764240677108D+00
    weight(4) =  0.412458794658703881567052971402D+00
    weight(5) =  0.412458794658703881567052971402D+00
    weight(6) =  0.341122692483504364764240677108D+00
    weight(7) =  0.210704227143506039382991065776D+00
    weight(8) =  0.357142857142857142857142857143D-01

  else if ( norder == 9 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.899757995411460157312345244418D+00
    xtab(3) =  - 0.677186279510737753445885427091D+00
    xtab(4) =  - 0.363117463826178158710752068709D+00
    xtab(5) =    0.0D+00
    xtab(6) =    0.363117463826178158710752068709D+00
    xtab(7) =    0.677186279510737753445885427091D+00
    xtab(8) =    0.899757995411460157312345244418D+00
    xtab(9) =    1.0D+00

    weight(1) =  0.277777777777777777777777777778D-01
    weight(2) =  0.165495361560805525046339720029D+00
    weight(3) =  0.274538712500161735280705618579D+00
    weight(4) =  0.346428510973046345115131532140D+00
    weight(5) =  0.371519274376417233560090702948D+00
    weight(6) =  0.346428510973046345115131532140D+00
    weight(7) =  0.274538712500161735280705618579D+00
    weight(8) =  0.165495361560805525046339720029D+00
    weight(9) =  0.277777777777777777777777777778D-01

  else if ( norder == 10 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.919533908166458813828932660822D+00
    xtab(3) =  - 0.738773865105505075003106174860D+00
    xtab(4) =  - 0.477924949810444495661175092731D+00
    xtab(5) =  - 0.165278957666387024626219765958D+00
    xtab(6) =    0.165278957666387024626219765958D+00
    xtab(7) =    0.477924949810444495661175092731D+00
    xtab(8) =    0.738773865105505075003106174860D+00
    xtab(9) =    0.919533908166458813828932660822D+00
    xtab(10) =   1.0D+00

    weight(1) =  0.222222222222222222222222222222D-01
    weight(2) =  0.133305990851070111126227170755D+00
    weight(3) =  0.224889342063126452119457821731D+00
    weight(4) =  0.292042683679683757875582257374D+00
    weight(5) =  0.327539761183897456656510527917D+00
    weight(6) =  0.327539761183897456656510527917D+00
    weight(7) =  0.292042683679683757875582257374D+00
    weight(8) =  0.224889342063126452119457821731D+00
    weight(9) =  0.133305990851070111126227170755D+00
    weight(10) = 0.222222222222222222222222222222D-01

  else if ( norder == 11 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.934001430408059134332274136099D+00
    xtab(3) =  - 0.784483473663144418622417816108D+00
    xtab(4) =  - 0.565235326996205006470963969478D+00
    xtab(5) =  - 0.295758135586939391431911515559D+00
    xtab(6) =    0.0D+00
    xtab(7) =    0.295758135586939391431911515559D+00
    xtab(8) =    0.565235326996205006470963969478D+00
    xtab(9) =    0.784483473663144418622417816108D+00
    xtab(10) =   0.934001430408059134332274136099D+00
    xtab(11) =   1.0D+00

    weight(1) =  0.181818181818181818181818181818D-01
    weight(2) =  0.109612273266994864461403449580D+00
    weight(3) =  0.187169881780305204108141521899D+00
    weight(4) =  0.248048104264028314040084866422D+00
    weight(5) =  0.286879124779008088679222403332D+00
    weight(6) =  0.300217595455690693785931881170D+00
    weight(7) =  0.286879124779008088679222403332D+00
    weight(8) =  0.248048104264028314040084866422D+00
    weight(9) =  0.187169881780305204108141521899D+00
    weight(10) = 0.109612273266994864461403449580D+00
    weight(11) = 0.181818181818181818181818181818D-01

  else if ( norder == 12 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.944899272222882223407580138303D+00
    xtab(3) =  - 0.819279321644006678348641581717D+00
    xtab(4) =  - 0.632876153031869677662404854444D+00
    xtab(5) =  - 0.399530940965348932264349791567D+00
    xtab(6) =  - 0.136552932854927554864061855740D+00
    xtab(7) =    0.136552932854927554864061855740D+00
    xtab(8) =    0.399530940965348932264349791567D+00
    xtab(9) =    0.632876153031869677662404854444D+00
    xtab(10) =   0.819279321644006678348641581717D+00
    xtab(11) =   0.944899272222882223407580138303D+00
    xtab(12) =   1.0D+00

    weight(1) =  0.151515151515151515151515151515D-01
    weight(2) =  0.916845174131961306683425941341D-01
    weight(3) =  0.157974705564370115164671062700D+00
    weight(4) =  0.212508417761021145358302077367D+00
    weight(5) =  0.251275603199201280293244412148D+00
    weight(6) =  0.271405240910696177000288338500D+00
    weight(7) =  0.271405240910696177000288338500D+00
    weight(8) =  0.251275603199201280293244412148D+00
    weight(9) =  0.212508417761021145358302077367D+00
    weight(10) = 0.157974705564370115164671062700D+00
    weight(11) = 0.916845174131961306683425941341D-01
    weight(12) = 0.151515151515151515151515151515D-01

  else if ( norder == 13 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.953309846642163911896905464755D+00
    xtab(3) =  - 0.846347564651872316865925607099D+00
    xtab(4) =  - 0.686188469081757426072759039566D+00
    xtab(5) =  - 0.482909821091336201746937233637D+00
    xtab(6) =  - 0.249286930106239992568673700374D+00
    xtab(7) =    0.0D+00
    xtab(8) =    0.249286930106239992568673700374D+00
    xtab(9) =    0.482909821091336201746937233637D+00
    xtab(10) =   0.686188469081757426072759039566D+00
    xtab(11) =   0.846347564651872316865925607099D+00
    xtab(12) =   0.953309846642163911896905464755D+00
    xtab(13) =   1.0D+00

    weight(1) =  0.128205128205128205128205128205D-01
    weight(2) =  0.778016867468189277935889883331D-01
    weight(3) =  0.134981926689608349119914762589D+00
    weight(4) =  0.183646865203550092007494258747D+00
    weight(5) =  0.220767793566110086085534008379D+00
    weight(6) =  0.244015790306676356458578148360D+00
    weight(7) =  0.251930849333446736044138641541D+00
    weight(8) =  0.244015790306676356458578148360D+00
    weight(9) =  0.220767793566110086085534008379D+00
    weight(10) = 0.183646865203550092007494258747D+00
    weight(11) = 0.134981926689608349119914762589D+00
    weight(12) = 0.778016867468189277935889883331D-01
    weight(13) = 0.128205128205128205128205128205D-01

  else if ( norder == 14 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.959935045267260901355100162015D+00
    xtab(3) =  - 0.867801053830347251000220202908D+00
    xtab(4) =  - 0.728868599091326140584672400521D+00
    xtab(5) =  - 0.550639402928647055316622705859D+00
    xtab(6) =  - 0.342724013342712845043903403642D+00
    xtab(7) =  - 0.116331868883703867658776709736D+00
    xtab(8) =    0.116331868883703867658776709736D+00
    xtab(9) =    0.342724013342712845043903403642D+00
    xtab(10) =   0.550639402928647055316622705859D+00
    xtab(11) =   0.728868599091326140584672400521D+00
    xtab(12) =   0.867801053830347251000220202908D+00
    xtab(13) =   0.959935045267260901355100162015D+00
    xtab(14) =   1.0D+00

    weight(1) =  0.109890109890109890109890109890D-01
    weight(2) =  0.668372844976812846340706607461D-01
    weight(3) =  0.116586655898711651540996670655D+00
    weight(4) =  0.160021851762952142412820997988D+00
    weight(5) =  0.194826149373416118640331778376D+00
    weight(6) =  0.219126253009770754871162523954D+00
    weight(7) =  0.231612794468457058889628357293D+00
    weight(8) =  0.231612794468457058889628357293D+00
    weight(9) =  0.219126253009770754871162523954D+00
    weight(10) = 0.194826149373416118640331778376D+00
    weight(11) = 0.160021851762952142412820997988D+00
    weight(12) = 0.116586655898711651540996670655D+00
    weight(13) = 0.668372844976812846340706607461D-01
    weight(14) = 0.109890109890109890109890109890D-01

  else if ( norder == 15 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.965245926503838572795851392070D+00
    xtab(3) =  - 0.885082044222976298825401631482D+00
    xtab(4) =  - 0.763519689951815200704118475976D+00
    xtab(5) =  - 0.606253205469845711123529938637D+00
    xtab(6) =  - 0.420638054713672480921896938739D+00
    xtab(7) =  - 0.215353955363794238225679446273D+00
    xtab(8) =    0.0D+00
    xtab(9) =    0.215353955363794238225679446273D+00
    xtab(10) =   0.420638054713672480921896938739D+00
    xtab(11) =   0.606253205469845711123529938637D+00
    xtab(12) =   0.763519689951815200704118475976D+00
    xtab(13) =   0.885082044222976298825401631482D+00
    xtab(14) =   0.965245926503838572795851392070D+00
    xtab(15) =   1.0D+00

    weight(1) =  0.952380952380952380952380952381D-02
    weight(2) =  0.580298930286012490968805840253D-01
    weight(3) =  0.101660070325718067603666170789D+00
    weight(4) =  0.140511699802428109460446805644D+00
    weight(5) =  0.172789647253600949052077099408D+00
    weight(6) =  0.196987235964613356092500346507D+00
    weight(7) =  0.211973585926820920127430076977D+00
    weight(8) =  0.217048116348815649514950214251D+00
    weight(9) =  0.211973585926820920127430076977D+00
    weight(10) = 0.196987235964613356092500346507D+00
    weight(11) = 0.172789647253600949052077099408D+00
    weight(12) = 0.140511699802428109460446805644D+00
    weight(13) = 0.101660070325718067603666170789D+00
    weight(14) = 0.580298930286012490968805840253D-01
    weight(15) = 0.952380952380952380952380952381D-02

  else if ( norder == 16 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.969568046270217932952242738367D+00
    xtab(3) =  - 0.899200533093472092994628261520D+00
    xtab(4) =  - 0.792008291861815063931088270963D+00
    xtab(5) =  - 0.652388702882493089467883219641D+00
    xtab(6) =  - 0.486059421887137611781890785847D+00
    xtab(7) =  - 0.299830468900763208098353454722D+00
    xtab(8) =  - 0.101326273521949447843033005046D+00
    xtab(9) =    0.101326273521949447843033005046D+00
    xtab(10) =   0.299830468900763208098353454722D+00
    xtab(11) =   0.486059421887137611781890785847D+00
    xtab(12) =   0.652388702882493089467883219641D+00
    xtab(13) =   0.792008291861815063931088270963D+00
    xtab(14) =   0.899200533093472092994628261520D+00
    xtab(15) =   0.969568046270217932952242738367D+00
    xtab(16) =   1.0D+00

    weight(1) =  0.833333333333333333333333333333D-02
    weight(2) =  0.508503610059199054032449195655D-01
    weight(3) =  0.893936973259308009910520801661D-01
    weight(4) =  0.124255382132514098349536332657D+00
    weight(5) =  0.154026980807164280815644940485D+00
    weight(6) =  0.177491913391704125301075669528D+00
    weight(7) =  0.193690023825203584316913598854D+00
    weight(8) =  0.201958308178229871489199125411D+00
    weight(9) =  0.201958308178229871489199125411D+00
    weight(10) = 0.193690023825203584316913598854D+00
    weight(11) = 0.177491913391704125301075669528D+00
    weight(12) = 0.154026980807164280815644940485D+00
    weight(13) = 0.124255382132514098349536332657D+00
    weight(14) = 0.893936973259308009910520801661D-01
    weight(15) = 0.508503610059199054032449195655D-01
    weight(16) = 0.833333333333333333333333333333D-02

  else if ( norder == 17 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.973132176631418314156979501874D+00
    xtab(3) =  - 0.910879995915573595623802506398D+00
    xtab(4) =  - 0.815696251221770307106750553238D+00
    xtab(5) =  - 0.691028980627684705394919357372D+00
    xtab(6) =  - 0.541385399330101539123733407504D+00
    xtab(7) =  - 0.372174433565477041907234680735D+00
    xtab(8) =  - 0.189511973518317388304263014753D+00
    xtab(9) =    0.0D+00
    xtab(10) =   0.189511973518317388304263014753D+00
    xtab(11) =   0.372174433565477041907234680735D+00
    xtab(12) =   0.541385399330101539123733407504D+00
    xtab(13) =   0.691028980627684705394919357372D+00
    xtab(14) =   0.815696251221770307106750553238D+00
    xtab(15) =   0.910879995915573595623802506398D+00
    xtab(16) =   0.973132176631418314156979501874D+00
    xtab(17) =   1.0D+00

    weight(1) =  0.735294117647058823529411764706D-02
    weight(2) =  0.449219405432542096474009546232D-01
    weight(3) =  0.791982705036871191902644299528D-01
    weight(4) =  0.110592909007028161375772705220D+00
    weight(5) =  0.137987746201926559056201574954D+00
    weight(6) =  0.160394661997621539516328365865D+00
    weight(7) =  0.177004253515657870436945745363D+00
    weight(8) =  0.187216339677619235892088482861D+00
    weight(9) =  0.190661874753469433299407247028D+00
    weight(10) = 0.187216339677619235892088482861D+00
    weight(11) = 0.177004253515657870436945745363D+00
    weight(12) = 0.160394661997621539516328365865D+00
    weight(13) = 0.137987746201926559056201574954D+00
    weight(14) = 0.110592909007028161375772705220D+00
    weight(15) = 0.791982705036871191902644299528D-01
    weight(16) = 0.449219405432542096474009546232D-01
    weight(17) = 0.735294117647058823529411764706D-02

  else if ( norder == 18 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.976105557412198542864518924342D+00
    xtab(3) =  - 0.920649185347533873837854625431D+00
    xtab(4) =  - 0.835593535218090213713646362328D+00
    xtab(5) =  - 0.723679329283242681306210365302D+00
    xtab(6) =  - 0.588504834318661761173535893194D+00
    xtab(7) =  - 0.434415036912123975342287136741D+00
    xtab(8) =  - 0.266362652878280984167665332026D+00
    xtab(9) =  - 0.897490934846521110226450100886D-01
    xtab(10) =   0.897490934846521110226450100886D-01
    xtab(11) =   0.266362652878280984167665332026D+00
    xtab(12) =   0.434415036912123975342287136741D+00
    xtab(13) =   0.588504834318661761173535893194D+00
    xtab(14) =   0.723679329283242681306210365302D+00
    xtab(15) =   0.835593535218090213713646362328D+00
    xtab(16) =   0.920649185347533873837854625431D+00
    xtab(17) =   0.976105557412198542864518924342D+00
    xtab(18) =   1.0D+00

    weight(1) =  0.653594771241830065359477124183D-02
    weight(2) =  0.399706288109140661375991764101D-01
    weight(3) =  0.706371668856336649992229601678D-01
    weight(4) =  0.990162717175028023944236053187D-01
    weight(5) =  0.124210533132967100263396358897D+00
    weight(6) =  0.145411961573802267983003210494D+00
    weight(7) =  0.161939517237602489264326706700D+00
    weight(8) =  0.173262109489456226010614403827D+00
    weight(9) =  0.179015863439703082293818806944D+00
    weight(10) = 0.179015863439703082293818806944D+00
    weight(11) = 0.173262109489456226010614403827D+00
    weight(12) = 0.161939517237602489264326706700D+00
    weight(13) = 0.145411961573802267983003210494D+00
    weight(14) = 0.124210533132967100263396358897D+00
    weight(15) = 0.990162717175028023944236053187D-01
    weight(16) = 0.706371668856336649992229601678D-01
    weight(17) = 0.399706288109140661375991764101D-01
    weight(18) = 0.653594771241830065359477124183D-02

  else if ( norder == 19 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.978611766222080095152634063110D+00
    xtab(3) =  - 0.928901528152586243717940258797D+00
    xtab(4) =  - 0.852460577796646093085955970041D+00
    xtab(5) =  - 0.751494202552613014163637489634D+00
    xtab(6) =  - 0.628908137265220497766832306229D+00
    xtab(7) =  - 0.488229285680713502777909637625D+00
    xtab(8) =  - 0.333504847824498610298500103845D+00
    xtab(9) =  - 0.169186023409281571375154153445D+00
    xtab(10) =   0.0D+00
    xtab(11) =   0.169186023409281571375154153445D+00
    xtab(12) =   0.333504847824498610298500103845D+00
    xtab(13) =   0.488229285680713502777909637625D+00
    xtab(14) =   0.628908137265220497766832306229D+00
    xtab(15) =   0.751494202552613014163637489634D+00
    xtab(16) =   0.852460577796646093085955970041D+00
    xtab(17) =   0.928901528152586243717940258797D+00
    xtab(18) =   0.978611766222080095152634063110D+00
    xtab(19) =   1.0D+00

    weight(1) =  0.584795321637426900584795321637D-02
    weight(2) =  0.357933651861764771154255690351D-01
    weight(3) =  0.633818917626297368516956904183D-01
    weight(4) =  0.891317570992070844480087905562D-01
    weight(5) =  0.112315341477305044070910015464D+00
    weight(6) =  0.132267280448750776926046733910D+00
    weight(7) =  0.148413942595938885009680643668D+00
    weight(8) =  0.160290924044061241979910968184D+00
    weight(9) =  0.167556584527142867270137277740D+00
    weight(10) = 0.170001919284827234644672715617D+00
    weight(11) = 0.167556584527142867270137277740D+00
    weight(12) = 0.160290924044061241979910968184D+00
    weight(13) = 0.148413942595938885009680643668D+00
    weight(14) = 0.132267280448750776926046733910D+00
    weight(15) = 0.112315341477305044070910015464D+00
    weight(16) = 0.891317570992070844480087905562D-01
    weight(17) = 0.633818917626297368516956904183D-01
    weight(18) = 0.357933651861764771154255690351D-01
    weight(19) = 0.584795321637426900584795321637D-02

  else if ( norder == 20 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.980743704893914171925446438584D+00
    xtab(3) =  - 0.935934498812665435716181584931D+00
    xtab(4) =  - 0.866877978089950141309847214616D+00
    xtab(5) =  - 0.775368260952055870414317527595D+00
    xtab(6) =  - 0.663776402290311289846403322971D+00
    xtab(7) =  - 0.534992864031886261648135961829D+00
    xtab(8) =  - 0.392353183713909299386474703816D+00
    xtab(9) =  - 0.239551705922986495182401356927D+00
    xtab(10) = - 0.805459372388218379759445181596D-01
    xtab(11) =   0.805459372388218379759445181596D-01
    xtab(12) =   0.239551705922986495182401356927D+00
    xtab(13) =   0.392353183713909299386474703816D+00
    xtab(14) =   0.534992864031886261648135961829D+00
    xtab(15) =   0.663776402290311289846403322971D+00
    xtab(16) =   0.775368260952055870414317527595D+00
    xtab(17) =   0.866877978089950141309847214616D+00
    xtab(18) =   0.935934498812665435716181584931D+00
    xtab(19) =   0.980743704893914171925446438584D+00
    xtab(20) =   1.0D+00

    weight(1) =  0.526315789473684210526315789474D-02
    weight(2) =  0.322371231884889414916050281173D-01
    weight(3) =  0.571818021275668260047536271732D-01
    weight(4) =  0.806317639961196031447768461137D-01
    weight(5) =  0.101991499699450815683781205733D+00
    weight(6) =  0.120709227628674725099429705002D+00
    weight(7) =  0.136300482358724184489780792989D+00
    weight(8) =  0.148361554070916825814713013734D+00
    weight(9) =  0.156580102647475487158169896794D+00
    weight(10) = 0.160743286387845749007726726449D+00
    weight(11) = 0.160743286387845749007726726449D+00
    weight(12) = 0.156580102647475487158169896794D+00
    weight(13) = 0.148361554070916825814713013734D+00
    weight(14) = 0.136300482358724184489780792989D+00
    weight(15) = 0.120709227628674725099429705002D+00
    weight(16) = 0.101991499699450815683781205733D+00
    weight(17) = 0.806317639961196031447768461137D-01
    weight(18) = 0.571818021275668260047536271732D-01
    weight(19) = 0.322371231884889414916050281173D-01
    weight(20) = 0.526315789473684210526315789474D-02

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOBATTO_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are between 2 and 20.'
    stop

  end if

  return
end
function log_gamma ( x )
!
!*******************************************************************************
!
!! LOG_GAMMA calculates the natural logarithm of GAMMA(X).
!
!
!  Discussion:
!
!    The method uses Stirling's approximation, and is accurate to about
!    12 decimal places.
!
!  Modified:
!
!    19 September 1998
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.  The routine
!    will fail if GAMMA(X) is not positive.  X should be greater than 0.
!
!    Output, real ( kind = 8 ) LOG_GAMMA, the natural logarithm of the
!    gamma function of X.
!
  implicit none
!
  real ( kind = 8 ) d_pi
  integer i
  integer k
  real ( kind = 8 ) log_gamma
  integer m
  real ( kind = 8 ) p
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
  if ( x < 0.5 ) then

    m = 1
    x2 = 1.0D+00 - x

  else

    m = 0
    x2 = x

  end if

  k = - 1

  do

    k = k + 1

    if ( 6.0D+00 < x2 + real ( k, kind = 8 ) ) then
      exit
    end if

  end do

  z = x2 + real ( k, kind = 8 )

  y = ( z - 0.5 ) * log ( z ) - z + 0.9189385332047D+00 + &
       ( ( ( ( ( &
       - 4146.0D+00 / z**2 &
       + 1820.0D+00 ) / z**2 &
       - 1287.0D+00 ) / z**2 &
       + 1716.0D+00 ) / z**2 &
       - 6006.0D+00 ) / z**2 &
       + 180180.0D+00 ) / z / 2162160.0D+00

  if ( 0 < k ) then

    do i = 1, k
      y = y - log ( x2 + real ( k - i, kind = 8 ) )
    end do

  end if

  if ( m /= 0 ) then

    p = d_pi ( ) / sin ( d_pi ( ) * ( 1.0D+00 - x2 ) )

    if ( p <= 0.0D+00 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOG_GAMMA - fatal error!'
      stop

    else

      y = log ( p ) - y

    end if

  end if

  log_gamma = y

  return
end
subroutine moulton_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! MOULTON_SET sets weights for Adams-Moulton quadrature.
!
!
!  Definition:
!
!    Adams-Moulton quadrature formulas are normally used in solving
!    ordinary differential equations, and are not suitable for general
!    quadrature computations.  However, an Adams-Moulton formula is
!    equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an implicit formula that relies on known values
!    of F(Y(X)) at X(M-N+1) through X(M), plus the unknown value at X(M+1).
!
!    Suppose the unknown function is denoted by Y(X), with derivative F(Y(X)),
!    and that approximate values of the function are known at a series of
!    X values, which we write as X(1), X(2), ..., X(M).  We write the value
!    Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y' = F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
!             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+2-I)) approximately.
!
!    Note that this formula is implicit, since the unknown value Y(M+1)
!    appears on the right hand side.  Hence, in ODE applications, this
!    equation must be solved via a nonlinear equation solver.  For
!    quadrature problems, where the function to be integrated is known
!    beforehand, this is not a problem, and the calculation is explicit.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!  Integration interval:
!
!    [ 0, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( 0 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( 2 - I )
!
!  Note:
!
!    The Adams-Moulton formulas require equally spaced data.
!
!    Here is how the formula is applied in the case with non-unit spacing:
!
!      Integral ( A <= X <= A+H ) F(X) dX =
!      H * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( A - (I-2)*H ),
!      approximately.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    page 915 ("Lagrangian Integration Coefficients").
!
!    Jean Lapidus and John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.  NORDER must be
!    between 1 and 10.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    WEIGHT(1) is the weight at X = 1, WEIGHT(2) the weight at X = 0, and so on.
!    The weights are rational.  The weights are not symmetric, and
!    some weights may be negative.  They should sum to 1.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 1 ) then

    weight(1) =  1.0D+00

  else if ( norder == 2 ) then

    d = 2.0D+00

    weight(1) =  1.0D+00 / d
    weight(2) =  1.0D+00 / d

  else if ( norder == 3 ) then

    d = 12.0D+00

    weight(1) =    5.0D+00 / d
    weight(2) =    8.0D+00 / d
    weight(3) =  - 1.0D+00 / d

  else if ( norder == 4 ) then

    d = 24.0D+00

    weight(1) =    9.0D+00 / d
    weight(2) =   19.0D+00 / d
    weight(3) =  - 5.0D+00 / d
    weight(4) =    1.0D+00 / d

  else if ( norder == 5 ) then

    d = 720.0D+00

    weight(1) =    251.0D+00 / d
    weight(2) =    646.0D+00 / d
    weight(3) =  - 264.0D+00 / d
    weight(4) =    106.0D+00 / d
    weight(5) =   - 19.0D+00 / d

  else if ( norder == 6 ) then

    d = 1440.0D+00

    weight(1) =    475.0D+00 / d
    weight(2) =   1427.0D+00 / d
    weight(3) =  - 798.0D+00 / d
    weight(4) =    482.0D+00 / d
    weight(5) =  - 173.0D+00 / d
    weight(6) =     27.0D+00 / d

  else if ( norder == 7 ) then

    d = 60480.0D+00

    weight(1) =    19087.0D+00 / d
    weight(2) =    65112.0D+00 / d
    weight(3) =  - 46461.0D+00 / d
    weight(4) =    37504.0D+00 / d
    weight(5) =  - 20211.0D+00 / d
    weight(6) =     6312.0D+00 / d
    weight(7) =    - 863.0D+00 / d

  else if ( norder == 8 ) then

    d = 120960.0D+00

    weight(1) =    36799.0D+00 / d
    weight(2) =   139849.0D+00 / d
    weight(3) = - 121797.0D+00 / d
    weight(4) =   123133.0D+00 / d
    weight(5) =  - 88547.0D+00 / d
    weight(6) =    41499.0D+00 / d
    weight(7) =  - 11351.0D+00 / d
    weight(8) =     1375.0D+00 / d

  else if ( norder == 9 ) then

    d = 3628800.0D+00

    weight(1) =   1070017.0D+00 / d
    weight(2) =   4467094.0D+00 / d
    weight(3) = - 4604594.0D+00 / d
    weight(4) =   5595358.0D+00 / d
    weight(5) = - 5033120.0D+00 / d
    weight(6) =   3146338.0D+00 / d
    weight(7) = - 1291214.0D+00 / d
    weight(8) =    312874.0D+00 / d
    weight(9) =   - 33953.0D+00 / d

  else if ( norder == 10 ) then

    d = 7257600.0D+00

    weight(1) =    2082753.0D+00 / d
    weight(2) =    9449717.0D+00 / d
    weight(3) = - 11271304.0D+00 / d
    weight(4) =   16002320.0D+00 / d
    weight(5) = - 17283646.0D+00 / d
    weight(6) =   13510082.0D+00 / d
    weight(7) =  - 7394032.0D+00 / d
    weight(8) =    2687864.0D+00 / d
    weight(9) =   - 583435.0D+00 / d
    weight(10) =     57281.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOULTON_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 through 10.'
    stop

  end if

  do i = 1, norder
    xtab(i) = real ( 2 - i, kind = 8 )
  end do

  return
end
subroutine nc_com ( norder, a, b, xtab, weight )
!
!*******************************************************************************
!
!! NC_COM computes a Newton-Cotes quadrature rule with given nodes.
!
!
!  Definition:
!
!    For the interval [A,B], the Newton-Cotes quadrature rule estimates
!
!      Integral ( A <= X <= B ) F(X) dX
!
!    using NORDER abscissas XTAB(I) and a weight vector WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the equally spaced abscissas include A and B.
!    For the OPEN rule, the equally spaced abscissas do not include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints of the interval
!    over which the quadrature rule is to be applied.
!
!    Input, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) diftab(norder)
  integer i
  integer j
  integer k
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) yvala
  real ( kind = 8 ) yvalb
!
  do i = 1, norder
!
!  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
!  and zero at the other nodes.
!
    diftab(1:norder) = 0.0D+00
    diftab(i) = 1.0D+00

    do j = 2, norder
      do k = j, norder
        diftab(norder+j-k) = ( diftab(norder+j-k-1) - diftab(norder+j-k) ) &
          / ( xtab(norder+1-k) - xtab(norder+j-k) )
      end do
    end do

    do j = 1, norder-1
      do k = 1, norder-j
        diftab(norder-k) = diftab(norder-k) - xtab(norder-k-j+1) * &
          diftab(norder-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of the polynomial at the left and
!  right endpoints.
!
    yvala = diftab(norder) / real ( norder, kind = 8 )
    do j = norder-1, 1, -1
      yvala = yvala * a + diftab(j) / real ( j, kind = 8 )
    end do
    yvala = yvala * a

    yvalb = diftab(norder) / real ( norder, kind = 8 )
    do j = norder-1, 1, -1
      yvalb = yvalb * b + diftab(j) / real ( j, kind = 8 )
    end do
    yvalb = yvalb * b

    weight(i) = yvalb - yvala

  end do

  return
end
subroutine ncc_com ( norder, xtab, weight )
!
!*******************************************************************************
!
!! NCC_COM computes a Newton-Cotes closed quadrature rule of given order.
!
!
!  Definition:
!
!    For the interval [-1,1], the Newton-Cotes open quadrature rule
!    estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the abscissas include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which should be
!    at least 2.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
!  Compute a closed quadrature rule.
!
  a = -1.0D+00
  b =  1.0D+00

  if ( norder == 1 ) then

    xtab(1) = 0.0D+00
    weight(1) = 2.0D+00
    return

  else

    do i = 1, norder
      xtab(i) = ( real ( norder - i,     kind = 8 ) * a   &
                + real (          i - 1, kind = 8 ) * b ) &
                / real ( norder     - 1, kind = 8 )
    end do

  end if

  call nc_com ( norder, a, b, xtab, weight )

  return
end
subroutine ncc_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! NCC_SET sets abscissas and weights for closed Newton-Cotes quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Note:
!
!    The closed Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with tabulated function data.
!
!    The rules are called "closed" because they include the endpoints.
!
!    The higher order rules involve negative weights.  These can produce
!    loss of accuracy due to the subtraction of large, nearly equal quantities.
!
!    NORDER = 2 is the trapezoidal rule.
!    NORDER = 3 is Simpson's rule.
!    NORDER = 4 is Simpson's 3/8 rule.
!    NORDER = 5 is Bode's rule.
!
!    The Kopal reference for NORDER = 12 lists
!      WEIGHT(6) = 15494566.0D+00 / 43545600.0D+00
!    but this results in a set of coeffients that don't add up to 2.
!    The correct value is
!      WEIGHT(6) = 15493566.0D+00 / 43545600.0.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Johnson,
!    Quarterly Journal of Mathematics,
!    Volume 46, Number 52, 1915.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 2 and 20.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!    The abscissas are uniformly spaced in the interval, and include
!    -1 and 1.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are symmetric, rational, and should sum to 2.0.
!    Some weights may be negative.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 2 ) then

    weight(1) = 1.0D+00
    weight(2) = 1.0D+00

  else if ( norder == 3 ) then

    d = 3.0D+00

    weight(1) = 1.0D+00 / d
    weight(2) = 4.0D+00 / d
    weight(3) = 1.0D+00 / d

  else if ( norder == 4 ) then

    d = 4.0D+00

    weight(1) = 1.0D+00 / d
    weight(2) = 3.0D+00 / d
    weight(3) = 3.0D+00 / d
    weight(4) = 1.0D+00 / d

  else if ( norder == 5 ) then

    d = 45.0D+00

    weight(1) =  7.0D+00 / d
    weight(2) = 32.0D+00 / d
    weight(3) = 12.0D+00 / d
    weight(4) = 32.0D+00 / d
    weight(5) =  7.0D+00 / d

  else if ( norder == 6 ) then

    d = 144.0D+00

    weight(1) = 19.0D+00 / d
    weight(2) = 75.0D+00 / d
    weight(3) = 50.0D+00 / d
    weight(4) = 50.0D+00 / d
    weight(5) = 75.0D+00 / d
    weight(6) = 19.0D+00 / d

  else if ( norder == 7 ) then

    d = 420.0D+00

    weight(1) =  41.0D+00 / d
    weight(2) = 216.0D+00 / d
    weight(3) =  27.0D+00 / d
    weight(4) = 272.0D+00 / d
    weight(5) =  27.0D+00 / d
    weight(6) = 216.0D+00 / d
    weight(7) =  41.0D+00 / d

  else if ( norder == 8 ) then

    d = 8640.0D+00

    weight(1) =  751.0D+00 / d
    weight(2) = 3577.0D+00 / d
    weight(3) = 1323.0D+00 / d
    weight(4) = 2989.0D+00 / d
    weight(5) = 2989.0D+00 / d
    weight(6) = 1323.0D+00 / d
    weight(7) = 3577.0D+00 / d
    weight(8) =  751.0D+00 / d

  else if ( norder == 9 ) then

    d = 14175.0D+00

    weight(1) =    989.0D+00 / d
    weight(2) =   5888.0D+00 / d
    weight(3) =  - 928.0D+00 / d
    weight(4) =  10496.0D+00 / d
    weight(5) = - 4540.0D+00 / d
    weight(6) =  10496.0D+00 / d
    weight(7) =  - 928.0D+00 / d
    weight(8) =   5888.0D+00 / d
    weight(9) =    989.0D+00 / d

  else if ( norder == 10 ) then

    d = 44800.0D+00

    weight(1) =   2857.0D+00 / d
    weight(2) =  15741.0D+00 / d
    weight(3) =   1080.0D+00 / d
    weight(4) =  19344.0D+00 / d
    weight(5) =   5778.0D+00 / d
    weight(6) =   5778.0D+00 / d
    weight(7) =  19344.0D+00 / d
    weight(8) =   1080.0D+00 / d
    weight(9) =  15741.0D+00 / d
    weight(10) =  2857.0D+00 / d

  else if ( norder == 11 ) then

    d = 299376.0D+00

    weight(1) =     16067.0D+00 / d
    weight(2) =    106300.0D+00 / d
    weight(3) =   - 48525.0D+00 / d
    weight(4) =    272400.0D+00 / d
    weight(5) =  - 260550.0D+00 / d
    weight(6) =    427368.0D+00 / d
    weight(7) =  - 260550.0D+00 / d
    weight(8) =    272400.0D+00 / d
    weight(9) =   - 48525.0D+00 / d
    weight(10) =   106300.0D+00 / d
    weight(11) =    16067.0D+00 / d

  else if ( norder == 12 ) then

    d = 43545600.0D+00

    weight(1) =     2171465.0D+00 / d
    weight(2) =    13486539.0D+00 / d
    weight(3) =   - 3237113.0D+00 / d
    weight(4) =    25226685.0D+00 / d
    weight(5) =   - 9595542.0D+00 / d
    weight(6) =    15493566.0D+00 / d
    weight(7) =    15493566.0D+00 / d
    weight(8) =   - 9595542.0D+00 / d
    weight(9) =    25226685.0D+00 / d
    weight(10) =  - 3237113.0D+00 / d
    weight(11) =   13486539.0D+00 / d
    weight(12) =    2171465.0D+00 / d

  else if ( norder == 13 ) then

    d = 31531500.0D+00

    weight(1) =      1364651.0D+00 / d
    weight(2) =      9903168.0D+00 / d
    weight(3) =    - 7587864.0D+00 / d
    weight(4) =     35725120.0D+00 / d
    weight(5) =   - 51491295.0D+00 / d
    weight(6) =     87516288.0D+00 / d
    weight(7) =   - 87797136.0D+00 / d
    weight(8) =     87516288.0D+00 / d
    weight(9) =   - 51491295.0D+00 / d
    weight(10) =    35725120.0D+00 / d
    weight(11) =   - 7587864.0D+00 / d
    weight(12) =     9903168.0D+00 / d
    weight(13) =     1364651.0D+00 / d

  else if ( norder == 14 ) then

    d = 150885504000.0D+00

    weight(1) =      6137698213.0D+00 / d
    weight(2) =     42194238652.0D+00 / d
    weight(3) =   - 23361540993.0D+00 / d
    weight(4) =    116778274403.0D+00 / d
    weight(5) =  - 113219777650.0D+00 / d
    weight(6) =    154424590209.0D+00 / d
    weight(7) =   - 32067978834.0D+00 / d
    weight(8) =   - 32067978834.0D+00 / d
    weight(9) =    154424590209.0D+00 / d
    weight(10) = - 113219777650.0D+00 / d
    weight(11) =   116778274403.0D+00 / d
    weight(12) =  - 23361540993.0D+00 / d
    weight(13) =    42194238652.0D+00 / d
    weight(14) =     6137698213.0D+00 / d

  else if ( norder == 15 ) then

    d = 2501928000.0D+00

    weight(1) =       90241897.0D+00 / d
    weight(2) =      710986864.0D+00 / d
    weight(3) =    - 770720657.0D+00 / d
    weight(4) =     3501442784.0D+00 / d
    weight(5) =   - 6625093363.0D+00 / d
    weight(6) =    12630121616.0D+00 / d
    weight(7) =  - 16802270373.0D+00 / d
    weight(8) =    19534438464.0D+00 / d
    weight(9) =  - 16802270373.0D+00 / d
    weight(10) =   12630121616.0D+00 / d
    weight(11) =  - 6625093363.0D+00 / d
    weight(12) =    3501442784.0D+00 / d
    weight(13) =   - 770720657.0D+00 / d
    weight(14) =     710986864.0D+00 / d
    weight(15) =      90241897.0D+00 / d

  else if ( norder == 16 ) then

    d = 3099672576.0D+00

    weight(1) =     105930069.0D+00 / d
    weight(2) =     796661595.0D+00 / d
    weight(3) =   - 698808195.0D+00 / d
    weight(4) =    3143332755.0D+00 / d
    weight(5) =  - 4688522055.0D+00 / d
    weight(6) =    7385654007.0D+00 / d
    weight(7) =  - 6000998415.0D+00 / d
    weight(8) =    3056422815.0D+00 / d
    weight(9) =    3056422815.0D+00 / d
    weight(10) = - 6000998415.0D+00 / d
    weight(11) =   7385654007.0D+00 / d
    weight(12) = - 4688522055.0D+00 / d
    weight(13) =   3143332755.0D+00 / d
    weight(14) =  - 698808195.0D+00 / d
    weight(15) =    796661595.0D+00 / d
    weight(16) =    105930069.0D+00 / d

  else if ( norder == 17 ) then

    d = 488462349375.0D+00

    weight(1) =       15043611773.0D+00 / d
    weight(2) =      127626606592.0D+00 / d
    weight(3) =    - 179731134720.0D+00 / d
    weight(4) =      832211855360.0D+00 / d
    weight(5) =   - 1929498607520.0D+00 / d
    weight(6) =     4177588893696.0D+00 / d
    weight(7) =   - 6806534407936.0D+00 / d
    weight(8) =     9368875018240.0D+00 / d
    weight(9) =  - 10234238972220.0D+00 / d
    weight(10) =    9368875018240.0D+00 / d
    weight(11) =  - 6806534407936.0D+00 / d
    weight(12) =    4177588893696.0D+00 / d
    weight(13) =  - 1929498607520.0D+00 / d
    weight(14) =     832211855360.0D+00 / d
    weight(15) =   - 179731134720.0D+00 / d
    weight(16) =     127626606592.0D+00 / d
    weight(17) =      15043611773.0D+00 / d

  else if ( norder == 18 ) then

    d = 1883051089920000.0D+00

    weight(1) =       55294720874657.0D+00 / d
    weight(2) =      450185515446285.0D+00 / d
    weight(3) =    - 542023437008852.0D+00 / d
    weight(4) =     2428636525764260.0D+00 / d
    weight(5) =   - 4768916800123440.0D+00 / d
    weight(6) =     8855416648684984.0D+00 / d
    weight(7) =  - 10905371859796660.0D+00 / d
    weight(8) =    10069615750132836.0D+00 / d
    weight(9) =   - 3759785974054070.0D+00 / d
    weight(10) =  - 3759785974054070.0D+00 / d
    weight(11) =   10069615750132836.0D+00 / d
    weight(12) = - 10905371859796660.0D+00 / d
    weight(13) =    8855416648684984.0D+00 / d
    weight(14) =  - 4768916800123440.0D+00 / d
    weight(15) =    2428636525764260.0D+00 / d
    weight(16) =   - 542023437008852.0D+00 / d
    weight(17) =     450185515446285.0D+00 / d
    weight(18) =      55294720874657.0D+00 / d

  else if ( norder == 19 ) then

    d = 7604556960000.0D+00

    weight(1) =       203732352169.0D+00 / d
    weight(2) =      1848730221900.0D+00 / d
    weight(3) =    - 3212744374395.0D+00 / d
    weight(4) =     15529830312096.0D+00 / d
    weight(5) =   - 42368630685840.0D+00 / d
    weight(6) =    103680563465808.0D+00 / d
    weight(7) =  - 198648429867720.0D+00 / d
    weight(8) =    319035784479840.0D+00 / d
    weight(9) =  - 419127951114198.0D+00 / d
    weight(10) =   461327344340680.0D+00 / d
    weight(11) = - 419127951114198.0D+00 / d
    weight(12) =   319035784479840.0D+00 / d
    weight(13) = - 198648429867720.0D+00 / d
    weight(14) =   103680563465808.0D+00 / d
    weight(15) =  - 42368630685840.0D+00 / d
    weight(16) =    15529830312096.0D+00 / d
    weight(17) =   - 3212744374395.0D+00 / d
    weight(18) =     1848730221900.0D+00 / d
    weight(19) =      203732352169.0D+00 / d

  else if ( norder == 20 ) then

    d = 2688996956405760000.0D+00

    weight(1) =       69028763155644023.0D+00 / d
    weight(2) =      603652082270808125.0D+00 / d
    weight(3) =    - 926840515700222955.0D+00 / d
    weight(4) =     4301581538450500095.0D+00 / d
    weight(5) =  - 10343692234243192788.0D+00 / d
    weight(6) =    22336420328479961316.0D+00 / d
    weight(7) =  - 35331888421114781580.0D+00 / d
    weight(8) =    43920768370565135580.0D+00 / d
    weight(9) =  - 37088370261379851390.0D+00 / d
    weight(10) =   15148337305921759574.0D+00 / d
    weight(11) =   15148337305921759574.0D+00 / d
    weight(12) = - 37088370261379851390.0D+00 / d
    weight(13) =   43920768370565135580.0D+00 / d
    weight(14) = - 35331888421114781580.0D+00 / d
    weight(15) =   22336420328479961316.0D+00 / d
    weight(16) = - 10343692234243192788.0D+00 / d
    weight(17) =    4301581538450500095.0D+00 / d
    weight(18) =   - 926840515700222955.0D+00 / d
    weight(19) =     603652082270808125.0D+00 / d
    weight(20) =      69028763155644023.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 2 through 20.'
    stop

  end if
!
!  The abscissas are uniformly spaced.
!
  do i = 1, norder
    xtab(i) = real ( 2 * i - 1 - norder, kind = 8 ) &
            / real ( norder - 1, kind = 8 )
  end do

  return
end
subroutine nco_com ( norder, xtab, weight )
!
!*******************************************************************************
!
!! NCO_COM computes a Newton-Cotes open quadrature rule of given order.
!
!
!  Definition:
!
!    For the interval [-1,1], the Newton-Cotes open quadrature rule
!    estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the OPEN rule, the abscissas do not include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the  rule.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  a = -1.0D+00
  b =  1.0D+00

  do i = 1, norder
    xtab(i) = ( real ( norder + 1 - i, kind = 8 ) * a   &
              + real (              i, kind = 8 ) * b ) &
              / real ( norder + 1,     kind = 8 )
  end do

  call nc_com ( norder, a, b, xtab, weight )

  return
end
subroutine nco_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! NCO_SET sets abscissas and weights for open Newton-Cotes quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Note:
!
!    The open Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with equally spaced data.
!
!    The rules are called "open" because they do not include the interval
!    endpoints.
!
!    Most of the rules involve negative weights.  These can produce loss
!    of accuracy due to the subtraction of large, nearly equal quantities.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 7, and 9.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are rational, symmetric, and should sum to 2.
!    Some weights may be negative.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) d
  integer i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder == 1 ) then

    weight(1) = 2.0D+00

  else if ( norder == 2 ) then

    weight(1) = 1.0D+00
    weight(2) = 1.0D+00

  else if ( norder == 3 ) then

    d = 3.0D+00

    weight(1) =   4.0D+00 / d
    weight(2) = - 2.0D+00 / d
    weight(3) =   4.0D+00 / d

  else if ( norder == 4 ) then

    d = 12.0D+00

    weight(1) = 11.0D+00 / d
    weight(2) =  1.0D+00 / d
    weight(3) =  1.0D+00 / d
    weight(4) = 11.0D+00 / d

  else if ( norder == 5 ) then

    d = 10.0D+00

    weight(1) =   11.0D+00 / d
    weight(2) = - 14.0D+00 / d
    weight(3) =   26.0D+00 / d
    weight(4) = - 14.0D+00 / d
    weight(5) =   11.0D+00 / d

  else if ( norder == 6 ) then

    d = 1440.0D+00

    weight(1) =  1222.0D+00 / d
    weight(2) = - 906.0D+00 / d
    weight(3) =  1124.0D+00 / d
    weight(4) =  1124.0D+00 / d
    weight(5) = - 906.0D+00 / d
    weight(6) =  1222.0D+00 / d

  else if ( norder == 7 ) then

    d = 945.0D+00

    weight(1) =    920.0D+00 / d
    weight(2) = - 1908.0D+00 / d
    weight(3) =   4392.0D+00 / d
    weight(4) = - 4918.0D+00 / d
    weight(5) =   4392.0D+00 / d
    weight(6) = - 1908.0D+00 / d
    weight(7) =    920.0D+00 / d

  else if ( norder == 9 ) then

    d = 4536.0D+00

    weight(1) =    4045.0D+00 / d
    weight(2) = - 11690.0D+00 / d
    weight(3) =   33340.0D+00 / d
    weight(4) = - 55070.0D+00 / d
    weight(5) =   67822.0D+00 / d
    weight(6) = - 55070.0D+00 / d
    weight(7) =   33340.0D+00 / d
    weight(8) = - 11690.0D+00 / d
    weight(9) =    4045.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCO_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 to 7, and 9.'
    stop

  end if
!
!  Set the abscissas.
!
  do i = 1, norder
    xtab(i) = real ( 2 * i - norder - 1, kind = 8 ) &
            / real ( norder + 1, kind = 8 )
  end do

  return
end
subroutine radau_set ( norder, xtab, weight )
!
!*******************************************************************************
!
!! RADAU_SET sets abscissas and weights for Radau quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-2).
!
!  Note:
!
!    The Radau rule is distinguished by the fact that the left endpoint
!    (-1) is always an abscissa.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 15.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive.  The weights are not symmetric.
!    The weights should sum to 2.  WEIGHT(1) should equal 2 / NORDER**2.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) weight(norder)
!
  if ( norder == 1 ) then

    xtab(1) =   - 1.0D+00
    weight(1) =   2.0D+00

  else if ( norder == 2 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =    1.0D+00 / 3.0D+00

    weight(1) =  0.5D+00
    weight(2) =  1.5D+00

  else if ( norder == 3 ) then

    xtab(1) =   - 1.0D+00
    xtab(2) =   - 0.289897948556635619639456814941D+00
    xtab(3) =     0.689897948556635619639456814941D+00

    weight(1) =  0.222222222222222222222222222222D+00
    weight(2) =  0.102497165237684322767762689304D+01
    weight(3) =  0.752806125400934550100150884739D+00

  else if ( norder == 4 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.575318923521694112050483779752D+00
    xtab(3) =    0.181066271118530578270147495862D+00
    xtab(4) =    0.822824080974592105208907712461D+00

    weight(1) =  0.125D+00
    weight(2) =  0.657688639960119487888578442146D+00
    weight(3) =  0.776386937686343761560464613780D+00
    weight(4) =  0.440924422353536750550956944074D+00

  else if ( norder == 5 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.720480271312438895695825837750D+00
    xtab(3) =  - 0.167180864737833640113395337326D+00
    xtab(4) =    0.446313972723752344639908004629D+00
    xtab(5) =    0.885791607770964635613757614892D+00

    weight(1) =  0.08D+00
    weight(2) =  0.446207802167141488805120436457D+00
    weight(3) =  0.623653045951482508163709823153D+00
    weight(4) =  0.562712030298924120384345300681D+00
    weight(5) =  0.287427121582451882646824439708D+00

  else if ( norder == 6 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.802929828402347147753002204224D+00
    xtab(3) =  - 0.390928546707272189029229647442D+00
    xtab(4) =    0.124050379505227711989974959990D+00
    xtab(5) =    0.603973164252783654928415726409D+00
    xtab(6) =    0.920380285897062515318386619813D+00

    weight(1) =  0.555555555555555555555555555556D-01
    weight(2) =  0.319640753220510966545779983796D+00
    weight(3) =  0.485387188468969916159827915587D+00
    weight(4) =  0.520926783189574982570229406570D+00
    weight(5) =  0.416901334311907738959406382743D+00
    weight(6) =  0.201588385253480840209200755749D+00

  else if ( norder == 7 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.853891342639482229703747931639D+00
    xtab(3) =  - 0.538467724060109001833766720231D+00
    xtab(4) =  - 0.117343037543100264162786683611D+00
    xtab(5) =    0.326030619437691401805894055838D+00
    xtab(6) =    0.703842800663031416300046295008D+00
    xtab(7) =    0.941367145680430216055899446174D+00

    weight(1) =  0.408163265306122448979591836735D-01
    weight(2) =  0.239227489225312405787077480770D+00
    weight(3) =  0.380949873644231153805938347876D+00
    weight(4) =  0.447109829014566469499348953642D+00
    weight(5) =  0.424703779005955608398308039150D+00
    weight(6) =  0.318204231467301481744870434470D+00
    weight(7) =  0.148988471112020635866497560418D+00

  else if ( norder == 8 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.887474878926155707068695617935D+00
    xtab(3) =  - 0.639518616526215270024840114382D+00
    xtab(4) =  - 0.294750565773660725252184459658D+00
    xtab(5) =    0.943072526611107660028971153047D-01
    xtab(6) =    0.468420354430821063046421216613D+00
    xtab(7) =    0.770641893678191536180719525865D+00
    xtab(8) =    0.955041227122575003782349000858D+00

    weight(1) =  0.03125D+00
    weight(2) =  0.185358154802979278540728972699D+00
    weight(3) =  0.304130620646785128975743291400D+00
    weight(4) =  0.376517545389118556572129261442D+00
    weight(5) =  0.391572167452493593082499534004D+00
    weight(6) =  0.347014795634501280228675918422D+00
    weight(7) =  0.249647901329864963257869293513D+00
    weight(8) =  0.114508814744257199342353728520D+00

  else if ( norder == 9 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.910732089420060298533757956283D+00
    xtab(3) =  - 0.711267485915708857029562959544D+00
    xtab(4) =  - 0.426350485711138962102627520502D+00
    xtab(5) =  - 0.903733696068532980645444599064D-01
    xtab(6) =    0.256135670833455395138292079035D+00
    xtab(7) =    0.571383041208738483284917464837D+00
    xtab(8) =    0.817352784200412087992517083851D+00
    xtab(9) =    0.964440169705273096373589797925D+00

    weight(1) =  0.246913580246913580246913580247D-01
    weight(2) =  0.147654019046315385819588499802D+00
    weight(3) =  0.247189378204593052361239794969D+00
    weight(4) =  0.316843775670437978338000849642D+00
    weight(5) =  0.348273002772966594071991031186D+00
    weight(6) =  0.337693966975929585803724239792D+00
    weight(7) =  0.286386696357231171146705637752D+00
    weight(8) =  0.200553298024551957421165090417D+00
    weight(9) =  0.907145049232829170128934984159D-01

  else if ( norder == 10 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.927484374233581078117671398464D+00
    xtab(3) =  - 0.763842042420002599615429776011D+00
    xtab(4) =  - 0.525646030370079229365386614293D+00
    xtab(5) =  - 0.236234469390588049278459503207D+00
    xtab(6) =    0.760591978379781302337137826389D-01
    xtab(7) =    0.380664840144724365880759065541D+00
    xtab(8) =    0.647766687674009436273648507855D+00
    xtab(9) =    0.851225220581607910728163628088D+00
    xtab(10) =   0.971175180702246902734346518378D+00

    weight(1) =  0.02D+00
    weight(2) =  0.120296670557481631517310522702D+00
    weight(3) =  0.204270131879000675555788672223D+00
    weight(4) =  0.268194837841178696058554475262D+00
    weight(5) =  0.305859287724422621016275475401D+00
    weight(6) =  0.313582457226938376695902847302D+00
    weight(7) =  0.290610164832918311146863077963D+00
    weight(8) =  0.239193431714379713376571966160D+00
    weight(9) =  0.164376012736921475701681668908D+00
    weight(10) = 0.736170054867584989310512940790D-01

  else if ( norder == 11 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.939941935677027005913871284731D+00
    xtab(3) =  - 0.803421975580293540697597956820D+00
    xtab(4) =  - 0.601957842073797690275892603234D+00
    xtab(5) =  - 0.351888923353330214714301017870D+00
    xtab(6) =  - 0.734775314313212657461903554238D-01
    xtab(7) =    0.210720306228426314076095789845D+00
    xtab(8) =    0.477680647983087519467896683890D+00
    xtab(9) =    0.705777100713859519144801128840D+00
    xtab(10) =   0.876535856245703748954741265611D+00
    xtab(11) =   0.976164773135168806180508826082D+00

    weight(1) =  0.165289256198347107438016528926D-01
    weight(2) =  0.998460819079680638957534695802D-01
    weight(3) =  0.171317619206659836486712649042D+00
    weight(4) =  0.228866123848976624401683231126D+00
    weight(5) =  0.267867086189684177806638163355D+00
    weight(6) =  0.285165563941007337460004408915D+00
    weight(7) =  0.279361333103383045188962195720D+00
    weight(8) =  0.250925377697128394649140267633D+00
    weight(9) =  0.202163108540024418349931754266D+00
    weight(10) = 0.137033682133202256310153880580D+00
    weight(11) = 0.609250978121311347072183268883D-01

  else if ( norder == 12 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.949452759204959300493337627077D+00
    xtab(3) =  - 0.833916773105189706586269254036D+00
    xtab(4) =  - 0.661649799245637148061133087811D+00
    xtab(5) =  - 0.444406569781935851126642615609D+00
    xtab(6) =  - 0.196994559534278366455441427346D+00
    xtab(7) =    0.637247738208319158337792384845D-01
    xtab(8) =    0.319983684170669623532789532206D+00
    xtab(9) =    0.554318785912324288984337093085D+00
    xtab(10) =   0.750761549711113852529400825472D+00
    xtab(11) =   0.895929097745638894832914608454D+00
    xtab(12) =   0.979963439076639188313950540264D+00

    weight(1) =  0.138888888888888888888888888888D-01
    weight(2) =  0.841721349386809762415796536813D-01
    weight(3) =  0.145563668853995128522547654706D+00
    weight(4) =  0.196998534826089634656049637969D+00
    weight(5) =  0.235003115144985839348633985940D+00
    weight(6) =  0.256991338152707776127974253598D+00
    weight(7) =  0.261465660552133103438074715743D+00
    weight(8) =  0.248121560804009959403073107079D+00
    weight(9) =  0.217868879026192438848747482023D+00
    weight(10) = 0.172770639313308564306065766966D+00
    weight(11) = 0.115907480291738392750341908272D+00
    weight(12) = 0.512480992072692974680229451351D-01

  else if ( norder == 13 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.956875873668299278183813833834D+00
    xtab(3) =  - 0.857884202528822035697620310269D+00
    xtab(4) =  - 0.709105087529871761580423832811D+00
    xtab(5) =  - 0.519197779050454107485205148087D+00
    xtab(6) =  - 0.299201300554509985532583446686D+00
    xtab(7) =  - 0.619016986256353412578604857936D-01
    xtab(8) =    0.178909837597084635021931298881D+00
    xtab(9) =    0.409238231474839556754166331248D+00
    xtab(10) =   0.615697890940291918017885487543D+00
    xtab(11) =   0.786291018233046684731786459135D+00
    xtab(12) =   0.911107073689184553949066402429D+00
    xtab(13) =   0.982921890023145161262671078244D+00

    weight(1) =  0.118343195266272189349112426036D-01
    weight(2) =  0.719024162924955289397537405641D-01
    weight(3) =  0.125103834331152358133769287976D+00
    weight(4) =  0.171003460470616642463758674512D+00
    weight(5) =  0.206960611455877074631132560829D+00
    weight(6) =  0.230888862886995434012203758668D+00
    weight(7) =  0.241398342287691148630866924129D+00
    weight(8) =  0.237878547660712031342685189180D+00
    weight(9) =  0.220534229288451464691077164199D+00
    weight(10) = 0.190373715559631732254759820746D+00
    weight(11) = 0.149150950090000205151491864242D+00
    weight(12) = 0.992678068818470859847363877478D-01
    weight(13) = 0.437029032679020748288533846051D-01

  else if ( norder == 14 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.962779269978024297120561244319D+00
    xtab(3) =  - 0.877048918201462024795266773531D+00
    xtab(4) =  - 0.747389642613378838735429134263D+00
    xtab(5) =  - 0.580314056546874971105726664999D+00
    xtab(6) =  - 0.384202003439203313794083903375D+00
    xtab(7) =  - 0.168887928042680911008441695622D+00
    xtab(8) =    0.548312279917645496498107146428D-01
    xtab(9) =    0.275737205435522399182637403545D+00
    xtab(10) =   0.482752918588474966820418534355D+00
    xtab(11) =   0.665497977216884537008955042481D+00
    xtab(12) =   0.814809550601994729434217249123D+00
    xtab(13) =   0.923203722520643299246334950272D+00
    xtab(14) =   0.985270697947821356698617003172D+00

    weight(1) =  0.102040816326530612244897959184D-01
    weight(2) =  0.621220169077714601661329164668D-01
    weight(3) =  0.108607722744362826826720935229D+00
    weight(4) =  0.149620539353121355950520836946D+00
    weight(5) =  0.183127002125729654123867302103D+00
    weight(6) =  0.207449763335175672668082886489D+00
    weight(7) =  0.221369811499570948931671683021D+00
    weight(8) =  0.224189348002707794238414632220D+00
    weight(9) =  0.215767100604618851381187446115D+00
    weight(10) = 0.196525518452982430324613091930D+00
    weight(11) = 0.167429727891086278990102277038D+00
    weight(12) = 0.129939668737342347807425737146D+00
    weight(13) = 0.859405354429804030893077310866D-01
    weight(14) = 0.377071632698969142774627282919D-01

  else if ( norder == 15 ) then

    xtab(1) =  - 1.0D+00
    xtab(2) =  - 0.967550468197200476562456018282D+00
    xtab(3) =  - 0.892605400120550767066811886849D+00
    xtab(4) =  - 0.778685617639031079381743321893D+00
    xtab(5) =  - 0.630779478886949283946148437224D+00
    xtab(6) =  - 0.455352905778529370872053455981D+00
    xtab(7) =  - 0.260073376740807915768961188263D+00
    xtab(8) =  - 0.534757226797460641074538896258D-01
    xtab(9) =    0.155410685384859484319182024964D+00
    xtab(10) =   0.357456512022127651195319205174D+00
    xtab(11) =   0.543831458701484016930711802760D+00
    xtab(12) =   0.706390264637572540152679669478D+00
    xtab(13) =   0.838029000636089631215097384520D+00
    xtab(14) =   0.932997190935973719928072142859D+00
    xtab(15) =   0.987166478414363086378359071811D+00

    weight(1) =  0.888888888888888888888888888889D-02
    weight(2) =  0.542027800486444943382142368018D-01
    weight(3) =  0.951295994604808992038477266346D-01
    weight(4) =  0.131875462504951632186262157944D+00
    weight(5) =  0.162854477303832629448732245828D+00
    weight(6) =  0.186715145839450908083795103799D+00
    weight(7) =  0.202415187030618429872703310435D+00
    weight(8) =  0.209268608147694581430889790306D+00
    weight(9) =  0.206975960249553755479027321787D+00
    weight(10) = 0.195637503045116116473556617575D+00
    weight(11) = 0.175748872642447685670310440476D+00
    weight(12) = 0.148179527003467253924682058743D+00
    weight(13) = 0.114135203489752753013075582569D+00
    weight(14) = 0.751083927605064397329716653914D-01
    weight(15) = 0.328643915845935322530428528231D-01

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RADAU_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of NORDER = ', norder
    write ( *, '(a)' ) '  Legal values are 1 to 15.'
    stop

  end if

  return
end
subroutine rule_adjust ( a, b, c, d, norder, x, w )
!
!*******************************************************************************
!
!! RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
!
!
!  Discussion:
!
!    Most quadrature rules are defined on a special interval, like
!    [-1,1] or [0,1].  To integrate over an interval, the abscissas
!    and weights must be adjusted.  This can be done on the fly,
!    or by calling this routine.
!
!    If the weight function W(X) is not 1, then the W vector will
!    require further adjustment by the user.
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the definition interval.
!
!    Input, real ( kind = 8 ) C, D, the endpoints of the integration interval.
!
!    Input, integer NORDER, the number of abscissas and weights.
!
!    Input/output, real ( kind = 8 ) X(NORDER), W(NORDER), the abscissas
!    and weights.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) w(norder)
  real ( kind = 8 ) x(norder)
!
  x(1:norder) = ( ( b - x(1:norder) ) * c + ( x(1:norder) - a ) * d ) &
    / ( b - a )

  w(1:norder) = ( ( d - c ) / ( b - a ) ) * w(1:norder)

  return
end
subroutine summer ( func, norder, xtab, weight, result )
!
!*******************************************************************************
!
!! SUMMER carries out a quadrature rule over a single interval.
!
!
!  Formula:
!
!    RESULT = SUM ( 1 <= I <= NORDER ) WEIGHT(I) * FUNC ( XTAB(I) )
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ), external :: func
  integer i
  real ( kind = 8 ) result
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUMMER - Fatal error!'
    write ( *, '(a)' ) '  NORDER must be at least 1.'
    write ( *, '(a,i6)' ) '  The input value was NORDER = ', norder
    stop
  end if

  result = 0.0D+00
  do i = 1, norder
    result = result + weight(i) * func ( xtab(i) )
  end do

  return
end
subroutine summer_gk ( func, norderg, weightg, resultg, norderk, xtabk, &
  weightk, resultk )
!
!*******************************************************************************
!
!! SUMMER_GK carries out Gauss-Kronrod quadrature over a single interval.
!
!
!  Note:
!
!    The abscissas for the Gauss-Legendre rule of order NORDERG are
!    not required, since they are assumed to be the even-indexed
!    entries of the corresponding Kronrod rule.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer NORDERG, the order of the Gauss-Legendre rule.
!
!    Input, real ( kind = 8 ) WEIGHTG(NORDERG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real ( kind = 8 ) RESULTG, the approximate value of the
!    integral, based on the Gauss-Legendre rule.
!
!    Input, integer NORDERK, the order of the Kronrod rule.  NORDERK
!    must equal 2 * NORDERG + 1.
!
!    Input, real ( kind = 8 ) XTABK(NORDERK), the abscissas of the Kronrod rule.
!
!    Input, real ( kind = 8 ) WEIGHTK(NORDERK), the weights of the Kronrod rule.
!
!    Output, real ( kind = 8 ) RESULTK, the approximate value of the integral,
!    based on the Kronrod rule.
!
  implicit none
!
  integer norderg
  integer norderk
!
  real ( kind = 8 ) fk
  real ( kind = 8 ), external :: func
  integer i
  real ( kind = 8 ) resultg
  real ( kind = 8 ) resultk
  real ( kind = 8 ) weightg(norderg)
  real ( kind = 8 ) weightk(norderk)
  real ( kind = 8 ) xtabk(norderk)
!
  if ( norderk /= 2 * norderg + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUMMER_GK - Fatal error!'
    write ( *, '(a)' ) '  NORDERK must equal 2 * NORDERG + 1.'
    write ( *, '(a,i6)' ) '  The input value was NORDERG = ', norderg
    write ( *, '(a,i6)' ) '  The input value was NORDERK = ', norderk
    stop
  end if

  resultg = 0.0D+00
  resultk = 0.0D+00

  do i = 1, norderk

    fk = func ( xtabk(i) )

    resultk = resultk + weightk(i) * fk

    if ( mod ( i, 2 ) == 0 )then
      resultg = resultg + weightg(i/2) * fk
    end if

  end do

  return
end
subroutine sum_sub ( func, a, b, nsub, norder, xlo, xhi, xtab, weight, result )
!
!*******************************************************************************
!
!! SUM_SUB carries out a composite quadrature rule.
!
!
!  Discussion:
!
!    SUM_SUB assumes the original rule was written for [XLO,XHI].
!
!  Integration interval:
!
!    [ A, B ]
!
!  Integral to approximate:
!
!    Integral ( A <= X <= B ) F(X) dX
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Input, integer NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be at least 1.
!
!    Input, real ( kind = 8 ) XLO, XHI, the left and right endpoints of
!    the interval over which the quadrature rule was defined.
!
!    Input, real ( kind = 8 ) XTAB(NORDER), the abscissas of a quadrature
!    rule for the interval [XLO,XHI].
!
!    Input, real ( kind = 8 ) WEIGHT(NORDER), the weights of the quadrature rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none
!
  integer norder
!
  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer i
  integer j
  integer nsub
  real ( kind = 8 ) quad_sub
  real ( kind = 8 ) result
  real ( kind = 8 ) result_sub
  real ( kind = 8 ) x
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xtab(norder)
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_sub
  real ( kind = 8 ) weight(norder)
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if

  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a,i6)' ) '  Nonpositive value of NORDER = ', norder
    stop
  end if

  if ( nsub < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a,i6)' ) '  Nonpositive value of NSUB = ', nsub
    stop
  end if

  if ( xlo == xhi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a)' ) '  XLO = XHI.'
    stop
  end if

  volume = 0.0D+00
  result = 0.0D+00

  do j = 1, nsub

    a_sub = ( real ( nsub - j + 1, kind = 8 ) * a   &
            + real (        j - 1, kind = 8 ) * b ) &
            / real ( nsub,         kind = 8 )

    b_sub = ( real ( nsub - j, kind = 8 )     * a   &
            + real (        j, kind = 8 )     * b ) &
            / real ( nsub,     kind = 8 )

    quad_sub = 0.0D+00
    do i = 1, norder
      x = ( ( xhi - xtab(i) ) * a_sub + ( xtab(i) - xlo ) * b_sub ) &
        / ( xhi - xlo )
      quad_sub = quad_sub + weight(i) * func ( x )
    end do

    volume_sub = ( b - a ) / ( ( xhi - xlo ) * real ( nsub, kind = 8 ) )
    result_sub = quad_sub * volume_sub

    volume = volume + volume_sub
    result = result + result_sub

  end do

  return
end
subroutine sum_sub_gk ( func, a, b, nsub, norderg, weightg, resultg, norderk, &
  xtabk, weightk, resultk, error )
!
!*******************************************************************************
!
!! SUM_SUB_GK carries out a composite Gauss-Kronrod rule.
!
!
!  Integration interval:
!
!    [ A, B ]
!
!  Integral to approximate:
!
!    Integral ( A <= X <= B ) F(X) dX
!
!  Approximate integral:
!
!    H = ( B - A ) / NSUB
!    XMID(J) = A + 0.5 * H * ( 2 * J - 1 )
!
!    Sum ( 1 <= J <= NSUB )
!      Sum ( 1 <= I <= NORDERK )
!        WEIGHTK(I) * F ( XMID(J) + 0.5 * H * XTABK(I) )
!
!  Note:
!
!    The Gauss-Legendre weights should be computed by LEGCOM or LEGSET.
!    The Kronrod abscissas and weights should be computed by KRONSET.
!
!    The orders of the Gauss-Legendre and Kronrod rules must satisfy
!    NORDERK = 2 * NORDERG + 1.
!
!    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
!    plus more points, resulting in an efficient and higher order estimate.
!
!    The difference between the Gauss-Legendre and Kronrod estimates
!    is taken as an estimate of the error in the approximation to the
!    integral.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Input, integer NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer NORDERG, the order of the Gauss-Legendre rule.
!    NORDERG must be at least 1.
!
!    Input, real ( kind = 8 ) WEIGHTG(NORDERG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real ( kind = 8 ) RESULTG, the approximate value of the
!    integral based on the Gauss-Legendre rule.
!
!    Input, integer NORDERK, the order of the Kronrod rule.
!    NORDERK must be at least 1.
!
!    Input, real ( kind = 8 ) XTABK(NORDERK), the abscissas of the
!    Kronrod rule.
!
!    Input, real ( kind = 8 ) WEIGHTK(NORDERK), the weights of the
!    Kronrod rule.
!
!    Output, real ( kind = 8 ) RESULTK, the approximate value of the
!    integral based on the Kronrod rule.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the approximation
!    error.  This is computed by taking the sum of the absolute values of
!    the differences between the Gauss-Legendre and Kronrod rules
!    over each subinterval.  This is usually a good estimate of
!    the error in the value RESULTG.  The error in the Kronrod
!    estimate RESULTK is usually much smaller.
!
  implicit none
!
  integer norderg
  integer norderk
!
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ) fk
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer i
  integer j
  integer nsub
  real ( kind = 8 ) partg
  real ( kind = 8 ) partk
  real ( kind = 8 ) resultg
  real ( kind = 8 ) resultk
  real ( kind = 8 ) x
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xtabk(norderk)
  real ( kind = 8 ) weightg(norderg)
  real ( kind = 8 ) weightk(norderk)
!
  resultg = 0.0D+00
  resultk = 0.0D+00
  error = 0

  if ( a == b ) then
    return
  end if

  if ( norderk /= 2 * norderg + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB_GK - Fatal error!'
    write ( *, '(a)' ) '  NORDERK must equal 2 * NORDERG + 1.'
    write ( *, '(a,i6)' ) '  The input value was NORDERG = ', norderg
    write ( *, '(a,i6)' ) '  The input value was NORDERK = ', norderk
    stop
  end if

  h = ( b - a ) / real ( nsub, kind = 8 )

  do j = 1, nsub

    xmid = a + 0.5D0 * h * real ( 2 * j - 1, kind = 8 )

    partg = 0.0D+00
    partk = 0.0D+00

    do i = 1, norderk

      x = xmid + 0.5D0 * h * xtabk(i)
      fk = func ( x )
      partk = partk + 0.5D0 * h * weightk(i) * fk

      if ( mod ( i, 2 ) == 0 ) then
        partg = partg + 0.5D0 * h * weightg(i/2) * fk
      end if

    end do

    resultg = resultg + partg
    resultk = resultk + partk
    error = error + abs ( partk - partg )

  end do

  return
end

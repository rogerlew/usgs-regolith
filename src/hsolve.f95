  subroutine h_solve(a,b,c,h0,h1,h,l_mode)
!  Solver for humped soil production model (see Peletier & Rasmussen, 2009)
!  Code updated Oct. 26 2021, RLB
  implicit none
! LOCAL VARIABLES
  integer:: i,mmax,itmax
  real (kind = 8):: hlb,hub,htb,hmax,htemp,fl,fu,ftemp   
  real (kind = 8):: x0,dx,eps 
  real :: c_mod 
! FORMAL ARGUMENTS
  real, intent(in):: a,b,h0,h1 ! a=sec_theta(i), b=dif_ratio(zo(i))
  real, intent(out):: h ! depth
  real, intent(in):: c ! c=transport function
  logical :: l_mode
!
  mmax=32 ! maximum iterations for bisection
  eps=0.001 ! smallest allowable depth increment (1 mm)
  x0=eps ! trial depth for bisection
  hmax=30. ! maximum allowable depth (in meters), h
  hlb=0.01 ! lower bound value of depth, h
  hub=0.0 ! upper bound value of depth, h
! Bracket the range in which the depth function changes sign.  
  if(l_mode) then ! Original mode
    if(c < 0.) then ! Apply formulas only to convex slopes
      fl=hlb+(h0*c)/(b*a*exp(-hlb/(h0*a))) ! Objective function at hlb
    else
      fl=hlb
    end if
    fu=fl
    do i=1,20
      hub=hub+1.5 !increase hub to widen bracket
      ftemp=fu;htemp=hub-1.5
      fu=hub+(h0*c)/(b*a*exp(-hub/(h0*a))) ! Objective function at hub
    if(fl*fu<0.) exit ! The function changes sign between hlb and hub
    end do
  else ! Modified mode 
    c_mod = c
    if(c > 0.) c_mod = -c ! Apply same formula to convex and concave slopes
    fl=hlb+(h0*c_mod*exp(-hlb/(h0*a)))/(b*a) ! Objective function at hlb
    fu=fl
    do i=1,20
      hub=hub+1.5 !increase hub
      ftemp=fu; htemp=hub-1.5
      fu=hub+(h0*c_mod*exp(-hlb/(h0*a)))/(b*a) ! Objective function at hub
      if(fl*fu<0.) exit ! The function changes sign between hlb and hub
    end do
  endif
  if(hub>=hmax) then ! Estimated soil depth range exceeds physically reasonable bounds
    h=0. ! Depth undefined, set to zero.
    return 
  end if
  fl=ftemp; if(hub>1.5) hlb=htemp
! Choose trial depth 
  if(fl<0.) then ! Orient bracket direction
    htb=hlb
    dx=hub-hlb ! bracket width
  else
    htb=hub
    dx=hlb-hub ! bracket width
  end if
! Solve for depth by bisection  
  call dbisect(mmax,eps,a,b,c,h0,htb,dx,x0,itmax,l_mode)
  h=htb; if (h >= hmax) h=0. ! Depth undefined, set to zero.
  end subroutine h_solve 

subroutine h_solve(a,b,c,h0,h1,h,l_test)
!  Solver for humped soil production model (see Peletier & Rasmussen, 2009)
!  Code updated 3/7/2019, RLB
  implicit none
! LOCAL VARIABLES
  integer:: i,mmax,itmax
  real (kind = 8):: hlb,hub,htb,hmax,htemp,fl,fu,ftemp   
  real (kind = 8):: x0,dx,eps 
! FORMAL ARGUMENTS
  real, intent(in):: a,b,h0,h1
  real, intent(out):: h
  real, intent(in):: c
  logical :: l_test
!
  mmax=20;eps=0.001;x0=eps 
  hmax=30. 
  hlb=0.01
  hub=0.0
  if(l_test) then ! test mode
    if(-(h0/hlb)*c/(b*a) >0.) then
      fl=hlb-h0*a*log(-(hlb/h0)*(b*a)/c) 
    else
      fl=hlb
    end if
    fu=fl
    do i=1,20
      hub=hub+1.5 !increase hub
      ftemp=fu;htemp=hub-1.5
      if(-(h0/hub)*c/(b*a) <0.) then
        cycle
      else
        fu=hub-h0*a*log(-(hub/h0)*(b*a)/c)
      endif
    if(fl*fu<0.) exit
    end do
  else ! production mode 
    if((h0/hlb)*c/(b*a) >0.) then
      fl=hlb-h0*a*log((h0/hlb)*c/(b*a)) 
    else
      fl=hlb
    end if
    fu=fl
    do i=1,20
      hub=hub+1.5 !increase hub
      ftemp=fu;htemp=hub-1.5
      if((h0/hub)*c/(b*a) <0.) then
        cycle 
      else
        fu=hub-h0*a*log((h0/hub)*c/(b*a)) 
      endif
    if(fl*fu<0.) exit
    end do
  endif
  if(hub>=hmax) then
    h=0.
    return 
  end if
  fl=ftemp;if(hub>1.5)hlb=htemp
!Choose trial depth 
  if(fl<0.) then
    htb=hlb
    dx=hub-hlb
  else
    htb=hub
    dx=hlb-hub
  end if
! Solve for depth by bisection  
  call dbisect(mmax,eps,a,b,c,h0,htb,dx,x0,itmax,l_test)
  h=htb; if (h>=hmax) h=0.
end subroutine h_solve 

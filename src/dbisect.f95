subroutine dbisect(mmax,eps,a,b,c,h0,rtb,dx,x0,itmax,l_test) 
! By Rex L. Baum, May 14, 2004, converted to F95 2/6/2015,RLB
! Format updated 3/7/2019, RLB
! Uses bisection method as described in Press and others, 1986, p. 246-247.
  implicit none
! LOCAL VARIABLES
  integer:: m
  real (kind = 8):: fm
! FORMAL ARGUMENTS
  integer:: itmax,mmax
  real:: a,b,c,h0
  real (kind = 8):: eps,x0 ! , c
  real (kind = 8):: rtb,dx, tol
  logical:: l_test
  itmax=0
  tol=1.0e-10
  if(l_test) then ! test mode
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb
      if(-(x0/h0)*(b*a)/c <0.) cycle
      fm=x0-h0*a*log(-(x0/h0)*(b*a)/c) !humped soil production function 
      if(fm <= 0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then
        itmax=m
        return
      end if 
    end do 
    return
  else ! Production mode
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb
      if((x0/h0)*c/(b*a) <0.) cycle
      fm=x0-h0*a*log((h0/x0)*c/(b*a)) !humped soil production function 
      if(fm.le.0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then
        itmax=m
        return
      end if 
    end do 
    return
  end if
end subroutine dbisect

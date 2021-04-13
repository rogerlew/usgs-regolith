subroutine dbisect(mmax,eps,a,b,c,h0,rtb,dx,x0,itmax,l_mode) 
! By Rex L. Baum, 2/6/2015
! Formula updated April 2021, RLB
! Uses bisection method as described in Press and others, 1986, p. 246-247.
  implicit none
! LOCAL VARIABLES
  integer:: m
  real (kind = 8):: fm
! FORMAL ARGUMENTS
  integer:: itmax, mmax
  real:: a, b, c, h0
  real (kind = 8):: eps, x0 
  real (kind = 8):: rtb, dx, tol
  logical:: l_mode
  itmax=0 ! iteration counter
  tol=1.0e-06 ! tolerace for testing convergence
  if(l_mode) then ! Original mode
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb ! trial depth for bisection
      if((x0/h0)*(b*a)/c <0.) cycle  ! Avoid negative argument of log()
      fm=x0-h0*a*log((x0/h0)*(b*a)/c) ! humped soil production function 
      if(fm <= 0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then ! Solution converged
        itmax=m ! number of iterations to achieve convergence
        return
      end if 
    end do 
    return
  else ! Modified mode
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb
      if((x0/h0)*c/(b*a) <0.) cycle ! Avoid negative argument of log()
      fm=x0-h0*a*log((h0/x0)*c/(b*a)) !humped soil production function 
      if(fm.le.0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then ! Solution converged
        itmax=m 
        return
      end if 
    end do 
    return
  end if
end subroutine dbisect

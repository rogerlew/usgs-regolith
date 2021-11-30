subroutine dbisect(mmax,eps,a,b,c,h0,rtb,dx,x0,itmax,l_mode) 
! By Rex L. Baum, 2/6/2015
! Formula updated April 2021, RLB
! Uses bisection method as described in Press and others, 1986, p. 246-247.
  implicit none
! LOCAL VARIABLES
  integer:: m
  real (kind = 8):: fm, tol
  real:: c_mod
! FORMAL ARGUMENTS
  integer:: itmax, mmax
  real, intent(in):: a, b, c, h0 ! a=sec_theta(i), b=dif_ratio(zo(i)),c=transport function
  real (kind = 8):: eps, x0 
  real (kind = 8):: rtb, dx
  logical:: l_mode
  itmax=0 ! iteration counter
  tol=1.0e-06 ! tolerace for testing convergence
  if(l_mode) then ! Original mode
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb ! trial depth for bisection
      fm=x0+(h0*c)/(b*a*exp(-x0/(h0*a))) ! humped soil production function 
      if(fm <= 0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then ! Solution converged
        itmax=m ! number of iterations to achieve convergence
        return
      end if 
    end do 
    return
  else ! Modified mode
    c_mod = c 
    do m=1,mmax
      dx=dx/2.
      x0=dx+rtb
      if(c > 0.) c_mod = -c ! Apply same formula to convex and concave slopes
      fm=x0+(h0*c_mod*exp(-x0/(h0*a)))/(b*a) !humped soil production function 
      if(fm.le.0.) rtb=x0
      if(abs(dx) <= eps .or. abs(fm) <= tol) then ! Solution converged
        itmax=m 
        return
      end if 
    end do 
    return
  end if
end subroutine dbisect

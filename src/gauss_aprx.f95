! Apply multiple instances of moving average filter to digital elevation grid
! Apply filter in y (N-S) and x (E-W) directions separately (see www.dspguide.com/ch24/3.htm)
! use mirror-image at boundaries and preserve (skip) no-data cells
! Rex L. Baum, USGS, January 2020, Latest revision 27 Jan. 2020.
 subroutine gauss_approx(z0,cta,imax,ncol,nrow,delx,dely,null,&
                        & no_data_int,zfiltered,ztemp,temp2,npoints)
  implicit none
!
! LOCAL VARIABLES
  integer:: i, j, k, k1, l !, xedge, yedge, corner
  real::cellsum,sum1,q
!
! FORMAL ARGUMENTS
  integer, intent(in):: ncol,nrow,no_data_int,imax,cta(ncol,nrow),npoints
  real (kind = 8), intent(in):: null,delx,dely
  real, intent(in):: z0(ncol,nrow) ! input array
  real:: ztemp(ncol,nrow),temp2(ncol,nrow)
  real, intent(out):: zfiltered(imax) 
!
  write(*,*) 'Entering subroutine gauss_approx()'
!  write(*,*) 'ncol, nrow ',ncol,nrow
! Search over entire grid.  Assume that upper left corner of grid corresponds to i=1,j=1.
! i=row number (y position), j=column number (x position)
! initialize filtered array with original data
  ztemp = z0
  do l = 1,4 ! Apply moving average in each direction 4x to approximate Gaussian.
  do i=1,nrow
     do j=1,ncol
! skip no-data cells
      if(cta(j,i)==no_data_int) then
        cycle 
      endif
!     Moving average in y direction
      cellsum = float(npoints)
      sum1=0.
      do k= -(npoints-1)/2, (npoints-1)/2
        k1 = k
!       Reflect phantom points across edges 
        if (i+k < 1 .or. i+k > nrow) k1 = -k
!       Skip null (no-data) cells
        if(cta(j,i+k1) /= no_data_int) then 
          sum1 = sum1 + ztemp(j,i+k1)
        else
          cellsum = cellsum - 1.
        endif
      end do
      q = sum1/cellsum
      zfiltered(cta(j,i)) = q
      temp2(j,i) = q
    end do
  end do
  ztemp = temp2 ! Update filtered array for next iteration.
!
  do i=1,nrow
     do j=1,ncol
! skip no-data cells
      if(cta(j,i)==no_data_int) then
        cycle 
      endif
!     Moving average in x direction
      cellsum = float(npoints)
      sum1=0.
      do k= -(npoints-1)/2, (npoints-1)/2
        k1 = k
!       Reflect phantom points across edges 
        if (j+k < 1 .or. j+k > ncol) k1 = -k
!       Skip null (no-data) cells
        if(cta(j+k1,i) /= no_data_int) then 
          sum1 = sum1 + ztemp(j+k1,i)
        else
          cellsum = cellsum - 1.
        endif
      end do
      q = sum1/cellsum
      zfiltered(cta(j,i)) = q
      temp2(j,i) = q
    end do
  end do
  ztemp = temp2 ! Update filtered array for next iteration.
  end do
  write(*,*) 'Approxmate Gaussian filtering completed'
 end subroutine gauss_approx

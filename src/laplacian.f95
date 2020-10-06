subroutine laplacian(z,z1,cta,imax,col,row,d2zdx2,d2zdy2,divgradz,delx,dely,null,nodati)
! Compute 2nd derivatives and laplacian of digital elevation grids 
! Rex L. Baum, USGS, May 26, 2017, Latest revision 29 May 2019.
  implicit none
!
! LOCAL VARIABLES
!!  integer,parameter:: double=kind(1d0)
  integer:: i,j,xedge,yedge
!
! FORMAL ARGUMENTS
  integer, intent(in):: col,row,nodati,imax,cta(col,row) 
  real (kind = 8), intent(in):: null,delx,dely
  real, intent(in):: z(imax),z1(col,row)
  real, intent(out):: d2zdx2(imax),d2zdy2(imax),divgradz(imax) 
!
  write(*,*) 'Entering subroutine laplacian'
! Search over entire grid.  Assume that upper left corner of grid corresponds to i=1,j=1.
! i=row number (y position), j=column number (x position)
  do i=1,row
    do j=1,col
! Compute slope in y-direction
      yedge=0
! Test for y-edge
      if(i==1) yedge=1
      if(i>1) then 
        if(z1(j,i-1)==null) yedge=1
      end if
      if(i==row) yedge=2
      if(i<row) then
        if(z1(j,i+1)==null) yedge=2
      endif
!            write(*,*) i,j,yedge,z1(j,i)
      select case (yedge)
! 3-point formulas compute slope components at edges
! Case formulas revised 29 May 2019 to handle narrow rows or columns of cells bounded by null (no-data) cells.
      case(1) ! north edge
        if(cta(j,i)/=nodati) then
          if(cta(j,i+1)==nodati .or. cta(j,i+2)==nodati) then
            d2zdy2(cta(j,i)) = 0.
          else
            d2zdy2(cta(j,i)) =&
            & (2.d0*z(cta(j,i+1))-z(cta(j,i+2))-z(cta(j,i)))/(3.d0*dely*dely)
          endif
        endif
      case(2) ! south edge
        if(cta(j,i)/=nodati) then
          if(cta(j,i-1)==nodati .or. cta(j,i-2)==nodati) then
            d2zdy2(cta(j,i)) = 0.
          else
            d2zdy2(cta(j,i)) =&
            & (2.d0*z(cta(j,i-1))-z(cta(j,i-2))-z(cta(j,i)))/(3.d0*dely*dely)
          endif
        endif
      case default
! Central difference formula computes slope components inside grid
      if(cta(j,i)/=nodati) d2zdy2(cta(j,i)) =&
      & (z(cta(j,i-1))+z(cta(j,i+1))-2.d0*z(cta(j,i)))/(dely*dely)
      end select
! Compute slope in x-direction
      xedge=0
! Test for x-edge
      if(j==col) xedge=3
      if(j<col) then
        if(z1(j+1,i)==null) xedge=3
      endif
      if(j==1) xedge=4
      if(j>1) then
        if(z1(j-1,i)==null) xedge=4
      end if
      select case (xedge)
      case(3) ! east edge
        if(cta(j,i)/=nodati) then
          if(cta(j-1,i)==nodati .or. cta(j-2,i)==nodati) then
            d2zdx2(cta(j,i)) = 0.
          else
            d2zdx2(cta(j,i)) =&
            & (2.d0*z(cta(j-1,i))-z(cta(j-2,i))-z(cta(j,i)))/(3.d0*delx*delx)
          endif
        endif
      case(4) ! west edge
        if(cta(j,i)/=nodati) then
          if(cta(j+1,i)==nodati .or. cta(j+2,i)==nodati) then
            d2zdx2(cta(j,i)) = 0.
          else
            d2zdx2(cta(j,i)) =&
            & (2.d0*z(cta(j+1,i))-z(cta(j+2,i))-z(cta(j,i)))/(3.d0*delx*delx)
          endif
        endif
      case default
! Central difference formula computes slope components inside grid
        if(cta(j,i)/=nodati) d2zdx2(cta(j,i)) =&
        & (z(cta(j+1,i))+z(cta(j-1,i))-2.d0*z(cta(j,i)))/(delx*delx)
      end select
! Compute Laplacian of z (divergence of gradient of z)   
      if(cta(j,i)/=nodati) divgradz(cta(j,i))=d2zdx2(cta(j,i))+d2zdy2(cta(j,i))
    end do
  end do
  write(*,*) '2nd derivatives completed'
end subroutine laplacian

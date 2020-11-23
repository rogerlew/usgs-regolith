subroutine xyslope(z,z1,cta,imax,ncol,nrow,dzdx,dzdy,delx,dely,null,nodati)
! Compute x and y components of slope for digital elevation grids, 
! Rex L. Baum, USGS, August 2011, Latest revision 23 Nov 2020.
  implicit none
!
! LOCAL VARIABLES
  integer:: i,j,xedge,yedge
!
! FORMAL ARGUMENTS
  integer, intent(in):: ncol,nrow,nodati,imax,cta(ncol,nrow)
  real (kind = 8), intent(in):: null,delx,dely
  real, intent(in):: z(imax),z1(ncol,nrow) 
  real, intent(inout):: dzdx(imax),dzdy(imax)
!
  write(*,*) 'Entering subroutine xyslope'
! Search over entire grid.  Assume that upper left corner of grid corresponds to i=1,j=1.
! i=row number (y position), j=column number (x position)
  do i=1,nrow
    do j=1,ncol
      if(nrow>2) then
! Compute slope in y-direction
        yedge=0
! Test for y-edge
        if(i==1) yedge=1
        if(i>1 .and. yedge==0) then 
          if(z1(j,i-1)==null) yedge=1
        end if
        if(i==nrow) yedge=2
        if(i<nrow .and. yedge==0) then
          if(z1(j,i+1)==null) yedge=2
        else if(i<nrow .and. yedge==1) then
          if(z1(j,i+1)==null) yedge=5        
        end if
        select case (yedge)
! 3-point formulas compute slope components at edges
! Case formulas revised 5/29/2019 to handle narrow rows/columns of data cells between null cells
        case(1) ! north edge
          if(cta(j,i)/=nodati) then
            if(cta(j,i+1)==nodati) then ! not enough ppoints to compute slope
              dzdy(cta(j,i)) = 0.
            else if(i<=nrow-2) then
              if(cta(j,i+2)==nodati) then ! Use 2-point forward difference formula
                dzdy(cta(j,i)) = (-z(cta(j,i+1))+z(cta(j,i)))/dely
              end if
            else if(i<nrow-2) then ! Use 3-point forward difference formula
              dzdy(cta(j,i)) =&
              & (-4.d0*z(cta(j,i+1))+z(cta(j,i+2))+3.d0*z(cta(j,i)))/(2.d0*dely)
            endif
          endif
        case(2) ! south edge
          if(cta(j,i)/=nodati) then
            if(cta(j,i-1)==nodati) then ! not enough ppoints to compute slope
              dzdy(cta(j,i)) = 0.
            else if(i>=2) then
              if(cta(j,i-2)==nodati) then ! Use 2-point backward difference formula
                dzdy(cta(j,i)) = (z(cta(j,i-1))-z(cta(j,i)))/dely
              end if
            else if(i>2) then  ! Use 3-point backward difference formula
              dzdy(cta(j,i)) =&
              & (4.d0*z(cta(j,i-1))-z(cta(j,i-2))-3.d0*z(cta(j,i)))/(2.d0*dely)
            endif
          endif
        case(5) ! isolated finger, not enough ppoints to compute slope
          if(cta(j,i)/=nodati) dzdy(cta(j,i)) = 0.
        case default
! Central difference formula computes slope components inside grid
          if(cta(j,i)/=nodati) dzdy(cta(j,i)) =&
        & (z(cta(j,i-1))-z(cta(j,i+1)))/(2.d0*dely) ! row number, i, increases in negative y-direction.
        end select
      end if
! Compute slope in x-direction
      if(ncol >2) then
        xedge=0
! Test for x-edge
        if(j==1) xedge=4
        if(j>1) then
          if(z1(j-1,i)==null) xedge=4
        end if
        if(j==ncol) xedge=3
        if(j<ncol .and. xedge==0) then
          if(z1(j+1,i)==null) xedge=3
        else if(j<ncol.and. xedge==4) then
          if(z1(j+1,i)==null) yedge=6        
        endif
        select case (xedge)
        case(3) ! east/right edge
          if(cta(j,i)/=nodati) then
            if(cta(j-1,i)==nodati) then
              dzdx(cta(j,i)) = 0.
            else if(j>=2) then
              if(cta(j-2,i)==nodati) then ! Use 2-point backward difference formula
                dzdx(cta(j,i)) = (z(cta(j,i))-z(cta(j-1,i)))/delx
              end if
            else if(j>2) then ! Use 3-point backward difference formula
              dzdx(cta(j,i)) =&
              & (-4.d0*z(cta(j-1,i))+z(cta(j-2,i))+3.d0*z(cta(j,i)))/(2.d0*delx)
            endif
          endif
        case(4) ! west/left edge
          if(cta(j,i)/=nodati) then
            if(cta(j+1,i)==nodati) then ! not enough ppoints to compute slope
              dzdx(cta(j,i)) = 0.
            else if (j<=ncol-2) then ! Use 2-point forward difference formula
              if (cta(j+2,i)==nodati) then
                dzdx(cta(j,i)) = (z(cta(j+1,i))-z(cta(j,i)))/delx
              end if
            else if(j<ncol-2) then ! Use 3-point forward difference formula
              dzdx(cta(j,i)) =&
              & (4.d0*z(cta(j+1,i))-z(cta(j+2,i))-3.d0*z(cta(j,i)))/(2.d0*delx)
            endif
          endif
        case(6) ! isolated finger, not enough ppoints to compute slope
          if(cta(j,i)/=nodati) dzdx(cta(j,i)) = 0.
        case default
! Central difference formula computes slope components inside grid
          if(cta(j,i)/=nodati) dzdx(cta(j,i)) =&
          & (z(cta(j+1,i))-z(cta(j-1,i)))/(2.d0*delx)
        end select
      end if
    end do
  end do
  write(*,*) '1st derivatives completed'
end subroutine xyslope

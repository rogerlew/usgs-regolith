subroutine laplacian(z,z1,cta,imax,col,row,d2zdx2,d2zdy2,divgradz,delx,dely,null,no_data_int)
! Compute 2nd derivatives and laplacian of digital elevation grids 
! Rex L. Baum, USGS, May 26, 2017, Latest revision 23 Nov 2020.
  implicit none
!
! LOCAL VARIABLES
!!  integer,parameter:: double=kind(1d0)
  integer:: i,j,xedge,yedge
!
! FORMAL ARGUMENTS
  integer, intent(in):: col,row,no_data_int,imax,cta(col,row) 
  real (kind = 8), intent(in):: null,delx,dely
  real, intent(in):: z(imax),z1(col,row)
  real, intent(out):: d2zdx2(imax),d2zdy2(imax),divgradz(imax) 
!
  write(*,*) 'Entering subroutine laplacian'
! Search over entire grid.  Assume that upper left corner of grid corresponds to i=1,j=1.
! i=row number (y position), j=column number (x position)
  do i=1,row
    do j=1,col
      if(row>2) then
! Compute 2nd derivatives in y-direction
        yedge=0
! Test for y-edge
        if(i==1) yedge=1
        if(i>1) then 
          if(z1(j,i-1)==null) yedge=1
        end if
        if(i==row) yedge=2
        if(i<row .and. yedge==0) then
          if(z1(j,i+1)==null) yedge=2
        else if(i<row .and. yedge==1) then
          if(z1(j,i+1)==null) yedge=5        
        end if
!
        select case (yedge)
! 3-point formulas compute derivatives components at edges
! Case formulas revised 29 May 2019 to handle narrow rows or columns of cells bounded by null (no-data) cells.
        case(1) ! north edge
          if(cta(j,i) /= no_data_int) then
            if(i>=row-1) then
              d2zdy2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivative
            else
              if(cta(j,i+1) == no_data_int .or. cta(j,i+2) == no_data_int) then
                d2zdy2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivative
                cycle
              else ! Apply forward-difference formula at edge
                d2zdy2(cta(j,i)) =&
                & (z(cta(j,i+2)) -2.d0*z(cta(j,i+1)) + z(cta(j,i)))/(dely*dely)
              end if
            end if
          endif
        case(2) ! south edge
          if(cta(j,i) /= no_data_int) then
            if(i<=2) then
              d2zdy2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
            else
              if(cta(j,i-1) == no_data_int .or. cta(j,i-2) == no_data_int) then
                d2zdy2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
              else! Apply backward-difference formula at edge
                d2zdy2(cta(j,i)) =&
                & (z(cta(j,i-2)) -2.d0*z(cta(j,i-1)) + z(cta(j,i)))/(dely*dely)
              endif
            end if
          endif
        case(5) ! isolated finger, not enough points to compute the 2nd derivatives
          if(cta(j,i)/=no_data_int) d2zdy2(cta(j,i)) = 0.
        case default
! Central difference formula computes 2nd derivative inside grid
        if(cta(j,i)/=no_data_int) d2zdy2(cta(j,i)) =&
        & (z(cta(j,i-1)) + z(cta(j,i+1)) - 2.d0*z(cta(j,i)))/(dely*dely)
        end select
      end if
!
! Compute 2nd derivatives in x-direction
      if(col >2) then
        xedge=0
! Test for x-edge
        if(j==1) xedge=4
        if(j>1) then
          if(z1(j-1,i)==null) xedge=4
        end if
        if(j==col) xedge=3
        if(j<col .and. xedge==0) then
          if(z1(j+1,i)==null) xedge=3
        else if(j<col.and. xedge==4) then
          if(z1(j+1,i)==null) yedge=6        
        endif
!
        select case (xedge)
        case(3) ! east edge
          if(cta(j,i) /= no_data_int) then
            if (j<=2) then
              d2zdx2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
            else
              if(cta(j-1,i) == no_data_int .or. cta(j-2,i) == no_data_int) then
                d2zdx2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
              else! Apply backward-difference formula at edge
                d2zdx2(cta(j,i)) =&
                & (z(cta(j-2,i)) -2.d0*z(cta(j-1,i)) + z(cta(j,i)))/(delx*delx)
              endif
            end if
          endif
        case(4) ! west edge
          if(cta(j,i) /= no_data_int) then
            if(j>=col-1)then
              d2zdx2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
            else
              if(cta(j+1,i) == no_data_int .or. cta(j+2,i) == no_data_int) then
                d2zdx2(cta(j,i)) = 0. ! not enough points to compute the 2nd derivatives
              else! Apply forward-difference formula at edge
                d2zdx2(cta(j,i)) =&
                & (z(cta(j+2,i)) -2.d0*z(cta(j+1,i)) + z(cta(j,i)))/(delx*delx)
              endif
            end if
          endif
        case(6) ! isolated finger, not enough points to compute the 2nd derivatives
          if(cta(j,i)/=no_data_int) d2zdx2(cta(j,i)) = 0.
        case default
! Central difference formula computes 2nd derivative inside grid
          if(cta(j,i)/=no_data_int) d2zdx2(cta(j,i)) =&
          & (z(cta(j+1,i))+z(cta(j-1,i))-2.d0*z(cta(j,i)))/(delx*delx)
        end select
      end if
! Compute Laplacian of z (divergence of gradient of z)   
      if(cta(j,i)/=no_data_int) divgradz(cta(j,i))=d2zdx2(cta(j,i))+d2zdy2(cta(j,i))
    end do
  end do
  write(*,*) '2nd derivatives completed'
end subroutine laplacian

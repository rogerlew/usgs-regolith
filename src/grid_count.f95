! procdure to create 2-d array that maps grid cell numbers (1-d array) to i-j (row-column )coordinates
! 2 Feb 2015, RLB
subroutine grid_count(ncol,nrow,imax,nodati,nodat,cell_row,cell_column,cta,pf1)
  implicit none
!!  integer,parameter:: double=kind(1d0)
! FORMAL ARGUMENTS
  integer, intent(in)::nrow,ncol,nodati,imax
  integer, intent(out)::cell_row(imax),cell_column(imax),cta(ncol,nrow)
  real, intent(in):: pf1(nrow*ncol)
  real (kind = 8), intent(in)::nodat
! LOCAL VARIABLES        
  integer::i,j,k,ctr
!
  cta=nodati
  ctr=0
  do i=1,nrow ! i=row number (y-position) 
    do j=1,ncol ! j=column number (x-position)
      k=j+ncol*(i-1)
       if(pf1(k)/=nodat) then
         ctr=ctr+1
!         write(*,*) k, ctr
         cta(j,i)=ctr
         cell_row(ctr)=i
         cell_column(ctr)=j
       end if
    end do
  end do
!  write(u(1),*) 'i,j,cta(i,j)'
end subroutine grid_count

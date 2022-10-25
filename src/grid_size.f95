! determine grid size parameters RLB 4/18/2011
  subroutine grid_size(elevfil,elfoldr,init, data_cells,row,col,nwf,ulog,usiz,uelev,&
   & celsiz,nodat)
  implicit none
! FORMAL ARGUMENTS
  character(len=255), intent(in):: elevfil,init
  character(len=224), intent(out):: elfoldr
  integer, intent(in):: ulog,usiz,uelev !uelev=u(12)
  integer, intent(out):: data_cells,row,col,nwf
  real(kind = 8), intent(inout):: celsiz, nodat
! LOCAL VARIABLES
  integer:: i,patlen
  logical:: ans,lnfil
  character (len=255):: infil, outfil, size_heading
  character (len=2):: pid(3)
  character (len=14):: header(6)
!
  write(ulog,*) 'Entering subroutine grid_size()'
  pid=(/'TI','GM','TR'/)
  
  patlen=scan(elevfil,'/\',.true.) ! find end of folder name
  elfoldr=elevfil(1:patlen) ! path to elevation grid
  elfoldr=adjustl(elfoldr)
  ans=.false.
  do i=1,3 ! check for presence of grid size file
    infil=trim(elfoldr)//pid(i)//'grid_size.txt'
    infil=adjustl(infil)
    inquire (file=trim(infil),exist=ans)
    write(*,*) trim(infil), ans
    if(ans) exit
  end do
  if(ans) then ! read existing grid size file
    open (usiz,file=trim(infil),status='unknown',err=420)
    read (usiz,'(a)') size_heading
    read (usiz,*) data_cells,row,col,nwf
    close (usiz)
    nodat=-9999.d0 ! set default value for nodat
  else ! read grid size data from elevation grid and create grid size file
    write(*,*) 'Obtaining grid-size parameters from elevation grid' 
    infil=elevfil; infil=adjustl(infil)
    inquire(file=infil,exist=lnfil)
    write(*,*) 'Status of elevfil:', lnfil, trim(elevfil)
    if(lnfil)then
      call ssizgrd(row,col,celsiz,nodat,data_cells,uelev,infil,header,ulog)
      outfil=trim(elfoldr)//pid(3)//'grid_size.txt'
      outfil=adjustl(outfil)
      open (usiz,file=trim(outfil),status='unknown',err=410)
      write (usiz,*) 'data_cells      row      col      nwf'
      nwf=1 ! dsctr is computed by TopoIndex; dsctr=1 is default value for no runoff routing.
      write (usiz,*) data_cells,row,col,nwf
      write (usiz,*) ''
      close (usiz)
      write (ulog,*) 'Grid size parameters from ', trim(infil)
      write (ulog,*) 'Created new size file, ', trim(outfil)
      write (ulog,*) 'data_cells      row      col      nwf'
      write (ulog,*) data_cells,row,col,nwf
      return
    else
      write(ulog,*) 'Elevation grid ', trim(infil), ' was not found'
      write(ulog,*) 'Edit path name of elevation grid in ', trim(init),&
                  &', then restart program.'
      write(*,*) 'Elevation grid ', trim(infil), ' was not found'
      write(*,*) 'Edit path name of elevation grid in ', trim(init),&
                  &', then restart program.'
      close(ulog)
      stop 'grid_size.f95, line 60'
    endif
  end if
  write (ulog,*) 'Grid size parameters from ', trim(infil)
  write (ulog,*) trim(adjustl(size_heading))
  write (ulog,*) data_cells,row,col,nwf
  return
! Error reporting
410 continue
  write (*,*) 'Error opening output file in subroutine grid_size()'
  write (*,*) '--> ',outfil
  write (*,*) 'Check file path and status'
  write (ulog,*) 'Error opening output file in subroutine grid_size()'
  write (ulog,*) '--> ',outfil
  write (ulog,*) 'Check file path and status'
  write(*,*) 'Press RETURN to exit'
  read*
  close (ulog)
  close (usiz)
  stop '410'
420 continue
  write (*,*) 'Error opening input file in subroutine grid_size()'
  write (*,*) '--> ',infil
  write (*,*) 'Check file path and status'
  write (ulog,*) 'Error opening input file in subroutine grid_size()'
  write (ulog,*) '--> ',infil
  write (ulog,*) 'Check file path and status'
  write(*,*) 'Press RETURN to exit'
  read*
  close (ulog)
  stop '420'
  end subroutine grid_size

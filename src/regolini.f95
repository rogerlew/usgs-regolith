Module read_inputs
implicit none

contains

! reads regolith initialization file and starts log file, R.L. Baum, USGS  
  subroutine regolini(ulog,uini,dg2rad,trans_model,chan_thresh,chan_depth,&
        & num_zones,max_zones,num_steps,hump_prod,sc_deg,suffix,folder,&
        & elevfil,slopefil,flo_accfil,pv_curvfil,ndxfil,zonfil,heading,lasc,&
        & h0,sc,dif_ratio,depth_max,depth_min,C0,C2,zon,vrsn,bldate,date,time,&
        & outfil,topoSmooth,soilSmooth,n_points,l_deriv,l_test,power)
  implicit none
! LOCAL VARIABLES
!!  integer,parameter:: double=kind(1d0)
  integer:: linct ,i, iz
  character (len=255):: init,title
  character (len=31):: scratch
  logical:: ans, reset_flag
! FORMAL ARGUMENTS 
  integer, intent(in):: ulog,uini
  integer, intent(out)::num_steps, max_zones, num_zones, n_points
  integer, allocatable, intent(out):: zon(:)
  real, allocatable, intent(out):: sc_deg(:),h0(:),sc(:),dif_ratio(:)
  real, allocatable, intent(out):: depth_max(:),depth_min(:),C0(:),C2(:)
  real, intent(out):: chan_thresh,chan_depth,power ! threshold upslope contributing area for channels
  real(kind = 8),intent(in)::dg2rad
  logical, allocatable, intent(out):: hump_prod(:)  
  logical, intent(out):: lasc, l_deriv, l_test
  logical, intent(out):: topoSmooth,soilSmooth
  character (*), intent(in):: vrsn,bldate,date,time
  character (*), intent(out):: trans_model,suffix,folder,heading(16)
  character (*), intent(out):: elevfil,slopefil,flo_accfil,pv_curvfil,ndxfil,zonfil 
  character (*), intent(inout):: outfil
  init='rg_in.txt'
  ans=.false.; reset_flag=.false.
  inquire (file=trim(init),exist=ans)
  if(ans) then
    open (uini,file=trim(init),status='old',err=201)
    write (*,*) 'Opening default initialization file'
  else
    write (*,*) 'Cannot locate default initialization file, <rg_in.txt>'
    write (*,*) 'Type name of initialization file and'
    write (*,*) 'press RETURN to continue'
    read (*,'(a)') init
    open (uini,file=trim(init),status='old',err=201)
  end if
! Read contents of initialization file.
  linct=1
! Project title  
    read (uini,'(a)',err=420) heading(1); linct=linct+1
    read (uini,'(a)',err=420) title; linct=linct+1
! Model Identification code
    read (uini,'(a)',err=420) heading(2); linct=linct+1
    read (uini,'(a)',err=420) trans_model; linct=linct+1
! Channel parameters   
    read (uini,'(a)',err=420) heading(3); linct=linct+1
    read (uini,*,err=420) chan_thresh,chan_depth 
        linct=linct+1
! Number of soil-depth zones, maximum zone number, number of steps for NDSD
    read (uini,'(a)',err=420) heading(4); linct=linct+1
    read (uini,*,err=420) num_zones, max_zones, num_steps
    linct=linct+1
    if(num_zones<=1) then 
      max_zones=1
      num_zones=1
    end if
! Allocate and initialize soil zone parameters
    allocate(h0(max_zones), sc_deg(max_zones), sc(max_zones), dif_ratio(max_zones))
    allocate(depth_max(max_zones), depth_min(max_zones), C0(max_zones), C2(max_zones))
    allocate(hump_prod(max_zones),zon(max_zones))
    h0=0.; sc_deg=0.; sc=0.; dif_ratio=0.; 
    depth_max=0.; depth_min=0.; C0=0.; C2=0.
    hump_prod=.false.; zon=0
! Soil depth zone parameters
  do i=1,num_zones
    read (uini,*,err=420) scratch,iz ! property zone number (NOTE, USING iz ALLOWS ZONE NUMBERS TO BE OUT OF ORDER, a maximum zone number needed to skip zones, would require storing empty zones used in larger area, but not present current extent)
    linct=linct+1
    read (uini,'(a)',err=420) heading(5); linct=linct+1
    read (uini,*,err=420) h0(iz),sc_deg(iz),dif_ratio(iz),depth_max(iz),&
        & depth_min(iz),C0(iz),C2(iz),hump_prod(iz) 
    linct=linct+1
    sc(iz)=tan(sc_deg(iz)*dg2rad) ! convert angle of stability to radians and compute tangent
    zon(i) = iz
  end do
! Model exponent
    read (uini,'(a)',err=420) heading(6); linct=linct+1
    read (uini,*,err=420) power; linct=linct+1
!  path names of input files
    read (uini,'(a)',err=420) heading(7); linct=linct+1
! File name of digital elevation grid (elevfil)
    read (uini,'(a)',err=420) elevfil; linct=linct+1
! File name of slope grid (slopefil)
    read (uini,'(a)',err=420) heading(8); linct=linct+1
    read (uini,'(a)',err=420) slopefil; linct=linct+1
! File name of flow-accumulation grid (flo_accfil)
    read (uini,'(a)',err=420) heading(9); linct=linct+1
    read (uini,'(a)',err=420) flo_accfil; linct=linct+1
! File name of elevation index grid file (ndxfil)
    read (uini,'(a)',err=420) heading(10); linct=linct+1
    read (uini,'(a)',err=420) ndxfil; linct=linct+1
! File name of plan-view curvature grid (pv_curvfil)
    read (uini,'(a)',err=420) heading(11); linct=linct+1
    read (uini,'(a)',err=420) pv_curvfil; linct=linct+1
! File name of soil zone grid (zonfil)
    read (uini,'(a)',err=420) heading(12); linct=linct+1
    read (uini,'(a)',err=420) zonfil; linct=linct+1
! Folder (directory) where output grid files will be stored  (folder)
    read (uini,'(a)',err=420) heading(13); linct=linct+1
    read (uini,'(a)',err=420) folder; linct=linct+1  
! Identification code to be added to names of output files (suffix)
    read (uini,'(a)',err=420) heading(14); linct=linct+1
    read (uini,'(a)',err=420) suffix; linct=linct+1 
! Specify grid file extension (.asc if true, default is .txt), Output derivatives?, Test mode?
    read (uini,'(a)',err=420) heading(15); linct=linct+1
    read (uini,*,err=420) lasc, l_deriv, l_test; linct=linct+1  
! Specify topographic smoothing, soil smoothing, width of running average
    read (uini,'(a)',err=420) heading(16); linct=linct+1
    read (uini,*,err=420) topoSmooth, soilSmooth, n_points; linct=linct+1  
    if(n_points < 3) then ! Minimum width of running average
      n_points = 3 ! 
      reset_flag = .true.
    else if(mod(n_points,2) == 0) then ! Even integer
      n_points = n_points+1
      reset_flag = .true.
    endif
!-----------------------------------------    
  close(uini)
!***********************************************
! Open log file in output folder and write copy of initialization file.
  folder=adjustl(folder)
  suffix=adjustl(suffix)
  trans_model=adjustl(trans_model)
  outfil=trim(folder)//'RG_Log_'//trim(trans_model)//'_'//trim(suffix)//'.txt'
  outfil=adjustl(outfil)
  open (ulog,file=trim(outfil),status='unknown',err=411)
!   write (ulog,*) 'initialization file -->',init
!   write (ulog,*) '-- LISTING OF INITIALIZATION FILE --'  
  write (ulog,*) ''
  write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
  write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
  write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  write (ulog,*) 'initialization file --> ',init
  write (ulog,*) '-- LISTING OF INITIALIZATION FILE --'  
! write copy of data to log file
!! FINISH REARRANGING THE OUTPUT TO MATCH NEW INPUT ORDER.
!! AFTER FINISH REVISING REGOLINI, THEN REVISE SOIL MODELS TO USE ZONES
  heading(1)=adjustl(heading(1))
  title=adjustl(title)
  write (ulog,*) trim(heading(1))
  write (ulog,*) trim(title)
  write (*,*) trim(title)
  heading(2)=adjustl(heading(2))
  trans_model=adjustl(trans_model)
  write (ulog,*) trim(heading(2))
  write (ulog,*) trim(trans_model)
  heading(3)=adjustl(heading(3))
  write (ulog,*) trim(heading(3))
  write (ulog,*) chan_thresh,chan_depth
  heading(4)=adjustl(heading(4))
  write (ulog,*) trim(heading(4))
  write (ulog,*) num_zones, max_zones, num_steps
  heading(5)=adjustl(heading(5))
  do i=1,num_zones
    write (ulog,*) trim(scratch),': ',zon(i)
    write (ulog,*) trim(heading(5))
    write (ulog,*) h0(zon(i)),sc_deg(zon(i)),dif_ratio(zon(i)),&
      & depth_max(zon(i)),depth_min(zon(i)),C0(zon(i)),C2(zon(i)),hump_prod(zon(i))
  end do
  heading(6)=adjustl(heading(6))
  write (ulog,*) trim(heading(6))
  write (ulog,*) power 
  heading(7)=adjustl(heading(7))
  elevfil=adjustl(elevfil)
  write (ulog,*) trim(heading(7))
  write (ulog,*) trim(elevfil)
  heading(8)=adjustl(heading(8))
  slopefil=adjustl(slopefil)
  write (ulog,*) trim(heading(8))
  write (ulog,*) trim(slopefil)
  heading(9)=adjustl(heading(9))
  flo_accfil=adjustl(flo_accfil)
  write (ulog,*) trim(heading(9))
  write (ulog,*) trim(flo_accfil)
  heading(10)=adjustl(heading(10))
  ndxfil=adjustl(ndxfil)
  write (ulog,*) trim(heading(10))
  write (ulog,*) trim(ndxfil)
  heading(11)=adjustl(heading(11))
  pv_curvfil=adjustl(pv_curvfil)
  write (ulog,*) trim(heading(11))
  write (ulog,*) trim(pv_curvfil)
  heading(12)=adjustl(heading(12))
  zonfil=adjustl(zonfil)
  write (ulog,*) trim(heading(12))
  write (ulog,*) trim(zonfil)
!  location of output files  
  heading(13)=adjustl(heading(13))
  folder=adjustl(folder)
  write (ulog,*) trim(heading(13))
  write (ulog,*) trim(folder)
!  output-file ID code
  heading(14)=adjustl(heading(14))
  suffix=adjustl(suffix)
  write (ulog,*) trim(heading(14))
  write (ulog,*) trim(suffix)
  heading(15)=adjustl(heading(15))
  write (ulog,*) trim(heading(15))
  write (ulog,*) lasc, l_deriv, l_test
  heading(16)=adjustl(heading(16))
  write (ulog,*) trim(heading(16))
  write (ulog,*) topoSmooth, soilSmooth, n_points
  if (reset_flag) write(ulog,*) 'Incorrect value of running average <n_points> was reset to: ', n_points
  write (ulog,*) '-- END OF INITIALIZATION DATA --'  
  write (ulog,*) ''
  write (ulog,*) title
  write (ulog,*) ''
  return
  201  continue
! Report error opening initialization file.
    write (*,*) '*** Error opening intialization file ***'
    write (*,*) '--> ',trim(init)
    write (*,*) 'Check file location and name'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=412)
    write (ulog,*) ''
    write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
    write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
    write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
    write (ulog,*) '*** Error opening intialization file ***'
    write (ulog,*) '--> ',trim(init)
    write (ulog,*) 'Check file location and name'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '201'
  420  continue
! Report error reading initialization file.
    write (*,*) 'Error reading initialization file'
    write (*,*) '--> ',trim(init), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=412)
    write (ulog,*) ''
    write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
    write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
    write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
    write (ulog,*) 'Error reading initialization file'
    write (ulog,*) '--> ',trim(init), 'at line ',linct
    write (ulog,*) 'Check file contents and organization'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '420'
  411  continue
    write (*,*) 'Error opening output file in subroutine regolini()'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file path and status'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=412)
    write (ulog,*) 'Error opening output file in subroutine regolini()'
    write (ulog,*) '--> ',outfil
    write (ulog,*) 'Check file path and status'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '411'
  412  continue
    write (*,*) 'Error opening output file in subroutine regolini()'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file path and status'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '412'
  end subroutine regolini
  
end module read_inputs

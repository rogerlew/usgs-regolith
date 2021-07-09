Module read_inputs
implicit none

contains

! reads regolith initialization file and starts log file, R.L. Baum, USGS  
  subroutine regolini(ulog,uini,init,dg2rad,trans_model,chan_thresh,chan_depth,&
        & num_zones,max_zones,num_steps,hump_prod,theta_c_deg,suffix,folder,&
        & elevfil,slopefil,flo_accfil,pv_curvfil,ndxfil,zonfil,heading,lasc,&
        & h0,sc,dif_ratio,depth_max,depth_min,C0,C1,C2,zon,vrsn,bldate,date,time,&
        & outfil,topoSmooth,soilSmooth,n_points,l_deriv,l_mode,power)
  implicit none
! LOCAL VARIABLES
  integer:: linct ,i, j, iz, patlen
  character (len=400):: msg(4)
  character (len=255):: title
  character(len=224):: elfoldr
  character (len=31):: scratch
  character (len=8):: model_dscrp
  character (len=4):: valid_model_id(14), model_type, lin_num
  logical:: ans, reset_flag, id_valid, out_of_range
! FORMAL ARGUMENTS 
  integer, intent(in):: ulog,uini
  integer, intent(out)::num_steps, max_zones, num_zones, n_points
  integer, allocatable, intent(out):: zon(:)
  real, allocatable, intent(out):: theta_c_deg(:),h0(:),sc(:),dif_ratio(:)
  real, allocatable, intent(out):: depth_max(:),depth_min(:),C0(:),C1(:),C2(:)
  real, intent(out):: chan_thresh,chan_depth,power ! threshold upslope contributing area for channels
  real(kind = 8),intent(in)::dg2rad
  logical, allocatable, intent(out):: hump_prod(:)  
  logical, intent(out):: lasc, l_deriv, l_mode
  logical, intent(out):: topoSmooth,soilSmooth
  character (len=4), intent(out):: trans_model
  character (len=7), intent(in):: vrsn
  character (len=8), intent(in):: date
  character (len=10), intent(in):: time
  character (len=11), intent(in):: bldate
  character (len=16), intent(out):: suffix
  character (len=224), intent(out):: folder
  character (len=255), intent(out):: init,elevfil,slopefil,flo_accfil,pv_curvfil,ndxfil,zonfil 
  character (len=255), intent(inout):: outfil,heading(18)
  init='rg_in.txt'; init = adjustl(init)
  ans=.false.; reset_flag=.false.
  valid_model_id=(/'DRS1','DRS2','DRS3','LCSD','NASD','NDSD','NSD ','NSDA',&
                &'WNDX','ESD ','CESD','PSD ','LASD', 'LRSC'/)
  inquire (file=trim(init),exist=ans)
  if(ans) then
    open (uini,file=trim(init),status='old',err=201)
    write (*,*) 'Opening default initialization file'
  else
    write (*,*) 'Cannot locate default initialization file, <rg_in.txt>'
    write (*,*) 'Type name of initialization file and'
    write (*,*) 'press RETURN to continue'
    read (*,'(a)') init; init = adjustl(init)
    open (uini,file=trim(init),status='old',err=201)
  end if
! Read contents of initialization file.
  linct=1
! Project title  
    read (uini,'(a)',err=202) heading(1); linct=linct+1
    read (uini,'(a)',err=202) title; linct=linct+1
! Model Identification code, Original mode (.true.) or Modified mode (.false.)?
    read (uini,'(a)',err=202) heading(2); linct=linct+1
    read (uini,*,err=202) trans_model, l_mode; linct=linct+1
    trans_model = adjustl(trans_model)
!    write(*,*) trans_model, l_mode
! Test for valid model ID code
  id_valid=.false.
  do i=1,14
    if(trim(trans_model) == trim(adjustl(valid_model_id(i)))) id_valid=.true.
  end do
  if(.not. id_valid) then
    write (*,*) ''
    write(*,*) '***###*** Invalid model ID detected. ***###***'
    write(*,*) trim(trans_model)
    write(*,*) 'Edit ', trim(init), ' to correct model ID code, then restart program.'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=212)
    write(ulog,*) '***###*** Invalid model ID detected. ***###***'
    write(ulog,*) trim(trans_model)
    write(ulog,*) 'Edit ', trim(init), ' to correct model ID code, then restart program.'
    close(uini); close(ulog)
    stop 'regolini.f95, lines 64 - 83'
  endif
! Translate alternate model codes to internal codes  
  ! ESD model is DRS1
  if(trim(trans_model) == trim(adjustl(valid_model_id(1)))) trans_model=valid_model_id(10)
  ! CESD model is DRS3
  if(trim(trans_model) == trim(adjustl(valid_model_id(3)))) trans_model=valid_model_id(11)
  ! PSD model is DRS2
  if(trim(trans_model) == trim(adjustl(valid_model_id(2)))) trans_model=valid_model_id(12)
  ! LASD model is WNDX
  if(trim(trans_model) == trim(adjustl(valid_model_id(9)))) trans_model=valid_model_id(13)
!  Identify model type
  if(trans_model(1:1) == 'N' .or. trans_model(1:2) == 'LC') then
    model_type = 'PROC' ! Process-based
  else
    model_type = 'EMPR' ! Empirical
  endif
! Channel parameters   
    read (uini,'(a)',err=202) heading(3); linct=linct+1
    read (uini,*,err=202) chan_thresh,chan_depth 
    if (chan_thresh < 0. .or. chan_depth < 0.) then
      write(*,*) ''
      write(*,*) 'Negative value of channel threshold or channel depth at line ',&
            &linct,' of ',trim(init)
      write(*,*) 'Edit ',trim(init), ', then restart program.'
      outfil='RegolithLog.txt'; outfil=adjustl(outfil)
      open (ulog,file=trim(outfil),status='unknown',err=212)
      write(ulog,*) ''
      write(ulog,*) 'Negative value of channel threshold or channel depth at line ',&
            &linct,' of ',trim(init)
      write(ulog,*) 'Edit ',trim(init), ', then restart program.'
      close(uini); close (ulog)
      stop 'regolini.f95, lines 101 - 115'
    endif
    linct=linct+1
! Number of soil-depth zones, maximum zone number
    read (uini,'(a)',err=202) heading(4); linct=linct+1
    read (uini,*,err=202) num_zones, max_zones  !, num_steps
    linct=linct+1
    if(num_zones < 1 .or. max_zones < 1 .or. max_zones < num_zones) then 
        write(*,*) ''
        write(*,*) 'Edit ', trim(init), ' to correct the number of zones or &
            &maximum number of zones, then restart program.'
        outfil='RegolithLog.txt'; outfil=adjustl(outfil)
        open (ulog,file=trim(outfil),status='unknown',err=212)
        write(ulog,*) ''
        write(ulog,*) 'Edit ', trim(init), ' to correct the number of zones or &
            &maximum number of zones, then restart program.'
        close(uini); close(ulog)
        stop 'regolini.f95, lines 122 - 132'
    end if
    if(num_zones==1) then 
      max_zones=1
    end if
! Number of steps for NDSD
    read (uini,'(a)',err=202) heading(5); linct=linct+1
    read (uini,*,err=202) num_steps
    linct=linct+1
    if(trim(trans_model) == 'NDSD') then
      l_mode = .true. ! Modified mode NOT available for NDSD model
      if(num_steps<1 .or. num_steps>10000) then 
        write(*,*) ''
        write(*,*) 'Number of steps for NDSD model out of permissible range, 1 - 10,000'
        write(*,*) 'Edit ', trim(init), ' to correct the number of steps'
        outfil='RegolithLog.txt'; outfil=adjustl(outfil)
        open (ulog,file=trim(outfil),status='unknown',err=212)
        write(ulog,*) ''
        write(ulog,*) 'Number of steps for NDSD model out of permissible range, 1 - 10,000'
        write(ulog,*) 'Edit ', trim(init), ' to correct the number of steps'
        close(uini); close(ulog)
        stop 'regolini.f95, line 137 - 153'
      end if
    end if
! Allocate and initialize soil zone parameters
    allocate(depth_max(max_zones), depth_min(max_zones), zon(max_zones))
    allocate(theta_c_deg(max_zones), sc(max_zones))
    allocate(C0(max_zones), C1(max_zones), C2(max_zones))
    allocate(h0(max_zones), dif_ratio(max_zones), hump_prod(max_zones))
    h0=0.; theta_c_deg=0.; sc=0.; dif_ratio=0.; 
    depth_max=0.; depth_min=0.; C0=0.; C1=0.; C2=0.
    hump_prod=.false.; zon=0
! Soil depth zone parameters
! (NOTE, USING iz ALLOWS ZONE NUMBERS TO BE OUT OF ORDER, a maximum zone number needed to skip zones, would require storing empty zones used in larger area, but not present current extent)
  select case(model_type)
  case('EMPR') 
    do i=1,num_zones
      out_of_range = .false.
      read (uini,*,err=202) scratch,iz ! property zone number 
      if (iz < 1 .or. iz > max_zones) then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(1) = 'Zone number out of range at line '//trim(lin_num)//' of '//trim(init)
      else
        zon(i) = iz
        msg(1) = ''
      endif 
      linct=linct+1
      read (uini,'(a)',err=202) heading(6); linct=linct+1
      read (uini,*,err=202) theta_c_deg(iz), depth_min(iz), depth_max(iz),&
          & C0(iz), C1(iz), C2(iz)
  !    Test input parameters
      if(theta_c_deg(iz) < 0. .or. theta_c_deg(iz) > 90.)then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(2)='Critical slope angle, theta_c_deg, out of range at line '//trim(lin_num)//' of '//trim(init)
      else
        msg(2) = ''
      end if
      if(depth_max (iz) < 0. .or. depth_min(iz) < 0. .or. depth_min(iz) > depth_max(iz)) then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(3)='Maximum or minimum depth out of range at line '//trim(lin_num)//' of '//trim(init)
      else
        msg(3) = ''
      end if
      if(C0(iz) < 0. .or. C1(iz) < 0. .or. C2(iz) < 0.)then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(4)='Empirical parameter C0, C1, or C2 < 0 at '//trim(lin_num)//' of '//trim(init)
      else
        msg(4) = ''
      end if
      if(out_of_range)then
        outfil='RegolithLog.txt'; outfil=adjustl(outfil)
        open (ulog,file=trim(outfil),status='unknown',err=212)
        write(*,*) ''
        write(ulog,*) ''
        do j=1,4
          write(*,*) trim(adjustl(msg(j)))
          write(ulog,*) trim(adjustl(msg(j)))
        end do
        write(*,*) 'Edit ',trim(init), ' to correct these errors, then restart program.'
        write(ulog,*) 'Edit ',trim(init), ' to correct these errors, then restart program.'
        close(uini); close(ulog)
        stop 'regolini.f95, line 167 - 217'
      end if
      linct=linct+1
      sc(iz)=tan(theta_c_deg(iz)*dg2rad) ! convert angle of stability to radians and compute tangent
    end do
  case('PROC')
    do i=1,num_zones
      out_of_range = .false.
      read (uini,*,err=202) scratch,iz ! property zone number 
      if (iz < 1 .or. iz > max_zones) then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(1)='Zone number out of range at line '//trim(lin_num)//' of '//trim(init)
      else
        zon(i) = iz
        msg(1) = ''
      endif 
      linct=linct+1
      read (uini,'(a)',err=202) heading(7); linct=linct+1
      read (uini,*,err=202) theta_c_deg(iz), depth_min(iz), depth_max(iz),&
          & h0(iz), dif_ratio(iz), hump_prod(iz) 
  !    Test input parameters
      if(theta_c_deg(iz) < 0. .or. theta_c_deg(iz) > 90.)then
        out_of_range = .true.
        write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
        msg(2)='Critical slope angle, theta_c_deg, out of range at line '//trim(lin_num)//' of '//trim(init)
      else
        msg(2) = ''
      end if
      if (l_mode) then ! Allow negative minimum depth in Original mode
        if(depth_max (iz) < 0. .or. depth_min(iz) > depth_max(iz))then
          out_of_range = .true.
          write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
          msg(3)='Maximum or minimum depth out of range at line '//trim(lin_num)//' of '//trim(init)
        else
          msg(3) = ''
        end if
        if(h0(iz) < 0.)then
          out_of_range = .true.
          write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
          msg(4)='Model parameter h0 < 0 at '//trim(lin_num)//' of '//trim(init)
        else
          msg(4) = ''
        end if
      else
        if(depth_max (iz) < 0. .or. depth_min(iz) < 0. .or. depth_min(iz) > depth_max(iz))then
          out_of_range = .true.
          write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
          msg(3)='Maximum or minimum depth out of range at line '//trim(lin_num)//' of '//trim(init)
        else
          msg(3) = ''
        end if
        if(h0(iz) < 0. .or. dif_ratio(iz) < 0.)then
          out_of_range = .true.
          write(lin_num,'(i4)') linct; lin_num = adjustl(lin_num)
          msg(4)='Model parameter h0, or dif_ratio < 0 at '//trim(lin_num)//' of '//trim(init)
        else
          msg(4) = ''
        end if
      endif
      if(out_of_range)then
        outfil='RegolithLog.txt'; outfil=adjustl(outfil)
        open (ulog,file=trim(outfil),status='unknown',err=212)
        write(*,*) ''
        write(ulog,*) ''
        do j=1,4
          write(*,*) trim(adjustl(msg(j)))
          write(ulog,*) trim(adjustl(msg(j)))
        end do
        write(*,*) 'Edit ',trim(init), ' to correct these errors, then restart program.'
        write(ulog,*) 'Edit ',trim(init), ' to correct these errors, then restart program.'
        close(uini); close(ulog)
        stop 'regolini.f95, line 222 - 289'
      end if
      linct=linct+1
      sc(iz)=tan(theta_c_deg(iz)*dg2rad) ! convert angle of stability to radians and compute tangent
    end do
  case default ! undefined model type
    write(*,*) 'Model type is undefined.'
    stop 'regolini.f95, line 296'  
  end select
! Model exponent
    read (uini,'(a)',err=202) heading(8); linct=linct+1
    read (uini,*,err=202) power; linct=linct+1
    if(abs(power) > 10.) then
      write(*,*) 'Absolute value of exponent too big, |power| > 10.0: ', power
      write(*,*) 'Revise value of <power> in rg_in.txt and restart'
      write(ulog,*) 'Absolute value of exponent too big, |power| > 10.0: ', power
      write(ulog,*) 'Revise value of <power> in rg_in.txt and restart'
      close(uini); close(ulog) 
      stop 'regolini.f95, lines 299 - 307'  
    endif
    if(abs(power) < 0.1) then
      write(*,*) 'Absolute value of exponent too small, |power| < 0.1: ', power
      write(*,*) 'Revise value of <power> in rg_in.txt and restart'
      write(ulog,*) 'Absolute value of exponent too small, |power| < 0.1: ', power
      write(ulog,*) 'Revise value of <power> in rg_in.txt and restart'
      close(uini); close(ulog) 
      stop 'regolini.f95, lines 309 - 315'  
    endif
!  path names of input files
    read (uini,'(a)',err=202) heading(9); linct=linct+1
! File name of digital elevation grid (elevfil)
    read (uini,'(a)',err=202) elevfil; linct=linct+1
! File name of slope grid (slopefil)
    read (uini,'(a)',err=202) heading(10); linct=linct+1
    read (uini,'(a)',err=202) slopefil; linct=linct+1
! File name of flow-accumulation grid (flo_accfil)
    read (uini,'(a)',err=202) heading(11); linct=linct+1
    read (uini,'(a)',err=202) flo_accfil; linct=linct+1
! File name of elevation index grid file (ndxfil)
    read (uini,'(a)',err=202) heading(12); linct=linct+1
    read (uini,'(a)',err=202) ndxfil; linct=linct+1
! File name of plan-view curvature grid (pv_curvfil)
    read (uini,'(a)',err=202) heading(13); linct=linct+1
    read (uini,'(a)',err=202) pv_curvfil; linct=linct+1
! File name of soil zone grid (zonfil)
    read (uini,'(a)',err=202) heading(14); linct=linct+1
    read (uini,'(a)',err=202) zonfil; linct=linct+1
! Folder (directory) where output grid files will be stored  (folder)
    read (uini,'(a)',err=202) heading(15); linct=linct+1
    read (uini,'(a)',err=202) folder; linct=linct+1  
! Identification code to be added to names of output files (suffix)
    read (uini,'(a)',err=202) heading(16); linct=linct+1
    read (uini,'(a)',err=202) suffix; linct=linct+1 
! Specify grid file extension (.asc if true, default is .txt), Output derivatives?
    read (uini,'(a)',err=202) heading(17); linct=linct+1
    read (uini,*,err=202) lasc, l_deriv; linct=linct+1  
! Specify topographic smoothing, soil smoothing, width of running average
    read (uini,'(a)',err=202) heading(18); linct=linct+1
    read (uini,*,err=202) topoSmooth, soilSmooth, n_points; linct=linct+1  
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
! Open log file in folder occupied by elevfil and write copy of initialization file.
  ans = .false.
  patlen=scan(elevfil,'/\',.true.) ! find end of folder name
  elfoldr=elevfil(1:patlen) ! path to elevation grid
  elfoldr=adjustl(elfoldr)
  inquire(file=trim(adjustl(elevfil)), exist=ans)
  if(ans) then
    suffix=adjustl(suffix)
    trans_model=adjustl(trans_model)
!    outfil=trim(folder)//'RG_Log_'//trim(trans_model)//'_'//trim(suffix)//'.txt'
    if(trans_model(1:1)=='N' .or. trans_model(1:2)=='LC') then
      if(l_mode) then
        model_dscrp=trim(trans_model)//'_anl' ! Original mode based on analytic formula
      else
        model_dscrp=trim(trans_model)//'_mdf' ! Modified mode
      endif
      outfil=trim(elfoldr)//'RG_Log_'//trim(model_dscrp)//'_'//trim(suffix)//'.txt'
    else
      outfil=trim(elfoldr)//'RG_Log_'//trim(trans_model)//'_'//trim(suffix)//'.txt'
    end if
    outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=211)
  else
    write (*,*) ''
    write(*,*) '*###* Output directory, "', trim(elfoldr), '" was not found *###*'
    write(*,*) 'Create directory ', trim(elfoldr)
    write(*,*) 'or edit ', trim(init), ' to correct the directory path name'
    write(*,*) '(input variable "elevfil"), then restart program.'
    stop 'regolini.f95, line 359 - 386'
  end if
  write (ulog,*) ''
  write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
  write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
  write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  write (ulog,*) 'initialization file --> ', trim(init)
  write (ulog,*) '-- LISTING OF INITIALIZATION FILE --'  
! write copy of data to log file
  heading(1)=adjustl(heading(1))
  title=adjustl(title)
  write (ulog,*) trim(heading(1))
  write (ulog,*) trim(title)
  write (*,*) trim(title)
  heading(2)=adjustl(heading(2))
  trans_model=adjustl(trans_model)
  write (ulog,*) trim(heading(2))
  write (ulog,*) trim(trans_model), l_mode
  heading(3)=adjustl(heading(3))
  write (ulog,*) trim(heading(3))
  write (ulog,*) chan_thresh,chan_depth
  heading(4)=adjustl(heading(4))
  write (ulog,*) trim(heading(4))
  write (ulog,*) num_zones, max_zones
  heading(5)=adjustl(heading(5))
  write (ulog,*) trim(heading(5))
  write (ulog,*) num_steps
  select case(model_type)
  case('EMPR')
    heading(6)=adjustl(heading(6))
    do i=1,num_zones
      write (ulog,*) trim(scratch),': ',zon(i)
      write (ulog,*) trim(heading(6))
      write (ulog,*) theta_c_deg(zon(i)), depth_min(zon(i)), depth_max(zon(i)),&
        & C0(zon(i)), C1(zon(i)), C2(zon(i))
    end do
  case('PROC')
    heading(7)=adjustl(heading(7))
    do i=1,num_zones
      write (ulog,*) trim(scratch),': ',zon(i)
      write (ulog,*) trim(heading(7))
      write (ulog,*) theta_c_deg(zon(i)), depth_min(zon(i)), depth_max(zon(i)),&
        & h0(zon(i)), dif_ratio(zon(i)), hump_prod(zon(i))
    end do
  case default
    continue
  end select
  heading(8)=adjustl(heading(8))
  write (ulog,*) trim(heading(8))
  write (ulog,*) power 
  heading(9)=adjustl(heading(9))
  elevfil=adjustl(elevfil)
  write (ulog,*) trim(heading(9))
  write (ulog,*) trim(elevfil)
  heading(10)=adjustl(heading(10))
  slopefil=adjustl(slopefil)
  write (ulog,*) trim(heading(10))
  write (ulog,*) trim(slopefil)
  heading(11)=adjustl(heading(11))
  flo_accfil=adjustl(flo_accfil)
  write (ulog,*) trim(heading(11))
  write (ulog,*) trim(flo_accfil)
  heading(12)=adjustl(heading(12))
  ndxfil=adjustl(ndxfil)
  write (ulog,*) trim(heading(12))
  write (ulog,*) trim(ndxfil)
  heading(13)=adjustl(heading(13))
  pv_curvfil=adjustl(pv_curvfil)
  write (ulog,*) trim(heading(13))
  write (ulog,*) trim(pv_curvfil)
  heading(14)=adjustl(heading(14))
  zonfil=adjustl(zonfil)
  write (ulog,*) trim(heading(14))
  write (ulog,*) trim(zonfil)
!  location of output files  
  heading(15)=adjustl(heading(15))
  folder=adjustl(folder)
  write (ulog,*) trim(heading(15))
  write (ulog,*) trim(folder)
!  output-file ID code
  heading(16)=adjustl(heading(16))
  suffix=adjustl(suffix)
  write (ulog,*) trim(heading(16))
  write (ulog,*) trim(suffix)
  heading(17)=adjustl(heading(17))
  write (ulog,*) trim(heading(17))
  write (ulog,*) lasc, l_deriv
  heading(18)=adjustl(heading(18))
  write (ulog,*) trim(heading(18))
  write (ulog,*) topoSmooth, soilSmooth, n_points
  if (reset_flag) write(ulog,*) 'Incorrect value of running average <n_points> was reset to: ', n_points
  write (ulog,*) '-- END OF INITIALIZATION DATA --'  
  write (ulog,*) ''
  write (ulog,*) trim(title)
  write (ulog,*) ''
  return
  201  continue
! Report error opening initialization file.
    write (*,*) ''
    write (*,*) '*** Error opening intialization file ***'
    write (*,*) '--> ',trim(init)
    write (*,*) 'Check file location and name'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=212)
    write (ulog,*) ''
    write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
    write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
    write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
    write (ulog,*) '*** Error opening intialization file ***'
    write (ulog,*) '--> ',trim(init)
    write (ulog,*) 'Check file location and name'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '201 in regolini.f95 (line 482 - 499)'
  202  continue
! Report error reading initialization file.
    write (*,*) ''
    write (*,*) '*** Error reading initialization file ***'
    write (*,*) '--> ',trim(init), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=212)
    write (ulog,*) ''
    write (ulog,*) 'Starting Regolith ', vrsn,' ',bldate
    write (ulog,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
    write (ulog,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
    write (ulog,*) 'Error reading initialization file'
    write (ulog,*) '--> ',trim(init), 'at line ',linct
    write (ulog,*) 'Check file contents and organization'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '202 in regolini.f95 (line 500 - 517)'
  211  continue
    write (*,*) ''
    write (*,*) '*** Error opening output file in subroutine regolini() ***'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file path and directory status'
    outfil='RegolithLog.txt'; outfil=adjustl(outfil)
    open (ulog,file=trim(outfil),status='unknown',err=212)
    write (ulog,*) '*** Error opening output file in subroutine regolini() ***'
    write (ulog,*) '--> ',outfil
    write (ulog,*) 'Check file path and directory status'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '211 in regolini.f95 (line 518 - 530)'
  212  continue
    write (*,*) ''
    write (*,*) '*** Error opening output file in subroutine regolini() ***'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file path and status'
    write(*,*) 'Press RETURN to exit'
    read*
  stop '212 in regolini.f95 (line 531 - 538)'
  end subroutine regolini
  
end module read_inputs

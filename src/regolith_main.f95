! program to compute soil depth using various topographically based models
! by Rex L. Baum, U.S. Geological Survey, January 2015 - June 2021
! VARIABLE DEFINITIONS
! celsiz -- grid cell spacing/width
! chan_thresh -- threshold value of upslope contributing area for channels
! chan_depth -- mean depth of alluvial deposits in steep channels (channel gradient > sc/5.)
! col -- number of columns in grid
! dif_ratio -- diffusivity ratio, (ro_b*P0)/(ro_s*D), where D is soil diffusivity, 
!              ro_b and ro_s are density of bedrock and soil, respectively, and 
!              P0 is maximum bedrock lowering rate on a flat surface.
! grd -- total number of grid cells (col*row)
! imax --  number of (non-null) grid cells
! h0 -- characteristic soil depth, typically 0.5 m.
! h1 -- initial estimate of depth, h, in iterative solutions
! hump_prod -- logical variable true if humped soil exponential production function used, 
!              false if ordinary exponential soil production function is used.
! init -- name of initialization file (default name is "rg_in.txt")
! no_data_64 -- no-data value (64-bit/ 8-byte real) in ASCII grids
! no_data_32 -- no-data value (32-bit/ 4-byte real) in ASCII grids
! no_data_int -- no-data value (integer) in ASCII grids
! num_zones -- number of property zones
! param(6) parami(6) -- grid parameters (8-byte real) from headings of ascii grid files from floating point and integer grids, respectively
! row -- number of rows in ascii grids
! sc -- tangent of angle of stability 
! theta_c_deg -- angle of stability (degrees)
! theta_c_rad -- angle of stability (radians)
! ti -- the smallest value that can be represented by a double precision variable.
! tis -- the smallest value that can be represented by a single precision variable.
! trans_model-- sediment transport model (or empirical soil depth model) used
!
! contrib_area(:) -- upslope contributing area per unit area (ESRI "Flow Accumulation" from ArcGIS).
!                   multiply by celsiz^2 to obtain total upslope contributing area.
! cell_row(:) -- row number of each grid cell--reference from upper left corner
! cell_column(:) -- column number of each grid cell--reference from upper left corner
! cn(:,:) -- data cell numbers referenced by grid row and column--reference from lower left corner
! cos_theta(:) -- cos(theta)=sqrt(1+|Del(z)|^2)
! cta(:,:) -- data cell numbers referenced by grid row and column--reference from upper left corner
! del2gs(:) -- laplacian of ground surface
! dzdxgs(:) -- slope of ground surface in x coordinate direction
! dzdygs(:) -- slope of ground surface in y coordinate direction
! d2zdx2gs(:) -- 2nd derivative of ground surface in x corrdinate direction 
! d2zdy2gs(:) -- 2nd derivative of ground surface in y corrdinate direction
! mag_del_z -- magnitude of gradient of z, |Del(z)|, (tan(theta)), where theta is slope angle of ground surface)
! mag_del_z_sq -- squared magnitude of gradient of z, |Del(z)|^2, (also tan^2(theta))
! nl_slope_fac -- 1-(|Del(z)|/Sc)^2, where Sc is the tangent of the "angle of stability"
! Del_dotDelZ_nlso(:) -- Divergence of the gradient of z divided by the nonlinear slope factor (nl_slope_fac), Del.(Del z)/nl_slope_fac)  
! unused(:) -- placeholder for output not needed in further computations
! elev(:) -- digital elevations of ground surface
! elev0(:,:) -- digital elevations of ground surface, with no-data values set to zero.
! elev_index_lkup(:) -- index that ranks elevations at grid cells from highest to lowest.
! filtered(:) -- smoothed array output by subroutine gauss_approx()
! indexed_cell_number(:) -- cell number, indexed from highest to lowest
! n_points -- an odd integer, greater than 1, that determines the width of the running average used in the smoothing algorithm, gauss_approx
! pf1(:) -- grid stored as 1-d array, used to track locations of no-data values
! power -- exponent of DRS2 polynomial, or upslope area in NASD, NSDA, and WNDX models
! soil_depth -- computed depth of soil, h
! sec_theta -- 1/cos(slope angle) = sqrt(1+ |Del(z)|^2)
! temp -- temporary array to hold one row (line) of grid data (floating point)
! itemp -- temporary array to hold one row (line) of grid data (integer)
! trans_x(:) = dzdxgs(i)/nl_slope_fac(i)
! trans_y(:) = dzdygs(i)/nl_slope_fac(i)
! zo(:) -- property zone grid
! zon(:) -- zone ID number
!
! temp1, temp2 -- Arrays to hold intermediate values during smoothing.
!
program regolith
  use read_inputs
  implicit none
  integer, parameter:: ulen=30 ! length of array u(:) for file unit numbers
  integer:: imax,row,col, chek_row, chek_col, chek_celsiz, chek_imax, zo_min
  integer grd,i,j,j1,mnd,imx1,nwf !,patlen
  integer:: no_data_int,sctr,num_steps, num_zones, max_zones, n_points ! Added n_points 8/20/2020 RLB
  integer:: ncol,nrow,u(ulen)
  integer,allocatable:: cta(:,:),zo(:),itemp(:),pf2(:)
  integer, allocatable:: nxt(:),cell_row(:),cell_column(:)
  integer,allocatable:: dsctr(:),zon(:)
  integer, allocatable:: indexed_cell_number(:),elev_index_lkup(:) 
  real:: x1,dipdr,dip,slpdr,chan_thresh,chan_depth, tis 
  real, allocatable:: h0(:),sc(:),dif_ratio(:),depth_max(:),depth_min(:),C0(:),C1(:),C2(:)
  real, allocatable:: theta_c_rad(:), theta_c_deg(:)
  real, allocatable:: elev(:),slope(:),soil_depth(:) 
  real, allocatable:: plan_view_curv(:),tfg(:),temp(:),pf1(:)
  real, allocatable:: dzdxgs(:),dzdygs(:),d2zdx2gs(:),d2zdy2gs(:),del2gs(:) 
  real, allocatable:: mag_del_z(:), mag_del_z_sq(:),sec_theta(:),nl_slope_fac(:),slope_rad(:)
  real, allocatable:: slopgs(:),dipgs(:),contrib_area(:)
  real, allocatable:: unused(:),filtered(:),elev0(:,:),depth0(:,:),temp1(:,:),temp2(:,:)
  real, allocatable:: trans_x(:),trans_y(:),d_trans_x_dx(:),d_trans_y_dy(:)
  real, allocatable:: Del_dotDelZ_nlso(:), aspect_gs(:)
  real:: no_data_32
  real::dzdxgs_max,dzdygs_max,power
  real(kind = 8)::pi,param(6),parami(6),ti,dg2rad
  real (kind = 8):: no_data_64,celsiz, chek_nodat
  logical::lpvc,ans,lasc,lnfil 
  logical:: topoSmooth,soilSmooth,hump_prod_any,l_deriv, l_mode
  logical, allocatable:: hump_prod(:) ! Switches between exponential (.false.) or humped soil production model (.true.)
  character (len=255):: outfil,init !,infil
  character (len=255):: elevfil, slopefil, flo_accfil, pv_curvfil, ndxfil, zonfil
  character (len=255):: heading(18) !,size_heading
  character (len=224):: folder ,elfoldr
  character (len=31):: scratch ! Holder for intermediate output strings
  character (len=16):: suffix ! Alphanumeric code appended to output file names to identify files from a particular model run.
  character (len=14):: header(6) ! Labels from the first six lines of ASCII grid files.
  character (len=14):: elevSmo_fil='RG_elev_Smooth'
  character (len=11):: bldate ! date of current version 
  character (len=10):: time
  character (len=8):: date
  character (len=8):: slopgs_fil='RG_slope' ! output file for computed ground-surface slope angle grid
  character (len=7):: vrsn ! version number
  character (len=4):: trans_model,grxt
!!!  character (len=2):: pid(3)
  character (len=1):: tb ! tab character
! first executable statement ............  
  call date_and_time(date,time)
  no_data_32=-9999.
  u=(/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,&
    &31,32,33,34,35,36,37,38,39,40/) ! Unit numbers for input and output files.
!!!  pid=(/'TI','GM','TR'/)
  pi=3.141592653589793
  dg2rad=pi/180.D0
  vrsn='0.9.02b'; bldate='07 Jun 2021'
  mnd=6 ! Default value assumed if no integer grid is read.
  tb=char(9)
  call rgbanner(vrsn,bldate)
!  read initialization file
!! ********** Change zoned input parameters from dif_ratio to dif_ratio_list, etc *****************
  call regolini(u(1),u(2),init,dg2rad,trans_model,chan_thresh,chan_depth,&
     & num_zones,max_zones,num_steps,hump_prod,theta_c_deg,suffix,folder,&
     & elevfil,slopefil,flo_accfil,pv_curvfil,ndxfil,zonfil,heading,lasc,&
     & h0,sc,dif_ratio,depth_max,depth_min,C0,C1,C2,zon,vrsn,bldate,date,time,&
     & outfil,topoSmooth,soilSmooth,n_points,l_deriv,l_mode,power)
  allocate(theta_c_rad(max_zones))
  theta_c_rad=theta_c_deg*dg2rad
  write(*,*) 'chan_thresh: ', chan_thresh
! determine grid size parameters RLB 4/18/2011
  call grid_size(elevfil,elfoldr,init, imax,row,col,nwf,u(1),u(3),u(12),&
   & celsiz,no_data_64)
! Allocate & initialize arrays 
  grd=row*col
  imx1=imax
  allocate(pf1(grd),pf2(grd),tfg(imax),nxt(imax))
  allocate(dsctr(imax+1),slope(imax),slope_rad(imax))
  allocate(temp(col),itemp(col))
  allocate(soil_depth(imax),zo(imax),contrib_area(imax))
  allocate(elev(imax)) 
  allocate(dzdxgs(imax),dzdygs(imax),dipgs(imax),slopgs(imax))
  allocate(cta(col,row),cell_row(imax),cell_column(imax))
  if(l_deriv .and. trim(slopefil)=='none' .or. topoSmooth) then
    allocate(aspect_gs(imax))
    aspect_gs=0.
  end if
  if(topoSmooth .or. soilSmooth) then
    allocate(filtered(imax))
    filtered=0.
  endif
  if(trans_model(1:3)=='DRS')then
    allocate(plan_view_curv(imax))
    plan_view_curv=0.
  end if
  pf1=0.;pf2=0; tfg=0.  !
  nxt=0; dsctr=0; zo=1
  slope=0.; slope_rad=0. 
  temp=0.;itemp=0 
  soil_depth=0.; contrib_area=0.; elev=0.;  
  dzdxgs=0.;dzdygs=0.;dipgs=0.;slopgs=0.
! Choose file extension for grid files, Added 4/14/2010
  grxt='.txt'
  if(lasc) grxt='.asc'
! *****************************************************************
!  read gridded data from GIS
  write (*,*) 'Reading input grids'
  write(u(1),*) 'Input file name,            Cell count'
!  read digital elevations, elev 
  inquire(file=elevfil,exist=lnfil)
  write(*,*) 'Status of elevfil:', lnfil, trim(elevfil)
  if(lnfil) then
    call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(12),elevfil,header,u(1))
    if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
      write(u(1),*) '***###*** Grid-size parameters do not match ***###***'
      write(u(1),*) 'Delete Grid-size file (TIgrid_size.txt or TRgrid_size.txt)&
        & from directory containng elevation grid, then restart program.'
      write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
      write(*,*) '***###*** Grid-size parameters do not match ***###***'
      write(*,*) 'Delete Grid-size file (TIgrid_size.txt or TRgrid_size.txt)&
          & from directory containng elevation grid, then restart program.'
      write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
      close(u(1))
      stop ' regolith_main.f95, line 174'
    end if
    call srdgrd(grd,col,ncol,nrow,celsiz,no_data_64,&
     &    elev,pf1,sctr,imax,temp,u(4),elevfil,param,header,u(1))
    write(u(1),*) 'Elevation grid'
    write(u(1),*) trim(elevfil),sctr,' data cells'
  else
    write(u(1),*) 'Elevation grid was not found'
    write(u(1),*) 'Edit path name of elevation grid in ', trim(init), ' then restart program.'
    write(*,*) 'Elevation grid was not found'
    write(*,*) 'Edit path name of elevation grid in ', trim(init), ' then restart program.'
    close(u(1))
    stop ' regolith_main.f95, line 174'
  end if
! *****************************************************************
  write(u(1),*) '---------------******---------------'
! Create 2-d array, cta(j,i) that maps grid cell numbers (1-d array) to j-i (column-row)coordinates
  no_data_int=int(no_data_64) ! null or no-data value for integer grids and arrays
  cta=no_data_int
  call grid_count(ncol,nrow,imax,no_data_int,no_data_64,cell_row,cell_column,cta,pf1) 
! (Optional) Smooth Elevation grid
  if(topoSmooth) then
    allocate(elev0(col,row),temp1(col,row),temp2(col,row))
    elev0=0. ! initialize elev0(:,:) and copy values from elev(:)
    temp1 = 0.; temp2 = 0.
    do i=1,nrow
      do j=1,ncol
        if(cta(j,i)/=no_data_int) elev0(j,i)=elev(cta(j,i))
      end do
    end do  
    call gauss_approx(elev0,cta,imax,ncol,nrow,celsiz,celsiz,no_data_64,&
                      &no_data_int,filtered,temp1,temp2,n_points)  ! last parameter is width of moving average window.
    elev=filtered
    scratch=trim(elevSmo_fil)
    scratch=adjustl(scratch)
    outfil=trim(folder)//trim(scratch)//grxt
    call ssvgrd(filtered,imax,pf1,row,col,u(10),no_data_32,param,u(1),&
      & outfil,ti,header)
    write(u(1),*) 'Applied optional smoothing to soil depth grid'
    deallocate(elev0,temp1,temp2)
  endif
!  read slope angles
  dzdxgs_max=0.; dzdygs_max=0.
  if(trim(slopefil)=='none' .or. topoSmooth)then
! compute east-west and north-south slope gradients
    write(*,*) 'Computing E-W & N-S elevation gradients'
    call xyslope(elev,pf1,cta,imax,ncol,nrow,dzdxgs,dzdygs,celsiz,celsiz,no_data_64,no_data_int)
    dzdxgs_max=maxval(dzdxgs); dzdygs_max=maxval(dzdygs)
    write(*,*) 'Max X and Y slopes:', dzdxgs_max, dzdygs_max
    if(l_deriv .and. trim(slopefil)=='none' .or. topoSmooth) then
      do i=1,imax
        call aspect(dzdxgs(i),dzdygs(i),dipdr,dip,slpdr)
        dipgs(i)=dip ! Angle of steepest descent in radians
        slopgs(i)=dip*180./pi  ! Angle of steepest descent in degrees
        if(l_deriv .and. trim(slopefil)=='none' .or. topoSmooth) aspect_gs(i) = dipdr*180./pi
      end do 
    else
      do i=1,imax
        call aspect(dzdxgs(i),dzdygs(i),dipdr,dip,slpdr)
        dipgs(i)=dip ! Angle of steepest descent in radians
        slopgs(i)=dip*180./pi  ! Angle of steepest descent in degrees
      end do 
    end if
    slope = slopgs ! Fix discrepancy between slope and slopgs.  
    slope_rad = dipgs 
! Save ground-surface slope file 
    scratch=trim(slopgs_fil)
    if(topoSmooth) scratch=trim(slopgs_fil)//'_smo'
    scratch=adjustl(scratch)
    outfil=trim(folder)//trim(scratch)//grxt
    call ssvgrd(slopgs,imax,pf1,row,col,u(13),no_data_32,param,u(1),&
      & outfil,ti,header)
  else
!  FOR THIS AND OTHER CASES WHERE INPUT GRIDS ARE IMPORTED, USE ssizgrd TO CHECK SIZE MATCHING BEFORE ATTEMPTING TO READ THE GRID FILE.
    inquire(file=slopefil,exist=lnfil)
    write(*,*) 'Status of slopefil:', lnfil, trim(slopefil)
    if(lnfil) then
      call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(30),slopefil,header,u(1))
      if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
        write(*,*) '***###*** Grid mismatch: "', trim(slopefil), '" ***###***'
        write(*,*) 'Compare numbers of rows, columns and no-data cells'
        write(*,*) 'of slope grid and elevation grid.'
        write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
        write(u(1),*) '***###*** Grid mismatch: "', trim(slopefil), '" ***###***'
        write(u(1),*) 'Compare numbers of rows, columns and no-data cells'
        write(u(1),*) 'of slope grid and elevation grid.'
        write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
        close(u(1))
        stop 'regolith_main, line 264'
      end if
      call srdgrd(grd,col,ncol,nrow,celsiz,no_data_64,&
       &  slope,pf1,sctr,imax,temp,u(5),slopefil,param,header,u(1))
      write(u(1),*) 'slope angle grid'
      write(u(1),*) trim(slopefil),sctr,' data cells'
      slope_rad=slope*dg2rad ! convert slope angles to radians
    else
      write(u(1),*) 'Slope grid ', trim(slopefil), ' was not found'
      write(u(1),*) 'Edit path name in ', trim(init), ', then restart program.'
      write(*,*) 'Slope grid ', trim(slopefil), ' was not found'
      write(*,*)  'Edit path name in ', trim(init), ', then restart program.'
      close(u(1))
      stop  'regolith_main, line 264'
    endif
  endif   
!  read flow-accumulation grid (upslope contributing area per unit grid cell)
  inquire(file=flo_accfil,exist=lnfil)
  write(*,*) 'Status of flo_accfil: ', lnfil, ' ', trim(flo_accfil)
  if(lnfil) then
    call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(30),flo_accfil,header,u(1))
    if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
      write(*,*) 'Grid mismatch: "', trim(flo_accfil), '" ***###***'
      write(*,*) 'Compare numbers of rows, columns and no-data cells'
      write(*,*) 'of flow accumulation grid against elevation grid.'
      write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
      write(u(1),*) 'Grid mismatch: "', trim(flo_accfil), '" ***###***'
      write(u(1),*) 'Compare numbers of rows, columns and no-data cells'
      write(u(1),*) 'of flow accumulation grid and elevation grid.'
      write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
      close(u(1))
      stop 'regolith_main, line 297'
    end if
    call srdgrd(grd,col,ncol,nrow,celsiz,no_data_64,&
    & contrib_area,pf1,sctr,imax,temp,u(6),flo_accfil,param,header,u(1))
    write(u(1),*) 'Flow accumulation grid'
    write(u(1),*) trim(flo_accfil),sctr,' data cells'
    write(*,*) 'chan_thresh, celsiz: ', chan_thresh, celsiz
    chan_thresh=chan_thresh/(celsiz*celsiz) ! normalize by grid cell area.
! Test for models that require upslope contributing area 
  elseif (trans_model == 'WNDX' .or. trans_model == 'NASD' .or. trans_model &
          &== 'NSDA') then 
    write(u(1),*) 'Flow accumulation grid ', trim(flo_accfil), ' was not found'
    write(u(1),*) 'Edit path name in ', trim(init), ', then restart program.'
    write(*,*) 'Flow accumulation grid ', trim(flo_accfil), ' was not found'
    write(*,*)  'Edit path name in ', trim(init), ', then restart program.'
    close(u(1))
    stop 'regolith_main, line 297'
  else ! *** Exception here models that do not require contrbuting area. ***
    write(u(1),*) 'Flow accumulation grid ', trim(flo_accfil), ' was not found'
    write(u(1),*) 'Regolith will compute soil depths without adjusting for channel depth.'
    write(*,*) 'Flow accumulation grid ', trim(flo_accfil), ' was not found'
    write(*,*)  'Regolith will compute soil depths without adjusting for channel depth.'
  end if
! read cell number index to determine order of computation for NDSD model
  if(trans_model(1:3) == 'NDS') then
    allocate(elev_index_lkup(imax),indexed_cell_number(imax))
    elev_index_lkup=0; indexed_cell_number=0 
    inquire(file=ndxfil,exist=lnfil)
    write(*,*) 'Status of ndxfil:', lnfil, trim(ndxfil)
    if(lnfil) then
      call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(30),ndxfil,header,u(1))
      if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
         write(*,*) '***###*** Grid mismatch: "', trim(ndxfil), '" ***###***'
         write(*,*) 'Compare numbers of rows, columns and no-data cells'
         write(*,*) 'of cell-index grid and elevation grid.'
         write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         write(u(1),*) 'Grid mismatch: "', trim(ndxfil), '" ***###***'
         write(u(1),*) 'Compare numbers of rows, columns and no-data cells'
         write(u(1),*) 'of cell-index grid and elevation grid.'
         write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         close(u(1))
         stop 'regolith_main, line 334'
      end if
      call irdgrd(grd,col,ncol,nrow,celsiz,no_data_int,mnd,&
      &elev_index_lkup,pf2,sctr,imax,itemp,u(7),ndxfil,parami,header,u(1))
      write(u(1),*) 'Cell index grid'
      write(u(1),*) trim(ndxfil),sctr,' data cells'
      do i=1,imax
        j1=elev_index_lkup(i) ! grid cell number of the ith cell in decending elevation
        indexed_cell_number(j1)=i ! elevaton index (highest to lowest) of j1th cell 
      end do
    else
      write(u(1),*) 'Cell index grid needed by NDSD model ', trim(ndxfil), ' was not found'
      write(u(1),*) 'Edit path name in ', trim(init), ', then restart program.'
      write(*,*) 'Cell index grid needed by NDSD model ', trim(ndxfil), ' was not found'
      write(*,*)  'Edit path name in ', trim(init), ', then restart program.'
      close(u(1))
      stop 'regolith_main, line 334'
    end if
  end if
! read plan-view curvature grid 
  if(trans_model=='DRS3')then
    inquire(file=pv_curvfil,exist=lpvc)
    write(*,*) 'Status of pv_curvfil:', lpvc, trim(pv_curvfil)
    if(lpvc) then
      call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(30),pv_curvfil,header,u(1))
      if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
         write(*,*) '***###*** Grid mismatch: "', trim(pv_curvfil), '" ***###***'
         write(*,*) 'Compare numbers of rows, columns and no-data cells'
         write(*,*) 'of Plan-view curvature grid and elevation grid.'
         write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         write(u(1),*) '***###*** Grid mismatch: "', trim(pv_curvfil), '" ***###***'
         write(u(1),*) 'Compare numbers of rows, columns and no-data cells'
         write(u(1),*) 'of Plan-view curvature grid and elevation grid.'
         write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         close(u(1))
         stop 'regolith_main, line 373'
      end if
      call srdgrd(grd,col,ncol,nrow,celsiz,no_data_64,&
      & plan_view_curv,pf1,sctr,imax,temp,u(11),pv_curvfil,param,header,u(1))
      write(u(1),*) 'Plan-view curvature grid'
      write(u(1),*) trim(pv_curvfil),sctr,' data cells'
    else
      write(*,*) ''
      write(u(1),*) '***###*** P-V curvature grid was not found ***###***'
      write(u(1),*) 'P-V curvature grid: "', trim(pv_curvfil), '"'
      write(u(1),*) 'Edit ', trim(init), ' to correct file name and directory location or select a different model.'
      write(*,*) ''
      write(*,*) '***###*** P-V curvature grid was not found ***###***'
      write(*,*) 'P-V curvature grid: "', trim(pv_curvfil), '"'
      write(*,*) 'Edit ', trim(init), ' to correct file name and directory location or select a different model.'
      close(u(1))
      stop 'regolith_main.f95, line 373'
    end if
  end if
! Import property zone grid
  if(num_zones > 1)then
    ans=.false.
    inquire (file=trim(zonfil),exist=ans)
    write(*,*) 'Status of zonfil:', ans, trim(zonfil)
    if(ans) then
      call ssizgrd(chek_row,chek_col,chek_celsiz,chek_nodat,chek_imax,u(30),zonfil,header,u(1))
      if(chek_row/=row .or. chek_col/=col .or. chek_imax/=imax) then
         write(*,*) '***###*** Grid mismatch: "', trim(zonfil), '" ***###***'
         write(*,*) 'Compare numbers of rows, columns and no-data cells'
         write(*,*) 'of zone grid and elevation grid.'
         write(*,*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         write(u(1),*) '***###*** Grid mismatch: "', trim(zonfil), '" ***###***'
         write(u(1),*) 'Compare numbers of rows, columns and no-data cells'
         write(u(1),*) 'of zone grid and elevation grid.'
         write(u(1),*) 'chek_row, row; chek_col, col; chek_imax, imax', chek_row, row,&
          & '; ', chek_col, col, '; ', chek_imax, imax
         close(u(1))
         stop 'regolith_main, line 410'
      end if
      call irdgrd(grd,col,ncol,nrow,celsiz,no_data_int,mnd,&
             &   zo,pf2,sctr,imax,itemp,u(11),zonfil,parami,header,u(1))
      write(u(1),*) 'Property zone grid'
      write(u(1),*) trim(zonfil),sctr,' data cells'  
      zo_min = minval(zo)
      if(zo_min <=0) then
        write(u(1),*) '***###*** Negative property zone found ***###***'
        write(u(1),*) 'Edit property zone grid and restart program.'
        write(*,*) '***###*** Negative property zone found ***###***'
        write(u(1),*) 'Edit property zone grid and restart program.'
        close(u(1))
        stop 'regolith_main.f95, line 435'
      end if
    else
      write(u(1),*) '***###*** Property zone grid was not found. ***###***'
      write(u(1),*) 'zone grid: "',trim(zonfil),'"'
      write(u(1),*) 'Edit ', trim(init), ' to check zone file name and directory location.'
      write(*,*) ''
      write(*,*) '***###*** Property zone grid was not found. ***###***'
      write(*,*) 'zone grid: "',trim(zonfil),'"'
      write(*,*) 'Edit ', trim(init), ' to check zone file name and directory location.'
      close(u(1))
      stop 'regolith_main.f95, line 410'
    endif
  else
    zo=1; zon=1 ! only one property zone.
    write(*,*) 'Property zone grid not required'
    write(u(1),*) 'Property zone grid not required'
  end if
! compute gradients and related quantities
  if(trans_model(1:1)=='N' .or. trans_model(1:1)=='L') then 
    allocate(mag_del_z(imax), mag_del_z_sq(imax),sec_theta(imax),nl_slope_fac(imax))
    allocate(unused(imax))
    allocate(trans_x(imax),trans_y(imax),d_trans_x_dx(imax),d_trans_y_dy(imax))
    mag_del_z=0.;mag_del_z_sq=0.;sec_theta=0.;nl_slope_fac=0.
    trans_x=0.;trans_y=0.;d_trans_x_dx=0.;d_trans_y_dy=0.;unused=0.
! compute east-west and north-south slope gradients
    if(dzdxgs_max*dzdygs_max==0.)then
      write(*,*) 'Computing E-W & N-S elevation gradients'
      call xyslope(elev,pf1,cta,imax,ncol,nrow,dzdxgs,dzdygs,celsiz,celsiz,no_data_64,no_data_int)
    endif 
! compute magnitude and squared magnitude of elevation gradient vector, Del z,
! non-linear slope factor, 1 - (|Del (z)|/Sc)^2, and 1/cos(theta), where theta is the slope angle of the ground surface
    write(*,*) 'Computing magnitude and related quantities for elevation gradients'
    do i=1,imax
      mag_del_z_sq(i)=dzdxgs(i)*dzdxgs(i)+dzdygs(i)*dzdygs(i)
      mag_del_z(i)=sqrt(mag_del_z_sq(i))
      nl_slope_fac(i)=1.d0-mag_del_z_sq(i)/(sc(zo(i))*sc(zo(i)))
      sec_theta(i)=sqrt(1+mag_del_z_sq(i))
    end do    
    tis=tiny(x1)
  end if
! compute soil_depth according to different soil-depth models
  write(*,*) 'Selecting depth model'
  select case(trans_model)
    case('DRS1'); call derose(u(1),imax,chan_thresh,chan_depth,theta_c_rad,slope,slope_rad,&
      & dg2rad,contrib_area,plan_view_curv,soil_depth,trans_model,depth_max,&
      & depth_min,C0,C1,C2,zo,max_zones,power) ! DeRose exponential formula 
    case('DRS2'); call derose(u(1),imax,chan_thresh,chan_depth,theta_c_rad,slope,slope_rad,&
      & dg2rad,contrib_area,plan_view_curv,soil_depth,trans_model,depth_max,&
      & depth_min,C0,C1,C2,zo,max_zones,power) ! DeRose polynomial formula 
    case('DRS3'); call derose(u(1),imax,chan_thresh,chan_depth,theta_c_rad,slope,slope_rad,&
      & dg2rad,contrib_area,plan_view_curv,soil_depth,trans_model,depth_max,&
      & depth_min,C0,C1,C2,zo,max_zones,power) ! DeRose exponential formula with curvature 
    case('WNDX'); call wetness_ndx(u(1),imax,chan_thresh,chan_depth,theta_c_rad,dg2rad,&
      & contrib_area,slope_rad,soil_depth,depth_min,depth_max,C0,zo,max_zones,power) ! Modified wetness index 
    case('LRSC') ! Linear regression slope- and curvature-based formula
      if (l_deriv .or. l_mode) then
        allocate(d2zdx2gs(imax),d2zdy2gs(imax),del2gs(imax))
        d2zdx2gs=0.;d2zdy2gs=0.;del2gs=0.
        call laplacian(elev,pf1,cta,imax,col,row,d2zdx2gs,d2zdy2gs,del2gs,&
           & celsiz,celsiz,no_data_64,no_data_int)
      end if
      if(l_mode) then ! Pass 2nd derivatives from laplacian subroutine.
        call lrsc_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,&
        & chan_thresh,chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,mag_del_z,slope_rad,&
        & contrib_area,soil_depth,C0,C1,C2,depth_max,depth_min,&
        & unused,del2gs,d2zdx2gs,d2zdy2gs,zo,max_zones,l_mode)
      else ! Compute 2nd derivatives from 1st derivatives for greater smoothing.
        allocate(del2gs(imax)); del2gs=0.
        call lrsc_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,&
        & chan_thresh,chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,mag_del_z,slope_rad,&
        & contrib_area,soil_depth,C0,C1,C2,depth_max,depth_min,&
        & unused,del2gs,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode)
      endif
    case('LCSD') ! Classic curvature-based formula
      if (l_deriv .or. l_mode) then
        allocate(d2zdx2gs(imax),d2zdy2gs(imax),del2gs(imax))
        d2zdx2gs=0.;d2zdy2gs=0.;del2gs=0.
        call laplacian(elev,pf1,cta,imax,col,row,d2zdx2gs,d2zdy2gs,del2gs,&
           & celsiz,celsiz,no_data_64,no_data_int)
      end if
      if(l_mode) then ! Pass 2nd derivatives from laplacian subroutine.
        call lcsd_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,chan_thresh,&
        & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,slope_rad,&
        & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
        & unused,trans_x,trans_y,d2zdx2gs,d2zdy2gs,zo,max_zones,l_mode)
      else ! Compute 2nd derivatives from 1st derivatives for greater smoothing.
        call lcsd_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,chan_thresh,&
        & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,slope_rad,&
        & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
        & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode)
      endif
    case('NSD');  call nsd_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,chan_thresh,&
      & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,nl_slope_fac,slope_rad,&
      & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
      & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode)
    case('NSDA'); call nsd_a_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,chan_thresh,&
      & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,nl_slope_fac,slope_rad,&
      & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
      & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode,power)
    case('NASD'); call nasd_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,cta,chan_thresh,&
      & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,nl_slope_fac,slope_rad,&
      & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
      & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode,power)
    case('NDSD') 
      allocate(d2zdx2gs(imax),d2zdy2gs(imax),del2gs(imax))
      d2zdx2gs=0.;d2zdy2gs=0.;del2gs=0.
      call laplacian(elev,pf1,cta,imax,col,row,d2zdx2gs,d2zdy2gs,del2gs,&
         & celsiz,celsiz,no_data_64,no_data_int)
      call ndsd_depth(u(1),imax,col,row,grd,celsiz,no_data_64,no_data_int,&
         & cell_row,cell_column,indexed_cell_number,elev_index_lkup,cta,pf1,dzdxgs,&
         & dzdygs,del2gs,nl_slope_fac,sec_theta,soil_depth,num_steps,chan_thresh,&
         & chan_depth,contrib_area,theta_c_rad,slope_rad,hump_prod,h0,dif_ratio,&
         & depth_max,depth_min,tis,unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,&
         & zo,max_zones)
    case default
      write(*,*) ''
      write(*,*) '***###*** Invalid soil depth model selected. ***###***'
      write(*,*) 'Edit ', trim(init), ' to correct model code, then restart program'
      write(u(1),*) '***###*** Invalid soil depth model selected. ***###***'
      write(u(1),*) 'Edit ', trim(init), ' to correct model code, then restart program'
      close(u(1))
      stop 'regolith_main.f95, line 536'
  end select 
  write(u(1),*) 'Range soil depth: ', minval(soil_depth), maxval(soil_depth)
! (Optional) Apply smoothing algorithm to computed soil depth.
  if(soilSmooth)then
    allocate(depth0(col,row),temp1(col,row),temp2(col,row))
    depth0=0. ! ! initialize depth0(:,:) and copy values from soil_depth(:)
    temp1=0.; temp2=0.
    do i=1,nrow
      do j=1,ncol
        if(cta(j,i)/=no_data_int) depth0(j,i)=soil_depth(cta(j,i))
      end do
    end do  
    call gauss_approx(depth0,cta,imax,ncol,nrow,celsiz,celsiz,no_data_64,&
                      &no_data_int,filtered,temp1,temp2,n_points) ! last parameter is width of moving average window.
    soil_depth=filtered
    write(u(1),*) 'Applied optional smoothing to soil depth grid'
    write(u(1),*) 'Range smoothed soil depth: ', minval(soil_depth), maxval(soil_depth)
    deallocate(depth0,temp1,temp2)
  endif
!
! Save results    
  write(*,*) 'Saving results'
  ti=tiny(param(1)) 
  do i=1,imx1
    tfg(i)=soil_depth(i)
  end do
! Save soil depth file  
  scratch='RG_'//trim(trans_model)//"_"
  scratch=adjustl(scratch)
  hump_prod_any=.false.
  do i=1,max_zones
    if(hump_prod(i)) hump_prod_any=.true.
  end do
  if(hump_prod_any) scratch=trim(scratch)//'hmp_'
  if(soilSmooth) scratch=trim(scratch)//'smo_'
  if(trans_model(1:1)=='N' .or. trans_model(1:1)=='L') then
    if(l_mode) then
      scratch=trim(scratch)//'anl_' ! Original mode based on analytic formula
    else
      scratch=trim(scratch)//'mdf_' ! Modified mode
    endif
  end if
  outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
  call ssvgrd(tfg,imax,pf1,row,col,u(9),no_data_32,param,u(1),&
   & outfil,ti,header)
! Save series of files with derivatives or related quantities
  if(l_deriv) then
  ! save dzdxgs
      scratch = 'RG_dzdxgs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(dzdxgs,imax,pf1,row,col,u(17),no_data_32,param,u(1),&
       & outfil,ti,header)    
  ! save dzdygs
      scratch = 'RG_dzdygs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(dzdygs,imax,pf1,row,col,u(18),no_data_32,param,u(1),&
       & outfil,ti,header)    
  ! save trans_nsd
    if(trans_model(1:3)=='NSD') then
      scratch = 'RG_trans_nsd_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(unused,imax,pf1,row,col,u(19),no_data_32,param,u(1),&
       & outfil,ti,header)    
    endif
  ! save trans_nsd
    if(trans_model(1:3)=='NAS') then
      scratch = 'RG_trans_nasd_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(unused,imax,pf1,row,col,u(19),no_data_32,param,u(1),&
       & outfil,ti,header)    
    endif
  ! save d2zdx2gs
    if(trans_model(1:3)=='NDS' .or. trans_model == 'LCSD') then
      scratch = 'RG_d2zdx2gs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(d2zdx2gs,imax,pf1,row,col,u(19),no_data_32,param,u(1),&
       & outfil,ti,header)
    endif
  ! save d2zdy2gs
    if(trans_model(1:3)=='NDS' .or. trans_model == 'LCSD') then
      scratch = 'RG_d2zdy2gs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(d2zdy2gs,imax,pf1,row,col,u(20),no_data_32,param,u(1),&
       & outfil,ti,header)    
    endif
  ! save laplacian   
    if(trans_model(1:3)=='NDS' .or. trans_model == 'LCSD') then
      scratch='RG_del2gs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(del2gs,imax,pf1,row,col,u(21),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save mag_del_z_sq   
    if(trans_model(1:1)=='N' .or. trans_model(1:1)=='L') then
      scratch='RG_mag_del_z_sq_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(mag_del_z_sq,imax,pf1,row,col,u(22),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save sec_theta   
    if(trans_model(1:1)=='N' .or. trans_model(1:1)=='L') then
      scratch='RG_sec_theta_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(sec_theta,imax,pf1,row,col,u(23),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save nl_slope_fac   
    if(trans_model(1:1)=='N') then
      scratch='RG_nl_slope_fac_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(nl_slope_fac,imax,pf1,row,col,u(24),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save trans_x   
    if(trans_model(1:1)=='N') then
      scratch='RG_trans_x_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(trans_x,imax,pf1,row,col,u(25),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save trans_y 
    if(trans_model(1:1)=='N') then
      scratch='RG_trans_y_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(trans_y,imax,pf1,row,col,u(26),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save Del_dotDelZ_nlso   
    if(trans_model=='NDSD') then
      allocate(Del_dotDelZ_nlso(imax))
      Del_dotDelZ_nlso=0.
      Del_dotDelZ_nlso=d_trans_x_dx+d_trans_y_dy
      scratch='RG_Del_dotDelZ_nlso_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(Del_dotDelZ_nlso,imax,pf1,row,col,u(27),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  ! save aspect_gs  
    if(trim(slopefil)=='none' .or. topoSmooth)then  
      scratch='RG_aspect_gs_'
      scratch=adjustl(scratch)
      if(soilSmooth) scratch=trim(scratch)//'smo_'
      outfil=trim(folder)//trim(scratch)//trim(suffix)//grxt
      call ssvgrd(aspect_gs,imax,pf1,row,col,u(26),no_data_32,param,u(1),&
       & outfil,ti,header)    
    end if
  end if
! Final
  write (u(1),*) 'Regolith finished normally'
  call date_and_time(date,time)
  write (u(1),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
  write (u(1),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  close (u(1))
  stop 'Regolith finished normally!'
end program regolith

! Procedure to compute soil depth based on Linear Regression Slope and Curvature  model
! 4 June 2021, RLB, Latest revision 7 Jun 2021, RLB.
  subroutine lrsc_depth(ulog,imax,ncol,nrow,grd,celsiz,no_data_64,no_data_int,cta,&
    & chan_thresh,chan_depth,theta_c_rad,pf1,dz_gs_dx,dz_gs_dy,mag_del_z,slope_rad,&
    & contrib_area,soil_depth,C0,C1,C2,depth_max,depth_min,&
    & unused,Laplacian,d2z_gs_dx2,d2z_gs_dy2,zo,max_zones,l_mode)
    implicit none
! LOCAL VARIABLES
    integer:: i, chan_ctr
    real ::soil_depth_min, soil_depth_max
    real::h1, h2, sc
! FORMAL ARGUMENTS
    integer, intent(in)::ulog,imax,ncol,nrow,grd,no_data_int,cta(ncol,nrow),max_zones,zo(imax)
    real, intent(in):: C0(max_zones), C1(max_zones), C2(max_zones)
    real, intent(in):: depth_max(max_zones),depth_min(max_zones)
    real, intent(in)::chan_thresh,chan_depth,theta_c_rad(max_zones),pf1(grd),contrib_area(imax)
    real, intent(in):: dz_gs_dx(imax),dz_gs_dy(imax),mag_del_z(imax),slope_rad(imax)
    real, intent(inout)::soil_depth(imax), unused(imax)
    real, intent(inout):: d2z_gs_dx2(imax),d2z_gs_dy2(imax),Laplacian(imax)
    real (kind = 8),intent(in):: no_data_64,celsiz !
    logical, intent(in):: l_mode
    write(*,*) 'Entering subroutine lrsc_depth'
    if(l_mode) then
      continue ! 2nd derivatives passed directly from main program.
    else  ! Modified mode, recompute 2nd derivatives from 1st derivatives
      call xyslope(dz_gs_dx,pf1,cta,imax,ncol,nrow,d2z_gs_dx2,unused,celsiz,celsiz,no_data_64,no_data_int)
      call xyslope(dz_gs_dy,pf1,cta,imax,ncol,nrow,unused,d2z_gs_dy2,celsiz,celsiz,no_data_64,no_data_int)
    endif
    Laplacian = d2z_gs_dx2 + d2z_gs_dy2
    chan_ctr=0
    do i=1,imax
      sc = tan(theta_c_rad(zo(i)))
      if (mag_del_z(i) > sc) then
        h1 = 0. ! no contribution above critical slope
        h2 = 0.
      else
        h1 = C2(zo(i)) * (sc - mag_del_z(i))
        h2 = C0(zo(i)) + C1(zo(i)) * Laplacian(i) 
      end if
      soil_depth(i)= h1 + h2
      if(soil_depth(i) < depth_min(zo(i))) soil_depth(i)=depth_min(zo(i))
      if(soil_depth(i) > depth_max(zo(i)) ) soil_depth(i)=depth_max(zo(i))
! Compare slope angle in channels with 0.1*critical slope angle, and reduce thickness accordingly.
      if(contrib_area(i)>chan_thresh) then 
         if(slope_rad(i)>0.1*theta_c_rad(zo(i))) then
            if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
              chan_ctr = chan_ctr + 1
         end if
      end if
    end do
    soil_depth_min=minval(soil_depth)
    soil_depth_max=maxval(soil_depth)
    write(*,*) 'Computed depth using LRSC model'
    write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
    write(ulog,*) 'Computed depth using LRSC model'
    write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
    write(ulog,*) 'Range dz_gs_dx', minval(dz_gs_dx), maxval(dz_gs_dx)
    write(ulog,*) 'Range dz_gs_dy', minval(dz_gs_dy), maxval(dz_gs_dy)
    write(ulog,*) 'Range d2z_gs_dx2', minval(d2z_gs_dx2), maxval(d2z_gs_dx2)
    write(ulog,*) 'Range d2z_gs_dy2', minval(d2z_gs_dy2), maxval(d2z_gs_dy2)
    write(ulog,*) 'Channel grid cell count where depth changed, grid-cell threshold: ', chan_ctr, chan_thresh
    return 
  end subroutine lrsc_depth

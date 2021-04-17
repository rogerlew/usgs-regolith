! procedure to compute soil depth based on NSD transport model
! 3 Feb 2015, RLB Latest revision April 2021, RLB.
subroutine nsd_depth(ulog,imax,ncol,nrow,grd,celsiz,nodat,no_data_int,cta,chan_thresh,&
  & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,nl_slope_fac,slope_rad,&
  & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
  & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode)
  implicit none
! LOCAL VARIABLES
  integer:: i, chan_ctr
  real ::trans_nsd,soil_depth_min,soil_depth_max 
  real::h1 
! FORMAL ARGUMENTS
  integer, intent(in)::ulog,imax,ncol,nrow,grd,no_data_int,cta(ncol,nrow),max_zones,zo(imax)
  real, intent(in):: h0(max_zones),dif_ratio(max_zones),tis
  real, intent(in):: depth_max(max_zones),depth_min(max_zones),theta_c_rad(max_zones)
  real, intent(in)::chan_thresh,chan_depth,pf1(grd),contrib_area(imax)
  real, intent(in):: dzdxgs(imax),dzdygs(imax),sec_theta(imax),nl_slope_fac(imax),slope_rad(imax)
  real, intent(inout)::soil_depth(imax)
  real, intent(inout)::unused(imax)
  real, intent(inout):: trans_x(imax),trans_y(imax),d_trans_x_dx(imax),d_trans_y_dy(imax)
  real (kind = 8),intent(in):: nodat,celsiz !
  logical, intent(in):: hump_prod(max_zones), l_mode
  write(*,*) 'Entering subroutine nsd_depth'
  do i=1,imax
    if(abs(nl_slope_fac(i))<=tis) cycle
    trans_x(i)=dzdxgs(i)/nl_slope_fac(i) ! x-component of transport factor
    trans_y(i)=dzdygs(i)/nl_slope_fac(i) ! y-component of transport factor
  end do
  call xyslope(trans_x,pf1,cta,imax,ncol,nrow,d_trans_x_dx,unused,celsiz,celsiz,nodat,no_data_int)
  call xyslope(trans_y,pf1,cta,imax,ncol,nrow,unused,d_trans_y_dy,celsiz,celsiz,nodat,no_data_int)
  chan_ctr=0
  if (l_mode) then  ! Original mode consistent with analytical solutions
    do i=1,imax
      trans_nsd=d_trans_x_dx(i)+d_trans_y_dy(i)
      unused(i) = trans_nsd
      if (trans_nsd < 0.) cycle  ! Avoid negative arguments of log()
      if (abs(trans_nsd) <= 0.0001) cycle ! Avoid division by zero and very small numbers
      if (hump_prod(zo(i))) then
        h1=h0(zo(i))*sec_theta(i)*log((dif_ratio(zo(i))*sec_theta(i))/trans_nsd) 
        call h_solve(sec_theta(i),dif_ratio(zo(i)),trans_nsd,h0(zo(i)),h1,soil_depth(i),l_mode)
      else
        soil_depth(i)=h0(zo(i))*sec_theta(i)*log((dif_ratio(zo(i))*sec_theta(i))/trans_nsd)  
      end if
        if(soil_depth(i) < depth_min(zo(i))) soil_depth(i) = depth_min(zo(i))
        if(soil_depth(i) > depth_max(zo(i))) soil_depth(i) = depth_max(zo(i))
    end do
  else ! Modified mode
    do i=1,imax
      if(nl_slope_fac(i) <= tis) cycle ! Avoid division by zero
      trans_nsd=d_trans_x_dx(i)+d_trans_y_dy(i)
      unused(i) = trans_nsd
      if (trans_nsd < 0.) trans_nsd=-trans_nsd ! Estimate soil depth for positive and negatively curved ground
      if (abs(trans_nsd) <= tis) cycle 
      if (hump_prod(zo(i))) then
        h1=h0(zo(i))*sec_theta(i)*log(trans_nsd/(dif_ratio(zo(i))*sec_theta(i)))
        call h_solve(sec_theta(i),dif_ratio(zo(i)),trans_nsd,h0(zo(i)),h1,soil_depth(i),l_mode)
      else
        soil_depth(i)=h0(zo(i))*sec_theta(i)*log((dif_ratio(zo(i))*sec_theta(i))/trans_nsd) 
      end if
        if(soil_depth(i) < depth_min(zo(i))) soil_depth(i) = depth_min(zo(i))
        if(soil_depth(i) > depth_max(zo(i))) soil_depth(i) = depth_max(zo(i))
! Compare slope angle in channels with 0.2*critical slope angle, and reduce thickness accordingly.
        if(contrib_area(i) > chan_thresh) then 
           if(slope_rad(i)>0.2*theta_c_rad(zo(i))) then
              if(soil_depth(i) > chan_depth) soil_depth(i) = chan_depth ! Set to average alluvium depth.
              chan_ctr = chan_ctr + 1
           end if
        end if
    end do
  endif
  soil_depth_min=minval(soil_depth)
  soil_depth_max=maxval(soil_depth)
  write(*,*) 'Computed depth using NSD model'
  write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Computed depth using NSD model'
  write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Range trans_x', minval(trans_x), maxval(trans_x)
  write(ulog,*) 'Range trans_y', minval(trans_y), maxval(trans_y)
  write(ulog,*) 'Range d_trans_x_dx', minval(d_trans_x_dx), maxval(d_trans_x_dx)
  write(ulog,*) 'Range d_trans_y_dy', minval(d_trans_y_dy), maxval(d_trans_y_dy)
  write(ulog,*) 'Channel grid cells where depth changed, grid-cell threshold: ', chan_ctr, chan_thresh
  return 
end subroutine nsd_depth

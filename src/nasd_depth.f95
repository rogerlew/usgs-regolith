! procedure to compute soil depth based on NASD transport model
! 3 Feb 2015, RLB, Latest revision 13 Aug 2020.
subroutine nasd_depth(ulog,imax,ncol,nrow,grd,celsiz,nodat,no_data_int,cta,chan_thresh,&
  & chan_depth,theta_c_rad,pf1,dzdxgs,dzdygs,sec_theta,nl_slope_fac,slope_rad,&
  & contrib_area,soil_depth,hump_prod,h0,dif_ratio,depth_max,depth_min,tis,&
  & unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,zo,max_zones,l_mode,power)
  implicit none
! LOCAL VARIABLES
  integer::i,ctr,ctr1,neg_ctr,cycl_ctr, chan_ctr
  real::h1,soil_depth_min,soil_depth_max
  real::trans_nasd,aexpn
! FORMAL ARGUMENTS
  integer, intent(in)::ulog,imax,ncol,nrow,grd,no_data_int,cta(ncol,nrow),max_zones,zo(imax)
  real, intent(in):: h0(max_zones),dif_ratio(max_zones),tis
  real, intent(in):: depth_max(max_zones),depth_min(max_zones),theta_c_rad(max_zones)
  real, intent(in):: chan_thresh,chan_depth,pf1(grd),contrib_area(imax),power
  real, intent(in):: dzdxgs(imax),dzdygs(imax),sec_theta(imax),nl_slope_fac(imax),slope_rad(imax)
  real, intent(inout)::soil_depth(imax)
  real, intent(inout):: unused(imax)
  real, intent(inout):: trans_x(imax),trans_y(imax),d_trans_x_dx(imax),d_trans_y_dy(imax)
  real (kind = 8),intent(in):: nodat,celsiz !
  logical, intent(in):: hump_prod(max_zones), l_mode
  write(*,*) 'Entering subroutine nasd_depth'
  ctr=0; ctr1=0
  do i=1,imax
    if(contrib_area(i)>0.) then
      aexpn=contrib_area(i) ! contrib_area^1
      if(power /= 1. .and. aexpn /= 0.) aexpn=aexpn**power
      if(ncol == 1 .or. nrow == 1) aexpn=sqrt(aexpn)
      ctr=ctr+1
    else
      aexpn=0.
    end if
    if(abs(nl_slope_fac(i))<=tis) cycle
    trans_x(i)=aexpn*dzdxgs(i)/nl_slope_fac(i) ! x-component of transport factor
    trans_y(i)=aexpn*dzdygs(i)/nl_slope_fac(i) ! y-component of transport factor
    ctr1=ctr1+1
  end do
  call xyslope(trans_x,pf1,cta,imax,ncol,nrow,d_trans_x_dx,unused,celsiz,celsiz,nodat,no_data_int)
  call xyslope(trans_y,pf1,cta,imax,ncol,nrow,unused,d_trans_y_dy,celsiz,celsiz,nodat,no_data_int)
  ctr=0;neg_ctr=0;cycl_ctr=0
  chan_ctr=0
  if(l_mode) then  ! Original mode consistent with analytical solutions
    do i=1,imax
      trans_nasd=d_trans_x_dx(i)+d_trans_y_dy(i)
      unused(i) = trans_nasd
      if (trans_nasd < 0.) cycle ! Avoid negative arguments of log() 
      if (abs(trans_nasd) <= 0.0001) cycle ! Avoid division by zero and very small numbers
      cycl_ctr=cycl_ctr+1
      if (hump_prod(zo(i))) then
        h1=h0(zo(i))*sec_theta(i)*log((dif_ratio(zo(i))*sec_theta(i))/trans_nasd) ! 
        call h_solve(sec_theta(i),dif_ratio(zo(i)),trans_nasd,h0(zo(i)),h1,soil_depth(i),l_mode)
      else
        soil_depth(i)=h0(zo(i))*sec_theta(i)*log((dif_ratio(zo(i))*sec_theta(i))/trans_nasd) 
      endif
      ctr=ctr+1
      if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
      if(soil_depth(i)<0.) then
        soil_depth(i)=depth_min(zo(i)) 
        neg_ctr=neg_ctr+1
      end if
    end do
  else ! Modified mode
    do i=1,imax
      if(nl_slope_fac(i)<=tis) cycle
      trans_nasd=d_trans_x_dx(i)+d_trans_y_dy(i)
      unused(i) = trans_nasd
      if (trans_nasd<0) trans_nasd=-trans_nasd ! Estimate soil depth for positive and negatively curved ground
      if (trans_nasd<=tis) cycle ! Avoid division by zero
      cycl_ctr=cycl_ctr+1
      if (hump_prod(zo(i))) then
        h1=h0(zo(i))*sec_theta(i)*log(trans_nasd/(dif_ratio(zo(i))*sec_theta(i)))
        call h_solve(sec_theta(i),dif_ratio(zo(i)),trans_nasd,h0(zo(i)),h1,soil_depth(i),l_mode)
      else
        soil_depth(i)=h0(zo(i))*sec_theta(i)*log(trans_nasd/(dif_ratio(zo(i))*sec_theta(i)))
      end if
      ctr=ctr+1
      if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
      if(soil_depth(i)<0.) then
        soil_depth(i)=depth_min(zo(i)) 
        neg_ctr=neg_ctr+1
      end if
! Compare slope angle in channels with 0.2*critical slope angle, and reduce thickness accordingly.
      if(contrib_area(i)>chan_thresh) then 
         if(slope_rad(i)>0.2*theta_c_rad(zo(i))) then
            if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
              chan_ctr = chan_ctr + 1
         end if
      end if
    end do
  end if
  soil_depth_min=minval(soil_depth)
  soil_depth_max=maxval(soil_depth)
  write(*,*) 'Counters:', ctr,neg_ctr,cycl_ctr
  write(*,*) 'Computed depth using NASD model'
  write(*,*) 'Exponent of upslope contributing area = ', power
  write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Computed depth using NASD model'
  write(ulog,*) 'Exponent of upslope contributing area = ', power
  write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Contrib. area counters:',ctr,ctr1
  write(ulog,*) 'Range trans_x', minval(trans_x), maxval(trans_x)
  write(ulog,*) 'Range trans_y', minval(trans_y), maxval(trans_y)
  write(ulog,*) 'Range d_trans_x_dx', minval(d_trans_x_dx), maxval(d_trans_x_dx)
  write(ulog,*) 'Range d_trans_y_dy', minval(d_trans_y_dy), maxval(d_trans_y_dy)
  write(ulog,*) 'Channel grid cells where depth changed, grid-cell threshold: ', chan_ctr, chan_thresh
  return 
end subroutine nasd_depth

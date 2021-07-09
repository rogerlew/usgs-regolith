! Subroutine to compute regolith depth using empirical methods of DeRose
! with modifications to allow for slope curvature and a minimum depth.
! Both polynomial and exponential forms supported.
! 2017, RLB, Latest revision 10 Apr 2019.
subroutine derose(ulog,imax,chan_thresh,chan_depth,theta_c_rad,slope,slope_rad,&
  & dg2rad,contrib_area,plan_view_curv,soil_depth,trans_model,&
  & depth_max,depth_min,C0,C1,C2,zo,max_zones,power)
  implicit none
! LOCAL VARIABLES
  integer:: i, chan_ctr
  real::temp0,soil_depth_min,soil_depth_max
! FORMAL ARGUMENTS
  integer, intent(in):: ulog,imax,max_zones,zo(imax)
  real, intent(in):: C1(max_zones),depth_max(max_zones),depth_min(max_zones)
  real, intent(in):: C0(max_zones),power, C2(max_zones)
  real, intent(in):: chan_thresh,chan_depth,theta_c_rad(max_zones)
  real, intent(in):: slope(imax),slope_rad(imax),contrib_area(imax),plan_view_curv(imax)
  real, intent(inout):: soil_depth(imax)
  real(kind = 8),intent(in)::dg2rad
  character(len=4), intent(in)::trans_model
  write(*,*) 'Entering subroutine DeRose'
  chan_ctr=0
!!  write(*,*) 'i, temp0 '
  select case(trans_model)
  case('PSD') ! or 'DRS2'
    write(*,*) trans_model
    write (*,*) 'Degree of slope polynomial = ', power
    write (ulog,*) 'Degree of slope polynomial = ', power 
    do i=1,imax
      if(slope_rad(i) >= depth_min(zo(i))*dg2rad ) then ! slope angle was computed to radians in main program
        temp0=(C0(zo(i))-C1(zo(i))*slope(i))**power ! (C0-C1*SlopeAngle)^3 from DeRose (1996)
      else
        temp0=0.
      end if
!!      write(*,*) i, temp0
      if(temp0>depth_min(zo(i))) then
        soil_depth(i)=temp0
      else
        soil_depth(i)=depth_min(zo(i))
      end if
      if(slope_rad(i)>theta_c_rad(zo(i))) soil_depth(i)=0.
      if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
! Compare slope angle in channels with 0.1*critical slope angle and reduce thickness accordingly.
      if(contrib_area(i)>chan_thresh) then 
         if(slope_rad(i)>0.1*theta_c_rad(zo(i))) then
            if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
         end if
      end if
    end do
  case('CESD') !same as DRS3
    write(*,*) trans_model
    do i=1,imax
      if(slope_rad(i) >= depth_min(zo(i))*dg2rad ) then ! slope angle was computed to radians in main program
        temp0=(C0(zo(i))-C2(zo(i))*sign(1.,plan_view_curv(i)))*exp(-C1(zo(i))*slope(i)) ! Need to compute or import plan_view_curv array.
      else
        temp0=0.
      end if
      if(temp0>depth_min(zo(i))) then
        soil_depth(i)=temp0
      else
        soil_depth(i)=depth_min(zo(i))
      end if
      if(slope_rad(i)>theta_c_rad(zo(i))) soil_depth(i)=0.
      if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
!!      write(*,*) i,soil_depth(i), depth_max(zo(i)), C0(zo(i)), plan_view_curv(i), &
!!        & exp(-C1(zo(i))*slope(i)), C1(zo(i)), slope(i)
! Compare slope angle in channels with 0.1*critical slope angle and reduce thickness accordingly.
      if(contrib_area(i)>chan_thresh) then 
         if(slope_rad(i)>0.1*theta_c_rad(zo(i))) then
            if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
         end if
      end if
    end do
  case default !('ESD' or 'DRS1')
    write(*,*) trans_model
    do i=1,imax
      if(slope_rad(i) >= depth_min(zo(i))*dg2rad) then ! slope angle was computed to radians in main program
        temp0=C0(zo(i))*exp(-C1(zo(i))*slope(i)) ! C0*exp(-C1*SlopeAngle) from DeRose et al. 1991
      else
        temp0=0.
      end if
!!      write(*,*) i, temp0, zo(i), C1(zo(i)), slope(i), depth_min(zo(i))
      if(temp0 > depth_min(zo(i))) then
        soil_depth(i)=temp0
!!        write(*,*) '+', soil_depth(i)
      else
        soil_depth(i)=depth_min(zo(i))
!!        write(*,*) '-', soil_depth(i)
      end if
      if(slope_rad(i) > theta_c_rad(zo(i))) soil_depth(i)=0.
!!      write (*,*) slope_rad(i), theta_c_rad(zo(i)), soil_depth(i)
! Compare slope angle in channels with 0.1*critical slope angle and reduce thickness accordingly.
      if(contrib_area(i)>chan_thresh) then 
         if(slope_rad(i)>0.1*theta_c_rad(zo(i))) then
            if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
              chan_ctr = chan_ctr + 1
         end if
      end if
    end do
  end select 
  soil_depth_min=minval(soil_depth)
  soil_depth_max=maxval(soil_depth)
  write(*,*) 'Computed depth using modified DeRose formula, ', trans_model
  write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Computed depth using modified DeRose formula, ', trans_model
  write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Channel grid cells where depth changed, grid-cell threshold: ', chan_ctr, chan_thresh
  return 
end subroutine derose

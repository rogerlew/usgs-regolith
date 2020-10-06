! subroutine to estimate soil depth using method of Ho et al. 2012
! with minor modifications to (1) prevent division by zero,  
! (2) specify a minimum soil depth, and (3) allow contributing area to be raised to a power.
! Feb 2015, RLB, Latest revision 8 Sep 2020.

! w0(:) -- inverse of wetness index at minimum depth
subroutine wetness_ndx(ulog,imax,chan_thresh,chan_depth,sc_rad,dg2rad,contrib_area,&
  & slope_rad,soil_depth,depth_min,depth_max,C0,zo,max_zones,power)
  implicit none
! LOCAL VARIABLES
  integer:: i
  real::soil_depth_min,soil_depth_max
  real(kind = 8)::w0,temp0 
! FORMAL ARGUMENTS
  integer, intent(in):: ulog,imax, max_zones,zo(imax)
	real, intent(in):: depth_min(max_zones),C0(max_zones),sc_rad(max_zones),depth_max(max_zones)
  real, intent(in):: chan_thresh,chan_depth,contrib_area(imax)
  real, intent(in):: slope_rad(imax),power
  real, intent(inout)::soil_depth(imax)
  real(kind = 8),intent(in)::dg2rad
  write(*,*) 'Entering subroutine wetness_ndx'
!!  power=2.
  write (*,*) 'Exponent of upslope contributing area',power 
  if (power == 1.0) then ! conventional wetness index
    do i=1,imax
      w0=exp(depth_min(zo(i))/C0(zo(i)))
      if(slope_rad(i)>depth_min(zo(i))*dg2rad ) then ! slope angle was computed to radians in main program
        temp0=contrib_area(i)/tan(slope_rad(i))
      else
        temp0=w0
      end if
      if(temp0>w0) then
        soil_depth(i)=C0(zo(i))*log(temp0) ! h = C0*Log(a/tan(slope)) from Ho et al. (2012)
        if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
! Compare slope angle in channels with 0.2*critical slope angle, and reduce thickness accordingly.
        if(contrib_area(i)>chan_thresh) then 
           if(slope_rad(i)>0.2*sc_rad(zo(i))) then
              if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
           end if
        end if
      else 
        soil_depth(i)=depth_min(zo(i))
      end if
    end do
  else ! power/= 1, modified wetness index
    do i=1,imax
      w0=exp(depth_min(zo(i))/C0(zo(i)))
      if(slope_rad(i)>depth_min(zo(i))*dg2rad ) then ! slope angle was computed to radians in main program
        temp0=contrib_area(i)**power/tan(slope_rad(i))
      else
        temp0=w0
      end if
      if(temp0>w0) then
        soil_depth(i)=C0(zo(i))*log(temp0) ! h = C0*Log(a/tan(slope)) from Ho et al. (2012)
        if(soil_depth(i)>depth_max(zo(i))) soil_depth(i)=depth_max(zo(i))
! Compare slope angle in channels with 0.2*critical slope angle, and reduce thickness accordingly.
        if(contrib_area(i)>chan_thresh) then 
           if(slope_rad(i)>0.2*sc_rad(zo(i))) then
              if(soil_depth(i)>chan_depth) soil_depth(i)=chan_depth ! Set to average alluvium depth.
           end if
        end if
      else 
        soil_depth(i)=depth_min(zo(i))
      end if
    end do
  endif
  soil_depth_min=minval(soil_depth)
  soil_depth_max=maxval(soil_depth)
  write(*,*) 'Computed depth using modified wetness index'
  write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Computed depth using modified wetness index'
  write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Exponent of upslope contributing area',power 
  return 
end subroutine wetness_ndx

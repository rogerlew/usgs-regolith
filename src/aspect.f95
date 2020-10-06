subroutine aspect(dzdx,dzdy,dipdr,dip,slpdr)
! Computes average dip and dip direction for trial landslide mass, Rex Baum, USGS, August 2010
! Reformatted and updated 3/7/2019
  implicit none
! FORMAL ARGUMENTS
  real, intent(in)::dzdx,dzdy
  real, intent(out)::dipdr,dip,slpdr
! LOCAL VARIABLES        
  integer::ix,iy,k
  real ::pi 
  pi=3.1415926535
  if(dzdx>0.) then
    ix=3
  else if(dzdx==0.) then
    ix=2
  else
    ix=1
  endif
  if(dzdy>0.) then
    iy=3
  else if(dzdy==0.) then
    iy=2
  else
    iy=1
  endif
  k=ix+3*(iy-1)
  select case(k)
! For dip direction (dipdr), angles are measured clockwise from north (positive y-axis), East (positive x-axis) is 90 degrees, or pi/2. radians clockwise from north
! For slope direction (slpdr), angles are measured from the positive x-axis (counterclockwise is positive direction), consistent with main branch of the inverse tangent function.
! dzdy<0, plane dips to the north
  case(1); dipdr=pi/2.-atan(dzdy/dzdx)
    slpdr=-(pi-atan(dzdy/dzdx))
  case(2); dipdr=0.
    slpdr=pi/2.
  case(3); dipdr=pi/2.-atan(dzdy/dzdx)+pi
    slpdr=atan(dzdy/dzdx)
! dzdy==0, horizontal north-south line
  case(4); dipdr=pi/2.
    slpdr=pi
  case(5); dipdr=2.*pi ! horizontal plane, dipdr undefined
    slpdr=0.
  case(6); dipdr=3.*pi/2.
    slpdr=0.
! dzdy>0., plane dips to the south
  case(7); dipdr=pi/2.-atan(dzdy/dzdx)
    slpdr=pi+atan(dzdy/dzdx)
  case(8); dipdr=pi
    slpdr=pi/2.
  case(9); dipdr=pi/2.-atan(dzdy/dzdx) +pi
    slpdr=atan(dzdy/dzdx)
  end select
  if(dipdr>2.*pi) dipdr=dipdr-2.*pi
  if(dipdr<0.) dipdr=dipdr+2.*pi
  dip=dzdx*cos(slpdr)+dzdy*sin(slpdr)
!  write(*,*)' k,dip,dzdx,cos(slpdr),dzdy,sin(slpdr) ',k,dip,dzdx,cos(slpdr),dzdy,sin(slpdr)
  dip=atan(dip)
       if(dip<0.) dip=-dip
!  if(dip>pi/2.) write(*,*) 'dip>pi/2.; dzdx,dzdy,slpdr,dipdr,dip', dzdx,dzdy,slpdr,dipdr,dip
end subroutine aspect

! Procedure to compute soil depth based on NDSD transport model (Pelletier & Rasmussen, 2009)
! 25 May 2017, RLB, Latest revision 13 Aug 2020 (removed nds2_depth)
  subroutine ndsd_depth(ulog,imax,ncol,nrow,grd,celsiz,no_data_64,no_data_int,&
     & cell_row,cell_column,indexed_cell_number,elev_index_lkup,cta,pf1,dzdxgs,&
     & dzdygs,del2gs,nl_slope_fac,sec_theta,soil_depth,num_steps,chan_thresh,&
     & chan_depth,contrib_area,theta_c_rad,slope_rad,hump_prod,h0,dif_ratio,&
     & depth_max,depth_min,tis,unused,trans_x,trans_y,d_trans_x_dx,d_trans_y_dy,&
     & zo,max_zones)
  implicit none
! LOCAL VARIABLES
  integer::i,i0,j0,l,iup,jup,m,n, chan_ctr
  integer::ileft,iright,jnorth,jsouth,iup_cn,jup_cn, edge_count
  integer::ctr(4),maxd(imax),big(4)
  real(kind = 8)::n100,imax_dble
  real::Del_dotDelZ_nlso
  real::resid,resid_cell,resid_max,resid_mean,large !  h1 is trial value of h at cell i
  real::h1,hmin,lhs,rhs,soil_depth_min,soil_depth_max,num_steps_flt 
  logical edge(4)
! FORMAL ARGUMENTS 
  integer,intent(in):: ulog, imax, ncol, nrow, grd, no_data_int
  integer,intent(in):: cell_row(imax), cell_column(imax)
  integer,intent(in):: indexed_cell_number(imax), elev_index_lkup(imax)
  integer,intent(in):: cta(ncol,nrow)
  integer,intent(in):: num_steps, max_zones, zo(imax)
  real, intent(in):: h0(max_zones), dif_ratio(max_zones), tis
  real, intent(in):: depth_max(max_zones), depth_min(max_zones)
  real, intent(in):: theta_c_rad(max_zones)
  real, intent(in):: pf1(grd), chan_thresh,chan_depth
  real, intent(in):: contrib_area(imax), slope_rad(imax), sec_theta(imax)
  real, intent(in):: dzdxgs(imax), dzdygs(imax), del2gs(imax), nl_slope_fac(imax)
  real, intent(inout):: soil_depth(imax), unused(imax)
  real, intent(inout):: trans_x(imax), trans_y(imax)
  real, intent(inout):: d_trans_x_dx(imax), d_trans_y_dy(imax)
  real (kind = 8),intent(in):: no_data_64,celsiz !
  logical:: hump_prod(max_zones) 
  write(*,*) 'Entering subroutine ndsd_depth'
! initialize local arrays  
  large=huge(h1)
  do i=1,imax
    if(abs(nl_slope_fac(i))<=tis) cycle ! Avoid division by zero
    trans_x(i)=abs(dzdxgs(i))/nl_slope_fac(i) ! x-component of transport factor
    trans_y(i)=abs(dzdygs(i))/nl_slope_fac(i) ! y-component of transport factor
  end do
  call xyslope(trans_x,pf1,cta,imax,ncol,nrow,d_trans_x_dx,unused,celsiz,celsiz,&
        & no_data_64,no_data_int)
  call xyslope(trans_y,pf1,cta,imax,ncol,nrow,unused,d_trans_y_dy,celsiz,celsiz,&
        & no_data_64,no_data_int)
  write(*,*) 'd_trans_x_dx max, min', maxval(d_trans_x_dx), minval(d_trans_x_dx)
  write(*,*) 'd_trans_y_dy max, min', maxval(d_trans_y_dy), minval(d_trans_y_dy)
  num_steps_flt=float(num_steps)
  if(max_zones == 1) then
    maxd=int(num_steps_flt*depth_max(1))
    write(*,*) 'num_steps, maxd= ', num_steps, maxd(1)
  else
    do i=1,imax
      maxd(i)=int(num_steps_flt*depth_max(zo(i)))
    end do
  end if
  num_steps_flt=float(num_steps)
  resid_max=0.; resid_mean=0.
  ctr=0; big=0
  write(*,*) 'Percent completed using NDSD model'
  do n=1,imax ! In this loop, n is the cell topographic index, where 1 is the highest cell and imax is the lowest cell.
! m is the cell number, starting at the upper left corner of the grid and counting left to right, row by row.
    m=indexed_cell_number(n)
    if(m > imax .or. m < 1) write(*,*) 'n, m = ', n, m
    n100=100.d0*dble(n);imax_dble=dble(imax)
    if(mod(n100,imax_dble)<100) write(*,fmt='(2x,i3,a1)',advance='no')&
                                & int(n100/imax_dble), '%'
    if(nl_slope_fac(m)<=tis) then  ! Use minimum soil depth at and above the angle of stability
       soil_depth(m)=depth_min(zo(m))
       resid_cell = 0.
       cycle
    end if
    j0=cell_row(m) ! row of cell m
    i0=cell_column(m) ! column of cell m
! get elevation index of surrounding cells, paying attention to edge conditions  
! Edge configuration:
!
!            North 3
!  Left 1               Right 2
!            South 4
!
    ileft=0; iright=0; jnorth=0; jsouth=0; edge=.true.; edge_count = 4
    if(i0>1) then
       if(cta(i0-1,j0)>0) then
          ileft=elev_index_lkup(cta(i0-1,j0))
          edge(1)=.false.
          edge_count = edge_count - 1
       end if
    end if
    if(i0<ncol) then
       if(cta(i0+1,j0)>0) then
          iright=elev_index_lkup(cta(i0+1,j0))
          edge(2)=.false.
          edge_count = edge_count - 1
       end if
    end if
    if(j0>1) then
       if(cta(i0,j0-1)>0) then
          jnorth=elev_index_lkup(cta(i0,j0-1))
          edge(3)=.false.
          edge_count = edge_count - 1
       end if
    end if
    if(j0<nrow) then
       if(cta(i0,j0+1)>0) then
          jsouth=elev_index_lkup(cta(i0,j0+1))
          edge(4)=.false.
          edge_count = edge_count - 1
       end if
    end if
    if(edge_count > 2) then ! Use average thickness at isolated and protruding cells.
      soil_depth(m) = (depth_min(zo(m)) + depth_max(zo(m)))/2.
      cycle
    end if
! compare elevation index of neighbors to each other and to cell m (of index n) to find iup and jup.
    iup=0; jup=0; iup_cn=0; jup_cn=0
    if(edge(1)) then
      if (iright<n) then
        iup=iright
        iup_cn=indexed_cell_number(iright)
      end if
    else if(edge(2)) then
       if (ileft<n) then
        iup=ileft
        iup_cn=indexed_cell_number(ileft)
      end if
    else 
      if (iright<ileft) then
        iup=iright
        iup_cn=indexed_cell_number(iright)
      endif
      if (ileft<iright) then
        iup=ileft
        iup_cn=indexed_cell_number(ileft)
      endif
      if (iup>n) iup=0
    end if
    if(edge(3)) then
      if (jsouth<n) then
        jup=jsouth
        jup_cn=indexed_cell_number(jsouth)
      end if
    else if(edge(4)) then
       if (jnorth<n) then
        jup=jnorth
        jup_cn=indexed_cell_number(jnorth)
      end if
    else 
      if (jnorth<jsouth) then
        jup=jnorth
        jup_cn =indexed_cell_number(jnorth)
      endif
      if (jsouth<jnorth) then
        jup=jsouth
        jup_cn =indexed_cell_number(jsouth)
      endif
      if (jup>n) jup=0
    end if
! The following code block assumes vertical thickness on RHS consistent with Pelletier & Rasmussen's equation 21.
! Compute regolith thickness with input from upslope cells, adjusting formula based upon how many upslope cells exist
    if (iup>0 .and. jup>0) then
      resid_cell=large; resid=0.
      ctr(1)=ctr(1)+1
      if(soil_depth(iup_cn)==0. .and. soil_depth(jup_cn)==0.)&
          & write(*,*) 'n,m,j0,i0 ',n,m,j0,i0
      do l=1,maxd(m)
        h1=float(l)/num_steps_flt
        if (hump_prod(zo(m))) then
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*(h1/h0(zo(m)))*&
              & exp(-h1/(h0(zo(m))*sec_theta(m)))
        else
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*exp(-h1/(h0(zo(m))*sec_theta(m)))
        end if
        rhs=soil_depth(iup_cn)*trans_x(iup_cn)/sec_theta(iup_cn) - h1*trans_x(m)/sec_theta(m) +&  !rhs=soil_depth(iup_cn)*trans_x(iup_cn)-h1*trans_x(m)+&
            & soil_depth(jup_cn)*trans_y(jup_cn)/sec_theta(jup_cn) - h1*trans_y(m)/sec_theta(m)  ! & soil_depth(jup_cn)*trans_y(jup_cn)-h1*trans_y(m)
        resid=abs(lhs-rhs)
        if(resid<resid_cell) then
          hmin=h1
          resid_cell=resid
        end if
      end do
          if(resid_cell>10.) then
            big(1)=big(1)+1 
            write(ulog,*) n,m,slope_rad(m),hmin,lhs,rhs,resid_cell
          end if
    else if (iup==0 .and. jup>0) then
      resid_cell=large; resid=0.
      ctr(2)=ctr(2)+1
      do l=1,maxd(m)
        h1=float(l)/num_steps_flt
        if (hump_prod(zo(m))) then
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*(h1/h0(zo(m)))*&
              & exp(-h1/(h0(zo(m))*sec_theta(m)))
        else
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*exp(-h1/(h0(zo(m))*sec_theta(m)))
        end if
! Replace iup terms with with x component of Del(Del(z))/nl_slope_fac)
        rhs=h1*d_trans_x_dx(m)/sec_theta(m) +soil_depth(jup_cn)*trans_y(jup_cn)/sec_theta(jup_cn) -h1*trans_y(m)/sec_theta(m) !rhs=h1*d_trans_x_dx(m)+soil_depth(jup_cn)*trans_y(jup_cn)-h1*trans_y(m) 
        resid=abs(lhs-rhs)
        if(resid<resid_cell) then
          hmin=h1
          resid_cell=resid
        end if
      end do
          if(resid_cell>10.) then
            big(2)=big(2)+1 
            write(ulog,*) n,m,slope_rad(m),hmin,lhs,rhs,resid_cell
          end if
    else if (iup>0 .and. jup==0) then
      resid_cell=large; resid=0.
      ctr(3)=ctr(3)+1
      do l=1,maxd(m)
        h1=float(l)/num_steps_flt
        if (hump_prod(zo(m))) then
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*(h1/h0(zo(m)))*&
              & exp(-h1/(h0(zo(m))*sec_theta(m)))
        else
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*exp(-h1/(h0(zo(m))*sec_theta(m)))
        end if
! Replace jup terms with y component of Del(Del(z))/nl_slope_fac)
        rhs=soil_depth(iup_cn)*trans_x(iup_cn)/sec_theta(iup_cn) - h1*trans_x(m)/sec_theta(m) + h1*d_trans_y_dy(m)/sec_theta(m) !rhs=soil_depth(iup_cn)*trans_x(iup_cn)-h1*trans_x(m)+h1*d_trans_y_dy(m)
        resid=abs(lhs-rhs)
        if(resid<resid_cell) then
          hmin=h1
          resid_cell=resid 
        end if
      end do
          if(resid_cell>10.) then
            big(3)=big(3)+1 
            write(ulog,*) n,m,slope_rad(m),hmin,lhs,rhs,resid_cell
          end if
    else if (iup==0 .and. jup==0) then
      resid_cell=large; resid=0.
      ctr(4)=ctr(4)+1
      Del_dotDelZ_nlso=d_trans_x_dx(m)+d_trans_y_dy(m)
      do l=1,maxd(m)
        h1=float(l)/num_steps_flt
        if (hump_prod(zo(m))) then
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*(h1/h0(zo(m)))*&
              & exp(-h1/(h0(zo(m))*sec_theta(m)))
        else
          lhs=-dif_ratio(zo(m))*celsiz*sec_theta(m)*exp(-h1/(h0(zo(m))*sec_theta(m)))
        end if
        rhs=h1*Del_dotDelZ_nlso ! Replace iup and jup terms with Del(Del(z))/nl_slope_fac)
        resid=abs(lhs-rhs)
        if(resid<resid_cell) then
          hmin=h1
          resid_cell=resid
        end if
      end do
          if(resid_cell>10.) then
            big(4)=big(4)+1 
            write(ulog,*) n,m,slope_rad(m),hmin,lhs,rhs,resid_cell
          end if
    end if
    soil_depth(m)=hmin
    if(hmin > (1.1/num_steps_flt)) then
      if(hmin < (depth_max(zo(m)) - 0.1/num_steps_flt)) then
        if(resid_cell > resid_max) resid_max = resid_cell
        if(resid_cell > 10000) write(ulog,*) n,m,slope_rad(m),hmin,lhs,rhs,resid_cell,'#'
        resid_mean = resid_mean + resid_cell
      end if
    end if
  end do
  resid_mean = resid_mean/float(imax - big(1) - big(2) - big(3) - big(4))
! Adjust depths in channel areas and bare slopes
  chan_ctr=0
  do i=1,imax
    if(soil_depth(i)<depth_min(zo(i))) soil_depth(i)=depth_min(zo(i))
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
  write(*,*) ''
  write(*,*) 'Computed depth using NDSD model'
  write(*,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(ulog,*) 'Computed depth using NDSD model'
  write(ulog,*) 'Range soil_depth: ', soil_depth_min,' - ', soil_depth_max
  write(*,*) 'Maximum & mean residuals ', resid_max, resid_mean
  write(*,*) 'Count of residuals > 10', big
  write(ulog,*) 'Maximum & mean residuals ', resid_max, resid_mean
  write(*,*) 'Counters, c, n-s, e-w, peak: ', ctr
  write(ulog,*) 'Channel grid cells where depth changed, grid-cell threshold: ', chan_ctr, chan_thresh
  return 
  end subroutine ndsd_depth

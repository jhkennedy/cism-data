

program interp_velo_to_bamber_500m
  use netcdf
  use polar_stereo
  implicit none
  integer, parameter :: nx_v = 3010
  integer, parameter :: ny_v = 5460
  integer, parameter :: nx_b = 5001
  integer, parameter :: ny_b = 6001
  integer, parameter :: iters = 50
  real(rp) :: x0_v , y0_v
  real(rp) :: x0_b , y0_b
  real(rp), allocatable :: x_v(:)
  real(rp), allocatable :: y_v(:)
  real(rp), allocatable :: x_b(:)
  real(rp), allocatable :: y_b(:)
  real(rp), allocatable :: x_v_b (:,:)
  real(rp), allocatable :: y_v_b (:,:)
  real(rp), allocatable :: velx_v(:,:)
  real(rp), allocatable :: vely_v(:,:)
  real(rp), allocatable :: errx_v(:,:)
  real(rp), allocatable :: erry_v(:,:)
  real(rp), allocatable :: velx_b(:,:)
  real(rp), allocatable :: vely_b(:,:)
  real(rp), allocatable :: errx_b(:,:)
  real(rp), allocatable :: erry_b(:,:)
  integer :: i , j , n
  real(rp) :: minx_v_b , maxx_v_b
  real(rp) :: miny_v_b , maxy_v_b
  real(rp) :: dx_v_b
  real(rp) :: dy_v_b
  integer :: indx , indy
  integer :: indx_old , indy_old
  integer :: totcount
  real(rp) :: dist(4)
  integer :: ncid , ierr , ydim , xdim , xvar , yvar , vxvar , vyvar , exvar , eyvar


  allocate(x_v(nx_v))
  allocate(y_v(ny_v))
  allocate(x_b(nx_b))
  allocate(y_b(ny_b))
  allocate(x_v_b (nx_v,ny_v))
  allocate(y_v_b (nx_v,ny_v))
  allocate(velx_v(nx_v,ny_v))
  allocate(vely_v(nx_v,ny_v))
  allocate(errx_v(nx_v,ny_v))
  allocate(erry_v(nx_v,ny_v))
  allocate(velx_b(nx_b,ny_b))
  allocate(vely_b(nx_b,ny_b))
  allocate(errx_b(nx_b,ny_b))
  allocate(erry_b(nx_b,ny_b))


  !Form the velocity cartesian grid
  x0_v = -645000
  y0_v = -3370000
  do i = 1 , nx_v
    x_v(i) = x0_v + (i-1)*500._rp
  enddo
  do j = 1 , ny_v
    y_v(j) = y0_v + (j-1)*500._rp
  enddo
  call nc_check( nf90_open( 'coords.nc' , NF90_NOWRITE , ncid ) , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'x' , xvar )            , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'y' , yvar )            , __LINE__ )
  call nc_check( nf90_get_var( ncid , xvar , x_v_b )            , __LINE__ )
  call nc_check( nf90_get_var( ncid , yvar , y_v_b )            , __LINE__ )
  call nc_check( nf90_close( ncid )                             , __LINE__ )

  !Compute minimum and maximum values of projected velocity x and y in order to get a pseudo grid spacing
  minx_v_b = minval(x_v_b)
  miny_v_b = minval(y_v_b)
  maxx_v_b = maxval(x_v_b)
  maxy_v_b = maxval(y_v_b)
  dx_v_b = ( maxx_v_b - minx_v_b ) / nx_v
  dy_v_b = ( maxy_v_b - miny_v_b ) / ny_v
  write(*,*) minval(x_v) , maxval(x_v) , 500.
  write(*,*) minx_v_b , maxx_v_b , dx_v_b
  write(*,*)
  write(*,*) minval(y_v) , maxval(y_v) , 500.
  write(*,*) miny_v_b , maxy_v_b , dy_v_b

  !Form the Bamber cartesian grid
  x0_b = -1300000
  y0_b = -3500000
  do i = 1 , nx_b
    x_b(i) = x0_b + (i-1)*500._rp
  enddo
  do j = 1 , ny_b
    y_b(j) = y0_b + (j-1)*500._rp
  enddo

  !Read in velocity data onto velocity grid
  call nc_check( nf90_open( 'greenland_vel_mosaic500.nc' , NF90_NOWRITE , ncid ) , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'vx' , vxvar )                           , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'vy' , vyvar )                           , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'ex' , exvar )                           , __LINE__ )
  call nc_check( nf90_inq_varid( ncid , 'ey' , eyvar )                           , __LINE__ )
  call nc_check( nf90_get_var( ncid , vxvar , velx_v )                           , __LINE__ )
  call nc_check( nf90_get_var( ncid , vyvar , vely_v )                           , __LINE__ )
  call nc_check( nf90_get_var( ncid , exvar , errx_v )                           , __LINE__ )
  call nc_check( nf90_get_var( ncid , eyvar , erry_v )                           , __LINE__ )
  call nc_check( nf90_close( ncid )                                              , __LINE__ )

  !Loop through Bamber points, interpolating each point from velocity data projected onto Bamber grid.
  totcount = 0
  !$OMP PARALLEL DO PRIVATE(i,j,indx,indy,dist)
  do j = 1 , ny_b
    !$OMP CRITICAL
    totcount = totcount + 1
    !$OMP END CRITICAL
    write(*,*) real(totcount) / ny_b * 100.
    do i = 1 , nx_b
      !If this Bamber point is outside the domain of the projected velocity domain, then continue the loop.
      !A lot of points will fall into this if statement, so it's more efficient to test for it and not interpolate.
      if ( x_b(i) < minx_v_b .or. x_b(i) > maxx_v_b .or. y_b(j) < miny_v_b .or. y_b(j) > maxy_v_b ) then
        velx_b(i,j) = -2.e9_rp
        vely_b(i,j) = -2.e9_rp
        errx_b(i,j) = -2.e9_rp
        erry_b(i,j) = -2.e9_rp
        cycle
      endif
      !Search for the closest projected velocity point
      indy = ny_v
      indx = nx_v
      do n = 1 , iters
        do while ( y_b(j) > y_v_b(indx,indy) )
          indy = indy + 1
          if (indy > ny_v) then
            indy = indy - 1
            exit
          endif
        enddo
        do while ( y_b(j) < y_v_b(indx,indy) )
          indy = indy - 1
          if (indy < 1) then
            indy = indy + 1
            exit
          endif
        enddo
        do while ( x_b(i) > x_v_b(indx,indy) )
          indx = indx + 1
          if (indx > nx_v ) then
            indx = indx - 1
            exit
          endif
        enddo
        do while ( x_b(i) < x_v_b(indx,indy) )
          indx = indx - 1
          if (indx < 1) then
            indx = indx + 1
            exit
          endif
        enddo
      enddo
      !If these conditions are true, then we have found a valid surrounding four projected velocity points from which to interpolate this point
      if ( y_b(j) >= y_v_b(indx,indy) .and. x_b(i) >= x_v_b(indx,indy) .and. indx /= nx_v .and. indy /= ny_v ) then
        dist(1) = sqrt( (y_b(j) - y_v_b(indx  ,indy  ))**2 + (x_b(i) - x_v_b(indx  ,indy  ))**2 )
        dist(2) = sqrt( (y_b(j) - y_v_b(indx+1,indy  ))**2 + (x_b(i) - x_v_b(indx+1,indy  ))**2 )
        dist(3) = sqrt( (y_b(j) - y_v_b(indx  ,indy+1))**2 + (x_b(i) - x_v_b(indx  ,indy+1))**2 )
        dist(4) = sqrt( (y_b(j) - y_v_b(indx+1,indy+1))**2 + (x_b(i) - x_v_b(indx+1,indy+1))**2 )
        velx_b(i,j) = interp_kernel    ( (/ velx_v(indx,indy) , velx_v(indx+1,indy) , velx_v(indx,indy+1) , velx_v(indx+1,indy+1) /) , dist )
        vely_b(i,j) = interp_kernel    ( (/ vely_v(indx,indy) , vely_v(indx+1,indy) , vely_v(indx,indy+1) , vely_v(indx+1,indy+1) /) , dist )
        errx_b(i,j) = interp_kernel_max( (/ errx_v(indx,indy) , errx_v(indx+1,indy) , errx_v(indx,indy+1) , errx_v(indx+1,indy+1) /)        )
        erry_b(i,j) = interp_kernel_max( (/ erry_v(indx,indy) , erry_v(indx+1,indy) , erry_v(indx,indy+1) , erry_v(indx+1,indy+1) /)        )
      else
        velx_b(i,j) = -2.e9
        vely_b(i,j) = -2.e9
        errx_b(i,j) = -2.e9
        erry_b(i,j) = -2.e9
      endif
    enddo
  enddo

  !Write the netcdf file
  call nc_check( nf90_create( 'greenland_vel_mosaic500_bambergrid.nc' , NF90_CLOBBER , ncid ) , __LINE__ )
  call nc_check( nf90_def_dim( ncid , 'y'  , ny_b , ydim )                                    , __LINE__ )
  call nc_check( nf90_def_dim( ncid , 'x'  , nx_b , xdim )                                    , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'y'  , NF90_FLOAT , (/ ydim /)        , yvar )          , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'x'  , NF90_FLOAT , (/ xdim /)        , xvar )          , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'vx' , NF90_FLOAT , (/ xdim , ydim /) , vxvar )         , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'vy' , NF90_FLOAT , (/ xdim , ydim /) , vyvar )         , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'ex' , NF90_FLOAT , (/ xdim , ydim /) , exvar )         , __LINE__ )
  call nc_check( nf90_def_var( ncid , 'ey' , NF90_FLOAT , (/ xdim , ydim /) , eyvar )         , __LINE__ )
  call nc_check( nf90_put_att( ncid , vxvar , 'missing_value' , -2.e9_rp )                    , __LINE__ )
  call nc_check( nf90_put_att( ncid , vyvar , 'missing_value' , -2.e9_rp )                    , __LINE__ )
  call nc_check( nf90_put_att( ncid , exvar , 'missing_value' , -2.e9_rp )                    , __LINE__ )
  call nc_check( nf90_put_att( ncid , eyvar , 'missing_value' , -2.e9_rp )                    , __LINE__ )
  call nc_check( nf90_enddef( ncid )                                                          , __LINE__ )
  call nc_check( nf90_put_var( ncid , xvar  , x_b    )                                        , __LINE__ )
  call nc_check( nf90_put_var( ncid , yvar  , y_b    )                                        , __LINE__ )
  call nc_check( nf90_put_var( ncid , vxvar , velx_b )                                        , __LINE__ )
  call nc_check( nf90_put_var( ncid , vyvar , vely_b )                                        , __LINE__ )
  call nc_check( nf90_put_var( ncid , exvar , errx_b )                                        , __LINE__ )
  call nc_check( nf90_put_var( ncid , eyvar , erry_b )                                        , __LINE__ )
  call nc_check( nf90_close( ncid )                                                           , __LINE__ )



contains



  subroutine nc_check( ierr , line )
    implicit none
    integer, intent(in   ) :: ierr, line
    if (ierr /= NF90_NOERR) then
      write(*,*) 'NETCDF ERROR AT LINE: ',line
      write(*,*) NF90_STRERROR( ierr )
      stop
    endif
  end subroutine nc_check



  function interp_kernel( vals , dists )   result( ival )
    implicit none
    real(rp), intent(in) :: vals(4)
    real(rp), intent(inout) :: dists(4)
    real(rp) :: ival
    real(rp) :: mindist
    integer :: count
    integer :: i
    count = 0
    mindist = 1.e12_rp
    do i = 1 , 4
      if ( vals(i) == -2.e9_rp ) then
        count = count + 1
      else
        if (dists(i) < mindist) then
          ival = vals(i)
          mindist = dists(i)
        endif
      endif
    enddo
    if (count > 2) then
      ival = -2.e9_rp
    endif
    if ( minval(dists) > 500._rp) then
      write(*,*) '           *** More than 500! ***:  ', minval(dists)
      stop
    endif
  end function interp_kernel



  function interp_kernel_max( vals )   result( ival )
    implicit none
    real(rp), intent(in) :: vals(4)
    real(rp) :: ival
    integer :: count
    integer :: i
    ival = 0._rp
    count = 0
    do i = 1 , 4
      if ( vals(i) == -2.e9_rp ) then
        count = count + 1
      else
        ival = max( ival , vals(i) )
      endif
    enddo
    if (count > 2) then
      ival = -2.e9_rp
    endif
  end function interp_kernel_max



  subroutine compute_velo_lonlat( x_v , y_v , x_b , y_b )
    implicit none
    real(rp), intent(in   ) :: x_v( nx_v )
    real(rp), intent(in   ) :: y_v( ny_v )
    real(rp), intent(  out) :: x_b( nx_v , ny_v )
    real(rp), intent(  out) :: y_b( nx_v , ny_v )
    real(rp) :: lambda0_v, lambda0_b, phi1_v, phi1_b
    real(rp) :: mylon, mylat
    integer :: i , j , totcount
    lambda0_v = d2r(-42._rp )
    phi1_v    = d2r( 71._rp )
    lambda0_b = d2r(-39._rp )
    phi1_b    = d2r( 71._rp )
    totcount = 0
    !$OMP PARALLEL DO PRIVATE(i,j,mylon,mylat)
    do j = 1 , ny_v
      !$OMP CRITICAL
      totcount = totcount + 1
      !$OMP END CRITICAL
      write(*,*) real(totcount) / ny_v * 100.
      do i = 1 , nx_v
        call stereo_to_latlon( lambda0_v , phi1_v , mylon , mylat , x_v(i)   , y_v(j)   )
        call latlon_to_stereo( lambda0_b , phi1_b , mylon , mylat , x_b(i,j) , y_b(i,j) )
      enddo
    enddo
  end subroutine compute_velo_lonlat



end program




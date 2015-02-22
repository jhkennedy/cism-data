

module refine_coarsen
  implicit none
  private

  public :: refine



contains



  !Refines a variable by an integer f times, assuming boundaries are the same.
  subroutine refine( nxc , nyc , f , dc , df )
    implicit none
    integer, intent(in   ) :: nxc
    integer, intent(in   ) :: nyc
    integer, intent(in   ) :: f
    real(8), intent(in   ) :: dc( nxc       , nyc       )
    real(8), intent(  out) :: df((nxc-1)*f+1,(nyc-1)*f+1)
    integer :: i, j, k, l
    real(8) :: wts2(2), wts(4)
    do j = 1 , nyc
      do i = 1 , nxc
        df( (i-1)*f+1 , (j-1)*f+1 ) = dc( i , j )
        if ( i /= nxc ) then
          do k = 2 , f
            wts2 = (/ dble(f)/dble(k-1) , dble(f)/dble(f-(k-1)) /)
            wts2 = wts2 / sum(wts2)
            df( (i-1)*f+k , (j-1)*f+1 ) = wts2(1) * dc(i,j) + wts2(2) * dc(i+1,j)
          enddo
        endif
        if ( j /= nyc ) then
          do l = 2 , f
            wts2 = (/ dble(f)/dble(l-1) , dble(f)/dble(f-(l-1)) /)
            wts2 = wts2 / sum(wts2)
            df( (i-1)*f+1 , (j-1)*f+l ) = wts2(1) * dc(i,j) + wts2(2) * dc(i,j+1)
          enddo
        endif
        if ( j /= nyc .and. i /= nxc ) then
          do l = 2 , f
            do k = 2 , f
              wts = (/ dble(f)/dble(   k-1 )*dble(f)/dble(   l-1 ) , &
                       dble(f)/dble(f-(k-1))*dble(f)/dble(   l-1 ) , &
                       dble(f)/dble(   k-1 )*dble(f)/dble(f-(l-1)) , &
                       dble(f)/dble(f-(k-1))*dble(f)/dble(f-(l-1)) /)
              wts = wts / sum(wts)
              df( (i-1)*f+k , (j-1)*f+l ) =  wts(1) * dc( i   , j   ) + &
                                             wts(2) * dc( i+1 , j   ) + &
                                             wts(3) * dc( i   , j+1 ) + &
                                             wts(4) * dc( i+1 , j+1 )
            enddo
          enddo
        endif
      enddo
    enddo
  end subroutine refine



end module refine_coarsen






program refine_coarsen_prog
  use refine_coarsen, only: refine
  use netcdf
  implicit none
  integer, parameter :: f = 5
  character(len=256), parameter :: fname_c = 'Greenland_5km_v1.1.nc'
  character(len=256), parameter :: fname_f = 'Greenland1km.nc'
  real(8), allocatable :: dc(:,:)
  real(8), allocatable :: df(:,:)
  integer :: nxc, nyc
  integer :: nxf, nyf
  integer :: i, j, ierr
  integer :: nc_c, nc_f
  integer :: xdim_c, ydim_c
  integer :: xdim_f, ydim_f, tdim_f
  integer :: var_c, var_f
  ierr = nf90_open( trim(fname_c) , NF90_NOWRITE , nc_c )
  ierr = nf90_open( trim(fname_f) , NF90_WRITE   , nc_f )
  ierr = nf90_inq_dimid( nc_c , 'x1' , xdim_c )
  ierr = nf90_inq_dimid( nc_c , 'y1' , ydim_c )
  ierr = nf90_inquire_dimension( nc_c , xdim_c , len=nxc )
  ierr = nf90_inquire_dimension( nc_c , ydim_c , len=nyc )
  nxf = (nxc-1)*f+1
  nyf = (nyc-1)*f+1
  allocate( dc(nxc,nyc) )
  allocate( df(nxf,nyf) )
  ierr = nf90_inq_varid( nc_c , 'smb' , var_c )
  ierr = nf90_get_var( nc_c , var_c , dc , start=(/1,1,1/) , count=(/nxc,nyc,1/) )
  call refine( nxc , nyc , f , dc , df )
  ierr = nf90_inq_dimid( nc_f , 'x' , xdim_f )
  ierr = nf90_inq_dimid( nc_f , 'y' , ydim_f )
  ierr = nf90_inq_dimid( nc_f , 't' , tdim_f )
  ierr = nf90_redef( nc_f )
  ierr = nf90_def_var( nc_f , 'smb' , NF90_FLOAT , (/xdim_f,ydim_f,tdim_f/) , var_f )
  ierr = nf90_enddef( nc_f )
  ierr = nf90_put_var( nc_f , var_f , df , start=(/1,1,1/) , count=(/nxf,nyf,1/) )
  ierr = nf90_close( nc_c )
  ierr = nf90_close( nc_f )
end program refine_coarsen_prog




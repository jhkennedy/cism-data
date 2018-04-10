

module polar_stereo
  implicit none
  private
  integer, parameter :: rp = selected_real_kind(12)
  real(rp) :: pi = 3.1415926535897932384626433832795028841971_rp
  real(rp) :: a = 6378137.0_rp    !equatorial radius or semimajor axis of the ellipsoid of reference.
  real(rp) :: e = 0.08181919084262_rp    !eccentricity of the ellipsoid. sqrt(1-b^2/a^2).
  real(rp) :: k0 = 1._rp  !central scale factor

  public :: d2r
  public :: r2d
  public :: latlon_to_stereo
  public :: stereo_to_latlon
  public :: rp

contains


  function d2r( deg )
    implicit none
    real(rp) :: deg
    real(rp) :: d2r
    d2r = deg / 180._rp * pi
  end function d2r
  
  
  function r2d( rad )
    implicit none
    real(rp) :: rad
    real(rp) :: r2d
    r2d = rad * 180._rp / pi
  end function r2d
  
  
  subroutine latlon_to_stereo( lambda0 , phi1 , lambda , phi ,  x , y )
    implicit none
    real(rp), intent(in   ) :: lambda0 !longitude east of Greenwich (for longitude west of Greenwich,use a minus sign).
    real(rp), intent(in   ) :: phi1    !latitude
    real(rp), intent(in   ) :: lambda  !longitude east of Greenwich (for longitude west of Greenwich,use a minus sign).
    real(rp), intent(in   ) :: phi     !latitude
    real(rp), intent(  out) :: x
    real(rp), intent(  out) :: y
    real(rp) :: k
    real(rp) :: biga
    real(rp) :: chi
    real(rp) :: m
    real(rp) :: chi1
    real(rp) :: m1
    m  = cos(phi ) / sqrt( 1._rp - e**2 * sin(phi )**2 )
    m1 = cos(phi1) / sqrt( 1._rp - e**2 * sin(phi1)**2 )
    chi  = 2._rp * atan( tan( pi/4._rp + phi /2._rp ) * ( ( 1._rp - e*sin(phi ) ) / ( 1._rp + e*sin(phi ) ) )**(e/2._rp) ) - pi/2._rp
    chi1 = 2._rp * atan( tan( pi/4._rp + phi1/2._rp ) * ( ( 1._rp - e*sin(phi1) ) / ( 1._rp + e*sin(phi1) ) )**(e/2._rp) ) - pi/2._rp
    biga = 2._rp*a*k0*m1 / ( cos(chi1) * ( 1._rp + sin(chi1)*sin(chi) + cos(chi1)*cos(chi)*cos(lambda - lambda0) ) )
    x = biga * cos(chi) * sin(lambda - lambda0)
    y = biga * ( cos(chi1)*sin(chi) - sin(chi1)*cos(chi)*cos(lambda - lambda0) )
  end subroutine latlon_to_stereo
  
  
  subroutine stereo_to_latlon( lambda0 , phi1 , lambda , phi , x , y )
    implicit none
    real(rp), intent(in   ) :: lambda0 !longitude east of Greenwich (for longitude west of Greenwich,use a minus sign).
    real(rp), intent(in   ) :: phi1    !latitude
    real(rp), intent(  out) :: lambda  !longitude east of Greenwich (for longitude west of Greenwich,use a minus sign).
    real(rp), intent(  out) :: phi     !latitude
    real(rp), intent(in   ) :: x
    real(rp), intent(in   ) :: y
    real(rp) :: k
    real(rp) :: biga
    real(rp) :: chi
    real(rp) :: rho
    real(rp) :: chi1
    real(rp) :: m1
    real(rp) :: ce
    real(rp) :: phi_old
    rho = sqrt(x**2+y**2)
    m1 = cos(phi1) / sqrt( 1._rp - e**2 * sin(phi1)**2 )
    chi1 = 2._rp * atan( tan( pi/4._rp + phi1/2._rp ) * ( ( 1._rp - e*sin(phi1) ) / ( 1._rp + e*sin(phi1) ) )**(e/2._rp) ) - pi/2._rp
    ce = 2._rp * atan2( rho * cos(chi1) , 2._rp * a * k0 * m1 )
    if (rho <= 1.e-20_rp) then
      chi = chi1
    else
      chi = asin( cos(ce) * sin(chi1) + (y*sin(ce)*cos(chi1)/rho) )
    endif
    lambda = lambda0 + atan2( x*sin(ce) , rho*cos(chi1)*cos(ce) - y*sin(chi1)*sin(ce) )
    phi_old = chi
    phi = 1.e6_rp
    do while ( abs(phi - phi_old) / phi_old > 1.e-12_rp )
      phi_old = phi
      phi = 2._rp * atan( tan( pi/4._rp + chi/2._rp ) * ( ( 1._rp + e*sin(phi) ) / ( 1._rp - e*sin(phi) ) )**(e/2._rp) ) - pi / 2._rp
    enddo
end subroutine stereo_to_latlon


end module polar_stereo



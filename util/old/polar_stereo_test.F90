
program polar_stereo_test
  use polar_stereo
  implicit none
  real(rp) :: phi1
  real(rp) :: lambda0
  real(rp) :: testx, testy
  real(rp) :: testlambda0, testphi0
  real(rp) :: testlambda, testphi

  phi1    = d2r(71._rp)
  lambda0 = d2r(45._rp)
  testphi0    = phi1 + d2r(1.75_rp)
  testlambda0 = lambda0 - d2r(1.5_rp)
  call latlon_to_stereo( lambda0 , phi1 , testlambda0 , testphi0 , testx , testy )
  call stereo_to_latlon( lambda0 , phi1 , testlambda  , testphi  , testx , testy )
  write(*,*) abs(testphi0 - testphi) / testphi0 , abs(testlambda0 - testlambda) / testlambda0

end program polar_stereo_test


double precision function runge_kuttaf(f1, g1, dt, i)
  integer                                          ::i
  double precision                                 ::f1, dt
  complex*16                                       ::g1
!  double precision, parameter                      ::gamma = 0.2d12
  double precision, parameter                      ::omega_r = 3.2*1.6d-19/6.63d-34
  double precision                                 ::t
  double precision                                 ::omegat
  complex*16                                       ::p
!  double precision funcf
!  funcf(p, t) = aimag(g1+f1)
  runge_kuttaf = 1

!common omegat(1:10)

!write(*,*) omegat
end function runge_kuttaf



program main
  use constants 
  implicit none
  complex*16                                     ::p0 = 0d0
!  double precision                               ::gamma = 0.2d12
  double precision, parameter                    ::f0 = 0d0, t_end = 1d-12, dt = 1d-13,  &
                                                   T = 0.1d-12 , tgrid = t_end/dt
  double precision                               ::omegat(tgrid+1,1) = sqrt(0.1)*1d12 
  double precision                               ::f(tgrid+1,1)
  complex*16                                     ::p(tgrid+1,1)
  double precision                               ::runge_kuttaf
  common omegat
  p0 = runge_kuttaf(1d0, cmplx(1d0,1d0), 1d0, 10)
  gamma = 1
  





!write(*,*) omegat

end program main

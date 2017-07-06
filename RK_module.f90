module RK_module
!07/03/2017 creation
!Runge-Kutta for: FFT of electrical field & polarization; time evolution of polarization and density
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine RK(Pfreq_in, E_in,  f_in, p_in, n, &
               Pfreq_out, E_out, f_out, p_out)
  !07/03/2017 creation
  !Runga_Kutta subroutine
    double precision                                    ::n                    !time step
    double precision, intent(in)                        ::f_in(Ny)
    complex*16, intent(in)                              ::p_in(Ny)
    complex*16, intent(in)                              ::E_in(N_freq)
    complex*16, intent(in)                              ::Pfreq_in(N_freq)
    double precision, intent(out)                       ::f_out(Ny)
    complex*16, intent(out)                             ::p_out(Ny)
    complex*16, intent(out)                             ::E_out(N_freq)
    complex*16, intent(out)                             ::Pfreq_out(N_freq)
    double precision, dimension(Ny)                     ::kf1, kf2, kf3, kf4
    complex*16, dimension(Ny)                           ::kp1, kp2, kp3, kp4
    complex*16, dimension(N_freq)                       ::kE1, kE2, kE3, kE4
    complex*16, dimension(N_freq)                       ::kPfreq1, kPfreq2, kPfreq3, kPfreq4,& 
                                                        kPfreq1_hlp, kPfreq2_hlp, kPfreq3_hlp, &
                                                        kPfreq4_hlp 

    integer                                             ::i1
    n = n - 1.0d0         !to fix the shift of one time step

    !calculate k1
    kp1 = dt * funcp(n,f_in, p_in)               !polarization
    kf1 = dt * funcf(n,p_in)                     !density 
    kE1 = dt * Etime(n) * Ebind * exp((0.0d0, 1.0d0)*(dt*n+tstart)*freqgrid)      
    !Fourier transform of electrical field
    kPfreq1_hlp = dy*sum(y*p_in)/(2.0d0*pi)   !calculate macroscopic polarization
    kPfreq1 = dt * kPfreq1_hlp * exp((0.0d0, 1.0d0)*(tstart+dt*n)*freqgrid)       
    !Fourier transform of polarization

    !calculate k2
    kp2 = dt * funcp(n+0.50d0, f_in + kf1/2.0d0, p_in + kp1/2.0d0)
    kf2 = dt * funcf(n+0.50d0, p_in + kp1/2.0d0)
    kE2 = dt * Etime(n+0.50d0) * Ebind * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*n + 0.5d0*dt))
    kPfreq2_hlp = dy*sum(y*(p_in + kp1/2.0d0))/(2.0d0*pi)
    kPfreq2 = dt * kPfreq2_hlp * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*n + 0.5d0*dt))

    !calculate k3
    kp3 = dt * funcp(n+0.5d0, f_in + kf2/2.0d0, p_in + kp2/2.0d0)
    kf3 = dt * funcf(n+0.5d0, p_in + kp2/2.0d0)
    kE3 = dt * Etime(n+0.5d0) * Ebind * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*n + 0.5d0*dt))
    kPfreq3_hlp = dy*sum(y*(p_in+kp2/2.0d0))/(2.0d0*pi)
    kPfreq3 = dt * kPfreq3_hlp * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*n + 0.5d0*dt))

    !calculate k4
    kp4 = dt * funcp(n+1.0d0, f_in + kf3, p_in + kp3)
    kf4 = dt * funcf(n+1.0d0, p_in + kp3)
    kE4 = dt * Etime(n+1.0d0) * Ebind * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+1.0d0)))
    kPfreq4_hlp = dy*sum(y*(p_in + kp3))/(2.0d0*pi)
    kPfreq4 = dt * kPfreq4_hlp * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+1.0d0)))


    !final step for Runge-Kutta
    p_out = p_in + kp1/6.0d0 +kp2/3.0d0 + kp3/3.0d0 + kp4/6.0d0
    f_out = f_in + kf1/6.0d0 +kf2/3.0d0 + kf3/3.0d0 + kf4/6.0d0
    E_out = E_in + kE1/6.0d0 +kE2/3.0d0 + kE3/3.0d0 + kE4/6.0d0
    Pfreq_out = Pfreq_in + kPfreq1/6.0d0 +kPfreq2/3.0d0 + kPfreq3/3.0d0 + kPfreq4/6.0d0

    !--------output------------
!    write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
!    do i1 = 1, Ny
!      write(701, format_V) p_out(i1)
!    end do
!    write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
!    write(702, format_V) Etime(n)
    !--------end output--------
    end subroutine RK
end module RK_module

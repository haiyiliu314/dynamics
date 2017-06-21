  module constants
  implicit none
  integer, parameter                           ::Nt = 700000
  double precision, parameter                  ::t_end = 35.0d0, &
                                               dt = t_end/dble(Nt)
  integer, parameter                           ::Ny = 400, N_fine = 50, Nphi = 100
  double precision, parameter                  ::ymax = 40.0d0, dy = ymax/dble(Ny)
  double precision, parameter                  ::pi = 4.0d0*datan(1.0d0), &
                                               hbar = 0.6582119514d0, e = 1.60217662d-19
  double precision, parameter                  ::Ebind = 4.18d0, gamma = 0.39d0
  double precision, parameter                  ::sigmat = 0.15d0, tstart = -5.0d0, & 
                                               E_excit = 1.0d-3, shift = 4.0d0
  double precision                             ::coul_mat(Ny, Ny) = 0.0d0, Et(Nt+1) = 0.0d0, &
                                               time(Nt+1) = 0.0d0, argument(Nt+1) = 0.0d0, &
                                               y(Ny)=0.0d0, dexp_argument(Nt+1) = 0.0d0, &
                                               exp_argument(Nt+1) = 0.0d0, &
                                               y_fine(N_fine) = 0.0d0, dy_fine, f(Ny) = 0.0d0, &
                                               ft(Nt+1) = 0.0d0
  complex*16                                   ::p(Ny) = 0.0d0, pt(Nt+1) = 0.0d0
  integer, parameter                           ::N_freq = 800
  complex*16                                   ::p_freq(N_freq) = 0.0d0, E_freq(N_freq) = 0.0d0
  double precision                             ::freqgrid(N_freq) = 0.0d0, test(Nphi) = 0.0d0
  integer                                      ::i1
  character(80)                                ::list_file
  character(len=100)                           ::format_V

  contains
  double precision function Etime(tvia)
    double precision                             ::tvia
    Etime = E_excit*dexp(-(tvia*dt+tstart)*(tvia*dt+tstart)/(sigmat*sigmat)) !no E_bind
  end function Etime

  subroutine readdata
  integer                                        ::i1
  character(len = 100)                           ::format_V
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  open(unit = 100, file = 'Efreq_Ben.dat', status = 'old', action = 'read')
  do i1 = 1,N_freq
    read(100, format_V) freqgrid(i1)
  end do
  close(100)
  freqgrid = (freqgrid*Ebind+4d0*Ebind)/hbar
  write(list_file, '(A)') 'freqgrid.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) freqgrid(i1)
  end do
  close(700)
  end subroutine readdata

  subroutine constant
    integer                                      ::Ndo
    do Ndo = 1, Nt+1
      time(Ndo) = (dble(Ndo)-1.0d0)*dt
      Et(Ndo) = Etime(dble(Ndo))
      argument(Ndo) = -(dble(Ndo)*dt+tstart)*(dble(Ndo)*dt+tstart)/(sigmat*sigmat)
      dexp_argument(Ndo) = dexp(-(dble(Ndo)*dt+tstart)*(dble(Ndo)*dt+tstart)/(sigmat*sigmat))
      exp_argument(Ndo) = exp(-(dble(Ndo)*dt+tstart)*(dble(Ndo)*dt+tstart)/(sigmat*sigmat))
    end do
    do Ndo = 1,Ny
      y(Ndo) = dy*(dble(Ndo) - 0.50d0)
    end do
    dy_fine = dy/dble(N_fine)             !y step for finer grid

  end subroutine constant
    
  subroutine coul_matrix
    implicit none

    double precision, dimension(Nphi)                 ::a, b, c, w1
    integer                                         ::Ndo, Ndo_in
    double precision, dimension(N_fine)              ::fine_grid
    double precision                                 ::b_mid1, b_mid2(Nphi), b_mid3(Nphi)
    call constant
    do Ndo = 1, Nphi
      a(Ndo) = dble(Ndo)
    end do
    do Ndo = 1, N_fine
      fine_grid(Ndo) = (dble(Ndo)-0.50d0)*dy_fine
    end do
    b_mid1 = ( dble(Nphi) + 1.0d0 ) / (2.0d0 * pi )
    b_mid2 = dsin ( (2.0d0 * pi * a ) / ( dble(Nphi)+ 1.0d0 ))
    b_mid3 =  ( a - b_mid1 * b_mid2)
    b = pi / ( dble(Nphi) + 1.0d0 ) * b_mid3
    w1 = pi / ( dble(Nphi) + 1.0d0 ) * ( 1.0d0 - dcos( 2.0d0 * pi * a / ( dble(Nphi) + 1.0d0))) 
    do Ndo = 1, Ny
      y_fine = y(Ndo) - dy/2.0d0 + fine_grid
      do Ndo_in = 1, Nphi
        coul_mat(:, Ndo)=coul_mat(:,Ndo)+(1.0d0/(dsqrt(y(Ndo)*y(Ndo)+y*y -2.0d0*y*y(Ndo)*dcos(b(Ndo_in)))))*2.0d0*y*w1(Ndo_in)
      end do
      coul_mat(Ndo, Ndo) = 0.0d0
      do Ndo_in = 1, N_fine    
        coul_mat(Ndo,Ndo)=coul_mat(Ndo,Ndo) + sum( 1.0d0/ dsqrt( y(Ndo) *y(Ndo) + y_fine(Ndo_in) *y_fine(Ndo_in) - 2.0d0 * y_fine(Ndo_in) * y(Ndo) * dcos( b ) ) / dble(N_fine) * 2.0d0 * y_fine(Ndo_in)  * w1)
      end do  
    end do
    coul_mat = coul_mat/pi*dy
  end subroutine coul_matrix

  function funcp(nt_via, f_via, p_via)
    double precision                                       ::nt_via 
    complex*16                                             ::p_via(Ny), p_sum(Ny), test(Ny)
    double precision                                       ::f_via(Ny), f_sum(Ny)
    complex*16                                             ::funcp(Ny)

    f_sum = matmul(TRANSPOSE(coul_mat), f_via)
    p_sum = matmul(TRANSPOSE(coul_mat), p_via)
    funcp = -(0.0d0,1.0d0)*(y*y*p_via - ((2.0d0*f_sum)*p_via) + shift*p_via - (0.0d0,1.0d0) * &
    gamma * p_via/Ebind - (1.0d0-2.0d0*f_via)*(p_sum+Etime(nt_via)))/hbar*Ebind
  end function funcp

  function funcf(nt_via, p_via)
    double precision                                    ::nt_via 
    complex*16, dimension(Ny)                           ::p_via    
    double precision, dimension(Ny)                     ::funcf
!    funcf = dble(aimag(conjg(Etime(nt_via)+matmul(TRANSPOSE(coul_mat),p_via))*p_via*2.0d0))/hbar*Ebind
     funcf = 0d0
  end function funcf

 subroutine RK(Pfreq_in, E_in,  f_in, p_in, n, &
               Pfreq_out, E_out, f_out, p_out)
    double precision                                    ::n                    !time step
    double precision, dimension(Ny), intent(in)         ::f_in
    double precision, dimension(Ny), intent(out)        ::f_out
    double precision, dimension(Ny)                     ::kf1, kf2, kf3, kf4
    complex*16, dimension(Ny), intent(in)               ::p_in
    complex*16, dimension(Ny), intent(out)              ::p_out
    complex*16, dimension(Ny)                           ::kp1, kp2, kp3, kp4

    complex*16, intent(in)                              ::E_in(N_freq)
    complex*16, intent(out)                             ::E_out(N_freq)
    complex*16                                          ::kE1(N_freq), kE2(N_freq), kE3(N_freq)&
                                                        , kE4(N_freq)

    complex*16, intent(in)                              ::Pfreq_in(N_freq)
    complex*16, intent(out)                             ::Pfreq_out(N_freq)
    complex*16                                          ::kPfreq1(N_freq), kPfreq2(N_freq), kPfreq3(N_freq)&
                                                        , kPfreq4(N_freq), kPfreq1_hlp(N_freq), &
                                                        kPfreq2_hlp(N_freq), kPfreq3_hlp(N_freq), &
                                                        kPfreq4_hlp(N_freq) 

    integer                                             ::num = 1 , i1

    kp1 = dt * funcp(n,f_in, p_in)
    kf1 = dt * funcf(n,p_in)
    kE1 = dt * Etime(n) * Ebind * zexp((0.0d0, 1.0d0)*(dt*n+tstart)*freqgrid)
    kPfreq1_hlp = dy*sum(y*p_in)/(2.0d0*pi)
    kPfreq1 = dt * kPfreq1_hlp * zexp((0.0d0, 1.0d0)*(tstart+dt*n)*freqgrid)

    kp2 = dt * funcp(n+0.50d0, f_in + kf1/2.0d0, p_in + kp1/2.0d0)
    kf2 = dt * funcf(n+0.5d0, p_in + kp1/2.0d0)
    kE2 = dt * Etime(n+0.5d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+0.5d0)))
    kPfreq2_hlp = dy*sum(y*(p_in + kp1/2.0d0))/(2.0d0*pi)
    kPfreq2 = dt * kPfreq2_hlp * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+0.5d0)))

    kp3 = dt * funcp(n+0.5d0, f_in + kf2/2.0d0, p_in + kp2/2.0d0)
    kf3 = dt * funcf(n+0.5d0, p_in + kp2/2.0d0)
    kE3 = dt * Etime(n+0.5d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+0.5d0)))
    kPfreq3_hlp = dy*sum(y*(p_in+kp2/2.0d0))/(2.0d0*pi)
    kPfreq3 = dt * kPfreq3_hlp * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+0.5d0)))

    kp4 = dt * funcp(n+1.0d0, f_in + kf3, p_in + kp3)
    kf4 = dt * funcf(n+1.0d0, p_in + kp3)
    kE4 = dt * Etime(n+1.0d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+1.0d0)))
    kPfreq4_hlp = dy*sum(y*(p_in + kp3))/(2.0d0*pi)
    kPfreq4 = dt * kPfreq4_hlp * zexp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(n+1.0d0)))

    p_out = p_in + kp1/6.0d0 +kp2/3.0d0 + kp3/3.0d0 + kp4/6.0d0
    f_out = f_in + kf1/6.0d0 +kf2/3.0d0 + kf3/3.0d0 + kf4/6.0d0
    E_out = E_in + kE1/6.0d0 +kE2/3.0d0 + kE3/3.0d0 + kE4/6.0d0
    Pfreq_out = Pfreq_in + kPfreq1/6.0d0 +kPfreq2/3.0d0 + kPfreq3/3.0d0 + kPfreq4/6.0d0
    if (int(n) == (num*2500)) then
      write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      do i1 = 1, 1
        write(700, format_V) int(n)
      end do
      num = num+1
    end if
end subroutine RK

end module constants

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
                                               time(Nt+1) = 0.0d0, y(Ny)=0.0d0, &
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
  freqgrid = (freqgrid*Ebind)/hbar
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
    integer                                                ::num = 1 , i1
    f_sum = matmul(TRANSPOSE(coul_mat), f_via)
    p_sum = matmul(TRANSPOSE(coul_mat), p_via)
    funcp = -(0.0d0,1.0d0)*(y*y*p_via - ((2.0d0*f_sum)*p_via) + shift*p_via - (0.0d0,1.0d0) * &
    gamma * p_via/Ebind - (1.0d0-2.0d0*f_via)*(p_sum+Etime(nt_via)))/hbar*Ebind
    if (nt_via == dble(num)*2500) then
      write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      do i1 = 1, Ny
        write(700, format_V) funcp(i1)
      end do
      num = num+1
    end if
  end function funcp

  function funcf(nt_via, p_via)
    double precision                                    ::nt_via 
    complex*16, dimension(Ny)                           ::p_via    
    double precision, dimension(Ny)                     ::funcf
!    funcf = dble(aimag(conjg(Etime(nt_via)+matmul(TRANSPOSE(coul_mat),p_via))*p_via*2.0d0))/hbar*Ebind
     funcf = 0d0
  end function funcf

  subroutine RK_pf(f, g, f1, g1, n)
    double precision                                    ::n                    !time step
    double precision, dimension(Ny), intent(in)         ::f1
    double precision, dimension(Ny), intent(out)        ::f
    double precision, dimension(Ny)                     ::kf1, kf2, kf3, kf4
    complex*16, dimension(Ny), intent(in)               ::g1
    complex*16, dimension(Ny), intent(out)              ::g
    complex*16, dimension(Ny)                           ::kg1, kg2, kg3, kg4
    kg1 = dt * funcp(n,f1, g1)
    kf1 = dt * funcf(n,g1)
    kg2 = dt * funcp(n+0.5d0, f1 + kf1/2.0d0, g1 + kg1/2.0d0)
    kf2 = dt * funcf(n+0.5d0, g1 + kg1/2.0d0)
    kg3 = dt * funcp(n+0.5d0, f1 + kf2/2.0d0, g1 + kg2/2.0d0)
    kf3 = dt * funcf(n+0.5d0, g1 + kg2/2.0d0)
    kg4 = dt * funcp(n+1.0d0, f1 + kf3, g1 + kg3)
    kf4 = dt * funcf(n+1.0d0, g1 + kg3)
    g = g1 + kg1/6.0d0 +kg2/3.0d0 + kg3/3.0d0 + kg4/6.0d0
    f = f1 + kf1/6.0d0 +kf2/3.0d0 + kf3/3.0d0 + kf4/6.0d0
  end subroutine RK_pf

  subroutine RK_E(E_out, E_in, n)
    double precision                                  ::n                    !time step
    complex*16, intent(in)                            ::E_in(N_freq)
    complex*16, intent(out)                           ::E_out(N_freq)
    complex*16                                        ::k1(N_freq), k2(N_freq), k3(N_freq)&
                                                      , k4(N_freq)
    k1 = dt * Etime(n) * Ebind * zexp((0.0d0, 1.0d0)*(dt*n)*freqgrid)
    k2 = dt * Etime(n+0.5d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+0.5d0)))
    k3 = dt * Etime(n+0.5d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+0.5d0)))
    k4 = dt * Etime(n+1.0d0) * Ebind * zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+1.0d0)))
    E_out = E_in + k1/6.0d0 +k2/3.0d0 + k3/3.0d0 + k4/6.0d0
  end subroutine RK_E

  subroutine RK_P(P_out, P_in, n)
    double precision                                  ::n                    !time step
    complex*16, intent(in)                            ::P_in(N_freq)
    complex*16, intent(out)                           ::P_out(N_freq)
    complex*16                                        ::k1(N_freq), k2(N_freq), k3(N_freq)&
                                                      , k4(N_freq)
    k1 = dt * pt(n) * Ebind* zexp((0.0d0, 1.0d0)*(dt*n)*freqgrid)
    k2 = dt * pt(n) * Ebind* zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+0.5d0)))
    k3 = dt * pt(n) * Ebind* zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+0.5d0)))
    k4 = dt * pt(n+1.0d0) * Ebind* zexp((0.0d0, 1.0d0)*freqgrid*(dt*(n+1.0d0)))
    P_out = P_in + k1/6.0d0 +k2/3.0d0 + k3/3.0d0 + k4/6.0d0
  end subroutine RK_P

end module constants

  module constants
  implicit none
  double precision                             ::omegat = sqrt(0.1d0)*1d12
  integer, parameter                           ::Nt = 5d0
  double precision, parameter                  ::f0 = 0d0, t_end = 5d0, &
                                               dt = t_end/dble(Nt)
  integer, parameter                           ::Ny = 100d0, N_fine = 50d0, Nphi = 100d0
  double precision, parameter                  ::ymax = 4d0, dy = ymax/Ny
  double precision, parameter                  ::pi = 4d0*atan(1.), hbar = 4.135667662d0/2d0/pi, &
                                               e = 1.60217662d-19
  double precision, parameter                  ::Ebind = 4.18d0*e, gamma = 0.38d0*e
  double precision, parameter                  ::sigmat = 0.15d0, tstart = -3d0*sigmat,& 
                                               E_excit = 1d-3, shift = 4d0
  double precision                             ::coul_mat(Ny, Ny), Et, time(Nt+1), y(Ny)=0, &
                                               y_fine(N_fine), dy_fine
  complex*16                                   ::y1(Ny) = 0
  contains
  subroutine constant
    integer                                      ::Ndo
    do Ndo = 1d0, Nt+1d0
      time(Ndo) = (Ndo-1d0)*dt
    end do
    do Ndo = 1d0,Ny
      y(Ndo) = dy*(Ndo - 0.5d0)
    end do
    dy_fine = dy/N_fine

  end subroutine constant
  double precision function Etime(tvia)
    double precision                             ::tvia
    Etime = E_excit*dexp(-(tvia*dt+tstart)**2/(sigmat)**2)
  end function Etime

  subroutine coul_matrix
    implicit none

    double precision, dimension(Ny)                 ::a, b, c, w1
    integer                                         ::Ndo, Ndo_in
    double precision, dimension(N_fine)              ::fine_grid
    
    call constant
    do Ndo = 1d0, Ny
      a(Ndo) = Ndo
    end do
    do Ndo = 1d0, N_fine
      fine_grid(Ndo) = (Ndo-0.5d0)*dy_fine
    end do
    b = pi / ( Nphi + 1d0 ) * ( a - ( Nphi + 1d0 ) / (2d0 * pi ) * sin ( (2d0 * pi * a ) / ( Nphi &
    + 1d0 )))
    w1 = pi / ( Nphi + 1d0 ) * ( 1d0 - cos( 2d0 * pi * a / ( Nphi + 1d0))) 

    do Ndo = 1d0, Ny
      y_fine = y(Ndo) - dy/2d0 + fine_grid
      do Ndo_in = 1d0, Nphi
        coul_mat(:, Ndo)=coul_mat(:,Ndo)+(1d0/(sqrt(y(Ndo)**2d0+y**2d0-2d0*y*y(Ndo)*cos (b(Ndo_in)))))*2d0*y*w1(Ndo_in)
      end do
      coul_mat(Ndo, Ndo) = 0d0
      do Ndo_in = 1d0, N_fine
        coul_mat(Ndo,Ndo)=coul_mat(Ndo,Ndo) + sum( 1d0 / sqrt( y(Ndo) ** 2d0 + y_fine(Ndo_in) ** 2d0 - 2d0 * y_fine(Ndo_in) * y(Ndo) * cos( b ) ) / N_fine * 2d0 * y_fine(Ndo_in)  * w1)
      end do
    end do
    coul_mat = coul_mat/pi*dy
  end subroutine coul_matrix
  function funcp1(nt_via, f_via, p_via)
    double precision                                       ::nt_via 
    complex*16                                             ::p_via(Ny), p_sum(Ny), test(Ny)
    double precision                                       ::f_via(Ny), f_sum(Ny)
    complex*16                                             ::funcp1(Ny)
    f_sum = matmul(coul_mat, f_via)
    p_sum = matmul(coul_mat, p_via)
    funcp1 = -(0d0,1d0)*(y**2*p_via - (((2d0*f_sum)*p_via) + shift*p_via - (0d0,1d0) * gamma * p_via/Ebind - (1d0-2d0*f_via)*(p_sum+Etime(nt_via))))/hbar*Ebind
  end function funcp1
  function funcf1(nt_via, p_via)
    double precision                                    ::nt_via 
    complex*16, dimension(Ny)                           ::p_via    
    double precision, dimension(Ny)                     ::funcf1
    funcf1 = aimag(conjg(Etime(nt_via)+matmul(coul_mat,p_via))*p_via*2d0)/hbar*Ebind
  end function funcf1
  subroutine RK_pf_main(f, g, f1, g1, n)
    double precision                                    ::n                    !time step
    double precision, dimension(Ny), intent(in)         ::f1
    double precision, dimension(Ny), intent(out)        ::f
    double precision, dimension(Ny)                     ::kf1, kf2, kf3, kf4
    complex*16, dimension(Ny), intent(in)               ::g1
    complex*16, dimension(Ny), intent(out)              ::g
    complex*16, dimension(Ny)                           ::kg1, kg2, kg3, kg4
    kg1 = dt * funcp1(n,f1, g1)
    kf1 = dt * funcf1(n,g1)
    kg2 = dt * funcp1(n+0.5d0, f1 + kf1/2d0, g1 + kg1/2d0)
    kf2 = dt * funcf1(n+0.5d0, g1 + kg1/2d0)
    kg3 = dt * funcp1(n+0.5d0, f1 + kf2/2d0, g1 + kg2/2d0)
    kf3 = dt * funcf1(n+0.5d0, g1 + kg2/2d0)
    kg4 = dt * funcp1(n+1d0, f1 + kf3, g1 + kg3)
    kf4 = dt * funcf1(n+1d0, g1 + kg3)
    g = g1 + kg1/6d0 +kg2/3d0 + kg3/3d0 + kg4/6d0
    f = f1 + kf1/6d0 +kf2/3d0 + kf3/3d0 + kf4/6d0
  end subroutine RK_pf_main


end module constants

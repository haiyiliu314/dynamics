module coul_matrix
  implicit none
    use constants
  contains
  subroutine coul_matrix_main
    implicit none

    double precision, dimension(Ny)                 ::a, b, c, w1
    double precision, dimension(Ny)                 ::d
    integer                                         ::Ndo, Ndo_in
    double precision, dimension(N_fine)              ::fine_grid
    
    call constant
    do Ndo = 1, Ny
      a(Ndo) = Ndo
    end do
    do Ndo = 1, N_fine
      fine_grid(Ndo) = (Ndo-0.5d0)*dy_fine
    end do
    b = pi / ( Nphi + 1d0 ) * ( a - ( Nphi + 1d0 ) / (2d0 * pi ) * sin ( (2d0 * pi * a ) / ( Nphi &
    + 1d0 )))
    w1 = pi / ( N + 1d0 ) * ( 1d0 - cos( 2d0 * pi * a / ( N + 1d0))) 

    do Ndo = 1, Ny
      y_fine = y(Ndo) - dy/2d0 + fine_grid
      do Ndo_in = 1, Nphi
        d(:, Ndo)=d(:,Ndo)+(1/(sqrt(y(Ndo)**2+y**2-2*y*y(Ndo_in)*cos (b(Ndo_in)))))&
        *2*y*w1(Ndo_in)
      end do
      d(Ndo, Ndo) = 0
      do Ndo_in = 1, N_fine
        d(Ndo,Ndo)=d(Ndo,Ndo) + sum( 1 / sqrt( y(Ndo) ** 2 + y2(Ndo_in) ** 2 - 2 * y2(Ndo_in) * y&
        (Ndo) * cos( b ) ) / N2 * 2 * y2(Ndo_in)  * w1)
      end do
    end do
    coul_mat = d/pi*dy
  
  end subroutine coul_matric_main  
end module coul_mat



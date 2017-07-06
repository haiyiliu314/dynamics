module coul_mat_module
  use constants
  use N_of_grids
  implicit none
  contains
  subroutine coul_matrix
  !07/03/2017 creation
  !calculate Coulomb matrix(non-symmetric, i.e. already contains the integration term, i.e. 2pi*dy*y*V_{y-y'})
    implicit none
    double precision, dimension(Nphi)                ::a, b, c, w1   
!a:integer series   b&w1: Gauss-Chebyshev grid
    integer                                          ::Ndo, Ndo_in   
!Ndo & Ndo_in:integer number for do loop 
    double precision, dimension(N_fine)              ::fine_grid     
!fine_grid: finer grid for singularity removal
    double precision                                 ::b_mid1, b_mid2(Nphi), b_mid3(Nphi)  
!help terms to help calculating Gauss-Chebyshev grid
    call constant
    do Ndo = 1, Nphi
      a(Ndo) = dble(Ndo)             
      !a:integer series from 1 to Nphi
    end do
    do Ndo = 1, N_fine
      fine_grid(Ndo) = (dble(Ndo)-0.50d0)*dy_fine  
      !fine_grid: finer grid for singularity removal
    end do
    b_mid1 = ( dble(Nphi) + 1.0d0 ) / (2.0d0 * pi )
    b_mid2 = dsin ( (2.0d0 * pi * a ) / ( dble(Nphi)+ 1.0d0 ))
    b_mid3 =  ( a - b_mid1 * b_mid2)
    b = pi / ( dble(Nphi) + 1.0d0 ) * b_mid3
    w1 = pi / ( dble(Nphi) + 1.0d0 ) * ( 1.0d0 - dcos( 2.0d0 * pi * a / ( dble(Nphi) + 1.0d0))) 
    !b&w1: Gauss-Chebyshev grid
    do Ndo = 1, Ny
      y_fine = y(Ndo) - dy/2.0d0 + fine_grid   
      !calculate the finer grid for diagonal terms
!------------calculate the Coulomb matrix without removal of singularity-------------
      do Ndo_in = 1, Nphi
        coul_mat(:, Ndo)=coul_mat(:,Ndo)+(1.0d0/(dsqrt(y(Ndo)*y(Ndo)+y*y -2.0d0*y*y(Ndo)*dcos(b(Ndo_in)))))*2.0d0*y*w1(Ndo_in)
      end do
!-----------end calculation----------------------------------------------------------
      coul_mat(Ndo, Ndo) = 0.0d0     
      !eliminate the diagonal terms and calculate them as followed
!-----------calculate the diagonal terms with removal of singularity---------------------
      do Ndo_in = 1, N_fine    
        coul_mat(Ndo,Ndo)=coul_mat(Ndo,Ndo) + sum( 1.0d0/ dsqrt( y(Ndo) *y(Ndo) + y_fine(Ndo_in) *y_fine(Ndo_in) - 2.0d0 * y_fine(Ndo_in) * y(Ndo) * dcos( b ) ) / dble(N_fine) * 2.0d0 * y_fine(Ndo_in)  * w1)
      end do  
!-----------end calculation-----------------------------------------------------------
    end do
    coul_mat = coul_mat/pi*dy   
    !multiply with constants that missed in the do loop
  end subroutine coul_matrix

end module coul_mat_module

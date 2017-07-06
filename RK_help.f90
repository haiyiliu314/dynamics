  module RK_help
  !07/03/2017 Creation
  !All the subroutines and functions that Runge-Kutta needs, including calculation of E(t), time evolution of polarization and density, import frequency grids used for FFT
  use constants
  use N_of_grids
  implicit none
  contains
  double precision function Etime(tvia)
  !07/03/2017 creation
  !function for calculating electrical field(unit: binding energy)
    double precision                             ::tvia                         
    !tvia: time step(converted into double precision)
    Etime = E_excit*exp(-((tvia)*dt+tstart)*((tvia)*dt+tstart)/(sigmat*sigmat)) !no E_bind
  end function Etime

  subroutine readdata
  !07/03/2017 creation
  !import the grid of frequency grid Ben used, to make sure we are using the same frequency grid
  integer                                        ::i1
  character(len = 100)                           ::format_V
  !---------------import---------------
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  open(unit = 100, file = 'Efreq_Ben.dat', status = 'old', action = 'read')
  do i1 = 1,N_freq
    read(100, format_V) freqgrid(i1)
  end do
  close(100)
  !--------------end import--------------
  !-----------output the frequency grid I'm using, to ensure again that Ben and I are using the same grid---------------
  freqgrid = (freqgrid*Ebind+4d0*Ebind)/hbar
  write(list_file, '(A)') 'freqgrid.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) freqgrid(i1)
  end do
  close(700)
!------------end output-----------------------------
  end subroutine readdata





  function funcp(nt_via, f_via, p_via)
  !07/03/2017 creation
  !function to calculate polarization
    double precision                                       ::nt_via   
    !time step(converted into double precision)
    complex*16                                             ::p_via(Ny)  
    !polarization of last time step
    double precision                                       ::f_via(Ny)   
    !density of last time step
    complex*16                                             ::funcp(Ny)
    complex*16                                             ::p_sum(Ny), test(Ny)
    double precision                                       ::f_sum(Ny)
    f_sum = matmul(transpose(coul_mat), f_via)                          
    !calculate \Sigma V_{k'-k}f_{k'}
    p_sum = matmul(transpose(coul_mat), p_via)                          
    !calculate \Sigma V_{k'-k}P_{k'}
    funcp = -(0.0d0,1.0d0)*(y*y*p_via - ((2.0d0*f_sum)*p_via) + shift*p_via - (0.0d0,1.0d0) * &
    gamma * p_via/Ebind - (1.0d0-2.0d0*f_via)*(p_sum+Etime(nt_via)))/hbar*Ebind
    test = y*y*p_via


  end function funcp

  function funcf(nt_via, p_via)
  !07/03/2017 creation
  !function to calculate density
    double precision                                    ::nt_via 
    !nt_via: time step(converted into double precision)
    complex*16, dimension(Ny)                           ::p_via, p_sum    
    !p_via: polarization of last time step
    double precision, dimension(Ny)                     ::funcf
    p_sum = matmul(transpose(coul_mat), p_via)
    !funcf = dble(aimag(conjg(Etime(nt_via)+matmul(TRANSPOSE(coul_mat), p_via))*p_via*2.0d0))/hbar*Ebind
    !calculate the density
    funcf = 0d0
    !set the density as 0
  end function funcf

  end module RK_help

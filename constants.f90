  module constants  
  !07/03/2017 Creation
  !Initialization and calculation of all the constants and global variables 
  !07/06/2017 Revision remove RWA
  use N_of_grids
  implicit none

  double precision, parameter                  ::t_end = 7.0d0            
  !t_end: whole time
  double precision, parameter                  ::dt = t_end/dble(Nt)           
  !length of one time step
  double precision, parameter                  ::ymax = 10.0d0, dy = ymax/dble(Ny)  
  !ymax: maximum value for y   dy: length of one step for y grid
  double precision, parameter                  ::pi = 4.0d0*datan(1.0d0)        
  !physical parameters
  double precision, parameter                  ::hbar = 0.6582119514d0, e = 1.60217662d-19  
  !physical parameters
  double precision, parameter                  ::Ebind = 4.18d0, gamma = 0.39d0, Eg = 1.49d0  
  !Ebind: binding energy(meV)   gamma: dephasing factor (meV) Eg: ground state energy
  double precision, parameter                  ::sigmat = 0.15d0, tstart = -5.0d0  
  !sigmat: sigma for gaussian pulse(E(t) = exp(-t^2/sigmat^2))(ps)   tstart: starting point for time (ps)
  double precision, parameter                  ::E_excit = 1.0d-3, shift = 3.0d0 
  !E_excit: excitation level(unit: binding energy)  shift: shift caused by rotation frame (unit: binding energy)
  double precision                             ::coul_mat(Ny, Ny) = 0.0d0, Et(Nt+1) = 0.0d0
  !coul_mat: Coulomb matrix(non-symmetric)(unit: binding energy) Et: electrical field for excitation
  double precision                             ::y(Ny)=0.0d0, omega_1s = (Ebind + Eg)/hbar     
  !y: grid for y
  double precision                             ::y_fine(N_fine) = 0.0d0, dy_fine, f(Ny) = 0.0d0
  !y_fine: finer grid for removal of singularity   dy_fine: length of one step of finer grid for removal of ringularity  f:density
  double precision                             ::ft(Nt+1) = 0.0d0                             
  !ft: macroscopic density as a function of time
  complex*16                                   ::p(Ny) = 0.0d0, pt(Nt+1) = 0.0d0                  
  !p: polarization   pt: macroscopic polarization as a function of time
  complex*16                                   ::p_freq(N_freq) = 0.0d0, E_freq(N_freq) = 0.0d0   
  !p_freq: Fourier transform of polarization   E_freq: Fourier transform of electrical field
  double precision                             ::freqgrid(N_freq) = 0.0d0, test(Nphi) = 0.0d0     
  !freqgrid: frequency grid used for Fourier transform   test: irrelevant with code
  character(80)                                ::list_file, list_file1, list_file2                
  !list_file, list_file1, list_file2: file name for output
  character(len=100)                           ::format_V                                         
  !format_V: format for output

  contains


  subroutine constant
  !07/03/2017 creation
  !calculate y grid and step length for finder y grid
    integer                                      ::Ndo   !number for do loop
    do Ndo = 1,Ny
      y(Ndo) = dy*(dble(Ndo) - 0.50d0)           !y grid
    end do
    dy_fine = dy/dble(N_fine)             !y step for finer grid

  end subroutine constant

end module constants

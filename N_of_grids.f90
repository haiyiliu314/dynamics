module N_of_grids
!07/03/2017 creation
!define the dimension of all the tensors
  implicit none
  double precision, parameter                  ::scale1 = 0.01d0
  integer, parameter                           ::Nt = 20000*scale1, Nt_RWA = 2000        
  !Nt: number of time step
  integer, parameter                           ::Ny = 350, N_fine = 50, Nphi = 100, Nm_o = 2
  !Ny: number of step for y; N_fine: number of step for fine grid to remove singularity  
  !Nphi: number of step for integrating the angle when calculating Coulomb matrix
  !Nm_o: order of m; 
  integer, parameter                           ::N_freq = 800                                     
  !N_freq: number of points used for frequency grid in Fourier transform
  integer                                      ::i1, num = 1                                      
  !i1, num: integer for do loop for output
end module N_of_grids




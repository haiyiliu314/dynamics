module params
!07/03/2017 creation
!output main parameters that used in the run
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine output_params
    write(list_file1, '(A)') 'params.dat'    
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'   
    open(unit=600,file=list_file1)  
    write(600, *) 'time length' 
    write(600, format_V) t_end
    write(600, *) 'time length for one time step' 
    write(600, format_V) dt   
    write(600, *) 'number of time steps' 
    write(600, format_V) Nt
    write(600, *) 'sigma for Gaussiang electrical field' 
    write(600, format_V) sigmat
    write(600, *) 'maximum value of y' 
    write(600, format_V) ymax
    write(600, *) 'number of finer y points' 
    write(600, format_V) N_fine
    write(600, *) 'number of points for integrating phi' 
    write(600, format_V) Nphi
    close(600)
  end subroutine output_params
end module params

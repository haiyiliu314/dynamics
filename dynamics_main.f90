program main
  !07/03/2017 creation
  !main program to calculate dynamics of SBE
  !result: susceptibility
  !zero-density
  use constants
  use RK_help
  use N_of_grids
  use coul_mat_module
  use RK_module
  use params
  implicit none
  integer                                                       ::Ndo
  double precision                                              ::suscept(N_freq)  !susceptibility
  call constant  
  call coul_matrix
  call readdata
!-------------------open files for output-----------------
      write(list_file, '(A)') 'pk.dat'            !p_y
      open(unit=700,file=list_file)
  write(list_file1, '(A)') 't1_test.dat'           
  open(unit=701,file=list_file1)  
  write(list_file2, '(A)') 'Et.dat'               !E(t)
  open(unit=702,file=list_file2)
!--------------------end open-----------------------------
  do Ndo = 1, Nt

    ft(Ndo) = dy*sum(y*f)/(2.0d0*pi)   !calculate macroscopic density
    pt(Ndo) = dy*sum(y*p)/(2.0d0*pi)   !calculate macroscopic polarization
    call RK(p_freq, E_freq,  f, p, dble(Ndo), &   !input
               p_freq, E_freq, f, p)              !output
!-------------------output----------------------------
      write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      do i1 = 1, Ny
        write(700, format_V) p(i1)
      end do

  end do
  close(702)
  close(701)
      close(700)
!---------------end output   -------------------------
  ft(Ndo) = dy*sum(y*f)/(2.0d0*pi)   !calculate macroscopic density
  pt(Ndo) = dy*sum(y*p)/(2.0d0*pi)   !calculate macroscopic polarization


  suscept = aimag(p_freq / E_freq)   !calculate susceptibility

!-------------------output----------------------------
  write(list_file, '(A)') 'pt.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) pt(i1)
  end do
  close(700)

  write(list_file, '(A)') 'p_freq.dat'           !p_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) p_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'E_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) E_freq(i1)
  end do
  close(700)


  write(list_file, '(A)') 'ft.dat'           !f(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) ft(i1)
  end do
  close(700)

  write(list_file, '(A)') 'coul_mat.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
    write(700, format_V) coul_mat(i1, :)
  end do
  close(700)
!---------------end output--------------------------------------
  call output_params
end program main

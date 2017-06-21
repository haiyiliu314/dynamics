program main
  use constants
  implicit none
  integer                                                       ::Ndo
  double precision                                              ::suscept(N_freq)
  call constant  
  call coul_matrix
  call readdata
  write(list_file, '(A)') 'pk.dat'           !p(t)
  open(unit=700,file=list_file)
  do Ndo = 1, Nt
    ft(Ndo) = dy*sum(y*f)/(2.0d0*pi)
    pt(Ndo) = dy*sum(y*p)/(2.0d0*pi)
    call RK(p_freq, E_freq,  f, p, dble(Ndo), &
               p_freq, E_freq, f, p)
  end do
  ft(Ndo) = dy*sum(y*f)/(2.0d0*pi)
  pt(Ndo) = dy*sum(y*p)/(2.0d0*pi)
  close(700)

  suscept = aimag(p_freq / E_freq)


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

  write(list_file, '(A)') 'Et.dat'           !E(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) Et(i1)
  end do
  close(700)

  write(list_file, '(A)') 'time.dat'           !E(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) time(i1)+tstart
  end do
  close(700)

  write(list_file, '(A)') 'argument.dat'           !E(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) argument(i1)
  end do
  close(700)

  write(list_file, '(A)') 'dexp_argument.dat'           !E(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) dexp_argument(i1)
  end do
  close(700)

  write(list_file, '(A)') 'exp_argument.dat'           !E(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) exp_argument(i1)
  end do
  close(700)

  write(list_file, '(A)') 'coul_mat.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
    write(700, format_V) coul_mat(i1, :)
  end do
  close(700)
end program main

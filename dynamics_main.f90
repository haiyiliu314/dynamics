program main
  use constants
  implicit none
  integer                                                       ::Ndo
  call constant  
  call coul_matrix

  write(list_file, '(A)') 'pk.dat'           !p(t)
  open(unit=700,file=list_file)
  do Ndo = 1, Nt
    ft(Ndo) = dy*sum(y*f)
    pt(Ndo) = dy*sum(y*p)
    call RK_pf(f, p, f, p, dble(Ndo))  
  end do
  close(700)
    ft(Ndo) = dy*sum(y*f)
    pt(Ndo) = dy*sum(y*p)
  call readdata
  do Ndo = 1, Nt
    call RK_P(p_freq, p_freq, dble(Ndo))
  end do
  
  do Ndo = 1, Nt+1
    call RK_E(E_freq, E_freq, dble(Ndo))
  end do

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


  write(list_file, '(A)') 'coul_mat.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
    write(700, format_V) coul_mat(i1, :)
  end do
  close(700)
end program main

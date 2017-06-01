program main
  use constants
  implicit none
  integer                                                       ::Ndo, i1
  character(80)                                                 ::list_file
  character(len=100)                                            ::format_V
  call constant
  call coul_matrix
  do Ndo = 1, Nt
    ft(Ndo) = dy*sum(y*f)
    pt(Ndo) = dy*sum(y*p)
    call RK_pf(f, p, f, p, dble(Ndo))
  end do
    ft(Ndo) = dy*sum(y*f)
    pt(Ndo) = dy*sum(y*p)
  write(list_file, '(A)') 'pt.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', Nt+1, '(", ",SE24.16e3))'
  do i1 = 1, Nt+1
    write(700, format_V) pt(i1)
  end do
  close(700)

  call coul_matrix
  write(list_file, '(A)') 'coul_mat.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
    write(700, format_V) coul_mat(i1, :)
  end do
  close(700)

end program main

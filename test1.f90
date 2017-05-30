
double precision function runge_kuttaf(f1, g1, i)
use constants
  integer                                          ::i
  double precision                                 ::f1
  complex*16                                       ::g1,p
  double precision                                 ::k1, k2, k3, k4
  double precision                                 ::funcf
  funcf(p) = aimag(omegat*p)*2d0            !function to be solved

!------------------Runge-Kutta method------------------
  k1 = dt * funcf(g1)
  k2 = dt * funcf(g1+k1/2d0)
  k3 = dt * funcf(g1+k2/2d0)
  k4 = dt * funcf(g1+k3)
  runge_kuttaf = f1 + k1/6d0 +k2/3d0 + k3/3d0 + k4/6d0
end function runge_kuttaf

complex*16 function runge_kuttap(g1, i)
use constants
  integer                                          ::i
  complex*16                                       ::g1,p
  complex*16                                       ::funcp, k1, k2, k3, k4
  funcp(p) = -(0d0,1d0)*(-((0d0,1d0)*gamma) * p - omegat)               !function to be solved
!------------------Runge-Kutta method------------------
  k1 = dt * funcp(g1)
  k2 = dt * funcp(g1+k1/2d0)
  k3 = dt * funcp(g1+k2/2d0)
  k4 = dt * funcp(g1+k3)
  runge_kuttap = g1 + k1/6d0 +k2/3d0 + k3/3d0 + k4/6d0
end function runge_kuttap

program main
  use constants 
  implicit none
  complex*16                                     ::p0 = (0d0, 0d0)

  double precision                               ::f(Nt+1,1), f1(1,Nt+1)
  complex*16                                     ::p(Nt+1,1)
  double precision                               ::runge_kuttaf

  double precision                               ::ptest(Nt+1, 1), diff
  integer                                        ::i, i1, i2
  complex*16                                     ::runge_kuttap
  character(80)                                  ::list_file
  character(len=100)                             ::format_V

!------------------initialization------------------
  p(1,1) = p0
  f(1,1) = f0
!------------------Runge-Kutta loop------------------
  do i = 1,Nt+1
  p(i+1,1) = runge_kuttap(p(i,1), i)
  f(i+1,1) = runge_kuttaf(f(i,1), p(i,1), i)
  end do
!------------------Test new codes--------------------
call constant
call coul_matrix
  write(list_file, '(A)') 'coul_mat.dat'           !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
    write(700, format_V) coul_mat(i1, :)
  end do
  close(700)
f1 = transpose(f)
write(*,*) y
!------------------Data processing------------------
  ptest = abs(p)**2;
!  diff = abs(abs(p(Nt+1,1)) - sqrt(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2)) / sqrt(abs(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2))
!  write (*,*) 'p(1001) =' , p(Nt+1,1)
!  write(*,*) 'theory peak p', sqrt(0.1d0*((1d0-exp(-0.2d12*1d-12)) / (0.2d12*1d-12))**2)
!  write(*,*) 'diff=', diff

!------------------export data------------------
  write(list_file, '(A)') 'psquared.dat'           !|p|^2
  open(unit=700,file=list_file)

  DO i = 1, Nt+1

    write(700,*)   ptest(i, 1)
  END DO
  close(700)

  write(list_file, '(A)') 'f.dat'                 !f
  open(unit=701,file=list_file)

  DO i = 1, Nt+1

    write(701,*)   f(i, 1)
  END DO
  close(701)

  write(list_file, '(A)') 'dt.dat'               !dt
  open(unit=702,file=list_file)

    write(702,*)  dt

  close(702)

  write(list_file, '(A)') 'Nt.dat'            !numbet of points
  open(unit=703,file=list_file)

    write(703,*)  Nt

  close(703)


!  call system('matlab -r test04202017')
end program main

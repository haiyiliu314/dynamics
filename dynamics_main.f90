program main
  !07/03/2017 creation
  !main program to calculate dynamics of SBE
  !result: susceptibility
  !zero-density
  !07/06/2017 Revision   
  !remove RWA
  !07/07/2017 Upgrade   
  !get rid of arrays with dimension(1,Nt), use elements to increase the speed
  !08/23/2017 Revision
  !calculate the circular harmonics expansion
  use constants
  use RK_help
  use N_of_grids
  use coul_mat_module
  use RK_module
  use params
  implicit none
  integer                                                       ::Ndo, Ndo_m, i2, Ndo_in, Ndo_in1
  call constant  
  call coul_matrix
  call readdata
  write(list_file, '(A)') 'pk2.dat'           !pk
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  write(list_file, '(A)') 'p.dat'             !p
  open(unit=701,file=list_file)
  write(list_file2, '(A)') 'Et.dat'               !E(t)
  open(unit=702,file=list_file2)
  write(list_file, '(A)') 'pt.dat'           !p(t)
  open(unit=703,file=list_file)
  write(list_file, '(A)') 'ft.dat'           !f(t)
  open(unit=704,file=list_file)
  write(list_file, '(A)') 'At.dat'           !A(t)
  open(unit=705,file=list_file)
  write(list_file, '(A)') 'Jt.dat'           !J(t)
  open(unit=706,file=list_file)
  write(format_V1, '(A12, I8, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
  write(format_V2, '(A12, I8, A18)')   '(SE24.16e3, ', 2*Nm_o+1, '(", ",SE24.16e3))'
  write(format_V3, '(A12, I8, A18)')   '(SE24.16e3, ', 2, '(", ",SE24.16e3))'
  do Ndo = 1, Nt
    do Ndo_in = 1, 2*Nm_o+1
      ft(Ndo_in) = dy*sum(y*f(Ndo_in, :))/(2.0d0*pi)   !calculate macroscopic density
      pt(Ndo_in) = dy*sum(y*p(Ndo_in, :))/(2.0d0*pi)   !calculate macroscopic polarization
    end do
    J_THZ_t = dy*sum(y*y*(f(Nm_o, :)+ f(Nm_o+2, :)))/(4.0d0*pi)   !calculate macroscopic density
    write(703, format_V2) abs(pt)
    write(704, format_V2) abs(ft)
    call RK(E_freq, p_freq, f, p, dble(Ndo), A_freq, J_THZ_freq, &              !input
            E_freq, p_freq, f, p, A_freq, J_THZ_freq)                           !output
    write(702, format_V1) real(Etime(dble(Ndo)))
    do i1 = 1, 2*Nm_o+1
      write(701, format_V) real(p(i1, :))
    end do
    write(705, format_V3) real(Atime(dble(Ndo-1))), aimag(Atime(dble(Ndo)))
    write(706, format_V3) real(J_THZ_t), aimag(J_THZ_t)
  end do

  close(700)
  close(701)
  close(704)
  close(703)
  close(705)
  close(706)

!-------------------eigen vectors --------------------
  A = -1d0*TRANSPOSE(coul_mat(1, :, :))
  do i2 = 1, Ny
    A(i2, i2) = A(i2, i2) + y(i2)*y(i2)
  end do 
  LWORK = -1
  CALL ZGEEV(YES, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )

  A1 = -1d0*TRANSPOSE(coul_mat(2, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A1(i2, i2) = A1(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(YES, YES, Ny, A1, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A1, LDA, W1, VL1, LDVL,&
             VR1, LDVR, WORK, LWORK, RWORK, INFO )

  do Ndo = 1, Ny
    VL(:, Ndo) = VL(:, Ndo)/sqrt(real(sum(y*VL(:, Ndo)*conjg(VL(:, Ndo))))*dy*2d0*pi)
    VL1(:, Ndo) = VL1(:, Ndo)/sqrt(real(sum(y*VL1(:, Ndo)*conjg(VL1(:, Ndo))))*dy*2d0*pi)
  end do
  do Ndo = 1, Ny
    VR(:, Ndo) = VR(:, Ndo)/sqrt(real(sum(y*VR(:, Ndo)*conjg(VR(:, Ndo))))*dy*2d0*pi)
    VR1(:, Ndo) = VR1(:, Ndo)/sqrt(real(sum(y*VR1(:, Ndo)*conjg(VR1(:, Ndo))))*dy*2d0*pi)
  end do
!---------------end eigen vectors -------------------- 

!---------------analytical susceptibility-------------
  !calculate p(0)
  p_proj = matmul(p(Nm_o+1, :)*y, VL)*dy*2d0*pi
  !calculate transition matrix of J
  do Ndo = 1, Ny
    do Ndo_in = 1, Ny
      J_tran(Ndo, Ndo_in) = (2d0*pi)*sum(VL1(:,Ndo)*y*y*VR(:,Ndo_in))*dy
      J_tran1(Ndo, Ndo_in) = (2d0*pi)*sum(VL(:,Ndo)*y*y*VR1(:,Ndo_in))*dy
    end do
  end do 
  !calculate S
  resp = 0d0
  resp1 = 0d0
  do Ndo = 1, Ny   !lambda
    do Ndo_in = 1, Ny   !nu
      do Ndo_in1 = 1, Ny   !beta
        resp(Ndo, Ndo_in, :) = resp(Ndo, Ndo_in, :)+(W(Ndo_in1) - W(Ndo))*Ebind&
                             *J_tran1(Ndo_in, Ndo_in1)&
                             *J_tran(Ndo_in1, Ndo)  /  ((W(Ndo_in1) - W(Ndo))*Ebind - freqgrid*hbar&
                             -2d0*(0d0, 1d0)*gamma)*p_proj(Ndo_in)*conjg(p_proj(Ndo))
        resp1(Ndo, Ndo_in, :) = resp1(Ndo, Ndo_in, :)+(W(Ndo_in1) - W(Ndo))*Ebind&
                              *J_tran1(Ndo_in, Ndo_in1)&
                              *J_tran(Ndo_in1, Ndo)  /  ((W(Ndo_in1) - W(Ndo))*Ebind + freqgrid*hbar&
                              -2d0*(0d0, 1d0)*gamma)*p_proj(Ndo_in)*conjg(p_proj(Ndo))
      end do
    end do
  end do
  !calculate summation in (4.21)
  summ = 0d0
  do Ndo = 1, Ny  !lambda
    do Ndo_in = 1, Ny  !nu
      summ = summ + resp(Ndo, Ndo_in, :) - conjg(resp1(Ndo, Ndo_in, :))
    end do
  end do
  !calculate susceptibility
  ana_ses = summ/(freqgrid*freqgrid*(hbar*freqgrid+(0d0, 1d0)*gamma))
!-----------end analytical susceptibility-------------

!-------------------output--------------------
  write(list_file, '(A)') 'coulomb.dat'       !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nm_o+1
    do i2 = 1, Ny
      write(700, format_V) real(coul_mat(i1, i2, :))
    end do
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

  write(list_file, '(A)') 'A_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) A_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'J_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) J_THZ_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'eigenvectors_R.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VR(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvectors_L.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VL(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvalues.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W)
  close(700)

  write(list_file, '(A)') 'eigenvectors_m1_R.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VR1(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvectors_m1_L.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VL1(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvalues_m1.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W1)
  close(700)

  write(list_file, '(A)') 'ana_ses.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
      write(700, format_V) aimag(ana_ses)
  close(700)

!----------------end output------------------
  call output_params
end program main

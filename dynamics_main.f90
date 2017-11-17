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
  integer                                                       ::Ndo, Ndo_m, i2, Ndo_in, i3
  call constant  
  call coul_matrix
  call readdata

  write(format_V1, '(A12, I8, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
  write(format_V2, '(A12, I8, A18)')   '(SE24.16e3, ', 2*Nm_o+1, '(", ",SE24.16e3))'
  write(format_V3, '(A12, I8, A18)')   '(SE24.16e3, ', 2, '(", ",SE24.16e3))'
write(*,*) A_freq_para
!-------------------eigen vectors --------------------
  A = -1d0*TRANSPOSE(coul_mat(1, :, :))
  do i2 = 1, Ny
    A(i2, i2) = A(i2, i2) + y(i2)*y(i2)
  end do 
  LWORK = -1
  CALL ZGEEV(NO, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(NO, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )

  A1 = -1d0*TRANSPOSE(coul_mat(2, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A1(i2, i2) = A1(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(NO, YES, Ny, A1, LDA, W1, VL1, LDVL,&
             VR1, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(NO, YES, Ny, A1, LDA, W1, VL1, LDVL,&
             VR1, LDVR, WORK, LWORK, RWORK, INFO )

  A2 = -1d0*TRANSPOSE(coul_mat(3, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A2(i2, i2) = A2(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(NO, YES, Ny, A2, LDA, W2, VL2, LDVL,&
             VR2, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(NO, YES, Ny, A2, LDA, W2, VL2, LDVL,&
             VR2, LDVR, WORK, LWORK, RWORK, INFO )
!---------------end eigen vectors -------------------- 

!---------------sorting-------------------------------
  Etemp = real(W)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W = Etemp1
  VR = VR_temp

  Etemp = real(W1)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR1(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W1 = Etemp1
  VR1 = VR_temp

  Etemp = real(W2)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR2(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W2 = Etemp1
  VR2 = VR_temp

  write(list_file, '(A)') 'eigenvectors_R.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VR(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvectors_m1_R.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VR1(:,i2))
    end do
  close(700)

  write(list_file, '(A)') 'eigenvectors_m2_R.dat'       !eigen vectors
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do i2 = 1, Ny
      write(700, format_V) real(VR2(:,i2))
    end do
  close(700)
!-----------end sorting-------------------------------
!---------------Normalization-------------------------
  do Ndo = 1, Ny
    VL(:, Ndo) = VL(:, Ndo)/sqrt(real(sum(y*VL(:, Ndo)*conjg(VL(:, Ndo))))*dy/2d0/pi)
    VL1(:, Ndo) = VL1(:, Ndo)/sqrt(real(sum(y*VL1(:, Ndo)*conjg(VL1(:, Ndo))))*dy/2d0/pi)
  end do
  do Ndo = 1, Ny
    VR(:, Ndo) = VR(:, Ndo)/sqrt(real(sum(y*VR(:, Ndo)*conjg(VR(:, Ndo))))*dy/2d0/pi)
    VR1(:, Ndo) = VR1(:, Ndo)/sqrt(real(sum(y*VR1(:, Ndo)*conjg(VR1(:, Ndo))))*dy/2d0/pi)
  end do

  VL1 = VR1
  VL =VR
  VL2 =VR2
!-----------end normalization-------------------------

  write(list_file, '(A)') 'pk2.dat'           !pk
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  write(list_file, '(A)') 'p_r.dat'             !p
  open(unit=701,file=list_file)
  write(list_file2, '(A)') 'Et.dat'               !E(t)
  open(unit=702,file=list_file2)
  write(list_file, '(A)') 'pt.dat'           !p(t)
  open(unit=703,file=list_file)
  write(list_file, '(A)') 'p_i.dat'           !p(t)
  open(unit=710,file=list_file)
  write(list_file, '(A)') 'ft.dat'           !f(t)
  open(unit=704,file=list_file)
  write(list_file, '(A)') 'At.dat'           !A(t)
  open(unit=705,file=list_file)
  write(list_file, '(A)') 'Jt.dat'           !J(t)
  open(unit=706,file=list_file)
  write(list_file, '(A)') 'f_r.dat'             !f
  open(unit=707,file=list_file)
  write(list_file, '(A)') 'f_i.dat'             !f
  open(unit=708,file=list_file)
  write(list_file, '(A)') 'sphere.dat'             !f
  open(unit=709,file=list_file)
  i2 = 1
  i3 = 1
  do Ndo = 1, Nt
    !calculate p(0)

    do Ndo_in = 1, 2*Nm_o+1
      ft(Ndo_in) = dy*sum(y*f(Ndo_in, :))/(2.0d0*pi)   !calculate macroscopic density
      pt(Ndo_in) = dy*sum(y*p(Ndo_in, :))/(2.0d0*pi)/(0.125d0)**2   !calculate macroscopic polarization
    end do
    J_THZ_t = dy*sum(y*y*(f(Nm_o, :)+ f(Nm_o+2, :)))/(4.0d0*pi)   !calculate macroscopic density
    if(Ndo == i2*Nt/Nt_RWA) then
      write(703, format_V2) abs(pt)
      write(704, format_V2) abs(ft)
    end if
    call RK(E_freq, p_freq, f, p, dble(Ndo), A_freq, J_THZ_freq, &              !input
            E_freq, p_freq, f, p, A_freq, J_THZ_freq)                           !output
    if(Ndo == i2*Nt/Nt_RWA) then
      write(702, format_V1) real(Etime(dble(Ndo)))
      write(705, format_V3) real(Atime(dble(Ndo))), aimag(Atime(dble(Ndo)))
      write(706, format_V3) real(J_THZ_t), aimag(J_THZ_t)
      p_proj(i2, :) = matmul(p(Nm_o+1, :)*y, VL)*dy/2d0/pi
      p_proj1(i2, :) = matmul(p(Nm_o+2, :)*y, VL1)*dy/2d0/pi
      p_proj11(i2, :) = matmul(p(Nm_o, :)*y, VL1)*dy/2d0/pi
      p_proj2(i2, :) = matmul(p(Nm_o+3, :)*y, VL2)*dy/2d0/pi
      p_proj22(i2, :) = matmul(p(Nm_o-1, :)*y, VL2)*dy/2d0/pi
      i2 = i2+1
    end if
    if(Ndo == i3*10*Nt/Nt_RWA) then
      do i1 = 1, 2*Nm_o+1
        write(701, format_V) real(p(i1, :))
        write(710, format_V) aimag(p(i1, :))
        write(707, format_V) real(f(i1, :))
        write(708, format_V) aimag(f(i1, :))
        write(709, format_V) abs((f(i1, :) - 0.5d0))*abs(f(i1, :) - 0.5d0)+abs(p(i1,:))*abs(p(i1,:))
      end do
      i3 = i3+1
    end if
  end do
  write(*,*) i2
  close(700)
  close(701)
  close(704)
  close(703)
  close(705)
  close(706)
  close(707)
  close(708)
  close(709)
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

  write(list_file, '(A)') 'p_proj.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nt_RWA
      write(700, format_V) real(p_proj(i1,:))
      write(700, format_V) aimag(p_proj(i1,:))
  end do
  close(700)

  write(list_file, '(A)') 'p_proj1.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nt_RWA
      write(700, format_V) real(p_proj1(i1,:))
      write(700, format_V) aimag(p_proj1(i1,:))
  end do
  close(700)

  write(list_file, '(A)') 'p_proj11.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nt_RWA
      write(700, format_V) real(p_proj11(i1,:))
      write(700, format_V) aimag(p_proj11(i1,:))
  end do
  close(700)

  write(list_file, '(A)') 'p_proj2.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nt_RWA
      write(700, format_V) real(p_proj2(i1,:))
      write(700, format_V) aimag(p_proj2(i1,:))
  end do
  close(700)

  write(list_file, '(A)') 'p_proj22.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nt_RWA
      write(700, format_V) real(p_proj22(i1,:))
      write(700, format_V) aimag(p_proj22(i1,:))
  end do
  close(700)

  write(list_file, '(A)') 'J_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) J_THZ_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'eigenvalues_m2.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W2)
  close(700)

  write(list_file, '(A)') 'y.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) y
  close(700)

  write(*,*) Ny
  write(list_file, '(A)') 'para_mat.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'  
    write(700, format_V) t_end
    write(700, format_V) dt*1d6   
    write(700, *) Nt
    write(700, format_V) ymax
    write(700, *) Ny
    write(700, format_V) A_excit*Ebind*1.129572800243224d3/50d0
    write(700, *) gamma
    write(700, *) 2d0*pi/A_freq_para*1000d0
    write(700, *) Nt_RWA
    write(700, *) tstart_A
    write(700, *) sigmat_A
  close(700)

  write(list_file, '(A)') 'eigenvalues.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W)
  close(700)
  write(list_file, '(A)') 'eigenvalues_m1.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W1)
  close(700)
!----------------end output------------------
  call output_params
end program main

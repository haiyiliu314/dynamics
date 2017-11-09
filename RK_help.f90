  module RK_help
  !07/03/2017 Creation
  !All the subroutines and functions that Runge-Kutta needs, including calculation of E(t), time evolution of polarization and density, import frequency grids used for FFT
  use constants
  use N_of_grids
  implicit none
  contains


  subroutine readdata
  !07/03/2017 creation
  !import the grid of frequency grid Ben used, to make sure we are using the same frequency grid
  !08/20/2017 revision
  !generate a frequency grid from 1 to 9 binding energy
  integer                                        ::i1
  character(len = 100)                           ::format_V
  !i1: do loop parameter
  !format_V: viarable for output format
  !--------set frequence grid, from 1 E_B to 9 E_B ---
  do i1 = 1,N_freq
    freqgrid(i1) =1d0+0.01d0*i1-0.005d0
  end do
!-----------output the frequency grid---------------
  freqgrid = (freqgrid*Ebind)/hbar
  write(list_file, '(A)') 'freqgrid.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) freqgrid(i1)
  end do
  close(700)
!------------end output-----------------------------
  end subroutine readdata

  complex*16 function Etime(tvia)
  !07/03/2017 creation
  !function for calculating electrical field(unit: binding energy)
    double precision                             ::tvia                         
    !tvia: time step(converted into double precision)
!    Etime = E_excit*exp(-((tvia)*dt+tstart)*((tvia)*dt+tstart)/(sigmat*sigmat)) * (1d0 + exp(-(0d0,1d0)*2d0*omega_1s*(tvia*dt+tstart)))!no E_bind non RWA
    Etime = E_excit*exp(-((tvia)*dt+tstart)*((tvia)*dt+tstart)/(sigmat*sigmat)) * (1d0 + 0d0*exp(-(0d0,1d0)*2d0*omega_1s*(tvia*dt+tstart)))!no E_bind   RWA
  end function Etime

  complex*16 function Atime(tvia)
  !09/03/2017 creation
  !function for calculating magnetic field(unit: binding energy)
    double precision                             ::tvia                         
    !tvia: time step(converted into double precision)
!    Atime = A_excit*exp(-((tvia)*dt+tstart_A)*((tvia)*dt+tstart_A)/(sigmat_A*sigmat_A))*cos(A_freq_para*(tvia*dt+tstart_A)/hbar)*Ebind*2d0
!    Atime = A_excit*exp(-((tvia)*dt+tstart_A)*((tvia)*dt+tstart_A)/(sigmat_A*sigmat_A))*(exp((0d0, 1d0)*A_freq_para*(tvia*dt+tstart_A)/hbar)+exp(-(0d0, 1d0)*A_freq_para*(tvia*dt+tstart_A)/hbar))*Ebind        ! unit:meV    
    Atime = A_excit*exp(-((tvia)*dt+tstart_A)*((tvia)*dt+tstart_A)/(sigmat_A*sigmat_A))*cos((A_freq_para)*(tvia*dt+tstart_A)/hbar)*Ebind        ! unit:meV                    
  end function Atime


  subroutine RHS(nt_via, f_via, p_via, &
                 p_out, f_out, kE_out, kPfreq_out, kA_out, kJ_out)
  !08/25/2017 creation
  !calculate the right hand side of RK
    double precision, intent(in)                           ::nt_via   
    !time step(converted into double precision)
    complex*16, intent(in)                                 ::p_via(2*Nm_o+1, Ny)  
    !polarization of last time step
    complex*16, intent(in)                                 ::f_via(2*Nm_o+1, Ny)   
    !density of last time step
    complex*16, intent(out)                                ::p_out(2*Nm_o+1, Ny)  
    !polarization output
    complex*16, intent(out)                                ::f_out(2*Nm_o+1, Ny) 
    !density output
    complex*16, intent(out)                                ::kE_out(N_freq)
    complex*16, intent(out)                                ::kPfreq_out(N_freq)
    complex*16, intent(out)                                ::kA_out(N_freq)
    complex*16, intent(out)                                ::kJ_out(N_freq)
    !kE_out, kPfreq_out, kA_out, kJ_out: output parameters for Runge-Kutta methods
    integer                                                ::Ndo_m, Ndo_m1, Ndo_phi
    !Ndo_m, Ndo_m1: arguments for do loops, indices of m and m'
    complex*16                                             ::p_sum_part(Ny), &
                                                             f_sum_part(Ny), &
                                                             pp_sum(Ny), &
                                                             pp_sum_plus(Ny), &
                                                             pf_sum(Ny), &
                                                             fp_sum(Ny), &
                                                             p_sum_part_m(Ny)
    !help parameters to compute all the summations
    !p_sum_part: \Sigma V_{k'-k}P_{k'}
    !f_sum_part: \Sigma V_{k'-k}f_{k'}
    !pp_sum: \Sigma_{m'} P^*_{k, -m+m'} \Sigma_{k'}V^{m'}_{k,k'} P_{k', m'}
    !pp_sum_plus: \Sigma_{m'} P_{k, m+m'} \Sigma_{k'}V^{-m'}_{k,k'} P^*_{k', m'}
    !pf_sum: \Sigma_{m'} P^*_{k, m-m'} \Sigma_{k'}V^{m'}_{k,k'} f_{k', m'}
    !fp_sum: \Sigma_{m'} f^*_{k, m-m'} \Sigma_{k'}V^{m'}_{k,k'} P_{k', m'}
    !p_sum_part_m: \Sigma V^m_{k, k'}P^m_{k'}
    complex*16                                             ::pt1, J_THZ_t1
    !pt1:\Sigma_k P_{k}
    !J_THZ_t1: \Sigma_k J_{k} f_k

    pt1 = 0.0d0
    !Initialization of P(t)
    do Ndo_m = 1, 2*Nm_o+1
      p_sum_part = 0.0d0
      pp_sum = 0.0d0
      pp_sum_plus = 0.0d0
      pf_sum = 0.0d0
      fp_sum = 0.0d0
      p_sum_part_m = 0.0d0
      f_sum_part = 0.0d0
      do Ndo_m1 = 1, 2*Nm_o+1

        f_sum_part = matmul(f_via(Ndo_m1, :), coul_mat(abs(Ndo_m1-Nm_o-1)+1, :, :))                                                                               
        p_sum_part = matmul(p_via(Ndo_m1, :), coul_mat(abs(Ndo_m1-Nm_o-1)+1, :, :))                       
        if(abs(Ndo_m-Ndo_m1)<Nm_o+1) then
          fp_sum = fp_sum+f_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*p_sum_part
          pf_sum = pf_sum+p_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*f_sum_part
          pp_sum = pp_sum+conjg(p_via(-Ndo_m+Ndo_m1+(Nm_o+1), :))*p_sum_part
        end if
        if(abs(Ndo_m+Ndo_m1 - 2*(Nm_o+1))<Nm_o+1) then
          pp_sum_plus = pp_sum_plus+p_via(Ndo_m+Ndo_m1-(Nm_o+1), :)*conjg(p_sum_part)
        end if
!        write(700, format_V) aimag(conjg(p_via(-Ndo_m+Ndo_m1+(Nm_o+1), :))*p_sum_part)
      end do

      !magnetism coupling, P^{m-1}+P^{m+1}, and truncate at Nm_o, -Nm_o
      if(abs(Ndo_m-Nm_o-1)<Nm_o) then
        coup = (p_via(Ndo_m-1, :)+p_via(Ndo_m+1, :))/2d0
        else if(Ndo_m == 1) then
        coup = (p_via(Ndo_m+1, :))/2d0
        else if(Ndo_m == 2*Nm_o+1) then
        coup = (p_via(Ndo_m-1, :))/2d0
      end if

      p_sum_part_m = matmul(p_via(Ndo_m, :), coul_mat(abs(Ndo_m-Nm_o-1)+1, :, :))
      p_out(Ndo_m, :) = -(0.0d0,1.0d0)*(y*y*p_via(Ndo_m, :) - 2.0d0*pf_sum - y*Atime(nt_via)*coup+&
                        shift*p_via(Ndo_m, :) - (0.0d0,1.0d0) * gamma * p_via(Ndo_m, :)/Ebind -&
                        ((abs(dble(Ndo_m == (Nm_o+1)))-2.0d0*f_via(Ndo_m, :)) *Etime(nt_via)+&
                        (p_sum_part_m - 2.0d0*fp_sum) ) )/hbar*Ebind*dt
      f_out(Ndo_m, :) = (Etime(nt_via)*(p_via(Ndo_m, :) - conjg(p_via((2*Nm_o+2)-Ndo_m, :))) &
                        +(pp_sum_plus-pp_sum)- (0.0d0,1.0d0) *(dble(Ndo_m==Nm_o+1)+1)*gamma * f_via(Ndo_m, :)/Ebind)&
                        /hbar*Ebind*dt/(0.0d0, 1.0d0)
    end do
    J_THZ_t1 = 0d0
    J_THZ_t1 = dy*sum(y*y*(f_via(Nm_o, :)+f_via(Nm_o+2, :)))/(2.0d0*pi)
    kE_out = dt * Etime(nt_via) * Ebind * exp((0.0d0, 1.0d0)*(dt*nt_via+tstart)*freqgrid)  
    kA_out = dt * Atime(nt_via) * exp((0.0d0, 1.0d0)*(dt*nt_via+tstart)*freqgrid) 
    pt1 = dy*sum(y*p_via(Nm_o+1, :))/(2.0d0*pi)   !calculate macroscopic polarization
    kPfreq_out = dt * pt1 * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(nt_via)))
    kJ_out = dt * J_THZ_t1 * exp((0.0d0, 1.0d0)*freqgrid*(tstart+dt*(nt_via)))
  end subroutine RHS
  end module RK_help

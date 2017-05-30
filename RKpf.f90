module RK_pf
  implicit none
  contains
  subroutine RK_pf_main(f, g, f1, g1, n, dimen)
  implicit none
  use constants
  double precision                                    ::n, nt_via                     !time step
  integer                                             ::dimen
  double precision, dimension(dimen), intent(in)      ::f1
  double precision, dimension(dimen), intent(out)     ::f
  double precision, dimension(dimen)                  ::f_via, kfi, kf2, kf3, kf4
  complex*16, dimension(dimen), intent(in)            ::p1
  complex*16, dimension(dimen), intent(out)           ::p
  complex*16, dimension(dimen)                        ::p_via, kg1, kg2, kg3, kg4
  double precision                                    ::funcf
  complex*16                                          ::funcp
  funcf(nt_via, p_via) = aimag(conjg(Et(nt_via)+matmul(coul_mat,& p_via))*p*2d0)/hbar*Ebind
  funcp(nt_via, f_via, p_via) = -(0d0,1d0)*(y**2*p - (((2d0*matmul(A,fvia))*pvia) + shift*p &
  - (0d0,1d0) * gamma * p/Ebind - (1d0-2d0*f)*(matmul(A, p)+Etime(nt_via))))/hbar*Ebind
  kg1 = dt * funcp(n,f1, g1)
  kf1 = dt * funcf(n,g1)
  kg2 = dt * funcp(n+0.5d0, f1 + kf1/2d0, g1 + kg1/2d0)
  kf2 = dt * funcf(n+0.5d0, g1 + kg1/2d0)
  kg3 = dt * funcp(n+0.5d0, f1 + kf2/2d0, g1 + kg2/2d0)
  kf3 = dt * funcf(n+0.5d0, g1 + kg2/2d0)
  kg4 = dt * funcp(n+1d0, f1 + kf3, g1 + kg3)
  kf4 = dt * funcf(n+1d0, g1 + kg3)
  g = g1 + kg1/6d0 +kg2/3d0 + kg3/3d0 + kg4/6d0
  f = f1 + kf1/6d0 +kf2/3d0 + kf3/3d0 + kf4/6d0
  
  end subroutine RK_pf_main
  end module RK_pf

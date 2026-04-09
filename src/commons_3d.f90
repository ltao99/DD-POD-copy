MODULE commons_3d

  USE sizes
  IMPLICIT NONE

  !Subdomain

  INTEGER :: sub, check

  REAL(DP) :: start_aux, end_aux, factorization_time, solve_time

  INTEGER :: auxz! nsz - 1 if nsub_proc > 1, 0 otherwise

  INTEGER :: FLAG_globy ! 1 if nsub_proc > 1, 0 otherwise

  INTEGER :: FLAG_globz ! 1 if nsub_proc > 1, 0 otherwise

  INTEGER :: FLAG_DD_Z 

  LOGICAL :: new_case =  .FALSE.
  ! Mode

  CHARACTER(2) ::  modo
  INTEGER :: mnodos

  LOGICAL :: is_factorized = .FALSE.

  ! Input

  REAL(DP) :: xmin, xmax, ymin, ymax, zmin, zmax, freq, &
       omega, e0,e1,e1o,e2,e3,e4,e5,sigma0,sigma1,sigma2, &
       sigma3,sigma4,sigma5,a1_fac

  Complex(DP) :: coef_i

  REAL(DP) :: hx(ngx), hy(ngy), hz(ngz), zsup

  REAL(DP) :: hsz, hsx, hsy

  REAL(DP) :: h3, hxz, hyz, hxy, hxy_z, hxz_y, hyz_x

  INTEGER :: nx, ny, nz, nmat, izsup, numfreq

  INTEGER :: matnum(0:ngx+1,0:ngy+1,0:ngz+1), matnumloc(0:ngx+1,0:(ngy/nsy)+1,0:(ngz/nsz)+1)

  ! Total number of nodes
  integer :: total_g_nod, total_l_nod

  ! Connectivity Matrix
  integer :: lc_nod(ngx*(ngy/nsy)*(ngz/nsz),12)
  integer :: gl_nod(ngx*ngy*ngz,12)

  ! Non zero elements of diagonal and upper triangular side of LSE
  integer :: nnz,nrow, nnz_rb
  integer, allocatable:: ip(:), ij(:)
  integer, allocatable:: ip_csr(:), ij_csr(:)
  complex(dp) , allocatable:: apj(:)
  complex(dp) , allocatable:: apj_csr(:)
  complex(dp) , allocatable:: apj_ser(:,:)

  ! RHS
  complex(DP), allocatable :: bglobal(:)
  complex(DP), allocatable :: bglobal_constant(:)
  complex(DP), allocatable :: bglobal_dd(:)
  complex(DP), allocatable :: bglobal_temp(:)
  complex(DP), allocatable:: F_TE(:),F_TM(:)

  !ROM

  INTEGER:: basis_map(nsub_proc,2)
  INTEGER:: rank
  complex(DP), allocatable:: bglobal_rb(:)
  complex(DP), allocatable:: bglobal_rb_constant(:)
  complex(DP), allocatable:: bglobal_rb_dd(:)
  complex(DP), allocatable:: sol_MUMPS_rb(:)
  complex(DP), allocatable:: init_sol(:)


  TYPE :: FO_solution
  COMPLEX(DP), ALLOCATABLE :: full_order(:)     ! The reduced POD basis matrix
  END TYPE FO_Solution
  TYPE(FO_Solution) :: test_case(nsub_proc, 2)  ! Stores one basis per subdomain × solve

  TYPE :: BasisSet
  COMPLEX(DP), ALLOCATABLE :: U(:,:)     ! The reduced POD basis matrix
  COMPLEX(DP), ALLOCATABLE :: UH(:,:)
  COMPLEX(DP), ALLOCATABLE :: U_interf(:,:)
  COMPLEX(DP), ALLOCATABLE :: UH_interf(:,:)
  COMPLEX(DP), ALLOCATABLE :: X_mean(:)  ! Snapshots mean
  COMPLEX(DP), ALLOCATABLE :: KX_mean(:) ! Center Correction
  END TYPE BasisSet
  TYPE(BasisSet), ALLOCATABLE :: basis(:,:)  ! Stores one basis per subdomain × solve

  TYPE :: ReducedK
  COMPLEX(DP), POINTER :: Krb(:,:)  ! The reduced POD basis matrix
  END TYPE ReducedK
  TYPE(ReducedK) :: reduced_stiffness(nsub_proc, 2)  ! Stores 2 Krb per subdomain

  TYPE :: ReducedKcoo
  integer, allocatable:: ip_rb(:)
  integer, allocatable:: ij_rb(:)
  complex(DP), allocatable :: apj_rb(:)  ! The reduced POD basis matrix
  END TYPE ReducedKcoo
  TYPE(ReducedKcoo) :: reduced_stiffnes_coo(nsub_proc, 2)  ! Stores 2 Krb per subdomain

  ! Solution of LSE
  complex(DP), allocatable:: sol_MUMPS(:)

  ! FEM Parameters
  
  real(DP) :: deltaf(0:mnx+1), &
       deltab(0:mnx+1), deltan(0:mnz+1), deltas(0:mnz+1), &
       deltae(0:mny+1), deltaw(0:mny+1)

  real(DP) :: deltabf(0:mnx+1), &
       deltabb(0:mnx+1), deltabn(0:mnz+1), deltabs(0:mnz+1), &
       deltabe(0:mny+1), deltabw(0:mny+1)
 
  Integer::  delta_y, delta_z
  INTEGER, allocatable:: interf_y(:),interf_z(:)
  INTEGER, allocatable:: dofs_interface(:)
  INTEGER:: ndofs_interface

  REAL(DP) :: sig(0:mnmat+1), sig1(0:mnmat+1), sig2(0:mnmat+1), sig3(0:mnmat+1), sigb(0:mnmat+1), &
       sigp(0:mnmat+1), xnod(0:ngx+1), ynod(0:ngy+1), znod(0:ngz+1), &
       alphae(0:mnx+1,0:mnz+1), alphaw(0:mnx+1,0:mnz+1), &
       alphan(0:mnx+1,0:mny+1), alphas(0:mnx+1,0:mny+1), &
       alphaf(0:mny+1,0:mnz+1), alphab(0:mny+1,0:mnz+1), frequency(20)

  COMPLEX(DP) :: source(0:mnx+1,0:mny+1,0:mnz+1)

  ! Unknowns epsilon and lambeda
  COMPLEX(DP) :: bet_e(0:mnx+1,0:mny+1,0:mnz+1), bet_w(0:mnx+1,0:mny+1,0:mnz+1), &
       bet_n(0:mnx+1,0:mny+1,0:mnz+1), bet_s(0:mnx+1,0:mny+1,0:mnz+1), &
       lam_ex(0:mnx+1,0:mny+1,0:mnz+1), lam_ez(0:mnx+1,0:mny+1,0:mnz+1), &
       lam_wx(0:mnx+1,0:mny+1,0:mnz+1), lam_wz(0:mnx+1,0:mny+1,0:mnz+1), &
       lam_nx(0:mnx+1,0:mny+1,0:mnz+1), lam_ny(0:mnx+1,0:mny+1,0:mnz+1), &
       lam_sx(0:mnx+1,0:mny+1,0:mnz+1), lam_sy(0:mnx+1,0:mny+1,0:mnz+1)

  COMPLEX(DP) :: e_ex(0:mnx+1,0:mny+1,0:mnz+1), &
       e_wx(0:mnx+1,0:mny+1,0:mnz+1), e_nx(0:mnx+1,0:mny+1,0:mnz+1), &
       e_sx(0:mnx+1,0:mny+1,0:mnz+1), e_fy(0:mnx+1,0:mny+1,0:mnz+1), &
       e_by(0:mnx+1,0:mny+1,0:mnz+1), e_ny(0:mnx+1,0:mny+1,0:mnz+1), &
       e_sy(0:mnx+1,0:mny+1,0:mnz+1), e_ez(0:mnx+1,0:mny+1,0:mnz+1), &
       e_wz(0:mnx+1,0:mny+1,0:mnz+1), e_fz(0:mnx+1,0:mny+1,0:mnz+1), &
       e_bz(0:mnx+1,0:mny+1,0:mnz+1)


  COMPLEX(DP) :: eo_ex(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_wx(0:mnx+1,0:mny+1,0:mnz+1), eo_nx(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_sx(0:mnx+1,0:mny+1,0:mnz+1), eo_fy(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_by(0:mnx+1,0:mny+1,0:mnz+1), eo_ny(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_sy(0:mnx+1,0:mny+1,0:mnz+1), eo_ez(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_wz(0:mnx+1,0:mny+1,0:mnz+1), eo_fz(0:mnx+1,0:mny+1,0:mnz+1), &
       eo_bz(0:mnx+1,0:mny+1,0:mnz+1)

  COMPLEX(DP) :: lamo_ex(0:mnx+1,0:mny+1,0:mnz+1), &
       lamo_ez(0:mnx+1,0:mny+1,0:mnz+1), lamo_wx(0:mnx+1,0:mny+1,0:mnz+1), &
       lamo_wz(0:mnx+1,0:mny+1,0:mnz+1), lamo_nx(0:mnx+1,0:mny+1,0:mnz+1), &
       lamo_ny(0:mnx+1,0:mny+1,0:mnz+1), lamo_sx(0:mnx+1,0:mny+1,0:mnz+1), &
       lamo_sy(0:mnx+1,0:mny+1,0:mnz+1)




  ! The paralel part
  
  INTEGER :: procid, blockid, vblockid(2,nsub_proc), blockty, ncontrol, nsize, nmax

  INTEGER :: nprocs

  INTEGER ::  freq_host,p_id

  INTEGER :: COLOR_COMM, COLOR_RANK, COLOR_SIZE, COMM_SUBD, COMM_MODE, mode_rank, mode_size

  ! Convergence

  REAL(DP) :: errn, den
  Real(DP) :: local_error
  REAL(DP) :: suml(2), swork(2)

  ! Time

  REAL(DP) :: t0, diff_time, end_time, start_time

 
  ! Send / Receive tensors

  COMPLEX(DP) :: sendvec(1:4*(mnmax+2)*ngx), recvec(1:4*(mnmax+2)*ngx)
  
END MODULE commons_3d

! mumps_handler.f90
MODULE mumps_handler
  USE mpi
  IMPLICIT NONE
  INCLUDE 'zmumps_struc.h'

  TYPE(ZMUMPS_STRUC), SAVE :: mumps_par
END MODULE mumps_handler







  
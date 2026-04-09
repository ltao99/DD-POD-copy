subroutine assemble_Krb_from_coo
  USE sizes
  USE commons_3d

  ! Allocate rows, columns and values in csr format
  IF (ALLOCATED(ip_csr)) DEALLOCATE(ip_csr)
  IF (ALLOCATED(ij_csr)) DEALLOCATE(ij_csr)
  IF (ALLOCATED(apj_csr)) DEALLOCATE(apj_csr)           
  ALLOCATE(ip_csr(total_l_nod+1), ij_csr(nnz), apj_csr(nnz))

  ! Convert COO to CSR
  CALL coo_to_csr(total_l_nod, nnz, ip, ij, apj, ip_csr, ij_csr, apj_csr)

  ! Compute KRB
  rank = basis_map(sub, mnodos)
  CALL compute_Krb(total_l_nod, rank, ip_csr, ij_csr, apj_csr, basis(sub, mnodos)%U)
    
  ! Store KRB in COO Format
  CALL store_Krb_coo(rank)

end subroutine assemble_Krb_from_coo


subroutine coo_to_csr(nrows, nnz, ip_coo, ij_coo, apj_coo, row_ptr, col_ind, val_csr)
  USE SIZES
  implicit none
  integer, intent(in) :: nrows, nnz
  integer, intent(in) :: ip_coo(nnz), ij_coo(nnz)
  complex(DP), intent(in) :: apj_coo(nnz)
  integer, intent(out) :: row_ptr(nrows + 1), col_ind(nnz)
  complex(DP), intent(out) :: val_csr(nnz)

  integer :: i, k, idx
  integer, allocatable :: row_counts(:), counter(:)
  integer, allocatable :: sort_order(:)
  integer :: row

  ! Allocate helper arrays
  allocate(row_counts(nrows))
  allocate(counter(nrows+1))
  allocate(sort_order(nnz))

  ! Initialize
  row_counts = 0

  ! Count non-zeros per row
  do k = 1, nnz
     row = ip_coo(k)
     row_counts(row) = row_counts(row) + 1
  end do

  ! Compute row_ptr using cumulative sum
  row_ptr(1) = 0
  do i = 1, nrows
     row_ptr(i+1) = row_ptr(i) + row_counts(i)
  end do

  ! Initialize counter for filling col_ind and val_csr
  counter = row_ptr

  ! Fill CSR structure
  do k = 1, nnz
     row = ip_coo(k)
     idx = counter(row) + 1
     col_ind(idx) = ij_coo(k)
     val_csr(idx) = apj_coo(k)
     counter(row) = idx
  end do

  row_ptr = row_ptr + 1
  
  ! Deallocate helper arrays
  deallocate(row_counts, counter, sort_order)

end subroutine coo_to_csr


subroutine compute_Krb(n, r, row_ptr, col_ind, val_csr, U)
  USE SIZES
  USE commons_3d
  implicit none

  integer, intent(in) :: n, r
  integer, intent(in) :: row_ptr(1:n+1), col_ind(nnz)
  complex(DP), intent(in) :: val_csr(nnz)
  complex(DP), intent(in) :: U(n, r)

  complex(DP) :: K_reduced(r, r)
  complex(DP), allocatable :: KU(:,:)
  character(1) :: transa, matdescra(6)
  complex(DP), parameter :: alpha = (1.0_dp, 0.0_dp), beta = (0.0_dp, 0.0_dp)
  integer :: ldb, ldc

  ! These arrays are needed for mkl_zcsrmm
  integer, allocatable :: pointerB(:), pointerE(:)

  ! Allocate result buffer
  allocate(KU(n, r))
  KU = CMPLX(0.0_DP, 0.0_DP, KIND=DP)
  K_reduced = CMPLX(0.0_DP, 0.0_DP, KIND=DP)

  !===============================
  ! Step 1: Compute KU = K * U using MKL ZCSRMM
  !===============================

  transa = 'N'  ! No transpose

  matdescra(1) = 'S'   ! Symmetric
  matdescra(2) = 'U'   ! Upper triangular (only upper + diagonal stored)
  matdescra(3) = 'N'   ! Non-unit diagonal
  matdescra(4) = 'F'   ! Fortran-style indexing (1-based)
  matdescra(5) = '0'   ! Ignored
  matdescra(6) = '0'   ! Ignored

  ldb = n
  ldc = n

  allocate(pointerB(n), pointerE(n))
  pointerB = row_ptr(1:n)
  pointerE = row_ptr(2:n+1)

  call mkl_zcsrmm(transa, n, r, n, alpha, matdescra,        &
                  val_csr, col_ind, pointerB, pointerE,     &
                  U, ldb, beta, KU, ldc)

  deallocate(pointerB, pointerE)

  !===============================
  ! Step 2: Compute K_reduced = U^H * KU
  !===============================

  call zgemm('N', 'N', r, r, n, alpha, basis(sub, mnodos)%UH, r, &
      KU, n, beta, K_reduced, r)


  !===============================
  ! Store result
  !===============================
  if (associated(reduced_stiffness(sub, mnodos)%Krb)) then
    deallocate(reduced_stiffness(sub, mnodos)%Krb)
  end if
  allocate(reduced_stiffness(sub, mnodos)%Krb(r, r))
  reduced_stiffness(sub, mnodos)%Krb = K_reduced

  deallocate(KU)

end subroutine compute_Krb




subroutine store_Krb_coo(r)
  USE SIZES
  USE commons_3d
  USE globalvars

  implicit none
  integer, intent(in) :: r
  integer :: i, j, count
  complex(DP), pointer :: Krb(:,:)

  ! Local COO arrays
  integer, allocatable :: ip_rb_local(:), ij_rb_local(:)
  complex(DP), allocatable :: apj_rb_local(:)

  ! Sorting support
  integer, allocatable :: perm(:)
  integer :: p, a, b

  ! Point to the reduced stiffness matrix
  Krb => reduced_stiffness(sub, mnodos)%Krb

  ! Count non-zero entries (entire matrix)
  nnz_rb = 0
  do i = 1, r
    do j = 1, r
      if (abs(Krb(i,j)) > ACCURMACHINE) then
        nnz_rb = nnz_rb + 1
      end if
    end do
  end do

  ! Allocate arrays
  if (allocated(ip_rb_local)) deallocate(ip_rb_local)
  if (allocated(ij_rb_local)) deallocate(ij_rb_local)
  if (allocated(apj_rb_local)) deallocate(apj_rb_local)
  if (allocated(perm)) deallocate(perm)
  allocate(ip_rb_local(nnz_rb), ij_rb_local(nnz_rb), apj_rb_local(nnz_rb), perm(nnz_rb))

  ! Fill COO
  count = 0
  do i = 1, r
    do j = 1, r
      if (abs(Krb(i,j)) > ACCURMACHINE) then
        count = count + 1
        ip_rb_local(count) = i
        ij_rb_local(count) = j
        apj_rb_local(count) = Krb(i,j)
        perm(count) = count
      end if
    end do
  end do

  ! Sort perm by (ip_rb, ij_rb) using insertion sort
  do i = 2, nnz_rb
    p = perm(i)
    a = ip_rb_local(p)
    b = ij_rb_local(p)
    j = i - 1
    do while (j >= 1 .and. (ip_rb_local(perm(j)) > a .or. (ip_rb_local(perm(j)) == a .and. ij_rb_local(perm(j)) > b)))
      perm(j+1) = perm(j)
      j = j - 1
    end do
    perm(j+1) = p
  end do

  ! Apply permutation
  reduced_stiffnes_coo(sub, mnodos)%ip_rb = ip_rb_local(perm)
  reduced_stiffnes_coo(sub, mnodos)%ij_rb = ij_rb_local(perm)
  reduced_stiffnes_coo(sub, mnodos)%apj_rb = apj_rb_local(perm)

  ! Deallocate temporaries
  deallocate(ip_rb_local, ij_rb_local, apj_rb_local, perm)
end subroutine store_Krb_coo

!--------------------------------------------------------------------
subroutine prepare_rhs_local_space(r)
   use sizes
   use commons_3d
   use globalvars
   implicit none
   integer, intent(in) :: r


    ! Allocate reduced RHS vector
    if (allocated(bglobal_rb)) deallocate(bglobal_rb)
    allocate(bglobal_rb(rank))

    ! Allocate reduced solution vector
    if (allocated(sol_MUMPS_rb)) deallocate(sol_MUMPS_rb)
    allocate(sol_MUMPS_rb(rank))

    if (use_pod) then
      bglobal_rb = bglobal_rb_constant + bglobal_rb_dd
   else
   ! Initialize bglobal arrays
      bglobal_rb = bglobal_constant + bglobal_dd
   end if
end subroutine prepare_rhs_local_space
!--------------------------------------------------------------------

subroutine recover_full_solution(r)
   use sizes
   use commons_3d
   use globalvars
   implicit none
   integer, intent(in) :: r
   integer :: i, j, dof
   complex(dp), allocatable :: sol_MUMPS_interface(:)

   IF (use_pod) then
      if (allocated(sol_MUMPS_interface)) then
      deallocate(sol_MUMPS_interface)
      end if
      allocate(sol_MUMPS_interface(ndofs_interface))
      sol_MUMPS = (0.0_dp, 0.0_dp) 

      call zgemv('N', ndofs_interface, r, (1.0_dp, 0.0_dp),         &
           basis(sub, mnodos)%U_interf, ndofs_interface,     &
           sol_MUMPS_rb, 1,                                  &
           (0.0_dp, 0.0_dp), sol_MUMPS_interface, 1)

      DO j = 1, ndofs_interface
         sol_MUMPS(dofs_interface(j)) = sol_MUMPS_interface(j)
      END DO
    ELSE
      sol_MUMPS = sol_MUMPS_rb
    ENDIF
end subroutine recover_full_solution
!--------------------------------------------------------------------


  subroutine mult_K_xmean(n, nnz, row_ptr, col_ind, val_csr, x_mean, KX_mean)
  USE SIZES
  implicit none

  integer, intent(in) :: n,nnz
  integer, intent(in) :: row_ptr(n + 1)
  integer, intent(in) :: col_ind(nnz)
  complex(dp), intent(in) :: val_csr(nnz)
  complex(dp), intent(in) :: x_mean(n)
  complex(dp), intent(out) :: KX_mean(n)

  integer :: i, j, k
  complex(dp) :: val

  ! Matrix-vector multiplication exploiting symmetry
  do i = 1, n
     do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_ind(k)
        val = val_csr(k)

        Kx_mean(i) = Kx_mean(i) + val * x_mean(j)
        if (j /= i) then
           Kx_mean(j) = Kx_mean(j) + val * x_mean(i)
        end if
     end do
  end do

end subroutine mult_K_xmean



 subroutine correct_rhs(n, F, KX, F_cor)
    USE sizes
    implicit none

    integer, intent(in) :: n
    complex(DP), intent(in) :: F(n)
    complex(DP), intent(in) :: KX(n)
    complex(DP), intent(out) :: F_cor(n)

    F_cor = F - KX

  end subroutine correct_rhs


  
  
  
  
  
  
  
  
SUBROUTINE prepare_pod_bases
  USE sizes
  USE commons_3d
  USE globalvars

  DO sub = 1,nsub_proc
        CALL compute_and_save_pod_basis
        IF (use_pod) THEN
          CALL register_pod_basis
        ELSE
          basis_map(sub, mnodos) = total_l_nod
        ENDIF
  ENDDO

END SUBROUTINE prepare_pod_bases


SUBROUTINE compute_and_save_pod_basis
  USE sizes
  USE globalvars
  USE commons_3d
  USE, INTRINSIC :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER :: i, j, m, info, lda, ldu, ldvt, lwork, unit_write, nmodes
  CHARACTER(256) :: filename_raw, filename_out
  REAL(DP), ALLOCATABLE :: raw(:,:)
  COMPLEX(DP), ALLOCATABLE :: rb_raw(:,:), VT(:,:), work(:), U_aux(:,:)
  REAL(DP), ALLOCATABLE :: S(:), rwork(:)
  REAL(DP) :: tol, eigen1, ur, ui
  INTEGER :: NROWS, NCOLS, line_count, nan_count
  COMPLEX(DP), ALLOCATABLE :: U_trunc(:,:)
  CHARACTER(200) :: base_path
  LOGICAL :: file_exists

  IF (.NOT.(use_pod)) RETURN

  ! Dimensions
  NROWS = total_l_nod
  NCOLS = Nsnapshots
  base_path = pod_path

  ! Filenames
  WRITE(filename_raw, '(A,A,A,A,I0,A)') TRIM(base_path), 'sol_', modo, '_', blockid, '.txt'
  WRITE(filename_out, '(A,A,A,A,I0,A)') TRIM(base_path), 'sol_', modo, '_', blockid, 'svd.txt'

  ! === If SVD already exists, skip everything ===
  INQUIRE(FILE=filename_out, EXIST=file_exists)
  IF (file_exists) THEN
      PRINT *, "POD basis exists. Skipping SVD for:", TRIM(filename_out)
      RETURN
  END IF

  ! === Read snapshot matrix ===
  ALLOCATE(raw(NROWS * NCOLS, 2), rb_raw(NROWS, NCOLS))

  OPEN(10, FILE=filename_raw, STATUS='OLD', ACTION='READ')
  DO i = 1, NROWS * NCOLS
    READ(10, *) raw(i,1), raw(i,2)
  END DO
  CLOSE(10)

  m = 1
  DO j = 1, NCOLS
    DO i = 1, NROWS
      rb_raw(i,j) = raw(m,1) + aim * raw(m,2)
      m = m + 1
    END DO
  END DO
  DEALLOCATE(raw)

  ! === Compute SVD ===
  lda = NROWS
  ldu = NROWS
  ldvt = NCOLS

  ALLOCATE(S(MIN(NROWS, NCOLS)), U_aux(NROWS, MIN(NROWS, NCOLS)))
  ALLOCATE(VT(MIN(NROWS, NCOLS), NCOLS), work(1), rwork(5 * MIN(NROWS, NCOLS)))

  CALL ZGESVD('S', 'S', NROWS, NCOLS, rb_raw, lda, S, U_aux, ldu, VT, ldvt, work, -1, rwork, info)
  lwork = INT(REAL(work(1)))
  DEALLOCATE(work)
  ALLOCATE(work(lwork))

  CALL ZGESVD('S', 'S', NROWS, NCOLS, rb_raw, lda, S, U_aux, ldu, VT, ldvt, work, lwork, rwork, info)
  IF (info /= 0) THEN
    PRINT *, 'SVD failed with info =', info
    STOP
  END IF

  eigen1 = S(1)
  tol = tol_svd * eigen1
  nmodes = COUNT(S >= tol)

  PRINT *, "Computed rank:", nmodes

  ALLOCATE(U_trunc(NROWS, nmodes))
  U_trunc = U_aux(:, 1:nmodes)

  DEALLOCATE(U_aux, S, VT, work, rwork)

  ! === Write SVD basis file ===
  unit_write = 20
  OPEN(unit_write, FILE=filename_out, STATUS='REPLACE', ACTION='WRITE')
  WRITE(unit_write, '(I8)') nmodes

  line_count = 1
  nan_count = 0

  DO j = 1, nmodes
    DO i = 1, NROWS
      ur = REAL(U_trunc(i,j))
      ui = AIMAG(U_trunc(i,j))
      WRITE(unit_write, '(E30.20, 1X, E30.20, 1X)') ur, ui
      line_count = line_count + 1
    END DO
  END DO

  CLOSE(unit_write)

  DEALLOCATE(U_trunc, rb_raw)

  PRINT *, 'Wrote POD basis to "', TRIM(filename_out), '"'

END SUBROUTINE compute_and_save_pod_basis


SUBROUTINE register_pod_basis
  USE sizes
  USE globalvars
  USE commons_3d
  IMPLICIT NONE

  INTEGER :: i, j, m, info, unit_check, dof
  CHARACTER(256) :: filename_out
  CHARACTER(200) :: base_path
  INTEGER :: NROWS, NCOLS, nmodes

  REAL(DP), ALLOCATABLE :: raw(:,:)
  COMPLEX(DP), ALLOCATABLE :: U_aux(:,:), UH_aux(:,:)

  ! === Dimensions ===
  NROWS = total_l_nod
  base_path = pod_path

  ! Construct filename
  WRITE(filename_out, '(A,A,A,A,I0,A)') TRIM(base_path), 'sol_', modo, '_', blockid, 'svd2.txt'

  ! === Open basis file ===
  unit_check = 11
  OPEN(unit_check, FILE=filename_out, STATUS='OLD', ACTION='READ', IOSTAT=info)
  IF (info /= 0) THEN
    PRINT *, 'ERROR: Cannot open POD basis file: ', filename_out
    STOP
  END IF

  ! === NEW: Read number of modes from the first line ===
  READ(unit_check, *) nmodes
  basis_map(sub, mnodos) = nmodes
  NCOLS = nmodes

  ! === Allocate U and raw arrays ===
  ALLOCATE(U_aux(NROWS, NCOLS))
  ALLOCATE(raw(NROWS * NCOLS, 2))

  ! === Read POD basis matrix ===
  DO m = 1, NROWS * NCOLS
    READ(unit_check, *) raw(m,1), raw(m,2)
  END DO

  m = 1
  DO j = 1, NCOLS
    DO i = 1, NROWS
      U_aux(i,j) = raw(m,1) + aim * raw(m,2)
      m = m + 1
    END DO
  END DO

  ! === Allocate and store U and UH ===
  allocate(basis(sub, mnodos)%U(NROWS, NCOLS))
  allocate(UH_aux(NCOLS, NROWS))
  allocate(basis(sub, mnodos)%UH(NCOLS, NROWS))

  basis(sub, mnodos)%U  = U_aux
  UH_aux                = transpose(conjg(U_aux))
  basis(sub, mnodos)%UH = UH_aux

  ! === Allocate and extract interface basis ===
  allocate(basis(sub, mnodos)%U_interf(ndofs_interface, NCOLS))
  allocate(basis(sub, mnodos)%UH_interf(NCOLS, ndofs_interface))

  DO j = 1, ndofs_interface
    dof = dofs_interface(j)
    DO i = 1, NCOLS
      basis(sub, mnodos)%U_interf(j, i)  = U_aux(dof, i)
      basis(sub, mnodos)%UH_interf(i, j) = conjg(U_aux(dof, i))
    END DO
  END DO

  ! Cleanup
  CLOSE(unit_check)
  DEALLOCATE(raw, U_aux, UH_aux)

END SUBROUTINE register_pod_basis



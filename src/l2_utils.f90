SUBROUTINE write_l2_iteration_fulltag(l2val, mode, index)
  USE sizes
  USE commons_3d
  USE globalvars
  USE string_utils
  IMPLICIT NONE

  ! Arguments
  CHARACTER(*), INTENT(IN) :: mode
  REAL(DP), INTENT(IN) :: l2val
  INTEGER, INTENT(IN) :: index

  ! Locals
  CHARACTER(LEN=20) :: c_tol, c_relx, c_penal, c_conv
  CHARACTER(LEN=500) :: suffix, filename
  INTEGER :: unit_id, ios

  ! === Format parameters ===
  WRITE(c_tol, '(ES10.1)') tol_svd
  WRITE(c_relx, '(F6.3)') relx
  WRITE(c_penal, '(F6.3)') penal
  WRITE(c_conv, '(ES10.1)') convergence_tol
  
  ! === Replace 'E' with 'e' for filename compatibility ===
  CALL replace_char(c_tol, 'E', 'e')
  CALL replace_char(c_conv, 'E', 'e')

  ! === Construct suffix ===
  suffix = ""
  IF (use_pod) THEN
  suffix = '_NS'    // TRIM(itoa(Nsnapshots)) // &
           '_tol'   // TRIM(ADJUSTL(c_tol))   // &
           '_relx'  // TRIM(ADJUSTL(c_relx))  // &
           '_penal' // TRIM(ADJUSTL(c_penal)) // &
           '_maxit' // TRIM(itoa(maxiter))    // &
           '_conv'  // TRIM(ADJUSTL(c_conv))
  ENDIF

  IF (use_pod) THEN
    suffix = TRIM(suffix)//'_svd'
  ELSE
    suffix = TRIM(suffix)//'_no_svd'
  END IF

  ! === Compose full filename ===
  filename = 'outputs/l2_0.001Hz_' // TRIM(mode) // TRIM(suffix) // '.txt'

  ! === Open file in append mode ===
  OPEN(NEWUNIT=unit_id, FILE=filename, STATUS='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
  IF (ios /= 0) THEN
    PRINT *, 'Error opening file: ', TRIM(filename)
    STOP
  END IF

  ! === Write iteration info ===
  WRITE(unit_id, '(I6,1X,I6,1X,ES22.14)') index, sub, l2val

  CLOSE(unit_id)
END SUBROUTINE write_l2_iteration_fulltag


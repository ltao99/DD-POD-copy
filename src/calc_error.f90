SUBROUTINE calc_error(errno, deno)

  USE sizes
  USE commons_3d
  IMPLICIT NONE

  REAL(DP) :: errno, deno
  INTEGER :: j, k, l, y, z
  COMPLEX(DP) :: tmp

  errno = 0.d0
  deno = 0.d0

  DO z = 1, size(interf_z)
       l = z
    DO k = 1, mny
      DO j = 1, mnx
        errno = errno + &
            (REAL(eo_nx(j,k,l) - e_nx(j,k,l))**2 + AIMAG(eo_nx(j,k,l) - e_nx(j,k,l))**2 + &
             REAL(eo_sx(j,k,l) - e_sx(j,k,l))**2 + AIMAG(eo_sx(j,k,l) - e_sx(j,k,l))**2 + &
             REAL(eo_ny(j,k,l) - e_ny(j,k,l))**2 + AIMAG(eo_ny(j,k,l) - e_ny(j,k,l))**2 + &
             REAL(eo_sy(j,k,l) - e_sy(j,k,l))**2 + AIMAG(eo_sy(j,k,l) - e_sy(j,k,l))**2)

        deno = deno + &
            (REAL(eo_nx(j,k,l) + e_nx(j,k,l))**2 + AIMAG(eo_nx(j,k,l) + e_nx(j,k,l))**2 + &
             REAL(eo_sx(j,k,l) + e_sx(j,k,l))**2 + AIMAG(eo_sx(j,k,l) + e_sx(j,k,l))**2 + &
             REAL(eo_ny(j,k,l) + e_ny(j,k,l))**2 + AIMAG(eo_ny(j,k,l) + e_ny(j,k,l))**2 + &
             REAL(eo_sy(j,k,l) + e_sy(j,k,l))**2 + AIMAG(eo_sy(j,k,l) + e_sy(j,k,l))**2)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE calc_error













    
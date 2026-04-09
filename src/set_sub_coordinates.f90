MODULE initialization
    USE sizes
    USE commons_3d
  IMPLICIT NONE
CONTAINS

  SUBROUTINE set_subdomain_coordinates

    IF (nsub_proc .GT. 1) THEN
      DO sub = 1, nsub_proc
        vblockid(1, sub) = MOD(sub - 1, nsy)
        vblockid(2, sub) = (sub - 1) / nsy
      END DO
      auxz = nsz - 1
      FLAG_globy = 1
      FLAG_globz = 1
    ELSE
      vblockid(1,1) = MOD(blockid, nsy)
      vblockid(2,1) = blockid / nsy
      auxz = 0
      FLAG_globy = 0
      FLAG_globz = 0
      ALLOCATE(basis(nsub_proc, 2))
    END IF
  END SUBROUTINE set_subdomain_coordinates

END MODULE initialization

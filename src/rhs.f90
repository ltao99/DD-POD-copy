SUBROUTINE rhs_glob_constant(r)
   USE sizes
   USE globalvars
   USE commons_3d
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: r
   INTEGER :: j, k, l, i, dy, dz
   INTEGER :: igx, igy, igz

   ! Deallocate and allocate bglobal_constant
   IF (ALLOCATED(bglobal_constant)) THEN
      DEALLOCATE(bglobal_constant)
   END IF
   ALLOCATE(bglobal_constant(total_l_nod))
   bglobal_constant = (0.0D0, 0.0D0)

   ! Assemble bglobal_constant from source values
   DO dz = 1, nz
      igz = dz + ((nsz - 1) - vblockid(2, sub)) * nz
      DO dy = 1, ny
         igy = dy + vblockid(1, sub) * ny
         DO j = 1, nx
            igx = j
            i = matnumloc(j, dy, dz)
            k = dy + FLAG_globy * (vblockid(1, sub) * ny)
            l = dz + FLAG_globz * ((nsz - 1 - vblockid(2, sub)) * nz)

            IF (mnodos .EQ. 1) THEN
               bglobal_constant(lc_nod(i, 5)) = bglobal_constant(lc_nod(i, 5)) + source(j, k, l)  ! by
               bglobal_constant(lc_nod(i, 6)) = bglobal_constant(lc_nod(i, 6)) + source(j, k, l)  ! ez
               bglobal_constant(lc_nod(i, 9)) = bglobal_constant(lc_nod(i, 9)) + source(j, k, l)  ! sy
               bglobal_constant(lc_nod(i, 10)) = bglobal_constant(lc_nod(i, 10)) + source(j, k, l) ! ny
               bglobal_constant(lc_nod(i, 11)) = bglobal_constant(lc_nod(i, 11))                   !

            ELSE
               bglobal_constant(lc_nod(i, 1)) = bglobal_constant(lc_nod(i, 1)) + source(j, k, l)  ! wx
               bglobal_constant(lc_nod(i, 2)) = bglobal_constant(lc_nod(i, 2)) + source(j, k, l)  ! ex
               bglobal_constant(lc_nod(i, 3)) = bglobal_constant(lc_nod(i, 3)) + source(j, k, l)  ! sx
               bglobal_constant(lc_nod(i, 4)) = bglobal_constant(lc_nod(i, 4)) + source(j, k, l)  ! nx
               bglobal_constant(lc_nod(i, 5)) = bglobal_constant(lc_nod(i, 5))                    !by
               bglobal_constant(lc_nod(i, 11)) = bglobal_constant(lc_nod(i, 11))
            END IF
         END DO
      END DO
   END DO

   ! POD projection if enabled
   if (use_pod) then
   if (allocated(bglobal_rb_constant)) then
      deallocate(bglobal_rb_constant)
   end if
   allocate(bglobal_rb_constant(r))

   call zgemv('N', r, total_l_nod, (1.0_dp, 0.0_dp),          &
              basis(sub, mnodos)%UH, r,                      &
              bglobal_constant, 1,                           &
              (0.0_dp, 0.0_dp), bglobal_rb_constant, 1)
end if


END SUBROUTINE rhs_glob_constant
 


SUBROUTINE rhs_glob_dd(r)

   USE sizes
   USE commons_3d
   USE globalvars

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: r
   INTEGER :: j, k, l, i, dy, dz, y, z,dof
   INTEGER :: igx, igz, igy
   COMPLEX(DP) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, &
                  tmp7, tmp8, tmp9, tmp10, tmp11, tmp12

IF (ALLOCATED(bglobal_dd)) THEN
   DEALLOCATE(bglobal_dd)
END IF

   ALLOCATE(bglobal_dd(total_l_nod))

   DO j = 1, total_l_nod
      bglobal_dd(j) = (0.d0, 0.d0)
   END DO

! Z-direction interface contribution
DO z = 1, SIZE(interf_z)
   dz = interf_z(z)
   igz = dz + ((nsz - 1) - vblockid(2, sub)) * nz

   DO dy = 1, ny
      igy = dy + vblockid(1, sub) * ny

      DO j = 1, nx
         igx = j
         i = matnumloc(j, dy, dz)
         k = dy + FLAG_globy * (vblockid(1, sub) * ny)
         l = dz + FLAG_globz * (((nsz - 1) - vblockid(2, sub)) * nz)

         ! Main logic
         hxz = hx(igx) * hz(igz)
         hxy = hx(igx) * hy(igy)
         hyz = hy(igy) * hz(igz)

         tmp3 = (1.0_DP - deltabs(l)) * hxy * (bet_s(j, k, l) * eo_nx(j, k, l - 1) - lamo_nx(j, k, l - 1))
         tmp4 = (1.0_DP - deltabn(l)) * hxy * (bet_n(j, k, l) * eo_sx(j, k, l + 1) - lamo_sx(j, k, l + 1))
         tmp9 = (1.0_DP - deltabs(l)) * hxy * (bet_s(j, k, l) * eo_ny(j, k, l - 1) - lamo_ny(j, k, l - 1))
         tmp10 = (1.0_DP - deltabn(l)) * hxy * (bet_n(j, k, l) * eo_sy(j, k, l + 1) - lamo_sy(j, k, l + 1))

         ! Safe writes
         bglobal_dd(lc_nod(i, 3)) = tmp3  ! sx
         bglobal_dd(lc_nod(i, 4)) = tmp4  ! nx
         bglobal_dd(lc_nod(i, 9)) = tmp9  ! sy
         bglobal_dd(lc_nod(i, 10)) = tmp10 ! ny
      END DO
   END DO
END DO

IF (use_pod) THEN
   IF (ALLOCATED(bglobal_rb_dd)) THEN
      DEALLOCATE(bglobal_rb_dd)
   END IF
   ALLOCATE(bglobal_rb_dd(r))
   bglobal_rb_dd = (0.d0, 0.d0)

call zgemv('N', r, ndofs_interface, (1.0_dp, 0.0_dp),       &
           basis(sub, mnodos)%UH_interf, r,                & ! UH_interf(j, i) = conjg(U(i,j))
           bglobal_dd(dofs_interface), 1,                  & ! only interface DoFs
           (0.0_dp, 0.0_dp), bglobal_rb_dd, 1)

END IF


END SUBROUTINE rhs_glob_dd




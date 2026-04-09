

!***********************************************************  


 SUBROUTINE change_sol_surf(sol, izsup_loc)
     USE sizes
     USE commons_3d
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: izsup_loc
     INTEGER :: dy, dz,j, k, l, i,z
     complex(dp) :: sol(1:total_l_nod)
   
     DO dz = izsup_loc, izsup_loc
      l = dz + FLAG_globz*(((nsz-1)-vblockid(2,sub))*nz)
        DO dy = 1, ny
         k = dy + FLAG_globy*(vblockid(1,sub)*ny)
           DO j = 1,nx
              i = matnumloc(j,dy,dz)
              e_wx(j, k, l) = sol(lc_nod(i, 1))
              e_ex(j, k, l) = sol(lc_nod(i, 2))
              e_sx(j, k, l) = sol(lc_nod(i, 3))
              e_nx(j, k, l) = sol(lc_nod(i, 4))
              e_by(j, k, l) = sol(lc_nod(i, 5)) 
              e_fy(j, k, l) = sol(lc_nod(i, 6))
              e_wz(j, k, l) = sol(lc_nod(i, 7))
              e_ez(j, k, l) = sol(lc_nod(i, 8))
              e_sy(j, k, l) = sol(lc_nod(i, 9))
              e_ny(j, k, l) = sol(lc_nod(i, 10))
              e_bz(j, k, l) = sol(lc_nod(i, 11))
              e_fz(j, k, l) = sol(lc_nod(i, 12))
           END DO
        END DO
     END DO

    

     RETURN
   END SUBROUTINE change_sol_surf

!***********************************************************  
   
 SUBROUTINE change_sol(sol)
     USE sizes
     USE commons_3d
     IMPLICIT NONE
     
     INTEGER :: dy, dz,j, k, l, i, y, z
     complex(dp) :: sol(1:total_l_nod)




    DO z = 1, SIZE(interf_z)
       dz = interf_z(z)
         l = dz +  FLAG_globz*(((nsz-1)-vblockid(2,sub))*nz)
         DO dy = 1, ny
            k = dy + FLAG_globy*(vblockid(1,sub)*ny)
            DO j = 1,nx
               i = matnumloc(j,dy,dz)

               e_sx(j, k, l) = sol(lc_nod(i, 3))
               e_nx(j, k, l) = sol(lc_nod(i, 4))
               e_sy(j, k, l) = sol(lc_nod(i, 9))
               e_ny(j, k, l) = sol(lc_nod(i, 10))
           END DO
        END DO
     END DO    

      RETURN
   END SUBROUTINE change_sol


!***********************************************************  



   SUBROUTINE update_e

     USE sizes 
     USE globalvars
     USE commons_3d
     IMPLICIT NONE
    
     INTEGER :: j, k, l, y, z
     
     DO z = 1, SIZE(interf_z)
      l = interf_z(z)
      DO k = 1, mny
         DO j = 1, mnx
               e_nx(j,k,l)=(1-relx)*e_nx(j,k,l)+relx*eo_nx(j,k,l)
               e_sx(j,k,l)=(1-relx)*e_sx(j,k,l)+relx*eo_sx(j,k,l)
               e_ny(j,k,l)=(1-relx)*e_ny(j,k,l)+relx*eo_ny(j,k,l)   
               e_sy(j,k,l)=(1-relx)*e_sy(j,k,l)+relx*eo_sy(j,k,l)
             END DO
          END DO
       END DO

   END SUBROUTINE update_e



   SUBROUTINE update_l
     USE sizes
     USE globalvars
     USE commons_3d
     IMPLICIT NONE
   
     INTEGER :: j, k, l, kp, km, lp, lm, y, z
   
      DO z = 1, SIZE(interf_z)
         l = interf_z(z)
         DO k = 1, mny
            DO j = 1, mnx

              kp = k + 1
              km = k - 1
   
              lp = l + 1
              lm = l - 1
   
              lam_nx(j, k, l) = -lamo_sx(j, k, lp) + bet_n(j, k, l) * &
                   (eo_sx(j, k, lp) - e_nx(j, k, l))
              lam_nx(j, k, l) = (1-relx) * lam_nx(j, k, l) + relx * lamo_nx(j, k, l)
   
              lam_ny(j, k, l) = -lamo_sy(j, k, lp) + bet_n(j, k, l) * &
                   (eo_sy(j, k, lp) - e_ny(j, k, l))
              lam_ny(j, k, l) = (1-relx) * lam_ny(j, k, l) + relx * lamo_ny(j, k, l)
   
              lam_sx(j, k, l) = -lamo_nx(j, k, lm) + bet_s(j, k, l) * &
                   (eo_nx(j, k, lm) - e_sx(j, k, l))
              lam_sx(j, k, l) = (1-relx) * lam_sx(j, k, l) + relx * lamo_sx(j, k, l)
  
              lam_sy(j, k, l) = -lamo_ny(j, k, lm) + bet_s(j, k, l) * &
                   (eo_ny(j, k, lm) - e_sy(j, k, l))
              lam_sy(j, k, l) = (1-relx) * lam_sy(j, k, l) + relx * lamo_sy(j, k, l)
   
           END DO
        END DO
     END DO
   
     RETURN
   END SUBROUTINE update_l

          
!***********************************************************  

   SUBROUTINE update_eo_lo
      USE sizes
      USE commons_3d
      IMPLICIT NONE
   
      INTEGER :: j, k, l, y, z
   
      DO z = 1, SIZE(interf_z)
         l = interf_z(z)
         DO k = 1, mny
            DO j = 1, mnx
   
               lamo_nx(j, k, l) = lam_nx(j, k, l)
               lamo_sx(j, k, l) = lam_sx(j, k, l)
               lamo_ny(j, k, l) = lam_ny(j, k, l)
               lamo_sy(j, k, l) = lam_sy(j, k, l)
      
               eo_nx(j, k, l) = e_nx(j, k, l)
               eo_sx(j, k, l) = e_sx(j, k, l)
               eo_ny(j, k, l) = e_ny(j, k, l)
               eo_sy(j, k, l) = e_sy(j, k, l)
            END DO
         END DO
      END DO
   
   END SUBROUTINE update_eo_lo




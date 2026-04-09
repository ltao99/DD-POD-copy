       SUBROUTINE inguess
       
         USE sizes
         USE commons_3d
         IMPLICIT NONE
         
         INTEGER :: j, k, l

         DO j= 0, mnx+1
            DO k= 0, mny+1
               DO l= 0, mnz+1

                  eo_ex(j,k,l)= (0.d0,0.d0)
                  eo_wx(j,k,l)= (0.d0,0.d0)
                  eo_nx(j,k,l)=  (0.d0,0.d0)
                  eo_sx(j,k,l)=  (0.d0,0.d0)
                  
                  eo_fy(j,k,l)= (0.d0,0.d0)
                  eo_by(j,k,l)= (0.d0,0.d0)
                  eo_ny(j,k,l)=  (0.d0,0.d0)
                  eo_sy(j,k,l)=  (0.d0,0.d0)
            
                  eo_ez(j,k,l)= (0.d0,0.d0)
                  eo_wz(j,k,l)= (0.d0,0.d0)
                  eo_fz(j,k,l)= (0.d0,0.d0)
                  eo_bz(j,k,l)= (0.d0,0.d0)

                  e_ex(j,k,l)= (0.d0,0.d0)
                  e_wx(j,k,l)= (0.d0,0.d0)
                  e_nx(j,k,l)= (0.d0,0.d0)
                  e_sx(j,k,l)= (0.d0,0.d0)
                  
                  e_fy(j,k,l)= (0.d0,0.d0)
                  e_by(j,k,l)= (0.d0,0.d0)
                  e_ny(j,k,l)= (0.d0,0.d0)
                  e_sy(j,k,l)= (0.d0,0.d0)
            
                  e_ez(j,k,l)= (0.d0,0.d0)
                  e_wz(j,k,l)= (0.d0,0.d0)
                  e_fz(j,k,l)= (0.d0,0.d0)
                  e_bz(j,k,l)= (0.d0,0.d0)
                  
                  lamo_ex(j,k,l)= (0.d0,0.d0)
                  lamo_ez(j,k,l)= (0.d0,0.d0)
                  lamo_wx(j,k,l)= (0.d0,0.d0)
                  lamo_wz(j,k,l)= (0.d0,0.d0)
                  
                  lamo_nx(j,k,l)= (0.d0,0.d0)
                  lamo_ny(j,k,l)= (0.d0,0.d0)
                  lamo_sx(j,k,l)= (0.d0,0.d0)
                  lamo_sy(j,k,l)= (0.d0,0.d0)

               ENDDO
            ENDDO
         ENDDO

         RETURN
       END SUBROUTINE inguess








      SUBROUTINE output
               
        USE sizes 
        USE commons_3d
        USE globalvars
        USE string_utils
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        
        INTEGER :: j,k,ll,m,flag,jj,izsup_loc,orig_group,ierr,grupo_local(nsy-FLAG_globy*(nsy-1)),l,i
        INTEGER :: new_group,new_comm,blockid_loc,k_sup,blockid_surf,nelm,idx,blk
        COMPLEX(DP) :: tmp1,tmp2,ep,hp,ep_s,hp_s,zyx,tmp3,tmp4,zxy
        COMPLEX(DP) :: envio_datos(12*ngx*mny),recibo_datos(12*ngx*ngy)
        REAL(DP) :: ryx,z,fase,rxy,hx_1,hy_1,hz_1
        COMPLEX(DP) :: eg_sy(ngx,ngy),eg_ny(ngx,ngy),eg_fy(ngx,ngy)
        COMPLEX(DP) :: eg_by(ngx,ngy),eg_fz(ngx,ngy),eg_bz(ngx,ngy)
        COMPLEX(DP) :: eg_ez(ngx,ngy),eg_wz(ngx,ngy),eg_ex(ngx,ngy)
        COMPLEX(DP) :: eg_wx(ngx,ngy),eg_nx(ngx,ngy),eg_sx(ngx,ngy)
        COMPLEX(DP) :: ce_x_xy(ngx,ngy),ce_y_xy(ngx,ngy)
        COMPLEX(DP) :: ce_x_yx(ngx,ngy),ce_y_yx(ngx,ngy)
        COMPLEX(DP) :: cm_x_xy(ngx,ngy),cm_y_xy(ngx,ngy)
        COMPLEX(DP) :: cm_x_yx(ngx,ngy),cm_y_yx(ngx,ngy)
        CHARACTER(LEN=50)  :: grid_info
        CHARACTER(8) :: dd
        CHARACTER(LEN=300) :: archivo, archivo1, archivo2, suffix
        CHARACTER(LEN=20) :: c_tol, c_relx, c_penal, c_conv
        CHARACTER(25) :: cfreq
        CHARACTER(80) :: filename
        INTEGER :: igx,igz,igy,indx
        COMPLEX(DP), ALLOCATABLE :: sendbuf(:), recvbuf(:)
        REAL(DP), ALLOCATABLE :: hx_inv(:), hy_inv(:), x_c(:), y_c(:)
        INTEGER :: status(MPI_STATUS_SIZE)

        blockid_surf = CEILING( real(ngz - (izsup - 1), dp) / real(nz, dp) + 0.0000000001 ) - 1

        IF(blockid .EQ. blockid_surf)THEN

        IF (mode_rank .EQ. 1) THEN

        WRITE(c_tol,'(ES10.1)') tol_svd
        WRITE(c_relx,'(F6.3)') relx
        WRITE(c_penal,'(F6.3)') penal
        WRITE(c_conv,'(ES10.1)') convergence_tol
        WRITE(cfreq,'(E12.5)') freq

        CALL replace_char(c_tol, 'E', 'e')
        CALL replace_char(c_conv, 'E', 'e')

        suffix = ""

        IF (use_pod) THEN
          suffix = '_NS'//TRIM(itoa(Nsnapshots))// &
                    '_tol'//TRIM(ADJUSTL(c_tol))// &
                    '_relx'//TRIM(ADJUSTL(c_relx))// &
                    '_penal'//TRIM(ADJUSTL(c_penal))// &
                    '_maxit'//TRIM(itoa(maxiter))// &
                    '_conv'//TRIM(ADJUSTL(c_conv))
          suffix = TRIM(suffix)//'_svd'
          ELSE
          suffix = TRIM(suffix)//'_no_svd'
        ENDIF

        CALL system('mkdir -p outputs')

        indx = INT(1.0_DP / freq)
        WRITE(grid_info, '(A,I0,A,I0,A,I0,A,I0)') '_ngx', ngx, '_ngy', ngy, '_ngz', ngz, '_np', nprocs

        archivo  = 'outputs/fields_for_freq=' // TRIM(cfreq) // 'Hz' // TRIM(grid_info) // TRIM(suffix) // '.txt'
        archivo1 = 'outputs/rxy_y_fase_freq=' // TRIM(cfreq) // 'Hz' // TRIM(grid_info) // TRIM(suffix) // '.txt'
        archivo2 = 'outputs/ryx_y_fase_freq=' // TRIM(cfreq) // 'Hz' // TRIM(grid_info) // TRIM(suffix) // '.txt'

        ! After receive, proceed with file output on this side
        OPEN(unit=100, file=archivo, status='replace')
        OPEN(unit=200+indx, file=archivo1, status='replace')
        OPEN(unit=300+indx, file=archivo2, status='replace')

        ALLOCATE(hx_inv(nx), x_c(nx))
        ALLOCATE(hy_inv(ngy), y_c(ngy))

        DO j = 1, nx
          hx_inv(j) = 1.0_dp / hx(j)
          x_c(j)    = ((xnod(j) + hx(j)/2.0_dp - (xnod(ngx+1) - xnod(1))/2.0_dp) / 1000.0_dp)
        END DO
        DO k = 1, ngy
          hy_inv(k) = 1.0_dp / hy(k)
          y_c(k)    = (ynod(k) + hy(k)/2.0_dp - (ynod(ngy+1) - ynod(1))/2.0_dp) / 1000.0_dp
       END DO
        ENDIF

           IF(nprocs.GT.1)THEN
              IF (mnz .EQ. 1 .AND. nsy .EQ.1) THEN
                 izsup_loc = 1
              ELSE
                 izsup_loc= nz - MOD(ngz-(izsup-1), nz)
              ENDIF
           ELSE
              izsup_loc= izsup-1
           ENDIF

              PRINT*,'la superficie es el nodo global', izsup
              PRINT*,'el valor local de izsup es (#elemento para computo'
              PRINT*,'de los campos en superficie)',izsup_loc
              PRINT*,'el aire ocupa',ngz-(izsup-1),'elementos en la vertical'
           
           IF (blockid .LT. ny) THEN
              IF (use_pod) THEN
                 sol_MUMPS = matmul(basis(1,mnodos)%U, sol_MUMPS_rb)
                 CALL change_sol_surf(sol_MUMPS(1:total_l_nod), izsup_loc)
              ELSE
                 CALL change_sol_surf(sol_MUMPS(1:total_l_nod), izsup_loc)
              ENDIF
           ENDIF
          
            DO k= 1,ngy            
               DO m= 1, ngx
                    eg_sy(m,k)= e_sy(m,k,izsup_loc)
                    eg_ny(m,k)= e_ny(m,k,izsup_loc)
                    eg_fy(m,k)= e_fy(m,k,izsup_loc)
                    eg_by(m,k)= e_by(m,k,izsup_loc)       
                    eg_fz(m,k)= e_fz(m,k,izsup_loc)     
                    eg_bz(m,k)= e_bz(m,k,izsup_loc) 
                    eg_ez(m,k)= e_ez(m,k,izsup_loc)
                    eg_wz(m,k)= e_wz(m,k,izsup_loc)  
                    eg_ex(m,k)= e_ex(m,k,izsup_loc)           
                    eg_wx(m,k)= e_wx(m,k,izsup_loc)       
                    eg_nx(m,k)= e_nx(m,k,izsup_loc)    
                    eg_sx(m,k)= e_sx(m,k,izsup_loc)            
                ENDDO
             ENDDO

                igz  = izsup - 1
                z    = znod(izsup)
                hz_1 = 1.d0 / hz(igz)
                ep_s = ep(z)
                hp_s = hp(z)
                blk = ngx * ngy

              PRINT*,'z de la sup. en output para campos primarios',z

                 DO j= 1, nx
                    hx_1= 1.d0/hx(j)
                    DO k= 1, ngy 
                       hy_1= 1.d0/hy(k)
                       IF (mode_rank .EQ. 0) THEN
            
                         tmp1= eg_ny(j,k)
                         tmp2= .125d0*(-8.d0*coef_i*(hz_1*(eg_ny(j,k)-eg_sy(j,k)) &
                            +hy_1*(-eg_ez(j,k)+eg_wz(j,k)))+ &
                            (-28.d0*coef_i*hz_1*(eg_ny(j,k)+eg_sy(j,k)-eg_fy(j,k)-eg_by(j,k))))
                         ce_y_yx(j,k)= ep_s+tmp1
                         cm_x_yx(j,k)= hp_s+tmp2
                         ce_x_yx(j,k)= eg_nx(j,k)
                         cm_y_yx(j,k)= .125d0*(-8.d0*coef_i*(hz_1*(eg_sx(j,k)-eg_nx(j,k)) &
                            +hx_1*(eg_fz(j,k)-eg_bz(j,k))) &
                            -(-28.d0*coef_i*hz_1*(eg_nx(j,k)+eg_sx(j,k)-eg_ex(j,k)-eg_wx(j,k))))
                     ELSE
                         ce_x_xy(j,k)= eg_nx(j,k)+ep_s
                         ce_y_xy(j,k)= eg_ny(j,k)
                         cm_x_xy(j,k)= .125d0*(-8.d0*coef_i*(hz_1*(eg_ny(j,k)-eg_sy(j,k)) &
                            +hy_1*(-eg_ez(j,k)+eg_wz(j,k)))+ &
                            (-28.d0*coef_i*hz_1*(eg_ny(j,k)+eg_sy(j,k)-eg_fy(j,k)-eg_by(j,k))))
                       
                         tmp4= .125d0*(-8.d0*coef_i*(hz_1*(eg_sx(j,k)-eg_nx(j,k)) &
                            +hx_1*(eg_fz(j,k)-eg_bz(j,k))) &
                            -(-28.d0*coef_i*hz_1*(eg_nx(j,k)+eg_sx(j,k)-eg_ex(j,k)-eg_wx(j,k))))
                         cm_y_xy(j,k)= tmp4-hp_s
                      ENDIF
                    ENDDO
                 ENDDO

                   ! -------- Single packed exchange of 4 YX matrices (TE->TM) --------
                nelm = 4*ngx*ngy
                ALLOCATE(sendbuf(nelm), recvbuf(nelm))
                IF (mode_rank .EQ. 0) THEN
                   ! -------- pack all four arrays at once --------
                k = 0
                DO j = 1, ngy
                   DO i = 1, ngx
                   k = k + 1
                   sendbuf(         k) = ce_y_yx(i,j)
                   sendbuf(   blk + k) = cm_x_yx(i,j)
                   sendbuf(2*blk + k) =  ce_x_yx(i,j)
                   sendbuf(3*blk + k) =  cm_y_yx(i,j)
                   END DO
                ENDDO

               CALL MPI_Send(sendbuf, nelm, MPI_DOUBLE_COMPLEX, 1, 77, COMM_SUBD, ierr)
               ELSE
               CALL MPI_Recv(recvbuf, nelm, MPI_DOUBLE_COMPLEX, 0, 77, COMM_SUBD, MPI_STATUS_IGNORE, ierr)

   ! -------- unpack all four arrays at once --------
               k = 0
               DO j = 1, ngy
                 DO i = 1, ngx
                 k = k + 1
                 ce_y_yx(i,j) = recvbuf(         k)
                 cm_x_yx(i,j) = recvbuf(   blk + k)
                 ce_x_yx(i,j) = recvbuf(2*blk + k)
                 cm_y_yx(i,j) = recvbuf(3*blk + k)
                 END DO
               END DO
             END IF

             IF (mode_rank .EQ. 1) THEN

                 DO j =1,nx
                    DO k = 1,ngy
!88                     FORMAT(18g16.8)
                       tmp1= ce_x_xy(j,k)*cm_x_yx(j,k)-ce_x_yx(j,k)*cm_x_xy(j,k)
                       tmp2= ce_y_xy(j,k)*cm_y_yx(j,k)-ce_y_yx(j,k)*cm_y_xy(j,k)
                       tmp4= -cm_x_xy(j,k)*cm_y_yx(j,k)+cm_x_yx(j,k)*cm_y_xy(j,k)
                       zxy  = tmp1 / tmp4
                       fase = MODULO( ATAN2(AIMAG(zxy), REAL(zxy))*360.d0/twopi + 360.d0, 360.d0 )
                       rxy  = ( REAL(zxy)**2 + AIMAG(zxy)**2 ) / (omega*amu)

                       WRITE(100,88) x_c(j), y_c(k), &
                       DREAL(ce_x_yx(j,k)), DIMAG(ce_x_yx(j,k)), &
                       DREAL(ce_x_xy(j,k)), DIMAG(ce_x_xy(j,k)), &
                       DREAL(ce_y_yx(j,k)), DIMAG(ce_y_yx(j,k)), &
                       DREAL(ce_y_xy(j,k)), DIMAG(ce_y_xy(j,k)), &
                       DREAL(cm_x_yx(j,k)), DIMAG(cm_x_yx(j,k)), &
                       DREAL(cm_x_xy(j,k)), DIMAG(cm_x_xy(j,k)), &
                       DREAL(cm_y_yx(j,k)), DIMAG(cm_y_yx(j,k)), &
                       DREAL(cm_y_xy(j,k)), DIMAG(cm_y_xy(j,k))

88                     FORMAT(18g16.8)

                       WRITE(200+indx,99) x_c(j), y_c(k), rxy, fase

                       zyx  = -tmp2 / tmp4
                       fase = MODULO( ATAN2(AIMAG(zyx), REAL(zyx))*360.d0/twopi + 360.d0, 360.d0 )
                       ryx  = ( REAL(zyx)**2 + AIMAG(zyx)**2 ) / (omega*amu)

                       WRITE(300+indx,99) x_c(j), y_c(k), ryx, fase

99                     FORMAT(4g16.8)

                    ENDDO
                 ENDDO

                 CLOSE(100)
                 CLOSE(200+indx)
                 CLOSE(300+indx)

           ENDIF ! si soy el procesado que escribe
        ENDIF ! si estoy en el grupo_local de procesadores

889     FORMAT(g16.8)

        e1= e1o
        RETURN
      END SUBROUTINE output

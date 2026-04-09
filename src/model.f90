      SUBROUTINE model

        USE sizes
        USE commons_3d
        IMPLICIT NONE

        REAL(DP) :: epsilon, prof(7)
        INTEGER :: i, j, k, l, mm

150     FORMAT(e24.8)

        epsilon= 1.d-5
!*******************************************************************
!        READ(5,170) 
        nmat=ngx*ngy*ngz
!******************************************************************
      !   OPEN(15,file='conduc1',status='unknown')
      !   READ(15,150)(sig1(i), i= 1, nmat)
      !   CLOSE(15)

      !   OPEN(15,file='conduc2',status='unknown')
      !   READ(15,150)(sig2(i), i= 1, nmat)
      !   CLOSE(15)

      !   OPEN(15,file='conduc3',status='unknown')
      !   READ(15,150)(sig3(i), i= 1, nmat)
      !   CLOSE(15)

      !   OPEN(15,file='conduc_input_back',status='unknown')
      !   READ(15,150)(sigb(i), i= 1, nmat)

        OPEN(15,file='conduc_input',status='unknown')
        READ(15,150)(sig(i), i= 1, nmat)

        mm= 1
        DO l= 1, ngz
           DO k= 1, ngy
              DO j= ngx, 1, -1
                 matnum(j,k,l)= mm
                 mm = mm + 1
              ENDDO
           ENDDO
        ENDDO


        mm= 1
        DO l= 1, nz
           DO k= 1, ny
              DO j= ngx, 1, -1
                 matnumloc(j,k,l)= mm
                 mm = mm + 1
              ENDDO
           ENDDO
        ENDDO

      ! Handle boundaries in the x-direction
        DO k = 0, ngy + 1
           DO l = 0, ngz + 1
               matnum(0, k, l) = matnum(1, k, l)         ! Left boundary
               matnum(ngx + 1, k, l) = matnum(ngx, k, l) ! Right boundary
           END DO
        END DO
     
       ! Handle boundaries in the y-direction
       DO j = 0, ngx + 1
           DO l = 0, ngz + 1
               matnum(j, 0, l) = matnum(j, 1, l)         ! Bottom boundary
               matnum(j, ngy + 1, l) = matnum(j, ngy, l) ! Top boundary
           END DO
       END DO
     
       ! Handle boundaries in the z-direction
       DO j = 0, ngx + 1
           DO k = 0, ngy + 1
               matnum(j, k, 0) = matnum(j, k, 1)         ! Front boundary
               matnum(j, k, ngz + 1) = matnum(j, k, ngz) ! Back boundary
           END DO
       END DO

     

!c  Reads data for primary fields. Six-layer earth model

        READ(15,150) e1, e2, e3, e4, e5
        READ(15,150) sigma0, sigma1, sigma2, sigma3, sigma4, sigma5
           
        e1o= e1
!*****************************************************************        
        close(15)
           
        IF(blockid.EQ.ncontrol .AND. ncontrol .EQ. mode_rank) THEN
           WRITE(6,353)
353        FORMAT(1x/,' Primary Fields (six layer model) ')
           
           WRITE(6,354) e1, e2, e3, e4, e5, sigma0, sigma1, sigma2, sigma3, &
                sigma4, sigma5
354        FORMAT(    ' thickness (layer 1 - Air) = ',e20.8,' m  ',/, &
                ' thickness (layer 2) = ',e20.8,' m  ',/, &
                ' thickness (layer 3) = ',e20.8,' m  ',/, &
                ' thickness (layer 4) = ',e20.8,' m  ',/, &
                ' thickness (layer 5) = ',e20.8,' m  ',/, &
                ' conductivity (layer 1 - Air) = ',e20.8,' 1/ohm*m ',/, &
                ' conductivity (layer 2) = ',e20.8,' 1/ohm*m ',/, &
                ' conductivity (layer 3) = ',e20.8,' 1/ohm*m ',/, &
                ' conductivity (layer 4) = ',e20.8,' 1/ohm*m ',/, &
                ' conductivity (layer 5) = ',e20.8,' 1/ohm*m ',/, &
                ' conductivity (semispace) = ',e20.8,' 1/ohm*m ')
           
        ENDIF
        
        zsup= znod(ngz+1)- e1
        !      zsup=e5+e4+e3+e2    !el eje z tiene el cero en el fondo del dominio en este programa
        
        l= 1
        DO WHILE(sig((matnum(1,1,l))).GT.sigma0*10.0)
           l= l+1
        ENDDO
        izsup= l
        
        prof(1)= znod(ngz+1)
        prof(2)= znod(ngz+1)-e1
        prof(3)= prof(2)-e2
        prof(4)= prof(3)-e3
        prof(5)= prof(4)-e4
        prof(6)= prof(5)-e5
        prof(7)= 0.d0
        
        IF(blockid.EQ.ncontrol .AND. ncontrol .EQ. mode_rank) THEN
           DO j= 1, 7
              WRITE(*,*)'indice ',j,' posicion de la interfaz en z ',prof(j)
           ENDDO
        ENDIF
        
        DO j= 2, 6
           l= 0
           DO WHILE(znod(l+1).LE.prof(j))
              l= l+1
           ENDDO
           IF(j.EQ.2)izsup= l
           prof(j)= znod(l)
        ENDDO
        
        e1= prof(1)-prof(2)
        e2= prof(2)-prof(3)
        e3= prof(3)-prof(4)
        e4= prof(4)-prof(5)
        e5= prof(5)-prof(6)
        
        
        IF(blockid.EQ.ncontrol .AND. ncontrol .EQ. mode_rank) THEN         
           WRITE(*,*)'e1',e1,'e2',e2,'e3',e3,'e4',e4,'e5',e5
           
           WRITE(6,*)'posicion de la superficie en el modelo',zsup, &
                znod(ngz+1),e1
           WRITE(6,*)'posicion de la superficie segun la grilla', &
                izsup,znod(izsup)
           CALL flush(6)
        ENDIF
        
        RETURN
      END SUBROUTINE model
                        

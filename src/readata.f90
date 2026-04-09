      SUBROUTINE readata

        USE sizes
        USE commons_3d
        USE globalvars
        IMPLICIT NONE
        
        INTEGER :: i, j, k, l, nreg, nn0, nn
        
        ncontrol= 0

150     FORMAT(e30.20)
170     FORMAT(i15)
155     FORMAT(5e15.5)

        open(5,file='input_mesh_MT',status='unknown')

        read(5,*)numfreq

        do j=1,numfreq
           read(5,*) frequency(j)
        enddo

        READ(5,*)
        
        READ(5,170) nreg
        
        nn0= 0
        DO i=1,nreg
           READ(5,150) xmin
           READ(5,150) xmax
           READ(5,170) nn
           DO j= nn0+1, nn0+nn+1 
              xnod(j)= xmin + (xmax - xmin)*(float(j-nn0-1)/float(nn))
           ENDDO
           nn0 = nn+nn0
        ENDDO

        
        READ(5,170) nreg
        
        nn0= 0
        
        DO i=1,nreg
           READ(5,150) ymin
           READ(5,150) ymax
           READ(5,170) nn
           DO k= nn0+1, nn0+nn+1
              ynod(k)= ymin + (ymax - ymin)*(float(k-nn0-1)/float(nn))
           ENDDO
           nn0= nn+nn0
        ENDDO
        READ(5,170) nreg
        nn0= 0
        DO i= 1, nreg
           READ(5,150) zmin
           READ(5,150) zmax
           READ(5,170) nn
           DO l= nn0+1, nn0+nn+1
              znod(l)= zmin + (zmax - zmin)*(float(l-nn0-1)/float(nn))
           ENDDO
           nn0= nn+nn0
        ENDDO
        
        DO j= 1, ngx
           hx(j)= xnod(j+1)-xnod(j)
        ENDDO
        
        DO k= 1, ngy
           hy(k)= ynod(k+1)-ynod(k)
        ENDDO
        
        DO l= 1, ngz
           hz(l)= znod(l+1)-znod(l)
        ENDDO

        close(5)

        IF (nsy == 1 .and. nsz == 1 .and. nsub_proc == 1) THEN
         maxiter = 1
         relx= 0.d0
        ENDIF
                
!c
!c   computes local grid sizes for each processor
!c
        nx= ngx          !       
        ny= ngy/nsy      !   OJO QUE NGY Y NGZ DEBEN SER SIEMPRE DIVISIBLES
        nz= ngz/nsz      !   POR EL MAXIMO NUMERO DE PROCS. EN UNA DIRECCION                             ! /
        nmax= MAX(nz,ny)
!c
!c     nsize es el tamanio del vector para message passing en setbuf.

        nsize= 4*(nmax+2)*ngx 

        IF(ncontrol.EQ.blockid .AND. ncontrol.EQ.mode_rank) THEN

           WRITE(6,103) nx, ny, nz
103        FORMAT(/1x,'# subdomains  in x-direction nx = ',i5, &
                /,1x,'# subdomains  in y-direction ny = ',i5, &
                /,1x,'# subdomains  in z-direction nz = ',i5)

           WRITE(6,1031) ngx, ngy, ngz
1031       FORMAT(/1x,'# subdomains  in x-direction ngx = ',i5, &
                /,1x,'# subdomains  in y-direction ngy = ',i5, &
                /,1x,'# subdomains  in z-direction ngz = ',i5)

!c      write(6,104) freq
104        FORMAT(/1x,' frequency = ',g20.8,' Hz ')
           
           WRITE(6,105) maxiter
105        FORMAT(/1x,'  maxiter =  max. # inner iterations ',i5)
           
           WRITE(6,106) convergence_tol
106        FORMAT(/1x,'  reduc =  relative error required  ',e20.8)

!c      WRITE(6,108) omega
108        FORMAT(/1x,'omega = ',g20.8,' Hz ')

           WRITE(6,109)
109        FORMAT(/1x,' vector xnod  in readata '//)
           WRITE(6,155)(xnod(j), j= 1, ngx+1)
           
           WRITE(6,111)
111        FORMAT(/1x,' vector ynod  in readata '//)
           WRITE(6,155)(ynod(j), j= 1, ngy+1)
           
           WRITE(6,110)
110        FORMAT(/1x,' vector znod  in readata '//)
           WRITE(6,155)(znod(j), j= 1, ngz+1)
           
           WRITE(6,121)
121        FORMAT(/1x,' vector hx  in readata '//)
           WRITE(6,155)(hx(j), j= 1, ngx)
           
           WRITE(6,122)
122        FORMAT(/1x,' vector hy  in readata '//)
           WRITE(6,155)(hy(j), j= 1, ngy)
           
           WRITE(6,123)
123        FORMAT(/1x,' vector hz  in readata '//)
           WRITE(6,155)(hz(j), j= 1, ngz)
             
        ENDIF
!******************************************************
        RETURN
      END SUBROUTINE readata
        

SUBROUTINE msource

   USE sizes 
   USE commons_3d
   IMPLICIT NONE

   REAL(DP) :: z, tmp, sigmap
   INTEGER :: i, j, k, l
   COMPLEX(DP) ::  ep
   INTEGER :: igz, igy
   !character(40) :: filename

   !WRITE(filename, '("msource_ddg_", I0,".txt")') blockid

   !OPEN(UNIT=10, FILE=filename, STATUS="REPLACE")

   DO j= 1, mnx
      DO k= 1, mny
         DO l= 1, mnz      
            source(j,k,l)= dcmplx(0.d0,0.d0)
         ENDDO
      ENDDO
   ENDDO

   DO l= 1, mnz    
      igz= ((nsz-1-auxz)-vblockid(2,1))*nz + l  

      IF(znod(ngz+1)-znod(igz+1).LT.e1) sigmap= sigma0

      IF(znod(ngz+1)-znod(igz+1).GE.e1.AND. &
           znod(ngz+1)-znod(igz+1).LT.e2+e1) &
      sigmap= sigma1

      IF(znod(ngz+1)-znod(igz+1).GE.e2+e1.AND. &
           znod(ngz+1)-znod(igz+1).LT.e1+e2+e3) &
           sigmap= sigma2

      IF(znod(ngz+1)-znod(igz+1).GE.e3+e1+e2.AND. &
           znod(ngz+1)-znod(igz+1).LT.e1+e2+e4+e3) &
           sigmap= sigma3

      IF(znod(ngz+1)-znod(igz+1).GE.e4+e1+e2+e3.AND. &
           znod(ngz+1)-znod(igz+1).LT.e1+e2+e5+e3+e4) &
           sigmap= sigma4

      IF(znod(ngz+1)-znod(igz+1).GE.e5+e1+e3+e2+e4) sigmap= sigma5

      DO k= 1, mny
         igy= vblockid(1,1)*ny + k
         DO j= mnx,1,-1
            i= matnum(j,igy,igz)
            tmp= (sig(i)-sigmap)
            IF(znod(igz+1).GT.znod(izsup)) tmp= 0.d0
            h3= hx(j)*hy(igy)*hz(igz)
            z= znod(igz) + hz(igz)/2.0d0
            source(j,k,l)= -.25d0*h3*tmp*ep(z)
            ! Write the data in three columns, handling complex numbers
         !WRITE(blockid+10,* ) j, k, l, i,  sig(i),real(source(j,k,l)), aimag(source(j,k,l))
         ENDDO
      ENDDO
   ENDDO


END SUBROUTINE msource

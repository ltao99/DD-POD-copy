      SUBROUTINE inbeta

        USE sizes
        USE commons_3d
        USE globalvars
        IMPLICIT NONE

        REAL(DP) :: alpha
        INTEGER :: j, k, i, l
        INTEGER :: igz, igy

        DO j= 1, mnx
           DO k= 1, mny
              igy= vblockid(1,1)*ny + k
              DO l= 1, mnz
                 igz= ((nsz-1-auxz)-vblockid(2,1))*nz+ l

                 i= matnum(j,igy-1,igz)
                 alpha= dsqrt(sig(i)/(2.0*omega*amu))
                 i= matnum(j,igy,igz)
                 alpha= alpha + dsqrt(sig(i)/(2.0*omega*amu))
                 bet_w(j,k,l)= .5d0*alpha*penal*(1.d0 -aim)
                 
                 i= matnum(j,igy+1,igz)
                 alpha= dsqrt(sig(i)/(2.0*omega*amu))
                 i= matnum(j,igy,igz)
                 alpha= alpha + dsqrt(sig(i)/(2.0*omega*amu))
                 bet_e(j,k,l)= .5d0*alpha*penal*(1.d0 -aim)
                 
                 i= matnum(j,igy,igz-1)
                 alpha= dsqrt(sig(i)/(2.0*omega*amu))
                 i= matnum(j,igy,igz)
                 alpha= alpha + dsqrt(sig(i)/(2.0*omega*amu))
                 bet_s(j,k,l)= .5d0*alpha*penal*(1.d0 -aim)
                 
                 i= matnum(j,igy,igz+1)
                 alpha= dsqrt(sig(i)/(2.0*omega*amu))
                 i= matnum(j,igy,igz)
                 alpha= alpha + dsqrt(sig(i)/(2.0*omega*amu))
                 bet_n(j,k,l)= .5d0*alpha*penal*(1.d0 -aim)
                 
              ENDDO
           ENDDO
        ENDDO
!1034   CONTINUE

        RETURN
      END SUBROUTINE inbeta

     
!**********************************************************************

     SUBROUTINE setdelta

       USE sizes
       USE commons_3d
       IMPLICIT NONE
       
       INTEGER :: j,k,l,y, z, i, subz, suby,dy, dz

       ! Deltas for Alphas

       DO j= 1, mnx
          deltaf(j)= 0.d0
          deltab(j)= 0.d0
       ENDDO
!232    CONTINUE

       DO j= 1, mny
          deltae(j)= 0.d0
          deltaw(j)= 0.d0
       ENDDO
!234    CONTINUE

       DO j= 1, mnz
          deltas(j)= 0.d0
          deltan(j)= 0.d0
       ENDDO
!233    CONTINUE

       IF(blockty.EQ.0)THEN
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 1.d0
          deltas(1)= 0.d0
          deltan(mnz)= 1.d0
       ELSEIF(blockty.EQ.1)THEN
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 0.d0
          deltas(1)= 0.d0
          deltan(mnz)= 1.d0
       ELSEIF(blockty.EQ.2)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 1.d0
          deltaw(1)= 0.d0
          deltas(1)= 0.d0
          deltan(mnz)= 1.d0 
       ELSEIF(blockty.EQ.3)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 1.d0
          deltas(1)= 0.d0
          deltan(mnz)= 0.d0   
       ELSEIF(blockty.EQ.4)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 0.d0
          deltas(1)= 0.d0
          deltan(mnz)= 0.d0   
       ELSEIF(blockty.EQ.5)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 1.d0
          deltaw(1)= 0.d0
          deltas(1)= 0.d0
          deltan(mnz)= 0.d0 
       ELSEIF(blockty.EQ.6)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 1.d0
          deltas(1)= 1.d0
          deltan(mnz)= 0.d0 
       ELSEIF(blockty.EQ.7)THEN 
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 0.d0
          deltaw(1)= 0.d0
          deltas(1)= 1.d0
          deltan(mnz)= 0.d0   
       ELSEIF(blockty.EQ.8) THEN
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 1.d0
          deltaw(1)= 0.d0
          deltas(1)= 1.d0
          deltan(mnz)= 0.d0   
      ELSEIF(blockty.EQ.10)THEN
         deltab(1)= 1.d0
         deltaf(mnx)= 1.d0
         deltae(mny)= 1.d0
         deltaw(1)= 1.d0
         deltas(1)= 0.d0
         deltan(mnz)= 1.d0
      ELSEIF(blockty.EQ.11)THEN
         deltab(1)= 1.d0
         deltaf(mnx)= 1.d0
         deltae(mny)= 1.d0
         deltaw(1)= 1.d0
         deltas(1)= 0.d0
         deltan(mnz)= 0.d0
      ELSEIF(blockty.EQ.12)THEN
         deltab(1)= 1.d0
         deltaf(mnx)= 1.d0
         deltae(mny)= 1.d0
         deltaw(1)= 1.d0
         deltas(1)= 1.d0
         deltan(mnz)= 0.d0       
      ENDIF

       IF(blockty.EQ.9)THEN
          deltab(1)= 1.d0
          deltaf(mnx)= 1.d0
          deltae(mny)= 1.d0
          deltaw(1)= 1.d0
          deltas(1)= 1.d0
          deltan(mnz)= 1.d0   
       ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Deltas for Betas
       DO j= 1, mnx
         deltabf(j)= 1.d0
         deltabb(j)= 1.d0
      ENDDO
!232    CONTINUE
      DO j= 1, mny
         deltabe(j)= 1.d0
         deltabw(j)= 1.d0
      ENDDO
!234    CONTINUE
      DO j= 1, mnz
         deltabs(j)= 1.d0
         deltabn(j)= 1.d0
      ENDDO

      IF(blockty.EQ.0)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 1.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 1.d0
      ELSEIF(blockty.EQ.1)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 0.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 1.d0
      ELSEIF(blockty.EQ.2)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 0.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 1.d0
      ELSEIF(blockty.EQ.3)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 1.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 0.d0 
      ELSEIF(blockty.EQ.4)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 0.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 0.d0    
      ELSEIF(blockty.EQ.5)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 0.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 0.d0 
      ELSEIF(blockty.EQ.6)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 1.d0
         deltabs(1)= 1.d0
         deltabn(mnz)= 0.d0 
      ELSEIF(blockty.EQ.7)THEN 
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 0.d0
         deltabw(1)= 0.d0
         deltabs(1)= 1.d0
         deltabn(mnz)= 0.d0    
      ELSEIF(blockty.EQ.8)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 0.d0
         deltabs(1)= 1.d0
         deltabn(mnz)= 0.d0
      ELSEIF(blockty.EQ.10)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 1.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 1.d0
      ELSEIF(blockty.EQ.11)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 1.d0
         deltabs(1)= 0.d0
         deltabn(mnz)= 0.d0
      ELSEIF(blockty.EQ.12)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 1.d0
         deltabs(1)= 1.d0
         deltabn(mnz)= 0.d0      
      ENDIF

      IF(blockty.EQ.9)THEN
         deltabb(1)= 1.d0
         deltabf(mnx)= 1.d0
         deltabe(mny)= 1.d0
         deltabw(1)= 1.d0
         deltabs(1)= 1.d0
         deltabn(mnz)= 1.d0 
      ENDIF


      IF(nsub_proc .GT.1) THEN
         DO suby = 1,nsy-1
            k = suby * ny
            deltabe(k) = 0.d0
            deltabw(k + 1)= 0.d0
         ENDDO
         DO subz = 1,nsz-1
            l = subz * nz
            deltabn(l) = 0.d0
            deltabs(l + 1) = 0.d0
         ENDDO
      ENDIF
      
       RETURN
 END SUBROUTINE setdelta
     

!**********************************************************************


   SUBROUTINE setalpha

       USE sizes
       USE commons_3d
       USE globalvars
       IMPLICIT NONE
       
       INTEGER :: i, j, k, l, igz, igy

       DO k= 1, mny
         igy= vblockid(1,1)*ny + k
          DO l= 1, mnz

            igz= ((nsz-1-auxz)-vblockid(2,1))*nz+ l
             
             i= matnum(1,igy,igz)
             
             alphab(k,l)= dsqrt(sig(i)/(2.0*omega*amu))
             
             i= matnum(mnx,igy,igz)
             
             alphaf(k,l)= dsqrt(sig(i)/(2.0*omega*amu))
          ENDDO
       ENDDO
!232    CONTINUE

       DO j= 1, mnx 
          DO k= 1, mny

             igy= vblockid(1,1)*ny + k
             i= matnum(j,igy,1)
             alphas(j,k)= dsqrt(sig(i)/(2.0*omega*amu))
             
             i= matnum(j,igy,ngz)
             alphan(j,k)= dsqrt(sig(i)/(2.0*omega*amu))
          ENDDO
       ENDDO
!234    CONTINUE

       DO j= 1, mnx
          DO l= 1, mnz

             igz= ((nsz-1-auxz)-vblockid(2,1))*nz+ l
             i= matnum(j,1,igz)
             alphaw(j,l)= dsqrt(sig(i)/(2.0*omega*amu))
             
             i= matnum(j,ngy,igz)
             alphae(j,l)= dsqrt(sig(i)/(2.0*omega*amu))
          ENDDO
       ENDDO
!236    CONTINUE

      !  DO j = 1,mnx
      !    DO l = 1,mnz
      !       WRITE(blockid,*) alphaw(j,l),alphae(j,l)
      !    ENDDO
      ! ENDDO

      !  DO j = 1,mnx
      !    DO l = 1,mny
      !       WRITE(blockid*10,*) alphas(j,l),alphan(j,l)
      !    ENDDO
      ! ENDDO

      ! DO j = 1,mny
      !    DO l = 1,mnz
      !       WRITE(blockid*100,*) alphas(j,l),alphan(j,l)
      !    ENDDO
      ! ENDDO


       RETURN
   END SUBROUTINE setalpha

         SUBROUTINE subtype
         USE sizes
         USE commons_3d


         IMPLICIT NONE

       IF (nsy == 1 .AND. nsz > 1) THEN
              IF (vblockid(2,1) == 0) THEN
                     blockty = 10
              ELSEIF (vblockid(2,1) == nsz - 1) THEN
                     blockty = 12
              ELSE
                     blockty = 11
              ENDIF
       ENDIF

!c
!c       In order to run serial code uncomment
!c       the following line
!c
         IF(nprocs.EQ.2)blockty= 9  
         
         RETURN
       END SUBROUTINE subtype
      
!************************************************************************

       SUBROUTINE setbuff
         USE sizes
         USE commons_3d

         IMPLICIT NONE        

         IF(blockty.EQ.10)THEN
              CALL bufdat10

         ELSEIF(blockty.EQ.11)THEN
              CALL bufdat11

         ELSEIF(blockty.EQ.12)THEN
              CALL bufdat12
         ENDIF
         
         RETURN
       END SUBROUTINE setbuff
    
!************************************************************************** 

       SUBROUTINE tbrytovec
         USE sizes
         USE commons_3d
         
         IMPLICIT NONE
         
         INTEGER :: j, k, ind, ind1, mm, ind2
         INTEGER :: n0, n1, n2, n3, n4
         
         n0= ny*nx
         
         DO  k= 1, ny
            n1= (k-1)*nx
            n2= n0+n1
            n3= 2*n0+n1
            n4= 3*n0+n1
            DO j= 1, nx
               sendvec(n1+j)=  eo_nx(j,k,nz)
               sendvec(n2+j)=  eo_ny(j,k,nz)
               sendvec(n3+j)=  lamo_nx(j,k,nz)
               sendvec(n4+j)=  lamo_ny(j,k,nz)
            ENDDO
         ENDDO
         
         RETURN 
       END SUBROUTINE tbrytovec
            
!************************************************************************** 


       SUBROUTINE bbrytovec
         USE sizes
         USE commons_3d

         IMPLICIT NONE
         
         INTEGER :: j, k, ind, ind1, mm, ind2
         INTEGER :: n0, n1, n2, n3, n4
         
         n0= ny*nx
         
         DO  k= 1, ny
            n1= (k-1)*nx
            n2= n0+n1
            n3= 2*n0+n1
            n4= 3*n0+n1
            DO j= 1, nx
               sendvec(n1+j)=  eo_sx(j,k,1)
               sendvec(n2+j)=  eo_sy(j,k,1)
               sendvec(n3+j)=  lamo_sx(j,k,1)
               sendvec(n4+j)=  lamo_sy(j,k,1)
            ENDDO
         ENDDO
         
         RETURN 
       END SUBROUTINE bbrytovec
       
     
!***************************************************************************
       SUBROUTINE vectotop
         USE sizes
         USE commons_3d

         IMPLICIT NONE
         
         INTEGER :: j, k, ind, ind1, mm, ind2
         INTEGER :: n0, n1, n2, n3, n4
         
         n0= ny*nx
         
         DO  k= 1, ny
            n1= (k-1)*nx
            n2= n0+n1
            n3= 2*n0+n1
            n4= 3*n0+n1
            DO  j= 1, nx
               eo_sx(j,k,nz+1)= recvec(n1+j)
               eo_sy(j,k,nz+1)= recvec(n2+j)
               lamo_sx(j,k,nz+1)= recvec(n3+j)
               lamo_sy(j,k,nz+1)= recvec(n4+j)                 
            ENDDO
         ENDDO
         
         RETURN 
       END SUBROUTINE vectotop

    
!************************************************************************** 

       
       SUBROUTINE vectobot
         USE sizes
         USE commons_3d

         IMPLICIT NONE
         
         INTEGER :: ind, ind1, j, k, mm, ind2
         INTEGER :: n0, n1, n2, n3, n4
         
         n0= ny*nx
         
         DO  k= 1, ny
            n1= (k-1)*nx
            n2= n0+n1
            n3= 2*n0+n1
            n4= 3*n0+n1
            DO  j= 1, nx
               eo_nx(j,k,0)= recvec(n1+j)
               eo_ny(j,k,0)= recvec(n2+j)
               lamo_nx(j,k,0)= recvec(n3+j)
               lamo_ny(j,k,0)= recvec(n4+j) 
            ENDDO
         ENDDO
         
         RETURN 
       END SUBROUTINE vectobot
           
!***************************************************************************


     SUBROUTINE bufdat10
       USE sizes
       USE commons_3d

       IMPLICIT NONE
       
!c      down shift to set up buffer

         CALL bbrytovec

!c       if (blockid.eq.ncontrol) then
!c         write(6,3520)
!c3520     format(1x/,' llame a bbrytovec ')
!c       endif  

!c       exchange info - send to bottom
         CALL cmsend(2*(blockid+nsy)+mode_rank,1,sendvec,nsize)

!c       if (blockid.eq.ncontrol) then
!c         write(6,3527)
!c3527     format(1x/,' llame a cmsend ')
!c       endif  
 
!c       exchange info - receive from botton 
         CALL cmrec(2*(blockid+nsy)+mode_rank,1,recvec,nsize)

!c       if (blockid.eq.ncontrol) then
!c         write(6,3528)
!c3528     format(1x/,' llame a cmrec ')
!c       endif  
	
!c       store information received 

         CALL vectobot

!c       if (blockid.eq.ncontrol) then
!c         write(6,3521)
!c3521     format(1x/,' bufdat0 ')
!c       endif  


       RETURN
     END SUBROUTINE bufdat10

!************************************************************************


     SUBROUTINE bufdat11
       USE sizes
       USE commons_3d

       IMPLICIT NONE
       
!c    /* down shift to set up buffer */

      CALL bbrytovec

!c	/* exchange info - receive from up send to botton */

      CALL cmsenrec(2*(blockid-nsy)+mode_rank,1,recvec,nsize, &
           2*(blockid+nsy)+mode_rank,1,sendvec,nsize)

!c	/* store information received */

      CALL vectotop

!c	/* up shift to set down buffer */

      CALL tbrytovec

!c	/* exchange info - receive from down send to up */

      CALL cmsenrec(2*(blockid+nsy)+mode_rank,1,recvec,nsize, &
           2*(blockid-nsy)+mode_rank,1,sendvec,nsize)

!c	/* store information received */

      CALL vectobot

!c         write(6,3528)
!c3528     format(1x/,' sali de bufdat4 ')


       RETURN
     END SUBROUTINE bufdat11

     
!************************************************************************


     SUBROUTINE bufdat12
       USE sizes
       USE commons_3d

       IMPLICIT NONE
       
!c	/* down shift to set up buffer */
!c	/* exchange info - receive from up */

       CALL cmrec(2*(blockid-nsy)+mode_rank,1,recvec,nsize)

!c	/* store information received */

       CALL vectotop

!c	/* up shift to set down buffer */

       CALL tbrytovec

!c	/* exchange info with node on top*/

       CALL cmsend(2*(blockid-nsy)+mode_rank,1,sendvec,nsize)

!c         write(6,3564)
!c3564     format(1x/,' sali de bufdat8 ')


       RETURN
     END SUBROUTINE bufdat12

!************************************************************************

     SUBROUTINE cmsend(b_id,dumb,sendd,size)
       USE sizes
       USE commons_3d
       
       INCLUDE 'mpif.h'
        
       
       INTEGER :: b_id, dumb, size, itag
       
!c        complex *16 sendd(0:4*(mnmax+2)**2)
       
       COMPLEX(DP) :: sendd(1:4*(mnmax+2)*ngx)

!csp2-----------------------------------------------------------------------
       itag= 1
       CALL MPI_SEND(sendd,size,MPI_DOUBLE_COMPLEX,b_id, &
            itag,MPI_COMM_WORLD,mpierr)

       IF(mpierr .NE. 0)THEN
          WRITE(6,*)' MPI error !!! (11)'
          STOP
       ENDIF
     
!csp2  -----------------------------------------------------------
       RETURN
     END SUBROUTINE cmsend

     
!***********************************************************************


     SUBROUTINE cmrec(b_id,dumb,recc,size)
       USE sizes
       USE commons_3d
      
       INCLUDE 'mpif.h'
      
       
       INTEGER :: status(MPI_STATUS_SIZE)

       INTEGER :: b_id, dumb, size, itag
       
!c        complex *16 recc(0:4*(mnmax+2)**2)

       COMPLEX(DP) :: recc(1:4*(mnmax+2)*ngx)
       
!csp2---------------------------------------------------------

       itag= 1
       CALL MPI_RECV(recc,size,MPI_DOUBLE_COMPLEX,b_id, &
            itag,MPI_COMM_WORLD,status,mpierr)
     
       IF(mpierr .NE. 0)THEN
          WRITE(6,*)' MPI error !!! (12)'
          STOP
       ENDIF
       
!csp2  ------------------------------------------------------
       RETURN
     END SUBROUTINE cmrec


!************************************************************************


     SUBROUTINE cmsenrec(b_id1,dumb1,recc,size1, &
          b_id2,dumb2,sendd,size2)
       USE sizes
       USE commons_3d
       
       INCLUDE 'mpif.h'
       
       INTEGER :: status(MPI_STATUS_SIZE)
       INTEGER :: b_id1, dumb1, b_id2, dumb2, size1, size2, itag

!c        complex *16 sendd(0:4*(mnmax+2)**2),recc(0:4*(mnmax+2)**2)

       COMPLEX(DP) :: sendd(1:4*(mnmax+2)*ngx), recc(1:4*(mnmax+2)*ngx)
        
!csp2-----------------------------------------------------------------------

        itag= 1
        CALL MPI_SEND(sendd,size2,MPI_DOUBLE_COMPLEX,b_id2, &
             itag,MPI_COMM_WORLD,mpierr)
        IF(mpierr.NE. 0)THEN
           WRITE(6,*)' MPI error !!! (13)'
           STOP
        ENDIF

        itag= 1
        
        CALL MPI_RECV(recc,size1,MPI_DOUBLE_COMPLEX,b_id1, &
             itag,MPI_COMM_WORLD,status,mpierr)
     
        IF(mpierr.NE. 0)THEN
           WRITE(6,*)' MPI error !!! (14)'
           STOP
        ENDIF
        
!csp2  ---------------------------------------------------------------
       RETURN
     END SUBROUTINE cmsenrec

!**************************************************************************     

 

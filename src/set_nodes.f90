subroutine set_nodes_loc
   !      modified for global calculations
   USE sizes
   USE commons_3d
   implicit none 
   integer :: j,jk,ie_g_left,ie_g_down,aux_g
   integer :: nx_g,ny_g,nz_g,iz_g,iy_g,ie_g,ie_t
   integer :: rn1,rn2,rn3
   nx_g=ngx
   ny_g=ngy/nsy
   nz_g=ngz/nsz
   ie_t=ngx*ny_g*nz_g
   
   !     Set local ( entire domain) nodes
   rn1=nxs
   iz_g=1
   iy_g=1
   ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g

   forall (j=1:12) lc_nod(ie_g,j)=j
   ie_g=ie_g-1
      
   lc_nod(ie_g,1:4)=lc_nod(ie_g+1,1:4)+12
   lc_nod(ie_g,5)=6
   lc_nod(ie_g,6:9)=lc_nod(ie_g+1,6:9)+11
   lc_nod(ie_g,10)=21
   lc_nod(ie_g,11)=12
   lc_nod(ie_g,12)=22

   do j=3,nx_g
      jk= (j-3)*10
      ie_g= ie_g-1
      lc_nod(ie_g,1)=23+jk
      lc_nod(ie_g,2)=24+jk
      lc_nod(ie_g,3)=25+jk
      lc_nod(ie_g,4)=26+jk
      lc_nod(ie_g,5)=17+jk
      lc_nod(ie_g,6)=27+jk
      lc_nod(ie_g,7)=28+jk
      lc_nod(ie_g,8)=29+jk
      lc_nod(ie_g,9)=30+jk
      lc_nod(ie_g,10)=31+jk
      lc_nod(ie_g,11)=22+jk
      lc_nod(ie_g,12)=32+jk
   enddo

   aux_g=lc_nod(ie_g,12)


   IF ((nsy * nsz) .EQ. (ngy * ngz)) THEN
      aux_g=lc_nod(ie_g,12)
      GO TO 50
   ENDIF


   ! do iz_g=2,nz_g
   !    iy_g=1
   !    ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
   !    ie_g_down= ie_g-nx_g*ny_g

   !    lc_nod(ie_g,1)=aux_g+1
   !    lc_nod(ie_g,2)=aux_g+2
   !    lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
   !    lc_nod(ie_g,4)=aux_g+3
   !    lc_nod(ie_g,5)=aux_g+4
   !    lc_nod(ie_g,6)=aux_g+5
   !    lc_nod(ie_g,7)=aux_g+6
   !    lc_nod(ie_g,8)=aux_g+7
   !    lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
   !    lc_nod(ie_g,10)=aux_g+8
   !    lc_nod(ie_g,11)=aux_g+9
   !    lc_nod(ie_g,12)=aux_g+10

   !    ie_g=ie_g-1
   !    ie_g_down=ie_g_down-1

   !    lc_nod(ie_g,1)=lc_nod(ie_g+1,1)+10
   !    lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+10
   !    lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
   !    lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 10
   !    lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
   !    lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 9
   !    lc_nod(ie_g,7)=lc_nod(ie_g +1,7)+ 9
   !    lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 9
   !    lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
   !    lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 9
   !    lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
   !    lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8

   !    do j=2,nx_g-1
   !       ie_g=ie_g-1
   !       ie_g_down=ie_g_down-1

   !       lc_nod(ie_g,1)=lc_nod(ie_g+1,1)+8
   !       lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+8
   !       lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
   !       lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 8
   !       lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
   !       lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 8
   !       lc_nod(ie_g,7)=lc_nod(ie_g +1,7)+ 8
   !       lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 8
   !       lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
   !       lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 8
   !       lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
   !       lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8
   !    enddo
   !    aux_g=lc_nod(ie_g,12)
   ! enddo

   iy_g=2
   ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
   ie_g_left= nx_g*ny_g*(iz_g-1)+(iy_g-1)*nx_g
   ! element at the same position in x as ie_g but at the left
   ! First element of stripe 2 at the bottom, need to conected with 1 first
   ! element of stripe 1

   lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
   lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn1-1   ! it adds nxs-1 instead of nxs because this elements shares the left node
   lc_nod(ie_g,3)=lc_nod(ie_g_left,3)+ rn1-1
   lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn1-1
   lc_nod(ie_g,5)=lc_nod(ie_g_left,5)+ rn1-1
   lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn1-1
   lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
   lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn1-2
   lc_nod(ie_g,9)=lc_nod(ie_g_left,9)+ rn1-2
   lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn1-2
   lc_nod(ie_g,11)=lc_nod(ie_g_left,11)+ rn1-2
   lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn1-2

   ! Second element of iy_g stripe at iz_g=1
   ie_g=ie_g-1
   ie_g_left=ie_g_left-1      ! element at the same position in x as ie_g but at the left
   lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
   lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn1-3
   lc_nod(ie_g,3)=lc_nod(ie_g_left,3)+ rn1-3
   lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn1-3
   lc_nod(ie_g,5)=lc_nod(ie_g +1,6)     ! ie_g_left shares this too, same rn1-3
   lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn1-3
   lc_nod(ie_g,7)=lc_nod(ie_g_left,8)   ! ie_g_left does't share this one node less
   lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn1-4
   lc_nod(ie_g,9)=lc_nod(ie_g_left,9)+ rn1-4
   lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn1-4
   lc_nod(ie_g,11)=lc_nod(ie_g+1,12)      ! ie_g_left shares this too, same rn1-3
   lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn1-4

   do j=2,nx_g-1
      ie_g=ie_g-1
      ie_g_left=ie_g_left-1
      lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
      lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+ 8
      lc_nod(ie_g,3)=lc_nod(ie_g +1,3)+ 8
      lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 8
      lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
      lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 8
      lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
      lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 8
      lc_nod(ie_g,9)=lc_nod(ie_g +1,9)+ 8
      lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 8
      lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
      lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8
   enddo
   aux_g=lc_nod(ie_g,12)
  
   ! IF ((nsy * nsz) .EQ. ((ngy * ngz)/2)) THEN
   !    aux_g=lc_nod(ie_g,12)
   !    GO TO 50
   ! ENDIF

   do iy_g=3,ny_g
      ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
      ie_g_left=nx_g*ny_g*(iz_g-1)+(iy_g-1)*nx_g
      rn2=rn1-(2*nx_g)

      ! element at the same position in x as ie_g but at the left
      ! First element of stripe 2 at the bottom, need to conected with 1 first
      ! element of stripe 1
      lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
      lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn2
      lc_nod(ie_g,3)=lc_nod(ie_g_left,3)+ rn2
      lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn2
      lc_nod(ie_g,5)=lc_nod(ie_g_left,5)+ rn2
      lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn2
      lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
      lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn2
      lc_nod(ie_g,9)=lc_nod(ie_g_left,9)+ rn2
      lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn2
      lc_nod(ie_g,11)=lc_nod(ie_g_left,11)+ rn2
      lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn2

      ie_g=ie_g-1
      ie_g_left=ie_g_left-1        ! element at the same position in x as ie_g but at the left

      lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
      lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn2
      lc_nod(ie_g,3)=lc_nod(ie_g_left,3)+ rn2
      lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn2
      lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
      lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn2
      lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
      lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn2
      lc_nod(ie_g,9)=lc_nod(ie_g_left,9)+ rn2
      lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn2
      lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
      lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn2

      do j=2,nx_g-1
         ie_g=ie_g-1
         ie_g_left=ie_g_left-1
         lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
         lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+ 8
         lc_nod(ie_g,3)=lc_nod(ie_g +1,3)+ 8
         lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 8
         lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
         lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 8
         lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
         lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 8
         lc_nod(ie_g,9)=lc_nod(ie_g +1,9)+ 8
         lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 8
         lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
         lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8
      enddo
      aux_g=lc_nod(ie_g,12)
   enddo

   do iz_g=2,nz_g
      iy_g=1
      ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
      ie_g_down= ie_g-nx_g*ny_g

      lc_nod(ie_g,1)=aux_g+1
      lc_nod(ie_g,2)=aux_g+2
      lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
      lc_nod(ie_g,4)=aux_g+3
      lc_nod(ie_g,5)=aux_g+4
      lc_nod(ie_g,6)=aux_g+5
      lc_nod(ie_g,7)=aux_g+6
      lc_nod(ie_g,8)=aux_g+7
      lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
      lc_nod(ie_g,10)=aux_g+8
      lc_nod(ie_g,11)=aux_g+9
      lc_nod(ie_g,12)=aux_g+10

      ie_g=ie_g-1
      ie_g_down=ie_g_down-1

      lc_nod(ie_g,1)=lc_nod(ie_g+1,1)+10
      lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+10
      lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
      lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 10
      lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
      lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 9
      lc_nod(ie_g,7)=lc_nod(ie_g +1,7)+ 9
      lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 9
      lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
      lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 9
      lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
      lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8

      do j=2,nx_g-1
         ie_g=ie_g-1
         ie_g_down=ie_g_down-1

         lc_nod(ie_g,1)=lc_nod(ie_g+1,1)+8
         lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+8
         lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
         lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+ 8
         lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
         lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 8
         lc_nod(ie_g,7)=lc_nod(ie_g +1,7)+ 8
         lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 8
         lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
         lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 8
         lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
         lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 8
      enddo

      iy_g=2
      ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
      ie_g_left= nx_g*ny_g*(iz_g-1)+(iy_g-1)*nx_g
      ie_g_down= ie_g-nx_g*ny_g

      lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
      lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn2-1
      lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
      lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn2-1
      lc_nod(ie_g,5)=lc_nod(ie_g_left,5)+ rn2-1
      lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn2-1
      lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
      lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn2-2
      lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
      lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn2-2
      lc_nod(ie_g,11)=lc_nod(ie_g_left,11)+ rn2-2
      lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn2-2


      ie_g=ie_g-1
      ie_g_down=ie_g_down-1
      ie_g_left= ie_g_left-1

      lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
      lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+8
      lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
      lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+8
      lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
      lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+ 7
      lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
      lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+ 7
      lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
      lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+ 7
      lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
      lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+ 6

      do j=2,nx_g-1
         ie_g=ie_g-1
         ie_g_down=ie_g_down-1
         ie_g_left= ie_g_left-1

         lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
         lc_nod(ie_g,2)=lc_nod(ie_g +1,2)+6
         lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
         lc_nod(ie_g,4)=lc_nod(ie_g +1,4)+6
         lc_nod(ie_g,5)=lc_nod(ie_g +1,6)
         lc_nod(ie_g,6)=lc_nod(ie_g +1,6)+6
         lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
         lc_nod(ie_g,8)=lc_nod(ie_g +1,8)+6
         lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
         lc_nod(ie_g,10)=lc_nod(ie_g +1,10)+6
         lc_nod(ie_g,11)=lc_nod(ie_g+1,12)
         lc_nod(ie_g,12)=lc_nod(ie_g +1,12)+6
      enddo

   !    print *, aux_g
   !    DO i = 1, 40
   !       WRITE(4, *) i, lc_nod(i, 1), lc_nod(i, 2), lc_nod(i, 3), lc_nod(i, 4), lc_nod(i, 5), lc_nod(i, 6), lc_nod(i, 7), lc_nod(i, 8), lc_nod(i, 9), lc_nod(i, 10), lc_nod(i, 11), lc_nod(i, 12)
   !   END DO
     
      do iy_g=3,ny_g
         ie_g= nx_g*ny_g*(iz_g-1)+iy_g*nx_g
         ie_g_left= nx_g*ny_g*(iz_g-1)+(iy_g-1)*nx_g
         ie_g_down= ie_g-nx_g*ny_g
         rn3= rn1-(4*nx_g)

         lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
         lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn3
         lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
         lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn3
         lc_nod(ie_g,5)=lc_nod(ie_g_left,5)+ rn3
         lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn3
         lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
         lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn3
         lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
         lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn3
         lc_nod(ie_g,11)=lc_nod(ie_g_left,11)+ rn3
         lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn3

         do j=1,nx_g-1
            ie_g=ie_g-1
            ie_g_down=ie_g_down-1
            ie_g_left= ie_g_left-1
            lc_nod(ie_g,1)=lc_nod(ie_g_left,2)
            lc_nod(ie_g,2)=lc_nod(ie_g_left,2)+ rn3
            lc_nod(ie_g,3)=lc_nod(ie_g_down,4)
            lc_nod(ie_g,4)=lc_nod(ie_g_left,4)+ rn3
            lc_nod(ie_g,5)=lc_nod(ie_g_left,5)+ rn3
            lc_nod(ie_g,6)=lc_nod(ie_g_left,6)+ rn3
            lc_nod(ie_g,7)=lc_nod(ie_g_left,8)
            lc_nod(ie_g,8)=lc_nod(ie_g_left,8)+ rn3
            lc_nod(ie_g,9)=lc_nod(ie_g_down,10)
            lc_nod(ie_g,10)=lc_nod(ie_g_left,10)+ rn3
            lc_nod(ie_g,11)=lc_nod(ie_g_left,11)+ rn3
            lc_nod(ie_g,12)=lc_nod(ie_g_left,12)+ rn3
         enddo
      enddo
      aux_g=lc_nod(ie_g,12)
   enddo
50 CONTINUE
      total_l_nod=aux_g

   ! Check if ip, ij, apj are allocated; if not, allocate them
   if (.not. allocated(sol_MUMPS)) then
      allocate(sol_MUMPS(1:total_l_nod))
   endif
   
end subroutine set_nodes_loc


Subroutine set_nodes_interface


       USE sizes
       USE commons_3d
       USE globalvars

       IMPLICIT NONE
       
       INTEGER :: j,k,l,y, z, i, subz, suby,dy, dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(blockty.EQ.0)THEN
         delta_y = 1
         delta_z = -1
      ELSEIF(blockty.EQ.1)THEN
         delta_y = 2
         delta_z = -1
      ELSEIF(blockty.EQ.2)THEN 
         delta_y = -1
         delta_z = -1
      ELSEIF(blockty.EQ.3)THEN 
         delta_y = 1
         delta_z = 2  
      ELSEIF(blockty.EQ.4)THEN   
         delta_y = 2
         delta_z = 2  
      ELSEIF(blockty.EQ.5)THEN 
         delta_y = -1
         delta_z = 2  
      ELSEIF(blockty.EQ.6)THEN 
         delta_y = 1
         delta_z = 1 
      ELSEIF(blockty.EQ.7)THEN 
         delta_y = 2
         delta_z = 1    
      ELSEIF(blockty.EQ.8)THEN
         delta_y = -1
         delta_z = 1 
      ELSEIF(blockty.EQ.10)THEN
         delta_y = 0
         delta_z = -1 
      ELSEIF(blockty.EQ.11)THEN
      IF (nz .EQ. 1) THEN
         delta_y = 0
         delta_z = 1
      ELSE
         delta_y = 0
         delta_z = 2 
      ENDIF
      ELSEIF(blockty.EQ.12)THEN
         delta_y = 0
         delta_z = 1        
      ENDIF

      IF(blockty.EQ.9)THEN
         delta_y = 0
         delta_z = 0   
      ENDIF

IF (nz == 1 .AND. blockty == 11) THEN
	ndofs_interface = 4 * (abs(delta_y) * nx * nz + abs(delta_z) * nx * ny)
ELSE
	ndofs_interface = 2 * (abs(delta_y) * nx * nz + abs(delta_z) * nx * ny)
ENDIF

ALLOCATE(dofs_interface(ndofs_interface))
k = 1

! Z-direction interface assignment
IF (delta_z == 0) THEN
   ALLOCATE(interf_z(0))
ELSEIF (delta_z == 1) THEN
   ALLOCATE(interf_z(1))
   interf_z(1) = nz

   DO z = 1, SIZE(interf_z)
      dz = interf_z(z)
      DO dy = 1, ny
         DO j = 1, nx
            i = matnumloc(j, dy, dz)
            IF (blockty == 11 .AND. nz == 1 .AND. use_pod) THEN
	       dofs_interface(k) = lc_nod(i, 3)
               dofs_interface(k + 1) = lc_nod(i, 4)
               dofs_interface(k + 2) = lc_nod(i, 9)
               dofs_interface(k + 3) = lc_nod(i, 10)
               k = k + 4
            ELSE
              dofs_interface(k) = lc_nod(i, 4)
              dofs_interface(k + 1) = lc_nod(i, 10)
              k = k + 2
          ENDIF
         END DO
      END DO
   END DO

ELSEIF (delta_z == -1) THEN
   ALLOCATE(interf_z(1))
   interf_z(1) = 1

   DO z = 1, SIZE(interf_z)
      dz = interf_z(z)
      DO dy = 1, ny
         DO j = 1, nx
            i = matnumloc(j, dy, dz)
            dofs_interface(k) = lc_nod(i, 3)
            dofs_interface(k + 1) = lc_nod(i, 9)
            k = k + 2
         END DO
      END DO
   END DO

ELSE
   ALLOCATE(interf_z(2))
   interf_z = [1, nz]

   DO z = 1, SIZE(interf_z)
      dz = interf_z(z)
      DO dy = 1, ny
         DO j = 1, nx
            i = matnumloc(j, dy, dz)
            IF (dz == 1) THEN
               dofs_interface(k) = lc_nod(i, 3)
               dofs_interface(k + 1) = lc_nod(i, 9)
            ELSE
               dofs_interface(k) = lc_nod(i, 4)
               dofs_interface(k + 1) = lc_nod(i, 10)
            END IF
            k = k + 2
         END DO
      END DO
   END DO
END IF


END SUBROUTINE set_nodes_interface




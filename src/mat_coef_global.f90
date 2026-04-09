subroutine mat_coef_global
   use quicksort, only: quick_sort
   USE commons_3d
   USE globalvars
   USE sizes
   integer :: i,j,k,l
   integer :: inod,jnod,ie_gl, nnaux
   integer  :: gli(total_l_nod*12), glj(total_l_nod*12)
   integer :: gl_aux(total_l_nod*12,3),sort_chico(12), k_chica(12) !new
   integer :: counti,count2,count3
   integer :: ksorti(12),sorti(12)
   complex(DP) :: mat_local(1:12,1:12)
   complex(DP) :: global_mat(total_l_nod*12)
   complex(DP) :: global_mat_aux(total_l_nod*12)
   complex(DP),allocatable  :: global_mat_aux2(:)
   integer,allocatable ::  ksort(:),sort_array(:),gl_aux2(:,:)

   nnaux=0
   do l=1,nz
      do k=1,ny
         do j = 1,nx
            call mat_coefl(j,k,l,mat_local,ie_gl)           
            do inod=1,12
               do jnod=1,12
                  if (abs(REAL(mat_local(inod,jnod)))+abs(AIMAG(mat_local(inod,jnod))) > ACCURMACHINE) then
                     nnaux=nnaux+1
                     gl_aux(nnaux,1)=lc_nod(ie_gl,inod)
                     gl_aux(nnaux,2)=lc_nod(ie_gl,jnod)
                     gl_aux(nnaux,3)= nnaux
                     global_mat_aux(nnaux)=  mat_local(inod,jnod)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo

   allocate( sort_array(nnaux),ksort(nnaux),gl_aux2(nnaux,3))
   sort_array=gl_aux(1:nnaux,1)
   ksort=gl_aux(1:nnaux,3)
   call quick_sort(sort_array,ksort)

   !Re-assign matrix with the order given by ksort
   gl_aux2(:,:)=gl_aux(ksort,:)
   ksorti=huge(0)
   sorti=huge(0)
   counti=1 
   count2=1
   count3=1
   do i=1,nnaux
      if ( gl_aux2(i,1) == counti) then
         sort_chico(count2)=gl_aux2(i,2)
         k_chica(count2)=gl_aux2(i,3)
         count2=count2+1
      else
         ksorti(1:count2-1)=k_chica(1:count2-1)
         sorti(1:count2-1)=sort_chico(1:count2-1)
         call quick_sort(sorti(1:count2-1),ksorti(1:count2-1))
         gl_aux2(count3:count3+count2-2,2)=sorti(1:count2-1)
         gl_aux2(count3:count3+count2-2,3)=ksorti(1:count2-1)

         count3=count3+count2-1
         counti= counti+1
         count2=1
         sort_chico(count2)=gl_aux2(i,2)
         k_chica(count2)=gl_aux2(i,3)
         count2=count2+1
      endif
   enddo

   allocate(global_mat_aux2(nnaux))
   global_mat_aux2(:)=global_mat_aux(gl_aux2(:,3))
   do i=2,nnaux
      if ( (gl_aux2(i,1) == gl_aux2(i-1,1)) .AND. (gl_aux2(i,2) == gl_aux2(i-1,2)) ) then
         global_mat_aux2(i-1)=global_mat_aux2(i-1)+global_mat_aux2(i)
         gl_aux2(i,1) =0            
         gl_aux2(i,2) =0        
      endif
   enddo


   ! Add all the elements that share the same  value
   nnz=0
   do i=1,nnaux
      if ( gl_aux2(i,1) == 0 .or. gl_aux2(i,2) == 0 ) cycle
      if (gl_aux2(i,2) < gl_aux2(i,1) ) cycle
      nnz=nnz+1
      global_mat(nnz)= global_mat_aux2(i)
      gli(nnz)=gl_aux2(i,1)
      glj(nnz)=gl_aux2(i,2)
   
   enddo
   nrow = total_l_nod

! Check if apj_ser is allocated; if not, allocate it
   if (.not. allocated(apj_ser)) then
      allocate(apj_ser(nnz, nsub_proc))
   endif
   
   ! Check if ip, ij, apj are allocated; if not, allocate them
   if (.not. allocated(ip) .or. .not. allocated(ij) .or. .not. allocated(apj)) then
      allocate(ip(nnz), ij(nnz), apj(nnz))
   endif

   ip(1:nnz)=gli(1:nnz)
   ij(1:nnz)=glj(1:nnz)
   apj(1:nnz)=global_mat(1:nnz)

   ! Assign values to a slice of apj_ser from apj
   !apj_ser(1:nnz, sub) = apj(1:nnz)

  IF (.NOT.(use_pod)) THEN
  ! K = KRB
  reduced_stiffnes_coo(sub, 1)%ip_rb = ip(1:nnz)
  reduced_stiffnes_coo(sub, 1)%ij_rb = ij(1:nnz)
  reduced_stiffnes_coo(sub, 1)%apj_rb = apj(1:nnz)

  reduced_stiffnes_coo(sub, 2)%ip_rb = ip(1:nnz)
  reduced_stiffnes_coo(sub, 2)%ij_rb = ij(1:nnz)
  reduced_stiffnes_coo(sub, 2)%apj_rb = apj(1:nnz)
  ENDIF

    
end subroutine mat_coef_global


subroutine mat_coefl(j,dy,dz,mat_loc,i)

      USE sizes
      USE commons_3d
      USE globalvars
        
            
      integer :: i,j,k,l,dy,dz
      complex(DP)  :: m(42)
      complex(DP)  :: temp, temp2, unity, alphawjl, alphaejl, alphasjk, alphanjk, alphabkl, alphafkl
      complex(DP) :: mat_loc(1:12,1:12)
  
  
      integer :: igx,igz,igy
  
      m(:) = (0.d0,0.d0)
      unity = (1.d0,0.d0)

      !  computes Global indices
   
      igx= j       
      igy= vblockid(1,sub)*ny + dy
      igz= ((nsz-1)-vblockid(2,sub))*nz+ dz
      
     
      i=matnumloc(igx,dy,dz)

      ! compute Global/Local indeces

      k = dy + FLAG_globy*(vblockid(1,sub)*ny)
      l = dz + FLAG_globz*(((nsz-1)-vblockid(2,sub))*nz)


      h3 = hx(igx)*hy(igy)*hz(igz)
      hxz = hx(igx)*hz(igz)
      hxy = hx(igx)*hy(igy)
      hyz = hy(igy)*hz(igz)
      hxz_y = hxz/hy(igy)
      hxy_z = hxy/hz(igz)
      hyz_x = hyz/hx(igx)
      alphawjl = hxz*(deltaw(k))*alphaw(j,l)*(unity-aim)
      alphaejl = hxz*(deltae(k))*alphae(j,l)*(unity-aim)
      alphasjk = hxy*(deltas(l))*alphas(j,k)*(unity-aim)
      alphanjk = hxy*(deltan(l))*alphan(j,k)*(unity-aim)
      alphabkl = hyz*(deltab(j))*alphab(k,l)*(unity-aim)
      alphafkl = hyz*(deltaf(j))*alphaf(k,l)*(unity-aim)
  
  !c    Fila 1
  
      temp2 = f1*h3*sig(matnum(igx,igy,igz))
  
      temp = temp2-coef_i*(f2*(hxy_z+hxz_y)+hxz_y)
      m(1)= temp + alphawjl + hxz*(1-deltabw(k))*bet_w(j,k,l)
      m(7)= temp + alphaejl + hxz*(1-deltabe(k))*bet_e(j,k,l) !c    Fila 2
  
      temp = f3*h3*sig(matnum(igx,igy,igz))
      m(2)  = temp -coef_i*(f2*(hxy_z+hxz_y)-hxz_y)
      m(13) = temp -coef_i*(f2*(hxy_z+hxz_y)-hxy_z)
      m(20) = temp -coef_i*(f2*(hxy_z+hyz_x)-hyz_x) !c    Fila 5
      m(27) = temp -coef_i*(f2*(hxz_y+hyz_x)-hxz_y) !c    Fila 7
      m(38) = temp -coef_i*(f2*(hxy_z+hyz_x)-hxy_z) !c    Fila 9
      m(41) = temp -coef_i*(f2*(hxz_y+hyz_x)-hyz_x) !c    Fila 11
  
      temp = f4*h3*sig(matnum(igx,igy,igz))
      m(3)  = temp +coef_i*f2*(hxy_z+hxz_y)
      m(21) = temp +coef_i*f2*(hxy_z+hyz_x) !c    Fila 5
      m(30) = temp +coef_i*f2*(hxz_y+hyz_x) !c    Fila 7
  
      m(4)=m(3)
  
      m(5)=coef_i*hz(igz)
  
      m(6)=-m(5)
  
  !c    Fila 2
  
      m(8)=m(3)
           
      m(9)=m(4)
  
      m(10)=m(6)
  
      m(11)=m(5)
  
  !c    Fila 3
  
      temp = temp2-coef_i*(f2*(hxy_z+hxz_y)+hxy_z)
      m(12)= temp + alphasjk + hxy*(1-deltabs(l))*bet_s(j,k,l)
      m(16)= temp + alphanjk + hxy*(1-deltabn(l))*bet_n(j,k,l)!c    Fila 4
  
      m(14) = coef_i*hy(igy)
  
      m(15) = -m(14)
  
  !c    Fila 4
  
      m(17) = m(15)
  
      m(18) = m(14)
  
  !c    Fila 5
  
      temp = temp2-coef_i*(f2*(hxy_z+hyz_x)+hyz_x)
      m(19)= temp + alphabkl
      m(23)= temp + alphafkl !c    Fila 6
  
      m(22) = m(21)
  
  !c    Fila 6
  
      m(24) = m(21)
  
      m(25) = m(21)
  
  !c    Fila 7
  
      temp = temp2-coef_i*(f2*(hxz_y+hyz_x)+hxz_y)
      m(26)= temp + alphawjl + hxz*(1-deltabw(k))*bet_w(j,k,l)
      m(32)= temp + alphaejl + hxz*(1-deltabe(k))*bet_e(j,k,l) !c    Fila 8
  
      m(28) = hx(igx)*coef_i
  
      m(29) = - m(28)
  
      m(31) = m(30)
  
  !c    Fila 8
  
      m(33) = m(29)
  
      m(34) = m(28)
  
      m(35) = m(30)
  
      m(36) = m(30)
  
      temp = temp2-coef_i*(f2*(hxy_z+hyz_x)+hxy_z)
      m(37)= temp + alphasjk + hxy*(1-deltabs(l))*bet_s(j,k,l) !c    Fila 9
      m(39)= temp + alphanjk + hxy*(1-deltabn(l))*bet_n(j,k,l) !c    Fila 10
  
      temp = temp2-coef_i*(f2*(hxz_y+hyz_x)+hyz_x)
      m(40)= temp + alphabkl !c    Fila 11
      m(42)= temp + alphafkl !c    Fila 12
  
  !    Armo la matriz local en un formato que me va a servir para armar la global
  
      mat_loc(:,:) = (0.d0,0.d0)  
      mat_loc(1,1)=m(1)
      mat_loc(1,2)=m(2)
      mat_loc(1,3)=m(3)
      mat_loc(1,4)=m(4)
      mat_loc(1,5)=m(5)
      mat_loc(1,6)=m(6)
      mat_loc(2,1)=m(2)
      mat_loc(2,2)=m(7)
      mat_loc(2,3)=m(8)
      mat_loc(2,4)=m(9)
      mat_loc(2,5)=m(10)
      mat_loc(2,6)=m(11)
      mat_loc(3,1)=m(3)
      mat_loc(3,2)=m(8)
      mat_loc(3,3)=m(12)
      mat_loc(3,4)=m(13)
      mat_loc(3,11)=m(14)
      mat_loc(3,12)=m(15)
      mat_loc(4,1)=m(4)
      mat_loc(4,2)=m(9)
      mat_loc(4,3)=m(13)
      mat_loc(4,4)=m(16)
      mat_loc(4,11)=m(17)
      mat_loc(4,12)=m(18)
      mat_loc(5,1)=m(5)
      mat_loc(5,2)=m(10)
      mat_loc(5,5)=m(19)
      mat_loc(5,6)=m(20)
      mat_loc(5,9)=m(21)
      mat_loc(5,10)=m(22)
      mat_loc(6,1)=m(6)
      mat_loc(6,2)=m(11)
      mat_loc(6,5)=m(20)
      mat_loc(6,6)=m(23)
      mat_loc(6,9)=m(24)
      mat_loc(6,10)=m(25)
      mat_loc(7,7)=m(26)
      mat_loc(7,8)=m(27)
      mat_loc(7,9)=m(28)
      mat_loc(7,10)=m(29)
      mat_loc(7,11)=m(30)
      mat_loc(7,12)=m(31)
      mat_loc(8,7)=m(27)
      mat_loc(8,8)=m(32)
      mat_loc(8,9)=m(33)
      mat_loc(8,10)=m(34)
      mat_loc(8,11)=m(35)
      mat_loc(8,12)=m(36)
      mat_loc(9,5)=m(21)
      mat_loc(9,6)=m(24)
      mat_loc(9,7)=m(28)
      mat_loc(9,8)=m(33)
      mat_loc(9,9)=m(37)
      mat_loc(9,10)=m(38)
      mat_loc(10,5)=m(22)
      mat_loc(10,6)=m(25)
      mat_loc(10,7)=m(29)
      mat_loc(10,8)=m(34)
      mat_loc(10,9)=m(38)
      mat_loc(10,10)=m(39)
      mat_loc(11,3)=m(14)
      mat_loc(11,4)=m(17)
      mat_loc(11,7)=m(30)
      mat_loc(11,8)=m(35)
      mat_loc(11,11)=m(40)
      mat_loc(11,12)=m(41)
      mat_loc(12,3)=m(15)
      mat_loc(12,4)=m(18)
      mat_loc(12,7)=m(31)
      mat_loc(12,8)=m(36)
      mat_loc(12,11)=m(41)
      mat_loc(12,12)=m(42)

 end subroutine mat_coefl
     FUNCTION ep(zq)

        USE sizes 
        USE commons_3d
        USE globalvars
        IMPLICIT NONE
        
        REAL(DP) :: z,a,zq,ww,www,hh1,hh2,hh3,hh4,hh5
        COMPLEX(DP) :: kk1,kk2,kk3,a0,a1,a2,a3,a4,a5,b1,b2,kk0,b0,b3,b4,ep
        COMPLEX(DP) :: ga45,eta45,ga345,eta345,ga2345,eta2345,ga12345
        COMPLEX(DP) :: eta12345,ga012345,kk4,kk5,tmp
        
!       z = zmax - zq
!       z=znod(ngz+1)-zq
        z=znod(izsup)-zq
!       z=(e1+e2+e3+e4+e5-zq)

   
      
        a=1.d0
! c     a0=(1.d0,0.d0)*a

        hh2= e2
        hh3= hh2+e3 
        hh4= hh3+e4
        hh5= hh4+e5

        kk0= aim*cdsqrt(aim*sigma0*amu*omega)
        kk1= aim*cdsqrt(aim*sigma1*amu*omega)
        kk2= aim*cdsqrt(aim*sigma2*amu*omega)
        kk3= aim*cdsqrt(aim*sigma3*amu*omega)
        kk4= aim*cdsqrt(aim*sigma4*amu*omega)
        kk5= aim*cdsqrt(aim*sigma5*amu*omega)

        ga45= cdexp(2.d0*aim*kk4*hh5)*(1.d0-(kk5/kk4))/(1.d0+(kk5/kk4))

        tmp= aim*kk4*hh4      

        eta45= (cdexp(tmp)-ga45*cdexp(-tmp))/(cdexp(tmp)+ga45*cdexp(-tmp))

        ga345= cdexp(2.d0*aim*kk3*hh4)*(1.d0-(kk4/kk3)*eta45)/ &
             (1.d0+(kk4/kk3)*eta45)

        tmp= aim*kk3*hh3      

        eta345= (cdexp(tmp)-ga345*cdexp(-tmp))/ &
             (cdexp(tmp)+ga345*cdexp(-tmp))

        ga2345= cdexp(2.d0*aim*kk2*hh3)*(1.d0-(kk3/kk2)*eta345)/ &
             (1.d0+(kk3/kk2)*eta345)


        tmp= aim*kk2*hh2      

        eta2345= (cdexp(tmp)-ga2345*cdexp(-tmp))/ &
             (cdexp(tmp)+ga2345*cdexp(-tmp))

        ga12345= cdexp(2.d0*aim*kk1*hh2)* &
             (1.d0-eta2345*(kk2/kk1))/(1.d0+eta2345*(kk2/kk1))


        eta12345= ((1.d0,0.d0)-ga12345)/((1.d0,0.d0)+ga12345) 

        ga012345= (1.d0-eta12345*(kk1/kk0))/(1.d0+eta12345*(kk1/kk0))   


        a1= a*(1.d0/(1.d0+ga12345))

        b1= ga12345*a1

        a0= a1*(1.d0+ga12345)/(1.d0+ga012345)         

        b0= ga012345*a0
      
      
        a2= a1*(cdexp(aim*kk1*hh2)+ga12345*cdexp(-aim*kk1*hh2))/ &
             (cdexp(aim*kk2*hh2)+ga2345*cdexp(-aim*kk2*hh2)) 

        b2= ga2345*a2

        a3= a2*(cdexp(aim*kk2*hh3)+ga2345*cdexp(-aim*kk2*hh3))/ &
             (cdexp(aim*kk3*hh3)+ga345*cdexp(-aim*kk3*hh3))

        b3= ga345*a3   
 
        a4= a3*(cdexp(aim*kk3*hh4)+ga345*cdexp(-aim*kk3*hh4))/ &
             (cdexp(aim*kk4*hh4)+ga45*cdexp(-aim*kk4*hh4))

        b4= ga45*a4   
 
        a5= a4*(cdexp(aim*(kk4-kk5)*hh5)+ga45*cdexp(-aim*(kk4+kk5)*hh5))

   
        IF(z.GE.(hh5))THEN
           ep = a5*cdexp(aim*kk5*z)
        ELSEIF (z.GE.hh4)THEN
           ep = a4*cdexp(aim*kk4*z) + b4*cdexp(-aim*kk4*z)
        ELSEIF(z.GE.hh3)THEN       
           ep = a3*cdexp(aim*kk3*z) + b3*cdexp(-aim*kk3*z)
        ELSEIF(z.GE.hh2)THEN       
           ep = a2*cdexp(aim*kk2*z) + b2*cdexp(-aim*kk2*z)
        ELSEIF(z.ge.0.d0.AND.z.LT.hh2)THEN       
           ep = a1*cdexp(aim*kk1*z) + b1*cdexp(-aim*kk1*z)
        ELSEIF(z.LT.0.d0)THEN
           ep = a0*cdexp(aim*kk0*z) + b0*cdexp(-aim*kk0*z)            
        ENDIF

      
        
        RETURN
      END FUNCTION ep


!**********************************************************************


    FUNCTION hp(zq)
       
       USE sizes 
       USE commons_3d
       USE globalvars
       IMPLICIT NONE
       
      REAL(DP) :: z,a,zq,ww,www,hh1,hh2,hh3,hh4,hh5
      COMPLEX(DP) :: kk1,kk2,kk3,a0,a1,a2,a3,a4,a5,b1,b2,kk0,b0,b3,b4
      COMPLEX(DP) :: ga45,eta45,ga345,eta345,ga2345,eta2345,ga12345
      COMPLEX(DP) :: eta12345,ga012345,kk4,kk5,tmp,hp

!      z = zmax - zq
!       z=znod(ngz+1)-zq
      z= znod(izsup)-zq

      a= 1.d0
      hh2= e2
      hh3= hh2+e3 
      hh4= hh3+e4
      hh5= hh4+e5

      kk0= aim*cdsqrt(aim*sigma0*amu*omega)
      kk1= aim*cdsqrt(aim*sigma1*amu*omega)
      kk2= aim*cdsqrt(aim*sigma2*amu*omega)
      kk3= aim*cdsqrt(aim*sigma3*amu*omega)
      kk4= aim*cdsqrt(aim*sigma4*amu*omega)
      kk5= aim*cdsqrt(aim*sigma5*amu*omega)

      ga45= cdexp(2.d0*aim*kk4*hh5)*(1.d0-(kk5/kk4))/(1.d0+(kk5/kk4))

      tmp= aim*kk4*hh4      

      eta45= (cdexp(tmp)-ga45*cdexp(-tmp))/(cdexp(tmp)+ga45*cdexp(-tmp))

      ga345= cdexp(2.d0*aim*kk3*hh4)*(1.d0-(kk4/kk3)*eta45)/ &
           (1.d0+(kk4/kk3)*eta45)

      tmp= aim*kk3*hh3      

      eta345= (cdexp(tmp)-ga345*cdexp(-tmp))/ &
           (cdexp(tmp)+ga345*cdexp(-tmp))

      ga2345= cdexp(2.d0*aim*kk2*hh3)*(1.d0-(kk3/kk2)*eta345)/ &
           (1.d0+(kk3/kk2)*eta345)


      tmp= aim*kk2*hh2      

      eta2345= (cdexp(tmp)-ga2345*cdexp(-tmp))/ &
           (cdexp(tmp)+ga2345*cdexp(-tmp))

      ga12345= cdexp(2.d0*aim*kk1*hh2)* &
           (1.d0-eta2345*(kk2/kk1))/(1.d0+eta2345*(kk2/kk1))


      eta12345= ((1.d0,0.d0)-ga12345)/((1.d0,0.d0)+ga12345) 

      ga012345= (1.d0-eta12345*(kk1/kk0))/(1.d0+eta12345*(kk1/kk0))   


      a1= a*(1.d0/(1.d0+ga12345))

      b1= ga12345*a1

      a0=a1*(1.d0+ga12345)/(1.d0+ga012345)         

      b0= ga012345*a0
      
      
      a2= a1*(cdexp(aim*kk1*hh2)+ga12345*cdexp(-aim*kk1*hh2))/ &
           (cdexp(aim*kk2*hh2)+ga2345*cdexp(-aim*kk2*hh2)) 

      b2= ga2345*a2

      a3= a2*(cdexp(aim*kk2*hh3)+ga2345*cdexp(-aim*kk2*hh3))/ &
           (cdexp(aim*kk3*hh3)+ga345*cdexp(-aim*kk3*hh3))

      b3= ga345*a3   
 
      a4= a3*(cdexp(aim*kk3*hh4)+ga345*cdexp(-aim*kk3*hh4))/ &
           (cdexp(aim*kk4*hh4)+ga45*cdexp(-aim*kk4*hh4))

      b4= ga45*a4   
 
      a5= a4*(cdexp(aim*(kk4-kk5)*hh5)+ga45*cdexp(-aim*(kk4+kk5)*hh5))


   
      IF(z.GE.(hh5))THEN
         hp = -(kk5/(amu*omega))*a5*cdexp(aim*kk5*z)
      ELSEIF(z.GE.hh4)THEN
         hp =-(kk4/(amu*omega))*(a4*cdexp(aim*kk4*z)-b4*cdexp(-aim*kk4*z))
      ELSEIF(z.GE.hh3)THEN       
         hp =-(kk3/(amu*omega))*(a3*cdexp(aim*kk3*z)-b3*cdexp(-aim*kk3*z))         
      ELSEIF(z.GE.hh2)THEN       
         hp =-(kk2/(amu*omega))*(a2*cdexp(aim*kk2*z)-b2*cdexp(-aim*kk2*z))        
      ELSEIF(z.GE.0.d0.AND.z.LT.hh2)THEN       
         hp =-(kk1/(amu*omega))*(a1*cdexp(aim*kk1*z)-b1*cdexp(-aim*kk1*z))      
      ELSEIF(z.LT.0.d0)THEN
         hp =-(kk0/(amu*omega))*(a0*cdexp(aim*kk0*z)-b0*cdexp(-aim*kk0*z))         
      ENDIF
 

      RETURN
    END FUNCTION hp

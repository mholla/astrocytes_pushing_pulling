! author: Karan Taneja
! unit system: mm-N-MPa
! dimension: plane strain
! documentation: 'Astrocytes in white matter respond to tensile cues during cortical folding: a numerical study'

!**********************************************************************
      module GlobalStorage
      ! number of nodes in the input file
      ! Important as simulation will crash without updating it.   
          real*8 inicoord(46000,3) 
      end module  
***********************************************************************

************************************************************************
      subroutine vumat (
c Read only -
     +     jblock, ndir, nshr, nstatev, nfieldv, nprops,lanneal, 
     +     stepTime, totalTime, dt, cmname, coordMp,charLength, 
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     +     charLength(*), strainInc(*),
     +     relSpinInc(*), tempOld(*),
     +     stretchOld(*),
     +     defgradOld(*),
     +     fieldOld(*), stressOld(*),
     +     stateOld(*), enerInternOld(*),
     +     enerInelasOld(*), tempNew(*),
     +     stretchNew(*),
     +     defgradNew(*),
     +     fieldNew(*),
     +     stressNew(*), stateNew(*),
     +     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname
      character*256 WHIT,GRAY

      parameter (     
     +     i_umt_nblock = 1,
     +     i_umt_npt    = 2,
     +     i_umt_layer  = 3,
     +     i_umt_kspt   = 4,
     +     i_umt_noel   = 5 )


      !--------------------------------------------------------
      ! 
      ! call particular user material to perform the analysis 
      ! 
      IF (CMNAME(1:4) .EQ. 'WHIT') THEN

      !
      ! this is white matter 
      !
      call  vumatXtrArg_white (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp,charLength, 
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
      !
      !
      ELSE IF(CMNAME(1:4) .EQ. 'GRAY') THEN
      !
      ! this is gray matter 
      !
      call  vumatXtrArg_gray (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp,charLength, 
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
      !
      !
      Endif

     
      end subroutine vumat
***********************************************************************
      subroutine vumatXtrArg_white (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops,lanneal, 
     +     stepTime, totalTime, dt, cmname, coordMp,charLength, 
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )

      use GlobalStorage

      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c
      character*80 cmname

      integer i,km

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_tau(3,3)
      real*8 sigma_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fe_tau(3,3)
      real*8 pwrinct,stress_power
      real*8 sigma_rot(3,3),rot_matrix(3,3),N_R(2,1)
      real*8 matProps(nprops),sigma_rad,sigma_tan
      real*8 theta_dot_1,f_2,zeta,ctheta,stheta
      real*8 coordx,coordy,coordz,thetag_t,thetag_tau,maj_axis,min_axis
      real*8 maj_min_ratio      



      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)

!! WM Ellipse params
      maj_axis = props(12) ! WM major axis (mm)
      min_axis = props(13) ! WM minor axis (mm)

      ! Pour initial coordinates into the global variable matrix 
      !
      if (totalTime.lt. 0.1) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if 



      ! Identity matrix for later use.
      ! 
      call onem(Iden)


      ! 
      ! START LOOP OVER MATERIAL POINTS:
      ! 
      do km=1,nblock
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! Dummy step, initalize state variables
            
            stateOld(km,1)   = one ! growth parameter at t=0
         endif

         ! Read old state variables
         
         thetag_t = stateOld(km,1) ! growth parameter at time t




          coordx = inicoord(nElement(km),1)
          coordy = inicoord(nElement(km),2)
          coordz = inicoord(nElement(km),3)


         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_white(matProps,nprops,F_tau,-1.0,sigma_tau,thetag_t,thetag_tau,
     +                       coordx,coordy,coordz,totalTime,
     +                       theta_dot_1,f_2)     

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_white(matProps,nprops,F_tau,dt,sigma_tau,thetag_t,thetag_tau,
     +                       coordx,coordy,coordz,totalTime,
     +                       theta_dot_1,f_2)     

         endif
         !---------------------------------------------------------------



         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         sigma_tau = matmul(transpose(R_tau),matmul(sigma_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = sigma_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = sigma_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = sigma_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = sigma_tau(1,3)
               endif
            endif
         endif

         !! Rotate stress tensor to get radial and tangential outputs
         ! S' = R.S.R^T


         ! Find normal vector

         maj_min_ratio = maj_axis/min_axis

         N_R(1,1) = 2.0*coordx/maj_axis**2.0
         N_R(2,1) = 2.0*coordy/min_axis**2.0
         zeta = atan(N_R(2,1)/N_R(1,1))

         ! Create rotation matrix
         ctheta = cos(zeta)
         stheta = sin(zeta)
         rot_matrix(1,1) = ctheta
         rot_matrix(1,2) = stheta
         rot_matrix(1,3) = 0.0
         rot_matrix(2,1) = -stheta
         rot_matrix(2,2) = ctheta
         rot_matrix(2,3) = 0.0
         rot_matrix(3,1) = 0.0
         rot_matrix(3,2) = 0.0
         rot_matrix(3,3) = 1.0

         ! Rotate the stress tensor
         sigma_rot = matmul(matmul(rot_matrix,sigma_tau),transpose(rot_matrix))

         ! Get radial and tangential components
         sigma_rad = sigma_rot(1,1)
         sigma_tan = sigma_rot(2,2)


         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx ! 
         stateNew(km,3) = coordy ! 
         stateNew(km,4) = coordz !
         call mdet(F_tau,detF)
         stateNew(km,5) = detF   
         stateNew(km,6) = (sigma_tau(1,1)+sigma_tau(2,2)+sigma_tau(3,3))/three
         stateNew(km,7) = theta_dot_1
         stateNew(km,8) = f_2   
         stateNew(km,9) = sigma_rad
         stateNew(km,10) = sigma_tan   

         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((stressOld(km,i)+stressNew(km,i))*
     +           strainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((stressOld(km,ndir+1)+stressNew(km,ndir+1))*
     +           strainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((stressOld(km,ndir+1) + stressNew(km,ndir+1))*
     +           strainInc(km,ndir+1)) +
     +           ((stressOld(km,ndir+2)+ stressNew(km,ndir+2)) *
     +           strainInc(km,ndir+2))+
     +           ((stressOld(km,ndir+3) + stressNew(km,ndir+3))*
     +           strainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_white
***********************************************************************
      subroutine vumatXtrArg_gray (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, 
     +     stepTime, totalTime, dt, cmname, coordMp, charLength, 
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )


      use GlobalStorage
      include 'vaba_param.inc'


      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)
c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c

      character*80 cmname

      integer i,km

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_tau(3,3)
      real*8 sigma_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fe_tau(3,3)
      real*8 pwrinct,stress_power
      real*8 matProps(nprops)
      real*8 thetag_t,thetag_tau
      real*8 coordx,coordy,coordz
      real*8 N_R(3,1)
 

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)


      ! pour initial coordinates into the global variable
      if (totalTime.lt. 0.1) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if


      ! Identity matrix for later use.
      !
      call onem(Iden)


      !
      ! START LOOP OVER MATERIAL POINTS:
      !
      do km=1,nblock

           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! Dummy step, initalize state variables
            
            stateOld(km,1)   = one ! growth parameter at t=0
         endif

         ! Read old state variables
         
         thetag_t = stateOld(km,1) ! growth parameter at time t




         ! reads in the original coordinates 
          coordx = inicoord(nElement(km),1)
          coordy = inicoord(nElement(km),2)
          coordz = inicoord(nElement(km),3)



         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_gray(matProps,nprops,F_tau,-1.0,sigma_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,N_R,totalTime)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_gray(matProps,nprops,F_tau,dt,sigma_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,N_R,totalTime)

         endif
         !---------------------------------------------------------------
         

         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz
         call mdet(F_tau,detF)

         stateNew(km,5) = detF 


         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         sigma_tau = matmul(transpose(R_tau),matmul(sigma_tau,R_tau))
         do i=1,ndir
            stressNew(km,i) = sigma_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = sigma_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = sigma_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = sigma_tau(1,3)
               endif
            endif
         endif
         stateNew(km,6) = (sigma_tau(1,1)+sigma_tau(2,2)+sigma_tau(3,3))/three

         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((stressOld(km,i)+stressNew(km,i))*
     +           strainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((stressOld(km,ndir+1)+stressNew(km,ndir+1))*
     +           strainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((stressOld(km,ndir+1) + stressNew(km,ndir+1))*
     +           strainInc(km,ndir+1)) +
     +           ((stressOld(km,ndir+2)+ stressNew(km,ndir+2)) *
     +           strainInc(km,ndir+2))+
     +           ((stressOld(km,ndir+3) + stressNew(km,ndir+3))*
     +           strainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_gray
***********************************************************************
      subroutine integ_white(Props,nprops,F_tau,dtime,sigma_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,totalTime,
     +                       theta_dot_1,f_2)
      implicit none

      character*256 ,fileName

      integer i,j,k,l,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_tau(3,3),sigma_tau(3,3)
      real*8 detF
      real*8 lambda,mu
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 alpha_bar, coordiff
      real*8 tmp
      real*8 Fginv(3,3)
      real*8 coordx,coordy,coordz,theta_dot_1,G_GM,gamma_1,f_phi
      real*8 T_1,totalTime
      real*8 majoraxis_reduced,minoraxis_reduced
      real*8 maj_min_ratio,maj_axis,min_axis
      real*8 T_2,f_H
      real*8 rad,psi,r_tilde,delta_bar,a_tilde
      real*8 f_2,N_gyri,gamma,alpha
      real*8 b_tilde,periods       



      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine,ten
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,ten=10.d0)

      ! standard deviation of gauss growth rate function
      alpha = 0.4d0

      ! Obtain WM material properties 
      !
       mu        = props(1)
       lambda    = props(2)
       G_GM      = props(3) ! grey matter growth rate
       delta_bar = props(4) ! scaled threshold for heaviside function
       alpha_bar = props(5) ! Smoothening for heaviside function used in Phase 3
       gamma_1   = props(6) ! G_wm/G_gm used in Phase 1
       T_1       = props(7) ! end of Phase 1
       N_gyri    = props(8) ! Number of proliferation zones
       gamma     = props(9) ! G_wm/G_gm in Phase 3 (pushing)
       b_tilde   = props(10) ! scaling factor used in calculating r_tilde
       T_2       = props(11) ! end of Phase 2
       maj_axis  = props(12) ! WM major axis (a)
       min_axis  = props(13) ! WM minor axis (b)

!! Setting up the growth rate calculations
       maj_min_ratio = maj_axis/min_axis

       a_tilde = b_tilde*maj_min_ratio

       ! For progenitor push effect
       majoraxis_reduced = maj_axis - a_tilde !white matter reduced to bring in the progenitor effect
       minoraxis_reduced = min_axis - b_tilde 


      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)

       psi = atan(coordy/coordx)
       rad = sqrt(coordx**2.0 + coordy**2.0)
       r_tilde = rad/sqrt((majoraxis_reduced*cos(psi))**2.0 + (minoraxis_reduced*sin(psi))**2.0)


       coordiff = (r_tilde - delta_bar)*one
       
       f_2 = sin(four*psi*(N_gyri - half)) + one

       call gauss(r_tilde,delta_bar,alpha,f_phi)
       
       theta_dot_1 = (G_GM*gamma_1)*half*f_phi*f_2 ! Scaled so that growth rate is highest at the grey matter layer




      if ((totalTime.le.T_1)) then !Push effect BEFORE grow time


     !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! iso 
      Fg_tau  = (thetag_tau**third)*Iden
      

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je




         return
      endif      
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       thetag_tau = thetag_t + (theta_dot_1)*dtime 


      ! update  kinematics 
      ! 

      Fg_tau  = (thetag_tau**third)*Iden
      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je

      else if ((totalTime.le.T_2).and.(totalTime.ge.T_1)) then ! Between the two push phases

      thetag_tau = thetag_t 


      ! update  kinematics 
      ! 

      Fg_tau  = (thetag_tau**third)*Iden
      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je      



      endif

      if ((totalTime.ge.T_2)) then !Push effect after grow time
      
C       ! KT: If the totalTime is > T_2, push effect from astrocytes is in play 


       call Hhat(coordiff,alpha_bar,f_H)
       
       theta_dot_1 = (G_GM*gamma)*half*f_H*f_2 ! Scaled so that growth rate is highest at the grey matter layer


       thetag_tau = thetag_t + (theta_dot_1)*dtime 


      Fg_tau  = (thetag_tau**third)*Iden
      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je
      endif


      end subroutine integ_white
****************************************************************************
      subroutine integ_gray(Props,nprops,F_tau,dtime,sigma_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,
     +                      N_R,totalTime)

      implicit none



      integer nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_tau(3,3),sigma_tau(3,3)
      real*8 detF
      real*8 N_R(3,1)
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 mu,lambda,G_GM
      real*8 coordx,coordy,coordz,tmp
      real*8 Fginv(3,3)
      real*8 T_1,totalTime
      real*8 maj_axis,min_axis

      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine,Pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,Pi=3.1415926d0)

      ! Obtain material properties
      !
       mu        = props(1)
       lambda    = props(2)
       G_GM      = props(3) ! grey matter growth rate
       T_1       = props(4) ! Time after which gray matter starts to grow in Phase 2.
       maj_axis  = props(5) ! WM major axis (a)
       min_axis  = props(6) ! WM minor axis (b)

      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)




      ! obtain referential surface outnormal of an elliptical surface
      N_R(1,1) = 2.0*coordx/maj_axis**2.0
      N_R(2,1) = 2.0*coordy/min_axis**2.0
      N_R(3,1) = 0.0

      tmp = sqrt(N_R(1,1)**2.0 + N_R(2,1)**2.0 + N_R(3,1)**2.0)

      N_R = N_R/tmp   
  

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(N_R,transpose(N_R))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je


         return
      endif   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!



      ! KT: GM grows only after T_1 to allow progenitors to grow in WM

       if (totalTime.le.T_1) then
           thetag_tau = one
       else
           thetag_tau = thetag_t + (G_GM)*dtime 
       endif


      ! update  kinematics 
      ! area 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(N_R,transpose(N_R))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      sigma_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je


      end subroutine integ_gray
C**********************************************************************
C   THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
****************************************************************************
      subroutine Hhat(x,alpha,y)

      implicit none

      real*8 x,alpha,y
      
      y = (exp(alpha*x))/(1.0 + exp(alpha*x))

      end subroutine Hhat
****************************************************************************
      subroutine gauss(x,mean,sigma,y)

      implicit none

      real*8 x,mean,sigma,y
      real*8 half,Pi,two
      parameter(half=0.5d0,two=2.d0,Pi=3.1415926d0)
      

      y = exp(-half * ((x - mean) / sigma)**2) / (sigma * dsqrt(two*Pi))

      end subroutine gauss      


C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END



C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END





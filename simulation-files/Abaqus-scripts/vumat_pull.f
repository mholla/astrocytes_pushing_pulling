************************************************************************
!  
!  
! A vumat for manuscript: 'Astrocytes in white matter respond to tensile
! cues during cortical folding: a numerical study'
! 

!**********************************************************************
      module GlobalStorage
      ! number of nodes in the input file   
          real*8 inicoord(46000,3) 
      end module  

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
c lanneal is a default parameter set to zero internally. Refer manual.
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
     +     ndir, nshr, nstatev, nfieldv, nprops,lanneal, 
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
     +     ndir, nshr, nstatev, nfieldv, nprops,lanneal, 
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


c$$$        implicit none ! This is used during compilation testing to make
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

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_tau(3,3)
      real*8 T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fe_tau(3,3)
      real*8 pwrinct,stress_power
      real*8 rot_stress(3,3),rot_matrix(3,3),normal_vector(2,1)
      real*8 matProps(nprops),rad_stress,tan_stress
      real*8 white_growth,ang_factor,rot_angle,ctheta,stheta
      real*8 coordx,coordy,coordz,thetag_t,thetag_tau,maj_axis,min_axis
      real*8 maj_min_ratio  



      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)
!! WM Ellipse params
      maj_axis = 34.5
      min_axis = 28.5

      ! Pour initial coordinates into the global variable matrix 
      !
      if (totaltime.lt. 0.1) then
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
            call integ_white(matProps,nprops,F_tau,-1.0,T_tau,thetag_t,thetag_tau,
     +                       coordx,coordy,coordz,totalTime,
     +                       white_growth,ang_factor)     

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_white(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                       coordx,coordy,coordz,totalTime,
     +                       white_growth,ang_factor)     

         endif
         !---------------------------------------------------------------




         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif
         !! Rotate stress tensor to get radial and tangential outputs
         ! S' = R.S.R^T


         ! Find normal vector

         maj_min_ratio = maj_axis/min_axis

         normal_vector(1,1) = 2.0*coordx/(maj_min_ratio**2.0)
         normal_vector(2,1) = 2.0*coordy
         rot_angle = atan(normal_vector(2,1)/normal_vector(1,1))

         ! Create rotation matrix
         ctheta = cos(rot_angle)
         stheta = sin(rot_angle)
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
         rot_stress = matmul(matmul(rot_matrix,T_tau),transpose(rot_matrix))

         ! Get radial and tangential components
         rad_stress = rot_stress(1,1)
         tan_stress = rot_stress(2,2)


         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx ! 
         stateNew(km,3) = coordy ! 
         stateNew(km,4) = coordz !
         call mdet(F_tau,detF)
         stateNew(km,5) = detF   
         stateNew(km,6) = (T_tau(1,1)+T_tau(2,2)+T_tau(3,3))/three
         stateNew(km,7) = white_growth
         stateNew(km,8) = ang_factor   
         stateNew(km,9) = rad_stress
         stateNew(km,10) = tan_stress   



         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
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

c$$$        implicit none ! This is used during compilation testing to make
      use GlobalStorage
      include 'vaba_param.inc'

c$$$!        When implicit none is used during compilation testing, all the following
c$$$!        variables need to be defined.
c$$$        integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
c$$$        real*8  stepTime, totalTime, dt, coordMp,  props,
c$$$     +          density, strainInc, relSpinInc, tempOld, stretchOld,
c$$$     +          defgradOld, fieldOld, stressOld, stateOld, 
c$$$     +          enerInternOld, enerInelasOld, tempNew, stretchNew,
c$$$     +          defgradNew, fieldNew, stressNew, stateNew, 
c$$$     +          enerInternNew, enerInelasNew

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

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail,CurrentElement

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_tau(3,3)
      real*8 T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fe_tau(3,3)
      real*8 pwrinct,stress_power
      real*8 matProps(nprops)
      real*8 thetag_t,thetag_tau
      real*8 coordx,coordy,coordz
      real*8 a0(3,1)
 

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)


      ! pour initial coordinates into the global variable
      if (totaltime.lt. 0.1) then
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


         ! calculate the current element 
         ! 
         CurrentElement  = nElement(km)-BaseElement


         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_gray(matProps,nprops,F_tau,-1.0,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,a0,totalTime)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_gray(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,a0,totalTime)

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
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif
         stateNew(km,6) = (T_tau(1,1)+T_tau(2,2)+T_tau(3,3))/three

         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_gray
***********************************************************************
      subroutine integ_white(Props,nprops,F_tau,dtime,T_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,totalTime,
     +                       white_growth,ang_factor)
      implicit none


      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,a,b,c,d,iterError,lenJobName,lenOutDir,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_tau(3,3),T_tau(3,3)
      real*8 detF
      real*8 Finv(3,3),B_tau(3,3)
      real*8 lambda,mu
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 Fginv(3,3)
      real*8 coordx,coordy,coordz,white_growth,Gctx,scale_fac,gauss_fac
      real*8 threshold,progen_grow_time,totalTime,pull_grow_time
      real*8 majoraxis_reduced,minoraxis_reduced,growth_crit
      real*8 nitl, thetag_dum, lnJe, trMe,reduction_major,reduction_minor
      real*8 res, dres, phig, dphig, xtol
      real*8 majorT,minorT,white_growth_tensile
      real*8 rad,phi,R_prime,scaled_r,scaled_thresh,maj_axis,min_axis
      real*8 ang_factor,N_gyri,scale_fac_pull
      real*8 periods,maj_min_ratio,gauss_std       



!!!!! Parameters
      !
      real*8 zero,one,two,half,three,third,nine,ten
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,ten=10.d0)


!! WM ellipse dimension
      maj_axis = 34.5d0 ! WM major axis
      min_axis = 28.5d0 ! WM minor axis


      xtol = 1.d-10 ! Tolerance for local newton iterations in growth parameter calculations 
      gauss_std = 0.4d0 ! Standard deviation of gauss function in Phase 1 growth rate

      growth_crit = zero ! critical growth criterion, Purely tensile growth in this study

      ! Obtain material properties 
      !
       mu        = props(1)
       lambda    = props(2)
       Gctx      = props(3) ! grey matter growth
       scaled_thresh = props(4) ! used to calculate parameter scaled_thresh \bar{\delta}
       scale_fac =  props(5) ! G_wm/G_gm used in Phase 1
       progen_grow_time = props(6) ! Time before which white matter grows in Phase1.
       N_gyri = props(7)! Number of proliferation zones
       scale_fac_pull = props(8) !scale_factor for pulling (gammahat)
       reduction_minor = props(9) ! Reducing from edge of WM (\tilde{b})
       pull_grow_time = props(10) !  Phase 3 pull effect starting after this time.
   


       reduction_major = reduction_minor*maj_min_ratio 

       majoraxis_reduced = maj_axis - reduction_major !white matter reduced to bring in the progenitor effect
       minoraxis_reduced = min_axis - reduction_minor 



      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)


      ! Compute the inverse of the deformation gradient
      !
      call m3inv(F_tau,Finv)



       phi = atan(coordy/coordx)
       rad = sqrt(coordx**2.0 + coordy**2.0)
       R_prime = sqrt((majoraxis_reduced*cos(phi))**2.0 + (minoraxis_reduced*sin(phi))**2.0)
       scaled_r = rad/R_prime


       periods = two*(two*N_gyri - one) !Periodic presence of progenitor growth

       ang_factor = sin(periods*phi)
       
       call gauss(scaled_r,scaled_thresh,gauss_std,gauss_fac)
       
       white_growth = half*(Gctx*scale_fac*one)*gauss_fac*(ang_factor+one) ! WM growth rate in Phase 1   
       white_growth_tensile = Gctx*scale_fac_pull*1000.0 !scaling with 1/mu_W, WM growth rate in Phase 3 


      if ((totalTime.le.progen_grow_time)) then !Push effect BEFORE grow time, Phase 1


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
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je




         return
      endif      
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      thetag_tau = thetag_t + (white_growth)*dtime 


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
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je

      else if ((totalTime.le.pull_grow_time).and.(totalTime.ge.progen_grow_time)) then 
      ! Between the progenitor push and astro pull phases,
      ! WM doesn't grow so theta_g doesn't change
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
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je

      endif

      if ((totalTime.ge.pull_grow_time)) then !Pull effect after grow time

      
      ! KT: If the totalTime is > pull_grow_time, tensile stress based growth 
      ! KT: Mandel stress derived growth criterion function
      ! KT: Taken from Ch 6 of Hitchhikers and associted code on Github

      ! Obtain dummy current growth var
      thetag_dum = thetag_t

      ! Start Newton Raphson Iteration
      nitl = 0
      
 200    continue      

        ! Get Growth tensor

        Fg_tau  = (thetag_dum**third)*Iden

        ! Get Fg^-1

        call m3inv(Fg_tau,Fginv) 

        ! Get Elastic F and associated deformation measures
      
        Fe_tau = matmul(F_tau,Fginv)

        Be_tau = matmul(Fe_tau,transpose(Fe_tau))

        B_tau = matmul(F_tau, transpose(F_tau))

        call mdet(Fg_tau,Jg)

        Je = detF/Jg

        lnJe = log(Je)

        ! Growth criterion based on trace of Kirchoff Stress

        trMe = 3*(lambda*lnJe - mu) + mu*(Be_tau(1,1) + Be_tau(2,2) + Be_tau(3,3))
        phig = trMe - growth_crit


        if (phig.le.zero) then ! no tensile growth
            ! print *,'No growth, trMe=',trME
            thetag_tau = thetag_dum
            res = zero
            dres = one

        else ! tensile growth
            ! print *,'Growth, trMe=',trME
            nitl = nitl + 1
            ! Growth criterion

            dphig = -1.d0/3.d0/thetag_dum*(2.d0*mu*(Be_tau(1,1) + Be_tau(2,2) + Be_tau(3,3)) + 9.d0*lambda)

            ! Non-linear residual
            res = thetag_dum - thetag_t - (white_growth_tensile)*phig*dtime
            dres = one - ((white_growth_tensile)*dphig)*dtime

            ! updated thetag
            thetag_tau = thetag_dum - res / dres

            ! Check for convergence

            if ((nitl.lt.50).and.(dabs(res).gt.xtol)) go to 200
            ! if (nitl.eq.50) print *, '>no convergence!!!! |r|=', dabs(res)

        endif
       ! ends local newton iteration 

      ! update  kinematics 


      Fg_tau  = (thetag_tau**third)*Iden

      call m3inv(Fg_tau,Fginv) 
      
      Fe_tau = matmul(F_tau,Fginv)

      Be_tau = matmul(Fe_tau,transpose(Fe_tau))

      B_tau = matmul(F_tau, transpose(F_tau))

      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je


      endif


      end subroutine integ_white
****************************************************************************
      subroutine integ_gray(Props,nprops,F_tau,dtime,T_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,
     +                      a0,totalTime)

      implicit none


      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_tau(3,3),T_tau(3,3)
      real*8 detF
      real*8 Finv(3,3),B_tau(3,3)
      real*8 a0(3,1)
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),fac,thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 mu_g,lambda_g,Gctx
      real*8 coordx,coordy,coordz,tmp
      real*8 Fginv(3,3)
      real*8 grow_time,totalTime

      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine,Pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,Pi=3.1415926d0)

      ! Obtain material properties 
      !
       mu_g      = props(1)
       lambda_g  = props(2)
       Gctx      = props(3)
       grow_time = props(4) ! Time after which gray matter starts to grow in Phase 2. 

      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)


      ! Compute the inverse of the deformation gradient
      !
      call m3inv(F_tau,Finv)


      ! obtain referential surface outnormal of an elliptical surface
      a0(1,1) = 2.0*coordx/1.2**2.0
      a0(2,1) = 2.0*coordy/1.0**2.0
      a0(3,1) = 0.0

      tmp = sqrt(a0(1,1)**2.0 + a0(2,1)**2.0 + a0(3,1)**2.0)

      a0 = a0/tmp   
  

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

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
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je


         return
      endif   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!



      ! KT: GM grows only after grow_time to allow progenitors to grow in WM

       if (totalTime.le.grow_time) then
           thetag_tau = one
       else
           thetag_tau = thetag_t + (Gctx)*dtime 
       endif


      ! update  kinematics 
      ! area 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

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
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je


      end subroutine integ_gray
****************************************************************************
****************************************************************************
C C**********************************************************************
C C   THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C C**********************************************************************
      subroutine gauss(x,mean,sigma,y)

      implicit none

      real*8 x,mean,sigma,y
      real*8 half,Pi,two
      parameter(half=0.5d0,two=2.d0,Pi=3.1415926d0)
      

      y = exp(-half * ((x - mean) / sigma)**2) / (sigma * dsqrt(two*Pi))

      end subroutine gauss      

C**********************************************************************
      SUBROUTINE ONEM(A)
      REAL*8 A(3,3)
      INTEGER I, J

      DO I = 1, 3
         DO J = 1, 3
            IF (I .EQ. J) THEN
               A(I,J) = 1.0D0
            ELSE
               A(I,J) = 0.0D0
            END IF
         END DO
      END DO

      RETURN
      END

C**********************************************************************
      SUBROUTINE MTRANS(A, ATRANS)
      REAL*8 A(3,3), ATRANS(3,3)
      INTEGER I, J

      DO I = 1, 3
         DO J = 1, 3
            ATRANS(J,I) = A(I,J)
         END DO
      END DO

      RETURN
      END


C**********************************************************************
      SUBROUTINE MDET(A, DET)
      REAL*8 A(3,3), DET

      DET =   A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     +       + A(1,3)*A(2,1)*A(3,2)
     +       - A(3,1)*A(2,2)*A(1,3) - A(3,2)*A(2,3)*A(1,1)
     +       - A(3,3)*A(2,1)*A(1,2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE M3INV(A, AINV)
      REAL*8 A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)
      INTEGER I, J

      CALL MDET(A, DET)
      IF (DET .EQ. 0.D0) THEN
         WRITE(*,10)
         STOP
      END IF

      CALL MCOFAC(A, ACOFAC)
      CALL MTRANS(ACOFAC, AADJ)

      DO I = 1, 3
         DO J = 1, 3
            AINV(I,J) = AADJ(I,J) / DET
         END DO
      END DO

10    FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR')

      RETURN
      END

C**********************************************************************
      SUBROUTINE MCOFAC(A, ACOFAC)
      REAL*8 A(3,3), ACOFAC(3,3)

      ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      ACOFAC(3,1) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
      ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

      RETURN
      END



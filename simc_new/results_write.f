	subroutine results_ntu_write(main,vertex,orig,recon,success)

	implicit none
	include 'radc.inc'
	include 'hbook.inc'
	include 'simulate.inc'

	integer i
	real*8	ntu(80)
	type(event_main):: main
	type(event):: vertex, orig, recon

!local (e,e'p) calculations:
	real*8 poftheta		!p as calculated from theta, assuming elastic.
	real*8 corrsing		!'corrected singles' for elastic
!	real*8 mm,mm2		!missing mass (assuming struck nucleon)
!	real*8 mmA,mmA2		!missing mass (assuming struck nucleus)
	real*8 Pm_Heepx,Pm_Heepy,Pm_Heepz	!Pm components for Heepcheck.
   	real*8 v1(4),v2(4),v3(4),ampi0/0.135/,amgam/0./,costh,phi
	real*8 x1,x2,y1,y2,dy,npsxmin,npsymin,npsemin,mchk
	real*8 grnd		!random # generator.

!local (e,e'pi/K) calculations:
!	real*8 t		!t
	real*8 dummy

	logical success

	if (debug(2)) write(6,*)'r_ntu_write: entering'

!If event did not make it thru spectrometers, return.  Later, we will want to
! add option to write event even if it failed when doing_phsp.

	if (.not.success) return

c pyb don't write zero weight events (for example from
c excl pi- on p or low energy rhos)
        if(main%weight.eq.0.) return

C Pyb added this section for NPS experiment.
c if doing pi0 in final state, have it decay into two photons
c with 98.55 b.r. 
	if(which_pion.eq.2.or.which_pion.eq.3) then
	 v1(1) = vertex%p%E / 1000.
	 v1(2) = vertex%p%P * vertex%up%x / 1000.
	 v1(3) = vertex%p%P * vertex%up%y / 1000.
	 v1(4) = vertex%p%P * 
     >    sqrt(1. - vertex%up%x**2 - vertex%up%y**2)  / 1000. 
	 costh = -1. + 2.0 * grnd()
	 phi = 2. * pi * *grnd()
	 mchk = sqrt(v1(1)**2 - v1(2)**2 - v1(3)**2 - v1(4)**2)
	 if(mchk-ampi0.gt.0.001) write(6,'(''mchk'',2f8.4)')
     >    mchk,ampi0
         call  DECAY (V1,V2,V3,ampi0,amgam,amgam,COSTH,PHI)
c	 write(6,'(/6f7.3)') v1,costh,phi
	 x1 = v2(2)/v2(1) * 300.
	 dy = atan(v2(3)/v2(4)) - spec%p%theta
	 y1 = dy * 300.
c 	 write(6,'(6f7.3)') v2,x1,y1
	 x2 = v3(2)/v3(1) * 300.
	 dy = atan(v3(3)/v3(4)) - spec%p%theta
	 y2 = dy * 300.
c         write(6,'(6f7.3)') v3,x2,y2
c would at least one photon hit the NPS at 300 cm and
c have an energy above 150 MeV?
	 npsxmin = 40.
	 npsymin = 40.
	 npsemin = 0.15
	 if((abs(x1).lt.npsxmin .and.
     >       abs(y1).lt.npsymin .and.
     >       v2(1).gt.npsemin) .or. 
     >      (abs(x2).lt.npsxmin .and.
     >       abs(y2).lt.npsymin .and.
     >      v3(1).gt.npsemin)) then
	  write(88,288) 
     >    main%weight,
     >    main%sigcc,				!d5sig
     >    ntup%sigcm,                          !pion sig_cm
     >    min(999.,max(-999.,recon%e%delta)), 
     >    min(999.,max(-999.,recon%e%yptar*100.)), 
     >    min(999.,max(-999.,recon%e%xptar*100.)), 
     >    min(999.,max(-999.,recon%e%z)),
     >    min(999.,max(-999.,main%FP%e%x)), 
     >    min(999.,max(-999.,main%FP%e%dx*100.)), 
     >    min(999.,max(-999.,main%FP%e%y)), 
     >    min(999.,max(-999.,main%FP%e%dy*100.)),
     >    min(999.,max(-999.,v2(2))), 
     >    min(999.,max(-999.,v2(3))), 
     >    min(999.,max(-999.,v2(4))), 
     >    min(999.,max(-999.,v3(2))), 
     >    min(999.,max(-999.,v3(3))), 
     >    min(999.,max(-999.,v3(4))), 
     >    min(999.,max(-999.,-main%target%z*spec%p%sin_th)),
     >    min(999.,max(-999.,-main%target%rastery)),
     >	  min(999.,max(-999.,main%thetacm)),
     >    min(999.,max(-999.,main%phicm)), 
     >    min(999.,max(-999.,main%t/1.e6)),
     >    min(999.,max(-999.,orig%e%delta)),
     >    min(999.,max(-999.,orig%e%xptar*100.)),			!mr
     >    min(999.,max(-999.,orig%e%yptar*100.)),			!mr
     >    min(999.,max(-999.,-main%target%z*spec%e%sin_th)),
     >    min(999.,max(-999.,recon%theta_pq)),
     >    min(999.,max(-999.,recon%phi_pq)),
     >    min(999.,max(-999.,main%phi_pq))

 288	  format(3e12.4,28f9.4)
	 endif
	 return
	endif

c pyb added this section
c If doing_dvcs, just write out the single photon
c xxx not done yet!

*	if (phot1.gt.lim1) write(6,*) 'phot1,lim1=',phot1,lim1
*	if (phot2.gt.lim2) then
*	   write(6,*) 'phot2,lim2=',phot2,lim2
*	   write(6,*) egamma_tot_max-egamma_used(1)
*	   write(6,*) vertex.e.E - edge.e.E.min
*	endif
*	if (phot3.gt.lim3) write(6,*) 'phot3,lim3=',phot3,lim3


!Calculate some proton/pion/kaon specific stuff:
!
!	if (doing_pion .or. doing_kaon) then
!	  mm2 = recon.Em**2 - recon.Pm**2
!	  mm  = sqrt(abs(mm2)) * abs(mm2)/mm2
!	  mmA2= (recon.nu + targ.M - recon.p.E)**2 - recon.Pm**2
!	  mmA = sqrt(abs(mmA2)) * abs(mmA2)/mmA2
!	  t = recon.Q2 - Mh2
!     >	    + 2*(recon.nu*recon.p.E - recon.p.P*recon.q*cos(recon.theta_pq))
!	endif

	if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  poftheta = Mp*Ebeam / (2*ebeam*sin(recon%e%theta/2.)**2 + Mp)
	  corrsing = recon%e%P - poftheta
	  Pm_Heepz = -(recon%Pmy*recon%uq%y+recon%Pmz*recon%uq%z)
     >		/ sqrt(recon%uq%y**2+recon%uq%z**2)
	  Pm_Heepy =  (recon%Pmz*recon%uq%y-recon%Pmy*recon%uq%z)
     >		/ sqrt(recon%uq%y**2+recon%uq%z**2)
	  Pm_Heepx =  -recon%Pmx
	endif

	if(electron_arm.eq.1 .or. electron_arm.eq.3.or. electron_arm.eq.7)then !electron = right side.
	  ntu(1) = recon%e%delta
	  ntu(2) = recon%e%yptar			!mr
	  ntu(3) = recon%e%xptar			!mr
	  ntu(4) = recon%e%z
	  ntu(5) = main%FP%e%x
	  ntu(6) = main%FP%e%dx				!mr
	  ntu(7) = main%FP%e%y
	  ntu(8) = main%FP%e%dy				!mr
	  ntu(9) = orig%e%delta
	  ntu(10) = orig%e%yptar			!mr
	  ntu(11) = orig%e%xptar			!mr
	  ntu(12) = main%target%z*spec%e%sin_th
	  ntu(13) = recon%p%delta
	  ntu(14) = recon%p%yptar			!mr
	  ntu(15) = recon%p%xptar			!mr
 	  ntu(16) = recon%p%z
 	  ntu(17) = main%FP%p%x
	  ntu(18) = main%FP%p%dx			!mr
 	  ntu(19) = main%FP%p%y
 	  ntu(20) = main%FP%p%dy			!mr
	  ntu(21) = orig%p%delta
	  ntu(22) = orig%p%yptar			!mr
	  ntu(23) = orig%p%xptar			!mr
 	  ntu(24) = -main%target%z*spec%p%sin_th
	else if (electron_arm.eq.2 .or. electron_arm.eq.4 .or.
     >		 electron_arm.eq.5 .or. electron_arm.eq.6.or. electron_arm.eq.8) then  !e- = left.
	  ntu(1) = recon%p%delta
	  ntu(2) = recon%p%yptar			!mr
	  ntu(3) = recon%p%xptar			!mr
	  ntu(4) = recon%p%z
	  ntu(5) = main%FP%p%x
	  ntu(6) = main%FP%p%dx				!mr
	  ntu(7) = main%FP%p%y
	  ntu(8) = main%FP%p%dy				!mr
c	  ntu(9) = vertex%p%delta
c	  ntu(10) = vertex%p%yptar			!mr
c	  ntu(11) = vertex%p%xptar			!mr
	  ntu(9) = orig%p%delta
	  ntu(10) = orig%p%yptar			!mr
	  ntu(11) = orig%p%xptar			!mr
	  ntu(12) = main%target%z*spec%p%sin_th
	  ntu(13) = recon%e%delta
	  ntu(14) = recon%e%yptar			!mr
	  ntu(15) = recon%e%xptar			!mr
	  ntu(16) = recon%e%z
 	  ntu(17) = main%FP%e%x
	  ntu(18) = main%FP%e%dx			!mr
 	  ntu(19) = main%FP%e%y
 	  ntu(20) = main%FP%e%dy			!mr
	  ntu(21) = orig%e%delta
	  ntu(22) = orig%e%yptar			!mr
	  ntu(23) = orig%e%xptar			!mr
	  ntu(24) = -main%target%z*spec%e%sin_th
	else
	  write (6,*) 'results_write not yet set up for your spectrometers.'
	endif
	ntu(25) = recon%q/1000.				!q - GeV/c	
	ntu(26) = recon%nu/1000.			!nu - GeV
	ntu(27) = recon%Q2/1.e6				!Q^2 - (GeV/c)^2
	ntu(28) = recon%W/1000.				!W - GeV/c
	ntu(29) = recon%epsilon				!epsilon
	ntu(30) = recon%Em/1000.			!GeV
	ntu(31) = recon%Pm/1000.			!GeV/c
	ntu(32) = recon%theta_pq			!theta_pq - radians
	ntu(33) = recon%phi_pq				!phi_pq - radians

	if (doing_pion .or. doing_kaon .or. doing_delta) then
	  ntu(34) = ntup%mm/1000.			!missmass (nucleon)
	  ntu(35) = ntup%mmA/1000.			!missmass (nucleus)
	  ntu(36) = recon%p%P/1000.			!ppi - GeV/c
	  ntu(37) = ntup%t/1.e6				!t - GeV^2
	  ntu(38) = recon%PmPar/1000.
	  ntu(39) = recon%PmPer/1000.
	  ntu(40) = recon%PmOop/1000.
	  ntu(41) = -main%target%rastery		!fry - cm
	  ntu(42) = ntup%radphot/1000.			!radphot - GeV
	  dummy = pferx*vertex%uq%x + pfery*vertex%uq%y + pferz*vertex%uq%z
	  if (dummy.eq.0) dummy=1.e-20
	  ntu(43) = pfer/1000.*abs(dummy)/dummy		!p_fermi - GeV/c
	  ntu(44) = main%sigcc				!d5sig
	  ntu(45) = ntup%sigcm				!pion sig_cm
	  ntu(46) = main%weight
	  ntu(47) = decdist				!decay distance (cm)
	  ntu(48) = sqrt(Mh2_final)
	  ntu(49) = pfer/1000.*dummy			!p_fermi along q.
	  ntu(50) = vertex%Q2/1.e6
	  ntu(51) = main%w/1.e3
	  ntu(52) = main%t/1.e6
	  ntu(53) = main%phi_pq
	  if(using_tgt_field) then
	     ntu(54) = recon%theta_tarq
	     ntu(55) = recon%phi_targ
	     ntu(56) = recon%beta
	     ntu(57) = recon%phi_s
	     ntu(58) = recon%phi_c
	     ntu(59) = main%beta
	     ntu(60) = vertex%phi_s
	     ntu(61) = vertex%phi_c	     
	     if (doing_kaon) then
		ntu(62) = ntup%sigcm1 !sigcm - saghai model
		ntu(63) = ntup%sigcm2 !sigcm - factorized.
	     endif
	  else
	     if (doing_kaon) then
		ntu(54) = ntup%sigcm1 !sigcm - saghai model
		ntu(55) = ntup%sigcm2 !sigcm - factorized.
	     endif
	  endif
	else if (doing_semi.or.doing_rho) then
	  ntu(34) = ntup%mm/1000.			!missmass (nucleon)
	  ntu(35) = recon%p%P/1000.			!ppi - GeV/c
	  ntu(36) = ntup%t/1.e6				!t - GeV^2
	  ntu(37) = -main%target%rastery		!fry - cm
	  ntu(38) = ntup%radphot/1000.			!radphot - GeV
	  ntu(39) = main%sigcc				!d5sig
	  ntu(40) = main%sigcent			!d5sig - central kin.
	  ntu(41) = main%weight
	  ntu(42) = decdist				!decay distance (cm)
	  ntu(43) = sqrt(Mh2_final)
	  ntu(44) = recon%zhad
	  ntu(45) = vertex%zhad
	  ntu(46) = recon%pt2/1.e06
	  ntu(47) = vertex%pt2/1.e06
	  ntu(48) = recon%xbj
	  ntu(49) = vertex%xbj
	  ntu(50) = acos(vertex%uq%z)
	  ntu(51) = ntup%sigcm
	  ntu(52) = main%davejac
	  ntu(53) = main%johnjac
	  dummy = pferx*vertex%uq%x + pfery*vertex%uq%y + pferz*vertex%uq%z
	  ntu(54) = pfer/1000.*abs(dummy)/dummy		!p_fermi - GeV/c
	  ntu(55) = ntup%xfermi
	  ntu(56) = main%phi_pq
	  if(using_tgt_field) then
	     ntu(57) = recon%theta_tarq
	     ntu(58) = recon%phi_targ
	     ntu(59) = recon%beta
	     ntu(60) = recon%phi_s
	     ntu(61) = recon%phi_c
	     ntu(62) = main%beta
	     ntu(63) = vertex%phi_s
	     ntu(64) = vertex%phi_c	     
	     if(doing_rho) then
		ntu(65) = ntup%rhomass
		ntu(66) = ntup%rhotheta
	     endif
	  else
	     if(doing_rho) then
		ntu(57) = ntup%rhomass
		ntu(58) = ntup%rhotheta
	     endif
	  endif
	else if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  ntu(34) = corrsing/1000.
	  ntu(35) = Pm_Heepx/1000.
	  ntu(36) = Pm_Heepy/1000.
	  ntu(37) = Pm_Heepz/1000.
	  ntu(38) = recon%PmPar/1000.
	  ntu(39) = recon%PmPer/1000.
	  ntu(40) = recon%PmOop/1000.
	  ntu(41) = -main%target%rastery		!fry - cm
	  ntu(42) = ntup%radphot/1000.			!radphot - GeV
	  ntu(43) = main%sigcc
	  ntu(44) = main%weight
	endif


c	call HFN(NtupleID,ntu)
	do i=1,NtupleSize
	   write(NtupleIO) ntu(i)
	enddo
	if (debug(2)) write(6,*)'r_ntu_write: ending'
	return
	end


	subroutine results_ntu_write1(vertex,recon,main,success)

	implicit none
	include 'hbook.inc'
	include 'simulate.inc'

	integer*4 nentries
	parameter (nentries = 9)

	real*8	ntu(nentries)
	logical success
	type(event_main):: main
	type(event):: vertex, recon

	if (debug(2)) write(6,*)'r_ntu_write: entering'
	ntu(1) = vertex%e%delta
	ntu(2) = vertex%e%yptar
	ntu(3) = -vertex%e%xptar
	ntu(4) = main%SP%e%z
	if(success)then
	  ntu(5) = recon%e%delta
	  ntu(6) = recon%e%yptar
	  ntu(7) = recon%e%xptar
	  ntu(8) = recon%e%z
	else
	  ntu(5) = 30.
	  ntu(6) = 0.1
	  ntu(7) = 0.1
	  ntu(8) = 4.
	endif
	ntu(9) = main%weight
c	call HFN(NtupleID,ntu)
	if (debug(2)) write(6,*)'r_ntu_write: ending'
	return
	end
C SUBROUTINE TO DECAY A PARTICLE INTO TWO OTHER PARTICLES
C AT angle COSTH in THE C.M. FRAME WHERE THE DIRECTION OF THE
C INITIAL PARTICLE NEED NOT BE ALONG THE Z-AXIS. PHI SHOULD BE 0 TO 2 PI
C vectors are (E, px, py, pz)
      SUBROUTINE DECAY (V1,V2,V3,M1,M2,M3,COSTH,PHI)
      IMPLICIT NONE
      REAL*8 V1(4),V2(4),V3(4),KV(4),M1,M2,M3,V1P,G,K,COSTH,SINTH,PHI
      REAL*8 CHK,CHK1,CHK2,KCHK,BETA,VT(4),VTT(4),VDM
! Find magnitude of momentum of initial particle
      V1P=SQRT(V1(2)**2+V1(3)**2+V1(4)**2)
      G=V1(1)/M1
      BETA=V1P/V1(1)

! Find c.m. momentum and pick decay isotropically
      K=SQRT(( ((M1**2-M2**2-M3**2)/2.)**2 - (M2*M3)**2 )/M1**2)
      KCHK=SQRT((M1**2-(M2+M3)**2)*(M1**2-(M2-M3)**2))/2./M1
      IF(ABS(K-KCHK).GT.0.001) WRITE(6,'(1X,''ERROR,K,KCHK='',2F10.4)')
     >  K,KCHK
      IF(ABS(SQRT(K*K+M2*M2)+SQRT(K*K+M3*M3)-M1).GT.0.001)
     >  WRITE(6,'(1X,''ERROR IN CM ENERGY'')')
      SINTH=SQRT(1.-COSTH**2)
! Put k into vector for first decay particle
      KV(1)=SQRT(K*K+M2*M2)
      KV(2)=K*SINTH*COS(PHI)
      KV(3)=K*SINTH*SIN(PHI)
      KV(4)=K*COSTH
! Transform k to lab frame for first decay particle
! First find vector as though decay were along Z axis
      V2(1)=G*(KV(1)+BETA*KV(4))
      VT(2)=KV(2)
      VT(3)=KV(3)
      VT(4)=G*(BETA*KV(1)+KV(4))
! Now rotate vector in xz plane
      VDM=SQRT(V1(2)**2+V1(4)**2)
      VTT(2)= VT(2)*V1(4)/VDM+VT(4)*V1(2)/VDM
      VTT(3)= VT(3)
      VTT(4)=-VT(2)*V1(2)/VDM+VT(4)*V1(4)/VDM
! Now rotate vector in yz plane
      VDM=SQRT(V1(3)**2+V1(4)**2)
      V2(2)= VTT(2)
      V2(3)= VTT(3)*V1(4)/VDM + VTT(4)*V1(3)/VDM
      V2(4)=-VTT(3)*V1(3)/VDM + VTT(4)*V1(4)/VDM
! Transform k to lab frame for first  decay particle
      CHK=V2(1)**2-V2(2)**2-V2(3)**2-V2(4)**2-M2**2
      IF(ABS(CHK).GT.0.005) WRITE(6,'(1X,''ERROR, CHK='',E9.2)') CHK
! Put k into vector for second decay particle
      KV(1)=SQRT(K*K+M3*M3)
      KV(2)=-K*SINTH*COS(PHI)
      KV(3)=-K*SINTH*SIN(PHI)
      KV(4)=-K*COSTH
! Transform k to lab frame for second decay particle
! First find vector as though decay were along Z axis
      V3(1)=G*(KV(1)+BETA*KV(4))
      VT(2)=KV(2)
      VT(3)=KV(3)
      VT(4)=G*(BETA*KV(1)+KV(4))
! Now rotate vector in xz plane
      VDM=SQRT(V1(2)**2+V1(4)**2)
      VTT(2)= VT(2)*V1(4)/VDM+VT(4)*V1(2)/VDM
      VTT(3)= VT(3)
      VTT(4)=-VT(2)*V1(2)/VDM+VT(4)*V1(4)/VDM
! Now rotate vector in yz plane
      VDM=SQRT(V1(3)**2+V1(4)**2)
      V3(2)= VTT(2)
      V3(3)= VTT(3)*V1(4)/VDM + VTT(4)*V1(3)/VDM
      V3(4)=-VTT(3)*V1(3)/VDM + VTT(4)*V1(4)/VDM
      CHK=V3(1)**2-V3(2)**2-V3(3)**2-V3(4)**2-M3**2
      IF(ABS(CHK).GT.0.005) WRITE(6,'(1X,''ERROR, CHK='',E9.2)') CHK
      CHK1=V1(1)-V2(1)-V3(1)
      CHK2=V1(2)-V2(2)-V3(2)+V1(3)-V2(3)-V3(3)+V1(4)-V2(4)-V3(4)
      IF(ABS(CHK1).GT.0.005.OR.ABS(CHK2).GT.0.05) WRITE(6,
     >  '(1X,''E12='',14F5.2)') CHK1,CHK2,V1,V2,V3
      RETURN
      END

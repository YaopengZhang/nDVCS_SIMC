      program makeinfile_pi0
c this is for the 2023-24 NPS runs
      implicit none
      integer i,j,k,it,ikin,itarg,iprocess
      integer ngen
      real*8 ebeam,phms,thhms,pipeak
      real*8 pnps,thnps
      real*8 hmp_x, hmp_y, pmp_x, pmp_y,z,zp
      real*8 tgt_A,tgt_z,tgt_amu,tgt_rho,tgt_t

! kinematic settings
      integer nkin/28/
      real*8 e0sv(28)/10.54, 10.54, 10.54, 10.54, 10.54, 
     >   10.54,10.54,10.54,
     >   8.457,8.457,8.457,8.457,8.457,8.457,8.457,
     >  10.538, 10.54, 10.54, 10.54, 10.54, 10.54,
     >  10.54, 10.54, 10.54, 10.54, 10.54, 10.54, 10.54/
      real*8 ep0sv(28)/6.67, 5.89, 6.11, 4.64, 4.64, 
     >   5.233,5.233,5.878,
     >   4.042,4.042,4.042,4.726,4.726,3.803,3.803,
     >   4.637, 6.667, 6.667,5.253, 5.253,4.637,
     >   5.878,5.878,5.038,5.038, 2.416, 2.416, 4.149 /
      real*8 the0sv(28)/12.49, 16.48, 12.37, 16.44, 16.433, 
     >   16.91,16.91,16.483,
     >   17.010,17.010,17.010,16.75,16.75,22.92,22.92,
     >    16.440, 12.48, 12.48, 16.93, 16.93, 16.44,
     >    16.48, 16.48, 19.35, 19.35, 26.800, 26.800, 15.20/
      real*8 thp0sv(28)/20.58, 18.72, 15.96, 12.30, 14.39, 
     >   15.45,14.39,17.720,
     >   12.460,14.360,16.600,17.080,18.980,16.58,12.45,
     >   12.200, 19.15, 21.70,13.43, 16.88, 14.00,
     >   16.72, 20.15, 14.08, 17.51, 7.500, 7.500, 11.20/
      character*3 tname(3)/'LH2','LD2','DUM'/
      character*5 pname(4)/' dvcs',' excl','delta','sidis'/
      character*80 fname,fname2

c make the json file. Each job does all targets and processes for
c a given kinematic setting

      open(unit=50,file='simcnps.json')
       write(50,198,advance="no")
 198   format('{"name": "simcnps", "jobs": [')

c loop over kinematics
      do ikin = 1, nkin
        write(crun,'(i2.2)') ikin
        fname='kin'//trim(crun)
        open(unit=51,file=fname)
        write(50,199,advance="no") ikin,ikin
 199    format('{"os": "general", "partition": "production", ', 
     >  '"command": ["source /u/group/nps/bosted/',
     >  'simc_new/scripts/kin',i2.2,
     >  '"], "diskBytes": 2000000000, "ramBytes": 1000000000,', 
     >  '"cpu_Cores": 1, "account":', 
     >  '"hallc", "name": "simc',i2.2,'", "shell": "/usr/',
     >  'bin/tcsh", "timeSecs": 28800}')
        write(51,1223) fname2
 1223   format(
     >  '#!/bin/tcsh'/
     >  'cd /u/group/c-sidis/bosted/simc'/
c loop over targets
       do itarg = 1,3
c loop over processes
        do iprocess=1,4

c open the input file
         write(ckin,'(i2.2)') ikin
         write(fname2,123) ikin,tname(itarg),pname(iprocess)
 123     format('simc_',i2.2,'_',a,'_',a,'.inp')
         fname = 'infiles'//fname2
         open(unit=15,file=fname)

c add to batch script
         write(51,1223) fname2
 1223    format(
     >    'setenv SIMCIN ',a/
     >    './simc')

c variables for each kinematic setting
         ebeam = e0sv(ikin)
         phms  = ep0sv(ikin)
         pnps = ebeam - phms
         thhms = the0sv(ikin)
         thnps = thp0sv(ikin)

c get mis-pointing of HMS spectrometer
         if(thhms.lt.40.) then
          hmp_y = 0.1 * (0.52 - 0.012*thhms + 0.002 * thhms**2)
         else 
          hmp_y = 0.1 * (0.52 - 0.012 * 40. + 0.002 * 40.**2)
         endif
         pmp_y = 0.1 * (-0.6)
         if (thhms < 50) then
          hmp_x = 0.1 * (2.37 - 0.086 * thhms + 0.0012 * thhms**2)
         else
          hmp_x = 0.1 * (2.37 - 0.086 * 50.   + 0.0012 * 50.**2)
         endif
         pmp_x = 0.1 * (-1.26);

c target vairables, A, Z, dentisy and tgt_t is density times lenth
         if(itarg.eq.1) then
          tgt_A = 1.
          tgt_z = 1.
          tgt_amu = 1.007276
          tgt_rho = 0.07231
          tgt_t = 717.8 ! NEED TO UPDATE!
         endif
         if(itARG.eq.2) then
          tgt_A = 2.
          tgt_z = 1.
          tgt_amu = 2.01355
          tgt_rho = 0.167
          tgt_t = 1662.3 ! NEED TO UPDATE
         endif
c Dummy. NEED TO UPDATE!
         if(itarg.eq.3) then
          tgt_A = 27.
          tgt_z = 13.
          tgt_amu = 26.98154
          tgt_rho = 2.700
          tgt_t = 362.000
          tgt_t = 362.000 / 2.
         endif

c start writing to the input file here
         write(15,'("; This is a CTP file")')
         write(15,'(" ")')
         write(15,'("begin parm experiment")')
c how many events to generate. Making it smallish for now
         if(iprocess.eq.1) ngen = 10000
         if(iprocess.eq.2) ngen = 10000
         if(iprocess.eq.3) ngen =  2000
         if(iprocess.eq.4) ngen =  20000
         write(15,3222) ngen * zfact
 3222    format('  ngen = ',i7,'	        ;  POS: # of successe')

         write(15,'("EXPER%charge = 1.0      ;  total charge (mC)")' )
         write(15,'("doing_phsp = 0		;  (ONE = TRUE)")' )
         write(15,'("doing_kaon = 0		;  (ONE = TRUE) ")' )

         if(iprocess.gt.1) then
           write(15,'("doing_pion = 1		;  (ONE = TRUE)" )' )
         else
           write(15,'("doing_pion = 0		;  (ONE = TRUE)" )' )
         endif
      
         write(15,'("which_pion = 2          ;  2 = pi0'')')

         if(iprocess.ne.3) then
          write(15,'("doing_delta = 0         ; H(e,ep)pi0 ")' ) 
         else
          write(15,'("doing_delta = 1         ; H(e,ep)pi0 ")' ) 
         endif

         write(15,'("doing_rho = 0           ; exclusive rho ")')

         if(iprocess.eq.4) then
          write(15,'("doing_semi = 1          ; doing semi-inclusive? Need 
     >    to set doing_pion or doing_kaon ")') 
         else
          write(15,'("doing_semi = 0          ; doing semi-inclusive? Need 
     >    to set doing_pion or doing_kaon ")') 
         endif

         write(15,'("doing_hplus = 0	        ; positive hadrons? (only 
     >    for semi or rho)" )' )

c for NPS, this means decay final state pi0 into 2 gammas or no
         if(iprocess.gt.1) then
          write(15,'("doing_decay = 1		;  1=decay ON, 0=decay OFF ")')
         else
           write(15,'("doing_decay = 0		;  1=decay ON, 0=decay OFF ")')
         endif

c this is not used in NPS
         write(15,'("ctau = 780.4	        ;  decay length (cm) ")' )


         write(15,'("end parm experiment")' )
         write(15,'(" ")')

         write(15,'("begin parm kinematics_main")' )
         write(15,'("Ebeam = ",f9.2)') ebeam*1000.0
         write(15,'("dEbeam = 0.05	         ;  beam energy variation 
     >  (%)")' )

        write(15,'("electron_arm = 1         ; 1=hms,2=sos,3=hrsr,4=hrsl,
     >  5=shms,7=calo (7=calo HMS side),")')

        write(15,'("hadron_arm = 8           ; 1=hms,2=sos,3=hrsr,4=hrsl,
     >  5=shms,7=calo (8=calo (SOS side)")' )

c this isn't used in NPS
        write(15,'("use_first_cer = 0        ; Use first cerenkov 
     >  in shms?")') 

        write(15,'("spec%e%P = ",f9.2)') abs(phms)*1000.0
c		 ;  e arm central momentum (MeV/c)" )' )
        write(15,'("spec%e%theta = ",f6.3)') thhms		
c                ;  e arm angle setting (degrees) ")' )
        write(15,'("spec%p%P = ",f9.2)') abs(pnps)*1000.0		
c                 ;  p arm central momentum (MeV/c) ")' )
        write(15,'("spec%p%theta = ",f6.3)') thnps		
c                ;  p arm angle setting (degrees)

        write(15,'("end parm kinematics_main")' )
        write(15,'("  ")' )

        write(15,'("begin parm target")' )
        write(15,'("targ%A = ",f4.1)') tgt_A			
c                 ;  target A
        write(15,'("targ%Z = ",f4.1)') tgt_z			
c                 ;  target Z
        write(15,'("targ%mass_amu = ",f9.6)') tgt_amu
c                 ;  target mass in amu
        write(15,'("targ%mrec_amu = 0.0      ;  recoil mass in amu (eep=A-1 
     >  system,pion=A-2)")' )
        write(15,'("targ%rho = ",f6.4)' ) tgt_rho 
c                 ;  target density (g/cm^3)
        write(15,'("targ%thick = ",f8.3)' ) tgt_t
c 		 ;  target thick (mg/cm^2) for 10 cm target
        write(15,'("targ%angle = 0.	         ; target angle (for solid 
     >  target) (degrees)")' )
        write(15,'("targ%abundancy = 100.0	 ;  target purity (%)")' )
        write(15,'("targ%can = 1		 ;  1=beer can (fpi), 
     >    2=pudding can (nucpi)" )' )
        write(15,'("end parm target")' )
        write(15,'("  ")' )

        write(15,'("end parm target")' )
        write(15,'("  ")' )

        write(15,'("begin parm debug         ; (ONES give helpful 
     >  debug info)")')
        write(15,'("debug(1) = 0	         ; turns on output from 
     >  brem.f")')
        write(15,'("debug(2) = 0		 ;  into/outa subs.")' )
        write(15,'("debug(3) = 0		 ;  spit out values (init. 
     >  and main loop).")' )
        write(15,'("debug(4) = 0		 ;  mostly comp_ev, gen_rad 
     >  diagnostics.")' )
        write(15,'("debug(5) = 0		 ;  a bit of everything.")' )
        write(15,'("end parm debug")' )
        write(15,'(" ")' )
        write(15,'(" ")' )

c these ranges are all a bit bigger than we will actually use
        write(15,'("begin parm e_arm_accept")' )
        write(15,'("SPedge%e%delta%min = -14.0   ;  delta min (SPECTROMETER 
     >  ACCEPTANCE !)")' )
        write(15,'("SPedge%e%delta%max =  14.0   ;  delta max")' )
        write(15,'("SPedge%e%yptar%min = -32.0   ; .yptar.min = 
     >  {TF} / 1000 (mrad)" )' )
        write(15,'("SPedge%e%yptar%max =  32.0   ; .yptar.max = 
     >  {TF} / 1000")')
        write(15,'("SPedge%e%xptar%min = -70.0  ; .xptar.min = 
     >  {TF} / 1000 (mrad)" )' )
        write(15,'("SPedge%e%xptar%max =  70.0  ; .xptar.max = 
     >  {TF} / 1000" )' )
        write(15,'("end parm e_arm_accept")' )
        write(15,'(" ")' )

c these params are for the NPS calorimeter
        write(15,'("begin parm p_arm_accept")' )
c set this to go down to 1.5 GeV pi0's or single photons
        delmin = -100. * (1 - 1.5 / pnps)
        write(15,'("SPedge%p%delta%min = '',f72,'' ; delta min 
     >  ACCEPTANCE !)")' ) delmin
        write(15,'("SPedge%p%delta%max =  5.    ;  delta max")' )
        write(15,'("SPedge%p%yptar%min = -130.0   ; .yptar.min = 
     >  {TF} / 1000 (mrad)" )' )
        write(15,'("SPedge%p%yptar%max =  130.0   ; .yptar.max = 
     >  {TF} / 1000" )' )
        write(15,'("SPedge%p%xptar%min = -130.0  ; .xptar.min = 
     >  {TF} / 1000 (mrad)" )' )
        write(15,'("SPedge%p%xptar%max =  130.0  ; .xptar.max = 
     > {TF} / 1000")' )
       write(15,'("end parm p_arm_accept")' )
       write(15,'(" ")')

c changed from 1 to 3, but not sure why!
       write(15,'("targ%fr_pattern = 3.	     ;  raster pattern: 1=square, 
     > 2=circular")' )
       write(15,'("targ%fr1 = 0.1	             ; horizontal size OR 
     > inner radius(2)" )' )
       write(15,'("targ%fr2 = 0.1	             ; vertical size OR 
     > outer radius(2)")' )
       write(15,'("targ%xoffset = 0.0	     ;  target x-offset (cm): +x 
     > = beam right")' )
c set  beam position to zero
       write(15,'("targ%yoffset =  0.0	     ;  target y-offset (cm): +y 
     > = up")' )
c for lh2 and ld2
       if(itarg.le.2) then
        write(15,'("targ%zoffset = 0.0	     ;  target z-offset (cm): +z 
     >  = downstream")' )
       else
c for dummy
        write(15,'("targ%zoffset = 5.0	     ;  target z-offset (cm): +z 
     >  = downstream")' )
       endif
       write(15,'("end parm beamandtergetinfo")' )
       write(15,'(" ")')

       write(15,'(";These are offsets applied before the call to 
     > the single arm montecarlos.")' )
       write(15,'("begin parm spect_offset")' )
       write(15,
     >  '("spec%e%offset%x = ",f5.2,10x,"; x offset (cm)")') hmp_x
       write(15,
     >  '("spec%e%offset%y = ",f5.2,10x," ; y offset (cm)")') hmp_y
       write(15,'("spec%e%offset%z = 0.	       ; z offset (cm)")')
       write(15,'("spec%e%offset%xptar = 0.       ; xptar 
     > offset (mr) ! x(y)ptar is slope, so")' )
       write(15,'("spec%e%offset%yptar = 0.       ; yptar 
     > offset (mr) ! it is really unitless.")' )
c here is where you can put in mis-pointing of NPS, if needed
       write(15,'("spec%p%offset%x = 0.0 	       ; x offset (cm)")' )
       write(15,'("spec%p%offset%y = 0.0 	       ; y offset (cm)")' )
       write(15,'("spec%p%offset%z = 0.	       ; z offset (cm)")' )
       write(15,'("spec%p%offset%xptar = 0.0      ; xptar 
     > offset (mr)")')
       write(15,'("spec%p%offset%yptar = 0.0      ; yptar 
     > offset (mr)")' )
       write(15,'("end parm spect_offset")' )
       write(15,'(" ")' )

       write(15,'("begin parm simulate")' )
       write(15,'("hard_cuts = 0               ;  (ONE = TRUE) SPedge 
     > and Em.max are hard cuts(ntuple)")' )
c alwasys do radiation
       write(15,'("using_rad = 1		    ;  (ONE = TRUE)")' )
       write(15,'("use_expon = 0		    ;  (LEAVE AT 0)")' )
       write(15,'("one_tail = 0		    ;  0=all, 1=e, 2=epr, 3=p, 
     > -3=all but p")' )
       write(15,'("intcor_mode = 1	            ;  (LEAVE AT 1)")' )
c NEED TO CHECK THIS1
       write(15,'("spect_mode = -2	            ;  0=e+p arms, -1=p arm, 
     > -2=e arm only, 1=none")' )
       write(15,'("cuts%Em%min = 0.	    ;  (Em.min=Em.max=0.0 gives 
     > wide open cuts)")' )
       write(15,'("cuts%Em%max = 0.	    ;  Must be wider than cuts 
     > in analysis(elast. or e,ep)" )' )
       write(15,'("using_Eloss = 1	            ;  (ONE = TRUE)")' )
       write(15,'("correct_Eloss = 1	    ;  ONE = correct reconstructed 
     > events for eloss.")' )
       write(15,'("correct_raster = 1	    ;  ONE=Reconstruct events 
     > using raster matrix elements")' )
       write(15,'("mc_smear = 1.	            ;  ONE = target & 
     > hut mult scatt AND DC smearing.")' )
       write(15,'("deForest_flag = 0	    ;  0=sigcc1, 1=sigcc2, 
     > -1=sigcc1 ONSHELL")' )
       write(15,'("rad_flag = 0		    ;  (radiative option #1...
     > see init.f)")' )
       write(15,'("extrad_flag = 2	            ;  (rad. option 
     > #2...see init.f)")' )
       write(15,'("lambda(1) = 0.0	            ;  if rad_flag.eq.4 
     > then lambda(1) = {TF}")' )
       write(15,'("lambda(2) = 0.0	            ;  if rad_flag.eq.4 
     > then lambda(2) = {TF}")' )
       write(15,'("lambda(3) = 0.0	            ;  if rad_flag.eq.4 
     > then lambda(3) = {TF}")' )
       write(15,'("Nntu = 1		    ;  ONE = generate ntuples")' )
       write(15,'("using_Coulomb = 0	    ;  (ONE = TRUE)")' )
       write(15,'("dE_edge_test = 0.	    ; (move around energy 
     > edges)")')
       write(15,'("use_offshell_rad = 1	    ;  (ONE = TRUE)")')
       egammax = (ebeam - 0.85*abs(phms)) * 1000. 
       write(15,'("Egamma_gen_max = ",f7.0,"         ;  Set >0 to hardwire the 
     >   Egamma limits. set big enough!")' ) egammax
       write(15,'("do_fermi = 0")' )
c used in SIDIS model sometimes
       write(15,'("pt_b_param = 3.8      ; was 4.9, now 3.8 
     > data for H and D +-")' )
       write(15,'("sigc_flag = 1         ; 0 = bin in z, 1 = bin in pt2, 
     > -1 = bin in xbj" )' )
c can set these to make central bin xsection tables
       write(15,'("sigc_nbin = 100.       ; number of bins for central 
     > cross section calc" )' )
       write(15,'("sigc_kin_min = 0.00   ; minumum z (or pt2) for central 
     > cross section calc")' ) 
       write(15,'("sigc_kin_max = 1.00   ; maximum z (or pt2) for central 
     > cross section calc")' )
       write(15,'("sigc_kin_ind = 0.55   ; value for independent variable 
     > (z or pt2 in GeV2)")' ) 
       rseed = ikin + 100*itarg + 10000*iprocess + int(100000000.*phms) 
       write(15,'("random_seed = ",i10,"   ; randm seed)")' ) rseed
       write(15,'(" ")' )

       write(15,'("end parm simulate")' )

       close(15)

       enddo ! end of loop over process
       enddo ! end of loop over target
      enddo ! loop over ikin

c close json file
      write(50,197)
 197  format(']}')

      stop
      end

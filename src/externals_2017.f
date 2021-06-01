      program EXTERNALS
c thetop                                                                        
C Program to calculate the measured cross section using Tsai's formula. 
C We use equivalent radiator approximation to estimate the internal     
C bremstrahlung contribution. However the integrals over radiation      
C length initial and final energy losses are done in full glory without 
C any approximations except at the edges of kinematics where the        
C divergences creep in.  Code by S. Dasu, modified by LWW.              
c Modified so list of input/output files is always looked
c for in extern.inp pb                                                       
c Modifed to allow multiple files to be looked at for Jan 05
c Modified to have Fermi smeaing of inelastic. see SECNUCLW line 2544
c Modified to have correct limits of integration for inelastic
c for nuclear target: see CONTINUUM and FUNC4. PEB 10/05
C=======================================================================
      IMPLICIT NONE
      COMMON /SET/   E0SET,EPSET,THSET      
      REAL  E0SET,EPSET,THSET   
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      COMMON /CON/   PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL 
      COMMON /OMG/   DELTA,R
      REAL   DELTA,R                                        
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL   
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                   
      COMMON /ERR/   REPS1,AEPS1
      REAL  REPS1,AEPS1
      COMMON /TRL/   TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL  TTARG,TWALL,TBEAM,TSPEC
      INTEGER NSEG,JTARG,ITARG,ipsv,ithsv,iesv,iqsv,iwsv,izz
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL  
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET  
      REAL TTOT,TLEN
      INTEGER ISEG154  
      COMMON/E154/ISEG154,TTOT,TLEN                       
      REAL  SIGMA_INEL(2),DELTA_INEL(2),                       
     >      SIGMA_QELA(2),DELTA_QELA(2),                       
     >      SIGMA_QEPK(2),DELTA_QEPK(2),                       
     >      SIGMA_QE_PK_COS(2),DELTA_QE_PK_COS(2),  
     >      SIGMA_ELPK(2),DELTA_ELPK(2),                       
     >      TOTAL(2),ZERO(18),SIGMA_CORR,Q2SET,SIGMA_BORN,RAT_SMR,
     >      VTRUN,TTRUN,VTSTOP,TTSTOP,VTSTART,TTSTART,X0
      REAL  SIGMA_INEL_154,SIGMA_QELA_154,SIGMA_QE_PK_COS_154,
     > SIGMA_ELPK_154,SIGMA_BORNPLSQE,SIGMA_elastic
      real e0sv(100000),epsv(100000),thsv(100000),xsv(100000)
      real w2sv(100000),sbsv(100000),srsv(100000),srqesv(100000)
      real srinsv(100000),srelsv(100000), ddnom
      INTEGER I,II,ISEG_NUM,j,npts
      REAL           CONTINUUM,QETAIL,QEPEAK,EPEAK,ATAILFL1,ATAILFL_QE
      EXTERNAL       CONTINUUM,QETAIL,QEPEAK,EPEAK,ATAILFL1,ATAILFL_QE
      logical usegrd
      common/testing/prttst,usegrd
      logical prttst,doeg1b,good
      common/experiment/ doeg1b
      real*4  xval(37),tmptmp,psf1(5),psf2(5),dr,w2last,ymax,ymin
      real*8 q28,w8,w18,w28,rcsim
      COMMON/ioana/xval

      real*8 jane0(23)/1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,
     >   2.3, 2.3, 2.3, 2.3, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4,
     >  4.6, 4.6, 4.6, 4.6/
      real*8 janth(23)/10.8, 13., 16., 19., 22., 28., 45., 55., 70.,
     >   20., 30., 45., 60., 14., 20., 28., 36., 40., 45.,
     >   10.65, 16.0, 20.0, 25.0/
c      logical dojan05/.true./
      logical dojan05/.false./
c      logical dofromin/.true./
      logical dofromin/.false./
c      logical dosimona/.true./
      logical dosimona/.false./
c      logical donucr/.true./
      logical donucr/.false./
c Put "Eg1b" on first line of extern.inpe for speical eg1b in/out format
      integer nfile, ifile, IAjan, itrgsim
      real*8 sigexp(100000),sigexper(100000),fact,Afromin,Zfromin
      real*8 qv,ams,amf2,a,bb,c,disc,costh,ysv(100000),xplt(100000)
      real*8 ysvp(100000)
      real*8 thfromin(6)/18.,22.,26.,32.,40.,50./
      real*8 thsimona(6)/38.,41.,45.,55.,60.,70./
      character*30 fname,fnamep
! used to tell quasiy8 if really smearing elastic peak
      logical doing_elas
      common/doelas/ doing_elas

c set ture if doing many points. 
c ONLY USE if Jlab kinematics up to 6 gev
c      usegrd=.true.
      usegrd=.false.

c set true if want proton elastic to be smeared and treated as q.e.
c if so, will smear by amound smrdEp
cc code currently modified so only does this for A=1
cc need to change back if want it for A>1 or q.e., inel. smearing
cc      smrelasp = .true.
      smrelasp = .false.
cc see below for smrdEp since depnds on kinematics

c set true if want to use Arenhoevel for q.e., threshold
cc      usearen = .true.
      usearen = .false.

cc      DELTA = 0.010                                                 
! xxx try making much smaller to get A(q2), q.e. at low Q2
      DELTA = 0.002
c      write(6,'(''A b='',e11.3)') b

      REPS1 = 0.0002  !0.001  changed 21 Nov 00
      AEPS1 = 1.E-19                                                    
      PI    = 3.1415927                                                 
      PM    = 0.93828                                                   
      PP    = 0.13957                                                   
      EM    = 0.000511                                                  
      AL    = 7.29735E-03                                               
                                                                        
      call ARENHOVEL_INIT()

c      write(6,'(1x,''got here'')')
c      open(unit=17,file='output/externals/kpp_shms_488.out')
c      open(unit=17,file='output/externals/single_point_short.out')
c      open(unit=17,file='output/externals/single_point_pol_3He.out')
c      open(unit=17,file='output/externals/test_A1nd2n_hms.out')
c      open(unit=17,file='output/externals/test_A1nd2n_shms.out')
c      open(unit=17,file='output/externals/test_dflay_hms.out')
c      open(unit=17,file='output/externals/test_dflay_shms.out')
c      open(unit=17,file='output/externals/N2_ref_hms_4.out')
c      open(unit=17,file='output/externals/N2_ref_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_test_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_test_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_O_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_O_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_Si_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_Si_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_Ba_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_Ba_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_Al_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_Al_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_Ca_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_Ca_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_Sr_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_Sr_shms_B.out')
c      open(unit=17,file='output/externals/Empty_up_win_H_hms_4.out')
c      open(unit=17,file='output/externals/Empty_up_win_H_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_O_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_O_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_Si_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_Si_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_Ba_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_Ba_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_Al_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_Al_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_Ca_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_Ca_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_Sr_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_Sr_shms_B.out')
c      open(unit=17,file='output/externals/Empty_down_win_H_hms_4.out')
c      open(unit=17,file='output/externals/Empty_down_win_H_shms_B.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_20deg.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_20deg.out')
      open(unit=17,file='output/externals/d2n_3He_hms_kin_C.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_kin_Z.out')
c      open(unit=17,file='output/externals/N2_ref_hms_C.out')
c      open(unit=17,file='output/externals/N2_ref_shms_Z.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_kin_A.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_kin_X.out')
c      open(unit=17,file='output/externals/N2_ref_hms_A.out')
c      open(unit=17,file='output/externals/N2_ref_shms_X.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_1pass_lo.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_1pass_lo.out')
c      open(unit=17,file='output/externals/N2_ref_hms_1pass_lo.out')
c      open(unit=17,file='output/externals/N2_ref_shms_1pass_lo.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_1pass_hi.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_1pass_hi.out')
c      open(unit=17,file='output/externals/N2_ref_hms_1pass_hi.out')
c      open(unit=17,file='output/externals/N2_ref_shms_1pass_hi.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_kin_3.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_kin_C.out')
c      open(unit=17,file='output/externals/N2_ref_hms_3.out')
c      open(unit=17,file='output/externals/N2_ref_shms_C.out')
c      open(unit=17,file='output/externals/d2n_3He_hms_kin_4.out')
c      open(unit=17,file='output/externals/d2n_3He_shms_kin_B.out')
c      open(unit=17,file='output/externals/N2_ref_hms_4.out')
c      open(unit=17,file='output/externals/N2_ref_shms_B.out')

c      open(unit=18,file='louk.out')
c      open(unit=19,file='ratioy.top')
c      open(unit=29,file='ratiow2.top')
c      write(19,'(1x,''set device postscript portrait'')')
c      write(29,'(1x,''set device postscript portrait'')')
c      write(18,1818)
c 1818 format(
c     >  1x,'                                   internal       ',
c     >  1x,'                     internal  + external'/
c     >  1x,'   E    Ep    Th     s_born      s_rad     s_elas',
c     >      '      s_dis      s_rad',
c     >     '     s_elas      s_dis')
CCCC      CALL READIN(LEN_FILENAME,FILENAME)                                                       

      do i=1,4
        do j=1,4
          x=0.1*i
          Q2=j
          CALL R1998(X,Q2,R,DR,GOOD)
c          write(6,'(1x,''x,q2,r,dr'',4f8.3)') x,q2,r,dr
        enddo
      enddo

c      write(6,'(1x,''got here'')')
c      write(6,'(''A b='',e11.3)') b
      CALL READIN_EXT()                                                       

c      write(6,'(1x,''got here'')')
c      write(6,'(''A b='',e11.3)') b
      CALL RADLENGTH(X0)                                                
c      write(6,'(''A b='',e11.3)') b
                                                                        
c read in params for d2model
c           write(6,*)'Reading parameters for Ioana''s model'
           open(unit=77,file
     >          ='input/d2model/parms.dat'
     >          ,status='old')
           do I=1,37
              read(77,*)j
              read(77,*)xval(i)
*              write(*,*)xval(i)
              read(77,*)tmptmp
              read(77,*)tmptmp
              read(77,*)tmptmp
           enddo	
           close(77)

C Constants for bremstrahlung and ionization spectra for target and for 
C aluminum:                                                             
c changed to uze iz=1 if neutron
      izz = max(iz,1)
      ETA  = LOG(1194./izz**0.6667)/LOG(184.15/izz**0.3333)               
      B    = 4./3.*(1.+1./9.*(izz+1.)/(izz+ETA)/LOG(184.15/izz**0.3333))   
      AX0  = 0.000154*izz/amuM*X0                                        
      BA   = 1.3667511                                                  
      AX0A =  .0017535                                                  
c      write(6,'(''b...'',5f10.3)') eta,b,ax0,ba,ax0a

C     TASUM=0.                                                          
C     TBSUM=0.                                                          
C Unless specified, set default number of target segments for integrals:
c                                                                        
cc took this out. nseg=0 means skip external
cc      IF (NSEG.LE.0) NSEG = 4                                           
                                                                        
      WRITE(10,'(                                                       
     >''    E1    TH    Q2 COS_K         SMR     NUC  '',     
     >''        INT-SMR'',/,                                               
     >''    EP     X    W2 QE-PK  INELA  Q-ELA  ELAST'',     
     >''   TOTAL  MODEL'',/)')                                              

c check out pauli supp.
      do ii=1,5,2
        do j=1,10
          e0set=float(ii)
          q2set=0.05*j
          do i=1,5
            call PAULI_SUPPRESSION(i-1,1,2,q2set,e0set,
     >        psf1(i),psf2(i))
          enddo
c          write(6,'(1x,12f5.2)') e0set,q2set,psf1,psf2
        enddo
      enddo

      
c Loop over kinematic points:                                           
      nfile=0
      if(dojan05) nfile = nfile + 1
      if(dofromin) nfile = nfile + 1
      if(dosimona) nfile = nfile + 1
      if(donucr) nfile = nfile + 1
      if(nfile.gt.1) then
        write(6,'(''cant do more than one special case at same time'')')
        stop
      endif
      nfile = 1
      if(dojan05) nfile=23
      if(dofromin) nfile=6
      if(dosimona) nfile=6
      if(donucr) nfile=28
c      if(donucr) nfile=2 ! for test
      if(donucr) open(unit=3,file='e99118/files.txt')
      do ifile=1,nfile
c      do ifile=20,20
      if(dojan05) then 
        close(unit=7)
        open(unit=7,file='jan05.dat')
      endif
      if(dofromin) then 
        close(unit=7)
c old        open(unit=7,file='fromin.newdat')
        open(unit=7,file='fomin.dat')
      endif
      if(dosimona) then 
        close(unit=7)
        open(unit=7,file='simona.dat')
      endif
      if(donucr) then 
        read(3,'(a)') fnamep
        write(fname,'(''e99118/h''a)') fnamep(2:20)
        if(IA.eq.2) write(fname,'(''e99118/d''a)') fnamep(2:20)
        close(unit=7)
        open(unit=7,file=fname)
      endif

      npts=0                                                                     
      w2last=-100.
      if(dofromin) w2last=100.
      DO 99 ii=1,100000                                                   
c           CALL TIME(VTSTART,TTSTART)                                   
        if(dojan05) then
          read(7,*,end=100) IAjan,E0SET,EPSET,THSET,
     >      sigexp(npts+1),sigexper(npts+1) 
          if(IAjan.ne.IA) goto 99
          if(abs(thset - janth(ifile)).gt.0.5) goto 99
          if(abs(e0set - jane0(ifile)).gt.0.5) goto 99
        elseif(dofromin) then
          e0set=5.77
          read(7,*,end=100) EPSET,THSET,Afromin,Zfromin,
     >      sigexp(npts+1),sigexper(npts+1) 
          if(abs(Afromin - float(IA)).gt.0.3) goto 99
          if(abs(Zfromin - float(IZ)).gt.0.3) goto 99
          if(abs(thset - thfromin(ifile)).gt.0.5) goto 99
c          write(6,'(1x,''fromin'',6f10.3)') EPSET,THSET,Afromin,Zfromin,
c     >      sigexp(npts+1),sigexper(npts+1) 
        elseif(dosimona) then
          read(7,'(i2,3f9.4,27x,3f11.4)',end=100) itrgsim,
     >      E0SET,THSET,EPSET,
     >      sigexp(npts+1),sigexper(npts+1),rcsim
          if(abs(thset - thsimona(ifile)).gt.0.5) goto 99
          if(itrgsim.eq.11.and.IA.eq.2) goto 99
          if(itrgsim.eq.13.and.IA.eq.1) goto 99
! for H2, values are born, so need to divide by rc used
          if(itrgsim.eq.11) then
            CALL SECNUCLW(E0SET,EPSET,TH,SIGMA_BORN)                     
            write(99,'(1x,3f6.2,f10.3)') E0SET,EPSET,TH,
     >        sigexp(npts+1)/SIGMA_BORN,
     >        sigexper(npts+1)/SIGMA_BORN
            sigexp(npts+1) = sigexp(npts+1) / rcsim
            sigexper(npts+1) = sigexper(npts+1) / rcsim
          endif
        elseif(donucr) then
          read(7,*,end=100) E0SET,THSET,EPSET,
     >      sigexp(npts+1),sigexper(npts+1) 
        elseif(doeg1b) then
c          read(7,'(1x,2i3,2f6.3,F7.3)') ipsv,ithsv,
c     >      E0SET,EPSET,THSET
c          read(7,'(3i3,3F8.3)',end=100) iesv,iqsv,iwsv,
c     >      E0SET,EPSET,THSET
          read(7,*,end=100) iesv,iqsv,iwsv,
     >      E0SET,EPSET,THSET
          if(iesv.eq.0) goto 100
        else
c           READ(7,'(1X,2F6.3,F7.3)',END=100) E0SET,EPSET,THSET                
c           READ(7,'(1X,F6.3,2F8.3)',END=100) E0SET,EPSET,THSET                
           READ(7,*,END=100) E0SET,EPSET,THSET                
c           write(6,'(1x,3f10.3)') E0SET,EPSET,THSET
           IF (E0SET.LE.0.) GOTO 100                                    
        endif
C          CALL INTERPOL(E0SET,EPSET,THSET)                             
                                                                        
           TH    = THSET                                                
           THR   = TH*PI/180.                                           
           SIN2  = SIN(THR/2.)**2                                       
           Q2SET = 4.*E0SET*EPSET*SIN2                                  
           ANU   = E0SET-EPSET                                          
           X     = Q2SET/2./PM/ANU                                      
           Y     = ANU/E0SET                                            
           EPS   = 1./( 1.+2.*SIN2/(1.-SIN2)*(1.+Q2SET/(2.*PM*X)**2))   
           W2    = PM**2+2.*PM*ANU-Q2SET                                
ccc        if(dojan05.and.(w2 - w2last).lt.0.2 .and.w2.gt.1.3) goto 99
           if(dojan05.and.(w2.lt.1.1.or.w2.gt.4.0)) goto 99
           if(dofromin.and.(w2last-w2).lt.0.1) goto 99
cc           if(iA.lt.1.5.and.w2.le.1.17.and.(.not.smrelasp)) goto 99
c           if(w2.le.0.) goto 99
           if(IA.gt.2.and.x.gt.float(ia).and.(.not.smrelasp)) goto 99
           if(x.gt.3.0) goto 99
           w2last=w2
           q28 = q2set
           w8 = sqrt(max(0.,w2))
           npts = npts+1
           e0sv(npts)=e0set
           epsv(npts)=epset
           thsv(npts)=thset
           xsv(npts)=x
           w2sv(npts)=w2
           costh = cos(thr)
           E0 = E0set
           EP = EPSEt
           QV = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH) 
           AMS   = avgM-PM 
           AMF2  = PM**2
           A = 4.*QV**2-4.*(E0-EP+avgM)**2 
           BB = -4.*QV*(AMS**2-AMF2-QV**2+(E0-EP+avgM)**2) 
           C = (E0-EP+avgM)**4+(AMS**2-AMF2-QV**2)**2
     >       -2.*(E0-EP+avgM)**2*(AMF2+QV**2+AMS**2)
           DISC  = BB**2-4.*A*C
           ysv(npts)=0.
           IF (A.ne.0..and.DISC.gt.0.) then
             Ysv(npts) = (-BB-SQRT(DISC))/2./A
             Ysvp(npts) = (-BB+SQRT(DISC))/2./A
           endif
! Default resolutions
! if using grid, then resolution smearing in W2 is
! hardwired: see line 3160!
! good for BARONS SOS with .2% and 4 mr
!*** Th resolution for eA elastic: should use A* 0.98 in future
c          smrdEp = sqrt((ep * 0.002)**2 +
c    >      (0.005 * (2. * e0 * ep * sin(thr)) /
c    >               (2. * 0.9383 + q2set / ep))**2) 
! this is for HMS: 0.15% and 1 mr
c          smrdEp = sqrt((ep * 0.0015)**2 +
c    >      (0.001 * (2. * e0 * ep * sin(thr)) /
c    >               (2. * 0.9383 + q2set / ep))**2)
! Make smaller for test
           smrdEp = sqrt((ep * 0.0005)**2 +
     >      (0.0001 * (2. * e0 * ep * sin(thr)) /
     >               (2. * 0.9383 + q2set / ep))**2)
C Resolution in W at W=M for particular experiment
           IF(INDEX(TARGET,'EG1b').GT.0) THEN 
             smrdEp = ep * 0.02 / sqrt(e0) 
           ENDIF
           IF(INDEX(TARGET,'EG4').GT.0) THEN 
             smrdEp = ep * 0.02 / sqrt(e0) 
           ENDIF
C Born cross section:                                    
           prttst=.true.
           CALL SECNUCLW(E0SET,EPSET,TH,SIGMA_BORN)                     
c           write(98,'(''sig born'',e11.4)') sigma_born

           doing_elas = .false.
           CALL QUASIY8(E0SET,EPSET,TH,SIGMA_CORR)                     

! added smeared elastic peak born
           doing_elas = .true.
           sigma_elastic = 0.
           if(smrelasp) CALL QUASIY8(E0SET,EPSET,TH,SIGMA_elastic)
c           write(6,'(''elastic'',3f8.3,e11.3)') 
c     >       E0SET,EPSET,TH,SIGMA_elastic
           doing_elas = .false.

           prttst=.false.
C          CALL  FYXSEC8(E0SET,EPSET,TH,SIGMA_CORR)                     
           write(*,'(1x,''e,ep,th,sigb,sigqe='',3f8.2,2e12.4)')
     >       E0SET,EPSET,TH,SIGMA_BORN,sigma_corr
C Radiation length integral. Loop over internal and internal+external:  
! If eg1b type target, replace index 1 with results with ttarg/2.
           IF(INDEX(TARGET,'EG1b').GT.0 .or.
     >        INDEX(TARGET,'EG4').GT.0) THEN 
             JTARG = 2               
             ttarg = ttarg / 2.0
             SIGMA_INEL(1) =0.
! inelastic
             if(w2.gt.1.17.or.IA.gt.1) CALL SIMP(
     >          0.,TTARG,NSEG,CONTINUUM,SIGMA_INEL(1))             
! quasielastic
             CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_QELA(1)) 
! elastic
             if(smrelasp)  then
               doing_elas = .true.
               CALL SIMP(0.,TTARG,NSEG,QETAIL,SIGMA_ELPK(1)) 
               doing_elas = .false.
             else
              IF(NUC_METHOD.EQ.0)  THEN
               CALL SIMP(0.,TTARG,NSEG,EPEAK,SIGMA_ELPK(1))
              ELSEIF(NUC_METHOD.EQ.1)  THEN
               CALL SIMP(0.,TTARG,NSEG,ATAILFL1,SIGMA_ELPK(1)) 
               SIGMA_ELPK(1) = 0.001 *SIGMA_ELPK(1) 
              ENDIF
             endif
             SIGMA_INEL(1) = SIGMA_INEL(1)/TTARG                       
             SIGMA_QELA(1) = SIGMA_QELA(1)/TTARG  
             SIGMA_ELPK(1) = SIGMA_ELPK(1)/TTARG
c             write(6,'(''tst 1'',f7.3,3e11.4)') ttarg,
c     >         sigma_elpk(1) , sigma_qela(1) , sigma_inel(1)
             ttarg = ttarg * 2.0
           ELSE
             JTARG = 1                                                    
             SIGMA_INEL(1) = 0.                               
! inelastic
             if(w2.gt.1.17.or.IA.gt.1) SIGMA_INEL(1) = CONTINUUM(0.)   
! quasielastic
             SIGMA_QELA(1) = QETAIL(0.)                                   
cc no longer used
             SIGMA_QEPK(1) = 0.
cc             SIGMA_QEPK(1) = QEPEAK(0.)
! Q-E using peak and cosk: no longer used
             SIGMA_QE_PK_COS(1) = 0.
cc             SIGMA_QE_PK_COS(1) = .001 *ATAILFL_QE(0.)
! elastic
             if(smrelasp)  then
               doing_elas = .true.
               SIGMA_ELPK(1) = QETAIL(0.)                                   
               doing_elas = .false.
             else
              IF(NUC_METHOD.EQ.0)  THEN
                SIGMA_ELPK(1) = EPEAK(0.)  !original external    
              ELSEIF(NUC_METHOD.EQ.1)  THEN
c              prttst=.true.
! E140 Mo-Tsai code (pb->nb)
               SIGMA_ELPK(1) = 0.001*ATAILFL1(0.) 
c              prttst=.false.
              ELSE
c               WRITE(6,'('' *** BAD parameter***  NUC_METHOD='',I2,
c     >           '' SHOULD BE 0 OR 1'')')NUC_METHOD 
                STOP
              ENDIF
             endif
           ENDIF

! now do with full radiation lengths
! (unless nseg=0)
           if(nseg.eq.0) then
             write(17,'(1x,5f7.3,6e11.4,3i3,4e11.4)') 
     >       e0set,epset,thset,xsv(npts),q2set,SIGMA_BORN,
     >       sigma_elpk(1)+sigma_qepk(1)+sigma_inel(1),
     >       sigma_elpk(1),sigma_qepk(1),sigma_inel(1)
             goto 99
           endif
           JTARG = 2               
! ineasltic
           SIGMA_INEL(2) =0.
           if(w2.gt.1.17.or.IA.gt.1) CALL SIMP(
     >         0.,TTARG,NSEG,CONTINUUM,SIGMA_INEL(2))             
! quasielastic
           CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_QELA(2)) 
! skip this
           SIGMA_QE_PK_COS(2) = 0.
cc           CALL SIMP(0.,TTARG,NSEG,ATAILFL_QE,SIGMA_QE_PK_COS(2))
C PERM.    CALL SIMP(0.,TTARG,NSEG,   QEPEAK,SIGMA_QEPK(2)) 
! elastic
           if(smrelasp)  then
             doing_elas = .true.
             CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_ELPK(2)) 
             doing_elas = .false.
           else
            IF(NUC_METHOD.EQ.0)  THEN
              CALL SIMP(0.,TTARG,NSEG,EPEAK,SIGMA_ELPK(2))  !original external
c              write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
            ELSEIF(NUC_METHOD.EQ.1)  THEN
              CALL SIMP(0.,TTARG,NSEG,ATAILFL1,SIGMA_ELPK(2)) 
              SIGMA_ELPK(2) = 0.001 *SIGMA_ELPK(2)   ! pb->nb 
c              write(*,*) '###########################################'
            ENDIF
           endif

cc no longer used
           SIGMA_QE_PK_COS(2) = 0.
cc           SIGMA_QE_PK_COS(2) = 0.001 *SIGMA_QE_PK_COS(2)  !pb->nb
           write(*,111) ttarg,SIGMA_INEL(2), SIGMA_QELA(2),
     >       SIGMA_QE_PK_COS(2), SIGMA_ELPK(2)
 111       format(1x,'ttarg,sig=',f8.4,4e11.3)
           SIGMA_INEL(2) = SIGMA_INEL(2)/TTARG                          
           SIGMA_QELA(2) = SIGMA_QELA(2)/TTARG  
cc           SIGMA_QE_PK_COS(2) = SIGMA_QE_PK_COS(2)/TTARG                          
           SIGMA_ELPK(2) = SIGMA_ELPK(2)/TTARG                          
c             write(6,'(''tst 2'',f7.3,3e11.4)') ttarg,
c     >         sigma_elpk(2) , sigma_qela(2) , sigma_inel(2)
! Total inelastic (DIS + Res. + q.e.) born xsection
           SIGMA_BORNPLSQE = SIGMA_BORN+SIGMA_CORR
           DO 98 I=1,2                                                  
!*** changed to divide by SIGMA_BORNPLSQE in all cases (before
! was only for QELA).
            DELTA_INEL(I) = (SIGMA_INEL(I)/ SIGMA_BORNPLSQE-1.)*100.         
            DELTA_QELA(I) = (SIGMA_QELA(I)/ SIGMA_BORNPLSQE)*100.
            DELTA_QEPK(I) = (SIGMA_QEPK(I)/ SIGMA_BORNPLSQE)*100. 
            DELTA_QE_PK_COS(I) = (SIGMA_QE_PK_COS(I)/ SIGMA_BORNPLSQE)*100.
            DELTA_ELPK(I) = (SIGMA_ELPK(I)/ SIGMA_BORNPLSQE)*100.            
            TOTAL(I)      =  DELTA_INEL(I)+DELTA_QELA(I)+DELTA_ELPK(I)  
 98        CONTINUE                                                     
                                                                        
           RAT_SMR = 1.                                                 
           IF(DELTA_QEPK(1).NE.0.0.and.w2.gt.0.95) 
     >       RAT_SMR = DELTA_QELA(1)/DELTA_QEPK(1)
                                                                        
c           CALL TIME(VTSTOP,TTSTOP)                                     
           VTRUN = VTSTOP-VTSTART                                       
           TTRUN = TTSTOP-TTSTART                                       
C          WRITE(11,'(1X,2F8.4)') TBSUM/12.,TASUM/12.                   
C          TASUM=0.                                                     
C          TBSUM=0.                                                     
c            WRITE(10,'(1X,2f7.3,f6.2,5F8.2,F7.3/
c    >          1x,2f7.3,F6.2,5f8.2,1PE10.3,/)')
c    >          E0SET, TH, Q2SET, DELTA_QE_PK_COS(1),  
c    >          DELTA_INEL(1),DELTA_QELA(1), DELTA_ELPK(1), 
c    >          TOTAL(1), RAT_SMR,EPSET, x, W2, 
c    >          DELTA_QE_PK_COS(2), DELTA_INEL(2),DELTA_QELA(2), 
c    >          DELTA_ELPK(2), TOTAL(2), SIGMA_BORNPLSQE  
             ddnom = SIGMA_BORN + SIGMA_CORR + sigma_elastic
             if(ddnom.ne.0.) WRITE(10,'(4f6.2,e10.3,9f6.2)')
     >          E0, EP, TH, w2, ddnom, 
     >          sigma_elastic/ddnom,  
     >          sigma_corr/ddnom,  
     >          sigma_born/ddnom,  
     >          sigma_elpk(2)/ddnom,
     >          sigma_qela(2)/ddnom,
     >          sigma_inel(2)/ddnom
      WRITE(*,'(1X,2f7.3,f6.2,5F8.2,F7.3/1x,2f7.3,F6.2,5f8.2,    
     >           1PE10.3,/)') 
     >    E0SET, TH, Q2SET, DELTA_QE_PK_COS(1),  DELTA_INEL(1),
     >       DELTA_QELA(1), DELTA_ELPK(1), TOTAL(1), RAT_SMR,
     >    EPSET,  x,    W2, DELTA_QE_PK_COS(2), DELTA_INEL(2),
     >       DELTA_QELA(2), DELTA_ELPK(2), TOTAL(2), SIGMA_BORNPLSQE  
      sbsv(npts)= SIGMA_BORNPLSQE
      srsv(npts)= sigma_qela(2) + sigma_elpk(2) + sigma_inel(2)
      srqesv(npts) = sigma_qela(2)
      srelsv(npts) = sigma_elpk(2)
      srinsv(npts) = sigma_inel(2)

c      write(17,'(36x,2e11.4)') sigma_born,sigma_corr 
      if(srsv(npts).ne.0.) then 
         write(17,'(1x,5f7.3,6e11.4,3i3,4e11.4)') 
     >        e0sv(npts),epsv(npts),thsv(npts),
     >        xsv(npts),q2set,
     >        max(0.,sbsv(npts)),
     >        max(0.,srsv(npts)),
     >        max(0.,srelsv(npts)),
     >        max(0.,srqesv(npts)),
     >        max(0.,srinsv(npts)),
     >        max(0.,SIGMA_BORN),iesv,iqsv,
     >        iwsv,
c     >        max(0.,sigma_elpk(1)),
     >        max(0.,sigma_qela(1) + sigma_elpk(1) + sigma_inel(1)),
     >        max(0.,sigma_elastic)
      endif
c      write(18,'(1x,4f10.3)') e0sv(npts),epsv(npts),thsv(npts),
c     >  1./(1+total(2)/100.)
c      write(18,'(1x,f4.1,f6.2,f6.1,7e11.4)') 
c     >  e0sv(npts),epsv(npts),thsv(npts),
c     >  sbsv(npts),SIGMA_BORNPLSQE*(1+total(1)/100.),
c     >  sigma_elpk(1),sigma_inel(1),srsv(npts),
c     >  srelsv(npts),srinsv(npts)
 99   CONTINUE                                                          
100   CONTINUE                                                          
      if(dojan05.or.dofromin.or.dosimona.or.donucr) then
       ymax=0.
       do i=1,npts
         ymax = max(ymax, sigexp(i))
! if Fromin, plot versus y instead of W2
         write(44,'(1x,''w2,y'',3f8.3)') w2sv(i),ysv(i),
     >     ysvp(i)
         xplt(i) = sqrt(max(0.,w2sv(i)))
         if(dofromin) xplt(i) = ysv(i)
       enddo
       if(dojan05) ymin=ymax/100.
       if(donucr) ymin=ymax/10.
       if(dofromin) ymin=ymax/100000.
       if(dosimona) ymin=ymax/300.
       if(dofromin) thset=thfromin(ifile)
       if(dosimona) thset=thsimona(ifile)
       if(dojan05) thset = janth(ifile)
       if(dojan05) E0set = jane0(ifile)
c       open(unit=16,file='extern.top')
c       write(16,'(1x,''set device postscript portrait'')')
c       write(16,1626) IA,IZ,E0set,THSET,ymin,ymax*1.1
 1626  format(1x,'new frame'/
     >  1x,'title top ',1h','A=',i3,' Z=',i2,' E0=',f5.3,
     >    ' Th=',f4.1,1h'/
     >  1x,'set bar size 0. ; set order x y dy '/
     >  1x,'set scale y log ; set limits y  ',2e12.5/
     >  1x,'set sym 9O size 0.5'/
     >  1x,'title left ',1h','dsigma (nb/sr/gev)',1h')
c       write(19,1926) IA,IZ,E0set,THSET
c       write(29,1926) IA,IZ,E0set,THSET
 1926  format(1x,'new frame'/
     >  1x,'title top ',1h','A=',i3,' Z=',i2,' E0=',f5.3,
     >    ' Th=',f4.1,1h'/
     >  1x,'set bar size 0. ; set order x y dy '/
     >  1x,'set scale y log ; set limits y  0.3 3.0'/
     >  1x,'set sym 9O size 2.0'/
     >  1x,'title left ',1h','data/model',1h')
       if(dofromin) then
c        write(16,1636)
 1636    format(1x,'title bottom ',1h','y',1h')
       else
c         write(16,1637)
 1637    format(1x,'title bottom ',1h','W2 (GeV)',1h')
       endif
c       write(19,1636)
c       write(29,1637)
c       write(16,'(1x,f8.4,2e12.4)') (xplt(i),
c     >   sigexp(i),sigexper(i),i=1,npts)
c       write(16,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(16,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(16,'(1x,''plot ; set symbol 9O size 2.0'')')
       fact = 1000. * float(IA)
       if(dofromin) fact=1000.
       if(donucr) fact=1000.
       if(dosimona) fact=1.
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srsv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 1O size 2.0'')')
       if(dofromin) write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srelsv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 2O size 2.0'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srqesv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 3O size 2.0'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srinsv(i)/fact,i=1,npts)
c      write(16,'(1x,''plot ; set symbol 4O size 2.0'')')
       do i=1,npts
         if(srsv(i).ne.0.) then
c           write(19,'(1x,f8.4,2e12.4)') ysv(i),
c     >       max(0.3,min(3.,
c     >       sigexp(i)/srsv(i)*fact)),
c     >       sigexper(i)/srsv(i)*fact
c           write(29,'(1x,f8.4,2e12.4)') w2sv(i),
c     >       max(0.3,min(3.,
c     >       sigexp(i)/srsv(i)*fact)),
c     >       sigexper(i)/srsv(i)*fact
         endif
       enddo
c       write(19,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(19,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(19,'(1x,''plot ; set symbol 9O size 2.0'')')
c       write(29,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(29,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(29,'(1x,''plot ; set symbol 9O size 2.0'')')
c       write(19,'(1x,''-10. 1. ; 100. 1. ; join dash'')')
c       write(29,'(1x,''-10. 1. ; 100. 1. ; join dash'')')
      else
c       write(16,1616)
c 1616  format(1x,'title bottom ',1h','W**2 (GeV)**2',1h'/
c     >  1x,'title left ',1h','dsigma (nb/sr/gev)',1h')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),sbsv(i),i=1,npts)
c       write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srsv(i),i=1,npts)
c       write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srelsv(i),i=1,npts)
c       write(16,'(1x,''join dot'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srqesv(i),i=1,npts)
c       write(16,'(1x,''join dash'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srinsv(i),i=1,npts)
c       write(16,'(1x,''join dotdash'')')
      endif
      WRITE(10,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,4F6.2,F6.3,/,25X,6F6.2,    
     >           1PE10.3,/)') (ZERO(i),i=1,18)                          

! loop over files
      enddo 
      CLOSE(10)
                                                                        
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION CONTINUUM(T)                                             
      IMPLICIT NONE
      REAL T                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /CON/ PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                              
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                        
      COMMON /ERR/ REPS1,AEPS1  
      REAL  REPS1,AEPS1                                        
      COMMON /SIG/ CSTYPE                                               
      CHARACTER*1  CSTYPE
      REAL E0LO,E0HI,ANS4,ANS3,ANS2,ANS1,EPHI,EPLO,QUADMO_R
      INTEGER NLVL
      REAL FUNC1,FUNC2,FUNC3,FUNC4
      EXTERNAL     FUNC1,FUNC2,FUNC3,FUNC4                              
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                          
                                                                        
      CSTYPE = 'I'                                                      
      CALL RADIATORS(T)                                                 
C Soft photon part; E0-R*DELTA ---> E0 and EP ---> EP+DELTA:            
                                                                        
      ANS1 = FUNC1(T)                                                   
                                                                        
C Soft photons along initial electron only; single integration.         
C E0-R*DELTA --> E0 analytically and EP+DELTA --> EPMAX(E0) numerically:
                                                                        
c PEB modified to use avgM for IA>1 10/05
      EPLO = EP+DELTA                                                   
      if(IA.le.1) then
        EPHI = (E0-PP-PP**2/2./PM)/(1.+E0*(1.-COS(THR))/PM)               
      else
        EPHI = (E0-PP-PP**2/2./avgM)/(1.+E0*(1.-COS(THR))/avgM)        
      endif
C$$      ANS2 = RGAUSS(FUNC2,EPLO,EPHI,REPS1) 
      if(ephi.le.eplo) then
cc      write(6,'(1x,''error, eplo,ephi='',2f8.3)') eplo,ephi
        ans2=0.
      else
        ANS2 = QUADMO_R(FUNC2,EPLO,EPHI,REPS1,NLVL) 
      endif

C Soft photons along final electron only; single integration.           
C EP --> EP+DELTA analytically and E0MIN(EP) --> E0-R*DELTA numerically:
                                                                        
c PEB modified to use avgM for IA>1 10/05
      if(IA.le.1) then
        E0LO = (EP+PP+PP**2/2./PM)/(1.-EP*(1.-COS(THR))/PM)               
      else
        E0LO = (EP+PP+PP**2/2./avgM)/(1.-EP*(1.-COS(THR))/avgM)          
      endif
      E0HI = E0-R*DELTA                                                 
C$$     ANS3 = RGAUSS(FUNC3,E0LO,E0HI,REPS1)                              
      if(ephi.le.eplo) then
c        write(6,'(1x,''error, e0lo,e0hi='',2f8.3)') e0lo,e0hi
        ans3=0.
      else
        ANS3 = QUADMO_R(FUNC3,E0LO,E0HI,REPS1,NLVL)    
      endif

C Hard photon part; double integration.  Set up the integration over the
C initial energy of the electrons.                                      
                                                                        
c PEB modified to use avgM for IA>1 10/05
      if(IA.le.1) then
        E0LO = (EPLO+PP+PP**2/2./PM)/(1.-EPLO*(1.-COS(THR))/PM)           
      else
        E0LO = (EPLO+PP+PP**2/2./avgM)/(1.-EPLO*(1.-COS(THR))/avgM)   
      endif
C$$      ANS4 = RGAUSS(FUNC4,E0LO,E0HI,REPS1)                              
      if(ephi.le.eplo) then
cc      write(6,'(1x,''error4, e0lo,e0hi='',2f8.3)') e0lo,e0hi
        ans4=0.
      else
        ANS4 = QUADMO_R(FUNC4,E0LO,E0HI,REPS1,NLVL)  
      endif
      CONTINUUM = ANS1+ANS2+ANS3+ANS4                                   
cc      if(CONTINUUM.gt.0.) 
c       write(6,'(1x,''continuum='',5e12.4)') ans1,ans2,ans3,ans4,
c     >  CONTINUUM                                                                       
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION QETAIL(T)                                                
      IMPLICIT NONE
      REAL T                                                            
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                    
      COMMON /CON/   PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL                                      
      COMMON /OMG/   DELTA,R 
      REAL   DELTA,R                                              
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL   
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                  
      COMMON /ERR/   REPS1,AEPS1
      REAL  REPS1,AEPS1                                        
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                          
      COMMON /SIG/   CSTYPE                                             
      CHARACTER*1    CSTYPE
      REAL*4 ANS,ANS1,ANS2,ANS3,ANS4,EPLO,EPCN,EPHI,E0LO,E0CN,E0HI
      INTEGER NLVL
      REAL*4   FUNC1,FUNC2,FUNC3,FUNC4Q,QUADMO_R
      EXTERNAL       FUNC1,FUNC2,FUNC3,FUNC4Q                      
      logical doing_elas
      common/doelas/ doing_elas

      QETAIL = 0.                                                    
      CSTYPE = 'Q'                                                 
      CALL RADIATORS(T)                                           
                                                                  
C Soft photon part; E0-R*DELTA ---> E0 and EP ---> EP+DELTA:       
                                                                     
      ANS1 = FUNC1(T)                                          
                                                                     
C Soft photons along initial electron only; single integration.     
C E0-R*DELTA ---> E0 analytically and EP Integral numerically.        
C lower limit defined by kinematic point and upper limit by x=M_targ.
                                                                      
      EPLO = EP+DELTA                                                 
      EPCN = E0/(1.+E0*(1.-COS(THR))/PM)                             
      EPHI = E0/(1.+E0*(1.-COS(THR))/avgM)                           
      EPCN = MAX(EPLO,EPCN)                                           
ccc pyb 04/07
ccc      EPHI = E0 - DELTA

      ANS  = QUADMO_R(FUNC2,EPLO,EPCN,REPS1,NLVL)
      ANS2 = QUADMO_R(FUNC2,EPCN,EPHI,REPS1,NLVL)                     
      ANS2 = ANS + ANS2                                               
                                                                 
C Soft photons along final electron only; single integration.         
C EP ---> EP+DELTA analytically and E0 Integral numerically.          
C upper limit defined by kinematic point and lower limit by x=M_targ. 
                                                                      
      E0LO = EP/(1.-EP*(1.-COS(THR))/avgM)                            
      E0CN = (EP+PP+PP**2/2./PM)/(1.-EP*(1.-COS(THR))/PM)             
      E0HI = E0-R*DELTA                                               
      E0CN = MIN(E0CN,E0HI)                                           
ccc pyb 04/07
ccc      E0HI = E0 - DELTA


      ANS  = QUADMO_R(FUNC3,E0LO,E0CN,REPS1,NLVL)  !NEW
      ANS3 = QUADMO_R(FUNC3,E0CN,E0HI,REPS1,NLVL)  !NEW
      ANS3 = ANS + ANS3                                               
                                                                      
C Hard photon part; double integration.  Set up the integration over t
C initial energy of the electrons. upper limit defined by delta cutof 
C and lower limit by x=M_target,above the W2=M_proton+M_pion for EP   
C integral.                                                           
                                                              
      E0LO = EPLO/(1.-EPLO*(1.-COS(THR))/avgM)                       
      ANS4 = QUADMO_R(FUNC4Q,E0LO,E0HI,0.01,NLVL)  
                                     
      QETAIL = ANS1 + ANS2 + ANS3 + ANS4                       

      if(doing_elas) then
c        write(6,'(''q'',7f6.3,4e10.3)') t,e0,ep,eplo,ephi,
c     >    e0lo,e0hi,ans1, ans2,
c     >    ans3, ans4
      endif
      RETURN                                                         
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION QEPEAK(T)                                                
                                                                        
C Function to integrate the quasi elastic radiative cross section in    
C the delta function approximation for the peak in the cross section.   
C The only integral is over the E0.                                     
      IMPLICIT NONE
      REAL*4 T                                                                        
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /CON/ PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                               
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL   
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                      
      COMMON /ERR/ REPS1,AEPS1
      REAL  REPS1,AEPS1                                                  
      COMMON /SIG/ CSTYPE                                               
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM !put in by SER 4/14/93  
      INTEGER IZ,IA 
      REAL AVGN,AVGA,AVGM,AMUM   
      REAL*4 EPP,EPLO,E0HI,E0LO,E0P,FAC,SIGEFF,AI11,AI1,AI2,
     >       ANS1,ANS2,ANS3,SIGMA
      REAL*4 FTSAI,BREMS,QUADMO_R
      INTEGER NLVL
      CHARACTER*1  CSTYPE                                               
      EXTERNAL     FUNCQE                                               
                                                                        
      QEPEAK = 0.                                                       
      IF (iA.LE.1) RETURN                                               
      CSTYPE = 'P'                                                      
      CALL RADIATORS(T)                                                 
                                                                        
C Soft photon part; E0-R*DELTA ---> E0 done analytically along incident 
C electron. EP for elastic peak for incident electron of energy E0:     
                                                                        
      EPP    = PM*E0/(PM+2.*E0*SIN(THR/2.)**2)                          
      CALL CROSS(E0,EPP,TH,SIGMA)   !EPP not used here
      SIGEFF = SIGMA*FTSAI(E0,EPP,THR)                                  
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/E0 << 1.                                                    
                                                                        
      AI11 = 1.+0.5772*BTBI-0.62*BTBI**2                                
      AI1  = AI11*(R*DELTA/E0)**BTB*(1.-ATBX0/(1.-BTBI)/R/DELTA)        
      AI2  = BREMS(EPP,EPP-EP,TA,2)                                     
      ANS1 = AI1*SIGEFF*AI2                                             
                                                                        
C Hard photon part; E0LO ---> E0-R*DELTA done numerically:              
                                                                        
      EPLO = EP+DELTA                                                   
      E0HI = E0-R*DELTA                                                 
      E0LO = PM*EPLO/(PM-2.*EPLO*SIN(THR/2.)**2)                        
C$$      ANS2 = RGAUSS(FUNCQE,E0LO,E0HI,REPS1)                             
      ANS2 = QUADMO_R(FUNCQE,E0LO,E0HI,REPS1,NLVL)
C Soft photon part; E0P(EP+DELTA) ---> E0P(EP) done analytically along  
C scattered electron. E0 for elastic peak for scattered electron of     
C energy EP:                                                            
                                                                        
      E0P    = PM*EP/(PM-2.*EP*SIN(THR/2.)**2)                          
      CALL CROSS(E0P,EP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0P,EP,THR)                                  
      FAC    = (PM+2.*E0P*SIN(THR/2.)**2)/(PM-2.*EP*SIN(THR/2.)**2)     
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/EP << 1.                                                    
                                                                        
      AI1  = BREMS(E0,E0-E0P,TB,1)                                      
      AI11 = 1.+0.5772*BTAI-0.62*BTAI**2                                
      AI2  = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)            
      ANS3 = AI1*SIGEFF*FAC*AI2                                         
                                                                        
      QEPEAK = ANS1+ANS2+ANS3                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION EPEAK(T)                                                 
                                                                        
C Function to integrate the quasi elastic radiative cross section in    
C the delta function approximation for the peak in the cross section.   
C The only integral is over the E0.                                     
      IMPLICIT NONE
      REAL*4 T                                                                  
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                     
      COMMON /CON/   PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                      
      COMMON /OMG/   DELTA,R 
      REAL   DELTA,R                                            
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL     
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      COMMON /ERR/   REPS1,AEPS1
      REAL  REPS1,AEPS1                                        
      COMMON /SIG/   CSTYPE
      CHARACTER*1    CSTYPE                                                                                          
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA 
      REAL AVGN,AVGA,AVGM,AMUM             
      COMMON /PRECISION/ DREPS1,DAEPS1
      REAL*8 DREPS1,DAEPS1 

      REAL*4 EPP,SIGMA,SIGEFF,FTSAI,FAC,AI11,AI1,AI2,EPLO,E0HI,E0LO,
     >       E0P,ANS1,ANS2,ANS3
      REAL*4 QUADMO_R,BREMS       
      INTEGER NLVL

      EXTERNAL       FUNCE
                                                                        
      CSTYPE = 'E'                                                      
      CALL RADIATORS(T)                                                 
      DAEPS1=1.D-17                       
                                                 
C Soft photon part; E0-R*DELTA ---> E0 done analytically along incident 
C electron. EP for elastic peak for incident electron of energy E0:     
                                                                        
      EPP    = avgM*E0/(avgM+2.*E0*SIN(THR/2.)**2)                      
      CALL CROSS(E0,EPP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0,EPP,THR)                                  
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/E0 << 1.                                                    
                                                                        
      AI11 = 1.+0.5772*BTBI-0.62*BTBI**2                                
      AI1  = AI11*(R*DELTA/E0)**BTB*(1.-ATBX0/(1.-BTBI)/R/DELTA)        
      AI2  = BREMS(EPP,EPP-EP,TA,2)                                     
      ANS1 = AI1*SIGEFF*AI2                                             
                                                                        
C Hard photon part; E0LO ---> E0-R*DELTA done numerically:              
                                                                        
      EPLO = EP+DELTA                                                   
      E0HI = E0-R*DELTA                                                 
      E0LO = avgM*EPLO/(avgM-2.*EPLO*SIN(THR/2.)**2)                    
C$$      ANS2 = RGAUSS(FUNCE,E0LO,E0HI,REPS1)                              
      ANS2 = QUADMO_R(FUNCE,E0LO,E0HI,REPS1,NLVL)       
C Soft photon part; E0P(EP+DELTA) ---> E0P(EP) done analytically along  
C scattered electron. E0 for elastic peak for scattered electron of     
C energy EP:                                                            
                                                                        
      E0P    = avgM*EP/(avgM-2.*EP*SIN(THR/2.)**2)                      
      CALL CROSS(E0P,EP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0P,EP,THR)                                  
      FAC    = (avgM+2.*E0P*SIN(THR/2.)**2)/(avgM-2.*EP*SIN(THR/2.)**2) 
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/EP << 1.                                                    
                                                                        
      AI1  = BREMS(E0,E0-E0P,TB,1)                                      
      AI11 = 1.+0.5772*BTAI-0.62*BTAI**2                                
      AI2  = AI11 * (DELTA/EP)**BTA * (1.- ATAX0/(1.-BTAI)/DELTA)       
      ANS3 = AI1*SIGEFF*FAC*AI2                                         
                                                                        
      EPEAK = ANS1+ANS2+ANS3                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC1(T)                                                 
      IMPLICIT NONE
      REAL*4 T                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                      
      COMMON /OMG/ DELTA,R  
      REAL   DELTA,R                                             
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL*4 SIGMA,SIGEFF,AI11,AI1,AI2,FTSAI                                                                             
C The constant part of the radiative cross section                      
                                                                        
      CALL CROSS(E0,EP,TH,SIGMA)                                        
      SIGEFF = SIGMA*FTSAI(E0,EP,THR)                                   
      AI11   = 1. + 0.5772 * BTBI - 0.62*BTBI**2                              
      AI1    = AI11 * (R*DELTA/E0)**BTB * (1.-ATBX0 / 
     >  (1.-BTBI) / R / DELTA)      
      AI11   = 1.+0.5772*BTAI-0.62*BTAI**2                              
      AI2    = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)          
      FUNC1  = SIGEFF*AI1*AI2                                           

c      write(6,'(''func1'',2e10.2,9f7.3)') sigma,sigeff,AI11,AI1,R,
c     >   delta,e0,btb,atbx0,BTBI,r

      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC2(EPP)                                               
      IMPLICIT NONE
      REAL*4 EPP
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                                  
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL*4 SIGMA,SIGEFF,AI11,AI1,AI2,FTSAI,BREMS
C The E0P integral is done analytically to obtain AI1                   
                                                                        
      CALL CROSS(E0,EPP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0,EPP,THR)                                  
      AI11   = 1.+0.5772*BTBI-0.62*BTBI**2                              
      AI1    = AI11*(R*DELTA/E0)**BTB*(1.-ATBX0/(1.-BTBI)/R/DELTA)      
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNC2  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC3(E0P)                                               
      IMPLICIT NONE
      REAL*4 E0P
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2  
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                            
      COMMON /OMG/ DELTA,R     
      REAL   DELTA,R                                         
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL    
      REAL*4 SIGMA,SIGEFF,FTSAI,AI1,AI11,AI2,BREMS

C The EPP integral is done analytically to obtain AI2                   
                                                                        
      CALL CROSS(E0P,EP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0P,EP,THR)                                  
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI11   = 1.+0.5772*BTAI-0.62*BTAI**2                              
      AI2    = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)          
      FUNC3  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4(E0P)                                               
      IMPLICIT NONE
      REAL E0P                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                        
      COMMON /CON/ PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                              
      COMMON /FUNCOM/ E0PP 
      REAL*4 E0PP                                                
      COMMON /ERR/ REPS1,AEPS1 
      REAL  REPS1,AEPS1
      REAL*4 EPLO,EPHI,QUADMO_U
      INTEGER NLVL
      EXTERNAL     FUNC4P                                               
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
                                                                        
C Set up the EPP integration                                            
                                                                        
      E0PP  = E0P                                                       
      EPLO  = EP+DELTA                                                  
c PEB fixed to use avgM if IA>1 10/05
      if(IA.le.1) then
        EPHI  = (E0P-PP-PP**2/2./PM)/(1.+E0P*(1.-COS(THR))/PM)            
      else
        EPHI  = (E0P-PP-PP**2/2./avgm)/(1.+E0P*(1.-COS(THR))/avgM)     
      endif
C$$      FUNC4 = UGAUSS(FUNC4P,EPLO,EPHI,REPS1)                            
      FUNC4 = QUADMO_U(FUNC4P,EPLO,EPHI,REPS1,NLVL)   
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4Q(E0P)
      IMPLICIT NONE
      REAL E0P
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                    
      COMMON /CON/   PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL                                      
      COMMON /OMG/   DELTA,R
      REAL   DELTA,R                                             
      COMMON /FUNCOM/   E0PP 
      REAL*4 E0PP                                                
      COMMON /ERR/   REPS1,AEPS1 
      REAL  REPS1,AEPS1                                       
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*4 EPLO,EPCN,EPHI, ANS1,ANS2,QUADMO_U
      INTEGER NLVL
      EXTERNAL     FUNC4P                                               
                                                                        
C Set up the EPP integration from kinematic limit to x=M_target         
C quasi elastic case                                                    
                                                                        
      E0PP = E0P                                                        
      EPLO = EP+DELTA                                                   
      EPCN = (E0P-PP-PP**2/2./PM)/(1.+E0P*(1.-COS(THR))/PM)             
      EPHI = E0P/(1.+E0P*(1.-COS(THR))/avgM)                            
      EPCN = MAX(EPLO,EPCN)                                             
                                                                        
C$$      CALL SIMP2(EPLO,EPCN,100,FUNC4P,ANS1)                             
C$$      CALL SIMP2(EPCN,EPHI,100,FUNC4P,ANS2)                             

      ANS1 = QUADMO_U(FUNC4P,EPLO,EPCN,REPS1,NLVL)  !NEW
      ANS2 = QUADMO_U(FUNC4P,EPCN,EPHI,REPS1,NLVL)  !NEW
      FUNC4Q = ANS1+ANS2                                                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4P(EPP)                                              
      IMPLICIT NONE
      REAL*4 EPP                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                      
      COMMON /FUNCOM/ E0P                                                  
      REAL*4 E0P
      REAL*4 SIGMA,SIGEFF,AI1,AI2,BREMS,FTSAI                                                                        
C The complete integrand in (C.1)                                       
C Sigma_r is of course approximated by equivalent radiator method       
                                                                        
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
c     if(E0-E0P.le.0.) write(6,'(1x,''ERROR,e0,e0p='',2f8.4)')
c    >  e0,e0p
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
c     if(EPP-EP.le.0.) write(*,'(1x,''ERROR,epp,ep='',2f8.4)')
c    >  ep,epp
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNC4P = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNCQE(E0P)                                              
      IMPLICIT NONE
      REAL E0P                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                      
      COMMON /CON/ PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                       
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL       
      REAL*4 EPP,SIGMA,SIGEFF,AI1,AI2,FTSAI,BREMS
C Integrand after the delta function killed the EPP integral. EPP is now
C constrained as below:                                                 
                                                                        
      EPP    = PM*E0P/(PM+2.*E0P*SIN(THR/2.)**2)                        
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNCQE = SIGEFF*AI1*AI2                                           
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNCE(E0P)              
      IMPLICIT NONE
      REAL E0P                                                                        
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                     
      COMMON /CON/   PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                     
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                       
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*4 EPP,SIGMA,SIGEFF,AI1,AI2,FTSAI,BREMS
                                                                        
C For Elastic electron-Nucleus collision as opposed to e-nucleon.       
C Integrand after the delta function killed the EPP integral            
C EPP is now constrained as below; Peak is at x=target_mass             
                                                                        
      EPP    = avgM*E0P/(avgM+2.*E0P*SIN(THR/2.)**2)                    
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNCE  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================


                                                             
      SUBROUTINE RADIATORS(T)                                           
      IMPLICIT NONE
      COMMON /SET/   E0SET,EPSET,THSET 
      REAL  E0SET,EPSET,THSET                                    
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                    
      COMMON /CON/   PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                    
      COMMON /OMG/   DELTA,R 
      REAL   DELTA,R                                             
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA     
      COMMON /TRL/TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG 
      REAL  TTARG,TWALL,TBEAM,TSPEC
      INTEGER NSEG,JTARG,ITARG,nprnt          
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM,thwin
      REAL DELP1,DEL01,DEL02,DELP2,DELP,TBEFOR,TAFTER,TAAL,TBAL,
     >  TEQUI,DEL0,T,x_end,ztarg,cryo_cm,wall_cm,cell_diam,
     >  cell_wall
      REAL ECIR,ECOR,ENTEC,TARGLEN,tcm,c_1,c_2,rho_tar,ang_b
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      save                                                                        
C If Internal corrections only, no target material only equiv. radiator:
                                                                        
      TBEFOR = T                                                        
      TAFTER = 0.                                                       
      TBAL   = 0.                                                       
      TAAL   = 0.                                                       
      ATBX0  = 0.                                                       
      ATAX0  = 0.                                                       
      BTBI   = 0.                                                       
      BTAI   = 0.                                                       
      DEL0   = 0.                                                       
      DELP   = 0.                                                       
      IF (JTARG.EQ.1) GO TO 1                                           
                                                                        
C Radiator quantities. Ionization spectrum should not have equivalent   
C radiator. Correct factors AX0 are used for the aluminium and target   
C parts:                                                                

C     IF (TARGET.EQ.'V_CYL') CALL V_CYL   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_087') CALL E_087   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_089') CALL E_089   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_140') CALL E_140   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E140XH1') CALL E140XH1 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XH2') CALL E140XH2 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XD1') CALL E140XD1 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XD2') CALL E140XD2 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'OTHER') CALL USERTARG(T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'SOLID') TAFTER = (TTARG-T)/COS(THR)                
C     TBAL  = TBAL+TBEAM                                                
C     TAAL  = TAAL+TSPEC                                                
C
c default for a simple target
      TBEFOR=t                        ! target
      TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
      TBAL= TBEAM                     ! material before the target
      TAAL= TSPEC                     ! material after the target

      IF(INDEX(TARGET,'E154').GT.0) THEN                    
       CALL RADLEN154(TARGET,T,TH,TAFTER,TBAL,TAAL,TBEFOR)  
      ELSEIF(INDEX(TARGET,'E142').GT.0) THEN                                
       CALL RADLEN42(TARGET,T,THR,TAFTER,TBAL,TAAL)              
      ELSEIF(INDEX(TARGET,'E143').GT.0) THEN                            
       CALL RADLEN43(TARGET,T,THR,TAFTER,TBAL,TAAL)              
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN                            
       CALL RADLEN49(TARGET,T,THR,TAFTER,TBAL,TAAL)              
c     ELSEIF(ITARG.LE.7) THEN                                           
c      CALL RADLENGT2(ITARG,T,THR,TAFTER,TBAL,TAAL)                     
      ENDIF
      IF(INDEX(TARGET,'EG1b').GT.0) THEN  ! Eg1btargets
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
        TBAL= 0.0024               ! Three 71 um Al windows (incl. banjo)
! get angle of exit window assuming radius of curvature is
! 3.5 times larger than distance to window (i.e. 35 inches
! versus 10 inches from target to window)
        thwin = thr / 3.5
        TAAL= 0.0014 / cos(thr)  + ! Banjo exit, 4K and 100K shields
     >    0.0031 / cos(thr + thwin) + ! 280 um exit widnow bowed inward
     >    0.0016 +                                    ! DC1
     >   (172. - th * 1.3 - 23.1 / cos(thr)) / 30420.  ! air
! extra 50 um of aluminum in 4K shield beyond 29 deg
        if(th.gt.29.) TAAL = TAAL + 0.00056
! extra 50 um of aluminum in 100K shield beyond 21 deg
        if(th.gt.21.) TAAL = TAAL + 0.00056
! 12 layers of superinsulation:
        if(th.gt.19.) TAAL = TAAL + 0.00029

        nprnt= nprnt+1
        if(nprnt.le.100) write(*,111) tbal,tbefor,tafter,taal
 111    format(1x,'EG1b rad len=',4f10.5)
      endif

      IF(INDEX(TARGET,'EG4').GT.0) THEN  ! Eg4 targets
        call eg4radlen(th,tbal,taal)
        nprnt= nprnt+1
        if(nprnt.le.100) write(*,122) tbal,tbefor,tafter,taal
 122    format(1x,'EG4 rad len=',4f10.5)
      endif

! New for E99-118
      IF(INDEX(TARGET,'E99S').GT.0) THEN  ! solid targets
        TBEFOR=t/cos(0.354)             ! target/cos(20.3 deg)
        TAFTER=(ttarg-t)/cos(thr-0.354) ! target(cos(th-20.3 deg)
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC                     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,112) tbal,tbefor,tafter,taal
 112    format(1x,'E99S rad len=',4f10.5)
      endif
      IF(INDEX(TARGET,'E99L').GT.0) THEN ! liquid in LD2 or LH2
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
        TBAL= TBEAM                     ! material before the target
        x_end = 4.095 * (ttarg-t)/ttarg * tan(thr)   ! position at endcap
        if(x_end.lt.2.0) then           ! went through endcap
          TAAL= TSPEC                     ! material after the target
        else
          taal=twall/sin(thr)+0.00606     ! went through side wall
        endif
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,113) tbal,tbefor,tafter,taal
 113    format(1x,'E99L rad len=',4f10.5)
      endif
      IF(INDEX(TARGET,'E99D').GT.0) THEN  ! Dummy aluminum 
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th): default
        x_end = 4.09 * tan(thr)         ! position at 2nd endcap from first
        if(x_end.gt.2.38) then             ! misses 2nd endcap +/- 2.38 cm   
          if(abs(t-ttarg/2.).lt.0.00001) then ! middle
            tafter = ttarg/2.0/2.0/cos(thr)  ! average
          else
            if(t.lt.ttarg/2.-0.00001) tafter=(ttarg/2.-t)/cos(thr)
            if(t.gt.ttarg/2.+0.00001) tafter=(ttarg-t)/cos(thr)
          endif
        endif
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC                     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,114) t/ttarg,tbal,tbefor,tafter,taal
 114    format(1x,'E99D rad len=',5f10.5)
      endif
      IF(INDEX(TARGET,'E99E').GT.0) THEN  ! for LH2, LD2 endcaps
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th): default
        x_end = 4.09 * tan(thr)         ! position at 2nd endcap from first
        if(x_end.gt.2.0) then           ! misses  2nd endcap
          if(abs(t-ttarg/2.).lt.0.0001) then ! middle
            tafter = ttarg/2.0/2.0/cos(thr)  ! average
          else
            if(t.lt.ttarg/2.-0.0001) tafter=(ttarg/2.-t)/cos(thr)
            if(t.gt.ttarg/2.+0.0001) tafter=(ttarg-t)/cos(thr)
          endif
        endif
        TBAL= TBEAM                      ! material before the target
        TAAL= TSPEC                      ! material after the target
! add liquid and walls to tbal
        if(t.lt.ttarg/2.-0.0001) then
          if(x_end.lt.2.0) then
            taal = taal + 0.0052/cos(thr)
          else
            taal = taal + (2.0/sin(thr)/890.+twall/sin(thr))
          endif
        endif
        if(t.gt.ttarg/2.+0.0001) then
          tbal = tbal + 0.0052
        endif
        if(abs(t-ttarg/2.).lt.0.0001) then ! middle
          tbal = tbal + 0.0052/2.0
          if(x_end.lt.2.0) then
            taal = taal + 0.0052/cos(thr)/2.0
          else
            taal = taal + (2.0/sin(thr)/890.+twall/sin(thr))/2.
          endif
        endif
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,115) t/ttarg,tbal,tbefor,tafter,taal
 115    format(1x,'E99E rad len=',5f10.5)
      endif
! For Tuna Can targets, make sure to use correct input values
! for ttarg, tbeam, and tspec. twall is not used.
      IF(INDEX(TARGET,'TUNA').GT.0) THEN ! liquid in LD2 or LH2
        TBEFOR=t                        ! target
c note: if change cell length, make sure ttarg matches for
c       actual density and material of target
        cell_diam = 3.96
        cell_wall = 0.0125
! z in cm relative to center 
        ztarg = cell_diam * (t / ttarg - 0.5) 
        call tunatarg(cell_diam, cell_wall, 
     >    ztarg, thr, cryo_cm, wall_cm)
        TAFTER = ttarg * cryo_cm / cell_diam
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC + wall_cm / 8.9     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,123) tbal,tbefor,tafter,taal,cryo_cm
 123    format(1x,'TUNA rad len=',5f10.5)
      endif

      IF(INDEX(TARGET,'CRYO17').GT.0) THEN     ! target geometry for
                                               ! 10 cm LH2 and LD2
        ECIR=1.315*2.54                        ! endcap inner radius (cm
        ECOR=1.320*2.54                        ! endcap outer radius (cm
        TARGLEN=3.942*2.54                     ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm
        
        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'SMALLRL').GT.0) THEN     ! Small radiation length
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T
        TAFTER=ttarg
        TAAL=TSPEC                               
        
      ENDIF

      IF(INDEX(TARGET,'A1ND2NPROD').GT.0) THEN     ! target geometry for
                                               ! pol 3He cell
        ECIR=(0.87/2)*2.54-0.15                ! endcap inner radius (cm
        ECOR=ECIR+0.015                        ! endcap outer radius (cm
        TARGLEN=40.0                           ! target length (cm
        c_1=1.38091E-5                         !up win thickness change
        c_2=0.110839                           !up win thickness change                                                
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=(t*targlen/ttarg)+(targlen/2.0)    ! t in cm
        if(tcm.lt.0) tcm=0
        if(tcm.gt.targlen) tcm=targlen

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! gas target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          rho_tar=(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2) 
     >    +(targlen-ECIR-tcm)*cos(thr))*sin(thr)!distance away tar center
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen !gas target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt((ECOR+c_1*exp(rho_tar/c_2))**2
     >    -((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    /7.04                                ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'A1ND2NREFN2').GT.0) THEN     ! target geometry for
                                               ! pol 3He cell
        ECIR=(0.87/2)*2.54-0.15                ! endcap inner radius (cm
        ECOR=ECIR+0.015                        ! endcap outer radius (cm
        TARGLEN=40.0                           ! target length (cm
        c_1=1.38091E-5                         !up win thickness change
        c_2=0.110839                           !up win thickness change                                                
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=(t*targlen/ttarg)+(targlen/2.0)    ! t in cm
        if(tcm.lt.0) tcm=0
        if(tcm.gt.targlen) tcm=targlen

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! gas target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          rho_tar=(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2) 
     >    +(targlen-ECIR-tcm)*cos(thr))*sin(thr)!distance away tar center
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen !gas target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt((ECOR+c_1*exp(rho_tar/c_2))**2
     >    -((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    /7.04                                ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'A1ND2NEMPTYUPWIN').GT.0) THEN     ! target geometry for
                                               ! pol 3He cell
        ECIR=(0.87/2)*2.54-0.15                ! endcap inner radius (cm
        ECOR=ECIR+0.015                        ! endcap outer radius (cm
        TARGLEN=0.015                          ! target length (cm
        c_1=1.38091E-5                         !up win thickness change
        c_2=0.110839                           !up win thickness change                                                
        ENTEC=40.0-ECIR                        ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=(t*targlen/ttarg)+(40.0/2.0)    ! t in cm
        if(tcm.lt.0) tcm=0
        if(tcm.gt.targlen) tcm=targlen

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=sin(PI-thr-asin((ECIR+targlen-tcm)*sin(thr)/ECIR))
     >    *(ECIR/sin(thr))*(ttarg/targlen)     ! glass target(SSA triangle)
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          rho_tar=(sqrt(ECOR**2-((ENTEC-tcm)*sin(thr))**2) 
     >    +(ENTEC-tcm)*cos(thr))*sin(thr)!distance away tar center
          TAFTER=                              ! e goes throught end cap
     >    sin(PI-thr-asin((ECIR+targlen-tcm)*sin(thr)/ECIR))
     >    *(ECIR/sin(thr))*(ttarg/targlen)     ! glass target(SSA triangle)
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt((ECOR+c_1*exp(rho_tar/c_2))**2
     >    -((ENTEC-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((ENTEC-tcm)*sin(thr))**2))
     >    /7.04                                ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'A1ND2NEMPTYDOWNWIN').GT.0) THEN     ! target geometry for
                                               ! pol 3He cell
        ECIR=(0.87/2)*2.54-0.15                ! endcap inner radius (cm
        ECOR=ECIR+0.015                        ! endcap outer radius (cm
        TARGLEN=0.015                          ! target length (cm
        c_1=1.38091E-5                         !up win thickness change
        c_2=0.110839                           !up win thickness change                                                
        ENTEC=40.0-ECIR                        ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T-40.0                          ! distance traveled bf
                                               ! vertex
        tcm=(t*targlen/ttarg)+(40.0/2.0)    ! t in cm
        if(tcm.lt.(40.0)) tcm=40.0
        if(tcm.gt.40.0+targlen) tcm=targlen+40.0

          rho_tar=ECOR*sin(thr-asin(sin(PI-thr)*(ECIR+tcm-40.0)/ECOR)) 
                                               !distance away tar center
          ang_b=asin((ECIR+tcm-40.0)*sin(PI-thr)/(ECOR+c_1*exp(rho_tar/c_2))) 
                                               !angle b for SSA triangle
          TAFTER=                              ! e goes throught end cap
     >    ECIR*sin(PI-ang_b)/sin(ang_b)/7.04   ! glass target(SSA triangle)
          TAAL=TSPEC                           ! material af the target

      ENDIF



      IF(INDEX(TARGET,'10OLD17').GT.0) THEN    ! target geometry for
                                               ! old 10 cm LH2 and LD2
        ECIR=0.795*2.54                        ! endcap inner radius (cm
        ECOR=0.800*2.54                        ! endcap outer radius (cm
        TARGLEN=10                             ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'4OLD17').GT.0) THEN     ! target geometry for
                                               ! old 4 cm LH2 and LD2
        ECIR=0.795*2.54                        ! endcap inner radius (cm
        ECOR=0.800*2.54                        ! endcap outer radius (cm
        TARGLEN=4                              ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      ATBX0 = AX0*TBEFOR+AX0A*TBAL                                      
      ATAX0 = AX0*TAFTER+AX0A*TAAL                                      
      BTBI  = B*TBEFOR+BA*TBAL                                          
      BTAI  = B*TAFTER+BA*TAAL                                          
C Redefine Initial and Final energies subtracTIng the most probable     
C ionization loss and calculate kinematic dependent quantities:         
                                                                        
      DEL01 = 0.                                                        
      DEL02 = 0.                                                        
      DELP1 = 0.                                                        
      DELP2 = 0.                                                        
      IF (TBEFOR.GT.0.) DEL01 = AX0*TBEFOR                              
     >        *(LOG(3.E9*AX0*TBEFOR*E0SET**2/EM**2/iZ**2)-0.5772)       
      IF (TBAL.GT.0.) DEL02 = AX0A*TBAL                                     
     >        *(LOG(3.E9*AX0A*TBAL*E0SET**2/EM**2/13.**2)-0.5772)       
      IF (TAFTER.GT.0.) DELP1 = AX0*TAFTER                              
     >        *(LOG(3.E9*AX0*TAFTER*EPSET**2/EM**2/iZ**2)-0.5772)       
      IF (TAAL.GT.0.) DELP2 = AX0A*TAAL                                 
     >        *(LOG(3.E9*AX0A*TAAL*EPSET**2/EM**2/13.**2)-0.5772)       
      DEL0  = DEL01+DEL02                                               
      DELP  = DELP1+DELP2                                               
                                                                        
ccc pyb TURNED OFF energy loss smearing!
      DEL0 = 0.
      DELP = 0.

1     E0 = E0SET-DEL0                                                   
      EP = EPSET+DELP                                                   
      Q2 = 4.*E0*EP*SIN(THR/2.)**2                                      
      R  = (PM+E0*(1-COS(THR)))/(PM-EP*(1-COS(THR)))                    
C Equivalent radiator for internal bremstrahlung correction             
      TEQUI = AL/PI*(LOG(Q2/EM**2)-1.)/B                                
                                                                        
      TB  = TBEFOR+TBAL+TEQUI                                           
      TA  = TAFTER+TAAL+TEQUI                                           
      BTB = B*(TBEFOR+TEQUI)+BA*TBAL                                    
      BTA = B*(TAFTER+TEQUI)+BA*TAAL                                    
c      write(6,'(''rad'',10f7.3)') b,tequi,tb,ta,btb,bta
      RETURN                                                            
      END                                                               
                                                                        
C=====================================================================  
                                                                        
      SUBROUTINE RADLEN42 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E142  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.       
!  Writes diagnostic and error messages to unit IU.                     
!                                                                       
!  Calculations are carried out in target coordinate system where the   
!  x axis is along the axis of the target and y axis perpendicular to   
!  x axis with positive values towards 8 GeV spectrometer.  The origin o
!  this coordinate system is at the in cap of the target cell.          
!  The out cap of the target cell is assumed to be ellipsoidal given by 
!  (B*(X-XC))**2 + (A*Y)**2 - (A*B)**2 = 0                              
!                                                                       
!  Changed r.l. of LH2 to 61.28 (1988 value) from 63.                   
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL LHE3/0./,LWALL/0./,LENDCAP/0./,XCM,TTARG,TWALL,TBEAM,TSPEC   
      REAL X,TR,TBAL,TAAL,XMAX/0./                                      
      REAL TAFTER         ! target material after scattering.           
      REAL RLENDCAP/0./   ! radiation length in end cap after scatter   
      REAL ECRL_F/0./, ECRL_R/0./ !rad length of end cap, front and rear
      INTEGER ITARG,NSEG,JTARG
      CHARACTER*7    TARGET                                             
      REAL      DAT(4,11)                                               
                                                                        
      REAL RADIUS/.825/                        
      LOGICAL FIRST/.TRUE./,FIRST1/.TRUE./                              
      LOGICAL INHE,INEC                                                 
      DATA    DAT/                                                      
! TARGPARAM: radiation lengths of various types of Al from P. 52        
! Info on wire arrayS: Jerry Davis says 5 mil diamter Al wire 25 mil spa
!              and there are two arrays.                                
! Info on Hymen: John Mark says it was 1 mil Al.                        
! Info on In-cap: from Hunter.  .01 to .012 inches total for both ends  
! Info on He3: from Hunter 2.32E20 atoms/cc.                            
!              30 cm long                                               
! Info on cell wall. Value of 4.4 mils is from p. 24                    
!              Actual value varies with x and is hardwired into program.
! Infor on endcap comes Hunter.                                         
! Info on insulation comes from J. Mark (10 layers of 1/4 mil mylar)    
!              it pretty much lines up with endcap of long targets,     
!              so for short targets need to take this into account.     
! Info on 1.6 gev window comes from p. 24 (sample)                      
! Thickness of 8 GeV exit window is 10.7 mils when not crinkled (Eislele
!                                                                       
!       Thickness     density      X0      Z    Name       Materl       
!         (cm)        (g/cm3)    (g/cm2)                                
     >    0.00400,    2.7,       24.0111, 13., 
     >    0.00254,    2.7,       24.0111, 13.,     
     >    0.01400,    2.52,      27.0000, 14.,       
     >    30.0280,    0.00116,   71.0000,  2.,         
     >    0.014,      2.52,      27.,     14.,       
     >    0.014,      2.52,      27.,     14.,       
     >    0.00635,    1.39,      39.95,   5.2,       
     >    0.00762,    2.68,      23.6311, 13.,     
     >    0.03048,    2.68,      23.6311, 13.,     
     >    16.0000,    0.001205,  36.975,  7.3,        
     >    0.02540,    2.68,      23.6311, 13./     
                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       DAT(2,4) = TTARG * DAT(3,4) / DAT(1,4) !Density = X0 *rl/len     
       DAT(1,5) = TWALL * DAT(3,5) / DAT(2,5) !len(cm)=  X0 *rl/den     
       XMAX = DAT(1,4) - DAT(1,3) - DAT(1,6)  !targ len minus end caps  
       ECRL_F = DAT(1,3) * DAT(2,3)/DAT(3,3)  !radiation length front ec
       ECRL_R = DAT(1,6) * DAT(2,6)/DAT(3,6)  !radiation length rear ec 
       WRITE( 6,'('' ****USING E142 TARGET MODEL '')')                  
      ENDIF                                                             
      INHE =.FALSE.                                                     
      IF(INDEX(TARGET,'E142F').GT.0) THEN  ! Gas in target goes thru wal
        INHE=.TRUE.                                                     
! Material after the scatter.                                           
! target model is a cylinder with flat end caps.                        
! Length in Helium                                                      
        LHE3 = RADIUS/SIN(TR)    ! in cm                                
        LWALL = DAT(1,5)/SIN(TR)                                        
        LENDCAP = 0.                                                    
      ELSEIF(INDEX(TARGET,'E142R').GT.0) THEN  ! exits thru endcap      
        INHE=.TRUE.                                                     
        XCM=X*DAT(3,4)/DAT(2,4)         ! changed to unit cm            
        LHE3 = (XMAX-XCM)/COS(TR)                                       
        LWALL =0.                                                       
        LENDCAP = DAT(1,6)/COS(TR)                                      
      ENDIF                                                             
      IF(INHE) THEN                                                     
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM + ! material before the target not including end ca 
     >        ECRL_F                       ! front end cap              
       TAAL= TSPEC +            !Non He3 material after scatter         
     >      LWALL  * DAT(2,5)/DAT(3,5) +                                
     >      LENDCAP* DAT(2,6)/DAT(3,6)                                  
       TAFTER = LHE3   * DAT(2,4)/DAT(3,4)  ! He3 after scatter         
      ENDIF                                                             
                                                                        
!End caps                                                               
      INEC=.FALSE.                                                      
      IF(INDEX(TARGET,'E142ECF').GT.0) THEN  !Front End Cap radiation co
        INEC=.TRUE.                                                     
        LHE3 = RADIUS/SIN(TR)                                           
        LWALL= DAT(1,5)/SIN(TR)                                         
        RLENDCAP = (ECRL_F -X)/COS(TR)    !front end cap after scatter  
        TBAL = TBEAM      !material before target                       
      ELSEIF(INDEX(TARGET,'E142ECR').GT.0) THEN  !rear end cap          
        INEC=.TRUE.                                                     
        LHE3 = 0.                                                       
        LWALL =0.                                                       
        RLENDCAP = (ECRL_F + ECRL_R -X)/COS(TR) !rear ec after scatter  
        TBAL = TBEAM +                                                  
     >        ECRL_F                       ! front end cap              
      ENDIF                                                             
                                                                        
      IF(INEC) THEN                                                     
       TAAL= TSPEC +            !material after target                  
     >      LWALL  * DAT(2,5)/DAT(3,5) +                                
     >      RLENDCAP                                                    
       TAFTER = LHE3   * DAT(2,4)/DAT(3,4)  ! He3 after scatter         
      ENDIF                                                             
                                                                        
      IF((.NOT.INEC) .AND. (.NOT.INHE))THEN !no target model            
       WRITE(6, '(''TARGET='',A8,'' NOT ALLOWED'')') TARGET             
       WRITE(10,'(''TARGET='',A8,'' NOT ALLOWED'')') TARGET             
       STOP                                                             
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C====================================================================== 
                                                                        
                                                                        
                                                                        
      SUBROUTINE RADLEN43 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E143  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the beginning of the target material.      
!  Steve Rock 4/9/93                                                    
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG ! total of amonia (nitrogen and H2 and He3 coolant     
      REAL TWALL,TBEAM,TSPEC                                            
      REAL X,TR,TBAL,TAAL                                         
      REAL TAFTER         ! target material after scattering.           
      INTEGER ITARG,NSEG,JTARG
      CHARACTER*7    TARGET 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                                             
      LOGICAL FIRST/.TRUE./
      CHARACTER*(*) STRING
      REAL TD                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       WRITE(6,'('' ENTERED RADLEN43'')')                               
      ENDIF                                                             
                                                                        
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM   ! material before the target amonia               
       TAAL= TSPEC    ! material after the target amonia.               
       TAFTER = TTARG -X   ! Amonia after scattering                    
       TD = 57.296*TR
! 10 deg spectrometer fit from test_targ2 6/2/98 with 302 cm of air in mags
       IF(TD.GE.8.3) THEN  
        IF((IA.EQ.14).OR.(IA.EQ.1)) THEN   ! N or p
         TAAL = -1.3171E-01 + 5.1404E-02*TD -5.7571E-03*TD**2 + 
     >     2.1474E-04*TD**3
        ELSEIF((IA.EQ.6).OR.(IA.EQ.2)) THEN  ! Li or d
         TAAL =  -1.3200E-01 + 5.1496E-02*TD -5.7681E-03*TD**2 + 
     >     2.1506E-04*TD**3
        ELSEIF((IA.EQ.9).or.(IA.EQ.12)) THEN   ! Long Be and long C
         TAAL =  -1.2472E-01  +4.9967E-02*TD -5.6091E-03*TD**2 +
     >     2.0969E-04*TD**3

        ENDIF
       ENDIF                                                                 
                                                                        
      RETURN                                                            
     
      ENTRY RADLEN43_INIT(STRING)
         STRING='** 10 deg T-Spec from fit in Targ_test2'
      RETURN

      END                                                               
                                                                        
                                                                        
C====================================================================== 
                                                                        
                                                                        
                                                                        
      SUBROUTINE RADLEN49 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E149  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the beginning of the target material.      
!  Steve Rock 4/26/93                                                   
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG ! total of amonia (nitrogen and H2 and He3 coolant     
      REAL TWALL,TBEAM,TSPEC                                            
      REAL X,TR,TBAL,TAAL
      REAL TAFTER         ! target material after scattering.           
      INTEGER ITARG,NSEG,JTARG                                                
      CHARACTER*7    TARGET                                             
!*** changed 4/8/03 pyb
      REAL RADIUS /3./                                                  
      LOGICAL FIRST/.TRUE./                                             
                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       WRITE(6,'('' ENTERED RADLEN49'')')                               
      ENDIF                                                             
                                                                        
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM   ! material before the target amonia               
       TAAL= TSPEC    ! material after the target amonia.               
       TAFTER = RADIUS/SIN(TR)/757. ! radiation lengths leaving target  
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                       
C=====================================================================  
                                                                        
      SUBROUTINE RADLEN154 (TARGET,X,THETAD,TAFTER,TBAL,TAAL,TBEFOR)
!_______________________________________________________________________
!  Calculates E154  radiation length (units of rl) for all materials   
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.                                             
      IMPLICIT NONE   
      CHARACTER*7    TARGET
      REAL X,THETAD,TBAL,TAAL !TH is scattering angle in degrees
      REAL TAFTER         ! target material after scattering.
      REAL TBEFOR         ! target material before scattering.  
      REAL TTOT,TLEN
      INTEGER ISEG154
      COMMON/E154/ISEG154,TTOT,TLEN    

C For Picard: from linda
       REAL*8 TB2(4)/0.00054D0, 0.00060D0, 0.00065D0, 0.00082D0/,
     >        TA2(4)/0.14631D0, 0.19431D0, 0.04434D0, 0.00067D0/,
     >        TLEN2(4)/3.D0, 5.D0, 2.D0, 20.D0/,
     >        TB5(4)/0.00062D0, 0.00075D0, 0.00081D0, 0.00090D0/,
     >        TA5(4)/0.07416D0, 0.09616D0, 0.04416D0, 0.00059D0/,
     >        TLEN5(4)/14.D0, 4.D0, 2.D0, 10.D0/

!use rho_ec =2.52, X0(rad length of glass)=28 gms/cm2 (guess from Wallet Cards
!end cap front = 69um =.00062 r.l.
!end cap rear  =61.6um .00055 r.l.
      REAL T_ECF/.00062/
      REAL T_ECR/.00055/

!assume all material is in end caps
      TAFTER=0.
      TBEFOR=0.
                            
         IF(THETAD.LT.3.5) THEN
           IF(INDEX(TARGET,'ECF').NE.0) THEN  ! front end cap
            TBAL= .5*TB2(1)
            TAAL = .5*TB2(1) +TA2(1)
            TLEN=1.
            TTOT=1.
           ELSEIF(INDEX(TARGET,'ECR').NE.0) THEN   !rear end cap
            TBAL =TB2(4) +.00014 +T_ECR/2.  !add last  10cm of 3He
            TAAL =TA2(4) -.00014 -T_ECR/2.
            TLEN=1.
            TTOT=1.
           ELSE
            TBAL = TB2(ISEG154)
            TAAL = TA2(ISEG154) 
            TLEN = TLEN2(ISEG154)
            TTOT =  TLEN2(1) +TLEN2(2) +TLEN2(3) +TLEN2(4) 
           ENDIF
         ELSEIF(THETAD.GT.4.0.AND.THETAD.LT.7.0) THEN
          IF(INDEX(TARGET,'ECF').NE.0) THEN  ! front end cap
            TBAL= .5*TB5(1) - .00010     !remove 7 cm of 3He
            TAAL = .5*TB5(1) +TA5(1) + .00010
            TLEN=1.
            TTOT=1.
           ELSEIF(INDEX(TARGET,'ECR').NE.0) THEN  !rear end cap
            TBAL =TB5(4) +.00007 +T_ECR/2.        !add last  5cm of 3He
            TAAL =TA5(4) -.00007 -T_ECR/2.
            TLEN=1.
            TTOT=1.
           ELSE
            TBAL = TB5(ISEG154)
            TAAL = TA5(ISEG154) 
            TLEN = TLEN5(ISEG154)
            TTOT =  TLEN5(1) +TLEN5(2) +TLEN5(3) +TLEN5(4) 
           ENDIF
         ENDIF

         RETURN
         END


                                                 
C=======================================================================
      SUBROUTINE RADLENGT2 (ITARG,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E140X radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TH (deg) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.       
!  Writes diagnostic and error messages to unit IU.                     
!                                                                       
!  Calculations are carried out in target coordinate system where the   
!  x axis is along the axis of the target and y axis perpendicular to   
!  x axis with positive values towards 8 GeV spectrometer.  The origin o
!  this coordinate system is at the in cap of the target cell.          
!  The out cap of the target cell is assumed to be ellipsoidal given by 
!  (B*(X-XC))**2 + (A*Y)**2 - (A*B)**2 = 0                              
!                                                                       
!  Changed r.l. of LH2 to 61.28 (1988 value) from 63.                   
!______________________________________________________________________ 
      COMMON /TRL/   TTARG,TWALL,TBEAM,TSPEC,KTARG,NSEG,JTARG!add 4/16/9
      REAL TTARG, TWALL, TBEAM, TSPEC                                   
      INTEGER KTARG !Dummy replacing ITARG which is in arguement        
      INTEGER JTARG ! 1 = internal rad only, 2= extern + intern         
      DIMENSION DAT(4,11)                                               
      DATA    DAT/                                                      
! TARGPARAM: radiation lengths of various types of Al from P. 52        
! Info on wire arrayS: Jerry Davis says 5 mil diamter Al wire 25 mil spa
!              and there are two arrays.                                
! Info on Hymen: John Mark says it was 1 mil Al.                        
! Info on In-cap is from p. 32 (sample in book, 3 mils thick)           
! Info on LH2: thickness is actually inner radius from 34, corrected    
!              for wall thickness and 0.4% shrinkage.                   
!              X0 value is from new 1988 particle                       
!              data book. (Used to be 63.3 in old books).               
! Values are overwritten if D2 is used                                  
! Info on cell wall. Value of 4.4 mils is from p. 24                    
!              Actual value varies with x and is hardwired into program.
! Infor on endcap comes from p. 24. Actuall thickness varies with y     
!              as is hardwired into program                             
! Info on insulation comes from J. Mark (10 layers of 1/4 mil mylar)    
!              it pretty much lines up with endcap of long targets,     
!              so for short targets need to take this into account.     
! Info on 1.6 gev window comes from p. 24 (sample)                      
! Thickness of 8 GeV exit window is 10.7 mils when not crinkled (Eislele
!                                                                       
!       Thickness     density      X0      Z    Name       Materl       
!         (cm)        (g/cm3)    (g/cm2)                                
     >    0.00400,    2.7,       24.0111, 13., ! Wire array  Al pure    
     >    0.00254,    2.7,       24.0111, 13., ! Hymen       Al pure    
     >    0.00762,    2.68,      23.6311, 13., ! In_cap      Al 5052    
     >    3.19600,    0.0707,    61.2800,  1., ! radius      LH2 or LD2 
     >    0.01143,    2.72,      23.6396, 13., ! Cell wall   Al 3004    
     >    0.01143,    2.72,      23.6396, 13., ! End_cap     Al 3004    
     >    0.00635,    1.39,      39.95,   5.2, ! Insulat     Mylar      
     >    0.00762,    2.68,      23.6311, 13., ! 1.6 window  Al 5052    
     >    0.03048,    2.68,      23.6311, 13., ! 8GEV windo  Al 5052    
     >    16.0000,    0.001205,  36.975,  7.3, ! Air /Vac gap VAC       
     >    0.02540,    2.68,      23.6311, 13./ ! quad_window Al 5052    
                                                                        
! Enter density and radiation length of liquid target                   
        IF(ITARG.EQ.1.OR.ITARG.EQ.2) THEN                               
          DAT(2,4)=0.0707                                               
          DAT(3,4)=61.28                                                
          IF(ITARG.EQ.1)  DAT(1,4)=15.836                               
          IF(ITARG.EQ.2)  DAT(1,4)=4.0470                               
        ELSE                                                            
          DAT(2,4)=0.1693                                               
          DAT(3,4)=122.60                                               
          IF(ITARG.EQ.4)  DAT(1,4)=15.763                               
          IF(ITARG.EQ.3)  DAT(1,4)=4.047                                
        ENDIF                                                           
        IF(ITARG.EQ.1) XMAX=15.818                                      
        IF(ITARG.EQ.2) XMAX=4.0292                                      
        IF(ITARG.EQ.3) XMAX=4.0292                                      
        IF(ITARG.EQ.4) XMAX=15.745                                      
                                                                        
                                                                        
        IF(ITARG.EQ.5.OR.ITARG.EQ.6) THEN                               
          TB=0.00045           ! wire array                             
     >      +0.00029           ! entrance window                        
            IF(ITARG.EQ.5) THEN ! 15 cm target                          
             IF(X.LE.0.0108)                                            
     >        TAFTER=(0.0108-X)/COS(TR)                                 
     >          +0.00011/COS(TR)                                        
             IF(X.GT.0.0108)                                            
     >        TAFTER=(0.0216-X)/COS(TR)                                 
     >        +0.00011+0.00011/COS(TR)                                  
            ELSE                                                        
              IF(TR.GT.0.671) THEN ! Doesn't hit downstream dummy       
               IF(X.LE.0.0108)                                          
     >          TAFTER=(0.0108-X)/COS(TR)                               
     >         +0.00011/COS(TR)                                         
               IF(X.GT.0.0108)                                          
     >          TAFTER=(0.0216-X)/COS(TR)                               
     >         +0.00011+0.00011/COS(TR)                                 
              ELSE                ! Does hit downstream dummy           
                TAFTER=(0.0216-X)/COS(TR) ! Target                      
     >         +0.00011+0.00011/COS(TR)                                 
              ENDIF                                                     
            ENDIF                                                       
            TA=  +0.00346      ! Chamber window                         
     >           +0.00053      ! air                                    
     >           +0.00288      ! spect window                           
          TAAL =TA                                                      
          TBAL =TB                                                      
          RETURN                                                        
          ENDIF                                                         
                                                                        
        IF(ITARG.EQ.7) THEN                                             
          TB=0.00045           ! wire array                             
     >      +0.00029           ! entrance window                        
          TAFTER=(TTARG-X)*COS(0.79)                                    
     >     /COS(ABS(0.79-TR))                                           
          TA=0.00022/SIN(TR)  !Mylar                                    
            TA=TA+0.00346      ! Chamber window                         
     >           +0.00053      ! air                                    
     >           +0.00288      ! spect window                           
          TAAL =TA                                                      
          TBAL =TB                                                      
          RETURN                                                        
          ENDIF                                                         
                                                                        
! Do liquild targets calculation                                        
                                                                        
      T1 = 0.0                                ! Clear sums              
      DO I = 1, 3                                                       
        TX = DAT(1,I)*DAT(2,I)/DAT(3,I)        ! Add wire array,        
        T1 = T1 + TX                           ! window, and in-cao     
      END DO                                                            
      TBAL=T1                                                           
      X=X*DAT(3,4)/DAT(2,4)         ! changed to unit cm                
! Find amount after scattering                                          
      TAAL=0.                                                           
! find horizontal shift for beam spot DY                                
      DY= 0.0                                                           
        DO I = 4,11                ! Do H2, wall, endcap, and mylar     
          IF(I.NE.8) THEN          ! exit window, air, and entr window  
            XEFF=DAT(1,I)          ! calculate geometry for first       
            IF(I.LE.7) CALL TARGET2(X,XMAX,TR,DY,I,ITARG,XEFF)          
            TX = XEFF*DAT(2,I)/DAT(3,I)                                 
            IF(I.EQ.4) TAFTER=TX                                        
            IF(I.GE.5) TAAL=TAAL+TX                                     
          ENDIF                                                         
        END DO                                                          
      X=X*DAT(2,4)/DAT(3,4)     ! changed back to r.l                   
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
      SUBROUTINE TARGET2 (X,XMAX,TR,DY,I,ITARG,XEFF)                    
!______________________________________________________________________ 
!  Calculates a thickness of material traversed by a scattered electron 
!  in a target contained in a target cell whose geometry is considered  
!  here.  Experiment = NE11. Only does 8 GeV side.                      
!  Input variables:                                                     
!  X    - x coordinate of scattering center,                            
!  XMAX - maximum x coordinate of LH2 = reference point (rear endcap)   
!  TR   - scattering angle alpha in rad.,                               
!  DY   - lateral beam offset in cm.  (+ = towards spectrometer)        
!  I    - indexes Ith element in the target cell.                       
!  ITARG - 1 or 4 for long cell, 2  or 3 for short cell                 
!  XEFF - thickness trversed in cm (0. if not in path)                  
!       - as input, it is nominal perpendicular thickness               
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      REAL*4  A /1.300/               !Ellipse x semiaxis (cm)          
      REAL*4  B /3.20000/             !Ellipse y semiaxis (cm)          
      REAL*4  PI/3.14159265/                                            
      REAL X,XMAX,TR,DY,XNOM,XEFF,XC,X1,X2,Y1,Y2,XWALL      
      INTEGER ITARG,I                                                   
      LOGICAL HIT                                                       
                                                                        
      XNOM=XEFF                       ! Save perpendicular thickness    
      XEFF=0.                         ! Set to 0 in case don't go thru  
      XC=XMAX-A                       ! x of ellipse center             
      IF(X .LT. 0.0 .OR. X .GT. XMAX) RETURN  ! Should never happen     
      CALL ENDCAP2(X,DY,XC,A,B,TR,X1,Y1,HIT) ! where his endcap or wal  
      IF(I.EQ.4) XEFF=SQRT((X1-X)**2+(Y1-DY)**2) ! Hydrogen             
      IF(HIT) THEN                                                      
        IF(I.EQ.6) THEN                ! Hit the endcap                 
          IF(Y1.LT.2.60) XNOM=0.0112   ! Put in y-dependence of endcap  
          IF(Y1.GE.2.60.AND.Y1.LT.2.9) XNOM=0.0178  ! based on Peter's m
          IF(Y1.GE.2.90) XNOM=0.0254   ! of an actual endcap.           
          CALL ENDCAP2(X,DY,XC,(A+XNOM),(B+XNOM),TR,X2,Y2,HIT)          
          IF(.NOT.HIT) THEN           ! Didn't hit other side of cap!   
            WRITE(6,'(1X,''ERROR IN PASSING ENDCAP'')')                 
          ELSE                                                          
            XEFF=SQRT((X2-X1)**2+(Y2-Y1)**2)                            
          ENDIF                                                         
        ENDIF                         ! If hit endcap, never hits wall b
        IF(I.EQ.7.AND.(ITARG.EQ.2.OR.ITARG.EQ.3)) THEN ! short target, c
          IF(X+(3.2-DY)*COS(TR)/SIN(TR).LT.9.5) THEN ! thru insultation 
            XEFF=XNOM/SIN(TR)  ! which ends at endcap of long target    
          ENDIF                ! (9.5 cm from begining of short target) 
        ENDIF                  ! where radius is 3.2 cm TARGPARAM:      
      ELSE                            ! Do wall and inslut.             
        IF(I.EQ.5) THEN        ! Put in rough x dependence of thickness 
          XWALL=XMAX-1.3-X     ! TARGPARAM: 1.3 cm is where wall begins 
          XWALL=XWALL-(3.2-DY)*COS(TR)/SIN(TR) ! Go from middle to edg  
          IF(XWALL.LT.0.50) XNOM=0.0254   ! based on Peter's meas.      
          IF(XWALL.GE.0.50.AND.XWALL.LT.1.00) XNOM=0.0238               
          IF(XWALL.GE.1.00.AND.XWALL.LT.1.50) XNOM=0.0152               
          IF(XWALL.GE.1.50.AND.XWALL.LT.2.50) XNOM=0.0127               
          IF(XWALL.GT.2.50) XNOM=0.0114                                 
        ENDIF                                                           
        IF(I.EQ.5.OR.I.EQ.7) XEFF=XNOM/SIN(TR)  ! Goes at angle TR      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
      SUBROUTINE  ENDCAP2 (X0,Y0,XC,A,B,TR,X,Y,HIT)                     
!______________________________________________________________________ 
!  Calculates an intersection point between a straight line and an      
!  ellipse (end cap of a target cell).  The line is given by:           
!  y = Y0 + (x - X0)*tan(TR),  and the ellipse by:                      
!  (B*(x - XC))**2 + (A*y)**2 = (A*B)**2.                               
!  The center of the target is taken as the y axis.                     
!  It set HIT .TRUE. if the solution exists, .FALSE. otherwise.         
!  If doesn't hit endcap, then give x,y of where hits target wall,      
!  assumed to be at radius B                                            
!  Output variables:                                                    
!  X,Y  - are the coordinates of the intersection.                      
!                                                                       
!  Assumes all x dimensions measured from in-cap of target, and         
!  positive Y is towards the 8 gev spectrometer (where electron goes).  
!  Note that the scattered electron can exit ONLY through that portion o
!  the ellipse for which XC<=x<=XC+A and y=>DY                          
!                                                                       
!  ZMS                  Jan. 1985...                                    
!  AFS                  2 May 85  simplified routine.                   
!  pyb  checked routine, added implicit none, provision for 90 deg      
!       and correct treatment of Y0                                     
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      REAL X0,Y0,XC,A,B,TR,X,Y                                          
      REAL TA,G,H,AA,BB,CC,DD                                           
      LOGICAL HIT                                                       
      REAL*4  PI/ 3.14159265/                                           
                                                                        
      HIT = .FALSE.                                                     
      X = 0.                                                            
      Y = 0.                                                            
                                                                        
! If TR very close to 90 degress, do special solution                   
      IF(ABS(TR-PI/2.).LT.0.017) THEN                                   
        IF(X0.GT.XC) THEN      ! Is it beyond center of endcap          
          HIT=.TRUE.                                                    
          X=X0                                                          
          Y=B*SQRT(1.0-(X0-XC)**2/A**2)                                 
        ELSE                   ! Must be upstream of endcap             
          X=X0                                                          
          Y=B                                                           
        ENDIF                                                           
        RETURN                                                          
      ENDIF                                                             
                                                                        
! Solve quadratic equation to see if hits endcap                        
      TA = TAN (TR)                                                     
      G  = A*TA                                                         
      H  = (X0 - XC)*TA - Y0                                            
      AA = 1.0 + G*G/(B*B)                                              
      BB = 2.0 * H                                                      
      CC = H*H - G*G                                                    
      DD = BB*BB - 4.0*AA*CC                                            
      IF(DD .GE. 0.0) THEN            ! See if radical postive          
        Y = 0.5*(-BB + SQRT (DD))/AA                                    
        X = X0 + (Y-Y0)/TA                                              
      END IF                                                            
      IF(X .GE. XC) THEN   ! Is solution in endcap region?              
        HIT = .TRUE.                                                    
      ELSE                 ! Otherwise, make it hit the wall            
        X=X0+(B-Y0)/TA                                                  
        Y=B                                                             
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C==============================================================         
                                                                        
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE CROSS(E0,EP,TH,SIGMA)                                  
                                                                        
C calculates cross section per NUCLEUS, not per nucleon.                
      IMPLICIT NONE 
      REAL E0,EP,TH,SIGMA
      COMMON /CON/ PI,PM,PP,EM,AL   
      REAL PI,PM,PP,EM,AL,W1,W2
      COMMON /SIG/ CSTYPE                                               
      CHARACTER*1  CSTYPE                                               
                                                                        
C inelastic tail case         = 'I'                                     
C quasielastic smeared case   = 'Q'                                     
C quasielastic unsmeared case = 'P'                                     
C elastic case                = 'E'                                     
                                                                        
      IF (CSTYPE.EQ.'I') CALL  SECNUCLW(E0,EP,TH,SIGMA)                  
      IF (CSTYPE.EQ.'Q') CALL  QUASIY8(E0,EP,TH,SIGMA)                  
C     IF (CSTYPE.EQ.'Q') CALL  FYXSEC8(E0,EP,TH,SIGMA)                 
      IF (CSTYPE.EQ.'P') CALL  QELASTIC(E0,TH,SIGMA,W1,W2)                  
      IF (CSTYPE.EQ.'E') CALL  ELASTIC(E0,EP,TH,SIGMA)                  
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4   FUNCTION FTSAI(E0,EP,TH)                                          
      IMPLICIT NONE
      REAL E0,EP,TH
                                                                  
      COMMON /CON/ PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL 
      REAL ALPI,Q2,VAC,VERTEX,ARG,SMALL,CORR,DELVAC,DSPEN
      PARAMETER    (ALPI=2.32281E-03)                                   
                                                                        
      Q2     = 4.*E0*EP*SIN(TH/2.)**2                                   
      VAC    = 2.*DELVAC(Q2)                                            
      VERTEX = 2.*ALPI*(-1.+0.75*ALOG(Q2/EM**2))                        
      ARG    = COS(TH*.5)**2                                            
      SMALL  = ALPI*(PI**2/6.-DSPEN(ARG))                               
      CORR   = -ALPI/2.*LOG(E0/EP)**2                                   
      FTSAI  = 1.+VAC+VERTEX+SMALL+CORR                                 
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION BREMS(E,EPS,T,ID)                                        
      IMPLICIT NONE
      REAL E,EPS,T,tmp
      INTEGER ID
C Bremstrahlung and ionization energy loss probablity function. E is the
C incoming energy of particle passing through T radiation lengths to    
C lose EPS energy. Formula from Tsai; Rev of Mod Physics 46,(1974)815,  
C equation 3.83. Complete screening.                                    
                                                                        
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL   
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                     
      COMMON /CON/ PI,PM,PP,EM,AL   
      REAL PI,PM,PP,EM,AL 
      REAL Y,WB,WI                                    
      REAL ATX0/0/, BT/0./  !initialization to make vm compiler happy   
                                                                        
      IF (T.EQ.0.) THEN                                                 
        BREMS = 0.                                                      
        RETURN                                                          
      ENDIF                                                             
                                                                        
      IF (ID.EQ.1) THEN                                                 
        ATX0 = ATBX0                                                    
        BT   = BTBI                                                     
      ELSE IF (ID.EQ.2) THEN                                            
        ATX0 = ATAX0                                                    
        BT   = BTAI                                                     
      ENDIF                                                             
                                                                        
C BREMS is Tsai eq.(3.83) normalized to be 1 at Y=0.  The higher order  
C order terms, with 2 approximations, result in the "-f" term of eq.    
C (3.66):                                                               
      Y     = EPS/E                                                     
      if(y.le.0.) then
c        write(6,'(1x,''ERROR, eps,e,atx0,bt='',4f10.3)')
c     >    eps,e,atx0,bt
        brems=0.
        return
      endif

      tmp = (1-Y+.75*Y**2)                                            
      WB    = B*T/EPS*tmp                                             
                                                                        
C Ionization:                                                           
      WI    = ATX0/MAX(EPS,10.*ATX0)**2*(1.+EPS**2/E/(E-EPS))**2        
      BREMS = (1.+0.5772*BT-0.62*BT**2) * Y**(B*T) * (WB+WI)                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION DSPEN(X)               

      IMPLICIT NONE                              
      REAL*4 X,FSPENS
      REAL  F1 /1.644934/                                       
                                                                        
      IF (X.LT.-1.) THEN                                                
           DSPEN = -.5*ALOG(1.-X)*ALOG(X**2/(1.-X))-F1+FSPENS(1./(1.-X))
      ELSE IF (X.LT.0.) THEN                                            
           DSPEN = -.5*ALOG(1.-X)**2-FSPENS(X/(X-1.))                   
      ELSE IF (X.LE..5) THEN                                            
           DSPEN = FSPENS(X)                                            
      ELSE IF (X.LE.1.) THEN                                            
           DSPEN = F1-ALOG(X)*ALOG(1.-X+1.E-10)-FSPENS(1.-X)            
      ELSE IF (X.LE.2.) THEN                                            
           DSPEN = F1-.5*ALOG(X)*ALOG((X-1.)**2/X)+FSPENS(1.-1./X)      
      ELSE                                                              
           DSPEN = 2*F1-.5*ALOG(X)**2-FSPENS(1./X)                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
      REAL*4 FUNCTION FSPENS(X)
      IMPLICIT NONE
      REAL*4 X,A,F,AN,TCH,B

      A   = 1.                                                          
      F   = .0                                                          
      AN  = .0                                                          
      TCH = 1.E-16                                                      
1     AN  = AN+1.                                                       
      A   = A*X                                                         
      B   = A/AN**2                                                     
      F   = F+B                                                         
      IF (B.GT.TCH) GO TO 1                                             
      FSPENS = F                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION DELVAC(T)                                                
                                                                        
C VACUUM POLARIZATION FOR LEPTONS AND HADRONS, includes terms for       
C electrons, muons, taus, and quarks.  The expression for small mass    
C reduces to Tsai's expression (2.10) in the slac pub.                  
      IMPLICIT NONE
      REAL*4 T, AMF2(9), COL(9), RTMF2,RT,ALT,AI34,SUML,A,B,C,SUMH
      INTEGER IC
      REAL*4  EPST / 1.E-3 /                                         
      REAL*4     ALPH_PI /.0023229 /                                    
      DATA       AMF2 / .26112E-6,.011164,3.1827,.0064,.0064,2.25,      
     +                  .09,1024.,20.25 /                               
      DATA       COL  / 3*1.,6*3. /                                     
                                                                        
C For newer treatment of vacuum correction, comment out the following   
C three lines.  These lines are exactly (2.10) of Tsai's Slac Pub '74:  
                                                                        
C     EM     = .000511                                                  
C     DELVAC = ALPH_PI*(-5./9.+ALOG(T/EM**2)/3.)                        
C     RETURN                                                            
                                                                        
      SUML = 0.                                                         
      DO 121 IC=1,3                                                     
           RTMF2 = T/AMF2(IC)                                           
           IF (RTMF2.LT.EPST) THEN                                      
                AI34 = RTMF2*(2./15.-1./70.*RTMF2)                      
           ELSE                                                         
                RT   = SQRT(1.+4./RTMF2)                                
                ALT  = RT*ALOG(RTMF2*(1.+RT)**2/4.)                     
                AI34 = -5./9.+4./3./RTMF2+(1./3.-2./3./RTMF2)*ALT       
           ENDIF                                                        
           SUML = SUML+COL(IC)*AI34                                     
121   CONTINUE                                                          
                                                                        
C HADRONIC VACUUM POLARIZATION TAKEN FROM BURKHARD, TASSO NOTE 192, 1982
                                                                        
      IF (T.LT.1.) THEN                                                 
           A = -1.345E-9                                                
           B = -2.302E-3                                                
           C =  4.091                                                   
      ELSE IF (T.LE.64.) THEN                                           
           A = -1.512E-3                                                
           B = -2.822E-3                                                
           C =  1.218                                                   
      ELSE                                                              
           A = -1.1344E-3                                               
           B = -3.0680E-3                                               
           C =  9.9992E-1                                               
      ENDIF                                                             
                                                                        
      SUMH   = -(A+B*ALOG(1.+C*T))                                      
      DELVAC = SUML*ALPH_PI+SUMH                                        
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE SIMP(A,B,N,FUNC,ANS)                                   
                                                                        
      REAL LOW                                                          
                                                                        
C Subroutine to do numerical integration using Simpson's rule           
C Input A lower limit                                                   
C Input B upper limit                                                   
C Input N number of bins; FUNC is called N+1 times; N HAS TO BE EVEN    
C Input FUNC Function subroutine address                                
C Output ANS estimate of the definite integral                          
                                                                        
      ANS = 0.                                                          
      IF (A.GE.B) RETURN                                                
      IF (N.LE.0) RETURN                                                
      IF (N.EQ.1) THEN                                                  
           ANS = FUNC((A+B)/2.)*(B-A)                                   
           RETURN                                                       
      ENDIF                                                             
      BINW = (B-A)/FLOAT(N)                                             
                                                                        
C FUNCTION AT LOWER LIMIT                                               
                                                                        
      LOW = FUNC(A)                                                     
                                                                        
C FUNCTION AT ODD BINS                                                  
                                                                        
      ODD = 0.                                                          
      DO 100 I = 1,N-1,2                                                
           ODD = ODD+FUNC(A+BINW*FLOAT(I))                              
100   CONTINUE                                                          
                                                                        
C FUNCTION AT EVEN BINS                                                 
                                                                        
      EVEN = 0.                                                         
      DO 200 I = 2,N-2,2                                                
           EVEN = EVEN+FUNC(A+BINW*FLOAT(I))                            
200   CONTINUE                                                          
                                                                        
C FUNCTION AT UPPER LIMIT                                               
                                                                        
      UP = FUNC(B)                                                      
                                                                        
C ANSWER                                                                
                                                                        
      ANS = (LOW+4.*ODD+2.*EVEN+UP)*BINW/3.                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE SIMP2(A,B,N,FUNC,ANS)                                  
                                                                        
      REAL LOW                                                          
      ANS = 0.                                                          
      IF (A.GE.B) RETURN                                                
      IF (N.LE.0) RETURN                                                
      IF (N.EQ.1) THEN                                                  
           ANS = FUNC((A+B)/2.)*(B-A)                                   
           RETURN                                                       
      ENDIF                                                             
      BINW = (B-A)/FLOAT(N)                                             
      LOW = FUNC(A)                                                     
      ODD = 0.                                                          
      DO 100 I = 1,N-1,2                                                
           ODD = ODD+FUNC(A+BINW*FLOAT(I))                              
100   CONTINUE                                                          
      EVEN = 0.                                                         
      DO 200 I = 2,N-2,2                                                
           EVEN = EVEN+FUNC(A+BINW*FLOAT(I))                            
200   CONTINUE                                                          
      UP = FUNC(B)                                                      
      ANS = (LOW+4.*ODD+2.*EVEN+UP)*BINW/3.                             
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE ELASTIC(E0,EP,TH,SIGMA)                                
                                                                        
C subroutine to calculate electron-nucleus electron scattering cross    
C section. Formulas from ASTRUC of Mo-Tsai program.                     
      IMPLICIT NONE
      REAL E0,EP,TH,SIGMA       
      REAL CHBAR,ALPHA,PM,PI,FDEL,FSHELL,FGAUSS
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828)                                       
      PARAMETER (PI    = 3.1415927)                                     

      COMMON    /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      REAL THR,QSQ,TAU
      REAL W1,W2,FF,CMOTT,CSMOTT,RECOIL
      REAL A_ULMAR,B_ULMAR,w2_old
      REAL*8 GE,GM,GEP,GEN,GMP,GMN
                                                                  
      sigma = 0.
c removed 3/25/08
c      if(iz.lt.1) return

      THR = TH*PI/180.                                                  
      QSQ = 4.*E0*EP*SIN(THR/2.)**2 
      CALL NUC_FORM_FACTOR(QSQ,W1,W2,FF)                        

c THIS SKIPS REST OF CODE!!! WHY???
      IF(1.EQ.1) GO TO 1111
      TAU = QSQ/4./PM**2                                                
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                                   
                                                                        
      IF (iA.EQ.1) THEN                                                 
           W1 = TAU*GMP**2                                              
           W2 = (GEP**2+W1)/(1.+TAU)                                    
      ELSEIF (iA.EQ.2) THEN  
        IF((IDUT.GE.11).AND.(IDUT.LE.14).AND.(QSQ.LE.3.5))THEN   !Tjon fits
           !Ulmar d2 elastic mode                                       
           CALL DEUT_U1(IDUT-10,QSQ,A_ULMAR,B_ULMAR)                     
           W1 = B_ULMAR/2.  ! sigma = sig_mot(A + B*tan...)
           W2 = A_ULMAR 
        ELSEIF(IDUT.EQ.1) THEN ! Linda Stuart's Model installed 5/30/96
           CALL FFD(DBLE(QSQ),GE,GM)   
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU) 
        ELSE   ! old  elastic deuterium from original code   
           FF  = FDEL(QSQ)                                              
           W1  = FF**2*TAU*.6667*(GMN+GMP)**2                           
           W2  = W1+(FF*(avgN*(GEN+TAU*GMN)+GEP+TAU*GMP)/(1.+TAU))**2   
        ENDIF
      ELSEIF (iA.EQ.3) THEN  !3HE  ! added 5/30/96  SER
        CALL FFHE3(DBLE(QSQ),GE,GM)
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU)  
           W2_old  = (iz*FSHELL(QSQ) )**2
      ELSEIF (iA.LE.20) THEN
           FF  = FSHELL(QSQ)                                            
           W1  = 0.                                                     
           W2  = (iZ*FF)**2                                             
      ELSE                                                              
           FF  = FGAUSS(QSQ)                                            
           W1  = 0.                                                     
           W2  = (iZ*FF)**2                                             
      ENDIF                                                             
 1111 CONTINUE

      CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
      CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
      RECOIL = avgM/(avgM+E0*(1.-COS(THR)))                             
      SIGMA  = (W2+2.*W1*TAN(THR/2.)**2)*CSMOTT*RECOIL                  
!Change elastic cross section by asymetry if doing asymetry run on prot 
c     IF((                                                              
c    >    ( (   (INDEX(TARGET,'E142')+INDEX(TARGET,'E143') ).GT.0)      
c    >           .AND.(IA.EQ.1)                                         
c    >    )   .OR.                                                      
c    >    ( (INDEX(TARGET,'E149') .GT.0 ).AND.(IA.EQ.2)!fake el d2 asym 
c    >    )                                                             
c    >   )                                                              
c    >      .AND.( INDEX(TARGET,'_P').GT.0)                             
c    >   )                                                              
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >               CSMOTT*RECOIL,SIGMA)                               
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE QELASTIC(E0,TH,SIGMA,W1,W2)                               
                                                                        
C Subroutine to calculate the quasi elastic cross section using NFORM.  
      IMPLICIT NONE
      REAL E0,TH,SIGMA                             
      REAL CHBAR,ALPHA,PM,PM24,PI
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828, PM24=3.5216)                                       
      PARAMETER (PI    = 3.1415927)                                     
      REAL FF/1./  !phony initialization to please compiler: SER 4/15/93
      COMMON/TARGT/ iZ,iA,avgN,avgA,avgM,amuM                       
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL THR,QSQ,TAU,W1,W2,CMOTT,CSMOTT,RECOIL,Pauli_sup1,PAULI_SUP2
      REAL EP,SINSQ,RAD,THR2
      REAL*8 GEP,GEN,GMP,GMN
      LOGICAL FIRST/.TRUE./

      IF(FIRST) THEN
         CMOTT  = CHBAR**2*0.001*ALPHA**2/4.  
         RAD = PI/180.
         FIRST=.FALSE.
      ENDIF
                                                                        
      SIGMA = 0.                                                        
      IF (iA.LE.1) RETURN                                               
                                                                        
      THR = TH *RAD
      THR2 = THR/2.
      SINSQ =  SIN(THR2)**2
      EP    = PM*E0/(PM+2.*E0*SINSQ)                                      
      QSQ = 4.*E0*EP*SINSQ                                     
      TAU = QSQ/PM24                                               
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)  

! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)

      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU) 
      CSMOTT = CMOTT*(1.-SINSQ)/(E0*SINSQ)**2              
      RECOIL = PM/(PM+E0*(1.-COS(THR)))                                 
      SIGMA  = (W2+2.0*TAN(THR2)**2*W1)*CSMOTT*RECOIL                 
! Are we changing quasi-elastic cross section due to Asymmetry?         
c     IF( ( (INDEX(TARGET,'E142')+INDEX(TARGET,'E143') +                
c    >       INDEX(TARGET,'E149') ).GT.0 )                              
c    >     .AND.(INDEX(TARGET,'_P').GT.0) )                             
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >  CSMOTT*RECOIL,SIGMA)                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE QUASIY8(E0,EP,TH,SIGMA)                                
                                                                        
C Calculates quasielastic cross section 
c using Donnelly/Sick super scaling model
c but Paris w.f. for deuteron
c for IA=1, does simple resolution smearing only
c changed to call F1F2QE 7/08 pyb

      IMPLICIT NONE     
      REAL E0,EP,TH,SIGMA
      REAL  CHBAR,ALPHA,PM,PI                                            
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828)                                       
      PARAMETER (PI    = 3.1415927)                                     
      REAL FF/1./  !phony initialization to please compiler: SER 4/15/93
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM 
      COMMON     /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      logical smrelasp, usearen
      real smrdEp,dep, dW
      common/smrproton/ smrelasp,smrdEp,usearen
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL PAULI_SUP1,PAULI_SUP2,FFF
      REAL*8 GEP,GEN,GMP,GMN
      real*8 qv,a,b,c,amf2,ams,disc,y,dydep,en,emp,denominator
      REAL THR,SINSQ,COSTH,QSQ,TAU,W1,W2,CMOTT,CSMOTT,RECOIL,DSIGE,FY,SD
      real*8 alfa,c1,c2,ARENHOVEL_SIG
      real*8 kappa,lam,lamp,taup,squigglef,psi,psip,nuL,nut
      real*8 kf,es,GM2bar,GE2bar,W1bar,W2bar,Delta,GL,GT
      common/testing/prttst,usegrd
      integer i,model,in_range
      logical prttst,usegrd
      real*8 z8,a8,q28,w28,f18,f28

! 100*prob(Pz) from Paris wave function in 10 MeV bins from -1 to 1 GeV
! integeral is 1. Integral is 100.
      integer izz
      real*8 depdpz,pz,qvp
! using just w*w for d-state (tail too big)
      real*8 fydr(200)/
     > 0.00000,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00000/
!*** NO D-sate (tail too small)
       real*8 fydn(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00030,
     > 0.00033,0.00038,0.00041,0.00047,0.00050,0.00057,0.00061,0.00068,
     > 0.00073,0.00081,0.00087,0.00095,0.00102,0.00111,0.00119,0.00128,
     > 0.00137,0.00147,0.00156,0.00167,0.00176,0.00189,0.00198,0.00210,
     > 0.00220,0.00232,0.00242,0.00254,0.00264,0.00276,0.00284,0.00296,
     > 0.00303,0.00313,0.00319,0.00327,0.00331,0.00338,0.00340,0.00344,
     > 0.00344,0.00347,0.00345,0.00348,0.00346,0.00349,0.00350,0.00358,
     > 0.00366,0.00386,0.00412,0.00454,0.00510,0.00594,0.00704,0.00859,
     > 0.01060,0.01334,0.01688,0.02156,0.02763,0.03544,0.04561,0.05878,
     > 0.07561,0.09741,0.12537,0.16167,0.20811,0.26892,0.34765,0.45096,
     > 0.58609,0.76539,1.00428,1.32177,1.75207,2.33167,3.11006,4.13194,
     > 5.42556,6.92172,8.36256,9.28786,9.28786,8.36256,6.92172,5.42556,
     > 4.13194,3.11006,2.33167,1.75207,1.32177,1.00428,0.76539,0.58609,
     > 0.45096,0.34765,0.26892,0.20811,0.16167,0.12537,0.09741,0.07561,
     > 0.05878,0.04561,0.03544,0.02763,0.02156,0.01688,0.01334,0.01060,
     > 0.00859,0.00704,0.00594,0.00510,0.00454,0.00412,0.00386,0.00366,
     > 0.00358,0.00350,0.00349,0.00346,0.00348,0.00345,0.00347,0.00344,
     > 0.00344,0.00340,0.00338,0.00331,0.00327,0.00319,0.00313,0.00303,
     > 0.00296,0.00284,0.00276,0.00264,0.00254,0.00242,0.00232,0.00220,
     > 0.00210,0.00198,0.00189,0.00176,0.00167,0.00156,0.00147,0.00137,
     > 0.00128,0.00119,0.00111,0.00102,0.00095,0.00087,0.00081,0.00073,
     > 0.00068,0.00061,0.00057,0.00050,0.00047,0.00041,0.00038,0.00033,
     > 0.00030,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/
! using 1.5 * (1-cos**2)**2 for dtate (better than either of above)
       real*8 fyd(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00031,
     > 0.00033,0.00038,0.00042,0.00048,0.00052,0.00058,0.00063,0.00070,
     > 0.00076,0.00084,0.00090,0.00099,0.00107,0.00117,0.00125,0.00136,
     > 0.00145,0.00158,0.00168,0.00181,0.00192,0.00207,0.00219,0.00235,
     > 0.00249,0.00265,0.00280,0.00297,0.00313,0.00332,0.00347,0.00368,
     > 0.00385,0.00406,0.00424,0.00447,0.00467,0.00491,0.00513,0.00541,
     > 0.00565,0.00597,0.00626,0.00665,0.00703,0.00751,0.00803,0.00868,
     > 0.00938,0.01030,0.01135,0.01265,0.01421,0.01616,0.01848,0.02143,
     > 0.02496,0.02942,0.03489,0.04170,0.05015,0.06060,0.07374,0.09018,
     > 0.11064,0.13651,0.16893,0.21017,0.26205,0.32886,0.41420,0.52472,
     > 0.66768,0.85544,1.10347,1.43058,1.87107,2.46124,3.25006,4.28232,
     > 5.58540,7.08967,8.53679,9.46493,9.46493,8.53679,7.08967,5.58540,
     > 4.28232,3.25006,2.46124,1.87107,1.43058,1.10347,0.85544,0.66768,
     > 0.52472,0.41420,0.32886,0.26205,0.21017,0.16893,0.13651,0.11064,
     > 0.09018,0.07374,0.06060,0.05015,0.04170,0.03489,0.02942,0.02496,
     > 0.02143,0.01848,0.01616,0.01421,0.01265,0.01135,0.01030,0.00938,
     > 0.00868,0.00803,0.00751,0.00703,0.00665,0.00626,0.00597,0.00565,
     > 0.00541,0.00513,0.00491,0.00467,0.00447,0.00424,0.00406,0.00385,
     > 0.00368,0.00347,0.00332,0.00313,0.00297,0.00280,0.00265,0.00249,
     > 0.00235,0.00219,0.00207,0.00192,0.00181,0.00168,0.00158,0.00145,
     > 0.00136,0.00125,0.00117,0.00107,0.00099,0.00090,0.00084,0.00076,
     > 0.00070,0.00063,0.00058,0.00052,0.00048,0.00042,0.00038,0.00033,
     > 0.00031,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/
      logical doing_elas,usef1f2qe
      common/doelas/ doing_elas

      SIGMA = 0.                                                        
      IF (iA.LE.1 .and. (.not.smrelasp)) RETURN
      IF (iA.LE.1 .and. (.not.doing_elas)) RETURN

      if(IA.eq.2 .and. usearen .and. (.not.doing_elas)) then
cc        model = 1
        model = 2
        sigma = ARENHOVEL_SIG(E0,EP,TH,MODEL,IN_RANGE)
        return
      endif

c xxx changed to use new routine unless arenhovel specified for D2
      usef1f2qe = .true.
      if(usef1f2qe .and. (.not.doing_elas)) then
       THR   = TH*PI/180.                                                
       SINSQ = SIN(THR/2.)**2                                            
       COSTH = COS(THR)                                                  
       QSQ   = 4.*E0*EP*SIN(THR/2.)**2                                   
       z8 = iz
       a8 = ia
       q28 = qsq
       w28 = pm**2 + 2. * pm * (e0 - ep) - qsq
 
c       CALL  F1F2QE07(z8, a8, q28, w28, F18, F28) 
c       CALL  F1F2QE09(z8, a8, q28, w28, F18, F28) 
       CALL  F1F2QE16(z8, a8, q28, w28, F18, F28) 
  
       W2 = F28 / (e0 - ep)
       W1 = F18 / PM
       CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
       CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
       sigma  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT
       if(abs(e0-5.157).lt.0.001.and.abs(ep-3.55).lt.0.01.and.
     >   abs(th-23.).lt.0.01)
     >  write(6,'(''dbgqe'',8f8.3)') e0,ep,th,q28,w28,F18,F28,sigma
       return
      endif

      DYDEP = 0.                                                        
      AMS   = avgM-PM                                                   
      THR   = TH*PI/180.                                                
      SINSQ = SIN(THR/2.)**2                                            
      COSTH = COS(THR)                                                  
      QSQ   = 4.*E0*EP*SIN(THR/2.)**2                                   
      TAU   = QSQ/4.0/PM**2                                             
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                          
c      if(prttst) write(8,'(1x,i3,11f7.3)') 
c     >  ig,e0,ep,th,qsq,gep,gmp,gmn,gen

                                                                        
! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     >  PAULI_SUP1,PAULI_SUP2)
      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)                           
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU)  
c     if(prttst) write(8,'(10f7.3)') pauli_sup2,gep**2,gen**2,w1,tau,w2
                                                                        
      CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
      CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
      RECOIL = PM/(PM+E0*(1.-COSTH))                                    
      DSIGE  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT*RECOIL                 

! new for proton. Use dw = e * dep / ep = dep / recoil
      if(doing_elas) then
        if(avgA.gt.1.5) then
         call NUC_FORM_FACTOR(QSQ,W1,W2,FFF)
         RECOIL = PM * avgA / (PM * avgA + E0 * (1.-COSTH)) 
         DSIGE  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT*RECOIL                 
cc         write(6,'(''dbg elas 1'',5f10.3)') QSQ,W1,W2,FFF,recoil
        endif
        dep = EP - E0 * recoil
        sigma = 0.
cc        write(6,'(''dbg elas 2'',5f10.3)') e0,ep,dep,smrdEp
        if(abs(dep/smrdEp).lt.5.) then
          sigma = dsige / smrdep  / sqrt(3.1416) * 
     >      exp(-dep**2 / smrdEp**2)
          if(sigma.lt.0.) write(6,'(''dbg elas 2'',7f8.3,3e11.4)') 
     >     e0,ep,dep,smrdEp,recoil,w1,w2,csmott,dsige,sigma
        endif
        return
      endif

! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
      if(IA.eq.3) kf=0.180
      if(iA.eq.3) Es=0.010 
c Use narrowe with for q.e.
cc      if(IA.eq.4) kf=0.200
      if(IA.eq.4) kf=0.180
      if(iA.eq.4) Es=0.015 
cc      if(iA.eq.4) Es=0.012 
      if(IA.gt.4) kf=0.220
      if(iA.gt.4) Es=0.015 
      if(IA.gt.7) kf=0.225
      if(iA.gt.7) Es=0.020 
      if(IA.gt.14) kf=0.225
      if(iA.gt.14) Es=0.020 
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
      if(ep.ge.e0) then 
cc        write(6,'(1x,''error,e0,ep='',2f8.3)') e0,ep
        return
      endif
      QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)                          
      kappa = qv / 2. / pm
      lam = (e0 - ep) / 2. / pm
      if(abs(kappa**2 - lam**2 - tau).gt.0.01) then
        write(6,'(1x,''error,tau='',3f8.3)') kappa**2, lam**2,tau
        return
      endif
      lamp = lam - Es / 2. / pm
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/pm)**2) -1.
      if(1.+lamp.le.0.) then
        write(6,'(1x,''error,lamp='',3f8.3)') lam,lamp
        return
      endif
      if(taup * (1. + taup).le.0.) then
        write(6,'(1x,''error,taup='',3f8.3)') kappa**2, lam**2,tau
        return
      endif
      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))
      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))
      nuL = (tau / kappa**2)**2
      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
      GM2bar = Pauli_sup1 * (iZ * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (iZ * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)
      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)
      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! Correction factors that make better agreement with world data
! on He and C from Day and NE5. See cn15.f
! BUT these don't work for Q2<0.3, so take them out again.
! actually, much better off without them as far as he/C concerned
      if(IA.gt.2.and.ia.lt.8) then
cc        GL = GL * (1.577 - 0.485 / qv)                
cc        GT = GT * (0.875 + 0.494 / qv)                   
      endif
      if(ia.ge.8) then
cc        GL = GL * (1.875 - 0.286 / qv)                
cc        GT = GT * (1.022 + 0.235 / qv)                   
      endif

! from Maria Barbaro: see superfy.m1
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
      sigma = csmott * FY * (nuL * GL + nuT * GT) / kf 
      if(sigma.lt.0.) then
cc        write(6,'(''ERROR, sigma...='',10f8.1)')
cc     >    sigma,csmott,fy,nuL,GL,nuT,GT
        sigma = 0.
      endif
      if(prttst) write(8,'(/1x,12f7.3)') e0,ep,th,W1,W2,csmott/10000.,
     >  fy,nul,gl,nut,gt,sigma/100.

! Use PWIA and Paris W.F. for deuteron
      if(IA.eq.2) then
! add extra effecive width by makeing qv wider for resolution
        qvp = qv
        if(smrelasp) qvp = qv * sqrt(1. + 
     >    (1.0 * smrdEp / 2. / 0.055 / qv)**2)
! qvp never used: THIS IS NOT IMPLEMENTED!!!!!
        pz = (2.*pm*(e0-ep) - qsq) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        depdpz = ep * qv / pm / e0
! Facotr of 100 and 0.01 cance to give:
        sigma = dsige * fyd(izz) / depdpz  
        if(prttst) write(8,'(i4,4f7.3)') izz,pz,depdpz,
     >    fyd(izz),sigma/100.
      endif

      return

! Are we changing quasi-elastic cross section due to Asymmetry?         
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0)                                
c    >     .AND.(INDEX(TARGET,'_P').GT.0) )                             
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >  CSMOTT*RECOIL,DSIGE)                                            
                                                                        
      QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)                          
      AMF2  = PM**2                                                     
      A     = 4.*QV**2-4.*(E0-EP+avgM)**2                               
      B     = -4.*QV*(AMS**2-AMF2-QV**2+(E0-EP+avgM)**2)                
      C     = (E0-EP+avgM)**4+(AMS**2-AMF2-QV**2)**2                    
     >       -2.*(E0-EP+avgM)**2*(AMF2+QV**2+AMS**2)                    
      DISC  = B**2-4.*A*C                                               
      IF (A.EQ.0..OR.DISC.LT.0.) RETURN                                 
      Y     = (-B-SQRT(DISC))/2./A                                      
                                                                        
C JACOBIAN DY/DEP                                                       
                                                                        
      EN    = SQRT(AMF2+(QV+Y)**2)                                      
      EMP   = SQRT(AMS**2+Y**2)                                         
      DENOMINATOR = (Y/EMP+(QV+Y)/EN)                                   
      IF (DENOMINATOR.LE.1.E-30) RETURN                                 
      DYDEP = (1.+(QV+Y)/EN*(EP-E0*COSTH)/QV)/DENOMINATOR               
                                                                        
C GAUSSIAN SMEARING FUNCTION                                            
                                                                        
      FY = 0.                                                           
      IF (iA.LE.4) THEN                                                 
! *** 034 way too small for 4He: swith to 100 
c**          SD = 0.034                                         
           SD = 0.100
      ELSE                                                              
           SD = 0.119                                                   
      END IF                                                            
c special case for C
c***      if(IA.eq.12) SD=0.09

      fy=0.
      IF ((Y**2/2./SD**2).LT.40.)                                       
     >     FY = 1./SD/SQRT(2.*PI)*EXP(-Y**2/2./SD**2)                   

! Modified to use fit to FY from nucl-th/9812078 (C. delgi Atti...)
! Modified 11/05 to get better agreement with Naida Fromin's data
      if(IA.eq.2) then
        a = 6.1 - 0.50
        alfa = .045 
        c1 = 0.018 
        c2 = 0.25 
      endif
      if(IA.eq.3) then
        a = 7.1 - 1.00
        alfa = .083
        c1 = 0.041 * 1.00
        c2 = 0.33 * 1.20
      endif
      if(IA.eq.4) then
        a = 6.8 - 1.00
        alfa = .167
        c1 = 0.106 * 1.20
        c2 = 0.65 * 1.20
      endif
      if(IA.gt.4.and.IA.le.30) then
        a = 5.1
        alfa = .166
        c1 = 0.083
        c2 = 0.57 * 1.20
      endif
      if(IA.gt.30) then
        a = 4.6
        alfa = .138
        c1 = 0.058
        c2 = 0.62 * 1.30
      endif
!*** peb pyb changed exponent from -6 to -8 10/05
!*** to try to get better agreement with data
! changed to 8 on 11/15/05
      fy = c1 * exp(-1.0 * a**2 * y**2) / (alfa**2 + y**2) +
     >   c2 * exp(-8.0 * abs(y))

      SIGMA = DSIGE*FY*DYDEP                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C====================================================================== 
      SUBROUTINE FYXSEC8(E0,E1,THT,XSEC)                                
C                                                                       
      DIMENSION x(1),a(4)                                               
      common/sigc /qx,qy,qz,q2,px,py,pz,p2,tpq                          
     >,pxf,pyf,pzf                                                      
      COMMON     /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                      
C initialize some values                                                
      nz = iz                                                           
      n = iA - iZ                                                       
      e0m=e0*1000.                                                      
      e1m=e1*1000.                                                      
!  iron coeff.                                                          
      nf=11                                                             
      nt=4                                                              
                                                                        
       if(iz.eq.26) then                                                
        a(1)=2.8316                                                     
        a(2)=44.624                                                     
        a(3)=0.3785                                                     
        a(4)=12.003                                                     
!   carbon coeff                                                        
       else                                                             
        a(1)=2.8757                                                     
        a(2)=41.922                                                     
        a(3)=0.3380                                                     
        a(4)=13.824                                                     
       endif                                                            
      XSEC = 0.0                                                        
      ep=e1m                                                            
      pz = pper                                                         
C     separation energy                                                 
      if(iz.eq.4) eps = 16.                                             
      if(iz.eq.26) eps=10.                                              
      hbc=0.1973289e+03                                                 
      fscnst=1./137.04                                                  
      rads=3.14159/180.                                                 
      rmp=.93828e+03                                                     
      tpm=2.*rmp                                                        
      tpm1=1./tpm                                                       
      thit=tht                                                          
      z=nz                                                              
      na=n+nz                                                           
      na1=na-1                                                          
      tma=2.*na*rmp                                                     
      tma1=2.*rmp*(na-1)                                                
      th2r=tht*rads/2.                                                  
      rma=na*rmp                                                        
      tn2=(tan(th2r))**2                                                
      thr=tht*rads                                                      
      tn2=(tan(th2r))**2.                                               
      Cost=COS(THT*RADS)                                                
      thtr=tht*rads                                                     
      sint=sin(thtr)                                                    
      cs2=(cos(th2r))**2.                                               
      s2=(sin(th2r))**2.                                                
      om=e0m-e1m                                                        
      ef=e0m-om                                                         
      q42=4.*e0m*ef*s2                                                  
      qx=e0m-ef*cost                                                    
      qy=-ef*sint                                                       
      qz=0.0                                                            
      qv2=q42+om**2                                                     
      qv=sqrt(qv2)                                                      
      q2=qv2                                                            
      qmp=q42/hbc/hbc                                                   
      tau=q42/4./rmp**2.                                                
      thpr=thp*rads                                                     
      csp=cos(thpr)                                                     
      tangle=qx/(-qy)                                                   
      thx=atan(tangle)                                                  
      e0gev = e0m/1000.                                                 
      epgev = e1m/1000.                                                 
      epsgev = eps/1000.                                                
      call yvalue(e0gev,epgev,tht,na,                                   
     >     epsgev,ygev)                                                 
      y = ygev                                                          
      p = ygev*1000.                                                    
      yp2=p                                                             
      px = p*sin(thx)                                                   
      py = -p*cos(thx)                                                  
      pxf=px+qx                                                         
      pyf=py+qy                                                         
      pzf=pz                                                            
      tpq=2.*e0m*px-(2.*px*cost+2.*py*sint)*ef                          
      p2=px*px+py*py+pz*pz                                              
      p=sqrt(p2)                                                        
      call sigcc(e0m,e1m,tht,sigp,sign)                                 
      sigt=z*sigp+n*sign                                                
      if (y .le. 0.0)x(1) = y                                           
      if (y .ge. 0.0)x(1) =-abs(y) ! sym                                
C     sum of two gaussians one located at zero                          
      functn = a(1)*exp(-a(2)*x(1)**2) +                                
     >    a(3)*exp(-a(4)*x(1)**2)                                       
C    ciofi's version                                                    
      top = qv                                                          
      bottom = sqrt(rmp**2 + qv2 + yp2**2 +                             
     > 2.*qv*yp2)                                                       
      domkdcos = top/bottom                                             
      dkcosdom = 1./domkdcos                                            
      dydom = dkcosdom                                                  
      xsec=fy*dydom*sigt*1000000.                                       
C      the om was in MeV                                                
22    return                                                            
      END                                                               
                                                                        
       subroutine yvalue(e0gev,epgev,theta,iatomic,eps,y)               
cccc       Implicit real*8 (a-h,o-z)                                        
!                                                                       
!subroutine to calulate the value of y given e,ep,theta,eps             
!                                                                       
       data rmp/0.93828/                                                 
       data rads/.01745329/                                             
       data thp/0.0/                                                    
       y = 0.0                                                          
       tpm=2.*rmp                                                       
       tpm1=1./tpm                                                      
!                                                                       
       na=iatomic                                                       
       na1=na-1                                                         
       tma=2.*na*rmp                                                    
       tma1=2.*rmp*(na-1)                                               
       th2r=theta*rads/2.                                               
       rma=na*rmp                                                       
       tn2=(tan(th2r))**2                                               
       thr=theta*rads                                                   
       tn2=(tan(th2r))**2.                                              
       cost=cos(theta*rads)                                             
       thtr=theta*rads                                                  
       sint=sin(thtr)                                                   
       cs2=(cos(th2r))**2.                                              
       s2=(sin(th2r))**2.                                               
                                                                        
       om=e0gev-epgev                                                   
       q42=4.*e0gev*epgev*s2                                            
       qx=e0gev-epgev*cost                                              
       qy=-epgev*sint                                                   
       qz=0.0                                                           
       qv2=q42+om**2                                                    
       qv=sqrt(qv2)                                                     
       q2=qv2                                                           
       thpr=thp*rads                                                    
       csp=cos(thpr)                                                    
       tangle=qx/(-qy)                                                  
       thx=atan(tangle)                                                 
                                                                        
       w=om - eps +na*rmp                                               
       wp=w**2+(na1*rmp)**2-rmp*rmp                                     
       c=4.*w*w*(na1*rmp)**2 +2.*wp*qv2-qv2*qv2-wp*wp                   
       b=qv*(4.*wp-4.*qv2)                                              
       a=4.*w*w-4.*qv2                                                  
       rad = b*b - 4.*a*c                                               
       if(rad .lt. 0.0)return                                           
       yp1 = (-b-sqrt(rad))/(2.*a)                                      
       yp2 = (-b+sqrt(rad))/(2.*a)                                      
       p = yp2                                                          
       if(p .ge. 0.0) p = abs(p)                                        
       y = p                                                            
                                                                        
       return                                                           
       end                                                              
                                                                        
C=========================================================              
       subroutine sigcc(e0m,e1m,tht,sigp,sign)                          
C gives the elementary sigm ep and en cross section for                 
C moving nucleon                                                        
      common/sigc /qx,qy,qz,q2,px,py,pz,p2,tpq,                         
     > pxf,pyf,pzf
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL*8 GEP,GEN,GMP,GMN
                                                                        
      pm=938.28                                                          
      cost=cos(tht*3.1416/180.)                                         
      sint=sin(tht*3.1416/180.)                                         
      fscnst = 1./137.04                                                
      hbc=197.3289                                                      
      pm2=pm*pm                                                         
      tpm1=1./2./pm                                                     
      teomc=2.*e0m*(1.-cost)                                            
      tan2=(1.-cost)/(1.+cost)                                          
      sig=5./137.04**2/(1.-cost)/e0m**2/tan2                            
      qm2=teomc*e1m                                                     
      ep=sqrt(pm2+p2)                                                   
      epf=sqrt(pm2+p2+q2+tpq)                                           
      qm2k=q2-(epf-ep)**2                                               
      rq2=qm2/q2                                                        
      v2=-qm2*tpm1**2                                                   
      omv2=1./(1.-v2)                                                   
      pt2=p2-tpq*tpq*0.25/q2                                            
      qmp=qm2/1000000.                                                  
      call nform(IG,DBLE(qmp),gep,gen,gmp,gmn)                                   
      f1p=gep-v2*gmp                                                    
      f1n=gen-v2*gmn                                                    
      cf2p=gmp-gep                                                      
      cf2n=gmn-gen                                                      
      pstcaq=qx*py-qy*px                                                
      ve=((ep+epf)*rq2*0.5-sqrt((rq2+tan2)/q2)*pstcaq)**2               
     >    +TAN2*PZ**2                                                   
      vm=((rq2*0.5+tan2)*qm2k-rq2*0.5*qm2)*0.5                          
c     efpf=pxf*cost+pyf*sint                                            
      sigm=sig*omv2**2/ep/epf !(epf-efpf) in ingo's version             
      sigp=sigm*(ve*(f1p**2+qm2k*(tpm1*cf2p)**2)+vm*(f1p+cf2p)**2)      
      sign=sigm*(ve*(f1n**2+qm2k*(tpm1*cf2n)**2)+vm*(f1n+cf2n)**2)      
      sigp=sigp*hbc*hbc                                                 
      sign=sign*hbc*hbc                                                 
      return                                                            
      end                                                               
C=======================================================================
                                                                        
      SUBROUTINE SECNUCLW(E0,EP,TH,SIGMA)                               
!--------------------------------------------------------------------
! SIGMA is cross section nb/GeV/str/nucleus.                                       
C Bastardization of Javier's SECNUC, inspired by patches, patches, and  
C patches splattered all over XSECFIT, INEFT, and this routine. Target  
C information is passed through common /targt/.                         
                                                                        
C SECNUCLW calculates the nuclear cross section in the deep inelastic   
C region using a variety of fits.                                  
C Added Fermi smearing Sept. 2005 P. Bosted
C Changed arguement of INELAST to WSQ rather than W
c Changed 8/06 to use better smearing pyb
!--------------------------------------------------------------------
      implicit none
      real e0,ep,th,sigma,sinsq,cossq,tansq,qsq,wsq,w,csmott,f1c,x
      real avgn,avga,avgM,amuM,r,dr,nu,eps,kappa,sigres,flux
      REAL*8 W1,W2,sigt,rc,q2,w1t,w2t
      INTEGER ISM
      logical goodfit
      integer iz,ia,i,iwp
      logical smrelasp, usearen
      real smrdEp,smrdW
      common/smrproton/ smrelasp,smrdEp,usearen
      COMMON      /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                     
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/
      logical usegrd,first/.true./
      common/testing/prttst,usegrd
      logical prttst
      real*8 xval1(50),xvall(50),temp(4)
      integer iq,iw
      real*8 w1sv(0:1000,0:100),w2sv(0:1000,0:100),w1tmp(0:1000)
      real*8 delw,w1a,w1b,delq2,w2a,w2b,q2lo,q2hi,w2tmp(0:1000),psm(7)

c make grid if first time through
      if(usegrd.and.first) then
        do iq=0,100
          do iw=0,1000
            wsq = -10. + 20. * (iw-0.5) / 1000.
cc            q2 = 0.1 * iq
c change to use log(q2) grid
            q2 = 10**(-2. + 3.*float(iq)/100.)
            nu = (wsq - pm**2 + q2) / 2. / pm
            x = q2 / 2. / pm / nu
            W1sv(iw,iq) = 0.
            W2sv(iw,iq) = 0.
            if(nu.gt.0.0 .and. x.lt.2.) then
             CALL INELAST(q2,DBLE(wsq),W1sv(iw,iq),W2sv(iw,iq)) 
c scale by dipole to make smoother
             W1sv(iw,iq) = W1sv(iw,iq) * (1. + q2/0.71)**2
             W2sv(iw,iq) = W2sv(iw,iq) * (1. + q2/0.71)**2
             if(w1sv(iw,iq).lt.0.0.or.w2sv(iw,iq).lt.0.) 
     >        write(6,'(''negative w1,w2'',2i4,2f10.3)') 
     >        iw,iq, W1sv(iw,iq),W2sv(iw,iq)
            endif
          enddo
! Add resolution smearing. 
          if(smrelasp) then
           do iw=0,1000
            w1tmp(iw) = W1sv(iw,iq)
            w2tmp(iw) = W2sv(iw,iq)
            W1sv(iw,iq) = 0.
            W2sv(iw,iq) = 0.
           enddo
! THIS IS HARD-WIRED for typical experiment
           smrdW = 0.030
           psm(1) = exp(-0.060**2 / smrdW**2) / sqrt(3.1416)
           psm(2) = exp(-0.040**2 / smrdW**2) / sqrt(3.1416)
           psm(3) = exp(-0.020**2 / smrdW**2) / sqrt(3.1416)
           psm(5) = exp(-0.020**2 / smrdW**2) / sqrt(3.1416)
           psm(6) = exp(-0.040**2 / smrdW**2) / sqrt(3.1416)
           psm(7) = exp(-0.060**2 / smrdW**2) / sqrt(3.1416)
           psm(4) = 1. - psm(1)- psm(2)- psm(3)- psm(5)- psm(6)- psm(7) 
           do iw=3,997
             do iwp = iw-3,iw+3
               W1sv(iwp,iq) = W1sv(iwp,iq) + w1tmp(iw) * psm(iwp+4-iw)
               W2sv(iwp,iq) = W2sv(iwp,iq) + w2tmp(iw) * psm(iwp+4-iw)
             enddo
           enddo
          endif
c          write(6,'(1x,''making grid...'',i4,f7.3,4f10.6)') iq,q2,
c     >      w1sv(600,iq),w1sv(900,iq),w2sv(600,iq),w2sv(900,iq)
        enddo
        first=.false.
      endif

C Calculate QSQ and WSQ for this kinematic point:                       
                                                                        
      SINSQ  = SIN(TH*3.1415927/180./2.)**2                             
      COSSQ  = 1.0-SINSQ                                                
      TANSQ  = SINSQ/COSSQ                                              
      QSQ    = 4.0*E0*EP*SINSQ                                          
      WSQ    = PM*PM+2*PM*(E0-EP)-QSQ                                   
      NU     = E0 - EP
      sigma = 0.

      CSMOTT = (.001)*(19732./(2.0*ALPHA*E0*SINSQ))**2*COSSQ            

! Moved Fermi smearing into inelast
      if(usegrd) then
        w1 = 0.
        w2 = 0.
c        iq = int(qsq*10.) 
c changed to log10 q2 bins
        iq = 0
        if(qsq.gt.0.01) iq = (alog10(qsq) + 2.0) * 100. / 3.0
        iw = -1
        if(wsq .gt.-10.0) iw = int((wsq + 10.)*50.)
        if(iw.ge.0.and.iw.lt.1000 .and.
     >     iq.ge.0.and.iq.lt.100) then
          delw = (wsq + 10.)*50. - iw
          if(delw.lt.0.0.or.delw.gt.1.) write(6,'(''error delw2'',
     >      i2,2f10.4)') iq,sqrt(wsq),delw
          W1a = w1sv(iw,iq)   + delw*(w1sv(iw+1,iq)  -w1sv(iw,iq))
          W1b = w1sv(iw,iq+1) + delw*(w1sv(iw+1,iq+1)-w1sv(iw,iq+1))
          W2a = w2sv(iw,iq)   + delw*(w2sv(iw+1,iq)  -w2sv(iw,iq))
          W2b = w2sv(iw,iq+1) + delw*(w2sv(iw+1,iq+1)-w2sv(iw,iq+1))
cc          delq2 = qsq*10. - iq
          q2lo = 10**(-2. + 3.*float(iq)/100.)
          q2hi = 10**(-2. + 3.*float(iq+1)/100.)
          delq2 = (qsq - q2lo) / (q2hi - q2lo)
          if(delq2.lt.0.0.or.delq2.gt.1.) then
c            write(6,'(''error delq2'',
c     >      i2,4f10.4)') iq,q2lo,qsq,q2hi,delq2
            delq2 = 0.
          endif
          w1 = w1a + delq2 * (w1b - w1a)
          w2 = w2a + delq2 * (w2b - w2a)
c divide out dipole
          w1 = w1 / (1. + qsq/0.71)**2
          w2 = w2 / (1. + qsq/0.71)**2
        endif
      else
        CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1,W2)!get s.f. /nucleon 
        if(w1.lt.0.0.or.w2.lt.0.0) then
          write(6,'(1x,''error, e,ep,th,q2,w,w1,w2='',6f8.3)') e0,ep,
     >      th,qsq,w,w1,w2
        endif
      endif

      SIGMA = (W2+2.*TANSQ*W1)*CSMOTT
C             (per nucleon, including emc effect (which includes neutron
C              excess correction), if any)                              
                                                                        
      SIGMA  = avgA*SIGMA  !per nucleus                                               
      if(sigma.lt.0.0.and.wsq.gt.1.0) write(6,'(1x,''error, sigma='',e10.3,2i4,
     >   2f6.3,4e10.3)') sigma,iw,iq,wsq,qsq,w1,w2,tansq,csmott
      sigma = max(0., sigma)                                                                  
! change cross section for polarized beam and target?                   
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0) )
c    >     .AND.( (INDEX(TARGET,'_P').GT.0).OR.                         
c    >            (INDEX(TARGET,'_A').GT.0) )   )  !inelastic only      
c    > CALL ASYM_INE(E0,EP,TANSQ,QSQ,TARGET,SIGMA)                      


cc      write(6,'(1x,3f7.3,e11.4)') E0,EP,TH,SIGMA
      RETURN                                                            
      END                                                               
                           

C ----------------------------------------------------------------------
                                                                        
      SUBROUTINE INELAST(QQ,WSQnom,W1sm,W2sm)                              
! Choose which inelastic model is to be used for resonance and for DIS  
! Returns structure functions per NUCLEON (11/3/95 SER)
C CHANGED TO USE WSQ RATHER THAN W peb 10/05
c changed to do smearing here now

      Implicit NONE
      REAL  avgN, avgA, avgM, amuM 
      COMMON /TARGT/ iZ, iA, avgN, avgA, avgM, amuM 
      INTEGER IZ,IA, ism, nsmr
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL   
      REAL*8 QQ,WSQ,W,X,W1,W2,F2_E665,R8,EMCFAC,F1C,RC,wsqnom,w1sm,w2sm
      REAL*8 MP2/0.8803/,MP/.93828/,F2,F1,F2nF2p,pi,alpha,sigL,sigT,Z,A
      REAL Q_RES_LIMIT/.2/                                              
      REAL QQ4,X4,F2_4,W1_4,W2_4,FITEMC_N,R4,DR4,NU,F2Allm,FNP_NMC
      REAL rnphi,rnp,rnplo,fact
      LOGICAL GD
! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/
! use w**2 for d-state
      real*8 fydr(15)/0.0192,0.0303,0.0538,0.1091,0.2544,0.6857,2.0061,
     >  3.5775,2.0061,0.6857,0.2544,0.1091,0.0538,0.0303,0.0192/
! *** with NO D-sate!
       real*8 fydn(15)/0.0040,0.0102,0.0277,0.0764,0.2150,0.6410,1.9589,
     >  3.5301,1.9589,0.6410,0.2150,0.0764,0.0277,0.0102,0.0040/
! using 1.5 * (1-cos**2)**2 for dtate
       real*8 fyd(15)/ 0.0094,0.0187,0.0411,0.0970,0.2462,0.6866,2.0207,
     >  3.6003,2.0207,0.6866,0.2462,0.0970,0.0411,0.0187,0.0094/
      logical first/.true./
      REAL*8 DW2DPF,pf,kf,es,dw2des,pz,qv,PM/0.93828/

      W1sm = 0.
      W2sm = 0.
      pi = 3.141593
      alpha = 1./137.036
      nu = (wsqnom - MP**2 + qq) / 2. / MP
      if(nu .le. 0.) return

! If model is 4, use new code for all values of A
      if(INEL_MODEL.EQ.4) THEN 
        Z = IZ
        A = IA
c        call F1F2IN07(Z, A, QQ, WSQnom, F1, F2)
        call F1F2IN09(Z, A, QQ, WSQnom, F1, F2)
c      call gsmearing(Z,A,WSQnom,QQ,F1,F2)
        W1sm = F1 / MP
        W2sm = F2 / nu
c need to make xsection per nucleon here!
c (they are converted back to per nucleus in SECNUCLW
        w1sm = w1sm / A
        w2sm = w2sm / A
        return
      endif

! Apply Fermi smearing to inelatic, except for H2
      if(first) then
        IF(IA .eq. 1) nsmr=1
        if(IA .ge. 2) nsmr=15
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
        if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
        if(IA.eq.3) kf=0.180
        if(iA.eq.3) Es=0.010 
cc       if(IA.eq.4) kf=0.210
        if(IA.eq.4) kf=0.180
        if(iA.eq.4) Es=0.015 
cc        if(iA.eq.4) Es=0.012 
        if(IA.gt.4) kf=0.165
        if(iA.gt.4) Es=0.015 
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
! adjust pf to give right width based on kf
        pf = 0.5 * kf 
        first = .false.
      endif

! assume this is 2 * pf * qv
      qv = sqrt(nu**2 + qq)
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + PM) 
      do ism = 1,nsmr
       if(IA.le.1) wsq = wsqnom
       if(IA.eq.2) wsq = wsqnom + (-0.3 + 0.6 * 
     >        (float(ism)-0.5)/15.) * dw2dpf
       if(IA.gt.2) WSQ = WSqnom + XX(ISM) * PF * DW2DPF - es * dw2des
       if(WSQ.gt.(0.93828 + 0.13957)**2) then
        W = SQRT(WSQ)

        IF(INEL_MODEL.EQ.0) THEN                                          
         CALL INEFT(QQ,W,W1,W2)   !Bodek fit                              
! 8/06 changed to use proton fit for all nuclei too!
! modified by simple n/p model
        ELSEIF(INEL_MODEL.EQ.3) THEN 
c        CALL CHRISTY705(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
c         CALL CHRISTY31606(WSQ,QQ,F1C,RC) ! Christy H2 resonance fi
         CALL CHRISTY0507(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
         NU = (WSQ +QQ -MP2)/(2.*MP)
         W1 = F1C / MP
         W2 = W1 /(1.0 + NU*NU/QQ) * (1.0 + RC)
! Simple model for n/p. Constnat in resonance region,
! and DIS at high W. Smoothly join from W=1.7 to 2.3
! Fixed to work for any Z,A 8/06 pyb
         if(IA.gt.1) then
          rnplo = 0.66
          X4 = QQ/(WSQ -MP2 +QQ) 
          qq4 = qq
          rnphi = FNP_NMC(X4,QQ4)
          fact = max(0., min(1., (wsq-1.7)/0.6))
          rnp = rnplo * (1.-fact) + fact * rnphi
          W1 = W1 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
          W2 = W2 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
         endif


        ELSEIF((INEL_MODEL.EQ.9).OR.(INEL_MODEL.EQ.12)) THEN              
         X = QQ/(WSQ -MP2 +QQ)                                           
         IF(QQ.GT.Q_RES_LIMIT) THEN
! Only H2 model. Calls F2GLOB for DIS.  
! Remember that F2GLOB returns structure functions/nucleon for H2 and D2
! Convert to Real*4 for this call
          QQ4=QQ
          X4=X                                       
!stuart for res, F2_Glob DI
          CALL INELSTU(QQ4, X4, F2_4, W1_4, W2_4, INEL_MODEL) 
          W2=W2_4
          W1=W1_4 
         ELSE  !out of range of inelstu.                                  
          CALL INEFT(QQ,W,W1,W2)   !Bodek fit                             
         ENDIF

!NMC model for DIS, INEFT for resonance
        ELSEIF(INEL_MODEL.EQ.1) THEN  
         X4 = QQ/(WSQ -MP2 +QQ) 
         QQ4=QQ
         CALL NMCDIS(QQ4,X4,W1_4,W2_4,amuM)
         W2=W2_4
         W1=W1_4
        ELSEIF(INEL_MODEL.EQ.2) THEN   ! E665 Fit (includes resonances
         X = QQ/(WSQ -MP2 +QQ) 
         IF(amuM.LT.1.5) THEN
          CALL  F2PGLO(QQ,X,F2_E665,R8) ! Prot E665
          EMCFAC=1.
         ELSE  ! Deuteron or heavier 
          CALL  F2DGLO(QQ,X,F2_E665,R8) ! Deut E665 
!with neutron excess 8/19/98
          EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD)) 
         ENDIF
         CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
         NU = (W**2 +QQ -MP2)/(2.*MP)
         W2 = F2_E665/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
       ELSEIF(INEL_MODEL.EQ.6) THEN   ! F2allm
         X4 = QQ/(WSQ -MP2 +QQ) 
         qq4 = qq
         F2_4 = F2ALLM(X4,QQ4)   ! proton fit
         F2 = F2_4
         x = x4
         IF(amuM.GT.1.5) THEN
! Changed to FNP_NMC 7/19/06: seems to work much better!
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))    
          F2 = (FNP_NMC(X4,QQ4) + 1.) / 2.  * F2
         endif
         if(amuM.gt.2.5) then
           EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD))
           F2 = F2 * EMCFAC
         endif
         CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
         NU = (W**2 +QQ -MP2)/(2.*MP)
         W2 = F2/NU 
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
        ELSEIF(INEL_MODEL.GE.100) THEN !joining a DIS model with resonance model
         X4 = QQ/(W**2 -MP2 +QQ) 
         QQ4=QQ
         CALL F2_JOIN(QQ4,X4,amuM,iZ,INEL_MODEL,W1_4,W2_4,iA)
         W2=W2_4
         W1=W1_4
        ELSE                                                              
         WRITE(6,'('' NOT GOOD H2 INELASTIC MODEL='',I3)')INEL_MODEL        
         STOP                                                             
        ENDIF                                                             
        if(IA.eq.1) then
          W1sm = W1
          W2sm = W2
        endif
        if(IA.eq.2) then
          W1sm = W1sm + W1 * FYD(ISM)/10.
          W2sm = W2sm + W2 * FYD(ISM)/10.
        else
          W1sm = W1sm + W1 * FY(ISM)
          W2sm = W2sm + W2 * FY(ISM)
c         if(prttst) write(8,'(1x,''insm'',10f8.3)') 
c    >          qsq, wsqp,w1,w1p,fy(ism)
        endif
       endif ! test of wsq
      enddo ! loop on ism
c     CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1P,W2P)!get s.f. /nucleon 
c     if(prttst) write(98,'(1x,''inelas sm'',6f10.3)') 
c    >     qsq, w,w1,w1p,w2,w2p
                                                                        
      RETURN                                                            
      END                                                               
                                                       

                                                                        
                                                                        
C=======================================================================
                                                                        
      subroutine radlength(X0)                                          
                                                                        
C calculates the radiation length based on Tsai74, table III.6, except  
C for iZ=1 values taken from Physics Letters 170B--Review of Particle   
C Properties. amuM is average atomic mass in C12-based atomic mass      
C units:                                                                
                                                                        
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
                                                                        
      if (iZ.le.4) then                                                 
           if (iA.eq.1) X0 =  61.28                                     
           if (iA.eq.2) X0 = 122.6                                      
           if (iA.eq.3.and.iZ.eq.1) X0 = 183.7                          
           if (iZ.eq.2) X0 =  94.322*amuM/4.0026                        
           if (iZ.eq.3) X0 =  82.756*amuM/6.9390                        
           if (iZ.eq.4) X0 =  65.190*amuM/9.0122                        
           return                                                       
      endif                                                             
                                                                        
      zz    = (iZ*.00729735)**2                                         
      f     = 1.202*zz - 1.0369*zz**2 + 1.008*zz**3/(1.+zz)             
      denom = (log(184.15/iZ**.3333)-f)*iZ**2 + log(1194/iZ**.6667)*iZ  
                                                                        
      X0 = 716.405*amuM/denom                                           
                                                                        
      return                                                            
      end                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE INTERPOL(E0,EP,TH)                                     
                                                                        
      CHARACTER*1 STAR                                                  
      DIMENSION ZERO(18)                                                
      DATA ZERO /18*0./                                                 
                                                                        
100   READ(7,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,F6.3,1X,A1)')                
     >         E0,EP,TH,X,Q2,W2,TASI,STAR                               
                                                                        
      IF (STAR.NE.'*'.AND.STAR.NE.'X'.AND.STAR.NE.'Q') THEN             
           CALL SECNUCLW(E0,EP,TH,SLAC)                                 
           WRITE(10,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,4F6.2,F6.3,/,25X,     
     >           6F6.2,1PE10.3,/)') E0,EP,TH,(ZERO(i),i=1,14),SLAC      
           GOTO 100                                                     
      ELSE                                                              
           RETURN                                                       
      ENDIF                                                             
      END                                                               
                                                                        
C=======================================================================
                                                                        
                                                                        
C=======================================================================
                                                                        
      subroutine weiz                                                   
                                                                        
C Calculates weizsacker mass formula from Segre iff avgM not supplied by
C READIN. Accuracy = +/-.002 GeV.                                       

      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
                                                                        
      avgA = iA                                                         
      if (iA.le.5) then                                                 
           if (iA.eq.1)             avgM =  .93828                      
           if (iA.eq.2)             avgM = 1.87537                      
ccc pyb to try to make elastic peak work
cc           if (iA.eq.2)             avgM = 2.0
                      
           if (iA.eq.3.and.iZ.eq.1) avgM = 2.8095                       
           if (iA.eq.3.and.iZ.eq.2) avgM = 2.8094                       
           if (iA.eq.4)             avgM = 3.7284                       
           return                                                       
      endif                                                             
                                                                        
      extra = 0.                                                        
      iN    = iA-iZ                                                     
      iN2   = iN/2*2                                                    
      iZ2   = iZ/2*2                                                    
      if (iN2.eq.iN.and.iZ2.eq.iZ) extra = -.012*(iA**(-.5))            
      if (iN2.ne.iN.and.iZ2.ne.iZ) extra =  .012*(iA**(-.5))            
                                                                        
      avgM = iN*.939573 + iZ*(.938280+.000511) - iA*.01567 +            
     +       .01723*(iA**(2./3.)) + .09315*((.5*iA-iZ)**2.)/iA +        
     +       .0006965*(iZ**2.)*(iA**(-1./3.)) + extra                   
                                                                        
      return                                                            
      end                                                               
                                                                        
C=======================================================================
                                                                        
      REAL Function Fgauss(T)                                                
      IMPLICIT NONE                                                                             
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL T, RADIUS,X2,CHAR
                                                                        
      Fgauss = 0.                                                       
      Radius = 1.07*avgA**(1./3.)  
! from H. de Vries: Nuclear Charge Density Distributions from Elastic Electron
!   Scattering. in Atomic Data and Nuclear Data Tables 36, 495(1987)
! 8/9/96
      IF(IA.EQ.205)RADIUS=5.470
      IF(IA.EQ.56) RADIUS=3.729    
      If(IA.EQ.28) RADIUS=3.085
      IF(IA.EQ.27) RADIUS=3.035                                                        
      x2     = (T/.197328**2)*Radius**2                                 
      char   = (T/.197328**2)*(2.4**2)/6.                               
      if (char.lt.80) Fgauss = exp(-char)/(1.+x2/6.)                    
      Return                                                            
      End                                                               
                                                                        
C=======================================================================
                                                                        
      REAL Function Fshell(T)                                                
      IMPLICIT NONE                                                                          
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM   
      REAL T, RADIUS,X2,ALP,CHAR
                                                                     
      Fshell = 0.                                                       
      Radius = 1.07*avgA**(1./3.)     
! from H. de Vries: Nuclear Charge Density Distributions from Elastic Electron
!   Scattering. in Atomic Data and Nuclear Data Tables 36, 495(1987)
!8/9/96
      IF(IA.EQ.16) RADIUS=2.737
      IF(IA.EQ.15) RADIUS=2.611
      IF(IA.EQ.14) RADIUS=2.540
      IF(IA.EQ.12) RADIUS=2.464  
      IF(IA.EQ. 9) RADIUS=2.519
      IF(IA.EQ. 7) RADIUS=2.39 
      IF(IA.EQ. 6) RADIUS=2.56 
      IF(IA.EQ. 4) RADIUS=1.676                                       
      x2     = (T/.197328**2)*Radius**2                                 
      alp    = (iZ-2.)/3.                                               
      char   = x2*(2.+3.*alp)/(12.+30.*alp)                             
      if (char.lt.80) Fshell = exp(-char)*(1.-alp*x2/(6.+15.*alp))      
      Return                                                            
      End                                                               
                                                                        
C=======================================================================
                                                                        
      REAL FUNCTION FDEL(T)   
      IMPLICIT NONE       
      REAL SQF,T
                                        
      SQF  = SQRT(T/.197328**2)                                         
      FDEL = 1.58/SQF*(ATAN(SQF/0.93)-2.*ATAN(SQF/3.19)+ATAN(SQF/5.45)) 
      IF (FDEL.LE.0.) FDEL = 0.                                         
      RETURN                                                            
      END                                                               
                                                                        
C ----------------------------------------------------------------------
c                                                                       
c     subroutine time(vt,tt)                                            
c                                                                       
c      call vttime(ivtime,ittime)                                       0
c     vt = float(ivtime)/100.0                                          
c     tt = float(ittime)/100.0                                          
c     return                                                            
c     end                                                               
                                                                        
C ====================================================================  
                                                                        
       SUBROUTINE ASYM_INE(E0,EP,TANSQ,QSQ,TARGET,SIGMA)                
!-------------------------------------------------------------------    
! Modifies polarized cross section because of E142 asymmetry.           
! A1 parameterization from Emlyn on 3/16/93                             
! Crude numbers for beam and target polarization and dilution factor    
!-------------------------------------------------------------------    
      IMPLICIT NONE                                                     
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER iZ,iA                                                     
      REAL avgN,avgA,avgM,amuM                                          
      REAL E0,EP,TANSQ,QSQ,SIGMA,X,EPS,D,NU,A1/0./,APAR,DELTA           
      REAL N_OVER_P,TARG_DILUTE/0./,A1N/0./,A1P/0./                     
      REAL R/.2/,W_D/.058/ !deuteron d state                            
      CHARACTER*7 TARGET                                                
      LOGICAL FIRST/.TRUE. /,GOOD_TARG/.FALSE./                         
                                                                        
      NU = E0 -EP                                                       
      X = QSQ/(2.*.93828 *NU)                                            
      EPS = 1./( 1. +2.*(1.+ NU**2/QSQ)*TANSQ  )                        
      D =( 1. -EP/E0 *EPS )/(1. +EPS*R)                                 
      N_OVER_P = 1.- .8*X     !ratio of cross sections                  
      IF(INDEX(TARGET,'E142').GT.0) THEN!  Neutron A1                   
       A1 = 1.58 * X**1.4 * (1.- X) + X * (2.5476 * X - 1.5476)         
! We are using He3 cross section, so asymmetry is diluted.              
       TARG_DILUTE = 1. + 2./N_OVER_P  ! ratio of He to neutron sigma   
       GOOD_TARG =.TRUE.                                                
       IF(FIRST) WRITE(6,'('' NEUT TARG'')')                            
      ELSE IF(INDEX(TARGET,'E143').GT.0)THEN ! NH3 target, proton asym  
       A1 = 1.025*(X**.12) * (1.-EXP(-2.7*X))  ! EMC parameterization   
       A1P=A1                                                           
       IF(IA.EQ.1) THEN ! proton target (with N in NH3 as radiator)     
        TARG_DILUTE =1.                                                 
       ELSEIF(IA.EQ.2) THEN !deuteron target                            
        A1N = 1.58 * X**1.4 * (1.- X) + X * (2.5476 * X - 1.5476)       
        A1 = (A1P/(1.+N_OVER_P) +A1N *(N_OVER_P)/(1.+N_OVER_P) )*       
     >    (1.- 1.5*W_D)                                                 
        TARG_DILUTE =1.                                                 
       ELSE  ! Perhaps weird case using NH3 cross section??             
        TARG_DILUTE = (10. + 8.*N_OVER_P)/3. ! NH3/H sigma ratio        
       ENDIF                                                            
       GOOD_TARG =.TRUE.                                                
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN !Parity violation          
c      IF(FIRST) THEN                                                   
c       WRITE(10,'('' Asym is 100 times parity violation'')')           
c       WRITE( 6,'('' Asym is 100 times parity violation'')')           
c      ENDIF                                                            
       GOOD_TARG=.TRUE.                                                 
       D=1.0
       A1 =  0.8E-4 * QSQ/D         !Fake a Sigma = delta sigma         
c      A1 = 100. * 0.80E-4 *QSQ/D   ! 100 Times parity violation        
       TARG_DILUTE =1.                                                  
      ENDIF                                                             
      APAR = D * A1/TARG_DILUTE !change in He or NH3 cross section      
! alter neutron xsection by differences.                                
      DELTA =     (1.+ APAR)/(1.-APAR)                                  
! If model is sigma parallel, then we calculate sigma anti-parallel     
! Measured asymmetry = antiparrallel - parallel                         
c      SIGMA = SIGMA* DELTA                                              
      SIGMA = SIGMA * (1. + apar)                                             
                                                                        
      IF(FIRST.AND.GOOD_TARG) THEN                                      
       FIRST =.FALSE.                                                   
       WRITE(6 ,'('' sig chngd dueto asym; D,A1,DELT,dilut='',          
     >   F5.2,1P,2E8.1,0P,F6.3,A8)')                                    
     >  D,A1,DELTA,TARG_DILUTE,TARGET                                   
       WRITE(10,'('' sig chngd dueto asym; D,A1,DELT,DILUT='',          
     >   f5.2,3f8.4,A8)')                                    
     >  D,A1,DELTA,TARG_DILUTE,TARGET                                   
       ENDIF                                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
!====================================================================== 
                                                                        
      SUBROUTINE ASYM_QEL(E0,EP,THETAR,QSQ,TARGET,GEN,GMN,GEP,GMP,      
     > MOTR,SIGMA)                                                      
!--------------------------------------------------------------         
! Calculate the Difference in cross sections in quasielastic            
! scattering from polarized Helium.to                                   
! From  Blankleider and Woloshyn, PR C29, 538 (1984)                    
!                                  -Steve Rock  4/13/93                 
!---------------------------------------------------------------        
      IMPLICIT NONE                                                     
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER iZ,iA                                                     
      REAL avgN,avgA,avgM,amuM                                          
      REAL E0,EP,THETAR,QSQ,GEN,GMN,GEP,GMP,SIGMA,GE/0./,GM/0./         
      REAL MOTR ! Mott * recoil factor                                  
      REAL ASYM/0./!Fractional change in cross section(phony initialize 
      CHARACTER*7 TARGET                                                
      REAL TARG_DILUTE/0./,TAN2,SIGN,SIGHE,SIGP,ASYM_P,ASYM_N           
      REAL G1/0./,G2/0./,TAU                                    
      REAL BETA /0./  ! angle between beam and Helium spins             
      REAL COSB/1./                                                     
      REAL MP2/.8803/,MP/.93828/                                           
      REAL GEX,GMX,ASYM_F,FDEL                                          
      INTEGER NTIMES/1/,TIMES/0/                                        
                                                                        
! M.J. Aluard,Phys. Rev. Lett 37, 1258(76) (slac experiment)            
! GE is only in numerator                                               
      ASYM_F(GEX,GMX) = TAU*GM *     !statement function                
     > (2.*MP*GEX/E0 +GMX*(2.*TAU*MP/E0 +2.*(1.+TAU)*TAN2))             
     >  /(GEX**2 + TAU*(GMX)**2 *(1. +2.*(1.+TAU)*TAN2))                
                                                                        
! If model is sigma parallel, then we calculate sigma anti-parallel     
! Measured asymmetry = antiparrallel - parallel                         
                                                                        
      TAU = QSQ/(4.*MP2)                                                
      TAN2 = (TAN(THETAR/2.))**2                                        
      IF(( INDEX(TARGET,'E142')+INDEX(TARGET,'E143')).GT.0) THEN        
       IF(INDEX(TARGET,'E142').GT.0) THEN !neutron asymetry             
C       G1 = - GMN *(GEN+ TAU* GMN)/(2.*(1.+TAU) )                      
C       G2 =   GMN * (GMN -GEN)/(4.*(1.+TAU) )                          
        GE = GEN                                                        
        GM = GMN                                                        
        SIGN = (GEN**2 +TAU*GMN**2)/(1.+TAU) +2.*TAU*GMN**2*TAN2        
        SIGHE=  (GEN**2 +2.*GEP**2 +TAU*(GMN**2+2.*GMP**2))/(1.+TAU) +  
     >   2.*TAU*(GMN**2 +2.*GMP**2)*TAN2                                
        TARG_DILUTE= SIGN/SIGHE                                         
        ASYM =ASYM_F(GE,GM)                                             
        SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)    
       ELSEIF(INDEX(TARGET,'E143').GT.0) THEN !proton asym, NH3 target  
        IF(IA.EQ.1) THEN ! proton                                       
         GE = GEP                                                       
         GM = GMP                                                       
         TARG_DILUTE =1.                                                
         ASYM =ASYM_F(GE,GM)                                            
         SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)   
        ELSEIF(IA.EQ.2) THEN ! deuterium                                
         SIGN = (GEN**2 +TAU*GMN**2)/(1.+TAU) +2.*TAU*GMN**2*TAN2       
         SIGP =((1.-FDEL(QSQ))*GEP**2 +TAU*GMP**2)/(1.+TAU)             
     >          +2.*TAU*GMP**2*TAN2                                     
         ASYM_P =ASYM_F(GEP,GMP)                                        
         ASYM_N =ASYM_F(GEN,GMN)                                        
         SIGMA =(SIGN *(1.+TARG_DILUTE*ASYM_N)/(1.-TARG_DILUTE*ASYM_N)  
     >        + SIGP *(1.+ TARG_DILUTE*ASYM_P)/(1.-TARG_DILUTE*ASYM_P) )
     >        * MOTR                                                    
        ENDIF                                                           
       ENDIF                                                            
C       ASYM=  TAU*GM *                                                 
C     > (2.*MP*GE/E0 +GM*(2.*TAU*MP/E0 +2.*(1.+TAU)*TAN2))              
C     >  /(GE**2 + TAU*(GM)**2 *(1. +2.*(1.+TAU)*TAN2))                 
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN  !parity violation         
       IF(TIMES.EQ.0) THEN                                              
        WRITE(10,'('' Asym is 100 times parity violtation=0.6E-4'')')   
        WRITE( 6,'('' Asymm is 100 times parity violation=0.6E-4'')')   
       ENDIF                                                            
C      ASYM = (0.6E-4 *QSQ -1.)/2.!Fake ASYM to make sigma = delta sig  
       ASYM = 100. *0.6E-4*QSQ    ! 100 Times parity violation.         
       TARG_DILUTE=1.                                                   
       SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)     
      ENDIF                                                             
                                                                        
      IF(TIMES.LT.NTIMES) THEN                                          
        TIMES = TIMES + 1                                               
        WRITE(6 ,'('' Quasi chnged  asym;A,DILUTE,THR,E0,EP,Q2='',      
     >   1P,E8.1,0P,F5.2,F4.1,2F5.1,F6.2)')                             
     >   ASYM,TARG_DILUTE,THETAR*57.296,E0,EP,QSQ                       
        WRITE(10,'('' Sig Quasi chnged dueto asym: M.J.Aluard formula''/
     >  '' A,DILUTE,THR,E0,EP,Q2='',                                    
     >   1P,E8.1,0P,F5.2,F4.1,2F5.1,F6.2)')                             
     >   ASYM,TARG_DILUTE,THETAR*57.296,E0,EP,QSQ                       
       ENDIF                                                            
                                                                        
      RETURN                                                            
      END                                                               

C $MEMBER=QUADMO, DATE=75061922, USER=KCE                               
      REAL*4 FUNCTION QUADMO_R(FUNCT,LOWER,UPPER,EPSLON,NLVL)  
      REAL*4 FUNCT,LOWER,UPPER,EPSLON                                   
      INTEGER NLVL                                                      
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*4  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS


         IF(LOWER.EQ.UPPER) THEN
          QUADMO_R =0.
          RETURN
         ENDIF

         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)
         FM = FUNCT(.5*(L+R)) 
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=ABS(ESTL)+ABS(ESTR)        
      ABAREA=AREA+ABAREA-ABS(EST)     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((ABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.     
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO_R = ROMBRG /12.0D0                                         
      RETURN                                                            
      END                                                               

C $MEMBER=QUADMO, DATE=75061922, USER=KCE                               
      REAL*4 FUNCTION QUADMO_U(FUNCT,LOWER,UPPER,EPSLON,NLVL)  
      REAL*4 FUNCT,LOWER,UPPER,EPSLON                                   
      INTEGER NLVL                                                      
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*4  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS
   
         IF(LOWER.EQ.UPPER) THEN
          QUADMO_U =0.
          RETURN
         ENDIF

         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)
         FM = FUNCT(.5*(L+R)) 
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=ABS(ESTL)+ABS(ESTR)        
      ABAREA=AREA+ABAREA-ABS(EST)     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((ABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.     
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO_U = ROMBRG /12.0D0                                         
      RETURN                                                            
      END                                                               

!----------------------------------------------------------------------------



      REAL*4 FUNCTION ATAILFL1(T)
!-----------------------------------------------------------------------------------
! Calls the E139/E140 code function ATAILFL_139 to calculate the 
!  "exact" nuclear tail.
!   SER 7/24/96
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL T
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL  E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL           ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
! TBEFOR, TBAL are the target material and windows before scattering point
! TAFTER, TAAL are the target material and windows after scattering point
! NOTE: TB and TA include the equivilent radiator: This is added on in ATAILFL
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*8 TB1,TA1, AtWt,ATAILFL_139
      common/testing/ prttst,usegrd
      logical prttst,usegrd
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen
                                                                        
      ATAILFL1 = 0.     
c removed 3/25/08
c      if(iz.lt.1) return
cxx      if(ia.eq.1.and.smrelasp) return

      AtWt= avgM/.93828 * 1.00797  ! atomic weight so that H=1.00797

      CALL RADIATORS(T)
      TB1 = TBEFOR +TBAL
      TA1 = TAFTER +TAAL
      ATAILFL1 =ATAILFL_139(DBLE(E0),DBLE(EP),DBLE(SIN2),TB1,TA1,
c 3/25/08 modifiled to ensure z>0
     >   DBLE(max(1,iZ)),DBLE(avgN),AtWt,0) !last arguement=0 for elastic
      if(prttst) write(*,111) t,tb1,ta1,atailfl1
 111  format(1x,'t,tb1,ta1,atailfl1=',3f8.5,e11.2) 
      RETURN
      END  
!==============================================================



      REAL*4 FUNCTION ATAILFL_QE(T)
!-----------------------------------------------------------------------------------
! Calls the E139/E140 code function ATAILFL_139 to calculate the 
!  quasi-elastic tail integrating over COSK.
!   SER 7/31/96
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL T
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL  E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL           ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 

! TBEFOR, TBAL are the target material and windows before scattering point
! TAFTER, TAAL are the target material and windows after scattering point
! NOTE: TB and TA include the equivilent radiator: This is added on in ATAILFL
      REAL*8 TB1,TA1,ATAILFL_139

      CALL RADIATORS(T)
      TB1 = TBEFOR +TBAL
      TA1 = TAFTER +TAAL
      ATAILFL_QE = ATAILFL_139(DBLE(E0),DBLE(EP),DBLE(SIN2),TB1,TA1,
     >  DBLE(iZ),DBLE(avgN),DBLE(AvgA),1) !last arguement=1 for Quasi Elastic
      RETURN
      END  
!==============================================================
      SUBROUTINE NUC_FORM_FACTOR(QSQ,W1,W2,FF)
!----------------------------------------------------------------------------
! Get Nuclear Form Factor from various models
!-----------------------------------------------------------

      IMPLICIT NONE
      REAL QSQ,W1,W2,FF
      COMMON    /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
                       
      INTEGER IZ,IA       
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL                                              
      REAL AVGN,AVGA,AVGM,AMUM
      REAL A_ULMAR,B_ULMAR,w2_old,PM,TAU,FDEL,FSHELL,FGAUSS,FF_BESSEL
ccc      REAL ffHe4, ffBe, ffC12, ffAl
      REAL*8 GE,GM,GEP,GEN,GMP,GMN
      LOGICAL OUT_OF_RANGE
      PARAMETER (PM    = 0.93828)  

    
      TAU = QSQ/4./PM**2                                                
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                                   
                                                                        
      IF (iA.EQ.1) THEN                                                 
c fixed 3/25/08
        if(iz.eq.1) then
           W1 = TAU*GMP**2                                              
           W2 = (GEP**2+W1)/(1.+TAU)                                    
        else
           W1 = TAU*GMN**2                                              
           W2 = (GEN**2+W1)/(1.+TAU)                                    
        endif
      ELSEIF (iA.EQ.2) THEN  
        IF((IDUT.GE.11).AND.(IDUT.LE.14).AND.(QSQ.LE.3.5))THEN   !Tjon fits
           !Ulmar d2 elastic mode                                       
           CALL DEUT_U1(IDUT-10,QSQ,A_ULMAR,B_ULMAR)                     
           W1 = B_ULMAR/2.  ! sigma = sig_mot(A + B*tan...)
           W2 = A_ULMAR 
        ELSEIF(IDUT.EQ.1) THEN ! Linda Stuart's Model installed 5/30/96
           CALL FFD(DBLE(QSQ),GE,GM)   
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU) 
        ELSE   ! old  elastic deuterium from original code   
           FF  = FDEL(QSQ)                                              
           W1  = FF**2*TAU*.6667*(GMN+GMP)**2                           
           W2  = W1+(FF*(avgN*(GEN+TAU*GMN)+GEP+TAU*GMP)/(1.+TAU))**2   
        ENDIF
      ELSEIF(iA.GE.3) THEN
       W1=0.
       OUT_OF_RANGE =.FALSE.
       IF(NUC_MODEL.EQ.1) THEN
         FF = FF_BESSEL(QSQ,OUT_OF_RANGE) !only for some Nuclei,feor limited Q2
         W2 = (iZ*FF)**2    
       ENDIF
       IF(OUT_OF_RANGE.OR.(NUC_MODEL.EQ.0)) THEN !use if FF_BESSEL out of range
        IF (iA.EQ.3) THEN  !3HE  ! added 5/30/96  SER
           CALL FFHE3(DBLE(QSQ),GE,GM)
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU)  
           W2_old  = (iz*FSHELL(QSQ) )**2
cc these don't seem to be working right: need more testing
cc to see if actually better than ff_bessel
cc      ELSEIF (iA.EQ.4) THEN
cc         W2  = ffHe4(QSQ)
cc      ELSEIF (iA.EQ.9) THEN
cc         W2  = ffBe(QSQ)
cc      ELSEIF (iA.EQ.12) THEN
cc         W2  = ffC12(QSQ)
cc      ELSEIF (iA.EQ.27) THEN
cc         W2  = ffAl(QSQ)
        ELSEIF (iA.LE.20) THEN
           FF  = FSHELL(QSQ)
           W2  = (iZ*FF)**2                                              
        ELSE     !ia >20
           FF  = FGAUSS(QSQ) 
           W2  = (iZ*FF)**2                                             
        ENDIF
       ENDIF                  
      ENDIF   ! iA>+3
      RETURN
      END
!---------------------------------------------------------------------

!=========================================================================

      SUBROUTINE QE_PEAK(QSQ,E0,W1,W2)
!------------------------------------------------------
! E0 is a guess of the incident energy used for the vanOrden calculation
!  of Pauli Suppression.  It is calculated from Es- OM
!   where OM is the photon energy.
! Calculate nucleon W1 and W2 from Q2
!--------------------------------------------------------
      IMPLICIT NONE
      REAL*4 QSQ,E0,W1,W2
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON/TARGT/ iZ,iA,avgN,avgA,avgM,amuM                       
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM
      REAL TAU,PAULI_SUP1,PAULI_SUP2,PM,PM24
      REAL*8 GEP,GEN,GMP,GMN
      PARAMETER (PM    = 0.93828, PM24=3.5216)    



      TAU = QSQ/PM24                                               
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)


! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)

      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)                           
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU)  

      RETURN
      END


      SUBROUTINE PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)
!-----------------------------------------------------------------------
! Gets Pauli Suppression factor for quasi-elastic scattering from
! Several different Models.
! The VanOrden Model is very Slow.
! Used by both INTERNAL and EXTERNAL
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*4 QSQ,E0,PAULI_SUP1,PAULI_SUP2,FDEL
      REAL Q,TAU,FSHELL,FGAUSS,FF_BESSEL,XX
      INTEGER PAULI_MODEL,NUC_MODEL,iA
      REAL FF
      REAL*4 MP24/3.52/   
      LOGICAL OUT_OF_RANGE
      REAL*4 P_FERMI/.25/

      PAULI_SUP1 = 1.
      PAULI_SUP2 = 1.
      if(IA.eq.1) return

      IF(PAULI_MODEL.EQ.0) THEN
        FF = 0.                                                           
        IF (iA.EQ.2.AND.QSQ.LE.8.)  FF = FDEL(QSQ) 
        IF(IA.GT.2) THEN
         OUT_OF_RANGE =.FALSE.
         IF(NUC_MODEL.EQ.1) THEN
          FF = FF_BESSEL(REAL(QSQ),OUT_OF_RANGE) !only for some Nuclei,feor limited Q2
         ENDIF
         IF(OUT_OF_RANGE.OR.(NUC_MODEL.EQ.0)) THEN !use if FF_BESSEL out of range
          IF (iA.LE.20) THEN
            FF  = FSHELL(QSQ)
          ELSE     !ia >20
            FF  = FGAUSS(QSQ) 
          ENDIF
         ENDIF  
        ENDIF                           
        Pauli_sup2 = (1.-FF**2)    !old model from Stein
        PAULI_sup1 =1.
      ELSE
       IF(PAULI_MODEL.EQ.1) THEN !Tsai RMP 46,816(74) eq.B54
           TAU   = QSQ/MP24     
           Q=SQRT(QSQ*TAU+QSQ)
           IF((Q .GT. 2.*P_FERMI).OR.(iA.EQ.1)) THEN
              PAULI_SUP2 =1.0
           ELSE
              PAULI_SUP2=3.0*Q*(1.0-0.08333*(Q/P_FERMI)**2)
           ENDIF
       ELSEIF(PAULI_MODEL.EQ.3) THEN
           CALL  Q_E_VANORDEN(QSQ,E0,PAULI_SUP2)
       ELSEIF(PAULI_MODEL.EQ.2) THEN ! no suppression
           PAULI_SUP2 =1.0   ! No suppression
       ELSE
           PAULI_SUP2 =1.0   ! No suppression
       ENDIF
        PAULI_SUP1= Pauli_sup2
      ENDIF
      IF(PAULI_MODEL.EQ.4) THEN ! Louk's formula for Deuterium
        TAU   = QSQ/MP24     
        Q=SQRT(QSQ*TAU+QSQ)
        xx = q * 1000. / 100.
        Pauli_sup2 = (1.  - 1. / (1. + xx * (0.0060 + xx * (0.3882 + 
     >    xx * 0.2477))))
        PAULI_SUP1= Pauli_sup2
      endif
           
      RETURN
      END
!-------------------------------------------------------------------------

c      SUBROUTINE READIN_EXT(LEN_FILENAME,FILENAME)     
      SUBROUTINE READIN_EXT()     
! used to be Called from cexternal.c
      IMPLICIT NONE
      COMMON /TARGT/  iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM                                   
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON /TTYPE/  TARGET                                            
      CHARACTER*7     TARGET
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG,TWALL,TBEAM,TSPEC
      INTEGER ITARG,NSEG,JTARG
      INTEGER I
      CHARACTER*72    COMMENT(6),RAD_STRING/'  '/                           
      CHARACTER*80    EXTERNAL_RUNPLAN,EXTERNAL_TARGET,EXTERNAL_OUT
      logical doeg1b
      common/experiment/ doeg1b
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen


c      OPEN(UNIT=20,FILE=FILENAME(1:LEN_FILENAME))
c      OPEN(UNIT=20,FILE='input/externals/kpp_shms_488.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_short.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_pol_3He.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_A1n_d2n_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_A1n_d2n_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_A1n_d2n_hms_dflay.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_A1n_d2n_shms_dflay.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_test_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_test_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_O_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_O_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Si_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Si_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Ba_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Ba_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Al_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Al_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Ca_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Ca_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Sr_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_Sr_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_H_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_up_win_H_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_O_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_O_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Si_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Si_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Ba_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Ba_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Al_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Al_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Ca_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Ca_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Sr_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_Sr_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_H_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_Empty_down_win_H_shms_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_20deg.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_20deg.inp')
      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_kin_C.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_kin_Z.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_C.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_Z.inp)
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_kin_A.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_kin_X.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_A.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_X.inp)
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_1pass_lo.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_1pass_lo.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_1pass_lo.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_1pass_lo.inp)
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_1pass_hi.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_1pass_hi.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_1pass_hi.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_1pass_hi.inp)
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_kin_3.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_kin_C.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_3.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_C.inp)
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_hms_kin_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_d2n_3He_shms_kin_B.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_hms_4.inp')
c      OPEN(UNIT=20,FILE='input/externals/test_N2_ref_shms_B.inp)


      READ(20,'(A)') COMMENT(1)
      WRITE(6,'(A)') COMMENT(1)
      READ(20,'(A)') EXTERNAL_RUNPLAN
      WRITE(6,'(A)') EXTERNAL_RUNPLAN
      READ(20,'(A)') EXTERNAL_TARGET
      WRITE(6,'(A)') EXTERNAL_TARGET
      READ(20,'(A)') EXTERNAL_OUT
      WRITE(6,'(A)') EXTERNAL_OUT
      close(unit=20)

      OPEN(UNIT=10,FILE=EXTERNAL_OUT)
      OPEN(UNIT=7,FILE=EXTERNAL_RUNPLAN)  
      OPEN(UNIT=5,FILE=EXTERNAL_TARGET)     

      do i=1,6
        READ(7,'(a)') COMMENT(i)             
        write(6,'(a)') comment(i)
      enddo
      doeg1b=.false.
      if(comment(1)(1:4).eq.'EG1b' .or.
     >   comment(1)(2:5).eq.'EG1b') then
        doeg1b=.true.
        write(6,'(''using eg1b format'')')
      endif
      READ(5,100) (COMMENT(i),i=4,6),                                   
     +            iZ,iA,                                                
     +            avgA,avgM,                                            
     +            target,                                               
     +            ttarg,twall,tbeam,tspec,                              
     +            NSEG,                                                 
     +    IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL

    
      IF(TARGET.EQ.'E140XH2') ITARG=1                                   
      IF(TARGET.EQ.'E140XH1') ITARG=2                                   
      IF(TARGET.EQ.'E140XD1') ITARG=3                                   
      IF(TARGET.EQ.'E140XD2') ITARG=4                                   
      IF(TARGET.EQ.'E140XAL') ITARG=5                                   
      IF(TARGET.EQ.'E140XAS') ITARG=6                                   
      IF(TARGET.EQ.'E140XBE') ITARG=7                                   
100   FORMAT      (3(A72,/),/,                                          
     +             2(15X,I5,/),/,                                       
     +             2(15X,F14.8,/),/,                                    
     +              (20X,A7,/),/,                                       
     +             4(15X,F14.8,/),                                      
     +             15X,I5,/,/,/,         ! NSEG, skip blank and ing   
     +             15X,I5,/,/,   ! IG,  elastic nucleon form factor model
     +             15X,I5,/,     ! IDUT                                      
     +             15X,I5/    ! INEL_MODEL
     +             15X,I5/     !Pauli suppression     
     +             15X,I5/     ! Nuclear Tail Method: NUC_METHOD
     +             15X,I5)     ! Nuclear Form Factor: NUC_MODEL

      call Q_E_VANORDEN_INIT(REAL(iZ),REAL(iA),.25,1)!Fermi Momentum=.25
      if (avgA.eq.0..or.avgM.eq.0) call weiz                            
      if(IA.eq.1) avgA = 1.
! pyb change avgM if smearing proton elastic
      if(smrelasp) avgM = max(avgM, 2.0)
      avgN = avgA-iZ                                                    
      amuM = avgM/.931501                                               
      IF(INDEX(TARGET,'E143').GT.0)CALL RADLEN43_INIT(RAD_STRING)
      WRITE(10,'(''EXTERNAL/TSAI:  VERSION OF 8/8/96''/
     > /A72/4(A72/)/,A72/A72/)')   
     +COMMENT(1),COMMENT(2),COMMENT(3),COMMENT(4),COMMENT(5),COMMENT(6),
     > RAD_STRING
      WRITE(10,200) iZ,iA,avgA,amuM,avgM,target,ttarg,twall,tbeam,tspec,
     + Nseg,IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
200   FORMAT('TARGET (iZ,iA) = (',I3,',',I3,');   avgA = ',F7.3,        
     +                 '; amuM=',F7.3,'; avgM=',F7.3,/                  
     +                 ' targ=',A7,'; ttarg(rl)=',F7.5,'; twall=',F7.5, 
     +                 '; tbeam=',F7.5,'; tspec=',F7.5/                 
     +       '         Nseg  = ',I3/' MODELS: P elastc(ikk)=',i2,
     >     '; deut elastic=',i2, '; inelastc=',I3, 
     >     '; Pauli Supresion=',i2 /
     >     ' NUCLEAR METHOD=',I2,';   NUCLEAR_MODEL=',I2/)
                                                                        
      RETURN                                                            
      END                                                               
                          
      	REAL function FF_BESSEL ( T ,OUT_OF_RANGE)
        IMPLICIT NONE
C	Calculates PWBA form factor  for Si-28 or O-16 using 
C	Fourier-Bessel coefficients given by H. de Vries et al., 
C	Atomic Data and Nuclear Data Tables (1986).
C
C	Note: Other nuclides can be entered, but at this time this
C        is only good for 
C       He-3, C-12, N-15, O-16, Al-27, Si-28, Fe-56, Cu-65,
!----------------------------------------------------------------
!3/28/03 Corrected for divide 0/0 when qm =0  - Steve Rock
        
        REAL T
        LOGICAL OUT_OF_RANGE
	common /targt/ iZ, iA, avgN, aavgA, avgM, amuM
        real avgN, aavgA, avgM, amuM
        integer iZ, iA
	real*8 a(20),a16(14),a28(14),a3(12),a12(16),a15(8),a27(12),
     >   a56(17), a65(17),a196(15)
        INTEGER N,J,I
        REAL*8 R_MAX,QMU,Q,FF,QMUX,QM,QP,SINQR

        REAL*8 PI4/12.56637062/ ! 4*PI
        REAL*8 PI2/ 6.28318531/ ! 2*PI
C
C     de Vries, de Vries & de Jager FB coefficients:


!3He
      data a3/0.20020E-01, 0.41934E-01, 0.36254E-01, 0.17941E-01,
     1        0.46608E-02, 0.46834E-02, 0.52042E-02, 0.38280E-02,
     2        0.25661E-02, 0.14182E-02, 0.61390E-03, 0.22929E-03/
!12C
      data a12/0.15721E-01, 0.38732E-01, 0.36808E-01, 0.14671E-01,
     1        -0.43277E-02,-0.97752E-02,-0.68908E-02,-0.27631E-02,
     2        -0.63568E-03, 0.71809E-04, 0.18441E-03, 0.75066E-04,
     3         0.51069E-04, 0.14308E-04, 0.23170E-05, 0.68465E-06/


!15N7
      data a15/0.25491E-01, 0.50618E-01, 0.29822E-01, -0.55196E-02,
     1        -0.15913E-01,-0.76184E-02, -.23992E-02, -0.47940E-03/
!16 O8
      data a16/0.20238E-01, 0.44793E-01, 0.33533E-01, 0.35030E-02,
     1          -0.12293E-01,-0.10329E-01,-0.34036E-02,-0.41627E-03,
     2          -0.94435E-03,-0.25571E-03, 0.23759E-03,-0.10603E-03,
     3           0.41480E-04, 0.00000E-03/

!27Al
      data a27/0.43418E-01, 0.60298E-01,  0.28950E-02, -0.23522E-01,
     1        -0.79791E-02, 0.23010E-02,  0.10794E-02,  0.12574E-03,
     2        -0.13021E-03, 0.56563E-04, -0.18011E-04,  0.42869E-05/

!28Si
      data a28/0.33495E-01, 0.59533E-01, 0.20979E-01,-0.16900E-01,
     1          -0.14998E-01,-0.93248E-03, 0.33266E-02, 0.59244E-03,
     2          -0.40013E-03, 0.12242E-03,-0.12994E-04,-0.92784E-05,
     3           0.72595E-05,-0.42096E-05/

!56Fe26 
      data a56/ 
     1  .42018E-1,  .62337E-1,  .23995e-3, -.32776E-1, -.79941E-2,
     2  .10844E-1,  .49123e-2, -.22144e-2, -.18146E-3,  .37261E-3,
     3 -.23296E-3,  .11494E-3, -.50596E-4,  .20652E-4, -.79428E-5,
     4  .28986E-5, -.10075E-5/         

!65Cu29 
        data a65/0.45444E-01, 0.59544E-01, -.94968E-02, -.31561e-01,
     1           0.22898E-03, 0.11189E-01, 0.37360E-02, -.64873E-03,
     2	         -.51133E-03, 0.43765E-03, -.24276E-03, 0.11507E-03, 
     3           -.49761E-04, 0.20140E-04, -.76945E-05, 0.28055E-05,
     4		 -.97411E-06 /


!196Pt78
         data a196/
     1     .50218e-1,  .53722e-1, -.35015e-1, -.34588e-1,  .23564e-1,
     2     .14340e-1, -.13270e-1, -.51212e-2,  .56088e-2,  .14890e-2,
     3    -.10928e-2,  .55662e-3, -.50557e-4, -.19708e-3,  .24016e-3/  

C	Change to units of fm**(-1)
	q = sqrt ( T ) / 0.197328
	if( q .le. 0.005 ) then
            OUT_OF_RANGE=.FALSE. ! line added 4/29/98
	    FF_BESSEL = 1.00000
	    return
	endif

        FF_BESSEL=0.
        OUT_OF_RANGE=.FALSE.
        R_max=0.
	if( iA .eq. 28 ) then
            if(q.gt.2.64) OUT_OF_RANGE=.TRUE.
	    R_max = 8.
	    n = 14
	    do i = 1,n
	        a(i) = a28(i)
	    enddo
        elseif( iA .eq. 16 ) then
            if(q.gt.2.77) OUT_OF_RANGE=.TRUE.
            R_max=8.
	    n = 13
	    do  i = 1,n
	      a(i) = a16(i)
	    enddo
        elseif( iA .eq. 3 ) then
            if(q.gt.10. ) OUT_OF_RANGE=.TRUE.
            R_max=5.
	    n = 12
	    do  i = 1,n
	      a(i) = a3(i)
	    enddo
        elseif( iA .eq. 12 ) then
            if(q.gt.4.01 ) OUT_OF_RANGE=.TRUE.
            R_max=8.
	    n = 16
	    do  i = 1,n
	      a(i) = a12(i)
	    enddo
        elseif( iA .eq. 15 ) then
            if(q.gt.3.17) OUT_OF_RANGE=.TRUE.
            R_max=7
	    n = 8
	    do  i = 1,n
	      a(i) = a15(i)
	    enddo
        elseif( iA .eq. 27 ) then
            if(q.gt.2.70) OUT_OF_RANGE=.TRUE.
            R_max=7
	    n = 12
	    do  i = 1,n
	      a(i) = a27(i)
	    enddo
        elseif( iA .eq. 56 ) then
            if(q.gt.2.22) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.51) OUT_OF_RANGE=.TRUE.
            R_max=9
	    n = 17
	    do  i = 1,n
	      a(i) = a56(i)
	    enddo

        elseif(( iA .eq. 64).or.(iA.eq.65) ) then
            if(q.gt.2.22) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.51) OUT_OF_RANGE=.TRUE.
            R_max=9
	    n = 17
	    do  i = 1,n
	      a(i) = a65(i)
	    enddo
        elseif( iA .eq. 196)  then
            if(q.gt.2.28) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.34) OUT_OF_RANGE=.TRUE.
            R_max=12
	    n = 15
	    do  i = 1,n
	      a(i) = a196(i)
	    enddo
	else      
            out_of_range=.true.
	endif
        if(out_of_range.or.r_max.eq.0.) then
          ff_bessel=0.
          Return
        endif 


	qmu = 3.14159265 / R_max

        ff=0.
        sinqR = sin(q*R_max)
        do j=1,n
	 qmux = qmu * float(j)
	 qm =  q - qmux
         qp =  q + qmux
         if(abs(qm).gt.1.e-6) then
          ff=ff+ a(j)*((-1.)**j)*sinqR/(qm*qp)
         else
          ff= ff +a(j)*R_max**2/(PI2*j)
         endif
        enddo
        if((q*R_max.gt.1.E-20).and.(ff.lt.1.E20)) then
         ff_bessel =PI4/FLOAT(IZ)/q *ff
        else
         ff_bessel=0.
        endif 
	Return
	End

                                                                       


                                                                       
                          
!-----------------------------------------------------------------------
C-----
C       SUBROUTINE deut_elastic
C
C       AUTHOR: S. Van Verst
C       DATE:   DEC-1991
C
C       PURPOSE:
C               Calculate cross section and polarizations for deuteron
C               elastic scattering.
C
C       MODIFICATIONS:
C               Steve Rock 10/93: modified for use in Radiative
C               Correction program.  Only input arguement is QSQ.
C       DEFINITIONS:
!          IMODEL is the model number: 1= impulse approx
!                                      2= mec
!                                      3= rsc
!                                      4= rsc+mec
!            QSQ in GeV
!            A and B are from Sigma = sig_mott * [A+B*Tan2(theta/s)]
C----------------------------------------------------------------------
C
      SUBROUTINE DEUT_U1(IMODEL,QSQ,A,B)
 
      IMPLICIT NONE
      CHARACTER*6 FCFILE, FQFILE, FMFILE
      CHARACTER*20 C_FILE,Q_FILE,M_FILE
      INTEGER J, IMODEL,ILUN
 
      REAL QSQ  ! QSQ in Gev/c
      DOUBLE PRECISION QMU2,ETA,QMU2_FM
      DOUBLE PRECISION RINTERPQ
      DOUBLE PRECISION FC(100),FQ(100),FM(100),Q_2(100)
      DOUBLE PRECISION G_C,G_M,G_Q,GC2,GM2,GQ2
      REAL A,B
      DOUBLE PRECISION ALPHA,HBARC,MD
      LOGICAL FIRST /.TRUE./
      INTEGER N_D_PTS
 
      PARAMETER (MD     = 1875.630)        !Deuteron mass in MeV
      PARAMETER (ALPHA  = 7.29735E-3)      !fine structure constant
      PARAMETER (HBARC   = 197.3286)       ! MeV-Fermi
 
      if(first) then
           DO J=1,100
              FC(J)  = 0.
              FM(J)  = 0.
              FQ(J)  = 0.
              Q_2(J) = 0.
           ENDDO
C
C ------------------------------------------------------------------------
C       Ask for deuteron form factor model (from Tjon), then read them
C       in from files in DAT$D
C ------------------------------------------------------------------------
C
           WRITE(10,
     >      '('' Deut Elastic Model: ia,iamec,rsc,rscmec (1,2,3,4)='',
     >       I3)')IMODEL
           IF(IMODEL .EQ. 1)THEN
              FCFILE =  'iactjn'
              FQFILE =  'iaqtjn'
              FMFILE =  'iamtjn'
           ELSEIF(IMODEL .EQ. 2) THEN
              FCFILE =  'mecctj'
              FQFILE =  'mecqtj'
              FMFILE =  'mecmtj'
           ELSEIF(IMODEL .EQ. 3) THEN
              FCFILE =  'rscctj'
              FQFILE =  'rscqtj'
              FMFILE =  'rscmtj'
           ELSE
              FCFILE =  'rscmct'
              FQFILE =  'rscmqt'
              FMFILE =  'rscmmt'
           ENDIF
           C_FILE = FCFILE//'.tjon_input'
           Q_FILE = FQFILE//'.tjon_input'
           M_FILE = FMFILE//'.tjon_input'
           WRITE(10,'('' TJON DEUT ELASTIC FILESS TO BE OPEN='',
     >      10A,10A,10A)')    C_FILE,Q_FILE,M_FILE
           OPEN(UNIT=20,FILE=M_FILE,STATUS='OLD')
           OPEN(UNIT=21,FILE=C_FILE,STATUS='OLD')
           OPEN(UNIT=22,FILE=Q_FILE,STATUS='OLD')
C
           DO J=1,100
              READ(20,*,END=9) Q_2(J),FM(J)
              READ(21,*) Q_2(J),FC(J)
              READ(22,*) Q_2(J),FQ(J)
              N_D_PTS = J
           ENDDO
C
9          DO ILUN = 20,22
              CLOSE(UNIT=ILUN)
           ENDDO
       first=.false.
      endif
 
 
C
C----------------------------------------------------------------------
C     Calculate some kinematical quantities: UNITS ARE MEV
C----------------------------------------------------------------------
 
      QMU2 =    1.E6 *QSQ ! change from GeV**2 to MeV**2
      ETA     = QMU2/(4.D0 * MD*MD)
 
C----------------------------------------------------------------------
C     Get deuteron form factors
C----------------------------------------------------------------------
 
      QMU2_FM  = QMU2/HBARC**2   ! change to inverse fermis.
      G_C = RINTERPQ(Q_2,FC,QMU2_FM,N_D_PTS)
      G_Q = RINTERPQ(Q_2,FQ,QMU2_FM,N_D_PTS)*(MD/HBARC)**2
      G_M = RINTERPQ(Q_2,FM,QMU2_FM,N_D_PTS)
 
      GC2 = G_C**2
      GM2 = G_M**2
      GQ2 = G_Q**2
      A   = GC2 + (2.D0*ETA/3.D0) * GM2 + (8.D0*ETA*ETA/9.D0) * GQ2
      B   = (4.D0*ETA/3.D0) * (1.D0+ETA) * GM2
      RETURN
      END
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Linear Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming a linear
C       dependence:  Y = AX + B
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERP(X,Y,X0,NBINS)
      IMPLICIT NONE
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B
      INTEGER NBINS,I1/0/,I2/0/,I
C
C$$      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
C$$     #    WRITE(10,'('' Extrapolating outside range: X='',G10.3)') X0
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-1
         I2 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
1     DET = X(I1)-X(I2)
      A = (Y(I1)-Y(I2))/DET
      B = (X(I1)*Y(I2)-X(I2)*Y(I1))/DET
C
      RINTERP = A*X0 + B
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Quadratic Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming a quadratic
C       dependence:  Y = AX^2 + BX + C
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPQ(X,Y,X0,NBINS)
      IMPLICIT NONE
C
      DOUBLE PRECISION X(100),Y(100),X0,DET,A,B,C
      INTEGER NBINS,I1/0/,I2/0/,I3/0/,I
C
C$$      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
c$$     #    WRITE(10,'('' Extrapolating outside range: X='',G10.3)') X0
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
1     DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y(I1)*(X(I2)-X(I3)) - X(I1)*(Y(I2)-Y(I3)) + Y(I2)*X(I3)
     #        - X(I2)*Y(I3) )/DET
      B = ( X(I1)**2*(Y(I2)-Y(I3)) - Y(I1)*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y(I3) - X(I3)**2*Y(I2) )/DET
      C = ( X(I1)**2*(X(I2)*Y(I3)-X(I3)*Y(I2))
     #        - X(I1)*(X(I2)**2*Y(I3)-X(I3)**2*Y(I2))
     #        + Y(I1)*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPQ = A*X0**2 + B*X0 + C
C
      RETURN
      END
 
 
 
 
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Exponential Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming the exponential
C       form:  Y = EXP(AX^2 + BX + C)
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPEXP(X,Y,X0,NBINS)
      IMPLICIT NONE
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B,C,Y1,Y2,Y3
      INTEGER NBINS,I1/0/,I2/0/,I3/0/,I
C
C$$      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
C$$     #    WRITE(10,'('' Extrapolating outside range: X='',G10.3)') X0
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C     If all three Y-values are > 0, perform quadratic interpolation
C     on their logarithms; otherwise return zero as the result.
C ----------------------------------------------------------------------
C
1     IF(Y(I1).GT.0.D0 .AND. Y(I2).GT.0.D0 .AND. Y(I3).GT.0.D0) THEN
         Y1 = LOG(Y(I1))
         Y2 = LOG(Y(I2))
         Y3 = LOG(Y(I3))
      ELSE
 
         WRITE(10,'('' RiNTERPEXP:non-positive y-value; set to 0'')')
         RINTERPEXP = 0.D0
         RETURN
      ENDIF
C
      DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y1*(X(I2)-X(I3)) - X(I1)*(Y2-Y3) + Y2*X(I3)
     #        - X(I2)*Y3 )/DET
      B = ( X(I1)**2*(Y2-Y3) - Y1*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y3 - X(I3)**2*Y2 )/DET
      C = ( X(I1)**2*(X(I2)*Y3-X(I3)*Y2)
     #        - X(I1)*(X(I2)**2*Y3-X(I3)**2*Y2)
     #        + Y1*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPEXP = EXP(A*X0**2 + B*X0 + C)
C
      RETURN
      END
 
 
 
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     2-Dimensional Linear Interpolation Routine:
C
C       Calculates F(X0,Y0) given array F(X,Y) by two successive
C       interpolations, first in X and then in Y.
C
C       Assumes uniform spacing of array in X and Y.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY,DELX,DELY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200),Y(200),F(200,200)
C
      I = DINT((X0-X(1))/DELX) + 1
      J = DINT((Y0-Y(1))/DELY) + 1
      IF((I+1.GT.NX).OR.(I.LT.1).OR.(J+1.GT.NY).OR.(J.LT.1))THEN
        RINTERP2D = 0.D0
        RETURN
      ENDIF
C
      RINTX1 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J)-F(I,J))
     #                  + F(I,J)
      RINTX2 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J+1)-F(I,J+1))
     #                  + F(I,J+1)
      RINTERP2D = ((Y0-Y(J))/(Y(J+1)-Y(J)))*(RINTX2-RINTX1) + RINTX1
C
      RETURN
      END

                          
!-----------------------------------------------------------------------

      REAL*8 FUNCTION ATAILFL_139(ES,EP,SNSQ,TB,TA,ZP,ZN,A,IFLG)  
      IMPLICIT NONE
      REAL*8 ES,EP,SNSQ,TB,TA,ZP,ZN,A
      INTEGER IFLG
      REAL*8 FUNXL_139,elastcl_139,PHIS,PHIP,BTR,PEAKS,PEAKP,COEFF,
     >       VAC_TERM_S,VAC_TERM_P,FACP,FACS,PART1,PART2,PART3,ABREM,
     >       ATAIL,COST,QUADMO,TM,Q2,XI,OMP,OMS,Q2S,Q2P,VS,VP,ETA,DLZ,
     >       B, VERTEX_TERM_S, VERTEX_TERM_P,TOT
      INTEGER NLVL
      REAL*4 BREMS_139,DELVAC
      real*8 DDILOG
      EXTERNAL FUNXL_139
      COMMON/TAL1/ESX,EPX,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1  
      REAL*8  ESX,EPX,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1  
      COMMON/TAL2/ZPX,ZNX,AX,IFLGX
      REAL*8 ZPX,ZNX,AX
      INTEGER IFLGX
      REAL*8 ETAIL,QTAIL,QE_SMEAR_PEAK,QE_DELTA_EXACT,QE_DELTA_PEAK, 
     >       TARGET_RADIATION,QE_SOFT                                           
      COMMON/WHICHT/ ETAIL, QTAIL, QE_SMEAR_PEAK, QE_DELTA_EXACT,               
     >               QE_DELTA_PEAK,TARGET_RADIATION,QE_SOFT                     
      REAL*8 EMSQ/.26113D-6/,PI/3.1415926535D0/                           
      REAL*8 ALPHA/7.2973515D-3/,TWO_ALPH_PI/.00464564D0/                                            
      common/testing/ prttst,usegrd
      logical prttst,usegrd

C------------------------------------------------------------------------
C     CALCULATES THE RADIATIVE TAIL OF ELASTIC OR QUASI-ELASTIC         
C     SCATTERING USING THE EXACT EXPRESSION (TSAI)                      
C                                                                       
C     GIVEN ES = INCIDENT ENERGY IN GEV                                 
C           EP = SCATTERED ENERGY IN GEV                                
C           SNSQ OF THE SCATTERING ANGLE                                
C           TB = RADIATOR "BEFORE" IN RADIATION LENGTHS                 
C           TA = RADIATOR "AFTER"  IN RADIATION LENGTHS                 
C           ZP = ATOMIC NUMBER (NO. OF PROTONS)                         
C           ZN = AVERAGE NUMBER OF NEUTRONS                             
C           A= ATOMIC WEIGHT (CARBON-12 SCALE) H=1.00797                
C           IFLG = 0 FOR ELASTIC, 1 FOR QUASI-ELASTIC                   
C                                                                       
C     ANSWER IS IN PB/(GEV-SR)                                          
! Version of ATAILFL, originally from E139/E140 code, 
! adopted for use in EXTERNAL 7/25/96
!------------------------------------------------------------------------------
      ATAILFL_139=0.D0                                                      
      IF (EP.LE.0.D0) RETURN                                            
C                                                                       
C     START FILLING THE COMMON BLOCK FOR THE INTEGRAND                  
      EPX=EP                                                            
      ESX=ES                                                            
      ZPX=ZP                                                            
      ZNX=ZN                                                            
      AX=A                                                              
      IFLGX=IFLG                                                        
C                                                                       
C     CALCULATE TARGET DEPENDENT PARAMETERS                             
C                                                                       
C *** NOTE FACTORS OF 1440 & 183 HAVE BEEN CHANGED TO 1194 & 184.15     
C     TO AGREE WITH TSAI, REV. MOD. PHYS. 46, 4 (1975), ON P. 828       
C                                              25 JUNE 83 AFS           
C                                                                       
      DLZ=DLOG(184.15D0/ZP**(1.D0/3.D0))                                
      ETA=DLOG(1194.D0/ZP**(2.D0/3.D0))/DLZ                             
      B=4.D0/3.D0*(1.D0+(ZP+1.D0)/(9.D0*(ZP+ETA)*DLZ))                  
      XI=.11D0*(TB+TA)/((ZP+ETA)*DLZ)                                   
      TM=.93828D0                                                      
      IF (IFLG.EQ.0) TM=A*TM/1.00797D0                                  
      TMSQ=TM**2                                                        
C                                                                       
C     CALCULATE KINEMATICS                                              
      Q2=4.D0*ES*EP*SNSQ                                                
      OMS=ES-EP/(1.D0-2.D0*EP*SNSQ/TM)                                  
      OMP=ES/(1.D0+2.D0*ES*SNSQ/TM)-EP                                  
      IF (OMP.LE.0.) RETURN                                             
      Q2S=4.D0*(ES-OMS)*EP*SNSQ                                         
      Q2P=4.D0*ES*(EP+OMP)*SNSQ                                         
      VS=OMS/ES                                                         
      VP=OMP/(EP+OMP)                                                   
C                                                                       
C     CALCULATE QED TERMS AND MULT. PHOTON NORMALIZATION                
      TOT=TB+TA                                                         
C Second order term in B*TOT added 5/28/87 5:07 PM- Steve Rock          
      FAC1=(1.D0+0.5772D0*B*TOT- .66*(B*TOT)**2)                        
     >   -ALPHA/(2.D0*PI)*DLOG(ES/EP)**2                                
Cxxx dropped this term temporarily to see if can get code to work
     *+ALPHA/PI*(PI**2/6.D0-DDILOG(1.D0-SNSQ))                           
C     FACS=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2S/EMSQ))  
C     FACP=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2P/EMSQ))  
                                                                        
                                                                        
                                                                        
C*******8/19/87**  Steve Rock ***************************************   
C Use the Bardin Calculation of Vertex terms instead of only electron   
C  and Tsai's Vertex term.                                              
C     FAC2=2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(-QSQ/EMSQ))      
                                                                        
C  Vacuum terms of Bardin including 3 leptons and quarks                
      VAC_TERM_S= 2. * DELVAC(SNGL(Q2S)) !DELVAC was implicitly real*8 before
      VAC_TERM_P= 2. * DELVAC(SNGL(Q2P))                                
                                                                        
C  Vertex term of Tsai (eq. 2.11)                                       
      VERTEX_TERM_S = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG(Q2S/EMSQ) )      
      VERTEX_TERM_P = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG(Q2P/EMSQ) )      
C7/29/96      FAC2 = VAC_TERM + VERTEX_TERM                                     
C     FACS=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2S/EMSQ))  
C     FACP=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2P/EMSQ))  
      FACS = FAC1 + VAC_TERM_S + VERTEX_TERM_S                          
      FACP = FAC1 + VAC_TERM_P + VERTEX_TERM_P                          
                                                                        
C*******************************************************************    
                                                                        
                                                                        
C     EQUIVALENT RADIATOR FOR INTERNAL                                  
      BTR=ALPHA/PI*(DLOG(Q2/EMSQ)-1.D0)                                 
C                                                                       
C     BREMSSTRAHLUNG SPECTRUM                                           
C$$   PHIS=1.D0-VS+0.75D0*VS**2                                         
C$$   PHIP=1.D0-VP+0.75D0*VP**2                                         
                                                                        
C Installed 6/4/87 to get better value of bremstrung spectrum using     
C  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C     This is the Old-E139.E140 version of BREMS which is different than
c     the one traditionally in External.  SER-7/25/96
      PHIS= BREMS_139(VS,ZP)                                                
      PHIP= BREMS_139(VP,ZP)                                                
C                                                                       
C     S-PEAK                                                            
      PEAKS=(TM+2.D0*(ES-OMS)*SNSQ)/(TM-2.D0*EP*SNSQ)*FACS*             
     * ELASTCL_139(ES-OMS,SNSQ,ZP,ZN,A,IFLG)*(B*TB/OMS*PHIS                 
     * +XI/(2.D0*DMAX1(OMS,10.D0*XI)**2))                               
C                                                                       
C     P-PEAK                                                            
      PEAKP=FACP*                                                       
     * ELASTCL_139(ES,SNSQ,ZP,ZN,A,IFLG)*(B*TA/OMP*PHIP                     
     * +XI/(2.D0*DMAX1(OMP,10.D0*XI)**2))                               
C                                                                       
C     CALCULATE QUANTITIES INDEPENDENT OF PHOTON ANGLE                  
      SVEC=DSQRT(ES**2-EMSQ)                                            
      PVEC=DSQRT(EP**2-EMSQ)                                            
      COST=1.D0-2.D0*SNSQ                                               
      SP=ES*EP-SVEC*PVEC*COST                                           
      USQ=2.D0*EMSQ+TMSQ+2.D0*TM*(ES-EP)-2.D0*SP                        
      U0=ES+TM-EP                                                       
      UVEC=DSQRT(U0**2-USQ)                                             
      COSP=(SVEC*COST-PVEC)/UVEC                                        
      COSS=(SVEC-PVEC*COST)/UVEC                                        
C                                                                       
C     INTEGRATE OVER COS(THETAK) IN THREE PIECES                        
C     MAKING CERTAIN THAT THE S AND P PEAKS ARE FOUND                   
                                                                        
      PART1=QUADMO(FUNXL_139,-1.D0,COSP,1.D-4,NLVL)                         
      PART2=QUADMO(FUNXL_139,COSP,COSS,1.D-4,NLVL)                          
      PART3=QUADMO(FUNXL_139,COSS,1.D0,1.D-4,NLVL)                          
                                                                        
C$$   PART1=DCADRE(FUNXL,-1.D0,COSP,0.D0,1.D-4,ERROR,IER)               
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART1 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
C$$   PART2=DCADRE(FUNXL,COSP,COSS,0.D0,1.D-4,ERROR,IER)                
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART2 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
C$$   PART3=DCADRE(FUNXL,COSS,1.D0,0.D0,1.D-4,ERROR,IER)                
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART3 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
                                                                        
                                                                        
      COEFF=389.44D6*ALPHA**3*TM/PI*EP/ES                               
      ATAIL=COEFF*(PART1+PART2+PART3)                                   
C                                                                       
C     PUT IT ALL TOGETHER                                               
      ABREM=PEAKS+PEAKP                                                 
      QE_SOFT=VS**(B*TB+BTR)*VP**(B*TA+BTR)                             
                                                                        
C     TARGET_RADIATION=TARCOR(ES,EP,SNSQ,ZP,A,IFLG)                     
      TARGET_RADIATION= 1.                                              
                                                                        
      ATAILFL_139=(ATAIL*TARGET_RADIATION+ABREM)*QE_SOFT 
      if(prttst) then
        write(*,111) tb,ta,COEFF,PART1,PART2,PART3
 111    format(/1x,' tb,ta,co,p1,p2,p3,at=',2f7.4,5e10.2)
        write(*,112) iflg,atail,ATAILFL_139,abrem,qe_soft,peaks,peakp
 112    format(1x,' iflg,at,ab,qe,ps,pb=',i2,6e11.3)
      endif
C                                                                       
C DEBUG:                                                                
C                                                                       
c      WRITE(*,1000) COEFF,PART1,PART2,PART3,ATAIL,ABREM,QE_SOFT  
c1000  FORMAT(' COEFF,PART1,PART2,PART3,ATAIL,ABREM,QE_SOFT:',    
c     >       / 7(1PE14.5) )                                             
c      WRITE(*,2000) ATAILFL_139                                             
2000  FORMAT(' ATAILFL =',1PE14.5)                                      
C                                                                       
      RETURN                                                            
      END                                                               

                          
!-----------------------------------------------------------------------
      SUBROUTINE Q_E_VANORDEN(Q2_ELAS_GEV,E_GEV,SUPPRESSION)
!---------------------------------------------------------------------------
C   This program compute the quasi-elastic cross section
C   based on the Van Orden calculation using the fermi gas model
!  input energy now in GeV
! It returns the Suppression factor for the quasi-elastic peak.
!-----------------------------------------------------------------------------

      IMPLICIT NONE
      REAL E_GEV,Z,A,KF_GEV,SUPPRESSION
      REAL SUM,AP,FMT,FACTOR,N,MT,TH,SINSQ,FMOTT,WMOTT,TANSQ,EP_ELAS,
     >  EF,RLR,RTR,SRL,RATIO,W,WBAR,QV,Q2BAR,F2,XETA,ZETA,T1,T2,E1,
     >  E,SL,ST,RS,RSSIG,XSECT,CROSS,QT,C,Q2_ELAS,E2,W0,W1,Q2,EINC,KF,
     >  KF3,Z1,A1,QV2,WSQ,MT4PI,MP2,Q2_ELAS_GEV
      INTEGER IOPT,I0,I00,I,IOPT1
      REAL MUP2/7.784/,MUN2/3.65/
      REAL EBAR/1./,WDEL/2./ !$$ are these OK?
      REAL  ALPH/.0072972/, PI/3.1415927/,AMASS/931.5/,MP/939./
      LOGICAL FIRST/.TRUE./
!-----------------------------------------------------------------------------------------
! In this notation, W=nu
!---------------------------------------------------------------------------------------     
      EINC = 1000.*E_GEV
      Q2_ELAS = 1.E6 * Q2_ELAS_GEV

C$$      OPEN(7,FILE='Q_E_VANORDEN2')
      SUM=0.

      IF(IOPT.EQ.1)THEN
 1     EP_ELAS = EINC - Q2_ELAS/(2.*MP)
       IF(EP_ELAS.GT.0) THEN
        SINSQ = Q2_ELAS/(4.*EINC*EP_ELAS)
        TH = 2.*ASIN(SQRT(SINSQ))
       ELSE
        EINC = EINC +5.
        GO TO 1
       ENDIF

       FMOTT=ALPH*ALPH*COS(TH/2.)**2/4./SINSQ/SINSQ*197.3**2*1.E-26
       WMOTT=FMOTT/EINC**2
       TANSQ=TAN(TH/2.)**2
       QT = sqrt(Q2_ELAS**2/(4.*MP2) + Q2_ELAS)
       IF(QT.GT. 2.*KF) THEN
        SUPPRESSION=1.
        RETURN
       ENDIF
       W0=MAX(EINC-(EP_ELAS+2.*KF),2.)
       W1=MAX(EINC-(EP_ELAS-2.*KF),0.)
C$$       WRITE(7,'(/''E_INC,THET,WDEL,FERMI_MOM='',1P9E10.2)')
C$$     >   EINC,THETA,WDEL,KF
C$$       WRITE(7,'(''EBAR,Z,A,Q2_ELAS,W0,W1='',1P10E9.2)')
C$$     >  EBAR,Z,A,Q2_ELAS,w0,w1
C$$       WRITE(7,'(/''        NU          SIGMA(nb)'')')
      ENDIF
      
      RLR=0.
      RTR=0.
      SRL=0.
C      QV=0.
      RATIO=1.
C
      W=W0
      I00=W0/WDEL
      I0=W1/WDEL
      DO 17 I=I00,I0
       WBAR=W-EBAR
       IF(WBAR.LE.0.)GO TO 15
       IF(IOPT.EQ.1)Q2=4.*EINC*(EINC-W)*SINSQ
c       WSQ = W**2
       QV2 = Q2+WSQ
       QV=SQRT(QV2)
       Q2BAR=QV2-WBAR**2
       E1=SQRT(QV2*MP2/Q2BAR+QV2/4.)-WBAR/2.
       IF(EF.LT.E1)GO TO 18     ! do not calculate 
       RATIO=Q2/(Q2+WSQ)
       F2=1./(1.+Q2/855.**2)**4
       XETA=Q2/(4.*MP2)
       ZETA=1.+XETA
!       T1=F2*Q2/2.*(((1.+2.79*XETA)/ZETA+(1.79/ZETA))**2+N/Z*3.65)
!         T1 = 2Mp**2 * DIPOLE *Tau*(MuP2 + N/Z * MuN2) = Tau*(Gmp**2 + Gmn**2 )
!        = 2Mp*Tau*(F1+(Mu-1)F2)**2 where
!          F1= DIPOLE*(1+TAU*Mu)/(1+Tau)    F2= DIPOLE/(1+Tau)
!       T2=2.*MP22*(((1.+2.79*XETA)/ZETA)**2+XETA*((1.79/ZETA)**2+
!     >  N/Z*1.91**2))*F2  !$$$*** I think the neutron term should be divided by ZETA
!         T2=2Mp**2 *DIPOLE* (Gep +Tau*Gmp)/(1+Tau)  + neutron
!           =2Mp**2 * (F1**2 +Tau*(MuP-1)F2**2)
     
! Below is Steve's Redoing
       T1 =F2*Q2/2.*(MUP2 +N/Z*MUN2)
       T2 =2.*MP2*F2*
     >   ( (1.+MUP2*XETA)/ZETA +N/Z* (0.+MUN2*XETA)/ZETA)
       E2=EF-WBAR
       E=E1
       IF(E2.GT.E1)E=E2

       RLR=(.75*Z/(KF3*QV))*(T2*((EF**3-E**3)/3.+W*(EF**2-E**2)/2.
     >  +WSQ*(EF-E)/4.)/MP2-QV2*T1*(EF-E)/Q2)

       RTR=(.75*Z/(KF3*QV))*(2.*T1*(EF-E)+T2*(Q2BAR*(EF**3-E**3)/
     > (3.*QV2)+Q2BAR*WBAR*(EF**2-E**2)/(2.*QV2)-
     > (Q2BAR**2/(4.*QV2)
     > +MP2)*(EF-E))/MP2)
15     CONTINUE
       SL=RLR*MT4PI
       ST=RTR*MT4PI

       RS=RATIO*RATIO*SL+(0.5*RATIO+TANSQ)*ST
       RSSIG=WMOTT*RS
       XSECT=FACTOR*WMOTT*RS
       IF(SL.EQ.0.AND.ST.EQ.0.)GO TO 18
C       SRL=SRL+RLR
       IF(IOPT.EQ.1)THEN
        SUM = SUM + XSECT* WDEL
C$$        WRITE(7,35)W,XSECT
       ELSE
C$$        WRITE(7,36)W,RLR,RTR
       ENDIF
35     FORMAT(1X,1P1E15.2,1X,1P1E15.2)
36     FORMAT(1X,F6.2,2E15.8)
 18    W=W+WDEL
 17   CONTINUE

      F2=1./(1.+Q2_ELAS/855.**2)**4
      XETA=Q2_ELAS/(4.*MP2)
      ZETA=1.+XETA
 
      CROSS =WMOTT* F2* 1.E33 *
     >         ( Z*( (1.+MUP2*XETA)/ZETA +2*TANSQ*MUP2 *XETA)+
     >           N*( (0.+MUN2*XETA)/ZETA +2*TANSQ*MUN2 *XETA))
C$$   WRITE(7,'(''SIGMA: SUM='',1PE9.2,'';  PEAK='',1PE9.2)')SUM,CROSS
      SUPPRESSION = SUM/CROSS
! Tsai
      QT = sqrt(Q2_ELAS**2/(2.*938.28)**2 + Q2_ELAS)
      IF(QT.GT. 2.*KF) THEN
       C=1.
      ELSE
       C = .75*QT/KF *(1.-.083333*(QT/KF)**2)
      ENDIF
C$$      WRITE(7,'('' SUPPRESSION: TSAI='',F5.3,'':   Meziani='',F5.3)')
C$$     >   C,sum/cross
      RETURN


      ENTRY  Q_E_VANORDEN_INIT(Z1,A1,KF_GEV,IOPT1)  
       Z = Z1
       A = A1
       KF= 1000.*KF_GEV
       KF3 = KF**3
       IOPT=IOPT1
       AP=ALPH/PI
       FMT=A*AMASS
       FACTOR= 4.0*PI/FMT *1.E33
       N=A-Z
       MP2 = MP**2
       MT=931.5*(Z+N)
       EF=SQRT(KF**2+MP2)
       MT4PI = MT/4./PI
      RETURN
      END
                          
!-----------------------------------------------------------------------
      SUBROUTINE F2_JOIN(QQ,X,amuM,iZ,INEL_MODEL,W1,W2,iA) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon by joining the resonance
!  model with the DIS model
! amuM=atomic weight
! INEL_MODEL =rrdd where rr is resonance model number and dd is DIS model 
!             It must be >=100
! DIS_MODEL= 1 ineft
!            2 f2nmc
!            3 f2nmc95
!            4 SMC
!            5 E665
!            6 f2allm (HERA fit, obtained from Antje Burrel
!            9 f2glob model 9
!           12 f2glob model 12
!           33 rock_res9
! RES_MODEL  1 ineft
!            2 H2Model
!            3 rock_res9
!            5 E665
!            
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
! EMC effect corrected for neutron excess 8/19/98
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,ERR_LO,ERR_HI,NU,WSQ,f2allm
!      real fnp_nmc,corr
      real w1model,w2model,w1in,w2in,fnp_nmc
      INTEGER IZ,iA
      REAL MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC_N,FRAC,amuM,W1_DIS,W2_DIS,DSLOPE,ST,SLOPE,SY
      real f2nf2p,emcfacp,ROCK_RES9
      REAL*8 W18,W28,QQ8,W8,amuM8,F28,RDUM
      INTEGER INEL_MODEL,IHD,RES_MODEL,DIS_MODEL
      LOGICAL GD
      CHARACTER*1 HD(2)/'H','D'/
!** changed from 0.3 to 0.0 to avoid discontinuities. H2model
!** is well enough behaved that this should be fine
      REAL Q_RES_LIMIT/0.0/  
      logical ReallyJoin
      common/testing/prttst,usegrd
      logical prttst,usegrd

      RES_MODEL =INEL_MODEL/100
      DIS_MODEL =INEL_MODEL -RES_MODEL*100

      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      QQ8=QQ
      W8= SQRT(max(0.,WSQ))
      amuM8 =amuM 
      IHD =MIN(IA,2)   !no more than deuteron

      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC_N(X,REAL(IA),REAL(IZ),GD)) !with neutron excess 8/19/98
!*** If using F2ALL, use R1990 because this used to get these F2 values
!*** changed 1/19/02
      IF(DIS_MODEL.EQ.6) THEN
        CALL R1990F(X,QQ,R,DR)
      ELSE
        CALL R1998(X,QQ,R,DR,GD)
      ENDIF

!* If A>2 and doing Anje fit, only use F2ALL
      ReallyJoin=.true.
      if(IA.gt.2.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 
!* If RES_MODEL and DIS_Model both 6, don't join
      if(RES_MODEL.eq.6.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 

!** Don't join if using rock_res9: should be good everywhere
      if(DIS_MODEL.EQ.33) ReallyJoin=.false.

      IF(WSQ.GT.3.5.or.(.not.ReallyJoin)) THEN        !DIS
	IF(DIS_MODEL.EQ.1) THEN   !Bodek fit  REAL*8      
          CALL INEFT(QQ8,W8,W18,W28) !already has EMC correction
          W1=W18
          W2=W28
        ENDIF
        IF((DIS_MODEL.GE.2.AND.DIS_MODEL.LE.6).or.
     >     DIS_MODEL.EQ.33 ) THEN
         IF(DIS_MODEL.EQ.2) CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
         IF(DIS_MODEL.EQ.3) CALL F2NMC_NEW(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.4) CALL F2SMC98(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.5)  THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          F2 = F28
         ENDIF
         IF(DIS_MODEL.EQ.33) F2=ROCK_RES9(ihd,x,qq)
         W2 = F2/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF
!** Changed 1/19/02 to use Anje's n/p fit and nuclear dep.
!** Changed 7/20/06 back to NMC for n/p (works MUCH better!)
        IF(DIS_MODEL.EQ.6) THEN  !(HERA fit, obtained from Antje Burrel
          F2 = F2ALLM(X,QQ)   ! proton fit
          IF(iHD.EQ.2)  F2 = (FNP_NMC(X,QQ)+1.)/2. *F2
! This is what Antje uses:
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))
c          IF(IHD.EQ.2)  F2 = (F2NF2P+1.)/2. *F2
!* For C, Al, Cu, Au, override EMCFAC
          emcfacp=emcfac
          if(ia.eq.12) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** use carbon for Al since Antje's seems to low at high x
          if(ia.eq.27) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** turn off Antje's Al corr., 
!         if(ia.eq.27) emcfacp=(0.825-0.46*x-0.19*exp(-21.8*x)+             !!!ALUMINUM  !!!
!    >                    (0.34*x**(-4.91))*x**5.0)
          if(ia.eq.65) emcfacp=(1.026-0.56*x-0.34*exp(-45.7*x)+             !!!COPPER   !!!
     >                    (0.26*x**(-4.41))*x**5.0) 
          if(ia.eq.197) emcfacp=(0.970-1.433*x-0.334*exp(-54.53*x)+          !!!GOLD    !!!
     >                  1.074*x**0.711)
!**** changed to use Rock emc, not Antje!
!***          W2 = F2/NU *EMCFACP
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
          if(prttst) write(*,'(1x,''W2,F2,R,emcfacp='',4f10.4)')
     >      W2,F2,R,emcfacp  
        ENDIF

        IF((DIS_MODEL.EQ.9).OR.(DIS_MODEL.EQ.12)) THEN
          CALL F2GLOB(X,QQ,HD(IHD),DIS_MODEL,F2,ST,SY,SLOPE,DSLOPE,GD)
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R)
        ENDIF
        if(prttst) write(*,'(1x,''model,x,q2,F2,r='',i3,5f10.4)')
     >    DIS_MODEL,x,qq,f2,r,emcfac  
      ENDIF

      IF(WSQ.LT.4.3.and.ReallyJoin) THEN !Resonance
        W1_DIS =W1
        W2_DIS =W2
        IF(RES_MODEL.EQ.1) THEN  !old Bodek resonance fit: EMC corrected
          CALL INEFT(QQ8,W8,W18,W28)
          W1=W18
          W2=W28
        ELSEIF(RES_MODEL.EQ.5) THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          W2 = F28/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF 
        IF(RES_MODEL.EQ.2) THEN !Stuart Fit to Resonance, valid forQ2>.75
          CALL INEFT(QQ8,W8,W18,W28) !already had EMC, neutron corrected
          W1=W18
          W2=W28
          w1in = w1
          w2in = w2
c old way          IF((QQ.LT.10.).AND.(QQ.GE.Q_RES_LIMIT).AND.iA.lt.3) THEN
          IF(QQ.GE.0.3.AND.iA.lt.3) THEN
            if(ia.eq.1) CALL H2MODEL(QQ,WSQ,W1model,W2model)
            if(ia.eq.2) CALL D2MODEL_IOANA(QQ,WSQ,W1model,W2model)
            if(prttst) WRITE(*,'('' H2/D2Model,q2,wsq,w1,w2='',
     >        5f8.3)') QQ,WSQ,W1model,W2model,EMCFAC 
            if(qq.gt.0.5.and.qq.lt.15.) then
              w1=w1model
              w2=w2model
            endif
            if(qq.ge.0.3.and.qq.le.0.5) then
              frac=(qq-0.3)/(0.5-0.3)
              w1 = w1model * frac + w1in *(1.0-frac)
              w2 = w2model * frac + w2in *(1.0-frac)
            endif
            if(qq.ge.10.0.and.qq.le.15.) then
              frac=(qq - 10.0)/(15.0 - 10.0)
              w1 = w1in * frac + w1model *(1.0-frac)
              w2 = w2in * frac + w2model *(1.0-frac)
            endif
          ENDIF
        ENDIF 

        IF(WSQ.GT.3.5) THEN   !interpolate
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           if(prttst) WRITE(*,'(1x,''joining w2'',3f10.3)') 
     >       w2_dis,w2,frac
           W1 = W1_DIS*FRAC + W1*(1.0 - FRAC)
           W2 = W2_DIS*FRAC + W2*(1.0 - FRAC)
        ENDIF
      ENDIF
!**** If using Antje model, then also use Blok correction factors
!* taken out 1/19/02
!      if(DIS_MODEL.EQ.6) THEN
!        if(ihd.eq.1) call f2corr1h(x,corr)
!        if(ihd.eq.2) call f2corr2h(x,corr)
!        w1=w1*corr
!        w2=w2*corr
!      endif
      RETURN
      END


                          
!-----------------------------------------------------------------------
      SUBROUTINE NMCDIS(QQ,X,W1,W2,amuM) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon. 
!  NMC fit for DIS and Old Bodek fit for Resonance.
! IA=atomic weight
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,NU,WSQ, MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC,FRAC,amuM
      REAL*8 W18,W28,QQ8,W8,amuM8
      INTEGER IA,IHD
      LOGICAL GOOD


      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      IF(WSQ.LT.1.16) RETURN
      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC(X,amuM,GOOD)) !isoscaler nucleus
      CALL R1998(X,QQ,R,DR,GOOD)
      IA = amuM+.5  !make into integer
      amuM8 =amuM
      IF(WSQ.GT.3.5) THEN
        
        IHD =MIN(IA,2)   !no more than deuteron
        CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
        W2 = F2/NU *EMCFAC
        W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
      ENDIF

      IF(WSQ.LT.4.3) THEN
         QQ8=QQ
         W8= SQRT(max(0.,WSQ))
         IF(WSQ.LE.3.5) THEN
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit: EMC corrected
           W1=W18
           W2=W28
         ELSE  !interpolate
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit       
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + W18*(1.0 - FRAC)
           W2 = W2*FRAC + W28*(1.0 - FRAC)
         ENDIF
      ENDIF
      F2 = NU*W2
      RETURN
      END


                          
!-----------------------------------------------------------------------
      REAL FUNCTION BREMS_139 (Y,Z)                                              
                                                                        
C=======================================================================
C                                                                       
C  Returns value of Bremstrung probability multiplied by k (photon energ
C     normalized to 1 at Y=0                                            
C                                                                       
C  Formula from Tsai; Rev of Mod Physics 46,(1974)815,  equation 3.83   
C      Complete screening                                               
C                                                                       
C  Y = k/E : fraction of beam energy carried by photon                  
C  Z =     : atomic number of target                                    
C                                                                       
C      6/4/87 Steve Rock                                                
C=======================================================================
      IMPLICIT NONE
                                                                        
      REAL*8 Z,Y                                                        
      REAL LRAD,LRADP,FC,ZZ                                          
      REAL ALPHA/7.2973515E-3/                                          
      LOGICAL FIRST /.TRUE./                                                                  
                                                                        
C TSAI page 486 Table B.2  and eq. 3.67, 3.68                           
      IF(FIRST) THEN                                                                  
       IF(Z.EQ.1)THEN                                                    
        LRAD = 5.31                                                      
        LRADP= 6.144                                                     
       ELSEIF(Z.EQ.2) THEN                                               
        LRAD = 4.79                                                      
        LRADP=5.62                                                       
       ELSEIF(Z.EQ.3) THEN                                               
        LRAD = 4.74                                                      
        LRADP= 5.805                                                     
       ELSEIF(Z.EQ.4) THEN                                               
        LRAD = 4.71                                                      
        LRADP= 5.924                                                     
       ELSE                                                              
        LRAD = ALOG(184.15) - DLOG(Z)/3.                                 
        LRADP= ALOG(1194.) -2.*DLOG(Z)/3.                                
       ENDIF                                                             
                                                                        
C Get coulomb correction from Tsai eq. 3.3                              
       ZZ=  (Z*ALPHA)**2
       FC = 1.202 * ZZ - 1.0369 * ZZ**2 +                  
     >  1.008 *ZZ**3/(1.+ ZZ)
       FIRST=.FALSE.
      ENDIF                                
                                                                        
                                                                        
C  Tsai eq.(3.83) normalized to be 1 at Y=0                             
      BREMS_139 = (1 -Y +.75 *Y**2) +                                       
     >         (1.-Y)/12. *(Z+1.)/(Z *(LRAD -FC) + LRADP)               
                                                                        
      RETURN                                                            
      END                                                               
!---------------------------------------------------------------------
C     $MEMBER=QUADMO, DATE=75061922, USER=KCE
      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   
      REAL*8 FUNCT,LOWER,UPPER,EPSLON                                   
         INTEGER NLVL                                                   
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         
         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)                                                  
         FM = FUNCT(0.5*(L+R))                                          
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=DABS(ESTL)+DABS(ESTR)                                        
      ABAREA=AREA+ABAREA-DABS(EST)                                      
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO = ROMBRG /12.0D0                                           
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ELASTCL_139(E,SNSQ,ZP,ZN,A,IFLG)
      IMPLICIT NONE
      REAL*8 E,SNSQ,ZP,ZN,A
      INTEGER IFLG            
      REAL*8 TM,RECOIL,EP,SIGM
      REAL QSQ_4,W1,W2,FF,E_4
      REAL*8 ALPHA/.00729735D0/                                           
                                                                       
C     GIVEN INCIDENT ELECTRON ENERGY E IN GEV                           
C           AND SNSQ OF SCATTERING ANGLE                                
C           ZP IS ATOMIC NUMBER (NO. OF PROTONS)                        
C           ZN IS AVERAGE NUMBER OF NEUTRONS                            
C           A IS ATOMIC WEIGHT (CARBON-12 SCALE) H=1.00797              
C           IFLG = 0 FOR ELASTIC, 1 FOR QUASI-ELASTIC                   
C     RETURNS THE CROSS SECTION FOR ELASTIC OR QUASI ELASTIC SCATTERING 
C     IN PB/SR (1.E-36 CM**2/SR)                                        

      TM=.93828D0                                                      
      IF (IFLG.EQ.0) TM=TM*A/1.00797D0                                  
      RECOIL=1.D0+2.D0*E/TM*SNSQ                                        
      EP=E/RECOIL                                                       
      QSQ_4=4.D0*E*EP*SNSQ                                                
      E_4 = E
      SIGM=(19732.D0*ALPHA/(2.D0*E*SNSQ))**2*(1.D0-SNSQ)/RECOIL         
C$$      CALL ASTRUL_139(QSQ,ZP,ZN,A,W1,W2,IFLG)     
      IF(IFLG.EQ.0)CALL NUC_FORM_FACTOR(QSQ_4,W1,W2,FF)
      IF(IFLG.EQ.1) CALL QE_PEAK(QSQ_4,E_4,W1,W2)   

      ELASTCL_139=SIGM*(W2+2.D0*W1*SNSQ/(1.D0-SNSQ)) 
c      write(6,'(''elastcl139'',i2,2f6.3,e11.3,2f8.3)') iflg,e,ep,
c     >  sigm,w2,w1

      RETURN                                                            
      END                                                               


!---------------------------------------------------------------------
       SUBROUTINE F2DGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(deuteron) for E665                          C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2DHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2DHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2DHIW(QQ,NU)

c.....low W parametrization for deuteron
      LOWWF2 = F2LOWW(2,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END
       SUBROUTINE F2PGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(proton) for E665                            C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2PHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2PHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2PHIW(QQ,NU)

c.....low W parametrization for proton
      LOWWF2 = F2LOWW(1,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END
      DOUBLE PRECISION FUNCTION F2HI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** hydrogen BCDMS-like function
      DATA DPEMC /-0.1011D0,2.562D0,0.4121D0,-0.518D0,5.957D0,
     &           -10.197D0,4.695D0,
     &             0.364D0,-2.764D0,0.015D0,0.0186D0,-1.179D0,
     &             8.24D0,-36.36D0,47.76D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2HI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION F2DI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** deuterium BCDMS-like function
      DATA DPEMC /-0.0996D0,2.489D0,0.4684D0,-1.924D0,8.159D0,
     &            -10.893D0,4.535D0,
     &              0.252D0,-2.713D0,0.0254D0,0.0299D0,-1.221D0,
     &              7.5D0,-30.49D0,40.23D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2DI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)
C
      RETURN
      END
      REAL FUNCTION F2DHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(deuteron) for W>5                                *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2DDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2DI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2DI15,F2DDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2DI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2DDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2DHIW = F2

      RETURN
      END
      REAL FUNCTION F2PHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) for W>5                                  *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2PDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2HI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2HI15,F2PDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2HI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2PDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2PHIW = F2

      RETURN
      END
      REAL FUNCTION F2DDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(deuteron)                                        *
* Uses n/p parametrization of E665 and model of Donnachie and*
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
* currently sets n/p to 1                                    *
*                                                            *
* 21may94 AVK                                                *
* set n/p = 0.94                                             *
* private communication of E665 measurement at low x and Q2  *
* from P. G. Spentzouris                                     *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)
      F2 = F2*(1.0+0.94)/2.0

      IF (F2.LT.0.0) F2 = 0.0
      F2DDL = F2

      RETURN
      END
      REAL FUNCTION F2PDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)

      IF (F2.LT.0.0) F2 = 0.0
      F2PDL = F2

      RETURN
      END
      REAL FUNCTION F2DOLA(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
*                                                            *
* The form of the function is a sum of two Regge powers, and *
* is fit to photoproduction data over a wide range of energy *
* and to NMC data for Q2<10 Gev2                             *
* Only the simple function represented by Eqns 4a and 4b are *
* used here, the refinements discussed in the paper are      *
* omitted. F2(neutron) is also omitted.                      *
*                                                            *
* units of Q2ext and Nuext are GeV^2, GeV                    *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,A,SMALLA,B,SMALLB,POWER1,POWER2,TWOM,F2
      PARAMETER (A=0.324,B=0.098,SMALLA=0.562,SMALLB=0.01113)
      PARAMETER (POWER1=0.0808,POWER2=-0.4525)
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      REAL TERM1,TERM2,TWOMNU
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

      TWOMNU = TWOM*NU

      TERM1= A * TWOMNU**POWER1 / (Q2+SMALLA)**(POWER1+1.0)

      TERM2= B * TWOMNU**POWER2 / (Q2+SMALLB)**(POWER2+1.0)

      F2 = Q2 * (TERM1 + TERM2)

      IF (F2.LT.0.0) F2 = 0.0
      F2DOLA = F2

      RETURN
      END
      REAL FUNCTION F2LOWW(ITARGE,Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*   A V Kotwal 30dec93                                       *
*                                                            *
*Return F2 for 1.1<W<5 over a large range of Q2 down to Q2=0 *
*                                                            *
* Based on parametrizations of SLAC,DESY and Daresbury data. *
* Calls two other routines for the resonance region 1.1<W<2  *
* and the inelastic region 2<W<5                             *
* Assume proton and deuteron same for resonance region       *
*                                                            *
* ITARG = 1 for H2, 2 for D2                                 *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
*                                                            *
*     PGS Feb 20, 94                                         *
* Proton and neutron are not the same, see Close chapter 7.2 *
* Modify deuteron resonances accordingly.                    *
* Simple scaling of the D resonance parametrization          *
* according to data (flat line fit in 8 points) by 0.88      *
* Data from J. Franz et al Z.Phys.C 96 P. and F. 10,         *
* 105-116(1981)                                              *
* fit a flat line in the H/D plots p. 110                    *
**************************************************************
      REAL M2,TWOM,F2,Q2,NU,W2,W,WMIN,WINTER,WMAX,NUMAX,WFUNC
      PARAMETER (M2=0.8803,TWOM=1.87654) ! protonmass**2, 2*protonmass
      PARAMETER (WMIN=1.12,WINTER=1.95,WMAX=5.0)
      INTEGER ITARGE,ITARG
      REAL F2LONU,F2PRES,LONUF2,PRESF2,WDIFF
      EXTERNAL F2LONU,F2PRES
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      ITARG = ITARGE

      W2 = M2 + TWOM*NU - Q2
      IF (W2.LT.M2)  W2=M2
      W  = SQRT( W2 )

      IF (NU.LT.0.001) NU = 0.001
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (ITARG.GT.2) ITARG = 2
      IF (ITARG.LT.1) ITARG = 1

      IF (W.LT.WMIN) THEN
c.......force to zero smoothly as W->M
        PRESF2 = F2PRES(Q2,WMIN)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        F2 = PRESF2 * ((W2-M2)/(WMIN**2-M2))**2
      ELSEIF (W.LE.WMAX) THEN
c.....inelastic region, joined smoothly to resonance region
        WDIFF = (W-WINTER)/0.03
        IF (WDIFF.GT.10.0) THEN
          WFUNC = 0.0
        ELSEIF (WDIFF.LT.-10.0) THEN
          WFUNC = 1.0
        ELSE
          WFUNC = 1.0/(1.0+EXP(WDIFF))
        ENDIF
        PRESF2 = F2PRES(Q2,W)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        LONUF2 = F2LONU(ITARG,Q2,NU)
        F2=WFUNC*PRESF2+(1.0-WFUNC)*LONUF2
      ELSE
C.......force to zero at high W
        NUMAX = (WMAX**2-M2+Q2)/TWOM
        LONUF2 = F2LONU(ITARG,Q2,NUMAX)
c        F2 = LONUF2 * EXP(10.0*(WMAX-W))
        F2 = LONUF2
      ENDIF

      IF (F2.LT.0.0) F2 = 0.0

      F2LOWW = F2

      RETURN
      END
       REAL FUNCTION F2LONU(ITARG,Q2EXT,NUEXT)
       IMPLICIT NONE
***************************************************************
* parametrizations for F2 at low nu    A V Kotwal 27dec93     *
* Original reference is F.W.Brasse et al, NP B39, 421 (1972)  *
* SLAC and DESY electroproduction and photoproduction         *
* data has been used in the fit. This fit can be used in the  *
* range 2<W<5 over a large Q2 range down to Q2=0              *
*                                                             *
* This parametrization is used in J. Franz et al,             *
* 'Inclusive Electron Scattering in the Low Q2 region'        *
* Z. Physics C 10 105-116 (1981)                              *
* in the range W>2, nu<6.5, 0.08<Q2<1.0                       *
*                                                             *
* ITARG=1 for H2, 2 for D2                                    *
* Q2EXT = Q2 in GeV^2                                         *
* NUEXT = Nu in GeV                                           *
***************************************************************
       REAL OMEGAW,NU,Q2,OMEGA,YOW,XPRIME,F2P,F2N,F2D,F2
       INTEGER ITARG
       REAL MW2,A2,B3,B4,B5,B6,B7
       PARAMETER (MW2=1.43,A2=0.42,B3=0.933,B4=-1.494,B5=9.021)
       PARAMETER (B6=-14.50,B7=6.453)
       REAL Q2EXT,NUEXT

       Q2 = Q2EXT
       NU = NUEXT
       IF (NU.LT.0.0) NU = 0.001
       IF (Q2.LT.0.0) Q2 = 0.001

       OMEGAW = (1.87654*NU+MW2)/(Q2+A2)
       OMEGA  =  1.87654*NU/Q2
       IF (OMEGA.LT.1.0) OMEGA = 1.0

       YOW    = 1.0 - 1.0/OMEGAW
       XPRIME = Q2/(1.87654*NU+0.8803)
       IF (XPRIME.GT.0.9999) XPRIME=0.9999

       F2P = B3+B4*YOW+B5*YOW**2+B6*YOW**3+B7*YOW**4
       F2P = OMEGAW*F2P*YOW**3/OMEGA
       IF (F2P.LT.0.0) F2P=0.0

       F2N = F2P*(1.0-XPRIME)
       IF (F2N.LT.0.0) F2N=0.0

       F2D = (F2N+F2P)/2.0

       IF (ITARG.EQ.1) THEN
         F2=F2P
       ELSEIF (ITARG.EQ.2) THEN
         F2=F2D
       ELSE
         F2=0.0
       ENDIF

       IF (F2.LE.0.0) F2 = 0.0

9999   CONTINUE
       F2LONU=F2
       RETURN
       END
       REAL FUNCTION F2PRES(Q2EXT,WEXT)
       IMPLICIT NONE
***********************************************************
* returns F2(proton) in the resonance region 1.1<W<2.0    *
* paremetrization of the cross-section is obtained from   *
* F.W.Brasse et al,                                       *
* ' Parametrization of the Q2 Dependence of the           *
* Gamma_v-P Total Cross-Sections in the Resonance Region' *
* Nuclear Physics B110, 413 (1976)                        *
*                                                         *
* for high energy muon scattering, the virtual photon     *
* polarization is close to 1 for low W, low Q2 scattering *
* such as resonance excitation.                           *
* Hence the parametrization for epsilon>0.9 is used from  *
* this reference. The average epsilon for this data is    *
* 0.957 . This is used together with R=sigmaL/sigmaT to   *
* convert the cross-section to F2                         *
*                                                         *
* A V Kotwal 28dec93                                      *
* inputs: Q2EXT = Q2 in GeV^2                             *
*         WEXT  = W  in GeV                               *
***********************************************************
       REAL M2,TWOM,F2P,Q2,NU,W2,W,Q,Q0,WLO,WHI,LNQOQ0,M,XX
       INTEGER NWBIN,IWBIN,N
       PARAMETER (NWBIN=56)
       REAL WMIN,WINTER,WMAX,DW1,DW2
       PARAMETER (WMIN=1.11,WINTER=1.77,WMAX=1.99)
       PARAMETER (DW1=0.015,DW2=0.02)
       PARAMETER (M=0.93828)    ! proton mass
       PARAMETER (M2=M*M,TWOM=2.0*M)
       REAL EPSILON,R,SIGMA,PI2AL4,Y,Y1,Y2,GD2,STRUW2,ERR
       PARAMETER (EPSILON=0.957)
       PARAMETER (PI2AL4=112.175)  ! 4.pi**2.alpha_em (microbarn.GeV**2)
       REAL Q2EXT,WEXT
       REAL A(NWBIN),B(NWBIN),C(NWBIN)
       DATA A /
     &  5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178
     & ,6.125,5.999,5.769,5.622,5.431,5.288,5.175,5.131
     & ,5.003,5.065,5.045,5.078,5.145,5.156,5.234,5.298
     & ,5.371,5.457,5.543,5.519,5.465,5.384,5.341,5.328
     & ,5.275,5.296,5.330,5.375,5.428,5.478,5.443,5.390
     & ,5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158
     & ,5.159,5.178,5.182,5.195,5.160,5.195,5.163,5.172 /
       DATA B /
     &  0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878
     & ,1.887,1.927,2.041,2.089,2.148,2.205,2.344,2.324
     & ,2.535,2.464,2.564,2.610,2.609,2.678,2.771,2.890
     & ,2.982,3.157,3.188,3.315,3.375,3.450,3.477,3.471
     & ,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519
     & ,4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322
     & ,5.546,5.623,5.775,5.894,6.138,6.151,6.301,6.542 /
       DATA C /
     &   0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080
     & ,-0.065,-0.056,-0.065,-0.056,-0.043,-0.034,-0.054,-0.018
     & ,-0.046,-0.015,-0.029,-0.048,-0.032,-0.046,-0.084,-0.115
     & ,-0.105,-0.159,-0.164,-0.181,-0.203,-0.220,-0.245,-0.264
     & ,-0.239,-0.302,-0.299,-0.318,-0.388,-0.393,-0.466,-0.588
     & ,-0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798
     & ,-1.048,-0.980,-1.021,-1.092,-1.313,-1.341,-1.266,-1.473 /

      Q2 = Q2EXT
      W  = WEXT
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (W.LT.M) W = M
      W2 = W*W
      NU = (W2+Q2-M2)/TWOM

      IF (NU.LT.0.001) NU = 0.001

* Q = three momentum transfer to the hadronic system in the lab frame
      Q = ABS(SQRT( Q2 + NU*NU ))
      IF (Q.LT.0.01) Q = 0.01

* Q0 is the value of Q at q2=0 for the same W
      Q0 = (W2 - M2)/TWOM
      IF (Q0.LT.0.01) Q0 = 0.01

c.....find the bin of W, and bin edges
      IF (W.GE.WMAX) THEN
        IWBIN = NWBIN
        WLO   = WMAX
        WHI   = WMAX
      ELSEIF (W.GT.WINTER) THEN
        N     = INT((W-WINTER)/DW2)
        IWBIN = N + 45
        WLO   = FLOAT(N)*DW2 + WINTER
        WHI   = WLO + DW2
      ELSEIF (W.GT.WMIN) THEN
        N     = INT((W-WMIN)/DW1)
        IWBIN = N + 1
        WLO   = FLOAT(N)*DW1 + WMIN
        WHI   = WLO + DW1
      ELSE
        IWBIN = 1
        WLO   = WMIN
        WHI   = WMIN
      ENDIF

      LNQOQ0 = ALOG(Q/Q0)

c.....eqn 3
c.....parameter d is fixed d=3 in section 3.1
      Y1 = A(IWBIN) + B(IWBIN)*LNQOQ0 + C(IWBIN)*(ABS(LNQOQ0))**3
      IF (Y1.LT.0.0) Y1 = 0.0

      IF (IWBIN.LT.NWBIN) THEN
       Y2=A(IWBIN+1)+B(IWBIN+1)*LNQOQ0+C(IWBIN+1)*(ABS(LNQOQ0))**3
c.....linear interpolation inside bin
       Y = Y1 + (W-WLO)*(Y2-Y1)/(WHI-WLO)
      ELSE
       Y = Y1
      ENDIF
      IF (Y.LT.0.0) Y = 0.0

* The parametrization returns log(sigma/GD**2) where GD is the
* nucleon dipole form factor, and sigma is the total virtual photoabsorption
* cross-section
      GD2   = 1.0 / (1.0 + Q2/0.71)**4
      SIGMA = GD2*EXP(Y)

c.....R=sigmaL/sigmaT
c.....use 1990 SLAC analysis result
      XX = Q2/(TWOM*NU)
      IF (XX.GT.0.99) XX=0.99
      IF (XX.LT.0.0) XX=0.0
      CALL R1990F(XX,Q2,R,ERR)

c.....'An introduction to Quarks and Leptons', F.E.Close
c......eqn 9.44,9.46
c......also eqn 1 from reference
      STRUW2 = SIGMA*(1.0+R)*Q2 / ((1.0+EPSILON*R)*(Q2+NU*NU))
c......the Hand convention for the virtual photon flux has been followed
c......(reference 3 in this reference)
      STRUW2 = STRUW2*Q0/PI2AL4
      F2P    = NU*STRUW2
      IF (F2P.LT.0.0) F2P = 0.0

9999  CONTINUE
      F2PRES = F2P

      RETURN
      END

      SUBROUTINE R1990F(X,QQ2,R,ERR)
C
C    Ref: L.W.Whitlow SLAC Report 357 (1990)         and
C         L.W.Whitlow et al.: PL B250 (1990) 193
C
C    for Q2 < 0.35 we extrapolate R as a constant with a rough error of 100 %
C
      REAL A(3)  / .06723, .46714, 1.89794 /,
     >     B(3)  / .06347, .57468, -.35342 /,
     >     C(3)  / .05992, .50885, 2.10807 /
C
      DATA QMAX /64./
      Q2=QQ2
      IF(Q2.LT.0.35) Q2=0.35
      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))
      RLOG  = FAC/LOG(Q2/.04)
      Q2THR = 5.*(1.-X)**5
      RA   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))
      RB   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)
      RC   = C(1)*RLOG + C(2)/SQRT((Q2-Q2THR)**2+C(3)**2)
      R     = (RA+RB+RC)/3.
      IF (Q2.GE.0.35) THEN
        Q = MIN(Q2,QMAX)
        S = .006+.03*X**2
        AA= MAX(.05,8.33*X-.66)
        XLOW  = .020+ABS(S*LOG(Q/AA))
        XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)
        D1SQUARE=XLOW**2+XHIGH**2
        D2SQUARE= ((RA-R)**2+(RB-R)**2+(RC-R)**2)/2.
        D3    = .023*(1.+.5*R)
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3
        ERR    = SQRT(D1SQUARE+D2SQUARE+D3**2)
      ELSE
        ERR = R
      ENDIF
      RETURN
      END


!---------------------------------------------------------------------

       Subroutine F2SMC98(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the SMC parametrisation
! in SMC: Adeva et al (SMC) Phys. Rev. D 58, 112001 
! in Tables 12 and 13.


      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi,errd_lo,errd_hi,FNP_NMC,F2D
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2SMC98_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       CALL F2SMC98_DO(2,x4,qsq4,F2D,errd_lo,errd_hi)  !deuteron
       F2 = 2.*F2D- F2
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END

      SUBROUTINE F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)/
     >  -0.24997, 2.3963, 0.22896, 0.08498, 3.8608, -7.4143, 3.4342,  !p      
     >  -0.28151, 1.0115, 0.08415,-0.72973, 2.8647, -2.5328, 0.47477 /
      real*8 b(4,2)/ 0.11411, -2.2356, 0.03115, 0.02135,    !p
     >               0.20040, -2.5154, 0.02599, 0.01859/ 

      real*8 c(4,2)/ -1.4517, 8.4745, -34.379, 45.888,    !p
     >               -1.3569, 7.8938, -29.117, 37.657 /
              
!lower limits
      real*8 al(7,2)/
     >  -0.25196, 2.4297, 0.21913,  0.21630, 3.4645, -6.9887, 3.2771, !p  
     >  -0.28178, 1.1694, 0.09973, -0.85884, 3.4541, -3.3995, 0.86034/
      real*8 bl(4,2)/ 0.13074, -2.2465, 0.02995, 0.02039, 
     >                0.20865, -2.5475, 0.02429, 0.01760/  
      real*8 cl(4,2)/ -1.4715, 8.9108, -35.714, 47.338,
     >                -1.3513, 8.3602, -31.710, 41.106/

!upper limits
      real*8 au(7,2)/
     > -0.24810, 2.3632, 0.23643, -0.03241, 4.2268, -7.8120, 3.5822,    !p
     > -0.28047, 0.8217, 0.06904, -0.60191, 2.2618, -1.6507, 0.08909 /
      real*8 bu(4,2)/0.09734, -2.2254, 0.03239, 0.02233,
     >               0.18711, -2.4711, 0.02802, 0.01973  /
      real*8 cu(4,2)/ -1.4361, 8.1084, -33.306, 44.717,
     >                -1.3762, 7.6113, -27.267, 35.100  /

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER      
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
       Subroutine F2NMC(T,x,qsq,F2,err_f2) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in Phy Lett. b295 159-168 Proton and Deuteron F2 Structure Functions  
! in deep inelastic muon scattering   P. Amaudruz et. al.               
!                                                                       
                                                                        
                                                                        
      real a(7,2)/-0.1011,2.562,0.4121,-0.518,5.967,-10.197,4.685,      
     >            -0.0996,2.489,0.4684,-1.924,8.159,-10.893,4.535/      
      real b(4,2)/0.364,-2.764,0.0150,0.0186,                           
     >            0.252,-2.713,0.0254,0.0299/                           
      real c(4,2)/-1.179,8.24,-36.36,47.76,                             
     >            -1.221,7.50,-30.49,40.23/                             
      real Lam/0.25/                                                    
      real Qsqo/20.0/                                                   
      real log_term                                                     
      real A_x,B_x,C_x,F2,err_f2,x,qsq                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
                                                                        
      A_x=x**a(1,t)*(1-x)**a(2,t)*(a(3,t)+a(4,t)*(1-x)+a(5,t)*(1-x)**2  
     >    +a(6,t)*(1-x)**3+a(7,t)*(1-x)**4)                             
      B_x=b(1,t)+B(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = alog(qsq/lam**2)/alog(qsqo/lam**2)                     
      F2=A_x*(log_term)**B_x*(1+C_x/qsq)                                
      return                                                            
      end                                                               

!---------------------------------------------------------------------
      REAL FUNCTION FITEMC_N(X,A,Z,GOODFIT)                                         
!---------------------------------------------------------------------  
! Modified FITEMC.F with Neutron excess correction and proton=1 added 8/19/98
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! Z= number of protons
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 8/19/98 **  If proton, return 1.
! 8/19/98 **  Add neutron excess correction.
! 11/05 Modified PYB to return value at x=0.7 if x>0.7 since
!       Fermi smearing in INELASTU will account for Fermi smearing rise
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL ALPHA, C,LN_C,X,A,Z ,X_U,SIG_N_P,F_IS
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points
!Term    Coeficient     Error
      REAL*8 ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                 
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /                                  

             !
!Chisq=         22.    for 30 points
!Term    Coeficient     Error 
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        

      IF(A.LT.1.5) THEN    ! Added 8/19/98
       FITEMC_N=1.
       GOODFIT=.TRUE.
       RETURN
      ENDIF                                                                                
c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c      IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC_N  =  C *A**ALPHA    !isoscaler
      SIG_N_P = 1.-0.8*X_U
      F_IS = .5*(1.+ SIG_N_P)/(Z/A +(1.-Z/A)*SIG_N_P)
      FITEMC_N = FITEMC_N/F_IS
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o inelstu.o -c inelstu.f"; -*-
       SUBROUTINE INELSTU(X,q2,F2,W1,W2,IMOD)
! New call to Linda Stuarts fit to resonance Larry's DIS fit
! Steve Rock 10/14/93

*******************************************************************************
* This subroutine returns W1 and W2, the inelastic structure functions
*
* 2/93, LMS.
******************************************************************************
       IMPLICIT NONE

       INTEGER IMOD
       REAL X, Q2, NU, W1, W2,  MP/0.93828/, MP2/.8803/,
     >        WSQ, WW1, WW2, FRAC
       REAL   F2,R,DR,ST,SY,SLOPE,DSLOPE
       LOGICAL GOOD

       W1 = 0.0
       W2 = 0.0

       WSQ = MP2   + Q2*(1./X - 1. )
       NU = (WSQ +Q2 -MP2)/(2.*MP)
       IF(WSQ.LT.1.16) RETURN

! Get proton structure functions
       CALL R1998(X,Q2,R,DR,GOOD)
       IF(WSQ.GT.3.5) THEN
         CALL F2GLOB(X,Q2,'H',IMOD,F2,ST,SY,SLOPE,DSLOPE,GOOD)
         W2 = F2/NU
         W1 = (1.0 + NU*NU/Q2)*W2/(1.0 + R)
       ENDIF
! H2 model only officially works above q2=.75
       IF(Q2.LT.10.0.AND.WSQ.LT.4.3) THEN
         IF(WSQ.LE.3.5) THEN
           CALL H2MODEL(Q2,WSQ,W1,W2)
           F2 = NU*W2
         ELSE
           CALL H2MODEL(Q2,WSQ,WW1,WW2)
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + WW1*(1.0 - FRAC)
           W2 = W2*FRAC + WW2*(1.0 - FRAC)
         ENDIF
       ENDIF
       F2 = NU*W2
       RETURN
       END

        SUBROUTINE H2MODELLINDA(QSQ,WSQ,W1,W2)
!**** NO LONGER USED 5/16/06
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
********************************************************************************
        IMPLICIT NONE
        logical goroper/.false./

        REAL*4  WSQ, QSQ, W1, W2, R_NRES, SIG_RES, SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, DIPOLE, K, NU,
     >          TAU, PI/3.14159265/, AM/.93828/, ALPHA/0.00729735/,
     >          CONV/0.0025767/,sigroper(2)

! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.16) return
        IF(WSQ.GT.4.3.OR.QSQ.GT.10.0) THEN
          WRITE(6,'(''  H2MODEL called outside of kinematic range'')')
          RETURN
        ENDIF

! H2MOD_FIT returns transverse cross sections in units of
! microbarns/(dipole FF)**2
!        CALL H2MOD_FIT_old(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2)
        CALL H2MOD_FIT(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)

        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.

        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
        SIG_RES  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1)) *
     >    DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES
        SIGL = R_NRES*SIG_NRES(1)

! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END


      SUBROUTINE H2MOD_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J
      REAL W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, 
     >     RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER

      REAL PI/3.14159265/, ALPHA/0.00729735/, AM/.93828/,
     >     MPPI/1.07783/, MDELTA/1.2340/, MPI/0.13957/, 
     >     GAM1WID/0.0800/, GAM2WID/0.0900/, MASS1/1.5045/,
     >     ROPERWID/0.0500/, MASSROPER/1.4000/, 
     >     MASS2/1.6850/, DELWID/0.1200/, FIT(30), SQRTWDIF,
     >     XR/0.1800/, MQDEP/3.40/

      REAL RCOEF(25)/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >  -1.5416E+02,  -1.4891E+02,   2.4102E+02,   2.5823E+00, 
     >   7.1004E+00,  -8.9771E+00,   1.3744E+00,  -1.2085E+00, 
     >   1.1218E-01/

      REAL COEF(22)/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/

      REAL RERR1(200)/
     >  2.6120E+02,-9.4211E+02, 4.0844E+03, 7.4994E+02,-3.5317E+03,
     >  3.1703E+03, 2.1238E-01,-6.1568E-01, 4.1479E-01, 1.9720E-02,
     >  1.2891E-01,-4.1615E+00, 4.7246E+00, 2.8090E-03, 6.1657E-02,
     >  1.3120E+00,-9.4379E+00, 9.0902E+00,-1.3636E-03, 2.8054E-02,
     >  9.6123E-02,-1.1465E+03, 3.9099E+03,-3.0097E+03,-1.0604E+01,
     > -1.2214E+00,-8.3549E-01, 1.3696E+04, 3.9085E+03,-1.5369E+04,
     >  1.2745E+04, 2.9942E+01, 7.7268E+00, 1.0113E+01,-4.3868E+04,
     >  1.5709E+05,-3.0207E+03, 1.2809E+04,-1.1075E+04,-2.0442E+01,
     > -7.5843E+00,-1.0773E+01, 3.2597E+04,-1.2457E+05, 1.0283E+05,
     > -1.6960E-01, 5.9410E-01,-4.5486E-01,-1.0715E-02,-2.6512E-03,
     >  1.0153E-03, 6.4074E+00,-1.9189E+01, 1.3612E+01, 9.2486E-03,
     >  2.7904E-01, 6.3576E+00,-7.8552E+00, 1.5302E-02,-1.1506E-01,
     > -4.7552E-02,-1.0171E+01,-1.5884E+00, 1.6223E+01,-1.1379E-04,
     >  4.9212E-01,-2.4354E+00, 1.7921E+01,-1.7223E+01, 4.0778E-03,
     > -4.5558E-02,-1.8539E-01, 7.9930E+00,-7.1588E+01, 7.1512E+01,
     > -2.1529E-03, 1.8337E-01, 7.7590E-01, 7.3007E+02,-2.5219E+03,
     >  1.9547E+03, 6.1102E+00, 1.2970E+00,-1.3084E+00,-9.4932E+03,
     >  3.1542E+04,-2.3894E+04,-5.9583E+00, 8.1005E-02, 3.6885E-01,
     >  9.5708E+03,-2.4911E+03, 9.4342E+03,-7.7120E+03,-1.8608E+01,
     > -1.1065E+00, 6.5015E+00, 3.1755E+04,-1.1529E+05, 9.1964E+04,
     >  1.8347E+01,-2.5899E+00, 7.1169E-01,-3.2268E+04, 1.1891E+05,
     >  1.9339E+03,-7.7737E+03, 6.6128E+03, 1.3392E+01,-7.3587E-02,
     > -4.9353E+00,-2.4223E+04, 9.2768E+04,-7.6712E+04,-1.3210E+01,
     >  1.2513E+00,-4.5156E+00, 2.4541E+04,-9.5131E+04, 7.8848E+04,
     >  1.9952E-02,-7.1332E-02, 5.5522E-02, 9.8804E-04, 2.3682E-04,
     > -7.9762E-05,-6.3638E-01, 1.9492E+00,-1.4036E+00,-9.9312E-04,
     > -7.8634E-05, 8.2617E-05, 6.8002E-01,-2.1138E+00, 1.5308E+00,
     >  1.3008E-04,-1.0389E+02, 3.5942E+02,-2.7883E+02,-6.0671E-01,
     > -1.3016E-01, 1.4621E-01, 1.2841E+03,-4.3361E+03, 3.3132E+03,
     >  7.0701E-01, 1.2805E-01, 1.3355E-01,-1.4645E+03, 4.9522E+03,
     > -3.7686E+03,-1.0047E-01, 2.7406E+02, 3.5483E+02,-1.3433E+03,
     >  1.0978E+03, 1.9033E+00, 5.3726E-02,-8.1621E-01,-4.3612E+03,
     >  1.6110E+04,-1.2957E+04,-2.2247E+00,-2.1299E-01,-5.8178E-01,
     >  4.9755E+03,-1.8393E+04, 1.4724E+04, 3.1774E-01,-9.2555E+02,
     >  3.4086E+03,-2.7508E+02, 1.1025E+03,-9.3653E+02,-1.4100E+00,
     >  7.3163E-02, 6.6492E-01, 3.3590E+03,-1.3073E+04, 1.0893E+04,
     >  1.6311E+00, 2.4826E-01, 8.3308E-01,-3.7999E+03, 1.4772E+04,
     > -1.2252E+04,-2.3255E-01, 7.0167E+02,-2.7162E+03, 2.2434E+03,
     >  3.0688E+00,-1.0328E+01, 7.8828E+00, 3.6601E-03, 1.3367E-03,
     > -2.9672E-03,-3.2441E+01, 1.0979E+02,-8.3795E+01,-6.6345E-03/


      real rerr2(125)/
     > 3.7074E-02,
     >-5.7300E-02, 1.5212E-02, 4.5952E-04,
     > 1.1568E-04,-2.9315E-04,-4.6018E-01,
     > 9.3624E-01,-4.5908E-01,-6.2914E-05,
     > 1.1699E-03, 2.0141E-03, 6.9968E-02,
     >-1.9348E-01, 1.2176E-01, 5.4214E-07,
     > 1.3845E-04, 2.5311E-03,-2.5396E-03,
     >-1.2757E-04, 2.4971E-04,-1.2737E-04,
     > 7.2023E-03,-4.1404E-03, 4.6704E-04,
     > -4.6388E-03,-5.2545E-03, 4.0159E+01,-1.3481E+02, 1.0186E+02,
     >  1.1796E-03,-9.1088E+00, 3.0200E+01,-2.2552E+01, 4.3562E-01,
     > -1.0404E+01, 3.8414E+01,-3.0978E+01,-1.4730E-02, 4.6327E-03,
     >  1.9716E-02, 1.1236E+02,-4.1952E+02, 3.3862E+02, 2.4150E-02,
     >  1.1098E-02, 2.0122E-02,-1.3812E+02, 5.1058E+02,-4.0773E+02,
     > -4.1791E-03, 3.0702E+01,-1.1132E+02, 8.7622E+01,-1.4199E+00,
     >  5.0230E+00, 8.0171E+00,-3.1384E+01, 2.6350E+01, 1.3147E-02,
     > -6.1508E-03,-1.6808E-02,-8.7538E+01, 3.4530E+02,-2.8922E+02,
     > -1.9581E-02,-1.0895E-02,-2.4705E-02, 1.0611E+02,-4.1369E+02,
     >  3.4296E+02, 3.2847E-03,-2.3191E+01, 8.8502E+01,-7.2288E+01,
     >  1.0469E+00,-3.8873E+00, 3.1142E+00,
     > 1.1348E+00,-1.7657E+00, 4.7686E-01,
     > 1.6653E-02, 4.3488E-04,-7.5168E-03,
     >-1.6696E+01, 3.4692E+01,-1.7470E+01,
     >-4.9697E-03, 4.4232E-02, 5.7617E-02,
     > 5.7800E+00,-1.3886E+01, 7.9819E+00,
     > 3.4744E-04,-5.4411E-01, 1.2683E+00,
     >-7.0771E-01, 1.1282E-02,-2.4800E-02,
     > 1.2909E-02, 1.5171E-01,-6.0417E-01,
     > 7.7405E-01,-5.8981E-02,-5.8502E-03,
     > 8.8611E-04, 5.8326E-03, 6.5418E+00,
     >-1.2978E+01, 6.1069E+00, 1.2462E-03,
     >-1.8442E-02,-2.7954E-02,-1.8335E+00,
     > 4.3674E+00,-2.4393E+00,-6.2354E-05,
     > 1.4746E-01,-3.4127E-01, 1.8285E-01,
     >-3.0479E-03, 6.8138E-03,-3.4673E-03,
     >-7.5270E-02, 4.0914E-02/

      REAL ERR1(200)/
     >  3.7797E+02,-1.2732E+03, 4.8470E+03, 9.7589E+02,-3.9592E+03,
     >  3.3447E+03, 1.9629E-01,-4.2402E-01, 1.9757E-01, 3.0613E-02,
     > -4.0257E-01,-2.0922E+00, 3.0126E+00, 3.8385E-03, 7.3553E-02,
     >  1.4084E+00,-8.4718E+00, 7.8586E+00,-1.6484E-03, 2.2185E-02,
     >  7.4896E-02,-1.5627E+03, 5.0106E+03,-3.7125E+03,-1.1701E+01,
     > -6.9186E-01,-1.4263E+00, 1.5792E+04, 5.0288E+03,-1.7793E+04,
     >  1.3974E+04, 3.1643E+01, 5.0040E+00, 9.9958E+00,-4.8540E+04,
     >  1.6247E+05,-3.7498E+03, 1.4066E+04,-1.1461E+04,-2.0806E+01,
     > -5.0428E+00,-9.7813E+00, 3.5056E+04,-1.2382E+05, 9.7850E+04,
     > -2.0038E-01, 5.9769E-01,-4.0397E-01,-1.5776E-02,-3.7509E-03,
     >  5.7496E-04, 7.2218E+00,-2.0335E+01, 1.3722E+01, 1.2562E-02,
     >  1.4708E+00, 1.8510E+00,-4.1856E+00, 1.9572E-02,-1.3469E-01,
     > -3.7791E-02,-1.5215E+01, 1.8843E+01,-9.9384E-01, 5.4133E-04,
     >  5.6775E-01,-2.4158E+00, 1.5245E+01,-1.4180E+01, 5.3668E-03,
     > -3.5419E-02,-1.4360E-01, 7.8707E+00,-5.7677E+01, 5.5406E+01,
     > -7.5727E-04, 1.4127E-01, 5.8964E-01, 1.0277E+03,-3.3407E+03,
     >  2.4943E+03, 6.1372E+00, 2.0731E+00,-1.0628E-01,-1.1445E+04,
     >  3.6033E+04,-2.6376E+04,-6.4849E+00,-1.5437E+00,-3.1093E+00,
     >  1.1966E+04,-3.3062E+03, 1.1473E+04,-8.9323E+03,-1.7658E+01,
     > -3.0298E+00, 2.4862E+00, 3.6140E+04,-1.2237E+05, 9.3797E+04,
     >  1.8377E+01, 2.4649E-01, 9.5713E+00,-3.7362E+04, 1.2613E+05,
     >  2.4733E+03,-8.9836E+03, 7.2301E+03, 1.2133E+01, 1.0120E+00,
     > -2.0972E+00,-2.6581E+04, 9.4364E+04,-7.4804E+04,-1.2397E+01,
     >  5.8276E-01,-9.1893E+00, 2.7145E+04,-9.6250E+04, 7.6086E+04,
     >  2.4070E-02,-7.3772E-02, 5.1165E-02, 1.4597E-03, 3.3977E-04,
     > -2.6275E-05,-7.2542E-01, 2.0676E+00,-1.4052E+00,-1.3577E-03,
     > -1.4477E-04,-8.5451E-05, 7.4811E-01,-2.1217E+00, 1.4288E+00,
     >  1.7439E-04,-1.6022E+02, 5.2231E+02,-3.9172E+02,-4.1771E-01,
     > -2.3133E-01,-1.9119E-02, 1.6931E+03,-5.4146E+03, 4.0099E+03,
     >  6.5228E-01, 4.5766E-01, 6.7254E-01,-2.0266E+03, 6.3551E+03,
     > -4.6404E+03,-9.4689E-02, 4.2768E+02, 5.1531E+02,-1.7829E+03,
     >  1.3890E+03, 1.1798E+00, 3.1335E-01,-2.5902E-01,-5.3955E+03,
     >  1.8502E+04,-1.4311E+04,-1.8045E+00,-9.6753E-01,-2.0260E+00,
     >  6.3626E+03,-2.1445E+04, 1.6387E+04, 2.6350E-01,-1.3456E+03,
     >  4.5055E+03,-3.8598E+02, 1.3911E+03,-1.1170E+03,-7.9328E-01,
     > -7.6840E-02, 2.5967E-01, 4.0005E+03,-1.4347E+04, 1.1455E+04,
     >  1.1795E+00, 6.2629E-01, 1.6961E+00,-4.6485E+03, 1.6399E+04,
     > -1.2954E+04,-1.7187E-01, 9.8638E+02,-3.4363E+03, 2.7002E+03,
     >  6.0266E+00,-1.9528E+01, 1.4686E+01,-1.7956E-02, 3.3364E-03,
     >  1.2080E-03,-5.5018E+01, 1.7933E+02,-1.3517E+02, 7.9955E-03/


      REAL ERR2(53)/
     > -2.1546E-02,-2.3493E-02, 7.4315E+01,-2.3518E+02, 1.7398E+02,
     > -6.4429E-04,-1.9950E+01, 6.3147E+01,-4.6881E+01, 1.2816E+00,
     > -1.9366E+01, 6.5755E+01,-5.0971E+01, 5.7005E-02, 3.3439E-04,
     >  5.5786E-03, 1.7715E+02,-6.1369E+02, 4.7999E+02,-2.9558E-02,
     >  5.5461E-02, 7.1075E-02,-2.3560E+02, 7.9154E+02,-6.0792E+02,
     >  2.7242E-03, 6.3265E+01,-2.0981E+02, 1.6050E+02,-4.0749E+00,
     >  1.3388E+01, 1.4562E+01,-5.1058E+01, 4.0557E+01,-4.3474E-02,
     > -4.4868E-03,-6.3041E-03,-1.3274E+02, 4.7814E+02,-3.8441E+02,
     >  2.5677E-02,-3.8538E-02,-5.8204E-02, 1.7424E+02,-6.0799E+02,
     >  4.8014E+02,-2.6425E-03,-4.6992E+01, 1.6058E+02,-1.2570E+02,
     >  3.0554E+00,-1.0258E+01, 7.9929E+00/

      LOGICAL FIRST/.TRUE./

! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
*       write(*,*)'dd1',mdelta,mpi,epircm,first
!Define error matrix:
        K = 0
        if (goroper) then 
          DO J = 1,25
            DO I = 1,J
              K = K + 1
              if (k.le.200) RERRMAT(I,J) = RERR1(K)
              if (k.le.325.and.k.gt.200) RERRMAT(I,J)=RERR2(K-200)
            ENDDO
          ENDDO
        endif
       if (.not.goroper) then  
          DO J = 1,22
            DO I = 1,J
              K = K + 1
              if (k.le.200) ERRMAT(I,J) = ERR1(K)
              if (k.le.253.and.k.gt.200) ERRMAT(I,J)=ERR2(K-200)
            ENDDO
          ENDDO
        endif
        if (goroper) then
          DO J = 1,25
              DO I = J+1,25
                RERRMAT(I,J) = RERRMAT(J,I)
              ENDDO
          ENDDO
        endif
        if (.not.goroper) then
          DO J = 1,22
              DO I = J+1,22
                ERRMAT(I,J) = ERRMAT(J,I)
              ENDDO
          ENDDO
        endif
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN
*      write(*,*)'sigdel=',sigdel,mdelta

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(23) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
      ENDDO
      if (goroper) FIT(24)  = FIT(23)/sqrt(Q2)
      if (goroper) FIT(25)  = FIT(23)/q2
      DO I = 1,4
        FIT(I + 12)  = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 16)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 19)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,25
          RSIGERR(J) = 0.0
          DO I = 1,25
            RSIGERR(J) = RSIGERR(J) + FIT(J)*FIT(I)*RERRMAT(I,J)
            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*RERRMAT(I,J)
          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
             SIG_RES2(2) = SIG_RES2(2) + RSIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
             SIG_RES1(2) = SIG_RES1(2) + RSIGERR(J)
          elseIF(j.ge.23.and.j.le.25) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
            SIGroper(2) = SIGroper(2) + RSIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
            SIG_NRES(2) = SIG_NRES(2) + RSIGERR(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
          SIGERR(J) = 0.0
          DO I = 1,22
            SIGERR(J) = SIGERR(J) + FIT(J)*FIT(I)*ERRMAT(I,J)
            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*ERRMAT(I,J)
          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
             SIG_RES2(2) = SIG_RES2(2) + SIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
            SIG_RES1(2) = SIG_RES1(2) + SIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
            SIG_NRES(2) = SIG_NRES(2) + SIGERR(J)
          ENDIF
        ENDDO
      endif

! ERRCHECK should agree with SIGTOTER.
C      ERRCHECK = SQRT(ABS(SIG_RES2(2) + SIG_RES1(2) + SIG_NRES(2)))
C      SIGTOTER = SQRT(SIGTOTER)
      SIG_RES2(2) = SQRT(ABS(SIG_RES2(2)))
      SIG_RES1(2) = SQRT(ABS(SIG_RES1(2)))
      SIG_NRES(2) = SQRT(ABS(SIG_NRES(2)))

      RETURN
      END

      SUBROUTINE H2MOD_FIT_old(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2)
********************************************************************************
* This is a 22 parameter model fit to the NE11 and E133 data and three
* QSQ points of generated BRASSE data to constrain fits at low QSQ. It has
* three background terms and three resonances. Each of these is multiplied
* by a three-term polynomial in Q2. The DELTA resonance has a four term
* polynomial.
*
* 8/91, LMS.
********************************************************************************
      IMPLICIT NONE

      INTEGER I, J
      REAL W2, Q2, SIG_NRES, SIG_RES1, SIG_RES2,  KLR,
     >     KRCM/0./,EPIRCM/0./, PPIRCM, GG, GPI, DEN, SIGDEL,
     >     W, KCM, K, EPICM, PPICM, WDIF

      REAL PI/3.14159265/, ALPHA/0.00729735/, AM/.93828/,
     >     MPPI/1.07783/, MDELTA/1.2355/, MPI/0.13957/,
     >     GAM1WID/0.0676/, GAM2WID/0.1025/, MASS1/1.5062/,
     >     MASS2/1.6800/, DELWID/0.1260/, FIT(30), SQRTWDIF,
     >     XR/0.2060/, MQDEP/2.07/

      REAL COEF(22)/
     >  -0.3982E+04,   0.1764E+05,  -0.1590E+05,   0.1391E+02,
     >  -0.1233E+02,  -0.8206E+01,   0.9006E+04,  -0.3581E+05,
     >   0.3094E+05,  -0.2668E+01,   0.3306E+02,   0.2652E+02,
     >  -0.3134E+04,   0.2059E+05,  -0.1817E+05,   0.1327E+00,
     >   0.6619E+03,  -0.3914E+04,   0.3579E+04,  -0.4032E+02,
     >   0.2233E+03,  -0.2066E+03/

      LOGICAL FIRST/.TRUE./

! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
      ENDIF
      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2)
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model.
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 +
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
      ENDDO
      DO I = 1,4
        FIT(I + 12)  = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 16)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 19)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO

! Find sig_t (in microbarns/gd**2).
      SIG_NRES = 0.0
      SIG_RES1 = 0.0
      SIG_RES2 = 0.0
      DO J = 1,22
        IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
          SIG_RES2 = SIG_RES2 + FIT(J)*COEF(J)
        ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
          SIG_RES1 = SIG_RES1 + FIT(J)*COEF(J)
        ELSE
          SIG_NRES = SIG_NRES + FIT(J)*COEF(J)
        ENDIF
      ENDDO

      RETURN
      END

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2glob.o -c f2glob.f"; -*-
!File E.13. F1990.FORTRN.                                               
!Reference:  L.W.Whitlow, SLAC-Report-357,                              
!            Ph.D. Thesis, Stanford University,                         
!            March 1990.                                                
!For details see file HELP.DOCUMENT.                                    
                                                                        
!Program contains 145 lines of Fortran code, of 72 characters each, with
!no subroutines.  Program requires File E.14 as input.                  
                                                                        
                                                                        
      SUBROUTINE F2GLOB(X,Q2,Target,MODEL,F2,ST,SY,slope,dslope,GOODFIT)
                                                                        
! Returns F2 and related quantities from the either the LAMBDA12 model  
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.           
!                                                                       
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2  
!                                                                       
! Further, program returns uncertainty in F2 based on both statistics   
! (ST) and systematic effects due to our choice of models (SY).         
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical
! uncertainty in this slope.                                            
!                                                                       
! Best model is LAMBDA12.  SY is estimated to be the difference between 
! the two models.                                                       
!                                                                       
! Errors due to overall normalization are not included in program output
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.      
! Systematic errors due to radiative corrections are shown in Reference 
! to be very kinematic independent, and are everywhere <.5%, and thus,  
! are ignored by this routine (also see documentation to dFRC in file   
! HELP.DOCUMENT).                                                       
!                                                                       
! Coefficients and correlation matrix elements are from File            
! E.13 F1990.MATRICES, dated 27Jan90.                                   
                                                                        
      IMPLICIT NONE                                                     
      LOGICAL GOODFIT,FIRST/.TRUE./                                     
      REAL   X, Q2, F2,ST,SY,SLOPE,DSLOPE                               
      REAL   XPP,XP,Y,POLY,F2B,POL1,STB,Q,QTH,DQ,F2TH,QUAD,SCB,F2L      
      REAL   BINDING,STL                                                
      INTEGER MODEL, I, J, K                                            
      CHARACTER*1 TARGET                                                
      REAL*8    B(2,9),L(2,9,9)                      ! OMEGA9 variables 
      REAL*8    C(2,12),M(2,12,12)                   ! LAMBDA12 variable
      REAL*8 V(9),Z(12),U(12),LIN                                       
                                                                        
                                                                        
! Model #9  27 Jan 90.                                                  
      DATA (B(1,J),J=1,9)/                        !  HYDROGEN Ci        
     >   0.7338659870D0,  11.0245522588D0,   2.6185804129D0,            
     >   4.0956321483D0,   0.1206495422D0,   1.9714128709D0,            
     >   3.8893348719D0, -14.0507358314D0,   8.8080576075D0/            
      DATA ((L(1,J,K),K=1,9),J=1,9)/              !  HYDROGEN MijDiDj   
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,             
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,             
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,             
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,             
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,             
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,             
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,             
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,             
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,             
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,             
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,             
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,             
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,             
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,             
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,             
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,             
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,             
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,             
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,             
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,             
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,             
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,             
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,             
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,             
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,             
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,             
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/             
                                                                        
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci       
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,             
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,             
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/             
                                                                        
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj  
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,             
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,             
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,             
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,             
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,             
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,             
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,             
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,             
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,             
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,             
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,             
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,             
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,             
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,             
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,             
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,             
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,             
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,             
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,             
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,             
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,             
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,             
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,             
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,             
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,             
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,             
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/             
                                                                        
!    MODEL #12:    27Jan90.                                             
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci            
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,             
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,             
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,             
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/             
      DATA ((M(1,J,K),K=1,12),J=1,12)/     !     HYDROGEN MijDiDj       
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,             
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,             
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,             
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,             
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,             
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,             
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,             
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,             
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,             
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,             
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,             
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,             
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,             
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,             
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,             
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,             
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,             
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,             
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,             
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,             
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,             
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,             
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,             
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0,             
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,             
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,             
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,             
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,             
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,             
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,             
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,             
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,             
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,             
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,             
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,             
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,             
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,             
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,             
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,             
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,             
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,             
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,             
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,             
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,             
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,             
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,             
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,             
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/             
                                                                        
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci           
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,             
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,             
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,             
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/             
      DATA ((M(2,J,K),K=1,12),J=1,12)/     !     DEUTERIUM MijDiDj      
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,             
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,             
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,             
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,             
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,             
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,             
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,             
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,             
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,             
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,             
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,             
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,             
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,             
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,             
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,             
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,             
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,             
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,             
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,             
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,             
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,             
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,             
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,             
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0,             
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,             
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,             
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,             
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,             
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,             
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,             
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,             
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,             
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,             
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,             
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,             
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,             
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,             
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,             
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,             
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,             
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,             
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,             
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,             
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,             
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,             
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,             
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,             
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/               
      !---------------------------------------------------------------- 
      !---------------------------------------------------------------- 
                                                                        
                                                                        
                                                                        
      i = 1                                                             
      IF (TARGET.EQ.'D') i = 2                                          
      BINDING = 1./(1.-EXP(-MIN(20.,7.7*(1./X+.93828**2/Q2-1.))))       
      IF (i.EQ.1) BINDING = 1.                                          
                                                                        
      !OMEGA9 MODEL FIRST:                                              
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                             
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                             
           Y    = 1.-XP                                                 
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+      
     >            B(i,9)*Y**7                                           
           F2B  = X/XPP*BINDING*POLY                                    
          !-----------------------------------------------------------  
          !V(k) is the derivative of F_2 with respect to parameter k.   
           V(1) = -F2B/XPP/(Q2/X+B(i,2))                                
           V(2) =  F2B/XPP/(Q2/X+B(i,2))**2*(Q2+B(i,1))                 
           POL1 =  3.*B(i,5)*Y**2+4.*B(i,6)*Y**3+5.*B(i,7)*Y**4+        
     >             6.*B(i,8)*Y**5+7.*B(i,9)*Y**6                        
           V(3) = -F2B*POL1/POLY/(Q2/X+B(i,4))                          
           V(4) =  F2B*POL1/POLY/(Q2/X+B(i,4))**2*(Q2+B(i,3))           
           DO 10 j = 5,9                                                
10         V(j) =  F2B/POLY*Y**(j-2)                                    
           STB = 0.                                                     
           DO 11 j = 1,9                                                
           DO 11 k = 1,9                                                
11         STB = STB + L(i,j,k)*V(j)*V(k)                               
           STB = SQRT(STB)*BINDING                                      
                                                                        
      !LAMBDA12 MODEL NEXT:                                             
           Y    = 1.-X                                                  
           q    = LOG(Q2)                                               
           qth  = .2+3.2*X                                              
           dq   = q-qth                                                 
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                  
     >            C(i,4)*Y**6+C(i,5)*Y**7                               
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                   
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq       
           IF (q.GT.qth) QUAD = 0.                                      
           SCB  = (1.+LIN+QUAD)                                         
           F2L  = F2th*SCB*BINDING                                      
          !-----------------------------------------------------------  
          !Z(k) is the derivative of F_2 with respect to parameter k.   
           DO 20 j = 1,5                                                
20         Z(j) = SCB*Y**(j+2)                                          
           Z(6) = 0.                                                    
           IF (q.LT.qth) Z(6) = F2th*dq**2                              
           Z(7) = Z(6)*X                                                
           Z(8) = Z(6)*X**2                                             
           DO 21 j = 9,12                                               
21         Z(j) = F2th*X**(j-9)*dq                                      
           STL = 0.                                                     
           DO 22 j = 1,12                                               
           DO 22 k = 1,12                                               
22         STL = STL + M(i,j,k)*Z(j)*Z(k)                               
           STL = SQRT(STL)*BINDING                                      
                                                                        
          !U(k) is the derivative of slope with respect to parameter k. 
           SLOPE= F2th*LIN/dq*BINDING                                   
           DO 30 j = 1,5                                                
30         U(j) = LIN/dq*Y**(j+2)                                       
           DO 31 j = 6,8                                                
31         U(j) = 0.                                                    
           DO 32 j = 9,12                                               
32         U(j) = Z(j)/dq                                               
           DSLOPE = 0.                                                  
           DO 33 j = 1,12                                               
           DO 33 k = 1,12                                               
33         DSLOPE = DSLOPE + M(i,j,k)*U(j)*U(k)                         
           DSLOPE = SQRT(DSLOPE)                                        
      !---------------------------------------------------------------- 
                                                                        
      F2 = 0.                                                           
      ST = 0.                                                           
      IF (MODEL.EQ. 9) THEN                                             
           F2 = F2B                                                     
           ST = STB                                                     
      ELSEIF (MODEL.EQ.12) THEN                                         
           F2 = F2L                                                     
           ST = STL                                                     
      ELSE                                                              
           WRITE(*,'('' F1990: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')    
      ENDIF                                                             
      SY = ABS(F2B-F2L)                                                 
                                                                        
      GOODFIT = .TRUE.                                                  
      !The following cuts define the region of applicability of F1990.  
      !In order they are:                                               
      !     [radiative corrections convergence criteria in Q2] .and.    
      !     [radiative corrections convergence criteria in x]  .and.    
      !     [stay out of resonance region, W2.ge.3.0]          .and.    
      !     [limitation imposed by maximum beam energy].                
      IF ((Q2.LT..566).OR.(X.LT..062).OR.(X.LT.Q2/(2.*.93828*21.))      
     >   .OR.(X.GT.1./((3.-.93828**2)/Q2+1.)))     THEN                 
                                                                        
C         WRITE(*,'('' WARNING[F1990]: OUTSIDE RECOMMENDED RANGE.'')')  
          GOODFIT=.FALSE.                                               
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               


!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
       Subroutine F2NMC_new(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in CERN_PPE/95-138  Sept 4, 1995  Proton and Deuteron F2 Structure Functions 
! in deep inelastic muon scattering   P. Arneodo et. al.               
!  Published in Phys.Lett.B364:107-115,1995 
!   e-Print Archive: hep-ph/9509406 
!                   
      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi ,FNP_NMC
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2NMC_NEW_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       F2 = F2 * FNP_NMC(X4,QSQ4)
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END

      SUBROUTINE F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)
     >  /-0.02778,  2.926,   1.0362, -1.840, 8.123, -13.074, 6.215, 
     >   -0.04858,  2.863,    .8367, -2.532, 9.145, -12.504, 5.473/ 
      real*8 b(4,2)/ 0.285,   -2.694,   0.0188,  0.0274,  
     >              -0.008,   -2.227,   0.0551,  0.0570/
      real*8 c(4,2)/-1.413,    9.366, -37.79,   47.10, 
     >              -1.509,    8.553, -31.20,   39.98/

!lower limits
      real*8 al(7,2)
     >  /-0.01705,  2.851,   0.8213, -1.156, 6.836, -11.681, 5.645, 
     >   -0.02732,  2.676,    .3966, -0.608, 4.946,  -7.994, 3.686/ 
      real*8 bl(4,2)/ 0.325,   -2.767,   0.0148,  0.0226,  
     >                0.141,   -2.464,   0.0299,  0.0396/
      real*8 cl(4,2)/-1.542,   10.549, -40.81,   49.12, 
     >               -2.128,   14.378, -47.76,   53.63/

!upper limits
      real*8 au(7,2)
     >  /-0.05711,  2.887,   0.9980, -1.758, 7.890, -12.696, 5.992, 
     >    -0.04715,  2.814,    .7286, -2.151, 8.662, -12.258, 5.452/ 
      real*8 bu(4,2)/ 0.247,   -2.611,   0.0243,  0.0307,  
     >               -0.048,  -2.114,   0.0672,  0.0677/
      real*8 cu(4,2)/-1.348,    8.548, -35.01,   44.43, 
     >               -1.517,    9.515, -34.94,   44.42/
      

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               

!---------------------------------------------------------------------

       SUBROUTINE FFHE3(Q2,GE,GM)
***************************************************************************
* This subroutine returns elastic charge and magnetic form factors for
* HE3. Low Q2 parameterization is from McCarthy, et al PRC 15, 1396 (1977).
* High Q2 parameterization for the charge form factor is from Arnold, 
* et al, PRL 40, 1429 (1978). The high Q2 parameterization is for a 
* measured combination of the charge and magnetic form factors. Here, it
* is assumed that the small magnetic form factor can be obtained from the
* low Q2 parameterization evaluated at high Q2.
* Errors are included for future model dependent studies.
*
* 2/94, LMS.
* Copied from ~stuart/e143rc/code/models.f 5/30/96
***************************************************************************
       IMPLICIT NONE
       LOGICAL ERROR/.FALSE./
       REAL*8 Q2, Q2FM, GE, GM, FC, FM, FCERR, FMERR, S1, S2,
     >        S3, S4, S5, S6, Q, TAU, MU/-3.2D0/, M/2.80833D0/, AQ2,
     >        AQ2ERR, FCHIGH, FCHIGHERR, FRAC, Z,
     >        HC2/0.0389379D0/        ! (GeV fm)**2
       REAL*8 AFC/0.675D0/, AFCERR/0.008D0/, BFC/0.366D0/, 
     >        BFCERR/0.025D0/, CFC/0.836D0/, CFCERR/0.032D0/,
     >        DFC/-0.00678D0/, DFCERR/0.00083D0/, PFC/0.9D0/,
     >        PFCERR/0.16D0/, Q0FC/3.98D0/, Q0FCERR/0.09D0/
       REAL*8 AFM/0.654D0/, AFMERR/0.024D0/, BFM/0.456D0/,
     >        BFMERR/0.029D0/, CFM/0.821D0/, CFMERR/0.053D0/
       REAL*8 AA/0.034D0/, AAERR/0.004D0/, BA/2.72D0/, BAERR/0.09D0/


       TAU = Q2/(4.D0*M*M)
       Q2FM = Q2/HC2
       Q = SQRT(Q2FM)
       IF(Q2.LT.0.8D0) THEN 
         FC = ABS(EXP(-AFC*AFC*Q2FM) - 
     >            BFC*BFC*Q2FM*EXP(-CFC*CFC*Q2FM)
     >          + DFC*EXP(-((Q - Q0FC)/PFC)**2))
         IF(ERROR) THEN
           S1 = 2.D0*AFC*Q2FM*AFCERR*EXP(-AFC*AFC*Q2FM)
           S2 = 2.D0*BFC*Q2FM*BFCERR*EXP(-CFC*CFC*Q2FM)
           S3 = 2.D0*CFC*BFC*BFC*Q2FM*Q2FM*CFCERR*EXP(-CFC*CFC*Q2FM)
           S4 = DFCERR*EXP(-((Q - Q0FC)/PFC)**2)
           S5 = 2.D0*DFC*(Q - Q0FC)/(PFC*PFC)*Q0FCERR*
     >              EXP(-((Q - Q0FC)/PFC)**2)
           S6 = S5*(Q - Q0FC)*PFCERR/(PFC*Q0FCERR)
           FCERR = SQRT(S1*S1 + S2*S2 + S3*S3 + S4*S4 + S5*S5 + S6*S6)
         ENDIF
       ENDIF

       FM = ABS(EXP(-AFM*AFM*Q2FM) - BFM*BFM*Q2FM*EXP(-CFM*CFM*Q2FM))
       IF(ERROR) THEN
         S1 = 2.D0*AFM*Q2FM*AFMERR*EXP(-AFM*AFM*Q2FM)
         S2 = 2.D0*BFM*Q2FM*BFMERR*EXP(-CFM*CFM*Q2FM)
         S3 = 2.D0*CFM*BFM*BFM*Q2FM*Q2FM*CFMERR*EXP(-CFM*CFM*Q2FM)
         FMERR = SQRT(S1*S1 + S2*S2 + S3*S3)
       ENDIF

       IF(Q2.GT.0.7D0) THEN
         AQ2 = AA*EXP(-BA*Q2)
         IF(ERROR) THEN
           S1 = AAERR*EXP(-BA*Q2)
           S2 = AQ2*Q2*BAERR
           AQ2ERR = SQRT(S1*S1 + S2*S2)
         ENDIF
         FCHIGH = SQRT(ABS(AQ2*AQ2*(1.D0 + TAU) - FM*FM*MU*MU*TAU))
         IF(ERROR) THEN
           S1 = AQ2*AQ2ERR*(1.D0 + TAU)/FCHIGH
           S2 = FM*FMERR*MU*MU*TAU/FCHIGH
           FCHIGHERR = SQRT(S1*S1 + S2*S2)
         ENDIF
         IF(Q2.GE.0.8D0) THEN
           FC = FCHIGH
           FCERR = FCHIGHERR
         ELSE                      ! Require continuity over overlap region. 
           FRAC = (Q2 - 0.7D0)/0.1D0
           FC = FRAC*FCHIGH + (1.D0 - FRAC)*FC
           IF(ERROR) THEN
             S1 = FRAC*FCHIGHERR
             S2 = (1.D0 - FRAC)*FCERR
             FCERR = SQRT(S1*S1 + S2*S2)
           ENDIF
         ENDIF
       ENDIF

! Absorb Z from Mott cross section here.
       Z = 2.D0
       GE =  Z*ABS(FC)
       GM =  Z*MU*ABS(FM)
       RETURN
       END


!---------------------------------------------------------------------

      FUNCTION f2allm(x,q2)             
!-----------------------------------------------------------------
c       allm97, NMC published measured points Q2>0.75 GeV2
c       for values Q<1 use data of E665!
c       parameterization of F2 , according to
c       H.Abramowicz and A.Levy, hep-ph/9712415
c
c       3*10-6 < x  < 0.85, W2>3GeV2
c       0.   < Q2 < 5000 GeV2, dof=0.97
c
! From Peter Bosted, Nov 00 from Antje Burrel.
!------------------------------------------------------------------
      IMPLICIT NONE
      REAL F2ALLM,X,Q2
      REAL SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      COMMON/ALLM/SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
C  POMERON
      REAL S11,S12,S13,A11,A12,A13,B11,B12,B13,M12
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      REAL S21,S22,S23,A21,A22,A23,B21,B22,B23,M22
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      REAL M02,LAM2,Q02,ALFA,XMP2
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      REAL W2,W,Z
C                                                                               
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
C
      IF(Q2.EQ.0.) THEN                                                        
       S=0.
       Z=1.           
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                                   
       AP=A11                                                            
       BP=B11                                                               
       SP=S11                                                             
       F2P=SP*XP**AP                                                
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))          
       AR=A21                                 
       BR=B21                                 
       SR=S21
       F2R=SR*XR**AR              
C                                                                               
      ELSE                                                                      
       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))                       
       Z=1.-X                                      
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)         
       BP=B11+B12*S**B13                                                  
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)         
       F2P=SP*XP**AP*Z**BP                            
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))                                          
       AR=A21+A22*S**A23                                                     
       BR=B21+B22*S**B23                                                
       SR=S21+S22*S**S23                                                     
       F2R=SR*XR**AR*Z**BR                                                   
    
C                                                                               
      ENDIF

c      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z              
c      SIGal=CIN*(F2P+F2R)                                             
c      f2allm=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)
      F2ALLM = q2/(q2+m02)*(F2P+F2R)
 

      RETURN                                                                    
      END                                                                       

!---------------------------------------------------------------------
      REAL FUNCTION FNP_NMC(X,QSQ)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL X,QSQ,A,B,X2,X3
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      FNP_NMC = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END

!---------------------------------------------------------------------

      REAL FUNCTION FITEMC(X,A,GOODFIT)                                         
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
!    beyond that assume Fermi motion is giving the rise, and we
!    already are taking that into account with the y-smearing of
!    the inelastic
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*4 ALPHA, C,LN_C,X,A ,X_U
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
      fitemc = 1.
      if(A .lt. 2.5) return

c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c       IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------

********************************************************************************
       SUBROUTINE FFD(Q2,GE,GM)
***************************************************************************
* This subroutine returns elastic charge and magnetic form factors for
* deuterium. 
* Errors are included for future model dependent studies.
* ***** Note that this routine is currently under construction and fits
* will all be updated in due time.
* 2/94, LMS.
* Copied from ~stuart/e143rc/code/models.f 5/30/96
***************************************************************************
       IMPLICIT NONE
       REAL*8 MD
       PARAMETER ( MD     = 1.87561339D0)    ! Deuteron mass (GeV).
!       INCLUDE 'RADCON.INC'
       LOGICAL ERROR/.FALSE./
       REAL*8 Q2, GE, GM, AQ, BQ, Q, S1, S2, S3, 
     >        BQERR, AQERR, TAU
       REAL*8 BQA/0.0046D0/, BQAERR/0.0006D0/, BQB/6.8D0/,
     >        BQBERR/0.24D0/, BQC/9.44D-09/, BQCERR/1.28D-08/,
     >        BQD/5.46D0/, BQE/1.6588D0/
       REAL*8 AQA/24.4819D0/, AQAERR/0.1913D0/, AQB/-75.0541D0/, 
     >        AQBERR/0.0425D0/, AQC/162.5866D0/, AQCERR/0.0437D0/,
     >        AQD/3.1238D0/, AQDERR/0.5446D0/, AQE/1.3093D0/,
     >        AQEERR/0.7254D0/

       AQ = 1.D0/(AQA*Q2**4 + AQB*Q2**3 + AQC*Q2**2 +
     >           AQD*Q2 + AQE)**2
       IF(ERROR) THEN
         S1 =  2.D0/(AQA*Q2**4 + AQB*Q2**3 + AQC*Q2**2 +
     >           AQD*Q2 + AQE)**3
         S2 = (AQAERR*Q2**4)**2 + (AQBERR*Q2**3)**2 +
     >        (AQCERR*Q2**2)**2 + (AQDERR*Q2)**2 +
     >         AQEERR**2
         AQERR = SQRT(S1*S1*S2)
       ENDIF

       Q = SQRT(Q2)
       BQ = BQA*EXP(-BQB*Q2) + BQC*EXP(-BQD*(Q - BQE)**2)
       IF(ERROR) THEN
         S1 = BQAERR*EXP(-BQB*Q2)
         S2 = BQA*Q2*BQBERR*EXP(-BQB*Q2)
         S3 = BQCERR*EXP(-BQD*(Q - BQE)**2)
         BQERR = SQRT(S1*S1 + S2*S2 + S3*S3)
       ENDIF
       TAU = Q2/(4.D0*MD*MD)

! Note that A(Q2) = (GC(Q2))**2 + (8/9)*TAU**2*(GQ(Q2))**2 +
! (2/3)*TAU*(1+TAU)(GM(Q2))**2 and 
! B(Q2) = (4/3)*TAU*(1+TAU)**2*(GM(Q2))**2 where
! GC is the charge form factor, GQ is the quadrupole form factor and
! GM is the magnetic form factor. Here, it is assumed that GE and GM
! follow the same functional form as given for elastic nucleon
! scattering.
       GM = SQRT(BQ/(2.D0*TAU))
       GE = SQRT( AQ*(1.D0+ TAU) - TAU*GM*GM)

       RETURN
       END

!---------------------------------------------------------------------

      SUBROUTINE R1998(X,Q2,R,DR,GOODFIT)                               
                                                                        
!----------------------------------------------------------------       
! X      : Bjorken x                                                    
! Q2     : Q squared in (GeV/c)**2                                      
! R      :                                                              
! DR     : Absolute error on R                                          
! GOODFIT:  = .TRUE. if the X,Q2 are within the range of the fit.       
!-------------------------------------------------------------------    
! Model for R, based on a fit to world R measurements. Fit performed by 
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details   
! see Reference.                                                        
!                                                                       
! Three models are used, each model has three free parameters.  The     
! functional forms of the models are phenomenological and somewhat      
! contrived.  Each model fits the data very well, and the average of    
! the fits is returned.  The standard deviation of the fit values is    
! used to estimate the systematic uncertainty due to model dependence.  
!                                                                       
! Statistical uncertainties due to fluctuations in measured values have 
! have been studied extensively.  A parametrization of the statistical  
! uncertainty of R1990 is presented in FUNCTION DR1990.                 
!                                                                       
!                                                                       
! Each model fits very well.  As each model has its own strong points   
! and drawbacks, R1998 returns the average of the models.  The          
! chisquare for each fit (237 points with 6 parameters) are:  
!          ALL DATA  #PTS=237         |     X<  0.07 #PTS= 28
! FIT  #PARAM CHISQ  CHISQ/DF PROB(%) |  CHISQ  CHISQ/DF   PROB(%)
! R1998         217.4   0.94    73.1        28.7   1.06    37.5
! R1998a   6    219.4   0.95    69.8        28.8   1.07    37.1
! R1998b   6    217.7   0.94    72.6        27.8   1.03    42.2
! R1998c   6    221.9   0.96    65.5        30.9   1.15    27.4                

!                                                                       
! This subroutine returns reasonable values for R for all x and for all 
! Q2 greater than or equal to .3 GeV.                                   
!                                                                       
! The uncertainty in R originates in three sources:                     
!                                                                       
!     D1 = uncertainty in R due to statistical fluctuations of the data 
!          as reflected in the error of the fit. 
!          It is parameterized in FUNCTION DR1990, for details see     
!          Reference.                                                   
!                                                                       
!     D2 = uncertainty in R due to possible model dependence, approxi-  
!          mated by the variance between the models.                    
!                                                                       
!     D3 = uncertainty in R due to possible epsilon dependent errors    
!          in the radiative corrections, taken to be +/- .025.  See     
!          theses (mine or Dasu's) for details.  This is copied from R1990                       
!                                                                       
! and the total error is returned by the program:                       
!                                                                       
!     DR = is the total uncertainty in R, DR = sqrt(D1+2D) 7
!          DR is my best estimate of how well we have measured R.  At   
!          high Q2, where R is small, DR is typically larger than R.  If
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,  
!          you might wish to assume DR = 0 at very high Q2.             
!                                                                       
! NOTE:    In many applications, for example the extraction of F2 from  
!          measured cross section, you do not want the full error in R  
!          given by DR.  Rather, you will want to use only the D1 and D2
!          contributions, and the D3 contribution from radiative        
!          corrections propogates complexely into F2.  For more informa-
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or     
!          for explicite detail, see Reference.                         
!                                                                       
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!***** modified 10/28/01 pyb to give R at  max(0.5,Q**2)

      IMPLICIT NONE                                                     
      REAL QP2,FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1998,X,Q2
      REAL Q2_SAVE

!!** Note, in S. Rock version, this is 0.2
      REAL Q2_LIMIT/.5/
      REAL A(6) /4.8520E-02,  5.4704E-01,  2.0621E+00,
     >          -3.8036E-01,  5.0896E-01, -2.8548E-02/
      REAL B(6) /4.8051E-02,  6.1130E-01, -3.5081E-01, 
     >          -4.6076E-01,  7.1697E-01, -3.1726E-02/
      REAL C(6) /5.7654E-02,  4.6441E-01,  1.8288E+00,
     >           1.2371E+01, -4.3104E+01,  4.1741E+01/


      LOGICAL GOODFIT                                                   
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      Q2_SAVE = Q2
      IF(Q2.LT.Q2_LIMIT) THEN   ! added 10/28/01 to match Blok et al
       Q2 = Q2_LIMIT
      ENDIF
                                 
      FAC = 1.+12.*Q2/(Q2+1.)*(.125**2 /(.125**2+X**2))
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!      


!Model A
      QP2 = (A(2)/(Q2**4+A(3)**4)**.25) * (1.+A(4)*X +A(5)*X**2) !1.07   .030
      R_A   = A(1)*RLOG +QP2*X**A(6)                                           
!Model B
      QP2 = (1.+B(4)*X+B(5)*X**2)*(B(2)/Q2 +B(3)/(Q2**2+.09)) !1.06   .042
      R_B   =  B(1)*RLOG+QP2*X**B(6)   
!Model C
      Q2thr =C(4)*(X) +c(5)*x**2 +c(6)*X**3  
      QP2 =  C(2)/SQRT((Q2-Q2thr)**2+C(3)**2) 
      R_C   =  C(1)*RLOG+QP2    

      R     = (R_A+R_B+R_C)/3.                                          

      D1    = DR1998(X,Q2)                                              

      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)               
      D3    = .023*(1.+.5*R)                                            
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                       


      DR    = SQRT(D1**2+D2**2+D3**2)                                   
                                                                        
      GOODFIT = .TRUE.                                                  
      IF ((X.LT.0.02).OR.(Q2.LT.0.3)) GOODFIT = .FALSE.                                   

! Restore Q2
      Q2 = Q2_SAVE

!*** commented out 10/28/01
!      IF(Q2.LT.Q2_LIMIT) THEN   ! added 11/15 to avoid low Q2 problems caused by RLOG
!       R = Q2/Q2_LIMIT * R
!      ENDIF
      

      RETURN                                                            
      END                                                               
                                                                        
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                                        
      REAL FUNCTION DR1998(X,Q2)                                             
                                                                        
! Parameterizes the uncertainty in R1990 due to the statistical         
! fluctuations in the data.  Values reflect an average of the R-values  
! about a neighborhood of the specific (x,Q2) value.  That neighborhood 
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see       
! Reference.                                                            
!                                                                       
! This subroutine is accurate over all (x,Q2), not only the SLAC deep   
! inelastic range.  Where there is no data, for example in the resonance
! region, it returns a realistic uncertainty, extrapolated from the deep
! inelastic region (suitably enlarged).  We similarly estimate the      
! uncertainty at very large Q2 by extrapolating from the highest Q2     
! measurments.  For extremely large Q2, R is expected to fall to zero,  
! so the uncertainty in R should not continue to grow.  For this reason 
! DR1990 uses the value at 64 GeV for all larger Q2.                    
!                                                                       
! XHIGH accounts for the rapidly diminishing statistical accuracy for   
! x>.8, and does not contribute for smaller x.                          
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      REAL Q2,X                   

      DR1998 = .0078 -.013*X +(.070 -.39*X+.70*X**2)/(1.7+Q2)
      RETURN                                                            
      END                                                               

!------------------------------------------------------------------------

c       SUBROUTINE INEFT(QQ,W,W1,W2,amuM)                      
       SUBROUTINE INEFT(QQ,W,W1,W2)
C Modified 6feb87 by lww to accept target information passed through    
C common block /targt/.                                                 
                                                                        
C This program takes the old slac structure function model (Atwood,     
C Bodek, et.al.) and outputs values for W1 and W2 at given kinematics. 
! As of 11/3/95 this version is per NEUCLEON   ! Steve Rock

! amuM is atomic number, ie. 1. 2.xxx etc.
! 11/05 changed to FITEMC_N (with neutron excess) instead of FITEMC
! 8/06 changed back to FITEMC and put in neutron excess from
! diff of d and p fits pyb peb
c this doesnt work: put in CONSTNAT n/p ratio for test
      Implicit None 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      integer iz,ia
      REAL avgN, avgA, avgM, amuM               
      REAL*8 QQ,W,W1,W2,WW,V,VV,OMEGAP,SP,UNIV,BRES,SLACF2,B,fneut
      REAL*8 VW2,X,EMCFAC,UNIVD,BRESD,UNIVP,BRESP,vW2p,vW2d,vW2n
      REAL*8    C(24),CF(11),CD(24),CFD(11)                          
      REAL*8    EF(7) 
      REAL FITEMC
c      REAL FITEMC_N
                                                  
      REAL*8         PM / .93828/,PMPM/.8803/,TPM/1.876512/            
      REAL*8         R /  .18/,ALPHAX/137.0388/,THCONST/0.0174533/
      LOGICAL GOODFIT
      common/testing/prttst,usegrd
      logical prttst,usegrd
      DATA         EF / -0.00136693,-.00510425,-.0375986,-.0946004,     
     +                  -.122435,-.0112751,0.406435/                    
                                                                        
C FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE TERMS)    
                                                                        
      DATA   C(1) / 0.10741163E 01/,  C(2) / 0.75531124E 00/,           
     *       C(3) / 0.33506491E 01/,  C(4) / 0.17447015E 01/,           
     *       C(5) / 0.35102405E 01/,  C(6) / 0.10400040E 01/,           
     *       C(7) / 0.12299128E 01/,  C(8) / 0.10625394E 00/,           
     *       C(9) / 0.48132786E 00/,  C(10)/ 0.15101467E 01/,           
     *       C(11)/ 0.81661975E-01/,  C(12)/ 0.65587179E 00/,           
     *       C(13)/ 0.17176216E 01/,  C(14)/ 0.12551987E 00/,           
     *       C(15)/ 0.74733793E 00/,  C(16)/ 0.19538129E 01/,           
     *       C(17)/ 0.19891522E 00/,  C(18)/-0.17498537E 00/,           
     *       C(19)/ 0.96701919E-02/,  C(20)/-0.35256748E-01/,           
     *       C(21)/ 0.35185207E 01/,  C(22)/-0.59993696E 00/,           
     *       C(23)/ 0.47615828E 01/,  C(24)/ 0.41167589E 00/            
                                                                        
C CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION) OMEGAW FIT   
                                                                        
      DATA CF(1) / 0.25615498E 00/,  CF(2) / 0.21784826E 01/,           
     *     CF(3) / 0.89783738E 00/,  CF(4) /-0.67162450E 01/,           
     *     CF(5) / 0.37557472E 01/,  CF(6) / 0.16421119E 01/,           
     *     CF(7) / 0.37635747E 00/,  CF(8) / 0.93825625E 00/,           
     *     CF(9) / 0.10000000E 01/,  CF(10)/ 0.0           /,           
     *     CF(11)/ 0.50000000E 00/                                      
                                                                        
C FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONANT TERMS)   
                                                                        
      DATA  CD(1) / 0.10521935E 01/, CD(2) / 0.76111537E 00/,           
     *      CD(3) / 0.41469897E 01/, CD(4) / 0.14218146E 01/,           
     *      CD(5) / 0.37119053E 01/, CD(6) / 0.74847487E 00/,           
     *      CD(7) / 0.12399742E 01/, CD(8) / 0.12114898E 00/,           
     *      CD(9) / 0.11497852E-01/, CD(10)/ 0.14772317E 01/,           
     *      CD(11)/ 0.69579815E-02/, CD(12)/ 0.12662466E 00/,           
     *      CD(13)/ 0.15233427E 01/, CD(14)/ 0.84094736E-01/,           
     *      CD(15)/ 0.74733793E 00/, CD(16)/ 0.19538129E 01/,           
     *      CD(17)/ 0.19891522E 00/, CD(18)/-0.24480414E 00/,           
     *      CD(19)/ 0.14502846E-01/, CD(20)/-0.35256748E-01/,           
     *      CD(21)/ 0.35185207E 01/, CD(22)/-0.21261862E 00/,           
     *      CD(23)/ 0.69690531E 01/, CD(24)/ 0.40314293E 00/            
                                                                        
C CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)        
C OMEGAW FIT                                                            
                                                                        
      DATA CFD(1) / 0.47708776E 00/, CFD(2) / 0.21601918E 01/,          
     *     CFD(3) / 0.36273894E 01/, CFD(4) /-0.10470367E 02/,          
     *     CFD(5) / 0.49271691E 01/, CFD(6) / 0.15120763E 01/,          
     *     CFD(7) / 0.35114723E 00/, CFD(8) / 0.93825625E 00/,          
     *     CFD(9) / 0.10000000E 01/, CFD(10)/ 0.0           /,          
     *     CFD(11)/ 0.50000000E 00/                                     
                                                                        
C COMPUTE SOME KINEMATIC QUANTITIES                                     
                                                                        
      WW     = W**2                                                     
      V      = (WW+QQ-PMPM)/2.D0/PM                                     
      VV     = V*V                                                      
      OMEGAP = TPM*V/QQ+PMPM/QQ                                         
                                                                        
C OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION                      
      OMEGAP = DMIN1(20.0D0,OMEGAP)                                       
                                                                        
      SP = 1.0-EXP(-7.7*(OMEGAP-1.0))                                   
C pyb *** modified to get proton too when IA ne IZ*2
      IF (IA. ne. 2*IZ) THEN 
C          UNIVERSAL AND RESONANCE FIT FOR HYDROGEN                     
           UNIVP = SLACF2(W,QQ,CF)                                      
           BRESP = B(W,QQ,C)                                            
      endif
      if(amuM.ge.1.5) then
C          UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM                    
           UNIVd = SLACF2(W,QQ,CFD)/SP
           BRESd = B(W,QQ,CD)          
      ENDIF                                                             
                                                                        
C COMPUTE VW2,W2,W1                                                     
                                                                        
      if(amuM.le.1.5) then
        vW2 = UNIVp * BRESp
      else
        VW2p    = UNIVp * BRESp 
        VW2d    = UNIVd * BRESd / 2. ! per nucleon
        vW2n    = 2.*vW2d - vW2p
        fneut = float(IA - 2*IZ) / float(IA)
c       if(fneut.lt.0.) write(6,'(''ERROR, fneut='',f8.3,2i3)') 
c    >    fneut,ia,iz
c *** set fneut back to zero: assume n/p=1!
        fneut=0.
        vW2 = (1. - fneut) * vW2d + fneut * vW2n
c        if(prttst) write(98,'(8f8.3)')  qq,w,vw2p,vw2d,vw2n,
c     >    vw2,vw2n/vw2d,fneut
      endif

      W2     = VW2/V                                                    
      W1     = (1.0D0+VV/QQ)/(V*(1.0D0+R))*VW2                          
!      if(prttst) write(*,'(1x,''univ...='',6f10.4)') sp,univ,bres,
!     >  vw2,w2,w1

      IF (amuM.LE.2.5) RETURN                                               
      X      = QQ/2./PM/V
c Modified 11/05 to include neutron excess pyb peb
c Modified 8/06 to do neutron excess as above, not in EMC
      EMCFAC= FITEMC(REAL(X),REAL(amuM),GOODFIT)
c      EMCFAC= FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GOODFIT)
                                                                        
      W2     = W2*EMCFAC                                                
      W1     = W1*EMCFAC                                                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION SLACF2(WM,QSQ,CF)                                        
                                                                        
C UNIVERSAL FUNCTION FOR ATWOOD'S FIT                                   

      Implicit none                                      
      REAL*8    WM,QSQ,CF(11)                                               
      REAL*8    PM2/1.876512/, PMSQ/.8803/, PHTR/.61993/
      REAL*8    V,OMEGA,XX,XPX,OMEGAW,ARG

                                                                        
C OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING                             
                                                                        
      V      = (WM**2+QSQ-PMSQ)/PM2                                     
      OMEGA  = 2.*CF(8)*V/QSQ                                           
      XX     = 1./OMEGA                                                 
      XPX    = CF(9)+CF(10)*(XX-CF(11))**2                              
      OMEGAW = (2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))                         
      ARG    = 1.-1./OMEGAW                                             
                                                                        
      SLACF2 = OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+                    
     >         CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)                  
      SLACF2 = SLACF2*XPX                                               
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION B(WM,QSQ,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(24),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG,BRES
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,ressv(4)
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.8803/, PM2/1.876512/, PM/.93828/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      common/testing/prttst,usegrd
      logical prttst,usegrd
      DATA      LSPIN/1,2,3,2/                                       
                                                                        
C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(22)+C(23)*(X-C(24))**2                                 
      PIEMSQ = (C(1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(1)) B1 = C(2)                                         
      EB1 = C(3)*(WM-C(1))                                              
      IF (EB1.LE.25.) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(4)) B2 = (1.-C(2))                                    
      EB2 = C(5)*(WSQ-C(4)**2)                                          
      IF (EB2.LE.25.0) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(INDEX)                                            
           IF (I.EQ.1) RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2              
           RMA    = C(INDEX+1)                                          
           IF (I.EQ.3) RMA=RMA*(1.D0+C(20)/(1.D0+C(21)*QSQ))            
           RWD    = C(INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           ressv(i)=res
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
c      if(prttst) write(*,'(1x,''w,q2,res='',6f7.3)') wm,qsq,
c     >  ressv
                                                                   
C FORM VW2/F2                                                           
                                                                        
      B = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
!      if(prttst) write(*,'(1x,''b...'',6f10.5)') b,bbkg,xpx,ressum                                                                  
      RETURN                                                            
      END                                                               
                  
!------------------------------------------------------------------------


                                                                        
      SUBROUTINE NFORM(IG,QQG,GEP,GEN,GMP,GMN)                             

!-----------------------------------------------------------------------
C CALCULATE NUCLEON FORM FACTORS                                        
! Modified by Steve Rock on 6/21/96 adding Peters's Phys Rev C. fit
!  and putting IG in arguements
!
C IG =1 - DIPOLE WITH FORM FACTOR SCALING AND GEN=0.0                
C     2 - IJL FIVE PARAMETER MODEL FIT                               
C     3 - GEP AND GMP FROM IJL, GMN=DIPOLE, GEN=GALSTER              
C            WE CALL THIS THE "BEST FIT" NUCLEON FORM FACTORS           
C     4 - BEST FIT EXCEPT GEN = 0.0                                  
C     5 - BLATNIK + ZOVKO VDM FIT                                    
C     6 - JANNSENS 1966 STANDARD FIT                                 
C     7 - DIPOLE + F1N = 0.0                                         
C     8 - GES = 0.0,  GMS = 1.0                                      
C     9 - GES = 1.0,  GMS = 0.0                                      
C    10 - HOHLER1 - PROTON AND NEUTRON FROM FIT 8.2                  
C    11 - HOHLER2 - PROTON FROM FIT 5.3, NEUTRON FROM 8.2            
C    12 - GARI + KRUMPELMANN, Z PHYS. A322,689(1985)                 
C    13 - KORNER + KURODA, PHYS. REV. D16,2165(1977)  
C    14 - GARI AND KRUMPELMANN WITH NE11 FIT (PETER')                  
C     15 - Peter Bosted's fit from SLAC-PUB-6651 (NE11 data + other) in
c     Phys Rev C
C    16 - Radyushkin,  Acta Physica Polonica B15 403,
c     (1984) 
c    
C QQG = INPUT Q SQUARED (GEV**2)  
!---------------------------------------------------------------------------


      IMPLICIT NONE
      INTEGER IG
      REAL*8 QQG,GEP,GEN,GMP,GMN,QQ,TAU
      REAL*8 GT,T1,T2,ALPH,TOP,BOT,RHO,F1S,F1V,F2S,F2V
      REAL*8 RS,RV,F1E,F2E,F1M,F2M
      REAL*8 F1,F2,F3,GES,GMS,GEV,GMV
      REAL*8 F1RHO,F2RHO,F1P,F2P
      REAL*8 QQP,C1,C2,C3,C4,F2VK,F2SK,F1N,F2N
      REAL*8 Q,Q3,Q4
      REAL*8 FRAC1,FRAC2,GD
      INTEGER I,IN

C IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS)  
C PHYS LETT. 43B, 191(1973)              
      REAL*8   GAM, BR, BW, BF, AF
      REAL*8   RMN2, RMW2, RMF2, RMR2, GAMR, PI   
      REAL*8   RMPI/ .139/,  RMPI2/ .019321  /                         
                                                                        
C PARAMETERS FOR BLATNIK AND ZOVKO VDM FIT                              
C ACTA PHYSICA AUSTRIACA 39, 62(1974)                                   
C VECTOR MESON MASSES SQUARED (GEV UNITS)  
      REAL*8     TRO, TROP, TROPP, TFI, TOM, TOMP     
                                                                        
C FITTED PARAMETERS        
      REAL*8     RMUS, RMUV, BS, BV         
                                                                        
      REAL*8     RMRHO2, CRHO, RKRHO, RKS, RMOMG2, COMG, RKOMG, RKV     
      REAL*8     RLAM12, RLAM22, RLAMQ2
                                                                        
C PARAMETERS FOR KORNER AND KURODA                                      
C VECTOR MESON MASSES SQUARED USING REGGE PARAMETER ALPHA=1.    
      REAL*8     VRH1, VRH2, VRH3, VOM1, VOM2, VOM3          
      REAL*8     RMUP/ 2.792782/ ,RMUN/ -1.913148 /                      
      common/testing/ prttst,usegrd
      logical prttst,usegrd

! Parameters for Radyuskin: From Table on p 414 of paper.               
      REAL R_QSQ(12)   / 1.,   2.,   3.,   4.,   5.,   6.,              
     >                   8.,  10.,  12.,  15.,  20.,  30./              
      REAL R_GMP_U_D(12)/ .91, 1.01, 1.05, 1.05, 1.04, 1.02,            
     >                   0.97, 0.91, 0.86, 0.78, 0.67, 0.53/            
      REAL R_GEP_D(12)  /1.00, 1.13, 1.16, 1.15, 1.11, 1.06,            
     >                   0.95, 0.86, 0.77, 0.67, 0.54, 0.38/            
      REAL R_GMN_U_D(12)/0.82, 0.80, 0.79, 0.77, 0.76, 0.74,            
     >                   0.70, 0.65, 0.61, 0.56, 0.49, 0.38/            
      REAL R_GEN_D(12)  /-.13, -.12, -.10, -.06, -.03, 0.00,            
     >                   0.05, 0.08, 0.11, 0.13, 0.14, 0.14/    
      DATA     RMN2, RMW2, RMF2, RMR2, GAMR, PI                         
     *                / 0.8817, .6146, 1.0384, 0.5852 , .112 , 3.14159/
      DATA     GAM, BR, BW, BF, AF  /0.25, 0.672, 1.102, 0.112, -0.052/
      DATA     TRO, TROP, TROPP, TFI, TOM, TOMP                         
     *                       / 0.585, 1.30,  2.10,  1.039, 0.614, 1.40 /
      DATA     RMUS, RMUV, BS, BV  / -0.060, 1.853, -0.91, -1.10 /      
      DATA     RMRHO2, CRHO, RKRHO, RKS, RMOMG2, COMG, RKOMG, RKV       
     *         /0.6022, 0.377, 6.62, -0.12, 0.6147, 0.411, 0.163, 3.706/
      DATA     RLAM12, RLAM22, RLAMQ2  /  0.632, 5.153, 0.0841 /        
      DATA       VRH1, VRH2, VRH3, VOM1, VOM2, VOM3                       
     *             / 0.593, 1.593, 2.593, 0.614, 1.614, 2.614 /
C=======================================================================
                                                                        
      QQ  = QQG/.197328**2                                              
      TAU = QQG/(4.*RMN2)                                               
      GO TO (110,120,120,120,150,160,170,180,190,200,200,               
     > 220,230,240,250,260) IG                                             
C DIPOLE                                                                
  110 GEP = 1./(1.+QQG/0.71)**2                                         
      GEN = 0.0                                                         
      GMP = RMUP*GEP                                                    
      GMN = RMUN*GEP                                                    
      RETURN                                                            
                                                                        
C IJL 5 PARAMTER JOB 
  120 GT  = 0.5/(1.+GAM*QQG)**2                                         
      T1  = SQRT(QQG+4.*RMPI2)                                          
      T2  = SQRT(QQG)                                                   
      ALPH= 2.*T1*LOG((T1+T2)/(2.*RMPI))/(T2*PI)                       
      TOP = RMR2+8.*GAMR*RMPI/PI                                        
      BOT = RMR2+QQG+(4.*RMPI2+QQG)*GAMR*ALPH/RMPI                      
      RHO = TOP/BOT                                                     
      F1S = GT*((1.-BW-BF)+BW/(1.+QQG/RMW2)+BF/(1.+QQG/RMF2))           
      F1V = GT*((1.-BR)+BR*RHO)                                         
      F2S = GT*((-0.12-AF)/(1.+QQG/RMW2)+AF/(1.+QQG/RMF2))              
      F2V = GT*(3.706*RHO)                                              
      GEP = F1V+F1S-TAU*(F2V+F2S)                                       
      GEN = F1S-F1V-TAU*(F2S-F2V)                                       
      GMP = F1V+F1S+F2V+F2S                                             
      GMN = F1S-F1V+F2S-F2V                                             
      IF (IG.EQ.2) RETURN                                               
      GD  = 1./(1.+QQG/.71)**2                                          
      GMN = RMUN*GD                                                     
      GEN = -RMUN*TAU*GD/(1.+5.6*TAU)                                   
      IF (IG.EQ.3) RETURN                                               
      GEN = 0.0                                                         
      RETURN                                                            
                                                                        
C BLATNIK AND ZOVKO                                                     
  150 RS  = 1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))                
      RV  = 1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP))              
      F1E = (0.5-TAU*(RMUS+2.*RMN2*BS))*RS                              
      F2E = (0.5-TAU*(RMUV+2.*RMN2*BV))*RV                              
      F1M = (0.5+RMUS-0.5*BS*QQG)*RS                                    
      F2M = (0.5+RMUV-0.5*BV*QQG)*RV                                    
      GEP = F1E+F2E                                                     
      GMP = F1M+F2M                                                     
      GEN = F1E-F2E                                                     
      GMN = F1M-F2M                                                     
      RETURN                                                            
                                                                        
C JANNSSENS                                                             
  160 F1  = 1.+QQ/15.7                                                  
      F2  = 1.+QQ/26.7                                                  
      F3  = 1.+QQ/8.19                                                  
      GES = 0.5  *(2.5 /F1-1.6 /F2+0.10)                                
      GMS = 0.44 *(3.33/F1-2.77/F2+0.44)                                
      GEV = 0.5  *(1.16/F3-0.16)                                        
      GMV = 2.353*(1.11/F3-0.11)                                        
      GEP = GES+GEV                                                     
      GMP = GMS+GMV                                                     
      GEN = GES-GEV                                                     
      GMN = GMS-GMV                                                     
      RETURN                                                            
                                                                        
C DIPOLE + F1N = 0.0                                                    
  170 GEP = 1./(1.+QQG/0.71)**2                                         
      GEN = -RMUN*TAU*GEP                                               
      GMP =  RMUP*GEP                                                   
      GMN =  RMUN*GEP                                                   
      RETURN                                                            
                                                                        
  180 GEP = 0.0                                                         
      GEN = 0.0                                                         
      GMP = 1.0                                                         
      GMN = 0.0                                                         
      RETURN                                                            
                                                                        
  190 GEP = 1.0                                                         
      GEN = 0.0                                                         
      GMP = 0.0                                                         
      GMN = 0.0                                                         
      RETURN                                                            
                                                                        
C HOHLER1 AND HOHLER2                                                   
  200 F1RHO = 0.5*(0.955+0.090/(1.+QQG/0.355)**2)/(1.+QQG/0.536)        
      F2RHO = 0.5*(5.335+0.962/(1.+QQG/0.268))   /(1.+QQG/0.603)        
      F1S   =  0.71/(0.6129+QQG)-0.64/(1.0404+QQG)-0.13/(3.240+QQG)     
      F2S   = -0.11/(0.6129+QQG)+0.13/(1.0404+QQG)-0.02/(3.240+QQG)     
      F1V   = F1RHO+0.05/(1.464+QQG)-0.52/(6.0025+QQG)+0.28/(8.7025+QQG)
      F2V   = F2RHO-1.99/(1.464+QQG)+0.20/(6.0025+QQG)+0.19/(8.7025+QQG)
      GEP = F1V+F1S-TAU*(F2V+F2S)                                       
      GEN = F1S-F1V-TAU*(F2S-F2V)                                       
      GMP = F1V+F1S+F2V+F2S                                             
      GMN = F1S-F1V+F2S-F2V                                             
      IF (IG.EQ.10) RETURN                                              
                                                                        
C HOHLER2 - USE PROTON FIT 5.3                                          
      F1P = F1RHO+0.67/(0.6129+QQG)-0.39/(0.9216+QQG)-0.54/( 2.7556+QQG)
      F2P = F2RHO+0.04/(0.6129+QQG)-1.88/(1.2996+QQG)+0.24/(10.1761+QQG)
      GEP = F1P-TAU*F2P                                                 
      GMP = F1P+F2P                                                     
      RETURN                                                            
                                                                        
C GARI AND KRUMPELMANN                                                  
 220  QQP  = QQG*LOG(((RLAM22+QQG)/RLAMQ2))/LOG(RLAM22/RLAMQ2)          
      C1   = RLAM12/(RLAM12+QQP)                                        
      C2   = RLAM22/(RLAM22+QQP)                                        
      F1   = C1*C2                                                      
      F2   = F1*C2                                                      
      C3   = RMRHO2/(RMRHO2+QQG)                                        
      C4   = RMOMG2/(RMOMG2+QQG)                                        
      F1V  = (C3*CRHO+(1-CRHO))*F1                                      
      F1S  = (C4*COMG+(1-COMG))*F1                                      
      F2VK = (C3*CRHO*RKRHO+(RKV-CRHO*RKRHO))*F2                        
      F2SK = (C4*COMG*RKOMG+(RKS-COMG*RKOMG))*F2                        
      F1P  = 0.5*( F1S+F1V)                                             
      F1N  = 0.5*( F1S-F1V)                                             
      F2P  = 0.5*(F2SK+F2VK)                                            
      F2N  = 0.5*(F2SK-F2VK)                                            
      GEP  = F1P-TAU*F2P                                                
      GMP  = F1P+F2P                                                    
      GEN  = F1N-TAU*F2N                                                
      GMN  = F1N+F2N                                                    
      RETURN                                                            
                                                                        
C KORNER AND KURODA                                                     
  230 F1S = (1/(1+QQG/VOM1))*(1/(1+QQG/VOM2))                           
      F1V = (1/(1+QQG/VRH1))*(1/(1+QQG/VRH2))                           
      F2S = F1S*(1/(1+QQG/VOM3))                                        
      F2V = F1V*(1/(1+QQG/VRH3))                                        
      F1P = 0.5*F1S+0.5*F1V                                             
      F1N = 0.5*F1S-0.5*F1V                                             
      F2P = (RMUP-1)*(-0.0335*F2S+1.0335*F2V)                           
      F2N =    -RMUN*(-0.0335*F2S-1.0335*F2V)                           
      GEP = F1P-TAU*F2P                                                 
      GMP = F1P+F2P                                                     
      GEN = F1N-TAU*F2N                                                 
      GMN = F1N+F2N                                                     
      RETURN                                                            
                                                                        
C GARI AND KRUMPELMANN WITH NE11 FIT (PETER')                           
 240  QQP  = QQG*LOG(((RLAM22+QQG)/RLAMQ2))/LOG(RLAM22/RLAMQ2)          
      C1   = RLAM12/(RLAM12+QQP)                                        
      C2   = RLAM22/(RLAM22+QQP)                                        
      F1   = C1*C2                                                      
      F2   = F1*C2                                                      
      C3   = RMRHO2/(RMRHO2+QQG)                                        
      C4   = RMOMG2/(RMOMG2+QQG)                                        
      F1V  = (C3*CRHO+(1-CRHO))*F1                                      
      F1S  = (C4*COMG+(1-COMG))*F1                                      
      F2VK = (C3*CRHO*RKRHO+(RKV-CRHO*RKRHO))*F2                        
      F2SK = (C4*COMG*RKOMG+(RKS-COMG*RKOMG))*F2                        
      F1P  = 0.5*( F1S+F1V)                                             
      F1N  = 0.5*( F1S-F1V)                                             
      F2P  = 0.5*(F2SK+F2VK)                                            
      F2N  = 0.5*(F2SK-F2VK)                                            
      GMP  = F1P+F2P                                                    
      GEP  = GMP/RMUP                                                   
      GEN  = 0.0                                                        
      GMN  = GMP/RMUP * RMUN                                            
      RETURN                                                            

! Peter Bosted's fit from SLAC-PUB-6651 (NE11 data + other) in Phys Rev C
 250  CONTINUE
      Q = SQRT(QQG)
      Q3= Q*QQG
      Q4 = QQG*QQG
      TAU = QQG/3.52
      GEP = 1./  (1.+0.14*Q +3.01*QQG + 0.02*Q3 +1.20*Q4 +.32*Q**5)
      GMP = RMUP*GEP
      GMN = RMUN/(1.-1.74*Q +9.29*QQG - 7.63*Q3 +4.63*Q4)
      GEN = 1.25* RMUN*TAU/(1.+18.3*TAU)/((1.+QQG/.71)**2)
c     if(prttst) write(8,'(1x,11f7.3)') qqg,q,q3,q4,gep,gmp,gmn,gen
c !!! this line was missing up until 8/06!!!!
      return

 260  CONTINUE
! Radyushkin:                                                           
      GD = 1./(1.+ QQG/.71)**2                                          
      IF(QQG.LT.R_QSQ(1)) QQG=1.                                        
      IF(QQG.GT.R_QSQ(12))QQG=12.                                       
      DO I=1,11                                                         
       IF(QQG.GE.R_QSQ(I) .AND. QQG.LE.R_QSQ(I+1) ) THEN                
        IN = I                                                          
        GO TO 241                                                       
       ENDIF                                                            
      ENDDO                                                             
! Out of range.                                                         
      GMP=0.                                                            
      GMN=0.                                                            
      GEP=0.                                                            
      GEN=0.                                                            
      RETURN                                                      
! Do linear interpolation                                               
  241 FRAC1 = (QQG - R_QSQ(IN) )/(R_QSQ(IN+1) -R_QSQ(IN) )              
      FRAC2 = (R_QSQ(IN+1) -QQG)/(R_QSQ(IN+1) -R_QSQ(IN) )              
      GMP = RMUP*GD* (R_GMP_U_D(IN) * FRAC2 + R_GMP_U_D(IN+1) *FRAC1)   
      GMN = RMUN*GD* (R_GMN_U_D(IN) * FRAC2 + R_GMN_U_D(IN+1) *FRAC1)   
      GEP =      GD* (R_GEP_D  (IN) * FRAC2 + R_GEP_D  (IN+1) *FRAC1)   
      GEN =      GD* (R_GEN_D  (IN) * FRAC2 + R_GEN_D  (IN+1) *FRAC1)
      RETURN        
                                                                
      END                                                               
!------------------------------------------------------------------------
C      DOUBLE PRECISION FUNCTION FUNCAL_139(ES1)                             
C                                                                                                                                                
CC---------------------------------------------------------------------  
CC 5/6/87: Allows OMS to be < DEELIM * XI instead of 10*XI as in         
CC          previous version called FUNCA                                
CC         Called by RDIATL                                              
CC 5/29/87: Calculates equivilent radiator at QSQ(Q2S) of interaction    
CC           instead of previous nominal QSQ (calculated elsewhere)      
CC  Steve Rock                                                           
CC 8/20/87: Use Bardin Vacuum term, use single precission logs  -Steve   
CC--------------------------------------------------------------------                                                                           FUN00130
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL BREMS_139                                         
C      COMMON/RAD/ES,EP,SNSQ,R,FAC1,BTA,BTB,ZP,ZN,A,XI,IFLG              
C      COMMON/IFIT/ IFIT                                                 
C      include '/u/ra/ser/radiative/e139code/func.common'
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      REAL*8  DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      DATA ALPHA/.00729735D0/,PI/3.1415926535D0/,EMSQ   /.26113D-6/     
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C                                                                        
CC     CALCULATES THE INTEGRAND FOR ES1 INTEGRATION IN RADIAT            
C                                                                        
CC    CALCULATE KINEMATICS                                               
C      Q2S   =4. *    (ES1)*    (EP)*    (SNSQ)                          
C      OMS=ES-ES1                                                        
C      VS=OMS/ES                                                         
CC                                                                       
CC     QED TERMS AND MULT. PHOTON NORMALIZATION                          
C      VAC_TERM_L = 2.D0 * DELVAC(SNGL(Q2S))                             
C      DLOG_Q = DLOG(Q2S   /EMSQ   )                                     
C      VERTEX_TERM_L = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG_Q )              
C                                                                        
C      FACS=FAC1 + VAC_TERM_L + VERTEX_TERM_L                            
C                                                                        
CC     BREMSSTRAHLUNG SPECTRUM                                           
CC$$   PHIS=1.D0-VS+0.75D0*VS**2                                         
C                                                                        
CC Installed 6/4/87 to get better value of bremstrung spectrum using     
CC  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C      PHIS= BREMS_139(VS,ZP)                                                
C                                                                        
CC     PUT IT TOGETHER                                                   
C      CROSS=0.D0                                                        
C       IF(IFIT.EQ.8) THEN                                               
C        CALL XINELS(ES1,EP,SNSQ,CROSS)                                  
C       ELSE                                                             
C        IF (IFLG.EQ.0) CROSS=QUASIL(ES1,EP,SNSQ,ZP,ZN,A)                
C        IF (IFLG.EQ.1) CALL STRUCT(ES1,EP,SNSQ,CROSS,SIGMOT,110+IFIT)   
C       ENDIF                                                            
C                                                                        
CC     Equivalent raditior at scattering Q**2   5/28/87                  
C      BTEQS = (ALPHA/PI) *(DLOG_Q -1.D0)                                
C                                                                        
C      BTAS=B*DTAFTR+BTEQS                                               
C      BTBS=B*DTBEFR+BTEQS                                               
C                                                                        
CC$$   BTAS =BTA                                                         
CC$$   BTBS =BTB                                                         
C                                                                        
C      FUNCAL_139=(OMS/(R*EP))**BTAS*VS**BTBS*(BTBS/OMS*PHIS+                
C     * XI/(2.D0*DMAX1(OMS,DEELIM*XI)**2))*FACS*CROSS                    
CC                                                                       
C      RETURN                                                            
C      END                                                               
CC---------------------------------------------------------------------  
C 5/6/87: Allows OMP to be < DEELIM * XI instead of 10*XI as in         
C          previous version called FUNCB                                
C         Called by RDIATL                                              
C 5/29/83: Calculates equivilent radiator at QSQ(Q2S) of interaction    
C           instead of previous nominal QSQ (calculated elsewhere)      
C  Steve Rock                                                           
C 8/20/87: Use Bardin Vacuum term, use single precission logs  -Steve   
C--------------------------------------------------------------------   
                                                                        
                                                                        
                                                                        
C      DOUBLE PRECISION FUNCTION FUNCBL_139(EP1)                             
C                                                                        
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL BREMS_139                                                        
C      COMMON/RAD/ES,EP,SNSQ,R,FAC1,BTA,BTB,ZP,ZN,A,XI,IFLG              
C      COMMON/IFIT/ IFIT                                                 
C      include '/u/ra/ser/radiative/e139code/func.common'
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP
C      REAL*8 DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP                                          
C      DATA ALPHA/.00729735D0/,PI/3.1415926535D0/,EMSQ   /.26113D-6/     
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C                                                                        
CC     CALCULATES THE INTEGRAND FOR EP1 INTEGRATION IN RADIAT            
C                                                                        
CC     CALCULATE KINEMATICS                                              
C      Q2P=4.D0 *(ES)*(EP1)*(SNSQ)                                       
C      OMP=EP1-EP                                                        
C      VP=OMP/EP1                                                        
C      DLOG_Q = DLOG(Q2P/EMSQ)                                           
C                                                                        
CC     QED TERMS AND MULT. PHOTON NORMALIZATION                          
C      VAC_TERM_L = 2.D0 * DELVAC(SNGL(Q2P))                             
C      VERTEX_TERM_L = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG_Q )              
C                                                                        
C      FACP=FAC1 + VAC_TERM_L + VERTEX_TERM_L                            
CC                                                                       
CC     BREMSSTRAHLUNG SPECTRUM                                           
CC$$   PHIP=1.D0-VP+0.75D0*VP**2                                         
C                                                                        
CC Installed 6/4/87 to get better value of bremstrung spectrum using     
CC  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C      PHIP= BREMS_139(VP,ZP)                                                
CC                                                                       
CC     PUT IT TOGETHER                                                   
C      CROSS=0.D0                                                        
C       IF(IFIT.EQ.8) THEN                                               
C        CALL XINELS(ES,EP1,SNSQ,CROSS)                                  
C       ELSE                                                             
C        IF (IFLG.EQ.0) CROSS=QUASIL(ES,EP1,SNSQ,ZP,ZN,A)                
C        IF (IFLG.EQ.1) CALL STRUCT(ES,EP1,SNSQ,CROSS,SIGMOT,110+IFIT)   
C       ENDIF                                                            
C                                                                        
C                                                                        
CC     Equivalenr raditior at scattering Q**2   5/28/87                  
C      BTEQP = (ALPHA/PI) *(DLOG_Q -1.D0)                                
C                                                                        
C      BTAP=B*DTAFTR+BTEQP                                               
C      BTBP=B*DTBEFR+BTEQP                                               
C                                                                        
CC$$   BTAP =BTA                                                         
CC$$   BTBP = BTB                                                        
C                                                                        
C      FUNCBL_139=VP**BTAP*(R*OMP/ES)**BTBP  *(BTAP/OMP*PHIP+                
C     * XI/(2.D0*DMAX1(OMP,DEELIM*XI)**2))*FACP*CROSS                    
CC                                                                       
C      RETURN                                                            
C      END                            

!------------------------------------------------------------------------
                                   
C $MEMBER=FUNX, DATE=75061922, USER=KCE                                 
      DOUBLE PRECISION FUNCTION FUNXL_139(COSK) 
      IMPLICIT NONE                            

      REAL*8 COSK
      COMMON/TAL1/ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1 
      REAL*8  ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1 
      COMMON/TAL2/ZP,ZN,AT,IFLG                                         
      REAL*8 ZP,ZN,AT
      INTEGER IFLG
      REAL*8 OM,QSQ,FAC2,VAC_TERM,VERTEX_TERM,FAC,FW1,FW2,A,AP,BP,GNU
      REAL*8 X,Y,TERM1,TERM21,TERM22,TERM23,TERM24,TERM25,TERM26,TERM2
      REAL*8 EMSQ/.26113D-6/,PI/3.1415926535D0/                           
      REAL*8 ALPHA/7.2973515D-3/
      REAL*8 TWO_ALPH_PI/.00464564D0/  
      REAL*4 QSQ_4,W1,W2,FF,ESOM_4,delvac
      common/testing/ prttst,usegrd
      logical prttst,usegrd

                                                                        
C -----------------------------------------------------------------     
C 8/28/87: Single preciession log - Steve Rock                          
C----------------------------------------------------------------- 
      OM=.5D0*(USQ-TMSQ)/(U0-UVEC*COSK)                                 
      QSQ=2.D0*EMSQ-2.D0*SP-2.D0*OM*(ES-EP-UVEC*COSK)  
      QSQ_4 = -QSQ
      ESOM_4 = ES-OM
      IF(IFLG.EQ.0)CALL NUC_FORM_FACTOR(QSQ_4,W1,W2,FF) 
      IF(IFLG.EQ.1) CALL QE_PEAK(QSQ_4,ESOM_4,W1,W2)   
C$$      CALL ASTRUL_139(-QSQ,ZP,ZN,AT,W1_8,W2_8,IFLG)

                                                                        
C*******8/13/87**  Steve Rock ***************************************   
C Use the Bardin Calculation of Vertex terms instead of only electron   
C  and Tsai's Vertex term.                                              
C     FAC2=2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(-QSQ/EMSQ))      
C                                                                       
C  Vacuum terms of Bardin including 3 leptons and quarks                
      VAC_TERM = 2. * DELVAC(SNGL(-QSQ))                                 
C  Vertex term of Tsai (eq. 2.11)                                       
      VERTEX_TERM =TWO_ALPH_PI *(-1.+.75 *ALOG(-SNGL(QSQ)/SNGL(EMSQ)) ) 
      FAC2 = VAC_TERM + VERTEX_TERM                                     
C*******************************************************************    
                                                                        
                                                                        
      FAC=FAC1+FAC2                                                     
      FW1=FAC*W1                                                        
      FW2=FAC*W2                                                        
      A=OM*(EP-PVEC*COSP*COSK)                                          
      AP=OM*(ES-SVEC*COSS*COSK)                                         
                                                                        
C      Tsai on page 40 of SLAC-PUB-898 (1971) says that an uncertainty  
C      of zero / zero occurs when a = a' in his equation A.24, the      
C      integrand of which is calculated by this function.  To avoid this
C                                                                       
      FUNXL_139= 0.D0                                                       
      IF (A.EQ.AP) RETURN                                               
C                                                                       
      BP=-OM*PVEC*DSQRT((1.D0-COSP**2)*(1.D0-COSK**2))                  
      GNU=1.D0/(AP-A)                                                   
      X=DSQRT(A**2-BP**2)                                               
      Y=DSQRT(AP**2-BP**2)                                              
      TERM1=(A/X**3+AP/Y**3)*EMSQ*(2.D0*EMSQ+QSQ)+4.D0                  
     * +4.D0*GNU*(1.D0/X-1.D0/Y)*SP*(SP-2.D0*EMSQ)                      
     * +(1.D0/X-1.D0/Y)*(2.D0*SP+2.D0*EMSQ-QSQ)                         
      TERM21=2.D0*ES*(EP+OM)+.5D0*QSQ                                   
      TERM22=2.D0*EP*(ES-OM)+.5D0*QSQ                                   
      TERM23=2.D0*ES*EP-SP+OM*(ES-EP)                                   
      TERM24=2.D0*(ES*EP+ES*OM+EP**2)+.5D0*QSQ-SP-EMSQ                  
      TERM25=2.D0*(ES*EP-EP*OM+ES**2)+.5D0*QSQ-SP-EMSQ                  
      TERM26=EMSQ*(SP-OM**2)+SP*TERM23                                  
      TERM2=-A*EMSQ/X**3*TERM21-AP*EMSQ/Y**3*TERM22                     
     * -2.D0+2.D0*GNU*(1.D0/X-1.D0/Y)*TERM26+1.D0/X*TERM24              
     * -1.D0/Y*TERM25                                                   
      FUNXL_139=(FW2*TERM2+FW1*TERM1)*OM/(QSQ**2*(U0-UVEC*COSK)) 
c      If(FUNXL_139.GT.1.D15) THEN
c       BP= BP*(1.+1.D-20           )
c      ENDIF
      if(prttst) write(44,'(1x,7e12.5,10f7.3)')cosk,FUNXL_139,
     >   term1,fac,fac1,fac2
      RETURN                                                            
      END                                                               
!------------------------------------------------------------------------------------------

C $MEMBER=FUNXI, DATE=75061922, USErqsq=KCE                                
C      DOUBLE PRECISION FUNCTION FUNCOSK(COSK)                           
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL*8 MASS_MAX,MASS_MIN/1.0783D0/                                
C      COMMON/TAL1/ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1    
C      COMMON/TAL2/ZP,ZN,AT,IFLG                                         
C      COMMON/SOFT/ OM,B,SNSQ                                            
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      REAL*8 DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP                
C      COMMON/COSK/COSK_PASS                                             
C      EXTERNAL FUNMASS                                                  
C      DATA EMSQ/.26113D-6/,PI/3.1415926535D0/                           
C      DATA TM/.93828D0/                                                
C      DATA ALPHA/7.2973515D-3/                                          
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C      COSK_PASS = COSK                                                  
C                                                                        
C      MASS_MAX = DSQRT(USQ -2.D0 *(U0-UVEC*COSK)*DEE )                  
C                                                                        
C      FUNCOSK = QUADMO1(FUNMASS,MASS_MIN,MASS_MAX,3.D-4,NLVL)           
C                                                                        
C      RETURN                                                            
C      END                                                               
! correction routine from Blok/Vlados for 99118 use only
	SUBROUTINE F2corr1H(x,fact_cor)
***** theta in degree
        implicit none
      COMMON /KIN/   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2

 	REAL DIF1,DIF2,DIF3,DIF4,x,fact_cor,th1,th2,diffac
	REAL*4 F1,F2,DIFYGL,DIFTHETA
	REAL*4 GIPOT,COSYGL,GIPOSMALL
	integer i

	real P1(2)
	real P2(2)
	real P3(2)
	real P4(2)
 
*============ "START" THETA=10.60 ========================
 
 	 p1(1)=0.97094 
	 p1(2)=-0.98026E-01 
         dif1=p1(1)+p1(2)*x 
         	
*============ "FINISH" THETA=10.60 ========================	
	
*============ "START" THETA=14.60 ========================
		 
	 p2(1)=0.99776 
	 p2(2)=-0.34140E-01 
         dif2=p2(1)+p2(2)*x     
         	
*============ "FINISH" THETA=14.60 ========================	
	
*============ "START" THETA=18.60 =========================

 	 p3(1)=0.99493 
	 p3(2)=0.30907E-01 
         dif3=p3(1)+p3(2)*x  
         
*============ "FINISH" THETA=18.60 ========================	
	
*============ "START" THETA=22.60 =========================

 	 p4(1)=1.0284 
	 p4(2)=-0.47091E-01 
         dif4=p4(1)+p4(2)*x     
     
      
*============ "FINISH" THETA=22.60 ========================

 
         if(theta.lt.14.60) then
	      f1=dif1
	      f2=dif2
	      th1=10.60
	      th2=14.60
	 endif
	 
         if(theta.ge.14.60.and.theta.lt.18.60) then
	      f1=dif2
	      f2=dif3
	      th1=14.60
	      th2=18.60
	 endif

         if(theta.ge.18.60) then
	      f1=dif3
	      f2=dif4
	      th1=18.60
	      th2=22.60
	 endif
	       
	 if(f1.le.f2) difygl=abs(theta-th1)
	 if(f1.gt.f2) difygl=abs(theta-th2)

	 diffac=abs(f1-f2)	 
	 diftheta=th2-th1	 
         gipot=sqrt(diftheta**2+diffac**2)
	 cosygl=diftheta/gipot	 	   
         giposmall=difygl/cosygl
	 fact_cor=sqrt(giposmall**2-difygl**2)	 

	 if(f1.lt.f2) THEN
	  fact_cor=fact_cor+f1
	   else
	  fact_cor=fact_cor+f2
	 ENDIF   
	  
	 	          
! 	 WRITE(*,*)'x=:Factor=',x,fact_cor
  	 
	 do i=1,2
	  P1(i)=0.
	  P2(i)=0.
	  P3(i)=0.
	  P4(i)=0.
	 enddo
	 
	 dif1=0.
	 dif2=0.
	 dif3=0.
	 dif4=0.
	 f1=0.
	 f2=0.
	 
	 END


	SUBROUTINE F2corr2H(x,fact_cor)
***** theta in degree
        implicit none

 	REAL DIF1,DIF2,DIF3,DIF4,x,fact_cor,th1,th2,diffac
	REAL*4 F1,F2,DIFYGL,DIFTHETA
	REAL*4 GIPOT,COSYGL,GIPOSMALL
      COMMON /KIN/   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
        integer i

	real P1(2)
	real P2(2)
	real P3(2)
	real P4(2)
 
*============ "START" THETA=10.60 ========================
      	 
	 p1(1)=0.94718 
	 p1(2)=-0.15652E-01 
         dif1=p1(1)+p1(2)*x 
          	
*============ "FINISH" THETA=10.60 ========================	
	
*============ "START" THETA=14.60 ========================
	 
	 p2(1)=1.0024  
	 p2(2)=-0.30468E-01  
         dif2=p2(1)+p2(2)*x   
   	
*============ "FINISH" THETA=14.60 ========================	
	
*============ "START" THETA=18.60 =========================

	 
	 p3(1)=1.0327 
	 p3(2)=-0.32571E-01 
         dif3=p3(1)+p3(2)*x   
      
       
*============ "FINISH" THETA=18.60 ========================	
	
*============ "START" THETA=22.60 =========================

 	 p4(1)=1.0846 
	 p4(2)=-0.13985 
         dif4=p4(1)+p4(2)*x   
  
*============ "FINISH" THETA=22.60 ========================

 
         if(theta.lt.14.60) then
	      f1=dif1
	      f2=dif2
	      th1=10.60
	      th2=14.60
	 endif
	 
         if(theta.ge.14.60.and.theta.lt.18.60) then
	      f1=dif2
	      f2=dif3
	      th1=14.60
	      th2=18.60
	 endif

         if(theta.ge.18.60) then
	      f1=dif3
	      f2=dif4
	      th1=18.60
	      th2=22.60
	 endif
	       
	 if(f1.le.f2) difygl=abs(theta-th1)
	 if(f1.gt.f2) difygl=abs(theta-th2)

	 diffac=abs(f1-f2)	 
	 diftheta=th2-th1	 
         gipot=sqrt(diftheta**2+diffac**2)
	 cosygl=diftheta/gipot	 	   
         giposmall=difygl/cosygl
	 fact_cor=sqrt(giposmall**2-difygl**2)	 

	 if(f1.lt.f2) THEN
	  fact_cor=fact_cor+f1
	   else
	  fact_cor=fact_cor+f2
	 ENDIF   
			  
!* 	 WRITE(*,*)'x=:Factor=',x,fact_cor
  	 
	 do i=1,2
	  P1(i)=0.
	  P2(i)=0.
	  P3(i)=0.
	  P4(i)=0.
	 enddo
	 
	 dif1=0.
	 dif2=0.
	 dif3=0.
	 dif4=0
	 f1=0.
	 f2=0. 
	    	 
	 END

       SUBROUTINE D2MODEL_IOANA(QSQ,WSQ,W1,W2)
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
*****************************************************************************
        IMPLICIT NONE

        REAL*4  WSQ, QSQ
        REAL*4  R_NRES, SIG_RES(2), SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, 
     >          DIPOLE, K, NU,TAU, PI, AM, ALPHA, W1, 
     >          W2,CONV,sig_roper(2)
	integer i
	real*4  xval(37)
	logical goroper/.false./
*
        COMMON/ioana/xval
*       
*
* N.I. ...
*

        DATA PI/3.14159265/, AM/.93828/, ALPHA/0.00729735/,
     >          CONV/0.0025767/

*
! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.17) return
        IF(WSQ.LT.1.17.OR.WSQ.GT.5.OR.
     >    QSQ.GT.10.0) THEN
           do i=1,2	
              SIG_NRES(i) = 0.0
              SIG_RES1(i) = 0.0
              SIG_RES2(i) = 0.0
              SIG_ROPER(i)= 0.0 
           enddo
!          WRITE(*,*)'H2MODEL_IOANA called outside of kinematic range'
          RETURN
        ENDIF

!  returns transverse cross sections in units of
! microbarns/(dipole FF)**2
        CALL i_d2_model(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     &                  SIG_ROPER,goroper,xval)
*        write(*,*)'1',SIG_NRES
        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.
        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
*        write(*,*)'1',SIG_NRES(1),DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1)+sig_roper(1))*
     &                 DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES(1)
        SIGL = R_NRES*SIG_NRES(1)
*	write(*,*) 'I am here'

!        write(33,*)SIG_NRES(1),SIG_RES(1),R_NRES,DIPOLE,sigt
! 33     format(F15.3,2x,F15.3,2x,F15.3,2x,F15.3,2x,F15.3)      

! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)
*        write(*,*)'from program',w1,w2 

!*** PYB DIVIDE BY 2 TO GET PER NUCLEON
        W1=W1/2.
        W2=W2/2.

        RETURN
        END

         SUBROUTINE i_d2_model(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper,XVAL)
********************************************************************************
*
* This subroutine calculates model cross sections for deuterium in resonance region. 
* Cross section is returned in nanobarns. 
* For all the resonance regions we use nonrelativistc
* Breit-Wigner shapes (L.Stuart). 
* This subroutine is based on h2model.f from SLAC
* 
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
********************************************************************************
      IMPLICIT NONE
*
      REAL*4 XVAL(37)
* 
      logical goroper
      INTEGER I, J
*
* modified by I.N. 04/09/97
*

       
      REAL*4 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     KRCM, EPIRCM, PPIRCM,   
     >     W, KCM, K, EPICM, PPICM, WDIF 
 
      REAL*4 PI, ALPHA, AM,
     >     MPPI, MPI,  
     >     XR  

      integer jj
      real*4 kr,meta
      real*4 mrho,mn,am2
      real*4 gamma_gamma,gamma_pi,denom_pi
      real*4 sigma_delta,sigma_second,sigma_third
*
      real*4 mass(3),width(3),qdep(3)
      real*4 coeff(3,0:3),nres_coeff(0:4,3)
* Define vectors : MASS(3), WIDTH(3), which will contain the masses and widths of
* the three resonances: 1=DELTA, 2=second res reg, 3=third res reg
*
      DATA PI/3.14159265/
      DATA ALPHA/0.00729735/
      DATA AM/.93828/
      DATA MPI/0.13957/
      DATA MPPI/1.07783/ 
      data Meta/0.54745/ 
      data Mrho/0.7681/
*      DATA XR/0.18/
	data Mn/0.938/
      Am2 = Am*Am
*
* The fit parameters are called  XVAL 
* Assign them to masses, widths, etc.
*
      xr=XVAL(1)
      do i=1,3
         mass(i)  = XVAL(i+1)
*         write(*,*)xval(i)
      enddo
      do i=1,3
         width(i)  = XVAL(i+4)
*         write(*,*)xval(i+3)
      enddo
      do i=1,3
         qdep(i)   = XVAL(i+7)
*         write(*,*)xval(i+6)
      enddo     
      k=11
      do i=1,3
         do j=0,3
            coeff(i,j) = XVAL(k)
*            write(*,*)xval(k)
            k=k+1
         enddo
      enddo
      do i=0,4
         do j=1,3
           nres_coeff(i,j) = XVAL(k)
***!!!            write(*,*)'3',xval(k)
           k=k+1
         enddo
      enddo
*      write(*,*)'Done assigning parameters'
*
*      write(*,*) 'w2=', w2
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
****************************************************************
* This part is not used for deuterium.
****************************************************************
* K equivalent real photon energy needed to excite resonance W
* K is in Lab frame, KCM is in CM frame
* 
         K = (W2 - Am2)/2./Am
         KCM = (W2 - Am2)/2./W
         KR = (Mass(1)*Mass(1)-Am2)/2./Am
         KRCM = (Mass(1)*Mass(1) - Am2)/2./Mass(1)
* Decay mode for Delta: N,pi:
********
      EPICM = 0.5*( W2 + Mpi*Mpi - Am2 )/W

      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      EPIRCM = 0.5*(Mass(1)*Mass(1) + Mpi*Mpi - Am2 )/Mass(1)   
      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      
********
* Now calculate the partial widths:
*  gamma_pi(i)
*      write(*,*)'mass(1),w2 ',mass(1),' ' ,w2
*      write(*,*)'width(1) ',width(1)
*      write(*,*)'PPICM, PPIRCM ',PPICM,' ',PPIRCM
*      write(*,*)'XR ' ,xr
      gamma_pi = width(1)*(PPICM/PPIRCM)**3*(PPIRCM*PPIRCM+XR*XR)/
     >          (PPICM*PPICM+XR*XR)
*
*      write(*,*)'iuhu 2'
      gamma_gamma = width(1)*(KCM/KRCM)**2*(KRCM*KRCM+XR*XR)/
     >          (KCM*KCM+XR*XR)

      denom_pi = (W2 - Mass(1)*Mass(1))**2 + (Mass(1)*gamma_pi)**2
*      write(*,*)'iuhu 3'
***********************************************************************
***********************************************************************
* Now calculate cross sections for Delta:
      sigma_delta = width(1)/((W-Mass(1)*(1. + Q2*qdep(1)))**2 
     >               + 0.25*width(1)*width(1))
* For the second and third resonance regions:
      sigma_second = width(2)/((W-Mass(2)*(1. + Q2*qdep(2)))**2 
     >               + 0.25*width(2)*width(2))

*      write(*,*)'width(2)',' ',width(2)
*      write(*,*)'mass(2)',' ',mass(2)

*      write(*,*)'sigma_second',' ',sigma_second

      sigma_third  = width(3)/((W-Mass(3)*(1. + Q2*qdep(3)))**2
     >               + 0.25*width(3)*width(3))

*      write(*,*)'sigma_third',' ',sigma_third

*      write(*,*)'width(3)',' ',width(3)
*      write(*,*)'mass(3)',' ',mass(3)

*      write(*,*)'iuhu 5'
*
* Put in the Q2 dependence of AH(Q2):
* AH(Q2) = coeff(i,j)*Q2**j
*
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(1) = 0.0
      do j=0,3  
         sig_res1(1) = sig_res1(1)+ sigma_delta*coeff(1,j)*Q2**j
         sig_res2(1) = sig_res2(1)+
     >               sigma_second*coeff(2,j)*Q2**j+
     >               sigma_third *coeff(3,j)*Q2**j
      enddo   
*      write(*,*)'sigma_d,sigma_23',' ',sig_res1(1),' ',sig_res2(1)
      do j=0,4
         do jj=1,3
            sig_nres(1) = sig_nres(1)+
     >          nres_coeff(j,jj)*q2**j*sqrt(wdif**(2*jj-1))
*          write(*,*)'2',sig_nres(1),nres_coeff(j,jj),j,wdif,jj
         enddo
      enddo     
*
      RETURN
      END

      REAL FUNCTION  ROCK_RES9(TARG,X,Q2)
! April 22 2003 Version of deuterium fit.
! No JLAB data

      IMPLICIT NONE
      INTEGER TARG ! 1=H, 2=D
      REAL Q2
      REAL X,F2INEL
      REAL*8 WM,B03
      CHARACTER*1 TARG_STR ! D or H      
      REAL MP2/.8803/,MP/.93828/ 
! DATE =2003/ 4/22  TIME=13: 3
! MODEL# 9
! HYDROGEN
      REAL*8 CH(40)/
     >  1.074100000, 0.279494473, 4.712278073, 1.739263857, 1.431155682,
     >  0.184080601, 1.225717385, 0.110613494, 0.214285318, 1.510392164,
     >  0.105595846, 0.309575075, 1.725500000, 0.138478762,-1.452323856,
     >  1.953800000, 0.342265136,-0.044517487, 0.093869839,-0.029404837,
     >  1.199536964, 0.671722156, 2.033898985, 0.024107348, 0.452992650,
     >  0.624721200, 0.054377925, 0.099384033, 0.897445331,13.205228916,
     >  6.284002261, 0.013794386, 0.030764850,-0.007500832, 0.025566716,
     >  0.016759263, 0.966305525, 1.058950737, 0.000000000, 0.000000000/
     > 
!deuterium
      REAL*8 CD(40)/
     >   1.05220000,  0.75530000,  3.30000000,  1.70000000,  3.50000000,
     >   1.04000000,  1.22000000,  0.10000000,  0.48000000,  1.50000000,
     >   0.08000000,  0.65000000,  1.70000000,  0.12000000,  0.74000000,
     >   1.90000000,  0.19000000, -0.17000000,  0.00960000,  0.03500000,
     >   3.50000000, -0.60000000,  4.70000000,  0.41100000,  0.41100000,
     >   0.10000000,  0.41100000,  0.10000000,  0.10000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000/
    
      IF(TARG.EQ.1) TARG_STR ='H'
      IF(TARG.EQ.2) TARG_STR ='D'
      WM = SQRT(MP2 +Q2*(1./X -1.) )

      CALL F2GLOB_FAST(X,Q2,TARG_STR,9,F2INEL)

      IF(TARG.EQ.1) THEN  ! H2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CH)   
      ELSE     ! D2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CD)   
      ENDIF
      RETURN
      END
!----------------------------------------------------------------

      REAL*8 FUNCTION B03(WM,QSQ,NC2,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(80),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,BRES
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.8803/, PM2/1.876512/, PM/.93828/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      INTEGER   NC2 !offset for coeficients (0 for h2)
      DATA      LSPIN/1,2,3,2/                                       

C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(NC2+22)+C(NC2+23)*(X-C(NC2+24))**2                                 
      PIEMSQ = (C(NC2+1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(NC2+1)) B1 = C(NC2+2)                                         
      EB1 = C(NC2+3)*(WM-C(NC2+1))                                              
      IF( (EB1.LE.25.).AND.(B1.NE.0.)) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(NC2+4)) B2 = (1.-C(NC2+2))                                    
      EB2 = C(NC2+5)*(WSQ-C(NC2+4)**2)                                          
      IF( (EB2.LE.25.0).AND.(B2.NE.0.)) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(NC2+2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(NC2+INDEX)                                            
           IF (I.EQ.1) RAM=C(NC2+INDEX)+C(NC2+18)*SQRT(QSQ)! +C(NC2+30)*QSQ**2
     >       + C(NC2+25)/(QSQ+C(NC2+26))
           RMA    = C(NC2+INDEX+1)
           IF (I.EQ.2) RAM =RAM + C(NC2+27)/(C(NC2+28)+QSQ)
           IF (I.EQ.3) THEN
              RMA=RMA*(1.D0 +C(NC2+20)/(1.D0+C(NC2+21)*QSQ))
              RAM =RAM + C(NC2+19)/(C(NC2+29)+QSQ) 
           ENDIF
           IF (I.EQ.4) RAM =RAM + C(NC2+30)/(C(NC2+31)+QSQ)
           RWD    = C(NC2+INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
                                                                        
C FORM VW2/F2                                                           
                                                                        
      B03 = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
                                                                        
      RETURN                                                            
      END             



! -*-Mode: Fortran; compile-command: "f77 -o f2glob.o -c f2glob.f"; -*-
!File E.13. F1990.FORTRN.                                               
!Reference:  L.W.Whitlow, SLAC-Report-357,                              
!            Ph.D. Thesis, Stanford University,                         
!            March 1990.                                                
!For details see file HELP.DOCUMENT.                                    
                                                                        
!Program contains 145 lines of Fortran code, of 72 characters each, with
!no subroutines.  Program requires File E.14 as input.                  
                                                                        
                                                                        
      SUBROUTINE F2GLOB_FAST(X,Q2,Target,MODEL,F2)
                                                                        
! Returns F2 and related quantities from the either the LAMBDA12 model  
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.           
!                                                                       
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2  
!                                                                       
! Further, program returns uncertainty in F2 based on both statistics   
! (ST) and systematic effects due to our choice of models (SY).         
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical
! uncertainty in this slope.                                            
!                                                                       
! Best model is LAMBDA12.  SY is estimated to be the difference between 
! the two models.                                                       
!                                                                       
! Errors due to overall normalization are not included in program output
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.      
! Systematic errors due to radiative corrections are shown in Reference 
! to be very kinematic independent, and are everywhere <.5%, and thus,  
! are ignored by this routine (also see documentation to dFRC in file   
! HELP.DOCUMENT).                                                       
!                                                                       
! Coefficients and correlation matrix elements are from File            
! E.13 F1990.MATRICES, dated 27Jan90.                                   
                                                                        
      IMPLICIT NONE                                                     
      LOGICAL FIRST/.TRUE./                                     
      REAL   X, Q2, F2,ST
      REAL   XPP,XP,Y,POLY,Q,QTH,DQ,F2TH,QUAD,SCB      
      REAL   BINDING                                               
      INTEGER MODEL, I, J, K                                            
      CHARACTER*1 TARGET                                                
      REAL*8    B(2,9),L(2,9,9)                      ! OMEGA9 variables 
      REAL*8    C(2,12),M(2,12,12)                   ! LAMBDA12 variable
      REAL*8 LIN                                       
                                                                        
                                                                        
! Model #9  27 Jan 90.                                                  
      DATA (B(1,J),J=1,9)/                        !  HYDROGEN Ci        
     >   0.7338659870D0,  11.0245522588D0,   2.6185804129D0,            
     >   4.0956321483D0,   0.1206495422D0,   1.9714128709D0,            
     >   3.8893348719D0, -14.0507358314D0,   8.8080576075D0/            
      DATA ((L(1,J,K),K=1,9),J=1,9)/              !  HYDROGEN MijDiDj   
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,             
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,             
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,             
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,             
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,             
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,             
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,             
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,             
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,             
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,             
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,             
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,             
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,             
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,             
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,             
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,             
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,             
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,             
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,             
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,             
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,             
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,             
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,             
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,             
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,             
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,             
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/             
                                                                        
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci       
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,             
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,             
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/             
                                                                        
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj  
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,             
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,             
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,             
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,             
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,             
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,             
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,             
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,             
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,             
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,             
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,             
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,             
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,             
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,             
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,             
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,             
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,             
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,             
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,             
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,             
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,             
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,             
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,             
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,             
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,             
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,             
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/             
                                                                        
!    MODEL #12:    27Jan90.                                             
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci            
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,             
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,             
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,             
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/             
      DATA ((M(1,J,K),K=1,12),J=1,12)/     !     HYDROGEN MijDiDj       
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,             
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,             
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,             
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,             
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,             
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,             
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,             
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,             
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,             
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,             
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,             
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,             
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,             
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,             
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,             
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,             
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,             
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,             
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,             
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,             
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,             
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,             
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,             
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0,             
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,             
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,             
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,             
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,             
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,             
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,             
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,             
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,             
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,             
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,             
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,             
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,             
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,             
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,             
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,             
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,             
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,             
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,             
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,             
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,             
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,             
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,             
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,             
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/             
                                                                        
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci           
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,             
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,             
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,             
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/             
      DATA ((M(2,J,K),K=1,12),J=1,12)/     !     DEUTERIUM MijDiDj      
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,             
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,             
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,             
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,             
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,             
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,             
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,             
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,             
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,             
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,             
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,             
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,             
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,             
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,             
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,             
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,             
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,             
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,             
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,             
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,             
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,             
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,             
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,             
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0,             
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,             
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,             
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,             
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,             
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,             
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,             
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,             
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,             
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,             
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,             
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,             
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,             
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,             
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,             
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,             
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,             
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,             
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,             
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,             
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,             
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,             
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,             
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,             
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/               
      !---------------------------------------------------------------- 
      !---------------------------------------------------------------- 
                                                                        
                                                                        
                                                                        
      i = 1                                                             
      IF (TARGET.EQ.'D') i = 2                                          
      BINDING = 1./(1.-EXP(-MIN(20.,7.7*(1./X+.93828**2/Q2-1.))))       
      IF (i.EQ.1) BINDING = 1.                                          
                                                                        
      !OMEGA9 MODEL FIRST:
      IF(MODEL.EQ.9) THEN                                        
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                             
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                             
           Y    = 1.-XP                                                 
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+      
     >            B(i,9)*Y**7                                           
           F2  = X/XPP*BINDING*POLY                                    

      ELSEIF(MODEL.EQ.12) THEN

      !LAMBDA12 MODEL NEXT:                                             
           Y    = 1.-X                                                  
           q    = LOG(Q2)                                               
           qth  = .2+3.2*X                                              
           dq   = q-qth                                                 
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                  
     >            C(i,4)*Y**6+C(i,5)*Y**7                               
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                   
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq       
           IF (q.GT.qth) QUAD = 0.                                      
           SCB  = (1.+LIN+QUAD)                                         
           F2  = F2th*SCB*BINDING                                      
      ELSE                                                              
           WRITE(*,'('' F2GLOB: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')  
           F2=0.
           ST=-1.  
           RETURN
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               


      SUBROUTINE tunatarg(diam,thick,ztar,trad,cryo_cm, wall_cm)

c inputes: diam: cell diameter in cm
c          thick: cell wall thickness in cm
c          ztar in cm (0 is center of target)
c          trad in radians
c output: cryo_cm cm of cryogen for scttered electron
c         wall_cm amount of Al wall in cm

C This subroutine calculates the amount of target material in cm
C that is transversed in the cryogen and cell wall of the Hall C 
C tuna can target for an electron detected at a scattering angle, theta.
C This doesn't currently include the beam on target position.
C E. Christy 6/05. Modified by P. Bosted 6/05.

      IMPLICIT none

      integer i
      real diam,R,R1,thick,trad,cryo_cm,wall_cm
      real ztar,zmax,zmin,zint,xint,zdiff,li,lf(2)
      real a,b,c,ctan,ctan2

      R1 = diam/2.
      
      zmax =  diam/2 - 1.e-4
      zmin = -diam/2 + 1.e-4
      if(ztar.GT.zmax) ztar = zmax  !!! don't go outside target wall  !!!
      if(ztar.lt.zmin) ztar = zmin  !!! don't go outside target wall  !!!

      ctan = 1/tan(trad)
      ctan2 = ctan*ctan
      a = 1. + ctan2
      b = 2.*ztar*ctan
     
      do i=1,2      !!!  Loop over inner and outer radius of cell wall  !!!
       R = R1
       if(i.eq.2) R = R1 + thick
       c = ztar*ztar - R*R
       xint = ((-1.)*b + sqrt(b*b - 4.*a*c))/2./a !!! z intersection of line and circle     !!!   
       zint = sqrt(R*R - xint*xint)  !!! x intersection of line and circle     !!! 
       if(xint.LT.0) then   !!!  choose other solution  !!!
         xint = ((-1.)*b - sqrt(b*b - 4.*a*c))/2./a !!! z intersection of line and circle   !!! 
         zint = sqrt(R*R - zint*zint)  !!! x intersection of line and circle     !!! 
       endif
       zdiff = zint - ztar           !!! z distance transversed                !!!  
       lf(i) = sqrt(xint*xint + zdiff*zdiff)  !!!  total distance transversed  !!!
      enddo
      li = ztar + R1

      wall_cm = lf(2) - lf(1)
      cryo_cm = lf(1)

      return
      end

      subroutine christy(W2,Q2,F1,R)
!------------------------------------------------------------------------
! subroutine to return proton structure function F1 and ratio R=sigl/sigt
!
! inputs are electron missing mass squared W2 (GeV**2) 
!            momentum trasfer squared Q2 (GeV**2)
! inputs and outputs are Real*8
! the file christy.dat is needed to use this subroutine
! Note: if W2<1.155 GeV**2, values of zero are returned (below threshold)
! Fit done by Eric Christy 11/04
! Reference this fit as ???       
!------------------------------------------------------------------------
      IMPLICIT NONE

      real*8 w2,q2,xval1(41),xvall(41),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F1,FL,R
      integer i
      logical first/.true./
 
      mp = .93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.

      if(first) then
        first=.false.
        open(unit=15,file='christy.dat',status='old')
        do i=1,41
          read(15,*) temp(1)  ! par #
          read(15,*) XVAL1(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        do i=1,41
          read(15,*) temp(1)  ! par #
          read(15,*) XVALL(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(15)
      endif

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2/(w2+q2-mp2)
      if(xb.le.0.0) return

      call resmod(1,w2,q2,XVAL1,F1)
      call resmod(2,w2,q2,XVALL,FL)
      if(F1.le.0.0) return

       if(F1.le.0.0) return

      F1 = F1/8.d0/pi/pi/alpha/1000.d0
      F1 = F1*(w2-mp2)
      FL = FL/8.d0/pi/pi/alpha/1000.d0
      FL = FL*(w2-mp2)*2.d0*xb
      R = FL/ (2.D0 * XB * F1)

      return
      end

      SUBROUTINE RESMOD(sf,w2,q2,xval,fn) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,fn,xval(41),mass(4),width(4)
      REAL*8 height(4),fn_del,fn_s11,fn_f15,rescoef(4,3)
      REAL*8 nr_coef(4,4),wdif,fn_nr,fn_4,w2temp,wtemp
      REAL*8 roper_mass,roper_width,roper_height
      REAL*8 roper_mparm,roper_exp,mq2(4)
      REAL*8 alpha,pi
      INTEGER i,j,num,sf


      mp = 0.93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036
      W = sqrt(w2)
      w2temp = w2 - .93828*.93828
      wtemp = sqrt(w2temp)
      wdif = w - (0.937272 + 0.137)
      xb = q2/(q2+w2-mp2)
      if(sf.EQ.1) xb = 1.d0
      fn_nr = 0.0d0
      do i=1,4
        height(i) = 0.
      enddo
      
      do i=1,4
        mass(i) = xval(i)
      enddo

      do i=1,4
        mq2(i) = xval(36 + i)
      enddo

      if(q2.GT.0.1) then

        mass(2) = mass(2)*exp(-(q2-0.1)/mq2(2)) 
     &            + mq2(1)*(1.-exp(-(q2-0.1)/mq2(2))) 

        mass(3) = mass(3)*exp(-(q2-0.1)/mq2(4)) 
     &            + mq2(3)*(1.-exp(-(q2-0.1)/mq2(4)))

      endif   
      do i=1,4
        width(i) = xval(4+i)
      enddo
      num = 0
      do i=1,4
       do j=1,3
         num = num + 1
         rescoef(i,j)=xval(8 + num)
c           write(6,*) i,j,num,rescoef(i,j)
c           height(i) = height(i)+rescoef(i,j)*q2**(-1.*(float(j-1)))
         enddo
         height(i) = rescoef(i,1)*
     &             (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
c         if(w2.LT.1.35) height(i) = height(i)/(w2-1.30)    
      enddo
      num = 0     
      do i=1,4
       do j=1,4
         num = num + 1
         nr_coef(i,j)=xval(20 + num)
c         write(6,*) i,j,num,nr_coef(i,j)         
       enddo
      enddo

      do i=1,5
       roper_mass = xval(37)
       roper_width = xval(38)
       roper_height = xval(39)
       roper_mparm = xval(40)
       roper_exp = xval(41)
      enddo
c      write(6,*) "constant coef's are:  ",height

CC   Calculate Breit-Wigners for the 3 resonance regions   CC

      fn_del = width(1)/((W-mass(1))**2 
     &               + 0.25*width(1)*width(1))
      fn_s11 = width(2)/((W-mass(2))**2 
     &               + 0.25*width(2)*width(2))
      fn_f15 = width(3)/((W-mass(3))**2 
     &               + 0.25*width(3)*width(3))
      fn_4   = width(4)/((W-mass(4))**2 
     &               + 0.25*width(4)*width(4))

      fn_del = height(1)*fn_del
      fn_s11 = height(2)*fn_s11
      fn_f15 = height(3)*fn_f15
      fn_4   = height(4)*fn_4

c      roper_height = roper_height*(1.+q2/roper_mparm)**(-1.*roper_exp)
c      fn_roper = roper_width/((W-roper_mass)**2 
c     &               + 0.25*roper_width*roper_width) 
c      fn_roper = roper_height*fn_roper


      do i=1,4
       do j=1,4
c         fn_nr = fn_nr+
c     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))

         fn_nr = fn_nr+
     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))
     &                  /w2temp/xb

c         write(6,*) sf

c         if(sf.EQ.2) fn_nr = fn_nr/xb
                         
       enddo
      enddo

      fn = fn_del + fn_s11 + fn_f15 + fn_4 + fn_nr
      if(sf.EQ.2) then
        fn = fn*(1.-exp(-q2/roper_mass))
      endif

c      write(6,*) "IN model:  ",sf,w2,q2,fn

      RETURN 
      END 

      subroutine christy705(W2,Q2,F1,R)
!------------------------------------------------------------------------
! subroutine to return proton structure function F1 and ratio R=sigl/sigt
!
! inputs are electron missing mass squared W2 (GeV**2) 
!            momentum trasfer squared Q2 (GeV**2)
! inputs and outputs are Real*8
! the file christy.dat is needed to use this subroutine
! Note: if W2<1.155 GeV**2, values of zero are returned (below threshold)
! Fit done by Eric Christy 7/05
! Reference this fit as ???       
!------------------------------------------------------------------------
      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F1,FL,R
      integer i
      logical first/.true./
 
      mp = .93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.

      if(first) then
        first=.false.
        open(unit=15,file='f1parms.dat',status='old')
        do i=1,50
          read(15,*) temp(1)  ! par #
          read(15,*) XVAL1(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(unit=15)
        open(unit=15,file='flparms.dat',status='old')
        do i=1,50
          read(15,*) temp(1)  ! par #
          read(15,*) XVALL(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(15)
      endif

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2/(w2+q2-mp2)
      if(xb.le.0.0) return

      call resmod705(1,w2,q2,XVAL1,F1)
      call resmod705(2,w2,q2,XVALL,FL)
      if(F1.le.0.0) return

       if(F1.le.0.0) return

      F1 = F1/8.d0/pi/pi/alpha/0.389d3
      F1 = F1*(w2-mp2)
      FL = FL/8.d0/pi/pi/alpha/0.389d3
      FL = FL*(w2-mp2)*2.d0*xb
      R = FL/ (2.D0 * XB * F1)

      return
      end

CCC  Version 072205  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMOD705(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.93828
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w


      if(sf.EQ.1) then
        q2low = 0.15
      else
        q2low = 0.05
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + 0.1
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

        enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             
          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.61)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        if(sf.EQ.1) then 
         br_21_1 = 0.5
         br_21_2 = 0.5
        else
          br_21_1 = 0.985
          br_21_2 = 1.-br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

          width(6)=intwidth(6)*
     &     (xval(48)*petacm/petacmr(5)+(1.-xval(48))*ppi2cm/ppi2cmr(5))

          sig_4 =  (xval(48)*petacm+(1.-xval(48))*ppi2cm)/2./kcm/
     &           ((W2 - mass(6)**2.)**2. + (mass(6)*width(6))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

c        if(sf.EQ.2) then
c          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
c        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.(mp+mpi)**2.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 


CCC  Version 031606  -  Author:  M.E. Christy                        CCC
C*** changed to 7/7/06 version (see below)
CCC  Subroutine to get Transvese and Longitudinal eP cross sections  CCC 
CCC  from fits to L/T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


      SUBROUTINE CHRISTY31606(W2,Q2,F1,R)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval1 / 
c    &  0.12291E+01,0.15000E+01,0.15117E+01,0.16926E+01,0.16200E+01,
c    &  0.14261E+01,0.14100E+00,0.23398E+00,0.81124E-01,0.96627E-01,
c    &  0.18885E+00,0.26660E+00,0.85509E+02,0.10378E+02,0.34206E+00,
c    &  0.38546E+01,0.73111E+01,0.57629E+00,0.63691E-01,0.24072E+01,
c    &  0.19496E+02,0.38962E+01,0.25672E+00,0.29348E+01,0.10138E+02,
c    &  0.86788E+01,0.25859E+00,0.29227E+01,0.28345E+02,0.48552E+01,
c    &  0.13768E+00,0.49746E+01,0.36645E+01,0.96196E+01,0.29507E+00,
c    &  0.23934E+01,0.28438E+03,0.10536E+00,0.12390E+01,0.14593E+00,
c    &  -.24992E+03,0.67379E+00,0.23922E+01,-.24997E+00,-.18421E-02,
c    &  0.63658E-01,0.19539E+01,0.25343E+00,0.14434E+02,0.10982E+00 /
C  Simona new - h2model 7/7/06
c    & 0.12293E+01,0.15000E+01,0.15134E+01,0.16932E+01,0.16150E+01,
c    & 0.14201E+01,0.14247E+00,0.23399E+00,0.81908E-01,0.98779E-01,
c    & 0.14974E+00,0.25805E+00,0.86182E+02,0.10540E+02,0.31745E+00,
c    & 0.38879E+01,0.13214E+02,0.64008E+00,0.83161E-01,0.26896E+01,
c    & 0.18654E+02,0.45807E+01,0.24799E+00,0.30223E+01,0.10005E+02,
c    & 0.10525E+02,0.27374E+00,0.30040E+01,0.21101E+02,0.50802E+01,
c    & 0.12346E+00,0.51132E+01,0.17793E+01,0.21694E+02,0.27851E+00,
c    & 0.24741E+01,0.27502E+03,0.92225E-01,0.12263E+01,0.13867E+00,
c    & -.29351E+03,0.75890E+00,0.26699E+01,-.35054E+00,-.12799E-02,
c    & 0.73801E-01,0.20070E+01,0.51441E+00,0.30894E+02,0.94466E-01 /
C  Simona new - Christy
c    & 0.12293E+01,0.15000E+01,0.15133E+01,0.16927E+01,0.16150E+01,
c    & 0.14141E+01,0.14185E+00,0.23400E+00,0.78118E-01,0.95281E-01,
c    & 0.16498E+00,0.24073E+00,0.86154E+02,0.10245E+02,0.35896E+00,
c    & 0.38298E+01,0.12709E+02,0.67233E+00,0.85973E-01,0.25172E+01,
c    & 0.17335E+02,0.45604E+01,0.23752E+00,0.30148E+01,0.98678E+01,
c    & 0.94126E+01,0.27697E+00,0.29372E+01,0.23114E+02,0.49613E+01,
c    & 0.20701E+00,0.48640E+01,0.17441E+01,0.19868E+02,0.27348E+00,
c    & 0.24772E+01,0.27779E+03,0.89228E-01,0.12219E+01,0.14308E+00,
c    & -.29694E+03,0.75738E+00,0.26442E+01,-.35376E+00,-.11288E-02,
c    & 0.78975E-01,0.20200E+01,0.58638E+00,0.34075E+02,0.89893E-01 /
c Iteration of July 26, 2006
     & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
     & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
     & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
     & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
     & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
     & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
     & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
     & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
     & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
     & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


      data xvalL/
c    &  0.12438E+01,0.14948E+01,0.13100E+01,0.17300E+01,0.16700E+01,
c    &  0.13950E+01,0.76008E-01,0.53212E-01,0.14228E+00,0.82749E-01,
c    &  0.60000E-01,0.42000E+00,0.28298E+02,0.36911E+04,0.98237E+04,
c    &  0.11796E-03,0.63566E+01,0.26476E+05,0.36530E+05,0.33431E-03,
c    &  0.95795E+05,0.68291E+04,0.35634E+04,0.90878E+05,0.12332E+03,
c    &  0.73314E+06,0.26074E+06,0.74001E+02,0.00000E+00,0.92972E+04,
c    &  0.87089E+04,0.46171E+05,0.10478E+02,0.73525E+04,0.95318E+04,
c    &  0.29686E+00,0.69642E+03,0.13502E+01,0.00000E+00,0.34422E+01,
c    &  0.19016E+01,0.20149E+00,0.62176E+03,0.61984E+01,0.16046E+00,
c    &  0.16765E+01,0.00000E+00,0.40092E-11,0.71919E+00,0.33674E+00 /
C  Simona new - h2model 7/7/06
c    & 0.12437E+01,0.14943E+01,0.12808E+01,0.17330E+01,0.16700E+01,
c    & 0.13921E+01,0.71753E-01,0.45000E-01,0.57526E-01,0.73573E-01,
c    & 0.60000E-01,0.47000E+00,0.28690E+02,0.73923E+04,0.20383E+05,
c    & 0.00000E+00,0.53160E+01,0.49592E+05,0.64624E+05,0.00000E+00,
c    & 0.12151E+05,0.10000E+05,0.31121E+04,0.81442E+05,0.41594E+04,
c    & 0.38606E+06,0.96177E+05,0.42609E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17969E+02,0.69923E+04,0.75455E+04,
c    & 0.79081E+00,0.69950E+03,0.13591E+01,0.00000E+00,0.34511E+01,
c    & 0.19061E+01,0.22893E+00,0.10577E+04,0.79123E+01,0.12661E+00,
c    & 0.15646E+01,0.00000E+00,0.00000E+00,0.66613E+00,0.38768E+00 /
C  Simona new - Christy
c    & 0.12441E+01,0.14943E+01,0.12820E+01,0.17304E+01,0.16700E+01,
c    & 0.13955E+01,0.73364E-01,0.45064E-01,0.68170E-01,0.76574E-01,
c    & 0.60000E-01,0.46067E+00,0.29174E+02,0.73983E+04,0.20326E+05,
c    & 0.00000E+00,0.53291E+01,0.49445E+05,0.64863E+05,0.00000E+00,
c    & 0.13790E+05,0.10000E+05,0.37341E+04,0.61390E+05,0.44336E+04,
c    & 0.37016E+06,0.10266E+06,0.39083E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.16753E+02,0.70252E+04,0.75339E+04,
c    & 0.77189E+00,0.70838E+03,0.13292E+01,0.00000E+00,0.34217E+01,
c    & 0.19038E+01,0.21016E+00,0.98783E+03,0.80689E+01,0.12300E+00,
c    & 0.15563E+01,0.00000E+00,0.00000E+00,0.65655E+00,0.36031E+00 /
C  Simona - July 28, 2006
     & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
     & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
     & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
     & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
     & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
     & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
     & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
     & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
     & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
     & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmod316(1,w2,q2,xval1,sigT)
      call resmod316(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

   
      end

         SUBROUTINE H2MODEL(QQ4,WW4,W1,W2)

********************************************************************************
*
* This subroutine calculates model cross sections for inclusive electron-proton
* scattering in the resonance region. The cross section is returned in 
* nanobarns/sr GeV. This fit is a modified version of Linda Stuart's 8/91 fit. 
* One major difference is that the coefficients represent fit results from a 
* substantially expanded data set, including all inclusive SLAC data in the 
* range 1.1 < W^2 < 5. There are other differences; for a complete discussion, 
* see Keppel's Ph.D. thesis. 2/94 CEK
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
* SIG1     = Cross section in nb/sr/GeV**2 (DSIG/DOMEGA/DW**2)
* SIG2     = Cross section in nb/sr/GeV (DSIG/DOMEGA/DEP)
* SIG_NRES = Non-resonant contribution to SIG1.
* SIG_RES  = Resonant contribution to SIG1.
* SIG_NR   = Non-resonant contribution to SIG2.
* SIG_R    = Resonant contribution to SIG2.
* SIGroper = Possible Roper contribution to SIG2.
* goroper  = Logical variable, set true if including possible Roper strength
*
* SIG1, SIG2, SIG_NRES, and SIG_RES are 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
        IMPLICIT NONE

        logical goroper 
        logical goodfit 
        INTEGER I
        REAL*4  QQ4,WW4,W1,W2
        REAL*8  SIN2, SIG1(2), SIG2(2), SIG_RES(2), SIG_R(2), 
     >          SIG_RES1(2), SIG_RES2(2), SIG_NRES(2), SIG_NR(2), 
     >          DIPOLE, COS2, TAU, EPS, K, DW2DEP, DEPCONV, 
     >          PI, AM, ALPHA,NU,CONV,SIGT,SIGL,qq,ww,
     >          RADCON, R_NRES, sigroper(2), x, r, dr,fac,rlog
       
        qq = qq4
        ww = ww4

        RADCON = 3.141592654/180.
        PI = 3.14159265
        AM = .9382727
        ALPHA = 0.00729735
        CONV = 0.0025767

C        SIN2 = SIN(TH*RADCON/2.0)*SIN(TH*RADCON/2.0)
C        COS2 = 1 - SIN2
C        Q2 = 4.0*E*EP*SIN2
C        W2 = AM*AM + 2.0*AM*(E - EP) - 4.0*E*EP*SIN2

         goroper = .true.

C        write(6,*) q2,w2

        IF(WW.LT.1.15) THEN            ! Below pion threshold
          DO I = 1,2
            SIG1(I) = 0.0
            SIG2(I) = 0.0
            SIG_NRES(I) = 0.0
            SIG_RES(I) = 0.0
           SIG_NR(I) = 0.0
            SIG_R(I) = 0.0
            sigroper(I) = 0.0
          ENDDO
          RETURN
        ENDIF

        K = (WW - AM*AM)/(2.0*AM)
        NU = K + QQ/(2.0*AM)
        TAU = NU*NU/QQ
        x = qq/2./AM/nu
C        EPS = 1.0/(1.0 + 2.0*(1.0 + TAU)*SIN2/COS2)
        DIPOLE = 1.0/(1.0 + QQ/0.71)**2
c        DEPCONV = ALPHA*K*EP/(4.0*PI*PI*Q2*E)*(2.0/(1.0 - EPS))*1000.
c        DW2DEP = 2.0*AM + 4.0*E*SIN2

! H2MOD_FIT returns cross sections in units of microbarns/(dipole FF)**2
        CALL H2MODEL_FIT(QQ,WW,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)
c        call r1990(x,qq,r,dr,goodfit)
c        if (goodfit) then
c         R_NRES = r
c        else 
         R_NRES = 0.25/SQRT(QQ)
c        endif

        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
C        SIG_NRES(2) = SIG_NRES(2)*DIPOLE*DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1))
     >                 *DIPOLE*DIPOLE
C        SIG_RES(2)  = SQRT(SIG_RES1(2)**2 + SIG_RES2(2)**2 + 
C     >                 sigroper(2)**2)*
C     >                 DIPOLE*DIPOLE
C        SIG1(1) = SIG_NRES(1)*(1.0 + EPS*R_NRES) + SIG_RES(1)
C        SIG1(2) = SQRT( (SIG_NRES(2)*(1.0 + EPS*R_NRES))**2 +
C     >                   SIG_RES(2)**2)
C        SIG2(1) = SIG1(1)*DW2DEP
C        SIG2(2) = SIG1(2)*DW2DEP
C        sig_nr(1) = sig_nres(1)*dw2dep
C        sig_nr(2) = sig_nres(2)*dw2dep
C        sig_r(1) = sig_res(1)*dw2dep
C        sig_r(2) = sig_res(2)*dw2dep
        
C        sige = sig2(1)

         SIGT = SIG_NRES(1)+SIG_RES(1)
c         sigl = r*sigt  
        SIGL = R_NRES*SIG_NRES(1)
         W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
         W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END

      SUBROUTINE H2MODEL_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J, KK
      REAL*8 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     WR, KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, GAM, 
     >     GR, RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER, ERRCHECK

      REAL*8 PI, ALPHA, AM, MPPI, MDELTA, MPI, GAM1WID, GAM2WID, MASS1,
     >     ROPERWID, MASSROPER, MASS2, DELWID, FIT(30), SQRTWDIF,
     >     XR, MQDEP

       REAL*8 RCOEF(27), COEF(22), RERR1(200),rerr2(125), ERR1(200), 
     >        ERR2(53)
    
      LOGICAL FIRST

      FIRST = .TRUE.

c      pi = 3.14159265
c      alpha = 0.00729735
c      am = .9382727
c      MPPI = 1.07783
c      MDELTA = 1.2340 
c      MPI = 0.13957 
c      GAM1WID = 0.0800 
c      GAM2WID = 0.0900
c      MASS1 = 1.5045 
c      ROPERWID = 0.0500 
c      MASSROPER = 1.4000 
c      MASS2 = 1.6850
c      DELWID = 0.1200
c      XR = 0.1800 
c      MQDEP = 3.40 
      
      pi = 3.14159265
      alpha = 0.00729735
      am = .9382727
      MPPI = 1.07783
      MDELTA = 1.229 
      MPI = 0.13957 
      GAM1WID = 0.080
      GAM2WID = 0.090
      MASS1 = 1.5062 
      ROPERWID = 0.0500 
      MASSROPER = 1.4000 
      MASS2 = 1.6810
      DELWID = 0.120
      XR = 0.1800 
      MQDEP = 3.40 

      DATA RCOEF/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >   4.6273E-01,   1.7844E-01,   -1.5416E+02, -1.4891E+02,
     >   2.4102E+02,   2.5823E+00,    7.1004E+00, -8.9771E+00,
     >   1.3744E+00,  -1.2085E+00,    1.1218E-01/

      DATA COEF/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/
 
c      DATA RERR1/
c     >  2.6120E+02,-9.4211E+02, 4.0844E+03, 7.4994E+02,-3.5317E+03,
c     >  3.1703E+03, 2.1238E-01,-6.1568E-01, 4.1479E-01, 1.9720E-02,
c     >  1.2891E-01,-4.1615E+00, 4.7246E+00, 2.8090E-03, 6.1657E-02,
c     >  1.3120E+00,-9.4379E+00, 9.0902E+00,-1.3636E-03, 2.8054E-02,
c     >  9.6123E-02,-1.1465E+03, 3.9099E+03,-3.0097E+03,-1.0604E+01,
c     > -1.2214E+00,-8.3549E-01, 1.3696E+04, 3.9085E+03,-1.5369E+04,
c     >  1.2745E+04, 2.9942E+01, 7.7268E+00, 1.0113E+01,-4.3868E+04,
c     >  1.5709E+05,-3.0207E+03, 1.2809E+04,-1.1075E+04,-2.0442E+01,
c     > -7.5843E+00,-1.0773E+01, 3.2597E+04,-1.2457E+05, 1.0283E+05,
c     > -1.6960E-01, 5.9410E-01,-4.5486E-01,-1.0715E-02,-2.6512E-03,
c     >  1.0153E-03, 6.4074E+00,-1.9189E+01, 1.3612E+01, 9.2486E-03,
c     >  2.7904E-01, 6.3576E+00,-7.8552E+00, 1.5302E-02,-1.1506E-01,
c     > -4.7552E-02,-1.0171E+01,-1.5884E+00, 1.6223E+01,-1.1379E-04,
c     >  4.9212E-01,-2.4354E+00, 1.7921E+01,-1.7223E+01, 4.0778E-03,
c     > -4.5558E-02,-1.8539E-01, 7.9930E+00,-7.1588E+01, 7.1512E+01,
c     > -2.1529E-03, 1.8337E-01, 7.7590E-01, 7.3007E+02,-2.5219E+03,
c     >  1.9547E+03, 6.1102E+00, 1.2970E+00,-1.3084E+00,-9.4932E+03,
c     >  3.1542E+04,-2.3894E+04,-5.9583E+00, 8.1005E-02, 3.6885E-01,
c     >  9.5708E+03,-2.4911E+03, 9.4342E+03,-7.7120E+03,-1.8608E+01,
c     > -1.1065E+00, 6.5015E+00, 3.1755E+04,-1.1529E+05, 9.1964E+04,
c     >  1.8347E+01,-2.5899E+00, 7.1169E-01,-3.2268E+04, 1.1891E+05,
c     >  1.9339E+03,-7.7737E+03, 6.6128E+03, 1.3392E+01,-7.3587E-02,
c     > -4.9353E+00,-2.4223E+04, 9.2768E+04,-7.6712E+04,-1.3210E+01,
c     >  1.2513E+00,-4.5156E+00, 2.4541E+04,-9.5131E+04, 7.8848E+04,
c     >  1.9952E-02,-7.1332E-02, 5.5522E-02, 9.8804E-04, 2.3682E-04,
c     > -7.9762E-05,-6.3638E-01, 1.9492E+00,-1.4036E+00,-9.9312E-04,
c     > -7.8634E-05, 8.2617E-05, 6.8002E-01,-2.1138E+00, 1.5308E+00,
c     >  1.3008E-04,-1.0389E+02, 3.5942E+02,-2.7883E+02,-6.0671E-01,
c     > -1.3016E-01, 1.4621E-01, 1.2841E+03,-4.3361E+03, 3.3132E+03,
c     >  7.0701E-01, 1.2805E-01, 1.3355E-01,-1.4645E+03, 4.9522E+03,
c     > -3.7686E+03,-1.0047E-01, 2.7406E+02, 3.5483E+02,-1.3433E+03,
c     >  1.0978E+03, 1.9033E+00, 5.3726E-02,-8.1621E-01,-4.3612E+03,
c     >  1.6110E+04,-1.2957E+04,-2.2247E+00,-2.1299E-01,-5.8178E-01,
c     >  4.9755E+03,-1.8393E+04, 1.4724E+04, 3.1774E-01,-9.2555E+02,
c     >  3.4086E+03,-2.7508E+02, 1.1025E+03,-9.3653E+02,-1.4100E+00,
c     >  7.3163E-02, 6.6492E-01, 3.3590E+03,-1.3073E+04, 1.0893E+04,
c     >  1.6311E+00, 2.4826E-01, 8.3308E-01,-3.7999E+03, 1.4772E+04,
c     > -1.2252E+04,-2.3255E-01, 7.0167E+02,-2.7162E+03, 2.2434E+03,
c     >  3.0688E+00,-1.0328E+01, 7.8828E+00, 3.6601E-03, 1.3367E-03,
c     > -2.9672E-03,-3.2441E+01, 1.0979E+02,-8.3795E+01,-6.6345E-03/

c      DATA  rerr2/
c     > 3.7074E-02,
c     >-5.7300E-02, 1.5212E-02, 4.5952E-04,
c     > 1.1568E-04,-2.9315E-04,-4.6018E-01,
c     > 9.3624E-01,-4.5908E-01,-6.2914E-05,
c     > 1.1699E-03, 2.0141E-03, 6.9968E-02,
c     >-1.9348E-01, 1.2176E-01, 5.4214E-07,
c     > 1.3845E-04, 2.5311E-03,-2.5396E-03,
c     >-1.2757E-04, 2.4971E-04,-1.2737E-04,
c     > 7.2023E-03,-4.1404E-03, 4.6704E-04,
c     > -4.6388E-03,-5.2545E-03, 4.0159E+01,-1.3481E+02, 1.0186E+02,
c     >  1.1796E-03,-9.1088E+00, 3.0200E+01,-2.2552E+01, 4.3562E-01,
c     > -1.0404E+01, 3.8414E+01,-3.0978E+01,-1.4730E-02, 4.6327E-03,
c     >  1.9716E-02, 1.1236E+02,-4.1952E+02, 3.3862E+02, 2.4150E-02,
c     >  1.1098E-02, 2.0122E-02,-1.3812E+02, 5.1058E+02,-4.0773E+02,
c     > -4.1791E-03, 3.0702E+01,-1.1132E+02, 8.7622E+01,-1.4199E+00,
c     >  5.0230E+00, 8.0171E+00,-3.1384E+01, 2.6350E+01, 1.3147E-02,
c     > -6.1508E-03,-1.6808E-02,-8.7538E+01, 3.4530E+02,-2.8922E+02,
c     > -1.9581E-02,-1.0895E-02,-2.4705E-02, 1.0611E+02,-4.1369E+02,
c     >  3.4296E+02, 3.2847E-03,-2.3191E+01, 8.8502E+01,-7.2288E+01,
c     >  1.0469E+00,-3.8873E+00, 3.1142E+00,
c     > 1.1348E+00,-1.7657E+00, 4.7686E-01,
c     > 1.6653E-02, 4.3488E-04,-7.5168E-03,
c     >-1.6696E+01, 3.4692E+01,-1.7470E+01,
c     >-4.9697E-03, 4.4232E-02, 5.7617E-02,
c     > 5.7800E+00,-1.3886E+01, 7.9819E+00,
c     > 3.4744E-04,-5.4411E-01, 1.2683E+00,
c     >-7.0771E-01, 1.1282E-02,-2.4800E-02,
c     > 1.2909E-02, 1.5171E-01,-6.0417E-01,
c     > 7.7405E-01,-5.8981E-02,-5.8502E-03,
c     > 8.8611E-04, 5.8326E-03, 6.5418E+00,
c     >-1.2978E+01, 6.1069E+00, 1.2462E-03,
c     >-1.8442E-02,-2.7954E-02,-1.8335E+00,
c     > 4.3674E+00,-2.4393E+00,-6.2354E-05,
c     > 1.4746E-01,-3.4127E-01, 1.8285E-01,
c     >-3.0479E-03, 6.8138E-03,-3.4673E-03,
c     >-7.5270E-02, 4.0914E-02/


c      DATA ERR1/
c     >  3.7797E+02,-1.2732E+03, 4.8470E+03, 9.7589E+02,-3.9592E+03,
c     >  3.3447E+03, 1.9629E-01,-4.2402E-01, 1.9757E-01, 3.0613E-02,
c     > -4.0257E-01,-2.0922E+00, 3.0126E+00, 3.8385E-03, 7.3553E-02,
c     >  1.4084E+00,-8.4718E+00, 7.8586E+00,-1.6484E-03, 2.2185E-02,
c     >  7.4896E-02,-1.5627E+03, 5.0106E+03,-3.7125E+03,-1.1701E+01,
c     > -6.9186E-01,-1.4263E+00, 1.5792E+04, 5.0288E+03,-1.7793E+04,
c     >  1.3974E+04, 3.1643E+01, 5.0040E+00, 9.9958E+00,-4.8540E+04,
c     >  1.6247E+05,-3.7498E+03, 1.4066E+04,-1.1461E+04,-2.0806E+01,
c     > -5.0428E+00,-9.7813E+00, 3.5056E+04,-1.2382E+05, 9.7850E+04,
c     > -2.0038E-01, 5.9769E-01,-4.0397E-01,-1.5776E-02,-3.7509E-03,
c     >  5.7496E-04, 7.2218E+00,-2.0335E+01, 1.3722E+01, 1.2562E-02,
c     >  1.4708E+00, 1.8510E+00,-4.1856E+00, 1.9572E-02,-1.3469E-01,
c     > -3.7791E-02,-1.5215E+01, 1.8843E+01,-9.9384E-01, 5.4133E-04,
c     >  5.6775E-01,-2.4158E+00, 1.5245E+01,-1.4180E+01, 5.3668E-03,
c     > -3.5419E-02,-1.4360E-01, 7.8707E+00,-5.7677E+01, 5.5406E+01,
c     > -7.5727E-04, 1.4127E-01, 5.8964E-01, 1.0277E+03,-3.3407E+03,
c     >  2.4943E+03, 6.1372E+00, 2.0731E+00,-1.0628E-01,-1.1445E+04,
c     >  3.6033E+04,-2.6376E+04,-6.4849E+00,-1.5437E+00,-3.1093E+00,
c     >  1.1966E+04,-3.3062E+03, 1.1473E+04,-8.9323E+03,-1.7658E+01,
c     > -3.0298E+00, 2.4862E+00, 3.6140E+04,-1.2237E+05, 9.3797E+04,
c     >  1.8377E+01, 2.4649E-01, 9.5713E+00,-3.7362E+04, 1.2613E+05,
c     >  2.4733E+03,-8.9836E+03, 7.2301E+03, 1.2133E+01, 1.0120E+00,
c     > -2.0972E+00,-2.6581E+04, 9.4364E+04,-7.4804E+04,-1.2397E+01,
c     >  5.8276E-01,-9.1893E+00, 2.7145E+04,-9.6250E+04, 7.6086E+04,
c     >  2.4070E-02,-7.3772E-02, 5.1165E-02, 1.4597E-03, 3.3977E-04,
c     > -2.6275E-05,-7.2542E-01, 2.0676E+00,-1.4052E+00,-1.3577E-03,
c     > -1.4477E-04,-8.5451E-05, 7.4811E-01,-2.1217E+00, 1.4288E+00,
c     >  1.7439E-04,-1.6022E+02, 5.2231E+02,-3.9172E+02,-4.1771E-01,
c     > -2.3133E-01,-1.9119E-02, 1.6931E+03,-5.4146E+03, 4.0099E+03,
c     >  6.5228E-01, 4.5766E-01, 6.7254E-01,-2.0266E+03, 6.3551E+03,
c     > -4.6404E+03,-9.4689E-02, 4.2768E+02, 5.1531E+02,-1.7829E+03,
c     >  1.3890E+03, 1.1798E+00, 3.1335E-01,-2.5902E-01,-5.3955E+03,
c     >  1.8502E+04,-1.4311E+04,-1.8045E+00,-9.6753E-01,-2.0260E+00,
c     >  6.3626E+03,-2.1445E+04, 1.6387E+04, 2.6350E-01,-1.3456E+03,
c     >  4.5055E+03,-3.8598E+02, 1.3911E+03,-1.1170E+03,-7.9328E-01,
c     > -7.6840E-02, 2.5967E-01, 4.0005E+03,-1.4347E+04, 1.1455E+04,
c     >  1.1795E+00, 6.2629E-01, 1.6961E+00,-4.6485E+03, 1.6399E+04,
c     > -1.2954E+04,-1.7187E-01, 9.8638E+02,-3.4363E+03, 2.7002E+03,
c     >  6.0266E+00,-1.9528E+01, 1.4686E+01,-1.7956E-02, 3.3364E-03,
c     >  1.2080E-03,-5.5018E+01, 1.7933E+02,-1.3517E+02, 7.9955E-03/


c       DATA ERR2/
c     > -2.1546E-02,-2.3493E-02, 7.4315E+01,-2.3518E+02, 1.7398E+02,
c     > -6.4429E-04,-1.9950E+01, 6.3147E+01,-4.6881E+01, 1.2816E+00,
c     > -1.9366E+01, 6.5755E+01,-5.0971E+01, 5.7005E-02, 3.3439E-04,
c     >  5.5786E-03, 1.7715E+02,-6.1369E+02, 4.7999E+02,-2.9558E-02,
c     >  5.5461E-02, 7.1075E-02,-2.3560E+02, 7.9154E+02,-6.0792E+02,
c     >  2.7242E-03, 6.3265E+01,-2.0981E+02, 1.6050E+02,-4.0749E+00,
c     >  1.3388E+01, 1.4562E+01,-5.1058E+01, 4.0557E+01,-4.3474E-02,
c     > -4.4868E-03,-6.3041E-03,-1.3274E+02, 4.7814E+02,-3.8441E+02,
c     >  2.5677E-02,-3.8538E-02,-5.8204E-02, 1.7424E+02,-6.0799E+02,
c     >  4.8014E+02,-2.6425E-03,-4.6992E+01, 1.6058E+02,-1.2570E+02,
c     >  3.0554E+00,-1.0258E+01, 7.9929E+00/


! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
!Define error matrix:
c        KK = 0
c        if (goroper) then 
c          DO J = 1,25
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) RERRMAT(I,J) = RERR1(KK)
c              if (kK.le.325.and.kK.gt.200) RERRMAT(I,J)=RERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c       if (.not.goroper) then  
c          DO J = 1,22
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) ERRMAT(I,J) = ERR1(KK)
c              if (kK.le.253.and.kK.gt.200) ERRMAT(I,J)=ERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c        if (goroper) then
c          DO J = 1,25
c              DO I = J+1,25
c                RERRMAT(I,J) = RERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
c        if (.not.goroper) then
c          DO J = 1,22
c              DO I = J+1,22
c                ERRMAT(I,J) = ERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(25) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
        FIT(I + 12) = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 18)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 21)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO
      if (goroper) FIT(26)  = FIT(25)/sqrt(Q2)
      if (goroper) FIT(27)  = FIT(25)/Q2

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,27
c          RSIGERR(J) = 0.0
c          DO I = 1,25
c            RSIGERR(J) = RSIGERR(J) + FIT(J)*FIT(I)*RERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*RERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12.OR.J.EQ.17
     > .OR.J.EQ.18 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + RSIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
c             SIG_RES1(2) = SIG_RES1(2) + RSIGERR(J)
          elseIF(j.ge.25.and.j.le.27) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
c            SIGroper(2) = SIGroper(2) + RSIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + RSIGERR(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
c          SIGERR(J) = 0.0
c          DO I = 1,22
c            SIGERR(J) = SIGERR(J) + FIT(J)*FIT(I)*ERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*ERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + SIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
c            SIG_RES1(2) = SIG_RES1(2) + SIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + SIGERR(J)
          ENDIF
        ENDDO
      endif

! ERRCHECK should agree with SIGTOTER.
C      ERRCHECK = SQRT(ABS(SIG_RES2(2) + SIG_RES1(2) + SIG_NRES(2)))
C      SIGTOTER = SQRT(SIGTOTER)
c      SIG_RES2(2) = SQRT(ABS(SIG_RES2(2)))
c      SIG_RES1(2) = SQRT(ABS(SIG_RES1(2)))
c      SIG_NRES(2) = SQRT(ABS(SIG_NRES(2)))

      RETURN
      END

CCC  Version 050521  -  Author:  M.E. Christy   modified by S. Malace CCC
CCC  Subroutine to get Transvese eD cross sections  CCC 
CCC  from fits to T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


C S. Malace 7/22/06 to get D2 

      SUBROUTINE rescssim(W2,Q2,F1)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval11(50),xval12(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1
      integer i,npts,sf
 
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11872E+01,0.15218E+01,0.14125E+01,0.21081E+01,0.17094E+01,
c    & 0.28094E+01,0.27909E+00,0.13229E+00,0.43753E+00,0.54004E+00,
c    & 0.30149E+00,0.97815E+00,0.491276E+04,0.188911E+04,0.352732E+04,
c    & -.29889E+03,0.253598E+04,0.96323E+03,0.881774E+05,0.344864E+05,
c    & 0.28063E+02,0.633321E+05,0.761637E+05,-.473462E+04,0.90734E+02,
c    & 0.8165041E+06,0.264528E+06,0.808716E+05,0.458931E+05,0.88573E+03,
c    & 0.12218E+02,0.37465E+01,0.90133E+03,0.1007504E+05,0.4736814E+05,
c    & -.463262E+04,0.356829E+03,0.17208E+01,0.13029E+01,0.65986E-01,
c    & 0.6864799E+05,0.375928E+02,0.342095E+01,-.91793E+00,0.23611E-02,
c    & 0.27051E+00,0.136559E+01,0.65689E+00,0.171501E+03,0.26701E+01 /
c
c
c     data xval12/
c
c
C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11996E+01,0.153641E+01,0.142073E+01,0.212016E+01,0.171082E+01,
c    & 0.287999E+01,0.29076E+00,0.192957E+00,0.439143E+00,0.557558E+00,
c    & 0.326104E+00,0.97999E+00,0.45174E+04,0.191460E+04,0.353009E+04,
c    & -.31843E+03,0.600311E+04,0.87710E+03,0.8694928E+05,0.3897957E+05,
c    & 0.28447E+02,0.6324006E+05,0.7590074E+05,-.48051E+04,0.99998E+02,
c    & 0.8355478E+06,0.2756187E+06,0.845187E+05,0.490728E+05,0.9423E+03,
c    & 0.13145E+02,0.36846E+01,0.105307E+04,0.96985E+04,0.472598E+05,
c    & -.46531E+04,0.35792E+03,0.1734E+01,0.13036E+01,0.65573E-01,
c    & 0.908018E+05,0.395258E+02,0.341936E+01,-.90659E+00,0.24026E-02,
c    & 0.27467E+00,0.136335E+01,0.66170E+00,0.16993E+03,0.26701E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

cc Iteration2 follows: out in simonaD2SimFig3.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.119098E+01,0.153249E+01,0.142649E+01,0.213592E+01,0.17119E+01,
c    & 0.297601E+01,0.271804E+00,0.173939E+00,0.394625E+00,0.568884E+00,
c    & 0.296322E+00,0.979981E+00,0.49315E+04,0.26699E+04,0.199866E+04,
c    & -.29483E+03,0.99585E+04,0.177447E+04,0.995236E+05,0.376651E+05,
c    & 0.32501E+02,0.631396E+05,0.772169E+05,-.457579E+04,0.989009E+02,
c    & 0.8603541E+06,0.2483892E+06,0.8244512E+05,0.47721E+05,0.9463E+03,
c    & 0.14093E+02,0.39168E+01,0.12428E+04,0.1017448E+05,0.476793E+05,
c    & -.46163E+04,0.358308E+03,0.17098E+01,0.13262E+01,0.65102E-01,
c    & 0.942561E+05,0.65499E+02,0.30852E+01,-.954945E+00,0.20332E-02,
c    & 0.29676E+00,0.135420E+01,0.68348E+00,0.16509E+03,0.2669E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c Following is for simonaD2SimFit4.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11888E+01,0.15245E+01,0.142628E+01,0.213139E+01,0.17127E+01,
c    & 0.29098E+01,0.30225E+00,0.141739E+00,0.41125E+00,0.550194E+00,
c    & 0.29422E+00,0.97980E+00,0.456537E+04,0.181743E+04,0.25647E+04,
c    & -.36697E+03,0.24858E+04,0.78326E+03,0.952766E+05,0.346543E+05,
c    & 0.30294E+02,0.664647E+05,0.787561E+05,-.43423E+04,0.9999E+02,
c    & 0.9195253E+06,0.3002124E+06,0.871368E+05,0.51180E+05,0.9766E+03,
c    & 0.1320E+02,0.3838E+01,0.10652E+04,0.100941E+05,0.4728109E+05,
c    & -.460501E+04,0.36059E+03,0.17467E+01,0.130749E+01,0.67122E-01,
c    & 0.965392E+05,0.36514E+02,0.336334E+01,-.86248E+00,0.21560E-02,
c    & 0.27810E+00,0.135956E+01,0.657735E+00,0.16344E+03,0.26694E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c following is iteration 3 starting with F2allm. /rc/simonaD2SimFit5.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11871E+01,0.153116E+01,0.143781E+01,0.21428E+01,0.171332E+01,
c    & 0.29792E+01,0.26959E+00,0.171029E+00,0.39371E+00,0.570214E+00,
c    & 0.29644E+00,0.973246E+00,0.503525E+04,0.263216E+04,0.19497E+04,
c    & -.30213E+03,0.950003E+04,0.17909E+04,0.9990789E+05,0.374659E+05,
c    & 0.32625E+02,0.633281E+05,0.775036E+05,-.450712E+04,0.99982E+02,
c    &0.8657796E+06,0.2497579E+06,0.82912E+05,0.481666E+05,0.95763E+03,
c    & 0.14368E+02,0.39312E+01,0.13487E+04,0.101466E+05,0.476961E+05,
c    & -.46158E+04,0.35876E+03,0.17086E+01,0.132586E+01,0.64952E-01,
c    & 0.927868E+05,0.65323E+02,0.30805E+01,-.956083E+00,0.211227E-02,
c    & 0.29716E+00,0.135334E+01,0.68561E+00,0.16425E+03,0.26698E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c following is iteration 3 starting with Bodek. /rc/simonaD2SimFit6.out
C  Simona D2 fit to (simona data) + (ioana-delta)
      data xval11/
     & 0.118651E+01,0.15247E+01,0.14354E+01,0.213435E+01,0.17131E+01,
     & 0.29094E+01,0.30152E+00,0.14436E+00,0.400969E+00,0.55492E+00,
     & 0.29446E+00,0.97992E+00,0.460104E+04,0.17992E+04,0.25114E+04,
     & -.37399E+03,0.24733E+04,0.75449E+03,0.8939804E+05,0.3523602E+05,
     & 0.30091E+02,0.665424E+05,0.790445E+05,-.42655E+04,0.99998E+02,
     & 0.9217486E+06,0.2988304E+06,0.869219E+05,0.512361E+05,0.9782E+03,
     & 0.13294E+02,0.38278E+01,0.10849E+04,0.100891E+05,0.4726618E+05,
     & -.46058E+04,0.36094E+03,0.17463E+01,0.130750E+01,0.670575E-01,
     & 0.969567E+05,0.361745E+02,0.335147E+01,-.86421E+00,0.21919E-02,
     & 0.28124E+00,0.135992E+01,0.65969E+00,0.16259E+03,0.26693E+01 /


      data xval12/


C  Simona D2 fit to ioana data
     & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
     & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
     & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
     & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
     & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
     & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
     & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
     & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
     & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
     & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

       if(q2.lt.4.5.and.w2.lt.2.and.q2.gt.1.7) then

         call resmodsim(1,w2,q2,xval12,sigt)

       else

         call resmodsim(1,w2,q2,xval11,sigt)

      endif  

      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3

      return
      end

CCC  Version 061105  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

C THIS IS FOR D2 from S. MALACE 7/22/07
             
      SUBROUTINE RESMODSIM(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.9382727
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (0.937272 + 0.136)
      wr = wdif/w


c      if(q2.LT.0.3.AND.sf.EQ.1) then
      if(q2.LT.0.15) then
        lowq2 = .true.
        lmax = 2
      endif

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = 0.15
        elseif(l.EQ.2.AND.lowq2) then
          q2 = 0.25
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

       enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.71)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        br_21_1 = 0.5
        br_21_2 = 0.5
        if(sf.EQ.2) then
          br_21_1 = xval(48)
          br_21_2 = 1.- br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +(nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.+nr_coef(i,3)*wdif**(float(3*i-1)/2)
     &       *(1.-xth(1))**3.)/(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

        if(sf.EQ.2) then
c          sig_res = sig_res + sig_4L
          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.1.16.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(0.15-q2)
        else
          slope = sig/0.15
          sig = sig - slope*(0.15-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 

      SUBROUTINE CHRISTYX(W2,Q2,F1,R,sigt,sigl)
      
      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

c     data xval1 / 
c Iteration of July 26, 2006
c    & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
c    & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
c    & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
c    & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
c    & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
c    & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
c    & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
c    & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
c    & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
c    & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


c     data xvalL/
C  Simona - July 28, 2006
c    & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
c    & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
c    & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
c    & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
c    & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
c    & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
c    & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
c    & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
c    & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

c copied from 806 version
c     data xval1 / 
c    & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
c    & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
c    & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
c    & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
c    & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
c    & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
c    & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
c    & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
c    & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
c    & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


c     data xvalL/
c    & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
c    & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
c    & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
c    & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
c    & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
c    & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
c    & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
c    & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
c    & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

c from original version
      data xval1 / 
     &  0.12291E+01,0.15000E+01,0.15117E+01,0.16926E+01,0.16200E+01,
     &  0.14261E+01,0.14100E+00,0.23398E+00,0.81124E-01,0.96627E-01,
     &  0.18885E+00,0.26660E+00,0.85509E+02,0.10378E+02,0.34206E+00,
     &  0.38546E+01,0.73111E+01,0.57629E+00,0.63691E-01,0.24072E+01,
     &  0.19496E+02,0.38962E+01,0.25672E+00,0.29348E+01,0.10138E+02,
     &  0.86788E+01,0.25859E+00,0.29227E+01,0.28345E+02,0.48552E+01,
     &  0.13768E+00,0.49746E+01,0.36645E+01,0.96196E+01,0.29507E+00,
     &  0.23934E+01,0.28438E+03,0.10536E+00,0.12390E+01,0.14593E+00,
     &  -.24992E+03,0.67379E+00,0.23922E+01,-.24997E+00,-.18421E-02,
     &  0.63658E-01,0.19539E+01,0.25343E+00,0.14434E+02,0.10982E+00 /


      data xvalL/
     &  0.12438E+01,0.14948E+01,0.13100E+01,0.17300E+01,0.16700E+01,
     &  0.13950E+01,0.76008E-01,0.53212E-01,0.14228E+00,0.82749E-01,
     &  0.60000E-01,0.42000E+00,0.28298E+02,0.36911E+04,0.98237E+04,
     &  0.11796E-03,0.63566E+01,0.26476E+05,0.36530E+05,0.33431E-03,
     &  0.95795E+05,0.68291E+04,0.35634E+04,0.90878E+05,0.12332E+03,
     &  0.73314E+06,0.26074E+06,0.74001E+02,0.00000E+00,0.92972E+04,
     &  0.87089E+04,0.46171E+05,0.10478E+02,0.73525E+04,0.95318E+04,
     &  0.29686E+00,0.69642E+03,0.13502E+01,0.00000E+00,0.34422E+01,
     &  0.19016E+01,0.20149E+00,0.62176E+03,0.61984E+01,0.16046E+00,
     &  0.16765E+01,0.00000E+00,0.40092E-11,0.71919E+00,0.33674E+00 /

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmodX(1,w2,q2,xval1,sigT)
      call resmodX(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

   
      end

CCC  Version 031606  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMODX(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low,dq2,t,xpr
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2
      common/tst1/sigr,sig_nr

      lowq2 = .false.
      lmax = 1
      q2temp = q2
      dq2 = 0.05

      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w

      br(1,1) = 1.0     !!! single pion branching ratios
      br(2,1) = 0.5
      br(3,1) = 0.65
      br(4,1) = 0.65
      br(5,1) = 0.4
      br(6,1) = 0.65
      br(7,1) = 0.6

      if(sf.EQ.2) then 
        br(6,1) = xval(48)
        br(2,1) = xval(49)
      endif 

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  ? 4th resonance region

      do i=1,7
        x0(i) = 0.165
      enddo
      x0(4) = 0.6

      do i=1,7
        br(i,2) = 1.-br(i,1)
      enddo
    

      if(sf.EQ.1) then
        q2low = 0.00
      else
        q2low = 0.1
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + dq2
        endif

        dip = 1./(1.+q2/0.71)**2             !!!  Dipole parameterization  !!!
        dip2 = dip*dip

        xb = q2/(q2+w2-mp2)
        xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
        xpr = 1./xpr
c        t = log(log((q2+xval(50))/0.330**2)/log(xval(50)/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6              !!!  Read in resonance masses     !!!
          num = num + 1
          mass(i) = xval(i)
        enddo
        do i=1,6              !!!  Read in resonance widths     !!!
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

        if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
          mass(7) = xval(41)
          intwidth(7) = xval(42)
          width(7) = intwidth(7)
        else
          mass(7) = xval(47)
          intwidth(7) = xval(48)
          width(7) = intwidth(7) 
        endif

        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)

          pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)

          if(i.EQ.2) then
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif 

          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

          pgam(i) = intwidth(i)*pgam(i)

          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)

        enddo
 

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo

          if(sf.EQ.1) then

            if(i.eq.6) height(i) = rescoef(i,1)/
     &        (1.+ q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)


             height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.71)**rescoef(i,4)

          else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif

        enddo

CCC    End resonance Q^2 dependence calculations   CCC

     
        do i=1,3               !!!  Non-Res coefficients  !!!
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo


        if(sf.EQ.2) then      !!!  4th resonance region  !!!
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
        else
          height(7) = xval(49)*dip2 
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC

        sig_res = 0.0

        do i=1,7
          sigr(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(i) = sigr(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
          sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
          sig_res = sig_res + sigr(i)   
        enddo


CCC    Finish resonances / start non-res background calculation   CCC

 
        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          enddo

          sig_nr = sig_nr*xpr


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xpr)**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif


        sig = sig_res + sig_nr

          
        if(L.EQ.1) sigtemp = sig  

      enddo
       

CCC   Now extrapolate sig_L linearly to zero for Q^2 less than q2min   CCC

      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/dq2
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif


 1000  format(9f12.5)

      RETURN 
      END 


! Christy fit to proton
      SUBROUTINE CHRISTY806(W2,Q2,F1,R,sigt,sigl)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r,W1p,W2p,nu
      real*8 noverp,fnp_nmc
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval1 / 
     & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
     & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
     & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
     & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
     & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
     & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
     & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
     & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
     & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
     & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


      data xvalL/
     & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
     & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
     & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
     & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
     & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
     & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
     & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
     & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
     & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
     & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

      W1p=0.
      W2p=0.
      R=0.
      F1=0.
      sigl=0.
      sigt=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmod316(1,w2,q2,xval1,sigT)
      call resmod316(2,w2,q2,xvalL,sigL)


      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

      NU = q2 / 2. / mp / xb
      W1p = F1 / MP 
      W2p = W1p /(1.0 + NU*NU/Q2) * (1.0 + R)
      return
   
      end

CCC  Version 031606  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMOD316(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low,dq2,t,xpr
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2
      common/tst1/sigr,sig_nr

      lowq2 = .false.
      lmax = 1
      q2temp = q2
      dq2 = 0.05

      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w

      br(1,1) = 1.0     !!! single pion branching ratios
      br(2,1) = 0.5
      br(3,1) = 0.65
      br(4,1) = 0.65
      br(5,1) = 0.4
      br(6,1) = 0.65
      br(7,1) = 0.6

      if(sf.EQ.2) then 
        br(6,1) = xval(48)
        br(2,1) = xval(49)
      endif 

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  ? 4th resonance region

      do i=1,7
        x0(i) = 0.165
      enddo
      x0(4) = 0.6

      do i=1,7
        br(i,2) = 1.-br(i,1)
      enddo
    

      if(sf.EQ.1) then
        q2low = 0.00
      else
        q2low = 0.1
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + dq2
        endif

        dip = 1./(1.+q2/0.71)**2             !!!  Dipole parameterization  !!!
        dip2 = dip*dip

        xb = q2/(q2+w2-mp2)
        xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
        xpr = 1./xpr
c        t = log(log((q2+xval(50))/0.330**2)/log(xval(50)/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6              !!!  Read in resonance masses     !!!
          num = num + 1
          mass(i) = xval(i)
        enddo
        do i=1,6              !!!  Read in resonance widths     !!!
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

        if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
          mass(7) = xval(41)
          intwidth(7) = xval(42)
          width(7) = intwidth(7)
        else
          mass(7) = xval(47)
          intwidth(7) = xval(48)
          width(7) = intwidth(7) 
        endif

        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)

          pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)

          if(i.EQ.2) then
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif 

          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

          pgam(i) = intwidth(i)*pgam(i)

          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)

        enddo
 

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo

          if(sf.EQ.1) then

            if(i.eq.6) height(i) = rescoef(i,1)/
     &        (1.+ q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)


             height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.71)**rescoef(i,4)

          else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif

        enddo

CCC    End resonance Q^2 dependence calculations   CCC

     
        do i=1,3               !!!  Non-Res coefficients  !!!
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo


        if(sf.EQ.2) then      !!!  4th resonance region  !!!
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
        else
          height(7) = xval(49)*dip2 
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC

        sig_res = 0.0

        do i=1,7
          sigr(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(i) = sigr(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
          sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
          sig_res = sig_res + sigr(i)   
        enddo


CCC    Finish resonances / start non-res background calculation   CCC

 
        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          enddo

          sig_nr = sig_nr*xpr


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xpr)**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif


        sig = sig_res + sig_nr

          
        if(L.EQ.1) sigtemp = sig  

      enddo
       

CCC   Now extrapolate sig_L linearly to zero for Q^2 less than q2min   CCC

      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/dq2
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif


 1000  format(9f12.5)

      RETURN 
      END 


      SUBROUTINE ARENHOVEL_INIT()
      IMPLICIT NONE
      REAL*8 QSQ_MAX,QSQ_GEV,NU,E_NP_FIND,Q2_CM_TEMP,QSQ_CM_MAX,WSQ,NU_T
      INTEGER NQ_,NE_,NMODEL
c      PARAMETER(NQ_=26,NE_=90,NMODEL=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL=2)
      INTEGER NE(0:NQ_,NMODEL)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL),Q2_CM(0:NE_,0:NQ_,NMODEL),
     >  FLT(2,0:NE_,0:NQ_,NMODEL),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE

c      CHARACTER IN_FILE(NMODEL)*100,QSQ_STR*4
      CHARACTER IN_FILE(NMODEL)*100,QSQ_STR*4,qsq_str2*2
      INTEGER IQ,IE,IM
      REAL*8 CHBAR,ALPHA,PI
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PI    = 3.1415927)   
      REAL*8 CONV/25.6819/ ! GeV^2 to fm_2
c      REAL*8 QSQ_V(0:NQ_)/0.,.125, .25, .5, .75, 1.0, 1.25, 1.5, 1.75,
c     >   2., 2.25, 2.5, 3.0, 3.5, 4.0,  4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
c     >   7.5, 8.0, 8.5, 9.0, 9.5, 10./
      REAL*8 QSQ_V(0:NQ_)/0.,.125, .25, .5, .75, 1.0, 1.25, 1.5, 1.75,
     >   2., 2.25, 2.5, 3.0, 3.5, 4.0,  4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
     >   7.5, 8.0, 8.5, 9.0, 9.5, 10., 11., 12., 13., 14., 15., 16.,
     >  17., 18., 19., 20., 22., 24., 26., 28., 30., 32., 34., 36.,
     >  38., 40. /
      REAL*8 E_NP_ADD/20./! extrapolate E_np this many MeV, to F=0
  

! cross section in nb/sr/GeV
      CMOTT  = 0.389E6 * ALPHA**2/4.  ! pb  
      DO IQ=1,NQ_
       QSQ_VAL(IQ) = QSQ_V(IQ) ! put into common block
       IF(QSQ_VAL(IQ).LT.0.2) THEN
        WRITE(QSQ_STR,'(F4.3)')QSQ_VAL(IQ)
       ELSEIF(QSQ_VAL(IQ).GE.10.0) THEN
        WRITE(QSQ_STR,'(F4.1)')QSQ_VAL(IQ)
       ELSE
        WRITE(QSQ_STR,'(F4.2)')QSQ_VAL(IQ)
       ENDIF
       IN_FILE(1)='input/cross-sec-data/arenhovel_v18_'
     >   //QSQ_STR//'_T.dat'
       IN_FILE(2) ='input/cross-sec-data/arenhovel_bonn_' 
     >   //QSQ_STR//'_T.dat'
       if(QSQ_VAL(IQ).GT.10.0) THEN
         WRITE(QSQ_STR2,'(i2.2)')int(QSQ_VAL(IQ))
         IN_FILE(1)='input/cross-sec-data/Bosted_v18_'
     >     //QSQ_STR2//'_T.dat'
         IN_FILE(2)='input/cross-sec-data/Bosted_bonn_'
     >     //QSQ_STR2//'_T.dat'
       endif
       DO IM=1,NMODEL
        OPEN(9,FILE=IN_FILE(IM))
        READ(9,'(/)') ! skip 2 header lines
        DO IE=1,NE_
         if(QSQ_VAL(IQ).GT.10.0) THEN
           E_NP(IE,IQ,IM) = E_NP(IE,25,IM)
           READ(9,'(27X,2E11.3)',END=99) FLT(1,IE,IQ,IM),FLT(2,IE,IQ,IM)
           NE(IQ,IM)=IE +1 ! extra 1 for extrapolation beyond maximum Enp
         else
          READ(9,'(2X,4E11.3)',END=99) E_NP(IE,IQ,IM), 
     >     Q2_CM(IE,IQ,IM),FLT(1,IE,IQ,IM),FLT(2,IE,IQ,IM)
          NE(IQ,IM)=IE +1 ! extra 1 for extrapolation beyond maximum Enp
         endif
        ENDDO ! ie
 99     close(unit=9)
       ENDDO ! im
      ENDDO ! iq
! Put in zero values for QSQ=0
      DO IM=1,NMODEL
       NE(0,IM) = NE(1,IM)     
       DO IE=0,NE(0,IM)
        FLT(1,IE,0,IM) = 0.
        FLT(2,IE,0,IM) = 0.
        E_NP(IE,0,IM) = E_NP(IE,1,IM)
       ENDDO ! ie

       DO IQ=1,NQ_
! Put in zero values for Enp=0
        FLT(1,0,IQ,IM) = 0.
        FLT(2,0,IQ,IM) = 0.
        E_NP(0,IQ,IM) =  0.
       ENDDO
       DO IQ=0,NQ_
! Put in zero values for Enp > Enp Max (beyond quasi-elastic peak)
        E_NP(NE(IQ,IM),IQ,IM) =  E_NP(NE(IQ,IM)-1,IQ,IM) + E_NP_ADD
        E_NP(NE(IQ,IM)+1,IQ,IM) = E_NP(NE(IQ,IM),IQ,IM) +  E_NP_ADD
        FLT(1,NE(IQ,IM),IQ,IM) = 0. 
        FLT(1,NE(IQ,IM)+1,IQ,IM) = 0.
        FLT(2,NE(IQ,IM),IQ,IM) = 0.
        FLT(2,NE(IQ,IM)+1,IQ,IM) = 0.

       ENDDO ! iq
      ENDDO ! im
      QSQ_CM_MAX = QSQ_VAL(NQ_)/CONV
      RETURN
      END

!================================================================
      real*8 function  ARENHOVEL_SIG(E,EP,TH,MODEL,IN_RANGE)
!-----------------------------------------------
! Get cross section from the structure functions
!-----------------------------------------------
      IMPLICIT NONE
      REAL*4 E,EP,TH
      INTEGER MODEL
      INTEGER IN_RANGE
      INTEGER NQ_,NE_,NMODEL_
c      PARAMETER(NQ_=26,NE_=90,NMODEL_=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL_=2)
      INTEGER NE(0:NQ_,NMODEL_)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL_),Q2_CM(0:NE_,0:NQ_,NMODEL_),
     >  FLT(2,0:NE_,0:NQ_,NMODEL_),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE
      REAL*8 QSQ,W1,W2,SIN2T2,CSMOTT,NU
     

      SIN2T2 =(SIN(.017453*TH/2))**2
      QSQ = 4.*E*EP*SIN2T2
      NU = E -EP

      CALL ARENHOVEL_W12(QSQ,NU,W1,W2,MODEL,IN_RANGE)

      CSMOTT =  CMOTT/(E*SIN2T2)**2 
      ARENHOVEL_SIG = CSMOTT *(W2*(1.-SIN2T2) +2.*W1 *SIN2T2)
      RETURN
      END

!=========================================================================
      SUBROUTINE ARENHOVEL_W12(QSQ,NU,W1,W2,MODEL,IN_RANGE)
!---------------------------------------------
! Get W1 and W2 from F_L and F_T
!--------------------------------------------
      IMPLICIT NONE     
      REAL*8 QSQ,NU,W1,W2 
      INTEGER MODEL
      INTEGER IN_RANGE
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/
      REAL*8 F_OUT(2),QSQ_L,Wnp2
      REAL*8 PI2SQ /19.7392/ ! 2*PI**2
      REAL*8 GEV_FM/5.07/ ! 1/(GeV-F 
      REAL*8 ALPHA
      PARAMETER (ALPHA = 7.29735E-03)    

      CALL ARENHOVEL_CALC(QSQ,NU,F_OUT,MODEL,IN_RANGE)
      QSQ_L = QSQ +NU**2  ! 3 momentum squared (lab)
      Wnp2 = MD*(MD+2.*NU) -QSQ
      W2 = (Wnp2/MDSQ*(QSQ/QSQ_L)**2 *F_OUT(1) +
     >      .5*QSQ/QSQ_L *F_OUT(2))/PI2SQ * GEV_FM/ALPHA
      W1 = .5*F_OUT(2)/PI2SQ * GEV_FM/ALPHA

      RETURN
      END

!==================================================================
!--------------------------------------------------------
! Get F_L [F_out(1)] and F_T [F_OUT(2)] by interpolating tables 
!  provided by Arenhovel.
! Lowest Q2 =.125 fm^-2 : Interpolate linearly to zero Q2 and zero F's
! Interpolate linearly from lowest Enp given to Enp=0 and zero F's
! ERROR codes in IN_RANGE
! = 0 OK
! = 1 E_NP < 0    sets output to 0 (correct)
! = 2 Enp > Max in Tables. extrapolate to zero at 20+Enp Max    (OK)
! = 3 QSQ > max in tables
! = 4 QSQ > max in tables AND E_NP > MAX in tables
! <3 is OK to use
! >=3 nead other fit for higher Q^2
!-----------------------------------------------------------------

      SUBROUTINE ARENHOVEL_CALC(QSQ,NU,F_OUT,MODEL,IN_RANGE)
! Get F_L [F_out(1)] and F_T [F_OUT(2)] by interpolating tables 
!-----------------------------------------------------------
      IMPLICIT NONE
      REAL*8 QSQ,NU,F_OUT(2) ! f_out(1) = F_L, f_out(2) = F_t
      INTEGER MODEL
      INTEGER IN_RANGE ! 0= OK, 1=Enp<0, 2= Enp high, 3=Q2 hi, 4=Q2>max && E_NP>MAX 
    
      INTEGER NQ_,NE_,NMODEL_,JQ1,JQ2
c      PARAMETER(NQ_=26,NE_=90,NMODEL_=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL_=2)
      INTEGER NE(0:NQ_,NMODEL_)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL_),Q2_CM(0:NE_,0:NQ_,NMODEL_),
     >  FLT(2,0:NE_,0:NQ_,NMODEL_),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE

      INTEGER IQ,JQ,je(2),KE,LT
      REAL*8 Q2_CM_GEV,E_NP_FIND,Q2_CM_FIND
      REAL*8  denomq,denome(2),w1q,w2q,w1e(2),w2e(2),fac1,fac2,w1ee,w2ee
      REAL*8 CONV/25.6819/ ! GeV^2 to fm_2

      IN_RANGE =0 

      CALL  KIN_TRANSFORM_NORM_ARENHOVEL(QSQ,NU,E_NP_FIND,Q2_CM_GEV)
      IF(E_NP_FIND.LT.0.) THEN
       IN_RANGE =1
       GOTO 999
      ENDIF
      Q2_CM_FIND = Q2_CM_GEV *CONV  ! convert to fm^{-2}
      E_NP_FIND = E_NP_FIND *1000.   ! convert to MeV from GeV

      IF(Q2_CM_FIND.GT.QSQ_VAL(NQ_))THEN ! QSQ > max in tables
       IN_RANGE=3
       IF(E_NP_FIND.GT.E_NP(NE(NQ_,MODEL),NQ_,MODEL)) THEN
        IN_RANGE=4
        GO TO 999
       ENDIF
       IQ = NQ_-1  ! QSQ > max in tables
       W1Q =1.
       W2Q =0.
       DENOMQ=1.
       JQ1=IQ
       JQ2=IQ
       GOTO 1100
      ENDIF
      DO iq=0, NQ_-1
       if((Q2_CM_FIND.GE.QSQ_VAL(iq)).AND.
     >    (Q2_CM_FIND.LE.QSQ_VAL(iq+1)))      THEN
        denomq =ABS(QSQ_VAL(iq+1) -QSQ_VAL(iq))
        w1q = abs(Q2_CM_FIND-QSQ_VAL(iq))
        w2q = abs(QSQ_VAL(iq+1)-Q2_CM_FIND)
        GO TO 1000
       ENDIF
      ENDDO  ! iq
1000  CONTINUE
      JQ1=IQ
      JQ2=IQ+1
      IF((E_NP_FIND.GT.E_NP(NE(IQ,MODEL),IQ,MODEL)).AND.
     >   (E_NP_FIND.GT.E_NP(NE(IQ+1,MODEL),IQ+1,MODEL))) THEN
       IN_RANGE=2 
       GOTO 999
      ENDIF
 1100 CONTINUE
! find E_NP factors for each of the flanking QSQ's
      DO JQ=JQ1,JQ2
       if(jq.gt.nq_ .or.jq.lt.0) then
         write(6,'(''ERROR!'',4i4,4f10.3)') iq,jq,jq1,jq2,
     >     Q2_CM_FIND.GE.QSQ_VAL(nq_),QSQ_VAL(0),E_NP_FIND
       endif
       if(E_NP_FIND.LE.E_NP(0,JQ,MODEL)) THEN !e_np(0,jq,model) =0 always
        w1e(JQ-IQ+1)=0.
        w2e(JQ-IQ+1)=1.
        denome(JQ-IQ+1)=1.
        je(jq-iq+1)=0
       elseif(E_NP_FIND.GE.E_NP(NE(JQ,MODEL),JQ,MODEL)) THEN !Enp > max
        w1e(JQ-IQ+1)=1.
        w2e(JQ-IQ+1)=0.
        denome(JQ-IQ+1)=1.
        je(jq-iq+1)=NE(JQ,MODEL)
       else 
        DO ke=0, NE(JQ,MODEL) -1 
         if((E_NP_FIND.GE.E_NP(ke,jq,MODEL)).AND.
     >     (E_NP_FIND.LE.E_NP(ke+1,jq,MODEL)))      THEN
          denome(JQ-IQ+1) =
     >     ABS(E_NP(ke+1,jq,MODEL) -E_NP(ke,jq,MODEL))
          w1e(JQ-IQ+1) = abs(E_NP_FIND -E_NP(ke,jq,MODEL))
          w2e(JQ-IQ+1) = abs(E_NP(ke+1,jq,MODEL)-E_NP_FIND)
          je(JQ-IQ+1) = ke
          go to 2000
         endif
        ENDDO ! ke
       ENDIF !e_np_find
 2000  continue
      ENDDO  ! jq
      IF(JQ1.EQ.JQ2) THEN  ! high edge of Enp_find range
        JE(2) = JE(1)
        W1E(2) = W1E(1)
        W2E(2) = W2E(1)
        DENOME(2) = DENOME(1)
      ENDIF
      w1ee = (w2q*w1e(1)/denome(1) + w1q*w1e(2)/denome(2))/denomq
      w2ee = (w2q*w2e(1)/denome(1) + w1q*w2e(2)/denome(2))/denomq
          
      DO LT =1,2
       fac1 = 
     >  (w1q*FLT(lt,je(2),iq+1,MODEL)+w2q*FLT(lt,je(1),iq,MODEL))/denomq
       fac2 = 
     >  (w1q*FLT(lt,je(2)+1,iq+1,MODEL)+w2q*FLT(lt,je(1)+1,iq,MODEL))/
     >   denomq
       F_OUT(LT) =w1ee*fac2 + w2ee*fac1
      ENDDO
      RETURN

 999  CONTINUE
      F_OUT(1)=0.
      F_OUT(2)=0.
      RETURN
      END

      SUBROUTINE KIN_TRANSFORM_ARENHOVEL_NORM (E_NP,Q2_CM,WSQ,QSQ,NU)
!------------------------------------------------------
! transform form Arenhovel's coordinates of E_NP, q^2_cm to 
! what I normally use.
!-------------------------------------------------------
      IMPLICIT NONE
      REAL*8 E_NP, Q2_CM, WSQ, QSQ,NU
      REAL*8 C,ARG
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/

      
      C = MDSQ  - (E_NP +MN+MP)**2  -Q2_CM*(E_NP +MN+MP)**2/MDSQ
      ARG = 4.*MDSQ -4.*C
      NU = (-MD2 + SQRT(ARG))/2.
      QSQ = MDSQ+MD2*NU -(E_NP+MN+MP)**2
      WSQ = MPSQ +2.*MP*NU -QSQ
      RETURN
      END
    
      SUBROUTINE KIN_TRANSFORM_NORM_ARENHOVEL(QSQ,NU,E_NP,Q2_CM)   
!------------------------------------------------------
! Everything in GeV
! transform form what I normally use.QSQ,NU
! to Arenhovel's coordinates of E_NP, q^2_cm
!-------------------------------------------------------
      IMPLICIT NONE
      REAL*8 E_NP, Q2_CM, QSQ,NU
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/
      REAL*8 Wnp2

      Wnp2  = MDSQ +MD2*NU -QSQ                 ! eq 5 Arenhovel notes
      E_NP = SQRT(Wnp2) -MP -MN   ! eq. 1. Arenhovel notes
      Q2_CM = MDSQ*(NU**2+QSQ)/Wnp2
      RETURN
      END

      real function ffBe(QSQP) !Mott*beff=Beryllium cross section  
C-----Beryllium Form Factor (done by K.Slifer.  03/26/03). From
C-----Stein et al. PRD 12, 1884 (1975). Eq. A14
C
      implicit none
      real Z,A,Q,XA,XK,XX,XF,W1,W2,QSQP    

      Z   = 4.                  ! Atomic charge
      A   = 9.012               ! Atomic number

      Q   = SQRT(QSQP)/0.1973    ! Convert fom MeV^2 to fm^-1
      XA  = (Z-2.)/3.         ! (A15)
      XK  = 3.*(2.+5.*XA)/2.*(2.+3.*XA)
      XX  = Q*1.07*A**(1./3.)

      XF  = 1. - XA*XX**2. / ( 2.*XK*(2.+3.*XA) )
     &  *EXP(-1.*XX**2./4./XK)

      W2  = (Z*XF)**2.             ! (A14)

      ffBe = (W2)

      RETURN
      END

C     Al. Form Factor. 
C     (from Stein et Al, PRD, 12, 1884)
      REAL FUNCTION ffAl(Q2)
      implicit none
      real Z,A,Q,XB,XC,XF,W2,Q2	


      Z   = 13.                  ! Atomic charge
      A   = 26.98                ! Atomic number


      Q   = SQRT(Q2)/0.1973      ! Convert fom GeV to fm^-1
      XB  = 2.4                  ! fm
      XC  = 1.07*A**(1./3.)      ! fm

      XF  = 1.0/(1.0 + 1.0/6.0 * Q**2 * XC**2)
      XF  = XF*EXP(-1./6.*Q**2*XB**2.)

      W2  = (Z * XF)**2             ! (A14)

      FFal=W2

      RETURN
      END

      real function ffC12(QSQP)
C
C     K. Slifer 10/04/02
C     12C Form Factor
C     See K. C. Stansfield et al. PRC 3, 1448 (1971)
C
      implicit none
      real qsqp,hbarc,zt,at,xalpha,qsqpm,q2fm,xa

      HBARC= 197.327053 ! MEV-FM
      ZT   = 6.         ! ATOMIC CHARGE
      AT   = 12.        ! ATOMIC NUMBER
      XALPHA=(ZT-2.)/3.
      QSQPm=QSQP*1000000. !now in MeV^2
      Q2FM= QSQP/HBARC**2                   ! fm^-2
      if (Q2FM.LE.3.2) THEN
         xa=1.64                            ! fm
      elseif(Q2FM.GE.3.5) then
         xa=1.68                            ! fm
      else
         xa=0                               ! fm
      endif
      xa=xa/HBARC                           ! 1/MeV

      FFC12  = 1.0 - XALPHA/(2.0*(2.+3.*XALPHA)) *QSQPM*xa**2
      FFC12  = FFC12*EXP(-(QSQPM*xa**2)/4.)
      if(Q2FM.GE.3.2.and.Q2FM.LE.3.5) then  ! DIFFRACTION MINIMUM
         FFC12  = 1.D-5
      endif

      ffc12 = ffc12 * zt**2

      RETURN
      end

C
C     Copper Form Factor. 
C     (from Stein et Al, PRD, 12, 1884) 
      REAL FUNCTION FFCu(Q2)
      implicit none
      real Z,A,Q,XB,XC,XF,W2,Q2

      Z   = 29.
      A   = 63.54


      Q   = SQRT(Q2)/0.1973      ! Convert fom GeV to fm^-1
      XB  = 2.4                  ! fm
      XC  = 1.07*A**(1./3.)      ! fm

      XF  = 1.0/(1.0 + 1.0/6.0 * Q**2 * XC**2)
      XF  = XF*EXP(-1./6.*Q**2*XB**2.)

      W2  = (Z * XF)**2             ! (A14)

      FFCu=W2

      RETURN
      END

      real function ffhe4(qsq)
      implicit none
      real a,b,qsq,qf,z,w2
      A=.316
      B=.675
      z = 2.
      QF=sqrt(Qsq)/0.197              !Q MUST BE IN FM**-1
      W2=((1.-(A**2*QF**2)**6)*EXP(-B**2*QF**2))**2
      ffhe4 = z**2 * w2
      return
      end

CCC  Version 051407  -  Author:  M.E. Christy                               CCC
CCC  This fit version includes data from E00-116 (see thesis of S. Malace)  CCC 
CCC  as well preliminary data from E00-002.                                 CCC  
CCC  Subroutine to get Transvese and Longitudinal eP cross sections         CCC 
CCC  from fits cross sections over a range of epsilon.  The subroutine      CCC
CCC  resmod.f is required.  Units are in ub/Sr/Gev.                         CCC


      SUBROUTINE christy0507(W2,Q2,F1,R)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end

CCC  Version 040207  -  Author:  M.E. Christy                       CCC
CCC  This routine returns proton photo-absorbtion cross sections     CCC
CCC  for either transverse or longitudinal photons in units of       CCC
CCC  ub/Sr/Gev.                                                      CCC
CCC  Fit form is empirical.  Interpret physics from it at your       CCC
CCC  own risk.                                                       CCC



             

      subroutine eg4radlen(theta,TBEFOR,TAFTER)
! materials before, after target excluding anything
! between Banjo windows. From Piotr Konczykowski, Mon, 10 Sep 2007 
! theta in degrees, tbefor and tafter are in r.l.
      implicit none

      real He,theta,th
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
      real TBEFOR,TAFTER

      th=abs(theta)

      TBEFOR=0.00167415733

         p1=9.037982e-03
         p2=1.624980e-08
         p3=5.874948e-03
         p4=9.681330e-04
         p5=1.040954e-02
         p6=3.256229e-04
         p7=1.276336e-02
         p8=1.231782e-04
         p9=1.362711e-02
         p10=1.022152e-04
         p11=-3.142974e-06
         p12=7.551468e-08

      if((th.ge.0.).and.(th.lt.3.25)) then
         TAFTER=p1+p2*th
      elseif((th.ge.3.25).and.(th.lt.7.)) then  
         TAFTER=p3+p4*th
      elseif((th.ge.7.).and.(th.lt.11.9)) then  
         TAFTER=p5+p6*th
      elseif((th.ge.11.9).and.(th.lt.16.)) then  
         TAFTER=p7+p8*th
      elseif((th.ge.16.)) then  
         TAFTER=p9+p10*th+p11*th*th+p12*th*th*th
      endif   

      end

C=======================================================================
                                                                        
      SUBROUTINE F1F2IN07(Z, A, QSQ, Wsq, F1, F2)                       
!--------------------------------------------------------------------
! Fit to inelastic cross sections for A(e,e')X
! valid for all W<3 GeV and all Q2<10 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 10/20/2006 P. Bosted
!--------------------------------------------------------------------
      implicit none
      real*8 Z,A,qsq,wsq,w,f1c,x
      real*8 avgn,r,dr,nu,eps,kappa,sigres,flux,siglp,sigtp,F1pp,F1dp
      REAL*8 W1,W2,sigt,rc,w1pb,w2pb,F1,F2,sigl,F1d, F1p,qv
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des,fyuse
      real a4, x4, fitemc, emcfac
      logical goodfit
      INTEGER ISM

! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/

! This is for exp(-xx**2/2.), from teste.f
       real*8 xxp(99)/
     > -3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
     > -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
     > -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
     > -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
     > -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
     > -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
     >  0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
     >  0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
     >  1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
     >  1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
     >  2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000/
! these are 100x bigger for convenience
       real*8 fyp(99)/
     > 0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
     > 0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
     > 0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
     > 0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
     > 1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
     > 2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
     > 2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
     > 1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
     > 0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
     > 0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
     > 0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272/

      integer iz,ia,i
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/

! deuteron fit parameters
       real*8 xvald0(50)/
     >  0.1964E+01, 0.1086E+01, 0.5313E-02, 0.1265E+01, 0.8000E+01,
     >  0.2979E+00, 0.1354E+00, 0.2200E+00, 0.8296E-01, 0.9578E-01,
     >  0.1094E+00, 0.3794E+00, 0.8122E+01, 0.5189E+01, 0.3290E+01,
     >  0.1870E+01, 0.6110E+01,-0.3464E+02, 0.9000E+03, 0.1717E+01,
     >  0.4335E-01, 0.1915E+03, 0.2232E+00, 0.2119E+01, 0.2088E+01,
     > -0.3029E+00, 0.2012E+00, 0.1104E-02, 0.2276E-01,-0.4562E+00,
     >  0.2397E+00, 0.1204E+01, 0.2321E-01, 0.5419E+03, 0.2247E+00,
     >  0.2168E+01, 0.2266E+03, 0.7649E-01, 0.1457E+01, 0.1318E+00,
     > -0.7534E+02, 0.1776E+00, 0.1636E+01, 0.1350E+00,-0.5596E-02,
     >  0.5883E-02, 0.1934E+01, 0.3800E+00, 0.3319E+01, 0.1446E+00/
                                                                        
      IA = int(A)
      avgN = A - Z
      nu = (wsq - pm**2 + qsq) / 2. / pm
      qv = sqrt(nu**2 + qsq)

! Cross section for proton or neutron
      W1 = 0.
      W2 = 0.
      IF(IA .lt. 2 .and. wsq.gt.1.155) THEN
        call CHRISTY507(Wsq,Qsq,F1p,Rc,sigt,sigl)
! If neutron, subtract proton from deuteron. Factor of two to
! convert from per nucleon to per deuteron
        if(Z .lt. 0.5) then
          call resmodd(wsq,qsq,xvald0,F1d)
          F1p = F1d * 2.0 - F1p
        endif
        W1 = F1p / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      ENDIF

! For deuteron
      if(IA .eq. 2) then
c get Fermi-smeared R from Erics proton fit
        call pind(Wsq, Qsq, F1c, Rc, sigt, sigl)
c get fit to F1 in deuteron, per nucleon
        call resd(qsq, wsq, xvald0, F1d)
c convert to W1 per deuteron
        W1 = F1d / PM * 2.0
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      endif

! For nuclei
      IF(IA.gt.2) then
        sigt = 0.
        sigl = 0.
        F1d = 0.
        F1p = 0.
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
       if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
! changed 9/24/07
        if(IA.eq.3) kf=0.100
        if(iA.eq.3) Es=0.010 
! changed 9/24/07
        if(IA.eq.4) kf=0.170
        if(iA.eq.4) Es=0.015 
        if(IA.gt.4) kf=0.165
        if(iA.gt.4) Es=0.015 
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
cc for test
        if(iA.gt.60) Es=0.018

! adjust pf to give right width based on kf
        pf = 0.5 * kf 
! assume this is 2 * pf * qv
        DW2DPF = 2. * qv
        dw2des = 2. * (nu + PM) 
! switched to using 99 bins!
cc        do ism = 1,15
cc          fyuse = fy(ism)
cc          WSQP = WSQ + XX(ISM) * PF * DW2DPF - es * dw2des
        do ism = 1,99
          fyuse = fyp(ism)/100.
          WSQP = WSQ + XXp(ISM) * PF * DW2DPF - es * dw2des
          IF(WSQP.GT. 1.159) THEN
            call CHRISTY507(Wsqp,Qsq,F1pp,Rc,sigtp,siglp)
            call resmodd(wsqp,qsq,xvald0,F1dp)
            F1d = F1d + F1dp * Fyuse
            F1p = F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse
          ENDIF
        ENDDO
        Rc = 0.
        if(sigt .gt. 0.) Rc = sigl / sigt
        W1 = (2. * Z * F1d + (A - 2. * Z) * (2. * F1d - F1p)) / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq) 
      ENDIF

      A4 = A
      x4 = qsq / 2. / pm / nu
      emcfac = fitemc(x4, a4, goodfit)

! this version is per nucleon
      F1 = pm * W1 * emcfac  / A
      F2 = nu * W2 * emcfac  / A

! Add MEC correction

      RETURN                                                            
      END                                                               

      SUBROUTINE christy507(W2,Q2,F1,R,sigt,sigl)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end

      SUBROUTINE RESMOD507(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      INTEGER i,j,l,num,sf


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)


      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       


 1000  format(8f12.5)

      RETURN 
      END 

      subroutine pind(W2,Q2,F1,R,sigt,sigl)
! Calculate proton with Fermi smearing of a deuteron 
      implicit none
      real*8 q2,w2,F1,R,sigt,sigl,am/0.9383/,nu,qv,F1p,Rp,sigtp,siglp
      real*8 amd/1.8756/,w2p,pz
      integer ism
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1=0.
      R=0.
      sigt=0.
      sigl=0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
         w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fyd(ism) / 10.
          sigl = sigl + siglp * fyd(ism) / 10.
          F1   = F1   + F1p   * fyd(ism) / 10.
        endif
       enddo
      else
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
c Need avp2f term to get right behavior x>1! 
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fydf(ism) / 100.
          sigl = sigl + siglp * fydf(ism) / 100.
          F1   = F1   + F1p   * fydf(ism) / 100.
        endif
       enddo
      endif

      if(sigt.ne.0.) R = sigl / sigt
      return
      end
      

      subroutine resd(q2,w2,xval,F1)
! Calculate dueteron F1 by Fermi smearing of proton plus neutron 
! Add MEC term (not smeared)
      implicit none
      real*8 q2,w2,xval(50),F1,am/0.9383/,nu,qv,dw2dpf,w2p,sigp,f1sv
      real*8 sigtst,amd/1.8756/, pz, f1m,f2m
      integer ism,i
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/
      logical usemec/.true./
c      common/mecc/usemec

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1 = 0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
        w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fyd(ism) / 10.
        endif
       enddo
      else
       f1 = 0.
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
! Need avp2f term to get right behavior x>1!
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fydf(ism) / 100. 
        endif
       enddo
      endif
! add MEC term
! took out again, then put back in again
      call mec(1.D0,2.D0,q2,w2,f1m,f2m,xval)
      if(usemec) f1 = f1 + f1m
      return
      end

CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC
CCC *** This is modified to give fit to deuteron for sigt
CCC *** Modified to be for sigt only (former sf=1)
CCC *** Modified to take out lowq2 business.
CCC *** params 1-12 are now hard-wired and used for n/p ratio
ccc *** instead
ccc added code to pre-calculate w-dependent parameters
ccc differences from Eric 507
c a) definition of A**2 (factor 2 M W)
c b) powers of L in width for 2-pion
c c) 0.71 -> 0.91 in dipole formula
             
! changed to use F1 instead of sigt
      SUBROUTINE RESMODD(w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(5000,7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,xpr,alpha,pi,F1
      INTEGER i,j,l,lmax,num,iw
      real*8 noverp,fnp_nmc,x,a,b,sig_mec,brp(7,3)
      real*8 xval0(12)/
c new 5/07 values. Values 1 and 7 will be overridden below.
     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00/
      real*8 xvalold(50),w2sv,q2sv,sigrsv(7),md,w2p,wp,wdifp,xprp,nu
      logical first/.true./
      common/tst2/sigrsv,sig_nr,sig_mec
      real br2(7),br3(7)
      save

      sig = 0.
      if(w2.lt.1.07327**2 .or. w2.gt.25 .or. 
     >  q2.lt.0.0 .or. q2.gt.11.0) then
        write(6,'(1x,''error, q2 or w2 out of range'',2f8.3)') w2,q2
        return
      endif

c do this if fitting masses or widths, else set first true in above
      first=.false.
      if(xvalold(1).ne.xval(1) .or.
     >   xvalold(7).ne.xval(7)) first=.true.
      xvalold(1)=xval(1)
      xvalold(7)=xval(7)

      if(first) then
       mp = 0.9382727
       mpi = 0.135
       mpi2 = mpi*mpi
       meta = 0.547
       mp2 = mp*mp
       pi = 3.141593
       alpha = 1./137.036

! branching ratios
       br(1,1) = 1.0     
       br(2,1) = 0.5
       br(3,1) = 0.65
       br(4,1) = 0.65
       br(5,1) = 0.4
       br(6,1) = 0.65
       br(7,1) = 0.6

! angular momenta
       ang(1) = 1.       !!!  P33(1232)
       ang(2) = 0.       !!!  S11(1535)
       ang(3) = 2.       !!!  D13(1520)
       ang(4) = 3.       !!!  F15(1680)
       ang(5) = 0.       !!!  S15(1650)
       ang(6) = 1.       !!!  P11(1440) roper   
       ang(7) = 3.       !!!  ? 4th resonance region

! x0 parameter
       do i=1,7
c 2006
c         x0(i) = 0.165
c 2007
         x0(i) = 0.215
       enddo
c 2006
c      x0(4) = 0.6
c 2007
       x0(1) = 0.1446

! out branching ratio
       do i=1,7
         br(i,2) = 1.-br(i,1)
       enddo
    
! remember w2
       w2sv = w2

! uses xvals of 1-12, 47, and 48
! move masses, wdiths into local variables
! pyb changed to be fixed
       num = 0
       do i=1,6              
         num = num + 1
         mass(i) = xval0(i)
       enddo
       do i=1,6             
         num = num + 1
         intwidth(i) = xval0(num)
       enddo
! changed to allow delta width, mass to vary
! taken out again since xval(1) used in MEC
c       mass(1) = xval(1)
c       intwidth(1) = xval(7)
c 2006
c       mass(7) = xval(47)
c       intwidth(7) = xval(48)
c 2007
       mass(7) = 1.9341
       intwidth(7) = 0.380

! precalculate w-dependent quantites in 0.1 MeV bins
       do iw=1073,5000
        w = 0.001 * (iw+0.5)
        w2 = w**2
        wdif = w - (mp + mpi)
        wr = wdif/w

! Calculate kinematics needed for threshold Relativistic B-W 
        k = (w2 - mp2) / 2. / mp
        kcm = (w2 - mp2) / 2. / w
        epicm = (W2 + mpi**2 -mp2 ) / 2. / w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 ) / 2. / w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 ) / 2. / w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))
        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))
! Calculate partial widths
          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
          if(i.ne.2) then
c            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
c     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)
c make same as resmod507
            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     >       **(ang(i)+2.)
     >       * W / mass(i)
          else
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif
          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)
          pgam(i) = intwidth(i)*pgam(i)
          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(iw,i) = width(i) * pgam(i) / ((W2 - mass(i)**2.)**2. 
     &            + (mass(i)*width(i))**2.) *
     >            kr(i) / k * kcmr(i) / kcm / intwidth(i)
        enddo ! loop on i
c        write(55,'(i5,f7.2,7f10.4)') iw,w,(sigr(iw,i),i=1,7)
       enddo ! loop on iw

       w2 = w2sv
       first = .false.
       do i=1,50
         xvalold(i) = xval(i)
       enddo
! write table for article
c       open(unit=9,file='nhd.tbl1')
       do i=1,7
         br2(i)=br(i,2)
         br3(i)=0.
         if(i.eq.2) br2(i)=0.
         if(i.eq.2) br3(i)=br(i,2)
         if(i.le.6) then
c          write(9,199) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),(xval(12 + (i-1)*4 + j),j=1,4)
c 199      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & ',f7.1,' & ',f6.3,' & ',f6.3,' \\\\')
         else
c          write(9,198) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),xval(49)
c 198      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & 0.0 & 0.0 & 4.0 \\\\')
         endif
       enddo
c       close(unit=9)
c       open(unit=9,file='nhd.tbl2')
c       write(9,197) (xval(i),i=2,6)
       do i=1,2
c         write(9,197) (xval(36 + (i-1)*4 + j),j=1,4),xval(44+i)
c 197     format(f7.1,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',
c     >     f7.4,' \\\\')
       enddo
c       write(9,'(''xval50='',f10.4)') xval(50)
       close(unit=9)
      endif ! if first
      
! get parameters into local variables
      num = 12
! resonance height coefficients. xvals of 13-36
      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo
      enddo
!  Non-Res coefficients xvals of 37-44
      do i=1,2               
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo

! Begin resonance Q^2 dependence calculations   CCC
! uses xvals 49
      do i=1,6
        height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2) * q2 / (1. + rescoef(i,3) * q2))/
ccc     &          (1. + q2/0.71)**rescoef(i,4)
c make same as resmod507
     &          (1. + q2/0.91)**rescoef(i,4)
      enddo
ccc      dip = 1./(1. + q2 / 0.71)**2  
c make same as resmod507
      dip = 1./(1. + q2 / 0.91)**2  
      dip2 = dip**2
ccc      height(7) = xval(49)*dip2 
c make same as resmod507
      height(7) = xval(49)*dip 
      iw = int(1000.*sqrt(w2))
      sig_res = 0.
      do i=1,7
ccc        sigrsv(i) =  height(i) * sigr(iw,i)
c make same as resmod507 by squaring height
        sigrsv(i) =  height(i)**2 * sigr(iw,i)
        sig_res = sig_res + sigrsv(i) 
      enddo

c make same as resmod507
      sig_res = sig_res * sqrt(w2)

! Begin non-resonant part uses xvals 45, 46, 50
! Depends on both W2 and Q2 so can't easily precalculate
      sig_nr = 0.
c      xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
! to make same as resmod507
      xpr = 1.+(w2-(mp+mpi)**2)/(q2+0.05)
      xpr = 1./xpr
      w = sqrt(w2)
      wdif = w - (mp + mpi)
      do i=1,2  
        sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &   /(q2+nr_coef(i,2))**
     &   (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
      enddo
c new section for sqrt(W)
c make same as resmod507 by turning this off
ccc        sig_nr = sig_nr +(xval(2)*(wdif)**0.5)
ccc     &   /(q2 + xval(3))**
ccc     &   (xval(4) + xval(5) * q2 + xval(6) * q2**2)
        
      sig_nr = sig_nr * xpr
     
! Add third term to try to describe MEC, now using Wdiff in 
! deuteron rather than proton
! ** Taken out 10/17/06 
c     md = 2.*mp
c     nu = (q2 + w2 - mp2) / 2. / mp
c     w2p = md**2 + 2. * md * nu - q2
c     Wp = sqrt(w2p)
c     wdifp = wp - md
c     sig_mec = 0.
c     if(wdifp .gt. 0.) then
c       xprp = 1. + (w2p - (md)**2) / (q2 + xval(50))
c       xprp = 1. / xprp
c       sig_mec = (xval(1) + xval(2)*wdifp**(1/2.) +
c    >     xval(3)*wdifp) /
c    &    (q2 + xval(4))**(xval(5) + xval(6) * q2) *
c    >    xprp
c      endif
c     sig = sig_res + sig_nr + sig_mec

      sig = sig_res + sig_nr 
c      write(6,'(1x,i4,2f7.2,4e10.3)') iw,q2,w2,height(1),
c     >  sigr(iw,1),sig_res,sig_nr

! changed to use F1 instead of sigt
      F1 = sig * (w2-mp2)/8./pi/pi/alpha/0.3894e3
      sig = F1

      RETURN 
      END 

      subroutine mec(z,a,q2,w2,f1,f2,xval)
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,f2,xval(50),am/0.9383/,w,nu

      w = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am

! changed to use max(0.3,q2)
      f1 = xval(1) * exp(-(w - xval(2))**2/xval(3)) /
     >   (1. + max(0.3,q2) / xval(4))**xval(5) * nu ** xval(6)
      
      f2 = 0.
      if(q2.gt.0.) f2 = nu * (f1/am) / (1. + nu**2 / q2) 
      return
      end

      SUBROUTINE F1F2QE07(Z, A, qsq, wsq, F1, F2) 
c
C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      integer IA, izz, izzmin, izp, izznom, izdif

c Look up tables for deuteron case
       real*8 fyd(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN                                      

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
! changed 9/24/07
      if(IA.eq.3) kf=0.120
      if(iA.eq.3) Es=0.010 
! changed 9/24/07
      if(IA.eq.4) kf=0.170
      if(iA.eq.4) Es=0.015 
      if(IA.gt.4) kf=0.165
      if(iA.gt.4) Es=0.015 
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
! for test
      if(iA.gt.60) Es=0.018

! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip)) / kf

! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izz)))**2 - 
c   >      qv**2 + 2. * qv * pz - avp2(izz)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'',3i4,6f7.3)')
     >     izznom,izzmin,izz,qsq,wsq,w2p,dwmin/1.e6,pzmin,pznom
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      return
      end


C=======================================================================
                                                                        
      SUBROUTINE F1F2IN09(Z, A, QSQ, Wsq, F1, F2)                       
!--------------------------------------------------------------------
! Fit to inelastic cross sections for A(e,e')X
! valid for all W<3 GeV and all Q2<10 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 10/20/2006 P. Bosted
!--------------------------------------------------------------------
      implicit none
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 Z,A,qsq,wsq,w,f1c,x
      real*8 avgn,r,dr,nu,eps,kappa,sigres,flux,siglp,sigtp,F1pp,F1dp
      REAL*8 W1,W2,sigt,rc,w1pb,w2pb,F1,F2,sigl,F1d, F1p,qv
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des,fyuse
      real*8 xp,f1cor
      real a4, x4, fitemc, emcfac
      logical goodfit
      INTEGER ISM

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/

! This is for exp(-xx**2/2.), from teste.f
       real*8 xxp(99)/
     > -3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
     > -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
     > -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
     > -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
     > -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
     > -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
     >  0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
     >  0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
     >  1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
     >  1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
     >  2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000/
! these are 100x bigger for convenience
       real*8 fyp(99)/
     > 0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
     > 0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
     > 0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
     > 0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
     > 1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
     > 2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
     > 2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
     > 1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
     > 0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
     > 0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
     > 0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272/

      integer iz,ia,i
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/

! deuteron fit parameters
       real*8 xvald0(50)/
     >  0.1964E+01, 0.1086E+01, 0.5313E-02, 0.1265E+01, 0.8000E+01,
     >  0.2979E+00, 0.1354E+00, 0.2200E+00, 0.8296E-01, 0.9578E-01,
     >  0.1094E+00, 0.3794E+00, 0.8122E+01, 0.5189E+01, 0.3290E+01,
     >  0.1870E+01, 0.6110E+01,-0.3464E+02, 0.9000E+03, 0.1717E+01,
     >  0.4335E-01, 0.1915E+03, 0.2232E+00, 0.2119E+01, 0.2088E+01,
     > -0.3029E+00, 0.2012E+00, 0.1104E-02, 0.2276E-01,-0.4562E+00,
     >  0.2397E+00, 0.1204E+01, 0.2321E-01, 0.5419E+03, 0.2247E+00,
     >  0.2168E+01, 0.2266E+03, 0.7649E-01, 0.1457E+01, 0.1318E+00,
     > -0.7534E+02, 0.1776E+00, 0.1636E+01, 0.1350E+00,-0.5596E-02,
     >  0.5883E-02, 0.1934E+01, 0.3800E+00, 0.3319E+01, 0.1446E+00/

cc     
       real*8 F1M
       logical DEBUG/.TRUE./

      IA = int(A)
      avgN = A - Z
      nu = (wsq - pm**2 + qsq) / 2. / pm
      qv = sqrt(nu**2 + qsq)

      if(Wsq.le.0.0) W = 0.0
      W  = sqrt(Wsq)
      x  = QSQ / (2.0 * pm * nu)
      if(Wsq.le.0.0) x = 0.0

! Cross section for proton or neutron
      W1 = 0.
      W2 = 0.
      IF(IA .lt. 2 .and. wsq.gt.1.155) THEN
        call CHRISTY507(Wsq,Qsq,F1p,Rc,sigt,sigl)
! If neutron, subtract proton from deuteron. Factor of two to
! convert from per nucleon to per deuteron
        if(Z .lt. 0.5) then
          call resmodd(wsq,qsq,xvald0,F1d)
          F1p = F1d * 2.0 - F1p
        endif
        W1 = F1p / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      ENDIF

! For deuteron
      if(IA .eq. 2) then
c get Fermi-smeared R from Erics proton fit
        call pind(Wsq, Qsq, F1c, Rc, sigt, sigl)
c get fit to F1 in deuteron, per nucleon
        call resd(qsq, wsq, xvald0, F1d)
c convert to W1 per deuteron
        W1 = F1d / PM * 2.0
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      endif

! For nuclei
      IF(IA.gt.2) then
        sigt = 0.
        sigl = 0.
        F1d = 0.
        F1p = 0.
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
       if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
! changed 4/09
        if(IA.eq.3) kf=0.115
        if(iA.eq.3) Es=0.001 
! changed 4/09
        if(IA.gt.3) kf=0.19
        if(iA.gt.3) Es=0.017
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
c changed 5/09
        if(iA.gt.7) Es=0.0165

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.230
c

        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
! changed 5/09 
        if(iA.gt.55) Es=0.018 
 

! adjust pf to give right width based on kf
        pf = 0.5 * kf 
! assume this is 2 * pf * qv
        DW2DPF = 2. * qv
        dw2des = 2. * (nu + PM) 
! switched to using 99 bins!
cc        do ism = 1,15
cc          fyuse = fy(ism)
cc          WSQP = WSQ + XX(ISM) * PF * DW2DPF - es * dw2des
        do ism = 1,99
          fyuse = fyp(ism)/100.
          WSQP = WSQ + XXp(ISM) * PF * DW2DPF - es * dw2des
          xp = qsq/(wsqp-pm*pm+qsq)

          IF(WSQP.GT. 1.159) THEN
            call CHRISTY507(Wsqp,Qsq,F1pp,Rc,sigtp,siglp)
            call resmodd(wsqp,qsq,xvald0,F1dp)
c            f1cor = 1.0+xvalc(26)*xp+xvalc(27)*xp**2+xvalc(28)*xp**3+
c     &          xvalc(29)*xp**4+xvalc(30)*xp**5

            f1cor = xvalc(26)*(1.0-xvalc(27)*xp*xp)**xvalc(28)*
     &               (1.0-xvalc(29)*exp(-1.0*xvalc(30)*xp))
           

c            f1cor = f1cor/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq))

c            write(6,*) qsq,xp,f1cor
            F1dp = F1dp*f1cor
            F1pp = F1pp*f1cor
            F1d = F1d + F1dp * Fyuse
            F1p = F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse   !!!  Not ready to use yet
c            call rescsp(Wsqp,Qsq,sigtp,siglp)
c            call rescsn(Wsqp,Qsq,sigtn,sigln)
c            F1d =  F1d+(F1n + F1pp) * Fyuse
c            F1p =  F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse
          ENDIF
        ENDDO
        Rc = 0.
        if(sigt .gt. 0.) Rc = sigl / sigt
        W1 = (2. * Z * F1d + (A - 2. * Z) * (2. * F1d - F1p)) / PM 


CC TEST BELOW if commented out CC


c        W1= W1*(1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)
c        write(6,*) P(13),P(14),P(15),P(16),P(17)
c        W1= W1*(1.0+xvalc(26)*x+xvalc(27)*x**2+xvalc(28)*x**3+
c     &          xvalc(29)*x**4+xvalc(30)*x**5)

cc        if(W .GT. 1.3) 

c         if(W .GT. 0.0) 
c     >       W1=W1*(1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ))**2

        CALL MEC2009( Z , A , qsq , wsq , F1M )

        W1 = W1 + F1M
c        if(Wsq .gt.0.0 ) Rc = Rc * ( 1.0 + P(6) + P(23)*A )
        if(Wsq .gt.0.0 ) Rc = Rc * (1.0 + xvalc(35)*x + 
     &          xvalc(36)*x*x/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq)))
c        Rc = Rc/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq))

        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)

        DEBUG=.FALSE.
        IF( W1 .LE. 0.0 .AND. DEBUG ) THEN 
           write(*,*) 'test  = ', Z,A,W,QSQ,x,F1M,W1
           write(*,*) 'test1 = ', 
     >          (1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ)),
     >        (1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)
        ENDIF

      ENDIF

      A4 = A
      x4 = qsq / 2. / pm / nu
      emcfac = fitemc(x4, a4, goodfit)

      F1 = pm * W1 
      F2 = nu * W2  
c      F1 = pm * W1 * emcfac 
c      F2 = nu * W2 * emcfac 



      RETURN                                                            
      END                                          

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MEC2009(z,a,q2,w2,f1)

! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,am/0.9383/,w,nu
      real*8 a1,b1,c1,t1,dw2,dw22,emc,w2min,damp
      real*4 x4,a4,fitemc
      logical goodfit
      integer i
      real*8 pb(20)/ 
     >     0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >     0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >    -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >     0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >     0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 p18

      real*8 x, f1corr 

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

      f1 = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am
      x  = q2 / (2.0 * am * nu )

      if(a.lt.2.5) return

      p18 = p(18)
! special case for 3He
      if(a.gt.2.5 .and. a.lt.3.5) p18 = 70
! special case for 4He
      if(a.gt.3.5 .and. a.lt.4.5) p18 = 170.
! new values for C, Al, Cu
      if(a.gt.4.5) p18 = 215.
      if(a.gt.20.) p18 = 235.
      if(a.gt.50.) p18 = 230.

c      write(6,*) "here2 ", xvalc(1),xvalc(2),xvalc(3),xvalc(5)

       
       f1corr = P(0)*exp(-((W-P(1))**2)/(P(2)))/ 
     >      ((1.0 + MAX( 0.3 , Q2 ) / P(3) ) ** P(4) )*nu**P(5)
     >      *( 1.0 + P18 * A ** ( 1.0 + P(19) * x ) )

       if(a.eq.12) then

        a1 = q2**2*xvalc(1)*exp(-1.0*q2/xvalc(2))/
     &                                (xvalc(3)+q2)**xvalc(4)

        b1 = xvalc(5)+xvalc(6)*q2

c        c1 = xvalc(7)/(1.+xvalc(8)*sqrt(q2))

        c1 = xvalc(33)+xvalc(34)*q2

c         c1 = 0.290   !!! Test  !!!

        t1 = (w2-b1)**2/c1**2/2.
c        dw2 = w2+q2/xvalc(10)-1.*am*am
c        dw22 = w2+q2/3.-1.*am*am
        dw2 = w2+q2*(1.-1./2.2)-1.0*am*am

        if(dw2.LT.0.0) dw2 = 0.0
c        if(dw22.LT.0.0) dw22 = 0.0
        f1corr = a1*(exp(-1.*t1)*sqrt(dw2))

        f1corr = f1corr+xvalc(38)*exp(-1.*((w2-1.25)/xvalc(39))**2.0)*
     &               q2*exp(-1.*xvalc(40)*q2)


        x4 = x
        a4 = a
        goodfit = .true.
c        emc = fitemc(x4,a4,goodfit)
c        f1corr = f1corr/am/emc
        f1corr = f1corr/am
       endif

       f1 = f1corr

       if(f1 .le.1.0E-9 ) f1=0.0
c       write(*,*) 'vahe1= ', A, W*W, Q2, f1corr

      return
      end

      SUBROUTINE F1F2QE09(Z, A, qsq, wsq, F1, F2)
c
C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      real*8 mcor,ecor
      integer IA, izz, izzmin, izp, izznom, izdif
      
      real*8 a1,a2,a3,b1,b2,b3,b4,b5,gd

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

c Look up tables for deuteron case
       real*8 fyd(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

c     Peter Bosted's correction params

       real*8 pb(20)/ 0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >      0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >     -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >      0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/ 
       real*8 y,R

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

CCCCCCCCCCCCCCCCCCC  New FormFactors  CCCCCCCCCCCCCCCCCCCCCCCCC

c        GMp =2.792782/(1.-0.29940*q+5.65622*qsq-5.57350*q*qsq+      !  with Hall A data  !
c     &         6.31574*qsq**2-1.77642*q*qsq**2+0.30100*qsq**3)

c        GEp = 1./(1.-0.31453*q+7.15389*qsq-14.42166*q*qsq+
c     &         23.33935*qsq**2-14.67906*q*qsq**2+4.03066*qsq**3)

c        GMp = 2.792782     !!!    
c     &       * (1.D0 -1.465D0*tau + 1.260D0*tau**2 + 0.262D0*tau**3)
c     &       / (1.D0 + 9.627D0*tau + 0.0D0*tau**2 + 0.0D0*tau**3
c     &               + 11.179D0*tau**4 + 13.245D0*tau**5)


        a1 = -1.651D0
        a2 = 1.287D0
        a3 = -0.185D0
        b1 = 9.531D0
        b2 = 0.591D0
        b3 = 0.D0
        b4 = 0.D0
        b5 = 4.994D0

        GEp = (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)


        a1 = -2.151D0
        a2 = 4.261D0
        a3 = 0.159D0
        b1 = 8.647D0
        b2 = 0.001D0
        b3 = 5.245D0
        b4 = 82.817D0
        b5 = 14.191D0

        GMp = RMUP
     &      * (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

        GD = (1./(1 + qsq/0.71))**2
        GMn = -1.913148           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &  /(1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)

*        GMn = -1.913148*GD

        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        mcor = 1.0-0.09*qsq/(0.1+qsq)
        ecor = 1.0+0.15*qsq/(0.1+qsq)
        GMn = mcor*GMn
        GEn = ecor*GEn


! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
! changed 4/09
      if(IA.eq.3) kf=0.115
      if(iA.eq.3) Es=0.001 
! changed 4/09
      if(IA.gt.3) kf=0.19
      if(iA.gt.3) Es=0.017 
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
c changed 5/09
      if(iA.gt.7) Es=0.0165
c Test MEC  10/13
c         if(iA.gt.7) Es=0.015
c      if(IA.gt.7) kf = 0.228
c      if(iA.gt.7) Es= 0.015

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.228
c

c
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
! changed 5/09 
      if(iA.gt.55) Es=0.018 


! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
c      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
c     >   (1. + exp(-2.4291 * psip)) / kf

CCMEC - testing below  CCCC

      FY = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip)) / kf


! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izz)))**2 - 
c   >      qv**2 + 2. * qv * pz - avp2(izz)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'',3i4,6f7.3)')
     >     izznom,izzmin,izz,qsq,wsq,w2p,dwmin/1.e6,pzmin,pznom
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0. .and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else
        r = 0.4/qsq
      endif


c apply correction factors
      if(A.gt.2) then
         y = (wsq -amp**2) / qv
c         F1 = F1 * (1. + pb(8) + pb(9) * y +
c     >        pb(10)*y**2 + pb(11)*y**3 + pb(12)*y**4 )
c         R = R * (1. + pb(13))
c         F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)

cc correction to correction Vahe
         if(wsq.gt.0.0) then

c            F1=F1*(1.0+P(7)+P(8)*y+P(9)*y**2 +P(10)*y**3 +P(11)*y**4)
c            R = R * ( 1.0 + P(12) )
            F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)
            if(F1.LT.0.0) F1=0.0

         endif
      endif

      return
      end
      BLOCK DATA CORRECTION
      REAL*8 P(0:23)
      COMMON/PARCORR/P
c     DATA P/
c    c       5.1141e-02,   9.9343e-01,   5.3317e-02,   1.3949e+00, 
c    c       5.9993e+00,  -2.9393e-01,   9.9316e-02,   1.0935e-02, 
c    c       3.7697e-01,   1.3542e+00,  -7.4618e+00,   4.3540e+00, 
c    c      -3.7911e-01,   3.5105e-01,  -1.8903e+00,   4.9139e+00, 
c    c      -5.9923e+00,   2.5021e+00,   1.9943e+01,  -5.5879e-02, 
c    c       3.1914e-01,  -1.8657e-01,   3.2746e+01,   4.9801e-03/
       DATA P/
     c       5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
     c       6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
     c       4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
     c      -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
     c      -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
     c       2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03/
       END


C=======================================================================
                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, F1, F2)
                       
!--------------------------------------------------------------------

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,f1mec,f2mec
      real*8 nu,x,mp,mp2,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at
      real*8 xxp(100),fytot,fytot2,norm

      real a4, x4
      real*8 fitemct,emcfac,emcfacL,xfr,xfl
      logical goodfit
      INTEGER ISM,drn,wfn
      external off_mKP_fit

      real*8 xvalc(40) /
     & 0.42590E+00,0.84169E+01,0.34460E+00,0.68082E+01,0.81354E+00,
     & 0.12912E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97981E+00,0.96986E+00,0.10395E+01,0.99471E+00,0.98433E+00,
     & 0.10000E+01,0.97860E+00,0.10060E+01,0.98912E+00,0.99476E+00,
     & 0.99330E+00,0.99866E+00,0.10000E+01,0.10105E+01,0.10077E+01,
     & 0.10780E+01,0.48919E+00,0.12220E+01,0.13706E+00,0.41167E+02,
     & 0.60174E+01,0.12541E-01,0.17868E+00,0.54704E-01,-.22587E+00,
     & 0.18227E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+0 /

      drn = 5
      xfr = 0.95
      xfl = 1.E-3
      wfn = 2

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)
      a4 = A                                                    
      x4 = x

      if(A.GE.12.) then
        deltae = 0.015  !!! energy shift !!!
        kf = 0.228      !!! fermi momentum  !!!
      elseif(A.GE.27) then
        deltae = 0.017
        kf = 0.238
      elseif(A.GE.56) then
        deltae = 0.023
        kf = 0.241
      endif

      norm = 20.471
      norm = norm*2.0
c      norm = 1.0

      f1p = 0.0D0
      f1n = 0.0D0
      f2p = 0.0D0
      f2n = 0.0D0
      fLp = 0.0D0
      fLn = 0.0D0


c        sigt = 0.0D0
c        sigL = 0.0D0
      fytot = 0.0D0
      fytot2 = 0.0D0


! adjust pf to give right width based on kf
      pf = 0.5 * kf 
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + mp) 

      DO ism = 1,99

CCC   
        xxp(ism) = -3.0+6.0*float(ism-1)/98.0

        fyuse = 1.0D0/norm*exp(-0.5*xxp(ism)*xxp(ism))    !!! Gaussian !!!

CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des

CCC

        fytot = fytot+fyuse

c           write(6,2000) w2,q2,ism,xxp(ism),fyuse, wsqp, fytot


        IF(WSQP.GT. 1.159) THEN
          xt = q2/(q2+Wsqp-mp2)
          delta = off_mKP_fit(xt,wfn,drn)
          if(xt.GT.xfr) delta =  off_mKP_fit(xfr,wfn,drn)
          if(xt.LT.xfl) delta =  off_mKP_fit(xfl,wfn,drn)
          if(q2.LT.0.01) delta = 0.0

          offshell = 1.0D0/(1.0D0-delta)  

          offshell = 1.0D0      !!!  test

          x4 = xt

CCC   Next is medium modification factor2  CCC


          emcfac = xvalc(26)*(1.0-xvalc(27)*xt*xt)**xvalc(28)*
     &               (1.0-xvalc(29)*exp(-1.0*xvalc(30)*xt))

          emcfacL = (1.0 + xvalc(35)*xt + 
     &          xvalc(36)*xt*xt)/(1.0+xvalc(31)*log(1.0+xvalc(32)*q2)) 



          call sf(wsqp,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)
    
          f1pp = f1pp*emcfac*offshell
          f1nn = f1nn*emcfac*offshell
          fLpp = fLpp*emcfac*emcfacL*offshell
          fLnn = fLnn*emcfac*emcfacL*offshell
          f2pp = (2.*xt*f1pp+fLpp)/(1.+4.*xt*xt*mp2/q2)
          f2nn = (2.*xt*f1nn+fLnn)/(1.+4.*xt*xt*mp2/q2)

          F1p = F1p + F1pp * Fyuse
          F1n = F1n + F1nn * Fyuse
          F2p = F2p + F2pp * Fyuse
          F2n = F2n + F2nn * Fyuse      
          FLp = FLp + FLpp * Fyuse
          FLn = FLn + FLnn * Fyuse

         ENDIF

      ENDDO

      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

      call MEC2016(Z,A,w2,q2,f1mec)
      F1 = F1 + f1mec
      f2mec = 2.*x*f1mec/(1.+4.*x*x*mp2/q2)
      F2 = F2 + f2mec

c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MEC2016(z,a,w2,q2,f1corr)

! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2,dw22,emc,w2min,damp
      integer i
      real*8 x, f1corr

      real*8 xvalc(40) / 
     & 0.42590E+00,0.84169E+01,0.34460E+00,0.68082E+01,0.81354E+00,
     & 0.12912E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97981E+00,0.96986E+00,0.10395E+01,0.99471E+00,0.98433E+00,
     & 0.10000E+01,0.97860E+00,0.10060E+01,0.98912E+00,0.99476E+00,
     & 0.99330E+00,0.99866E+00,0.10000E+01,0.10105E+01,0.10077E+01,
     & 0.10780E+01,0.48919E+00,0.12220E+01,0.13706E+00,0.41167E+02,
     & 0.60174E+01,0.12541E-01,0.17868E+00,0.54704E-01,-.22587E+00,
     & 0.18227E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+0 /

      mp2 = mp*mp
      f1corr = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(a.lt.2.5) return

      a1 = A*q2**2*xvalc(1)*exp(-1.0*q2/xvalc(2))/
     &                                (xvalc(3)+q2)**xvalc(4)

      b1 = xvalc(5)+xvalc(6)*q2

c        c1 = xvalc(7)/(1.+xvalc(8)*sqrt(q2))

      c1 = xvalc(33)+xvalc(34)*q2

c         c1 = 0.290   !!! Test  !!!


      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1corr = a1*(exp(-1.*t1)*sqrt(dw2))

      f1corr = f1corr+xvalc(38)*exp(-1.*((w2-1.25)/xvalc(39))**2.0)*
     &               q2*exp(-1.*xvalc(40)*q2)


c       f1corr = f1corr

       if(f1corr.LE.1.0E-9 ) f1corr=0.0


      return
      end


      SUBROUTINE F1F2QE16(Z, A, qsq, wsq,F1, F2)

C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
c      REAL*8 P(0:23)
c      COMMON/PARCORR/P
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      real*8 mcor,ecor
      integer IA, izz, izzmin, izp, izznom, izdif
      
      real*8 a1,a2,a3,b1,b2,b3,b4,b5,gd

      real*8 xvalc(40) / 
     & 0.42590E+00,0.84169E+01,0.34460E+00,0.68082E+01,0.81354E+00,
     & 0.12912E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97981E+00,0.96986E+00,0.10395E+01,0.99471E+00,0.98433E+00,
     & 0.10000E+01,0.97860E+00,0.10060E+01,0.98912E+00,0.99476E+00,
     & 0.99330E+00,0.99866E+00,0.10000E+01,0.10105E+01,0.10077E+01,
     & 0.10780E+01,0.48919E+00,0.12220E+01,0.13706E+00,0.41167E+02,
     & 0.60174E+01,0.12541E-01,0.17868E+00,0.54704E-01,-.22587E+00,
     & 0.18227E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+0 /

c Look up tables for deuteron case
       real*8 fyd(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

c     Peter Bosted's correction params

       real*8 pb(20)/ 0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >      0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >     -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >      0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/ 
       real*8 y,R

c       write(6,*) A,Z

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

CCCCCCCCCCCCCCCCCCC  New FormFactors  CCCCCCCCCCCCCCCCCCCCCCCCC

c        GMp =2.792782/(1.-0.29940*q+5.65622*qsq-5.57350*q*qsq+      !  with Hall A data  !
c     &         6.31574*qsq**2-1.77642*q*qsq**2+0.30100*qsq**3)

c        GEp = 1./(1.-0.31453*q+7.15389*qsq-14.42166*q*qsq+
c     &         23.33935*qsq**2-14.67906*q*qsq**2+4.03066*qsq**3)

c        GMp = 2.792782     !!!    
c     &       * (1.D0 -1.465D0*tau + 1.260D0*tau**2 + 0.262D0*tau**3)
c     &       / (1.D0 + 9.627D0*tau + 0.0D0*tau**2 + 0.0D0*tau**3
c     &               + 11.179D0*tau**4 + 13.245D0*tau**5)


        a1 = -1.651D0
        a2 = 1.287D0
        a3 = -0.185D0
        b1 = 9.531D0
        b2 = 0.591D0
        b3 = 0.D0
        b4 = 0.D0
        b5 = 4.994D0

        GEp = (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)


        a1 = -2.151D0
        a2 = 4.261D0
        a3 = 0.159D0
        b1 = 8.647D0
        b2 = 0.001D0
        b3 = 5.245D0
        b4 = 82.817D0
        b5 = 14.191D0

        GMp = RMUP
     &      * (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

        GD = (1./(1 + qsq/0.71))**2
        GMn = -1.913148           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &  /(1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)

*        GMn = -1.913148*GD

        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        mcor = 1.0-0.09*qsq/(0.1+qsq)
        ecor = 1.0+0.15*qsq/(0.1+qsq)
        GMn = mcor*GMn
        GEn = ecor*GEn


! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
! changed 4/09
      if(IA.eq.3) kf=0.115
      if(iA.eq.3) Es=0.001 
! changed 4/09
      if(IA.gt.3) kf=0.19
      if(iA.gt.3) Es=0.017 
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
c changed 5/09
      if(iA.gt.7) Es=0.0165
c Test MEC  10/13
c         if(iA.gt.7) Es=0.015
c      if(IA.gt.7) kf = 0.228
c      if(iA.gt.7) Es= 0.015

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.228
c

c
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
! changed 5/09 
      if(iA.gt.55) Es=0.018 


! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip)) / kf

CCMEC - testing below  CCCC

      FY = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip)) / kf


! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izz)))**2 - 
c   >      qv**2 + 2. * qv * pz - avp2(izz)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'',3i4,6f7.3)')
     >     izznom,izzmin,izz,qsq,wsq,w2p,dwmin/1.e6,pzmin,pznom
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0. .and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else
        r = 0.4/qsq
      endif


c Make sure structure functions are positive
      if(A.gt.2) then
         if(wsq.gt.0.0) then
          F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)
          if(F1.LT.0.0) F1=0.0
         endif
      endif

c      write(6,*) "f1f2qe16:  ", wsq,qsq,f1,f2

      return
      end

      
C=======================================================================
C ***********************************************************************
        FUNCTION off_mKP_fit (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]), constrained to vanish at x=0
C
C  Defined such that F2d = F2d(conv) + del^off F2d 
C       with off_mKP_fit = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Fit by Eric Christy (6/14/12).
C  Fortran code by W. Melnitchouk (6/17/12).
C
C  wfn: 1 (AV18)
C	2 (CD-Bonn)
C	3 (WJC-1)
C	4 (WJC-2)
C
C  dRN: 0 (0.0%)	[% change in nucleon radius in the deuteron]
C	1 (0.3%)
C	2 (0.6%)
C	3 (0.9%)
C	4 (1.2%)
C	5 (1.5%)
C	6 (1.8%)
C	7 (2.1%)
C	8 (2.4%)
C	9 (2.7%)
C      10 (3.0%)
C
C ***********************************************************************
	IMPLICIT NONE
	REAL*8	off_mKP_fit, x
	INTEGER	wfn, dRN
	REAL*8  p(0:9)

	off_mKP_fit = 0.D0
	IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
	IF (dRN.LT.0 .OR. dRN.GT.10) RETURN
! .......................................................................
	IF (wfn.EQ.1) THEN		! AV18
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02628D0
	    p(1) = 1.20029D0
	    p(2) = 7.49503D0
	    p(3) = 2.01901D0
	    p(4) = 0.00789D0
	    p(5) = 0.46739D0
	    p(6) = 0.73242D0
	    p(7) = 0.00328D0
	    p(8) = 0.87228D0
	    p(9) = 0.06400D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = 0.03638D0
	    p(1) = 0.38307D0
	    p(2) = 8.01156D0
	    p(3) = 2.30992D0
	    p(4) = 0.09027D0
	    p(5) = 0.69521D0
	    p(6) = 0.75973D0
	    p(7) = -0.05098D0
	    p(8) = 1.18963D0
	    p(9) = -0.19192D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02260D0
	    p(1) = 1.45377D0
	    p(2) = 0.50628D0
	    p(3) = 13.92200D0
	    p(4) = 0.03558D0
	    p(5) = 0.75147D0
	    p(6) = 0.86335D0
	    p(7) = -0.01383D0
	    p(8) = 1.04749D0
	    p(9) = 0.42099D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06410D0
	    p(1) = 1.18883D0
	    p(2) = 6.96799D0
	    p(3) = 8.87113D0
	    p(4) = 0.02603D0
	    p(5) = 0.70504D0
	    p(6) = 1.44139D0
	    p(7) = 0.00004D0
	    p(8) = -1.14305D0
	    p(9) = 0.73785D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06237D0
	    p(1) = 2.03192D0
	    p(2) = 4.01755D0
	    p(3) = 6.83741D0
	    p(4) = 0.04701D0
	    p(5) = -0.00457D0
	    p(6) = 1.30967D0
	    p(7) = -0.00996D0
	    p(8) = -0.42418D0
	    p(9) = 0.27524D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.06759D0
	    p(1) = 1.95103D0
	    p(2) = 3.54215D0
	    p(3) = 11.77533D0
	    p(4) = 0.09269D0
	    p(5) = 0.56534D0
	    p(6) = 0.98398D0
	    p(7) = -0.03031D0
	    p(8) = 3.26913D0
	    p(9) = -0.45923D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.07007D0
	    p(1) = 2.30938D0
	    p(2) = 4.94226D0
	    p(3) = 8.95701D0
	    p(4) = 0.06933D0
	    p(5) = 0.07145D0
	    p(6) = 1.94887D0
	    p(7) = -0.01210D0
	    p(8) = 5.92311D0
	    p(9) = 0.14312D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11965D0
	    p(1) = 2.06149D0
	    p(2) = 5.38881D0
	    p(3) = 12.08265D0
	    p(4) = 0.19668D0
	    p(5) = 0.61820D0
	    p(6) = 0.80489D0
	    p(7) = -0.08735D0
	    p(8) = 3.74802D0
	    p(9) = -0.70773D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.14735D0
	    p(1) = 2.27109D0
	    p(2) = 8.23092D0
	    p(3) = 7.31581D0
	    p(4) = 0.11953D0
	    p(5) = 0.67459D0
	    p(6) = 1.59118D0
	    p(7) = -0.02700D0
	    p(8) = 4.52840D0
	    p(9) = -1.77765D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.27194D0
	    p(1) = 2.01340D0
	    p(2) = 10.71380D0
	    p(3) = 8.84886D0
	    p(4) = 0.09345D0
	    p(5) = 0.49802D0
	    p(6) = 1.28523D0
	    p(7) = -0.00474D0
	    p(8) = 0.58703D0
	    p(9) = 0.88354D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.69848D0
	    p(1) = 1.48173D0
	    p(2) = 17.44991D0
	    p(3) = 12.73730D0
	    p(4) = 0.13118D0
	    p(5) = 0.34598D0
	    p(6) = 1.65884D0
	    p(7) = -0.02215D0
	    p(8) = 1.21306D0
	    p(9) = 0.96399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.2) THEN		! CD-Bonn
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02820D0
	    p(1) = 0.85879D0
	    p(2) = 9.48856D0
	    p(3) = 2.18885D0
	    p(4) = 0.00070D0
	    p(5) = -5.61817D0
	    p(6) = 14.80512D0
	    p(7) = 0.00348D0
	    p(8) = -1.30292D0
	    p(9) = -0.73075D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.02996D0
	    p(1) = 0.35717D0
	    p(2) = 6.53843D0
	    p(3) = 3.88389D0
	    p(4) = 0.00758D0
	    p(5) = -16.50399D0
	    p(6) = 77.60083D0
	    p(7) = 0.00320D0
	    p(8) = 0.42334D0
	    p(9) = 0.28545D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03261D0
	    p(1) = 0.91185D0
	    p(2) = 8.49348D0
	    p(3) = 10.19681D0
	    p(4) = 0.01598D0
	    p(5) = 0.83748D0
	    p(6) = 1.55960D0
	    p(7) = 0.00085D0
	    p(8) = -0.63447D0
	    p(9) = 0.65632D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.03034D0
	    p(1) = 1.58677D0
	    p(2) = 3.21753D0
	    p(3) = 11.66572D0
	    p(4) = 0.04999D0
	    p(5) = 0.56688D0
	    p(6) = 0.94941D0
	    p(7) = -0.01453D0
	    p(8) = -0.89157D0
	    p(9) = 0.27160D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.04831D0
	    p(1) = 1.75241D0
	    p(2) = 4.74662D0
	    p(3) = 8.29052D0
	    p(4) = 0.04730D0
	    p(5) = 0.33550D0
	    p(6) = 1.18790D0
	    p(7) = -0.00678D0
	    p(8) = -0.42800D0
	    p(9) = 0.36573D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.09019D0
	    p(1) = 1.22091D0
	    p(2) = 1.30114D0
	    p(3) = 17.58701D0
	    p(4) = 0.08312D0
	    p(5) = 0.66902D0
	    p(6) = 0.60767D0
	    p(7) = -0.02035D0
	    p(8) = 0.95978D0
	    p(9) = 1.11322D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.17060D0
	    p(1) = 1.42115D0
	    p(2) = 7.24672D0
	    p(3) = 5.80680D0
	    p(4) = 0.09200D0
	    p(5) = 0.43367D0
	    p(6) = 1.56378D0
	    p(7) = -0.02338D0
	    p(8) = 0.44968D0
	    p(9) = 0.29678D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11026D0
	    p(1) = 1.85213D0
	    p(2) = 6.74413D0
	    p(3) = 7.74362D0
	    p(4) = 0.08467D0
	    p(5) = 0.24708D0
	    p(6) = 1.12274D0
	    p(7) = -0.01505D0
	    p(8) = 0.44209D0
	    p(9) = 0.36126D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.15291D0
	    p(1) = 1.83333D0
	    p(2) = 7.76495D0
	    p(3) = 7.04783D0
	    p(4) = 0.09206D0
	    p(5) = 0.08655D0
	    p(6) = 1.27460D0
	    p(7) = -0.01659D0
	    p(8) = 0.45536D0
	    p(9) = 0.29407D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.24143D0
	    p(1) = 1.50401D0
	    p(2) = 9.33393D0
	    p(3) = 11.62779D0
	    p(4) = 0.09454D0
	    p(5) = 0.36361D0
	    p(6) = 0.82058D0
	    p(7) = -0.00802D0
	    p(8) = 0.34851D0
	    p(9) = 0.50844D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.22196D0
	    p(1) = 1.87228D0
	    p(2) = 10.18898D0
	    p(3) = 9.21038D0
	    p(4) = 0.11850D0
	    p(5) = 0.34360D0
	    p(6) = 1.28278D0
	    p(7) = -0.01754D0
	    p(8) = 0.54540D0
	    p(9) = 0.53457D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.3) THEN		! WJC-1
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.02322D0
	    p(1) = 0.11213D0
	    p(2) = 3.71079D0
	    p(3) = 5.51496D0
	    p(4) = 0.00877D0
	    p(5) = 0.84639D0
	    p(6) = 0.66227D0
	    p(7) = -0.00621D0
	    p(8) = -0.39896D0
	    p(9) = 0.32012D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.00058D0
	    p(1) = 2.33827D0
	    p(2) = 2.35664D0
	    p(3) = 36.75823D0
	    p(4) = -0.00752D0
	    p(5) = 0.05286D0
	    p(6) = 1.27262D0
	    p(7) = 0.01269D0
	    p(8) = 1.72720D0
	    p(9) = 0.20652D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03373D0
	    p(1) = 0.93858D0
	    p(2) = 0.15704D0
	    p(3) = 10.71630D0
	    p(4) = -0.00235D0
	    p(5) = -0.11937D0
	    p(6) = 0.74925D0
	    p(7) = 0.00452D0
	    p(8) = 2.96830D0
	    p(9) = -2.89070D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.08982D0
	    p(1) = 0.73060D0
	    p(2) = -0.16543D0
	    p(3) = 12.37035D0
	    p(4) = 0.04407D0
	    p(5) = 0.47361D0
	    p(6) = 0.74570D0
	    p(7) = -0.00933D0
	    p(8) = 0.53186D0
	    p(9) = 0.26943D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.11990D0
	    p(1) = 1.19824D0
	    p(2) = 3.06386D0
	    p(3) = 8.55017D0
	    p(4) = 0.05815D0
	    p(5) = 0.06123D0
	    p(6) = 1.45024D0
	    p(7) = -0.01414D0
	    p(8) = 0.48172D0
	    p(9) = 0.25171D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.15292D0
	    p(1) = 1.01991D0
	    p(2) = 1.20661D0
	    p(3) = 13.31860D0
	    p(4) = 0.02571D0
	    p(5) = -1.56438D0
	    p(6) = 2.69042D0
	    p(7) = -0.00000D0
	    p(8) = 0.29759D0
	    p(9) = 0.97967D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.35935D0
	    p(1) = 0.44637D0
	    p(2) = -0.25510D0
	    p(3) = 16.70057D0
	    p(4) = 0.10634D0
	    p(5) = 0.61659D0
	    p(6) = 0.58524D0
	    p(7) = -0.03335D0
	    p(8) = 0.93904D0
	    p(9) = 0.89819D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.97384D0
	    p(1) = -0.24934D0
	    p(2) = -0.61349D0
	    p(3) = 18.43254D0
	    p(4) = 0.18772D0
	    p(5) = 0.49599D0
	    p(6) = 0.61366D0
	    p(7) = -0.08116D0
	    p(8) = 0.87175D0
	    p(9) = 0.24026D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.21641D0
	    p(1) = 1.74710D0
	    p(2) = 5.19387D0
	    p(3) = 10.61285D0
	    p(4) = 0.06655D0
	    p(5) = 0.01300D0
	    p(6) = 0.94503D0
	    p(7) = -0.00642D0
	    p(8) = 0.48859D0
	    p(9) = 0.16331D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.32283D0
	    p(1) = 1.71708D0
	    p(2) = 7.51556D0
	    p(3) = 9.68202D0
	    p(4) = 0.09871D0
	    p(5) = 0.18788D0
	    p(6) = 0.80490D0
	    p(7) = -0.01673D0
	    p(8) = 0.48879D0
	    p(9) = 0.21016D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.32064D0
	    p(1) = 2.07514D0
	    p(2) = 9.34847D0
	    p(3) = 8.17225D0
	    p(4) = 0.10772D0
	    p(5) = 0.50272D0
	    p(6) = 1.30663D0
	    p(7) = -0.01215D0
	    p(8) = 0.59432D0
	    p(9) = 0.65399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.4) THEN		! WJC-2
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.03490D0
	    p(1) = 0.78902D0
	    p(2) = -0.25256D0
	    p(3) = 7.98679D0
	    p(4) = 0.00913D0
	    p(5) = 0.74835D0
	    p(6) = 0.60145D0
	    p(7) = -0.00464D0
	    p(8) = 0.41358D0
	    p(9) = 0.22524D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.01119D0
	    p(1) = 0.50514D0
	    p(2) = 19.35710D0
	    p(3) = 3.32395D0
	    p(4) = 0.00670D0
	    p(5) = 1.38279D0
	    p(6) = 1.24216D0
	    p(7) = 0.00049D0
	    p(8) = 0.38623D0
	    p(9) = 0.23497D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02653D0
	    p(1) = 1.27315D0
	    p(2) = -0.53410D0
	    p(3) = 14.08029D0
	    p(4) = 0.01474D0
	    p(5) = 1.82129D0
	    p(6) = 1.99455D0
	    p(7) = -0.00090D0
	    p(8) = 3.96583D0
	    p(9) = 4.61316D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06301D0
	    p(1) = 1.10373D0
	    p(2) = -0.26356D0
	    p(3) = 15.04038D0
	    p(4) = 0.02428D0
	    p(5) = -0.15349D0
	    p(6) = 3.03168D0
	    p(7) = 0.00127D0
	    p(8) = -0.73818D0
	    p(9) = 0.07474D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06150D0
	    p(1) = 2.15792D0
	    p(2) = 2.18241D0
	    p(3) = 9.84713D0
	    p(4) = 0.03608D0
	    p(5) = -0.13604D0
	    p(6) = 1.12241D0
	    p(7) = -0.00695D0
	    p(8) = -0.35646D0
	    p(9) = 0.31793D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.07179D0
	    p(1) = 1.97917D0
	    p(2) = 3.47662D0
	    p(3) = 10.00224D0
	    p(4) = 0.04587D0
	    p(5) = 0.06416D0
	    p(6) = 1.10677D0
	    p(7) = -0.00391D0
	    p(8) = -0.42677D0
	    p(9) = 0.26619D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.09883D0
	    p(1) = 1.96788D0
	    p(2) = 5.19182D0
	    p(3) = 8.82173D0
	    p(4) = 0.06468D0
	    p(5) = 0.11297D0
	    p(6) = 1.63850D0
	    p(7) = -0.00872D0
	    p(8) = 0.52753D0
	    p(9) = 0.41794D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.14258D0
	    p(1) = 2.00822D0
	    p(2) = 6.23508D0
	    p(3) = 7.81846D0
	    p(4) = 0.07064D0
	    p(5) = -0.05869D0
	    p(6) = 1.24848D0
	    p(7) = -0.01160D0
	    p(8) = 0.48932D0
	    p(9) = 0.22001D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.16184D0
	    p(1) = 2.16963D0
	    p(2) = 7.62378D0
	    p(3) = 7.33369D0
	    p(4) = 0.09197D0
	    p(5) = 0.15692D0
	    p(6) = 1.80734D0
	    p(7) = -0.01561D0
	    p(8) = 0.53224D0
	    p(9) = 0.39357D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.20205D0
	    p(1) = 2.28733D0
	    p(2) = 9.10375D0
	    p(3) = 7.24877D0
	    p(4) = 0.08325D0
	    p(5) = 0.36941D0
	    p(6) = 2.39131D0
	    p(7) = -0.00057D0
	    p(8) = 0.41640D0
	    p(9) = 0.90531D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.95664D0
	    p(1) = 1.11409D0
	    p(2) = 19.00631D0
	    p(3) = 15.97282D0
	    p(4) = 0.15616D0
	    p(5) = 0.40229D0
	    p(6) = 0.85878D0
	    p(7) = -0.03123D0
	    p(8) = 6.75437D0
	    p(9) = -3.83159D0
	  ENDIF

	ENDIF

	off_mKP_fit = -( p(0) * x**p(3) * DEXP(p(1) * x**p(2))
     &                + p(4) * x*DEXP(((x-p(5))/p(6))**2)
     &                + x**0.5D0 * p(7) * DEXP(((x-p(9))/p(8))**2) )

	RETURN
	END

C ***********************************************************************
CCCC   Converts reduced cross sections to structure functions   CCCCC
CCCC   for protons and neutrons                                 CCCCC

      SUBROUTINE SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2
      Integer i
 

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp(w2,q2,sigTp,sigLp)
      call rescsn(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)

c      fLp = 2*q2/abs(w2-mp2+q2)*sigLp/0.3894e3/pi2/alpha/8.0*
c     &        abs(w2-mp2) 
c      fLn = 2*q2/abs(w2-mp2+q2)*sigLn/0.3894e3/pi2/alpha/8.0*
c     &        abs(w2-mp2) 
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)


      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end

      
      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis

      data xvalp / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)
c      call disp(w2,q2,sigTdis,sigLdis)
      if(w2.GT.8.0) then 
c        write(6,*) "resmod: ",w2,q2,sigTp,sigLp
        call disp(w2,q2,sigTdis,sigLdis)
c        write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
        if(w2.LE.10.0) then
         sigTp = (10.0-w2)*sigTp+(w2-8.0)*sigTdis
         sigLp = (10.0-w2)*sigLp+(w2-8.0)*sigLdis
         sigTp = sigTp/2.0
         sigLp = sigLp/2.0
         else
          sigTp = sigTdis
          sigLp = sigLdis
        endif
      endif

      return
      end



      SUBROUTINE RESCSN(w2,q2,sigtn,sigln)
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i

      data xvaln / 
     & 0.12300E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13223E+00,0.21779E+00,0.10504E+00,0.19322E+00,
     & 0.16957E+00,0.38000E+00,0.76000E+01,0.44484E+01,0.22412E+01,
     & 0.18607E+01,0.13787E-02,0.98347E+04,0.18160E+00,0.26770E+01,
     & 0.34272E+00,0.18609E-05,0.47226E+01,0.63725E-06,0.23236E-01,
     & 0.84561E+03,0.29815E+00,0.27082E+01,0.00000E+00,0.10000E+01,
     & 0.10000E+01,0.20000E+01,0.67680E+01,0.56785E+01,0.43946E+00,
     & 0.56784E+01,0.20954E+03,0.18496E+00,0.16769E+01,0.15454E+00,
     & -.83108E+02,0.30574E+01,0.10000E-02,0.93889E+00,-.89302E-02,
     & -.62215E-01,0.19800E+01,0.45000E+00,0.40853E+01,0.14462E+00,
     & 0.10046E+01,0.99612E+00,0.99930E+00,0.99357E+00,0.10192E+01,
     & 0.10000E+01,0.99669E+00,0.10002E+01,0.10000E+01,0.10000E+01,
     & 0.10094E+01,0.10064E+01,0.16270E+03,0.86884E+01,0.57864E+01,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.30131E+02,0.55479E+01,0.93622E+00,0.00000E+00,0.76450E+02,
     & 0.66412E-01,0.85428E+01,0.00000E+00,0.00000E+00,0.20000E+01,
     & 0.10000E+01,0.00000E+00,0.23539E+01,0.10360E-02,0.67452E+00,
     & 0.00000E+00,0.23312E+03,0.35946E+00,0.23000E+01,-.30155E-01,
     & 0.59653E+01,-.10616E+01,0.00000E+00,0.00000E+00,0.13868E+04,
     & 0.29402E+03,0.39655E+00,0.00000E+00,0.26453E+00,0.00000E+00 /


      do i=1,50
        xval1(i) = xvaln(i)
        xvalL(i) = xvaln(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodn(1,w2,q2,xval1,sigtn)
      call resmodn(2,w2,q2,xvalL,sigLn)

      return
      end

      SUBROUTINE RESMODP(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       
      if(w2.LT.1.159) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 
  


      SUBROUTINE RESMODN(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.939565

c      mp = 0.938272

      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)
      if(w.LT.(mp+mpi)) wdif(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) wdif(2) = 0.0

c      write(6,*) "here"

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+2.*mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LT.(mp+mpi)) xpr(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) xpr(2) = 0.0

      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)*dip

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
c     &      (1.-xb)**(nr_coef(i,2)+0.0*t )
     &      (1.-xpr(1))**(nr_coef(i,2)+nr_coef(i,3)*log(q2+m0))
c     &        *xpr(1)**(xval(41)+xval(42)*t)
     &        *xb**(xval(41)+xval(42)*log(q2+m0)+
     &                  nr_coef(i,4)*q2*sqrt(q2))
c     &        *(1./(1.+q2/nr_coef(i,3)))**nr_coef(i,4)
        enddo

CCC  next line is test CCC nonres3 directory ******

c        sig_nr = sig_nr + nr_coef(i,3)/(1.+q2/nr_coef(i,4))**2*
c     &               (1.-xpr(1))**2

        

      endif


      sig_res = abs(sig_res)
      sig_nr = abs(sig_nr)

      sig = sig_res + sig_nr
       
      if(w.LT.(mp+mpi)) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 

      SUBROUTINE DISP(w2,q2,sigt,sigl)
      IMPLICIT none

      real*8 w2,q2,q2t,x,sigt,sigl,f1,f2,fL,r,dr,r1,r2
      Integer i,modt
      character*1 targ
      real*8 mp2,pi,pi2,alpha,t1,t2
      logical goodfit


      targ = 'P'
      modt = 12
      mp2 = 0.938272
      mp2 = mp2*mp2
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.0359

      x = q2/(w2+q2-mp2)

      call f2allmer(x,q2,f2)
      call f2glober(x,q2,targ,modt,f2)
      call r1998er(x,q2,r,dr,goodfit)
      if(q2.LE.0.15) then
        q2t = 0.15
        call r1998er(x,q2t,r1,dr,goodfit)
        r = r1/0.25*q2
c        write(6,*) r1,r
        call f2allmer(x,q2,f2)
      endif


      f1 = f2/2./x/(r+1.0)*(1.0+4.0*mp2*x*x/q2)
      fL = 2.*x*r*f1    

      sigt = 0.3894e3*f1*pi2*alpha*8.0/abs(w2-mp2)
      sigL = r*sigt

      return
      end

!Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!            Ph.D. Thesis, Stanford University,                                 
!            March 1990.                                                        
!For details see file HELP.DOCUMENT.                                            
                                                                                
!Program contains 145 lines of Fortran code, of 72 characters each, with        
!no subroutines.  Program requires File E.14 as input.                          
                                                                                
                                                                                
      SUBROUTINE F2GLOBER(X,Q2,Target,MODEL,F2)        
                                                                                
! Returns F2 and related quantities from the either the LAMBDA12 model          
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.                   
!                                                                               
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2          
!                                                                               
! Further, program returns uncertainty in F2 based on both statistics           
! (ST) and systematic effects due to our choice of models (SY).                 
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical        
! uncertainty in this slope.                                                    
!                                                                               
! Best model is LAMBDA12.  SY is estimated to be the difference between         
! the two models.                                                               
!                                                                               
! Errors due to overall normalization are not included in program output        
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.              
! Systematic errors due to radiative corrections are shown in Reference         
! to be very kinematic independent, and are everywhere <.5%, and thus,          
! are ignored by this routine (also see documentation to dFRC in file           
! HELP.DOCUMENT).                                                               
!                                                                               
! Coefficients and correlation matrix elements are from File                    
! E.13 F1990.MATRICES, dated 27Jan90.                                           
                                                                                
      IMPLICIT NONE                                                             
      LOGICAL GOODFIT,FIRST                                             
      REAL*8   X, Q2, F2,ST,SY,SLOPE,DSLOPE                               
      REAL*8   XPP,XP,Y,POLY,F2B,POL1,STB,Q,QTH,DQ,F2TH,QUAD,SCB,F2L      
      REAL*8   BINDING,STL                                                
      INTEGER MODEL, I, J, K                                                    
      CHARACTER*1 TARGET                                                        
      REAL*8    B(2,9),L(2,9,9)                      ! OMEGA9 variables         
      REAL*8    C(2,12),M(2,12,12)                   ! LAMBDA12 variable        
      REAL*8 V(9),Z(12),U(12),LIN                                               
                             
      FIRST = .true.                                                   
                                                                                
! Model #9  27 Jan 90.                            !  HYDROGEN Ci            
      DATA (B(1,J),J=1,9)/                        
     >   0.7338659870D0,  11.0245522588D0,   2.6185804129D0,                    
     >   4.0956321483D0,   0.1206495422D0,   1.9714128709D0,                    
     >   3.8893348719D0, -14.0507358314D0,   8.8080576075D0/
!                                                 !  HYDROGEN MijDiDj
      DATA ((L(1,J,K),K=1,9),J=1,9)/              
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,                     
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,                     
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,                     
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,                     
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,                     
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,                     
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,                     
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,                     
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,                     
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,                     
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,                     
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,                     
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,                     
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,                     
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,                     
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,                     
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,                     
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,                     
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,                     
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,                     
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,                     
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,                     
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,                     
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,                     
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,                     
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,                     
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/                     
                                                                                
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci               
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,                     
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,                     
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/                     
                                                                                
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj          
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,                     
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,                     
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,                     
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,                     
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,                     
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,                     
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,                     
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,                     
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,                     
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,                     
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,                     
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,                     
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,                     
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,                     
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,                     
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,                     
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,                     
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,                     
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,                     
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,                     
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,                     
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,                     
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,                     
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,                     
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,                     
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,                     
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/                     
                                                                                
!    MODEL #12:    27Jan90.                                                     
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci                    
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,                     
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,                     
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,                     
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/                     
      DATA ((M(1,J,K),K=1,12),J=1,6)/     !     HYDROGEN MijDiDj               
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,                     
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,                     
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,                     
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,                     
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,                     
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,                     
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,                     
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,                     
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,                     
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,                     
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,                     
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,                     
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,                     
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,                     
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,                     
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,                     
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,                     
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,                     
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,                     
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,                     
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,                     
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,                     
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,                     
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0/                     
        DATA ((M(1,J,K),K=1,12),J=7,12)/             ! hydrogen
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,                     
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,                     
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,                     
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,                     
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,                     
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,                     
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,                     
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,                     
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,                     
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,                     
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,                     
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,                     
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,                     
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,                     
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,                     
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,                     
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,                     
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,                     
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,                     
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,                     
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,                     
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,                     
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,                     
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/                     
                                                                                
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci                   
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,                     
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,                     
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,                     
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/                     
      DATA ((M(2,J,K),K=1,12),J=1,6)/     !     DEUTERIUM MijDiDj              
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,                     
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,                     
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,                     
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,                     
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,                     
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,                     
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,                     
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,                     
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,                     
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,                     
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,                     
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,                     
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,                     
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,                     
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,                     
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,                     
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,                     
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,                     
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,                     
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,                     
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,                     
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,                     
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,                     
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0/                     
        DATA ((M(2,J,K),K=1,12),J=7,12)/             ! deuterium
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,                     
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,                     
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,                     
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,                     
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,                     
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,                     
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,                     
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,                     
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,                     
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,                     
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,                     
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,                     
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,                     
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,                     
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,                     
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,                     
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,                     
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,                     
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,                     
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,                     
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,                     
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,                     
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,                     
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/                       
      !----------------------------------------------------------------         
      !----------------------------------------------------------------         
                                                                                
                                                                                                                          
      i = 1                                                                     
      IF (TARGET.EQ.'D') i = 2                                                  
      BINDING = 1./(1.-EXP(-MIN(20.0,7.7*(1./X+.93828**2/Q2-1.))))               
      IF (i.EQ.1) BINDING = 1.                                                  
                                                                                
      !OMEGA9 MODEL FIRST:                                                      
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                                     
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                                     
           Y    = 1.-XP                                                         
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+              
     >            B(i,9)*Y**7                                                   
           F2B  = X/XPP*BINDING*POLY                                            
          !-----------------------------------------------------------          
          !V(k) is the derivative of F_2 with respect to parameter k.           
           V(1) = -F2B/XPP/(Q2/X+B(i,2))                                        
           V(2) =  F2B/XPP/(Q2/X+B(i,2))**2*(Q2+B(i,1))                         
           POL1 =  3.*B(i,5)*Y**2+4.*B(i,6)*Y**3+5.*B(i,7)*Y**4+                
     >             6.*B(i,8)*Y**5+7.*B(i,9)*Y**6                                
           V(3) = -F2B*POL1/POLY/(Q2/X+B(i,4))                                  
           V(4) =  F2B*POL1/POLY/(Q2/X+B(i,4))**2*(Q2+B(i,3))                   
           DO 10 j = 5,9                                                        
10         V(j) =  F2B/POLY*Y**(j-2)                                            
           STB = 0.                                                             
           DO 11 j = 1,9                                                        
           DO 11 k = 1,9                                                        
11         STB = STB + L(i,j,k)*V(j)*V(k)                                       
           STB = SQRT(STB)*BINDING                                              
                                                                                
      !LAMBDA12 MODEL NEXT:                                                     
           Y    = 1.-X                                                          
           q    = LOG(Q2)                                                       
           qth  = .2+3.2*X                                                      
           dq   = q-qth                                                         
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                          
     >            C(i,4)*Y**6+C(i,5)*Y**7                                       
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                           
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq               
           IF (q.GT.qth) QUAD = 0.                                              
           SCB  = (1.+LIN+QUAD)                                                 
           F2L  = F2th*SCB*BINDING                                              
          !-----------------------------------------------------------          
          !Z(k) is the derivative of F_2 with respect to parameter k.           
           DO 20 j = 1,5                                                        
20         Z(j) = SCB*Y**(j+2)                                                  
           Z(6) = 0.                                                            
           IF (q.LT.qth) Z(6) = F2th*dq**2                                      
           Z(7) = Z(6)*X                                                        
           Z(8) = Z(6)*X**2                                                     
           DO 21 j = 9,12                                                       
21         Z(j) = F2th*X**(j-9)*dq                                              
           STL = 0.                                                             
           DO 22 j = 1,12                                                       
           DO 22 k = 1,12                                                       
22         STL = STL + M(i,j,k)*Z(j)*Z(k)                                       
           STL = SQRT(STL)*BINDING                                              
                                                                                
          !U(k) is the derivative of slope with respect to parameter k.         
           SLOPE= F2th*LIN/dq*BINDING                                           
           DO 30 j = 1,5                                                        
30         U(j) = LIN/dq*Y**(j+2)                                               
           DO 31 j = 6,8                                                        
31         U(j) = 0.                                                            
           DO 32 j = 9,12                                                       
32         U(j) = Z(j)/dq                                                       
           DSLOPE = 0.                                                          
           DO 33 j = 1,12                                                       
           DO 33 k = 1,12                                                       
33         DSLOPE = DSLOPE + M(i,j,k)*U(j)*U(k)                                 
           DSLOPE = SQRT(DSLOPE)                                                
      !----------------------------------------------------------------         
                                                                                
      F2 = 0.                                                                   
      ST = 0.                                                                   
      IF (MODEL.EQ. 9) THEN                                                     
           F2 = F2B                                                             
           ST = STB                                                             
      ELSEIF (MODEL.EQ.12) THEN                                                 
           F2 = F2L                                                             
           ST = STL                                                             
      ELSE                                                                      
           WRITE(*,'('' F1990: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')            
      ENDIF                                                                     
      SY = ABS(F2B-F2L)                                                         
                                                                                
      GOODFIT = .TRUE.                                                          
      !The following cuts define the region of applicability of F1990.          
      !In order they are:                                                       
      !     [radiative corrections convergence criteria in Q2] .and.            
      !     [radiative corrections convergence criteria in x]  .and.            
      !     [stay out of resonance region, W2.ge.3.0]          .and.            
      !     [limitation imposed by maximum beam energy].                        
      IF ((Q2.LT..566).OR.(X.LT..062).OR.(X.LT.Q2/(2.*.93828*21.))              
     >   .OR.(X.GT.1./((3.-.93828**2)/Q2+1.)))     THEN                         
                                                                                
C         WRITE(*,'('' WARNING[F1990]: OUTSIDE RECOMMENDED RANGE.'')')          
          GOODFIT=.FALSE.                                                       
      ENDIF                                                                     
                                                                                
      RETURN                                                                    
      END                                                                       

! Modified by Yongguang Liang to include the R1998 fitting parameters
! 
!File R1990.FORTRN.                                                       
!Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!            Ph.D. Thesis, Stanford University,                                 
!            March 1990.                                                        
!For details see file HELP.DOCUMENT.                                            
                                                                                
!Program contains 135 lines of Fortran code, of 72 characters each, with        
!one subroutine.                                                                
                                                                                
                                                                                
      SUBROUTINE R1998er(X,Q2,R,DR,GOODFIT)                                       
                                                                                
! Model for R, based on a fit to world R measurements. Fit performed by         
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details           
! see Reference.                                                                
!                                                                               
! Three models are used, each model has three free parameters.  The             
! functional forms of the models are phenomenological and somewhat              
! contrived.  Each model fits the data very well, and the average of            
! the fits is returned.  The standard deviation of the fit values is            
! used to estimate the systematic uncertainty due to model dependence.          
!                                                                               
! Statistical uncertainties due to fluctuations in measured values have         
! have been studied extensively.  A parametrization of the statistical          
! uncertainty of R1990 is presented in FUNCTION DR1990.                         
!                                                                               
! The three model forms are given by:                                           
!                                                                               
!     R_A = A(1)/LOG(Q2/.04)*FAC + A(2)/[Q2\E34+A(3)\E34]\E3.25 ;                     
!     R_B = B(1)/LOG(Q2/.04)*FAC + B(2)/Q2 + B(3)/(Q2**2+.3**2) ;               
!     R_C = C(1)/LOG(Q2/.04)*FAC + C(2)/[(Q2-Q2thr)\E32+C(3)\E32]\E3.5 ,              
!                               ...where Q2thr = 5(1-X)\E35 ;                     
!           where FAC = 1+12[Q2/(1+Q2)][.125\E32/(.125\E32+x\E32)] gives the          
!           x-dependence of the logarithmic part in order to match Rqcd         
!           at high Q2.                                                         
!                                                                               
! Each model fits very well.  As each model has its own strong points           
! and drawbacks, R1990 returns the average of the models.  The                  
! chisquare for each fit (124 degrees of freedom) are:                          
!     R_A: 110,    R_B: 110,    R_C: 114,    R1990(=avg): 108                   
!                                                                               
! This subroutine returns reasonable values for R for all x and for all         
! Q2 greater than or equal to .3 GeV.                                           
!                                                                               
! The uncertainty in R originates in three sources:                             
!                                                                               
!     D1 = uncertainty in R due to statistical fluctuations of the data         
!          and is parameterized in FUNCTION DR1990, for details see             
!          Reference.                                                           
!                                                                               
!     D2 = uncertainty in R due to possible model dependence, approxi-          
!          mated by the variance between the models.                            
!                                                                               
!     D3 = uncertainty in R due to possible epsilon dependent errors            
!          in the radiative corrections, taken to be +/- .025.  See             
!          theses (mine or Dasu's) for details.                                 
!                                                                               
! and the total error is returned by the program:                               
!                                                                               
!     DR = is the total uncertainty in R, DR = sqrt(D1\E32+D2\E32+D3\E32).            
!          DR is my best estimate of how well we have measured R.  At           
!          high Q2, where R is small, DR is typically larger than R.  If        
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,          
!          you might wish to assume DR = 0 at very high Q2.                     
!                                                                               
! NOTE:    In many applications, for example the extraction of F2 from          
!          measured cross section, you do not want the full error in R          
!          given by DR.  Rather, you will want to use only the D1 and D2        
!          contributions, and the D3 contribution from radiative                
!          corrections propogates complexely into F2.  For more informa-        
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or             
!          for explicite detail, see Reference.                                 
!                                                                               
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      IMPLICIT NONE                                                             
      REAL*8 FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1990,X,Q2,W2                
      REAL*8 A(3), B(3), C(3), MP,MP2
      LOGICAL GOODFIT                                                           
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      DATA A/ .06723, .46714, 1.89794 /                                  
      DATA B/ .06347, .57468, -.35342 /                                   
      DATA C/ .05992, .50885, 2.10807 /

      MP = .9382727
      MP2 = MP*MP                                                                        

      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))                       
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!           
c      Q2thr = 5.*(1.-X)**5                                                      
                                                                                
c      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))                        
c      R_B   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)                          
c      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)                      


       W2 = Q2/X-Q2+MP2

c       IF (W2 .GT. 3.9 ) THEN 
       q2thr = 12.3708*x-43.1043*x**2  !! r1998
     & +41.7415*x**3

      r_a = 0.0485*rlog
     &+0.5470/sqrt(sqrt(q2**4+2.0621**4))
     &*(1.-0.3804*x+0.5090*x**2)*x**(-0.0285) 
      r_b = 0.0481*rlog
     &+(0.6114/q2-0.3509/(q2**2+.3**2))
     &*(1.-0.4611*x+0.7172*x**2)*x**(-0.0317)                    
      r_c = 0.0577*rlog
     & + 0.4644/sqrt((q2-q2thr)**2+1.8288**2)
c      ELSE
c       q2thr = -5.8119*x+0.9259*x**2  !! r1990
c     & +9.4095*x**3

c      r_a = 0.2839*rlog
c     &+0.0789/sqrt(sqrt(q2**4+0.0**4))
c     &*(1.+32.2507*x-28.5938*x**2)*x**(2.8009) 
c      r_b = 0.2518*rlog
c     &+(0.0110/q2-0.0069/(q2**2+.3**2))
c     &*(1.+57.0509*x-1.4183*x**2)*x**(-0.0934)                    
c      r_c = 0.0621*rlog
c     & + 1.4431/sqrt((q2-q2thr)**2+7.1946**2)
c       write(*,*) r_a, r_b, r_c 
c      ENDIF
      
      R     = (R_A+R_B+R_C)/3.                                                  
                                                                                
      D1    = DR1990(X,Q2)                                                      
      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)                       
      D3    = .023*(1.+.5*R)                                                    
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                               
      DR    = SQRT(D1**2+D2**2+D3**2)                                           
                                                                                
      GOODFIT = .TRUE.                                                          
      IF (Q2.LT..3) GOODFIT = .FALSE.                                           
      RETURN                                                                    
      END                                                                       
                                                                                
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
                                                                                
      FUNCTION DR1990(X,Q2)                                                     
                                                                                
! Parameterizes the uncertainty in R1990 due to the statistical                 
! fluctuations in the data.  Values reflect an average of the R-values          
! about a neighborhood of the specific (x,Q2) value.  That neighborhood         
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see               
! Reference.                                                                    
!                                                                               
! This subroutine is accurate over all (x,Q2), not only the SLAC deep           
! inelastic range.  Where there is no data, for example in the resonance        
! region, it returns a realistic uncertainty, extrapolated from the deep        
! inelastic region (suitably enlarged).  We similarly estimate the              
! uncertainty at very large Q2 by extrapolating from the highest Q2             
! measurments.  For extremely large Q2, R is expected to fall to zero,          
! so the uncertainty in R should not continue to grow.  For this reason         
! DR1990 uses the value at 64 GeV for all larger Q2.                            
!                                                                               
! XHIGH accounts for the rapidly diminishing statistical accuracy for           
! x>.8, and does not contribute for smaller x.                                  
                                                                                
                                                                                
      IMPLICIT NONE                                                             
      REAL*8 U(10,10),DR1990,QMAX,Q,S,A,XLOW,XHIGH,X,Q2                           
                                                                                
                                                                                
      QMAX = 64.                                                                
                                                                                
      Q = MIN(Q2,QMAX)                                                          
      S = .006+.03*X**2                                                         
      A = MAX(.05,8.33*X-.66)                                                   
                                                                                
      XLOW  = .020+ABS(S*LOG(Q/A))                                              
      XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)                          
                                                                                
      DR1990 = SQRT(XLOW**2+XHIGH**2)                                           
      RETURN                                                                    
      END                                                                       
!                                                                               
!End of file R1990.FORTRN.  135 Fortran lines.                                  





c       allm97, NMC published measured points Q2>0.75 GeV2
c       for values Q<1 use data of E665!
c       parameterization of F2 , according to
c       H.Abramowicz and A.Levy, hep-ph/9712415
c
c       3*10-6 < x  < 0.85, W2>3GeV2
c       0.   < Q2 < 5000 GeV2, dof=0.97
c
 
      SUBROUTINE f2allmer(x,q2,f2a)

      IMPLICIT NONE
 
      REAL*8 x,q2,M22,f2a
      REAL*8 SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      REAL*8 S11,A11,B11,M12,S21,A21,B21,M02,LAM2,Q02
      REAL*8 S12,S13,A12,A13,B12,B13,S22,S23,A22,A23
      REAL*8 B22,B23,w2,w,z
      REAL*8 ALFA,XMP2
C  POMERON
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )         
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      
C
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
        
C
      IF(Q2.EQ.0.) THEN
       S=0.
       Z=1.

C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11
       BP=B11
       SP=S11
       F2P=SP*XP**AP
C                                               
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21
       BR=B21
       SR=S21
       F2R=SR*XR**AR
C
      ELSE
       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))
       Z=1.-X   
C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)
       BP=B11+B12*S**B13
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)
       F2P=SP*XP**AP*Z**BP
C
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21+A22*S**A23
       BR=B21+B22*S**B23
       SR=S21+S22*S**S23
       F2R=SR*XR**AR*Z**BR
 
C
      ENDIF                                     
      
 
      f2a = q2/(q2+m02)*(F2P+F2R)
 
 
      RETURN
      END                                  


      

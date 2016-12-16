c     **********************************************************
c     ** Automaton S+D (phase transition + dislocations)      **
c     **                                                      **
c     ** Transforming elements: 2D square lattice             **
c     **                                                      **
c     ** Open boundaries                                      **
c     **                                                      **
c     **      FJPR 26/06/2006                                 **
c     **                                                      **
c     **      12/07/2006                                      **
c     **      General interaction along axis and diagonals    **
c     **                                                      **
c     ** Adapted to use long-range Gaussian kernel (SK model) **
c     **********************************************************

c     **
c	BRUTE FORCE

c     Driven by temperature
c     

c	Th0(s(i),change): threshold corresponding to the update of a spin s(i)
c	zsh(s(i),change): change in z associated to a certain update
c	sch(s(i),change):    "          s                       "
c	uns(): unstable sites at time t

!     Spanning avalanches analysis partially implemented for cooling runs
!     Unit 120: Output whether there is a 1D or 2D spanning avalanche in a series of runs
!     Size distributions of non-spanning and spanning avalanches.

      parameter(long=200)
      parameter(nTmax=10000)
      parameter(big=1.d15)
      parameter(small1=1.d-3)
      integer L,NL,nruns,z1,z2,limav,ncycle,nava
      integer xini,yini,nini,site

      integer nup,ndw

      integer s(long**2),sch(-1:1,2)
      double precision zchnn(long**2)
      double precision Th0(-1:1,2),zch(-1:1,2),signh(-1:1,2),Th1,Th2
      double precision z(long**2),h,overh,fD,feed,h1,h2,h3,h4,haux

      integer ngrid,igrid
      double precision xgrid(500),qini(500),q
      double precision J1mat(500),J2mat(500),J3mat(500)
      character*32 extension(500)

      integer dis(long**2)
      double precision fDTemp
      double precision f,J1,J2,J3,zeff,zabs
      integer uns(long**2),nuns
      integer nunsPT  ! Number of unstable with respect to phase transition.
      double precision strain,density
      real sigma
      double precision r,raux,Di,sgi,Dini
      double precision rhograd
      real cgrad(4,2)

      integer ProxG,ProxE,Fincapa
      integer size,t,sPT,tPT,sD,nshell,stress,ndis,absndis
      integer sums,sumt,sumsPT,sumsD,nflip,maxA,sumA,MaxsPT
      integer sumsCool,sumtCool,sumsPTCool
      integer sumsCoolNS,sumsCool1D,sumsCool2D
      integer sumdis,sumdis2,sumdis3
      double precision skew
      integer nshapes
      integer Hs(2*long**2),Ht(long**2),HsPT(2*long**2),HsD(long**2)
      integer HsCool(2*long**2),HtCool(long**2),HsPTCool(2*long**2)
      integer HsCoolNS(2*long**2),HsCool1D(2*long**2),HsCool2D(long**2)
      integer HA(long**2)
      integer shape(nTmax),shapesum(nTmax,50),sumshape(50),tshape(50)
      integer shapePT(nTmax)
      integer shapeD(nTmax)
      integer j(long**2,0:8) 
      integer i,n,n11,n22,nx,ny,nz
      integer coorx,coory
      integer nthermal
      integer x0,y0,x,y 

      character*12 filename
      character*5 shapename(50)
      character*1 increase,mapST
      character*5 randomtype

      integer initseed
      real*4 rvec,rv
      dimension rvec(108024)
      
      double precision sdi,sdi2,sti,sti2,stidi,ni,rho,dii,ti,nt
      character*1 zout
      character*3 BC !Boundary conditions DBC: dissipative (by default). FBC: fixed (number of grains conserved) 
      character*1 stress_out
      character*1 sdout

      integer aux
      double precision daux
      character*4 Buffer
      character us,ds,cs,ms,ts
      character*1 stab_out
      dimension hstab(0:20000,1:3)
      dimension sum_stab(1:3)
      dimension hdstab(0:20000,1:3)
      dimension sum_dstab(1:3)
      integer tglob

!      dimension DJij(long**2,long**2)
      dimension DJij(1,1)
      character*3 Lchar
      character*7 c_name(500)
      dimension cvalue(500)

      integer xmin,xman

      integer id
      double precision wtime
 
      dimension h_branch(0:1000) !histogram for branching rate
      dimension h_branch_g(0:1000,0:1000)
      dimension h_branchPT_g(0:1000,0:1000)
      dimension h_branchD_g(0:1000,0:1000)
      dimension obs_branch_g_c(0:1000),obs_branch_g_h(0:1000)
      dimension obs_branchPT_g_c(0:1000)
      dimension obs_branchD_g_c(0:1000)
      dimension Sum1_branch_g_c(0:1000),Sum1_branch_g_h(0:1000)
      dimension Sum1_branchPT_g_c(0:1000),Sum2_branchPT_g_c(0:1000)
      dimension Sum1_branchD_g_c(0:1000)

      dimension h_rf_g(0:1000,0:1000)
      dimension sum_rf_g(0:1000)

      dimension hstab_g(0:1000,0:1000,1:3)
      dimension sum_stab_g(0:1000,1:3)

      dimension Deltai(long**2) !Change of elastic strain after each configurational change.

      dimension hDeltai_g(0:1000,0:1000)
      dimension sum_Deltai_g(0:1000)

      dimension hDeltai(0:1000)

      dimension hdnflip_g(0:1000,0:1000)
      dimension sum_dnflip_g(0:1000)

      integer Hs_g_c(2*long**2,0:1000)
      integer HsPT_g(2*long**2,0:1000)

      integer HsNS_g_c(2*long**2,0:1000)
      integer Hs1D_g_c(2*long**2,0:1000)
      integer Hs2D_g_c(2*long**2,0:1000)

      integer Hs_g_h(2*long**2,0:1000)

      integer sums_g

      double precision rf(long**2) !effective random fields
      double precision E_rf(0:1000,4)
      integer num_rf(0:1000)
      double precision E_rf_h(0:1000,4)
      integer num_rf_h(0:1000)

      double precision E_Js_c(0:1000,4)

      character*1 Dis_cyc_out

      real*8 g0_h(0:2*long)
      integer N0_h(0:2*long)
      real*8 g_h(0:2*long)
      integer N_h(0:2*long)
      integer ix0,iy0

      character*1 forced_MA
      integer i_snap_cool,i_snap_heat,size_snap

      integer xWmin,xWmax,yWmin,yWmax
      integer xSmin,xSmax,ySmin,ySmax
      integer sox(0:long),soy(0:long)

      character*2 SR_LR
      integer SRtype
      dimension DJSR(0:8)
      
      character*1 temp_out

      integer s_M_prev(long**2),dis_M_prev(long**2) ! Martensite configuration fields in previous cycle 

      double precision f_s,f_d,f_w
      integer nflip_abs
      double precision cumm_f  !Cummulative fraction of martensite.

! output of a list with random fields for specific temperature
      character*1 rf_raw_out
      double precision tau_rf_out,tol_rf,hmin_g2,hmax_g2
      integer nc_min,nc_max,nb_g2

      integer sums_g_c,sumsNS_g_c,sums1D_g_c,sums2D_g_c

!     ---- Fixed parameters for histograms of branching, stabilities and temperature

!     -- Branching ratio
      nb_rf=500
      rf_min=-2.d0
      rf_max=2.d0
      aa_rf=(rf_max-rf_min)/dble(nb_rf)

!     -- Branching ratio
      nb_branch=500
      hmin_branch=0.d0
      hmax_branch=5.d0
      aa_branch=(hmax_branch-hmin_branch)/dble(nb_branch)

!     -- Stability
      nb_stab=1000
      ddmin=-1.d0
      ddmax=0.5d0
      aa=(ddmax-ddmin)/dble(nb_stab)
      nb_delta=1000
      ddeltamin=-0.5d0
      ddeltamax=0.5d0
      ad=(ddeltamax-ddeltamin)/dble(nb_delta)

!     -- Deltai - total change in local strain
      nb_Deltai=1000
      Deltai_min=-0.7d0
      Deltai_max=0.7d0
      aa_Deltai=(Deltai_max-Deltai_min)/dble(nb_Deltai)

!     -- Hist-onflip vs g
      nb_dnflip=100
      dnflip_min=0.d0
      dnflip_max=1.d0
      aa_dnflip=(dnflip_max-dnflip_min)/dble(nb_dnflip)

!     -- Temperature 
      nb_g=200
! -- Wide interval
!      hmin_g=-1.d0/8.d0
!      hmax_g=1.d0/8.d0

      open(20,file='SDNNNMulti_KernelGlobal_And_ShortRange.ini')
      read(20,*) L,Lchar
      read(20,*) SR_LR,SRtype
      read(20,*) n_c
      do i_c=1,n_c
         read(20,*) cvalue(i_c),c_name(i_c)
      enddo
      read(20,*) BC
      read(20,*) f
      read(20,*) fD
      read(20,*) fDTemp
      read(20,*) feed
      read(20,*) temp_out
      read(20,*) forced_MA
      read(20,*) randomtype
      read(20,*) r_Di
      read(20,*) coorx
      read(20,*) coory
      read(20,*) Dini
c      read(20,*) overh
      read(20,*) iseed
      read(20,*) ncycle
      read(20,*) limav
      read(20,*) size_snap
      read(20,*) filename
      read(20,*) mapST
      read(20,*) increase
      read(20,*) zout
      read(20,*) stress_out
      read(20,*) sdout
      read(20,*) Dis_cyc_out
      read(20,*) stab_out
      read(20,*) ngrid
      do i=1,ngrid
         read(20,*) J1mat(i),J2mat(i),J3mat(i),xgrid(i),
     &        qini(i),extension(i)
      enddo      
      read(20,*) goutmin,goutmax
      read(20,*) nshapes
      if(nshapes.gt.0)then
      do i=1,nshapes
         read(20,*) tshape(i),shapename(i)
      enddo
      endif
      read(20,*) rf_raw_out
      read(20,*) tau_rf_out,tol_rf,nc_min,nc_max,hmin_g2,hmax_g2,nb_g2
      close(20)
      
      
      
      print*,fDTemp
      open(20,file='results/'//filename//'/'//filename
     &     //extension(1)//'_'//BC//'_'//c_name(1)//'.info')
      write(20,*) 'L=',L,Lchar
      write(20,*) 'Interaction range, type (for SR):', SR_LR,SRtype
      
      write(20,*) 'Number of values of c:',n_c
      do i_c=1,n_c
         write(20,*) cvalue(i_c),c_name(i_c)
      enddo
      write(20,*) 'Boundary conditions:', BC
      write(20,*) 'f=',f
      write(20,*) 'J1,J2,J3=',J1,J2,J3
      write(20,*) 'fD=',fD
      write(20,*) 'fDTemp=',fDTemp
      write(20,*) 'feed=',feed
      write(20,*) 'temp_out=',temp_out
      write(20,*) 'randomtype=',randomtype
      if(randomtype.ne.'DiDET') write(20,*) 'r_Di=',r_Di
      if(randomtype.eq.'DiDET'.or.randomtype.eq.'sgGAU'
     &     .or.randomtype.eq.'runif')then
         write(20,*) 'coorx, coory =',coorx,coory
         write(20,*) 'Dini=',Dini
      endif
      write(20,*) 'Seed for random number generator=',iseed
      write(20,*) 'ncycle=',ncycle
      write(20,*) 'limav=',limav
      write(20,*) 'size_snap=',size_snap
      write(20,*) 'mapST=',mapST
      write(20,*) 'increase=',increase
      write(20,*) 'ngrid=',ngrid
      write(20,*) 'zout=',zout
      write(20,*) 'stress_out=',stress_out
      write(20,*) 'sdout=',sdout
      write(20,*) 'Dis_cyc_out=',Dis_cyc_out
      write(20,*) 'stab_out=',stab_out
      do i=1,ngrid
         write(20,*)J1mat(i),J2mat(i),J3mat(i),xgrid(i),
     &        qini(i),extension(i)
      enddo
      write(20,*) 'goutmin,goutmax=',goutmin,goutmax

      write(20,*) nshapes
      if(nshapes.gt.0)then
      do i=1,nshapes
         write(20,*) tshape(i),shapename(i)
      enddo
      endif

      write(20,*) 'rf_raw_out=',rf_raw_out
      write(20,*) 'tau_rf_out,tol_rf,nc_min,nc_max,hmi_g,hma_g,nb_g2=',
     &     tau_rf_out,tol_rf,nc_min,nc_max,hmin_g2,hmax_g2,nb_g2

!      close(20)

!      limav2=max(0,ncycle-10)
      limav2=ncycle+2
      limav3=max(0,ncycle-500)

!Output of parameters on the screen
      write(*,*) 'L=',L
      write(*,*) 'Number of values of c:',n_c
      do i_c=1,n_c
         write(*,*) cvalue(i_c),c_name(i_c)
      enddo
      write(*,*) 'Boundary conditions:', BC
      write(*,*) 'f=',f
      write(*,*) 'J1,J2,J3=',J1,J2,J3
      write(*,*) 'fD=',fD
      write(*,*) 'fDTemp=',fDTemp
      write(*,*) 'feed=',feed
      write(*,*) 'temp_out=',temp_out
      write(*,*) 'randomtype=',randomtype
      if(randomtype.ne.'DiDET') write(*,*) 'r_Di=',r_Di
      if(randomtype.eq.'DiDET'.or.randomtype.eq.'sgGAU')then
         write(*,*) 'coorx, coory =',coorx,coory
         write(*,*) 'Dini=',Dini
      endif
      write(*,*) 'Seed for random number generator=',iseed
      write(*,*) 'ncycle=',ncycle
      write(*,*) 'limav=',limav
      write(*,*) 'size_snap=',size_snap
      write(*,*) 'mapST=',mapST
      write(*,*) 'ngrid=',ngrid
      write(*,*) 'zout=',zout
      write(*,*) 'stress_out=',stress_out
      write(*,*) 'sdout=',sdout
      write(*,*) 'Stabilities out?:',stab_out
      do i=1,ngrid
         write(*,*) J1mat(i),J2mat(i),J3mat(i),xgrid(i),
     &        qini(i),extension(i)
      enddo

      NL=L*(L-1)
      z1=4
      z2=4
      zp=z1
      xmin=-1 !0
      xmax=L+1 !L-1

!----- Calculate the spanning window: [xWmin,xWmax]x[yWmin,yWmax]
      xWmin=L
      xWmax=0
      yWmin=L
      yWmax=0     
      do i=1,NL
         y=int(i/L)
         if(mod(i,L).eq.0) y=y-1
         x=i-y*L-1
         if(x.gt.xWmax)xWmax=x
         if(x.lt.xWmin)xWmin=x
         if(y.gt.yWmax)yWmax=y
         if(y.lt.yWmin)yWmin=y
      enddo

      open(118,file='results/'//filename//'/'//filename
     &     //extension(1)//'_'//BC//'_e_plasticity-delta_max.dat')

      do i_c=1,n_c

      call srand(iseed)

!      zeff=J1*float(z1)+J2*float(z2)
!      zabs=abs(J1)*float(z1)+abs(J2)*float(z2)
c     Initial configuration -----
c      do i=1,NL
c         z(i)=0.0d0
c	 s(i)=0
c      enddo

c     ----------------------------------------------------------------
c     ----------------------------------------------------------------
c     Nearest neighbours matrix j(n,c), {n=1,2,...NL}, {c=1,2,...z1}
c     -------------------------
c	OPEN BOUNDARIES

      do n=1,NL
         j(n,0)=n
!     Nearest neighbours
         j(n,1)=n-1
         j(n,3)=n+1
         j(n,2)=n-L
         j(n,4)=n+L
!     Next-to-nearest neighbours
         j(n,5)=n-L-1
         j(n,6)=n-L+1
         j(n,7)=n+L+1
         j(n,8)=n+L-1
      enddo 

c     Outgoing Boundary neighbours
c
c     x=0 and x=Lsite-1 planes
      do ny=0,L-2
         n11=1+ny*L
         n22=L+ny*L
         j(n11,1)=-1
         j(n11,5)=-1
         j(n11,8)=-1
         j(n22,3)=-1
         j(n22,6)=-1
         j(n22,7)=-1
c         j(n11,1)=n22
c         j(n22,3)=n11
      enddo

c     y=0 and y=Lsite-1 planes
      do nx=1,L
         n11=nx
         n22=nx+(L-2)*L
         j(n11,2)=-1
         j(n11,5)=-1
         j(n11,6)=-1
         j(n22,4)=-1
         j(n22,7)=-1
         j(n22,8)=-1
c         j(n11,2)=n22
c         j(n22,4)=n11
      enddo

c     ------------------------------------------------------------
c     ------------------------------------------------------------



!-------------------------------------------------
!----- Kernel J - Long range  --------------------
!-------------------------------------------------
      fJ=1.0d0
      print*,'Jmax=',fJ
      print*,'ATTENTION: Jmax set to 1'
      write(20,*)'ATTENTION: Jmax set to 1, i.e. disorder '
      close(20)

      if(SR_LR.eq.'LR')then
      open(11,file='Kernel/Jc_Global_'//c_name(i_c)//'_'
     &     //Lchar//'.dat')
      
      do i = 1,NL
         read(11,*)(DJij(i,ii), ii=1,NL)
      enddo 
      close(11)

! Statistics of the random field generated by a slip at site "it". The random field is obtained after deviding by e^2 (qini^2 in the notation of the code)
      it=coorx+(coory-1)*L
      avgJ_1=0.d0
      avgJ_2=0.d0
      do i=1,NL
            avgJ_1=avgJ_1+DJij(i,it)
            avgJ_2=avgJ_2+DJij(i,it)**2.d0
      enddo
      avgJ=avgJ_1/dble(NL)
      varJ=avgJ_2/dble(NL)-avgJ**2.d0
      sdJ=sqrt(varJ)

      print*,'avgJ,varJ,sdJ=',avgJ,varJ,sdJ

      zeff=0.5d0 
      zabs=0.d0
      print*,'zeff,zabs=',zeff,zabs
      endif

      open(104,file='results/'//filename//'/'//filename
     &     //'_'//BC//'_'
     &     //c_name(i_c)//'_'//extension(1)//
     &     '_Mean_Branching-Rate_g_vs_bare.dat')

      open(115,file='results/'//filename//'/'//filename
     &     //'_'//BC//'_'
     &     //c_name(i_c)//'_'//extension(1)//
     &     '_Deltah_Cooling_Min_Max_vs_e.dat')

      open(119,file='results/'//filename//'/'//filename
     &     //extension(1)//'_'//BC//'_'
     &     //c_name(i_c)//'_SNAP-POP.dat')

      open(120,file='results/'//filename//'/'//filename
     &     //extension(1)//'_'//BC//'_'
     &     //c_name(i_c)//'_0D_1D_2D_Spanning_Cooling.dat')

      open(121,file='results/'//filename//'/'//filename
     &     //extension(1)//'_'//BC//'_'
     &     //c_name(i_c)//'Branching_Realizations.dat')


      do igrid=1,ngrid

         print*,'igrid=',igrid,'of',ngrid
         print*,extension(igrid)

      if(zout.eq.'y') open(90,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//c_name(i_c)//'z.dat')
      if(mapST.eq.'y') open(91,file='results/'//filename
     &     //extension(igrid)//'/'//filename//'_'//BC//'_'
     &     //c_name(i_c)//'_mapST.dat')
      open(92,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_micre.dat')
      open(116,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.rf')
      open(93,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.s+')
      open(94,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.s-')
      open(95,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.d+')          ! positive slip
      open(96,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.d-')          ! negative slip
      open(97,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.gd')          ! gradient of slip
      open(98,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_rho.dat')          
      open(72,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_STdis.dat')
      if(sdout.eq.'y')open(73,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_map.sd')
      open(82,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Branching-Rate.dat')
!      open(83,file='results/'//filename//'/'//filename
!     &     //extension(igrid)//'_'//BC//'_'
!     &     //c_name(i_c)//'_Branching-Rate_All.dat')
      open(85,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Branching-Rate_Histo.dat')
      open(86,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Branching-Rate_Histo_g.dat')
      
      open(124,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Branching-Rate_PT_Histo_g.dat')

      open(87,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_sigma0_g.dat')
      open(88,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_sigma1_g.dat')
      open(89,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_sigmaS_g.dat')
      
      open(100,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_Deltai_g.dat')

      open(101,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_Deltai.dat')

      open(102,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Histo_nflip_g.dat')

      open(103,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_s_g_c.dat')

      open(34,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_s_NS_g_c.dat')

      open(35,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_s_1D_g_c.dat')

      open(36,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_s_2D_g_c.dat')

      open(105,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_sPT_g.dat')

      open(106,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_stat_rf_g_c.dat')

      open(107,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_stat_rf_cycles_Martensitic.dat')

      open(108,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_stat_rf_cycles_Austenite.dat')

      open(109,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_rf_g.dat')

      open(110,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_stat_rf_g_h.dat')

      open(111,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_s_g_h.dat')

      if(Dis_cyc_out.eq.'y')
     &     open(112,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_rf_cycles.dat')

      open(99,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Stat_branch_c.dat')

      open(127,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Stat_branchPT_c.dat')
      
      open(128,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Stat_branchD_c.dat')

      open(114,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Stat_branch_h.dat') 

      open(117,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Correlation_rf.dat')

      open(122,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_Active_Sites.dat')

      open(123,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//'_stat_Js_g_c.dat')

      if(rf_raw_out.eq.'y')then
         open(125,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'
     &        //c_name(i_c)//'_rf_raw_One-tau.dat')
         open(126,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'
     &        //c_name(i_c)//'_rf_raw_All.dat')
      endif

      q=qini(igrid)
      J1=J1mat(igrid)
      J2=J2mat(igrid)
      J3=J3mat(igrid)
      r=xgrid(igrid)

      
!-------------------------------------------------
!----- Kernel J - Short range  -------------------
!-------------------------------------------------
      if(SR_LR.eq.'SR')then
         if(SRtype.eq.1)then
            DJSR(1)=J1
            DJSR(2)=J2
            DJSR(3)=J1
            DJSR(4)=J2
            do i=5,8
               DJSR(i)=J3
            enddo
         endif
         if(SRtype.eq.2)then
            do i=1,4
               DJSR(i)=J1
            enddo
            do i=5,8
               DJSR(i)=J2
            enddo
         endif
c$$$         open(11,file='skel_JSR.in')
c$$$         do i=1,8
c$$$            read(11,*)aux1
c$$$            DJSR(i)=DJtemp(i)*aux1
c$$$         enddo
         aux1=0.d0
         aux2=0.d0
         do i=1,8
            aux1=aux1+DJSR(i)
            aux2=aux2+abs(DJSR(i))
         enddo
         DJSR(0)=-aux1
         zeff=aux1
         zabs=aux2
      endif

! -- Interval for T
      hmin_g=-(q**2.d0)/(2.d0) !-(q**2.d0)/(2.d0*(zeff+zabs))
      hmax_g=(q**2.d0)/(2.d0) !(q**2.d0)/(2.d0*(zeff+zabs))
      aa_g=(hmax_g-hmin_g)/dble(nb_g)
      hmin_g2=hmin_g2*(q**2.d0)
      hmax_g2=hmax_g2*(q**2.d0)
      aa_g2=(hmax_g2-hmin_g2)/dble(nb_g2)

      print*,'hmin_g2,hmax_g2,aa_g2=',hmin_g2,hmax_g2,aa_g2

!      stop
c     ----- Thresholds for dislocations ---------------
c
c      Authomatic
      Th0(-1,2)=-fD*(1.0d0-2.0d0*q)*(zeff+zabs)
      Th0(1,1)=fD*(1.0d0-2.0d0*q)*(zeff+zabs)
      print*,'Th0(1,1),Th0(-1,2)=',Th0(1,1),Th0(-1,2)
!      stop
c
c     By hand
c      Th0(-1,2)=-fD*2.0d0*float(z1)*J1
c      Th0(1,1)=fD*2.0d0*float(z1)*J1
c     -------------------------------------------
      
      zch(0,1)=q
      zch(0,2)=-q
      zch(-1,1)=q
      zch(-1,2)=-(1.0d0-2.0d0*q)
      zch(1,1)=(1.0d0-2.0d0*q)
      zch(1,2)=-q
      
c      signh(0,1)=1.0d0
c      signh(0,2)=-1.0d0
c      signh(-1,1)=0.0d0
c      signh(-1,2)=-1.0d0
c      signh(1,1)=0.0d0
c      signh(1,2)=1.0d0
      
      sch(0,1)=1
      sch(0,2)=-1
      sch(-1,1)=0
      sch(-1,2)=1
      sch(1,1)=-1
      sch(1,2)=0
c     ---------------------------------------------

c     ----------- matrix cgrad for coordinates of 'dual' lattice -------
      cgrad(1,1)=-0.5
      cgrad(1,2)=0.0
      cgrad(2,1)=0.0
      cgrad(2,2)=-0.5
      cgrad(3,1)=0.5
      cgrad(3,2)=0.0
      cgrad(4,1)=0.0
      cgrad(4,2)=0.5
c     ------------------------------------------------------------------

!--- vectors for correlations of slip field ---
      g_h(:)=0.d0
      N_h(:)=0.d0
!----------------------------------------------

      ivec=1
      raux=r/2.0d0
      raux_Di=r_Di/2.0d0
      
      do i=1,NL
         z(i)=0.0d0
         dis(i)=0
         s(i)=0
         rf(i)=0
      enddo



c     --------- Random numbers ------------
      initseed=REPLACESEED
      lcarry=108000
c     -------------------------------------
      call rcarin(initseed,rvec,lcarry)
      call rcarry(rvec,lcarry)

!      randomtype='DiDET'

      if(randomtype.eq.'rf_Di')then
         
         it=coorx+(coory-1)*L
         Di=Dini
         dis(it)=int(Dini)

         z(it)=z(it)-f*(J1*float(z1)+J2*float(z2))*Di
         rf(it)=rf(it)-f*(J1*float(z1)+J2*float(z2))*Di
         do k=1,z1
            if(j(it,k).gt.0)then
               z(j(it,k))=z(j(it,k))+J1*Di
               rf(j(it,k))=rf(j(it,k))+J1*Di
c     density=density+zchnn(site)
            endif
c     --- Open boundaries ---
            if(BC.eq.'FBC'.and.j(it,k).lt.0)then
               z(it)=z(it)+J1*Di
               rf(it)=rf(it)+J1*Di
            endif
c     if(j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------------------------
         enddo
         do k=1,z1-1
            if(j(j(it,k),k+1).gt.0.and.j(it,k).gt.0)then
               z(j(j(it,k),k+1))=z(j(j(it,k),k+1))
     &              +J2*Di
               rf(j(j(it,k),k+1))=rf(j(j(it,k),k+1))
     &              +J2*Di
            endif
c     --- Open boundaries ---
            if(BC.eq.'FBC'.and.j(j(it,k),k+1).lt.0)then
               z(it)=z(it)+J2*Di
               rf(it)=rf(it)+J2*Di
            endif
c     if(j(j(i,k),k+1).lt.0) z(i)=z(i)+J2*zchnn(i)
c     ---------------------------
         enddo
         if(j(j(it,1),4).gt.0.and.j(it,1).gt.0)then
            z(j(j(it,1),4))= 
     &           z(j(j(it,1),4))+J2*Di
            rf(j(j(it,1),4))= 
     &           rf(j(j(it,1),4))+J2*Di
         endif
         if(BC.eq.'FBC'.and.j(j(it,1),4).gt.0.and.j(it,1).gt.0)then
            z(it)=z(it)+J2*Di
            rf(it)=rf(it)+J2*Di
         endif
      endif
      
      if(randomtype.eq.'runif'.or.randomtype.eq.'sgGAU')then

         i=coorx+(coory-1)*L
!--- Initial heterogeneity: One dislocation loop
         Di=Dini
         dis(i)=int(Dini)
!--- Initial heterogeneity: one site in phase s=1
!         Di=q
!         s(i)=1
         if(SR_LR.eq.'LR')then
         do i1=1,NL
            z(i1)=z(i1)+DJij(i,i1)*Di
            z(i1)=z(i1)+r*(rand()-0.5d0)
            rf(i1)=rf(i1)+DJij(i,i1)*Di
         enddo
         endif

         if(SR_LR.eq.'SR')then
         do i1=0,8
            if(j(i,i1).gt.0)then
               z(j(i,i1))=z(j(i,i1))+DJSR(i1)*Di
               rf(j(i,i1))=rf(j(i,i1))+DJSR(i1)*Di
            endif
         enddo
         endif

!         i=coorx+(coory-1)*L
!         z(i)=r*(rand()-0.5d0)
         
!     z(15)=1.d0
      endif

      do i=1,NL

!         if(randomtype.eq.'runif')then
!            sgi=r*(rand()-0.5d0)
!            z(i)=sgi
!         endif

c     Initial randomness in sigma
         if(randomtype.eq.'sgGAU'.and.r.ne.0.0d0)then
 2101       continue
            zr1=rvec(ivec)
            zr2=rvec(ivec+1)
            ivec=ivec+2
            if (ivec.ge.lcarry) then
               call rcarry(rvec,lcarry)
               ivec=1
            endif
            xr= 2.0d0*zr1 -1.0d0
            yr= 2.0d0*zr2 -1.0d0
            rho_gauss=xr*xr+yr*yr
            if((rho_gauss.gt.1.0d0).or.(rho_gauss.eq.0.0d0))goto 2101
            sgi=xr*dsqrt(-2.0d0*log(rho_gauss)/rho_gauss)*r*fJ*q

            z(i)=sgi
            rf(i)=sgi
         endif


         if(randomtype.eq.'rf_Di'.and.r_Di.ne.0.0d0)then

c            goto 1212

            rv=rvec(ivec)
            if(rv.lt.raux_Di) Di=Dini
            if(rv.ge.raux_Di.and.rv.lt.r_Di) Di=-Dini
            if(rv.ge.r_Di) Di=0.0d0

            ivec=ivec+1
            if (ivec.ge.lcarry) then
               call rcarry(rvec,lcarry)
               ivec=1
            endif

            if(Di.ne.0.0d0)then

c               print*,'Index dislocation=',i
               dis(i)=int(Di)

               z(i)=z(i)-f*(J1*float(z1)+J2*float(z2))*Di
               rf(i)=rf(i)-f*(J1*float(z1)+J2*float(z2))*Di
               do k=1,z1
                  if(j(i,k).gt.0)then
                     z(j(i,k))=z(j(i,k))+J1*Di
                     rf(j(i,k))=rf(j(i,k))+J1*Di
c     density=density+zchnn(site)
                  endif
c     --- Open boundaries ---
                  if(BC.eq.'FBC'.and.j(i,k).lt.0)then
                     z(i)=z(i)+J1*Di
                     rf(i)=rf(i)+J1*Di
                  endif
c     if(j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------------------------
               enddo
               do k=1,z1-1
                  if(j(j(i,k),k+1).gt.0.and.j(i,k).gt.0)then
                     z(j(j(i,k),k+1))=z(j(j(i,k),k+1))
     &                    +J2*Di
                     rf(j(j(i,k),k+1))=rf(j(j(i,k),k+1))
     &                    +J2*Di
                  endif
c     --- Open boundaries ---
                  if(BC.eq.'FBC'.and.j(j(i,k),k+1).lt.0)then
                     z(i)=z(i)+J2*Di
                     rf(i)=rf(i)+J2*Di
                  endif
c     if(j(j(i,k),k+1).lt.0) z(i)=z(i)+J2*zchnn(i)
c     ---------------------------
               enddo
               if(j(j(i,1),4).gt.0.and.j(i,1).gt.0)then
                  z(j(j(i,1),4))= 
     &                 z(j(j(i,1),4))+J2*Di
               endif
            endif
 1212       continue

 2102       continue
            zr1=rvec(ivec)
            zr2=rvec(ivec+1)
            ivec=ivec+2
            if (ivec.ge.lcarry) then
               call rcarry(rvec,lcarry)
               ivec=1
            endif
            xr= 2.0d0*zr1 -1.0d0
            yr= 2.0d0*zr2 -1.0d0
            rho_gauss=xr*xr+yr*yr
            if((rho_gauss.gt.1.0d0).or.(rho_gauss.eq.0.0d0))goto 2102
            sgi=xr*dsqrt(-2.0d0*log(rho_gauss)/rho_gauss)*r*q

            z(i)=z(i)+sgi
         endif


         if(randomtype.eq.'sgRND'.and.r.ne.0.0d0)then
            sgi=r*(rvec(ivec)-0.5d0)
c     rv=rvec(ivec)
c     if(rv.lt.0.05) z(i)=1.0d0
c     if(rv.ge.0.05.and.rv.lt.0.1) z(i)=-1.0d0
c     if(rv.ge.0.1) z(i)=0.0d0
            
            ivec=ivec+1
            if (ivec.ge.lcarry) then
               call rcarry(rvec,lcarry)
               ivec=1
            endif

            z(i)=z(i)+sgi*(1.0-(J1*float(z1)+J2*float(z2)))

            do k=1,z1
               if(j(i,k).gt.0)then
                  z(j(i,k))=z(j(i,k))+J1*sgi
c     density=density+zchnn(site)
               endif
c     --- Open boundaries ---
c     if(j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------------------------
            enddo
            do k=1,z1-1
               if(j(j(i,k),k+1).gt.0.and.j(i,k).gt.0)then
                  z(j(j(i,k),k+1))=z(j(j(i,k),k+1))
     &                 +J2*sgi
               endif
c     --- Open boundaries ---
c     if(j(j(i,k),k+1).lt.0) z(i)=z(i)+J2*zchnn(i)
c     ---------------------------
            enddo
            if(j(j(i,1),4).gt.0.and.j(i,1).gt.0) 
     &           z(j(j(i,1),4))= 
     &           z(j(j(i,1),4))+J2*sgi
         endif

c     Initial randomness in Di
         if(randomtype.eq.'DiRND'.and.r_Di.ne.0.0d0)then
            rv=rvec(ivec)
            if(rv.lt.raux_Di) Di=Dini
            if(rv.ge.raux_Di.and.rv.lt.r_Di) Di=-Dini
            if(rv.ge.r_Di) Di=0.0d0

            ivec=ivec+1
            if (ivec.ge.lcarry) then
               call rcarry(rvec,lcarry)
               ivec=1
            endif

            if(Di.ne.0.0d0)then

c               print*,'Index dislocation=',i
               dis(i)=int(Di)

               z(i)=z(i)-f*(J1*float(z1)+J2*float(z2))*Di
               rf(i)=rf(i)-f*(J1*float(z1)+J2*float(z2))*Di
               
               do k=1,z1
                  if(j(i,k).gt.0)then
                     z(j(i,k))=z(j(i,k))+J1*Di
                     rf(j(i,k))=rf(j(i,k))+J1*Di
c     density=density+zchnn(site)
                  endif
c     --- Open boundaries ---
c     if(j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------------------------
               enddo
               do k=1,z1-1
                  if(j(j(i,k),k+1).gt.0.and.j(i,k).gt.0)then
                     z(j(j(i,k),k+1))=z(j(j(i,k),k+1))
     &                    +J2*Di
                     rf(j(j(i,k),k+1))=rf(j(j(i,k),k+1))
     &                    +J2*Di
                  endif
c     --- Open boundaries ---
c     if(j(j(i,k),k+1).lt.0) z(i)=z(i)+J2*zchnn(i)
c     ---------------------------
               enddo
               if(j(j(i,1),4).gt.0.and.j(i,1).gt.0)then
                  z(j(j(i,1),4))= 
     &                 z(j(j(i,1),4))+J2*Di
                  rf(j(j(i,1),4))= 
     &                 rf(j(j(i,1),4))+J2*Di
               endif
            endif

         endif
        
c	 if(r.eq.0.0d0) z(i)=0.0d0
!         print*,z(i)
         s(i)=0
c         dis(i)=0
c	z(i)=-0.1d0
c	z(i+1)=0.1d0
c      	z(i)=5.0d0*float(z1)
c	z(i+5)=-5.0d0*float(z1)
      enddo


      if(randomtype.eq.'DiDET')then
         
         i=coorx+(coory-1)*L
!--- Initial heterogeneity: One dislocation loop
         Di=Dini
         dis(i)=int(Dini)
!--- Initial heterogeneity: one site in phase s=1
!         Di=q
!         s(i)=1
         if(SR_LR.eq.'LR')then
         do i1=1,NL
            z(i1)=z(i1)+DJij(i,i1)*Di
            rf(i1)=rf(i1)+DJij(i,i1)*Di
!            print*,i,i1,DJij(i,i1),Di,L,z(i1)
         enddo
         endif

         if(SR_LR.eq.'SR')then
            do i1=0,8
               if(j(i,i1).gt.0)then
                  z(j(i,i1))=z(j(i,i1))+DJSR(i1)*Di
                  rf(j(i,i1))=rf(j(i,i1))+DJSR(i1)*Di
               endif
            enddo
         endif
      endif
      
!      do i=1,NL
!         write(82,*)z(i)
!      enddo
!      close(82)

!      stop
c      print*,'hola'
c      print*,'hola'

      do i=1,NL
         Hs(i)=0
         Ht(i)=0
         HsPT(i)=0
         HsD(i)=0
	 HA(i)=0
         HsCool(i)=0
         HsCoolNS(i)=0
         HsCool1D(i)=0
         HsCool2D(i)=0
         HtCool(i)=0
         HsPTCool(i)=0
         do ii=0,nb_g
            Hs_g_c(i,ii)=0
            HsNS_g_c(i,ii)=0
            Hs1D_g_c(i,ii)=0
            Hs2D_g_c(i,ii)=0
            Hs_g_h(i,ii)=0
            HsPT_g(i,ii)=0
         enddo
      enddo
      do i=NL+1,2*NL
         Hs(i)=0
         HsCool(i)=0
         HsPT(i)=0
         HsPTCool(i)=0
         do ii=0,nb_g
            Hs_g_c(i,ii)=0
            Hs_g_h(i,ii)=0
            HsPT_g(i,ii)=0
            HsNS_g_c(i,ii)=0
            Hs1D_g_c(i,ii)=0
            Hs2D_g_c(i,ii)=0
         enddo
      enddo

      if(nshapes.gt.0)then
      do k=1,nshapes
         do i=1,tshape(k)
            shapesum(i,k)=0
         enddo
         sumshape(k)=0
      enddo
      endif

      do i=0,nb_branch
         h_branch(i)=0.d0
         do ii=0,nb_g
            h_branch_g(i,ii)=0.d0
            h_branchPT_g(i,ii)=0.d0
         enddo
      enddo
      sum_branch=0.d0
      sum_branchPT=0.d0

      do i=0,nb_rf
         do ii=0,nb_g
            h_rf_g(i,ii)=0.d0
            sum_rf_g(ii)=0.d0
         enddo
      enddo

      do ii=0,nb_g
         do i=0,nb_stab
            hstab_g(i,ii,1)=0.d0
            hstab_g(i,ii,2)=0.d0
            hstab_g(i,ii,3)=0.d0
         enddo
         do i=0,nb_Deltai
            hDeltai_g(i,ii)=0.d0
         enddo
         do i=0,nb_dnflip
            hdnflip_g(i,ii)=0.d0
         enddo
         sum_stab_g(ii,1)=0.d0
         sum_stab_g(ii,2)=0.d0
         sum_stab_g(ii,3)=0.d0
         sum_Deltai_g(ii)=0.d0
         sum_dnflip_g(ii)=0.d0
         do ik=1,4
         E_rf(ii,ik)=0.d0
         E_rf_h(ii,ik)=0.d0
         E_Js_c(ii,ik)=0.d0
         enddo
         num_rf(ii)=0
         num_rf_h(ii)=0
         obs_branch_g_c(ii)=0.d0
         obs_branchPT_g_c(ii)=0.d0
         obs_branchD_g_c(ii)=0.d0
         obs_branch_g_h(ii)=0.d0
         Sum1_branch_g_c(ii)=0.d0
         Sum1_branchPT_g_c(ii)=0.d0
         Sum2_branchPT_g_c(ii)=0.d0
         Sum1_branchD_g_c(ii)=0.d0
         Sum1_branch_g_h(ii)=0.d0
      enddo
      hDeltai(0:nb_Deltai)=0.d0
      sum_Deltai=0.d0

      nflip=0
      nflip_abs=0
c     -----------------------------

      if(increase.eq.'n') ncycle=1
      
         
      if(stress_out.eq.'y')
     &     open(60,file='results/'//filename//'/'//filename
     &     //extension(igrid)//'_'//BC//'_'
     &     //c_name(i_c)//
     &     '_Stress.dat')

      if(stab_out.eq.'y')then
         open(81,file='results/'//filename//
     &        '/stability/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//
     &        '_StabilityFile_g.dat')
         write(81,*)'File  ','g'
      endif

      strain=0.0d0
      stress=0
      density=0.0d0
      ndis=0
      iaval_tot=0
      i_snap_cool=0
      i_snap_heat=0
      i_1D_span=0
      i_2D_span=0
      Ns1_glob=0
      Ns2_glob=0
      
      do 1000 icycle=1,ncycle
!         if(mod(icycle,100).eq.0)print*,'icycle=',icycle
         print*,'icycle=',icycle
! Condicional output of last cycle.
         if(icycle.eq.ncycle.and.stress_out.eq.'l')then
            stress_out='y'
            open(60,file='results/'//filename//'/'//filename
     &           //extension(igrid)//'_'//BC//'_'
     &           //c_name(i_c)//
     &           '_Stress.dat')
         endif

         if(icycle.eq.ncycle.and.Dis_cyc_out.eq.'l')then
            Dis_cyc_out='y'
            open(112,file='results/'//filename//'/'//filename
     &           //extension(igrid)//'_'//BC//'_'
     &           //c_name(i_c)//'_rf_cycles.dat')
         endif


!         print*,'icycle=',icycle
c     ------ DECREASING h (~Temperature) -------------
         iaval=0

         tglob=0
         delta_max1=-1.d0
         branch_min0=big
         branch_max0=-big
         branch_minTot=big
         branch_maxTot=-big
         branchPT_minTot=big
         branchPT_maxTot=-big

         if(randomtype.eq.'sgGAU'.and.r.ne.0.0d0.and.fD.gt.1.d0)then
            i=coorx+(coory-1)*L
!---  Initial heterogeneity: One dislocation loop
            Di=Dini
            dis(i)=int(Dini)
!---  Initial heterogeneity: one site in phase s=1
!     Di=q
!     s(i)=1
            if(SR_LR.eq.'LR')then
               do i1=1,NL
                  z(i1)=z(i1)+DJij(i,i1)*Di
                  z(i1)=z(i1)+r*(rand()-0.5d0)
                  rf(i1)=rf(i1)+DJij(i,i1)*Di
               enddo
            endif
            
            if(SR_LR.eq.'SR')then
               do i1=0,8
                  if(j(i,i1).gt.0)then
                     z(j(i,i1))=z(j(i,i1))+DJSR(i1)*Di
                     rf(j(i,i1))=rf(j(i,i1))+DJSR(i1)*Di
                  endif
               enddo
            endif

            do i=1,NL
 2111          continue
               zr1=rvec(ivec)
               zr2=rvec(ivec+1)
               ivec=ivec+2
               if (ivec.ge.lcarry) then
                  call rcarry(rvec,lcarry)
                  ivec=1
               endif
               xr= 2.0d0*zr1 -1.0d0
               yr= 2.0d0*zr2 -1.0d0
               rho_gauss=xr*xr+yr*yr
               if((rho_gauss.gt.1.0d0).or.(rho_gauss.eq.0.0d0))goto 2111
               sgi=xr*dsqrt(-2.0d0*log(rho_gauss)/rho_gauss)*r*fJ*q
               
               z(i)=sgi
               rf(i)=sgi
            enddo
         endif

      Th0(0,1)=big
      Th0(-1,1)=-big
      Th0(0,2)=-big
      Th0(1,2)=big

      irun=0
      MaxsPT=0
      nthermal=0
      avg_branch=0.d0 !mean branching ratio over avalanches in a cycle
      avg_branchPT=0.d0 !mean branching ratio for phase transition events
      avg_branchD=0.d0 !mean branching ratio over avalanches in a cycle

!Calculation of the statistics of the random slip field in austenite
      rf1=0.d0
      rf2=0.d0
      do i=1,NL
         rf1=rf1+rf(i)
         rf2=rf2+rf(i)**2.d0
      enddo

      avg_rftot=rf1/dble(NL)
      var_rftot=rf2/dble(NL)-avg_rftot**2.d0
!      print*,rf1,rf2,avg_rftot,var_rftot,fJ,qini(igrid),
!     &     var_rftot/((qini(igrid)*fJ)**2.d0)


      write(108,*)icycle,avg_rftot/(qini(igrid)*fJ),
     &     var_rftot/((qini(igrid)*fJ)**2.d0)

!      if(Dis_cyc_out.eq.'y')write(112,*)q*h/(2.d0*zeff),
!     &              avg_rftot/(qini(igrid)*fJ),
!     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
!     &     iaval_tot

      if(Dis_cyc_out.eq.'y'.and.iaval_tot.eq.0)
     &        write(112,*)-100.d0,
     &              avg_rftot/(qini(igrid)*fJ),
     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &     0,0,0,-1,cumm_f
!      stop


      do 1001
      irun=irun+1
c         print*,'irun, nflip =',irun,nflip


c         write(80,*) irun,nflip

c     - computation of h -
         h=-big
         h1=-big
         h2=-big
         h3=big
         h4=big
         do i=1,NL
            y=int(i/L)
            if(mod(i,L).eq.0) y=y-1
            x=i-y*L-1

            if(x.gt.xmin.and.x.lt.xmax)then
!               print*,x,y
            if(s(i).eq.0.and.z(i).lt.Th0(0,1).and.z(i).gt.h1)then
               h1=z(i)
               haux=-fDTemp*q*(zeff+zabs)+h1
               if(haux.gt.h) h=haux
            endif

            if(s(i).eq.-1.and.z(i).gt.Th0(-1,1).and.z(i).gt.h2)then
               h2=z(i)
               haux=fDTemp*q*(zeff+zabs)-h2
               if(haux.gt.h) h=haux
            endif

            if(s(i).eq.0.and.z(i).gt.Th0(0,2).and.z(i).lt.h3)then
               h3=z(i)
               haux=-fDTemp*q*(zeff+zabs)-h3
               if(haux.gt.h) h=haux
            endif

            if(s(i).eq.1.and.z(i).lt.Th0(1,2).and.z(i).lt.h4)then
               h4=z(i)
               haux=fDTemp*q*(zeff+zabs)+h4
               if(haux.gt.h) h=haux
            endif

            endif !endif for excluding vertical boundaries
         enddo

c         print*,'hi=',h1,h2,h3,h4,z(1)

c         print*,'hi=',h1,h2,h3,h4,z(1)

!         print*,'h=',h,'zmax=',-q*fDTemp*float(z1),'h1=',h1

c         stop
	 if(h.eq.-big)then
!            print*,'hola'
            goto 2001
         endif

	 
         h=h-feed
         Th0(0,1)=fDTemp*q*(zeff+zabs)+h
         Th0(1,2)=Th0(0,1)-2.0d0*fDTemp*q*(zeff+zabs)
         Th0(0,2)=-Th0(0,1)
         Th0(-1,1)=-Th0(0,1)+2.0d0*fDTemp*q*(zeff+zabs)
         g_new=q*h/(2.d0*(zeff+zabs))

         iaval=iaval+1
         iaval_tot=iaval_tot+1
         if(iaval_tot.ge.2)aval_dens=dble(size)/abs(g_new-g_old)
!         if(stress_out.eq.'y'.and.icycle.ge.limav)

      cumm_f=float(nflip_abs)/float(NL)

         if(stress_out.eq.'y'.and.iaval_tot.eq.1)
     &              write(60,*) q*h/(2.d0*(zeff+zabs)),strain/float(NL),
     &              float(nflip)/float(NL),float(ndis)/float(NL),icycle,
     &     iaval_tot,-1,cumm_f


         if(stress_out.eq.'y'.and.iaval_tot.ge.2)
     &              write(60,*) q*h/(2.d0*(zeff+zabs)),strain/float(NL),
     &              float(nflip)/float(NL),float(ndis)/float(NL),icycle,
     &     iaval_tot,-1,cumm_f

      if(Dis_cyc_out.eq.'y'.and.iaval_tot.eq.1)
     &        write(112,*)q*h/(2.d0*(zeff+zabs)),
     &              avg_rftot/(qini(igrid)*fJ),
     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &        iaval_tot,0,0,-1,cumm_f

         if(Dis_cyc_out.eq.'y'.and.iaval_tot.ge.2)
     &        write(112,*)q*h/(2.d0*(zeff+zabs)),
     &              avg_rftot/(qini(igrid)*fJ),
     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &     iaval_tot,aval_dens,size,-1,cumm_f

         g_old=q*h/(2.d0*(zeff+zabs))

         tau_aux1=g_old/(qini(igrid)**2.d0)
         dif_aux=abs(tau_aux1-tau_rf_out)
!         print*,g_old,tau_aux1,dif_aux,aa_g
         if(rf_raw_out.eq.'y'
     &        .and.dif_aux.lt.aa_g
     &        .and.icycle.ge.nc_min
     &        .and.icycle.le.nc_max)then
            do i1=1,NL
               write(125,*)g_old/(q**2.d0),rf(i1)/(qini(igrid)*fJ)
            enddo
         endif
         
         daux=(g_old-hmin_g2)/aa_g2
         dif_aux=daux-dble(int(daux))
!         print*,daux,dif_aux,tol_rf
         xh2=hmin_g2+aa_g2*dble(int(daux))
         if(rf_raw_out.eq.'y'
     &        .and.dif_aux.lt.tol_rf
     &        .and.hmin_g2.le.xh2
     &        .and.hmax_g2.ge.xh2
     &        .and.icycle.ge.nc_min
     &        .and.icycle.le.nc_max)then
            do i1=1,NL
               write(126,*)
     &              xh2/(q**2.d0),
     &              rf(i1)/(qini(igrid)*fJ)
            enddo
         endif

         if(temp_out.eq.'y')then
         if(icycle.ge.limav)then
            gaux=q*h/(2.d0*(zeff+zabs))
            ig=int((gaux-hmin_g)/aa_g)
            if(ig.gt.nb_g)ig=nb_g
            if(ig.lt.0)ig=0

            dnflip_aux=dble(nflip)/dble(NL)
            idnflip=int((dnflip_aux-dnflip_min)/aa_dnflip)
            if(idnflip.gt.nb_dnflip)idnflip=nb_dnflip

            hdnflip_g(idnflip,ig)=hdnflip_g(idnflip,ig)+1.d0
            sum_dnflip_g(ig)=sum_dnflip_g(ig)+1.d0
         endif
         if(icycle.ge.limav)then            
            do i=1,NL
               gaux=q*h/(2.d0*(zeff+zabs))
               ig=int((gaux-hmin_g)/aa_g)
               if(ig.gt.nb_g)ig=nb_g
               if(ig.lt.0)ig=0

               E_rf(ig,1)=E_rf(ig,1)+rf(i)
               E_rf(ig,2)=E_rf(ig,2)+rf(i)**2.d0
               E_rf(ig,3)=E_rf(ig,3)+rf(i)**3.d0
               E_rf(ig,4)=E_rf(ig,4)+rf(i)**4.d0
               num_rf(ig)=num_rf(ig)+1

               dJs=z(i)-rf(i)
               E_Js_c(ig,1)=E_Js_c(ig,1)+dJs
               E_Js_c(ig,2)=E_Js_c(ig,2)+dJs**2.d0
               E_Js_c(ig,3)=E_Js_c(ig,3)+dJs**3.d0
               E_Js_c(ig,4)=E_Js_c(ig,4)+dJs**4.d0

               irf=int((rf(i)-rf_min)/aa_rf)
               if(irf.lt.0)irf=0
               if(irf.gt.nb_rf)irf=nb_rf
               h_rf_g(irf,ig)=h_rf_g(irf,ig)+1.d0
               sum_rf_g(ig)=sum_rf_g(ig)+1.d0

               if(stab_out.eq.'y')then
               if(s(i).eq.0)then
                  s1=(Th0(0,1)-abs(z(i)))/(2.d0*(zeff+zabs))
                  istab=int((s1-ddmin)/aa)
!     print*,'s1,istab=',s1,istab
                  if(istab.gt.nb_stab)istab=nb_stab
                  if(istab.lt.0)istab=0
                  hstab_g(istab,ig,1)=hstab_g(istab,ig,1)+1.d0
                  sum_stab_g(ig,1)=sum_stab_g(ig,1)+1.d0
               endif
               if(s(i).eq.1)then
                  s2=(z(i)-Th0(1,2))/(2.d0*(zeff+zabs))
                  istab=int((s2-ddmin)/aa)
                  if(istab.gt.nb_stab)istab=nb_stab
                  if(istab.lt.0)istab=0
                  hstab_g(istab,ig,2)=hstab_g(istab,ig,2)+1.d0

                  s3=(Th0(1,1)-z(i))/(2.d0*(zeff+zabs))
                  istab=int((s3-ddmin)/aa)
                  if(istab.gt.nb_stab)istab=nb_stab
                  if(istab.lt.0)istab=0
                  hstab_g(istab,ig,3)=hstab_g(istab,ig,3)+1.d0

                  sum_stab_g(ig,2)=sum_stab_g(ig,2)+1.d0
                  sum_stab_g(ig,3)=sum_stab_g(ig,3)+1.d0
!                        print*,'s2,s3,istab=',s2,s3,istab
               endif
               if(s(i).eq.-1)then
                  s2=(Th0(-1,1)-z(i))/(2.d0*(zeff+zabs))
                  istab=int((s2-ddmin)/aa)
                  if(istab.gt.nb_stab)istab=nb_stab
                  if(istab.lt.0)istab=0
                  hstab_g(istab,ig,2)=hstab_g(istab,ig,2)+1.d0

                  s3=(z(i)-Th0(-1,2))/(2.d0*(zeff+zabs))
                  istab=int((s3-ddmin)/aa)
                  if(istab.gt.nb_stab)istab=nb_stab
                  if(istab.lt.0)istab=0
                  hstab_g(istab,ig,3)=hstab_g(istab,ig,3)+1.d0
                  sum_stab_g(ig,2)=sum_stab_g(ig,2)+1.d0
                  sum_stab_g(ig,3)=sum_stab_g(ig,3)+1.d0
               endif
               endif
            enddo
         endif
         endif !endif for temp_out

!--------- Stability pdfs just before the avalanche --------
               if(stab_out.eq.'y'.and.icycle.eq.ncycle)then
                  aux = iaval
!                  ts = char(ichar('0') + aux/10000)
!                  aux = mod(aux,10000)
                  ms = char(ichar('0') + aux/1000)
                  aux = mod(aux,1000)
                  cs = char(ichar('0') + aux/100)
                  aux = mod(aux,100)
                  ds = char(ichar('0') + aux/10)
                  aux = mod(aux,10)
                  us = char(ichar('0') + aux)
                  Buffer=ms//cs//ds//us

                  open(74,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_sigma0.dat')
                  open(75,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_sigma1.dat')
                  open(76,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_sigmaS.dat')
                  open(78,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_Pdelta_s0.dat')
                  open(79,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_Pdelta_s+1.dat')
                  open(80,file='results/'//filename//
     &                 '/stability/'//filename
     &     //extension(igrid)//'_'//BC//'_'//Buffer//
     &                 '_'//c_name(i_c)//'_Pdelta_s-1.dat')


                  write(81,*)Buffer,q*h/(2.d0*(zeff+zabs))

!                  print*,'aa=',aa
                  do i=0,nb_stab
                     hstab(i,1)=0.d0
                     hstab(i,2)=0.d0
                     hstab(i,3)=0.d0
                     hdstab(i,1)=0.d0
                     hdstab(i,2)=0.d0
                     hdstab(i,3)=0.d0
                  enddo
                  sum_stab(1)=0.d0
                  sum_stab(2)=0.d0
                  sum_stab(3)=0.d0
                  sum_dstab(1)=0.d0
                  sum_dstab(2)=0.d0
                  sum_dstab(3)=0.d0

                  do i=1,NL
                     if(s(i).eq.0)then
                        s1=(Th0(0,1)-abs(z(i)))/(2.d0*(zeff+zabs))
                        istab=int((s1-ddmin)/aa)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
!                        print*,'s1,istab=',s1,istab
                        if(istab.gt.nb_stab)print*,istab,s1,aa
                        hstab(istab,1)=hstab(istab,1)+1.d0
                        sum_stab(1)=sum_stab(1)+1.d0

                        s1=z(i)/(2.d0*(zeff+zabs))
                        istab=int((s1-ddeltamin)/ad)
                        hdstab(istab,1)=hdstab(istab,1)+1.d0
                        sum_dstab(1)=sum_dstab(1)+1.d0
                     endif
                     if(s(i).eq.1)then
                        s2=(z(i)-Th0(1,2))/(2.d0*(zeff+zabs))
                        istab=int((s2-ddmin)/aa)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hstab(istab,2)=hstab(istab,2)+1.d0
                        s3=(Th0(1,1)-z(i))/(2.d0*(zeff+zabs))
!                  if(iaval.eq.10)
!                        write(77,*)s(i),s2,z(i),Th0(1,2)
!                        write(77,*)s3,z(i),Th0(1,1),Th0(1,2)
                        istab=int((s3-ddmin)/aa)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hstab(istab,3)=hstab(istab,3)+1.d0
                        sum_stab(2)=sum_stab(2)+1.d0
                        sum_stab(3)=sum_stab(3)+1.d0
!                        print*,'s2,s3,istab=',s2,s3,istab

                        s1=z(i)/(2.d0*(zeff+zabs))
                        istab=int((s1-ddeltamin)/ad)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hdstab(istab,2)=hdstab(istab,2)+1.d0
                        sum_dstab(2)=sum_dstab(2)+1.d0
                     endif
                     if(s(i).eq.-1)then
                        s2=(Th0(-1,1)-z(i))/(2.d0*(zeff+zabs))
                        istab=int((s2-ddmin)/aa)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hstab(istab,2)=hstab(istab,2)+1.d0
                        s3=(z(i)-Th0(-1,2))/(2.d0*(zeff+zabs))
!                        write(77,*)s(i),s2,z(i),Th0(-1,1)
!                  if(iaval.eq.10)write(77,*)z(i),Th0(1,1),Th0(1,2)
                        istab=int((s3-ddmin)/aa)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hstab(istab,3)=hstab(istab,3)+1.d0
                        sum_stab(2)=sum_stab(2)+1.d0
                        sum_stab(3)=sum_stab(3)+1.d0
!                        print*,'s2,s3,istab,aa=',s2,s3,istab,aa

                        s1=z(i)/(2.d0*(zeff+zabs))
                        istab=int((s1-ddeltamin)/ad)
                        if(istab.gt.nb_stab)istab=nb_stab
                        if(istab.lt.0)istab=0
                        hdstab(istab,3)=hdstab(istab,3)+1.d0
                        sum_dstab(3)=sum_dstab(3)+1.d0
!                        print*,s1,z(i),istab,ad
                     endif
!                     print*,'istab,z(i),s(i),aa=',istab,z(i),s(i),aa
                  enddo
                  do i=0,nb_stab
                     delta=ddmin+aa*dble(i)
                     if(hstab(i,1).gt.0.d0)write(74,*)delta,
     &                    hstab(i,1)/sum_stab(1),
     &                    q*h/(2.d0*(zeff+zabs))
                     if(hstab(i,2).gt.0.d0)write(75,*)delta,
     &                    hstab(i,2)/sum_stab(2),
     &                    q*h/(2.d0*(zeff+zabs))
                     if(hstab(i,3).gt.0.d0)write(76,*)delta,
     &                    hstab(i,3)/sum_stab(3),
     &                    q*h/(2.d0*(zeff+zabs))
                  enddo
                  do i=0,nb_delta
                     ddelta=ddeltamin+ad*dble(i)
                     if(hdstab(i,1).gt.0.d0)write(78,*)ddelta,
     &                    hdstab(i,1)/sum_dstab(1),
     &                    q*h/(2.d0*(zeff+zabs))
                     if(hdstab(i,2).gt.0.d0)write(79,*)ddelta,
     &                    hdstab(i,2)/sum_dstab(2),
     &                    q*h/(2.d0*(zeff+zabs))
                     if(hdstab(i,3).gt.0.d0)write(80,*)ddelta,
     &                    hdstab(i,3)/sum_dstab(3),
     &                    q*h/(2.d0*(zeff+zabs))
                  enddo
                  close(74)
                  close(75)
                  close(76)
                  close(78)
                  close(79)
                  close(80)
               endif
!--------------------------------------------------------------



!	 print*,Th0(0,1),Th0(1,2),Th0(0,2),Th0(-1,1),Th0(-1,2),Th0(1,1)


c         ProxE=1
c        do i=1,NL
c            if(z(i).gt.Th0(s(i),1).or.
c     &           z(i).lt.Th0(s(i),2))then
c               Qe(ProxE)=i
cc     print*,j(nini,i)
c               ProxE=ProxE+1
cc               goto 1002
cc               print*,i
c            endif
c         enddo
c
c	 print*,'ProxE=',ProxE,Qe(1)
c	 if(ProxE.eq.1)goto 2001
	 
c     ------ Propagating the avalanche ----
 1002    stress=stress+1
         size=0
         t=0
         sPT=0
         tPT=0
         sD=0
         ProxG=1
         Fincapa=ProxE
         nshell=0
	 maxA=0
c         print*,'hola'

         nava=0
         bar_branch=0.d0
         bar_branchPT=0.d0
         bar_branchD=0.d0
!     Coordinates for the shadow window 
        xSmin=L+1
        xSmax=0
        ySmin=L
        ySmax=0
        sox(:)=0
        soy(:)=0

         do 1003
            nava=nava+1

	    nuns=0
            nunsPT=0
            nup=0
            ndw=0

!            if(size.gt.NL)print*,'size=',size

	    do i=1,NL
c               if(z(i).gt.Th0(s(i),1).or.
c     &              z(i).lt.Th0(s(i),2))then
c		  nuns=nuns+1
c		  uns(nuns)=i
               site=i
                  y=int(i/L)
                  if(mod(i,L).eq.0) y=y-1
                  x=i-y*L-1

                  if(x.gt.xmin.and.x.lt.xmax)then
c     casos indecisos
c                  if(z(site).gt.Th0(s(site),1).and.
c     &                 z(site).lt.Th0(s(site),2).and.
c     &                 s(site).eq.0.and.z(site).eq.0.0d0)then
c                  if(z(site).gt.Th0(s(site),1).and.
c     &                 z(site).lt.Th0(s(site),2).and.
c     &                 z(site).eq.0.0d0.and.s(site).eq.0)then
                  if(z(site).gt.Th0(s(site),1).and.
     &                 z(site).lt.Th0(s(site),2))then
                     nuns=nuns+1
                     nunsPT=nunsPT+1
                     uns(nuns)=site

                     if(z(site).eq.0.0d0)then
                        nthermal=nthermal+1
                        ivec=ivec+1
                        if (ivec.ge.lcarry) then
                           call rcarry(rvec,lcarry)
                           ivec=1
                        endif
                        if(rvec(ivec).le.0.5)then
                           size=size+1
                           strain=strain+zch(s(site),1)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.-1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.1)then
                              ndis=ndis+1
                              dis(site)=dis(site)+1
                              nunsPT=nunsPT-1
                           endif
                           
!                           z(site)=z(site)
!     &                          +DJij(site,site)*zch(s(site),1)
                           zchnn(site)=zch(s(site),1)
                           s(site)=sch(s(site),1)
                           density=density-float(z1)*zch(s(site),1)
                           ndw=ndw+1
                        else
                           size=size+1
                           strain=strain+zch(s(site),2)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.-1)then
                              ndis=ndis-1
                              dis(site)=dis(site)-1
                              nunsPT=nunsPT-1
                           endif
                           
!                           z(site)=z(site)
!     &                          +DJij(site,site)*zch(s(site),2)
                           zchnn(site)=zch(s(site),2)
                           s(site)=sch(s(site),2)
c     density=density-float(z1)*zch(s(site),2)
                           nup=nup+1
                        endif
                        
                     else
                        if(z(site).gt.0)then
                           size=size+1
                           strain=strain+zch(s(site),1)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.-1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.1)then
                              ndis=ndis+1
                              dis(site)=dis(site)+1
                              nunsPT=nunsPT-1
                           endif
                           
!                           z(site)=z(site)
!     &                          +DJij(site,site)*zch(s(site),1)
                           zchnn(site)=zch(s(site),1)
                           s(site)=sch(s(site),1)
                           density=density-float(z1)*zch(s(site),1)
                           ndw=ndw+1
                        else
                           if(z(site).lt.0)then
                              size=size+1
                              strain=strain+zch(s(site),2)
                              
                              if(s(site).eq.0)then
                                 sPT=sPT+1
                                 nflip=nflip+1
                                 nflip_abs=nflip_abs+1
                              endif
                              
                              if(s(site).eq.1)then
                                 sPT=sPT+1
                                 nflip=nflip-1
                                 nflip_abs=nflip_abs-1
                              endif
                              
                              if(s(site).eq.-1)then
                                 ndis=ndis-1
                                 dis(site)=dis(site)-1
                                 nunsPT=nunsPT-1
                              endif
                              
!                              z(site)=z(site)
!     &                             +DJij(site,site)*zch(s(site),2)
                              zchnn(site)=zch(s(site),2)
                              s(site)=sch(s(site),2)
c     density=density-float(z1)*zch(s(site),2)
                              nup=nup+1
                           endif
                        endif
                     endif
                     goto 1430
                  endif
                  
                  if(z(site).gt.Th0(s(site),1))then
                     size=size+1
                     strain=strain+zch(s(site),1)
                     
                     if(s(site).eq.0)then
                        sPT=sPT+1
                        nflip=nflip+1
                        nflip_abs=nflip_abs+1
                     endif
                     
                     if(s(site).eq.-1)then
                        sPT=sPT+1
                        nflip=nflip-1
                        nflip_abs=nflip_abs-1
                     endif
                     
                     if(s(site).eq.1)then
                        ndis=ndis+1
                        dis(site)=dis(site)+1
                        nunsPT=nunsPT-1
c                        print*,'hola'
                     endif
                     
!                     z(site)=z(site)
!     &                    +DJij(site,site)*zch(s(site),1)
                     zchnn(site)=zch(s(site),1)
                     s(site)=sch(s(site),1)
c                     density=density-float(z1)*zch(s(site),1)
                     nuns=nuns+1
                     nunsPT=nunsPT+1
                     uns(nuns)=site
                  else
                     
                     if(z(site).lt.Th0(s(site),2))then
                        size=size+1
                        strain=strain+zch(s(site),2)
                        
                        if(s(site).eq.0)then
                           sPT=sPT+1
                           nflip=nflip+1
                           nflip_abs=nflip_abs+1
                        endif
                        
                        if(s(site).eq.1)then
                           sPT=sPT+1
                           nflip=nflip-1
                           nflip_abs=nflip_abs-1
                        endif
                        
                        if(s(site).eq.-1)then
                           ndis=ndis-1
                           dis(site)=dis(site)-1
                           nunsPT=nunsPT-1
c                          print*,'hola2'
                        endif
                        
!                        z(site)=z(site)
!     &                       +DJij(site,site)*zch(s(site),2)
                        zchnn(site)=zch(s(site),2)
                        s(site)=sch(s(site),2)
c                        density=density-float(z1)*zch(s(site),2)

                        nuns=nuns+1
                        nunsPT=nunsPT+1
                        uns(nuns)=site
                     endif 
                  endif

                  endif ! endif for excluding vertical boundaries

c               endif
 1430          continue
               
	    enddo
            
c$$$            print*,'nuns,ndis,nflip,icycle cooling=',
c$$$     &           nuns,ndis,nflip,icycle   

            if(nuns.gt.NL)print*,'Problem: nuns > NL'
c            stop
!            if(nflip.eq.NL)nuns=0

            if(nuns.eq.0.and.nava.eq.1)goto 2001

            if(icycle.gt.limav2)Deltai(1:NL)=0.d0

            do k=1,nuns
               site=uns(k)
!               print*,zchnn(site)
               rf1=0.d0
               rf2=0.d0
               if(SR_LR.eq.'LR')then
               do i=1,NL
                     z(i)=z(i)+DJij(i,site)*zchnn(site)
                     if(zchnn(site).eq.zch(-1,2))rf(i)=rf(i)
     &                    -DJij(i,site)
                     if(zchnn(site).eq.zch(1,1))rf(i)=rf(i)
     &                    +DJij(i,site)

                     rf1=rf1+rf(i)
                     rf2=rf2+rf(i)**2.d0

                     if(icycle.gt.limav2)Deltai(i)=Deltai(i)+
     &                    DJij(i,site)*zchnn(site)
               enddo
               endif

               if(SR_LR.eq.'SR')then
                  do i1=0,8
                     i=j(site,i1)
                     if(j(site,i1).gt.0)then
                        z(i)=z(i)+DJSR(i1)*zchnn(site)
                        if(zchnn(site).eq.zch(-1,2))rf(i)=rf(i)
     &                       -DJSR(i1)
                        if(zchnn(site).eq.zch(1,1))rf(i)=rf(i)
     &                       +DJSR(i1)
                        
                        rf1=rf1+rf(i)
                        rf2=rf2+rf(i)**2.d0
                        if(icycle.gt.limav2)Deltai(i)=Deltai(i)+
     &                       DJij(i,site)*zchnn(site)
                     endif
                  enddo
               endif

               y=int(site/L)
               if(mod(site,L).eq.0) y=y-1
               x=site-y*L-1
               sox(x)=1
               soy(y)=1

               if(x.lt.xSmin)xSmin=x
               if(x.gt.xSmax)xSmax=x
               if(y.lt.ySmin)ySmin=y
               if(y.gt.ySmax)ySmax=y

               if(sdout.eq.'y'.and.icycle.eq.ncycle)then 
                  write(73,*) x,y,tglob,iaval,s(site),dis(site),
     &                 q*h/(2.d0*(zeff+zabs)),rf(site)/(qini(igrid)*fJ),
     &                 z(site)/q,0.5d0/(1.d0+dble(s(site))*z(site)/q),
     &                 z(site)/q-rf(site)/q
               endif
               abs_z1=abs(z(site))/q
               if(abs_z1.gt.delta_max1)delta_max1=abs_z1
            enddo
            
!            e_s=0.5d0/(1.d0+abs_z)
!            if(icycle.eq.1)print*,'Maximum delta =',abs_z
!            if(icycle.eq.1)print*,'e_s =',e_s

            if(icycle.gt.limav2)then
               gaux=q*h/(2.d0*(zeff+zabs))
               ig=int((gaux-hmin_g)/aa_g)
               if(ig.gt.nb_g)ig=nb_g
               if(ig.lt.0)ig=0

               do i=1,NL
                  Deltai_aux=Deltai(i)
                  iDeltai=int((Deltai_aux-Deltai_min)/aa_Deltai)
                  if(iDeltai.gt.nb_Deltai)iDeltai=nb_Deltai
                  if(temp_out.eq.'y')then
                  hDeltai_g(iDeltai,ig)=hDeltai_g(iDeltai,ig)+1.d0
                  sum_Deltai_g(ig)=sum_Deltai_g(ig)+1.d0
                  endif
                  hDeltai(iDeltai)=hDeltai(iDeltai)+1.d0
                  sum_Deltai=sum_Deltai+1.d0
!                  write(100,*)Deltai(i),q*h/(2.d0*(zeff+zabs))
               enddo

            endif
            
            nshell=nuns
            if(nshell.ne.0)then
               t=t+1
               tglob=tglob+1
!               if(t.le.nTmax.and.icycle.ge.limav)then
               if(t.le.nTmax)then
                  shape(t)=nshell
                  shapePT(t)=nunsPT
                  shapeD(t)=nshell-nunsPT
               endif
	       if(nshell.gt.maxA) maxA=nshell
               ! Double black in the file 73 for plotting time evolution               
               if(sdout.eq.'y'.and.icycle.eq.ncycle)then
                  write(73,*)
                  write(73,*)
               endif
               if(t.gt.1.and.shape(t).gt.0)then
                  bar_branch=bar_branch+
     &                 dble(shape(t))/dble(shape(t-1))
                  branch_loc=dble(shape(t))/dble(shape(t-1))
                  gaux=q*h/(2.d0*(zeff+zabs))
                  if(branch_loc.lt.branch_minTot)then
                     branch_minTot=branch_loc
                     g_branch_minTot=gaux
                  endif
                  if(branch_loc.gt.branch_maxTot)then
                     branch_maxTot=branch_loc
                     g_branch_maxTot=gaux
                  endif
               endif
               if(t.gt.1.and.shapePT(t).gt.0.and.shapePT(t-1).gt.0)
     &              bar_branchPT=bar_branchPT+
     &              dble(shapePT(t))/dble(shapePT(t-1))
               if(t.gt.1.and.shapeD(t).gt.0.and.shapeD(t-1).gt.0)
     &              bar_branchD=bar_branchD+
     &              dble(shapeD(t))/dble(shapeD(t-1))
               if(icycle.gt.100.and.icycle.lt.110.
     &              and.shape(t-1).gt.0)
     &              write(121,*)q*h/(2.d0*(zeff+zabs)),
     &              dble(shape(t))/dble(shape(t-1))
!               print*,'bar_branch=',bar_branch,shape(t-1),shape(t),
!     &              nshell

               
            endif
            
!            if(nflip.eq.NL) goto 1004
            
            
            if(nshell.eq.0)then  ! end of the avalanche
               if(size.gt.size_snap)i_snap_cool=1
               if(size.le.2*NL.and.icycle.ge.limav)then
                  Hs(size)=Hs(size)+1
                  HsCool(size)=HsCool(size)+1

                  if(temp_out.eq.'y')then
                     gaux=q*h/(2.d0*(zeff+zabs))
                     ig=int((gaux-hmin_g)/aa_g)
                     if(ig.gt.nb_g)ig=nb_g
                     if(ig.lt.0)ig=0
                     Hs_g_c(size,ig)=Hs_g_c(size,ig)+1
                  endif

               endif
               if(t.le.NL.and.icycle.ge.limav)then
                  Ht(t)=Ht(t)+1
                  HtCool(t)=HtCool(t)+1
               endif

               if(sPT.gt.2*NL)sPT=2*NL
               if(icycle.ge.limav)then
                  HsPT(sPT)=HsPT(sPT)+1                  
                  HsPTCool(sPT)=HsPTCool(sPT)+1

                  if(temp_out.eq.'y')then
                     gaux=q*h/(2.d0*(zeff+zabs))
                     ig=int((gaux-hmin_g)/aa_g)
                     if(ig.gt.nb_g)ig=nb_g
                     if(ig.lt.0)ig=0
                     HsPT_g(sPT,ig)=HsPT_g(sPT,ig)+1
                  endif
               endif
	       if(maxA.le.NL.and.icycle.ge.limav) HA(maxA)=HA(maxA)+1
               
	       if(sPT.gt.MaxsPT) MaxsPT=sPT

! Output of hysteresis cycles with the published driving parameter: g = q*h/(4.d0*J1*z1)
!	       if(stress_out.eq.'y'.and.icycle.ge.limav)
               cumm_f=float(nflip_abs)/float(NL)

               if(stress_out.eq.'y'.and.iaval_tot.ge.2)
     &              write(60,*) q*h/(2.d0*(zeff+zabs)),strain/float(NL),
     &              float(nflip)/float(NL),float(ndis)/float(NL),icycle,
     &     iaval_tot,-1,cumm_f

!Output of Disorder after avalanche
               avg_rftot=rf1/dble(NL)
               var_rftot=rf2/dble(NL)-avg_rftot**2.d0
               
               if(Dis_cyc_out.eq.'y'.and.iaval_tot.ge.2)
     &              write(112,*)q*h/(2.d0*(zeff+zabs)),
     &              avg_rftot/(qini(igrid)*fJ),
     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &              iaval_tot,aval_dens,size,-1,cumm_f
     
	       if(mapST.eq.'y') write(91,*) sPT,size,t

               if(t.gt.1)then
!     Only phase transition counted
                  duration=dble(t)
               avg_branch=avg_branch+bar_branch/duration
               avg_branchPT=avg_branchPT+bar_branchPT/duration
               avg_branchD=avg_branchD+bar_branchD/duration

!                bar_branch=bar_branchPT
!               write(83,*)icycle,iaval,bar_branch/dble(t),
!     &              q*h/(2.d0*(zeff+zabs))
               auxbranch=bar_branch/duration
               auxbranchPT=bar_branchPT/duration
               gaux=q*h/(2.d0*(zeff+zabs))
               if(auxbranch.gt.branch_max0)then
                  branch_max0=auxbranch
                  g_branch_max0=gaux
               endif
               if(auxbranch.lt.branch_min0)then
                  branch_min0=auxbranch
                  g_branch_min0=gaux
               endif

               if(icycle.ge.limav)then
!                  if(bar_branch.gt.0.d0)then
                  duration=dble(t)
                  baux=bar_branch/duration
                  ibranch=int((baux-hmin_branch)/aa_branch)
                  gaux=q*h/(2.d0*(zeff+zabs))
                  ig=int((gaux-hmin_g)/aa_g)
                  if(ibranch.gt.nb_branch)ibranch=nb_branch
                  if(ibranch.lt.0)ibranch=0
                  if(ig.gt.nb_g)ig=nb_g
                  if(ig.lt.0)ig=0
                  h_branch(ibranch)=h_branch(ibranch)+1.d0
                  if(temp_out.eq.'y')then
                  h_branch_g(ibranch,ig)=
     &                 h_branch_g(ibranch,ig)+1.d0
                  sum_branch=sum_branch+1.d0
                  obs_branch_g_c(ig)=obs_branch_g_c(ig)+1.d0
                  Sum1_branch_g_c(ig)=Sum1_branch_g_c(ig)+baux
!                  if(icycle.lt.20)write(121,*)gaux,baux
                  endif
!                  endif
               endif
               if(icycle.ge.limav)then
!                  if(bar_branchPT.gt.0.d0)then
                  duration=dble(t)
                  baux=bar_branchPT/duration
                  ibranch=int((baux-hmin_branch)/aa_branch)
                  gaux=q*h/(2.d0*(zeff+zabs))
                  ig=int((gaux-hmin_g)/aa_g)
                  if(ibranch.gt.nb_branch)ibranch=nb_branch
                  if(ibranch.lt.0)ibranch=0
                  if(ig.gt.nb_g)ig=nb_g
                  if(ig.lt.0)ig=0
!                  h_branch(ibranch)=h_branch(ibranch)+1.d0
                  if(temp_out.eq.'y')then
                  h_branchPT_g(ibranch,ig)=
     &                 h_branchPT_g(ibranch,ig)+1.d0
                  sum_branchPT=sum_branchPT+1.d0
                  obs_branchPT_g_c(ig)=obs_branchPT_g_c(ig)+1.d0
                  Sum1_branchPT_g_c(ig)=Sum1_branchPT_g_c(ig)+baux
                  Sum2_branchPT_g_c(ig)=Sum2_branchPT_g_c(ig)+baux**2.d0
!                  if(icycle.lt.20)write(121,*)gaux,baux
                  endif
!                  endif
               endif
               if(icycle.ge.limav)then
!                  if(bar_branchPT.gt.0.d0)then
                  duration=dble(t)
                  baux=bar_branchD/duration
                  ibranch=int((baux-hmin_branch)/aa_branch)
                  gaux=q*h/(2.d0*(zeff+zabs))
                  ig=int((gaux-hmin_g)/aa_g)
                  if(ibranch.gt.nb_branch)ibranch=nb_branch
                  if(ibranch.lt.0)ibranch=0
                  if(ig.gt.nb_g)ig=nb_g
                  if(ig.lt.0)ig=0
!                  h_branch(ibranch)=h_branch(ibranch)+1.d0
                  if(temp_out.eq.'y')then
                  h_branchD_g(ibranch,ig)=
     &                 h_branchD_g(ibranch,ig)+1.d0
                  sum_branchD=sum_branchD+1.d0
                  obs_branchD_g_c(ig)=obs_branchD_g_c(ig)+1.d0
                  Sum1_branchD_g_c(ig)=Sum1_branchD_g_c(ig)+baux
!                  if(icycle.lt.20)write(121,*)gaux,baux
                  endif
!                  endif
               endif
               endif

! Check whether the cluster is spanning or not
               isox=1
	       isoy=1
               do i=xWmin,xWmax
                  isox=isox*sox(i)
               enddo
               do i=yWmin,yWmax
                  isoy=isoy*soy(i)
               enddo

!	       if(xSmin.le.xWmin.and.xSmax.ge.xWmax)isox=1
!	       if(ySmin.le.yWmin.and.ySmax.ge.yWmax)isoy=1
c$$$               print*,'xWmin,xSmin=',xWmin,xSmin
c$$$               print*,'xWmax,xSmax=',xWmax,xSmax
c$$$               print*,'yWmin,ySmin=',yWmin,ySmin
c$$$               print*,'yWmax,ySmax=',yWmax,ySmax
c$$$               print*,'isox,isoy=',isox,isoy
               if(isox.eq.1.or.isoy.eq.1)then
		  if(isox.eq.1.and.isoy.eq.1)then
		     Ns2_hi=Ns2_hi+1
                     Ns2_glob=Ns2_glob+1
		     size2_hi=size2_hi+dfloat(size)
		     ntype=2
                     i_2D_span=1
                      if(size.le.2*NL.and.icycle.ge.limav)then
                         HsCool2D(size)=HsCool2D(size)+1
c$$$                         print*,'size,HsCool2D(size)=',size,
c$$$     &                        HsCool2D(size)
                         if(temp_out.eq.'y')then
                            gaux=q*h/(2.d0*(zeff+zabs))
                            ig=int((gaux-hmin_g)/aa_g)
                            if(ig.gt.nb_g)ig=nb_g
                            if(ig.lt.0)ig=0
                            Hs2D_g_c(size,ig)=Hs2D_g_c(size,ig)+1
                         endif
                      endif
		  else
		     Ns1_hi=Ns1_hi+1
                     Ns1_glob=Ns1_glob+1
		     size1_hi=size1_hi+dfloat(size)
		     ntype=1
                     i_1D_span=1
                     if(size.le.2*NL.and.icycle.ge.limav)then
                        HsCool1D(size)=HsCool1D(size)+1
                        if(temp_out.eq.'y')then
                            gaux=q*h/(2.d0*(zeff+zabs))
                            ig=int((gaux-hmin_g)/aa_g)
                            if(ig.gt.nb_g)ig=nb_g
                            if(ig.lt.0)ig=0
                            Hs1D_g_c(size,ig)=Hs1D_g_c(size,ig)+1
                         endif
                     endif
		  endif
	       else
		  Nns_hi=Nns_hi+1.0
		  sizens_hi=sizens_hi+dfloat(size)
		  ntype=0
                  if(size.le.2*NL.and.icycle.ge.limav)then
                     HsCoolNS(size)=HsCoolNS(size)+1
                     if(temp_out.eq.'y')then
                        gaux=q*h/(2.d0*(zeff+zabs))
                        ig=int((gaux-hmin_g)/aa_g)
                        if(ig.gt.nb_g)ig=nb_g
                        if(ig.lt.0)ig=0
                        HsNS_g_c(size,ig)=HsNS_g_c(size,ig)+1
                     endif
                  endif
	       endif
c$$$               print*,'i_1D_span,i_2D_span=',i_1D_span,i_2D_span

               icycDif=ncycle-icycle

               gaux=q*h/(2.d0*(zeff+zabs))
               igtrue=0
               if(gaux.gt.goutmin.and.gaux.lt.goutmax)igtrue=1
               if(nshapes.gt.0.and.igtrue.eq.1)then
               do i=1,nshapes
                  if(t.eq.tshape(i).and.icycle.ge.limav)then
                     sumshape(i)=sumshape(i)+1
                     do k=1,tshape(i)
                        shapesum(k,i)=shapesum(k,i)+shape(k)
                     enddo
                     goto 501
                  endif
               enddo
               endif

               goto 501
            endif
            
!     Output of instabilities. Preliminary 21/02/2013
!            hlim_l=-8.d0*q+0.01
!            hlim_h=hlim_h+0.1
!            if(nshell.eq.0.and.h.lt.hlim_h.and.h.gt.hlim_l)then
!               do i=1,NL
!                  if(s(i).eq.0)then
!                     dlambda_0=abs(Th0(0,1)-abs(z(i)))
!                     print*,'dlambda_0=',dlambda_0
!                  endif
!                  if(abs(s(i)).eq.1)then
!                     dlambda_1=abs(abs(z(i))-Th0(1,2))
!                     dlambda_s=abs(Th0(1,1)-abs(z(i)))
!                     print*,'dlambda_1,dlambda_s=',dlambda_1,dlambda_s
!                  endif
!               enddo
!            endif

!            if(nflip.eq.NL)goto 2001
 1003    enddo
 501     continue



 1001 enddo
     

 2001 continue
	
c	stop
c      do i=1,NL
c         print*,i,'increasing',z(i)
c      enddo

c      stop
      if(icycle.ge.2)then
         f_s=0.d0
         f_d=0.d0
         f_w=0.d0
         do i=1,NL
            if(s(i).ne.s_M_prev(i))f_s=f_s+1.d0
            if(dis(i).ne.dis_M_prev(i))f_d=f_d+1.d0
            wprev=dble(dis_M_prev(i))+q*dble(s_M_prev(i))
            wcurrent=dble(dis(i))+q*dble(s(i))
            if(wcurrent.ne.wprev)f_w=f_w+1.d0
         enddo
         write(122,*)icycle,f_s/dble(NL),f_d/dble(NL),f_w/dble(NL)
      endif
      s_M_prev(:)=s(:)
      dis_M_prev(:)=dis(:)

      if(icycle.eq.ncycle)then
c         do i=1,NL
c            write(92,*) i,s(i)
c         enddo
         do y=1,L-1
            do x=1,L
               n=x+(y-1)*L
c     Spin configuration
               write(92,*) n,s(n),z(n)
               if(s(n).gt.0) write(93,*) x,y
               if(s(n).lt.0) write(94,*) x,y
               write(116,*)x,y,rf(n)/(qini(igrid)*fJ)
            enddo
         enddo
         print*,'tglob,iaval=',tglob,iaval-1
      endif

      sdi=0.0d0
      sdi2=0.0d0
      sti=0.0d0
      sti2=0.0d0
      stidi=0.0d0
      ni=0.0d0
      nt=0.0d0
      rf1=0.d0
      rf2=0.d0

c      print*,'hola2'
      do i=1,NL
         dii=0.0d0
         ti=0.0d0

         if(j(i,3).gt.0)then

            if(s(i).ne.s(j(i,3)))nt=nt+1.0d0
            if(dis(i).ne.dis(j(i,3)))then
               dii=1.0d0
               if(s(i).ne.s(j(i,3)))ti=1.0d0
               sdi=sdi+dii
               sdi2=sdi2+dii**2.0d0
               sti=sti+ti
               sti2=sti2+ti**2.0d0
               stidi=stidi+ti*dii
               ni=ni+1.0d0
            endif
         endif
         dii=0.0d0
         ti=0.0d0
c     ---- Conditional correlation -------
         if(j(i,4).gt.0)then
            if(s(i).ne.s(j(i,4)))nt=nt+1.0d0

            if(dis(i).ne.dis(j(i,4)))then
               dii=1.0d0
               if(s(i).ne.s(j(i,4)))ti=1.0d0
               sdi=sdi+dii
               sdi2=sdi2+dii**2.0d0
               sti=sti+ti
               sti2=sti2+ti**2.0d0
               stidi=stidi+ti*dii
               ni=ni+1.0d0
            endif
         endif

!Calculation of the statistics of the random slip field in martensite      
         rf1=rf1+rf(i)
         rf2=rf2+rf(i)**2.d0
         if(icycle.eq.ncycle)then
            rftemp=0.d0
            if(SR_LR.eq.'LR')then
            do ii=1,NL
               rftemp=rftemp+DJij(i,ii)*dis(ii)
            enddo
            endif
            if(SR_LR.eq.'SR')then
               do ii=0,8
                  if(j(i,ii).gt.0)then
                     rftemp=rftemp+DJSR(ii)*dis(j(i,ii))
                  endif
               enddo
            endif
            write(200,*)rftemp,rf(i)
         endif
      enddo

      if(ni.gt.0.0d0)then

      sdi=sdi/ni
      sdi2=sdi2/ni
      sti=sti/ni
      sti2=sti2/ni
      stidi=stidi/ni

      rho=stidi-sdi*sti
      rho=stidi
c      rho=rho/sqrt((sdi2-sdi**2.0d0)*(sti2-sti**2.0d0))
      write(98,*)icycle,rho,nt/float(2*(L-1)**2),ni/float(2*(L-1)**2)
!      print*,icycle,rho,nt/float(2*(L-1)**2),ni/float(2*(L-1)**2)

     

      avg_rftot=rf1/dble(NL)
      var_rftot=rf2/dble(NL)-avg_rftot**2.d0
      write(107,*)icycle,avg_rftot/(qini(igrid)*fJ),
     &     var_rftot/((qini(igrid)*fJ)**2.d0)
      endif

      avg_branch=avg_branch/dble(iaval)
      avg_branchPT=avg_branchPT/dble(iaval)
      avg_branchD=avg_branchD/dble(iaval)

      write(82,*)r,avg_branch,avg_branchPT,g_branch_min0,branch_min0,
     &     g_branch_max0,branch_max0,
     &     g_branch_minTot,branch_minTot,
     &     g_branch_maxTot,branch_maxTot

      !Spatial correlation of the slip fields if icycle>limav3
      if(icycle.ge.limav3)then
         do ix0=1,L
         do iy0=1,L-1
            
            call Corr_0_2D(ix0,iy0,L,rf,g0_h,N0_h)
            do idist=0,2*L
               if(N0_h(idist).gt.0)then
                  g_h(idist)=g_h(idist)+g0_h(idist)
                  N_h(idist)=N_h(idist)+N0_h(idist)
               endif
            enddo
         enddo
         enddo
      endif


      if(icycle.eq.1)then
         delta_max=-1.d0
         avg_delta=0.d0
         avg_e_s=0.d0
         do site=1,NL
            abs_z=abs(z(site))/q
            avg_delta=avg_delta+abs_z
            if(abs_z.gt.delta_max)delta_max=abs_z
            avg_e_s=avg_e_s+0.5d0/(1.d0+abs_z)
         enddo
         avg_delta=avg_delta/dble(NL)
         e_s=0.5d0/(1.d0+delta_max)
         e_s1=0.5d0/(1.d0+delta_max1)
         e_s2=0.5d0/(1.d0+avg_delta)
         avg_e_s=avg_e_s/dble(NL)
         write(118,*)cvalue(i_c),e_s,e_s1,avg_e_s,e_s2
      endif


c     ------ INCREASING h -------------
      if(increase.eq.'y')then
         if(mod(icycle,100).eq.0)print*,'increasing',icycle
!         print*,'increasing,h=',h         
!         print*,'NL=',NL

!      do i=1,NL
!         if(s(i).eq.-1)write(84,*)z(i) !,s(i)
!      enddo
!      close(84)
!      stop

      tglob=0

      Th0(0,1)=-big
      Th0(-1,1)=big
      Th0(0,2)=big
      Th0(1,2)=-big

      irun=0
      do 1011
         irun=irun+1
c         print*,'irun, nflip =',irun,nflip

c         write(80,*) irun,nflip

         h=big
         h1=-big
         h2=-big
         h3=big
         h4=big
         do i=1,NL
            y=int(i/L)
            if(mod(i,L).eq.0) y=y-1
            x=i-y*L-1

            if(x.gt.xmin.and.x.lt.xmax)then
            if(s(i).eq.0.and.z(i).gt.Th0(0,1).and.z(i).gt.h1)then
               h1=z(i)
               haux=-fDTemp*q*(zeff+zabs)+h1
               if(haux.lt.h) h=haux
            endif

            if(s(i).eq.-1.and.z(i).lt.Th0(-1,1).and.z(i).gt.h2)then
               h2=z(i)
               haux=fDTemp*q*(zeff+zabs)-h2
               if(haux.lt.h) h=haux
            endif

            if(s(i).eq.0.and.z(i).lt.Th0(0,2).and.z(i).lt.h3)then
               h3=z(i)
               haux=-fDTemp*q*(zeff+zabs)-h3
               if(haux.lt.h) h=haux
            endif

            if(s(i).eq.1.and.z(i).gt.Th0(1,2).and.z(i).lt.h4)then
               h4=z(i)
               haux=fDTemp*q*(zeff+zabs)+h4
               if(haux.lt.h) h=haux
            endif
            endif
         enddo

	 if(h.eq.big)goto 1432
	 h=h+feed
         Th0(0,1)=fDTemp*q*(zeff+zabs)+h
         Th0(1,2)=Th0(0,1)-2.0d0*fDTemp*q*(zeff+zabs)
         Th0(0,2)=-Th0(0,1)
         Th0(-1,1)=-Th0(0,1)+2.0d0*fDTemp*q*(zeff+zabs)
         g_new=q*h/(2.d0*(zeff+zabs))

         iaval=iaval+1
         iaval_tot=iaval_tot+1
         if(iaval_tot.ge.2)aval_dens=dble(size)/abs(g_new-g_old)
!         if(stress_out.eq.'y'.and.icycle.ge.limav)               

         cumm_f=float(nflip_abs)/float(NL)

         if(stress_out.eq.'y'.and.iaval_tot.ge.2)
     &        write(60,*) q*h/(2.d0*(zeff+zabs)),strain/float(NL),
     &        float(nflip)/float(NL),float(ndis)/float(NL),icycle,
     &     iaval_tot,1,cumm_f

         if(Dis_cyc_out.eq.'y'.and.iaval_tot.ge.2)
     &        write(112,*)q*h/(2.d0*(zeff+zabs)),
     &        avg_rftot/(qini(igrid)*fJ),
     &        var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &        iaval_tot,aval_dens,size,1,cumm_f
         
         g_old=q*h/(2.d0*(zeff+zabs))

         if(temp_out.eq.'y')then
         if(icycle.ge.limav)then            
            do i=1,NL
               gaux=q*h/(2.d0*(zeff+zabs))
               ig=int((gaux-hmin_g)/aa_g)
               if(ig.gt.nb_g)ig=nb_g
               if(ig.lt.0)ig=0

               E_rf_h(ig,1)=E_rf_h(ig,1)+rf(i)
               E_rf_h(ig,2)=E_rf_h(ig,2)+rf(i)**2.d0
               E_rf_h(ig,3)=E_rf_h(ig,3)+rf(i)**3.d0
               E_rf_h(ig,4)=E_rf_h(ig,4)+rf(i)**4.d0
               num_rf_h(ig)=num_rf_h(ig)+1
            enddo
         endif
         endif
!! Transforming sites artificially to the austenite phase
         if(forced_MA.eq.'y')then
            if(Th0(1,2).le.0.d0)then
               do i=1,NL
                  if(s(i).ne.0)s(i)=0
               enddo
               goto 1432
            endif
         endif

!         print*,'Th(0,1),Th(1,2),Th(0,2),Th(-1,1),Th(1,1),Th(-1,2)=',
!     &        Th0(0,1),Th0(1,2),Th0(0,2),Th0(-1,1),Th0(1,1),Th0(-1,2)

!         stop
c     ------ Propagating the avalanche ----
 1012    stress=stress-1
         size=0
         t=0
         sPT=0
         tPT=0
         sD=0
         ProxG=1
         Fincapa=ProxE
         nshell=0
	 maxA=0
c         print*,'hola'

         nava=0
         bar_branch=0.d0
         do 1013
            nava=nava+1
            nuns=0
	    do i=1,NL
c               if(z(i).gt.Th0(s(i),1).or.
c     &              z(i).lt.Th0(s(i),2))then
c		  nuns=nuns+1
c		  uns(nuns)=i
                  site=i
                  y=int(site/L)
                  if(mod(site,L).eq.0) y=y-1
                  x=site-y*L-1

                  if(x.gt.xmin.and.x.lt.xmax)then                  
c     casos indecisos
c                  if(z(site).gt.Th0(s(site),1).and.
c     &                 z(site).lt.Th0(s(site),2).and.
c     &                 z(site).eq.0.0d0)then
c                  if(z(site).gt.Th0(s(site),1).and.
c     &                 z(site).lt.Th0(s(site),2).and.
c     &                 s(site).eq.0.and.z(site).eq.0.0d0)then
                  if(z(site).gt.Th0(s(site),1).and.
     &                 z(site).lt.Th0(s(site),2))then
                     nuns=nuns+1
                     uns(nuns)=site
                     
                     if(z(site).eq.0.0d0)then
                        nthermal=nthermal+1
                        ivec=ivec+1
                        if (ivec.ge.lcarry) then
                           call rcarry(rvec,lcarry)
                           ivec=1
                        endif
                        if(rvec(ivec).le.0.5)then
                           size=size+1
                           strain=strain+zch(s(site),1)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.-1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.1)then
                              ndis=ndis+1
                              dis(site)=dis(site)+1
                           endif
                           
!                           z(site)=z(site)-f*zeff*zch(s(site),1)
                           zchnn(site)=zch(s(site),1)
                           s(site)=sch(s(site),1)
                           density=density-float(z1)*zch(s(site),1)
                           ndw=ndw+1
                        else
                           size=size+1
                           strain=strain+zch(s(site),2)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.-1)then
                              ndis=ndis-1
                              dis(site)=dis(site)-1
                           endif
                           
!                           z(site)=z(site)-f*zeff*zch(s(site),2)
                           zchnn(site)=zch(s(site),2)
                           s(site)=sch(s(site),2)
c     density=density-float(z1)*zch(s(site),2)
                           nup=nup+1
                        endif
                        
                     else
                        if(z(site).gt.0)then
                           size=size+1
                           strain=strain+zch(s(site),1)
                           
                           if(s(site).eq.0)then
                              sPT=sPT+1
                              nflip=nflip+1
                              nflip_abs=nflip_abs-1
                           endif
                           
                           if(s(site).eq.-1)then
                              sPT=sPT+1
                              nflip=nflip-1
                              nflip_abs=nflip_abs+1
                           endif
                           
                           if(s(site).eq.1)then
                              ndis=ndis+1
                              dis(site)=dis(site)+1
                           endif
                           
!                           z(site)=z(site)-f*zeff*zch(s(site),1)
                           zchnn(site)=zch(s(site),1)
                           s(site)=sch(s(site),1)
                           density=density-float(z1)*zch(s(site),1)
                           ndw=ndw+1
                        else
                           if(z(site).lt.0)then
                              size=size+1
                              strain=strain+zch(s(site),2)
                              
                              if(s(site).eq.0)then
                                 sPT=sPT+1
                                 nflip=nflip+1
                                 nflip_abs=nflip_abs-1
                              endif
                              
                              if(s(site).eq.1)then
                                 sPT=sPT+1
                                 nflip=nflip-1
                                 nflip_abs=nflip_abs+1
                              endif
                              
                              if(s(site).eq.-1)then
                                 ndis=ndis-1
                                 dis(site)=dis(site)-1
                              endif
                              
!                              z(site)=z(site)-f*zeff*zch(s(site),2)
                              zchnn(site)=zch(s(site),2)
                              s(site)=sch(s(site),2)
c     density=density-float(z1)*zch(s(site),2)
                              nup=nup+1
                           endif
                        endif
                     endif
                     goto 1431
                  endif
                  
                  if(z(site).gt.Th0(s(site),1))then
                     size=size+1
                     strain=strain+zch(s(site),1)
                     
                     if(s(site).eq.0)then
                        sPT=sPT+1
                        nflip=nflip+1
                        nflip_abs=nflip_abs-1
                     endif
                     
                     if(s(site).eq.-1)then
                        sPT=sPT+1
                        nflip=nflip-1
                        nflip_abs=nflip_abs+1
                     endif
                     
                     if(s(site).eq.1)then
                        ndis=ndis+1
                        dis(site)=dis(site)+1
                     endif
                     
!                     z(site)=z(site)-f*zeff*zch(s(site),1)
                     zchnn(site)=zch(s(site),1)
                     s(site)=sch(s(site),1)
                     density=density-float(z1)*zch(s(site),1)
                     nuns=nuns+1
                     uns(nuns)=site
                  else
                     
                     if(z(site).lt.Th0(s(site),2))then
                        size=size+1
                        strain=strain+zch(s(site),2)
                        
                        if(s(site).eq.0)then
                           sPT=sPT+1
                           nflip=nflip+1
                           nflip_abs=nflip_abs-1
                        endif
                        
                        if(s(site).eq.1)then
                           sPT=sPT+1
                           nflip=nflip-1
                           nflip_abs=nflip_abs+1                           
                        endif
                        
                        if(s(site).eq.-1)then
                           ndis=ndis-1
                           dis(site)=dis(site)-1
                        endif
                        
!                        z(site)=z(site)-f*zeff*zch(s(site),2)
                        zchnn(site)=zch(s(site),2)
                        s(site)=sch(s(site),2)
                        density=density-float(z1)*zch(s(site),2)

                        nuns=nuns+1
                        uns(nuns)=site
                     endif   
                  endif
                  endif
 1431             continue                              
	    enddo

!            print*,'nuns,ndis,nflip,icycle heating=',
!     &           nuns,ndis,nflip,icycle

!            if(nflip.eq.NL)nuns=0

            if(nuns.gt.NL)print*,'Problem: nuns > NL'

	    if(nava.eq.1.and.nuns.eq.0)goto 1432

!	    print*,'increasing nuns,nflip,ndis,h=',nuns,nflip,ndis,h
!            print*,'Th0(0,1),Th0(1,2),Th0(1,1)=',
!     &           Th0(0,1),Th0(1,2),Th0(1,1)
c            nshell=nuns
	  
c     update of the neighbors
            do k=1,nuns
               site=uns(k)
               rf1=0.d0
               rf2=0.d0
               if(SR_LR.eq.'LR')then
               do i=1,NL
!     if(i.ne.site)z(i)=z(i)+DJij(i,site)*zchnn(site)
                  z(i)=z(i)+DJij(i,site)*zchnn(site)
                  if(zchnn(site).eq.zch(-1,2))rf(i)=rf(i)-DJij(i,site)
                  if(zchnn(site).eq.zch(1,1))rf(i)=rf(i)+DJij(i,site)
                  rf1=rf1+rf(i)
                  rf2=rf2+rf(i)**2.d0
               enddo
               endif

               if(SR_LR.eq.'SR')then
                  do i1=0,8
                     i=j(site,i1)
                     if(j(site,i1).gt.0)then
                        z(i)=z(i)+DJSR(i1)*zchnn(site)
                        if(zchnn(site).eq.zch(-1,2))rf(i)=rf(i)
     &                       -DJSR(i1)
                        if(zchnn(site).eq.zch(1,1))rf(i)=rf(i)
     &                       +DJSR(i1)
                        
                        rf1=rf1+rf(i)
                        rf2=rf2+rf(i)**2.d0
                     endif
                  enddo
               endif

               y=int(site/L)
               if(mod(site,L).eq.0) y=y-1
               x=site-y*L-1
               sox(x)=1
               soy(y)=1

               if(x.lt.xSmin)xSmin=x
               if(x.gt.xSmax)xSmax=x
               if(y.lt.ySmin)ySmin=y
               if(y.gt.ySmax)ySmax=y

               if(sdout.eq.'y'.and.icycle.eq.ncycle)then 
                  write(73,*) x,y,tglob,iaval,s(site),dis(site),
     &                 q*h/(2.d0*(zeff+zabs)),rf(site)/(qini(igrid)*fJ),
     &                 z(site)/q,0.5d0/(1.d0+dble(s(site))*z(site)/q),
     &                 z(site)/q-rf(site)/q
               endif

	    enddo

            nshell=nuns
            
	    if(nshell.ne.0)then
               t=t+1
               tglob=tglob+1
!               if(t.le.nTmax.and.icycle.ge.limav)then
               if(t.le.nTmax)then
                  shape(t)=nshell
               endif
	       if(nshell.gt.maxA) maxA=nshell
               if(t.gt.1)bar_branch=bar_branch+
     &              dble(shape(t))/dble(shape(t-1))
               if(sdout.eq.'y'.and.icycle.eq.ncycle)then
                  write(73,*)
                  write(73,*)
               endif
            endif

c            if(nflip.eq.NL) goto 1014

            if(nshell.eq.0)then
               if(size.gt.size_snap)i_snap_heat=1
               if(size.le.2*NL.and.icycle.ge.limav)then
                  Hs(size)=Hs(size)+1
                  if(temp_out.eq.'y')then
                     gaux=q*h/(2.d0*(zeff+zabs))
                     ig=int((gaux-hmin_g)/aa_g)
                     if(ig.gt.nb_g)ig=nb_g
                     if(ig.lt.0)ig=0
                     Hs_g_h(size,ig)=Hs_g_h(size,ig)+1
                  endif
               endif
               if(t.le.NL.and.icycle.ge.limav) Ht(t)=Ht(t)+1
               if(sPT.le.2*NL.and.icycle.ge.limav) HsPT(sPT)=HsPT(sPT)+1
	       if(maxA.le.NL.and.icycle.ge.limav) HA(maxA)=HA(maxA)+1
	       
	       if(sPT.gt.MaxsPT) MaxsPT=sPT

! Output of hysteresis cycles with the published driving parameter: g = q*h/(4.d0*J1*z1)
!	       if(stress_out.eq.'y'.and.icycle.ge.limav)

               cumm_f=float(nflip_abs)/float(NL)

               if(stress_out.eq.'y'.and.iaval_tot.ge.2)
     &              write(60,*) q*h/(2.d0*(zeff+zabs)),strain/float(NL),
     &              float(nflip)/float(NL),float(ndis)/float(NL),icycle,
     &     iaval_tot,1,cumm_f
     
!     Output of Disorder after avalanche
               avg_rftot=rf1/dble(NL)
               var_rftot=rf2/dble(NL)-avg_rftot**2.d0
               
               if(Dis_cyc_out.eq.'y'.and.iaval_tot.ge.2)
     &              write(112,*)q*h/(2.d0*(zeff+zabs)),
     &              avg_rftot/(qini(igrid)*fJ),
     &              var_rftot/((qini(igrid)*fJ)**2.d0),icycle,
     &              iaval_tot,aval_dens,size,1,cumm_f

	       if(mapST.eq.'y') write(91,*) sPT,size,t
     
               if(icycle.ge.limav)then
                  if(bar_branch.gt.0.d0)then
                     duration=dble(t)
                  baux=bar_branch/duration
                  if(temp_out.eq.'y')then
                     gaux=q*h/(2.d0*(zeff+zabs))
                     ig=int((gaux-hmin_g)/aa_g)
                     if(ig.gt.nb_g)ig=nb_g
                     if(ig.lt.0)ig=0
                     obs_branch_g_h(ig)=obs_branch_g_h(ig)+1.d0
                     Sum1_branch_g_h(ig)=Sum1_branch_g_h(ig)+baux
                  endif
                  endif
               endif

               gaux=q*h/(2.d0*(zeff+zabs))
               igtrue=0
               if(gaux.gt.goutmin.and.gaux.lt.goutmax)igtrue=1
               if(nshapes.gt.0.and.igtrue.eq.1)then
               do i=1,nshapes
                  if(t.eq.tshape(i).and.icycle.ge.limav)then
                     sumshape(i)=sumshape(i)+1
                     do k=1,tshape(i)
                        shapesum(k,i)=shapesum(k,i)+shape(k)
                     enddo
                     goto 511
                  endif
               enddo
               endif

               goto 511
            endif

 1013    enddo

 511     continue
 1011 enddo

c	endif of for the increasing
      endif 	

 1432 continue
!      if(stress_out.eq.'y'.and.icycle.ge.limav)then
      if(stress_out.eq.'y')then
!         write(60,*)
      endif
      if(Dis_cyc_out.eq.'y')then
!         write(112,*)
      endif
      
      if(zout.eq.'y')then
!         if(icycle.eq.ncycle.or.icycle.eq.1)then
         if(icycle.eq.ncycle)then
            do i=1,NL
c     print*,z(i)
               write(90,*) i,s(i),z(i),dis(i)
            enddo
            write(90,*) ' '
         endif
      endif

      if(icycle.eq.ncycle)then
c         do i=1,NL
c         write(92,*) i,s(i)
c         enddo
         do y=1,L-1
            do x=1,L
               n=x+(y-1)*L
c               write(92,*) n,s(n),z(n)
c               if(s(n).gt.0) write(93,*) nx-1,ny
c               if(s(n).lt.0) write(94,*) nx-1,ny
	       if(dis(n).gt.0) write(95,*) x,y,dis(n)
               if(dis(n).lt.0) write(96,*) x,y,dis(n)
!               if(sdout.eq.'y') write(73,*) x,y,s(n),dis(n)
            enddo
!            if(sdout.eq.'y') write(73,*)
         enddo


c     ------------ configuration of D-gradient links -------------
         do y0=3,L-1,2
            do y=y0,L-1
               x=y-y0+1
               i=x+(y-1)*L
               do k=1,z1
                  if(j(i,k).gt.0.and.dis(i).ne.dis(j(i,k)))then
                     write(97,*) float(x)+cgrad(k,1),
     &                    float(y)+cgrad(k,2),
     &                    min(dis(i),dis(j(i,k))),
     &                    abs(dis(j(i,k))-dis(i))
c                     write(97,*) float(x),float(y)
                  endif
c     --- Open boundaries ---
      if(BC.eq.'FBC'.and.j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------               
               enddo
            enddo
         enddo
         do x0=1,L,2
            do x=x0,L
               y=x-x0+1
               i=x+(y-1)*L
               do k=1,z1
                  if(j(i,k).gt.0.and.dis(i).ne.dis(j(i,k)))then
                     write(97,*) float(x)+cgrad(k,1),
     &                    float(y)+cgrad(k,2),
     &                    min(dis(i),dis(j(i,k))),
     &                    abs(dis(j(i,k))-dis(i))
c                     write(97,*) float(x),float(y)
                  endif
c     --- Open boundaries ---
      if(BC.eq.'FBC'.and.j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------               
               enddo
            enddo
         enddo
c     -----------------------------------------------------------

      endif

      sumdis=0
      sumdis2=0
      sumdis3=0
      absndis=0
      rhograd=0.0d0
      do i=1,NL
         sumdis=sumdis+dis(i)
         sumdis2=sumdis2+dis(i)**2
         sumdis3=sumdis3+dis(i)**3
         if(dis(i).ne.0) absndis=absndis+1

         do k=1,z1
            if(j(i,k).gt.0.and.dis(i).ne.dis(j(i,k)))then
               rhograd=rhograd+1.0d0
            endif
c     --- Open boundaries ---
      if(BC.eq.'FBC'.and.j(i,k).lt.0) z(i)=z(i)+J1*Di
c     ---------------------------
         enddo
      enddo

      skew=float(sumdis3)/float(absndis)
     &     -3.0d0*float(sumdis2)*float(sumdis)/float(absndis**2)+
     &     2.0d0*(float(sumdis)/float(absndis))**3
      skew=skew/(float(sumdis2)/float(absndis)
     &     -(float(sumdis)/float(absndis))**2.0)**(1.5d0)

      rhograd=rhograd/2.0d0
      write(72,*) icycle,MaxsPT,float(absndis)/float(NL),
     & rhograd/float(NL)
c      write(*,*) icycle,MaxsPT,float(absndis)/float(NL),
c     & rhograd/float(NL)
c,
c     & float(sumdis)/float(NL),
c     & sqrt(float(sumdis2)/float(NL)-(float(sumdis)/float(NL))**2.0d0)
c     & ,nthermal
c     ,     & skew
      
 1000 enddo

      
 1004 sums=0
      sumt=0
      sumsPT=0
      sumA=0
      sumsCool=0
      sumsCoolNS=0
      sumsCool1D=0
      sumsCool2D=0
      sumtCool=0
      sumsPTCool=0
      do i=1,NL
         sums=sums+Hs(i)
         sumt=sumt+Ht(i)
         sumsPT=sumsPT+HsPT(i)
	 sumA=sumA+HA(i)
         sumsCool=sumsCool+HsCool(i)
         sumtCool=sumtCool+HtCool(i)
         sumsPTCool=sumsPTCool+HsPTCool(i)
         sumsCoolNS=sumsCoolNS+HsCoolNS(i)
         sumsCool1D=sumsCool1D+HsCool1D(i)
         sumsCool2D=sumsCool2D+HsCool2D(i)
      enddo
      do i=NL+1,2*NL
         sums=sums+Hs(i)
         sumsCool=sumsCool+HsCool(i)
         sumsPT=sumsPT+HsPT(i)
         sumsPTCool=sumsPTCool+HsPTCool(i)
         sumsCoolNS=sumsCoolNS+HsCoolNS(i)
         sumsCool1D=sumsCool1D+HsCool1D(i)
         sumsCool2D=sumsCool2D+HsCool2D(i)
      enddo

      write(119,*)r,q,cvalue(i_c),i_snap_cool,i_snap_heat

      frac1=dble(Ns1_glob)/dble(ncycle)
      frac2=dble(Ns2_glob)/dble(ncycle)
      i_1D_span_f=0
      i_2D_span_f=0
      if(frac1.ge.1.d0)i_1D_span_f=1
      if(frac2.ge.1.d0)i_2D_span_f=1
      
      write(120,*)r,q,cvalue(i_c),i_1D_span,i_2D_span,
     &     i_1D_span_f,i_2D_span_f,frac1,frac2

      print*,r,q,cvalue(i_c),i_1D_span,i_2D_span,
     &     i_1D_span_f,i_2D_span_f,frac1,frac2

      open(20,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//
     &     '_'//c_name(i_c)//'.stat_info')
      write(20,*) 'sums=',sums
      write(20,*) 'sumsPT=',sumSPT
      write(20,*) 'sumt=',sumt
      write(20,*) 'sumA=',sumA

      open(30,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_s.dat')
      open(31,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_sNS.dat')
      open(32,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_s1D.dat')
      open(33,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_s2D.dat')
      open(40,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_t.dat')
      open(70,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_sPT.dat')
      open(71,file='results/'//filename//'/'//filename
     &        //extension(igrid)//'_'//BC//'_'//c_name(i_c)//'_A.dat')



      do i=1,NL
         if(Hs(i).gt.0.or.HsCool(i).gt.0) 
     &        write(30,*) i,float(Hs(i))/float(sums),
     &        float(HsCool(i))/float(sumsCool),sums,sumsCool
         if(HsCoolNS(i).gt.0)      
     &        write(31,*) i,float(HsCoolNS(i))/float(sumsCoolNS),
     &        sumsCoolNS
         if(HsCool1D(i).gt.0)      
     &        write(32,*) i,float(HsCool1D(i))/float(sumsCool1D),
     &        sumsCool1D
         if(HsCool2D(i).gt.0)      
     &        write(33,*) i,float(HsCool2D(i))/float(sumsCool2D),
     &        sumsCool2D
         if(Ht(i).gt.0.or.HtCool(i).gt.0)
     &        write(40,*) i,float(Ht(i))/float(sumt),
     &        float(HtCool(i))/float(sumtCool),sumt,sumtCool
         if(HsPT(i).gt.0.or.HsPTCool(i).gt.0)
     &        write(70,*) i,float(HsPT(i))/float(sumsPT),
     &        float(HsPTCool(i))/float(sumsPTCool),sumsPT,sumsPTCool
	 if(HA(i).gt.0) write(71,*) i,float(HA(i))/float(sumA)
      enddo
      do i=NL+1,2*NL
         if(Hs(i).gt.0.or.HsCool(i).gt.0) 
     &        write(30,*) i,float(Hs(i))/float(sums),
     &        float(HsCool(i))/float(sumsCool),sums,sumsCool
         if(HsPT(i).gt.0.or.HsPTCool(i).gt.0)
     &        write(70,*) i,float(HsPT(i))/float(sumsPT),
     &        float(HsPTCool(i))/float(sumsPTCool),sumsPT,sumsPTCool
         if(HsCoolNS(i).gt.0)      
     &        write(31,*) i,float(HsCoolNS(i))/float(sumsCoolNS),
     &        sumsCoolNS
         if(HsCool1D(i).gt.0)      
     &        write(32,*) i,float(HsCool1D(i))/float(sumsCool1D),
     &        sumsCool1D
         if(HsCool2D(i).gt.0)      
     &        write(33,*) i,float(HsCool2D(i))/float(sumsCool2D),
     &        sumsCool2D
      enddo
      
      if(nshapes.gt.0)then
      do k=1,nshapes
         if(sumshape(k).gt.0)then
            open(50,file='results/shapes/'//filename//shapename(k)
     &           //'_'//extension(igrid)//'_'//BC
     &           //'_'//c_name(i_c)//'.dat')

            do i=1,tshape(k)
               shape_norm=dble(shapesum(i,k))/dble(sumshape(k))
!               shape_norm=shape_norm/dble(tshape(k))
               write(50,*) dble(i)/dble(tshape(k)),shape_norm
            enddo
            close(50)
         endif
      enddo
      endif

      do i=0,nb_branch
         xbranch=hmin_branch+aa_branch*dble(i)
         if(h_branch(i).gt.0.d0)write(85,*)xbranch,
     &        h_branch(i)/sum_branch
      enddo
      
      avg_branch_max=-1.d0
      do ii=0,nb_g
         xg=hmin_g+aa_g*dble(ii)
         sum1=0.d0
         sum0=0.d0
         branch_max=-big
         branch_min=big
         do i=0,nb_branch
            xbranch=hmin_branch+aa_branch*dble(i)
            aux_branch_g=h_branch_g(i,ii)/(aa_branch*obs_branch_g_c(ii))
            if(h_branch_g(i,ii).gt.0.d0)then
               write(86,*)xg,xbranch,
     &           aux_branch_g
               sum0=sum0+aux_branch_g
               sum1=sum1+aux_branch_g*xbranch
               if(xbranch.gt.branch_max)branch_max=xbranch
               if(xbranch.lt.branch_min)branch_min=xbranch
            endif
         enddo
         write(86,*)
         do i=0,nb_branch
            xbranch=hmin_branch+aa_branch*dble(i)
            aux_branchPT_g=
     &           h_branchPT_g(i,ii)/(aa_branch*obs_branchPT_g_c(ii))
            if(h_branchPT_g(i,ii).gt.0.d0)then
               write(124,*)xg,xbranch,
     &           aux_branchPT_g
            endif
         enddo
         write(124,*)
!         write(86,*)
         avg_R=sum1/sum0
         if(avg_R.gt.avg_branch_max)then
            avg_branch_max=avg_R
            g_branch_max=xg
         endif
         if(obs_branch_g_c(ii).gt.0.d0)
     &        write(99,*)xg,Sum1_branch_g_c(ii)/obs_branch_g_c(ii)
         if(obs_branch_g_h(ii).gt.0.d0)
     &        write(114,*)xg,Sum1_branch_g_h(ii)/obs_branch_g_h(ii)

         if(obs_branchPT_g_c(ii).gt.0.d0)then
            avgbranchPT=Sum1_branchPT_g_c(ii)/obs_branchPT_g_c(ii)
            sdbranchPT2=Sum2_branchPT_g_c(ii)/obs_branchPT_g_c(ii)
            sdbranchPT2=sqrt(sdbranchPT2-avgbranchPT**2.d0)
            upperbranchPT=avgbranchPT+sdbranchPT2
            lowerbranchPT=avgbranchPT-sdbranchPT2
            write(127,*)xg,avgbranchPT,sdbranchPT2,
     &           lowerbranchPT,upperbranchPT
         endif
         if(obs_branchD_g_c(ii).gt.0.d0)
     &        write(128,*)xg,Sum1_branchD_g_c(ii)/obs_branchD_g_c(ii)
      enddo

      write(104,*)r,q,cvalue(i_c),g_branch_max,avg_branch_max

!      if(stab_out.eq.'y')then
         Deltah_cool_min=big
         Deltah_cool_max=-big
      do ii=0,nb_g
         xg=hmin_g+aa_g*dble(ii)
         if(num_rf(ii).gt.0)then
            E_rf_1=E_rf(ii,1)/(dble(num_rf(ii))*(qini(igrid)*fJ))
            E_rf_2=E_rf(ii,2)/(dble(num_rf(ii))*(qini(igrid)*fJ)**2.d0)
            E_rf_3=E_rf(ii,3)/(dble(num_rf(ii))*(qini(igrid)*fJ)**3.d0)
            E_rf_4=E_rf(ii,4)/(dble(num_rf(ii))*(qini(igrid)*fJ)**4.d0)
            dnu1=E_rf_1
            dnu2=E_rf_2-dnu1**2.d0
            dnu3=E_rf_3-3.d0*E_rf_2*dnu1+2.d0*dnu1**3.d0
            dnu4=E_rf_4-4.d0*E_rf_3*dnu1+
     &           6.d0*E_rf_2*dnu1**2.d0-3*dnu1**4.d0
            skew_rf=dnu3/dnu2**(3.d0/2.d0)
            exc_kurt_rf=dnu4/dnu2**2.d0-3.d0
            write(106,*)xg,dnu1,dnu2,skew_rf,exc_kurt_rf

            E_Js_1=E_Js_c(ii,1)/(dble(num_rf(ii))*(qini(igrid)*fJ))
            E_Js_2=E_Js_c(ii,2)/
     &           (dble(num_rf(ii))*(qini(igrid)*fJ)**2.d0)
            E_Js_3=E_Js_c(ii,3)/
     &           (dble(num_rf(ii))*(qini(igrid)*fJ)**3.d0)
            E_Js_4=E_Js_c(ii,4)/
     &           (dble(num_rf(ii))*(qini(igrid)*fJ)**4.d0)
            dnu1=E_Js_1
            dnu2=E_Js_2-dnu1**2.d0
            dnu3=E_Js_3-3.d0*E_Js_2*dnu1+2.d0*dnu1**3.d0
            dnu4=E_Js_4-4.d0*E_Js_3*dnu1+
     &           6.d0*E_Js_2*dnu1**2.d0-3*dnu1**4.d0
            skew_rf=dnu3/dnu2**(3.d0/2.d0)
            exc_kurt_rf=dnu4/dnu2**2.d0-3.d0
            write(123,*)xg,dnu1,dnu2,skew_rf,exc_kurt_rf

            Deltah_cool=sqrt(dnu2)
            if(Deltah_cool.lt.Deltah_cool_min)
     &           Deltah_cool_min=Deltah_cool
            if(Deltah_cool.gt.Deltah_cool_max)
     &           Deltah_cool_max=Deltah_cool
         endif
         if(num_rf_h(ii).gt.0)then
            E_rf_1=E_rf_h(ii,1)/(dble(num_rf_h(ii))*(qini(igrid)*fJ))
            E_rf_2=
     &         E_rf_h(ii,2)/(dble(num_rf_h(ii))*(qini(igrid)*fJ)**2.d0)
            E_rf_3=
     &         E_rf_h(ii,3)/(dble(num_rf_h(ii))*(qini(igrid)*fJ)**3.d0)
            E_rf_4=
     &         E_rf_h(ii,4)/(dble(num_rf_h(ii))*(qini(igrid)*fJ)**4.d0)
            dnu1=E_rf_1
            dnu2=E_rf_2-dnu1**2.d0
            dnu3=E_rf_3-3.d0*E_rf_2*dnu1+2.d0*dnu1**3.d0
            dnu4=E_rf_4-4.d0*E_rf_3*dnu1+
     &           6.d0*E_rf_2*dnu1**2.d0-3.d0*dnu1**4.d0            
            skew_rf=dnu3/dnu2**(3.d0/2.d0)
            exc_kurt_rf=dnu4/dnu2**2.d0-3.d0
            write(110,*)xg,dnu1,dnu2,skew_rf,exc_kurt_rf
         endif
         do i=0,nb_rf
            xrf=(rf_min+aa_rf*dble(i))/(qini(igrid)*fJ)
            if(h_rf_g(i,ii).gt.0.d0)write(109,*)xg,xrf,
     &           qini(igrid)*fJ*h_rf_g(i,ii)/(aa_rf*sum_rf_g(ii)),
     &           sum_rf_g(ii),aa_rf/(qini(igrid)*fJ)
         enddo
         write(109,*)
         do i=0,nb_stab
            delta=ddmin+aa*dble(i)

!            print*,ii,i,xg,delta,hstab_g(i,ii,1),sum_stab_g(ii,1)
            if(hstab_g(i,ii,1).gt.0.d0)write(87,*)xg,delta,
     &           hstab_g(i,ii,1)/(aa*sum_stab_g(ii,1)),
     &	         sum_stab_g(ii,1),aa
!            write(87,*)xg,delta,
!     &           hstab_g(i,ii,1)/sum_stab_g(ii,1)
            if(hstab_g(i,ii,2).gt.0.d0)write(88,*)xg,delta,
     &           hstab_g(i,ii,2)/(aa*sum_stab_g(ii,2)),
     &	         sum_stab_g(ii,2),aa
            if(hstab_g(i,ii,3).gt.0.d0)write(89,*)xg,delta,
     &           hstab_g(i,ii,3)/(aa*sum_stab_g(ii,3)),
     &           sum_stab_g(ii,3),aa
         enddo
         
         do i=0,nb_Deltai
            xD=Deltai_min+aa_Deltai*dble(i)
            if(hDeltai_g(i,ii).gt.0.d0)write(100,*)xg,xD,
     &           hDeltai_g(i,ii)/(aa_Deltai*sum_Deltai_g(ii)),
     &           sum_Deltai_g(ii),aa
         enddo

         do i=0,nb_dnflip
            xnflip=dnflip_min+aa_dnflip*dble(i)
            if(hdnflip_g(i,ii).gt.0.d0)write(102,*)xg,xnflip,
     &           hdnflip_g(i,ii)/(aa_dnflip*sum_dnflip_g(ii)),
     &           sum_dnflip_g(ii),aa
         enddo
!         write(87,*)
!         write(87,*)
!         write(88,*)
!         write(88,*)
!         write(89,*)
!         write(89,*)
         sums_g_c=0
         sumsNS_g_c=0
         sums1D_g_c=0
         sums2D_g_c=0
         sumsPT_g=0
         sums_g_h=0
         do i=1,2*NL
            sums_g_c=sums_g_c+Hs_g_c(i,ii)
            sumsNS_g_c=sumsNS_g_c+HsNS_g_c(i,ii)
            sums1D_g_c=sums1D_g_c+Hs1D_g_c(i,ii)
            sums2D_g_c=sums2D_g_c+Hs2D_g_c(i,ii)
            sumsPT_g=sumsPT_g+HsPT_g(i,ii)
            sums_g_h=sums_g_h+Hs_g_h(i,ii)
         enddo
         do i=1,2*NL
            if(i.le.NL.and.HsPT_g(i,ii).gt.0)write(105,*)xg,i,
     &           dble(HsPT_g(i,ii))/dble(sumsPT_g),sumsPT_g
            if(Hs_g_c(i,ii).gt.0)write(103,*)xg,i,
     &           dble(Hs_g_c(i,ii))/dble(sums_g_c),sums_g_c
            if(HsNS_g_c(i,ii).gt.0)write(34,*)xg,i,
     &           dble(HsNS_g_c(i,ii))/dble(sumsNS_g_c),sumsNS_g_c
            if(Hs1D_g_c(i,ii).gt.0)write(35,*)xg,i,
     &           dble(Hs1D_g_c(i,ii))/dble(sums1D_g_c),sums1D_g_c
            if(Hs2D_g_c(i,ii).gt.0)write(36,*)xg,i,
     &           dble(Hs2D_g_c(i,ii))/dble(sums2D_g_c),sums2D_g_c
            if(Hs_g_h(i,ii).gt.0)write(111,*)xg,i,
     &           dble(Hs_g_h(i,ii))/dble(sums_g_h),sums_g_h
         enddo
         write(103,*)
         write(34,*)
         write(35,*)
         write(36,*)
         write(105,*)
         write(111,*)
      enddo

      diff_aux=Deltah_cool_min-Deltah_cool_max
      nDeltah_width=1
      if(diff_aux.eq.0.d0)nDeltah_width=0
      write(115,*)cvalue(i_c),qini(igrid),
     &     Deltah_cool_min,Deltah_cool_max
     &     ,nDeltah_width

      do i=0,nb_Deltai
         xD=Deltai_min+aa_Deltai*dble(i)
         if(hDeltai(i).gt.0.d0)write(101,*)xD,
     &        hDeltai(i)/(aa_Deltai*sum_Deltai),
     &        sum_Deltai,aa_Deltai
      enddo
!      endif

! Output of slip field correlations
      do idist=0,2*L
         if(N_h(idist).gt.0)write(117,*)idist,
     &        g_h(idist)/(dble(N_h(idist))*(qini(igrid)*fJ)**2.d0)
      enddo
! ---------------------------------

      close(20)
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
      close(36)
      close(40)
      if(stress_out.eq.'y')close(60)
      close(70)
      close(71)
      close(72)
      if(sdout.eq.'y') close(73)
      if(zout.eq.'y') close(90)
      if(mapST.eq.'y') close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(96)
      close(97)
      close(98)

      close(82)
!      close(83)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)
      close(99)
      close(100)
      close(101)
      close(102)
      close(103)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      if(Dis_cyc_out.eq.'y')close(112)
      close(114)
      close(116)
      close(117)
      close(122)
      close(123)
      close(124)
      close(125)
      close(127)
      enddo

      close(104)
      close(115)
      close(119)
      close(120)

      enddo !loop for i_c. Values of c

      close(118)

      end

! --- Correlation function of slip disorder
      subroutine Corr_0_2D(ix0,iy0,L,rf,g0_h,N0_h)
      implicit real*8 (a-h,o-z)
      parameter(long=111)
      parameter(nTmax=10000)
      parameter(big=1.d15)   
      dimension g0_h(0:2*long)
      dimension N0_h(0:2*long)
      double precision rf(long**2) !effective random fields

      g0_h(:)=0.d0
      N0_h(:)=0
      n0=ix0+(iy0-1)*L
      do ix=1,L
      do iy=1,L-1
         n=ix+(iy-1)*L
         Dx=dble(ix-ix0)
         Dy=dble(iy-iy0)
         dist=sqrt(Dx**2.d0+Dy**2.d0)
         ih=int(dist)
         gaux=rf(n0)*rf(n)
!         print*,ix,iy,iz,gaux
         g0_h(ih)=g0_h(ih)+gaux
         N0_h(ih)=N0_h(ih)+1
      enddo
      enddo
      
      return
      end

c     ******************************************************************
c     *                     SUBRUTINA RCARIN                           *
c     ******************************************************************
c

        SUBROUTINE RCARIN(IJKL,RVEC,LENV)

C----------------------------------------------------------------------
C Inicializa valores antes de llamar a la subrutina RCARRY.
C IJKL debe estar en el rango 0<IJKL<900 000 000.
C Para conseguir los valores standar usados por Marsaglia y Zaman en su
C articulo poner IJKL = 54217137 (I=12, J=34, K=56, L=78)
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
        COMMON /RAN1/ CARRY
        DIMENSION RVEC(LENV+24)

        IJ = IJKL/30082
        KL = IJKL - 30082*IJ
        I = MOD(IJ/177,177) + 2
        J = MOD(IJ,177)     + 2
        K = MOD(KL/169,178) + 1
        L = MOD(KL,169)

        DO 2 II=24,1,-1
          S = 0.0
          T = 0.5
          DO 3 JJ=1,24
            M = MOD(MOD(I*J,179)*K,179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1,169)
            IF (MOD(L*M,64).GE.32) S = S+T
            T = 0.5*T
3         CONTINUE
          RVEC(II) = S
2       CONTINUE

        CARRY = 0.0

        RETURN
        END
c
c     ******************************************************************
c     *                     SUBRUTINA RCARRY                           *
c     ******************************************************************
c
        SUBROUTINE RCARRY(RVEC,LENV)
C----------------------------------------------------------------------
C Generador de numeros pseudo-aleatorios. Algoritmo de G. Marsaglia y
C A. Zaman. Genera numeros reales de 32-bits con mantisas de 24 bits,
C comprendidos entre 0 y 1 (1, explicitamente excluido).
C Periodo aproximado : 10**171.
C Admite la generacion de subsecuencias disjuntas.
C                   F. James, 1989
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
        DIMENSION RVEC(LENV+24)
        COMMON /RAN1/ CARRY
        PARAMETER (TWOM24=1.0/16777216.0)
C
        DO 100 IVEC=25,LENV+24
          UNI = RVEC(IVEC-24) - RVEC(IVEC-10) - CARRY
          IF (UNI.LT.0.) THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
          ELSE
            CARRY = 0.0
          ENDIF

          IF(UNI.EQ.0.)THEN
            UNI=RVEC(IVEC-24)*TWOM24
            in48=-48
            IF(UNI.EQ.0.)UNI=2**(in48)
          ENDIF

          RVEC(IVEC) = UNI
100     CONTINUE

        DO 200 I=1,24
           RVEC(I)=RVEC(LENV+I)
200     continue

        RETURN
        END

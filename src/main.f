! ************************************************
! Forest Light Environmental Simulator (FLiES)    
!
!     Version 2.4.2
!     Hideki Kobayashi
! 
!    as of 2013.11.19
!
!     H. Kobayashi and H. Iwabuchi (2008),
!   Remote Sensing of Environment,112(1), 173-185.
! *************************************************

      implicit none

c      include 'mpif.h'
      include 'common.inc'
      include 'math.inc'
    
!     --- only in main
      integer i,up,ip,iwl,iz, num
      integer knmix,knang,knzext
      parameter (knmix=10, knang=2000, knzext=200)
      real wkd(knkd), ang(knang)
      real re(knmix),Qext_ref(knmix),Qext(knmix),G(knmix),Qabs(knmix)
      real omg(knmix), phs(knmix, knang)
      real rnd, chi, ftau, conv, cv, lai, plai(100), spn

!     --- photon data  ---
      integer nscat, nscata, ichi, ikd
      real x,y,z,ux,uy,uz,w

!     --- function parapeters ----
      real fsin, fcos, facos, r_acos
      real*8 frnd

!     --- set by iparam
      integer flg,cflg,rtype,atype,nkd,nwl
      integer cbnz,ctnz
      real wl0,wls,span(100)
      real taur,ctaur,d
      real th0,ph0,tm,hfov

      real*8 tflx,bflx,dflx,obflx,odflx
      real*8 tpfd,bpfd,dpfd,obpfd,odpfd
      real*8 rflx, rbflx, rdflx
      real*8 scmpf(3,100),scmpp(3,100)
      
      real RF,RQ
      real fpc,fpf,fpv(100)

      real ext(knmix, knzext)  

!     --- set by rparam
      integer np
      integer npl(200),nmix
      integer amode,imode,cmode
      integer stype
      real alb(100)
      real dif,phi,th,ph, tgx, tgy
      real sinf0,cosf0,cosq0,sinq0
      real wq,spcf(400),spcq(400)

      real lr(5,100), lt(5,100), ulr(100), ult(100), str(5,100),sor(100)
      real tlr(5),tlt(5),tstr(5)

      character*81 fname(knmix),rfname

!     MPI related parapeters
      integer Nid, Nprod, ierr
      integer npproc
      
      Nid = 0
      Nprod = 1

!     **** MPI functions  *******************
c      call MPI_init(ierr)
c      call MPI_comm_size(MPI_COMM_WORLD,Nprod,ierr)
c      call MPI_comm_rank(MPI_COMM_WORLD,Nid,ierr)
!     ***************************************

      if(Nid .eq. 0) then
      write(*,*) "*********************************************"
      write(*,*) "3D canopy-atmosphere radiative transfer model" 
      write(*,*) "Forest Light Environmental Simulator (FLiES) "
      write(*,*) "                         by Hideki Kobayashi "
      write(*,*) "        Special thanks to Hironobu Iwabuchi  "
      write(*,*) "               Version 2.4.* Since 2008/5/1  "
      write(*,*) "*********************************************"
      end if

!     ---- Preparation of the initial condition -------

!     Initialize parameters
      write(*,*) "iparam"
      write(*,*) "Parameters initialization..."
      call iparam (
     &     flg,cflg,rtype,atype,nkd,
     &     wl0,wls,nwl,span,
     &     taur,ctaur,d,
     &     th0,ph0,tm,hfov,
     &     tflx,bflx,dflx,obflx,odflx,
     &     tpfd,bpfd,dpfd,obpfd,odpfd,
     &     rflx,rbflx,rdflx,scmpf,scmpp,
     &     cbnz,ctnz,RF,RQ,
     &     fpc,fpf,fpv,ext,Nid
     &  )

!     Initialize math function
      write(*,*) "imath"
      call imath

!     number the initial file read devide (5 for standard input)
      num = 5
c      num = 51
c      open(51, file = "init.inp")
      
!     Read required parameters
      write(*,*) "rparam"
      call rparam(
     &     np,amode,cmode,imode,nwl,npl,nmix,
     &     atype,rtype,stype,cbnz,ctnz,cflg,
     &     dif,th0,ph0,phi,th,ph,
     &     sinf0,cosf0,cosq0,sinq0,
     &     wq,span,wl0,wls,RF,RQ,
     &     d,spcf,spcq,taur, 
     &     ctaur,rfname,fname,alb,
     &     lr,lt,ulr,ult,str,sor,
     &     npproc,Nprod, num)

      if(stype .eq. 2) then
         write(*,*) "G-function" !     G-function LUT
         call igtbl
         write(*,*) "idivspace" !     Initialize space division
         call idivspace
         write(*,*) "ipf"       !   LUT for phase function
         call ipf
      end if

!     fish eye image at the forest floor
!     call local estimation for fish eye image

      if(np .eq. -4) then
         write(*,*) "Input the x, y position of view"
         read(num,*) tgx, tgy
         w = 1.

         write(*,*) "start fish eye"
         call vegfeye(w, tgx, tgy, ichi, ikd)              

         write(*,*) "wdata"
         call wdata(np, nwl, stype, imode, cmode, amode, span, 
     &     RF, RQ, cosq0, wl0, wls, tflx, bflx, dflx, 
     &     tpfd, bpfd, dpfd, scmpf, scmpp, rflx, rbflx, rdflx)
         stop
      end if

!     LAI calculation
      if(np .eq. -5) then
         write(*,*) "input sampling grid (m)"
         write(*,*) "suggested value is 0.1 (m)"
         read(num,*) spn
         call clai(lai, plai, spn, cv) 
         write(*,*) "LAI = ", lai, "; Crown cover (%) = ", cv * 100.
         write(*,*) "H(m) "," LAI_profile"
         do i = int(zmax + 1), 1, -1
            write(*,'(I5,F10.6)') i,plai(i)
         end do
         stop
      end if

c      close(51)

!     Interval to show the progress in standard I/O 
      if(cmode.eq.1)then
         up = 1000
      else
         up = 10000
      end if

      conv = 1.d-8

!     ---- start simulation --------
!     ************MPI function****************
c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
!     ****************************************

      write(*,*) "start simulation"

      if(amode .eq. 2)then      !     without atmosphere 
         
         do ip = 1, npproc

            if(mod(ip, up) .eq. 0)  write(*,*) ip," of ", npproc
            w = 1.0             ! initial weight of photon
            ikd = 0             ! initialization of CDK (Correlated K-dist)

!     selection of beam or diffuse flux
            if(real(frnd()) .gt. dif)then ! beam

               ux = sinq0 * cosf0
               uy = sinq0 * sinf0
               uz = cosq0                  
               if(abs(uz) .lt. 0.0174524) uz = sign(0.0174524, uz)
               nscat = 0
               nscata = nscat

            else                !diffuse

               th = 0.5 * pi + 0.5 * acos(1. - 2. * real(frnd()))
               ph = 2.0 * pi * real(frnd())
                   
               ux = fsin(th) * fcos(ph)
               uy = fsin(th) * fsin(ph)
               uz = fcos(th)
               if(abs(uz) .lt. 0.0174524) uz = sign(0.0174524, uz)
               nscat = 1
               nscata = nscat
            end if

!     initial position (x, y)
            x = xmax * real(frnd())
            y = ymax * real(frnd())

            tflx = tflx + w
            tpfd = tpfd + w * wq

            bflx = bflx + w * (1. - min(real(nscat), 1.))
            bpfd = bpfd + w * wq * (1. - min(real(nscat), 1.))                 
            dflx = dflx + w * min(real(nscat), 1. )
            dpfd = dpfd + w * wq * min(real(nscat), 1.)

!     surface reflection
            if(stype .eq. 1) then !     lambertian

               w = w * alb(1)
               th = 0.5 * facos(1.-2.* real(frnd()))
               ph = 2 * pi * real(frnd())

               ux = fsin(th) * fcos(ph)
               uy = fsin(th) * fsin(ph)
               uz = fcos(th)
               if(abs(uz) .lt. 0.0174524) uz = sign(0.0174524, uz)

            else                !     3-D surface

               do i = 1, nts
                  tlr(i) = lr(i,1)
                  tlt(i) = lt(i,1)
                  tstr(i) = str(i,1)
               end do

!     call the canopy radiation transfer module 

               call canort(x, y, ux, uy, uz, w, wq, cmode, nscat,
     &        tlr, tlt, ulr(1), ult(1), tstr, sor(1), ichi, ikd)

               rflx = rflx + w               
               rbflx = rbflx + w * (1. - min(real(nscata), 1.))                 
               rdflx = rdflx + w * min(real(nscata), 1. )
            end if

         end do
         
      else                      !     With atmosphere
         
         write(*,*) "# of sampling wavelength"
         do iwl = 1, nwl
            write(*, *) iwl, "  of", nwl, " (", npl(iwl),") "
            
!     preparation of the atmospheric optical parameters         
            call  prepatm(imode, wq, RQ, span, spcq, np, npproc, 
     &           wl0, wls, iwl, npl, nkd, wkd, ext, rfname, fname, 
     &           nmix, cflg, re, G, Qabs, Qext, Qext_ref, d, 
     &           cbnz, ctnz, taur, ctaur)

            if(npl(iwl) .eq. 0) then
               write(*, *) "npl error !"
               exit 
            end if

!     loop for single wavelength
            do ip = 1, npl(iwl)

               if(mod(ip, up) .eq. 0)  write(*,*) ip," of ", npl(iwl)
 
               w = 1.0
               x = xmax * real(frnd())
               y = ymax * real(frnd())
               z = zgrd(nz)
               ux = sinq0 * cosf0
               uy = sinq0 * sinf0
               uz = cosq0
               
               ftau = -log(max(1.0e-35, real(frnd())))
               chi = 1.0
               ichi = 1
               iz = nz
               nscat = 0
!     determination of the k-term
               rnd = real(frnd())
               ikd = 1
               do while(rnd .gt. wkd(ikd))
                  ikd = ikd + 1
               end do

               do               ! Monte Carlo core loop
                  call mc1dtrace(w, wq, x, y, z, ux, uy, uz, ftau,
     &                 chi, iz, ikd, nscat, ichi, iwl)

                  nscata = nscat
                  if (w .le. 0.0 .or. iz .gt. nz) exit

                  scmpf(1, iwl) = scmpf(1, iwl) + w
                  scmpf(2, iwl) = scmpf(2, iwl)
     &                 + w * (1. - min(real(nscat), 1.))
                  scmpf(3, iwl) = scmpf(3, iwl)
     &                 + w * min(real(nscat), 1. )

                  scmpp(1, iwl) = scmpp(1, iwl) + w * wq
                  scmpp(2, iwl) = scmpp(2, iwl) 
     &                 + w * wq * (1. - min(real(nscat), 1.))                 
                  scmpp(3, iwl) = scmpp(3, iwl) 
     &                 + w * wq * min(real(nscat), 1.)
                  
!     surface reflection
                  if(stype .eq. 1) then !     lambertian
                     w = w * alb(1)
                     if(w .lt. conv) exit
                     iz = 1

                     th = 0.5 * facos(1.-2.* real(frnd()))
                     ph = 2 * pi * real(frnd())
                     
                     ux = fsin(th) * fcos(ph)
                     uy = fsin(th) * fsin(ph)
                     uz = fcos(th)
                     if(abs(uz) .lt. 0.0174524) uz = sign(0.0174524, uz)
                    
                     rflx = rflx + w               
                     rbflx = rbflx + w * (1. - min(real(nscata), 1.))      
                     rdflx = rdflx + w * min(real(nscata), 1. )
                   
                  else          !     3-D surface                     
                     do i = 1, nts
                        tlr(i) = lr(i,iwl)
                        tlt(i) = lt(i,iwl)
                        tstr(i) = str(i,iwl)
                     end do
      
!     call the canopy radiation transfer module 
                     call canort(x, y, ux, uy, uz, w, wq, cmode, nscat,
     &                    tlr, tlt, ulr(iwl), ult(iwl), tstr, 
     &                    sor(iwl), ichi, ikd)
                     if(w .lt. conv) exit
                     iz = 1
                     rflx = rflx + w               
                     rbflx = rbflx + w * (1. - min(real(nscata), 1.))      
                     rdflx = rdflx + w * min(real(nscata), 1. )
                  end if
               end do           ! single photon loop
            end do              ! photon loop

!     summary of spectral flx
            tflx = tflx + scmpf(1,iwl)
            bflx = bflx + scmpf(2,iwl)
            dflx = dflx + scmpf(3,iwl)

            tpfd = tpfd + scmpp(1,iwl)
            bpfd = bpfd + scmpp(2,iwl)                 
            dpfd = dpfd + scmpp(3,iwl)

         end do                 ! loop for nwl (spectral loop)
      end if

!     ---- End simulation ------
!     ************MPI function****************
c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
!     ****************************************

!     ****************************************
c      call mpisum(tflx, bflx, dflx, 
c     &     tpfd, bpfd, dpfd, scmpf, scmpp, rflx, rbflx, rdflx)
!     ***************************************
      

!     Write the simulated dataset
      if(Nid .eq. 0) then
         write(*,*) "wdata"
         call wdata(np, nwl, stype, imode, cmode, amode, span, 
     &        RF, RQ, cosq0, wl0, wls, tflx, bflx, dflx, 
     &        tpfd, bpfd, dpfd, scmpf, scmpp, rflx, rbflx, rdflx)
      end if

!     ****************************************
c         call MPI_Finalize(ierr)
!     ****************************************

      stop
      end

! space definition
      integer size,ndivs(720),divs(720,300),mdiv
      integer ixmax,iymax,izmax,nzmin	
      parameter(size=300)
      real xmax,ymax,zmax,res,epsi,zmin
      common /space/xmax,ymax,zmax,ndivs,divs,res,epsi,mdiv,
     &      ixmax,iymax,izmax,zmin,nzmin 

! atmospheric common
      integer knz, knkd, knchi, knxr, knyr, knrdc, nlut, nlut1, nang
      parameter (knz=200, knkd=3, knchi=10)
      parameter (knxr=300, knyr=300, knrdc=700)
      parameter (nlut=2000, nlut1=nlut+1)
      integer klayer(0:knz)
      real zgrd(0:knz), zgrd_back(0:knz), chigrd(100) 
      real zgrdm(0:knz), trulut(6, knchi, knz)
      real absg1d(knz, knkd), extt1d(knz, knchi), abst1d(knz)
      real pflut(0:nlut1, knz), salut(0:nlut1, knz)
      real uxrtab(knrdc), uyrtab(knrdc), uzrtab(knrdc)
      double precision prdcF(knxr, knyr, knrdc),prdcQ(knxr, knyr, knrdc)
      real wrr, fsang
      integer nz, nchi, nrdc
      common /catmosppi/prdcF,prdcQ, pflut, salut, fsang, trulut, 
     &	absg1d, extt1d, abst1d, zgrd, chigrd, uxrtab, uyrtab, uzrtab,
     &  wrr, nz, nchi, nrdc, nang, klayer, zgrd_back, zgrdm
      save /catmosppi/

! canopy common

      integer nangc, nth,nph,nobj,iobj(6000),sobj(6000)  
      integer diftype,dirflag,bound
      integer irefl(2,size,size)
      integer mc,mb,mf,nts
      real obj(6000,5), sbar(5)
      real angt(100),angp(100)
      real brf(2,700)
      real brfc(2,700),brfs(2,700),brff(2,700)
      real refl(2,size,size), feye(0:90,0:360), rfeye(0:90,0:360)
      real ap(size,size,100),apf(size,size),aps(size,size)
      real apd(size,size,100),apb(size,size,100),apfd(size,size)
      real apnp(100)
      real tfpr, cfpr, bfpr, ffpr, sfpr 
      real ffdir(size,size),ffdif(size,size)
      real sfdir(size,size),sfdif(size,size)
      real gLAI,u(6),BAD(5),ir,fe,bp1,bp2
      real gtblc(0:180),gtblb(0:180),gtblf(0:180)
      real uxrc(knrdc), uyrc(knrdc), uzrc(knrdc)	
      real rb
      parameter (rb = 0.5)
      common /canopy/diftype,irefl,brf,sbar,
     &  brfc,brfs,brff,refl,dirflag,ap,apf,aps,apd,apb,apfd,apnp,
     &  tfpr, cfpr, bfpr, ffpr, sfpr, feye, rfeye, 
     &  ffdir,ffdif,sfdir,sfdif,bound,nth,nph,angt,angp,
     &  gLAI,u,ir,BAD,fe,nobj,iobj,sobj,obj,gtblc,gtblb,gtblf,
     &  bp1,bp2,
     &  uxrc, uyrc, uzrc, nangc, nts, mc, mb, mf

c scattering angle parameters

      real gmrc(19, 19, 37), gmtc(19, 19, 37)
      real gmrb(19, 19, 37), gmtb(19, 19, 37)
      real gmrf(19, 19, 37), gmtf(19, 19, 37)
      common /scatpahse/gmrc, gmtc, gmrb, gmtb, gmrf, gmtf

c for back trace
      integer nb,bg

      common /bug/nb, bg
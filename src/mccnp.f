!*******************************************************
!     This subroutine calculates photon trajectry
!     in the tubid medium by using Monte Calro technique.
!     
!     Written by H. Kobayashi
!     Last modified 08/04/02
!*******************************************************

      subroutine mccnp(w, wq, x, y, z, ux, uy, uz, cmode, nscat,
     &   tobj, inobj, io, face, lr, lt, str, ichi, ikd)

      implicit none

      include 'common.inc'
      include 'math.inc'

      integer m, inobj, nscat, tsobj
      integer ichi, ikd, cmode
      integer ith,i, cb
      integer io, io1, io2, io12, face
      integer ix, iy, iz

      real x, y, z, ux, uy, uz, uxo, uyo, uzo
      real th, ref, tr, mgn, conv
      real bp, lr, lt, str
      real sgm, tobj(5),tobj12(5), tobjb(5)
      real s, s1, s2, cf, ba, la, par, s12, cf12, rb12
      real w, wq, fd, psh, ssa
      real facos, uzm
      real rnd
      real*8 frnd
    
      conv = 1.d-8
      mgn = 1.d-2
      uzm = 0.0174524
      cf = sbar(iobj(inobj))
      ba = BAD(iobj(inobj))
      la = u(iobj(inobj))
      tsobj = sobj(inobj)

      cf12 = sbar(iobj(inobj))
      rb12 = 1.0
c      cf12 = 0.25
c      rb12 = 0.8

!     define the second canopy area
      tobj12(1) = tobj(1)
      tobj12(2) = tobj(2)
      tobj12(3) = tobj(3) - tobj(4) * (1.- rb12)* min(1, abs(tsobj - 5))
      tobj12(4) = tobj(4) * rb12
      tobj12(5) = tobj(5) * rb12      

!     define the branch dominant region
      tobjb(1) = tobj(1)
      tobjb(2) = tobj(2)
      tobjb(3) = tobj(3) - tobj(4) * (1. - rb) * min(1, abs(tsobj - 5))
      tobjb(4) = tobj(4) * rb
      tobjb(5) = tobj(5) * rb

!     check status:leaf dominant region, branch dominant or irregulary outside
      do 

         if(tsobj .eq. 1) then            
            call cones(s1, x, y, z, ux, uy, uz, tobj, face, io1)
            call cones(s12, x, y, z, ux, uy, uz, tobj12, face, io12)
            call cones(s2, x, y, z, ux, uy, uz, tobjb, face, io2)         
         else if(tsobj .eq. 2) then
            call cyls(s1, x, y, z, ux, uy, uz, tobj, face, io1)
            call cyls(s12, x, y, z, ux, uy, uz, tobj12, face, io12)
            call cyls(s2, x, y, z, ux, uy, uz, tobjb, face, io2)
         else if(tsobj .eq. 3) then
            call elpss(s1, x ,y ,z ,ux ,uy ,uz ,tobj, io1)
            call elpss(s12, x ,y ,z ,ux ,uy ,uz ,tobj12, io12)
            call elpss(s2, x ,y ,z ,ux ,uy ,uz ,tobjb, io2)    
         else if(tsobj .eq. 5) then
            call helpss(s1 ,x ,y ,z ,ux ,uy ,uz ,tobj, face, io1)
            call helpss(s12,x ,y ,z ,ux ,uy ,uz ,tobj12, face, io12)
            call helpss(s2 ,x ,y ,z ,ux ,uy ,uz ,tobjb, face, io2)
         end if
         
         if(io1 .eq. 1) return
      
!     if photon inside branch dominant region
         if(io2 .eq. 0) then

            do 
               rnd = real(frnd()) 
               
!     branch scattering or leaf scatetring
               th = facos(uz)
               ith = int(th * 180. / pi)
               bp = 0.5 + sign(0.5, bp2 - rnd)
               
               sgm = gtblb(ith) * ba * bp 
               sgm = sgm + 4. * cf * gtblc(ith) * la * (1. - bp) 
               sgm = sgm * bp + (sgm / fe) * (1. - bp)
               sgm = max(1.d-5, sgm)               
               ref = str * bp + lr * (1. - bp)
               tr = lt * (1. - bp)
               
               rnd = real(frnd())
               s = -log(max(rnd,1.e-10)) / sgm           
               s = min(s, 0.9d5)

               if(s .lt. s2) then
                  x = x + (s + mgn) * ux
                  y = y + (s + mgn) * uy 
                  z = z + (s + mgn) * uz 

                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

!     recollision loop for shoot clumping effect
                  do 
                     
!     collision forcing parameters
                     ssa = 1.0 - (1.0 - ref - tr) * fe
                     ssa = ssa * (1. - bp)
                     ssa = ssa + (ref + tr) * bp
                     fd = (1. - fe) / ssa
                     fd = fd * (1. - bp)

!     fpar samping (leave or branch)
                     bfpr = bfpr + w * wq * (1. - ssa) * bp 
                     cfpr = cfpr + w * wq * (1. - ssa) * (1. - bp) 
                     
                     ix = int(x * res) + 1
                     iy = int(y * res) + 1 
                     iz = int(z) + 1
 
                     ix = min(ix, size)
                     iy = min(iy, size)

                     par = w * wq * (1. - ssa)
                     ap(ix, iy, iz) = ap(ix, iy, iz) + par * (1. - bp)
                     apb(ix, iy, iz) = apb(ix, iy, iz) + (1. - bp)
     &                    * (1. - min(nscat, 1))
                     apd(ix, iy, iz) = apd(ix, iy, iz) + par * (1. - bp) 
     &                    * min(nscat, 1) 
                     apnp(iz) = apnp(iz) + par * bp 
  
                     w = w * ssa
                     nscat = nscat + 1
                    
                     call vegrroulette(w, epsi) 
                     if(w .lt. conv) return
                     
                     psh = (1. - 4. * cf) * (1. - bp)
                     if(real(frnd()) .ge. psh) exit

                  end do

!     radiance sampling
                  cb = int((1. - bp) + 2. * bp)
                  call vegrad(w, x, y, z, ux, uy, uz, 
     &                 ref, tr, cb, 1.0, fd, ichi, ikd)
                  
!     new photon direction after scattering
                  if(real(frnd()) .ge. fd) then
                     call scatvec(ref,tr, ux, uy, uz, uxo, uyo, uzo, mb)
                     ux = uxo
                     uy = uyo
                     uz = uzo
                     if(abs(uz) .lt. uzm) uz = sign(uzm, uz)
                  end if

!     check status
                  if(tsobj .eq. 1) then
                     call cones
     &                    (s2, x, y, z, ux, uy, uz, tobjb, face, io2)
                  else if(tsobj .eq. 2) then
                     call cyls(
     &                    s2, x, y, z, ux, uy, uz, tobjb, face, io2)  
                  else if(tsobj .eq. 3) then
                     call elpss(s2, x ,y ,z ,ux ,uy ,uz ,tobjb, io2)
                  else if(tsobj .eq. 5) then 
                     call helpss(s2,x,y,z,ux,uy,uz,tobjb,face,io2)
                  end if
                  if(io2. eq. 1) exit
 
               else
                  x = x + (s2 + mgn) * ux
                  y = y + (s2 + mgn) * uy 
                  z = z + (s2 + mgn) * uz

                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

                  exit
               end if
            end do
            
!     if photon inside the canopy dominant region
         else if(io12 .eq. 0) then

            do 

               rnd = real(frnd())            
!     branch scattering or leaf scatetring
               th = facos(uz)
               ith = int(th * 180. / pi)
               bp = 0.5 + sign(0.5, bp1 - rnd)

               sgm = gtblb(ith) * ba * bp 
               sgm = sgm + 4. * cf * gtblc(ith) * la * (1. - bp)  
               sgm = sgm * bp + (sgm / fe) * (1. - bp)
               sgm = max(1.d-5, sgm)
              
               ref = str * bp +lr * (1. - bp)
               tr = lt * (1. - bp)
               
               rnd = real(frnd())
               s = -log(max(rnd,1.e-10)) / sgm
               s = min(s, 0.9d5)
  
               if(s .gt. s2) then
                  x = x + (s2 + mgn) * ux
                  y = y + (s2 + mgn) * uy 
                  z = z + (s2 + mgn) * uz                
                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
                  exit
               else if(s .gt. s12) then
                  x = x + (s12 + mgn) * ux
                  y = y + (s12 + mgn) * uy 
                  z = z + (s12 + mgn) * uz 
                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
                  exit
               else
                  x = x + (s + mgn) * ux
                  y = y + (s + mgn) * uy 
                  z = z + (s + mgn) * uz
                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

                  do 
!     collision forcing parameters
                     ssa = 1.0 - (1.0 - ref - tr) * fe
                     ssa = ssa * (1. - bp)
                     ssa = ssa + (ref + tr) * bp
                     fd = (1. - fe) / ssa
                     fd = fd * (1. - bp)
 
!     fpar samping (leave or branch)                    
                     bfpr = bfpr + w * wq * (1. - ssa) * bp 
                     cfpr = cfpr + w * wq * (1. - ssa) * (1. - bp) 
                     
                     ix = int(x * res) + 1
                     iy = int(y * res) + 1 
                     iz = int(z) + 1

                     ix = min(ix, size)
                     iy = min(iy, size)

                     par = w * wq * (1. - ssa)
                     ap(ix, iy, iz) = ap(ix, iy, iz) + par * (1. - bp) 
                     apb(ix, iy, iz) = apb(ix, iy, iz) 
     &                    + (1. - bp) * (1. - min(nscat, 1))
                     apd(ix, iy, iz) = apd(ix, iy, iz) + par * (1. - bp) 
     &                    * min(nscat, 1) 
                     apnp(iz) = apnp(iz) + par * bp 

                     w = w * ssa
                     nscat = nscat + 1
                  
                     psh = (1. - 4. * cf) * (1. - bp)
                     if(real(frnd()) .ge. psh) exit

                  end do

                  call vegrroulette(w, epsi) 
                  if(w .lt. conv) return

                  cb = int((1. - bp) + 2. * bp)
                  call vegrad(w, x, y, z, ux, uy, uz, 
     &                 lr, lt, cb, 1.0, fd, ichi, ikd)

                  if(real(frnd()) .ge. fd) then
                     call scatvec(ref,tr, ux, uy, uz, uxo, uyo, uzo, mc)
                     ux = uxo
                     uy = uyo
                     uz = uzo 
                     if(abs(uz) .lt. uzm) uz = sign(uzm, uz)
                  end if
 
!     check status
                  if(tsobj .eq. 1) then
                     call cones
     &                    (s12, x, y, z, ux, uy, uz, tobj12, face, io12)
                    call cones
     &                    (s2, x, y, z, ux, uy, uz, tobjb, face, io2)
                  else if(tsobj .eq. 2) then
                     call cyls
     &                    (s12, x, y, z, ux, uy, uz, tobj12, face, io12)
                    call  cyls
     &                    (s2, x, y, z, ux, uy, uz, tobjb, face, io2)   
                  else if(tsobj .eq. 3) then
                     call elpss(s12, x ,y ,z ,ux ,uy ,uz ,tobj12, io12)
                     call elpss(s2, x ,y ,z ,ux ,uy ,uz ,tobjb, io2)
                  else if(tsobj .eq. 5) then 
                     call helpss
     &                    (s12,x ,y ,z ,ux ,uy ,uz,tobj12,face,io12)
                     call helpss
     &                    (s2,x ,y ,z ,ux ,uy ,uz,tobjb,face,io2)
                  end if
!     for unexpected case (photon is outside of canopy)
                    if(io12 .eq. 1) exit
               end if 

            end do
            

!     if photon inside the canopy dominant region
         else if(io1 .eq. 0) then

            do 
               rnd = real(frnd())            
!     branch scattering or leaf scatetring
               th = facos(uz)
               ith = int(th * 180. / pi)
               bp = 0.5 + sign(0.5, bp1 - rnd)

               sgm = gtblb(ith) * ba * bp 
               sgm = sgm + 4. * cf * gtblc(ith) * la * (1. - bp)  
               sgm = sgm * bp + (sgm / fe) * (1. - bp)
               sgm = max(1.d-5, sgm)
              
               ref = str * bp +lr * (1. - bp)
               tr = lt * (1. - bp)
               
               rnd = real(frnd())
               s = -log(max(rnd,1.e-10)) / sgm           
               s = min(s, 0.9d5)
  
               if(s .gt. s12) then

                  x = x + (s12 + mgn) * ux
                  y = y + (s12 + mgn) * uy 
                  z = z + (s12 + mgn) * uz
                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

                  exit
               else if(s .gt. s1) then
                  x = x + (s1 + mgn) * ux
                  y = y + (s1 + mgn) * uy 
                  z = z + (s1 + mgn) * uz 
                  return
               else

                  x = x + (s + mgn) * ux
                  y = y + (s + mgn) * uy 
                  z = z + (s + mgn) * uz
                  x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
                  y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

                  do 
!     collision forcing parameters
                     ssa = 1.0 - (1.0 - ref - tr) * fe
                     ssa = ssa * (1. - bp)
                     ssa = ssa + (ref + tr) * bp
                     fd = (1. - fe) / ssa
                     fd = fd * (1. - bp)
 
!     fpar samping (leave or branch)                    
                     bfpr = bfpr + w * wq * (1. - ssa) * bp 
                     cfpr = cfpr + w * wq * (1. - ssa) * (1. - bp) 
                     
                     ix = int(x * res) + 1
                     iy = int(y * res) + 1 
                     iz = int(z) + 1
 
                     ix = min(ix, size)
                     iy = min(iy, size)

                     par = w * wq * (1. - ssa)
                     ap(ix, iy, iz) = ap(ix, iy, iz) + par * (1. - bp) 
                     apb(ix, iy, iz) = apb(ix, iy, iz) 
     &                    + (1. - bp) * (1. - min(nscat, 1))
                     apd(ix, iy, iz) = apd(ix, iy, iz) + par * (1. - bp) 
     &                    * min(nscat, 1) 
                     apnp(iz) = apnp(iz) + par * bp 

                     w = w * ssa
                     nscat = nscat + 1
                  
                     psh = (1. - 4. * cf) * (1. - bp)
                     if(real(frnd()) .ge. psh) exit

                  end do

                  call vegrroulette(w, epsi) 
                  if(w .lt. conv) return

                  cb = int((1. - bp) + 2. * bp)
                  call vegrad(w, x, y, z, ux, uy, uz, 
     &                 lr, lt, cb, 1.0, fd, ichi, ikd)

                  if(real(frnd()) .ge. fd) then
                     call scatvec(ref,tr, ux, uy, uz, uxo, uyo, uzo, mc)
                     ux = uxo
                     uy = uyo
                     uz = uzo 
                     if(abs(uz) .lt. uzm) uz = sign(uzm, uz)
                  end if
 
!     check status
                  if(tsobj .eq. 1) then
                     call cones
     &                    (s1, x, y, z, ux, uy, uz, tobj, face, io1)
                     call cones
     &                    (s12, x, y, z, ux, uy, uz, tobj12, face, io12)
                  else if(tsobj .eq. 2) then
                     call cyls
     &                    (s1, x, y, z, ux, uy, uz, tobj, face, io1) 
                     call cyls
     &                    (s12, x, y, z, ux, uy, uz, tobj12, face, io12)  
                  else if(tsobj .eq. 3) then
                     call elpss(s1, x ,y ,z ,ux ,uy ,uz ,tobj, io1)
                     call elpss(s12, x ,y ,z ,ux ,uy ,uz ,tobj12, io12)
                  else if(tsobj .eq. 5) then 
                     call helpss(s1 ,x ,y ,z ,ux ,uy ,uz,tobj,face,io1)
                     call helpss
     &                    (s12 ,x ,y ,z ,ux ,uy ,uz,tobj12,face,io12)
                  end if
!     for unexpected case (photon is outside of canopy)
                    if(io1 .eq. 1) return
               end if 
            end do
         end if
      end do
      return
      end

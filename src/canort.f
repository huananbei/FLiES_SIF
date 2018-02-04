!**********************************************
!     Simulator for forest light environment 
!     start MC simulation
!     
!     by H. Kobayashi
!     modified 08/03/26
!     
!**********************************************
! This is the main suvroutine to do a photon tracing.
      
      subroutine canort(x, y, ux, uy, uz, w, wq, cmode, nscat,
     &     lr, lt, ulr, ult, str, sor, ichi, ikd)
      
      implicit none
      include 'common.inc'
      include 'math.inc'
      
!     ivox : id for the current big-voxel
      integer ivox,ichi,ikd,l,idiv,nscat, cmode
      integer io, tio, oface, pface, inobj,i, id
      integer ix, iy

      real x,y,z,ux,uy,uz,x1,y1,z1
      real w,wq
      real lr(*),lt(*),ulr,ult,str(*),sor
      real intv(3),s, so, sp, ts
      real tobj(5)
      real mgn,conv
      
!     wieght minimum limit
      conv = 1.d-8
      
!     initial distance from object
      so = 1.d5
      
      do i = 1, 3
         intv(i) = 50. / res
      end do

!     merginal value
      mgn = 1.d-2
      
      z = zmax
      z = z - mgn
   
!     do while poton exit from canopy space
      do
     
!     determinatin of first input voxel
         x = min(x,xmax)
         x = max(x,0.0)
         y = min(y,ymax)
         y = max(y,0.0)
         z = min(z,zmax)

         x1 = aint(x / intv(1)) 
         y1 = aint(y / intv(2))
         z1 = aint(z / intv(3))  
         ivox = ixmax * iymax * int(z1)
         ivox = ivox + int(y1) * iymax
         ivox = ivox + int(x1) + 1

         x1 = x1 * intv(1)
         y1 = y1 * intv(2)
         z1 = z1 * intv(3)
         
         io = 1

!     check the photon intersection with big-voxel walls
         call planes(sp, x, y, z, ux, uy, uz, x1, y1, z1, pface, intv)         

!     check the photon intersection with objects
         if(ndivs(ivox) .ne. 0)then
            
            so = 1.d5
            do idiv = 1, ndivs(ivox)
               
!     selected object number
               i = divs(ivox, idiv)

               do l = 1, 5
                  tobj(l) = obj(i,l)
               end do
               
               if(sobj(i) .eq. 1)then
                  call cones(ts, x, y, z, ux, uy, uz, tobj, oface, tio)       
               else if((sobj(i) .eq. 2) .or. (sobj(i) .eq. 4))then
                  call cyls(ts, x, y, z, ux, uy, uz, tobj, oface, tio)
               else if(sobj(i) .eq. 3) then
                  call elpss(ts ,x ,y ,z ,ux ,uy ,uz ,tobj, tio)
               else if(sobj(i) .eq. 5) then
                  call helpss(ts ,x ,y ,z ,ux ,uy ,uz ,tobj,oface,tio)
               end if

               if(ts .lt. so) then
                  so = ts
                  inobj = i
                  io = tio
               end if
               if(io .eq. 0) exit

            end do
         end if
               
!     canopy interaction
         if(so .le. sp .or. io .eq. 0) then
            
!     stem interaction
            
            if(sobj(inobj) .eq. 4) then
!     if the photon is in the stem objct by mistake, exit from stem
!     in other case, photon go to the stem surface
               x = x + so * ux 
               y = y + so * uy
               z = z + so * uz 
               
!     Monte Carlo on stem surface
               do l = 1, 5
                  tobj(l) = obj(inobj,l)
               end do

               call mcstm (w, wq, x, y, z, ux, uy, uz, nscat,
     &              cmode, tobj, oface, str, ichi, ikd)

               if(w .lt. conv) return
               x = x + mgn * ux
               y = y + mgn * uy
               z = z + mgn * uz
            else
               
               x = x + ((so + mgn) * ux) * real(io)
               y = y + ((so + mgn) * uy) * real(io)
               z = z + ((so + mgn) * uz) * real(io)               

!     Monte Carlo in canopy media
               id = iobj(inobj)

!         write(*,*) "canort9",x,y,z,x1,y1,z1,ivox,tobj,inobj
               call mccnp(w, wq, x, y, z, ux, uy, uz, cmode, nscat,
     &              tobj, inobj, io, oface, lr(id), lt(id), str(id),
     $              ichi, ikd)
!         write(*,*) "canort10",x,y,z,x1,y1,z1,ivox

               if(w .lt. conv) return
               x = x + mgn * ux
               y = y + mgn * uy
               z = z + mgn * uz      
               x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
               y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
               if(z .ge. zmax) return
            end if
            
!     big-voxel wall interaction
         else
            x = x + (sp + mgn) * ux
            y = y + (sp + mgn) * uy 
            z = z + (sp + mgn) * uz 

            x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
            y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
            
            if(z .le. 0.0) then
!     Monte Carlo in forest floor            
!     forest floor downward flux

               ix = int(x * res) + 1
               iy = int(y * res) + 1

               ix = min(ix, size)
               iy = min(iy, size)

               ffdir(ix, iy) = ffdir(ix, iy) 
     &              + w * wq * (1. - min(nscat, 1))
               ffdif(ix, iy) = ffdif(ix, iy) 
     &              + w * wq * min(nscat, 1)

               call  mcflr(w, wq, x, y, z, ux, uy, uz, nscat, 
     &     cmode, ulr, ult, sor, ichi, ikd) 

               if(w .lt. conv) return
               x = x + mgn * ux
               y = y + mgn * uy
               z = z + mgn * uz
               
            else if(z .ge. zmax) then
!     sky (exit from canopy space)
               return
            else
!     refresh the x, y positin using the boudary condition
               x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
               y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
            end if
         end if
      end do   
      return
      end
      

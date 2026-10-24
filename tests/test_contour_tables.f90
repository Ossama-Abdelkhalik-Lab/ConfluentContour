! Multi-revolution feasibility is decided from compiled-in t_min tables
! (contour_table_multi_rev.inc), with the exact minimum-time search only
! within 1e-4 of the tabulated value. This test computes t_min independently
! (golden-section search on the real time of flight over the revolution
! interval) and requires the right feasible/infeasible decision at
! t = t_min*(1 +- delta), delta = 1e-1 ... 1e-6, for random geometries and
! revolution counts on and between the table slices, both branches.
program test_contour_tables
   use iso_fortran_env, only: real64
   use confluentcontour, only: contour_one, prepare_geometry, transfer_time, geometry_data, pi
   implicit none
   integer,parameter::dp=real64
   integer,parameter::ms(9)=[1,2,3,5,7,12,20,50,200]
   real(dp),parameter::deltas(6)=[1e-1_dp,1e-2_dp,1e-3_dp,1e-4_dp,1e-5_dp,1e-6_dp]
   real(dp),parameter::gr=(sqrt(5._dp)-1)/2
   real(dp)::r1(3),r2(3),u(3),th,rr,tmin,t,va(3),vb(3),a,b,c,d,fc,fd
   type(geometry_data)::g
   integer::k,im,id,sgn,br,dr,st,hf,i,checks
   call random_seed(put=[(97,k=1,64)])
   checks=0
   do k=1,150
      call random_number(u)
      rr=0.05_dp+0.95_dp*u(1);th=1e-3_dp+(pi-2e-3_dp)*u(2);dr=merge(1,-1,u(3)>.5_dp)
      r1=[1._dp,0._dp,0._dp];r2=rr*[cos(th),sin(th),0._dp]
      call prepare_geometry(r1,r2,dr,g,st,hf)
      do im=1,size(ms)
         ! Golden-section minimum of T over the open revolution interval
         ! (normalized units; T is unimodal there).
         a=(2*ms(im)*pi)**2*(1+1e-9_dp);b=(2*(ms(im)+1)*pi)**2*(1-1e-9_dp)
         c=b-gr*(b-a);d=a+gr*(b-a);fc=tm(c);fd=tm(d)
         do i=1,200
            if(fc<fd)then
               b=d;d=c;fd=fc;c=b-gr*(b-a);fc=tm(c)
            else
               a=c;c=d;fc=fd;d=a+gr*(b-a);fd=tm(d)
            endif
            if(b-a<1e-13_dp*b)exit
         enddo
         tmin=min(fc,fd)*g%scale*sqrt(g%scale)
         do id=1,size(deltas)
            do sgn=-1,1,2
               t=tmin*(1+sgn*deltas(id))
               do br=-1,1,2
                  call contour_one(r1,r2,t,dr,va,vb,st,hf,ms(im),br)
                  checks=checks+1
                  if((sgn<0).neqv.(st==5))then
                     print '(a,i4,a,es10.2,a,i2,a,i3,a,es12.5,a,i2)','feasibility error m=',ms(im),' delta=',sgn*deltas(id), &
                        ' branch=',br,' status=',st,' r2=',rr,' d=',dr
                     error stop 'contour t_min classification'
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   print '(a,i0,a)','PASS: ',checks,' near-t_min feasibility decisions match an independent minimum time'
contains
   real(dp) function tm(z)
      real(dp),intent(in)::z
      tm=real(transfer_time(cmplx(z,0._dp,dp),g,ms(im)),dp)
   end function
end program

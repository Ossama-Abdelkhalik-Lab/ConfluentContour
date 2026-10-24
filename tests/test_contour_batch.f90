! Contour: the staged batch path (contour_batch) must reproduce contour_one
! exactly (status, half flag and velocities, bit for bit) on mixed cases:
! zero-revolution elliptic and hyperbolic, multi-revolution both branches
! (feasible, near-minimum and infeasible), invalid inputs, and a batch size
! that is not a multiple of the internal chunk.
program test_contour_batch
   use iso_fortran_env, only: real64
   use confluentcontour, only: contour_one, contour_batch, pi
   implicit none
   integer,parameter::dp=real64,n=5003
   real(dp)::r1(3,n),r2(3,n),tof(n),a1(3,n),a2(3,n),b1(3,n),b2(3,n),u(5),th,rr
   integer::dir(n),nrev(n),br(n),sa(n),ha(n),sb(n),hb(n),i,bad
   call random_seed(put=[(31,i=1,64)])
   do i=1,n
      call random_number(u)
      th=1e-3_dp+(pi-2e-3_dp)*u(1);rr=0.05_dp+0.95_dp*u(2)
      r1(:,i)=[1._dp,0._dp,0._dp];r2(:,i)=rr*[cos(th),sin(th),0.1_dp*(u(3)-.5_dp)]
      dir(i)=merge(1,-1,u(3)>.5_dp);br(i)=merge(1,-1,u(4)>.5_dp)
      nrev(i)=int(4*u(5));if(mod(i,7)==0)nrev(i)=int(20*u(5))
      ! flight times from very short (hyperbolic) to long multi-revolution
      tof(i)=10._dp**(-1.5_dp+3.5_dp*u(4))*(1+nrev(i))
   enddo
   tof(5)=-1;dir(6)=0;nrev(7)=-2;br(8)=3;r2(:,9)=r1(:,9)
   do i=1,n
      call contour_one(r1(:,i),r2(:,i),tof(i),dir(i),a1(:,i),a2(:,i),sa(i),ha(i),nrev(i),br(i))
   enddo
   call contour_batch(n,r1,r2,tof,dir,b1,b2,sb,hb,nrev,br)
   bad=0
   do i=1,n
      if(sa(i)/=sb(i).or.ha(i)/=hb(i))then
         bad=bad+1
      elseif(sa(i)==0)then
         if(any(a1(:,i)/=b1(:,i)).or.any(a2(:,i)/=b2(:,i)))bad=bad+1
      endif
   enddo
   print '(a,i0,a,i0,a,i0,a,i0)','cases ',n,'  status 0: ',count(sa==0),'  status 5: ',count(sa==5),'  mismatches: ',bad
   if(bad/=0)error stop 'contour_batch differs from contour_one'
   print '(a)','PASS: contour_batch matches contour_one bit for bit'
end program

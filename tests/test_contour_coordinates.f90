program test_contour_coordinates
   use confluentcontour
   use ieee_arithmetic
   implicit none
   type(geometry_data)::g
   type(contour_diagnostics)::diag
   real(dp)::p1(3),p2(3),v1(3),v2(3),theta,q,z,t,err,best,va(3),vb(3),left,right,z1,z2
   integer::status,half,m,branch,k,mode
   logical::fast
   p1=[1._dp,0._dp,0._dp]
   ! Direct entry point initializes its cache without an explicit init call.
   do mode=1,2
      fast=mode==1
      do k=2,12
         theta=10._dp**(-k)
         p2=[cos(theta),sin(theta),0._dp]
         call contour_one(p1,p2,2*pi-theta,-1,v1,v2,status,half,fast=fast,diagnostics=diag)
         if(status/=0.or..not.ieee_is_finite(diag%q_upper))error stop 'upper coordinate missing'
         if(diag%relative_tof_residual>5e-13_dp)error stop 'upper residual'
         err=max(norm2(v1-[0._dp,-1._dp,0._dp]),norm2(v2-[sin(theta),-cos(theta),0._dp]))
         if(err>2e-10_dp)error stop 'upper velocity'
         ! Compare the small boundary distance, not the rounded difference of z.
         q=theta*(4*pi-theta)
         if(abs(diag%q_upper/q-1)>2e-10_dp)error stop 'upper distance accuracy'
         call contour_one(p1,p2,theta,1,v1,v2,status,half,fast=fast,diagnostics=diag)
         if(status/=0.or.diag%relative_tof_residual>5e-13_dp)then
            print *, 'lower failure',mode,k,status,diag
            error stop 'lower coordinate'
         endif
         if(norm2(v1-[0._dp,1._dp,0._dp])>2e-10_dp)error stop 'lower velocity'
      enddo
   enddo
   ! Both roots must exist; one is the independent circular-orbit truth.
   theta=1.2_dp;p2=[cos(theta),sin(theta),0._dp]
   do m=1,10
      best=huge(1._dp)
      do branch=-1,1,2
         call contour_one(p1,p2,theta+2*m*pi,1,v1,v2,status,half,m,branch,diagnostics=diag)
         if(status/=0.and.status/=1)error stop 'multi-revolution status'
         if(.not.ieee_is_finite(diag%relative_tof_residual))error stop 'multi-revolution residual finite'
         if(diag%relative_tof_residual>2e-12_dp)error stop 'multi-revolution residual'
         best=min(best,max(norm2(v1-[0._dp,1._dp,0._dp]),norm2(v2-[-sin(theta),cos(theta),0._dp])))
      enddo
      if(best>2e-10_dp)error stop 'multi-revolution circular truth'
      call contour_one(p1,p2,.1_dp,1,v1,v2,status,half,m,-1)
      if(status/=5)error stop 'infeasible revolution count'
   enddo
   ! Tangent root and a just-infeasible request. Locate the minimum using
   ! real TOF samples, independently of the solver's derivative search.
   call prepare_geometry(p1,p2,1,g,status,half)
   left=4*pi*pi+.01_dp;right=16*pi*pi-.01_dp
   do k=1,120
      z1=left+(right-left)/3;z2=right-(right-left)/3
      if(real(transfer_time(cmplx(z1,0._dp,dp),g,1),dp)< &
         real(transfer_time(cmplx(z2,0._dp,dp),g,1),dp))then
         right=z2
      else
         left=z1
      endif
   enddo
   t=real(transfer_time(cmplx((left+right)/2,0._dp,dp),g,1),dp)
   do branch=-1,1,2
      call contour_one(p1,p2,t,1,v1,v2,status,half,1,branch)
      if(status/=1)then
         print *, 'merged status',status
         error stop 'merged root status'
      endif
      call contour_one(p1,p2,t*(1-1e-8_dp),1,v1,v2,status,half,1,branch)
      if(status/=5)error stop 'near minimum infeasible'
   enddo
   ! Exercise both signs through the parabolic boundary and strong hyperbolas.
   call prepare_geometry(p1,p2,1,g,status,half)
   do k=1,6
      z=(-1._dp)**k*10._dp**(-k)
      t=real(transfer_time(cmplx(z,0._dp,dp),g,0),dp)
      call reconstruct(g,z,0,va,vb)
      call contour_one(p1,p2,t,1,v1,v2,status,half,diagnostics=diag)
      if(status/=0.or.diag%relative_tof_residual>2e-12_dp)error stop 'parabolic continuity'
      if(max(norm2(v1-va),norm2(v2-vb))>2e-10_dp)error stop 'parabolic velocities'
   enddo
   call prepare_geometry(p1,p2,-1,g,status,half)
   do k=1,10
      z=(-1._dp)**k*10._dp**(-k)
      t=real(transfer_time(cmplx(z,0._dp,dp),g,0),dp)
      call reconstruct(g,z,0,va,vb)
      call contour_one(p1,p2,t,-1,v1,v2,status,half,diagnostics=diag)
      if(status/=0.or.diag%relative_tof_residual>2e-12_dp)error stop 'long-way parabolic continuity'
      if(max(norm2(v1-va),norm2(v2-vb))>2e-10_dp)error stop 'long-way parabolic velocities'
   enddo
   do k=1,12
      t=10._dp**(-k)
      call contour_one(p1,p2,t,1,v1,v2,status,half,diagnostics=diag)
      if(status/=0.or.diag%relative_tof_residual>2e-12_dp)error stop 'short hyperbolic arc'
      if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))error stop 'hyperbolic finite'
   enddo
   call contour_one(p1,p2,-1._dp,1,v1,v2,status,half,diagnostics=diag)
   if(status/=-1.or.any(ieee_is_finite(v1)))error stop 'invalid time'
   print '(a)', 'PASS: upper/lower coordinates, both modes, 1:10 revolutions, parabolic and short hyperbolic arcs'
end program

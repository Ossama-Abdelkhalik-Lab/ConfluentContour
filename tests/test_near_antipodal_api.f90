program test_near_antipodal_api
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_finite
   use confluentcontour, only: contour_one, pi
   use lambert_api, only: initialize_lambert, solve_lambert, lambert_success, &
      lambert_merged
   use two_body_reference, only: propagate_two_body
   implicit none
   integer, parameter :: dp=real64
   real(dp) :: r1(3),r2(3),v1(3),v2(3),scaled1(3),scaled2(3),contour1(3),contour2(3),arrival(6)
   real(dp) :: theta,tof,delta,position_error,velocity_error,max_position,max_velocity
   real(dp) :: length_unit,mu,time_unit,velocity_unit
   integer :: i,j,m,branch,direction,status,scaled_status,contour_status,half,solved,contour_solved

   call initialize_lambert()
   r1=[1._dp,0._dp,0._dp]
   max_position=0;max_velocity=0;solved=0;contour_solved=0
   do m=0,10
      if(m/=0.and.m/=1.and.m/=2.and.m/=10)cycle
      do i=1,12
         delta=10._dp**(-2._dp-10._dp*real(i-1,dp)/11._dp)
         do j=-1,1,2
            theta=pi+j*delta
            r2=[cos(theta),sin(theta),0._dp]
            tof=2*m*pi+theta
            direction=1
            if(j>0)direction=-1
            do branch=-1,1,2
               if(m==0.and.branch==1)cycle
               call solve_lambert(r1,r2,tof,1._dp,direction,m,branch,v1,v2,status)
               if(.not.accepted(status))then
                  print '(a,3(i0,1x),a,es16.8,a,i0)', &
                     'near-half rejected: m i branch=',m,i,branch,' delta=',delta,' status=',status
                  error stop 'near-half branch status'
               endif
               call propagate_two_body(r1,v1,tof,arrival)
               position_error=norm2(arrival(1:3)-r2)/(1+norm2(r2))
               velocity_error=norm2(arrival(4:6)-v2)/(1+norm2(v2))
               max_position=max(max_position,position_error)
               max_velocity=max(max_velocity,velocity_error)
               if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp)then
                  print '(a,3(i0,1x),a,3(es16.8,1x))', &
                     'near-half miss: m i branch=',m,i,branch, &
                     ' delta position velocity=',delta,position_error,velocity_error
                  error stop 'near-half terminal state'
               endif
               solved=solved+1
               if(i/=2.and.i/=8.and.i/=12)cycle
               length_unit=1.e150_dp;mu=1.e300_dp
               if(i==8)then
                  length_unit=1.e-150_dp;mu=1.e-300_dp
               endif
               if(i==12)then
                  length_unit=1.e6_dp;mu=398600._dp
               endif
               time_unit=exp((3*log(length_unit)-log(mu))/2)
               velocity_unit=exp((log(mu)-log(length_unit))/2)
               call solve_lambert(length_unit*r1,length_unit*r2,tof*time_unit,mu, &
                  direction,m,branch,scaled1,scaled2,scaled_status)
               if(.not.accepted(scaled_status).or..not.all(ieee_is_finite(scaled1)).or. &
                  .not.all(ieee_is_finite(scaled2)))error stop 'scaled near-half status'
               if(norm2(scaled1/velocity_unit-v1)>1.e-10_dp.or. &
                  norm2(scaled2/velocity_unit-v2)>1.e-10_dp)then
                  print '(a,3(i0,1x),a,es16.8)', &
                     'scaled near-half miss: m i branch=',m,i,branch,' delta=',delta
                  error stop 'scaled near-half velocity'
               endif
               call contour_one(r1,r2,tof,direction,contour1,contour2, &
                  contour_status,half,m,branch)
               if(.not.accepted(contour_status).or.half<0) &
                  error stop 'contour near-half status'
               call propagate_two_body(r1,contour1,tof,arrival)
               position_error=norm2(arrival(1:3)-r2)/(1+norm2(r2))
               velocity_error=norm2(arrival(4:6)-contour2)/(1+norm2(contour2))
               if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp)then
                  print '(a,3(i0,1x),a,3(es16.8,1x))', &
                     'contour near-half miss: m i branch=',m,i,branch, &
                     ' delta position velocity=',delta,position_error,velocity_error
                  error stop 'contour near-half terminal state'
               endif
               contour_solved=contour_solved+1
            enddo
         enddo
      enddo
   enddo
   r2=[-1._dp,0._dp,0._dp]
   do m=0,10
      if(m/=0.and.m/=1.and.m/=10)cycle
      tof=(2*m+1)*pi
      do branch=-1,1,2
         if(m==0.and.branch==1)cycle
         call solve_lambert(r1,r2,tof,1._dp,1,m,branch,v1,v2,status,[0._dp,0._dp,1._dp])
         if(.not.accepted(status))error stop 'exact antipodal branch status'
         call propagate_two_body(r1,v1,tof,arrival)
         position_error=norm2(arrival(1:3)-r2)/(1+norm2(r2))
         velocity_error=norm2(arrival(4:6)-v2)/(1+norm2(v2))
         if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp) &
            error stop 'exact antipodal terminal state'
         solved=solved+1
      enddo
   enddo
   print '(a,i0,a,i0,a,es12.4,a,es12.4)', 'PASS: near-antipodal API branches=',solved, &
      ' contour branches=',contour_solved, &
      ' max API position=',max_position,' max API velocity=',max_velocity
contains
   pure logical function accepted(code)
      integer, intent(in) :: code
      accepted=code==lambert_success.or.code==lambert_merged
   end function
end program

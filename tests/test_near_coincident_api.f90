program test_near_coincident_api
   use iso_fortran_env, only: real64
   use confluentcontour, only: contour_one, pi
   use lambert_api, only: initialize_lambert, solve_lambert, lambert_success
   use two_body_reference, only: propagate_two_body
   implicit none
   integer, parameter :: dp=real64
   real(dp) :: r1(3),r2(3),v1(3),v2(3),w1(3),w2(3),arrival(6)
   real(dp) :: delta,angle,tof,position_error,velocity_error,max_position,max_velocity
   real(dp) :: vt,vr,semimajor,u
   integer :: edge,i,direction,status,half,solved

   call initialize_lambert()
   r1=[1._dp,0._dp,0._dp]
   max_position=0;max_velocity=0;solved=0
   do edge=1,2
      do i=1,40
         delta=10._dp**(-2._dp-10._dp*real(i-1,dp)/39._dp)
         direction=1;angle=delta
         if(edge==2)then
            direction=-1;angle=2*pi-delta
         endif
         r2=[cos(angle),sin(angle),0._dp];tof=angle
         call solve_lambert(r1,r2,tof,1._dp,direction,0,-1,v1,v2,status)
         if(.not.accepted(status))then
            print '(a,2(i0,1x),a,es16.8,a,i0)', &
               'API near-coincident status: edge i=',edge,i,' delta=',delta,' status=',status
            error stop 'API near-coincident status'
         endif
         call check_terminal(r1,r2,tof,v1,v2,position_error,velocity_error)
         max_position=max(max_position,position_error)
         max_velocity=max(max_velocity,velocity_error)
         if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp) &
            error stop 'API near-coincident terminal state'
         solved=solved+1
         call contour_one(r1,r2,tof,direction,w1,w2,status,half,0,-1)
         if(.not.accepted(status).or.half<0)error stop 'contour near-coincident status'
         call check_terminal(r1,r2,tof,w1,w2,position_error,velocity_error)
         max_position=max(max_position,position_error)
         max_velocity=max(max_velocity,velocity_error)
         if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp) &
            error stop 'contour near-coincident terminal state'
         solved=solved+1
      enddo
   enddo
   do i=1,20
      u=real(i-1,dp)/19
      delta=10._dp**(-2._dp-10*u)
      vt=0.8_dp+0.4_dp*u;vr=0.2_dp*(2*real(mod(7*i,23),dp)/23-1)
      semimajor=1/(2-vt*vt-vr*vr)
      tof=2*pi*semimajor**1.5_dp-delta
      call propagate_two_body(r1,[vr,vt,0._dp],tof,arrival)
      r2=arrival(1:3)
      if(r2(2)>=0)error stop 'eccentric near-period orientation'
      call solve_lambert(r1,r2,tof,1._dp,-1,0,-1,v1,v2,status)
      if(.not.accepted(status))error stop 'API eccentric near-period status'
      call check_terminal(r1,r2,tof,v1,v2,position_error,velocity_error)
      max_position=max(max_position,position_error)
      max_velocity=max(max_velocity,velocity_error)
      if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp) &
         error stop 'API eccentric near-period terminal state'
      solved=solved+1
      call contour_one(r1,r2,tof,-1,w1,w2,status,half,0,-1)
      if(.not.accepted(status).or.half<0)error stop 'contour eccentric near-period status'
      call check_terminal(r1,r2,tof,w1,w2,position_error,velocity_error)
      max_position=max(max_position,position_error)
      max_velocity=max(max_velocity,velocity_error)
      if(position_error>1.e-11_dp.or.velocity_error>1.e-11_dp) &
         error stop 'contour eccentric near-period terminal state'
      solved=solved+1
   enddo
   print '(a,i0,a,es12.4,a,es12.4)', 'PASS: near-coincident solutions=',solved, &
      ' max position=',max_position,' max velocity=',max_velocity
contains
   pure logical function accepted(code)
      integer,intent(in)::code
      accepted=code==lambert_success
   end function
   subroutine check_terminal(start,finish,time,departure,arrival_velocity,er,ev)
      real(dp),intent(in)::start(3),finish(3),time,departure(3),arrival_velocity(3)
      real(dp),intent(out)::er,ev
      call propagate_two_body(start,departure,time,arrival)
      er=norm2(arrival(1:3)-finish)/(1+norm2(finish))
      ev=norm2(arrival(4:6)-arrival_velocity)/(1+norm2(arrival_velocity))
   end subroutine
end program

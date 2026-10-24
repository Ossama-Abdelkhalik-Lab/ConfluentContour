program test_lambert_api
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_nan
   use confluentcontour, only: pi
   use lambert_api
   implicit none
   integer, parameter :: dp=real64
   real(dp) :: r1(3),r2(3),v1(3),v2(3),theta,large,mu,time_scale
   real(dp) :: all_v1(3,3),all_v2(3,3)
   integer :: status,all_status(3),count

   theta=1.2_dp;r1=[1._dp,0._dp,0._dp];r2=[cos(theta),sin(theta),0._dp]
   call solve_lambert(r1,r2,theta,1._dp,1,0,-1,v1,v2,status)
   if(status/=lambert_not_initialized)error stop 'initialization status'
   call initialize_lambert()
   call solve_lambert(r1,r2,theta,1._dp,1,0,-1,v1,v2,status)
   if(.not.accepted(status))error stop 'ordinary status'
   if(norm2(v1-[0._dp,1._dp,0._dp])>2.e-11_dp)error stop 'ordinary departure'
   if(norm2(v2-[-sin(theta),cos(theta),0._dp])>2.e-11_dp)error stop 'ordinary arrival'

   call solve_lambert(r1,r2,theta/3,9._dp,1,0,-1,v1,v2,status)
   if(.not.accepted(status))error stop 'mu status'
   if(norm2(v1-[0._dp,3._dp,0._dp])>2.e-10_dp)error stop 'mu departure'
   if(norm2(v2-3*[-sin(theta),cos(theta),0._dp])>2.e-10_dp)error stop 'mu arrival'

   large=1.e150_dp;mu=1.e300_dp;time_scale=1.e75_dp
   call solve_lambert(large*r1,large*r2,theta*time_scale,mu,1,0,-1,v1,v2,status)
   if(.not.accepted(status))error stop 'large-scale status'
   if(norm2(v1/time_scale-[0._dp,1._dp,0._dp])>2.e-10_dp)error stop 'large-scale departure'
   if(norm2(v2/time_scale-[-sin(theta),cos(theta),0._dp])>2.e-10_dp)error stop 'large-scale arrival'

   large=1.e-150_dp;mu=1.e-300_dp;time_scale=1.e-75_dp
   call solve_lambert(large*r1,large*r2,theta*time_scale,mu,1,0,-1,v1,v2,status)
   if(.not.accepted(status))error stop 'small-scale status'
   if(norm2(v1/time_scale-[0._dp,1._dp,0._dp])>2.e-10_dp)error stop 'small-scale departure'
   if(norm2(v2/time_scale-[-sin(theta),cos(theta),0._dp])>2.e-10_dp)error stop 'small-scale arrival'

   r2=[-1._dp,0._dp,0._dp]
   call solve_lambert(r1,r2,pi,1._dp,1,0,-1,v1,v2,status)
   if(status/=lambert_plane_required.or..not.all(ieee_is_nan(v1)))error stop 'missing plane status'
   call solve_lambert(r1,r2,pi,1._dp,1,0,-1,v1,v2,status,[0._dp,0._dp,1._dp])
   if(.not.accepted(status))then
      print '(a,i0)', 'half-revolution status ',status
      error stop 'half-revolution status'
   endif
   if(norm2(v1-[0._dp,1._dp,0._dp])>2.e-10_dp)error stop 'half-revolution departure'
   if(norm2(v2-[0._dp,-1._dp,0._dp])>2.e-10_dp)error stop 'half-revolution arrival'
   call solve_lambert(r1,r2,pi,1._dp,-1,0,-1,v1,v2,status,[0._dp,0._dp,1._dp])
   if(.not.accepted(status))error stop 'reverse half-revolution status'
   if(norm2(v1-[0._dp,-1._dp,0._dp])>2.e-10_dp)error stop 'reverse half-revolution departure'

   r2=[-cos(1.e-6_dp),sin(1.e-6_dp),0._dp]
   call solve_lambert(r1,r2,pi-1.e-6_dp,1._dp,1,0,-1,v1,v2,status)
   if(.not.accepted(status))error stop 'near-antipodal status'
   if(norm2(v1-[0._dp,1._dp,0._dp])>1.e-10_dp)error stop 'near-antipodal departure'
   if(norm2(v2-[-sin(pi-1.e-6_dp),cos(pi-1.e-6_dp),0._dp])>1.e-10_dp) &
      error stop 'near-antipodal arrival'
   r2=[-1._dp,0._dp,0._dp]

   call solve_lambert(r1,r1,pi,1._dp,1,0,-1,v1,v2,status)
   if(status/=lambert_nonunique)error stop 'coincident status'
   call solve_lambert(r1,2*r1,pi,1._dp,1,0,-1,v1,v2,status,[0._dp,0._dp,1._dp])
   if(status/=lambert_unsupported_geometry)error stop 'same-ray status'
   call solve_lambert(r1,r2,pi,0._dp,1,0,-1,v1,v2,status,[0._dp,0._dp,1._dp])
   if(status/=lambert_invalid_input)error stop 'invalid mu status'
   call solve_lambert([0._dp,0._dp,0._dp],r2,pi,1._dp,1,0,-1,v1,v2,status)
   if(status/=lambert_invalid_input)error stop 'zero radius status'
   call solve_lambert(r1,r2,pi,1._dp,-huge(1)-1,0,-1,v1,v2,status)
   if(status/=lambert_invalid_input)error stop 'invalid direction status'

   r2=[cos(theta),sin(theta),0._dp]
   call solve_all_lambert(r1,r2,10._dp,1._dp,1,1,all_v1,all_v2,all_status,count)
   if(count/=3.or..not.all([accepted(all_status(1)),accepted(all_status(2)),accepted(all_status(3))]))then
      print '(a,4(i0,1x))', 'enumeration count/status ',count,all_status
      error stop 'branch enumeration'
   endif
   if(norm2(all_v1(:,2)-all_v1(:,3))<1.e-3_dp)error stop 'distinct branches collapsed'
   print '(a)', 'PASS: dimensional scaling, half-revolution plane, statuses, and branch enumeration'
contains
   pure logical function accepted(code)
      integer, intent(in) :: code
      accepted=code==lambert_success.or.code==lambert_merged
   end function
end program

! Public dimensional interface around the contour solver.
module lambert_api
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use confluentcontour, only: contour_initialize, contour_one, geometry_data, prepare_geometry, solve_z, reconstruct
   implicit none
   private
   integer, parameter :: dp=real64
   integer, parameter, public :: lambert_success=0
   integer, parameter, public :: lambert_merged=1
   integer, parameter, public :: lambert_no_solution=5
   integer, parameter, public :: lambert_invalid_input=-1
   integer, parameter, public :: lambert_plane_required=-2
   integer, parameter, public :: lambert_nonunique=-3
   integer, parameter, public :: lambert_numerical_failure=-4
   integer, parameter, public :: lambert_unsupported_geometry=-5
   integer, parameter, public :: lambert_not_initialized=-6
   integer, parameter, public :: lambert_accuracy_unresolved=-7
   real(dp), parameter :: near_antipodal_sine_floor=1.e-4_dp
   logical, save :: initialized=.false.
   public :: initialize_lambert, solve_lambert, solve_all_lambert
contains
   subroutine initialize_lambert()
      call contour_initialize()
      initialized=.true.
   end subroutine

   subroutine solve_lambert(r1,r2,tof,mu,direction,revolutions,branch,v1,v2,status,plane_normal)
      real(dp), intent(in) :: r1(3),r2(3),tof,mu
      integer, intent(in) :: direction,revolutions,branch
      real(dp), intent(out) :: v1(3),v2(3)
      integer, intent(out) :: status
      real(dp), optional, intent(in) :: plane_normal(3)
      real(dp) :: p1(3),p2(3),e1(3),e2(3),normal(3),cross(3)
      real(dp) :: scale,semi,log_time,log_velocity,velocity_unit,flight_time
      real(dp) :: w1(3),w2(3),root,normal_size,alignment,q_upper
      type(geometry_data) :: geometry
      integer :: legacy_status,half

      v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1
      if(.not.initialized)then
         status=lambert_not_initialized;return
      endif
      status=lambert_invalid_input
      if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)))return
      if(.not.ieee_is_finite(tof).or..not.ieee_is_finite(mu))return
      if(tof<=0.or.mu<=0.or.(direction/=1.and.direction/=-1).or.revolutions<0.or. &
         (branch/=1.and.branch/=-1))return
      scale=max(norm2(r1),norm2(r2))
      if(.not.ieee_is_finite(scale).or.scale<=0)return
      p1=r1/scale;p2=r2/scale
      if(norm2(p1)<=0.or.norm2(p2)<=0)then
         status=lambert_invalid_input;return
      endif
      if(all(r1==r2))then
         status=lambert_nonunique;return
      endif
      semi=(norm2(p1)+norm2(p2)+norm2(p2-p1))/2
      if(.not.ieee_is_finite(semi).or.semi<=0)then
         status=lambert_numerical_failure;return
      endif
      p1=p1/semi;p2=p2/semi
      log_time=log(tof)+0.5_dp*log(mu)-1.5_dp*(log(scale)+log(semi))
      log_velocity=0.5_dp*(log(mu)-log(scale)-log(semi))
      status=lambert_numerical_failure
      if(log_time>=log(huge(tof)).or.log_time<=log(tiny(tof)))return
      if(log_velocity>=log(huge(tof)).or.log_velocity<=log(tiny(tof)))return
      flight_time=exp(log_time);velocity_unit=exp(log_velocity)
      if(.not.ieee_is_finite(flight_time).or.flight_time<=0)return
      if(.not.ieee_is_finite(velocity_unit).or.velocity_unit<=0)return
      e1=p1/norm2(p1);e2=p2/norm2(p2)
      cross=cross3(e1,e2)
      if(norm2(cross)<=32*epsilon(1._dp))then
         if(.not.present(plane_normal))then
            status=lambert_plane_required;return
         endif
         normal=plane_normal
         status=lambert_invalid_input
         if(.not.all(ieee_is_finite(normal)))return
         normal_size=norm2(normal)
         if(normal_size<=0.or..not.ieee_is_finite(normal_size))return
         normal=normal/normal_size
         alignment=max(abs(dot_product(normal,e1)),abs(dot_product(normal,e2)))
         if(alignment>64*epsilon(1._dp))return
         if(norm2(e1+e2)>32*epsilon(1._dp))then
            status=lambert_unsupported_geometry;return
         endif
         geometry%r1=norm2(p1);geometry%r2=norm2(p2);geometry%scale=1
         geometry%e1=e1;geometry%e2=e2
         geometry%t1=cross3(direction*normal,e1)
         geometry%t2=cross3(direction*normal,e2)
         geometry%cs=0;geometry%ss=1;geometry%a=0;geometry%b=0
         geometry%y0=geometry%r1+geometry%r2;geometry%yend=geometry%y0
         call solve_z(geometry,flight_time,revolutions,branch,root,legacy_status,q_upper)
         if(legacy_status==0.or.legacy_status==1)then
            call reconstruct(geometry,root,revolutions,w1,w2,q_upper=q_upper)
         else
            call translate_status(legacy_status,status)
            return
         endif
      else
         if(dot_product(e1,e2)<0.and.norm2(cross)<near_antipodal_sine_floor)then
            call prepare_geometry(p1,p2,direction,geometry,legacy_status,half)
            if(legacy_status/=0.or.half<0)then
               status=lambert_accuracy_unresolved;return
            endif
            call solve_z(geometry,flight_time/geometry%scale/sqrt(geometry%scale), &
               revolutions,branch,root,legacy_status,q_upper)
            if(legacy_status==0.or.legacy_status==1) &
               call reconstruct(geometry,root,revolutions,w1,w2,stable_dc=.true.,q_upper=q_upper)
         else
            call contour_one(p1,p2,flight_time,direction,w1,w2,legacy_status,half,revolutions,branch)
         endif
         if(half<0)then
            status=lambert_numerical_failure;return
         endif
      endif
      call translate_status(legacy_status,status)
      if(status/=lambert_success.and.status/=lambert_merged)return
      if(.not.all(ieee_is_finite(w1)).or..not.all(ieee_is_finite(w2)))then
         status=lambert_numerical_failure;return
      endif
      v1=w1*velocity_unit;v2=w2*velocity_unit
      if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))then
         v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1
         status=lambert_numerical_failure
      endif
   end subroutine

   subroutine solve_all_lambert(r1,r2,tof,mu,direction,max_revolutions,v1,v2,status,count,plane_normal)
      real(dp), intent(in) :: r1(3),r2(3),tof,mu
      integer, intent(in) :: direction,max_revolutions
      real(dp), intent(out) :: v1(:,:),v2(:,:)
      integer, intent(out) :: status(:),count
      real(dp), optional, intent(in) :: plane_normal(3)
      integer :: m,slot
      v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1
      status=lambert_invalid_input;count=0
      if(max_revolutions<0.or.max_revolutions>(huge(max_revolutions)-1)/2)return
      if(size(v1,1)/=3.or.size(v2,1)/=3)return
      if(size(v1,2)<2*max_revolutions+1.or.size(v2,2)<2*max_revolutions+1)return
      if(size(status)<2*max_revolutions+1)return
      call solve_lambert(r1,r2,tof,mu,direction,0,-1,v1(:,1),v2(:,1),status(1),plane_normal)
      if(is_solution(status(1)))count=count+1
      do m=1,max_revolutions
         slot=2*m
         call solve_lambert(r1,r2,tof,mu,direction,m,-1,v1(:,slot),v2(:,slot),status(slot),plane_normal)
         call solve_lambert(r1,r2,tof,mu,direction,m,1,v1(:,slot+1),v2(:,slot+1),status(slot+1),plane_normal)
         if(is_solution(status(slot)))count=count+1
         if(is_solution(status(slot+1)).and.status(slot+1)/=lambert_merged)count=count+1
      enddo
   end subroutine

   pure logical function is_solution(status)
      integer, intent(in) :: status
      is_solution=status==lambert_success.or.status==lambert_merged
   end function

   subroutine translate_status(legacy,status)
      integer, intent(in) :: legacy
      integer, intent(out) :: status
      select case(legacy)
       case(0);status=lambert_success
       case(1);status=lambert_merged
       case(5);status=lambert_no_solution
       case default;status=lambert_numerical_failure
      end select
   end subroutine

   pure function cross3(a,b) result(c)
      real(dp), intent(in) :: a(3),b(3)
      real(dp) :: c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function
end module

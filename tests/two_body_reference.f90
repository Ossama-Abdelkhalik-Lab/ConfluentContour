! Test-only propagator: quad working precision reduces long-arc reference drift.
module two_body_reference
   use iso_fortran_env, only: real64, real128
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   integer, parameter :: dp=real64, qp=real128
   public :: propagate_two_body
contains
   pure function dynamics(y) result(dy)
      real(qp),intent(in) :: y(6)
      real(qp) :: dy(6),radius
      radius=norm2(y(1:3))
      dy(1:3)=y(4:6)
      dy(4:6)=-y(1:3)/radius**3
   end function

   subroutine propagate_two_body(r,v,tof,arrival)
      real(dp),intent(in) :: r(3),v(3),tof
      real(dp),intent(out) :: arrival(6)
      real(qp) :: y(6),candidate(6),lower(6),k1(6),k2(6),k3(6),k4(6),k5(6),k6(6),k7(6)
      real(qp) :: t,h,err,scale,factor,flight_time
      integer :: attempt
      y=[real(r,qp),real(v,qp)];t=0;flight_time=real(tof,qp);h=flight_time/100
      do attempt=1,1000000
         h=min(h,flight_time-t)
         if(h<=0)exit
         k1=dynamics(y)
         k2=dynamics(y+h*(k1/5))
         k3=dynamics(y+h*(3*k1/40+9*k2/40))
         k4=dynamics(y+h*(44*k1/45-56*k2/15+32*k3/9))
         k5=dynamics(y+h*(19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729))
         k6=dynamics(y+h*(9017*k1/3168-355*k2/33+46732*k3/5247+49*k4/176-5103*k5/18656))
         candidate=y+h*(35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84)
         k7=dynamics(candidate)
         lower=y+h*(5179*k1/57600+7571*k3/16695+393*k4/640-92097*k5/339200+187*k6/2100+k7/40)
         scale=maxval(abs(candidate-lower)/(1.e-18_qp+1.e-18_qp*max(abs(y),abs(candidate))))
         if(.not.ieee_is_finite(scale)) error stop 'nonfinite propagation step'
         err=max(scale,1.e-24_qp)
         if(err<=1)then
            y=candidate;t=t+h
            if(t>=flight_time)exit
         endif
         factor=max(0.2_qp,min(5._qp,0.9_qp*err**(-0.2_qp)))
         h=h*factor
         if(h<=spacing(t)) error stop 'propagation step underflow'
      enddo
      if(t<flight_time) error stop 'propagation step limit'
      arrival=real(y,dp)
   end subroutine
end module

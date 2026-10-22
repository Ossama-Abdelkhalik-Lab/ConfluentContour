!> Simple contour Lambert solver: the core of confluentcontour in one short routine.
!!
!! REFERENCE
!!   Negrete and Abdelkhalik (2024), DOI 10.2514/1.G008499: the time-of-flight
!!   root is the pole of 1/(T-t), located by a ratio of contour integrals.
!!   solvers/confluentcontour/ConfluentContourSolver.md, Sections 3.1, 3.3(d), 3.4, 4.9.
!!
!! ALGORITHM
!!   1. Half-anomaly time function (one sin/cos pair, one square root):
!!        T(v) = sqrt(y)*(y*F + a),  y = y0 + b*(1-cos v),
!!        F    = (v + m*pi - sin v cos v)/(sqrt2*sin(v)**3),
!!      v = sqrt(z)/2 - m*pi in (0,pi) on elliptic arcs, w = sqrt(-z)/2 on
!!      hyperbolic arcs (sinh/cosh in place of sin/cos).
!!   2. T', T'', T''' from the exact ODE  sin(v)*F' + 3*cos(v)*F = sqrt2.
!!   3. Closed-form starters: cubic Hermite models between the parabolic
!!      limit and the exact pole (elliptic) or branch point (hyperbolic).
!!   4. Confluent contour steps: the zero-radius limit of the paper's 4-node
!!      contour ratio (Householder order 3), kept inside the admissible
!!      interval by bisection.
!!   5. Multiple revolutions: Halley on T'(v) for the minimum time, then each
!!      branch on its own side of the minimum.
!!   6. Very short short-way hyperbolas (t < 0.1*tp): Newton in sqrt(y),
!!      since y is lost to cancellation near the branch point y = 0.
!!   No tables, sampled contours or safeguards; see confluentcontour for those.
!!
!! UNITS AND CONVENTIONS
!!   Same as contour_one: mu = 1, direction = +1 (along r1 x r2) or -1,
!!   nrev full revolutions, branch = -1 (low z) or +1 (high z) when nrev > 0.
!!
!! STATUS
!!    0: solution found
!!    5: flight time below the multi-revolution minimum
!!   -1: invalid inputs
!!   -2: nearly collinear geometry (undefined transfer plane)
!!   -4: root not converged or nonfinite result
!!   half: -1 degenerate, 0 ordinary, +1 near-collinear.
module contour_simple
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: contour_simple_one
   integer,parameter :: dp=real64
   real(dp),parameter :: pi=acos(-1.0_dp),sqrt2=sqrt(2.0_dp)
   real(dp),parameter :: step_tol=1e-7_dp   ! relative step that ends the iteration
   integer,parameter  :: max_steps=12
   ! Parabolic series of (v-sin(v)cos(v))/sin(v)**3 in v**2 (highest first),
   ! used for |v|<0.3 on zero-revolution arcs where the closed form cancels.
   real(dp),parameter :: pser(10)=[115665628927._dp/2601853536522240000._dp, &
      3821869001._dp/10751460894720000._dp,133541._dp/48117888000._dp,763967._dp/36324288000._dp, &
      1393481._dp/9081072000._dp,1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5,2._dp/3]

contains

   subroutine contour_simple_one(r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
   real(dp),intent(in)          :: r1(3),r2(3),tof
   integer,intent(in)           :: direction
   real(dp),intent(out)         :: v1(3),v2(3)
   integer,intent(out)          :: status,half
   integer,optional,intent(in)  :: nrev,branch
   real(dp) :: scale,rn1,rn2,e1(3),e2(3),h(3),hn,t1v(3),t2v(3),cs,ss,p,q,a,b,y0,yend, &
               t,tp,tv2,x,lo,hi,phi,sig,tmin,d(0:3),f(0:3),k0,kpi,c,s,y,dc,ratio, &
               vr1,vr2,vt1,vt2,wb,psi,yy,sw,dt
   integer  :: m,br,k,it

   v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1;half=0;status=-1
   m=0;br=-1
   if(present(nrev))m=nrev
   if(present(branch))br=branch
   if(abs(direction)/=1.or.m<0.or.abs(br)/=1.or..not.(tof>0).or..not.ieee_is_finite(tof))return

   !======= Geometry, normalized by the larger radius =======
   rn1=norm2(r1);rn2=norm2(r2)
   if(.not.(rn1>0.and.rn2>0).or..not.ieee_is_finite(rn1+rn2))return
   scale=max(rn1,rn2);e1=r1/rn1;e2=r2/rn2;rn1=rn1/scale;rn2=rn2/scale
   h=cross(e1,e2);hn=norm2(h)
   status=-2
   if(hn<32*epsilon(hn))then
      half=-1;return
   endif
   if(hn<1e-4_dp)half=1
   h=h*(direction/hn);t1v=cross(h,e1);t2v=cross(h,e2)
   cs=norm2(e1+e2)/2;ss=norm2(e2-e1)/2        ! |cos(theta/2)|, |sin(theta/2)|
   p=sqrt(rn1*rn2);a=direction*sqrt2*p*cs;b=sqrt2*a
   q=(rn1-rn2)/(sqrt(rn1)+sqrt(rn2))
   yend=q*q+2*p*ss*ss/(1+cs)                  ! y at v=pi (elliptic, any direction)
   y0=yend;if(direction<0)y0=q*q+2*p*(1+cs)   ! y at v=0
   t=tof/(scale*sqrt(scale))
   status=-4

   !======= Starters and admissible interval =======
   k=1;lo=0;hi=pi;yy=-1
   if(m==0)then
      tp=sqrt(y0)*(y0*sqrt2/3+a)                                 ! parabolic time
      tv2=y0*sqrt(y0)/(5*sqrt2)+b*(sqrt(y0)/sqrt2+a/(2*sqrt(y0)))/2 ! dT/d(v**2) at 0
      if(abs(t-tp)<=4*epsilon(t)*tp)then
         x=0                                                     ! parabola
      elseif(t>tp)then
         ! Ellipse: cubic Hermite for v**2 in phi=(tp/t)**(1/3), between
         ! v**2=(t-tp)/tv2 at phi=1 and pi-v=phi*(K/tp)**(1/3) at phi=0,
         ! with the exact upper pole T ~ pi*yend**1.5/(sqrt2*(pi-v)**3).
         phi=(tp/t)**(1._dp/3)
         x=hermite(phi,pi*pi,-2*pi*sqrt(yend)*(pi/sqrt2/tp)**(1._dp/3),0._dp,-3*tp/tv2)
         if(.not.(x>0))x=(t-tp)/tv2         ! model undershoot: parabolic limit
         x=min(sqrt(x),pi*(1-1e-3_dp))
      else
         ! Hyperbola in w=sqrt(-z)/2: w**2=(tp-t)/tv2 near the parabola.
         k=-1;hi=huge(hi)
         if(b>0)then
            ! Short way: T -> 0 at the branch point y(wb)=0.
            wb=2*asinh(sqrt(y0/(2*b)));hi=wb
            if(t<0.1_dp*tp)then
               ! Newton in Y=sqrt(y): T=Y*(Y**2*F+a) with the exact
               ! w(Y)=2*asinh(sqrt((y0-Y**2)/(2b))). Either term of
               ! T = a*Y + F*Y**3 alone overestimates Y; start from the smaller.
               call f_derivs(wb,-1,0,f)
               yy=min(t/a,(t/f(0))**(1._dp/3),0.9_dp*sqrt(y0))
               do it=1,max_steps
                  sw=(y0-yy*yy)/(2*b)                            ! sinh(w/2)**2
                  x=2*asinh(sqrt(sw))
                  call f_derivs(x,-1,0,f)
                  dt=3*yy*yy*f(0)+a-yy**3*f(1)*2*yy/(b*2*sqrt(sw*(1+sw)))
                  p=(t-yy*(yy*yy*f(0)+a))/dt
                  yy=min(max(yy+p,yy/2),(yy+sqrt(y0))/2)
                  if(abs(p)<=1e-10_dp*yy)exit
               enddo
               if(abs(p)>1e-10_dp*yy)return
               sw=(y0-yy*yy)/(2*b);x=2*asinh(sqrt(sw))
            else
               ! Cubic Hermite for w**2 in (t/tp)**2 between the parabola and
               ! wb, where T**2 ~ a**2*y has slope -a**2*b*sinh(wb).
               x=hermite((t/tp)**2,wb*wb,-2*wb*tp*tp/(a*a*b*sinh(wb)),0._dp,-tp/(2*tv2))
               if(.not.(x>0))x=(tp-t)/tv2
               x=min(sqrt(x),wb*(1-1e-6_dp))
            endif
         else
            ! Long way: T ~ sqrt(-b)*(y0+b)*exp(-w/2) as w -> infinity.
            psi=log(tp/t)
            x=sqrt(tp*psi/tv2)
            if(psi>1)x=max(x,2*log(sqrt(-b)*(y0+b)/t))
         endif
      endif
   else
      ! Multiple revolutions: minimum time by Halley on T'(v) = 0, started
      ! at the large-m minimum of y**1.5/sin(v)**3, cos v = b/(y0+b+sqrt(y0*yend)).
      x=acos(b/(y0+b+sqrt(y0*yend)))
      do it=1,20
         call time_derivs(x,1,m,y0,yend,a,b,d)
         p=-2*d(1)*d(2)/(2*d(2)**2-d(1)*d(3))
         p=min(max(x+p,x/2),(x+pi)/2)
         if(abs(p-x)<=1e-12_dp)exit
         x=p
      enddo
      x=p;call time_derivs(x,1,m,y0,yend,a,b,d)
      tmin=d(0)
      status=5
      if(t<tmin.or..not.(d(2)>0))return
      status=-4
      ! Each branch: cubic Hermite in sig=sqrt(1-(tmin/t)**(1/3)) between the
      ! minimum (v*-v ~ sqrt(6*tmin/T'')*sig) and the pole T ~ K/v**3 at sig=1.
      sig=sqrt(max(1-(tmin/t)**(1._dp/3),0._dp))
      k0=y0*sqrt(y0)*m*pi/sqrt2;kpi=yend*sqrt(yend)*(m+1)*pi/sqrt2
      if(br<0)then
         hi=x
         x=hermite(sig,x,-sqrt(6*tmin/d(2)),0._dp,-2*(k0/tmin)**(1._dp/3))
      else
         lo=x
         x=pi-hermite(sig,pi-x,-sqrt(6*tmin/d(2)),0._dp,-2*(kpi/tmin)**(1._dp/3))
      endif
      x=min(max(x,lo+(hi-lo)*1e-3_dp),hi-(hi-lo)*1e-3_dp)
   endif

   !======= Confluent contour steps =======
   if(x>0.and.yy<0)then
      if(.not.converge())return
   endif

   !======= Velocities (radial/transverse, cancellation-free) =======
   p=sign(1.0_dp,a)
   if(yy>=0)then
      y=yy*yy;dc=-2*sw-ss*ss/(1+cs)             ! cs-cosh(w), from y
   else
      if(k>0)then
         c=cos(x);s=sin(x)**2
         if(b>=0)then
            y=y0+2*b*sin(x/2)**2
         else
            y=yend-2*b*cos(x/2)**2
         endif
      else
         c=cosh(x);s=-sinh(x)**2;y=y0-2*b*sinh(x/2)**2
      endif
      dc=p*cs-c
      if(p*cs*c>0)dc=(s-ss*ss)/(p*cs+c)
   endif
   ratio=sqrt(rn2/rn1)
   vr1=sqrt2*(dc+p*cs*(rn2-rn1)/(rn1*(ratio+1)))/sqrt(y)
   vr2=-sqrt2*(dc+p*cs*(rn1-rn2)/(rn2*(1/ratio+1)))/sqrt(y)
   vt1=sqrt2*ratio*ss/sqrt(y);vt2=sqrt2*ss/(ratio*sqrt(y))
   v1=(vr1*e1+vt1*t1v)/sqrt(scale)
   v2=(vr2*e2+vt2*t2v)/sqrt(scale)
   status=0
   if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4

   contains

      !> Householder order-3 steps on T(x)=t, i.e. the r->0 limit of the
      !! 4-node trapezoidal ratio for the pole of 1/(T-t). A step against
      !! Newton's direction falls back to Newton; steps leaving (lo,hi) are
      !! replaced by bisection toward the violated end.
      logical function converge()
      real(dp) :: g0,c1,c2,c3,b1,b2,b3,dx,xn
      integer  :: j
      converge=.false.
      do j=1,max_steps
         call time_derivs(x,k,m,y0,yend,a,b,d)
         g0=d(0)-t;c1=d(1);c2=d(2)/2;c3=d(3)/6
         if(abs(g0)<=16*epsilon(t)*t)then
            converge=.true.;return          ! residual at rounding level
         endif
         b1=-c1;b2=-(c2*g0+c1*b1);b3=-(c3*g0*g0+c2*g0*b1+c1*b2)
         dx=g0*b2/b3
         if(dx*g0*c1>0)dx=-g0/c1
         if(.not.ieee_is_finite(dx))return
         xn=x+dx
         if(xn<=lo)then
            x=(x+lo)/2
         elseif(xn>=hi)then
            x=(x+hi)/2
         else
            x=xn
            if(abs(dx)<=step_tol*x.and.abs(g0)<=1e-3_dp*t)then
               converge=.true.;return
            endif
         endif
      enddo
      end function converge

   end subroutine contour_simple_one

   !> F and its first three derivatives in x: k=+1, x=v (elliptic);
   !! k=-1, x=w (hyperbolic, sinh/cosh).
   pure subroutine f_derivs(x,k,m,f)
   real(dp),intent(in)  :: x
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: f(0:3)
   real(dp) :: s,c,r,g(0:3),xx
   integer  :: j
   if(m==0.and.x<0.3_dp)then
      ! Parabolic series G(X), X=k*x**2, with Taylor coefficients g(j)=G^(j)/j!.
      xx=k*x*x;g=0;g(0)=pser(1)
      do j=2,10
         g(3)=g(3)*xx+g(2);g(2)=g(2)*xx+g(1);g(1)=g(1)*xx+g(0);g(0)=g(0)*xx+pser(j)
      enddo
      g=g/sqrt2
      f(0)=g(0);f(1)=2*k*x*g(1);f(2)=8*x*x*g(2)+2*k*g(1);f(3)=48*k*x**3*g(3)+24*x*g(2)
   else
      ! Closed form and the ODE  s*F' + 3*c*F = sqrt2  (s'=c, c'=-k*s).
      if(k>0)then
         s=sin(x);c=cos(x)
      else
         s=sinh(x);c=cosh(x)
      endif
      r=1/s
      f(0)=k*(x+m*pi-s*c)/sqrt2*r**3
      f(1)=(sqrt2-3*c*f(0))*r
      f(2)=(3*k*s*f(0)-4*c*f(1))*r
      f(3)=(3*k*c*f(0)+7*k*s*f(1)-5*c*f(2))*r
   endif
   end subroutine f_derivs

   !> T and its first three derivatives in the anomaly coordinate x (k as above).
   pure subroutine time_derivs(x,k,m,y0,yend,a,b,d)
   real(dp),intent(in)  :: x,y0,yend,a,b
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: d(0:3)
   real(dp) :: sh,ch,s,c,y,y1,y2,y3,f(0:3),q,q1,q2,q3,p,p1,p2,p3
   if(k>0)then
      sh=sin(x/2);ch=cos(x/2);c=(ch-sh)*(ch+sh)
      if(b>=0)then
         y=y0+2*b*sh*sh
      else
         y=yend-2*b*ch*ch
      endif
   else
      sh=sinh(x/2);ch=cosh(x/2);c=1+2*sh*sh
      y=y0-2*b*sh*sh
   endif
   s=2*sh*ch
   y1=k*b*s;y2=k*b*c;y3=-b*s
   call f_derivs(x,k,m,f)
   q=sqrt(y);q1=y1/(2*q);q2=(y2/2-q1*q1)/q;q3=(y3/2-3*q1*q2)/q
   p=y*f(0)+a;p1=y1*f(0)+y*f(1);p2=y2*f(0)+2*y1*f(1)+y*f(2)
   p3=y3*f(0)+3*y2*f(1)+3*y1*f(2)+y*f(3)
   d(0)=q*p
   d(1)=q1*p+q*p1
   d(2)=q2*p+2*q1*p1+q*p2
   d(3)=q3*p+3*q2*p1+3*q1*p2+q*p3
   end subroutine time_derivs

   !> Cubic Hermite on [0,1]: values p0,p1 and slopes m0,m1 at 0 and 1.
   pure real(dp) function hermite(u,p0,m0,p1,m1)
   real(dp),intent(in) :: u,p0,m0,p1,m1
   hermite=p0*(2*u**3-3*u*u+1)+m0*(u**3-2*u*u+u)+p1*(3*u*u-2*u**3)+m1*(u**3-u*u)
   end function hermite

   pure function cross(u,v) result(w)
   real(dp),intent(in) :: u(3),v(3)
   real(dp) :: w(3)
   w=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)]
   end function cross

end module contour_simple

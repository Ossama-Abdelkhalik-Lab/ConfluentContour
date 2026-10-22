!> Contour Lambert solver (solver ID 0), zero- and multiple-revolution arcs.
!!
!! REFERENCE
!!   Negrete and Abdelkhalik (2024), DOI 10.2514/1.G008499: the root of the
!!   time-of-flight equation as the pole of 1/(T-t), from a ratio of contour
!!   integrals (Eqs. 9, 14, 24-25).
!!
!! ALGORITHM (details: solvers/confluentcontour/ConfluentContourSolver.md)
!!   1. Normalize the endpoint geometry by its largest radius.
!!   2. Estimate the root: closed-form models refined by compiled-in tables
!!      in the Lancaster-Blanchard parameter (contour_table_*.inc).
!!   3. Extract the root with confluent contour steps: the exact zero-radius
!!      limit of the paper's contour ratio (Householder order 4 on elliptic
!!      arcs, order 3 elsewhere), with derivatives from exact ODEs of the
!!      time function. A converged step counts as verified.
!!   4. Otherwise fall back to sign-bracketed small trapezoidal contours, and
!!      beyond them to interval contours and a certified bisection safeguard
!!      (safeguarded_z).
!!   5. Recover the radial and transverse velocity components.
!!   contour_batch runs the same algorithm in stages over chunks of cases;
!!   its results equal contour_one bit for bit.
!!
!! UNITS AND CONVENTIONS
!!   The routines assume mu = 1 in the input distance/time units. For
!!   dimensional states, nondimensionalize first, or rescale time and
!!   velocity consistently with sqrt(mu). direction = +1 or -1 chooses the
!!   transfer angular momentum orientation; nrev is the number of full
!!   revolutions (default 0); branch = -1 or +1 selects the low-/high-z root
!!   when nrev > 0.
!!
!! STATUS
!!    0: solution found
!!    1: coincident multi-revolution roots at minimum time
!!    5: flight time below the multi-revolution minimum
!!   -1: invalid inputs
!!   -2: nearly collinear geometry (undefined transfer plane)
!!   -4: no sufficiently accurate root obtained
!!   The separate 'half' output is a geometry-conditioning flag:
!!   -1 = degenerate, 0 = ordinary, +1 = near-collinear.
module confluentcontour
   use iso_fortran_env, only: real64, int64
   use ieee_arithmetic
   use iso_c_binding, only: c_double
   implicit none
   private
   integer,parameter,public :: dp=real64
   real(dp),parameter,public :: pi=acos(-1.0_dp),sqrt2=sqrt(2.0_dp)
   real(dp),parameter :: ln2=log(2.0_dp)
   public :: contour_initialize,contour_one,contour_batch
   public :: geometry_data,prepare_geometry,reconstruct,transfer_time,solve_z,transfer_time_upper
   public :: contour_diagnostics
   !> Optional output for inspecting convergence and selected coordinates.
   type contour_diagnostics
      real(dp) :: z                    ! Root in the standard z coordinate
      real(dp) :: q_upper              ! Local upper-angle distance 4*w*(2*pi-w), if used
      real(dp) :: relative_tof_residual ! |T(z)-t|/t (or equivalent upper coordinate)
      logical  :: used_safeguard=.false.
   end type

   !> Precomputed, scale-free endpoint geometry and local velocity directions.
   type geometry_data
      real(dp) :: r1, r2, scale      ! Normalized endpoint radii and length scale
      real(dp) :: a, b, y0, yend    ! Coefficients/end values for y(anomaly)
      real(dp) :: cs, ss            ! |cos(theta/2)| and |sin(theta/2)|
      real(dp) :: e1(3), e2(3)      ! Unit radial vectors
      real(dp) :: t1(3), t2(3)      ! Oriented transverse unit vectors
   end type
   ! Parabolic series of (v-sin(v)cos(v))/sin(v)**3 in v**2 (highest first);
   ! ten terms keep truncation below 3e-19 for |v|<0.3, where the closed form
   ! loses at most eps/0.09 to cancellation.
   real(dp),parameter::pser(10)=[115665628927._dp/2601853536522240000._dp, &
      3821869001._dp/10751460894720000._dp,133541._dp/48117888000._dp,763967._dp/36324288000._dp, &
      1393481._dp/9081072000._dp,1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5,2._dp/3]
   real(dp),parameter::series_limit=.3_dp
   ! min over v of (M*pi+v-sin(v)cos(v))/sin(v)**3, M=1..64 (rounded down);
   ! T = y**1.5*h_M(v)/sqrt2 + a*sqrt(y) exactly, giving a t_min lower bound.
   real(dp),parameter::hmin(64)=[ &
      4.5666632942963430e+00_dp,7.7682145387745551e+00_dp,1.0934624772096207e+01_dp,1.4089860421122442e+01_dp, &
      1.7240094866659767e+01_dp,2.0387655636768560e+01_dp,2.3533618527909585e+01_dp,2.6678549991590746e+01_dp, &
      2.9822776861304746e+01_dp,3.2966501003349848e+01_dp,3.6109853859458525e+01_dp,3.9252924705678339e+01_dp, &
      4.2395776310615162e+01_dp,4.5538454096483910e+01_dp,4.8680991745548120e+01_dp,5.1823414763671224e+01_dp, &
      5.4965742819885243e+01_dp,5.8107991325977864e+01_dp,6.1250172529367759e+01_dp,6.4392296285710259e+01_dp, &
      6.7534370615647134e+01_dp,7.0676402112936970e+01_dp,7.3818396248284444e+01_dp,7.6960357598698877e+01_dp, &
      8.0102290022843604e+01_dp,8.3244196796655018e+01_dp,8.6386080719352165e+01_dp,8.9527944197114792e+01_dp, &
      9.2669789309732323e+01_dp,9.5811617864134561e+01_dp,9.8953431437720923e+01_dp,1.0209523141368673e+02_dp, &
      1.0523701901001964e+02_dp,1.0837879530345108e+02_dp,1.1152056124935804e+02_dp,1.1466231769839203e+02_dp, &
      1.1780406541044671e+02_dp,1.2094580506644772e+02_dp,1.2408753727835138e+02_dp,1.2722926259766147e+02_dp, &
      1.3037098152271446e+02_dp,1.3351269450493578e+02_dp,1.3665440195423324e+02_dp,1.3979610424366268e+02_dp, &
      1.4293780171347859e+02_dp,1.4607949467466176e+02_dp,1.4922118341200115e+02_dp,1.5236286818679480e+02_dp, &
      1.5550454923922274e+02_dp,1.5864622679043845e+02_dp,1.6178790104441626e+02_dp,1.6492957218958773e+02_dp, &
      1.6807124040029447e+02_dp,1.7121290583808087e+02_dp,1.7435456865284735e+02_dp,1.7749622898388071e+02_dp, &
      1.8063788696077734e+02_dp,1.8377954270427136e+02_dp,1.8692119632697973e+02_dp,1.9006284793407283e+02_dp, &
      1.9320449762388017e+02_dp,1.9634614548843768e+02_dp,1.9948779161398315e+02_dp,2.0262943608140600e+02_dp]
   ! Zero-revolution elliptic root tables (estimate refinement, see
   ! table_estimate): 64x64 per half, 64 KB, generated by
   ! tools/gen_table.py and tools/gen_mrev.py.
   include 'contour_table_zero_rev.inc'
   include 'contour_table_multi_rev.inc'
   ! Near-parabolic confluent steps in x=z/4 (kind 0) for -1/4<=x<=0.1225
   ! (-1<=z<=0.49). F uses 13 terms of its exact series in x (truncation
   ! below 3e-19 relative over that range).
   real(dp),parameter :: zstep_xmin=-.25_dp,zstep_xmax=.1225_dp
   real(dp),parameter :: zser(13)=[9203614011768182e-27_dp,7.826388945464640e-11_dp, &
      6.577162419890202e-10_dp,5.451608923746391e-09_dp, &
      115665628927._dp/2601853536522240000._dp,3821869001._dp/10751460894720000._dp, &
      133541._dp/48117888000._dp,763967._dp/36324288000._dp,1393481._dp/9081072000._dp, &
      1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5]
   ! A converged confluent step (agreement of consecutive orders at the
   ! rounding level) counts as verified, so solve_z skips its separate
   ! residual evaluation for it (over 220 000 benchmark cases that check
   ! never rejected a confluent root; worst |T-t|/t = 1.4e-14 vs 1e-10).
   ! Coordinates for fast bracketed contours (see coordinate_inverse):
   !   0: universal z (zero-revolution, default)
   !   1: upper elliptic angle w=pi-v
   !   2: lower elliptic angle v
   !   3: multi-revolution minimum-time critical function K
   !   4: hyperbolic anomaly w=sqrt(-z)/2
   !   5: Y=sqrt(y), for very short hyperbolic transfers
   type coordinate
      integer :: kind=0,m=0
      complex(dp) :: branch=(0._dp,0._dp) ! sqrt(y) branch point in this coordinate
      real(dp) :: inner=huge(1._dp)        ! domain edge facing the other root (M>0)
   end type
   ! Shared upper-semicircle quadrature nodes; filled once by initialize().
   complex(dp),save :: circle(0:128),nodes(0:768)
   real(dp),save :: nodes_r(0:768),nodes_i(0:768)
#ifdef CONTOUR_APPLE_SINCOS
   type,bind(C)::sincos_result
      real(c_double)::s,c
   end type
#endif
   interface
#ifdef CONTOUR_APPLE_SINCOS
      pure type(sincos_result) function sincos_stret(x) bind(C,name='__sincos_stret')
         import sincos_result,c_double
         real(c_double),value,intent(in)::x
      end function
#else
      pure subroutine sincos_c(x,s,c) bind(C,name='sincos')
         import c_double
         real(c_double),value,intent(in)::x
         real(c_double),intent(out)::s,c
      end subroutine
#endif
      pure real(c_double) function expm1(x) bind(C,name='expm1')
         import c_double
         real(c_double),value,intent(in)::x
      end function
      pure real(c_double) function cbrt(x) bind(C,name='cbrt')
         import c_double
         real(c_double),value,intent(in)::x
      end function
   end interface
   logical,save :: initialized=.false.
   ! Contour order targets (r/rho)**(2N-2)~10**(-2*node_digits).
   real(dp),parameter :: node_digits=7
   ! Cases per contour_batch chunk (staged together; the OpenMP work unit).
   integer,parameter :: batch_chunk=256
contains
   !> Lazily initialize the shared semicircle nodes, safely under OpenMP.
   !! The inexpensive atomic-read fast path avoids repeated critical sections.
   subroutine contour_initialize()
      integer j
      logical ready
      ! Fast path without the lock once the tables exist (set last, below).
      !$omp atomic read
      ready=initialized
      if(ready)return
      !$omp critical(contour_node_initialization)
      if(.not.initialized)then
         do j=0,768
            nodes(j)=cmplx(cos(pi*j/768),sin(pi*j/768),dp)
         enddo
         nodes(0)=cmplx(1._dp,0._dp,dp);nodes(768)=cmplx(-1._dp,0._dp,dp)
         nodes_r=real(nodes,dp);nodes_i=aimag(nodes)
         do j=0,128
            circle(j)=nodes(j*6)
         enddo
         !$omp flush
         !$omp atomic write
         initialized=.true.
      endif
      !$omp end critical(contour_node_initialization)
   end subroutine
   !> Cross product for fixed-length Cartesian 3-vectors.
   pure function cross3(a,b) result(c)
      real(dp),intent(in)::a(3),b(3)
      real(dp)::c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function
   !> Validate the endpoints and form dimensionless Lambert geometry.
   !! Invalid/collinear cases return before the geometry fields are used.
   subroutine prepare_geometry(r1,r2,d,g,status,half)
      real(dp),intent(in)::r1(3),r2(3)
      integer,intent(in)::d
      type(geometry_data),intent(out)::g
      integer,intent(out)::status,half
      real(dp)::a,b,h(3),hn,p,q
      status=-1;half=0
      if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)).or.abs(d)/=1) return
      a=sqrt(dot_product(r1,r1));b=sqrt(dot_product(r2,r2))
      if(.not.ieee_is_finite(a+b).or.a<tiny(a)*1e20_dp.or.b<tiny(b)*1e20_dp)then
         a=norm2(r1);b=norm2(r2)
      endif
      if(a<=0.or.b<=0.or..not.ieee_is_finite(a+b))return
      g%scale=max(a,b);g%r1=a/g%scale;g%r2=b/g%scale
      g%e1=r1/a;g%e2=r2/b;h=cross3(g%e1,g%e2);hn=sqrt(dot_product(h,h))
      status=-2
      if(hn<32*epsilon(a))then
         half=-1;return
      endif
      if(hn<1e-4_dp)half=1
      h=h*(d/hn);g%t1=cross3(h,g%e1);g%t2=cross3(h,g%e2)
      ! Norms of sums/differences avoid loss of 1 +/- cos(theta).
      g%cs=sqrt(dot_product(g%e1+g%e2,g%e1+g%e2))/2;g%ss=sqrt(dot_product(g%e2-g%e1,g%e2-g%e1))/2
      p=sqrt(g%r1*g%r2);g%a=d*sqrt2*p*g%cs;g%b=sqrt2*g%a
      q=(g%r1-g%r2)/(sqrt(g%r1)+sqrt(g%r2))
      if(d>0)then
         g%y0=q*q+2*p*g%ss*g%ss/(1+g%cs)
      else
         g%y0=q*q+2*p*(1+g%cs)
      endif
      g%yend=q*q+2*p*g%ss*g%ss/(1+g%cs)
      status=0
   end subroutine
   ! One libm call for sin and cos (GCC otherwise lowers the pair to cexp).
   pure subroutine sincos(x,s,c)
      real(dp),intent(in)::x
      real(dp),intent(out)::s,c
#ifdef CONTOUR_APPLE_SINCOS
      type(sincos_result)::p
      p=sincos_stret(x);s=p%s;c=p%c
#else
      call sincos_c(x,s,c)
#endif
   end subroutine
   ! Parabolic series /sqrt2 for real x=v**2 (negative for hyperbolic v).
   pure function pseries(x) result(f)
      real(dp),intent(in)::x
      real(dp)::f
      integer k
      f=pser(1)
      do k=2,10
         f=f*x+pser(k)
      enddo
      f=f/sqrt2
   end function
   ! Principal complex square root, using separate real and imaginary parts.
   ! Avoids a generic complex-square-root library call on hot paths.
   pure subroutine csqrt_parts(x,y,a,b)
      real(dp),intent(in)::x,y
      real(dp),intent(out)::a,b
      real(dp)::m
      if(y==0)then
         if(x>=0)then
            a=sqrt(x);b=0
         else
            a=0;b=sign(sqrt(-x),y)
         endif
         return
      endif
      m=sqrt(x*x+y*y)
      if(x>=0)then
         a=sqrt((m+x)/2);b=y/(2*a)
      else
         b=sign(sqrt((m-x)/2),y);a=y/(2*b)
      endif
   end subroutine
   ! Shared TOF core in real arithmetic. v is the angle whose half angle feeds
   ! sin/cos (one real sincos and expm1); u is the anomaly numerator; upper
   ! selects the w=pi-v form. Returns T-t as the fraction P/D so that a
   ! contour node needs one complex division.
   pure subroutine frac_core(ur,ui,vr,vi,upper,series_ok,g,t,pr,pim,dr,di)
      real(dp),intent(in)::ur,ui,vr,vi,t
      logical,intent(in)::upper,series_ok
      type(geometry_data),intent(in)::g
      real(dp),intent(out)::pr,pim,dr,di
      real(dp)::sr,cr,em,ep,rp,sy,cy,hr,hi,kr,ki,sn_r,sn_i,s2r,s2i,k2r,k2i,cr2,ci2
      real(dp)::yr,yi,fr,fi,nr,ni,qr,qi,x2r,x2i,ar,ai,t0,cr3,ci3,s3r,s3i
      integer j
      call sincos(vr/2,sr,cr)
      if(vi==0)then
         hr=sr;hi=0;kr=cr;ki=0
      else
         em=expm1(vi/2);ep=1+em;rp=1/ep;sy=.5_dp*em*(1+rp);cy=.5_dp*(ep+rp)
         hr=sr*cy;hi=cr*sy;kr=cr*cy;ki=-sr*sy
      endif
      ! s=sin v=2 sin(v/2)cos(v/2); sin^2(v/2), cos^2(v/2)
      sn_r=2*(hr*kr-hi*ki);sn_i=2*(hr*ki+hi*kr)
      s2r=hr*hr-hi*hi;s2i=2*hr*hi
      k2r=kr*kr-ki*ki;k2i=2*kr*ki
      if(upper)then
         ! v is the upper distance w; cos(pi-w)=-(1-2 sin^2(w/2)).
         if(g%b<0)then
            yr=g%yend-2*g%b*s2r;yi=-2*g%b*s2i
         else
            yr=g%y0+2*g%b*k2r;yi=2*g%b*k2i
         endif
         cr2=-(1-2*s2r);ci2=2*s2i
         x2r=ur*ur-ui*ui;x2i=2*ur*ui
      else
         if(g%b<0)then
            yr=g%yend-2*g%b*k2r;yi=-2*g%b*k2i
         else
            yr=g%y0+2*g%b*s2r;yi=2*g%b*s2i
         endif
         cr2=1-2*s2r;ci2=-2*s2i
         x2r=vr*vr-vi*vi;x2i=2*vr*vi
      endif
      call csqrt_parts(yr,yi,qr,qi)
      if(series_ok.and.x2r*x2r+x2i*x2i<series_limit**4)then
         ! Parabolic series (2/3+x/5+...)/sqrt2 in x=u**2: T-t over D=1.
         fr=pser(1);fi=0
         do j=2,10
            t0=fr*x2r-fi*x2i+pser(j);fi=fr*x2i+fi*x2r;fr=t0
         enddo
         fr=fr/sqrt2;fi=fi/sqrt2
         ar=yr*fr-yi*fi+g%a;ai=yr*fi+yi*fr
         pr=qr*ar-qi*ai-t;pim=qr*ai+qi*ar;dr=1;di=0
      else
         ! T-t = [sqrt(y)*(y*(u-s*c)+sqrt2*a*s**3) - sqrt2*t*s**3]/(sqrt2*s**3)
         nr=ur-(sn_r*cr2-sn_i*ci2);ni=ui-(sn_r*ci2+sn_i*cr2)
         cr3=sn_r*sn_r-sn_i*sn_i;ci3=2*sn_r*sn_i
         s3r=sqrt2*(cr3*sn_r-ci3*sn_i);s3i=sqrt2*(cr3*sn_i+ci3*sn_r)
         ar=yr*nr-yi*ni+g%a*s3r;ai=yr*ni+yi*nr+g%a*s3i
         pr=qr*ar-qi*ai-t*s3r;pim=qr*ai+qi*ar-t*s3i;dr=s3r;di=s3i
      endif
   end subroutine
   ! Critical-point function for M>0: sign(dT/dv)=sign(K) with
   ! K = s*(y'*P+2*y*P') - 6*c*y*P, P = y*(u-s*c)+sqrt2*a*s**3, y'=b*s.
   ! K has no square root, so it is entire in v.
   pure subroutine critical_core(vr,vi,m,g,kr,ki)
      real(dp),intent(in)::vr,vi
      integer,intent(in)::m
      type(geometry_data),intent(in)::g
      real(dp),intent(out)::kr,ki
      complex(dp)::sn,cs,y,n,s2,p,dp1,yd,v
      real(dp)::sr,cr,em,ep,rp,sy,cy
      complex(dp)::h,k
      call sincos(vr/2,sr,cr)
      if(vi==0)then
         h=cmplx(sr,0._dp,dp);k=cmplx(cr,0._dp,dp)
      else
         em=expm1(vi/2);ep=1+em;rp=1/ep;sy=.5_dp*em*(1+rp);cy=.5_dp*(ep+rp)
         h=cmplx(sr*cy,cr*sy,dp);k=cmplx(cr*cy,-sr*sy,dp)
      endif
      v=cmplx(vr,vi,dp)
      sn=2*h*k;cs=1-2*h*h
      if(g%b<0)then
         y=g%yend-2*g%b*k*k
      else
         y=g%y0+2*g%b*h*h
      endif
      n=v+m*pi-sn*cs;s2=sn*sn;yd=g%b*sn
      p=y*n+sqrt2*g%a*s2*sn
      dp1=yd*n+2*y*s2+3*sqrt2*g%a*s2*cs
      h=sn*(yd*p+2*y*dp1)-6*cs*y*p
      kr=real(h,dp);ki=aimag(h)
   end subroutine
   ! Short-way hyperbolic time in Y=sqrt(y) (kind 5). With
   ! S=sinh^2(w/2)=(y0-Y^2)/(2b): sinh(w)=2*sqrt(S(1+S)), cosh(w)=1+2S and
   ! w=2*asinh(sqrt(S)), so y is never formed by cancellation. Small S uses
   ! the parabolic series in w^2=4*S*(asinh(sqrt(S))/sqrt(S))^2.
   pure function sqrty_time(yr,yi,g) result(t)
      real(dp),intent(in)::yr,yi
      type(geometry_data),intent(in)::g
      complex(dp)::t,y,sv,a,b,w,w2,f,sh,ch,p
      ! asinh(sqrt(S))/sqrt(S) in S, highest power first.
      real(dp),parameter::as(10)=[-12155._dp/1245184,6435._dp/557056,-143._dp/10240, &
         231._dp/13312,-63._dp/2816,35._dp/1152,-5._dp/112,3._dp/40,-1._dp/6,1._dp]
      integer k
      y=cmplx(yr,yi,dp)
      sv=(g%y0-y*y)/(2*g%b)
      if(abs(sv)<.0225_dp)then
         p=as(1)
         do k=2,10
            p=p*sv+as(k)
         enddo
         w2=4*sv*p*p
         f=pser(1)
         do k=2,10
            f=f*(-w2)+pser(k)
         enddo
         f=f/sqrt2
      else
         a=sqrt(sv);b=sqrt(1+sv);w=2*log(a+b)
         sh=2*a*b;ch=1+2*sv
         f=(sh*ch-w)/(sqrt2*sh*sh*sh)
      endif
      t=y*(y*y*f+g%a)
   end function
   !> Evaluate complex time of flight via the shared P/D representation.
   pure subroutine tof_core(ur,ui,vr,vi,upper,series_ok,g,tr,ti)
      real(dp),intent(in)::ur,ui,vr,vi
      logical,intent(in)::upper,series_ok
      type(geometry_data),intent(in)::g
      real(dp),intent(out)::tr,ti
      real(dp)::pr,pim,dr,di,d2
      call frac_core(ur,ui,vr,vi,upper,series_ok,g,0._dp,pr,pim,dr,di)
      d2=1/(dr*dr+di*di)
      tr=(pr*dr+pim*di)*d2;ti=(pim*dr-pr*di)*d2
   end subroutine
   !> Evaluate (possibly complex) nondimensional flight time T(z) at rev m.
   !! Primarily used in contour integrals and residual validation.
   pure function transfer_time(z,g,m) result(t)
      complex(dp),intent(in)::z
      type(geometry_data),intent(in)::g
      integer,intent(in)::m
      complex(dp)::t
      real(dp)::ur,ui,tr,ti
      call csqrt_parts(real(z,dp),aimag(z),ur,ui)
      ur=ur/2;ui=ui/2
      call tof_core(ur,ui,ur-m*pi,ui,.false.,m==0,g,tr,ti)
      t=cmplx(tr,ti,dp)
   end function
   ! Real-axis transfer time (z<0 uses the hyperbolic form). Returns huge
   ! where the time is not real (y<0), matching the complex-test behavior.
   pure function real_time(z,g,m) result(t)
      real(dp),intent(in)::z
      type(geometry_data),intent(in)::g
      integer,intent(in)::m
      real(dp)::t,u,v,y,s,c,fun,sh,ch,w,e
      if(z>=0)then
         u=sqrt(z)/2;v=u-m*pi
         call sincos(v/2,sh,ch);s=2*sh*ch;c=1-2*sh*sh
         if(g%b<0)then
            y=g%yend-2*g%b*ch*ch
         else
            y=g%y0+2*g%b*sh*sh
         endif
         if(m==0.and.abs(v)<series_limit)then
            fun=pseries(v*v)
         else
            fun=(u-s*c)/(sqrt2*s*s*s)
         endif
      else
         if(m/=0)then
            t=huge(t);return
         endif
         w=sqrt(-z)/2
         e=expm1(w/2);sh=.5_dp*(e+e/(1+e))
         y=g%y0-2*g%b*sh*sh
         if(w<series_limit)then
            fun=pseries(-w*w)
         else
            ! sinh w=2*sh*ch, cosh w=1+2*sh**2 from the half angle.
            ch=.5_dp*(1+e+1/(1+e));s=2*sh*ch;c=1+2*sh*sh
            fun=(s*c-w)/(sqrt2*s*s*s)
         endif
      endif
      if(y<0.or..not.ieee_is_finite(y))then
         t=huge(t);return
      endif
      t=sqrt(y)*(y*fun+g%a)
   end function
   !> Real-axis derivative dT/dz, used to locate the multi-revolution minimum.
   pure function slope(z,g,m) result(dt)
      real(dp),intent(in)::z
      type(geometry_data),intent(in)::g
      integer,intent(in)::m
      real(dp)::dt,u,v,s,c,y,f,df
      u=sqrt(z)/2;v=u-m*pi;s=sin(v);c=cos(v)
      y=g%y0+2*g%b*sin(v/2)**2
      if(g%b<0)y=g%yend-2*g%b*cos(v/2)**2
      f=(u-s*c)/(sqrt2*s**3);df=sqrt2/s-3*f*c/s
      dt=(sqrt(y)*y*df+g%b*s*(1.5_dp*sqrt(y)*f+g%a/(2*sqrt(y))))/(8*u)
   end function
   !> Recover both roots from the first four Cauchy contour moments.
   !! The coordinate is centered and scaled to improve conditioning.
   subroutine pair_moments(g,tof,m,lo,hi,roots,ok)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::tof,lo,hi
      integer,intent(in)::m
      real(dp),intent(out)::roots(2)
      logical,intent(out)::ok
      real(dp)::mom(0:3),c,r,w,qa,qb,qc,disc
      complex(dp)::e,term,p
      integer j,k
      c=(lo+hi)/2;r=(hi-lo)/2;mom=0
      ! Centered/scaled moments have much better conditioning than raw z**k.
      do j=0,128
         e=circle(j);w=2;if(j==0.or.j==128)w=1
         term=w*e/(transfer_time(c+ r*e,g,m)-tof);p=(1.0_dp,0.0_dp)
         do k=0,3
            mom(k)=mom(k)+real(term*p,dp);p=p*e
         enddo
      enddo
      qa=mom(1)**2-mom(0)*mom(2);qb=mom(3)*mom(0)-mom(1)*mom(2);qc=mom(2)**2-mom(1)*mom(3)
      disc=qb*qb-4*qa*qc;ok=ieee_is_finite(disc).and.disc>=0.and.abs(qa)>tiny(qa)
      roots=0
      if(ok)then
         roots=c+r*[(-qb-sqrt(disc))/(2*qa),(-qb+sqrt(disc))/(2*qa)]
         if(roots(1)>roots(2))roots=roots(2:1:-1)
         ok=all(ieee_is_finite(roots)).and.roots(1)>lo.and.roots(2)<hi
      endif
   end subroutine
   !> Certified-sign-bracket fallback with contour updates and bisection.
   !! Handles zero-revolution and multi-revolution branches without Newton.
   subroutine safeguarded_z(g,tof,m,branch,z,status)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::tof
      integer,intent(in)::m,branch ! -1 lower z; +1 upper z
      real(dp),intent(out)::z
      integer,intent(out)::status
      real(dp)::lo,hi,mid,tmin,tp,zm,c,r,w,num,den,guess,t,f,inc,prev,roots(2),left,right
      complex(dp)::e,term
      integer i,j,n,pass
      logical ok
      status=-4;z=ieee_value(0.0_dp,ieee_quiet_nan)
      if(tof<=0.or..not.ieee_is_finite(tof).or.m<0.or.abs(branch)/=1)then
         status=-1;return
      endif
      if(m==0)then
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         if(abs(tof-tp)<=8*epsilon(tp)*tp)then
            z=0;status=0;return
         endif
         if(tof<tp)then
            hi=0
            if(g%b>0)then
               lo=-16*asinh(sqrt(g%y0/(2*g%b)))**2*(1-8*epsilon(tp))
            else
               lo=-1
               do i=1,16
                  if(real_time(lo,g,0)<tof)exit
                  lo=lo*2
               enddo
            endif
            if(real_time(lo,g,0)>tof)return
         else
            lo=0;hi=4*pi*pi*(1-1e-12_dp)
            if(real_time(hi,g,0)<tof)return
         endif
         inc=1
      else
         lo=(2*m*pi)**2*(1+1e-12_dp);hi=(2*(m+1)*pi)**2*(1-1e-12_dp)
         left=lo;right=hi
         do i=1,60
            mid=(left+right)/2
            if(slope(mid,g,m)<0)then
               left=mid
            else
               right=mid
            endif
         enddo
         zm=(left+right)/2;tmin=real_time(zm,g,m)
         if(tof<tmin-16*epsilon(tof)*tmin)then
            status=5;return
         endif
         if(abs(tof-tmin)<=16*epsilon(tof)*tmin)then
            z=zm;status=1;return
         endif
         call pair_moments(g,tof,m,lo,hi,roots,ok)
         if(branch<0)then
            hi=zm;inc=-1;guess=roots(1)
         else
            lo=zm;inc=1;guess=roots(2)
         endif
         if(ok.and.guess>lo.and.guess<hi)then
            f=real_time(guess,g,m)-tof
            if(inc*f<0)then
               lo=guess
            else
               hi=guess
            endif
         endif
      endif
      prev=huge(1.0_dp)
      do pass=1,96
         c=(lo+hi)/2;r=(hi-lo)/2
         if(r<=4*epsilon(c)*max(1.0_dp,abs(c)))then
            z=c
            if(abs(real_time(z,g,m)-tof)<=2e-12_dp*tof)status=0
            return
         endif
         guess=c
         do n=16,128,16
            if(n/=16.and.n/=32.and.n/=64.and.n/=128)cycle
            num=0;den=0
            do j=0,n
               e=circle(j*(128/n));w=2;if(j==0.or.j==n)w=1
               term=w*e/(transfer_time(c+r*e,g,m)-tof)
               den=den+real(term,dp);num=num+real(term*e,dp)
            enddo
            if(abs(den)<tiny(den).or..not.ieee_is_finite(den))cycle
            guess=c+r*num/den
            if(.not.ieee_is_finite(guess).or.guess<=lo.or.guess>=hi)then
               guess=c;cycle
            endif
            t=real_time(guess,g,m);f=t-tof
            if(abs(f)<=32*epsilon(tof)*tof.and.abs(guess-prev)<=64*epsilon(guess)*max(1.0_dp,abs(guess)))then
               z=guess;status=0;return
            endif
            prev=guess
            if(n>=32)exit
         enddo
         if(mod(pass,3)==0)guess=c
         f=real_time(guess,g,m)-tof
         if(inc*f<0)then
            lo=guess
         else
            hi=guess
         endif
      enddo
   end subroutine
   !> Flight time as a function of q=4*pi**2-z for upper elliptic arcs.
   ! Direct distance to the upper elliptic boundary. Never recover q by
   ! subtracting a rounded z when reconstructing velocities.
   pure function transfer_time_upper(q,g) result(t)
      complex(dp),intent(in)::q
      type(geometry_data),intent(in)::g
      complex(dp)::t
      real(dp)::ar,ai,dr,di,d2,wr,wi,tr,ti
      ! w=q/(2*(2*pi+sqrt(4*pi^2-q))) avoids cancellation near q=0.
      call csqrt_parts(4*pi*pi-real(q,dp),-aimag(q),ar,ai)
      dr=2*(2*pi+ar);di=2*ai;d2=dr*dr+di*di
      wr=(real(q,dp)*dr+aimag(q)*di)/d2;wi=(aimag(q)*dr-real(q,dp)*di)/d2
      call tof_core(pi-wr,-wi,wr,wi,.true.,.true.,g,tr,ti)
      t=cmplx(tr,ti,dp)
   end function

   ! Real upper-coordinate time, q=4*pi^2-z. Huge where not real.
   pure function real_time_upper(q,g) result(t)
      real(dp),intent(in)::q
      type(geometry_data),intent(in)::g
      real(dp)::t,w,s,y,u,f,sh,ch
      if(q<0.or.q>4*pi*pi)then
         t=real(transfer_time_upper(cmplx(q,0._dp,dp),g),dp);return
      endif
      w=q/(2*(2*pi+sqrt(4*pi*pi-q)))
      call sincos(w/2,sh,ch);s=2*sh*ch
      if(g%b<0)then
         y=g%yend-2*g%b*sh*sh
      else
         y=g%y0+2*g%b*ch*ch
      endif
      u=pi-w
      if(abs(u)<series_limit)then
         f=pseries(u*u)
      else
         f=(u+s*(1-2*sh*sh))/(sqrt2*s*s*s)
      endif
      if(y<0)then
         t=huge(t);return
      endif
      t=sqrt(y)*(y*f+g%a)
   end function

   ! Coordinates for the fast bracketed contour. kind 0: z (M=0); kind 1:
   ! upper angle w=pi-v; kind 2: angle v. For M revolutions u=v+M*pi.

   ! Returns 1/(T-t) at a complex node.
   pure function coordinate_inverse(cd,x,g,t) result(f)
      type(coordinate),intent(in)::cd
      complex(dp),intent(in)::x
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t
      complex(dp)::f
      real(dp)::ur,ui,pr,pim,dr,di,d2
      select case(cd%kind)
       case(1)
         call frac_core((cd%m+1)*pi-real(x,dp),-aimag(x),real(x,dp),aimag(x),.true.,cd%m==0,g,t,pr,pim,dr,di)
       case(2)
         call frac_core(real(x,dp)+cd%m*pi,aimag(x),real(x,dp),aimag(x),.false.,cd%m==0,g,t,pr,pim,dr,di)
       case(3)
         call critical_core(real(x,dp),aimag(x),cd%m,g,pr,pim);dr=1;di=0
       case(4)
         call frac_core(-aimag(x),real(x,dp),-aimag(x),real(x,dp),.false.,.true.,g,t,pr,pim,dr,di)
       case(5)
         f=1/(sqrty_time(real(x,dp),aimag(x),g)-t);return
       case default
         call csqrt_parts(real(x,dp),aimag(x),ur,ui)
         call frac_core(ur/2,ui/2,ur/2,ui/2,.false.,.true.,g,t,pr,pim,dr,di)
      end select
      d2=1/(pr*pr+pim*pim)
      f=cmplx((dr*pr+di*pim)*d2,(di*pr-dr*pim)*d2,dp)
   end function
   !> Real residual kernel in the coordinate selected for this contour.
   pure function coordinate_time(cd,x,g) result(t)
      type(coordinate),intent(in)::cd
      real(dp),intent(in)::x
      type(geometry_data),intent(in)::g
      real(dp)::t,ki
      select case(cd%kind)
       case(1)
         t=real_core((cd%m+1)*pi-x,x,.true.,cd%m==0,g)
       case(2)
         t=real_core(x+cd%m*pi,x,.false.,cd%m==0,g)
       case(3)
         call critical_core(x,0._dp,cd%m,g,t,ki)
       case(4)
         t=huge(t)
         if(x>0)t=real_time(-4*x*x,g,0)
       case(5)
         t=real(sqrty_time(x,0._dp,g),dp)
       case default
         t=real_time(x,g,0)
      end select
   end function
   ! coordinate_time for the angle kinds 1,2 from a known half-angle pair.
   pure function coordinate_time_h(cd,x,g,sh,ch) result(t)
      type(coordinate),intent(in)::cd
      real(dp),intent(in)::x,sh,ch
      type(geometry_data),intent(in)::g
      real(dp)::t
      if(cd%kind==1)then
         t=real_core_h((cd%m+1)*pi-x,x,.true.,cd%m==0,g,sh,ch)
      else
         t=real_core_h(x+cd%m*pi,x,.false.,cd%m==0,g,sh,ch)
      endif
   end function
   ! Real TOF from anomaly numerator u and angle v (upper: v is w=pi-v).
   ! Huge where the time is not real (y<0) or the angle is out of range.
   pure function real_core(u,v,upper,series_ok,g) result(t)
      real(dp),intent(in)::u,v
      logical,intent(in)::upper,series_ok
      type(geometry_data),intent(in)::g
      real(dp)::t,s,c,sh,ch
      call real_core_sc(u,v,upper,series_ok,g,t,s,c,sh,ch)
   end function
   ! real_core that also returns s=sin(v), c=cos(v) and the half-angle
   ! sh=sin(v/2), ch=cos(v/2) of the angle argument.
   pure subroutine real_core_sc(u,v,upper,series_ok,g,t,s,cv,sh,ch)
      real(dp),intent(in)::u,v
      logical,intent(in)::upper,series_ok
      type(geometry_data),intent(in)::g
      real(dp),intent(out)::t,s,cv,sh,ch
      call sincos(v/2,sh,ch);s=2*sh*ch;cv=(ch-sh)*(ch+sh)
      t=real_core_h(u,v,upper,series_ok,g,sh,ch)
   end subroutine
   ! real_core from a known half-angle pair sh=sin(v/2), ch=cos(v/2).
   pure function real_core_h(u,v,upper,series_ok,g,sh,ch) result(t)
      real(dp),intent(in)::u,v,sh,ch
      logical,intent(in)::upper,series_ok
      type(geometry_data),intent(in)::g
      real(dp)::t,s,y,f,c
      if(v<=0.or.v>=pi)then
         t=huge(t);return
      endif
      s=2*sh*ch;c=(ch-sh)*(ch+sh)
      if(upper)then
         c=-c
         if(g%b<0)then
            y=g%yend-2*g%b*sh*sh
         else
            y=g%y0+2*g%b*ch*ch
         endif
      else
         if(g%b<0)then
            y=g%yend-2*g%b*ch*ch
         else
            y=g%y0+2*g%b*sh*sh
         endif
      endif
      if(series_ok.and.abs(u)<series_limit)then
         f=pseries(u*u)
      else
         f=(u-s*c)/(sqrt2*s*s*s)
      endif
      if(y<0)then
         t=huge(t);return
      endif
      t=sqrt(y)*(y*f+g%a)
   end function
   ! Half-angle pair at v+d from the pair at v: rotation by d/2 with the
   ! fdlibm sin/cos kernels (|d/2|<=0.78, no range reduction), else libm.
   ! Its latency is a fraction of a libm sincos on the critical path. The
   ! rotation is only absolutely accurate, so libm is used when v+d moves
   ! within half the distance to 0 or pi (sin or cos of the half angle would
   ! lose relative accuracy, which T~1/sin(v)**3 needs).
   pure subroutine rotate_half(v,sh,ch,d,shn,chn)
      real(dp),intent(in)::v,sh,ch,d
      real(dp),intent(out)::shn,chn
      real(dp),parameter::s1=-1.66666666666666324348e-01_dp,s2=8.33333333332248946124e-03_dp, &
         s3=-1.98412698298579493134e-04_dp,s4=2.75573137070700676789e-06_dp,s5=-2.50507602534068634195e-08_dp, &
         s6=1.58969099521155010221e-10_dp
      real(dp),parameter::c1=4.16666666666666019037e-02_dp,c2=-1.38888888888741095749e-03_dp, &
         c3=2.48015872894767294178e-05_dp,c4=-2.75573143513906633035e-07_dp,c5=2.08757232129817482790e-09_dp, &
         c6=-1.13596475577881948265e-11_dp
      real(dp)::y,z,z2,z4,sn,cm
      y=d/2
      if(abs(y)>.78_dp.or.abs(v+d)<.5_dp*abs(v).or.abs(pi-(v+d))<.5_dp*abs(pi-v))then
         call sincos((v+d)/2,shn,chn);return
      endif
      z=y*y;z2=z*z;z4=z2*z2
      sn=y+z*y*((s1+z*s2)+z2*(s3+z*s4)+z4*(s5+z*s6))
      ! cm=cos(y)-1 keeps the small rotation exact to rounding.
      cm=z*(-.5_dp+z*((c1+z*c2)+z2*(c3+z*c4)+z4*(c5+z*c6)))
      shn=sh+(sh*cm+ch*sn);chn=ch+(ch*cm-sh*sn)
   end subroutine
   !> Distance to the closest modeled singularity, setting contour order.
   ! Nearest-singularity distance used to choose the node count.
   pure real(dp) function bdist(c,br)
      real(dp),intent(in)::c
      complex(dp),intent(in)::br
      bdist=sqrt((c-real(br,dp))**2+aimag(br)**2)
   end function
   pure function coordinate_rho(cd,c) result(rho)
      type(coordinate),intent(in)::cd
      real(dp),intent(in)::c
      real(dp)::rho
      real(dp),parameter::r3=sqrt(3._dp)
      select case(cd%kind)
       case(1)
         rho=min(r3*c,bdist(c,cd%branch))
         if(cd%m==0)rho=min(rho,2*(pi-c))
       case(2)
         rho=min(r3*(pi-c),bdist(c,cd%branch))
         if(cd%m==0)then
            rho=min(rho,2*c)
         else
            rho=min(rho,r3*c)
         endif
       case(3)
         rho=min(1._dp,c,pi-c)
       case(4)
         rho=min(2*c,bdist(c,cd%branch),pi)
       case(5)
         ! sqrt(1+S) branch at Y=sqrt(y0+2b); T is odd in Y.
         rho=min(abs(real(cd%branch,dp)-c),c+real(cd%branch,dp))
       case default
         rho=min(r3*(4*pi*pi-c),bdist(c,cd%branch))
      end select
      if(cd%m>0.and.cd%kind/=3)rho=min(rho,2*abs(cd%inner-c))
   end function

   ! Evaluates 1/(T-t) at the interior semicircle nodes j=1..n-1 of the
   ! contour c+r*exp(i*pi*j/n) for coordinate kinds 0,1,2,4. Library calls run
   ! in a scalar loop; the arithmetic runs branch-free over all nodes so the
   ! compiler can vectorize it. Returns node directions and 2*e/(T-t).
   subroutine contour_nodes(cd,g,t,c,r,n,er,ei,fr,fi,hc)
      type(coordinate),intent(in)::cd
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,c,r
      integer,intent(in)::n
      real(dp),intent(in),optional::hc(2)
      real(dp),intent(out)::er(16),ei(16),fr(16),fi(16)
      real(dp),dimension(16)::ur,ui,vr,vi,sr,cr,sy,cy
      real(dp)::vmax,imax,umax
      real(dp)::ybase,ycoef,csign,syj,cyj,hrj,hij,krj,kij,snrj,snij,s2r,s2i,k2r,k2i
      real(dp)::qrj,qij,yrj,yij,cr2j,ci2j,mmj,saj,sbj,sqr,sqi,nrj,nij,c3r,c3i,s3r,s3i
      real(dp)::x2r,x2i,tr,ti,t0,arj,aij,prj,pij,d2j
      real(dp)::xa,kk,ya,za,z2,z4,ra,sa1,ca1,hza,wa,ke,p2,q2,wo,we,em
      integer(int64)::iq
      real(dp),parameter::rmagic=6755399441055744._dp
      real(dp),parameter::pio2a=1.57079632673412561417e+00_dp,pio2b=6.07710050650619224932e-11_dp
      real(dp),parameter::ln2a=6.93147180369123816490e-01_dp,ln2b=1.90821492927058770002e-10_dp
      real(dp),parameter::ks1=-1.66666666666666324348e-01_dp,ks2=8.33333333332248946124e-03_dp, &
         ks3=-1.98412698298579493134e-04_dp,ks4=2.75573137070700676789e-06_dp,ks5=-2.50507602534068634195e-08_dp, &
         ks6=1.58969099521155010221e-10_dp
      real(dp),parameter::kc1=4.16666666666666019037e-02_dp,kc2=-1.38888888888741095749e-03_dp, &
         kc3=2.48015872894767294178e-05_dp,kc4=-2.75573143513906633035e-07_dp,kc5=2.08757232129817482790e-09_dp, &
         kc6=-1.13596475577881948265e-11_dp
      ! sinh(r)/r and (cosh(r)-1)/r**2 in w=r**2, |r|<=ln2/2 (Taylor, next
      ! terms below 4e-19 relative).
      real(dp),parameter::ho(0:6)=[1._dp,1._dp/6,1._dp/120,1._dp/5040,1._dp/362880,1._dp/39916800,1._dp/6227020800._dp]
      real(dp),parameter::he(0:6)=[1._dp/2,1._dp/24,1._dp/720,1._dp/40320,1._dp/3628800,1._dp/479001600, &
         1._dp/87178291200._dp]
      ! Taylor coefficients of sin, cos-1 (and sinh, cosh-1 with alternating
      ! signs removed) for the center-offset node angles.
      real(dp),parameter::ts(5)=[-1._dp/6,1._dp/120,-1._dp/5040,1._dp/362880,-1._dp/39916800]
      real(dp),parameter::tc(5)=[-1._dp/2,1._dp/24,-1._dp/720,1._dp/40320,-1._dp/3628800]
      real(dp)::hcs,hcc,h,wb
      logical::upper,use_series,use_k2
      integer j,nn,stride
      nn=n-1;stride=768/n;upper=cd%kind==1
      ! Node angles v (half angle feeds sin/cos) and anomaly numerators u.
      select case(cd%kind)
       case(1)
         do j=1,nn
            er(j)=nodes_r(j*stride);ei(j)=nodes_i(j*stride)
            vr(j)=c+r*er(j);vi(j)=r*ei(j)
            ur(j)=(cd%m+1)*pi-vr(j);ui(j)=-vi(j)
         enddo
       case(2,3)
         do j=1,nn
            er(j)=nodes_r(j*stride);ei(j)=nodes_i(j*stride)
            vr(j)=c+r*er(j);vi(j)=r*ei(j)
            ur(j)=vr(j)+cd%m*pi;ui(j)=vi(j)
         enddo
       case(4)
         do j=1,nn
            er(j)=nodes_r(j*stride);ei(j)=nodes_i(j*stride)
            vr(j)=-r*ei(j);vi(j)=c+r*er(j)
            ur(j)=vr(j);ui(j)=vi(j)
         enddo
       case default
         ! u=sqrt(z)/2, branchless principal root (Im z/=0 on interior nodes).
         do j=1,nn
            er(j)=nodes_r(j*stride);ei(j)=nodes_i(j*stride)
            xa=c+r*er(j);ya=r*ei(j)
            za=sqrt(xa*xa+ya*ya);sa1=sqrt((za+abs(xa))/2);ca1=sa1*(ya/(za+abs(xa)))
            ur(j)=merge(sa1,abs(ca1),xa>=0)/2;ui(j)=merge(ca1,sign(sa1,ya),xa>=0)/2
            vr(j)=ur(j);vi(j)=ui(j)
         enddo
      end select
      vmax=0;imax=0;umax=0
      do j=1,nn
         vmax=max(vmax,abs(vr(j)));imax=max(imax,abs(vi(j)));umax=max(umax,ur(j)**2+ui(j)**2)
      enddo
      h=r/2
      if(present(hc).and.cd%kind>=1.and.cd%kind<=3.and.h<=.1_dp.and.r<c)then
         ! Angle kinds: v/2=c/2+a+i*b, |a|,|b|<=h. a=(vr-c)/2 is exact
         ! (Sterbenz, r<c), so the trigonometry is that of the rounded node
         ! angle vr that also enters u, as with direct evaluation.
         ! From the exact center pair (sin, cos of c/2), sin(v/2) and cos(v/2)
         ! need only sin a, cos a-1, sinh b, cosh b-1 of |a|,|b|<=h: short
         ! Taylor polynomials with no range reduction (truncation below
         ! 3e-21 relative for h<=0.02 with degree 8, h<=0.1 with degree 12).
         hcs=hc(1);hcc=hc(2)
         if(h<=.02_dp)then
            do j=1,nn
               xa=(vr(j)-c)/2;ya=vi(j)/2;za=xa*xa;wb=ya*ya
               sa1=xa+xa*za*(ts(1)+za*(ts(2)+za*ts(3)))
               ca1=za*(tc(1)+za*(tc(2)+za*(tc(3)+za*tc(4))))
               sy(j)=ya+ya*wb*(-ts(1)+wb*(ts(2)-wb*ts(3)))
               cy(j)=1+wb*(-tc(1)+wb*(tc(2)+wb*(-tc(3)+wb*tc(4))))
               sr(j)=hcs+(hcs*ca1+hcc*sa1);cr(j)=hcc+(hcc*ca1-hcs*sa1)
            enddo
         else
            do j=1,nn
               xa=(vr(j)-c)/2;ya=vi(j)/2;za=xa*xa;wb=ya*ya;z2=za*za;ra=wb*wb
               sa1=xa+xa*za*((ts(1)+za*ts(2))+z2*((ts(3)+za*ts(4))+z2*ts(5)))
               ca1=za*((tc(1)+za*tc(2))+z2*((tc(3)+za*tc(4))+z2*tc(5)))
               sy(j)=ya+ya*wb*((-ts(1)+wb*ts(2))+ra*((-ts(3)+wb*ts(4))-ra*ts(5)))
               cy(j)=1+wb*((-tc(1)+wb*tc(2))+ra*((-tc(3)+wb*tc(4))-ra*tc(5)))
               sr(j)=hcs+(hcs*ca1+hcc*sa1);cr(j)=hcc+(hcc*ca1-hcs*sa1)
            enddo
         endif
      elseif(cd%kind==4.and.h<=.1_dp.and.r<=c/2)then
         ! Hyperbolic kind: v/2=-b'+i*(c/2+a), a=(vi-c)/2 exact, |a|,|b'|<=h.
         ! sinh, cosh of the imaginary part from the center pair (one expm1)
         ! and the offset series (r<=c/2 keeps sinh free of cancellation);
         ! sin, cos of the small real part directly from the series.
         em=expm1(c/2);hcs=.5_dp*(em+em/(1+em));hcc=.5_dp*(1+em+1/(1+em))
         do j=1,nn
            xa=(vi(j)-c)/2;ya=vr(j)/2;za=xa*xa;wb=ya*ya;z2=za*za;ra=wb*wb
            sa1=xa+xa*za*((-ts(1)+za*ts(2))+z2*((-ts(3)+za*ts(4))-z2*ts(5)))
            ca1=za*((-tc(1)+za*tc(2))+z2*((-tc(3)+za*tc(4))-z2*tc(5)))
            sy(j)=hcs+(hcs*ca1+hcc*sa1);cy(j)=hcc+(hcc*ca1+hcs*sa1)
            sr(j)=ya+ya*wb*((ts(1)+wb*ts(2))+ra*((ts(3)+wb*ts(4))+ra*ts(5)))
            cr(j)=1+wb*((tc(1)+wb*tc(2))+ra*((tc(3)+wb*tc(4))+ra*tc(5)))
         enddo
      elseif(imax<=80.and.vmax<=16)then
         ! In-line branch-free sin/cos (fdlibm kernels, Estrin order, Cody-Waite
         ! reduction by pi/2 with a 33-bit leading part) and sinh/cosh; the
         ! 1.5*2**52 shift rounds to nearest without a library call. sinh/cosh
         ! reduce x=k*ln2+r, |r|<=ln2/2, and take e^(+-r)=1+E+-O from the
         ! odd/even series: no division.
         do j=1,nn
            xa=vr(j)/2
            kk=(xa*(2/pi)+rmagic)-rmagic
            ya=(xa-kk*pio2a)-kk*pio2b
            za=ya*ya;z2=za*za;z4=z2*z2
            sa1=ya+za*ya*((ks1+za*ks2)+z2*(ks3+za*ks4)+z4*(ks5+za*ks6))
            ra=za*((kc1+za*kc2)+z2*(kc3+za*kc4)+z4*(kc5+za*kc6))
            hza=.5_dp*za;wa=1-hza
            ca1=wa+(((1-wa)-hza)+za*ra)
            iq=iand(int(kk,int64),3_int64)
            sr(j)=merge(sa1,merge(ca1,merge(-sa1,-ca1,iq==2),iq==1),iq==0)
            cr(j)=merge(ca1,merge(-sa1,merge(-ca1,sa1,iq==2),iq==1),iq==0)
            ke=(vi(j)/2*(1/ln2)+rmagic)-rmagic
            xa=(vi(j)/2-ke*ln2a)-ke*ln2b
            za=xa*xa;z2=za*za;z4=z2*z2
            wo=xa*((ho(0)+za*ho(1))+z2*(ho(2)+za*ho(3))+z4*((ho(4)+za*ho(5))+z2*ho(6)))
            we=za*((he(0)+za*he(1))+z2*(he(2)+za*he(3))+z4*((he(4)+za*he(5))+z2*he(6)))
            p2=transfer(ishft(int(ke,int64)+1023_int64,52),1._dp)
            q2=transfer(ishft(1023_int64-int(ke,int64),52),1._dp)
            ! k=0 gives p2-q2=0 exactly, so sinh=O with no cancellation.
            sy(j)=.5_dp*((p2-q2)*(1+we)+(p2+q2)*wo)
            cy(j)=.5_dp*((p2+q2)*(1+we)+(p2-q2)*wo)
         enddo
      else
         do j=1,nn
            call sincos(vr(j)/2,sr(j),cr(j))
            em=expm1(vi(j)/2)
            sy(j)=.5_dp*em*(1+1/(1+em));cy(j)=.5_dp*(1+em+1/(1+em))
         enddo
      endif
      ! Per-contour choices hoisted out of the node loop.
      if(upper)then
         use_k2=g%b>=0;csign=-1
      else
         use_k2=g%b<0;csign=1
      endif
      if(use_k2)then
         ybase=merge(g%yend,g%y0,g%b<0);ycoef=-2*g%b
         if(g%b>=0)ycoef=2*g%b
      else
         ybase=merge(g%yend,g%y0,g%b<0);ycoef=2*g%b
         if(g%b<0)ycoef=-2*g%b
      endif
      use_series=cd%m==0.and.umax<series_limit**2
      if(cd%kind==3)then
         ! Critical-point function K (see critical_core); returns 2*e/K.
         do j=1,nn
            syj=sy(j);cyj=cy(j)
            hrj=sr(j)*cyj;hij=cr(j)*syj;krj=cr(j)*cyj;kij=-sr(j)*syj
            snrj=2*(hrj*krj-hij*kij);snij=2*(hrj*kij+hij*krj)
            s2r=hrj*hrj-hij*hij;s2i=2*hrj*hij
            k2r=krj*krj-kij*kij;k2i=2*krj*kij
            qrj=merge(k2r,s2r,g%b<0);qij=merge(k2i,s2i,g%b<0)
            yrj=merge(g%yend,g%y0,g%b<0)+merge(-2*g%b,2*g%b,g%b<0)*qrj
            yij=merge(-2*g%b,2*g%b,g%b<0)*qij
            cr2j=1-2*s2r;ci2j=-2*s2i
            ! n=u-s*c, s2=s*s, P=y*n+sqrt2*a*s^3, P'=b*s*n+2*y*s2+3*sqrt2*a*s2*c
            nrj=ur(j)-(snrj*cr2j-snij*ci2j);nij=ui(j)-(snrj*ci2j+snij*cr2j)
            c3r=snrj*snrj-snij*snij;c3i=2*snrj*snij
            s3r=c3r*snrj-c3i*snij;s3i=c3r*snij+c3i*snrj
            arj=yrj*nrj-yij*nij+sqrt2*g%a*s3r;aij=yrj*nij+yij*nrj+sqrt2*g%a*s3i
            prj=g%b*(snrj*nrj-snij*nij)+2*(yrj*c3r-yij*c3i)+3*sqrt2*g%a*(c3r*cr2j-c3i*ci2j)
            pij=g%b*(snrj*nij+snij*nrj)+2*(yrj*c3i+yij*c3r)+3*sqrt2*g%a*(c3r*ci2j+c3i*cr2j)
            ! K = s*(y'*P+2*y*P') - 6*c*y*P with y'=b*s
            x2r=g%b*(snrj*arj-snij*aij)+2*(yrj*prj-yij*pij)
            x2i=g%b*(snrj*aij+snij*arj)+2*(yrj*pij+yij*prj)
            tr=snrj*x2r-snij*x2i;ti=snrj*x2i+snij*x2r
            t0=yrj*arj-yij*aij;d2j=yrj*aij+yij*arj
            prj=tr-6*(cr2j*t0-ci2j*d2j);pij=ti-6*(cr2j*d2j+ci2j*t0)
            d2j=2/(prj*prj+pij*pij)
            arj=prj*d2j;aij=-pij*d2j
            fr(j)=er(j)*arj-ei(j)*aij;fi(j)=er(j)*aij+ei(j)*arj
         enddo
         return
      endif
      ! Fused, branch-free node arithmetic (registers only). T-t is formed as
      ! [sqrt(y)*(y*n+a*S)-t*S]/S with S=sqrt2*s**3 and n=u-s*c, or S=1 and
      ! n=parabolic series (M=0 near v=0, chosen per contour).
      if(use_series)then
         do j=1,nn
            syj=sy(j);cyj=cy(j)
            hrj=sr(j)*cyj;hij=cr(j)*syj;krj=cr(j)*cyj;kij=-sr(j)*syj
            s2r=hrj*hrj-hij*hij;s2i=2*hrj*hij
            k2r=krj*krj-kij*kij;k2i=2*krj*kij
            qrj=merge(k2r,s2r,use_k2);qij=merge(k2i,s2i,use_k2)
            yrj=ybase+ycoef*qrj;yij=ycoef*qij
            mmj=sqrt(yrj*yrj+yij*yij);saj=sqrt((mmj+abs(yrj))/2);sbj=saj*(yij/(mmj+abs(yrj)))
            sqr=merge(saj,abs(sbj),yrj>=0);sqi=merge(sbj,sign(saj,yij),yrj>=0)
            x2r=ur(j)*ur(j)-ui(j)*ui(j);x2i=2*ur(j)*ui(j)
            tr=pser(1)*x2r+pser(2);ti=pser(1)*x2i
            t0=tr*x2r-ti*x2i+pser(3);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(4);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(5);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(6);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(7);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(8);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(9);ti=tr*x2i+ti*x2r;tr=t0
            t0=tr*x2r-ti*x2i+pser(10);ti=tr*x2i+ti*x2r;tr=t0
            nrj=tr/sqrt2;nij=ti/sqrt2
            arj=yrj*nrj-yij*nij+g%a;aij=yrj*nij+yij*nrj
            prj=sqr*arj-sqi*aij-t;pij=sqr*aij+sqi*arj
            d2j=2/(prj*prj+pij*pij)
            arj=prj*d2j;aij=-pij*d2j
            fr(j)=er(j)*arj-ei(j)*aij;fi(j)=er(j)*aij+ei(j)*arj
         enddo
      else
         do j=1,nn
            syj=sy(j);cyj=cy(j)
            hrj=sr(j)*cyj;hij=cr(j)*syj;krj=cr(j)*cyj;kij=-sr(j)*syj
            snrj=2*(hrj*krj-hij*kij);snij=2*(hrj*kij+hij*krj)
            s2r=hrj*hrj-hij*hij;s2i=2*hrj*hij
            k2r=krj*krj-kij*kij;k2i=2*krj*kij
            qrj=merge(k2r,s2r,use_k2);qij=merge(k2i,s2i,use_k2)
            yrj=ybase+ycoef*qrj;yij=ycoef*qij
            cr2j=csign*(1-2*s2r);ci2j=-csign*2*s2i
            mmj=sqrt(yrj*yrj+yij*yij);saj=sqrt((mmj+abs(yrj))/2);sbj=saj*(yij/(mmj+abs(yrj)))
            sqr=merge(saj,abs(sbj),yrj>=0);sqi=merge(sbj,sign(saj,yij),yrj>=0)
            nrj=ur(j)-(snrj*cr2j-snij*ci2j);nij=ui(j)-(snrj*ci2j+snij*cr2j)
            c3r=snrj*snrj-snij*snij;c3i=2*snrj*snij
            s3r=sqrt2*(c3r*snrj-c3i*snij);s3i=sqrt2*(c3r*snij+c3i*snrj)
            arj=yrj*nrj-yij*nij+g%a*s3r;aij=yrj*nij+yij*nrj+g%a*s3i
            prj=sqr*arj-sqi*aij-t*s3r;pij=sqr*aij+sqi*arj-t*s3i
            d2j=2/(prj*prj+pij*pij)
            arj=(s3r*prj+s3i*pij)*d2j;aij=(s3i*prj-s3r*pij)*d2j
            fr(j)=er(j)*arj-ei(j)*aij;fi(j)=er(j)*aij+ei(j)*arj
         enddo
      endif
   end subroutine

   ! One contour over a verified bracket [c-r,c+r]; fm/fp are the known endpoint
   ! residuals. Returns the full estimate, the N/2 subsample estimate and ok.
   subroutine bracket_pole(g,t,cd,c,r,n,fm,fp,x,xh,ok,residue,hc)
      type(geometry_data),intent(in)::g
      type(coordinate),intent(in)::cd
      real(dp),intent(in)::t,c,r,fm,fp
      real(dp),intent(in),optional::hc(2)
      integer,intent(in)::n
      real(dp),intent(out)::x,xh
      logical,intent(out)::ok
      real(dp),intent(out),optional::residue
      real(dp)::num,den,numh,denh,a,b,er(16),ei(16),fr(16),fi(16)
      complex(dp)::e,term
      integer j
      ! Endpoints: e=+1 at c+r and e=-1 at c-r with half weight.
      num=1/fp+1/fm;den=1/fp-1/fm;numh=num;denh=den
      if(cd%kind/=5.and.n<=16)then
         call contour_nodes(cd,g,t,c,r,n,er,ei,fr,fi,hc)
         do j=1,n-1
            a=er(j)*fr(j)-ei(j)*fi(j);b=fr(j)
            num=num+a;den=den+b
            if(mod(j,2)==0)then
               numh=numh+a;denh=denh+b
            endif
         enddo
      else
         do j=1,n-1
            e=nodes(j*(768/n))
            term=2*e*coordinate_inverse(cd,c+r*e,g,t)
            a=real(term*e,dp);b=real(term,dp)
            num=num+a;den=den+b
            if(mod(j,2)==0)then
               numh=numh+a;denh=denh+b
            endif
         enddo
      endif
      x=c+r*num/den;xh=x
      if(mod(n,2)==0)xh=c+r*numh/denh
      ! The sign bracket contains the root: an estimate a rounding hair
      ! outside is clamped to the endpoint.
      if(abs(x-c)>=r.and.abs(x-c)<=1.01_dp*r)x=c+sign(r,x-c)
      if(abs(xh-c)>=r.and.abs(xh-c)<=1.01_dp*r)xh=c+sign(r,xh-c)
      ! Trapezoid residue of 1/(T-t): r*den/(2n) = 1/T'(root).
      if(present(residue))residue=r*den/(2*n)
      ok=ieee_is_finite(x).and.ieee_is_finite(xh).and.abs(x-c)<=r
   end subroutine

   ! Estimate-driven contour: residual at the estimate sets a bracket whose
   ! far endpoint is confirmed by sign; one contour then locates the root.
   ! inc is the sign of dT/dx; [lo,hi] is the open coordinate domain.
   ! dxdt: optional estimate of |dx/dT|, giving the smallest radius whose
   ! nodes resolve T-t above rounding. Later passes size the bracket from the
   ! contour residue 1/T' and the residual at the new center.
   ! h0: sin/cos of x0/2 for the angle kinds 1,2. Later half-angle pairs are
   ! then rotated from known ones (rotate_half) instead of calling sincos, and
   ! hx returns the pair at the root (NaN when unavailable).
   subroutine fast_root(g,t,cd,x0,d0,inc,lo,hi,x,status,residue,dxdt,fa0,h0,hx,checked)
      type(geometry_data),intent(in)::g
      type(coordinate),intent(in)::cd
      real(dp),intent(in)::t,x0,d0,lo,hi
      integer,intent(in)::inc
      real(dp),intent(out)::x
      real(dp),intent(out),optional::residue,hx(2)
      real(dp),intent(in),optional::dxdt,fa0,h0(2)
      integer,intent(out)::status
      logical,intent(out),optional::checked
      real(dp)::xa,fa,xb,fb,c,r,rho,xh,e,step,edge,lg,res,slope_t,floor,sha,cha,shb,chb,sn,cn
      real(dp)::digits
      integer pass,expand,n
      logical ok,have
      xa=x0;step=d0;slope_t=0;res=0
      digits=node_digits
      if(present(checked))checked=.false.
      if(present(dxdt))slope_t=abs(dxdt)
      status=-4;x=ieee_value(t,ieee_quiet_nan)
      if(present(hx))hx=x
      have=.false.
      if(present(h0))have=cd%kind==1.or.cd%kind==2
      if(have)then
         sha=h0(1);cha=h0(2)
      endif
      do pass=1,4
         if(pass==1.and.present(fa0))then
            fa=fa0
         elseif(have)then
            fa=coordinate_time_h(cd,xa,g,sha,cha)-t
         else
            fa=coordinate_time(cd,xa,g)-t
         endif
         if(.not.ieee_is_finite(fa))return
         if(fa==0)then
            x=xa;status=0
            if(have.and.present(hx))hx=[sha,cha]
            return
         endif
         ! Confluent (r->0) contour steps from the estimate; the contour
         ! below remains the fallback.
         if(have.and.t>0.and.pass==1)then
            call confluent_step(cd,g,t,xa,sha,cha,lo,hi,x,sn,cn,ok)
            if(ok)then
               status=0
               if(present(hx))hx=[sn,cn]
               return
            endif
         endif
         if(cd%kind==4.and.t>0.and.pass==1)then
            call confluent_step_h(g,t,xa,lo,hi,x,ok)
            if(ok)then
               status=0;return
            endif
         endif
         if(pass>1)step=max(step,2*abs(fa)*slope_t)
         floor=max(1e-11_dp*abs(xa),2e3_dp*epsilon(t)*abs(t)*slope_t,tiny(t))
         do expand=1,6
            step=max(step,floor)
            if(-fa*inc>0)then
               edge=hi-xa
            else
               edge=xa-lo
            endif
            if(edge<=0)return
            step=min(step,.9_dp*edge)
            xb=xa-sign(step,fa*inc)
            if(have)then
               call rotate_half(xa,sha,cha,xb-xa,shb,chb)
               fb=coordinate_time_h(cd,xb,g,shb,chb)-t
            else
               fb=coordinate_time(cd,xb,g)-t
            endif
            if(.not.ieee_is_finite(fb))return
            if(fb==0)then
               x=xb;status=0
               if(have.and.present(hx))hx=[shb,chb]
               return
            endif
            if((fa>0).neqv.(fb>0))exit
            ! Same sign: first aim past the secant root through (xa,fa),
            ! (xb,fb) when it points onward; afterwards expand geometrically
            ! (a convex T makes the secant undershoot repeatedly).
            e=fb/(fa-fb)
            if(expand==1.and.e>0.and.e<4)then
               step=step*(e+.25_dp*e+.01_dp)
            else
               step=4*step
            endif
            xa=xb;fa=fb
            if(have)then
               sha=shb;cha=chb
            endif
         enddo
         if((fa>0).eqv.(fb>0))return
         c=(xa+xb)/2;r=abs(xb-xa)/2
         rho=coordinate_rho(cd,c)
         ! Smallest N with (r/rho)**(2N-2) negligible: N>=1+7/log10(rho/r).
         lg=rho/r
         if(lg>=10**(digits/3))then
            n=4
         elseif(lg>=10**(digits/5))then
            n=6
         elseif(lg>=10**(digits/7))then
            n=8
         elseif(lg>=10**(digits/11))then
            n=12
         else
            n=16
         endif
         if(have)then
            ! Half-angle pair at the center for the node trigonometry.
            ! A libm pair: exact, and off the rotation chain's latency.
            call sincos(c/2,sn,cn)
            if(xb>xa)then
               call bracket_pole(g,t,cd,c,r,n,fa,fb,x,xh,ok,res,[sn,cn])
            else
               call bracket_pole(g,t,cd,c,r,n,fb,fa,x,xh,ok,res,[sn,cn])
            endif
         elseif(xb>xa)then
            call bracket_pole(g,t,cd,c,r,n,fa,fb,x,xh,ok,res)
         else
            call bracket_pole(g,t,cd,c,r,n,fb,fa,x,xh,ok,res)
         endif
         if(present(residue))residue=res
         if(.not.ok)return
         ! Full-rule error estimate from the N/2 rule: e*e/r for fast
         ! convergence, but only e*(r/rho)**N when r/rho is not small.
         e=abs(x-xh)
         e=e*max(e/r,min(1._dp,r/rho)**n)
         ! Half-angle pair at x from the nearer bracket end.
         if(have)then
            if(abs(x-xa)<=abs(x-xb))then
               call rotate_half(xa,sha,cha,x-xa,sn,cn)
            else
               call rotate_half(xb,shb,chb,x-xb,sn,cn)
            endif
            sha=sn;cha=cn
         endif
         if(e<=2*epsilon(t)*max(abs(x),1e-300_dp))then
            status=0
            if(have.and.present(hx))hx=[sha,cha]
            return
         endif
         xa=x;step=min(4*e+4*epsilon(t)*abs(x),r)
         if(ieee_is_finite(res).and.res/=0)slope_t=abs(res)
      enddo
   end subroutine

   ! Elliptic confluent step of order 4 (Householder/Koenig from the Taylor
   ! coefficients of 1/(T-t): the r->0 limit of the 5-node trapezoidal
   ! contour ratio). Closed-form v-derivatives from the exact ODE
   ! s*F1 + 3*c*F = sqrt2 (any revolution count):
   !   F2=(3sF-4cF1)/s, F3=(3cF+7sF1-5cF2)/s, y1=b*s, y2=b*c, y3=-b*s.
   ! Iterates; accepts when the order-4 minus order-3 difference D predicts
   ! an error D*(D/|step|)**(1/4) below eps*|x|.
   subroutine confluent_step(cd,g,t,xin,shin,chin,lo,hi,x,shn,chn,ok)
      type(coordinate),intent(in)::cd
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,xin,shin,chin,lo,hi
      real(dp),intent(out)::x,shn,chn
      logical,intent(out)::ok
      real(dp)::s,c,u,y,sd,is,f,f1,f2,f3,y1,y2,y3,q,iq,q1,q2,q3,p,p1,p2,p3
      real(dp)::g0,t1,t2,t3,b1,b2,b3,x3,dlt,dd,x0,sh,ch
      real(dp)::f4,q4,p4,t4,b4,x4
      integer it
      integer,parameter::confluent_iters=4
      ok=.false.;x0=xin;sh=shin;ch=chin;x=x0;shn=sh;chn=ch
      do it=1,confluent_iters
         s=2*sh*ch;c=(ch-sh)*(ch+sh)
         if(cd%kind==1)then
            u=(cd%m+1)*pi-x0;c=-c;sd=-1
            if(g%b<0)then
               y=g%yend-2*g%b*sh*sh
            else
               y=g%y0+2*g%b*ch*ch
            endif
         else
            u=x0+cd%m*pi;sd=1
            if(g%b<0)then
               y=g%yend-2*g%b*ch*ch
            else
               y=g%y0+2*g%b*sh*sh
            endif
         endif
         if(.not.(y>0).or.(cd%m==0.and.abs(u)<series_limit))return
         is=1/s
         f=(u-s*c)*is*is*is/sqrt2
         f1=(sqrt2-3*c*f)*is
         f2=(3*s*f-4*c*f1)*is
         f3=(3*c*f+7*s*f1-5*c*f2)*is
         y1=g%b*s;y2=g%b*c;y3=-y1
         q=sqrt(y);iq=1/q
         q1=.5_dp*y1*iq
         q2=(.5_dp*y2-q1*q1)*iq
         q3=(.5_dp*y3-3*q1*q2)*iq
         p=y*f+g%a;p1=y1*f+y*f1;p2=y2*f+2*y1*f1+y*f2;p3=y3*f+3*y2*f1+3*y1*f2+y*f3
         g0=q*p-t
         t1=sd*(q1*p+q*p1)
         t2=(q2*p+2*q1*p1+q*p2)/2
         t3=sd*(q3*p+3*q2*p1+3*q1*p2+q*p3)/6
         b1=-t1;b2=-(t2*g0+t1*b1);b3=-(t3*g0*g0+t2*g0*b1+t1*b2)
         ! Order 4 (5-node confluent limit): one more derivative from the
         ! same ODE, s F4 = -3sF + 10cF1 + 12sF2 - 6cF3, y4 = -b c.
         f4=(-3*s*f+10*c*f1+12*s*f2-6*c*f3)*is
         q4=(-.5_dp*y2-4*q1*q3-3*q2*q2)*iq
         p4=-y2*f+4*y3*f1+6*y2*f2+4*y1*f3+y*f4
         t4=(q4*p+4*q3*p1+6*q2*p2+4*q1*p3+q*p4)/24
         b4=-(t4*g0*g0*g0+t3*g0*g0*b1+t2*g0*b2+t1*b3)
         x3=x0+g0*b2/b3;x4=x0+g0*b3/b4
         if(.not.(x4>lo.and.x4<hi))return
         call sincos(x4/2,shn,chn)
         dd=abs(x4-x3);dlt=abs(x4-x0)
         ! Order-4 error ~ D*(D/|step|)**(1/4), D = |x4-x3|.
         if(dd<=4*epsilon(x4)*abs(x4).or.(dd<dlt.and.dd**5<=dlt*(epsilon(x4)*abs(x4))**4))then
            x=x4;ok=.true.;return
         endif
         if(.not.(dd<dlt))return
         x0=x4;sh=shn;ch=chn
      enddo
   end subroutine
   ! Zero-revolution elliptic root from the compiled-in tables. The root
   ! depends only on the Lancaster-Blanchard parameter
   ! lam = d*sqrt(r1*r2)*cos(theta/2)/s_L and s = (tp/t)**(1/3). Tables use
   ! tau = lam/(1+sqrt(1-lam**2)) (grid stretched toward |lam| = 1) and
   ! store A = (pi-v)/s for s <= 1/2 and B = z/(1-s) for s > 1/2, which are
   ! smooth through t -> infinity and t -> tp. Local 4x4 cubic Lagrange
   ! interpolation. Leaves zg, qg unchanged outside the tables.
   pure subroutine table_estimate(g,s,zg,qg)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::s
      real(dp),intent(inout)::zg,qg
      real(dp),parameter::bq=4*pi*pi,smax=.99_dp
      real(dp)::tau,u,p,wt(4),ws(4),acc,w,zt,qt
      integer i,j,a,b
      if(.not.(s>0.and.s<=smax))return
      call lb_coord(g,tau,w)
      if(.not.(abs(tau)<=tab_tmax))return
      u=(tau+tab_tmax)*((tab_n-1)/(2*tab_tmax))
      i=min(max(int(u)-1,0),tab_n-4);p=u-i
      call lagrange4(p,wt)
      if(s<=.5_dp)then
         u=s*(2*(tab_n-1))
      else
         u=(s-.5_dp)*(2*(tab_n-1))
      endif
      j=min(max(int(u)-1,0),tab_n-4);p=u-j
      call lagrange4(p,ws)
      acc=0
      if(s<=.5_dp)then
         do b=1,4
            do a=1,4
               acc=acc+wt(a)*ws(b)*tab_a(i+a,j+b)
            enddo
         enddo
         w=acc*s
         qt=4*w*(2*pi-w);zt=bq-qt
      else
         do b=1,4
            do a=1,4
               acc=acc+wt(a)*ws(b)*tab_b(i+a,j+b)
            enddo
         enddo
         zt=acc*(1-s);qt=bq-zt
      endif
      if(zt>0.and.qt>0)then
         zg=zt;qg=qt
      endif
   end subroutine
   pure subroutine lagrange4(p,w)
      real(dp),intent(in)::p
      real(dp),intent(out)::w(4)
      w(1)=-(p-1)*(p-2)*(p-3)/6;w(2)=p*(p-2)*(p-3)/2
      w(3)=-p*(p-1)*(p-3)/2;w(4)=p*(p-1)*(p-2)/6
   end subroutine
   ! Hyperbolic confluent step (coordinate kind 4, x = w = sqrt(-z)/2). With
   ! S = sinh w, C = cosh w, F = (S*C - w)/(sqrt2 S**3) satisfies the same
   ! exact ODE as the elliptic case, S F1 + 3 C F = sqrt2, so
   !   F2 = -(4C F1 + 3S F)/S,  F3 = -(5C F2 + 7S F1 + 3C F)/S,
   ! and y = y0 - 2b sinh(w/2)**2: y1 = -b S, y2 = -b C, y3 = -b S.
   ! Order-3 step; accepts when the order-3 minus order-2 (Halley)
   ! difference D predicts an error D*(D/|step|)**(1/3) below eps*|x|.
   ! sinh/cosh are recomputed from one expm1 per step.
   subroutine confluent_step_h(g,t,xin,lo,hi,x,ok)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,xin,lo,hi
      real(dp),intent(out)::x
      logical,intent(out)::ok
      real(dp)::e,hs,hc,s,c,y,is,f,f1,f2,f3,y1,y2,y3,q,iq,q1,q2,q3,p,p1,p2,p3
      real(dp)::g0,t1,t2,t3,b1,b2,b3,x2,x3,dlt,dd,x0
      integer it
      integer,parameter::confluent_iters=4
      ok=.false.;x0=xin;x=x0
      do it=1,confluent_iters
         if(.not.(x0>=series_limit))return
         e=expm1(x0/2);hs=.5_dp*(e+e/(1+e));hc=.5_dp*(1+e+1/(1+e))
         s=2*hs*hc;c=1+2*hs*hs
         y=g%y0-2*g%b*hs*hs
         if(.not.(y>0))return
         is=1/s
         f=(s*c-x0)*is*is*is/sqrt2
         f1=(sqrt2-3*c*f)*is
         f2=-(4*c*f1+3*s*f)*is
         f3=-(5*c*f2+7*s*f1+3*c*f)*is
         y1=-g%b*s;y2=-g%b*c;y3=y1
         q=sqrt(y);iq=1/q
         q1=.5_dp*y1*iq
         q2=(.5_dp*y2-q1*q1)*iq
         q3=(.5_dp*y3-3*q1*q2)*iq
         p=y*f+g%a;p1=y1*f+y*f1;p2=y2*f+2*y1*f1+y*f2;p3=y3*f+3*y2*f1+3*y1*f2+y*f3
         g0=q*p-t
         t1=q1*p+q*p1
         t2=(q2*p+2*q1*p1+q*p2)/2
         t3=(q3*p+3*q2*p1+3*q1*p2+q*p3)/6
         b1=-t1;b2=-(t2*g0+t1*b1);b3=-(t3*g0*g0+t2*g0*b1+t1*b2)
         x2=x0+g0*b1/b2;x3=x0+g0*b2/b3
         if(.not.(x3>lo.and.x3<hi))return
         dd=abs(x3-x2);dlt=abs(x3-x0)
         if(dd<=4*epsilon(x3)*abs(x3).or.(dd<dlt.and.dd*dd*dd*dd<=dlt*(epsilon(x3)*abs(x3))**3))then
            x=x3;ok=.true.;return
         endif
         if(.not.(dd<dlt))return
         x0=x3
      enddo
   end subroutine
   ! Lancaster-Blanchard parameter in the stretched table coordinate:
   ! lam = d*sqrt(r1*r2)*cos(theta/2)/s_L, tau = lam/(1+sqrt(1-lam**2)),
   ! with 1-lam**2 = c/s_L (chord over semiperimeter).
   ! Also returns the semiperimeter s_L (times scale with s_L**1.5).
   pure subroutine lb_coord(g,tau,sl)
      type(geometry_data),intent(in)::g
      real(dp),intent(out)::tau,sl
      real(dp)::c
      c=sqrt((g%r1-g%r2)**2+4*g%r1*g%r2*g%ss*g%ss);sl=(g%r1+g%r2+c)/2
      tau=(g%a/(sqrt2*sl))/(1+sqrt(c/sl))
   end subroutine
   ! Multi-revolution table estimate (see contour_table_multi_rev.inc).
   ! stat = 5: t is below the tabulated minimum time by a safe margin
   ! (infeasible); stat = 1: estimate v of the requested branch and the
   ! minimum-time location vstar; stat = 0: not available (outside the
   ! tables or too close to t_min), use the model path.
   ! The root depends on lam, sigma = branch*sqrt(1-(t_min/t)**(1/3)) and
   ! mu = 1/m; slices in mu are interpolated by 4-point Lagrange.
   pure subroutine mrev_estimate(g,t,m,branch,v,vstar,stat)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t
      integer,intent(in)::m,branch
      real(dp),intent(out)::v,vstar
      integer,intent(out)::stat
      real(dp),parameter::margin=1e-4_dp
      real(dp)::tau,u,p,wt(4),wc(4),ws(4),wm(4),mu,cm,r,sig,acc,sl
      integer i,ic,j,k,km(4),nk,a,b,kk
      stat=0;v=0;vstar=0
      call lb_coord(g,tau,sl)
      if(.not.(abs(tau)<=mr_tmax))return
      ! slices in mu = 1/m (mr_mu decreasing)
      mu=1._dp/m;nk=0
      do k=1,mr_nm
         if(mr_m(k)==m)then
            nk=1;km(1)=k;wm(1)=1
         endif
      enddo
      if(nk==0)then
         kk=mr_nm-1
         do k=1,mr_nm-1
            if(mu>=mr_mu(k+1))then
               kk=k;exit
            endif
         enddo
         kk=min(max(kk-1,1),mr_nm-3);nk=4
         do a=1,4
            km(a)=kk+a-1;wm(a)=1
            do b=1,4
               if(b/=a)wm(a)=wm(a)*(mu-mr_mu(kk+b-1))/(mr_mu(kk+a-1)-mr_mu(kk+b-1))
            enddo
         enddo
      endif
      ! fine tau grid: cbrt(t_min/m) and v*
      u=(tau+mr_tmax)*((mr_nc-1)/(2*mr_tmax))
      ic=min(max(int(u)-1,0),mr_nc-4);p=u-ic
      call lagrange4(p,wc)
      cm=0;vstar=0
      do k=1,nk
         do a=1,4
            cm=cm+wm(k)*wc(a)*mr_c(ic+a,km(k))
            vstar=vstar+wm(k)*wc(a)*mr_s(ic+a,km(k))
         enddo
      enddo
      ! r = (t_min/t)**(1/3); the table holds cbrt(t_min/m)/sqrt(s_L).
      r=cm*sqrt(sl)*cbrt(m/t)
      if(r>1+margin)then
         stat=5;return
      endif
      if(r>=1-margin)return
      sig=branch*sqrt(1-r)
      u=(tau+mr_tmax)*((mr_nt-1)/(2*mr_tmax))
      i=min(max(int(u)-1,0),mr_nt-4);p=u-i
      call lagrange4(p,wt)
      u=(sig+1)*((mr_ns-1)/2._dp)
      j=min(max(int(u)-1,0),mr_ns-4);p=u-j
      call lagrange4(p,ws)
      acc=0
      do k=1,nk
         do b=1,4
            do a=1,4
               acc=acc+wm(k)*wt(a)*ws(b)*mr_v(i+a,j+b,km(k))
            enddo
         enddo
      enddo
      v=acc
      if(v>0.and.v<pi.and.vstar>0.and.vstar<pi.and.(v-vstar)*branch>0)stat=1
   end subroutine
   ! Near-parabolic confluent step in x = z/4 = v**2 (coordinate kind 0, both
   ! signs of z). Here T is a polynomial and one square root: F is the
   ! parabolic series in x and 1-cos(v) = sum (-1)**(k+1) x**k/(2k)!, so T,
   ! T', T'', T''' follow from one Horner pass with derivatives. Valid for
   ! |x| <= zstep_xmax, where the ten-term series is below 1e-17 relative.
   ! Order-3 step and acceptance test as in confluent_step_h.
   subroutine confluent_step_z(g,t,zin,zlo,zhi,z,ok)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,zin,zlo,zhi
      real(dp),intent(out)::z
      logical,intent(out)::ok
      ! 1-cos(v) in x=v**2: coefficients of x**1..x**8 (highest first)
      real(dp),parameter::cser(8)=[-1._dp/20922789888000._dp,1._dp/87178291200._dp,-1._dp/479001600._dp, &
         1._dp/3628800._dp,-1._dp/40320._dp,1._dp/720._dp,-1._dp/24._dp,.5_dp]
      real(dp)::x0,f,f1,f2,f3,c,c1,c2,c3,y,y1,y2,y3,q,iq,q1,q2,q3,p,p1,p2,p3
      real(dp)::g0,t1,t2,t3,b1,b2,b3,x2,x3,dlt,dd,lo,hi
      integer it,j
      integer,parameter::confluent_iters=4
      ok=.false.;z=zin;x0=zin/4;lo=zlo/4;hi=zhi/4
      do it=1,confluent_iters
         if(.not.(x0>=zstep_xmin.and.x0<=zstep_xmax))then
            z=4*x0;return
         endif
         f=zser(1);f1=0;f2=0;f3=0
         do j=2,13
            f3=f3*x0+f2;f2=f2*x0+f1;f1=f1*x0+f;f=f*x0+zser(j)
         enddo
         f3=f3*x0+f2;f2=f2*x0+f1;f1=f1*x0+f;f=f*x0+2._dp/3
         f=f/sqrt2;f1=f1/sqrt2;f2=2*f2/sqrt2;f3=6*f3/sqrt2
         c=cser(1);c1=0;c2=0;c3=0
         do j=2,8
            c3=c3*x0+c2;c2=c2*x0+c1;c1=c1*x0+c;c=c*x0+cser(j)
         enddo
         ! c(x) above lacks the overall factor x: 1-cos v = x*c(x)
         c3=6*c3;c2=2*c2
         y=g%y0+g%b*(x0*c);y1=g%b*(c+x0*c1);y2=g%b*(2*c1+x0*c2);y3=g%b*(3*c2+x0*c3)
         if(.not.(y>0))return
         q=sqrt(y);iq=1/q
         q1=.5_dp*y1*iq
         q2=(.5_dp*y2-q1*q1)*iq
         q3=(.5_dp*y3-3*q1*q2)*iq
         p=y*f+g%a;p1=y1*f+y*f1;p2=y2*f+2*y1*f1+y*f2;p3=y3*f+3*y2*f1+3*y1*f2+y*f3
         g0=q*p-t
         t1=q1*p+q*p1
         t2=(q2*p+2*q1*p1+q*p2)/2
         t3=(q3*p+3*q2*p1+3*q1*p2+q*p3)/6
         b1=-t1;b2=-(t2*g0+t1*b1);b3=-(t3*g0*g0+t2*g0*b1+t1*b2)
         x2=x0+g0*b1/b2;x3=x0+g0*b2/b3
         if(.not.(x3>lo.and.x3<hi))return
         dd=abs(x3-x2);dlt=abs(x3-x0)
         if(dd<=4*epsilon(x3)*abs(x3).or.(dd<dlt.and.dd*dd*dd*dd<=dlt*(epsilon(x3)*abs(x3))**3))then
            z=4*x3;ok=.true.;return
         endif
         z=4*x3
         if(.not.(dd<dlt))return
         x0=x3
      enddo
   end subroutine
   ! Very short hyperbolic arcs: confluent step in Y = sqrt(y) (kind 5).
   ! With S = sinh(w/2)**2 = (y0 - Y**2)/(2b), G(S) = F(w(S)) satisfies
   !   2S(1+S) G' + 3(1+2S) G = sqrt2,
   ! so G'' = -(5(1+2S)G' + 6G)/(2S(1+S)), G''' = -(7(1+2S)G'' + 16G')/(2S(1+S)),
   ! and T = Y**3 G(S(Y)) + a Y with S' = -Y/b, S'' = -1/b. Declines where the
   ! closed form for G cancels (S < 0.0225, the solver's series region).
   subroutine confluent_step_y(g,t,xin,lo,hi,x,ok)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,xin,lo,hi
      real(dp),intent(out)::x
      logical,intent(out)::ok
      ! G(S) = (sqrt2/3) 2F1(1,3;5/2;-S): coefficients of S**11..S**0, used for
      ! S < 0.0225 where the closed form cancels (truncation below 1e-19).
      real(dp),parameter::gser(12)=[ &
         -1.14063450214713069e+00_dp, &
         1.09676394437224101e+00_dp, &
         -1.05106544669006419e+00_dp, &
         1.00328974456778841e+00_dp, &
         -9.53125257339399012e-01_dp, &
         9.00173854153876851e-01_dp, &
         -8.43912988269259534e-01_dp,7.83633489107169456e-01_dp,-7.18330698348238705e-01_dp,6.46497628513414857e-01_dp,-5.65685424949238014e-01_dp,4.71404520791031734e-01_dp]
      real(dp)::yy,s,r,sh,c,w,gg,g1,g2,g3,id,s1,s2,h1,h2,h3,y2
      integer j
      real(dp)::g0,t1,t2,t3,b1,b2,b3,x2,x3,dlt,dd,x0
      integer it
      integer,parameter::confluent_iters=4
      ok=.false.;x0=xin;x=x0
      do it=1,confluent_iters
         yy=x0;s=(g%y0-yy*yy)/(2*g%b)
         if(.not.(s>=0))return
         if(s<.0225_dp)then
            gg=gser(1);g1=0;g2=0;g3=0
            do j=2,12
               g3=g3*s+g2;g2=g2*s+g1;g1=g1*s+gg;gg=gg*s+gser(j)
            enddo
            g2=2*g2;g3=6*g3
         else
            r=sqrt(s);c=sqrt(1+s);sh=2*r*c;w=2*log(r+c);c=1+2*s
            gg=(sh*c-w)/(sqrt2*sh*sh*sh)
            id=1/(2*s*(1+s))
            g1=(sqrt2-3*c*gg)*id
            g2=-(5*c*g1+6*gg)*id
            g3=-(7*c*g2+16*g1)*id
         endif
         s1=-yy/g%b;s2=-1/g%b
         h1=g1*s1;h2=g2*s1*s1+g1*s2;h3=(g3*s1*s1+3*g2*s2)*s1
         y2=yy*yy
         g0=yy*(y2*gg+g%a)-t
         t1=3*y2*gg+y2*yy*h1+g%a
         t2=(6*yy*gg+6*y2*h1+y2*yy*h2)/2
         t3=(6*gg+18*yy*h1+9*y2*h2+y2*yy*h3)/6
         b1=-t1;b2=-(t2*g0+t1*b1);b3=-(t3*g0*g0+t2*g0*b1+t1*b2)
         x2=x0+g0*b1/b2;x3=x0+g0*b2/b3
         if(.not.(x3>lo.and.x3<hi))return
         dd=abs(x3-x2);dlt=abs(x3-x0)
         if(dd<=4*epsilon(x3)*abs(x3).or.(dd<dlt.and.dd*dd*dd*dd<=dlt*(epsilon(x3)*abs(x3))**3))then
            x=x3;ok=.true.;return
         endif
         if(.not.(dd<dlt))return
         x0=x3
      enddo
   end subroutine
   ! Zero-revolution elliptic estimate: inverse cubic Hermite of the upper
   ! distance q(phi), phi=T**(-1/3), matching q=0 (T~K/q**3) and q=4*pi^2
   ! (parabolic time and slope). Returns the estimate in z and in q and a
   ! local slope |dx/dphi|.
   subroutine elliptic_estimate(g,t,tp,der,zg,qg,slope_phi,ct,sout)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,tp,der
      real(dp),intent(out)::zg,qg,slope_phi,ct
      real(dp),intent(out),optional::sout
      real(dp),parameter::bq=4*pi*pi,ktop=8*pi**(4._dp/3)/2**(1._dp/6)
      real(dp)::cp,s,eps1,m0,m1
      ct=cbrt(t);cp=cbrt(tp)
      s=cp/ct;eps1=(t-tp)/(ct*ct+ct*cp+cp*cp)/ct
      if(present(sout))sout=s
      m0=ktop*sqrt(g%y0+2*g%b)/cp;m1=3*tp/der
      if(s>.5_dp)then
         zg=bq*eps1*eps1*(3-2*eps1)-s*eps1*eps1*m0+s*s*eps1*m1
         qg=bq-zg
      else
         qg=s*(1-s)**2*m0+s*s*(3-2*s)*bq+s*s*(s-1)*m1
         ! A steep parabolic slope can overshoot; keep the top asymptote.
         if(qg<.25_dp*s*m0)qg=s*m0
         zg=bq-qg
      endif
      slope_phi=abs(((3*s-4)*s+1)*m0+6*s*(1-s)*bq+(3*s-2)*s*m1)*cp
   end subroutine

   subroutine elliptic_fast(g,t,tp,der,x,upper,status,hx,checked)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,tp,der
      real(dp),intent(out)::x,hx(2)
      logical,intent(out)::upper,checked
      integer,intent(out)::status
      real(dp),parameter::bq=4*pi*pi
      real(dp)::zg,qg,sp,d0,f,gam,ct,sh0,ch0,sr,w,sn,cn
      complex(dp)::zb,sb
      type(coordinate)::cd
      logical::ok
      status=-4;x=ieee_value(t,ieee_quiet_nan);upper=.false.;hx=x;checked=.false.
      if(.not.(der>0))return
      ! Fast path: table estimate and confluent steps straight away. The
      ! closed-form estimate, branch point and bracket/contour set-up below
      ! run only when the table does not apply or the step declines.
      zg=ieee_value(t,ieee_quiet_nan);qg=zg
      call table_estimate(g,cbrt(tp/t),zg,qg)
      if(zg>.49_dp)then
         upper=zg>bq/2
         if(upper)then
            cd%kind=1
            w=min(max(qg,tiny(qg)),bq/2);w=w/(2*(2*pi+sqrt(bq-w)))
         else
            cd%kind=2
            w=sqrt(zg)/2
         endif
         call sincos(w/2,sh0,ch0)
         call confluent_step(cd,g,t,w,sh0,ch0,0._dp,pi,x,sn,cn,ok)
         if(ok)then
            status=0;hx=[sn,cn];checked=.true.
            if(upper)then
               x=4*x*(2*pi-x)
            else
               x=4*x*x
            endif
            return
         endif
         x=ieee_value(t,ieee_quiet_nan)
      endif
      call elliptic_estimate(g,t,tp,der,zg,qg,sp,ct,sr)
      call table_estimate(g,sr,zg,qg)
      if(.not.ieee_is_finite(zg).or..not.ieee_is_finite(qg))return
      upper=zg>bq/2
      ! Steps below use f**(-1/3)-t**(-1/3)=(t-f)/(cbrt(t)*(t+2f)), the
      ! (1,1) Pade form of u**(-1/3)-1 (within 3% for 1/2<f/t<2): no cbrt
      ! on the critical path.
      ! Branch point of sqrt(y) in z (only sizes the node count; log forms
      ! of asinh/acosh suffice).
      if(g%b>0)then
         gam=sqrt(g%y0/(2*g%b))
         gam=log(gam+sqrt(1+gam*gam))
         zb=cmplx(-16*gam*gam,0._dp,dp);sb=cmplx(0._dp,4*gam,dp)
      else
         gam=max(1._dp,sqrt(g%y0/(-2*g%b)))
         gam=2*log(gam+sqrt(gam*gam-1))
         sb=2*cmplx(pi,gam,dp);zb=sb*sb
      endif
      if(upper)then
         ! Upper angle w=pi-v from q without cancellation; dw/dq=1/(8*pi-4*w).
         cd%kind=1;cd%branch=pi-sb/2
         qg=min(max(qg,tiny(qg)),bq/2)
         zg=qg/(2*(2*pi+sqrt(bq-qg)))
         call sincos(zg/2,sh0,ch0);f=coordinate_time_h(cd,zg,g,sh0,ch0)
         d0=abs((t-f)/(ct*(t+2*f)))*sp/(8*pi-4*zg)
         call fast_root(g,t,cd,zg,1.25_dp*d0,-1,0._dp,pi,x,status, &
            dxdt=sp/(8*pi-4*zg)/(3*t*ct),fa0=f-t,h0=[sh0,ch0],hx=hx,checked=checked)
         if(status==0)x=4*x*(2*pi-x)
      elseif(zg>.49_dp)then
         ! Angle v=sqrt(z)/2 avoids a complex sqrt per node; the even image
         ! root at -v is kept out of the contour by the node-count rule.
         cd%kind=2;cd%branch=sb/2
         zg=sqrt(min(zg,bq/2))/2
         call sincos(zg/2,sh0,ch0);f=coordinate_time_h(cd,zg,g,sh0,ch0)
         d0=abs((t-f)/(ct*(t+2*f)))*sp/(8*zg)
         call fast_root(g,t,cd,zg,1.25_dp*d0,1,0._dp,pi,x,status,dxdt=sp/(8*zg)/(3*t*ct),fa0=f-t, &
            h0=[sh0,ch0],hx=hx,checked=checked)
         if(status==0)x=4*x*x
      else
         cd%kind=0;cd%branch=zb
         zg=max(zg,tiny(zg))
         call confluent_step_z(g,t,zg,0._dp,bq,x,ok)
         if(ok)then
            status=0;checked=.true.;return
         endif
         f=coordinate_time(cd,zg,g)
         d0=abs((t-f)/(ct*(t+2*f)))*sp
         call fast_root(g,t,cd,zg,1.25_dp*d0,1,0._dp,bq,x,status,dxdt=sp/(3*t*ct),fa0=f-t,checked=checked)
      endif
   end subroutine


   ! Multi-revolution branch estimate and bracketed contour. With
   ! H=sqrt(y)*(y*(u-s*c)+sqrt2*a*s**3), T=H/(sqrt2*s**3), so
   ! phi=T**(-1/3)=sin(v)*g(v) with smooth g=(sqrt2/H)**(1/3). log g is
   ! modeled linearly from the exact end values and one exact value at the
   ! model minimum vm; T(vm)<t proves one root on each side of vm.
   ! Returns v (lower branch) or w=pi-v (upper branch).
   subroutine multirev_fast(g,t,m,branch,x,upper,status,hx,checked)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t
      integer,intent(in)::m,branch
      real(dp),intent(out)::x,hx(2)
      logical,intent(out)::upper,checked
      integer,intent(out)::status
      real(dp)::lg1,lgm,vm,tm,v1,t1,phit,qa,qb,v,sv,cv,dq,d0,fx,gam,h,dh,vm0,vguess,vs,res,tmin,yv,tpp,ymin,ymax
      real(dp)::yp,c3,sh0,ch0,w0
      complex(dp)::vb
      type(coordinate)::cd
      integer i
      logical near_min,ok
      status=-4;x=ieee_value(t,ieee_quiet_nan);upper=branch>0;vguess=x;hx=x;checked=.false.
      ! T = y**1.5*h_M(v)/sqrt2+a*sqrt(y) exactly, with y monotone between
      ! y0 and y0+2b and h_(M+1)>=h_M+pi: a rigorous lower bound on t_min.
      ymin=min(g%y0,g%y0+2*g%b);ymax=max(g%y0,g%y0+2*g%b)
      h=hmin(min(m,64))+max(m-64,0)*pi
      if(t<(ymin*sqrt(ymin)*h/sqrt2+min(g%a*sqrt(ymin),g%a*sqrt(ymax)))*(1-64*epsilon(t)))then
         status=5;return
      endif
      ! Table path: classify against the tabulated minimum time and start
      ! the confluent steps from the tabulated branch root.
      call mrev_estimate(g,t,m,branch,v,vs,i)
      if(i==5)then
         status=5;return
      elseif(i==1)then
         cd%m=m
         ! Confluent steps straight from the table root; set-up for the
         ! bracket/contour fallback only if they decline.
         if(.not.upper)then
            cd%kind=2;w0=v
         else
            cd%kind=1;w0=pi-v
         endif
         call sincos(w0/2,sh0,ch0)
         if(.not.upper)then
            call confluent_step(cd,g,t,w0,sh0,ch0,0._dp,vs,x,sv,cv,ok)
         else
            call confluent_step(cd,g,t,w0,sh0,ch0,0._dp,pi-vs,x,sv,cv,ok)
         endif
         if(ok)then
            status=0;hx=[sv,cv];checked=.true.;return
         endif
         if(g%b>0)then
            gam=sqrt(g%y0/(2*g%b))
            vb=cmplx(0._dp,2*log(gam+sqrt(1+gam*gam)),dp)
         else
            gam=max(1._dp,sqrt(g%y0/(-2*g%b)))
            gam=2*log(gam+sqrt(gam*gam-1))
            vb=cmplx(pi,gam,dp)
         endif
         if(.not.upper)then
            cd%kind=2;cd%inner=vs;cd%branch=vb
            call real_core_sc(v+m*pi,v,.false.,.false.,g,fx,sv,cv,sh0,ch0)
            call fast_root(g,t,cd,v,1e-3_dp*min(v,vs-v),-1,0._dp,vs,x,status,fa0=fx-t, &
               h0=[sh0,ch0],hx=hx,checked=checked)
         else
            v=pi-v
            cd%kind=1;cd%inner=pi-vs;cd%branch=pi-vb
            call real_core_sc((m+1)*pi-v,v,.true.,.false.,g,fx,sv,cv,sh0,ch0)
            call fast_root(g,t,cd,v,1e-3_dp*min(v,pi-vs-v),-1,0._dp,pi-vs,x,status,fa0=fx-t, &
               h0=[sh0,ch0],hx=hx,checked=checked)
         endif
         if(status==0)return
         status=-4;checked=.false.;hx=ieee_value(t,ieee_quiet_nan)
      endif
      ! Model in l(v)=log(g(v)/g(0)): one logarithm per anchor.
      ! g(0)=(sqrt2/(M*pi))**(1/3)/sqrt(y0), g(pi) likewise with M+1, y(pi).
      yp=g%y0+2*g%b
      lg1=log(m/real(m+1,dp)*(g%y0/yp)*sqrt(g%y0/yp))/3
      ! Minimum of sin(v)*g(v) for log-linear g, then one exact value there.
      ! vm need only lie in (0,pi): x/(1+0.28*x**2) approximates atan(x)
      ! to 5e-3 and stays below pi/2 in magnitude.
      vm=lg1/pi;vm=pi/2+vm/(1+.28_dp*vm*vm)
      call real_core_sc(vm+m*pi,vm,.false.,.false.,g,tm,sv,cv,sh0,ch0)
      if(.not.ieee_is_finite(tm))return
      ! c3=(phi*g(0))**3 scale: phi/g(0)=cbrt(c3/T)/sin(v) without cbrt(T).
      c3=m*pi/sqrt2*g%y0*sqrt(g%y0)
      lgm=log(c3/(tm*sv**3))/3;vm0=vm
      ! Quadratic l through v=0, vm, pi; its hump maximum solves
      ! cot(v)+qa+2*qb*v=0 (model algebra only).
      call log_g_quadratic()
      ! Refining the split point only matters when t is near T(vm).
      v1=vm
      if(t<1.1_dp*tm)then
         do i=1,2
            call sincos(v1,sv,cv)
            h=cv/sv+qa+2*qb*v1;dh=-1/(sv*sv)+2*qb
            v1=min(max(v1-h/dh,.05_dp),pi-.05_dp)
         enddo
      endif
      if(abs(v1-vm)>1e-3_dp)then
         call real_core_sc(v1+m*pi,v1,.false.,.false.,g,t1,sv,cv,sh0,ch0)
         if(ieee_is_finite(t1).and.t1<tm)then
            vm=v1;tm=t1;lgm=log(c3/(tm*sv**3))/3
            call log_g_quadratic()
         endif
      endif
      vguess=ieee_value(t,ieee_quiet_nan);near_min=.false.;tmin=tm
      if(.not.(t>tm))then
         ! Undecided: locate the TOF minimum exactly (zero of the entire
         ! function K), then classify with t_min=T(v*).
         cd%kind=3;cd%m=m
         call fast_root(g,0._dp,cd,vm,1.25_dp*abs(vm-vm0)+1e-2_dp,1,0._dp,pi,vs,status,res)
         if(status/=0)then
            status=-4;return
         endif
         call real_core_sc(vs+m*pi,vs,.false.,.false.,g,tmin,sv,cv,sh0,ch0)
         if(t<tmin*(1-16*epsilon(t)))then
            status=5;return
         endif
         if(t<=tmin*(1+16*epsilon(t)))then
            x=vs;status=1;upper=.false.;return
         endif
         ! T''(v*)=K'/(2*sqrt2*sqrt(y)*s**4) with 1/K' the contour residue.
         yv=g%y0+2*g%b*sin(vs/2)**2
         tpp=1/(res*2*sqrt2*sqrt(yv)*sv**4)
         ! Quadratic model near v* only while it stays well inside the branch.
         d0=sqrt(2*(t-tmin)/abs(tpp))
         if(.not.upper.and.d0<.3_dp*vs)vguess=vs-d0
         if(upper.and.d0<.3_dp*(pi-vs))vguess=vs+d0
         vm=vs;tm=tmin;lgm=log(c3/(tmin*sv**3))/3
         call log_g_quadratic()
         status=-4;near_min=.true.
      endif
      phit=cbrt(c3/t)
      ! sqrt(y) branch point in v (node count only).
      if(g%b>0)then
         gam=sqrt(g%y0/(2*g%b))
         vb=cmplx(0._dp,2*log(gam+sqrt(1+gam*gam)),dp)
      else
         gam=max(1._dp,sqrt(g%y0/(-2*g%b)))
         gam=2*log(gam+sqrt(gam*gam-1))
         vb=cmplx(pi,gam,dp)
      endif
      cd%m=m
      ! phi/g(0) = sin(v)*exp(l(v)); estimates are model-only fixed points.
      if(.not.upper)then
         v=vm/2
         do i=1,2
            v=min(asin(min(1._dp,phit*exp(-(qa+qb*v)*v))),vm*(1-1e-3_dp))
         enddo
         if(ieee_is_finite(vguess))v=vguess
         cd%kind=2;cd%inner=vm;cd%branch=vb
         call real_core_sc(v+m*pi,v,.false.,.false.,g,fx,sv,cv,sh0,ch0);dq=qa+2*qb*v
         dq=1/abs(exp((qa+qb*v)*v)*(cv+sv*dq))
         d0=abs(phit*(t-fx)/(t+2*fx))*dq;dq=dq*phit/(3*t)
         if(near_min)call quadratic_step(vm)
         call fast_root(g,t,cd,v,1.25_dp*d0,-1,0._dp,vm,x,status,dxdt=dq,fa0=fx-t, &
            h0=[sh0,ch0],hx=hx,checked=checked)
      else
         ! Upper root in w=pi-v, w in (0,pi-vm).
         v=(pi-vm)/2
         do i=1,2
            v=min(asin(min(1._dp,phit*exp(-(qa+qb*(pi-v))*(pi-v)))),(pi-vm)*(1-1e-3_dp))
         enddo
         if(ieee_is_finite(vguess))v=pi-vguess
         cd%kind=1;cd%inner=pi-vm;cd%branch=pi-vb
         call real_core_sc((m+1)*pi-v,v,.true.,.false.,g,fx,sv,cv,sh0,ch0);dq=qa+2*qb*(pi-v)
         dq=1/abs(exp((qa+qb*(pi-v))*(pi-v))*(cv-sv*dq))
         d0=abs(phit*(t-fx)/(t+2*fx))*dq;dq=dq*phit/(3*t)
         if(near_min)call quadratic_step(pi-vm)
         call fast_root(g,t,cd,v,1.25_dp*d0,-1,0._dp,pi-vm,x,status,dxdt=dq,fa0=fx-t, &
            h0=[sh0,ch0],hx=hx,checked=checked)
      endif
   contains
      ! Near the located minimum (xs,tmin) the log-g slope is unreliable:
      ! the quadratic with zero slope at xs through (v,fx) gives the step,
      ! exact for quadratic T whatever T'' is.
      subroutine quadratic_step(xs)
         real(dp),intent(in)::xs
         real(dp)::vn
         if(.not.(fx>tmin))return
         vn=xs+(v-xs)*sqrt((t-tmin)/(fx-tmin))
         if(.not.ieee_is_finite(vn))return
         d0=abs(vn-v);dq=d0/abs(fx-t)
      end subroutine
      subroutine log_g_quadratic()
         ! l = qa*v + qb*v**2 through (vm,lgm) and (pi,lg1).
         qb=(lg1/pi-lgm/vm)/(pi-vm)
         qa=lgm/vm-qb*vm
      end subroutine
   end subroutine


   ! Zero-revolution hyperbolic transfers (t<tp). Short way: inverse cubic
   ! Hermite of z(tau), tau=T**2, between the y=0 branch point (T~a*sqrt(y))
   ! and the parabolic point. Long way: parabolic estimate in ln(tp/t), or the
   ! exp(-w/2) asymptote for t<tp/e. Contours in z near the parabola, else in
   ! w=sqrt(-z)/2.
   subroutine hyperbolic_fast(g,t,tp,der,z,status,checked)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,tp,der
      real(dp),intent(out)::z
      integer,intent(out)::status
      logical,intent(out)::checked
      real(dp)::wb,zb,yd,s,m0,m1,zg,slope,f,d0,w,lo
      type(coordinate)::cd
      logical::ok
      status=-4;z=ieee_value(t,ieee_quiet_nan);checked=.false.
      if(.not.(der>0))return
      if(g%b>0)then
         ! sinh(2*asinh(q))=2*q*sqrt(1+q*q).
         f=sqrt(g%y0/(2*g%b))
         wb=2*asinh(f);zb=-4*wb*wb
         yd=g%b*2*f*sqrt(1+f*f)/(8*wb)
         s=(t/tp)**2
         m0=tp*tp/(g%a*g%a*yd);m1=tp/(2*der)
         zg=(2*s**3-3*s*s+1)*zb+(s**3-2*s*s+s)*m0+(s**3-s*s)*m1
         slope=abs((6*s*s-6*s)*zb+(3*s*s-4*s+1)*m0+(3*s*s-2*s)*m1)/(tp*tp)*2*t
         zg=min(max(zg,zb*(1-1e-12_dp)),-tiny(t))
         lo=zb
      else
         ! Long way: w**2=tp*ln(tp/t)/(4*der) near the parabola, and the
         ! asymptote T ~ sqrt(-b)*(y0+b)*exp(-w/2) once t < tp/e
         ! (dz/dT ~ 16*w/T there).
         wb=huge(t);zb=-huge(t);lo=-huge(t)
         s=log(tp/t);w=sqrt(tp*s/(4*der));slope=1/der
         if(s>1)then
            w=max(w,2*log(sqrt(-g%b)*(g%y0+g%b)/t));slope=max(slope,16*w/t)
         endif
         zg=-4*w*w
      endif
      if(zg>-1)then
         cd%kind=0;cd%branch=cmplx(zb,0._dp,dp)
         call confluent_step_z(g,t,zg,lo,0._dp,z,ok)
         if(ok)then
            status=0;checked=.true.;return
         endif
         ! The z-step left its series range on the hyperbolic side: continue
         ! with the w-coordinate step from its last iterate.
         if(z<4*zstep_xmin.and.z>lo)then
            call confluent_step_h(g,t,sqrt(-z)/2,0._dp,wb,z,ok)
            if(ok)then
               status=0;z=-4*z*z;checked=.true.;return
            endif
         endif
         if(g%b<0)cd%branch=4*cmplx(pi,2*acosh(max(1._dp,sqrt(g%y0/(-2*g%b)))),dp)**2
         f=real_time(zg,g,0)
         d0=abs(f-t)*slope
         call fast_root(g,t,cd,zg,1.25_dp*d0+1e-3_dp*abs(zg),1,lo,0._dp,z,status,dxdt=slope,fa0=f-t,checked=checked)
      else
         cd%kind=4;cd%branch=cmplx(wb,0._dp,dp)
         w=sqrt(-zg)/2
         call confluent_step_h(g,t,w,0._dp,wb,z,ok)
         if(ok)then
            status=0;z=-4*z*z;checked=.true.;return
         endif
         f=real_time(zg,g,0)
         d0=abs(f-t)*slope/(8*w)
         call fast_root(g,t,cd,w,1.25_dp*d0+1e-3_dp*w,-1,0._dp,wb,z,status,dxdt=slope/(8*w),fa0=f-t,checked=checked)
         if(status==0)z=-4*z*z
      endif
   end subroutine


   ! Short-way hyperbolic arcs with t<0.1*tp in Y=sqrt(y): inverse cubic
   ! Hermite of Y(T) through T=0 (dY/dT=1/a) and the parabolic point
   ! (Y=sqrt(y0), dY/dT=b/(16*sqrt(y0)*T'(0))), then a bracketed contour.
   subroutine hyperbolic_short(g,t,tp,der,z,y,status)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,tp,der
      real(dp),intent(out)::z,y
      integer,intent(out)::status
      real(dp)::s,m0,m1,yp,yg,slope,f,d0,x
      type(coordinate)::cd
      logical::ok
      status=-4;z=ieee_value(t,ieee_quiet_nan);y=z
      if(.not.(der>0).or..not.(g%a>0).or..not.(g%b>0))return
      yp=sqrt(g%y0);s=t/tp
      m0=tp/g%a;m1=tp*g%b/(16*yp*der)
      yg=(s**3-2*s*s+s)*m0+(3-2*s)*s*s*yp+(s**3-s*s)*m1
      slope=abs((3*s*s-4*s+1)*m0+6*s*(1-s)*yp+(3*s*s-2*s)*m1)/tp
      yg=min(max(yg,1e-3_dp*t/g%a),yp*(1-1e-12_dp))
      call confluent_step_y(g,t,yg,0._dp,yp,x,ok)
      if(ok)then
         status=0;y=x*x;z=-16*asinh(sqrt((g%y0-y)/(2*g%b)))**2;return
      endif
      cd%kind=5;cd%branch=cmplx(sqrt(g%y0+2*g%b),0._dp,dp)
      f=real(sqrty_time(yg,0._dp,g),dp)
      d0=abs(f-t)*slope
      call fast_root(g,t,cd,yg,1.25_dp*d0+1e-6_dp*yg,1,0._dp,yp,x,status,dxdt=slope,fa0=f-t)
      if(status/=0)return
      if(abs(real(sqrty_time(x,0._dp,g),dp)-t)>1e-10_dp*t)then
         status=-4;return
      endif
      y=x*x;z=-16*asinh(sqrt((g%y0-y)/(2*g%b)))**2
   end subroutine

   !> Shared TOF evaluator for the slower, higher-order fallback contours.
   pure function legacy_time(x,g,upper) result(t)
      complex(dp),intent(in)::x
      type(geometry_data),intent(in)::g
      logical,intent(in)::upper
      complex(dp)::t
      if(upper)then
         t=transfer_time_upper(x,g)
      else
         t=transfer_time(x,g,0)
      endif
   end function

   !> Relative time-of-flight residual for a fallback candidate root.
   pure function residual(x,g,t,upper) result(r)
      real(dp),intent(in)::x,t
      type(geometry_data),intent(in)::g
      logical,intent(in)::upper
      real(dp)::r
      complex(dp)::value
      value=legacy_time(cmplx(x,0._dp,dp),g,upper)
      r=huge(t)
      if(ieee_is_finite(real(value)).and.ieee_is_finite(aimag(value)))r=abs(value-t)/t
   end function

   !> Estimate one enclosed simple root via the contour residue ratio.
   subroutine pole(g,t,c,r,n,upper,x,ok)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,c,r
      integer,intent(in)::n
      logical,intent(in)::upper
      real(dp),intent(out)::x
      logical,intent(out)::ok
      real(dp)::num,den,scale,w
      complex(dp)::e,f,term
      integer j
      ok=.false.;x=c;num=0;den=0;scale=0
      if(r<=0)return
      do j=0,n
         e=nodes(j*(768/n));w=2
         if(j==0.or.j==n)w=1
         f=legacy_time(c+r*e,g,upper)-t
         if(.not.ieee_is_finite(real(f)).or..not.ieee_is_finite(aimag(f)))return
         if(abs(f)<=tiny(t))then
            ! An endpoint may land exactly on a real root.
            if(j==0.or.j==n)then
               x=real(c+r*e,dp);ok=.true.
            endif
            return
         endif
         term=w*e/f
         den=den+real(term,dp);num=num+real(term*e,dp);scale=scale+abs(real(term,dp))
      enddo
      ! Relative cancellation gate: unlike max(1,scale), this is invariant
      ! under rescaling the TOF equation on very long transfers.
      if(abs(den)<=128*epsilon(t)*scale.or..not.ieee_is_finite(den))return
      x=c+r*num/den
      ok=ieee_is_finite(x).and.x>c-r.and.x<c+r
   end subroutine

   !> Recenter smaller contours around a candidate until residual passes.
   subroutine refine_candidate(g,t,upper,counts,tol,x,err)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,tol
      logical,intent(in)::upper
      integer,intent(in)::counts(:)
      real(dp),intent(inout)::x,err
      real(dp)::c,r,base,candidate,rr,dist
      real(dp),parameter::fq(3)=[.5_dp,.15_dp,.04_dp],fz(3)=[1._dp,.25_dp,.0625_dp]
      integer i,j
      logical ok
      c=x;dist=min(c,4*pi*pi-c);base=min(2._dp,.75_dp*dist)
      do i=1,3
         if(err<=tol)return
         if(upper)then
            r=fq(i)*min(c,4*pi*pi-c)
         else
            r=base*fz(i)
         endif
         if(r<=0)return
         do j=1,size(counts)
            call pole(g,t,c,r,counts(j),upper,candidate,ok)
            if(.not.ok)cycle
            rr=residual(candidate,g,t,upper)
            if(rr<err)then
               x=candidate;err=rr;c=candidate
            endif
            if(err<=tol)return
         enddo
      enddo
   end subroutine

   !> Try several quadrature orders, retaining the best residual.
   subroutine try_contour(g,t,c,r,counts,upper,tol,x,err)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t,c,r,tol
      integer,intent(in)::counts(:)
      logical,intent(in)::upper
      real(dp),intent(inout)::x,err
      real(dp)::candidate,rr
      integer j
      logical ok
      do j=1,size(counts)
         call pole(g,t,c,r,counts(j),upper,candidate,ok)
         if(.not.ok)cycle
         rr=residual(candidate,g,t,upper)
         call refine_candidate(g,t,upper,counts,tol,candidate,rr)
         if(rr<err)then
            x=candidate;err=rr
         endif
         if(err<=tol)return
      enddo
   end subroutine

   !> Older contour-search fallback for zero-revolution elliptic roots.
   !! Retained because fast estimates can fail on difficult geometries.
   subroutine elliptic_root(g,t,upper,fast,x,status)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t
      logical,intent(in)::upper,fast
      real(dp),intent(out)::x
      integer,intent(out)::status
      logical ok_global
      integer::local(4),global(3),nl,i,j,k
      real(dp)::tol,c,r,err,guesses(2),base,distance
      real(dp),parameter::mult(7)=[1._dp,.5_dp,2._dp,.25_dp,4._dp,.125_dp,8._dp]
      real(dp),parameter::exclusions(4)=[.5_dp,.05_dp,.005_dp,.0005_dp]
      real(dp),parameter::epsz(3)=[.5_dp,.1_dp,.02_dp]
      real(dp),parameter::boundary=4*pi*pi
      local=[8,16,32,64];nl=3;global=[6,12,24];tol=1e-11_dp
      if(.not.fast)then
         local=[16,32,64,128];nl=4;global=[16,32,64];tol=5e-13_dp
      endif
      ! A tighter internal gate protects velocities when the TOF is poorly
      ! conditioned. The public fast target remains 1e-11.
      tol=min(tol,5e-13_dp)
      status=-4;err=huge(t);x=ieee_value(t,ieee_quiet_nan)
      if(upper)then
         guesses=[(512*pi**4/sqrt2/t)**(1._dp/3)*sqrt(g%yend),boundary/4]
         do k=1,2
            do i=1,7
               c=guesses(k)*mult(i)
               if(c<=0.or.c>=boundary.or..not.ieee_is_finite(c))cycle
               r=.85_dp*min(c,boundary-c)
               call try_contour(g,t,c,r,local(:nl),upper,tol,x,err)
               if(err<=tol)then
                  status=0;return
               endif
            enddo
         enddo
         do i=1,4
            call try_contour(g,t,boundary/2,boundary/2-exclusions(i),local,upper,tol,x,err)
            if(err<=tol)then
               status=0;return
            endif
         enddo
      else
         do j=1,3
            do i=1,3
               c=boundary/2;r=c-epsz(i)
               call pole(g,t,c,r,global(j),upper,base,ok_global)
               if(.not.ok_global)cycle
               distance=residual(base,g,t,upper)
               call refine_candidate(g,t,upper,local(:nl),tol,base,distance)
               if(distance<err)then
                  x=base;err=distance
               endif
               if(err<=tol)then
                  status=0;return
               endif
            enddo
         enddo
      endif
   end subroutine
   !> Contour-search fallback for zero-revolution hyperbolic roots.
   subroutine hyperbolic_root(g,t,fast,z,status)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::t
      logical,intent(in)::fast
      real(dp),intent(out)::z
      integer,intent(out)::status
      real(dp)::lo,c,r,x,rr,local_r,refined
      integer::counts(4),a,b,j,n
      real(dp),parameter::exclusions(3)=[.01_dp,.001_dp,.0001_dp]
      logical ok
      status=-4;z=ieee_value(t,ieee_quiet_nan)
      ! The physical y=0 boundary is exact for short-way geometry; the
      ! negative-axis search supplies a bound for the long-way branch.
      if(g%b>0)then
         lo=-16*asinh(sqrt(g%y0/(2*g%b)))**2*(1-32*epsilon(t))
      else
         lo=-1
         do a=1,18
            if(real_time(lo,g,0)<t)exit
            lo=2*lo
         enddo
      endif
      if(real_time(lo,g,0)>t)return
      counts=[8,16,32,64];n=3
      if(.not.fast)then
         counts=[16,32,64,128];n=4
      endif
      c=lo/2
      do a=1,3
         r=-lo/2*(1-exclusions(a))
         do b=1,n
            call pole(g,t,c,r,counts(b),.false.,x,ok)
            if(.not.ok)cycle
            rr=residual(x,g,t,.false.)
            if(rr<=5e-13_dp)then
               z=x;status=0;return
            endif
            local_r=.5_dp*min(x-lo,-x)
            do j=1,n
               call pole(g,t,x,local_r,counts(j),.false.,refined,ok)
               if(.not.ok)cycle
               if(residual(refined,g,t,.false.)<=5e-13_dp)then
                  z=refined;status=0;return
               endif
            enddo
         enddo
      enddo
   end subroutine
   !> Solve the normalized TOF equation and return z (and optionally q).
   !! Selects elliptic, hyperbolic, or multi-revolution fast routes, validates
   !! their TOF residuals, and falls back to safeguarded contours if needed.
   !! 'half' optionally carries the half-angle pair for reconstruction; its
   !! values are NaN if the chosen solver path did not calculate them.
   subroutine solve_z(g,tof,m,branch,z,status,q_upper,fast,used_safeguard,half)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::tof
      integer,intent(in)::m,branch
      real(dp),intent(out)::z
      integer,intent(out)::status
      real(dp),optional,intent(out)::q_upper,half(2)
      logical,optional,intent(in)::fast
      logical,optional,intent(out)::used_safeguard
      real(dp)::tp,q,der,hx(2)
      logical mode,upper,checked
      call contour_initialize()
      z=ieee_value(0._dp,ieee_quiet_nan);status=-1
      if(present(q_upper))q_upper=z
      if(present(half))half=z
      if(present(used_safeguard))used_safeguard=.false.
      if(tof<=0.or..not.ieee_is_finite(tof).or.m<0.or.abs(branch)/=1)return
      mode=.true.;if(present(fast))mode=fast
      if(m==0)then
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         if(abs(tof-tp)<=8*epsilon(tp)*tp)then
            z=0;status=0;return
         endif
         ! Within rounding of the parabola no contour can resolve T-t; the
         ! linear term of the parabolic expansion is exact to O(z**2).
         der=sqrt(g%y0)*g%y0/(20*sqrt2)+g%b/8*(sqrt(g%y0)/sqrt2+g%a/(2*sqrt(g%y0)))
         if(der>0.and.abs(tof-tp)<=1e3_dp*epsilon(tp)*sqrt(g%y0)*(g%y0*sqrt2/3+abs(g%a)))then
            z=(tof-tp)/der;status=0;return
         endif
         if(tof>tp)then
            call elliptic_fast(g,tof,tp,der,q,upper,status,hx,checked)
            if(status==0.and..not.checked)then
               if(upper)then
                  if(abs(real_time_upper(q,g)-tof)>1e-10_dp*tof)status=-4
               else
                  if(abs(real_time(q,g,0)-tof)>1e-10_dp*tof)status=-4
               endif
            endif
            if(status==0)then
               z=q
               if(upper)then
                  z=4*pi*pi-q
                  if(present(q_upper))q_upper=q
               endif
               if(present(half))half=hx
               return
            endif
            call elliptic_root(g,tof,g%b<0,mode,q,status)
            if(status==0)then
               z=q
               if(g%b<0)then
                  z=4*pi*pi-q
                  if(present(q_upper))q_upper=q
               endif
               return
            endif
         else
            call hyperbolic_fast(g,tof,tp,der,z,status,checked)
            if(status==0.and..not.checked)then
               if(abs(real_time(z,g,0)-tof)>1e-10_dp*tof)status=-4
            endif
            if(status==0)return
            call hyperbolic_root(g,tof,mode,z,status)
            if(status==0)return
         endif
      else
         call multirev_fast(g,tof,m,branch,q,upper,status,hx,checked)
         if(status==5)return
         if(status==1)then
            z=4*(q+m*pi)**2;return
         endif
         if(status==0.and..not.checked)then
            if(upper)then
               if(abs(real_core((m+1)*pi-q,q,.true.,.false.,g)-tof)>1e-10_dp*tof)status=-4
            else
               if(abs(real_core(q+m*pi,q,.false.,.false.,g)-tof)>1e-10_dp*tof)status=-4
            endif
         endif
         if(status==0)then
            if(upper)then
               z=4*((m+1)*pi-q)**2
               if(present(q_upper))q_upper=4*q*(2*pi-q)
            else
               z=4*(q+m*pi)**2
            endif
            if(present(half))half=hx
            return
         endif
      endif
      ! Last resort for any case the routes above leave unresolved: the
      ! certified sign-bracket contour with bisection (safeguarded_z).
      if(present(used_safeguard))used_safeguard=.true.
      call safeguarded_z(g,tof,m,branch,z,status)
      if(.not.mode.and.(status==0.or.status==1))then
         if(abs(transfer_time(cmplx(z,0._dp,dp),g,m)-tof)>5e-13_dp*tof)status=-4
      endif
   end subroutine

   !> Short-way hyperbolic TOF in q=sqrt(y), avoiding cancellation at y=0.
   pure function time_q(q,g) result(t)
      complex(dp),intent(in)::q
      type(geometry_data),intent(in)::g
      complex(dp)::t,u,u2,f,s
      integer k
      ! q=sqrt(y) resolves extremely short hyperbolic arcs without subtracting
      ! two nearly equal numbers at the y=0 universal-variable boundary.
      u=2*asinh(sqrt((g%y0-q*q)/(2*g%b)));s=sinh(u)
      if(abs(u)<series_limit)then
         u2=-u*u;f=pser(1)
         do k=2,10
            f=f*u2+pser(k)
         enddo
         f=f/sqrt2
      else
         f=(s*cosh(u)-u)/(sqrt2*s**3)
      endif
      t=q*(q*q*f+g%a)
   end function
   !> Recover physical-space endpoint velocities from a solved anomaly z.
   !! The optional upper-boundary coordinate q_upper and cached half-angle
   !! avoid subtractive cancellation and repeated trigonometry. y_override
   !! supplies a directly solved y value for extremely short hyperbolic arcs.
   subroutine reconstruct(g,z,m,v1,v2,y_override,stable_dc,q_upper,half)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::z
      integer,intent(in)::m
      real(dp),intent(out)::v1(3),v2(3)
      real(dp),optional,intent(in)::y_override,q_upper,half(2)
      logical,optional,intent(in)::stable_dc
      real(dp)::cv,sv2,y,ratio,dc,d,vr1,vr2,vt1,vt2,u,sh,ch,e
      logical::use_quotient,have_q,have_h
      ! One half-angle evaluation: sin u=2*sh*ch and cos u=(ch-sh)*(ch+sh).
      ! Near u=pi/2 (antipodal) dc=cs-cos(u) needs the relative accuracy of
      ! cos(u) itself and dc's quotient form needs sin(u)**2=1-cos(u)**2 to
      ! the last bit, so |cos u|<1/4 recomputes both from u.
      have_q=.false.;have_h=.false.
      if(present(q_upper))have_q=ieee_is_finite(q_upper)
      if(present(half))have_h=ieee_is_finite(half(1)).and.ieee_is_finite(half(2))
      if(have_h.and.(have_q.or.z>=0))then
         sh=half(1);ch=half(2)
      elseif(have_q)then
         u=q_upper/(2*(2*pi+sqrt(4*pi*pi-q_upper)))
         call sincos(u/2,sh,ch)
      elseif(z>=0)then
         u=sqrt(z)/2-m*pi
         call sincos(u/2,sh,ch)
      endif
      if(have_q)then
         cv=-(ch-sh)*(ch+sh);sv2=(2*sh*ch)**2
         if(abs(cv)<.25_dp)then
            cv=-cos(q_upper/(2*(2*pi+sqrt(4*pi*pi-q_upper))));sv2=1-cv*cv
         endif
         if(g%b<0)then
            y=g%yend-2*g%b*sh*sh
         else
            y=g%y0+2*g%b*ch*ch
         endif
      elseif(z>=0)then
         cv=(ch-sh)*(ch+sh);sv2=(2*sh*ch)**2
         if(abs(cv)<.25_dp)then
            cv=cos(sqrt(z)/2-m*pi);sv2=1-cv*cv
         endif
         y=g%y0+2*g%b*sh*sh
         if(g%b<0)y=g%yend-2*g%b*ch*ch
      else
         ! sinh/cosh of u/2 from one expm1; cosh u=1+2*sh**2.
         u=sqrt(-z)/2;e=expm1(u/2)
         sh=.5_dp*(e+e/(1+e));ch=.5_dp*(1+e+1/(1+e))
         cv=1+2*sh*sh;sv2=-(2*sh*ch)**2;y=g%y0-2*g%b*sh*sh
      endif
      d=sign(1.0_dp,g%a);dc=d*g%cs-cv
      use_quotient=d*g%cs*cv>0
      if(present(stable_dc))then
         if(stable_dc.and.abs(d*g%cs+cv)<0.5_dp)use_quotient=.false.
      endif
      if(use_quotient)dc=(sv2-g%ss*g%ss)/(d*g%cs+cv)
      if(present(y_override))then
         y=y_override;dc=(y-g%y0)/g%b-g%ss**2/(1+g%cs)
      endif
      ratio=sqrt(g%r2/g%r1)
      vr1=sqrt2*(dc+d*g%cs*(g%r2-g%r1)/(g%r1*(ratio+1)))/sqrt(y)
      vr2=-sqrt2*(dc+d*g%cs*(g%r1-g%r2)/(g%r2*(1/ratio+1)))/sqrt(y)
      vt1=sqrt2*ratio*g%ss/sqrt(y);vt2=sqrt2*g%ss/(ratio*sqrt(y))
      v1=(vr1*g%e1+vt1*g%t1)/sqrt(g%scale)
      v2=(vr2*g%e2+vt2*g%t2)/sqrt(g%scale)
   end subroutine
   !> Solve one Lambert problem for the requested direction/revolution/branch.
   !! r1/r2 are Cartesian endpoint positions, tof is elapsed transfer time,
   !! and v1/v2 are the output velocities. The optional diagnostics report
   !! the root, TOF residual, and whether a safeguard was needed.
   subroutine contour_one(r1,r2,tof,direction,v1,v2,status,half,nrev,branch,fast,diagnostics)
      real(dp),intent(in)::r1(3),r2(3),tof
      integer,intent(in)::direction
      real(dp),intent(out)::v1(3),v2(3)
      integer,intent(out)::status,half
      integer,optional,intent(in)::nrev,branch
      logical,optional,intent(in)::fast
      type(contour_diagnostics),optional,intent(out)::diagnostics
      type(geometry_data)::g
      real(dp)::z,t,y,tp,q_upper,hs(2)
      logical safeguard
      integer m,b
      call contour_initialize()
      m=0;b=-1;if(present(nrev))m=nrev;if(present(branch))b=branch
      v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1
      if(present(diagnostics))then
         diagnostics%z=v1(1);diagnostics%q_upper=v1(1)
         diagnostics%relative_tof_residual=v1(1);diagnostics%used_safeguard=.false.
      endif
      call prepare_geometry(r1,r2,direction,g,status,half)
      if(status/=0)return
      if(m<0.or.abs(b)/=1)then
         status=-1;return
      endif
      t=tof/g%scale/sqrt(g%scale)
      if(m==0.and.g%a>0.and.ieee_is_finite(t).and.t>0)then
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         if(t<.1_dp*tp)then
            call hyperbolic_short(g,t,tp,sqrt(g%y0)*g%y0/(20*sqrt2)+ &
               g%b/8*(sqrt(g%y0)/sqrt2+g%a/(2*sqrt(g%y0))),z,y,status)
            if(status==0)then
               call reconstruct(g,z,m,v1,v2,y)
               if(present(diagnostics))then
                  diagnostics%z=z;diagnostics%used_safeguard=.true.
                  diagnostics%relative_tof_residual=abs(time_q(cmplx(sqrt(y),0._dp,dp),g)-t)/t
               endif
               if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4
               return
            endif
         endif
      endif
      call solve_z(g,t,m,b,z,status,q_upper,fast,safeguard,hs)
      if(present(diagnostics))then
         diagnostics%z=z;diagnostics%q_upper=q_upper;diagnostics%used_safeguard=safeguard
      endif
      if(status<0.or.status==5)return
      if(present(diagnostics))then
         if(ieee_is_finite(q_upper).and.m==0)then
            diagnostics%relative_tof_residual=residual(q_upper,g,t,.true.)
         else
            diagnostics%relative_tof_residual=abs(transfer_time(cmplx(z,0._dp,dp),g,m)-t)/t
         endif
      endif
      call reconstruct(g,z,m,v1,v2,q_upper=q_upper,half=hs)
      if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4
   end subroutine
   !> Solve n independent transfers (column-major 3-by-n vectors).
   !! The outer loop is OpenMP parallel; shared quadrature tables are read-only
   !! after contour_one performs thread-safe lazy initialization.
   !> Batch solve. Cases are processed in chunks of batch_chunk, and each
   !! chunk runs in stages over all of its cases (geometry and table
   !! estimate, half-angle pairs, confluent step, velocities), so that
   !! independent cases overlap in the CPU instead of each solve being one
   !! long dependent chain. Zero-revolution elliptic cases inside the
   !! tables and multi-revolution cases decided by the tables take the
   !! staged path; every other case, and any case whose step declines, is
   !! solved by contour_one. Results are identical to calling contour_one
   !! case by case. Chunks run in parallel under OpenMP.
   ! contour_one after prepare_geometry, without diagnostics: lets the
   ! staged batch reuse geometry it has already prepared.
   subroutine solve_prepared(g,tof,m,b,v1,v2,status)
      type(geometry_data),intent(in)::g
      real(dp),intent(in)::tof
      integer,intent(in)::m,b
      real(dp),intent(out)::v1(3),v2(3)
      integer,intent(out)::status
      real(dp)::z,t,y,tp,q_upper,hs(2)
      v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1
      t=tof/g%scale/sqrt(g%scale)
      if(m==0.and.g%a>0.and.ieee_is_finite(t).and.t>0)then
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         if(t<.1_dp*tp)then
            call hyperbolic_short(g,t,tp,sqrt(g%y0)*g%y0/(20*sqrt2)+ &
               g%b/8*(sqrt(g%y0)/sqrt2+g%a/(2*sqrt(g%y0))),z,y,status)
            if(status==0)then
               call reconstruct(g,z,m,v1,v2,y)
               if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4
               return
            endif
         endif
      endif
      call solve_z(g,t,m,b,z,status,q_upper,half=hs)
      if(status<0.or.status==5)return
      call reconstruct(g,z,m,v1,v2,q_upper=q_upper,half=hs)
      if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4
   end subroutine
   subroutine contour_batch(n,r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
      integer,intent(in)::n,direction(n)
      real(dp),intent(in)::r1(3,n),r2(3,n),tof(n)
      real(dp),intent(out)::v1(3,n),v2(3,n)
      integer,intent(out)::status(n),half(n)
      integer,optional,intent(in)::nrev(n),branch(n)
      integer::j
      call contour_initialize()
      ! A single chunk skips the OpenMP region: even a serialized region
      ! costs about 180 ns, and forking a team 10-50 us. Chunk times vary
      ! with the case mix (and performance/efficiency cores), hence dynamic.
      if(n<=batch_chunk)then
         call run_chunk(1)
      else
         !$omp parallel do schedule(dynamic)
         do j=1,n,batch_chunk
            call run_chunk(j)
         enddo
      endif
   contains
      subroutine run_chunk(j)
         integer,intent(in)::j
         integer::k,i,mm(batch_chunk),bb(batch_chunk)
         k=min(n,j+batch_chunk-1)
         do i=j,k
            mm(i-j+1)=0;bb(i-j+1)=-1
            if(present(nrev))mm(i-j+1)=nrev(i)
            if(present(branch))bb(i-j+1)=branch(i)
         enddo
         call staged_chunk(k-j+1,r1(:,j:k),r2(:,j:k),tof(j:k),direction(j:k),v1(:,j:k),v2(:,j:k), &
            status(j:k),half(j:k),mm,bb)
      end subroutine
   end subroutine
   subroutine staged_chunk(n,r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
      integer,intent(in)::n,direction(n),nrev(*),branch(*)
      real(dp),intent(in)::r1(3,n),r2(3,n),tof(n)
      real(dp),intent(out)::v1(3,n),v2(3,n)
      integer,intent(out)::status(n),half(n)
      real(dp),parameter::bq=4*pi*pi
      ! Fixed-size work arrays (n<=batch_chunk): automatic arrays would be
      ! heap-allocated by gfortran on every call (about 200 ns).
      type(geometry_data)::g(batch_chunk)
      real(dp)::t(batch_chunk),x(batch_chunk),sh(batch_chunk),ch(batch_chunk),hi(batch_chunk)
      integer::kd(batch_chunk)
      logical::gok(batch_chunk)
      real(dp)::zg,qg,w,sn,cn,z,q,tp,v,vs,ymin,ymax,h,nan
      type(coordinate)::cd
      integer i,m,stat
      logical ok
      nan=ieee_value(nan,ieee_quiet_nan)
      ! Stage 1: geometry, classification and table estimate.
      ! kd: 0 scalar solve, -5 infeasible, 1 upper coordinate w=pi-v, 2 v.
      do i=1,n
         kd(i)=0;gok(i)=.false.;m=nrev(i)
         if(m<0.or.abs(branch(i))/=1)cycle
         call prepare_geometry(r1(:,i),r2(:,i),direction(i),g(i),status(i),half(i))
         if(status(i)/=0)cycle
         gok(i)=.true.
         t(i)=tof(i)/g(i)%scale/sqrt(g(i)%scale)
         if(.not.(t(i)>0.and.ieee_is_finite(t(i))))cycle
         if(m==0)then
            tp=sqrt(g(i)%y0)*(g(i)%y0*sqrt2/3+g(i)%a)
            if(.not.(t(i)>tp*(1+1e-9_dp)))cycle
            if(.not.(sqrt(g(i)%y0)*g(i)%y0/(20*sqrt2)+g(i)%b/8*(sqrt(g(i)%y0)/sqrt2+ &
               g(i)%a/(2*sqrt(g(i)%y0)))>0))cycle
            zg=nan;qg=nan
            call table_estimate(g(i),cbrt(tp/t(i)),zg,qg)
            if(.not.(zg>.49_dp))cycle
            hi(i)=pi
            if(zg>bq/2)then
               w=min(max(qg,tiny(qg)),bq/2);x(i)=w/(2*(2*pi+sqrt(bq-w)));kd(i)=1
            else
               x(i)=sqrt(zg)/2;kd(i)=2
            endif
         else
            ymin=min(g(i)%y0,g(i)%y0+2*g(i)%b);ymax=max(g(i)%y0,g(i)%y0+2*g(i)%b)
            h=hmin(min(m,64))+max(m-64,0)*pi
            if(t(i)<(ymin*sqrt(ymin)*h/sqrt2+min(g(i)%a*sqrt(ymin),g(i)%a*sqrt(ymax)))*(1-64*epsilon(t)))then
               kd(i)=-5;cycle
            endif
            call mrev_estimate(g(i),t(i),m,branch(i),v,vs,stat)
            if(stat==5)then
               kd(i)=-5
            elseif(stat==1)then
               if(branch(i)>0)then
                  x(i)=pi-v;hi(i)=pi-vs;kd(i)=1
               else
                  x(i)=v;hi(i)=vs;kd(i)=2
               endif
            endif
         endif
      enddo
      ! Stage 2: half-angle pairs at the estimates.
      do i=1,n
         if(kd(i)>0)call sincos(x(i)/2,sh(i),ch(i))
      enddo
      ! Stage 3: confluent steps.
      do i=1,n
         if(kd(i)<=0)cycle
         cd%kind=kd(i);cd%m=nrev(i)
         call confluent_step(cd,g(i),t(i),x(i),sh(i),ch(i),0._dp,hi(i),w,sn,cn,ok)
         if(ok)then
            x(i)=w;sh(i)=sn;ch(i)=cn
         else
            kd(i)=0
         endif
      enddo
      ! Stage 4: velocities (as contour_one/solve_z would form them).
      do i=1,n
         m=nrev(i)
         if(kd(i)==-5)then
            v1(:,i)=nan;v2(:,i)=nan;status(i)=5
         elseif(kd(i)==1)then
            q=4*x(i)*(2*pi-x(i))
            if(m==0)then
               z=bq-q
            else
               z=4*((m+1)*pi-x(i))**2
            endif
            call reconstruct(g(i),z,m,v1(:,i),v2(:,i),q_upper=q,half=[sh(i),ch(i)])
            status(i)=0
         elseif(kd(i)==2)then
            z=4*(x(i)+m*pi)**2
            call reconstruct(g(i),z,m,v1(:,i),v2(:,i),half=[sh(i),ch(i)])
            status(i)=0
         endif
      enddo
      ! Stage 5: every other case, reusing the prepared geometry.
      do i=1,n
         if(kd(i)/=0)cycle
         if(gok(i))then
            call solve_prepared(g(i),tof(i),nrev(i),branch(i),v1(:,i),v2(:,i),status(i))
         else
            call contour_one(r1(:,i),r2(:,i),tof(i),direction(i),v1(:,i),v2(:,i),status(i),half(i),nrev(i),branch(i))
         endif
      enddo
   end subroutine
end module confluentcontour

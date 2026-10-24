! Simple contour solver (contour_simple_one) against contour_one on random
! zero-revolution elliptic, short- and long-way hyperbolic, near-parabolic
! and multi-revolution transfers, including feasibility (status 5).
program test_contour_simple
   use iso_fortran_env, only: real64
   use confluentcontour, only: contour_one, contour_initialize, pi
   use contour_simple, only: contour_simple_one
   implicit none
   integer,parameter::dp=real64
   character(*),parameter::names(5)=['elliptic    ','short hyp   ','long-way hyp','parabolic   ','multi-rev   ']
   ! Long-way hyperbolas with tof << tp are conditioned at about 1e-10 for
   ! every solver (Izzo agrees with contour_one only to 3e-10 there).
   real(dp),parameter::limit(5)=[1e-11_dp,1e-11_dp,2e-9_dp,1e-11_dp,1e-11_dp]
   real(dp)::r1(3),r2(3),tof,u(9),v1(3),v2(3),w1(3),w2(3),err,worst(5),th,cs,p,q,a,y0,rn1,rn2,sc
   integer::i,s,sr,h,dir,m,b,c,bad(5),cases(5)
   call contour_initialize()
   call random_seed(put=[(7,i=1,64)])
   worst=0;bad=0;cases=0
   do i=1,50000
      call random_number(u)
      th=6.2_dp*u(3)+0.04_dp
      r1=[1+u(1),0._dp,0._dp]
      r2=(0.3_dp+3*u(2))*[cos(th),sin(th)*cos(u(4)),sin(th)*sin(u(4))]
      dir=merge(1,-1,u(5)<0.5_dp);m=0;b=-1;c=1+int(5*u(6))
      select case(c)
      case(1); tof=0.3_dp+20*u(7)**2
      case(2); tof=10**(-3+3*u(7)); dir=1
      case(3); tof=10**(-2+2.5_dp*u(7)); dir=-1
      case(4)
         ! Within 1e-12..1 relative of the parabolic time.
         rn1=norm2(r1);rn2=norm2(r2);sc=max(rn1,rn2)
         cs=norm2(r1/rn1+r2/rn2)/2;p=sqrt(rn1*rn2)/sc
         a=dir*sqrt(2._dp)*p*cs;q=(rn1-rn2)/sc/(sqrt(rn1/sc)+sqrt(rn2/sc))
         y0=q*q+2*p*(1-cs*cs)/(1+cs);if(dir<0)y0=q*q+2*p*(1+cs)
         tof=sqrt(y0)*(y0*sqrt(2._dp)/3+a)*sc**1.5_dp*(1+(2*u(7)-1)*10**(-12*u(8)))
      case(5); m=1+int(10*u(7)**2);b=merge(1,-1,u(8)<0.5_dp);tof=(2*m+1)*pi*(0.6_dp+3*u(9))
      end select
      call contour_one(r1,r2,tof,dir,w1,w2,sr,h,m,b)
      if(sr/=0.and.sr/=5)cycle
      cases(c)=cases(c)+1
      call contour_simple_one(r1,r2,tof,dir,v1,v2,s,h,m,b)
      err=0
      if(s/=sr)then
         err=huge(1._dp)
      elseif(s==0)then
         err=max(norm2(v1-w1),norm2(v2-w2))/norm2(w1)
      endif
      worst(c)=max(worst(c),err)
      if(err>limit(c))bad(c)=bad(c)+1
   enddo
   do c=1,5
      print '(a12,a,i6,a,es9.2,a,i0)',names(c),'  cases ',cases(c),'  worst ',worst(c),'  bad ',bad(c)
   enddo
   if(any(bad>0)) error stop 'contour_simple disagrees with contour_one'
end program

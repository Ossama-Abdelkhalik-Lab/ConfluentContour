! Independent circular truth, including both singular elliptic boundaries.
! Optional output file is consumed by test_contour_matlab_parity.m.
program test_contour_port
   use iso_fortran_env, only: int64
   use confluentcontour
   use ieee_arithmetic
   implicit none
   integer,parameter::n=10000
   real(dp)::r1(3,n),r2(3,n),tof(n),v1(3,n),v2(3,n),truth1(3,n),truth2(3,n)
   real(dp)::theta,radius,u,error,maxerr,tolerance,samples(5),tmp
   integer::direction(n),status(n),half(n),i,j,k,failed,unit,ios
   integer(int64)::start,finish,rate
   character(1024)::output
   character(32)::mode
   mode="extreme"
   if(command_argument_count()>1)call get_command_argument(2,mode)
   do i=1,n
      u=real(mod(7919*i,10001),dp)/10001
      select case(mod(i,5))
       case(0)
         theta=.01_dp+(pi-.02_dp)*u
       case(1)
         theta=10._dp**(-9+7*u)
       case(2)
         theta=pi-10._dp**(-9+7*u)
       case(3)
         theta=pi+.01_dp+(pi-.02_dp)*u
       case(4)
         theta=2*pi-10._dp**(-9+7*u)
      end select
      if(trim(mode)=="ordinary")theta=.1_dp+(pi-.2_dp)*u
      radius=10._dp**(-6+12*real(mod(3571*i,10001),dp)/10001)
      r1(:,i)=[radius,0._dp,0._dp]
      r2(:,i)=radius*[cos(theta),sin(theta),0._dp]
      tof(i)=theta*radius*sqrt(radius)
      direction(i)=1;if(theta>pi)direction(i)=-1
      truth1(:,i)=[0._dp,1._dp,0._dp]/sqrt(radius)
      truth2(:,i)=[-sin(theta),cos(theta),0._dp]/sqrt(radius)
   enddo
   call contour_initialize()
   do k=1,5
      call system_clock(start,rate)
      call contour_batch(n,r1,r2,tof,direction,v1,v2,status,half)
      call system_clock(finish)
      samples(k)=real(finish-start,dp)/real(rate,dp)
   enddo
   do i=2,5
      tmp=samples(i);j=i-1
      do while(j>=1)
         if(samples(j)<=tmp)exit
         samples(j+1)=samples(j);j=j-1
      enddo
      samples(j+1)=tmp
   enddo
   maxerr=0;failed=0
   do i=1,n
      error=huge(1._dp)
      if(status(i)==0.and.all(ieee_is_finite(v1(:,i))).and.all(ieee_is_finite(v2(:,i)))) &
         error=max(norm2(v1(:,i)-truth1(:,i)),norm2(v2(:,i)-truth2(:,i)))/norm2(truth1(:,i))
      maxerr=max(maxerr,error)
      tolerance=2e-10_dp
      ! Normalizing nearly parallel endpoints can perturb their relative
      ! radii by one ulp; velocity sensitivity grows as 1/sin(theta).
      if(abs(r2(2,i))/norm2(r2(:,i))<1e-5_dp)tolerance=3e-8_dp
      if(error>tolerance)then
         failed=failed+1
         if(failed<=10)print '(a,i0,a,i0,a,es12.4)', 'case=',i,' status=',status(i),' error=',error
      endif
   enddo
   print '(a,i0,a,i0,a,es12.4)', 'circular truth: ',n-failed,'/',n,' max relative endpoint error=',maxerr
   print '(a,i0)', 'nonzero statuses: ',count(status/=0)
   print '(a,3es14.5)', 'full batch min/median/max seconds: ',samples(1),samples(3),samples(5)
   if(command_argument_count()>0)then
      call get_command_argument(1,output)
      open(newunit=unit,file=trim(output),status='replace',iostat=ios)
      if(ios/=0)error stop 'cannot open parity data'
      do i=1,n
         write(unit,'(13(es25.16,1x),2(i4,1x))')r1(:,i),r2(:,i),tof(i),v1(:,i),v2(:,i),direction(i),status(i)
      enddo
      close(unit)
   endif
   if(failed/=0)error stop 'contour circular regression'
end program


!-----Romberg积分,适合光顺函数的积分
subroutine Slam_Romberg(RomK,a,b,Ct,Num0,Point1,Intab  )

  

  implicit none
  integer(4)::RomK   !---阶次
  real(8)::a,b   !---下限，上限
  real(8)::Ct
  integer(4)::Num0    
  real(8),dimension(Num0,2)::Point1   !---轮廓线

  real(8)::Intab   !---积分结果

  real(8),dimension(RomK,RomK)::RomTkj  !---系数
  integer(4)::Mk

  

  integer(4)::i,j,k
  real(8)::s,t,s1,t1
  real(8)::h0,h1

  RomTkj=0.0;

  h1=b-a
  s=Ct*sin(a)   !----temx
  !----插值计算f(s)
  call SlamPoint_Interpolation(Num0,Point1(:,1),Point1(:,2),s,s1 )

  t=Ct*sin(b)   !----temx2
  call SlamPoint_Interpolation(Num0,Point1(:,1),Point1(:,2),t,t1 )

  RomTkj(1,1)=h1/2.0*(s1+t1)

  h0=h1;
  do i=2,RomK
    h1=h0*0.5
    Mk=2**(i-1)
    Mk=Mk/2

    t1=0.0;
    do j=1,Mk
      s=a+real(2*j-1)*h1
      s=Ct*sin(s)    !---temx

      call SlamPoint_Interpolation(Num0,Point1(:,1),Point1(:,2),s,s1 )
      t1=t1+s1
    end do
    t1=t1*h0;

    RomTkj(i,1)=(RomTkj(i-1,1)+t1 )/2.0

    do j=2,i
      RomTkj(i,j)=RomTkj(i,j-1)+(RomTkj(i,j-1)-RomTkj(i-1,j-1)  )/( 4.0**real(j-1)-1.0 )
    end do

    h0=h1;
  end do

  Intab=RomTkj(RomK,RomK)

  return
end subroutine Slam_Romberg


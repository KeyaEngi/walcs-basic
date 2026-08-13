



subroutine spline3(num,point,ds1,dsn,SSS )

   implicit none

   integer(4)::num
   real(8),dimension(num,2)::point
   real(8)::ds1,dsn   !---边界条件/一阶导数
   real(8),dimension(num-1,4)::SSS      !------系数

   real(8),dimension(num,num)::AA,invAA
   real(8),dimension(num,1)::BB
   real(8),dimension(num,1)::MM     !---------二次导数

   real(8)::modifycoef

   integer(4)::i


   !-----构建系数矩阵
   call spline3_coefficient(num,point,ds1,dsn,AA,BB(:,1)  ) 

   !-----求解线性代数方程组
   call Gauss_Jordan(num,num,AA,invAA )
   MM=matmul(invAA,BB   )

   modifycoef=1.0         !------------------修正系数（修改刚度）
   do i=1,num
       MM(i,1)=MM(i,1)*modifycoef
   end do

   !----得到各分段插值函数系数
   sss=0.0;
   call spline3_slove_function(num,point,MM,SSS )
    

   return
end subroutine spline3


subroutine spline3_coefficient(num,point,ds1,dsn,AA,d)
     
     implicit none

     integer(4)::num
     real(8),dimension(num,2)::point
     real(8)::ds1,dsn

     integer(4)::i,j,k
     real(8),dimension(1:2)::point1,point2,point3
     real(8),dimension(num)::u,v,d
     real(8),dimension(num,num)::AA
     real(8)::h1,h2,hh1,hh2
  

     do i=2,num-1
       do j=1,2
         point1(j)=point(i-1,j)
       end do
       do j=1,2
         point2(j)=point(i,j)
       end do
       do j=1,2
         point3(j)=point(i+1,j)
       end do
       h1=point2(1)-point1(1)
       h2=point3(1)-point2(1)
       u(i)=h1/(h1+h2)
       v(i)=h2/(h1+h2)
       hh1=point2(2)-point1(2)
       hh2=point3(2)-point2(2)
       d(i)=(hh2/h2-hh1/h1)/(h2+h1)*6
     end do 
  
     v(1)=1
     u(num)=1

     !-----自然边界条件，以首尾段直线斜浪作为边界一阶导
!!     d(1)=0
!!     d(num)=0               !首尾两端用直线
     
     !---守节点边界条件
     point2(:)=point(1,:);
     point3(:)=point(2,:);
     h2=point3(1)-point2(1);
     hh2=point3(2)-point2(2);

     d(1)=(hh2/h2-ds1)/h2*6.0

     !---尾节点边界条件
     point1(:)=point(num-1,:);
     point2(:)=point(num,:);
     h1=point2(1)-point1(1);
     hh1=point2(2)-point1(2);

     d(num)=(dsn-hh1/h1)/h1*6.0


     AA=0

     do i=2,num-1
       AA(i,i)=2
       AA(i,i+1)=v(i)
       AA(i,i-1)=u(i)
     end do
     AA(1,1)=2
     AA(1,2)=v(1)
     AA(num,num)=2
     AA(num,num-1)=u(num)

     return
end subroutine spline3_coefficient


subroutine spline3_slove_function(num,point,MM,SS)
    implicit none

    integer(4)::num
    real(8),dimension(num,2)::point
    real(8),dimension(num)::MM
    real(8),dimension(num-1,4)::SS

    integer(4)::i
    real(8)::a1,a2,b1,b2,c1,c2,c3,c4,d1,d2,d3,d4
    real(8)::h

    do i=1,num-1
      h=point(i+1,1)-point(i,1)
      a1=-MM(i)/(h*6)
      a2=MM(i+1)/(h*6)
      SS(i,1)=a1+a2                    !---------------三次项系数 

      b1=MM(i)*3*point(i+1,1)/(h*6)
      b2=-MM(i+1)*3*point(i,1)/(h*6)
      SS(i,2)=b1+b2                     !----------------二次项系数 

      c1=-MM(i)*3*point(i+1,1)*point(i+1,1)/(h*6)
      c2=MM(i+1)*3*point(i,1)*point(i,1)/(h*6)
      c3=-(-MM(i)*h*h/6+point(i,2))/h
      c4=(-MM(i+1)*h*h/6+point(i+1,2))/h
      SS(i,3)=c1+c2+c3+c4                     !--------------一次项系数

      d1=MM(i)*point(i+1,1)*point(i+1,1)*point(i+1,1)/(h*6)
      d2=-MM(i+1)*point(i,1)*point(i,1)*point(i,1)/(h*6)
      d3=(-MM(i)*h*h/6+point(i,2))*point(i+1,1)/h
      d4=-(-MM(i+1)*h*h/6+point(i+1,2))*point(i,1)/h
      SS(i,4)=d1+d2+d3+d4                        !-----------常数项
    end do

    return
end subroutine spline3_slove_function










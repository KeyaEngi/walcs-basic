MODULE ArrayOperations
!========================================================================================   
!Function... 
!    1.basic algrithm for the vectors and matrixs
!    
!CopyRight... 
!    1.Origional code by Tong XiaoWang, modified by Li Zhi-fu,2015.06.25 
!
!=========================================================================================
IMPLICIT NONE
!
!任何两个实型或复数型向量的点乘。
!用法：A(1:3).dot.B(1:3)
!
Interface Operator (.dot.)
  Module Procedure R_RDot
end interface
!
!任何两个实型或复数型向量的叉乘。
!用法：A(1:3).X.B(1:3)
!
Interface Operator (.X.)
  Module Procedure R_RCross
end interface
!

Interface Operator (.trans.)
  Module Procedure R_Transfer
End Interface
!!求矢量的模。
!用法：.mod.V(1:3)
!
Interface Operator (.mod.)
  Module Procedure VectorRModulus
end interface
!

!求沿某矢量的单位向量。
!用法：.unit.V(1:3)
!
Interface Operator (.unit.)
  Module Procedure Unit_VectorR
end interface
!


!矢量由于物体的转动运动所产生的旋转（到全局坐标）。
!这里是通常的线性计算。用法：R(1:3).rota.Zeta(4:6)
!
Interface Operator (.rota.)
  Module Procedure Rotate_VectorR
End Interface
!
!矢量由于物体的转动运动所产生的旋转。矢量由全局坐标转到局部坐标。
!这里是非线性计算。用法：R(1:3).rotGL.Zeta(4:6)
!
Interface Operator (.rotGL.)
  Module Procedure VectorG2L 
End Interface

!矢量由于物体的转动运动所产生的旋转。矢量由局部坐标转到全局坐标。
!这里是非线性计算。用法：R(1:3).rotLG.Zeta(4:6)
!
Interface Operator (.rotLG.)
  Module Procedure VectorL2G
End Interface



CONTAINS
!
! Dot product of two vectors   
!点乘：.dot.
!任何两个实型或复数型向量的点乘。用法：A(1:3).dot.B(1:3)
!
   FUNCTION R_RDot(v1,v2) 
     real*8 R_RDot
	 real*8, intent(in):: v1(3), v2(3)
     R_RDot = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
   END FUNCTION R_RDot  
   

!
! Cross product of two vectors, Matrix(3*3) times vector   
!叉乘：.X.
!任何两个实型或复数型向量的叉乘。用法：A(1:3).X.B(1:3)
!
   FUNCTION R_RCross(v1,v2) 
     real*8 R_RCross(3)
	 real*8, intent(in):: v1(3), v2(3)
      R_RCross(1) = v1(2) * v2(3) - v1(3) * v2(2)
	 R_RCross(2) = v1(3) * v2(1) - v1(1) * v2(3)
	 R_RCross(3) = v1(1) * v2(2) - v1(2) * v2(1)
   END FUNCTION R_RCross  



!Motion转换：.trans.
!把实型或复型的6自由度Motion从一个点转换到另一点，得到该点3个平动自由度的Motion。用法：R(1:3).trans.Zeta(1:6)，R是两点间的矢量。
!                                                                         w(1:3)+w(4:6).x.r(1:3)
   Function R_Transfer(P,zeta)
     real*8, intent(in):: P(3), Zeta(6)
     real*8 R_Transfer(3)          
     R_Transfer(1) = Zeta(1) + Zeta(5)*P(3) - Zeta(6)*P(2)
	 R_Transfer(2) = Zeta(2) + Zeta(6)*P(1) - Zeta(4)*P(3)
	 R_Transfer(3) = Zeta(3) + Zeta(4)*P(2) - Zeta(5)*P(1)   
  end function R_Transfer  

!
! Calculate modulus of a vector
!求模：.mod.
!求矢量的模。用法：.mod.V(1:3)
!
   Function VectorRModulus(v)
     real*8,intent(in):: v(3)
     real*8:: VectorRModulus
     VectorRModulus = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))      
   end function VectorRModulus
!
! Calculate the unit vector of a vector
!单位向量：.unit.
!求沿某矢量的单位向量。用法：.unit.V(1:3)
!
   FUNCTION Unit_VectorR(v) 
     real*8:: Unit_VectorR(3)
     real*8, intent(in):: v(3)
     Unit_VectorR = v/sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
   END FUNCTION Unit_VectorR
!

  
   
   
!
!矢量由于物体的转动运动所产生的旋转（到全局坐标）。
!这里是通常的线性计算。用法：R(1:3).rota.Zeta(4:6)
! 
   
   FUNCTION Rotate_VectorR(R,theta) 
     real*8:: Rotate_VectorR(3)
     real*8, intent(in):: R(3), theta(3)
     Rotate_VectorR(1) = R(1)+ theta(2) * R(3) - theta(3) * R(2)
	 Rotate_VectorR(2) = R(2)+ theta(3) * R(1) - theta(1) * R(3)
	 Rotate_VectorR(3) = R(3)+ theta(1) * R(2) - theta(2) * R(1)
   END FUNCTION Rotate_VectorR   
!
!矢量由于物体的转动运动所产生的旋转。矢量由局部坐标转到全局坐标。
!这里是非线性计算。用法：R(1:3).rotGL.Zeta(4:6)
!   
   FUNCTION VectorL2G(R,A)
     real*8:: VectorL2G(3)
     real*8, intent(in):: R(3), A(3)
     real*8:: T(3,3)
      
     T(1,1)= cos(A(2))*cos(A(3))
     T(2,1)= sin(A(1))*sin(A(2))*cos(A(3))+cos(A(1))*sin(A(3))
     T(3,1)=-cos(A(1))*sin(A(2))*cos(A(3))+sin(A(1))*sin(A(3))

     T(1,2)=-cos(A(2))*sin(A(3))
     T(2,2)=-sin(A(1))*sin(A(2))*sin(A(3))+cos(A(1))*cos(A(3))
     T(3,2)= cos(A(1))*sin(A(2))*sin(A(3))+sin(A(1))*cos(A(3))

     T(1,3)= sin(A(2))
     T(2,3)=-sin(A(1))*cos(A(2))
     T(3,3)= cos(A(1))*cos(A(2))
     VectorL2G = matmul(T,R)
   END FUNCTION VectorL2G

!
!矢量由于物体的转动运动所产生的旋转。矢量由全局坐标转到局部坐标。
!这里是非线性计算。用法：R(1:3).rotLG.Zeta(4:6)
!
   FUNCTION VectorG2L(R,A)
     real*8 VectorG2L(3)
     real*8, intent(in):: R(3), A(3)
     real*8 T(3,3)
      
     T(1,1)= cos(A(2))*cos(A(3))
     T(1,2)= sin(A(1))*sin(A(2))*cos(A(3))+cos(A(1))*sin(A(3))
     T(1,3)=-cos(A(1))*sin(A(2))*cos(A(3))+sin(A(1))*sin(A(3))

     T(2,1)=-cos(A(2))*sin(A(3))
     T(2,2)=-sin(A(1))*sin(A(2))*sin(A(3))+cos(A(1))*cos(A(3))
     T(2,3)= cos(A(1))*sin(A(2))*sin(A(3))+sin(A(1))*cos(A(3))

     T(3,1)= sin(A(2))
     T(3,2)=-sin(A(1))*cos(A(2))
     T(3,3)= cos(A(1))*cos(A(2))
     VectorG2L = matmul(T,R)
   END FUNCTION VectorG2L
   
   
     FUNCTION RotateL2G(R,A)
     real*8 RotateL2G(3)
     real*8, intent(in):: R(3), A(3)
     real*8 T(3,3)
      
     T(1,1)= cos(A(3))/cos(A(2))
     T(1,2)= -sin(A(3))/cos(A(2))
     T(1,3)=0.

     T(2,1)=sin(A(3))
     T(2,2)=cos(A(3))
     T(2,3)=0.

     T(3,1)= -cos(A(3))*tan(A(2))
     T(3,2)=sin(A(3))*tan(A(2))
     T(3,3)=1.
     RotateL2G = matmul(T,R)
   END FUNCTION RotateL2G
!
!矢量由于物体的转动运动所产生的旋转。矢量由全局坐标转到局部坐标。
!这里是非线性计算。用法：R(1:3).rotLG.Zeta(4:6)
!
   FUNCTION RotateG2L(R,A)
     real*8 RotateG2L(3)
     real*8, intent(in):: R(3), A(3)
     real*8 T(3,3)
      
     T(1,1)= cos(A(2))*cos(A(3))
     T(1,2)= sin(A(3))
     T(1,3)= 0.

     T(2,1)=-cos(A(2))*sin(A(3))
     T(2,2)=cos(A(3))
     T(2,3)= 0.

     T(3,1)=sin(A(2))
     T(3,2)=0.
     T(3,3)=1.
     RotateG2L = matmul(T,R)
   END FUNCTION RotateG2L
        
END MODULE ArrayOperations

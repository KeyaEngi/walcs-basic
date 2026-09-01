
!功能：计算梁的振动运动，用于弹性的考虑
    
subroutine SectMotion(y,t)
    use Slamming,only:SlamNumLine
    use ShipHullVar,only:MEsecLoad,MEloadr
    use Constant,only:NR,coefT
    
    implicit none
    
    real(8),dimension(2*NR)::y
    
    real(8)::t
 
    integer(4)::oj,i
    
    !*****以上是变量定义******
    
    MEsecLoad=0.0;
    do oj=1,SlamNumLine
        do i=7,NR
            MEsecLoad(1:6,oj)=MEsecLoad(1:6,oj)+MEloadr(oj,i,1:6)*y(i)
    
        end do
    end do
    !**人工赋值的梁振动系数，防止计算发散20230427**
    MEsecLoad=MEsecLoad*coefT
    
    
    
    
    
end subroutine SectMotion
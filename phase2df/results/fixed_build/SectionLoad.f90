subroutine SectionLoad( y,t,smtf,it )
    use ShipHullVar,only:Eloadr,NBSECT
    use Constant,only:NR,Nramp
    
    
    implicit none
    
    real(8),dimension(1:6,1:Nbsect)::EsecLoad   !---弹性体模态对应载荷
    
   
    
    real(8),dimension(2*NR)::y
    
    real(8)::t
    
    integer(4)::oj,i
    
    !20230522尝试添加缓载减小第一个周期剖面载荷震荡
    real(8)::smtf    
    integer(4)::it
    
    
    !******以上是变量定义
    
    EsecLoad=0.0;
    
    if( it<Nramp )then
    do oj=1,Nbsect
        do i=7,NR
            EsecLoad(1:6,oj)=EsecLoad(1:6,oj)+Eloadr(oj,i,1:6)*y(i)*smtf

        end do
    end do
    
    else
    do oj=1,Nbsect
        do i=7,NR
            EsecLoad(1:6,oj)=EsecLoad(1:6,oj)+Eloadr(oj,i,1:6)*y(i)

        end do
    end do    
 
    end if
    
   
    EsecLoad=EsecLoad
    
    

    
    !**输出剖面载荷**
    write(24,"(f9.3,3x,\ )",advance='NO') t
    write(25,"(f9.3,3x,\ )",advance='NO') t
    write(26,"(f9.3,3x,\ )",advance='NO') t
    write(27,"(f9.3,3x,\ )",advance='NO') t
    write(28,"(f9.3,3x,\ )",advance='NO') t
    write(29,"(f9.3,3x,\ )",advance='NO') t
    
    do i=1,NBSECT
        
    
        write(24,"(e12.4,3x,\ )",advance='NO') eSecLoad(1,i)
        write(25,"(e12.4,3x,\ )",advance='NO') eSecLoad(2,i)
        write(26,"(e12.4,3x,\ )",advance='NO') eSecLoad(3,i)
        write(27,"(e12.4,3x,\ )",advance='NO') eSecLoad(4,i)
        write(28,"(e12.4,3x,\ )",advance='NO') eSecLoad(5,i)
        write(29,"(e12.4,3x,\ )",advance='NO') eSecLoad(6,i)
    
    end do
    

    
    write(24,"(/)"); write(25,"(/)"); write(26,"(/)")
    write(27,"(/)"); write(28,"(/)"); write(29,"(/)")
    
    
    
    
    
    
end subroutine SectionLoad
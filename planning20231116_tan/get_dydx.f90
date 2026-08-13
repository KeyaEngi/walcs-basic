subroutine get_dydx( it,ca,air_num,air_point,dydx )
    use lift,only:dydx
    use Constant,only:pi
    implicit none
    
    integer(4)::it
    integer(4)::air_num
    real(8)::air_point(air_num,2)
    
    real(8)::dydx(it)
    
    integer(4)::i,j,k
    
    real(8)::xx(it),yy(it)
    
    real(8)::kexl
    
    real(8)::ca
    
    !先找baita对应的x
    do i=1,it        
        xx(i)=ca*( 1-cos(i/real(it)*pi) )/2        
    end do
    
    !插值找对应的y,并计算dydx
    dydx=0.0;
    yy=0.0;
    
    do i=1,it
        do j=1,air_num-1
            if( xx(i)>=air_point(j,1).and.xx(i)<air_point(j+1,1) )then
                
                kexl=( xx(i)-air_point(j,1) )/( air_point(j+1,1)-air_point(j,1) )
                
                yy(i)=air_point(j,2)*(1.0-kexl)+air_point(j+1,2)*kexl
                
                dydx(i)=( yy(i)-air_point(j,2) )/( xx(i)-air_point(j,1) )
                
            end if   
        end do
    end do
    
 
end subroutine get_dydx
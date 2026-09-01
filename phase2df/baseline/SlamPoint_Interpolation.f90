
!---½Úµã²åÖµ

subroutine SlamPoint_Interpolation(Num,Px,Py,x,y   )

  implicit none
  integer(4)::Num
  real(8),dimension(Num)::Px,Py
  real(8)::x,y

  integer(4)::i,j,k
  real(8)::s,t,s1,s2

  if(x<Px(1) ) then
    y=Py(1)
  elseif(x>Px(Num) ) then
    y=Py(Num)
  else
    do i=1,Num-1
      s1=Px(i); s2=Px(i+1);

      if( (x-s1)*(x-s2)<=0.0 ) then
        t=(x-s1)/(s2-s1)
        y=Py(i)*(1.0-t)+Py(i+1)*t
        exit
      end if

    end do
  end if



  return
end subroutine SlamPoint_Interpolation



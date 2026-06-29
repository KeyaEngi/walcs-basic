!---------------------------------------------------------
!用高斯消去法解方程组子程序,修改：参数h改为hh /by sun
!---------------------------------------------------------
subroutine Gauss(c,n,nm)
!列选主元Gauss消去法解代数方程组
!最后c(1:n,n+1:nm)中存放方程组的解
implicit none
integer(4)::n,nm	!n--增广矩阵的行数    nm--增广矩阵的列数
real(8),dimension(n,nm)::c	!增广矩阵
integer(4)::i,j,k,l,m,n1	!局部工作单元
real(8)::hh,g				!局部工作单元
n1=n+1
do k=1,n
	g=0.0
	do i=k,n
		hh=c(i,k)
		if(abs(hh)<=abs(g))cycle
		g=hh
		m=i
	end do			  !以上为选主元过程
	if(m/=k)then
		do j=k,nm
			hh=c(k,j)
			c(k,j)=c(m,j)
			c(m,j)=hh
		end do
	end if			  !以上为行交换
	g=1.0/g
	l=k+1
	do j=l,nm
		hh=c(k,j)*g
		c(k,j)=hh
		if(l>n)cycle
		do i=l,n
			c(i,j)=c(i,j)-c(i,k)*hh
		end do
	end do
end do				 !以上为消元
do l=1,n
	i=n1-l
	do j=i+1,n
		hh=c(i,j)
		do k=n1,nm
			c(i,k)=c(i,k)-c(j,k)*hh
		end do
	end do
end do			   !回代
end subroutine Gauss
!--------------------------------------------------------

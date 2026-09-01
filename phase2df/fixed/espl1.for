!同spline.for,只是参数列表有变化
!=============================================================================================================================
!=============================================================================================================================
!程序功能:三次样条插值微分
!变量说明:(X,y)给定点坐标;n,给定点数目;n,给定点数目;dy1,dyn,两端点一阶导数值;xx,指定插值点纵坐标;m,指定插值点数目;
!         DY,DDY,给定结点处的一阶导数与二阶导数;S,DS,DDS,指定插值点处的函数值、一阶和二阶导数值;T,[x1,xn]区间的积分值;
!         H,本子程序的工作数组.
!根据本程序需要改变参数列表 2005.4.18
!=============================================================================================================================
      subroutine espl1(x,y,n,dy1,dyn,xx,m,s,ds,dds)
	dimension x(n),y(n),xx(m),dy(n),ddy(n)
	dimension s(m),ds(m),dds(m),h(n)
	double precision x,y,xx,dy,ddy,s,ds,dds,h,dy1,dyn,t,h0,h1,beta,
     +                 alpha 
	dy(1)=-0.5
	h0=x(2)-x(1)
	h(1)=3.0*(y(2)-y(1))/(2.0*h0)-dy1*h0/4.0
	do 10 j=2,n-1
		h1=x(j+1)-x(j)
		alpha=h0/(h0+h1)
		beta=(1.0-alpha)*(y(j)-y(j-1))/h0
		beta=3.0*(beta+alpha*(y(j+1)-y(j))/h1)
		dy(j)=-alpha/(2.0+(1.0-alpha)*dy(j-1))
		h(j)=(beta-(1.0-alpha)*h(j-1))
		h(j)=h(j)/(2.0+(1.0-alpha)*dy(j-1))
		h0=h1
10	continue
	dy(n)=(3.0*(y(n)-y(n-1))/h1+dyn*h1/2.0-h(n-1))/(2.0+dy(n-1))
	do 20 j=n-1,1,-1
20	dy(j)=dy(j)*dy(j+1)+h(j)
	do 30 j=1,n-1
30	h(j)=x(j+1)-x(j)
	do 40 j=1,n-1
		h1=h(j)*h(j)
		ddy(j)=6.0*(y(j+1)-y(j))/h1-2.0*(2.0*dy(j)+dy(j+1))/h(j)
40	continue	
	h1=h(n-1)*h(n-1)
	ddy(n)=6.0*(y(n-1)-y(n))/h1+2.0*(2.0*dy(n)+dy(n-1))/h(n-1)
	t=0.0
	do 50 i=1,n-1
		h1=0.5*h(i)*(y(i)+y(i+1))
		h1=h1-h(i)*h(i)*h(i)*(ddy(i)+ddy(i+1))/24.0
		t=t+h1
50	continue
	do 70 j=1,m
		if(xx(j).ge.x(n))then
			i=n-1
		else
			i=1
60			if(xx(j).gt.x(i+1))then
				i=i+1
				goto 60
			endif
		endif
		h1=(x(i+1)-xx(j))/h(i)
		s(j)=(3.0*h1*h1-2.0*h1*h1*h1)*y(i)
		s(j)=s(j)+h(i)*(h1*h1-h1*h1*h1)*dy(i)
		ds(j)=6.0*(h1*h1-h1)*y(i)/h(i)
		ds(j)=ds(j)+(3.0*h1*h1-2.0*h1)*dy(i)
		dds(j)=(6.0-12.0*h1)*y(i)/(h(i)*h(i))
		dds(j)=dds(j)+(2.0-6.0*h1)*dy(i)/h(i)
		h1=(xx(j)-x(i))/h(i)
		s(j)=s(j)+(3.0*h1*h1-2.0*h1*h1*h1)*y(i+1)
		s(j)=s(j)-h(i)*(h1*h1-h1*h1*h1)*dy(i+1)
		ds(j)=ds(j)-6.0*(h1*h1-h1)*y(i+1)/h(i)
		ds(j)=ds(j)+(3.0*h1*h1-2.0*h1)*dy(i+1)
		dds(j)=dds(j)+(6.0-12.0*h1)*y(i+1)/(h(i)*h(i))
		dds(j)=dds(j)-(2.0-6.0*h1)*dy(i+1)/h(i)
70	continue
	return
	end
							

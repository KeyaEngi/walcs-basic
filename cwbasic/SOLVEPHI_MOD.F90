!-------------------------------------------------------------------
!用列主元GAUSS消去法或GMRES法求解速度势PHI(1:2,1:4,1:7,1:NAPNEL)
!模块编制：孙葳 2011.12.9
!-------------------------------------------------------------------
MODULE SOLVEPHI_MOD
USE MESH_MOD,ONLY:NPanel,NB,xav,e,Factor,TarNum,Space
USE ENVIRONMENT_MOD,ONLY:ZG,H,G0
USE CAL_MOD,ONLY:WATERDEPTH
USE GMRES_MOD
USE INOUTACCESS_MOD

CONTAINS
!*****************************************************************************
!程序功能：1.形成用于计算扰动流场源强的线性方程组的系数矩阵A和用来计算速度势
!          及其梯度的系数矩阵B. 
!          2.生成右端项并求解扰动流场的源强，计算各面元中心速度势梯度及速度势  
!输入参数：GcD(i,j),GsD(i,j)为j面元对i面元格林函数法向诱导速度的实和虚部
!          GcS(i,j,1:4)，GsS(i,j,1:4)为j面元对i面元格林函数诱导速度及其本身的
!          实部和虚部。 
!输出参数：Phi(1:2,1:4,1:7,1:NPANEL)为各面元中心速度势梯度及速度势        
!调用子程序：Gauss           
!******************************************************************************
subroutine solve(GcD,GsD,GcS,GsS,omeI0,HEADI,k0,Pnu,PHI,SSTRENGTH)
implicit none

! integer(4),intent(in)::id
real(8),intent(in)::omeI0,HEADI,k0,Pnu
real(8),INTENT(IN),DIMENSION(1:NPANEL,1:NPANEL)::GcD,GsD     
real(8),INTENT(IN),DIMENSION(1:NPANEL,1:NPANEL,1:4)::GcS,GsS 
REAL(8),INTENT(OUT),DIMENSION(1:2,1:4,1:7,1:NPANEL) ::PHI
REAL(8),INTENT(OUT),DIMENSION(1:2*NPANEL,1:7) ::SSTRENGTH

real(8)::xh,c1,c2,c3,c4,c5,c6
integer(4)::i,j,k,L

real(8),allocatable,dimension(:,:)::A
real(8),allocatable,dimension(:,:,:)::B
real(8),allocatable,dimension(:,:,:)::D

ALLOCATE(A(1:2*NPanel,1:2*NPanel+7),B(1:2*NPanel,1:2*NPanel,1:4))
allocate(D(1:2*NPanel,1:4,1:7))

A=0.0;B=0.0
do i=1,NPanel
	A(i,1:NPanel)=GcD(i,1:NPanel)
	A(i,NPanel+1:2*NPanel)=-GsD(i,1:NPanel)
	A(NPanel+i,1:NPanel)=GsD(i,1:NPanel)
	A(NPanel+i,NPanel+1:2*NPanel)=GcD(i,1:NPanel)
	do j=1,4
		B(i,1:NPanel,j)=GcS(i,1:NPanel,j)
		B(i,NPanel+1:2*NPanel,j)=-GsS(i,1:NPanel,j)
		B(NPanel+i,1:NPanel,j)=GsS(i,1:NPanel,j)
		B(NPanel+i,NPanel+1:2*NPanel,j)=GcS(i,1:NPanel,j)
	enddo
enddo

! A(:,2*NPanel+1:2*NPanel+7)=0.0

do i=1,NB
	A(i,2*NPanel+1:2*NPanel+3)=e(1:3,3,i)
	A(i,2*NPanel+4)=xav(2,i)*e(3,3,i)-(xav(3,i)-ZG)*e(2,3,i)
	A(i,2*NPanel+5)=(xav(3,i)-ZG)*e(1,3,i)-xav(1,i)*e(3,3,i)

	A(i,2*NPanel+6)=xav(1,i)*e(2,3,i)-xav(2,i)*e(1,3,i)
	if(abs(omeI0)<1e-10)then
		c1=0.0
	else
		c1=g0*k0/omeI0
	endif
	xh=Pnu*H
	if((xh>7).or.(trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		c2=exp(Pnu*(xav(3,i)))
		c3=c2
	else
		c2=cosh(k0*(xav(3,i)+H))/cosh(k0*H)
		c3=sinh(k0*(xav(3,i)+H))/cosh(k0*H)
	endif
	c4=e(1,3,i)*cos(HEADI)-e(2,3,i)*sin(HEADI)
	c5=xav(1,i)*cos(HEADI)-xav(2,i)*sin(HEADI)
	c6=k0*c5
	A(i,2*NPanel+7)=-(-c1*c2*c4*cos(c6)-c1*c3*e(3,3,i)*sin(c6))
	A(i+NPanel,2*NPanel+7)=-(-c1*c2*c4*sin(c6)+c1*c3*e(3,3,i)*cos(c6))
enddo

call Gauss(A,2*NPanel,2*NPanel+7)

SSTRENGTH(:,1:7)=A(:,2*NPANEL+1:2*NPANEL+7)  !输出源强 by sun

D=0.0
do L=1,7
	do K=1,4	
			do i=1,2*NPanel
				do j=1,2*NPanel
					D(i,K,L)=D(i,K,L)+B(i,j,k)*A(j,2*NPanel+L)
				enddo
			enddo
	enddo
enddo

DO I=1,NPANEL
   do L=1,7
	   do K=1,4
	   	   Phi(1,K,L,I)=D(I,K,L)
		   Phi(2,K,L,I)=D(NPanel+I,K,L)
	   enddo
   enddo
END DO

deallocate(A,B,D)

end subroutine solve
!-----------------------------------------------------------------------




!*****************************************************************************************
!程序功能：应用GMRES方法生成右端项并求解扰动流场的源强，计算各面元中心速度势梯度及速度势     
!注意：m-重启动数（取50至100间，越大越容易收敛，同时计算效率变慢）
!调用子程序：BGMRES              
!*****************************************************************************************
subroutine solve_gmres(G0Dc,G0Ds,omeI0,HEADI,k0,Pnu,Phi,SSTRENGTH)
implicit none
CHARACTER(LEN=500)::BASS
real(8),intent(in),dimension(1:NPanel,1:NPanel,1:3,1:Factor)::G0Dc,G0Ds
real(8),intent(in)::omeI0,HEADI,k0,Pnu
real(8),intent(out),dimension(1:2,1:2,1:7,1:Factor,1:NB)::Phi !(Phi(实部、虚部，？？，辐射1~6绕射7，对称性因子，湿表面序号）
REAL(8),INTENT(OUT),DIMENSION(1:2,1:7,1:Factor,1:NPANEL) ::SSTRENGTH !源强（实部、虚部，1~6辐射、7绕射，对称因子，总面元数）
complex,allocatable,dimension(:,:)::BB,KK
complex,allocatable,dimension(:,:)::RK
complex,allocatable,dimension(:,:,:)::A1
complex,allocatable,dimension(:)::N1,N2,N3,N4 !对称性
complex,allocatable,dimension(:,:)::Q
real(8)::xh,c0,c01,c1,c2,c3,c4
real(8)::c40,c50,c60 !仅仅供无对称性使用
complex::c5,c6,c7,e1,e2,e3,e4,m6
integer(4)::i,j,k,L,jj,k2,k3,k4
integer::M,ITR !重启动数和迭代次数
allocate(A1(1:NPanel,1:NPanel,1:Factor))
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if(factor==1)then !不对称
        allocate(KK(1:NPanel,1:NPanel),BB(1:NPanel,1:7),Q(1:NPanel,1:7))

                    BB=0.0
                    M=50
                    A1=0.0
                    A1(:,:,1)=cmplx(G0Dc(:,:,3,1),G0Ds(:,:,3,1))
                    phi=0.0

        do i=1,NB
	                if(abs(omeI0)<1e-10)then
		                c0=0.0
	                else
		                c0=g0*k0/omeI0
	                endif
	                xh=Pnu*H
	                if(xh>7.or.(trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		                c1=exp(Pnu*(xav(3,i)))
		                c2=c1
		                c01=c0*c1	
                    else
                    c1=cosh(k0*(xav(3,i)+H))/cosh(k0*H)
                    c2=sinh(k0*(xav(3,i)+H))/cosh(k0*H)
                    c01=c0*c1
                    endif
	                c40=e(1,3,i)*cos(HEADI)-e(2,3,i)*sin(HEADI)
	                c50=xav(1,i)*cos(HEADI)-xav(2,i)*sin(HEADI)
	                c60=k0*c50
	                BB(i,7)=cmplx(-(-c0*c1*c40*cos(c60)-c0*c2*e(3,3,i)*sin(c60)),-(-c0*c1*c40*sin(c60)+c0*c2*e(3,3,i)*cos(c60)))
	                BB(i,1:3)=cmplx(e(1:3,3,i),0.0)
	                BB(i,4)=cmplx(xav(2,i)*e(3,3,i)-(xav(3,i)-ZG)*e(2,3,i),0.0)
	                BB(i,5)=cmplx((xav(3,i)-ZG)*e(1,3,i)-xav(1,i)*e(3,3,i),0.0)
	                BB(i,6)=cmplx(xav(1,i)*e(2,3,i)-xav(2,i)*e(1,3,i),0.0)
        enddo
        !A为系数矩阵，BB（i）为b矩阵，Q为求得解，Npanel为系数矩阵的维数。
                    KK=A1(:,:,1)
                    call bgmres(KK,BB(:,1),Q(:,1),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,2),Q(:,2),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,3),Q(:,3),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,4),Q(:,4),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,5),Q(:,5),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,6),Q(:,6),NPanel,1,M,ITR)
                    call bgmres(KK,BB(:,7),Q(:,7),NPanel,1,M,ITR)


                     !源强对称性输出：
        do i=1,7
               SSTRENGTH(1,i,1,1:NPANEL)=REAL(Q(1:NPANEL,i))    !源强实部
               SSTRENGTH(2,i,1,1:NPANEL)=AIMAG(Q(1:NPANEL,i))   !源强虚部
        end do
            
        Phi=0.0
        do L=1,7
            do jj=1,2
	            do i=1,NB
				          do j=1,NPanel
				            Phi(1,jj,L,1,I)=Phi(1,jj,L,1,I)+G0Dc(i,j,jj,1)*SSTRENGTH(1,L,1,j)-G0Ds(i,j,jj,1)*SSTRENGTH(2,L,1,j) !(实部，？？，辐射1~6绕射7，对称性因子，湿表面序号
				            Phi(2,jj,L,1,I)=Phi(2,jj,L,1,I)+G0Dc(i,j,jj,1)*SSTRENGTH(2,L,1,j)+G0Ds(i,j,jj,1)*SSTRENGTH(1,L,1,j)!(虚部，？？，1~6绕射7，对称性因子，湿表面序号
			        enddo
		        enddo
	        enddo
        enddo



        allocate(RK(1:1,1),N1(1:1),N2(1:1),N3(1:1),N4(1:1))
        deallocate(BB,Q,RK,N1,N2,N3,N4)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
else if(Factor==2)then
        allocate(KK(1:NPanel,1:NPanel),BB(1:NPanel,1:6),Q(1:NPanel,1:6))
 
                    BB=0.0
                    M=100
                    A1=0.0
                    A1(:,:,1)=cmplx(G0Dc(:,:,3,1),G0Ds(:,:,3,1))
                    A1(:,:,2)=cmplx(G0Dc(:,:,3,2),G0Ds(:,:,3,2))
                    phi=0.0
        do i=1,NB
	                BB(i,1:3)=cmplx(e(1:3,3,i),0.0)
	                BB(i,4)=cmplx(xav(2,i)*e(3,3,i)-(xav(3,i)-ZG)*e(2,3,i),0.0)
	                BB(i,5)=cmplx((xav(3,i)-ZG)*e(1,3,i)-xav(1,i)*e(3,3,i),0.0)
	                BB(i,6)=cmplx(xav(1,i)*e(2,3,i)-xav(2,i)*e(1,3,i),0.0)

        enddo
        !A为系数矩阵，BB（i）为b矩阵，Q为求得解，Npanel为系数矩阵的维数。
        KK=A1(:,:,1)+A1(:,:,2)
        !call LSACG(KK,BB(:,1),Q(:,1))
        call bgmres(KK,BB(:,1),Q(:,1),NPanel,1,M,ITR)

        call bgmres(KK,BB(:,3),Q(:,3),NPanel,1,M,ITR)

        call bgmres(KK,BB(:,5),Q(:,5),NPanel,1,M,ITR)

        KK=A1(:,:,1)-A1(:,:,2)

        call bgmres(KK,BB(:,2),Q(:,2),NPanel,1,M,ITR)

        call bgmres(KK,BB(:,4),Q(:,4),NPanel,1,M,ITR)

        call bgmres(KK,BB(:,6),Q(:,6),NPanel,1,M,ITR)

         !源强对称性输出：
 
         !!surge    + - - +
               SSTRENGTH(1,1,1,1:NPANEL)=REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,1,1:NPANEL)=AIMAG(Q(1:NPANEL,1))

         
               SSTRENGTH(1,1,2,1:NPANEL)=REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,2,1:NPANEL)=AIMAG(Q(1:NPANEL,1))   
       
        !!sway    +  +  -  -
               SSTRENGTH(1,2,1,1:NPANEL)=REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,1,1:NPANEL)=AIMAG(Q(1:NPANEL,2))
       
         
               SSTRENGTH(1,2,2,1:NPANEL)=-REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,2))   
       
       
         !!heave   +   +   +   +
               SSTRENGTH(1,3,1,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,1,1:NPANEL)=AIMAG(Q(1:NPANEL,3))
        
         
               SSTRENGTH(1,3,2,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,2,1:NPANEL)=AIMAG(Q(1:NPANEL,3))   
       
       
        !!roll   +   +   -   -
               SSTRENGTH(1,4,1,1:NPANEL)=REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,1,1:NPANEL)=AIMAG(Q(1:NPANEL,4))
         
               SSTRENGTH(1,4,2,1:NPANEL)=-REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,4))   
       
       
        !!pitch   +   -   -   +
               SSTRENGTH(1,5,1,1:NPANEL)=REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,1,1:NPANEL)=AIMAG(Q(1:NPANEL,5))

               SSTRENGTH(1,5,2,1:NPANEL)=REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,2,1:NPANEL)=AIMAG(Q(1:NPANEL,5))   
       
       
         !!yaw   +   -   -   +
               SSTRENGTH(1,6,1,1:NPANEL)=REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,1,1:NPANEL)=AIMAG(Q(1:NPANEL,6))
         
               SSTRENGTH(1,6,2,1:NPANEL)=-REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,6))   
       
        Phi=0.0
        do L=1,6
            do jj=1,2
	            do i=1,NB
			        do K=1,2

				          do j=1,NPanel
				            Phi(1,jj,L,1,I)=Phi(1,jj,L,1,I)+G0Dc(i,j,jj,k)*SSTRENGTH(1,L,k,j)-G0Ds(i,j,jj,k)*SSTRENGTH(2,L,k,j)
				            Phi(2,jj,L,1,I)=Phi(2,jj,L,1,I)+G0Dc(i,j,jj,k)*SSTRENGTH(2,L,k,j)+G0Ds(i,j,jj,k)*SSTRENGTH(1,L,k,j)
				    				    
				        enddo
			        enddo
		        enddo
	        enddo
        enddo

        Phi(:,:,1,2,:)=Phi(:,:,1,1,:)
        !sway
        Phi(:,:,2,2,:)=-Phi(:,:,2,1,:)
        !heave
        Phi(:,:,3,2,:)=Phi(:,:,3,1,:)
        !roll
        Phi(:,:,4,2,:)=-Phi(:,:,4,1,:)
        !pitch
        Phi(:,:,5,2,:)=Phi(:,:,5,1,:)
        !yaw
        Phi(:,:,6,2,:)=-Phi(:,:,6,1,:)
        deallocate(BB,Q)
        allocate(RK(1:Npanel,1:2),N1(1:Npanel),N4(1:Npanel),N3(1:Npanel),N2(1:Npanel))
        rk=0
        N1=0;N2=0;N3=0;N4=0
        do i=1,NB

        c0=0.0;c1=0.0;c2=0.0;c3=0.0
        c4=0.0;c5=0.0;c6=0.0;c7=0.0
        c01=0.0;e1=0.0;e2=0.0;e3=0.0;e4=0.0
	        if(abs(omeI0)<1e-10)then
		        c0=0.0
	        else
		        c0=g0/omeI0
	        endif
	        xh=Pnu*H
	        if(xh>7.or.(trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		        c1=exp(Pnu*(xav(3,i)))
		        c01=c0*c1
	        c3=k0*xav(1,i)*cos(HEADI)   !k0*x*cos(beta)
	        c4=k0*xav(2,i)*sin(HEADI)    !k0*y*sin(beta)
	        c5=cmplx(k0*e(3,3,i),0)
	        c6=cmplx(0,k0*cos(Headi)*e(1,3,i))
	        c7=cmplx(0,k0*sin(Headi)*e(2,3,i))
	        e1=cmplx(0,cos(c3)*cos(c4)*c01)
	        e2=cmplx(0,sin(c3)*sin(c4)*c01)
	        e3=cmplx(-sin(c3)*cos(c4)*c01,0)
	        e4=cmplx(cos(c3)*sin(c4)*c01,0)
	        N1(i)=e1*c5+e3*c6-e4*c7
	        N2(i)=e2*c5+e4*c6-e3*c7
	        N3(i)=e3*c5+e1*c6-e2*c7
	        N4(i)=e4*c5+e2*c6-e1*c7	
		
		
	        else
            c1=cosh(k0*(xav(3,i)+H))/cosh(k0*H)
	        c2=sinh(k0*(xav(3,i)+H))/cosh(k0*H)
	        c01=c0*c1
	        c3=k0*xav(1,i)*cos(HEADI)   !k0*x*cos(beta)
	        c4=k0*xav(2,i)*sin(HEADI)    !k0*y*sin(beta)
	        c5=cmplx(k0*e(3,3,i),0)
	        c6=cmplx(0,k0*cos(Headi)*e(1,3,i))
	        c7=cmplx(0,k0*sin(Headi)*e(2,3,i))
	        e1=cmplx(0,cos(c3)*cos(c4)*c01)
	        e2=cmplx(0,sin(c3)*sin(c4)*c01)
	        e3=cmplx(-sin(c3)*cos(c4)*c01,0)
	        e4=cmplx(cos(c3)*sin(c4)*c01,0)
	        N1(i)=e1*c2*c5/c1+e3*c6-e4*c7
	        N2(i)=e2*c2*c5/c1+e4*c6-e3*c7
	        N3(i)=e3*c2*c5/c1+e1*c6-e2*c7
	        N4(i)=e4*c2*c5/c1+e2*c6-e1*c7
	        endif
        end do
            KK=A1(:,:,1)+A1(:,:,2)
            call bgmres(KK,-N1-N3,RK(:,1),NPanel,1,M,ITR)
            KK=A1(:,:,1)-A1(:,:,2)
            call bgmres(KK,-N2-N4,RK(:,2),NPanel,1,M,ITR)


        !输出绕射源强

               SSTRENGTH(1,7,1,1:NPANEL)=REAL(RK(1:NPANEL,1)+RK(1:NPANEL,2))
               SSTRENGTH(2,7,1,1:NPANEL)=AIMAG(RK(1:NPANEL,1)+RK(1:NPANEL,2))
              
               SSTRENGTH(1,7,2,1:NPANEL)=REAL(RK(1:NPANEL,1)-RK(1:NPANEL,2))
               SSTRENGTH(2,7,2,1:NPANEL)=AIMAG(RK(1:NPANEL,1)-RK(1:NPANEL,2))

        !输出绕射源强
			        do i=1,NB
			            do k=1,2
          select case(k)
			          case(1)
			             k4=2
			             case(2)
			             k4=1
			             end select
				        do j=1,NPanel
				        !速度势的实部和虚部
					        PHI(1,1,7,1,i)=PHI(1,1,7,1,i)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
					        PHI(2,1,7,1,i)=PHI(2,1,7,1,i)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
				    
				            Phi(1,1,7,2,I)=Phi(1,1,7,2,I)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				            Phi(2,1,7,2,I)=Phi(2,1,7,2,I)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				        !速度势的x方向导数实部和虚部
					        PHI(1,2,7,1,i)=PHI(1,2,7,1,i)+real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
					        PHI(2,2,7,1,i)=PHI(2,2,7,1,i)+AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
				    
				            Phi(1,2,7,2,I)=Phi(1,2,7,2,I)+real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				            Phi(2,2,7,2,I)=Phi(2,2,7,2,I)+AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))				    			    
				        enddo
				        end do
	        enddo
else if(factor==4)then
        allocate(KK(1:NPanel,1:NPanel),BB(1:NPanel,1:7),Q(1:NPanel,1:7))
        !m-重启动数（取50至100间）
        BB=0.0
        M=100
        A1=0.0
        A1(:,:,1)=cmplx(G0Dc(:,:,3,1),G0Ds(:,:,3,1))
        A1(:,:,2)=cmplx(G0Dc(:,:,3,2),G0Ds(:,:,3,2))
        A1(:,:,3)=cmplx(G0Dc(:,:,3,3),G0Ds(:,:,3,3))
        A1(:,:,4)=cmplx(G0Dc(:,:,3,4),G0Ds(:,:,3,4))
        phi=0.0
        do i=1,NB

	        BB(i,1:3)=cmplx(e(1:3,3,i),0.0)
	        BB(i,4)=cmplx(xav(2,i)*e(3,3,i)-(xav(3,i)-ZG)*e(2,3,i),0.0)
	        BB(i,5)=cmplx((xav(3,i)-ZG)*e(1,3,i)-xav(1,i)*e(3,3,i),0.0)
	        BB(i,6)=cmplx(xav(1,i)*e(2,3,i)-xav(2,i)*e(1,3,i),0.0)

        enddo
        !A为系数矩阵，BB（i）为b矩阵，Q为求得解，Npanel为系数矩阵的维数。
        KK=A1(:,:,1)-A1(:,:,2)-A1(:,:,3)+A1(:,:,4)
        !call LSACG(KK,BB(:,1),Q(:,1))
        call bgmres(KK,BB(:,1),Q(:,1),NPanel,1,M,ITR)
        KK=A1(:,:,1)+A1(:,:,2)-A1(:,:,3)-A1(:,:,4)
        call bgmres(KK,BB(:,2),Q(:,2),NPanel,1,M,ITR)
        KK=A1(:,:,1)+A1(:,:,2)+A1(:,:,3)+A1(:,:,4)
        call bgmres(KK,BB(:,3),Q(:,3),NPanel,1,M,ITR)
        KK=A1(:,:,1)+A1(:,:,2)-A1(:,:,3)-A1(:,:,4)
        call bgmres(KK,BB(:,4),Q(:,4),NPanel,1,M,ITR)
        KK=A1(:,:,1)-A1(:,:,2)-A1(:,:,3)+A1(:,:,4)
        call bgmres(KK,BB(:,5),Q(:,5),NPanel,1,M,ITR)
        KK=A1(:,:,1)-A1(:,:,2)+A1(:,:,3)-A1(:,:,4)
        call bgmres(KK,BB(:,6),Q(:,6),NPanel,1,M,ITR)

         !源强对称性输出：
 
         !!surge    + - - +
               SSTRENGTH(1,1,1,1:NPANEL)=REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,1,1:NPANEL)=AIMAG(Q(1:NPANEL,1))
       
               SSTRENGTH(1,1,2,1:NPANEL)=-REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,1))  
            
               SSTRENGTH(1,1,3,1:NPANEL)=-REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,3,1:NPANEL)=-AIMAG(Q(1:NPANEL,1)) 
         
               SSTRENGTH(1,1,4,1:NPANEL)=REAL(Q(1:NPANEL,1))
               SSTRENGTH(2,1,4,1:NPANEL)=AIMAG(Q(1:NPANEL,1))   
       
        !!sway    +  +  -  -
               SSTRENGTH(1,2,1,1:NPANEL)=REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,1,1:NPANEL)=AIMAG(Q(1:NPANEL,2))
       
               SSTRENGTH(1,2,2,1:NPANEL)=REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,2,1:NPANEL)=AIMAG(Q(1:NPANEL,2))  
            
               SSTRENGTH(1,2,3,1:NPANEL)=-REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,3,1:NPANEL)=-AIMAG(Q(1:NPANEL,2)) 
         
               SSTRENGTH(1,2,4,1:NPANEL)=-REAL(Q(1:NPANEL,2))
               SSTRENGTH(2,2,4,1:NPANEL)=-AIMAG(Q(1:NPANEL,2))   
       
       
         !!heave   +   +   +   +
               SSTRENGTH(1,3,1,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,1,1:NPANEL)=AIMAG(Q(1:NPANEL,3))
       
               SSTRENGTH(1,3,2,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,2,1:NPANEL)=AIMAG(Q(1:NPANEL,3))  
            
               SSTRENGTH(1,3,3,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,3,1:NPANEL)=AIMAG(Q(1:NPANEL,3)) 
         
               SSTRENGTH(1,3,4,1:NPANEL)=REAL(Q(1:NPANEL,3))
               SSTRENGTH(2,3,4,1:NPANEL)=AIMAG(Q(1:NPANEL,3))   
       
       
        !!roll   +   +   -   -
               SSTRENGTH(1,4,1,1:NPANEL)=REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,1,1:NPANEL)=AIMAG(Q(1:NPANEL,4))
       
               SSTRENGTH(1,4,2,1:NPANEL)=REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,2,1:NPANEL)=AIMAG(Q(1:NPANEL,4))  
            
               SSTRENGTH(1,4,3,1:NPANEL)=-REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,3,1:NPANEL)=-AIMAG(Q(1:NPANEL,4)) 
         
               SSTRENGTH(1,4,4,1:NPANEL)=-REAL(Q(1:NPANEL,4))
               SSTRENGTH(2,4,4,1:NPANEL)=-AIMAG(Q(1:NPANEL,4))   
       
       
        !!pitch   +   -   -   +
               SSTRENGTH(1,5,1,1:NPANEL)=REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,1,1:NPANEL)=AIMAG(Q(1:NPANEL,5))
       
               SSTRENGTH(1,5,2,1:NPANEL)=-REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,5))  
            
               SSTRENGTH(1,5,3,1:NPANEL)=-REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,3,1:NPANEL)=-AIMAG(Q(1:NPANEL,5)) 
         
               SSTRENGTH(1,5,4,1:NPANEL)=REAL(Q(1:NPANEL,5))
               SSTRENGTH(2,5,4,1:NPANEL)=AIMAG(Q(1:NPANEL,5))   
       
       
         !!yaw   +   -   -   +
               SSTRENGTH(1,6,1,1:NPANEL)=REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,1,1:NPANEL)=AIMAG(Q(1:NPANEL,6))
       
               SSTRENGTH(1,6,2,1:NPANEL)=-REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,2,1:NPANEL)=-AIMAG(Q(1:NPANEL,6))  
            
               SSTRENGTH(1,6,3,1:NPANEL)=REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,3,1:NPANEL)=AIMAG(Q(1:NPANEL,6)) 
         
               SSTRENGTH(1,6,4,1:NPANEL)=-REAL(Q(1:NPANEL,6))
               SSTRENGTH(2,6,4,1:NPANEL)=-AIMAG(Q(1:NPANEL,6))   
       
        Phi=0.0
        do L=1,6
            do jj=1,2
	            do i=1,NB
			        do K=1,4

				          do j=1,NPanel
				            Phi(1,jj,L,1,I)=Phi(1,jj,L,1,I)+G0Dc(i,j,jj,k)*SSTRENGTH(1,L,k,j)-G0Ds(i,j,jj,k)*SSTRENGTH(2,L,k,j)
				            Phi(2,jj,L,1,I)=Phi(2,jj,L,1,I)+G0Dc(i,j,jj,k)*SSTRENGTH(2,L,k,j)+G0Ds(i,j,jj,k)*SSTRENGTH(1,L,k,j)
				    

				    
				        enddo
			        enddo
		        enddo
	        enddo
        enddo

        Phi(:,1,1,2,:)=-Phi(:,1,1,1,:)
        Phi(:,2,1,2,:)=Phi(:,2,1,1,:)
        Phi(:,1,1,3,:)=-Phi(:,1,1,1,:)
        Phi(:,2,1,3,:)=Phi(:,2,1,1,:)
        Phi(:,:,1,4,:)=Phi(:,:,1,1,:)

        !sway
        Phi(:,1,2,2,:)=Phi(:,1,2,1,:)
        Phi(:,2,2,2,:)=-Phi(:,2,2,1,:)
        Phi(:,1,2,3,:)=-Phi(:,1,2,1,:)
        Phi(:,2,2,3,:)=Phi(:,2,2,1,:)
        Phi(:,:,2,4,:)=-Phi(:,:,2,1,:)
        !heave
        Phi(:,1,3,2,:)=Phi(:,1,3,1,:)
        Phi(:,1,3,3,:)=Phi(:,1,3,1,:)
        Phi(:,2,3,2,:)=-Phi(:,2,3,1,:)
        Phi(:,2,3,3,:)=-Phi(:,2,3,1,:)
        Phi(:,:,3,4,:)=Phi(:,:,3,1,:)
        !roll
        Phi(:,1,4,2,:)=Phi(:,1,4,1,:)
        Phi(:,1,4,3,:)=-Phi(:,1,4,1,:)
        Phi(:,2,4,2,:)=-Phi(:,2,4,1,:)
        Phi(:,2,4,3,:)=Phi(:,2,4,1,:)
        Phi(:,:,4,4,:)=-Phi(:,:,4,1,:)
        !pitch
        Phi(:,1,5,2,:)=-Phi(:,1,5,1,:)
        Phi(:,1,5,3,:)=-Phi(:,1,5,1,:)
        Phi(:,2,5,2,:)=Phi(:,2,5,1,:)
        Phi(:,2,5,3,:)=Phi(:,2,5,1,:)
        Phi(:,:,5,4,:)=Phi(:,:,5,1,:)
        !yaw
        Phi(:,1,6,2,:)=-Phi(:,1,6,1,:)
        Phi(:,1,6,3,:)=Phi(:,1,6,1,:)
        Phi(:,2,6,2,:)=Phi(:,2,6,1,:)
        Phi(:,2,6,3,:)=-Phi(:,2,6,1,:)
        Phi(:,:,6,4,:)=-Phi(:,:,6,1,:)
        deallocate(BB,Q)
        !allocate(RK(1:Npanel,1:4),N1(1:Npanel),N2(1:Npanel),N3(1:Npanel),N4(1:Npanel),PHI1(1:2,1:4,1:NB))
        allocate(RK(1:Npanel,1:4),N1(1:Npanel),N2(1:Npanel),N3(1:Npanel),N4(1:Npanel))
        N1=0;N2=0;N3=0;N4=0
        do i=1,Npanel

        c0=0.0;c1=0.0;c2=0.0;c3=0.0
        c4=0.0;c5=0.0;c6=0.0;c7=0.0
        c01=0.0;e1=0.0;e2=0.0;e3=0.0;e4=0.0
	        if(abs(omeI0)<1e-10)then
		        c0=0.0
	        else
		        c0=g0/omeI0
	        endif
	        xh=Pnu*H
	        if(xh>7.or.(trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		        c1=exp(Pnu*(xav(3,i)))
		        c01=c0*c1
	        c3=k0*xav(1,i)*cos(HEADI)   !k0*x*cos(beta)
	        c4=k0*xav(2,i)*sin(HEADI)    !k0*y*sin(beta)
	        c5=cmplx(k0*e(3,3,i),0)
	        c6=cmplx(0,k0*cos(Headi)*e(1,3,i))
	        c7=cmplx(0,k0*sin(Headi)*e(2,3,i))
	        e1=cmplx(0,cos(c3)*cos(c4)*c01)
	        e2=cmplx(0,sin(c3)*sin(c4)*c01)
	        e3=cmplx(-sin(c3)*cos(c4)*c01,0)
	        e4=cmplx(cos(c3)*sin(c4)*c01,0)
	        N1(i)=e1*c5+e3*c6-e4*c7
	        N2(i)=e2*c5+e4*c6-e3*c7
	        N3(i)=e3*c5+e1*c6-e2*c7
	        N4(i)=e4*c5+e2*c6-e1*c7	
		
		
	        else
            c1=cosh(k0*(xav(3,i)+H))/cosh(k0*H)
	        c2=sinh(k0*(xav(3,i)+H))/cosh(k0*H)
	        c01=c0*c1
	        c3=k0*xav(1,i)*cos(HEADI)   !k0*x*cos(beta)
	        c4=k0*xav(2,i)*sin(HEADI)    !k0*y*sin(beta)
	        c5=cmplx(k0*e(3,3,i),0)
	        c6=cmplx(0,k0*cos(Headi)*e(1,3,i))
	        c7=cmplx(0,k0*sin(Headi)*e(2,3,i))
	        e1=cmplx(0,cos(c3)*cos(c4)*c01)
	        e2=cmplx(0,sin(c3)*sin(c4)*c01)
	        e3=cmplx(-sin(c3)*cos(c4)*c01,0)
	        e4=cmplx(cos(c3)*sin(c4)*c01,0)
	        N1(i)=e1*c2*c5/c1+e3*c6-e4*c7
	        N2(i)=e2*c2*c5/c1+e4*c6-e3*c7
	        N3(i)=e3*c2*c5/c1+e1*c6-e2*c7
	        N4(i)=e4*c2*c5/c1+e2*c6-e1*c7
	        endif


        end do

        KK=A1(:,:,1)+A1(:,:,2)+A1(:,:,3)+A1(:,:,4)
        call bgmres(KK,-N1,RK(:,1),NPanel,1,M,ITR)
        KK=A1(:,:,1)-A1(:,:,2)+A1(:,:,3)-A1(:,:,4)
        call bgmres(KK,-N2,RK(:,2),NPanel,1,M,ITR)
        KK=A1(:,:,1)-A1(:,:,2)-A1(:,:,3)+A1(:,:,4)
        call bgmres(KK,-N3,RK(:,3),NPanel,1,M,ITR)
        KK=A1(:,:,1)+A1(:,:,2)-A1(:,:,3)-A1(:,:,4)
        call bgmres(KK,-N4,RK(:,4),NPanel,1,M,ITR)

        !输出绕射源强

               SSTRENGTH(1,7,1,1:NPANEL)=REAL(RK(1:NPANEL,1)+RK(1:NPANEL,2)+RK(1:NPANEL,3)+RK(1:NPANEL,4))
               SSTRENGTH(2,7,1,1:NPANEL)=AIMAG(RK(1:NPANEL,1)+RK(1:NPANEL,2)+RK(1:NPANEL,3)+RK(1:NPANEL,4))
       
               SSTRENGTH(1,7,2,1:NPANEL)=REAL(RK(1:NPANEL,1)-RK(1:NPANEL,2)-RK(1:NPANEL,3)+RK(1:NPANEL,4))
               SSTRENGTH(2,7,2,1:NPANEL)=AIMAG(RK(1:NPANEL,1)-RK(1:NPANEL,2)-RK(1:NPANEL,3)+RK(1:NPANEL,4))
       
               SSTRENGTH(1,7,3,1:NPANEL)=REAL(RK(1:NPANEL,1)+RK(1:NPANEL,2)-RK(1:NPANEL,3)-RK(1:NPANEL,4))
               SSTRENGTH(2,7,3,1:NPANEL)=AIMAG(RK(1:NPANEL,1)+RK(1:NPANEL,2)-RK(1:NPANEL,3)-RK(1:NPANEL,4))
       
               SSTRENGTH(1,7,4,1:NPANEL)=REAL(RK(1:NPANEL,1)-RK(1:NPANEL,2)+RK(1:NPANEL,3)-RK(1:NPANEL,4))
               SSTRENGTH(2,7,4,1:NPANEL)=AIMAG(RK(1:NPANEL,1)-RK(1:NPANEL,2)+RK(1:NPANEL,3)-RK(1:NPANEL,4))

        !输出绕射源强
			        do i=1,NB
			            do k=1,4
        select case(k)
			          case(1)
			             k2=2
			             k3=3
			             k4=4
			             case(2)
			             k2=1
			             k3=4
			             k4=3
			             case(3)
			             k2=4
			             k3=1
			             k4=2
			             case(4)
			             k2=3
			             k3=2
			             k4=1
			             end select
				        do j=1,NPanel
				        !速度势的实部和虚部
					        PHI(1,1,7,1,i)=PHI(1,1,7,1,i)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
					        PHI(2,1,7,1,i)=PHI(2,1,7,1,i)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
	
		                    Phi(1,1,7,2,I)=Phi(1,1,7,2,I)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k2,j),SSTRENGTH(2,7,k2,j)))
				            Phi(2,1,7,2,I)=Phi(2,1,7,2,I)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k2,j),SSTRENGTH(2,7,k2,j)))
				    
				            Phi(1,1,7,3,I)=Phi(1,1,7,3,I)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k3,j),SSTRENGTH(2,7,k3,j)))
				            Phi(2,1,7,3,I)=Phi(2,1,7,3,I)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k3,j),SSTRENGTH(2,7,k3,j)))
				    
				            Phi(1,1,7,4,I)=Phi(1,1,7,4,I)+real(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				            Phi(2,1,7,4,I)=Phi(2,1,7,4,I)+AIMAG(CMPLX(G0Dc(i,j,1,k),G0Ds(i,j,1,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				        !速度势的x方向导数实部和虚部
					        PHI(1,2,7,1,i)=PHI(1,2,7,1,i)+real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
					        PHI(2,2,7,1,i)=PHI(2,2,7,1,i)+AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k,j),SSTRENGTH(2,7,k,j)))
	
		                    Phi(1,2,7,2,I)=Phi(1,2,7,2,I)-real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k2,j),SSTRENGTH(2,7,k2,j)))
				            Phi(2,2,7,2,I)=Phi(2,2,7,2,I)-AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k2,j),SSTRENGTH(2,7,k2,j)))
				    
				            Phi(1,2,7,3,I)=Phi(1,2,7,3,I)-real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k3,j),SSTRENGTH(2,7,k3,j)))
				            Phi(2,2,7,3,I)=Phi(2,2,7,3,I)-AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k3,j),SSTRENGTH(2,7,k3,j)))
				    
				            Phi(1,2,7,4,I)=Phi(1,2,7,4,I)+real(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))
				            Phi(2,2,7,4,I)=Phi(2,2,7,4,I)+AIMAG(CMPLX(G0Dc(i,j,2,k),G0Ds(i,j,2,k))*cmplx(SSTRENGTH(1,7,k4,j),SSTRENGTH(2,7,k4,j)))				    
				    
				        enddo
				        end do
	        enddo
	end if
	
	
!BASS=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.BASS'
!OPEN (30021,FILE=BASS)
!do i=1,NB  !p
!if(xav(2,i)>=0)then
!do j=1,7
!
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),SSTRENGTH(1,j,1,I),SSTRENGTH(2,j,1,I)
!end do
!endif
!END DO
!do i=1,NB  !p
!if(xav(2,i)<=0)then
!do j=1,7
!
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),SSTRENGTH(1,j,1,I),SSTRENGTH(2,j,1,I)
!end do
!end if
!END DO

end subroutine solve_gmres
END MODULE SOLVEPHI_MOD
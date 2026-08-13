!----积分瞬时湿表面网格得到合外力
subroutine SurfaceIntegral1(Nwh,Nl,Node0,dtbphi,Force,InstUr )

  use ArrayOperations
!  use PanelGeometry,only:Nwh,Nl
  use verification
  !use ShipHullVar,only:NR
  use GauLegCoe
  use Constant,only: rou,NR
  Implicit None

  integer(4)::Nwh,Nl
  real(8),dimension(1:Nwh,1:Nl,1:3 )::Node0   !---重新划分的网格
  real*8 :: dtbphi(1:Nwh,1:Nl)
  real*8 :: Force(1:NR)
  !real(8),dimension(3)::Cog1
  real(8),dimension(Nwh,Nl,NR,3)::InstUr   !----瞬时节点各模态位移

  real(8),dimension(4)::pasiq
  real*8 :: dtphiq
  real(8),dimension(4,NR,3)::eur
  real(8),dimension(NR,3)::tempur

  real*8 :: vnaq(1:6)
  integer :: i,j,jj,j1,k

  real(8),dimension(1:3)::vel


     Force=0.0
  	Do i=1,Nwh-1
	Do j=1,Nl
        if (j<Nl) then
          x1=Node0(i,j,1);y1=Node0(i,j,2);z1=Node0(i,j,3)
          x2=Node0(i+1,j,1);y2=Node0(i+1,j,2);z2=Node0(i+1,j,3)
          x3=Node0(i+1,j+1,1);y3=Node0(i+1,j+1,2);z3=Node0(i+1,j+1,3)
          x4=Node0(i,j+1,1);y4=Node0(i,j+1,2);z4=Node0(i,j+1,3)

          pasiq(1)=dtbphi(i,j);      pasiq(2)=dtbphi(i+1,j);
          pasiq(3)=dtbphi(i+1,j+1);  pasiq(4)=dtbphi(i,j+1);
          eur(1,1:NR,1:3)=InstUr(i,j,1:NR,1:3);       eur(2,1:NR,1:3)=InstUr(i+1,j,1:NR,1:3);
          eur(3,1:NR,1:3)=InstUr(i+1,j+1,1:NR,1:3);   eur(4,1:NR,1:3)=InstUr(i,j+1,1:NR,1:3);
        else
          x1=Node0(i,j,1);y1=Node0(i,j,2);z1=Node0(i,j,3)
          x2=Node0(i+1,j,1);y2=Node0(i+1,j,2);z2=Node0(i+1,j,3)
          x3=Node0(i+1,1,1);y3=Node0(i+1,1,2);z3=Node0(i+1,1,3)
          x4=Node0(i,1,1);y4=Node0(i,1,2);z4=Node0(i,1,3)

          pasiq(1)=dtbphi(i,j);      pasiq(2)=dtbphi(i+1,j);
          pasiq(3)=dtbphi(i+1,1);    pasiq(4)=dtbphi(i,1);
          eur(1,1:NR,1:3)=InstUr(i,j,1:NR,1:3);     eur(2,1:NR,1:3)=InstUr(i+1,j,1:NR,1:3);
          eur(3,1:NR,1:3)=InstUr(i+1,1,1:NR,1:3);   eur(4,1:NR,1:3)=InstUr(i,1,1:NR,1:3);
        end if
        
        DxDxi=x2-x1;DyDxi=y2-y1;DzDxi=z2-z1
        DxDet=x4-x1;DyDet=y4-y1;DzDet=z4-z1
        hs1_x(1) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(1) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(1) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(1)=(hs1_x(1)**2+hs1_y(1)**2+hs1_z(1)**2)**0.5
        
        DxDxi=x3-x2;DyDxi=y3-y2;DzDxi=z3-z2
        DxDet=x1-x2;DyDet=y1-y2;DzDet=z1-z2
        hs1_x(2) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(2) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(2) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(2)=(hs1_x(2)**2+hs1_y(2)**2+hs1_z(2)**2)**0.5
          
        DxDxi=x4-x3;DyDxi=y4-y3;DzDxi=z4-z3
        DxDet=x2-x3;DyDet=y2-y3;DzDet=z2-z3
        hs1_x(3) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(3) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(3) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(3)=(hs1_x(3)**2+hs1_y(3)**2+hs1_z(3)**2)**0.5
        
        DxDxi=x1-x4;DyDxi=y1-y4;DzDxi=z1-z4
        DxDet=x3-x4;DyDet=y3-y4;DzDet=z3-z4
        hs1_x(4) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(4) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(4) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(4)=(hs1_x(4)**2+hs1_y(4)**2+hs1_z(4)**2)**0.5

    
	   Do jj=1,2
	   Do j1=1,2

		ph1=(1.+cor2(j1))*(1.+cor2(jj))/4.
		ph2=(1.-cor2(j1))*(1.+cor2(jj))/4.
		ph3=(1.-cor2(j1))*(1.-cor2(jj))/4.
		ph4=(1.+cor2(j1))*(1.-cor2(jj))/4.

          xq=ph1*x1+ph2*x2+ph3*x3+ph4*x4
	      yq=ph1*y1+ph2*y2+ph3*y3+ph4*y4
          zq=ph1*z1+ph2*z2+ph3*z3+ph4*z4

          hs_x=(ph1*hs1_x(1)+ph2*hs1_x(2)+ph3*hs1_x(3)+ph4*hs1_x(4))/4.        
          hs_y=(ph1*hs1_y(1)+ph2*hs1_y(2)+ph3*hs1_y(3)+ph4*hs1_y(4))/4.         
          hs_z=(ph1*hs1_z(1)+ph2*hs1_z(2)+ph3*hs1_z(3)+ph4*hs1_z(4))/4. 
          hs=(ph1*hs1(1)+ph2*hs1(2)+ph3*hs1(3)+ph4*hs1(4))/4.
          
          !-----压力积分节点
          dtphiq=ph1*pasiq(1)+ph2*pasiq(2)+ph3*pasiq(3)+ph4*pasiq(4)
          !----节点各模态位移积分点
          tempur(1:NR,1:3)=ph1*eur(1,1:NR,1:3)+ph2*eur(2,1:NR,1:3)+ph3*eur(3,1:NR,1:3)+ph4*eur(4,1:NR,1:3)


          !-----此处与FNL有出入,正负号问题(注意核查)//负号表示投影向量指向船壳内部
          vnaq(1)=-hs_x;
          vnaq(2)=-hs_y;
          vnaq(3)=-hs_z;
          vel(1)=xq; vel(2)=yq; vel(3)=zq;  !---关于重心
          vnaq(4:6)=R_Rcross( vel(1:3),vnaq(1:3) )
	    	    	    	
          Do k=1,NR
              Force(k)=Force(k)-ww(j1)*ww(jj)*( tempur(k,1)*vnaq(1)+tempur(k,2)*vnaq(2)+tempur(k,3)*vnaq(3) )*dtphiq
          Enddo
	   end do
	   end do

     end do
     end do
     Force=rou*Force

     return
end subroutine SurfaceIntegral1




!输入：smtf（光滑函数）；tmpx（第一个砰击典型节点平动坐标系下的瞬时位置）；time（time0时刻）
!输出：fetaZ（与time0对应的入射波面升高）；DerPhiI（与time0对应的波面的瞬时绝对速度）    
    

subroutine Slam_WaveElevation(smtf,tmpx,time,fetaZ,DerPhiI  )
  
 
  use Constant
  !20231113修改
  use IrreWaveVar
  
  implicit none

  real(8)::smtf    !----平滑函数
  real(8),dimension(3)::tmpx   !---计算点位置(平动坐标系)
  real(8)::time
  real(8)::tmpt                    

  real(8)::fetaZ                  !---点位置(x,y)对应的波面起伏Z
  real(8)::fetaVx,fetaVy,fetaVz   !---与(x,y)对应的波面上的点的运动速度,本质上是入射势的梯度

  real(8),dimension(3)::DerPhiI
  integer(4)::i,j,k,jj

  !-----默认z=0.0
  !20231113修改
    if( IrreCtrl==0 )then
        tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*time   
        fetaZ=amp*cos(tmpt )
        fetaVx=amp*ome*cos(head)* cos( tmpt ) 
        fetaVy=amp*ome*sin(head)* cos( tmpt )
        fetaVz=amp*ome*sin(tmpt )                      !---对应点入射波垂向速度
        
    else
        fetaZ=0.0;   fetaVx=0.0;
        fetaVy=0.0;  fetaVz=0.0;
        do jj=1,IrreNum
            tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*time + Irrepha(jj)
            fetaZ=fetaZ+Irreamp(jj)*cos( tmpt )
            fetaVx=fetaVx+Irreamp(jj)*Irreome(jj)*cos(head)* cos( tmpt )
            fetaVy=fetaVy+Irreamp(jj)*Irreome(jj)*sin(head)* cos( tmpt )
            fetaVz=fetaVz+Irreamp(jj)*Irreome(jj)*sin(tmpt )
        end do
   
    end if
      

   fetaZ=fetaZ*smtf

   DerPhiI(1)=fetaVx*smtf;
   DerPhiI(2)=fetaVy*smtf;
   DerPhiI(3)=fetaVz*smtf;

  return
end subroutine Slam_WaveElevation
    
    
    
    
    
    
subroutine Slam_WaveAngle(smtf,tmpx,time,DerPhiI  )
  
  
  use Constant
  !20231113修改
  use IrreWaveVar
  implicit none

  real(8)::smtf    !----平滑函数
  real(8),dimension(3)::tmpx   !---计算点位置(平动坐标系)
  real(8)::time
  real(8)::tmpt                    

  real(8)::fetaZ                  !---点位置(x,y)对应的波面起伏Z
  real(8)::fetadx,fetady          !---与(x,y)对应波面的导数dz/dx;dz/dy

  real(8),dimension(2)::DerPhiI
  integer(4)::i,j,k,jj

  !-----默认z=0.0
  !20231113修改
  
    if( IrreCtrl==0 )then
        tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*time   
        fetaZ=amp*cos(tmpt )
        fetadx=-amp*wavek*cos(head)*sin( tmpt )
        
        fetady=-amp*wavek*sin(head)*sin( tmpt )
        
        !!fetaVx=amp*ome*cos(head)* cos( tmpt ) 
        !!fetaVy=amp*ome*sin(head)* cos( tmpt )
        !!fetaVz=amp*ome*sin(tmpt )                      !---对应点入射波垂向速度
   
    else
        fetaZ=0.0;   
        fetadx=0.0; fetady=0.0;
        
        do jj=1,IrreNum
            
            tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*time + Irrepha(jj)
            fetaZ=fetaZ+Irreamp(jj)*cos( tmpt )
            
            fetadx=fetadx-Irreamp(jj)*Irrek(jj)*cos(head)*sin( tmpt )
            fetady=fetady-Irreamp(jj)*Irrek(jj)*sin(head)*sin( tmpt )
            
            
        end do
  
    end if
   
    
   

   fetaZ=fetaZ*smtf
   DerPhiI(1)=fetadx*smtf;
   DerPhiI(2)=fetady*smtf;
   
   
   DerPhiI(1)=atan( DerPhiI(1) );
   
   DerPhiI(2)=atan( DerPhiI(2) );
   

   !!DerPhiI(1)=fetaVx*smtf;
   !!DerPhiI(2)=fetaVy*smtf;
   !!DerPhiI(3)=fetaVz*smtf;

  return
end subroutine Slam_WaveAngle  
    
    
    
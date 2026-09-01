Module IrreWaveVar
!=========================================================================================   
!Function... 
!   1.Irregular wave generation
!    
!CopyRight... 
!   1.Origional code by Li Zhi-fu,2015.12.02
! 
!=========================================================================================    
    implicit none 
    
    real*8 :: IrreEle
    !integer*4:: IrreCtrl
    integer*4:: IrreNum        !×Ó²¨¸öÊý
    real*8,allocatable,dimension(:):: IrreAmp
    real*8,allocatable,dimension(:):: IrreOme,IrreOmee
    real*8,allocatable,dimension(:):: Irrek
    real*8,allocatable,dimension(:):: Irrepha
    
Endmodule IrreWaveVar
subroutine PreIrreWave(hs,tp)
!=========================================================================================   
!Function... 
!   1.Irregular wave generation
!    
!CopyRight... 
!   1.Origional code by Wang Siyu,2018
!=========================================================================================
    
    use IrreWaveVar
    
    implicit none
    
    real(8):: hs,t1,tp
    real(8) :: wmin,wmax,deltw
    real(8) :: xx,aa,bb
    real(8),allocatable,dimension(:):: ww,sw
    integer :: i


    t1=0.772*tp

    wmin=0.05; wmax=3.5*2*3.1415926/1.294/t1  

    aa=173*hs**2/t1**4; bb=691/t1**4

    allocate(ww(1:IrreNum),sw(1:IrreNum))

    deltw=(wmax-wmin)/(IrreNum-1)

    do i=1,IrreNum
        ww(i)=wmin+(i-1)*deltw
    enddo
    Irreome(1)=wmin;   Irreome(IrreNum)=wmax

    call random_seed () 

    do i=1,IrreNum-1
        call random_number (xx) 
        Irreome(i+1)=xx*(ww(i+1)-ww(i))+ww(i)
    enddo


    do i=1,IrreNum
        sw(i)=aa/Irreome(i)**5*exp(-bb/Irreome(i)**4)
    enddo

    call random_seed ()
    do i=1,IrreNum 
        if(i/=1)then
            IrreAmp(i)=((Irreome(i)-Irreome(i-1))*2.0*sw(i))**0.5
        else
            IrreAmp(i)=((Irreome(i+1)-Irreome(i))*2.0*sw(i))**0.5
        endif
        call random_number (xx) 
        IrrePha(i)=2*3.1415926*xx
    enddo

end subroutine PreIrreWave
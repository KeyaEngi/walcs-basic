!x+iy-->(amp,phase)»¡¶È
subroutine comp_to_ampha(x,y,amp,phase)
implicit none
real(8)::x,y
real(8),intent(out)::amp,phase
real(8)::PI
real(8)::esp
esp=1.0e-8
! PI=4.0*atan(1.0)
PI=3.141592654

if(abs(x)<=esp)	x=0.0
if(abs(y)<=esp)	y=0.0
amp=sqrt(x*x+y*y)
if(abs(x)<1e-6)then
	if(abs(y)<1e-6)then
		phase=0.0
	elseif(y>0.0)then
		!phase=PI/2
		phase=90.0
	else
		!phase=-PI/2
		phase=270.0
	endif
elseif(x>0.0)then
	if(y>0.0)then
		phase=atan(y/x)/PI*180.0 
	else
		phase=atan(y/x)/PI*180.0+360.0
	endif
elseif(x<0.0)then
	if(y>0.0)then
		phase=(PI+atan(y/x))/PI*180.0
	else
		phase=(-PI+atan(y/x)+2*PI)/PI*180.0
	endif
endif
end subroutine comp_to_ampha
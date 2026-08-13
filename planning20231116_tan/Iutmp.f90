subroutine iutmp(iu)
!=========================================================================================
!Function... 
!   1.用于略去输入数据文件中的#注释内容:#必须在第一列，注释内容以字符"#"打头,该行在#后的
!     内容都被认为是注释内容
!    
!CopyRight... 
!   1.Origional code by Shan Penghao
!-----------------------------------------------------------------------------------------
!Notes...
!   1.默认读取数据的列数不能大于200，否则就要加大字符串C的声明长度
! 
!-----------------------------------------------------------------------------------------
!variable ...
!   1.iu : 输入文件号
!   2.iu_tmp : 临时文件号
!   3.rf_tmp : 临时文件名    
!
!=========================================================================================
    
    implicit none 
    
    integer::iu               
    integer::iu_tmp           
    logical::ex,exf
    character(len=14)::rf_tmp 
    character(len=2)::no
    character(len=2000)::c
    integer::pp
    
!   以下为判断文件号iu_tmp是否已经被使用的代码
iu_tmp=1
do
	inquire(iu_tmp,opened=ex)    !判断文件号iu_tmp是否已经被使用
	call qinttostr(iu_tmp,no,2)  !此时no为长度为2的字符串
    rf_tmp='SWIFT'//NO//'.tmp'
	inquire(file=rf_tmp,exist=exf)
	if(ex==.false..and.exf==.false.)exit
	iu_tmp=iu_tmp+1
end do

open(iu_tmp,file=rf_tmp)        !rf_tmp用于临时保存不含注释的文件内容
do while( .true. )
	    read(iu,'(a)',iostat=pp)c   
	    if  (pp/=0) exit    !判断文件是否正常读取内容，否，则退出
	    if (c /=" ") then    !如果读到的不是空行
               if(c(1:1)=='#')then      
                    cycle
               else
                   write(iu_tmp,'(a)')c
               end if
	    else                                        !空行
    	
	    end if
end do
close(iu);close(iu_tmp)

open(iu,file=rf_tmp,status='old',dispose='delete')

end subroutine iutmp

      
subroutine qinttostr(m,str,n)
!=========================================================================================
!Function... 
!   1.利用文件号创建临时的文件名:将一个整型数据转换为包含n个字符的字符串(取末尾的几个数字)，
!     1--01；2--02；...8--08
!    
!CopyRight... 
!   1.Origional code by Li Zhi-fu,2015.09.06
!-----------------------------------------------------------------------------------------
!Notes...
!   1.n不能超过10(若需要更长的字符串可通过调整程序中数组inc的大小实现)
! 
!-----------------------------------------------------------------------------------------
!variable ...
!   1.iu : 输入文件号
!   2.iu_tmp : 临时文件号
!   3.rf_tmp : 临时文件名    
!
!=========================================================================================

    implicit none
    
    integer::m,n
    character(len=n)::str
    integer::n1,i
    integer,dimension(1:10)::inc

n1=m
do i=1,n
	inc(i)=n1-int(n1/10)*10
	n1=int(n1/10)
end do

str=char(48+inc(n))

do i=2,n
	str=trim(str)//char(48+inc(n-i+1))
end do

end subroutine qinttostr


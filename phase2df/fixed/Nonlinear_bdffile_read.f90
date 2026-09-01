
!-----非线性计算的第二种方式
!-----读取bdf文件

subroutine Nonlinear_bdffile_read(tempNonbdf)
  
  use Constant,only:InAccess,OutAccess
  use PanelGeometry,only:bdfnum_ele,bdfnum_node,bdfele_kind,bdfele_node,bdfnode
  implicit none

  character(len=100)::tempNonbdf

    character(len=500)::char1,char2,char3
    character(len=8)::temchar,temchar1,temchar2
    character(len=9)::temchar3

    character(len=16)::tem2char,tem2char1,tem2char2
    character(len=17)::tem2char3

    integer(8)::ielem,inode    !----壳单元个数，总的节点数
    integer(8)::telement(1:5000000),tnode(1:5000000,1:4) !---壳单元，节点编号
    integer(8)::tv_ele(1:5000000)   !---单元类型记录
    integer(8)::fx(1:5000000)   !---壳单元属性编号

    integer(8),allocatable,dimension(:)::bdfnode_ID,bdfele_ID


    integer(8)::length1,Nlap

    integer(8)::i,j,k,ii,jj,kk



!*******************************************
!*******************************************

  open(unit=12,file=trim(adjustl(InAccess))//'\'//trim(adjustl(tempNonbdf)) )

    !----中间文件
    open(unit=21,file='tem_Cquad4.txt')   !---四边形单元存储
    open(unit=22,file='tem_Ctria3.txt')   !---三角形单元存储
    open(unit=23,file='tem_Node.txt')     !---节点存储


    ielem=0;
    inode=0;


    telement=0;
    tnode=0;

    tv_ele=0;


    read(12,'(A)')	char2
    read(12,'(A)')	char3
    do
        char1=char2    !----当前行
        char2=char3    !----前进一行
        if(char2(1:7)=='ENDDATA')	then
		    exit     !读到文件末尾，跳出循环
	    endif
        read(12,'(A)')	char3    !----再下一行


        if( char1(1:8)=='CQUAD4  ') then     !---四边形单元
            ielem=ielem+1;
            tv_ele(ielem)=4;

            read(char1,"(A8,i8,i8,i8,i8,i8,i8 )" ) temchar,telement(ielem),fx(ielem),tnode(ielem,1),tnode(ielem,2),tnode(ielem,3),tnode(ielem,4)
            write(21,"(6(i8,1x))") telement(ielem),fx(ielem),(tnode(ielem,j),j=1,4)

        elseif(char1(1:8)=='CQUAD4* '  ) then   !----四边形单元(16位存储)
            ielem=ielem+1;
            tv_ele(ielem)=4;
            
            read(char1,"(A8,i16,i16,i16,i16 )" )   temchar,telement(ielem),fx(ielem),tnode(ielem,1),tnode(ielem,2)  
            read(char2,"(A8,i16,i16)") temchar,tnode(ielem,3),tnode(ielem,4)

            write(21,"(6(i16,1x))") telement(ielem),fx(ielem),(tnode(ielem,j),j=1,4)
        elseif(char1(1:8)=='CTRIA3  ' ) then   !---三角形单元
            ielem=ielem+1;
            tv_ele( ielem)=3;

            read(char1,"(A8,i8,i8,i8,i8,i8 )" ) temchar,telement(ielem),fx(ielem),tnode(ielem,1),tnode(ielem,2),tnode(ielem,3)
            write(22,"(5(i8,1x))") telement(ielem),fx(ielem),(tnode(ielem,j),j=1,3)
        elseif(char1(1:8)=='CTRIA3* '  ) then   !--三角形单元(16位存储)
            ielem=ielem+1;
            tv_ele( ielem)=3;

            read(char1,"(A8,i16,i16,i16,i16 )" )   temchar,telement(ielem),fx(ielem),tnode(ielem,1),tnode(ielem,2) 
            read(char2,"(A8,i16)") temchar,tnode(ielem,3)

            write(22,"(5(i16,1x))") telement(ielem),fx(ielem),(tnode(ielem,j),j=1,3)


        elseif(char1(1:8)=='GRID    ' ) then   !---处理网格节点信息
            inode=inode+1;

            length1=len(trim(adjustl(char1)))
		    Nlap=int(length1/8)+1

            jj=0
            do i=2,Nlap
                 if( char1((i-1)*8+1:i*8 )/='        ' ) then
                     do j=2,8
                         if( char1((i-1)*8+j:(i-1)*8+j )=='-'.or.char1((i-1)*8+j:(i-1)*8+j )=='+' ) exit
                     end do
                     if(j<=8) then
                         temchar1(1:j-1)=char1( (i-1)*8+1:(i-1)*8+j-1 )
                         temchar2(j:8)=char1( (i-1)*8+j:(i-1)*8+8 )
                         temchar3=trim(adjustl(temchar1(1:j-1)))//'E'//trim(adjustl(temchar2(j:8) ))

                         write(23,"(a10,1x,\)",advance='NO') temchar3
                     else
                         temchar1(1:8)=char1( (i-1)*8+1:(i-1)*8+8 )
                         write(23,"(a10,1x,\)",advance='NO') temchar1
                     end if
                 end if

                 jj=jj+1
                 if(jj==6) exit
            end do
            write(23,"(/)")

        end if

    end do

    close(21,status='delete')
    close(22,status='delete')
    close(23)

    write(*,*) 'bdf文件读取完毕'      
    close(12)


    !--------处理单元信息
    bdfnum_ele=ielem      !----壳单元个数
    bdfnum_node=inode     !----总的节点个数

    allocate(bdfele_ID(1:bdfnum_ele),bdfnode_ID(1:bdfnum_node),bdfele_kind(1:bdfnum_ele) )
    allocate(bdfele_node(1:bdfnum_ele,1:4),bdfnode(1:bdfnum_node,1:3) )

    do i=1,bdfnum_ele
        bdfele_ID(i)=telement(i);   !---单元编号
        bdfele_node(i,1:4)=tnode(i,1:4)   !---单元对应节点编号
        bdfele_kind(i)=tv_ele(i)    !---单元种类编号
    end do

    open(unit=23,file='tem_Node.txt')    !----节点信息
    do i=1,bdfnum_node
        read(23,*) bdfnode_ID(i),(bdfnode(i,j),j=1,3)
    end do
    close(23,status='delete')

!----对单元对应节点ID进一步处理
!----节点node按读入顺序排列，顺序号即为ID号，需要修改相应单元对应节点ID号
    write(*,*) '处理壳单元节点对应信息'
    i=1;
    do 
        
        if(i>bdfnum_ele) exit

        kk=0;
        do j=1,bdfele_kind(i)
            k=bdfele_node(i,j)

            do ii=1,bdfnum_node
                if(k==bdfnode_ID(ii) ) then
                    bdfele_node(i,j)=ii
                    kk=kk+1;
                    exit
                end if
            end do
        end do

        if( kk/=bdfele_kind(i) ) then
             do j=i+1,bdfnum_ele
                 bdfele_ID(j-1)=bdfele_ID(j);
                 bdfele_node(j-1,:)=bdfele_node(j,:);
                 bdfele_kind(j-1)=bdfele_kind(j);
             end do

             bdfnum_ele=bdfnum_ele-1
             cycle
        end if

        i=i+1;
    end do

!*******此时有用的数据包括
!*******  bdfnum_ele  壳单元个数
!*******  bdfele_kind 壳单元类型(三角形或四边形)
!*******  bdfele_node 壳单元对应节点编号
!*******  bdfnum_node 节点个数
!*******  bdfnode(i,j) 节点坐标(当前还是在patran内部坐标系下)






  return
end subroutine Nonlinear_bdffile_read





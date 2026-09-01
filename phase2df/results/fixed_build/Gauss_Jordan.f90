

subroutine Gauss_Jordan(row,col,matrixIni,matrix1 )

  implicit none

        integer(4)::row    !-----行数
        integer(4)::col    !-----列数       行列数应该相等
        real(8),dimension(row,col)::matrixIni
        
        real(8),dimension(row,col)::matrix
        real(8),dimension(row,col)::matrix1
        real(8)::max
        integer(4),allocatable,dimension(:)::is,js     !-----记录主元位置信息（is为行数，js为列数）

        integer(4)::i,j,k


        allocate(is(row),js(col))

        
        matrix=matrixIni
        !----------全选主元

        do k=1,row
            max=matrix(k,k)
            is(k)=k; js(k)=k

            do i=k,row     !-----确定主元位置
                do j=k,col
                    if(abs(max)<abs(matrix(i,j)) ) then
                        max=matrix(i,j)
                        is(k)=i; js(k)=j
                    end if
                end do
            end do

!!            write(*,"(e15.6)") max
            if(abs(max)==0) then
                write(*,"(a30,i12)") '无逆矩阵',k
                exit
            end if

            !-----进行行列交换
            if(is(k)/=k) then
                call matrix_swap_row(row,col,matrix,is(k),k)   !----行交换
            end if
            if(js(k)/=k) then
                call matrix_swap_col(row,col,matrix,js(k),k)   !----列交换
            end if

            matrix(k,k)=1.0/matrix(k,k)   !---处理主元

            do j=1,col
                if(j/=k) then
                    matrix(k,j)=matrix(k,j)*matrix(k,k)
                end if
            end do
            
            do i=1,row
                if(i/=k) then
                    do j=1,col
                        if(j/=k) then
                            matrix(i,j)=matrix(i,j)-matrix(i,k)*matrix(k,j)
                        end if
                    end do
                end if
            end do

            do i=1,row
                if(i/=k) then
                    matrix(i,k)=-matrix(i,k)*matrix(k,k)
                end if
            end do

        end do

        !-------恢复
        do k=row,1,-1
            if(js(k)/=k) then
                call matrix_swap_row(row,col,matrix,js(k),k)
            end if
            if(is(k)/=k) then
                call matrix_swap_col(row,col,matrix,is(k),k)
            end if
        end do

        matrix1=matrix


        return
end subroutine Gauss_Jordan


    subroutine matrix_swap_row(row,col,matrix,k1,k2)   !-----对于矩阵的k1和k2行进行互换
        implicit none
        integer(4)::row,col
        real(8),dimension(row,col)::matrix
        integer(4)::k1,k2

        integer(4)::i
        real(8)::s

        do i=1,col
            s=matrix(k1,i)
            matrix(k1,i)=matrix(k2,i)
            matrix(k2,i)=s
        end do

        return
    end subroutine matrix_swap_row

    subroutine matrix_swap_col(row,col,matrix,k1,k2)   !-----对于矩阵的k1和k2列进行互换
        implicit none
        integer(4)::row,col
        real(8),dimension(row,col)::matrix
        integer(4)::k1,k2

        integer(4)::i
        real(8)::s

        do i=1,row
            s=matrix(i,k1)
            matrix(i,k1)=matrix(i,k2)
            matrix(i,k2)=s
        end do

        return
    end subroutine matrix_swap_col




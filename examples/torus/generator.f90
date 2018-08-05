
program generate_torus
    implicit none
    
    integer, parameter :: num_rows_points = 21
    integer, parameter :: num_cols_points = 21
    integer, parameter :: num_rows_blocks = 1
    integer, parameter :: num_cols_blocks = 1
    real*8, parameter :: pi = 4.d0*datan(1.d0)
    real*8, parameter :: a = 2.0
    real*8, parameter :: b = 1.0
    
    integer :: i, j, k, l
    real*8 :: u, v
    real*8, dimension(num_rows_points, num_cols_points) :: x, y, z
    integer, dimension(num_rows_blocks, num_cols_blocks) :: num_rows_per_block, num_cols_per_block
    integer, dimension(num_rows_blocks, num_cols_blocks) :: start_u, end_u, start_v, end_v
    
    
    do i = 1, num_rows_points
        do j = 1, num_cols_points
            u = real(j-1)/(num_cols_points-1)
            v = real(i-1)/(num_rows_points-1)
            x(i, j) = (a + b * cos(v)) * cos(u) 
            y(i, j) = (a + b * cos(v)) * sin(u) 
            z(i, j) = b * sin(v)
        end do
    end do
    
    
    do i = 1, num_rows_blocks
        do j = 1, num_cols_blocks
            
            if (i .eq. 1) then
                start_u(i,j) = 1
            else 
                start_u(i,j) = (i-1)*(real(num_rows_points, 8)/(num_rows_blocks))
            end if
            end_u(i,j) = i*(real(num_rows_points, 8)/(num_rows_blocks))
            
            if (j .eq. 1) then
                start_v(i,j) = 1
            else 
                start_v(i,j) = (j-1)*(real(num_cols_points, 8)/(num_cols_blocks))
            end if
            end_v(i,j) = j*(real(num_cols_points, 8)/(num_cols_blocks))
        
            
            print *, start_u(i,j), end_u(i,j), start_v(i,j), end_v(i,j)
        end do
    end do
    
    
    
    
    open(2,file='mesh.xyz',form='unformatted')
    write (2) num_rows_blocks, num_cols_blocks
    
        
    do i = 1, num_rows_blocks
        do j = 1, num_cols_blocks

            write (2) end_u(i,j)-start_u(i,j)+1, end_v(i,j)-start_v(i,j)+1
            
        end do
    end do
    
    
    
    do i = 1, num_rows_blocks
        do j = 1, num_cols_blocks
            
            write (2) ((x(k, l), l = start_v(i,j), end_v(i,j)), k = start_u(i,j), end_u(i,j)), &
                ((y(k, l), l = start_v(i,j), end_v(i,j)), k = start_u(i,j), end_u(i,j)), &
                ((z(k, l), l = start_v(i,j), end_v(i,j)), k = start_u(i,j), end_u(i,j))

        end do
    end do
    
    
    

end program generate_torus



      
     

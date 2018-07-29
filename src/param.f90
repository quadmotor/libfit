
! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains procedures related to parameterization and 
!    the creation of parameter distributions.

module param_m

    use grid_m
    implicit none


contains

    !******************************************************************
    !***** Set Parameter Values -- Uniform
    !******************************************************************
    subroutine param_uniform(num_rows, num_cols, u, v)
        integer, intent(in) :: num_rows, num_cols
        real*8, dimension(num_rows*num_cols), intent(out) :: u, v
        
        integer :: i, j, ip
        
        do i = 1, num_rows
            do j = 1, num_cols
                ip = (i-1)*num_cols + j
                u(ip) = real(j-1, 8) / (num_cols-1)
                v(ip) = real(i-1, 8) / (num_rows-1)
            end do
        end do
        
    end subroutine param_uniform
    
    !******************************************************************
    !***** Set Parameter Values -- Arc Length
    !******************************************************************
    subroutine param_arclength(num_rows, num_cols, u, v, x, y, z)
        integer, intent(in) :: num_rows, num_cols
        real*8, dimension(num_rows*num_cols), intent(inout) :: u, v, x, y, z
        
        
        integer :: ip, i, j
        real*8 :: curvelen
        
        do i = 1, num_rows
        
            curvelen = 0.d0
            u((i-1)*num_cols + 1) = 0.d0
            do j = 2, num_cols
                ip = (i-1)*num_cols + j
                curvelen = curvelen + sqrt( &
                    (x(ip) - x(ip-1))**2 + &
                    (y(ip) - y(ip-1))**2 + &
                    (z(ip) - z(ip-1))**2)
                    
                u(ip) = curvelen
            enddo
            
            do j = 1, num_cols
                ip = (i-1)*num_cols + j
                u(ip) = u(ip) / curvelen
            end do
            
        enddo
        
        do j = 1, num_cols
            
            curvelen = 0.d0
            v(j) = 0.d0
            do i = 2, num_rows
                ip = (i-1)*num_cols+j
                curvelen = curvelen + sqrt( &
                    (x(ip) - x(ip - num_cols))**2 + &
                    (y(ip) - y(ip - num_cols))**2 + &
                    (z(ip) - z(ip - num_cols))**2)

                v(ip) = curvelen
            enddo
            do i = 1, num_rows
                ip = (i-1)*num_cols + j
                v(ip) = v(ip) / curvelen
            end do
        enddo
        
    end subroutine param_arclength
    
    !******************************************************************
    !***** Normalize Parameter Values to [0, 1]
    !******************************************************************
    subroutine normalize_param(num_rows, num_cols, u, v)
        implicit none
        integer, intent(in) :: num_rows, num_cols
        real*8, dimension(num_rows*num_cols), intent(inout) :: u, v
        
        
        integer :: i, j, ip
        
        do i = 1, num_rows
            do j = 1, num_cols
                ip = (i-1)*num_cols+j
                    
                   u(ip) = u(ip) - u(ip - j + 1)
                u(ip) = u(ip) / u(i*num_cols)
                
                v(ip) = v(ip) - v(j)
                v(ip) = v(ip) / v((num_rows-1)*num_cols+j)

            end do
        end do
    
    end subroutine normalize_param
    
    !******************************************************************
    !***** Parameterize Grid
    !******************************************************************
    subroutine param_grid(grid) 
        type(grid_t), intent(inout) :: grid
        
        integer :: i, j
        
        do i = 1, grid%num_rows
            do j = 1, grid%num_cols
            
                call param_arclength(grid%blocks(i, j)%num_rows, grid%blocks(i,j)%num_cols, &
                    grid%blocks(i,j)%u, grid%blocks(i,j)%v, grid%blocks(i,j)%x, grid%blocks(i,j)%y, &
                    grid%blocks(i,j)%z)
            
                call normalize_param(grid%blocks(i,j)%num_rows, grid%blocks(i,j)%num_cols, &
                    grid%blocks(i,j)%u, grid%blocks(i,j)%v)
            end do
        end do
        
        
    
    end subroutine param_grid



end module param_m

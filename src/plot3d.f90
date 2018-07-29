! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains a data type representing Plot3D files and
!    procedures read/write from that file.

module plot3d_m
    use grid_m
    implicit none

    type, public :: plot3d_t
        character(len=:), allocatable :: file_name
        type(grid_t), pointer :: grid
            
    contains
        procedure, public :: init => plot3d_init
        procedure, public :: destroy => plot3d_destroy 
        procedure, public :: read_xyz => plot3d_read_xyz
        procedure, public :: write_xyz => plot3d_write_xyz
    
    end type plot3d_t


contains

    !******************************************************************
    !***** Initialize Plot3D Type
    !******************************************************************
    subroutine plot3d_init(self, grid, file_name)
        class(plot3d_t), intent(inout) :: self
        type(grid_t), target, intent(inout) :: grid
        character(len=*), intent(in) :: file_name
        
        call self%destroy()
        
        self%grid => grid
        self%file_name = file_name
        
    
    end subroutine plot3d_init
    
    !******************************************************************
    !***** Destroy Plot3D Type
    !******************************************************************
    subroutine plot3d_destroy(self)
        class(plot3d_t), intent(inout) :: self
        
        if (allocated(self%file_name)) deallocate(self%file_name)
        
        nullify(self%grid)
    
    end subroutine plot3d_destroy


    !******************************************************************
    !***** Read Grid
    !******************************************************************
    subroutine plot3d_read_xyz(self)
        class(plot3d_t), intent(inout) :: self
        
        integer :: i, j, k, ip, num_rows, num_cols
        
        open(2, file=self%file_name, form='unformatted')

        ! Read number of rows of blocks and number of cols of blocks
        read (2) num_rows, num_cols
        call self%grid%init(num_rows, num_cols)
        
        ! Read each block's number of rows of points and number of cols of points
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols
                read (2) num_rows, num_cols
                call self%grid%blocks(i,j)%init(num_rows, num_cols)
                                
            end do
        end do
        
        ! Read each block's real-space data
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols
                read (2) (self%grid%blocks(i,j)%x(k), k = 1, self%grid%blocks(i,j)%num_points), &
                    (self%grid%blocks(i,j)%y(k), k = 1, self%grid%blocks(i,j)%num_points), &
                    (self%grid%blocks(i,j)%z(k), k = 1, self%grid%blocks(i,j)%num_points)
             
            end do
        end do
            
    end subroutine plot3d_read_xyz
    
    
    !******************************************************************
    !***** Write Grid
    !******************************************************************
    subroutine plot3d_write_xyz(self)
        class(plot3d_t), intent(inout) :: self
        
        integer :: i, j, k, num_rows, num_cols, num_blocks
        
        open(2, file='test.xyz', form='unformatted')

        ! Write number of rows of blocks and number of cols of blocks
        write (2) self%grid%num_rows, self%grid%num_cols
    
    
        ! Write each block's number of rows of points and number of cols of points
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols
                write (2) self%grid%blocks(i,j)%num_rows, self%grid%blocks(i,j)%num_cols
            end do
        end do
        
        ! Write each block's real-space data
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols
                write (2) self%grid%blocks(i,j)%x, self%grid%blocks(i,j)%y, self%grid%blocks(i,j)%z
            end do
        end do
            
    end subroutine plot3d_write_xyz

    

end module plot3d_m










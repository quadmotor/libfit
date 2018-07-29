! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains data types representng multiblock 
!    structured meshes.


module grid_m
    use block_m
    implicit none

    type, public :: grid_t
        ! Structured grids blocks
        type(block_t), dimension(:, :), allocatable :: blocks
        ! Number of rows of blocks, number of columns of blocks, and number of blocks
        integer :: num_rows, num_cols, num_blocks
    contains
        procedure, public :: init => grid_init
        procedure, public :: destroy => grid_destroy
    end type grid_t
            

contains

    !******************************************************************
    !***** Initialize Grid Type
    !******************************************************************
    subroutine grid_init(self, num_rows, num_cols)
        class(grid_t), intent(inout) :: self
        integer, intent(in) :: num_rows, num_cols
        
        self%num_rows = num_rows
        self%num_cols = num_cols
        self%num_blocks = self%num_rows * self%num_cols
        
        allocate(self%blocks(num_rows, num_cols))
    
    end subroutine grid_init
    
    !******************************************************************
    !***** Destroy Grid Type
    !******************************************************************
    subroutine grid_destroy(self)
        class(grid_t), intent(inout) :: self
        
        integer :: i, j
        
        if (allocated(self%blocks)) then
            
            do i = 1, self%num_rows
                do j = 1, self%num_cols
                    call self%blocks(i,j)%destroy()
                end do
            end do
        end if    
        
        if (allocated(self%blocks)) deallocate(self%blocks)
    
    end subroutine grid_destroy



end module grid_m

! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains a data type representing a numpy file
!    and procedure to read/write from that file. "Numpy files" are used for 
!    matplotlib visualization purposes as the python Numpy module can easily 
!    read matrices from file in this form. X, Y, and Z data are stored in separate
!     files.


module numpy_m
    use grid_m
    use bezier_m
    implicit none

    type, public :: numpy_t
        character(len=:), allocatable :: file_name_x, file_name_y, file_name_z
    
    contains
        procedure, public :: init => numpy_init
        procedure, public :: destroy => numpy_destroy
        procedure, public :: write_xyz => numpy_write_xyz
        procedure, public :: write_cp => numpy_write_cp
    
    end type numpy_t

contains
    
    !******************************************************************
    !***** Initialize Numpy Type
    !******************************************************************
    subroutine numpy_init(self, file_name_x, file_name_y, file_name_z)
        class(numpy_t), intent(inout) :: self
        character(len=*), intent(in) :: file_name_x, file_name_y, file_name_z
        
        call self%destroy()
        
        self%file_name_x = file_name_x
        self%file_name_y = file_name_y
        self%file_name_z = file_name_z
    
    end subroutine numpy_init
    
    !******************************************************************
    !***** Destroy Numpy Type
    !******************************************************************
    subroutine numpy_destroy(self)
        class(numpy_t), intent(inout) :: self
        
        if (allocated(self%file_name_x)) deallocate(self%file_name_x)
        if (allocated(self%file_name_y)) deallocate(self%file_name_y)
        if (allocated(self%file_name_z)) deallocate(self%file_name_z)
        

    
    end subroutine numpy_destroy
    

    
    !******************************************************************
    !***** Write Grid
    !******************************************************************
    subroutine numpy_write_xyz(self, grid)
        class(numpy_t), intent(inout) :: self
        type(grid_t), intent(inout) :: grid
    
        integer :: i, j, k, l, ib, ip
        
        open(2, file=self%file_name_x)
        open(3, file=self%file_name_y)
        open(4, file=self%file_name_z)


        do i = 1, grid%num_rows
            do j = 1, grid%num_cols
                ! write block delimiter
                write (2, '(A)') "HEAD"
                write (3, '(A)') "HEAD"
                write (4, '(A)') "HEAD"
                do k = 1, grid%blocks(i,j)%num_rows
                    do l = 1, grid%blocks(i,j)%num_cols
                        ip = (k-1)*grid%blocks(i,j)%num_cols+l
                        ! write block's real-space coordinates
                        write (2, '(f15.8)', advance='no') grid%blocks(i,j)%x(ip)
                        write (3, '(f15.8)', advance='no') grid%blocks(i,j)%y(ip)
                        write (4, '(f15.8)', advance='no') grid%blocks(i,j)%z(ip)
                    end do
                    write (2, *)
                    write (3, *)
                    write (4, *)
                end do
                write (2, '(A)') "END"
                write (3, '(A)') "END"
                write (4, '(A)') "END"
            end do
        end do
        
        close(2)
        close(3)
        close(4)
    
    end subroutine numpy_write_xyz
    
    
    !******************************************************************
    !***** Write Control Points
    !******************************************************************
    subroutine numpy_write_cp(self, bezier_grid)
        class(numpy_t), intent(inout) :: self
        type(bezier_grid_t), intent(inout) :: bezier_grid
        
        integer :: i, j, k, l, ip
        
        open(2, file=self%file_name_x)
        open(3, file=self%file_name_y)
        open(4, file=self%file_name_z)

        do i = 1, bezier_grid%grid%num_rows
            do j = 1, bezier_grid%grid%num_cols
            
                write (2, '(A)') "HEAD"
                write (3, '(A)') "HEAD"
                write (4, '(A)') "HEAD"
                do k = 1, bezier_grid%num_v_cp
                    do l = 1, bezier_grid%num_u_cp
                        ip = (l-1)*bezier_grid%num_v_cp+k
                        
                        write (2, '(f15.8)', advance='no') bezier_grid%blocks(i,j)%x_cp(ip)
                        write (3, '(f15.8)', advance='no') bezier_grid%blocks(i,j)%y_cp(ip)
                        write (4, '(f15.8)', advance='no') bezier_grid%blocks(i,j)%z_cp(ip)
                    end do
                    write (2, *)
                    write (3, *)
                    write (4, *)
                end do
                write (2, '(A)') "END"
                write (3, '(A)') "END"
                write (4, '(A)') "END"
            end do
        end do
        
        close(2)
        close(3)
        close(4)
        
    
    end subroutine numpy_write_cp
    


    

end module numpy_m




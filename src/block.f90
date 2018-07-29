
module block_m
    implicit none
    
    type, public :: block_t
        integer :: num_rows, num_cols, num_points
        real*8, dimension(:), allocatable :: x, y, z, u, v
    contains
        procedure, public :: init => block_init
        procedure, public :: destroy => block_destroy
    end type block_t
    
contains
    
    !******************************************************************
    !***** Initialize Block Type
    !******************************************************************
    subroutine block_init(self, num_rows, num_cols)
        class(block_t), intent(inout) :: self
        integer, intent(in) :: num_rows, num_cols
        
        call self%destroy()
        
        self%num_rows = num_rows
        self%num_cols = num_cols
        self%num_points = self%num_rows * self%num_cols
        
        allocate(self%x(self%num_points))
        allocate(self%y(self%num_points))
        allocate(self%z(self%num_points))
        
        allocate(self%u(self%num_points))
        allocate(self%v(self%num_points))
    
    end subroutine block_init
    
    !******************************************************************
    !***** Destroy Block Type
    !******************************************************************
    subroutine block_destroy(self)
        class(block_t), intent(inout) :: self 
        
        if (allocated(self%x)) deallocate(self%x)
        if (allocated(self%y)) deallocate(self%y)
        if (allocated(self%z)) deallocate(self%z)
        
        if (allocated(self%u)) deallocate(self%u)
        if (allocated(self%v)) deallocate(self%v)
    
    end subroutine block_destroy
    
end module block_m

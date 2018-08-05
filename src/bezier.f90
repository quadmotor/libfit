! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains a data type representing an Bezier representation
!    of a structured mesh and wrappers procedures relating to Bezier surfaces.

module bezier_m
    use grid_m
    use bezier_helper_m
    use param_m
    implicit none

    ! Bezier representation of structured mesh's block interfaces
    type, private :: bezier_interf_t
   
        real*8, dimension(:, :), allocatable :: xyz_bernstein
        real*8, dimension(:), allocatable :: x_cp, y_cp, z_cp

    end type bezier_interf_t

    ! Bezier representation of structured mesh's block
    type, private :: bezier_block_t
        real*8, dimension(:, :), allocatable :: xyz_bernstein
        real*8, dimension(:), allocatable :: x_cp, y_cp, z_cp

        real*8, dimension(:, :), allocatable :: table_bernstein
        real*8, dimension(:), allocatable :: u_cp, v_cp
        real*8, dimension(:, :), allocatable :: u_table, v_table
        real*8, dimension(:), allocatable :: u_totals, v_totals
    end type bezier_block_t

    ! Bezier representation of structured mesh
    type, public :: bezier_grid_t
        type(grid_t), pointer :: grid

        integer :: u_order, v_order
        integer :: num_u_cp, num_v_cp, num_cp
        integer :: num_rows_table, num_cols_table, num_points_table

        type(bezier_block_t), dimension(:, :), allocatable :: blocks
        type(bezier_interf_t), dimension(:, :), allocatable :: u_interfs, v_interfs

    contains
        procedure, public :: init => bezier_grid_init
        procedure, public :: destroy => bezier_grid_destroy
        procedure, public :: generate_new_grid => bezier_grid_generate_new_grid
        procedure, public :: set_xyz_cp_wrapper => bezier_grid_set_xyz_cp_wrapper
    end type bezier_grid_t


contains

    !******************************************************************
    !***** Initialize Bezier Grid Type
    !******************************************************************
    subroutine bezier_grid_init(self, grid, u_order, v_order, num_rows_table, num_cols_table)
        class(bezier_grid_t), intent(inout) :: self
        type(grid_t), target, intent(in) :: grid
        integer, intent(in) :: u_order, v_order, num_rows_table, num_cols_table
    
        real*8, dimension(:), allocatable :: u_temp, v_temp, x_temp, y_temp, z_temp
        integer :: i, j, k, l, ip

        call self%destroy()

        self%u_order = u_order
        self%v_order = v_order
        self%num_u_cp = u_order + 1
        self%num_v_cp = v_order + 1
        self%num_cp = self%num_u_cp * self%num_v_cp

        self%num_rows_table = num_rows_table
        self%num_cols_table = num_cols_table
        self%num_points_table = num_rows_table * num_cols_table

        self%grid => grid

        allocate(self%blocks(self%grid%num_rows, self%grid%num_cols))
        ! Interfaces running in v direction
        allocate(self%v_interfs(self%grid%num_rows, self%grid%num_cols-1))
        ! Interfaces running in u direction
        allocate(self%u_interfs(self%grid%num_rows-1, self%grid%num_cols))    
        
        ! Initialize each Bezier block u interface representation
        do i = 1, self%grid%num_rows-1
            do j = 1, self%grid%num_cols
                allocate(self%u_interfs(i, j)%x_cp(self%num_u_cp))
                allocate(self%u_interfs(i, j)%y_cp(self%num_u_cp))
                allocate(self%u_interfs(i, j)%z_cp(self%num_u_cp))
                allocate(self%u_interfs(i, j)%xyz_bernstein(self%grid%blocks(i, j)%num_cols, self%num_u_cp))
                
                allocate(u_temp(self%grid%blocks(i,j)%num_cols))
                allocate(x_temp(self%grid%blocks(i,j)%num_cols))
                allocate(y_temp(self%grid%blocks(i,j)%num_cols))
                allocate(z_temp(self%grid%blocks(i,j)%num_cols))
                
                k = self%grid%blocks(i,j)%num_rows
                do l = 1, self%grid%blocks(i,j)%num_cols
                    ip = (k-1)*self%grid%blocks(i,j)%num_cols+l
                    u_temp(l) = self%grid%blocks(i,j)%u(ip)
                    x_temp(l) = self%grid%blocks(i,j)%x(ip)
                    y_temp(l) = self%grid%blocks(i,j)%y(ip)
                    z_temp(l) = self%grid%blocks(i,j)%z(ip)
                end do
                
                
                call set_xyz_bernstein_curve(self%grid%blocks(i,j)%num_cols, self%num_u_cp, &
                    u_temp, self%u_interfs(i,j)%xyz_bernstein)
                    
                call set_xyz_cp(self%grid%blocks(i,j)%num_cols, self%num_u_cp, &
                    x_temp, y_temp, z_temp, self%u_interfs(i,j)%x_cp, self%u_interfs(i,j)%y_cp, &
                    self%u_interfs(i,j)%z_cp, self%u_interfs(i,j)%xyz_bernstein)
            
                deallocate(u_temp)
                deallocate(x_temp)
                deallocate(y_temp)
                deallocate(z_temp)
    
            end do
        end do

        ! Initialize each Bezier block v interface representation
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols-1
                allocate(self%v_interfs(i, j)%x_cp(self%num_v_cp))
                allocate(self%v_interfs(i, j)%y_cp(self%num_v_cp))
                allocate(self%v_interfs(i, j)%z_cp(self%num_v_cp))
                allocate(self%v_interfs(i, j)%xyz_bernstein(self%grid%blocks(i,j)%num_rows, self%num_v_cp))
                
                allocate(v_temp(self%grid%blocks(i,j)%num_rows))
                allocate(x_temp(self%grid%blocks(i,j)%num_rows))
                allocate(y_temp(self%grid%blocks(i,j)%num_rows))
                allocate(z_temp(self%grid%blocks(i,j)%num_rows))
                
                l = self%grid%blocks(i,j)%num_cols
                do k = 1, self%grid%blocks(i,j)%num_rows
                    ip = (k-1)*self%grid%blocks(i,j)%num_cols+l
                    v_temp(k) = self%grid%blocks(i,j)%v(ip)
                    x_temp(k) = self%grid%blocks(i,j)%x(ip)
                    y_temp(k) = self%grid%blocks(i,j)%y(ip)
                    z_temp(k) = self%grid%blocks(i,j)%z(ip)
                end do
                
                call set_xyz_bernstein_curve(self%grid%blocks(i,j)%num_rows, self%num_v_cp, &
                    v_temp, self%v_interfs(i,j)%xyz_bernstein)
                    
                call set_xyz_cp(self%grid%blocks(i,j)%num_rows, self%num_v_cp, &
                    x_temp, y_temp, z_temp, self%v_interfs(i,j)%x_cp, self%v_interfs(i,j)%y_cp, &
                    self%v_interfs(i,j)%z_cp, self%v_interfs(i,j)%xyz_bernstein)

                deallocate(v_temp)
                deallocate(x_temp)
                deallocate(y_temp)
                deallocate(z_temp)
    
            end do
        end do
 
        ! Fix adjacent control point curves
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols-1
                if (i .ne. self%grid%num_rows) then
                    self%v_interfs(i,j)%x_cp(self%num_v_cp) = self%v_interfs(i+1,j)%x_cp(1)
                    self%v_interfs(i,j)%y_cp(self%num_v_cp) = self%v_interfs(i+1,j)%y_cp(1)
                    self%v_interfs(i,j)%z_cp(self%num_v_cp) = self%v_interfs(i+1,j)%z_cp(1)
                end if
            end do
        end do
 
        ! Fix intersecting control points matrices, snap u's to v's
        do i = 1, self%grid%num_rows-1
            do j = 1, self%grid%num_cols
                
                if (j .ne. self%grid%num_cols) then
                    self%u_interfs(i,j)%x_cp(self%num_u_cp) = self%v_interfs(i,j)%x_cp(self%num_v_cp)
                    self%u_interfs(i,j)%y_cp(self%num_u_cp) = self%v_interfs(i,j)%y_cp(self%num_v_cp)
                    self%u_interfs(i,j)%z_cp(self%num_u_cp) = self%v_interfs(i,j)%z_cp(self%num_v_cp)
                end if
                
                if (j .ne. 1) then
                    self%u_interfs(i,j)%x_cp(1) = self%v_interfs(i,j-1)%x_cp(self%num_v_cp)
                    self%u_interfs(i,j)%y_cp(1) = self%v_interfs(i,j-1)%y_cp(self%num_v_cp)
                    self%u_interfs(i,j)%z_cp(1) = self%v_interfs(i,j-1)%z_cp(self%num_v_cp)
                end if
                
            end do
        end do
       
        ! Initialize each Bezier block representation
        do i = 1, self%grid%num_rows
            do j = 1, self%grid%num_cols
        
                allocate(self%blocks(i, j)%x_cp(self%num_cp))
                allocate(self%blocks(i, j)%y_cp(self%num_cp))
                allocate(self%blocks(i, j)%z_cp(self%num_cp))
                allocate(self%blocks(i, j)%xyz_bernstein(self%grid%blocks(i, j)%num_points, self%num_cp))
        
                allocate(self%blocks(i, j)%u_cp(self%num_cp))
                allocate(self%blocks(i, j)%v_cp(self%num_cp))
                allocate(self%blocks(i, j)%table_bernstein(self%num_points_table, self%num_cp))
                allocate(self%blocks(i, j)%u_totals(self%num_rows_table))
                allocate(self%blocks(i, j)%v_totals(self%num_cols_table))
                allocate(self%blocks(i, j)%u_table(self%num_rows_table, self%num_cols_table))
                allocate(self%blocks(i, j)%v_table(self%num_cols_table, self%num_rows_table))
            
                call set_xyz_bernstein(self%grid%blocks(i,j)%num_points, self%num_u_cp, &
                    self%num_v_cp, self%grid%blocks(i,j)%u, self%grid%blocks(i,j)%v, self%blocks(i,j)%xyz_bernstein)
                    
            
                call self%set_xyz_cp_wrapper(i, j)
                
                call set_table(self%num_rows_table, self%num_cols_table, self%num_u_cp, self%num_v_cp, &
                    self%blocks(i, j)%x_cp, self%blocks(i, j)%y_cp, self%blocks(i, j)%z_cp, self%blocks(i, j)%u_table, &
                    self%blocks(i, j)%v_table, self%blocks(i, j)%u_totals, self%blocks(i, j)%v_totals)
                
                call set_table_bernstein(self%num_rows_table, self%num_cols_table, self%num_u_cp, self%num_v_cp, &
                    self%blocks(i, j)%u_table, self%blocks(i, j)%v_table, self%blocks(i, j)%u_totals, self%blocks(i, j)%v_totals, &
                    self%blocks(i, j)%table_bernstein)
                
                call set_table_cp(self%num_rows_table, self%num_cols_table, self%num_cp, self%blocks(i, j)%u_cp, &
                    self%blocks(i, j)%v_cp, self%blocks(i, j)%table_bernstein)
                
            end do
        end do
        
                
    end subroutine bezier_grid_init
    
        
    !******************************************************************
    !***** Subtract Known Control Points
    !******************************************************************
    subroutine bezier_grid_set_xyz_cp_wrapper(self, row, col)
        class(bezier_grid_t), intent(inout) :: self
        integer, intent(in) :: row, col
        
        logical :: flag
        integer :: i, j, k, ip, jp
        integer :: startu, endu, startv, endv
        integer :: num_points, num_u_cp, num_v_cp, num_cp
        real*8, dimension(:), allocatable :: x, y, z, x_cp, y_cp, z_cp
        real*8, dimension(:,:), allocatable :: bernstein
        
        call set_xyz_bernstein(self%grid%blocks(row,col)%num_points, &
                self%num_u_cp, self%num_v_cp, self%grid%blocks(row,col)%u, &
                self%grid%blocks(row,col)%v, self%blocks(row,col)%xyz_bernstein)
        
                
        ! Calculate Bernstein subsection indices
        if (col .eq. 1) then
            startu = 1
        else
            startu = 2
        end if
        
        if (col .eq. self%grid%num_cols) then
            endu = self%num_u_cp
        else
            endu = self%num_u_cp-1
        endif
        
        if (row .eq. 1) then
            startv = 1
        else
            startv = 2
        end if
        
        if (row .eq. self%grid%num_rows) then
            endv = self%num_v_cp
        else
            endv = self%num_v_cp-1
        end if

        num_u_cp = endu - startu + 1
        num_v_cp = endv - startv + 1
        num_cp = num_u_cp * num_v_cp
        num_points = self%grid%blocks(row,col)%num_points
        allocate(bernstein(num_points, num_cp))
        allocate(x_cp(num_cp))
        allocate(y_cp(num_cp))
        allocate(z_cp(num_cp))
        allocate(x(num_points))
        allocate(y(num_points))
        allocate(z(num_points))
        
        ! Subsect Bernstein matrix
        do i = 1, num_points
            do j = startv, endv
                do k = startu, endu
                    ip = (k-1)*self%num_v_cp+j
                    jp = (k-startu)*num_v_cp+(j-startv+1)
                    bernstein(i, jp) = self%blocks(row,col)%xyz_bernstein(i, ip)    
                end do
            end do
        end do
    
        
        ! Subtract from xyz vectors
        do i = 1, num_points
            x(i) = self%grid%blocks(row,col)%x(i)
            y(i) = self%grid%blocks(row,col)%y(i)
            z(i) = self%grid%blocks(row,col)%z(i)
        

            do j = 1, self%num_v_cp
                do k = 1, self%num_u_cp
                    ip = (k-1)*self%num_v_cp+j
                    flag = .false.
                    
                    if (k .lt. startu) then
                        x(i) = x(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col-1)%x_cp(j)
                        y(i) = y(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col-1)%y_cp(j)
                        z(i) = z(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col-1)%z_cp(j)
                        flag = .true.
                    end if
                    
                    
                    if (k .gt. endu) then
                        x(i) = x(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col)%x_cp(j)
                        y(i) = y(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col)%y_cp(j)
                        z(i) = z(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%v_interfs(row, col)%z_cp(j)
                        flag = .true.
                    end if
                    
                    if ((j .lt. startv) .and. (flag .eqv. .false.)) then
                        x(i) = x(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row-1, col)%x_cp(k)
                        y(i) = y(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row-1, col)%y_cp(k)
                        z(i) = z(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row-1, col)%z_cp(k)
                    end if
                    
                    
                    if ((j .gt. endv) .and. (flag .eqv. .false.)) then
                        x(i) = x(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row, col)%x_cp(k)
                        y(i) = y(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row, col)%y_cp(k)
                        z(i) = z(i) - self%blocks(row,col)%xyz_bernstein(i, ip)*self%u_interfs(row, col)%z_cp(k)
                    end if 
                     
                
                end do
            end do

        
        end do
        
        ! Solve for new control points
        call set_xyz_cp(num_points, num_cp, x, y, z, x_cp, y_cp, z_cp, bernstein)
        
        
        ! Copy new control points to final control points
        do i = startv, endv
            do j = startu, endu
                ip = (j-1)*self%num_v_cp+i    
                jp = (j-startu)*num_v_cp+(i-startv+1)
                
                if (((j .ge. startu) .or. (j .le. endu)) .and. &
                    ((i .ge. startv) .or. (i .le. endv))) then
                    
                    self%blocks(row,col)%x_cp(ip) = x_cp(jp)
                    self%blocks(row,col)%y_cp(ip) = y_cp(jp)
                    self%blocks(row,col)%z_cp(ip) = z_cp(jp)
                end if

            end do
        end do
        
                
        ! Copy top interface control points to final control points
        if (row .ne. 1) then
            i = 1
            do j = 1, self%num_u_cp
                ip = (j-1)*self%num_v_cp+i
            
                self%blocks(row,col)%x_cp(ip) = self%u_interfs(row-1, col)%x_cp(j)
                self%blocks(row,col)%y_cp(ip) = self%u_interfs(row-1, col)%y_cp(j)
                self%blocks(row,col)%z_cp(ip) = self%u_interfs(row-1, col)%z_cp(j)
            end do    
        end if
        
        ! Copy bottom interface control points to final control points
        if (row .ne. self%grid%num_rows) then
            i = self%num_v_cp
            do j = 1, self%num_u_cp
                ip = (j-1)*self%num_v_cp+i
                
                self%blocks(row,col)%x_cp(ip) = self%u_interfs(row, col)%x_cp(j)
                self%blocks(row,col)%y_cp(ip) = self%u_interfs(row, col)%y_cp(j)
                self%blocks(row,col)%z_cp(ip) = self%u_interfs(row, col)%z_cp(j)
            end do        
        end if
        
        ! Copy left interface control points to final control points
        if (col .ne. 1) then
            j = 1
            do i = 1, self%num_v_cp
                ip = (j-1)*self%num_v_cp+i
            
                self%blocks(row,col)%x_cp(ip) = self%v_interfs(row, col-1)%x_cp(i)
                self%blocks(row,col)%y_cp(ip) = self%v_interfs(row, col-1)%y_cp(i)
                self%blocks(row,col)%z_cp(ip) = self%v_interfs(row, col-1)%z_cp(i)
            end do
        end if
        
        ! Copy right interface control points to final control points
        if (col .ne. self%grid%num_cols) then
            j = self%num_u_cp
            do i = 1, self%num_v_cp
                ip = (j-1)*self%num_v_cp+i
            
                self%blocks(row,col)%x_cp(ip) = self%v_interfs(row, col)%x_cp(i)
                self%blocks(row,col)%y_cp(ip) = self%v_interfs(row, col)%y_cp(i)
                self%blocks(row,col)%z_cp(ip) = self%v_interfs(row, col)%z_cp(i)
            end do
            
        end if
        
        
        deallocate(bernstein)
        deallocate(x_cp)
        deallocate(y_cp)
        deallocate(z_cp)
        deallocate(x)
        deallocate(y)
        deallocate(z)

    end subroutine bezier_grid_set_xyz_cp_wrapper
    
    !******************************************************************
    !***** Destroy Bezier Grid Type
    !******************************************************************
    subroutine bezier_grid_destroy(self)
        class(bezier_grid_t), intent(inout) :: self
                
        integer :: i, j
        
        
        if (associated(self%grid) .and. allocated(self%v_interfs)) then
            do i = 1, self%grid%num_rows
                do j = 1, self%grid%num_cols-1
                    if (allocated(self%v_interfs(i, j)%xyz_bernstein)) deallocate(self%v_interfs(i,j)%xyz_bernstein)
                    if (allocated(self%v_interfs(i, j)%x_cp)) deallocate(self%v_interfs(i, j)%x_cp)
                    if (allocated(self%v_interfs(i, j)%y_cp)) deallocate(self%v_interfs(i, j)%y_cp)
                    if (allocated(self%v_interfs(i, j)%z_cp)) deallocate(self%v_interfs(i, j)%z_cp)
                end do
            end do
        end if
        
        if (associated(self%grid) .and. allocated(self%u_interfs)) then
            do i = 1, self%grid%num_rows-1
                do j = 1, self%grid%num_cols
                    if (allocated(self%u_interfs(i,j)%xyz_bernstein)) deallocate(self%u_interfs(i,j)%xyz_bernstein)
                    if (allocated(self%u_interfs(i,j)%x_cp)) deallocate(self%u_interfs(i,j)%x_cp)
                    if (allocated(self%u_interfs(i,j)%y_cp)) deallocate(self%u_interfs(i,j)%y_cp)
                    if (allocated(self%u_interfs(i,j)%z_cp)) deallocate(self%u_interfs(i,j)%z_cp)
                end do
            end do
        end if

        if (associated(self%grid) .and. allocated(self%blocks)) then
            do i = 1, self%grid%num_rows
                do j = 1, self%grid%num_cols
                    if (allocated(self%blocks(i,j)%xyz_bernstein)) deallocate(self%blocks(i,j)%xyz_bernstein)
                    if (allocated(self%blocks(i,j)%x_cp)) deallocate(self%blocks(i,j)%x_cp)
                    if (allocated(self%blocks(i,j)%y_cp)) deallocate(self%blocks(i,j)%y_cp)
                    if (allocated(self%blocks(i,j)%z_cp)) deallocate(self%blocks(i,j)%z_cp)
                    
        
                    if (allocated(self%blocks(i,j)%table_bernstein)) deallocate(self%blocks(i,j)%table_bernstein)
                    if (allocated(self%blocks(i,j)%u_cp)) deallocate(self%blocks(i,j)%u_cp)
                    if (allocated(self%blocks(i,j)%v_cp)) deallocate(self%blocks(i,j)%v_cp)
                    if (allocated(self%blocks(i,j)%u_totals)) deallocate(self%blocks(i,j)%u_totals)
                    if (allocated(self%blocks(i,j)%v_totals)) deallocate(self%blocks(i,j)%v_totals)
                    if (allocated(self%blocks(i,j)%u_table)) deallocate(self%blocks(i,j)%u_table)
                    if (allocated(self%blocks(i,j)%v_table)) deallocate(self%blocks(i,j)%v_table)
                end do
            end do
        end if
        
        if (allocated(self%blocks)) deallocate(self%blocks)
        if (allocated(self%u_interfs)) deallocate(self%u_interfs)
        if (allocated(self%v_interfs)) deallocate(self%v_interfs)
         
        nullify(self%grid)
        
    end subroutine bezier_grid_destroy
    

    
    
    
    !******************************************************************
    !***** Generate New Grid
    !******************************************************************
    subroutine bezier_grid_generate_new_grid(self, grid, num_rows, num_cols, u_gen, v_gen)
        class(bezier_grid_t), intent(inout) :: self
        type(grid_t), intent(inout) :: grid
        integer, intent(in) :: num_rows, num_cols
        real*8, dimension(num_rows*num_cols), intent(inout) :: u_gen, v_gen

        integer :: i, j, k, l, ip
        real*8, dimension(:), allocatable :: x_c, y_c, z_c, u_act_c, v_act_c
        real*8, dimension(num_rows*num_cols) :: u_act, v_act
        
        call grid%init(self%grid%num_rows, self%grid%num_cols)
        
        ! Generate new mesh's xyz values
        do i = 1, grid%num_rows
            do j = 1, grid%num_cols
        
                call grid%blocks(i, j)%init(num_rows, num_cols)
        
                call uv_casteljau(num_rows*num_cols, self%num_u_cp, self%num_v_cp, &
                    self%blocks(i, j)%u_cp, self%blocks(i, j)%v_cp, u_gen, v_gen, u_act, v_act)
                
                call normalize_param(num_rows, num_cols, u_act, v_act)
                    
                call xyz_casteljau(num_rows*num_cols, self%num_u_cp, self%num_v_cp, &
                    self%blocks(i, j)%x_cp, self%blocks(i, j)%y_cp, self%blocks(i, j)%z_cp, &
                    u_act, v_act, grid%blocks(i, j)%x, grid%blocks(i, j)%y, grid%blocks(i, j)%z)
                
            end do
        end do
        
        ! Generate new mesh's u curves
        do i = 1, grid%num_rows-1
            do j = 1, grid%num_cols
        
                ! Allocate memory for new xyz data
                allocate(x_c(grid%blocks(i,j)%num_cols))
                allocate(y_c(grid%blocks(i,j)%num_cols))
                allocate(z_c(grid%blocks(i,j)%num_cols))
                allocate(u_act_c(grid%blocks(i,j)%num_cols))
                allocate(v_act_c(grid%blocks(i,j)%num_cols))
                
                ! copy uv data from block                        REPLACE THIS WITH ARC LENGTH TABLE ALONG CURVE
                l = grid%blocks(i,j)%num_rows
                do k = 1, grid%blocks(i,j)%num_cols
                    ip = (l-1)*grid%blocks(i,j)%num_cols+k
                    u_act_c(k) = u_act(ip)
                    v_act_c(k) = v_act(ip)
                end do
                
                call xyz_casteljau_curve(grid%blocks(i,j)%num_cols, self%num_u_cp, self%u_interfs(i,j)%x_cp, &
                    self%u_interfs(i,j)%y_cp, self%u_interfs(i,j)%z_cp, u_act_c, x_c, y_c, z_c)
                
                    
                ! Copy into blocks
                l = grid%blocks(i,j)%num_rows
                do k = 1, grid%blocks(i,j)%num_cols
                    ip = (l-1)*grid%blocks(i,j)%num_cols+k
                    grid%blocks(i,j)%x(ip) = x_c(k)
                    grid%blocks(i,j)%y(ip) = y_c(k)
                    grid%blocks(i,j)%z(ip) = z_c(k)
                end do
                l = 1
                do k = 1, grid%blocks(i+1,j)%num_cols
                    ip = (l-1)*grid%blocks(i+1,j)%num_cols+k
                    grid%blocks(i+1,j)%x(ip) = x_c(k)
                    grid%blocks(i+1,j)%y(ip) = y_c(k)
                    grid%blocks(i+1,j)%z(ip) = z_c(k)
                end do
                
                deallocate(x_c)
                deallocate(y_c)
                deallocate(z_c)
                deallocate(u_act_c)
                deallocate(v_act_c)

            end do
        end do
        
        ! Generate new mesh's v curves
        do i = 1, grid%num_rows
            do j = 1, grid%num_cols-1
        
                allocate(x_c(grid%blocks(i,j)%num_rows))
                allocate(y_c(grid%blocks(i,j)%num_rows))
                allocate(z_c(grid%blocks(i,j)%num_rows))
                allocate(u_act_c(grid%blocks(i,j)%num_rows))
                allocate(v_act_c(grid%blocks(i,j)%num_rows))
                
                ! Copy uv data from block                        REPLACE THIS WITH ARC LENGTH TABLE ALONG CURVE
                k = grid%blocks(i,j)%num_cols
                do l = 1, grid%blocks(i,j)%num_rows
                    ip = (l-1)*grid%blocks(i,j)%num_cols+k
                    u_act_c(l) = u_act(ip)
                    v_act_c(l) = v_act(ip)
                end do
                
        
                call xyz_casteljau_curve(grid%blocks(i,j)%num_rows, self%num_v_cp, self%v_interfs(i,j)%x_cp, &
                    self%v_interfs(i,j)%y_cp, self%v_interfs(i,j)%z_cp, v_act_c, x_c, y_c, z_c)
                
                    
                ! Copy into blocks
                k = grid%blocks(i,j)%num_cols
                do l = 1, grid%blocks(i,j)%num_rows
                    ip = (l-1)*grid%blocks(i,j)%num_cols+k
                    grid%blocks(i,j)%x(ip) = x_c(l)
                    grid%blocks(i,j)%y(ip) = y_c(l)
                    grid%blocks(i,j)%z(ip) = z_c(l)
                end do
                k = 1
                do l = 1, grid%blocks(i,j+1)%num_rows
                    ip = (l-1)*grid%blocks(i,j+1)%num_cols+k
                    grid%blocks(i,j+1)%x(ip) = x_c(l)
                    grid%blocks(i,j+1)%y(ip) = y_c(l)
                    grid%blocks(i,j+1)%z(ip) = z_c(l)
                end do
                
                deallocate(x_c)
                deallocate(y_c)
                deallocate(z_c)
                deallocate(u_act_c)
                deallocate(v_act_c)

            end do
        end do
        
            
    end subroutine bezier_grid_generate_new_grid
    


end module bezier_m




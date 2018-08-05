
! Project: libfit
! File: examples/torus/main.f90 
! Author: Cameron Mackintosh (2017)
! Description: This program demonstrates an example usage of the
!    libfit library.


program torus
    use grid_m
    use plot3d_m
    use numpy_m
    use bezier_m
    use param_m
    implicit none

    integer, parameter :: num_rows_new = 21, num_cols_new = 21
    type(grid_t) :: grid_in, grid_out
    type(plot3d_t) :: plot3d_in
    type(numpy_t) :: numpy_in, numpy_out, numpy_cp
    type(bezier_grid_t) :: bezier
    real*8, dimension(num_rows_new*num_cols_new) :: u_new, v_new
    
    ! Initialize readers
    call plot3d_in%init(grid_in, "mesh.xyz")
    call numpy_in%init('output/x_in.dat', 'output/y_in.dat', 'output/z_in.dat')
    call numpy_out%init('output/x_out.dat', 'output/y_out.dat', 'output/z_out.dat')
    call numpy_cp%init('output/x_cp.dat', 'output/y_cp.dat', 'output/z_cp.dat')

    ! Read coarse mesh from plot3d file
    call plot3d_in%read_xyz()
    ! Write coarse mesh to numpy file -- used for easy plotting with matplotlib
    call numpy_in%write_xyz(grid_in)
    ! Parameterize coarse mesh by arc length
    call param_grid(grid_in)
    ! Initialize bezier class with uv orders 5x5 and arc length table dimensions 50x50
    call bezier%init(grid_in, 5, 5, 50, 50)
  
    ! Create a generating computational-space distribution
    call param_uniform(num_rows_new, num_cols_new, u_new, v_new)
    ! Generate fine mesh using computational-space distribution
    call bezier%generate_new_grid(grid_out, num_rows_new, num_cols_new, u_new, v_new)
    
    ! Write fine mesh to numpy file -- used for easy plotting with matplotlib
    call numpy_out%write_xyz(grid_out)
    ! Write fine mesh's control points to numpy file -- used for easy plotting with matplotlib
    call numpy_cp%write_cp(bezier)
    
    call plot3d_in%destroy()
    call numpy_in%destroy()
    call numpy_out%destroy()
    call bezier%destroy()
    call grid_in%destroy()
    call grid_out%destroy()
    
end program torus


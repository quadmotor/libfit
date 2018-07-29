
! o888  o88   oooo       ooooooooooo o88    o8  
!  888  oooo   888ooooo   888    88  oooo o888oo
!  888   888   888    888 888ooo8     888  888  
!  888   888   888    888 888         888  888  
! o888o o888o o888ooo88  o888o       o888o  888o
!
! Created: June 2017
! Author: Cameron Mackintosh
! Description: This module contains procedures relating to Bezier surfaces.


module bezier_helper_m

    implicit none


contains

!******************************************************************    
!******************************************************************
!***** XYZ Subroutines
!******************************************************************    
!******************************************************************    

    !******************************************************************
    !***** Set XYZ Bernstein Curve
    !******************************************************************
    subroutine set_xyz_bernstein_curve(num_points, num_t_cp, t, bernstein)
        integer, intent(in) :: num_points, num_t_cp
        real*8, dimension(num_points) :: t
        real*8, dimension(num_points, num_t_cp), intent(out) :: bernstein
        
        real*8 :: bi, bj
        integer :: ip, ic, t_order
        
        t_order = num_t_cp - 1
        
        do ip = 1, num_points
            do ic = 1, num_t_cp
                bernstein(ip, ic) = binomial(t_order, ic-1) * ((1.d0-t(ip))**(t_order-(ic-1)) * t(ip)**(ic-1))
            end do
        end do

    end subroutine set_xyz_bernstein_curve

    !******************************************************************
    !***** Set XYZ Bernstein
    !******************************************************************
    subroutine set_xyz_bernstein(num_points, num_u_cp, num_v_cp, u, v, bernstein)
        integer, intent(in) :: num_points, num_u_cp, num_v_cp
        real*8, dimension(num_points) :: u, v
        real*8, dimension(num_points, num_u_cp*num_v_cp), intent(out) :: bernstein
        
        real*8 :: bi, bj
        integer :: i, j, ip, ic, u_order, v_order
        
        u_order = num_u_cp - 1
        v_order = num_v_cp - 1
        
        do ip = 1, num_points
            do i = 1, num_u_cp
                do j = 1, num_v_cp
                    ic = (i-1)*num_v_cp + j
                    bi = binomial(u_order, i-1) * ((1.d0-u(ip))**(u_order-(i-1)) * u(ip)**(i-1))
                    bj = binomial(v_order, j-1) * ((1.d0-v(ip))**(v_order-(j-1)) * v(ip)**(j-1))
                    bernstein(ip, ic) = bi*bj
                end do
            end do
        end do

    end subroutine set_xyz_bernstein
    

    !******************************************************************
    !***** Set XYZ Control Points
    !******************************************************************
    subroutine set_xyz_cp(num_points, num_cp, x, y, z, x_cp, y_cp, z_cp, bernstein)
        use lsqr_m, only: lsqr
        integer, intent(in) :: num_points, num_cp
        real*8, dimension(num_points) :: x, y, z
        real*8, dimension(num_cp) :: x_cp, y_cp, z_cp
        real*8, dimension(num_points, num_cp) :: bernstein
        
        real*8, dimension(num_cp) :: sex, sey, sez
        real*8 :: anorm, acond, rnorm, arnorm, xnorm, eps
        integer :: info, iters, itlim
        
        eps = 1.d-16
        itlim = 4 * num_cp
        
        
        call lsqr(num_points, num_cp, bernstein, x, 0.d0, .true., x_cp, &
            sex, eps, eps, 1.d6, itlim, 0, info, iters, anorm, acond, rnorm, arnorm, xnorm)

        call lsqr(num_points, num_cp, bernstein, y, 0.d0, .true., y_cp, &
            sey, eps, eps, 1.d6, itlim, 0, info, iters, anorm, acond, rnorm, arnorm, xnorm)
            
        call lsqr(num_points, num_cp, bernstein, z, 0.d0, .true., z_cp, &
            sez, eps, eps, 1.d6, itlim, 0, info, iters, anorm, acond, rnorm, arnorm, xnorm)
        
    end subroutine set_xyz_cp
    
    
    !******************************************************************
    !***** XYZ Casteljau Curve
    !******************************************************************
    subroutine xyz_casteljau_curve(num_points, num_t_cp, x_cp, y_cp, z_cp, t, x, y, z)
        integer, intent(in) :: num_points, num_t_cp
        real*8, dimension(num_t_cp), intent(in) :: x_cp, y_cp, z_cp
        real*8, dimension(num_points), intent(inout) :: t, x, y, z
        
        integer :: i, j, ip, ic, num_t_cp_t
        real*8, dimension(num_t_cp) :: x_cp_t, y_cp_t, z_cp_t
        
        
        do i = 1, num_points
            
            do j = 1, num_t_cp
                x_cp_t(j) = x_cp(j)
                y_cp_t(j) = y_cp(j)
                z_cp_t(j) = z_cp(j)
            end do
            
            num_t_cp_t = num_t_cp
            do
                if (num_t_cp_t .lt. 1) exit
                do j = 1, num_t_cp_t-1
                    x_cp_t(j) = (1.d0-t(i))*x_cp_t(j) + t(i)*x_cp_t(j+1)
                    y_cp_t(j) = (1.d0-t(i))*y_cp_t(j) + t(i)*y_cp_t(j+1)
                    z_cp_t(j) = (1.d0-t(i))*z_cp_t(j) + t(i)*z_cp_t(j+1)
                end do
                num_t_cp_t = num_t_cp_t - 1
            end do
            
            x(i) = x_cp_t(1)
            y(i) = y_cp_t(1)
            z(i) = z_cp_t(1)
        end do
    
    
    end subroutine xyz_casteljau_curve
    
    !******************************************************************
    !***** XYZ Casteljau 
    !******************************************************************
    subroutine xyz_casteljau(num_points, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u, v, x, y, z)
        integer, intent(in) :: num_points, num_u_cp, num_v_cp
        real*8, dimension(num_u_cp*num_v_cp), intent(in) :: x_cp, y_cp, z_cp
        real*8, dimension(num_points), intent(inout) :: x, y, z, u, v
        
        integer :: i, j, ip, ic, t_num_u_cp, t_num_v_cp
        real*8, dimension(num_u_cp) :: x_cp_u, y_cp_u, z_cp_u
        real*8, dimension(num_v_cp) :: x_cp_v, y_cp_v, z_cp_v
    
        do ip = 1, num_points
            
            do i = 1, num_u_cp
                do j = 1, num_v_cp
                    ic = (i-1)*num_v_cp + j
                    x_cp_v(j) = x_cp(ic)
                    y_cp_v(j) = y_cp(ic)
                    z_cp_v(j) = z_cp(ic)
                end do
                
                t_num_v_cp = num_v_cp
                do
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        x_cp_v(j) = (1.d0-v(ip))*x_cp_v(j) + v(ip)*x_cp_v(j+1)
                        y_cp_v(j) = (1.d0-v(ip))*y_cp_v(j) + v(ip)*y_cp_v(j+1)
                        z_cp_v(j) = (1.d0-v(ip))*z_cp_v(j) + v(ip)*z_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp - 1
                end do
                
                x_cp_u(i) = x_cp_v(1)
                y_cp_u(i) = y_cp_v(1)
                z_cp_u(i) = z_cp_v(1)
            end do
            
            
            t_num_u_cp = num_u_cp
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    x_cp_u(i) = (1.d0-u(ip))*x_cp_u(i) + u(ip)*x_cp_u(i+1)
                    y_cp_u(i) = (1.d0-u(ip))*y_cp_u(i) + u(ip)*y_cp_u(i+1)
                    z_cp_u(i) = (1.d0-u(ip))*z_cp_u(i) + u(ip)*z_cp_u(i+1)
                end do
                t_num_u_cp = t_num_u_cp-1
            end do
            
            x(ip) = x_cp_u(1)
            y(ip) = y_cp_u(1)
            z(ip) = z_cp_u(1)
            
        end do
    
    
    end subroutine xyz_casteljau
    
    !******************************************************************
    !***** XYZ Casteljau Deriv U
    !******************************************************************
    subroutine xyz_casteljau_du(num_points, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u, v_fixed, x, y, z)
        integer, intent(in) :: num_points, num_u_cp, num_v_cp
        real*8, dimension(num_u_cp*num_v_cp) :: x_cp, y_cp, z_cp
        real*8, dimension(num_points), intent(in) :: u
        real*8, intent(in) :: v_fixed
        real*8, dimension(num_points), intent(inout) :: x, y, z
        
        integer :: i, j, ip, ic, t_num_u_cp, t_num_v_cp, u_order
        real*8 :: x1, y1, z1, x2, y2, z2
        real*8, dimension(num_u_cp) :: x_cp_u, y_cp_u, z_cp_u
        real*8, dimension(num_v_cp) :: x_cp_v, y_cp_v, z_cp_v
        
        u_order = num_u_cp - 1
        
        do ip = 1, num_points
            
            do i = 1, num_u_cp-1
                do j = 1, num_v_cp
                    ic = (i-1)*num_v_cp + j
                    x_cp_v(j) = x_cp(ic+num_v_cp)
                    y_cp_v(j) = y_cp(ic+num_v_cp)
                    z_cp_v(j) = z_cp(ic+num_v_cp)
                end do            
                t_num_v_cp = num_v_cp
                do
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        x_cp_v(j) = (1.d0-v_fixed)*x_cp_v(j) + v_fixed*x_cp_v(j+1)
                        y_cp_v(j) = (1.d0-v_fixed)*y_cp_v(j) + v_fixed*y_cp_v(j+1)
                        z_cp_v(j) = (1.d0-v_fixed)*z_cp_v(j) + v_fixed*z_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp-1
                end do
                
                x_cp_u(i) = x_cp_v(1)
                y_cp_u(i) = y_cp_v(1)
                z_cp_u(i) = z_cp_v(1)
            end do
            t_num_u_cp = num_u_cp-1
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    x_cp_u(i) = (1.d0-u(ip))*x_cp_u(i) + u(ip)*x_cp_u(i+1)
                    y_cp_u(i) = (1.d0-u(ip))*y_cp_u(i) + u(ip)*y_cp_u(i+1)
                    z_cp_u(i) = (1.d0-u(ip))*z_cp_u(i) + u(ip)*z_cp_u(i+1)
                end do
                t_num_u_cp = t_num_u_cp-1            
            end do
            x1 = x_cp_u(1)
            y1 = y_cp_u(1)
            z1 = z_cp_u(1)
            
            
            do i = 1, num_u_cp-1
                do j = 1, num_v_cp
                    ic = (i-1)*num_v_cp + j 
                    x_cp_v(j) = x_cp(ic)
                    y_cp_v(j) = y_cp(ic)
                    z_cp_v(j) = z_cp(ic)
                end do
                
                t_num_v_cp = num_v_cp
                do 
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        x_cp_v(j) = (1.d0-v_fixed)*x_cp_v(j) + v_fixed*x_cp_v(j+1)
                        y_cp_v(j) = (1.d0-v_fixed)*y_cp_v(j) + v_fixed*y_cp_v(j+1)
                        z_cp_v(j) = (1.d0-v_fixed)*z_cp_v(j) + v_fixed*z_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp-1
                end do
                x_cp_u(i) = x_cp_v(1)
                y_cp_u(i) = y_cp_v(1)
                z_cp_u(i) = z_cp_v(1)
            end do
            t_num_u_cp = num_u_cp-1
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    x_cp_u(i) = (1.d0-u(ip))*x_cp_u(i) + u(ip)*x_cp_u(i+1)
                    y_cp_u(i) = (1.d0-u(ip))*y_cp_u(i) + u(ip)*y_cp_u(i+1)
                    z_cp_u(i) = (1.d0-u(ip))*z_cp_u(i) + u(ip)*z_cp_u(i+1) 
                end do
                t_num_u_cp = t_num_u_cp-1
            end do
            x2 = x_cp_u(1)
            y2 = y_cp_u(1)
            z2 = z_cp_u(1)
            
            x(ip) = u_order * (x1-x2)
            y(ip) = u_order * (y1-y2)
            z(ip) = u_order * (z1-z2)
            
        end do
        
    end subroutine xyz_casteljau_du
    
    !******************************************************************
    !***** XYZ Casteljau Deriv V
    !******************************************************************
    subroutine xyz_casteljau_dv(num_points, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u_fixed, v, x, y, z)
        integer, intent(in) :: num_points, num_u_cp, num_v_cp
        real*8, dimension(num_u_cp*num_v_cp) :: x_cp, y_cp, z_cp
        real*8, intent(in) :: u_fixed
        real*8, dimension(num_points), intent(in) :: v
        real*8, dimension(num_points), intent(inout) :: x, y, z
        
        integer :: i, j, ip, ic, t_num_u_cp, t_num_v_cp, v_order
        real*8 :: x1, y1, z1, x2, y2, z2
        real*8, dimension(num_u_cp) :: x_cp_u, y_cp_u, z_cp_u
        real*8, dimension(num_v_cp) :: x_cp_v, y_cp_v, z_cp_v
        
        v_order = num_v_cp - 1
        
        do ip = 1, num_points
            
            ic = 1
            do i = 1, num_u_cp
                do j = 1, num_v_cp-1
                    x_cp_v(j) = x_cp(ic+1)
                    y_cp_v(j) = y_cp(ic+1)
                    z_cp_v(j) = z_cp(ic+1)
                    ic = ic + 1
                end do            
                ic = ic + 1

                t_num_v_cp = num_v_cp-1
                do
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        x_cp_v(j) = (1.d0-v(ip))*x_cp_v(j) + v(ip)*x_cp_v(j+1)
                        y_cp_v(j) = (1.d0-v(ip))*y_cp_v(j) + v(ip)*y_cp_v(j+1)
                        z_cp_v(j) = (1.d0-v(ip))*z_cp_v(j) + v(ip)*z_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp-1
                end do
                x_cp_u(i) = x_cp_v(1)
                y_cp_u(i) = y_cp_v(1)
                z_cp_u(i) = z_cp_v(1)
            end do
            t_num_u_cp = num_u_cp
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    x_cp_u(i) = (1.d0-u_fixed)*x_cp_u(i) + u_fixed*x_cp_u(i+1)
                    y_cp_u(i) = (1.d0-u_fixed)*y_cp_u(i) + u_fixed*y_cp_u(i+1)
                    z_cp_u(i) = (1.d0-u_fixed)*z_cp_u(i) + u_fixed*z_cp_u(i+1)
                end do
                t_num_u_cp = t_num_u_cp-1            
            end do
            x1 = x_cp_u(1)
            y1 = y_cp_u(1)
            z1 = z_cp_u(1)
            
            ic = 1
            do i = 1, num_u_cp
                do j = 1, num_v_cp-1
                    x_cp_v(j) = x_cp(ic)
                    y_cp_v(j) = y_cp(ic)
                    z_cp_v(j) = z_cp(ic)
                    ic = ic + 1
                end do
                ic = ic + 1
                
                t_num_v_cp = num_v_cp-1
                do 
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        x_cp_v(j) = (1.d0-v(ip))*x_cp_v(j) + v(ip)*x_cp_v(j+1)
                        y_cp_v(j) = (1.d0-v(ip))*y_cp_v(j) + v(ip)*y_cp_v(j+1)
                        z_cp_v(j) = (1.d0-v(ip))*z_cp_v(j) + v(ip)*z_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp-1
                end do
                x_cp_u(i) = x_cp_v(1)
                y_cp_u(i) = y_cp_v(1)
                z_cp_u(i) = z_cp_v(1)
            end do
            t_num_u_cp = num_u_cp
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    x_cp_u(i) = (1.d0-u_fixed)*x_cp_u(i) + u_fixed*x_cp_u(i+1)
                    y_cp_u(i) = (1.d0-u_fixed)*y_cp_u(i) + u_fixed*y_cp_u(i+1)
                    z_cp_u(i) = (1.d0-u_fixed)*z_cp_u(i) + u_fixed*z_cp_u(i+1)
                end do
                t_num_u_cp = t_num_u_cp-1
            end do
            x2 = x_cp_u(1)
            y2 = y_cp_u(1)
            z2 = z_cp_u(1)
            
            x(ip) = v_order * (x1-x2)
            y(ip) = v_order * (y1-y2)
            z(ip) = v_order * (z1-z2)
            
        end do
    
    end subroutine xyz_casteljau_dv

!******************************************************************    
!******************************************************************
!***** UV Subroutines
!******************************************************************    
!******************************************************************    
    
    !******************************************************************
    !***** Set UV Bernstein
    !******************************************************************
    subroutine set_table_bernstein(num_rows, num_cols, num_u_cp, num_v_cp, u_table, v_table, u_totals, v_totals, bernstein)
        integer, intent(in) :: num_rows, num_cols, num_u_cp, num_v_cp
        real*8, dimension(num_rows, num_cols), intent(in) :: u_table, v_table 
        real*8, dimension(num_rows), intent(in) :: u_totals
        real*8, dimension(num_cols), intent(in) :: v_totals
        real*8, dimension(num_rows*num_cols, num_u_cp*num_v_cp), intent(inout) :: bernstein 
    
        
        integer :: i, j, k, l, ip, ic, u_order, v_order
        real*8 :: bi, bj
        real*8, dimension(num_rows*num_cols) :: u, v
        
        u_order = num_u_cp - 1
        v_order = num_v_cp - 1
        
        ip = 1
        do i = 1, num_rows
            do j = 1, num_cols
                u(ip) = real(u_table(i, j)) / u_totals(i)
                v(ip) = real(v_table(j, i)) / v_totals(j)
                ip = ip + 1
            enddo
        enddo

        ip = 1
        do ip = 1, num_rows*num_cols
            ic = 1
            do i = 1, num_u_cp
                do j = 1, num_v_cp
                    bi = binomial(u_order, i-1) * (1.d0 - u(ip))**(u_order-(i-1)) * u(ip)**(i-1)
                    bj = binomial(v_order, j-1) * (1.d0 - v(ip))**(v_order-(j-1)) * v(ip)**(j-1)
                    bernstein(ip, ic) = bi*bj
                    ic = ic + 1
                end do
            enddo
        end do
    end subroutine set_table_bernstein
    

    !******************************************************************
    !***** Set UV Control Points
    !******************************************************************
    subroutine set_table_cp(num_rows, num_cols, num_cp, u_cp, v_cp, bernstein)
        use lsqr_m, only: lsqr
        integer, intent(inout) :: num_rows, num_cols, num_cp
        real*8, dimension(num_cp), intent(inout) :: u_cp, v_cp
        real*8, dimension(num_rows*num_cols, num_cp), intent(in) :: bernstein
    
        real*8, dimension(num_cp) :: seu, sev
        real*8 :: anorm, acond, rnorm, arnorm, xnorm, eps
        integer :: info, iters, itlim
        
        real*8, dimension(num_rows*num_cols) :: u, v
        integer :: i, j, ip

        eps = 1.d-16
        itlim = 10 * num_cp
        
        ip = 1
        do i = 1, num_rows
            do j = 1, num_cols
                u(ip) = real(j-1) / (num_cols-1)
                v(ip) = real(i-1) / (num_rows-1)
                ip = ip + 1
            enddo
        enddo
    
        
        call lsqr(num_rows*num_cols, num_cp, bernstein, u, 0.d0, .true., u_cp, &
            seu, eps, eps, 1.d6, itlim, 0, info, iters, anorm, acond, rnorm, arnorm, xnorm)
            
        call lsqr(num_rows*num_cols, num_cp, bernstein, v, 0.d0, .true., v_cp, &
            sev, eps, eps, 1.d6, itlim, 0, info, iters, anorm, acond, rnorm, arnorm, xnorm)
    
    end subroutine set_table_cp
    
    
    !******************************************************************
    !***** Set Arc Length Table
    !******************************************************************
    subroutine set_table(num_rows, num_cols, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u_table, v_table, u_totals, v_totals)
        integer, intent(in) :: num_rows, num_cols, num_u_cp, num_v_cp
        real*8, dimension(num_u_cp*num_v_cp) :: x_cp, y_cp, z_cp
        real*8, dimension(num_rows, num_cols), intent(inout) :: u_table, v_table
        real*8, dimension(num_rows), intent(inout) :: u_totals
        real*8, dimension(num_cols), intent(inout) :: v_totals 
        
        real*8, dimension(num_cols) :: u, xu, yu, zu
        real*8, dimension(num_rows) :: v, xv, yv, zv
        real*8 :: arc_len
        integer :: i, j, ip
        
        do i = 1, num_cols
            u(i) = real(i-1) / (num_cols-1)
        end do
        do i = 1, num_rows
            v(i) = real(i-1) / (num_rows-1)
        end do
        
        do i = 1, num_rows
            call xyz_casteljau_du(num_cols, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u, v(i), xu, yu, zu)
            
            do j = 1, num_cols
                if (j == 1) then
                    arc_len = 0.d0
                else 
                    arc_len = arc_len + sqrt(xu(j)**2 + yu(j)**2 + zu(j)**2) * (1.d0 / num_cols)
                end if
                u_table(i, j) = arc_len
            end do
            u_totals(i) = arc_len
        end do
        
        do i = 1, num_cols
        
            call xyz_casteljau_dv(num_rows, num_u_cp, num_v_cp, x_cp, y_cp, z_cp, u(i), v, xv, yv, zv)
            
            do j = 1, num_rows
                if (j == 1) then
                    arc_len = 0.d0
                else
                    arc_len = arc_len + sqrt(xv(j)**2 + yv(j)**2 + zv(j)**2) * (1.d0 / num_rows)
                end if
                v_table(i, j) = arc_len
            end do
            v_totals(i) = arc_len
        end do
        

    end subroutine set_table


    !******************************************************************
    !***** Table Casteljau
    !******************************************************************
     subroutine uv_casteljau(num_points, num_u_cp, num_v_cp, u_cp, v_cp, u, v, u_new, v_new)
        integer, intent(in) :: num_points, num_u_cp, num_v_cp
        real*8, dimension(num_u_cp*num_v_cp), intent(in) :: u_cp, v_cp
        real*8, dimension(num_points), intent(inout) :: u, v
        real*8, dimension(num_points), intent(out) :: u_new, v_new
        
        integer :: i, j, ip, ic, t_num_u_cp, t_num_v_cp
        real*8, dimension(num_u_cp) :: u_cp_u, v_cp_u
        real*8, dimension(num_v_cp) :: u_cp_v, v_cp_v
    
        do ip = 1, num_points
            
            do i = 1, num_u_cp
                do j = 1, num_v_cp
                    ic = (i-1)*num_v_cp + j
                    u_cp_v(j) = u_cp(ic)
                    v_cp_v(j) = v_cp(ic)
                end do
                
                t_num_v_cp = num_v_cp
                do
                    if (t_num_v_cp .lt. 1) exit
                    do j = 1, t_num_v_cp-1
                        u_cp_v(j) = (1.d0-v(ip))*u_cp_v(j) + v(ip)*u_cp_v(j+1)
                        v_cp_v(j) = (1.d0-v(ip))*v_cp_v(j) + v(ip)*v_cp_v(j+1)
                    end do
                    t_num_v_cp = t_num_v_cp - 1
                end do
                
                u_cp_u(i) = u_cp_v(1)
                v_cp_u(i) = v_cp_v(1)
            end do
            
            t_num_u_cp = num_u_cp
            do
                if (t_num_u_cp .lt. 1) exit
                do i = 1, t_num_u_cp-1
                    u_cp_u(i) = (1.d0-u(ip))*u_cp_u(i) + u(ip)*u_cp_u(i+1)
                    v_cp_u(i) = (1.d0-u(ip))*v_cp_u(i) + u(ip)*v_cp_u(i+1)
                end do
                t_num_u_cp = t_num_u_cp-1
            end do
            
            u_new(ip) = u_cp_u(1)
            v_new(ip) = v_cp_u(1)
            
        end do
    
    
    end subroutine uv_casteljau
    
    
    !******************************************************************
    !***** Binomial function
    !******************************************************************
    function binomial(n, t) result(res)
        integer, intent(in) :: n, t
        
        integer :: i, res, k

        k = t

        if (k .gt. n - k) k = n - k
        
        res = 1
        do i = 0, k-1
            res = res * (n - i)
            res = res / (i + 1)
        end do 
        
    end function binomial


end module bezier_helper_m

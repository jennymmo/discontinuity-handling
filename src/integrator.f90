module integrator_module
    use interpolator_module, only: interpolator
    use parameters, only: WP, DP, DI

    implicit none
    private
    public :: integrate_normal
    public :: integrate_discontinuity_handling
    public :: rk4, rk3_kutta, rk3_heun, rk2, rk1

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Bisection and Hermite to find time of boundary crossing !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function bisection(X0, X1, ks, boundary, h) result(theta)
        ! Find the dimensionless time theta that corresponds to the
        ! boundary crossing
        ! crossing happens for t_i <= t <= t_i+1
        ! theta = 0 means t = t_i
        ! theta = 1 means t = t_i+1

        implicit none
        real(WP), intent(in) :: X0, X1, boundary, h
        real(WP), dimension(:), intent(in) :: ks
        real(WP) :: theta

        ! Local variables
        real(WP) :: a, b, c, fa, fb, fc, tol, X_

        ! Initialize
        tol = 1e-9
        theta = 1e12_WP
        a = 0.0_WP
        b = 2.0_WP
        fa = X0 - boundary
        fb = X1 - boundary
        ! Initial checks
        if (fa .eq. 0.0) then
            theta = 0
        else if (fb .eq. 0.0) then
            theta = 1
        end if

        if (theta .gt. 1) then
            do while ((b-a) .gt. tol)
                ! Find midpoint
                c = (a+b)/2
                ! Note: this call to Hermite assumes that ks(1) is the derivative at point X0
                ! and that ks(size(ks)) is the derivative at point X1
                ! This is automatically the case for FSAL integrators,
                ! but must be set up manually outside this subroutine for other integrators
                X_ = hermite(X0, X1, ks(1), ks(size(ks)), h, c)
                fc = X_ - boundary
                ! If fa and fc have opposite sign, the boundary was crossed between them
                if (fa*fc .lt. 0) then
                    b = c
                ! If not, then the boundary must have been crossed between b and c
                else
                    a = c
                end if
            end do
            c = (a+b)/2
            theta = c
        end if
    end function bisection

    function hermite(X0, X1, k0, k1, h, t) result(X)
        implicit none
        real(WP), intent(in) :: X0, X1, k0, k1, h, t
        real(WP) :: X
        X = (1-t)*X0 + t*X1 + t*(t-1)*( (1-2*t)*(X1-X0) + (t-1)*h*k0 + t*h*k1)
    end function hermite


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Helper routines to find locate cell boundaries etc. !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_cell_indices(X, resolution, Xref, cell_indices)
        ! Determines which cell in the spatial grid the particle is currently in.
        ! Takes in the position of the particle, the resolution of the grid (in m),
        ! And a reference position. The reference position is assumed to be the 
        ! coordinates of  the point that would be at the origin if the coordinate  
        ! system was shifted so that the entire particle trajectory is in the first
        ! quadrant. For example, if X = {[x,y] : x,y > 0} then Xref could be [0,0]. 
        ! The reference position must be at a grid corner.
        implicit none
        ! IO
        real(WP), intent(in),  dimension(:)          :: X, Xref
        integer,  intent(in)                         :: resolution
        integer, dimension(size(X)), intent(out) :: cell_indices

        ! Calculate which cell X is in
        ! Adding 1 to get 1-indexing to be consistent with fortran convention,
        ! although this probably doesn't matter in practice since this is only used
        ! to check if a particle has moved to a different cell
        cell_indices = floor( (X - Xref) / real(resolution) )
    end subroutine get_cell_indices

    subroutine get_num_crossed_boundaries(X0, X1, Xref, resolution, crossed)
        implicit none
        ! IO
        real(WP),  dimension(:),        intent(in)  :: X0, X1, Xref
        integer,                        intent(in)  :: resolution
        integer,   dimension(size(X0)), intent(out) :: crossed
        ! Local
        integer, dimension(size(X0)) :: cell, next_cell

        ! Find out which grid cell the particle is in at this step ...
        call get_cell_indices(X0, resolution, Xref, cell)
        ! ... And after the next step
        call get_cell_indices(X1, resolution, Xref, next_cell)
        crossed = next_cell - cell
    end subroutine get_num_crossed_boundaries

    subroutine get_crossed_boundaries(X0, X1, Xref, resolution, boundaries)
        implicit none
        ! IO
        real(WP),  dimension(:),        intent(in)  :: X0, X1, Xref
        integer,                        intent(in)  :: resolution
        real(WP), dimension(:,:), allocatable, intent(out) :: boundaries
        integer,   dimension(size(X0)) :: crossed
        ! Local
        integer :: i, j, Nbounds
        integer, dimension(size(X0)) :: cell, next_cell

        ! Find out which grid cell the particle is in at this step ...
        call get_cell_indices(X0, resolution, Xref, cell)
        ! ... And after the next step
        call get_cell_indices(X1, resolution, Xref, next_cell)
        crossed = next_cell - cell

        ! Allocate space for the required number of crossings
        Nbounds = maxval(abs(crossed))
        allocate(boundaries(size(X0), Nbounds))
        boundaries = 0.0

        do i = 1, size(X0)
            if (crossed(i) .ne. 0) then
                do j = 1, abs(crossed(i))
                    if (crossed(i) > 0) then
                        boundaries(i,j) = Xref(i) + resolution*(cell(i) + j)
                    else
                        boundaries(i,j) = Xref(i) + resolution*(cell(i) - (j-1))
                    endif
                enddo
            endif
        end do
    end subroutine get_crossed_boundaries

    subroutine get_first_crossing(X, X5, ks, h, crossed, boundaries, theta, boundary, direction)
        implicit none
        ! IO
        real(WP), dimension(:),    intent(in) :: X, X5
        integer,  dimension(:),    intent(in) :: crossed
        real(WP), dimension(:,:),  intent(in) :: ks, boundaries
        real(WP), intent(in) :: h
        ! Output:
        ! Direction in which the earliest boundary is crossed
        ! Time to crossing, and position of boundary
        integer, intent(out)  :: direction
        real(WP), intent(out) :: theta, boundary

        ! Internal variables to store thetas
        ! (one for each direction)
        real(WP), dimension(:,:), allocatable :: thetas
        integer :: n_boundary, i
        integer, dimension(2) :: loc
        ! Allocate maximal needed size
        allocate(thetas(size(X), 1))
        ! Initialise to large value, as we will later look for minimum
        thetas = 1e12_WP

        ! Find time to first crossed boundary in any direction
        do i = 1, size(X)
            if (crossed(i) .ne. 0) then
                thetas(i,1) = bisection(X(i), X5(i), ks(:,i), boundaries(i,1), h)
            end if
        end do

        ! Find value and location of smallest theta
        theta = minval(thetas)
        loc = minloc(thetas)
        ! Find direction and boundary number from location in thetas
        direction = loc(1)
        n_boundary = loc(2)
        ! Find the boundary in the given direction
        boundary = boundaries(direction, n_boundary)

    end subroutine get_first_crossing


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Subroutines to integrate ODEs   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine integrate_normal(X, t0, tmax, h0, f, method, Nsteps)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(t=0), given by X at input.
        ! Solution is found by repeatedly calling a fixed-step
        ! Runge-Kutta method, with timestep h0.
        ! The routine adjusts the last timestep to stop exactly at tmax,
        ! and returns the last position.
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(in)                    :: t0, tmax, h0
        type(interpolator), intent(inout)       :: f
        integer(DP), intent(out)                :: Nsteps
        interface
            subroutine method(X, t, h, f, ks)
                import :: WP, interpolator
                real(WP), intent(inout), dimension(:)   :: X
                real(WP), intent(inout)                 :: t
                real(WP), intent(in)                    :: h
                type(interpolator), intent(inout)       :: f
                real(WP), allocatable,  intent(out),   dimension(:,:) :: ks
            end subroutine method
        end interface
        ! local variables
        real(WP) :: t, h
        real(WP), dimension(:,:), allocatable :: ks

        t = t0
        h = h0
        Nsteps = 0
        ! Loop over timesteps until t == tmax.
        do while (t < tmax)
            ! If remaining time until tmax is smaller than timestep,
            ! adjust h to stop exactly at tmax.
            h = min(h, tmax - t)
            call method(X, t, h, f, ks)
            ! Count number of steps, for convenient comparison to other methods
            Nsteps = Nsteps + 1
        end do
    end subroutine

    subroutine integrate_discontinuity_handling(X, t0, tmax, h0, f, method, Xref, resolution, Nsteps, callfactor)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x,t),
        ! with initial value, X(t=t0), given by X at input.
        ! Solution is found by repeatedly calling a fixed-step 
        ! Runge-Kutta method.
        ! The routine adjusts the last timestep to stop exatly at tmax,
        ! and returns the last position.
        ! Additionally, the routine keeps track of the locations at which 
        ! there are discontinuities in the data. With the help of a dense
        ! output procedure the integration is stopped and restarted at 
        ! these locations.
        implicit none
        ! IO
        real(WP),           intent(inout), dimension(:)       :: X
        real(WP),           intent(in)                        :: t0, tmax, h0
        type(interpolator), intent(inout)                     :: f
        real(WP),           intent(in),    dimension(size(X)) :: Xref
        integer,            intent(in)                        :: resolution
        real(DP),           intent(out)                       :: Nsteps
        real(DP),           intent(in)                        :: callfactor
        interface
            subroutine method(X, t, h, f, ks)
                import :: interpolator, WP
                real(WP),               intent(inout), dimension(:) :: X
                real(WP),               intent(inout)               :: t
                real(WP),               intent(in)                  :: h
                type(interpolator),     intent(inout)               :: f
                real(WP), allocatable,  intent(out),   dimension(:,:) :: ks
            end subroutine method
        end interface

        ! Local variables:
        real(WP) :: t, h

        t = t0
        Nsteps = 0._DP
        ! Loop over timesteps until t == tmax.
        do while ( t < tmax)
            ! If remaining time until tmax is smaller than timestep,
            ! adjust h to stop exactly at tmax.
            h = min(h0, tmax-t)
            ! Here we use a wrapper routine that will call the integrator one or more times
            ! (depending on boundary crossings)
            ! in order to complete a step of duration h
            call onestep_fixed(X, t, h, f, method, Xref, resolution, Nsteps, callfactor)
        enddo
    end subroutine integrate_discontinuity_handling

    recursive subroutine onestep_fixed(X, t, h, f, method, Xref, resolution, Nsteps, callfactor)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x,t),
        ! with initial value, X(t=t0), given by X at input.
        ! Solution is found by repeatedly calling a fixed-step 
        ! Runge-Kutta method.
        ! The routine adjusts the last timestep to stop exatly at tmax,
        ! and returns the last position.
        ! Additionally, the routine keeps track of the locations at which 
        ! there are discontinuities in the data. With the help of a dense
        ! output procedure the integration is stopped and restarted at 
        ! these locations.
        implicit none
        ! IO
        real(WP),           intent(inout), dimension(:)       :: X
        real(WP),           intent(inout)                     :: t
        real(WP),           intent(in)                        :: h
        type(interpolator), intent(inout)                     :: f
        real(WP),           intent(in),    dimension(size(X)) :: Xref
        integer,            intent(in)                        :: resolution
        real(DP),           intent(inout)                     :: Nsteps
        real(DP),           intent(in)                        :: callfactor
        interface
            subroutine method(X, t, h, f, ks)
                import :: interpolator, WP
                real(WP),               intent(inout), dimension(:)   :: X
                real(WP),               intent(inout)                 :: t
                real(WP),               intent(in)                    :: h
                type(interpolator),     intent(inout)                 :: f
                real(WP), allocatable,  intent(out),   dimension(:,:) :: ks
            end subroutine method
        end interface

        ! local variables
        real(WP), dimension(size(X)) :: X_
        integer, dimension(size(X)) :: crossed
        real(WP), dimension(:,:), allocatable :: ks, boundaries
        real(WP) :: t_, h_, theta, boundary
        integer :: direction
        !real :: speed
        X_ = X
        t_ = t

        ! Make one step with the chosen integrator
        call method(X_, t_, h, f, ks)
        ! Determine if the particle is still in the same cell
        call get_num_crossed_boundaries(X, X_, Xref, resolution, crossed)

        if (any(crossed .ne. 0)) then
            ! One or more cell border is crossed - handle discontinuities
            ! Update derivative at endpoint, needed for Hermite interpolation
            ! and bisection that happens below (including in the get_first_crossing subroutine)
            ks(size(ks, 1),:) = f%eval(X_, t_)
            ! Add a fractional step to account for the one extra eval
            Nsteps = Nsteps + callfactor

            ! Using the bisection method, the theta corresponding to the
            ! first border crossing is found.
            ! See e.g. chapter II.6 of 'Solving Ordinary Differential
            ! Equations I Nonstiff problems' by Hairer et al. (2008)
            call get_crossed_boundaries(X, X_, Xref, resolution, boundaries)
            call get_first_crossing(X, X_, ks, h, crossed, boundaries, theta, boundary, direction)

            !print*, 'theta before = ', theta
            theta = bisection(X(direction), X_(direction), ks(:,direction), boundary, h)
            !print*, 'theta after  = ', theta

            if (theta < 1e-3) then
                ! We will never hit the boundary _exactly_, so if we are sufficiently close, just move on.
                X = X_
                t = t_
                Nsteps = Nsteps + 1._DP
            else
                ! Adjust the timestep and take a trial step, deliberately too short
                h_ = theta * h * 0.9
                X_ = X
                t_ = t
                call method(X_, t_, h_, f, ks)
                ! Update derivative at endpoint of trial step, needed below
                ks(size(ks, 1),:) = f%eval(X_, t_)
                ! Extrapolate to find the duration of a step that will hit the boundary,
                ! and take that step
                theta = bisection(X(direction), X_(direction), ks(:, direction), boundary, h_)
                h_ = theta * h_
                X_ = X
                t_ = t
                call method(X_, t_, h_, f, ks)

                ! Add two steps to account for the two calls to method,
                ! plus fractional step for the one extra eval
                Nsteps = Nsteps + 2._DP + callfactor

                ! Then take another step to complete one timestep
                call onestep_fixed(X_, t_, h-h_, f, method, Xref, resolution, Nsteps, callfactor)
                X = X_
                t = t_
            endif
        else
            ! No cell borders are crossed, continue integration
            X = X_
            t = t_
            ! Increment step count by one
            Nsteps = Nsteps + 1._DP
        end if
    end subroutine onestep_fixed


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Fixed step Runge-Kutta methods !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine rk4( X, t, h, f, ks)
        ! Make one step with the classic 4th-order Runge-Kutta method.
        ! Calculates next position using timestep h.
        ! See, e.g., Griffiths (2010, p.131)
        implicit none
        ! IO
        real(WP),               intent(inout), dimension(:)   :: X
        real(WP),               intent(inout)                 :: t
        real(WP),               intent(in)                    :: h
        type(interpolator),     intent(inout)                 :: f
        real(WP), allocatable,  intent(out),   dimension(:,:) :: ks

        allocate(ks(4, size(X)))
        ! Evaluations of f(x,t)
        ks(1,:) = f%eval(X,          t      )
        ks(2,:) = f%eval(X + ks(1,:)*h/2, t + h/2)
        ks(3,:) = f%eval(X + ks(2,:)*h/2, t + h/2)
        ks(4,:) = f%eval(X + ks(3,:)*h,   t + h  )

        ! Next position
        X = X + h*(ks(1,:) +2*ks(2,:) + 2*ks(3,:) + ks(4,:)) / 6
        t = t + h
    end subroutine

    subroutine rk3_kutta( X, t, h, f, ks)
        ! Kutta's third-order method
        implicit none
        ! IO
        real(WP),               intent(inout), dimension(:)   :: X
        real(WP),               intent(inout)                 :: t
        real(WP),               intent(in)                    :: h
        type(interpolator),     intent(inout)                 :: f
        real(WP), allocatable,  intent(out),   dimension(:,:) :: ks

        allocate(ks(3, size(X)))
        ! Evaluations of f(x,t)
        ks(1,:) = f%eval(X,                            t      )
        ks(2,:) = f%eval(X + ks(1,:)*h/2,              t + h/2)
        ks(3,:) = f%eval(X - ks(1,:)*h +  ks(2,:)*2*h, t + h  )

        ! Next position
        X = X + h*(ks(1,:) + 4*ks(2,:) + ks(3,:)) / 6
        t = t + h
    end subroutine

    subroutine rk3_heun( X, t, h, f, ks)
        ! Heun's third-order method
        implicit none
        ! IO
        real(WP),               intent(inout), dimension(:)   :: X
        real(WP),               intent(inout)                 :: t
        real(WP),               intent(in)                    :: h
        type(interpolator),     intent(inout)                 :: f
        real(WP), allocatable,  intent(out),   dimension(:,:) :: ks

        allocate(ks(3, size(X)))
        ! Evaluations of f(x,t)
        ks(1,:) = f%eval(X,                            t          )
        ks(2,:) = f%eval(X + ks(1,:)*h/3,              t +   h/3  )
        ks(3,:) = f%eval(X + ks(2,:)*2*h/3, t + 2*h/3  )

        ! Next position
        X = X + h*(ks(1,:) + 3*ks(3,:)) / 4
        t = t + h
    end subroutine

    subroutine rk2( X, t, h, f, ks)
        ! Explicit trapezoid second-order method
        implicit none
        ! IO
        real(WP),               intent(inout), dimension(:)   :: X
        real(WP),               intent(inout)                 :: t
        real(WP),               intent(in)                    :: h
        type(interpolator),     intent(inout)                 :: f
        real(WP), allocatable,  intent(out),   dimension(:,:) :: ks

        allocate(ks(2, size(X)))
        ! Evaluations of f(x,t)
        ks(1,:) = f%eval(X,             t    )
        ks(2,:) = f%eval(X + ks(1,:)*h, t + h)

        ! Next position
        X = X + h*(ks(1,:) + ks(2,:)) / 2
        t = t + h
    end subroutine

    subroutine rk1( X, t, h, f, ks)
        ! Explicit trapezoid second-order method
        implicit none
        ! IO
        real(WP),               intent(inout), dimension(:)   :: X
        real(WP),               intent(inout)                 :: t
        real(WP),               intent(in)                    :: h
        type(interpolator),     intent(inout)                 :: f
        real(WP), allocatable,  intent(out),   dimension(:,:) :: ks

        allocate(ks(1, size(X)))
        ! Evaluations of f(x,t)
        ks(1,:) = f%eval(X, t)

        ! Next position
        X = X + h*ks(1,:)
        t = t + h
    end subroutine

end module

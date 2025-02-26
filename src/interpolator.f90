module interpolator_module

use bspline_module, only: bspline_3d
use parameters, only: SP, WP

implicit none

private
public :: interpolator


type interpolator
    type(bspline_3d) :: fvx, fvy
    real(WP), dimension(:), allocatable :: xc, yc, tc
    contains
        private
        procedure, public  :: eval => evaluate
        procedure, public  :: init => initialize
        procedure, public  :: destroy => destroy
end type interpolator

contains

    subroutine destroy(this)
        implicit none
        class(interpolator), intent(inout)     :: this
        call this%fvx%destroy()
        call this%fvy%destroy()
    end subroutine

    subroutine initialize(this, xc, yc, tc, gvx, gvy, order)
        implicit none
        class(interpolator), intent(inout)     :: this
        real(WP), intent(in), dimension(:)     :: xc, yc, tc
        real(WP), intent(in), dimension(:,:,:) :: gvx, gvy
        integer,  intent(in)                   :: order
        integer                                :: iflag
        ! allocate space for, and store, coordinate arrays
        allocate(this % xc(size(xc)))
        allocate(this % yc(size(yc)))
        allocate(this % tc(size(tc)))
        this % xc = xc
        this % yc = yc
        this % tc = tc
        ! initialise interpolation objects from bspline-fortran
        call this%fvx%initialize(xc, yc, tc, gvx, order, order, order, iflag, extrap = .true.)
        if (iflag /= 0) then
            print*, this%fvx%status_message(iflag)
        endif
        call this%fvy%initialize(xc, yc, tc, gvy, order, order, order, iflag, extrap = .true.)
        if (iflag /= 0) then
            print*, this%fvy%status_message(iflag)
        endif
    end subroutine

    function evaluate(this, X, t) result(V)
        implicit none
        !!!! inputs !!!!
        class(interpolator),    intent(inout) :: this
        real(WP), dimension(2), intent(in)    :: X
        real(WP),               intent(in)    :: t

        !!!! output !!!!
        ! velocity vector to return
        real(WP), dimension(2)                :: V

        !!!! local variables !!!!
        ! status flag
        integer                               :: iflag
        ! order of derivative to return (0 returns the function value itself)
        integer, parameter                    :: d = 0

        ! Check that we are within bounds
        if (.not. ((this%xc(1) <= X(1)) .and. (X(1) <= this%xc(size(this%xc))))) then
            print*, 'x-position is outside region covered by data'
        endif
        if (.not. ((this%yc(1) <= X(2)) .and. (X(2) <= this%yc(size(this%yc))))) then
            print*, 'y-position is outside region covered by data'
        endif
        if (.not. ((this%tc(1) <= t) .and. (t <= this%tc(size(this%tc))))) then
            print*, 'time is outside interval covered by data'
        endif

        ! Interpolate x-component of velocity
        call this%fvx%evaluate(X(1), X(2), t, d, d, d, V(1), iflag)
        if (iflag /= 0) then
            print*, this%fvx%status_message(iflag)
            print*, 'X = ', X, 't = ', t
        endif

        ! Interpolate y-component of velocity
        call this%fvy%evaluate(X(1), X(2), t, d, d, d, V(2), iflag)
        if (iflag /= 0) then
            print*, this%fvx%status_message(iflag)
            print*, 'X = ', X, 't = ', t
        endif
    end function

end module interpolator_module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8

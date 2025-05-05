program run

use parameters,             only: SP, DP, WP, input_folder
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use experiment_module,      only: experiment_discontinuity_handling
use interpolator_module,    only: interpolator
use integrator_module,      only: rk4, rk3_kutta, rk3_heun, rk2, rk1
use output_module,          only: get_output_filename

implicit none
! Arrays for particles
real(wp), dimension(:,:), allocatable :: X0      ! two-component vectors
! Coordinate arrays
real(WP), dimension(:),     allocatable :: xc, yc, tc
! Velocity x and y components
real(WP), dimension(:,:,:), allocatable :: u, v
! Derived type to evaluate interpolated current data
type(interpolator) :: f
! Time
real(wp) :: t0, tmax
! List of timesteps to test
real(WP), dimension(10), parameter :: timesteps = (/ 60, 90, 120, 180, 300, 600, 900, 1200, 1800, 3600 /)
! Order of interpolation (comes from command line argument)
integer :: order
! Loop variables
integer :: i
! Variables for filenames
character(len=256) :: currentdata_filename
character(len=256) :: initial_position_filename
character(len=256) :: output_filename

! Name identifying the dataset
character(len=16) :: dataset_name
integer :: resolution ! Resolution of the dataset in meters
real(WP), dimension(2) :: xref ! Position of bottom left corner of the dataset

! Meta-variables for command line arguments
integer :: num_args, length, info
character(16) :: arg

real(DP) :: hackfactor ! Used to count number of evaluations in the event detection scheme



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prepare current data, interpolator, initial positions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions, hardcoded for nordic4
currentdata_filename = trim(input_folder) // 'nordic4_larger.nc'
initial_position_filename   = trim(input_folder) // 'initial_positions_nordic4_backwards_discontinuity_handling.txt'

! Set resolution and order of interpolation
dataset_name = 'nordic4'
resolution = 4000
order = 2

! Read initial particle positions from file
print*, 'Reading initial positions from file'
call read_initial_positions(initial_position_filename, X0)
! Get data from NetCDF file
print*, 'Reading current data from file'
call get_current(currentdata_filename, u, v, xc, yc, tc)

! Set reference position for detection of boundary crossings
xref = (/ xc(1), yc(1) /) ! Position of bottom left corner of the dataset

! Duration of integration (taken from time coordinates of data)
t0   = tc(6)
tmax = tc(6 + 25*24)


!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run simulation !!!!
!!!!!!!!!!!!!!!!!!!!!!!!

! Create interpolator of desired order from discrete data
print*, 'Creating interpolator'
call f%init(xc, yc, tc, u, v, order)

! Run simulations with rk4
output_filename = get_output_filename('backward_run_discontinuity_handling', 'rk4', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
print*, 'Xref = ', Xref
hackfactor=1._DP / 4._DP
call experiment_discontinuity_handling(X0, t0, tmax, timesteps, f, rk4, output_filename, Xref, resolution, hackfactor, .true.)


! Deallocate interpolator
call f%destroy()

deallocate(u)
deallocate(v)

end program

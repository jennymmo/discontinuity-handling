program run

use parameters,             only: SP, DP, WP, input_folder, timesteps
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use experiment_module,      only: experiment_normal
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
! Order of interpolation (comes from command line argument)
integer :: order

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

! Check for presence of required command line arguments
num_args = command_argument_count()
if (num_args .ne. 2) then
    print*, "run.f90"
    print*, "Usage:"
    print*, "run dataset_name order"
    print*, "    dataset_name : String, name of input dataset [norkyst800, nordic4, arctic20]"
    print*, "    order        : Integer, order of bspline interpolation (order = degree + 1)"
    print*, " "
    print*, "Incorrect number of command arguments received, stopping"
    stop
else
    ! Read dataset name
    call get_command_argument(1, arg, length, info)
    read(arg, *) dataset_name
    ! Read interpolation order
    call get_command_argument(2, arg, length, info)
    read(arg, *) order
    ! Print arguments for checking
    print*, "Running simulation with these parameters:"
    print*, "    dataset_name = ", trim(dataset_name)
    print*, "    order        = ", order
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prepare current data, interpolator, initial positions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions
currentdata_filename = trim(input_folder) // trim(dataset_name) // '.nc'
initial_position_filename   = trim(input_folder) // 'initial_positions_' // trim(dataset_name) // '.txt'

! Set resolution
if (trim(dataset_name) == "norkyst800") then
    resolution = 800
else if (trim(dataset_name) == "nordic4") then
    resolution = 4000
else if (trim(dataset_name) == "arctic20") then
    resolution = 20000
else
    print*, "Invalid dataset name specified, stopping"
    stop
endif


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
tmax = tc(6 + 72)


!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run simulation !!!!
!!!!!!!!!!!!!!!!!!!!!!!!

! Create interpolator of desired order from discrete data
print*, 'Creating interpolator'
call f%init(xc, yc, tc, u, v, order)

! Run simulations with rk4
output_filename = get_output_filename('runs_normal', 'rk4', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
call experiment_normal(X0, t0, tmax, timesteps, f, rk4, output_filename)

! Run simulations with rk3 (Heun)
output_filename = get_output_filename('runs_normal', 'rk3_heun', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
call experiment_normal(X0, t0, tmax, timesteps, f, rk3_heun, output_filename)

! Run simulations with rk3 (Kutta)
output_filename = get_output_filename('runs_normal', 'rk3_kutta', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
call experiment_normal(X0, t0, tmax, timesteps, f, rk3_kutta, output_filename)

! Run simulations with rk2
output_filename = get_output_filename('runs_normal', 'rk2', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
call experiment_normal(X0, t0, tmax, timesteps, f, rk2, output_filename)

! Run simulations with rk1 (Euler's method)
output_filename = get_output_filename('runs_normal', 'rk1', dataset_name, order)
print*, 'Running simulation with timesteps ', timesteps
print*, 'Storing results to ', output_filename
call experiment_normal(X0, t0, tmax, timesteps, f, rk1, output_filename)

! Deallocate interpolator
call f%destroy()

deallocate(u)
deallocate(v)

end program

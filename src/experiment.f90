module experiment_module

    use hdf5
    use h5lt
    use parameters,             only: SP, DP, WP, output_folder
    use interpolator_module,    only: interpolator
    use integrator_module,      only: integrate_normal
    use integrator_module,      only: integrate_discontinuity_handling
    use output_module,          only: write_to_hdf5, create_hdf5_file, close_hdf5_file

    implicit none
    private
    public :: experiment_normal
    public :: experiment_discontinuity_handling

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Subroutine to run experiment with fixed-step Runge-Kutta method !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine experiment_normal(X0, t0, tmax, timesteps, f, method, output_filename)
        ! This subroutine takes the initial positions at time t0, for a number,
        ! Np, of particles, along with a list of timesteps, an interpolator
        ! object to evaluate the velocity field, and a numerical ODE integrator.
        ! It then obtains the position at time tmax for each
        ! particle, writes to file, and repeats for each timestep.
        implicit none
        ! Input variables
        real(WP), dimension(:,:), intent(in)    :: X0
        real(WP), dimension(:),   intent(in)    :: timesteps
        real(WP),                 intent(in)    :: t0, tmax
        character(len=*),         intent(in)    :: output_filename
        type(interpolator),       intent(inout) :: f
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
        real(WP), dimension(:,:), allocatable   :: X
        integer(hid_t)                          :: file_id
        integer                                 :: n, idt, Np
        integer(DP), dimension(size(X0,2))      :: Nsteps
        real(WP)                                :: h, tic, toc

        ! Assuming that X0 has shape (2, Np) or (3, Np)
        ! depending on two or three dimensions
        Np = size(X0, 2)
        allocate(X( size(X0,1), size(X0,2) ))

        ! Create file for output
        call create_hdf5_file(output_filename, file_id)
        ! Open text file for writing number of calls and runtime
        open(42, file = trim(output_filename) // '.txt')
        write(42,*) '# Timestep               Runtime                       Naccepted         Nrejected'

        ! Scan through timesteps
        do idt = 1, size(timesteps)
            h = timesteps(idt)
            print*, 'Running simulation with timestep ', h
            ! Reset initial positions
            X = X0
            ! Measure computational time
            call cpu_time(tic)
            ! Transport each particle from time t0 to tmax
            do n = 1, Np
                call integrate_normal(X(:,n), t0, tmax, h, f, method, Nsteps(n))
            end do
            ! Measure computational time
            call cpu_time(toc)
            ! Write end positions to hdf5 file
            call write_to_hdf5(file_id, X, timesteps(idt))
            ! Print information on timing and number of steps:
            ! timestep, runtime, accepted steps, rejected steps
            write(42,*) timesteps(idt), toc - tic, sum(Nsteps), 0
            flush(42)
        end do
        ! Clean up
        deallocate(X)
        call close_hdf5_file(file_id)
        close(42)
    end subroutine

    subroutine experiment_discontinuity_handling(X0, t0, tmax, timesteps, f, method, output_filename, Xref, resolution, callfactor)
        implicit none
        !Input variables 
        real(WP), dimension(:,:), intent(in) :: X0
        real(WP),                 intent(in) :: t0, tmax
        real(WP), dimension(:),   intent(in) :: timesteps
        character(len=*),         intent(in) :: output_filename
        real(WP), dimension(2),   intent(in) :: Xref
        integer,                  intent(in) :: resolution
        real(DP),                 intent(in) :: callfactor
        type(interpolator),       intent(inout) :: f 

        interface
            subroutine method(X, t, h, f, ks)
                import :: WP, interpolator
                real(WP), intent(inout), dimension(:) :: X 
                real(WP), intent(inout)               :: t 
                real(WP), intent(in)                  :: h 
                type(interpolator), intent(inout)     :: f 
                real(WP), allocatable,  intent(out),   dimension(:,:) :: ks
            end subroutine method 
        end interface 
        ! Local variables
        real(WP), dimension(:,:), allocatable :: X 
        integer(hid_t)                        :: file_id 
        integer                               :: n, idt, Np 
        real(DP), dimension(size(X0,2))    :: Nsteps 
        real(WP)                              :: h, tic, toc

        ! Assuming that X0 has shape (2, Np)
        Np = size(X0, 2)
        allocate(X( size(X0,1), size(X0,2)))

        ! Create file for output
        call create_hdf5_file(output_filename, file_id)
        ! Open text file for writing number of calls and runtime
        open(42, file = trim(output_filename) // '.txt')
        write(42,*) '# Timestep               Runtime                       Naccepted         Nrejected'

        ! Scan through timesteps
        do idt = 1, size(timesteps)
            h = timesteps(idt)
            print*, 'Running simulation with timestep ', h
            ! Reset initial positions
            X = X0
            ! Measure computational time
            call cpu_time(tic)
            ! Transport each particle from time t0 to tmax
            do n = 1, Np
                call integrate_discontinuity_handling(X(:,n), t0, tmax, h, f, method, Xref, resolution, Nsteps(n), callfactor)
            end do
            ! Measure computational time
            call cpu_time(toc)
            ! Write end positions to hdf5 file
            call write_to_hdf5(file_id, X, timesteps(idt))
            ! Print information on timing and number of steps:
            ! timestep, runtime, accepted steps, rejected steps
            write(42,*) timesteps(idt), toc - tic, sum(Nsteps), 0
            flush(42)
        end do
        ! Clean up
        deallocate(X)
        call close_hdf5_file(file_id)
        close(42)
    end subroutine

    end module

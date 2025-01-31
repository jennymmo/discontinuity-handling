module parameters
    implicit none
    ! Precision parameters
    integer,  parameter :: SP = kind(1.0E0)
    integer,  parameter :: DP = kind(1.0D0)
    integer,  parameter :: WP = DP
    integer,  parameter :: WI = kind(1)
    ! Double integer
    integer,  parameter :: DI = SELECTED_INT_KIND(17)

    ! Mathematical constants
    real(DP), parameter :: PI = 3.1415926535897932384626433_DP

    ! Folders for input and output
    character(len =  11), parameter :: input_folder = '../data/'
    character(len =  12), parameter :: output_folder = '../results/'

    ! List of timesteps to test
    real(WP), dimension(14), parameter :: timesteps = (/10, 20, 30, 60, 90, 120, 180, 240, 300, 600, 900, 1200, 1800, 3600 /)

end module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8

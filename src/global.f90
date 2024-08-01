Module global

    implicit none
    double precision, parameter :: PI = acos(-1.d0)

    interface
      function intensity(x_int)
        double precision :: intensity
        double precision, intent(in) :: x_int
      end function intensity
    end interface

    ! Computational domain
    integer :: ny, nz, nym, nzm
    double precision, allocatable, dimension(:), public :: y, z
	double precision :: dy, dz
	
	! Computational time
	integer :: ntime
	double precision :: dt

    ! SEM parameter
    double precision :: sigma   !< eddy size
    double precision :: eps     !<
    double precision :: neddy   !< the number of symthetic eddies

    double precision :: v_b     !< volume of the virtual box
    double precision :: ystart, yend, zstart, zend

	integer :: int_x, int_y, int_z
	double precision, allocatable, dimension(:), public :: eddy_len, expos, eypos, ezpos
	double precision, allocatable, dimension(:), public :: euint, evint, ewint, etint
	double precision, allocatable, dimension(:,:), public :: u_inlet, v_inlet, w_inlet, t_inlet
	double precision, allocatable, dimension(:,:), public :: u_comb, v_comb, w_comb, t_comb

End module global

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
    integer :: ny, nz
    double precision, allocatable, dimension(:), public :: y, z

    ! SEM parameter
    double precision :: sigma   !< eddy size
    double precision :: eps     !<
    double precision :: Neddy   !< the number of symthetic eddies

    double precision :: V_b     !< volume of the virtual box
    double precision :: ystart, yend, zstart, zend


End module global

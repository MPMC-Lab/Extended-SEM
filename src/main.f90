Program main

    use global

    implicit none

    call read_input

    write(*,*) ny, nz

End program main

Subroutine read_input

    use global

    ny=101
    nz=201

End subroutine read_input

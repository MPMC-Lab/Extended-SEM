Program main
    use global
    implicit none
	integer :: it

    call inputs
	call arrs_allocate
	call set_parameters
	
	do it=1,ntime
		call xsem_process
	end do
	
	call arrs_deallocate

    write(*,*) ny, nz

End program main


Subroutine inputs
    use global
	implicit none

    ny=101
    nz=201
	nym=ny-1
	nzm=nz-1
	
	ntime=1000
	dt=0.01

	sigma = 0.5
	eps   = 1e-05
	neddy = 1000
	
End subroutine inputs


Subroutine arrs_allocate
	use global
	implicit none
	integer :: i, j, k
	
	allocate(y(0:ny), z(0:nz))
	allocate(eddy_len(0:neddy), expos(0:neddy), eypos(0:neddy), ezpos(0:neddy))
	allocate(u_inlet(0:ny,0:nz), v_inlet(0:ny,0:nz), w_inlet(0:ny,0:nz), t_inlet(0:ny,0:nz))
	allocate(euint(0:ny,0:nz), evint(0:ny,0:nz), ewint(0:ny,0:nz), etint(0:ny,0:nz))
	allocate(u_comb(0:ny,0:nz), v_comb(0:ny,0:nz), w_comb(0:ny,0:nz), t_comb(0:ny,0:nz))
	
	y(0) = 0.d0
	y(ny)= 5.d0
	dy= (y(0)-y(ny))/nym
	do j=1, nym
		y(j) = real(j-1)*dy
	end do
	
	z(0) = 0.d0
	z(nz)= 5.d0
	dz= (z(0)-z(nz))/nzm
	do k=1, nzm
		z(k) = real(k-1)*dz
	end do
End Subroutine arrs_allocate

Subroutine arrs_deallocate
	use global
	implicit none
	deallocate(y, z)
	deallocate(eddy_len, expos, eypos, ezpos)
	deallocate(u_inlet, v_inlet, w_inlet, t_inlet)
	deallocate(euint, evint, ewint, etint)
	deallocate(u_comb, v_comb, w_comb, t_comb)
End Subroutine arrs_deallocate


Subroutine set_parameters
    use global
	implicit none
	double precision :: Stmp

	ystart= y(1) - sigma
	yend  = y(ny) + sigma
	zstart= z(1) - sigma
	zend  = z(nz) + sigma
	
	Stmp= (yend-ystart)*(zend-zstart)
	v_b = Stmp*2*sigma
	
End subroutine set_parameters


Subroutine xsem_process
	use global
	implicit none
	integer :: it, i, j, k
	integer :: zero
	double precision :: tmp(1:3),tmp2(1:6)
	double precision :: x0, y0, z0

	if(ntime.EQ.1) then
		write(*,*) '  '
		write(*,*) 'Initialize eddy positions'
		write(*,*) '  '
		
		call random_seed
        call random_number(int_x)
        call random_number(int_y)
        call random_number(int_z)

		zero=int(0,4)
		
		do it = 1,neddy
			eddy_len(it) = sigma
			!eddy_len(it,1) = sigma
			!eddy_len(it,2) = sigma
			!eddy_len(it,3) = sigma
			
			call random_number(tmp)
			expos(it) = -sigma + 2*sigma*tmp(1)
			eypos(it) = ystart + (yend-ystart)*tmp(2)
			ezpos(it) = zstart + (zend-zstart)*tmp(3)
			
			euint(it) = intensity_det(int_x(it)-0.5)
			evint(it) = intensity_det(int_y(it)-0.5)
			ewint(it) = intensity_det(int_z(it)-0.5)
			etint(it) = intensity_det(int_z(it)-0.5)
		end do

	end if

	!>-----------------------<!
	!> Generate fluctuations <!
	!>-----------------------<!
	u_inlet= 0.d0
	v_inlet= 0.d0
	w_inlet= 0.d0
	t_inlet= 0.d0
	
	do k = 1,nz
	do j = 1,ny
	do it = 1,neddy

		x0 = (0.0  - expos(it))/eddy_len(it) 
		y0 = (y(j) - eypos(it))/eddy_len(it)
		z0 = (z(k) - ezpos(it))/eddy_len(it)

		!> Shape function
		if ( abs(x0).LE.1 .AND. abs(y0).LE.1 .AND. abs(z0).LE.1) THEN
		f = sqrt(1.5) * (1.- abs(x0)) * &
			sqrt(1.5) * (1.- abs(y0)) * &
			sqrt(1.5) * (1.- abs(z0))

		u_inlet(j,k) = u_inlet(j,k) + sqrt(v_b/eddy_len(it)**3) * exint(it)*f
		v_inlet(j,k) = v_inlet(j,k) + sqrt(v_b/eddy_len(it)**3) * eyint(it)*f
		w_inlet(j,k) = w_inlet(j,k) + sqrt(v_b/eddy_len(it)**3) * ezint(it)*f
		t_inlet(j,k) = t_inlet(j,k) + sqrt(v_b/eddy_len(it)**3) * etint(it)*f
		end if
	
	end do
	end do
	end do
	
	u_inlet(1:ny,1:nz) = u_inlet(1:ny,1:nz) / sqrt(real(neddy,8))
	v_inlet(1:ny,1:nz) = v_inlet(1:ny,1:nz) / sqrt(real(neddy,8))
	w_inlet(1:ny,1:nz) = w_inlet(1:ny,1:nz) / sqrt(real(neddy,8))
	t_inlet(1:ny,1:nz) = t_inlet(1:ny,1:nz) / sqrt(real(neddy,8))
	
	!>-------------------------------<!
	!> Combine mean and fluctuations <!
	!>-------------------------------<!
	u_comb = 0.d0
	v_comb = 0.d0
	w_comb = 0.d0
	t_comb = 0.d0
	
	do k = 1,n3
	do j = 1,n2
		a(1:4,1:4)     = 0.0
		r_loc(1:4,1:4) = 0.0
		u_ins(1:4,1)   = 0.0
		u_fluc(1:4,1)  = 0.0

		u_mean(1:4,1) = (/u_read(j,k),v_read(j,k),
		>                      w_read(j,k),t_read(j,k)/)
		u_tmp(1:4,1)  = (/u_inlet(j,k),v_inlet(j,k),   
		>                      w_inlet(j,k),t_inlet(j,k)/)
		r_loc(1,1:4) = (/rs(1,j,k),rs(4,j,k),rs(5,j,k),ths(2,j,k)/)
		r_loc(2,1:4) = (/rs(4,j,k),rs(2,j,k),rs(6,j,k),ths(3,j,k)/)
		r_loc(3,1:4) = (/rs(5,j,k),rs(6,j,k),rs(3,j,k),ths(4,j,k)/)
		r_loc(4,1:4) = (/ths(2,j,k),ths(3,j,k),ths(4,j,k),ths(1,j,k)/)

	call chol(a,r_loc,4,eps)
	! call mat_mul(a,u_tmp,u_fluc,4,1,4)
	do ii = 1,4
	do jj = 1,1
		do kk = 1,4
			u_fluc(ii,jj) = u_fluc(ii,jj) + a(ii,kk) * u_tmp(kk,jj)
		end do
	end do
	end do

	u_ins(1:4,1) = u_mean(1:4,1) + u_fluc(1:4,1)

	u_comb(j,k) = u_ins(1,1)
	v_comb(j,k) = u_ins(2,1)
	w_comb(j,k) = u_ins(3,1)
	t_comb(j,k) = u_ins(4,1)

	end do
	end do
	
	
	!>-------------------<!
	!> Convecting eddies <!
	!>-------------------<!
	call random_seed
	u_conv = sum(u_comb)/real(nym*nzm)
	v_conv = sum(v_comb)/real(nym*nzm)
	w_conv = sum(w_comb)/real(nym*nzm)

	!write(103,*) U_conv, V_conv, W_conv, DT
	!write(*,"(1A,e12.5)") 'XSEM CONVECTIVE VELOCITY U=', U_conv

	do it = 1,neddy

	!write(104,*) expos(it)*10,eypos(it),ezpos(it)

	expos(it) = expos(it) + u_conv*dt
	eypos(it) = eypos(it) + v_conv*dt
	ezpos(it) = ezpos(it) + w_conv*dt

		if ( (expos(it)-(-sigma))*(expos(it)-(sigma)) .gt. 0.0 .or. &
			 (eypos(it) - y_start)*(eypos(it) - y_end) .gt. 0.0 ) then

		eddy_len(it) = sigma
		!eddy_len(it,1) = sigma
		!eddy_len(it,2) = sigma
		!eddy_len(it,3) = sigma
		!expos(it) = - sigma

		call random_number(tmp2)
		expos(it) = 0.0 - sigma
		eypos(it) = y_start + (y_end-y_start)*tmp2(1)
		ezpos(it) = z_start + (z_end-z_start)*tmp2(2)

		exint(it) = intensity_det(tmp2(3)-0.5)
		eyint(it) = intensity_det(tmp2(4)-0.5)
		ezint(it) = intensity_det(tmp2(5)-0.5)
		etint(it) = intensity_det(tmp2(6)-0.5)

		end if

		!periodic boundary conditions
		if ( ezpos(it) .lt. z_start ) then
		tmp_z = z_start - ezpos(it)
		ezpos(it) = z_end - tmp_z
		end if

		if ( ezpos(it) .gt. z_end ) then
		tmp_z = ezpos(it) - z_end
		ezpos(it) = z_start + tmp_z
		end if

	END DO
	
	
	

End Subroutine xsem_process

Real*8 function intensity_det(x_int)
	implicit none
	real, intent(in) :: x_int

	if ( x_int .gt. 0.0 ) then
		intensity_det = 1
	else
		intensity_det = -1
	end if	

	return
end
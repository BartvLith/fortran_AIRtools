module system_generation
	use sparse_matrices
	implicit none
	real(kind=8),parameter :: pi = 4d0*atan(1d0)
	
	type data
		integer :: n,p,nth
		real(kind=8) :: r,dw,sd
		real(kind=8),allocatable :: T(:),theta(:)
	contains
		procedure :: read => read_data
	end type data
	
	
	
contains
	recursive subroutine quicksort(a,IND, first, last)
	  implicit none
	  real(kind=8) ::  a(:)
	  integer :: IND(:)
	  integer, optional :: first, last
	  
	  !local variables
	  real*8 :: x, t
	  integer :: i, j,f,l,m

	  if (present(first)) then
		  f=first
	  else
		  f = 1
	  endif
	  if (present(last)) then
		  l=last
	  else
		  l = size(a)
	  endif

	  x = a( (f+l) / 2 )
	  i = f
	  j = l
	  do
	     do while (a(i) < x)
	        i=i+1
	     end do
	     do while (x < a(j))
	        j=j-1
	     end do
	     if (i >= j) exit
	     t = a(i);  a(i) = a(j);  a(j) = t
		 m = IND(i); IND(i) = IND(j); IND(j)=m
		 
	     i=i+1
	     j=j-1
	  end do
	  if (f < i-1) call quicksort(a,IND, f, i-1)
	  if (j+1 < l)  call quicksort(a,IND, j+1, l)
	end subroutine quicksort

	function fanlineartomo(dat) result(Amat)
		class(data),intent(in) :: dat
		type(csr_matrix) :: Amat
		
		integer :: i,j,k,m,n_int,count
		real(kind=8) :: R,dw,sd
		real(kind=8) :: x0,y0,a,b,x0th,y0th,ath,bth,dew,xdth,ydth,xm,ym,test_val
		real(kind=8),dimension(dat%n+1) :: x,y,tx,yx,ty,xy
		real(kind=8),dimension(2*dat%n+2) :: t,xxy,yxy,aval,temp_xxy,temp_yxy
		integer,dimension(2*dat%n+2) :: indices,col
		real(kind=8) :: omega(dat%p),deposall(dat%p)
		real(kind=8),allocatable,dimension(:) :: depos
	
	
		!sparse matrix stuff
		integer :: nnz
		real(kind=8),allocatable,dimension(:) :: temp
		integer,allocatable,dimension(:) :: jtemp
		
		R = dat%R*dat%n
		dw = dat%dw*dat%n
		sd = dat%sd*dat%n
		
		!write(*,*) R,dw,sd
		
		x0 = 0d0
		y0 = R
		dew = dw/dat%p

		do i=1,dat%n+1
			x(i) = i - .5d0*dat%n - 1d0
			y(i) = i - .5d0*dat%n - 1d0
		enddo

		

		!Set angles to match linear detector
		if (mod(dat%p,2)==1) then
			allocate(depos((dat%p-1)/2))
			do i=1,(dat%p-1)/2
				depos(i) = i*dew
			enddo
			deposall = (/-depos((dat%p-1)/2:1:-1),0d0,depos/)
			omega = (/atan(sd/depos((dat%p-1)/2:1:-1))-.5d0,0d0,.5d0*pi - atan(sd/depos)  /)
		else
			allocate(depos(dat%p/2))
			do i=1,dat%p/2
				depos(i) = (i-.5d0)*dew
			enddo
			deposall = (/-depos(dat%p/2:1:-1),depos/)
			omega = (/atan(sd/depos(dat%p/2:1:-1))-.5d0*pi,.5d0*pi - atan(sd/depos)  /)
		endif
		

		!initiate CSR arrays
		nnz = 0
		allocate(Amat%IA(0:dat%nth*dat%p))
		Amat%IA(0) = 0
		
		!Perform the entire computation once without recording values to determine nnz.
		!Once nnz is determined, the computation is repeated but now values are recorded.
		!This is waaaaay faster than growing the arrays on the fly. It easily saves a
		!factor 10 in computing time.
		
		do i=1,dat%nth!loop over all angles

			!rotated starting position
			x0th = cos(dat%theta(i))*x0 - sin(dat%theta(i))*y0
			y0th = sin(dat%theta(i))*x0 + cos(dat%theta(i))*y0
			
			
		
		
			!rotated ray parameters
			ath = -x0th/r
			bth = -y0th/r
			
			do j=1,dat%p!loop over all rays
				!find ray parameters
				a = cos(omega(j))*ath - sin(omega(j))*bth
				b = sin(omega(j))*ath + cos(omega(j))*bth
			
				xdth = cos(dat%theta(i))*deposall(j) - sin(dat%theta(i))*(r-sd)
				ydth = sin(dat%theta(i))*deposall(j) + cos(dat%theta(i))*(r-sd)
				
				!use the parametrisation of the line to get the y-coordinates of the intersections of constant x
				tx = (x-x0th)/a
				yx = b*tx + y0th


				!use the parametrisation of the line to get the x-coordinates of the intersections of constant y
				ty = (y-y0th)/b
				xy = a*ty+x0th

				!collect the intersection times and coordinates
				t = (/tx,ty/)
				xxy = (/x,xy/)
				yxy = (/yx,y/)
			
				!sort the intersection times and sort along the coordinates
				indices = (/( k, k = 1, 2*dat%n+2 )/)
				call quicksort(t,indices)
				xxy = xxy(indices)
				yxy = yxy(indices)
			
			
				!not all intersection points are valid, sift out those that fall outside
				n_int = 0
				temp_xxy = 0d0
				temp_yxy = 0d0
				do k=1,2*dat%n+2
					if (xxy(k) >= -.5d0*dat%n .and. xxy(k) <= .5d0*dat%n .and. yxy(k) >= -.5d0*dat%n .and. yxy(k) <= .5d0*dat%n) then
						n_int = n_int+1
						temp_xxy(n_int) = xxy(k)
						temp_yxy(n_int) = yxy(k)
					endif
				enddo
				
				!remove double points
				m = n_int
				n_int = 0
				do k=1,m-1
					if ( sqrt( ( temp_xxy(k) - temp_xxy(k+1) )**2 + ( temp_yxy(k) - temp_yxy(k+1)  )**2   ) > 1d-10 ) then
						n_int = n_int+1
					
					endif
				enddo
			
				!number of elements in the row is one less than the number of intersections
				
				nnz = nnz + max(n_int,0)
				
				!record the number of nonzeros up to this point
				Amat%IA( (i-1)*dat%p+j ) = Amat%IA( (i-1)*dat%p+j-1 ) + max(n_int,0)
			enddo
		enddo
		
		!Now that nnz is known, we can allocate IA and vals.
		Amat%m = dat%nth*dat%p
		Amat%n = dat%n**2
		Amat%nnz = nnz
		allocate(Amat%vA(nnz))
		allocate(Amat%JA(nnz))
		
		count = 0
		do i=1,dat%nth!loop over all angles

			!rotated starting position
			x0th = cos(dat%theta(i))*x0 - sin(dat%theta(i))*y0
			y0th = sin(dat%theta(i))*x0 + cos(dat%theta(i))*y0
			
			
		
		
			!rotated ray parameters
			ath = -x0th/r
			bth = -y0th/r
			
			do j=1,dat%p!loop over all rays
				!find ray parameters
				a = cos(omega(j))*ath - sin(omega(j))*bth
				b = sin(omega(j))*ath + cos(omega(j))*bth
			
				xdth = cos(dat%theta(i))*deposall(j) - sin(dat%theta(i))*(r-sd)
				ydth = sin(dat%theta(i))*deposall(j) + cos(dat%theta(i))*(r-sd)
				
				!use the parametrisation of the line to get the y-coordinates of the intersections of constant x
				tx = (x-x0th)/a
				yx = b*tx + y0th


				!use the parametrisation of the line to get the x-coordinates of the intersections of constant y
				ty = (y-y0th)/b
				xy = a*ty+x0th

				!collect the intersection times and coordinates
				t = (/tx,ty/)
				xxy = (/x,xy/)
				yxy = (/yx,y/)
			
				!sort the intersection times and sort along the coordinates
				indices = (/( k, k = 1, 2*dat%n+2 )/)
				call quicksort(t,indices)
				xxy = xxy(indices)
				yxy = yxy(indices)
			
			
				!not all intersection points are valid, sift out those that fall outside
				n_int = 0
				temp_xxy = 0d0
				temp_yxy = 0d0
				do k=1,2*dat%n+2
					if (xxy(k) >= -.5d0*dat%n .and. xxy(k) <= .5d0*dat%n .and. yxy(k) >= -.5d0*dat%n .and. yxy(k) <= .5d0*dat%n) then
						n_int = n_int+1
						temp_xxy(n_int) = xxy(k)
						temp_yxy(n_int) = yxy(k)
					endif
				enddo



				!remove double points
				m = 0
				do k=1,n_int-1
					test_val = sqrt( ( temp_xxy(k) - temp_xxy(k+1) )**2 + ( temp_yxy(k) - temp_yxy(k+1)  )**2 )
					if ( test_val  > 1d-10 ) then
						m = m+1
						count = count+1
						Amat%vA(count) = test_val
						
						xm = .5d0*( temp_xxy(k) + temp_xxy(k+1) + dat%n )
						ym = .5d0*( temp_yxy(k) + temp_yxy(k+1) + dat%n )

						!compute column numbers
						Amat%JA(count) = floor(xm)*dat%n + (dat%n - floor(ym))
						
					
					endif
				enddo
				
				!write(*,*) count, Amat%IA((i-1)*dat%p + j)
				
			enddo
		enddo
		
		if (count<nnz) stop "Didn't hit all elements... somehow"
		
		
	end function fanlineartomo
	
	subroutine read_data(dat,inputfile,verbose)
		implicit none
		
		character(30),intent(in) :: inputfile
		class(data),intent(out) :: dat
		logical,optional :: verbose
		 
		integer :: l,i
		character(30) :: trash
		logical :: verb
		
		verb = .false.
		if (present(verbose)) verb = verbose
	
		if (verb) write(*,*) "Opening file: ",trim(inputfile)
	
		open (unit = 10, file = trim(inputfile))
	
		read(10,*) trash, dat%n
		read(10,*) trash, dat%p
		read(10,*) trash, dat%nth
		read(10,*) trash, dat%r
		read(10,*) trash, dat%dw
		read(10,*) trash, dat%sd
	
		read(10,*) trash,l
		if (l /= dat%nth) stop "Angle vector has the wrong length, nth is incorrect."
		
		allocate(dat%theta(l))
		do i=1,l
			read(10,*) dat%theta(i)
		enddo
		dat%theta = dat%theta/180d0*pi
		
		
		read(10,*) trash,l
		if (l /= dat%nth*dat%p) stop "Data has the wrong size, or p and nth are incorrect."
	
		allocate(dat%T(l))
		do i=1,l
			read(10,*) dat%T(i)
		enddo
	
		close(10)
	end subroutine read_data
	
	subroutine read_ordering(ordering,inputfile,verbose)
		implicit none
		
		character(30),intent(in) :: inputfile
		integer, allocatable :: ordering(:)
		logical,optional :: verbose
		 
		integer :: l,i
		character(30) :: trash
		logical :: verb
		
		verb = .false.
		if (present(verbose)) verb = verbose
	
		if (verb) write(*,*) "Opening file: ",trim(inputfile)
	
		open (unit = 11, file = trim(inputfile))
		
		read(11,*) trash,l
		
		allocate(ordering(l))
		do i=1,l
			read(11,*) ordering(i)
		enddo
		
		close(11)
	end subroutine read_ordering
	
	subroutine add_white_noise(bt,e,b,eta)
		implicit none
		real(kind=8),intent(in) :: b(:),eta
		
		real(kind=8),dimension(size(b)) :: bt,e
		
		!local
		real(kind=8),parameter :: tol = 1d-14
		real(kind=8) :: u,x,f,df
		integer :: i,k,max_iter
		
		!need the inverse of erf
		
		max_iter = 100
		
		
		do k = 1,size(b)
		
		call random_number(u)
			x = 0d0
			do i=1,max_iter
				f = .5d0*(1+erf(x/sqrt(2d0))) - u
				df = sqrt(2d0/pi)*exp(-.5d0*x**2)
			
				if (abs(f/df) < tol) exit
				
				if (i==max_iter .and. abs(f/df)>tol) stop "Newton method did not converge."
			
				x = x - f/df
			enddo
			e(k) = x
		enddo
		
		e = eta*norm2(b)*e/norm2(e)
		bt = b + e
		
	end subroutine add_white_noise
	
end module system_generation
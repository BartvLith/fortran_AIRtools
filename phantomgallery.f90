module phantomgallery
	use system_generation, only : pi
	implicit none
	
	! All functions in this module create a phantom with the of size N x N.
	!
	! Input:
	!    N    Scalar denoting the nubmer of discretization intervals in each
	!          dimesion, such that the phantom head consists of N^2 cells.
	!
	! Output:
	!    X    The modified phantom head as a 2-dimensional array
	!
	!Some functions have optional inputs to modify the phantoms somewhat.
contains
	
	function linspace(a,b,n) result(x)
		!Creates an array of size n with values linearly distributed between a and b.
		!The endpoints, a and b, are included.
		implicit none
		real(kind=8),intent(in) :: a,b
		integer,intent(in) :: n
		real(kind=8) :: x(n)
		
		!local
		integer :: i
		real(kind=8) :: c,d
		
		c = (b-a)/(n-1)
		d = a-c
		
		do i=1,n
			x(i) = c*i+d
		enddo
	end function linspace
	
	subroutine meshgrid(x,y,N,corners)
		!create a rectangular meshgrid of N x N on the unit square [-1,1]^2
		!optional input: the four corners of rectangle.
		!corners = (/x0,x1,y0,y1/)
		
		implicit none
		real(kind=8),dimension(0:N-1,0:N-1),intent(out) :: x,y
		integer,intent(in) :: N
		real(kind=8),optional :: corners(4)
		
		!local
		real(kind=8) :: x0,x1,y0,y1
		integer :: i
		
		if (present(corners)) then
			x0 = corners(1)
			x1 = corners(2)
			y0 = corners(3)
			y1 = corners(4)
		else
			x0 = -1
			x1 = 1
			y0 = -1
			y1 = 1
		endif
			
		
		
		!fill out a coordinate mesh for x
		x(0,:) = linspace(x0,x1,N)
		do i=1,N-1
			x(i,:) = x(0,:)
		enddo
		
		!coordinate mesh y is x rotated by 90 degrees
		y(:,0) = linspace(y1,y0,N)
		do i=1,N-1
			y(:,i) = y(:,0)
		enddo
	end subroutine meshgrid
	
	function choose_phantom(phantom,N) result(x)
		!Just a wrapper function that picks the right phantom based on the input string phantom
		implicit none
		character(*),intent(in) :: phantom
		integer,intent(in) :: N
		real(kind=8) :: x(N,N)
		
		select case(phantom)
			case('shepplogan')
				x = shepplogan(N)
			case('smooth')
				x = smooth(N)
			case('threephases')
				x = threephases(N)
			case default
				stop 'Not a recognised phantom.'
		end select
		
	end function choose_phantom
	
	function shepplogan(N) result(x)
		! This head phantom is the same as the Shepp-Logan except the intensities
		! are changed to yield higher contrast in the image.
		!
		! Peter Toft, "The Radon Transform - Theory and Implementation", PhD
		! thesis, DTU Informatics, Technical University of Denmark, June 1996.
		
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		
		!local variables
		real(kind=8) :: e(10,6)
		real(kind=8),dimension(0:N-1,0:N-1) :: xc,yc
		real(kind=8) :: a2,b2,x0,y0,phi,a
		integer :: i
		
		!information on all the elipses
		!              A          a          b         x0            y0        phi
		!         -----------------------------------------------------------------
		e(1,:) = (/   1d0  , .6900d0  , .9200d0  ,    0d0  ,     0000d0  ,     0d0 /)
		e(2,:) = (/ -.8d0  , .6624d0  , .8740d0  ,    0d0  ,   -.0184d0  ,     0d0 /)
		e(3,:) = (/ -.2d0  , .1100d0  , .3100d0  ,  .22d0  ,        0d0  ,   -18d0 /)
		e(4,:) = (/ -.2d0  , .1600d0  , .4100d0  , -.22d0  ,        0d0  ,    18d0 /)
		e(5,:) = (/  .1d0  , .2100d0  , .2500d0  ,    0d0  ,      .35d0  ,     0d0 /)
		e(6,:) = (/  .1d0  , .0460d0  , .0460d0  ,    0d0  ,       .1d0  ,     0d0 /)
		e(7,:) = (/  .1d0  , .0460d0  , .0460d0  ,    0d0  ,      -.1d0  ,     0d0 /)
		e(8,:) = (/  .1d0  , .0460d0  , .0230d0  , -.08d0  ,    -.605d0  ,     0d0 /)
		e(9,:) = (/  .1d0  , .0230d0  , .0230d0  ,    0d0  ,    -.606d0  ,     0d0 /)
		e(10,:)= (/  .1d0  , .0230d0  , .0460d0  ,  .06d0  ,    -.605d0  ,     0d0 /)
		
		
		call meshgrid(xc,yc,N)
		
		x = 0d0
		!for every elipse 
		do i=1,size(e,1)
			A = e(i,1)
			a2 = e(i,2)**2
			b2 = e(i,3)**2
			x0 = e(i,4)
			y0 = e(i,5)
			phi = e(i,6)*pi/180d0
			
			
			where (((xc-x0)*cos(phi) + (yc-y0)*sin(phi))**2/a2 + (( (yc-y0)*cos(phi) - (xc-x0)*sin(phi)))**2/b2 <= 1d0 ) x=x+A
		enddo
		
	end function shepplogan
	
	function smooth(N,pin) result(x)
		!SMOOTH Creates a 2D test image of a smooth function
		!Per Christian Hansen, May 8, 2012, DTU Compute.
		
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		integer,optional :: pin
		
		real(kind=8),dimension(0:N-1,0:N-1) :: xc,yc
		integer :: i,p
		real(kind=8) :: sigma,c(4,2),a(4)
		
		if (present(pin)) then
			p = pin
		else
			p = 4
		endif
		
		call meshgrid(xc,yc,N,(/1,N,N,1/)*1d0)
 		sigma = 0.25d0*N
 		c(1,:) = (/0.6d0*N, 0.6d0*N/)
		c(2,:) = (/0.5*N, 0.3*N/)
		c(3,:) = (/ 0.2*N ,0.7*N/)
		c(4,:) =  (/0.8*N , 0.2*N/)
 		a = (/1d0, 0.5d0, 0.7d0, 0.9d0/)
		x = 0d0
 		do i=1,p
			x = x + a(i)*exp( - (xc-c(i,1))**2/(1.2d0*sigma)**2 - (yc-c(i,2))**2/sigma**2  )
! 		    im = im + a(i)*exp( - (I-c(i,1)).^2/(1.2*sigma)^2 - (J-c(i,2)).^2/sigma^2);
 		enddo
 		x = x/maxval(x)
		
	end function smooth
	
	function threephases(N,pin,seed) result(x)
		!THREEPHASES Creates a 2D test image with three different phases
		!Per Christian Hansen, Sept. 30, 2014, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		integer,optional :: pin,seed(:)
		
		!local
		integer :: i,p
		real(kind=8),dimension(0:N-1,0:N-1) :: xc,yc
		real(kind=8) :: sigma1,sigma2,t1,t2,im1(N,N),im2(N,N)
		real(kind=8),allocatable :: c(:,:)
		
		if (present(pin)) then
			p = pin
		else
			p = 100
		endif
		if (present(seed)) call random_seed(put=seed)
		
		call meshgrid(xc,yc,N,(/1,N,N,1/)*1d0)
		
		sigma1 = 0.025d0*N
		allocate(c(p,2))
		call random_number(c)
		c = c*N
		im1 = 0d0
		
		!generate first image
		do i=1,p
			im1 = im1 + exp(-abs(xc-c(i,1))**3/(2.5d0*sigma1)**3 - abs(yc-c(i,2))**3/(2.5d0*sigma1)**3  )
		enddo
		t1 = .35d0
		where( im1 <  t1) im1 = 0
		where( im1 >= t1) im1 = 2d0
		
		!generate second image
		sigma2 = 0.025d0*N
		call random_number(c)
		c = c*N
		im2 = 0d0
		do i=1,p
			im2 = im2 + exp(-abs(xc-c(i,1))**3/(2.5d0*sigma2)**2 - abs(yc-c(i,2))**3/(2.5d0*sigma2)**2  )
		enddo
		t2 = 0.55d0
		where( im2 < t2) im2 = 0d0
		where( im2 >= t2) im2 = 1d0
		
		!combine the two images
		x = im1 + im2
		where( x == 3d0 ) x = 1d0
		x = x/maxval(x);
	end function threephases
	
	function vectorise(x) result(v)
		!Vectorise an N-by-N matrix into a vector according to columns.
		!This is the same way MATLAB vectorises things, conveniently.
		implicit none
		real(kind=8),intent(in) :: x(:,:)
		real(kind=8) :: v(size(x))
		
		integer :: i,j,k
		
		v = 0d0
		k = 0
		do j=1,size(x,2)
			do i=1,size(x,1)
				k = k+1
				v(k) = x(i,j)
			enddo
		enddo
	
	end function vectorise
end module phantomgallery
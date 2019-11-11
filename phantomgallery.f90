module phantomgallery
	use system_generation, only : pi
	implicit none
	
	
contains
	
	function linspace(a,b,n) result(x)
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
	
	function shepplogan(N) result(x)
		
		! This function creates the modifed Shepp-Logan phantom with the
		! discretization N x N, and returns it as a vector.
		!
		! Input:
		!    N    Scalar denoting the nubmer of discretization intervals in each
		!          dimesion, such that the phantom head consists of N^2 cells.
		!
		! Output:
		!    X    The modified phantom head reshaped as a vector
		!
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
		
		
		!fill out a coordinate mesh for x
		xc(0,:) = ( (/ (i,i=0,N-1) /) -(N-1)/2d0 )  /((N-1)/2d0)
		do i=1,N-1
			xc(i,:) = xc(0,:)
		enddo
		
		!coordinate mesh y is x rotated by 90 degrees
		yc = transpose(xc(:,N-1:0:-1))
		
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
	
	
	function vectorise(x) result(v)
		!vectorise an N-by-N matrix into a vector according to how matlab reshapes things
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
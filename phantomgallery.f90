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
	
	subroutine intmeshgrid(x,y,N)
		!create a rectangular meshgrid of N x N on the unit square [-1,1]^2
		!optional input: the four corners of rectangle.
		!corners = (/x0,x1,y0,y1/)
		
		implicit none
		integer,dimension(0:N-1,0:N-1),intent(out) :: x,y
		integer,intent(in) :: N
		
		!local
		integer :: i
		
			
		!fill out a coordinate mesh for x
		x(0,:) = [(i, i=1,N)]
		do i=1,N-1
			x(i,:) = x(0,:)
		enddo
		
		!coordinate mesh y is x rotated by 90 degrees
		y(:,0) = [(i, i=1,N)]
		do i=1,N-1
			y(:,i) = y(:,0)
		enddo
	end subroutine intmeshgrid
	
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
			case('3phases')
				x = threephases(N)
			case('3phasesmooth')
				x = threephases_smooth(N)
			case('binary')
				x = binary(N)
			case('4phases')
				x = fourphases(N)
			case('mst')
				x = mst(N)
			case('grains')
				x = grains(N)
			case('bubbles')
				x = bubbles(N)
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
		
		integer,dimension(0:N-1,0:N-1) :: xc,yc
		integer :: i,p
		real(kind=8) :: sigma,c(4,2),a(4)
		
		if (present(pin)) then
			p = pin
		else
			p = 4
		endif
		
		call intmeshgrid(xc,yc,N)
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
	
	function threephases_smooth(N,pin,vin,seed) result(x)
		!THREEPHASES Creates a 2D test image with three different phases
		!Per Christian Hansen, Sept. 30, 2014, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		integer,optional :: pin,seed(:)
		real(kind=8),optional :: vin
		
		!local
		integer :: i,p
		real(kind=8),dimension(0:N-1,0:N-1) :: xc,yc
		real(kind=8) :: sigma1,sigma2,sigma_bg,t1,t2,im1(N,N),im2(N,N),bg(N,N),v
		real(kind=8),allocatable :: c(:,:)
		
		if (present(pin)) then
			p = pin
		else
			p = 100
		endif
		if (present(seed)) call random_seed(put=seed)
		if (present(vin)) then
			v = vin
		else
			v = 1.8d0
		endif
		
		call meshgrid(xc,yc,N,(/1,N,N,1/)*1d0)
		
		!nullify images
		im1 = 0d0
		im2 = 0d0
		
		sigma1 = 0.025d0*N
		allocate(c(p,2))
		call random_number(c)
		c = c*N
				
		!generate first image
		do i=1,p
			im1 = im1 + exp(-abs(xc-c(i,1))**3/(2.5d0*sigma1)**3 - abs(yc-c(i,2))**3/(2.5d0*sigma1)**3  )
		enddo
		t1 = .35d0
		where( im1 <  t1) im1 = 0
		where( im1 >= t1) im1 = (im1 - minval(im1))/maxval(im1)*v + .8d0
		
		!generate second image
		sigma2 = 0.025d0*N
		call random_number(c)
		c = c*N
		do i=1,p
			im2 = im2 + exp(-abs(xc-c(i,1))**3/(2.5d0*sigma2)**2 - abs(yc-c(i,2))**3/(2.5d0*sigma2)**2  )
		enddo
		t2 = 0.55d0
		where( im2 < t2) im2 = 0d0
		where( im2 >= t2) im2 = (im2 - minval(im2))/maxval(im2)*v + .3d0
			
		!generate smooth background
		call random_number(c)
		c = c*N
		bg = 0d0
		do i=1,p
			call random_number(sigma_bg)
			sigma_bg = N*(5d0*sigma_bg + .5d0)
			bg = bg + (1d-1/p)*exp(-.5d0*(xc-c(i,1))**2/(sigma_bg)**2 - .5d0*(yc-c(i,2))**2/(sigma_bg)**2  )
		enddo
		
		!combine the two images
		x = bg
		where (im1>0) x = im1
		where (im2>0) x = im2
		x = x/maxval(x);
	end function threephases_smooth
	
	function binary(N,pin,stepsin,seed) result(x)
		!binary creates a 2D test image with two phases that look like tunnels
		!Bart van Lith, Nov. 13, 2019, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		integer,optional :: pin,stepsin,seed(:)
		
		!local
		integer :: i,j,p,steps
		real(kind=8),dimension(0:N-1,0:N-1) :: xc,yc,r2
		real(kind=8) :: psize,stepsize,vert(3),horz(3),prob,pos(2),sigma
		
		if (present(pin)) then
			p = pin
		else
			p = 100
		endif
		if (present(stepsin)) then
			steps = stepsin
		else
			steps = N
		endif
		if (present(seed)) call random_seed(put=seed)
		
		x = 0d0
		!the idea is to generate p few random walks that little ball around
		!the current position at every step.
		!All the specifics can be randomised like step size, probability, etc.
		
		call meshgrid(xc,yc,N)
		psize = xc(1,2)-xc(1,1)
		
		!vert = (/.3d0,.7d0,.1d0/) !probabilities for up, stay, down
		!horz = (/.3d0,.4d0,.3d0/) !go left, stay, right
		
		do j=1,p
		
			
			call random_number(stepsize)
			
			call random_number(vert)
			vert(2) = (10d0 + 10d0/stepsize)*vert(2)
			vert = (/vert(1),vert(1)+vert(2),sum(vert)/)
			vert = vert/maxval(vert)
			call random_number(horz)
			!horz(1) = (3d0 + 5d0*stepsize)*horz(1)
			horz(2) = 0d0
			horz = (/horz(1),horz(1)+horz(2),sum(horz)/)
			horz = horz/maxval(horz)
			
			stepsize = .75d0*psize*(stepsize+1)
		
			call random_number(pos)
			pos = 2d0*pos-1d0

			sigma = .025d0
		
			do i=1,steps
				call random_number(prob)
				if (prob<vert(1)) then
					pos(2) = pos(2)-stepsize
				elseif (prob>vert(2)) then
					pos(2) = pos(2)+stepsize
				endif
			
				call random_number(prob)
				if (prob<horz(1)) then
					pos(1) = pos(1)-stepsize
				elseif (prob>horz(2)) then
					pos(1) = pos(1)+stepsize
				endif
				
				if (abs(pos(1))>1d0 .or. abs(pos(2))>1d0 ) exit
				
				r2 = (xc-pos(1))**2 + (yc-pos(2))**2
! 				x = x + exp(-.5d0*((xc-pos(1))**2)/sigma**2 )*exp(-.5d0*((yc-pos(2))**2 )/sigma**2 )*&
! 				( (yc-pos(2))**4 - 6d0*sigma**2*(yc-pos(2))**2 + 3d0*sigma**4  )/(3d0*sigma**4)
! 				where (x<.07d0) x = 0d0
! 				where (x>=.07d0) x = 1d0
				where(r2 <= stepsize**2) x = 1d0
			enddo
			
		enddo
		
		
	end function binary
	
	function fourphases(N,seed) result(x)
		!FOURPHASES Creates a 2D test image with four different phases.
		!Bart van Lith, Nov. 13, 2019, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: x(N,N)
		integer,optional :: seed(:)
		
		integer :: i
		
		if (present(seed)) call random_seed(put=seed)
		
		x = 0d0
		do i=1,8
			x = x + (-1d0)**i*binary(N,80)
		enddo
		where (x < 0d0) x = 0d0
		x = x/maxval(x)
		
	end function fourphases
	
	function mst(N,pin,seed) result(img)
		!MST creates a 2D test image with a minimum spanning tree over normal random points.
		!Bart van Lith, Nov. 13, 2019, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: img(N,N)
		integer,optional :: seed(:),pin
		
		integer :: i,j,ii,jj,k,p,pow,t,l,addt,addl,lind
		integer,allocatable :: tree(:),loose(:)
		real(kind=8) :: dmin,d,xm(2),r,th,u
		integer,dimension(1:N,1:N) :: x,y
		real(kind=8),allocatable :: dist(:,:),x0(:,:)
		
		if (present(pin)) then
			p = pin
		else
			p = 30*sqrt(N*1d0)
		endif
		if (present(seed)) call random_seed(put=seed)
		
		
		pow = 1
		call intmeshgrid(x,y,N)
		
		allocate(x0(2,p))
		do i=1,2
			do j=1,p
				x0(i,j) = 0d0
				!do while (.not. (1d0 < x0(i,j) .and. x0(i,j) < N*1d0))
					x0(i,j) = .2d0*N*randn() + .5d0*N
				!enddo
			enddo
		enddo
		
		allocate(dist(p,p))
		do i=1,p
			do j=1,p
				dist(i,j) = norm2(x0(:,i)-x0(:,j))
			enddo
		enddo
		
		
		!use Prim's algorithm to determine the minimum spanning tree
		!there are faster algorithms, but those are harder to implement
		img = 0d0
		allocate(tree(p))
		allocate(loose(p-1))
		l = p-1
		t = 1
		tree(t) = 1
		loose = [(i,i=2,p)]
		do i = 2,p
			
			dmin = huge(dmin)
			do j=1,t
				do k=1,l
					d = dist(tree(j),loose(k))
					if (d<dmin) then
						dmin = d
						addt = tree(j)
						addl = k
					endif
				enddo
			enddo
			
			t = t+1
			tree(t) = loose(addl)
			ii = tree(t)
			jj = addt
			loose(1:l-1) = (/loose(1:addl-1),loose(addl+1:l)/)
			l = l-1
			
			xm = .5d0*(x0(:,ii)+x0(:,jj))
			th = atan2(x0(2,ii)-x0(2,jj),x0(1,ii)-x0(1,jj))
			r = sum(abs(xm-x0(:,jj))**pow )
			call random_number(u)
			!call draw_line(img,x0(:,addt),x0(:,tree(t)))
			where ( abs(  cos(th)*(x-xm(1)) + sin(th)*(y-xm(2))  )**pow&
				 +abs( -sin(th)*(x-xm(1)) + cos(th)*(y-xm(2)) )**pow < r )&
				  img = 1d0+u!N - (abs(xm(1)-.5d0*N) + abs(xm(2)-.5d0*N))
			 !img( nint( x0(2,ii) ) , nint( x0(1,ii) )  ) = 1d0
			 !img( nint( x0(2,jj) ) , nint( x0(1,jj) )  ) = 1d0
		enddo
		
		img = img/maxval(img)
		
		
	end function mst
	
	function grains(N,randomtype,pin,seed) result(img)
		!GRAINS Creates a test image of Voronoi cells

		!Jakob Sauer Jorgensen, October 9, 2012, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: img(N,N)
		character(8),optional :: randomtype
		integer,optional :: seed(:),pin
		
		integer :: i,j,k,kmin,p
		real(kind=8) :: dmin,d,xloc(2)
		real(kind=8),dimension(1:N,1:N) :: x,y
		real(kind=8),allocatable :: x0(:,:),vals(:)
		character(8) :: rt
		
		if (present(pin)) then
			p = pin
		else
			p = 3*sqrt(N*1d0)
		endif
		if (present(randomtype)) then
			rt = randomtype
		else
			rt = 'uniform'
		endif
		if (present(seed)) call random_seed(put=seed)
		
		
		!create random points
		allocate(x0(2,p))
		if (rt == 'gaussian') then
			do i=1,2
				do j=1,p
					x0(i,j) = -2d0
					do while (.not. (-1d0 < x0(i,j) .and. x0(i,j) < 1d0))
						x0(i,j) = randn()
					enddo
				enddo
			enddo
		elseif (rt == 'uniform') then
			call random_number(x0)
			x0 = 1.5d0*(2d0*x0-1d0)
		else
			stop "Not a recognised input for randomtype. Use 'gaussian' or 'uniform'."
		endif
		
		allocate(vals(p))
		call random_number(vals)
		
		call meshgrid(x,y,N)
		
		do i=1,N
			do j=1,N
				xloc = (/x(i,j),y(i,j)/)
				dmin = huge(dmin)
				do k=1,p
					d = norm2(xloc - x0(:,k))
					
					if (d<dmin) then
						dmin = d
						kmin = k
					endif
				enddo
				
				img(i,j) = kmin
			enddo
		enddo
		
		img = img/maxval(img)
		
	end function grains
	
	function bubbles(N,randomtype,pin,seed) result(img)
		!BUBBLES Creates a test image of a random circle packing.

		!Bart van Lith, November 17, 2019, DTU Compute.
		implicit none
		integer, intent(in) :: N
		real(kind=8) :: img(N,N)
		character(8),optional :: randomtype
		integer,optional :: pin
		integer,optional :: seed(:)
		
		integer :: i,j,kmin,p
		real(kind=8) :: dmin,d,xloc(2),u
		real(kind=8),dimension(1:N,1:N) :: x,y
		real(kind=8),allocatable :: x0(:,:),vals(:)
		character(8) :: rt
		
		if (present(pin)) then
			p = pin
		else
			p = 15*sqrt(N*1d0)
		endif
		if (present(randomtype)) then
			rt = randomtype
		else
			rt = 'uniform'
		endif
		if (present(seed)) call random_seed(put=seed)
		
		call meshgrid(x,y,N)
		
		!create random points
		allocate(x0(2,p))
		if (rt == 'gaussian') then
			do i=1,2
				do j=1,p
					x0(i,j) = -2d0
					do while (.not. (-1d0 < x0(i,j) .and. x0(i,j) < 1d0))
						x0(i,j) = .5d0*randn()
					enddo
				enddo
			enddo
		elseif (rt == 'uniform') then
			call random_number(x0)
			x0 = 1.5d0*(2d0*x0-1d0)
		else
			stop "Not a recognised input for randomtype. Use 'gaussian' or 'uniform'."
		endif
		img = 0d0
		do i=1,p
			dmin = huge(dmin)
			do j=1,p
				if (i==j) cycle
				
				
				d = norm2(x0(:,i)-x0(:,j))
				if (d<dmin) then
					dmin = d
					kmin = j
				endif
			enddo
			
			call random_number(u)
			where ((x-x0(1,i))**2 + (y-x0(2,i))**2 <= (1d0*dmin)**2  ) img = kmin
		enddo
		
		img = img/maxval(img)
		
	end function bubbles
	
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
	
	subroutine swapcoord(p1, p2)
	    real(kind=8), intent(inout) :: p1, p2
	    real(kind=8) :: t

	    t = p2
	    p2 = p1
	    p1 = t
	end subroutine swapcoord
 
	subroutine draw_line(img, x0,x1)
		implicit none
		real(kind=8), intent(inout) :: img(:,:)
		real(kind=8), intent(in) :: x0(2),x1(2)

		real(kind=8),dimension(2) :: rfrom, rto
		real(kind=8) :: gradient,intery,xend,yend,xgap,dx, dy
		integer :: error, ystep, x, y
		integer :: xpx11,ypx11,xpx12,ypx12,N
		logical :: steep
		
		N = size(img,1)

		rfrom = x0
		rto = x1
		steep = (abs(x1(2) - x0(2)) > abs(x1(1) -x0(1)))
		if ( steep ) then
		   call swapcoord(rfrom(1), rfrom(2))
		   call swapcoord(rto(1), rto(2))
		end if
		if ( rfrom(1) > rto(1) ) then
		   call swapcoord(rfrom(1), rto(1))
		   call swapcoord(rfrom(2), rto(2))
		end if

		dx = rto(1) - rfrom(1)
		dy = rto(2) - rfrom(2)
		gradient = dy/dx


		xend = nint(rfrom(1))
		yend = rfrom(2) + gradient*(xend-rfrom(1))
		xgap = 1d0 - (rfrom(1) + .5d0 - floor(rfrom(1) + .5d0))
		xpx11 = max(1d0,min(xend*1d0,N*1d0))
		ypx11 = max(1d0,1d0*floor(yend*1d0))
		if (steep) then
			!write(*,*) xpx11,ypx11
			img(xpx11,ypx11) = (1d0 - (yend - floor(yend)))*xgap
			if (ypx11+1 <= N ) img(xpx11,ypx11+1) = (yend-floor(yend)) * xgap
		else
			!write(*,*) xpx11,ypx11
			img(ypx11,xpx11) = (1d0 - (yend - floor(yend)))*xgap
			if (xpx11<=N) img(ypx11,xpx11+1) = (yend-floor(yend)) * xgap
		endif
		intery = yend + gradient


		xend = nint(rto(1))
		yend = rto(2) + gradient*(xend-rto(1))
		xgap = 1d0 - (rto(1) + .5d0 - floor(rto(1) + .5d0))
		xpx12 = max(1d0,min(xend*1d0,N*1d0))
		ypx12 = max(1d0,1d0*floor(yend*1d0))
		if (steep) then
			img(xpx12,ypx12) = (1d0 - (yend - floor(yend)))*xgap
			if (ypx12+1 <= N ) img(xpx12,ypx12+1) = (yend-floor(yend)) * xgap
		else
			img(ypx12,xpx12) = (1d0 - (yend - floor(yend)))*xgap
			if (xpx12<=N) img(ypx12,xpx12+1) = (yend-floor(yend)) * xgap
		endif

		if (steep) then
			do x = xpx11+1,xpx12-1
				if (floor(intery)+1 > N .or. floor(intery)<1) exit
				img(x,floor(intery)) = 1d0 - (intery - floor(intery))
				img(x,floor(intery)+1) = intery - floor(intery)
		        intery = intery + gradient
			enddo
		else
			do x = xpx11+1,xpx12-1
				if (floor(intery)+1 > N .or. floor(intery)<1) exit
				img(floor(intery),x) = 1d0 - (intery - floor(intery))
				img(floor(intery)+1,x) = intery - floor(intery)
		        intery = intery + gradient
			enddo
		endif

	end subroutine draw_line
	
	function randn(maxits,tolerance) result(x)
		implicit none
		integer,optional :: maxits
		real(kind=8),optional :: tolerance
		real(kind=8) :: x
		
		integer :: i,max_iter
		real(kind=8) ::f,df,tol,u
		
		if (present(maxits)) then
			max_iter = maxits
		else
			max_iter = 100
		endif
		if (present(tolerance)) then
			tol = tolerance
		else
			tol = 1d-14
		endif
		
		x = 0d0
		call random_number(u)
		do i=1,max_iter
			f = .5d0*(1+erf(x/sqrt(2d0))) - u
			df = sqrt(2d0/pi)*exp(-.5d0*x**2)
		
			if (abs(f/df) < tol) exit
			
			if (i==max_iter .and. abs(f/df)>tol) stop "Newton method did not converge."
		
			x = x - f/df
		enddo
	end function randn
	  
end module phantomgallery
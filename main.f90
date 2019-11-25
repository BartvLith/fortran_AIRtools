program main
	use system_generation
	use sparse_matrices
	use airms
	use phantomgallery
	!$ use omp_lib

	implicit none
	
	
	!real(kind=8),parameter :: pi = 4d0*atan(1d0)
	integer :: max_its
	
	
	integer :: i,j,k
	real(kind=8) :: sigma,test
	
	!OOP sparse
	type(csr_matrix) :: A
	real(kind=8) :: grid_size
	
	!Kaczmarz stuff
	type(opts) :: options
	real(kind=8),allocatable :: Xstore(:),Ystore(:),xout(:),EG(:)
	integer :: final_it
	
	!Shepp-Logan
	real(kind=8),allocatable :: phantom(:,:),xref(:),b(:),bt(:),noise(:)
	
	!timing
	integer :: t0,t1,clock_rate,clock_max
	
	!input
	integer :: cac,skip,cull
	character(32) :: flagchar
	character(12) :: recognised_methods(10),recognised_sr(6),recognised_phantom(10)
	character(2) :: flag
	character(30) :: inputfile,outputfile,inputmethod,inputstoprule,inputphantom,inputorder,matrixfile,oraclefile
	logical :: test_mode,verbose_mode,stop_rule,customord,storematrix,readmatrix,useoracle
	type(data) :: dat
	real(kind=8) :: inputlb,inputub,inputrelaxpar,inputnoise
	real(kind=8),allocatable :: data_temp(:,:),theta_temp(:)
	integer,allocatable :: ordering(:)
	
	
	!======================================
	!
	! INPUT
	!
	!======================================
	
	skip = 0
	max_its = 100
	verbose_mode = .false.
	test_mode = .false.
	stop_rule = .false.
	customord = .false.
	storematrix = .false.
	useoracle = .false.
	readmatrix = .false.
	inputnoise = 1d-3
	outputfile = 'sol.txt'
	inputmethod = 'kaczmarz'
	inputphantom = ' '
	recognised_methods = (/'kaczmarz    ','randkaczmarz','rand        ','symkaczmarz ','sym         '&
						  ,'graddescent ','gd          ','cimmino     ','sart        ','sirt        '/)
	recognised_sr =      (/'errorgauge  ','mutualstep  ','eg          ','ms          ','twin        '&
						  ,'oracle      '/)
	recognised_phantom = (/'shepplogan  ','            ','smooth      ','3phases     ','3phasesmooth'&
						  ,'binary      ','4phases     ','mst         ','grains      ','bubbles     '/)
	call options%initialise()
	cac = COMMAND_ARGUMENT_COUNT()
	if (cac == 0) stop "Supply input file or activate test mode."
	if (cac > 9) stop "Too many input arguments."
	
	do i=1,cac
		CALL GET_COMMAND_ARGUMENT(i,flagchar)
		flagchar = adjustl(flagchar)
		flag = flagchar(1:2)
		
		select case (flag)
			case('-h')
				write(*,*) "Main CT AIR tools input flags:"
				write(*,*) "-i[inputfile]    Provide suitable inputfile."
				write(*,*) "-o[outputfile]   Choose output file name."
				write(*,*) "-M[method]       Choose which method to use. Standard is Kaczmarz."
				write(*,*) "                 Choose from: ",("'"//trim(recognised_methods(j))//"', ", j=1,8),&
											 "and '",trim(recognised_methods(9)),"'."
				write(*,*) "-S[stoprule]     Choose what stopping rule to use. Standard is none."
				write(*,*) "                 Choose from 'errorgauge' or 'mutualstep'. Omit if no stopping rule is required."
				write(*,*) "-K[iterations]   Set maximum number of iterations, standard is 100."
				write(*,*) "-O[orderfile]    Ordering of rows, supplied as a list of integers in a file. "
				write(*,*) "-L[lowerbound]   Set lower bound."
				write(*,*) "-U[upperbound]   Set upper bound."
				write(*,*) "-R[relaxpar]     Set relaxation parameter."
				write(*,*) "-t[phantom]      Test mode, runs a test phantom. Standard is Shepp-Logan."
				write(*,*) "                 Possible inputs: 'shepplogan'."
				write(*,*) "-c[cull]         Cull factor. Only use every cull angles, cull must be an integer."
				write(*,*) "-v               Verbose mode."
				write(*,*) "-h               Help."
				stop
			case('-i')
				inputfile = flagchar(3:32)
			case('-o')
				outputfile = flagchar(3:32)
			case('-v')
				write(*,*) "Verbose mode on!"
				verbose_mode = .true.
				options%verbose = .true.
			case('-t')
				test_mode = .true.
				if (flagchar(3:32) /= " ") READ(flagchar(3:32),*) inputphantom
				if ( .not. any(inputphantom == recognised_phantom)) stop "Not a recognised phantom in -t. Use -h for help."
				if (inputphantom == " ") inputphantom = 'shepplogan'
				write(*,*) "Test mode active. Ignoring input file."
			case('-K')
				READ(flagchar(3:32),*)max_its
			case('-c') !set row skip
				READ(flagchar(3:32),*) skip
			case('-M')
				READ(flagchar(3:32),*) inputmethod
				if ( .not. any(inputmethod == recognised_methods) ) then
					stop "Not a recognised stopping rule. Use -h for help."
				endif
				if (trim(inputmethod) == 'rand') inputmethod = 'randkaczmarz'
				if (trim(inputmethod) == 'sym') inputmethod = 'symkaczmarz'
				if (trim(inputmethod) == 'gd') inputmethod = 'graddescent'
				if (trim(inputmethod) == 'sirt') inputmethod = 'sart'
			case('-n')
				READ(flagchar(3:32),*) inputnoise
				if (inputnoise<0d0) stop "Provide positive noise level."
			case('-S')
				stop_rule = .true.
				READ(flagchar(3:32),*) inputstoprule
				if ( .not. any(inputstoprule == recognised_sr)) stop "Not a recognised stopping rule. Choose from 'errorgauge' or 'mutualstep'."
				if (trim(inputstoprule) == 'eg' .or. trim(inputstoprule) == 'twin') inputstoprule = 'errorgauge'
				if (trim(inputstoprule) == 'ms') inputstoprule = 'mutualstep'
			case('-L')
				READ(flagchar(3:32),*) inputlb
				call options%set_lbound(inputlb)
			case('-U')
				READ(flagchar(3:32),*) inputub
				call options%set_ubound(inputub)
			case('-R')
				READ(flagchar(3:32),*) inputrelaxpar
				options%relaxpar = inputrelaxpar
			case('-O')
				READ(flagchar(3:32),*) inputorder
				customord = .true.
				if (verbose_mode) write(*,*) "Reading custom ordering."
				call read_ordering(ordering,inputorder,verbose_mode)
			case('-r')
				READ(flagchar(3:32),*) matrixfile
				readmatrix = .true.
			case('-s')
				READ(flagchar(3:32),*) matrixfile
				storematrix = .true.
			case('-f')
				READ(flagchar(3:32),*) oraclefile
				useoracle = .true.
			case('-H')
				write(*,*) "Oh noes, the hamster! RUN!!!"
				stop
			case default
				stop "Not a recognised flag, use -h for help."
		end select
	enddo
	
	if (useoracle .and. .not. (trim(inputstoprule)=='oracle')) then
		write(*,*) "Overriding stopping rule. Using oracle instead."
		inputstoprule = 'oracle'
		stop_rule = .true.
	endif
	
	if (storematrix .and. readmatrix) stop "Choose either storing or reading matrix."
	
	if (verbose_mode) then
		!$ write(*,*) "Parallel processing enabled."
		!$ write(*,*) "Number of threads:",omp_get_max_threads()
	endif
	
	
	!======================================
	!
	! SETUP
	!
	!======================================
	
	if (test_mode) then
		dat%n=256
		dat%nth = 400
		dat%p = 1.5d0*dat%n
		dat%r  = sqrt(2d0)
		dat%dw = 1d0!2d0*sqrt(2d0)
		dat%sd = 1d0
		
		allocate(dat%theta(dat%nth))
		allocate(phantom(dat%n,dat%n))
		allocate(xref(dat%n**2))
		allocate(b(dat%nth*dat%p))
		allocate(bt(dat%nth*dat%p))
		allocate(noise(dat%nth*dat%p))
		
		dat%theta = linspace(0d0,pi,dat%nth)
		
		call system_clock ( t0, clock_rate, clock_max )
		if (readmatrix) then
			call read_csr_matrix(A,matrixfile,verbose_mode)
		else
			if (verbose_mode) write(*,*) "Constructing matrix..."
			A = paralleltomo(dat)
		endif
		
		phantom = choose_phantom(trim(inputphantom),dat%n)
		xref = vectorise(phantom)
		b = A%matvecmul(xref)
		call add_white_noise(bt,noise,b,inputnoise)
		
		call system_clock ( t1, clock_rate, clock_max )

	else !read in inputfile
		if (verbose_mode) write(*,*) "Reading in data..."
		
		call dat%read(inputfile,verbose_mode)
		
		if (verbose_mode) write(*,*) "done."
		
		if (skip>0) then
			if (verbose_mode) write(*,*) "Culling data..."
		
			allocate(data_temp(dat%p,dat%nth))
			data_temp = reshape(dat%T,(/dat%p,dat%nth/))
		

			cull = size(dat%theta(1:dat%nth:skip))
			allocate(theta_temp(cull))
			theta_temp = dat%theta(1:dat%nth:skip)
			call move_alloc(theta_temp,dat%theta)
			deallocate(dat%T)
			allocate(dat%T(cull*dat%p))
			dat%T = vectorise(data_temp(:,1:dat%nth:skip))
			dat%nth = cull
			if (verbose_mode) write(*,*) "done."
		endif
		
		if (useoracle) then
			allocate(xref(dat%n**2))
			call read_ref(xref,oraclefile,dat%n**2,verbose_mode)
		endif
		
		call system_clock ( t0, clock_rate, clock_max )
		
		
		!This is somehow needed, see the demo in matlab.
		grid_size = dat%dw*dat%R/dat%sd
		dat%R = dat%R/grid_size
		dat%dw = dat%dw/grid_size
		dat%sd = dat%sd/grid_size
		
		if (readmatrix) then
			call read_csr_matrix(A,matrixfile,verbose_mode)
		else
			if (verbose_mode) write(*,*) "Constructing matrix..."
			A = fanlineartomo(dat)
			call A%multiply_by_scalar(grid_size/dat%N)		
		endif
		
		allocate(bt(A%m))
		bt = -log(dat%T)
		
		
		call system_clock ( t1, clock_rate, clock_max )
		
	endif
	
	if (storematrix) then
	 	call store_matrix(A,matrixfile,verbose_mode)
	endif
	
	
	
	
	
	if (verbose_mode) then
		write(*,*) "done in ",real(t1-t0)/real(clock_rate)," seconds!"
		write(*,*) " "
		write(*,*) "Matrix size: ",A%m," by ",A%n
		write(*,*) "nnz: ",A%nnz
		write(*,*) "sparsity level: ",(1d0*A%nnz)/(1d0*A%m*A%n)
		write(*,*) "space saved: ", (1d0-(2d0*A%nnz + A%m + 1)/(1d0*A%m*A%n))*100d0 ,"%"
		write(*,*) " "
		if (test_mode) write(*,*) "Noise level: ",inputnoise
	endif	
	
	
	!======================================
	!
	! EXECUTION
	!
	!======================================
	
	call system_clock ( t0, clock_rate, clock_max )
	
	allocate(EG(max_its))
	allocate(Xstore(dat%n**2))
	allocate(Ystore(dat%n**2))
	allocate(xout(dat%n**2))
	
	if (trim(inputmethod) == 'graddescent') then
		call gradient_descent(xout,A,bt,max_its,options = options)
	elseif (trim(inputmethod) == 'cimmino' .or. trim(inputmethod) == 'sart') then
		call sirt(xout,inputmethod,A,bt,(/max_its/),options = options)
	else
		if (customord) then
			if (stop_rule) then
				if (trim(inputstoprule)=='oracle') then
					call art(Xstore,inputmethod,A,bt,(/max_its/),xref=xref,order=ordering,options = options)
					xout = Xstore
				else
					call twin_art(xout,final_it,Xstore,Ystore,EG,inputmethod,inputstoprule,A,bt,(/max_its/),order=ordering,options = options)
				endif
			else
				call art(Xstore,inputmethod,A,bt,(/max_its/),order=ordering,options = options)
				xout = Xstore
			endif
		else
			if (stop_rule) then
				if (trim(inputstoprule)=='oracle') then
					call art(Xstore,inputmethod,A,bt,(/max_its/),xref=xref,options = options)
					xout = Xstore
				else
					call twin_art(xout,final_it,Xstore,Ystore,EG,inputmethod,inputstoprule,A,bt,(/max_its/),options = options)
				endif
			else
				call art(Xstore,inputmethod,A,bt,(/max_its/),options = options)
				xout = Xstore
			endif
		endif
	endif
	
	if (verbose_mode .and. test_mode) write(*,*) "Final error: ",norm2(xref-xout)/norm2(xref)
	
	call system_clock ( t1, clock_rate, clock_max )
	if (verbose_mode) write(*,*) "Done in ",real(t1-t0)/real(clock_rate)," seconds."
	
	
	!======================================
	!
	! WRITE AWAY DATA
	!
	!======================================
	
	open(unit=23,file=trim(outputfile),action="write",status="replace")
	do k=1,A%n
		write(23,*) xout(k)
	enddo
	close(23)
	
	if (stop_rule) then
		open(unit=24,file="errorgauge.txt",action="write",status="replace")
		do k=1,max_its
			write(24,*) EG(k)
		enddo
		close(24)
	endif

contains
	
	
	
	
	
end program main





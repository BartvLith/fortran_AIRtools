module AIRMs
	use sparse_matrices
		
	private :: initialise_options,set_upperbound,set_lowerbound
	
	type opts
		real(kind=8) :: relaxpar,ubound,lbound,damp
		logical :: verbose,upb,lwrb,rowprobability
	contains
		procedure :: initialise => initialise_options
		procedure :: set_ubound => set_upperbound
		procedure :: set_lbound => set_lowerbound
	end type opts
		
contains
	
	subroutine initialise_options(options)
		!Initialises standard choices of the options
		!Use this subroutine if you only want to set a few options and otherwise standard.
		!For instance, you want vanilla Kaczarmz to be verbose, use:
		
		!call initialise_options(options)
		!options%verbose = .true.
		
		
		implicit none
		class(opts) :: options
		
		options%upb = .false.
		options%lwrb = .false.
		options%ubound = 0d0
		options%lbound = 1d0
		options%damp = 0d0
		options%relaxpar = 1d0
		options%verbose = .false.
		options%rowprobability = .false.
	end subroutine initialise_options
	
	subroutine set_upperbound(options,ubound)
		implicit none
		class(opts),intent(inout) :: options
		real(kind=8) :: ubound
		
		options%upb = .true.
		options%ubound = ubound
		
	end subroutine set_upperbound
	
	subroutine set_lowerbound(options,lbound)
		implicit none
		class(opts),intent(inout) :: options
		real(kind=8) :: lbound
		
		options%lwrb = .true.
		options%lbound = lbound
		
	end subroutine set_lowerbound
	
	subroutine art(Xout,art_method,A,b,Kin,x0,order,options)
		!ART  General interface for all Kaczmarz/ART methods
		!
		!   call art(Xout,art_method,A,b,Kin)
		!   call art(Xout,art_method,A,b,Kin,x0)
		!   call art(Xout,art_method,A,b,Kin,x0,order)
		!   call art(Xout,art_method,A,b,Kin,x0,order,options)
		!
		! Implements a general ART method for the linear system Ax = b:
		!       
		!       x^{k+1} = x^k + relaxpar*(b_i - a_i'*x^k)/(||a_i||_2^2)*a_i
		!
		! where a_i' is the i-th row of A, and the order for i is chosen in 
		! different ways by the provided methods (see kaczmarz, symkaczmarz, 
		! randkaczmarz for details) or a custom order specified by the user.
		! One iteration consists of m such steps (m-1 for symkaczmarz).
		!
		! Input:
		!   art_method  Either one of the strings 'kaczmarz', 'symkaczmarz', or 
		!               'randkaczmarz' to specify one of the provided methods;
		!               default is 'kaczmarz'.
		!               Or a row index vector of length m (for matrix A m-by-n)
		!               with a desired fixed order in which to step through all 
		!               rows of A. Please see demo_custom for an example.
		!   A           m times n matrix, or a function that implements matrix-
		!               vector multiplication with A and A'; see explanation below.
		!   b           m times 1 vector containing the right-hand side.
		!   Kin         Number of iterations. If K is a scalar, then K is the 
		!               maximum number of iterations and only the last iterate is 
		!               returned. If K is a vector, then max(K) is the maximum
		!               number of iterations and only iterates corresponding to the
		!               values in K are returned, together with the last iterate.
		!   x0          n times 1 starting vector. Default: x0 = 0.
		!	order		A specific order to go through the rows. Default order is
		!				cyclic 1,2,...,m.
		!   options     Struct with the following fields:
		!      relaxpar The relaxation parameter. If relaxpar is a scalar < 2 then
		!               it is used in each iteration; default value is 1.
		!               Alternatively, relaxpar can be a function with a diminishing
		!               parameter, e.g., @(j) 1/sqrt(j), where j counts the total
		!               number of row updates.
		!      lbound   Lower bound in box constraint [lbound,ubound]. If scalar,
		!               this value is enforced on all elements of x in each 
		!               iteration. If vector, it must have same size as x and 
		!               then enforces elementwise lower bounds on x. If empty, no
		!               bound is enforced. +/-Inf can be used.
		!      ubound   Upper bound in box constraint [lbound,ubound]. If scalar,
		!               this value is enforced on all elements of x in each 
		!               iteration. If vector, it must have same size as x and 
		!               then enforces elementwise lower bounds on x. If empty, no
		!               bound is enforced. +/-Inf can be used.
		!      damp     A parameter damp to avoid division by very small row norms
		!               by adding damp*max_i{||a_i||_2^2} to ||a_i||_2^2.
		!      verbose  Nonnegative integer specifying whether progress is printed
		!               to screen during iterations. Default=0: no info printed.
		!               1: Print in every iteration. Larger than 1: Print every
		!               verbose'th iteration and first and last.
		!
		! Output:
		!   Xout        Matrix containing the saved iterations as the columns.
		
		implicit none
		!inputs
		character(*) :: art_method
		class(csr_matrix),intent(in) :: A
		real(kind=8),intent(in) :: b(A%m)
		integer,intent(in) :: Kin(:)
		real(kind=8),intent(in),optional :: x0(A%n)
		integer,optional :: order(:)
		class(opts),intent(in),optional :: options
		
		!outputs
		real(kind=8),intent(out) :: Xout(A%n,size(Kin))
		
		!local variables
		integer :: ii,i,j,k,l,iter,maxiter,next_record,nord,ri
		integer :: Kuse(size(Kin))
		real(kind=8),dimension(A%n) :: x
		real(kind=8) :: ub,lb,dmp,omg,cumsum
		logical :: verb,upb,lwrb,isrand,rowprob
		integer :: ord(A%m),tempord(A%m)
		real(kind=8),dimension(A%m) :: nrA,cumul
		logical :: perctdone(9)
		
		if (present(options)) then
			upb = options%upb
			lwrb = options%lwrb
			ub = options%ubound
			lb = options%lbound
			dmp = options%damp
			omg = options%relaxpar
			rowprob = options%rowprobability
			verb = options%verbose
			
			if (verb) then
				write(*,*) "================================"
				write(*,*) " "
				write(*,*) "ART"
				write(*,*) " "
				write(*,*) "================================"
				write(*,*) "Options detected:"
				write(*,*) "Verbose mode activated. I'll do a lot of talking."
				if (upb) then
					write(*,*) "Using upper bound: ",ub
				else
					write(*,*) "No upper bound present."
				endif
				if (lwrb) then
					write(*,*) "Using lower bound: ",lb
				else
					write(*,*) "No lower bound present."
				endif
				if (dmp==0) then
					write(*,*) "No damping coefficient present."
				else
					write(*,*) "Damping coefficient: ",dmp
				endif
				if (omg<=0 .or. omg >=2) stop "Relaxation parameter needs to be in (0,2)."
				if (omg==1d0) then
					write(*,*) "Using standard relaxation parameter: 1."
				else
					write(*,*) "Using relaxation parameter: ",omg
				endif
				write(*,*) "================================"
			endif
		else
			!standard options means no bounds applied, no damping and a relaxation parameter of 1.
			upb = .false.
			lwrb = .false.
			ub = 0d0
			lb = 1d0
			dmp = 0d0
			omg = 1d0
			rowprob = .false.
			verb = .false.
		endif
		
		if (present(x0)) then
			if (size(x0) /= A%n) stop "Size of x0 should be equal to n (#columns)."
			x = x0
			if (verb) write(*,*) "Using custom initial condition."
		else
			x = 0d0
			if (verb) write(*,*) "Using all-zero initial condition."
		endif
		
		if (size(b) /= A%m) stop "Size of b should be equal to m (#rows)."
		
		isrand = .false.
		Kuse = Kin
		select case (trim(art_method))
			case ('kaczmarz')
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using Kaczmarz and custom ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				else
					nord = A%m
					ord(1:A%m) = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using Kaczmarz and standard sequential ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				endif
			case ('symkaczmarz')
				if (any(mod(Kin,2)==1)) stop "Kin should all be multiples of 2."
				Kuse = Kuse/2
				
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using symmetric Kaczmarz and custom ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				else
					nord = A%m
					ord = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using symmetric Kaczmarz and standard sequential ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				endif
			case ('randkaczmarz')
				isrand = .true.
				!Ord needs to be present as it implements the skipping of zero-norm rows.
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using random Kaczmarz; ignoring order."
				else
					nord = A%m
					ord = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using random Kaczmarz."
				endif
				
				if (verb .and. rowprob) then
					write(*,*) "Using fancy row probability selection,",&
					" i.e., probability to select rows is proportional to row norm."
				 else
					 write(*,*) "Using uniform row probabilities."
				 endif
			case default
				stop "Unknown ART method. Use kaczmarz, symkaczmarz or randkaczmarz"
		end select
		
		if (verb) write(*,*) "================================"
		
		!Apply damping.
		nrA = A%row_norms(2)
		nrA = nrA**2
		nrA = nrA + dmp*maxval(nrA)
		
		!select only nonzero rows
		if (verb) write(*,*) "Selecting nonzero rows..."
		l = nord
		nord = 0
		do i=1,l
			if (nrA(ord(i))>0) then !perhaps some tolerance should be used
				nord = nord+1
				tempord(nord) = ord(i)
			endif
		enddo
		ord(1:nord) = tempord(1:nord)
		

		
		if (verb) write(*,*) "done!"
		if (verb) write(*,*) "================================"
		
		
		!if random kaczmarz was chosen as the method
		if (isrand) then
			ord = scramble(nord)
			if (rowprob) then
				if (verb) write(*,*) "Calcultating row probabilities..."
				!row selection probability is proportional to row norms
				
				cumsum = 0
				do i=1,nord
					cumsum = cumsum+nrA(i)
					cumul(i) = cumsum
				enddo
				cumul = cumul/sum(nrA)
				cumul(nord+1:A%m) = 1d0
				if (verb) write(*,*) "done!"
				if (verb) write(*,*) "================================"
			endif
		endif
		
		
		maxiter = maxval(Kuse)
		l = 1
		do iter=1,maxiter
			
			!scramble the rows if random
			if (isrand) ord(1:nord) = scramble(nord)
			
			!go through the nonzero rows
			call basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(x,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			if (iter==Kuse(l)) then
				Xout(:,l) = x
				l = l+1
				if (verb) then
					if (trim(art_method)=='symkaczmarz') then
						write(*,*) "Storing iteration: ",2*iter," out of ",2*maxiter
					else
						write(*,*) "Storing iteration: ",iter," out of ",maxiter
					endif
				endif
			else
				if (verb) then
					if (trim(art_method)=='symkaczmarz') then
						write(*,*) "Iteration: ",2*iter," out of ",2*maxiter
					else
						write(*,*) "Iteration: ",iter," out of ",maxiter
					endif
				endif
			endif
			
		enddo
		
		if (verb) then
			write(*,*) "================================"
			write(*,*) "Verbose signing off."
			write(*,*) "================================"
			write(*,*) " "
			write(*,*) "ART"
			write(*,*) " "
			write(*,*) "================================"
		endif
	end subroutine art
	
	subroutine twin_art(xout,final_it,Xstore,Ystore,EG,art_method,stoprule,A,b,Kin,x0,y0,order,slack,options)
		!ART  General interface for all Kaczmarz/ART methods
		!
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin)
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin,x0)
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin,x0,y0)
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin,x0,y0,order)
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin,x0,y0,order,slack)
		!   call art(xout,Xstore,Ystore,art_method,stoprule,A,b,Kin,x0,y0,order,slack,options)
		!
		! Implements a general ART method for the linear system Ax = b:
		!       
		!       x^{k+1} = x^k + relaxpar*(b_i - a_i'*x^k)/(||a_i||_2^2)*a_i
		!
		! where a_i' is the i-th row of A, and the order for i is chosen in 
		! different ways by the provided methods (see kaczmarz, symkaczmarz, 
		! randkaczmarz for details) or a custom order specified by the user.
		! One iteration consists of m such steps (m-1 for symkaczmarz).
		!
		! Input:
		!   art_method  Either one of the strings 'kaczmarz', 'symkaczmarz', or 
		!               'randkaczmarz' to specify one of the provided methods;
		!               default is 'kaczmarz'.
		!	stoprule	Either one of the strings 'errorgauge' or 'mutualstep';
		!				default is 'errorgauge'.
		!   A           m times n matrix, or a function that implements matrix-
		!               vector multiplication with A and A'; see explanation below.
		!   b           m times 1 vector containing the right-hand side.
		!   Kin         Number of iterations. If K is a scalar, then K is the 
		!               maximum number of iterations and only the last iterate is 
		!               returned. If K is a vector, then max(K) is the maximum
		!               number of iterations and only iterates corresponding to the
		!               values in K are returned, together with the last iterate.
		!   x0          n times 1 starting vector. Default: x0 = 0.
		!   y0          n times 1 starting vector. Default: y0 = 0.
		!	order		A specific order to go through the rows. Default order is
		!				cyclic 1,2,...,m.
		!	slack		Number of iterations the twin algorithm keeps going after
		!				finding a minimum in the error gauge.
		!   options     Struct with the following fields:
		!      relaxpar The relaxation parameter. If relaxpar is a scalar < 2 then
		!               it is used in each iteration; default value is 1.
		!               Alternatively, relaxpar can be a function with a diminishing
		!               parameter, e.g., @(j) 1/sqrt(j), where j counts the total
		!               number of row updates.
		!      lbound   Lower bound in box constraint [lbound,ubound]. If scalar,
		!               this value is enforced on all elements of x in each 
		!               iteration. If vector, it must have same size as x and 
		!               then enforces elementwise lower bounds on x. If empty, no
		!               bound is enforced. +/-Inf can be used.
		!      ubound   Upper bound in box constraint [lbound,ubound]. If scalar,
		!               this value is enforced on all elements of x in each 
		!               iteration. If vector, it must have same size as x and 
		!               then enforces elementwise lower bounds on x. If empty, no
		!               bound is enforced. +/-Inf can be used.
		!      damp     A parameter damp to avoid division by very small row norms
		!               by adding damp*max_i{||a_i||_2^2} to ||a_i||_2^2.
		!      verbose  Nonnegative integer specifying whether progress is printed
		!               to screen during iterations. Default=0: no info printed.
		!               1: Print in every iteration. Larger than 1: Print every
		!               verbose'th iteration and first and last.
		!
		! Output:
		!	xout		Final output of the algorithms
		!	final_it	Final iteration number
		! Xstore,Ystore Matrices containing the saved iterations as the columns.
		!	EG			Error gauge for each iteration.
		
		implicit none
		!inputs
		character(*) :: art_method,stoprule
		class(csr_matrix),intent(in) :: A
		real(kind=8),intent(in) :: b(A%m)
		integer,intent(in) :: Kin(:)
		real(kind=8),intent(in),optional :: x0(A%n),y0(A%n)
		integer,optional :: order(:),slack
		class(opts),intent(in),optional :: options
		
		!outputs
		real(kind=8),intent(out) :: xout(A%n)
		integer :: final_it
		real(kind=8),dimension(A%n,size(Kin)),intent(out) :: Xstore,Ystore
		real(kind=8),intent(out) :: EG(maxval(Kin))
		
		!local variables
		integer :: ii,i,j,k,l,iter,maxiter,startiter,next_record,nord,ri
		integer :: Kuse(size(Kin)),slck,trck
		real(kind=8),dimension(A%n) :: x,y,v,w,u,z
		real(kind=8) :: ub,lb,dmp,omg,cumsum,minEG
		logical :: verb,upb,lwrb,isrand,rowprob
		integer :: ord(A%m),tempord(A%m)
		real(kind=8),dimension(A%m) :: nrA,cumul
		real(kind=8) :: rhs(2),det,mat(3),alph,bet
		logical :: perctdone(9)
		
		if (present(slack)) then
			slck = slack
		else
			slck = 7
		endif
		
		if (present(options)) then
			upb = options%upb
			lwrb = options%lwrb
			ub = options%ubound
			lb = options%lbound
			dmp = options%damp
			omg = options%relaxpar
			rowprob = options%rowprobability
			verb = options%verbose
			
			if (verb) then
				write(*,*) "================================"
				write(*,*) " "
				write(*,*) "TWIN ART"
				write(*,*) " "
				write(*,*) "================================"
				write(*,*) "Options detected:"
				write(*,*) "Verbose mode activated. I'll do a lot of talking."
				if (upb) then
					write(*,*) "Using upper bound: ",ub
				else
					write(*,*) "No upper bound present."
				endif
				if (lwrb) then
					write(*,*) "Using lower bound: ",lb
				else
					write(*,*) "No lower bound present."
				endif
				if (dmp==0) then
					write(*,*) "No damping coefficient present."
				else
					write(*,*) "Damping coefficient: ",dmp
				endif
				if (omg<=0 .or. omg >=2) stop "Relaxation parameter needs to be in (0,2)."
				if (omg==1d0) then
					write(*,*) "Using standard relaxation parameter: 1."
				else
					write(*,*) "Using relaxation parameter: ",omg
				endif
				write(*,*) "================================"
			endif
		else
			!standard options means no bounds applied, no damping and a relaxation parameter of 1.
			upb = .false.
			lwrb = .false.
			ub = 0d0
			lb = 1d0
			dmp = 0d0
			omg = 1d0
			rowprob = .false.
			verb = .false.
		endif
		
		isrand = .false.
		Kuse = Kin
		select case (trim(art_method))
			case ('kaczmarz')
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using Kaczmarz and custom ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				else
					nord = A%m
					ord(1:A%m) = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using Kaczmarz and standard sequential ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				endif
			case ('symkaczmarz')
				if (any(mod(Kin,2)==1)) stop "Kin should all be multiples of 2."
				Kuse = Kuse/2
				
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using symmetric Kaczmarz and custom ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				else
					nord = A%m
					ord = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using symmetric Kaczmarz and standard sequential ordering."
					if (present(options) .and. rowprob) write(*,*) "Ignoring row probability selection."
				endif
			case ('randkaczmarz')
				isrand = .true.
				!Ord needs to be present as it implements the skipping of zero-norm rows.
				if (present(order)) then
					nord = size(order)
					ord(1:nord) = order
					if (verb) write(*,*) "Using random Kaczmarz; ignoring order."
				else
					nord = A%m
					ord = (/ (i,i=1,A%m) /)
					if (verb) write(*,*) "Using random Kaczmarz."
				endif
				
				if (verb .and. rowprob) then
					write(*,*) "Using fancy row probability selection,",&
					" i.e., probability to select rows is proportional to row norm."
				 else
					 write(*,*) "Using uniform row probabilities."
				 endif
			case default
				stop "Unknown ART method. Use kaczmarz, symkaczmarz or randkaczmarz."
		end select
		
		select case (trim(stoprule))
			case ('errorgauge')
				if (verb) write(*,*) "Using errorgauge as a stopping rule."
			case ('mutualstep')
				if (verb) write(*,*) "Using mutualstep as a stopping rule."
			case default
				stop "Unknown stopping rule "//trim(stoprule)//". Use errorgauge or mutualstep."
		end select
		
		if (verb) write(*,*) "================================"
		
		!Apply damping.
		nrA = A%row_norms(2)
		nrA = nrA**2
		nrA = nrA + dmp*maxval(nrA)
		
		!select only nonzero rows
		if (verb) write(*,*) "Selecting nonzero rows..."
		l = nord
		nord = 0
		do i=1,l
			if (nrA(ord(i))>0) then !perhaps some tolerance should be used
				nord = nord+1
				tempord(nord) = ord(i)
			endif
		enddo
		ord(1:nord) = tempord(1:nord)

		
		if (verb) write(*,*) "done!"
		if (verb) write(*,*) "================================"
		
		
		!if random kaczmarz was chosen as the method
		if (isrand) then
			ord = scramble(nord)
			if (rowprob) then
				if (verb) write(*,*) "Calcultating row probabilities..."
				!row selection probability is proportional to row norms
				
				cumsum = 0
				do i=1,nord
					cumsum = cumsum+nrA(i)
					cumul(i) = cumsum
				enddo
				cumul = cumul/sum(nrA)
				cumul(nord+1:A%m) = 1d0
				if (verb) write(*,*) "done!"
				if (verb) write(*,*) "================================"
			endif
		endif
		
		if (present(x0) .and. .not. present(y0)) then
			x = x0
			y = x0
			
			if (verb) write(*,*) "Single initial condition provided, performing sweeps first."
			
			call basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(x,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			call basic_sweep(y,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(y,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			EG(1) = norm2(x-y)
			
			startiter = 2
		elseif (.not. present(x0) .and. present(y0)) then
			x = y0
			y = y0
			
			if (verb) write(*,*) "Single initial condition provided, performing sweeps first."
			
			call basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(x,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			call basic_sweep(y,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(y,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			EG(1) = norm2(x-y)
			
			startiter = 2
		elseif (present(x0) .and. present(y0)) then
			x = x0
			y = y0
			if (norm2(x0-y0) < 1d-10 ) stop "Initial conditions x0 and y0 too close."
			if (verb) write(*,*) "Using custom initial condition for both x and y."
			startiter = 1
		else
			x = 0
			y = 0 
			
			if (verb) write(*,*) "No initial conditions provided, performing sweeps first."
			
			call basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(x,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			call basic_sweep(y,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(y,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			
			EG(1) = norm2(x-y)
			
			startiter = 2
		endif
		
		maxiter = maxval(Kuse)
		l = 1
		minEG = huge(minEG)
		do iter=startiter,maxiter
			
			!scramble the rows if random
			if (isrand) ord(1:nord) = scramble(nord)
			
			!go through the nonzero rows
			if (trim(stoprule) == 'mutualstep') v = x !record previous iteration
			call basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(x,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			if (trim(stoprule) == 'mutualstep') then
				v = x-v !Record search direction
				x = x-v !Revert to previous iteration
			endif
			
			if (trim(stoprule) == 'mutualstep') w = y
			call basic_sweep(y,A,b,ord(nord:1:-1),nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			if (trim(art_method) == 'symkaczmarz') then
				call basic_sweep(y,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
			endif
			if (trim(stoprule) == 'mutualstep') then
				w = y-w
				y = y-w
				
			endif
			
			if (trim(stoprule) == 'errorgauge') then
				EG(iter) = norm2(x-y)
			
				if (EG(iter) < minEG) then
					minEG = EG(iter)
					v = x
					w = y
					trck = 0
					final_it = iter
				else
					trck = trck +1
					if (trck >=slck) exit
				endif
			elseif (trim(stoprule) == 'mutualstep') then
				w = w/norm2(w)
				v = v/norm2(v)
				
				!v = u
				!w = z - dot_product(z,u)*u
				
				
				mat = (/ dot_product(v,v), dot_product(w,w) , -dot_product(v,w) /)
				rhs = (/ -dot_product(v,x-y) , dot_product(w,x-y)  /)
				det = mat(1)*mat(2)-mat(3)**2
				alph = (rhs(1)*mat(2) - rhs(2)*mat(3) )/det
				bet = ( mat(1)*rhs(2) - mat(3)*rhs(1) )/det
				
				x = x + alph*v
				y = y + bet*w
				EG(iter) = norm2(x-y)
				
				if (verb) write(*,*) "(alpha,beta): (",alph,bet,")"
				
				if (abs(alph) + abs(bet) < 1d-5) exit
				
			endif
				
			
			if (iter==Kuse(l)) then
				Xstore(:,l) = x
				Ystore(:,l) = y
				l = l+1
				if (verb) then
					if (trim(art_method)=='symkaczmarz') then
						write(*,*) "Storing iteration: ",2*iter," out of ",2*maxiter
					else
						write(*,*) "Storing iteration: ",iter," out of ",maxiter
					endif
				endif
			else
				if (verb) then
					if (trim(art_method)=='symkaczmarz') then
						write(*,*) "Iteration: ",2*iter," out of ",2*maxiter
					else
						write(*,*) "Iteration: ",iter," out of ",maxiter
					endif
				endif
			endif
			
		enddo
		
		if (verb) then
			write(*,*) "================================"
			write(*,*) "Computing output solution..."
		endif
		
		if (trim(stoprule) == 'errorgauge') then
			xout = .5d0*(v+w)
		elseif (trim(stoprule) == 'mutualstep') then
			xout = .5d0*(x+y)
			final_it = iter
		endif
		
		
		if (verb) then
			write(*,*) "done!"
			write(*,*) "================================"
			write(*,*) "Final iteration: ",final_it
			write(*,*) "================================"
			write(*,*) "Verbose signing off."
			write(*,*) "================================"
			write(*,*) " "
			write(*,*) "ART"
			write(*,*) " "
			write(*,*) "================================"
		endif
	end subroutine twin_art
	
	subroutine basic_sweep(x,A,b,ord,nord,isrand,rowprob,omg,nrA,cumul,upb,lwrb,ub,lb)
		implicit none
		class(csr_matrix),intent(in) :: A
		real(kind=8),intent(inout) :: x(A%n)
		real(kind=8),intent(in) :: b(A%m),cumul(A%m),nrA(A%m),omg,ub,lb
		integer,intent(in) :: ord(2*A%m),nord
		logical,intent(in) :: isrand,rowprob,upb,lwrb
		
		!local
		real(kind=8) :: u,inp,coeff
		integer :: ri,ii,i
		
		do ii=1,nord
			i = ord(ii)

			if (isrand .and. rowprob) then
				!select fancy row
				call random_number(u)
				ri = 1
				!this could be sped up by using binary search
				do while (cumul(ri) < u)
					ri = ri+1
				enddo
			else
				ri = i
			endif

			!The actual ART step
			inp = A%rowvecmult(ri,x)
			coeff = ( omg*( b(ri) - inp )/nrA(ri) )
		
			!The bounds application is built into row addition the subroutine.
			!The rationale is as follows: when using a sparse matrix,
			!only the elements in the row that are nonzero are added,
			!so also only those will exceed the bounds.
			if (upb .and. .not. lwrb) then
				 call A%add_scalar_mult_row(ri,coeff,x,upb=ub)
			 elseif (.not. upb .and. lwrb) then
				 call A%add_scalar_mult_row(ri,coeff,x,lwrb=lb)
			 elseif (upb .and. lwrb) then
				 call A%add_scalar_mult_row(ri,coeff,x,lwrb=lb,upb=ub)
			 else
				 call A%add_scalar_mult_row(ri,coeff,x)
			 endif
		enddo
	end subroutine basic_sweep
	
	subroutine gradient_descent(xout,A,b,max_iter,x0,rho_in,tolerance,options)
		implicit none
		class(sparse_matrix),intent(in) :: A
		real(kind=8),intent(in) :: b(A%m)
		real(kind=8),optional :: x0(A%n),rho_in,tolerance
		integer,optional :: max_iter
		class(opts),optional :: options
		real(kind=8),intent(out) :: xout(A%n)
		
		!local
		integer :: it,maxiter
		real(kind=8) :: sigma,rho
		real(kind=8) :: r(A%m),xn(A%n)
		logical :: upb,lwrb,verb
		real(kind=8) :: ub,lb,nb,tol
				
		if (present(options)) then
			upb = options%upb
			lwrb = options%lwrb
			ub = options%ubound
			lb = options%lbound
			verb = options%verbose
			
			if (verb) then
				write(*,*) "================================"
				write(*,*) " "
				write(*,*) "GRADIENT DESCENT"
				write(*,*) " "
				write(*,*) "================================"
				write(*,*) "Options detected:"
				write(*,*) "Verbose mode activated. I'll do a lot of talking."
				if (upb) then
					write(*,*) "Using upper bound: ",ub
				else
					write(*,*) "No upper bound present."
				endif
				if (lwrb) then
					write(*,*) "Using lower bound: ",lb
				else
					write(*,*) "No lower bound present."
				endif
			endif
		else
			upb = .false.
			lwrb = .false.
			verb = .false.
		endif
		
		if (present(max_iter)) then
			maxiter = max_iter
		else
			maxiter = 100
		endif
		
		if (present(rho_in)) then
			rho = rho_in
		else
			rho = 1.5d0
		endif
		
		if (present(x0)) then
			xout = x0
			if (verb) write(*,*) "Using custom initial conditions."
		else
			xout = 0d0
		endif
		
		if (present(tolerance)) then
			tol = tolerance
			if (tol<0d0) stop "Must use positive tolerance."
		else
			tol = 1d-6
		endif
		
		if (verb) then
			write(*,*) "Using rho: ",rho
			if (tol==0) then
				write(*,*) "No tolerance set."
			else
				write(*,*) "Using tolerance: ",tol
			endif
			write(*,*) "================================"
			write(*,*) "Estimating matrix 2-norm..."
		endif
		
		sigma = normest(A)
		
		if (verb) then
			write(*,*) "done!"
			write(*,*) "Matrix norm: ",sigma
			write(*,*) "================================"
		endif
		
		sigma = sigma**2
		r = -b
		nb = norm2(b)**2
		
		if (verb) write(*,*) "Initial objective value: ",norm2(r)**2/nb
		
		do it=1,maxiter			
			xn = xout - rho/sigma*A%transvecmul(r)
			
			!apply bounds
			if (upb .and. .not. lwrb) then
				xn = min(xn,ub)
			elseif (.not. upb .and. lwrb) then
				xn = max(xn,lb)
			elseif (upb .and. lwrb) then
				xn = min(xn,ub)
				xn = max(xn,lb)
			endif
			
			if (norm2(xn-xout) < tol) then
				if (verb) write(*,*) "Tolerance reached, exiting."
				exit
			endif
			
			xout = xn
			
			r = A%matvecmul(xout)-b
			
			if (verb) write(*,*) "Objective value at iteration ",it," of ",maxiter,": ",norm2(r)**2/nb
		enddo
		
		if (verb) then
			write(*,*) "done!"
			write(*,*) "================================"
			write(*,*) "Verbose signing off."
			write(*,*) "================================"
			write(*,*) " "
			write(*,*) "GRADIENT DESCENT"
			write(*,*) " "
			write(*,*) "================================"
		endif
		
	end subroutine gradient_descent
	
	subroutine sirt(xout,sirt_method,A,b,Kin,x0,options)
		implicit none
		!inputs
		character(*) :: sirt_method
		class(sparse_matrix) :: A
		real(kind=8),intent(in) :: b(A%m)
		integer,intent(in) :: Kin(:)
		real(kind=8),intent(in),optional :: x0(A%n)
		class(opts),intent(in),optional :: options
		
		!local
		integer :: i,maxits
		real(kind=8),dimension(A%m) :: Mdiag,r
		real(kind=8),dimension(A%n)::x,xn,xout,Ddiag
		real(kind=8) :: ub,lb
		logical :: verb,upb,lwrb
		
		maxits = maxval(Kin)
		if (present(options)) then
			upb = options%upb
			lwrb = options%lwrb
			ub = options%ubound
			lb = options%lbound
			!dmp = options%damp
			!omg = options%relaxpar
			verb = options%verbose
			
			if (verb) then
				write(*,*) "================================"
				write(*,*) " "
				write(*,*) "SIRT"
				write(*,*) " "
				write(*,*) "================================"
				write(*,*) "Options detected:"
				write(*,*) "Verbose mode activated. I'll do a lot of talking."
				if (upb) then
					write(*,*) "Using upper bound: ",ub
				else
					write(*,*) "No upper bound present."
				endif
				if (lwrb) then
					write(*,*) "Using lower bound: ",lb
				else
					write(*,*) "No lower bound present."
				endif
				write(*,*) "================================"
			endif
		else
			!standard options means no bounds applied, no damping and a relaxation parameter of 1.
			upb = .false.
			lwrb = .false.
			ub = 0d0
			lb = 1d0
			verb = .false.
		endif
		
		if (present(x0)) then
			if (size(x0) /= A%n) stop "Size of x0 should be equal to n (#columns)."
			x = x0
			if (verb) write(*,*) "Using custom initial condition."
		else
			x = 0d0
			if (verb) write(*,*) "Using all-zero initial condition."
		endif
		
		
		
		select case(sirt_method)
			case ('cimmino')
				Mdiag = A%m * (A%row_norms(2))**2
				Ddiag = 1d0
				if (verb) write(*,*) "Using Cimmino's method."
			case ('sirt')
				Mdiag = A%row_norms(1)
				call A%transpose()
				Ddiag = A%row_norms(1)
				call A%transpose()
				if (verb) write(*,*) "Using SIRT method."
			case default
				stop "Not a recognised SIRT method."
		end select
		
		
		if(verb) write(*,*) "================================"
		
		r = 0d0
		do i = 1,maxits
			xn = x + A%transvecmul( (b - A%matvecmul(x))/Mdiag )/Ddiag
			
			!apply bounds
			if (upb .and. .not. lwrb) then
				xn = min(xn,ub)
			elseif (.not. upb .and. lwrb) then
				xn = max(xn,lb)
			elseif (upb .and. lwrb) then
				xn = min(xn,ub)
				xn = max(xn,lb)
			endif
			
			x = xn
			
			
			if (verb) write(*,*) "Iteration ",i," of ",maxits
			
			
			!TO ADD: record the iterations specified in Kin
			
		enddo
		
		xout = x
		
		
		if (verb) then
			write(*,*) "================================"
			write(*,*) "Verbose signing off."
			write(*,*) "================================"
			write(*,*) " "
			write(*,*) "SIRT"
			write(*,*) " "
			write(*,*) "================================"
		endif
		
		
		
	end subroutine sirt
	
	function normest(A,tolerance,max_its) result(sigma)
		!Use the power method to find the largest singular value of A.
		!Brute force, but it works.
		!Also finds the singular vector, but that one is discarded.

		implicit none
		class(sparse_matrix) :: A
		real(kind=8),optional :: tolerance
		integer, optional :: max_its
		real(kind=8) :: sigma

		!local variables
		real(kind=8) :: vn(A%n),v(A%n),w(A%m)
		real(kind=8) :: tol
		integer :: maxit,i

		if (present(tolerance)) then
			tol = tolerance
		else
			tol = 1d-5 !need only a crude estimation
		endif

		if (present(max_its)) then
			maxit = max_its
		else
			maxit = 100
		endif

		call random_number(v) !set up a random vector

		do i=1,maxit
			w = A%matvecmul(v) !w = A*v
			vn = A%transvecmul(w) !vn = A'*A*v
			vn = vn/norm2(vn) !normalise

			if ( norm2(vn-v) < tol ) exit

			v = vn
		enddo

		w = A%matvecmul(v)

		sigma = sqrt(dot_product(w,w)/dot_product(v,v))

	end function normest
	
	function scramble( number_of_values ) result(array)
		!Got this one from the internet.
		!return integer array of random values 1 to N.
		implicit none
		integer,intent(in)    :: number_of_values
		integer,allocatable   :: array(:)
		integer               :: i, j, k, m, n
		integer               :: temp
		real(kind=8)          :: u

		array=[(i,i=1,number_of_values)]

		! The intrinsic RANDOM_NUMBER(3f) returns a real number (or an array
		! of such) from the uniform distribution over the interval [0,1). (ie.
		! it includes 0 but not 1.).
		!
		! To have a discrete uniform distribution on
		! the integers {n, n+1, ..., m-1, m} carve the continuous distribution
		! up into m+1-n equal sized chunks, mapping each chunk to an integer.
		!
		! One way is:
		!   call random_number(u)
		!   j = n + FLOOR((m+1-n)*u)  ! choose one from m-n+1 integers

		n=1
		m=number_of_values
		do k=1,2
		 do i=1,m
		    call random_number(u)
		    j = n + FLOOR((m+1-n)*u)
		    ! switch values
		    temp=array(j)
		    array(j)=array(i)
		    array(i)=temp
		 enddo
		enddo

	end function scramble
	
	
end module AIRMs
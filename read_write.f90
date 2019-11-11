module read_write
	implicit none
	
	type data
		integer :: n,p,nth
		real(kind=8) :: r,dw,sd
		real(kind=8),allocatable :: T(:),theta(:)
	contains
		procedure :: read => read_data
	end type data
	
contains
	
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
		
		
		read(10,*) trash,l
		if (l /= dat%nth*dat%p) stop "Data has the wrong size, or p and nth are incorrect."
	
		allocate(dat%T(l))
		do i=1,l
			read(10,*) dat%T(i)
		enddo
	
		close(10)
	end subroutine read_data
	
end module read_write
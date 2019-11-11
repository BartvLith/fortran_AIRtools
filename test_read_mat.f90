program test_read_mat
	implicit none
	integer :: n,p,nth,l,i
	real(kind=8) :: r,dw,sd
	real(kind=8),allocatable :: T(:)
	character(30) :: trash
	character(100) :: num1char
	
	if (COMMAND_ARGUMENT_COUNT()==0) stop "Supply input file."
	if (COMMAND_ARGUMENT_COUNT()>1) stop "Only 1 input file can be used."
	
	
	CALL GET_COMMAND_ARGUMENT(1,num1char)
	num1char = adjustl(num1char)
	
	if (num1char(1:2) == '-i') then
		write(*,*) "Opening file: ",trim(num1char(3:100))
	else
		stop "Input flag not recognised. Only -i is available."
	endif
	
	open (unit = 10, file = trim(num1char(3:100)))
	
	read(10,*) trash, n
	read(10,*) trash, p
	read(10,*) trash, nth
	read(10,*) trash, r
	read(10,*) trash, dw
	read(10,*) trash, sd
	
	read(10,*) trash,l
	
	if (l /= p*nth) stop "Measurements have the wrong size, or p and nth are incorrect."
	
	allocate(T(l))
	do i=1,l
		read(10,*) T(i)
	enddo
	
	close(10)
	
	
	
	write(*,*) n,p,nth,r,dw,sd,T(10)
	
end program test_read_mat
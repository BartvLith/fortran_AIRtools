module sparse_matrices
	implicit none
	
	real(kind=8),parameter :: TOL = 1d-15
	
	!This module implements some aspects of sparse matrices.
	
	
	!Compressed Sparse Row matrix format:
	!Say A is an m-by-n matrix with nnz nonzeros.
	!CSR stores A as three one-dimensional arrays.
	!vA is the component values scanned row-by-row starting at the top and from left to right.
	!vA has length nnz.
	!IA is an array of length m+1, with IA(0) = 0 and IA(k) is the number of nonzeros encountered
	!so far in the 	scanning process per row.
	!JA records the column indices of each component of vals. JA also has length nnz.
	
	!The Compressed Sparse Column matrix format is similar but runs through the columns.
	!As a consequence, the subroutine for CSC matrix-vector multiplication applied to a CSR matrix
	!results in an application of the transpose.
	
	
	!to add:
	!- infinity-norms for row norms
	!- column norms
	
	public sparse_matrix,CSR_matrix,CSC_matrix
	private print_line,csr_matrowvec,sparsematmul,sparsetransmul,MMRS,MMCS,&
			transposition,conv_to_CSC,conv_to_CSR,reorder,SM_set_element,SM_get_element,&
			csr_get_element,SM_delete_element,compute_row_norms,add_mult_row,scalar_mult,&
			compute_column_norms
	
	type sparse_matrix
		!This type is only a helper, it contains all the basic information that defines CSR and CSC matrices:
		!- vA, list of nonzero values
		!- IA, cumulative list of nonzeros per row/column (CSR/CSC)
		!- JA, index of the non-scan direction
		!- m,n size of the matrix
		!-nnz number of nonzeros
		
		!This type is private, and as such cannot be used outside this module.
		!The reason is that both CSR and CSC matrices submit to the basic format,
		!moreover, the matrix-vector multiplication algorithm of CSC matrices applied
		!to a CSR matrix is equivalent to computing the transpose matrix-vector multiplication.
		!Hence, for polymorphism reasons, this class carries all the defining attributes and procedures.
		
		integer :: m,n,nnz
		real(kind=8),allocatable :: vA(:)
		integer,allocatable :: IA(:),JA(:)
	contains
		!matrix-vector multiplication
		!	w = A%matvecmul(v)
		!w in R^m is A*v
		!v in R^n
		procedure :: matvecmul => sparsematmul 
		
		
		!multiply input vector by transpose
		!	v = A%transvecmul(w)
		!v in R^n is A'*v
		!w in R^m
		procedure :: transvecmul => sparsetransmul 
		
		
		!multiply a vector by row k (CSC) or column k (CSR) of A
		!	phi = A%rowvecmult(k,v), phi is scalar, A is type(csr_matrix), k in {1,A%m} and v in R^n
		!	phi = A%rowvecmult(k,v), phi is scalar, A is type(csc_matrix), k in {1,A%n} and v in R^m
		procedure :: rowvecmult => csr_matrowvec 
		
		
		!print a row (CSC) or a column (CSR)
		!	call A%print_part()
		procedure :: print_part => print_line
		
		
		!Transpose the matrix A
		!	call A%transpose()
		procedure :: transpose => transposition
		
		
		!Convert a matrix to CSR type
		!	B = A%convert_to_CSR()
		!A is any sparse matrix
		!B is type(csr_matrix)
		procedure :: convert_to_CSR => conv_to_CSR
		
		
		!Convert a matrix to CSC type
		!	B = A%convert_to_CSC()
		!A is any sparse matrix
		!B is type(csc_matrix)
		procedure :: convert_to_CSC => conv_to_CSC
		
		
		!Set a single element (i,j) to a certain value phi
		!	call A%set_element(i,j,phi)
		!A is type(csr_matrix) or type(csc_matrix)
		!i in {1,A%m},j in {1,A%n} are integers
		!phi is a scalar
		
		!If phi is set to 0d0, the element is deleted.
		!In fact, delete_element calls this subroutine.
		
		!NOTE:
		!==========================================================================
		!Never use this procedure unless absolutely unavoidable, it's very slow.
		!The vA and JA arrays have to be reallocated to extend by one element.
		!==========================================================================
		procedure :: set_element => SM_set_element
		
		
		!Get a single element (i,j)
		!	phi = A%get_element(i,j)
		!A is type(csr_matrix) or type(csc_matrix)
		!i in {1,A%m},j in {1,A%n} are integers
		!phi is a scalar
		procedure :: get_element => SM_get_element
		
		
		!Delete a single element (i,j)
		!	call A%delete_element(i,j)
		!A is a type(csr_matrix) or type(csc_matrix)
		!i in {1,A%m},j in {1,A%n} are integers
		
		!NOTE:
		!==========================================================================
		!Never use this procedure unless absolutely unavoidable, it's very slow.
		!The vA and JA arrays have to be reallocated to extend by one element.
		!==========================================================================
		procedure :: delete_element => SM_delete_element
		
		!Compute the p-norms of each row
		!	nrA = A%row_norms(p)
		!A is type(csr_matrix) or type(csc_matrix)
		!nrA is in R^m; the vector of row norms
		!p is an integer, cannot be infinite
		procedure :: row_norms => compute_row_norms
		
		!Compute the p-norms of each column
		!	ncA = A%column_norms(p)
		!A is type(csr_matrix) or type(csc_matrix)
		!ncA is in R^n; the vector of column norms
		!p is an integer, cannot be infinite
		procedure :: column_norms => compute_column_norms
		
		!Multiply the whole matrix by a scalar phi
		!	call A%multiply_by_scalar(phi)
		!A is type(sparse_matrix)
		!phi is a sclaar
		procedure :: multiply_by_scalar => scalar_mult
		
	end type sparse_matrix
	
	!The following types are basically a flag to be able to tell what type of matrix it is.
	!As explained above, CSR and CSC matrices only differ in interpretation.
	
	type, extends(sparse_matrix) :: CSR_matrix
	contains
		!Add a scalar, phi, times the multiple of row k to an input vector x and bound
		!	call A%add_scalar_mult_row(k,phi,x)
		!	call A%add_scalar_mult_row(k,phi,x,lwrb = lb)
		!	call A%add_scalar_mult_row(k,phi,x,upb = ub)
		!	call A%add_scalar_mult_row(k,phi,x,lb,ub)
		!A is type(csr_matrix)
		!k is an integer in {1,A%m}
		!phi is a scalar
		!x is in R^n, and is intent(inout)
		
		!optional inputs:
		!lb is a scalar lower bound such that x = max(x,lb)
		!ub is an scalar upper bound such that x = min(x,ub)
		procedure :: add_scalar_mult_row => add_mult_row
	end type CSR_matrix

	type, extends(sparse_matrix) :: CSC_matrix
	end type CSC_matrix

contains
	
	subroutine construct_csr_matrix(m,n,vA,IA,JA,A)
		implicit none
		integer,intent(in) :: m,n
		real(kind=8),intent(in) :: vA(:)
		integer,intent(in) :: IA(:),JA(:)
		type(csr_matrix),intent(out) :: A

		A%n = n
		A%m = m
		A%nnz = IA(ubound(IA,1))
		allocate(A%vA(size(vA)))
		allocate(A%IA(0:m))
		allocate(A%JA(size(JA)))

		A%vA = vA
		A%IA = IA
		A%JA = JA
	end subroutine construct_csr_matrix
	
	subroutine construct_csc_matrix(m,n,vA,IA,JA,A)
		implicit none
		integer,intent(in) :: m,n
		real(kind=8),intent(in) :: vA(:)
		integer,intent(in) :: IA(:),JA(:)
		type(csc_matrix),intent(out) :: A

		A%n = n
		A%m = m
		A%nnz = IA(ubound(IA,1))
		allocate(A%vA(size(vA)))
		allocate(A%IA(0:m))
		allocate(A%JA(size(JA)))

		A%vA = vA
		A%IA = IA
		A%JA = JA
	end subroutine construct_csc_matrix
	
	subroutine print_line(A,k)
		!prints row k of a CSR matrix
		!used for testing
		implicit none
		class(sparse_matrix),intent(in) :: A
		integer, intent(in) :: k

		!local variables
		integer :: i

		select type (A)
		    type is (sparse_matrix)
				stop "Storage format must be defined"
			class is (csr_matrix)
				if (A%IA(k)>A%IA(k-1)) then
					write(*,*) "row ",k
					do i=A%IA(k-1)+1,A%IA(k)
						write(*,*) "( ",k," ,", A%JA(i),"  ) ",A%vA(i)
					enddo
				else
					write(*,*) "empty row"
				endif
			class is (csc_matrix)
				if (A%IA(k)>A%IA(k-1)) then
					write(*,*) "column ",k
					do i=A%IA(k-1)+1,A%IA(k)
						write(*,*) "( ",A%JA(i)," ,",k,"  ) ",A%vA(i)
					enddo
				else
					write(*,*) "empty column"
				endif
		end select
	end subroutine print_line

	function csr_matrowvec(A,k,v) result(f)
		!multiply vector v by a single row of a sparse matrix in CSR format
		implicit none
		class(sparse_matrix),intent(in) :: A
		real(kind=8),intent(in) :: v(:)
		integer, intent(in) :: k
		real(kind=8) :: f
		integer :: i
		
		f = 0d0
		if (A%IA(k)>A%IA(k-1)) then !only if the row contains any nonzeros
			do i=A%IA(k-1)+1,A%IA(k)!loop through nonzeros
				!multiply only the nonzero components
				f = f + v(A%JA(i))*A%vA(i)

				!write(*,*) "( ",k," ,", A%JA(i),"  ) ",A%vAi),v(A%JA(i))
			enddo
		endif
	end function csr_matrowvec
	
	function sparsematmul(A,v) result(w)
		!multiply vector v by a sparse matrix
		implicit none
		class(sparse_matrix),intent(in) :: A
		real(kind=8), intent(in) :: v(:)
		real(kind=8) :: w(A%m)
		
		select type (A)
		    type is (sparse_matrix)
				stop "Matrix type must be either CSR or CSC."
		    class is (csr_matrix)
				!write(*,*) "CSR matrix detected"
				if (size(v)==A%n) then
					w = MMRS(A,v)
				else
					stop "Vector input size incorrect."
				endif
			class is (csc_matrix)
				!write(*,*) "CSC matrix detected"
				if (size(v)==A%n) then
					w = MMCS(A,v)
				else
					stop "Vector input size incorrect."
				endif
		end select
		
	end function sparsematmul
	
	function sparsetransmul(A,v) result(w)
		!multiply vector v by a sparse matrix
		implicit none
		class(sparse_matrix) :: A
		real(kind=8), intent(in) :: v(:)
		real(kind=8) :: w(A%n)
		
		!local
		integer :: exec,temp
		
		select type (A)
		    type is (sparse_matrix)
				stop "Matrix type must be either CSR or CSC."
		    class is (csr_matrix)
				!write(*,*) "CSR matrix detected"
				if (size(v)==A%m) then
					exec = 1
				else
					stop "Vector input size incorrect."
				endif
			class is (csc_matrix)
				!write(*,*) "CSC matrix detected"
				if (size(v)==A%m) then
					exec = 2
				else
					stop "Vector input size incorrect."
				endif
		end select
		
		if (exec == 1) then
			temp = A%n
			A%n = A%m
			A%m = temp
			w = MMCS(A,v)
			temp = A%n
			A%n = A%m
			A%m = temp
		elseif (exec == 2) then
			temp = A%n
			A%n = A%m
			A%m = temp
			w = MMRS(A,v)
			temp = A%n
			A%n = A%m
			A%m = temp
		endif
		
	end function sparsetransmul

	function MMRS(A,v) result(w)
		!multiply vector v by a sparse matrix
		implicit none
		class(sparse_matrix),intent(in) :: A
		real(kind=8), intent(in) :: v(A%n)
		real(kind=8) :: w(A%m)

		!local variables
		integer :: k,i

		w = 0d0 !nullify everything
		
		!$omp parallel shared(w)
		!$OMP do
		do k=1,A%m !loop through the rows
			w(k) = A%rowvecmult(k,v)
		enddo
		!$OMP end do
		!$omp end parallel
	end function MMRS
	
	function MMCS(A,v) result(w)
		!multiply vector v by a sparse matrix
		implicit none
		class(sparse_matrix),intent(in) :: A
		real(kind=8), intent(in) :: v(A%n)
		real(kind=8) :: w(A%m)

		!local variables
		integer :: k,i

		w = 0d0 !nullify everything
		
		!$omp parallel shared(w)
		!$omp do
		do k=1,A%n !run through the components of v
			if (A%IA(k)>A%IA(k-1)) then !only if the column contains any nonzeros
				do i=A%IA(k-1)+1,A%IA(k)!loop through nonzeros
					w(A%JA(i)) = w(A%JA(i)) + v(k)*A%vA(i)
				enddo
			endif
		enddo
		!$omp end do
		!$omp end parallel
	end function MMCS	
	
	subroutine transposition(A)
		implicit none
		class(sparse_matrix), intent(inout) :: A
		
		integer :: temp
		type(sparse_matrix) :: B
		
		B = reorder(A)
		temp = A%m
		A%m = A%n
		A%n = temp
		call move_alloc(B%IA,A%IA)
		call move_alloc(B%JA,A%JA)
		call move_alloc(B%vA,A%vA)
		A%nnz = B%nnz
		
	end subroutine transposition
	
	function conv_to_CSC(A) result(B)
		implicit none
		class(sparse_matrix) :: A
		type(CSC_matrix) :: B
			
		logical :: reord
		type(sparse_matrix) :: C
			
		select type (A)
		    type is (sparse_matrix)
				reord = .false.
		    class is (csr_matrix)
				reord = .true.
			class is (csc_matrix)
				reord = .false.
				B = A
				return
		end select
		
		if (reord) then 
			C = reorder(A)
			B%m = C%m
			B%n = C%n
			call move_alloc(C%IA,B%IA)
			call move_alloc(C%JA,B%JA)
			allocate(B%vA(size(A%va)))
			B%vA = C%va
			B%nnz = C%nnz
		else
			B%m = A%m
			B%n = A%n
			call move_alloc(A%IA,B%IA)
			call move_alloc(A%JA,B%JA)
			allocate(A%vA(size(A%va)))
			B%vA = A%va
			B%nnz = A%nnz
		endif
		
	end function conv_to_CSC
	
	function conv_to_CSR(A) result(B)
		implicit none
		class(sparse_matrix) :: A
		type(CSR_matrix) :: B
			
		type(sparse_matrix) :: C
		logical :: reord
			
		select type (A)
		    type is (sparse_matrix)
				reord = .false.
		    class is (csr_matrix)
				reord = .false.
				B = A
				return
			class is (csc_matrix)
				reord = .true.
		end select
		
		if (reord) then 
			C = reorder(A)
			B%m = C%m
			B%n = C%n
			call move_alloc(C%IA,B%IA)
			call move_alloc(C%JA,B%JA)
			allocate(C%vA(size(A%va)))
			B%vA = C%va
			B%nnz = C%nnz
		else
			B%m = A%m
			B%n = A%n
			call move_alloc(A%IA,B%IA)
			call move_alloc(A%JA,B%JA)
			allocate(A%vA(size(A%va)))
			B%vA = A%va
			B%nnz = A%nnz
		endif
	end function conv_to_CSR
	
	function reorder(A) result(B)
		!This subroutine converts matrices between types CSR and CSC.
		!However, given the relation between CSR and CSC types, we have the following:

		!If A is a CSR matrix, then B is:
		! 1) the CSC form of A
		! 2) the CSR form of A^T

		!If A is a CSC matrix, then B is:
		! 1) the CSR form of A
		! 2) the CSC form of A^T

		implicit none
		class(sparse_matrix),intent(in) :: A
		type(sparse_matrix) :: B

		!local variables
		integer :: i,j,k,nnz,nr,nc,t,cumsum,col,row,dest,last

		!find the size of the matrix
		nr = A%m
		nc = A%n
		nnz = A%nnz
		
		
		B = A
		
		deallocate(B%IA)
		
		if (size(A%IA)==(A%m+1)) then
			allocate(B%IA(0:A%n))
		elseif (size(A%IA)==(A%n+1)) then
			allocate(B%IA(0:A%m))
		else
			stop "Wrong array sizes in reorder."
		endif
			
		
		
		!count the number of nonzeros in each column
		B%IA = 0
		do i = 1,nnz
			B%IA(A%JA(i)) = B%IA(A%JA(i))+1
		enddo

		!IB is the cumulative sum of the row nnz's
		cumsum = 0
		do i=0,nc
			t = B%IA(i)
			cumsum = cumsum + t
			B%IA(i) = cumsum
		enddo


		B%vA = 0
		B%JA = 0
		do i=1,nr
			do j=A%IA(i-1)+1,A%IA(i)
				col = A%JA(j)-1 !find correct column
				dest = B%IA(col)+1 !find correct place in array

				B%JA(dest) = i !record
				B%vA(dest) = A%vA(j)

				B%IA(col) = B%IA(col)+1 !update IB so that every value gets a unique spot
			enddo
		enddo


		!At this point, IB is shifted one location to the left, while the last element is still nnz.
		!We have to shift each element one to the right, setting the first element to 0.
		last = 0
		do col=0,nc-1
			t = B%IA(col)
			B%IA(col) = last
			last = t
		enddo

	end function reorder
	
	subroutine SM_set_element(A,i,j,aij)
		implicit none
		class(sparse_matrix),intent(inout) ::A
		integer,intent(in) :: i,j
		real(kind=8),intent(in) :: aij
		
		integer :: local_type
		
		select type (A)
		    type is (sparse_matrix)
				stop "Storage format must be known, use only on CSR or CSC matrices."
		    class is (csr_matrix)
				local_type = 1
			class is (csc_matrix)
				local_type = 2
		end select
		
		if (local_type == 1) then
			call csr_set_element(A,i,j,aij)
		elseif (local_type == 2) then
			call csr_set_element(A,j,i,aij)
		endif
	end subroutine SM_set_element
	
	subroutine SM_delete_element(A,i,j)
		implicit none
		class(sparse_matrix),intent(inout) ::A
		integer,intent(in) :: i,j
		
		call A%set_element(i,j,0d0)
	end subroutine SM_delete_element
	
	subroutine csr_set_element(A,i,j,aij)
		implicit none
		class(sparse_matrix),intent(inout) ::A
		integer,intent(in) :: i,j
		real(kind=8),intent(in) :: aij

		!local variables
		logical :: present
		integer :: k,l,nnz
		integer, allocatable ::jtemp(:)
		real(kind=8),allocatable :: atemp(:)

		nnz = A%nnz
		
		!NOTE: this function lacks the possibility that elements are set to zero
		

		!check if a nonzero element already exists in position i,j
		present = .false.
		do k=A%IA(i-1)+1,A%IA(i)
			if (A%JA(k)==j) then
				present = .true.
				exit
			endif
		enddo

		!if so, overwrite it
		if (present) then
			
			if (abs(aij)>TOL) then
				A%vA(k) = aij
			else !remove element k
				allocate(jtemp(1:nnz-1))
				allocate(atemp(1:nnz-1))
				if (k ==1 ) then
					!remove first element
					jtemp = A%JA(2:nnz)
					atemp = A%vA(2:nnz)
				elseif (k==nnz) then
					!remove last element
					jtemp = A%JA(1:nnz-1)
					atemp = A%vA(1:nnz-1)
				else
					!remove an element from the middle somewhere
					jtemp = A%JA([ (l,l=1,k-1) , (l,l=k+1,nnz) ])
					atemp = A%vA([ (l,l=1,k-1) , (l,l=k+1,nnz) ])
				endif
				call move_alloc(jtemp,A%JA )
				call move_alloc(atemp,A%vA )
				A%IA(i:ubound(A%IA,1)) = A%IA(i:ubound(A%IA,1))-1
				A%nnz = nnz-1
			endif
		else
			!if the element is not present and nonzero
			if (abs(aij)>TOL) then
				!This involves putting j in the correct place in arrays IA and valsA.
				!The array IA needs to be updated, add 1 to every entry starting from i


				k=A%IA(i-1)+1 !new position for a in the valsA and IA arrays

				!allocate a temporary array of size nnz+1
				allocate(jtemp(1:nnz+1))

				!If JA is nonempty copy the old contents,
				!leaving space for j.
				if (allocated(A%JA)) then
					jtemp(1:k-1) = A%JA(1:k-1)
					jtemp(k+1:nnz+1) = A%JA(k:nnz)
				endif
				call move_alloc(jtemp,A%JA) !reallocate JA and copy jtemp
				A%JA(k) = j !add the index j

				allocate(atemp(1:nnz+1))
				if (allocated(A%vA)) then
					atemp(1:k-1) = A%vA(1:k-1)
					atemp(k+1:nnz+1) = A%vA(k:nnz)
				endif
				call move_alloc(atemp,A%vA)
				A%vA(k) = aij

				A%IA(i:ubound(A%IA,1)) = A%IA(i:ubound(A%IA,1))+1
			
				A%nnz = A%nnz+1
			endif

		endif

	end subroutine csr_set_element

	function SM_get_element(A,i,j) result(aij)
		implicit none
		class(sparse_matrix),intent(in) :: A
		integer, intent(in) :: i,j
		real(kind=8) :: aij
		
		select type (A)
		    type is (sparse_matrix)
				stop "Storage format must be known, use only on CSR or CSC matrices."
		    class is (csr_matrix)
				aij = csr_get_element(A,i,j)
			class is (csc_matrix)
				aij = csr_get_element(A,j,i)
		end select
	end function SM_get_element

	function csr_get_element(A,i,j) result(aij)
		implicit none
		class(sparse_matrix) :: A
		integer,intent(in) :: i,j
		real(kind=8):: aij

		!local variables
		logical :: present
		integer :: k

		!Look through row i for a record of element j.
		!Since the list of column entries is potentially unordered,
		!a simple linear search is used.

		present = .false.
		do k=A%IA(i-1)+1,A%IA(i)
			if (A%JA(k)==j) then
				present = .true.
				exit
			endif
		enddo

		if (present) then
			aij = A%vA(k)
		else
			aij = 0d0
		endif

	end function csr_get_element
	
	function compute_row_norms(A,p) result(nrA)
		implicit none
		class(sparse_matrix) :: A
		integer,intent(in) :: p
		
		integer :: i,k
		real(kind=8) :: nrA(A%m)
		logical :: trans
		
		select type (A)
		    type is (sparse_matrix)
				stop "Storage format must be defined"
			class is (csr_matrix)
				trans = .false.
			class is (csc_matrix)
				trans = .true.
		end select
		
		if (trans) call A%transpose()
		
		
		nrA = 0d0
		!$omp parallel shared(nrA)
		!$omp do
		do k=1,A%m
			if (A%IA(k)>A%IA(k-1)) then !only if the row contains any nonzeros
				do i=A%IA(k-1)+1,A%IA(k)!loop through nonzeros
					!multiply only the nonzero components
					nrA(k) = nrA(k) + abs(A%vA(i))**p
				enddo
			endif
		enddo
		!$omp end do
		!$omp end parallel
		
		nrA = nrA**(1d0/p)
		
		if (trans) call A%transpose()
		
		
	end function compute_row_norms
	
	function compute_column_norms(A,p) result(ncA)
		implicit none
		class(sparse_matrix) :: A
		integer,intent(in) :: p
		real(kind=8) :: ncA(A%n)
		type(sparse_matrix) :: B
		
		call A%transpose()
		ncA = A%row_norms(p)
		call A%transpose()
		
	end function compute_column_norms
	
	subroutine add_mult_row(A,k,c,w,lwrb,upb)
		!add a scalar multiple, c, of the kth row to w
		!the row being very sparse, we should only adjust those elements
		!that need adjusting
		implicit none
		class(csr_matrix),intent(in) :: A
		integer, intent(in) :: k
		real(kind=8),intent(inout) :: w(A%n),c
		real(kind=8),intent(in),optional :: lwrb,upb

		!local variables
		integer :: i

		if (A%IA(k)>A%IA(k-1)) then !only if the row contains any nonzeros
			do i=A%IA(k-1)+1,A%IA(k)!loop through nonzeros
				!record nonzero components
				w(A%JA(i)) = w(A%JA(i)) + c*A%vA(i)
				if (present(lwrb)) then
					w(A%JA(i)) = max(w(A%JA(i)),lwrb)
				endif

				if (present(upb)) then
					w(A%JA(i)) = min(w(A%JA(i)),upb)
				endif
			enddo
		endif
		
	end subroutine add_mult_row
	
	subroutine scalar_mult(A,c)
		!multiply the whole matrix with scalar c
		implicit none
		class(sparse_matrix),intent(inout) :: A
		real(kind=8) :: c
		
		A%vA = c*A%vA
	end subroutine scalar_mult
	
end module sparse_matrices
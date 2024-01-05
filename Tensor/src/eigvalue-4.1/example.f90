program Name
	use Tools
	use mpi
	use Tensor_Tools
	use eigen_value_Tools
	implicit none
	type(eigenvalue)::Solver
	type(Tensor)::Res(2),parameter(2),phi, TMP
	integer::ncv ! number of lanzcos basis
	integer::num_ouput
	integer::print_num
	integer::N
	integer::Legi,dimi


	ncv=10
	print_num=10
	num_ouput=3
	N=50

	call Solver%set_ncv(ncv)
	call Solver%set_Hphi(H_phi)
	call Solver%set_print_num(print_num)
	call initialData(parameter,N)
	Res=Solver%eig(parameter,'SR',num_ouput,.true.)

	call Res(1)%print

	Legi=2 ! the col of the tensor
	dimi=1 ! the first col
	phi=Res(2)%getSubTensor(Legi,dimi) ! the col ith of Res(2) is the eigen state for eigen value ith


	TMP=(.H.phi)*parameter(1)*phi
	call TMP%print

	!************ use lapack subroutine to get all the eigen value


	Res=parameter(1)%eig()
	TMP=dble(Res(1))
	call TMP%sort()
	call TMP%print


contains

	subroutine H_phi(parameter,phi) !return H*|phi>, if phi is empty, generate phi
		type(Tensor)::parameter(:) ! It can store any thing, for example the data of H
		type(Tensor)::phi   ! input phi, output H*Phi
		type(Tensor):: H, phi0
		H=parameter(1)
		phi0=parameter(2)
		if(.not.phi%getFlag())then
			phi=phi0
		else
			phi=H*phi
		end if
		return
	end subroutine

	subroutine initialData(parameter,N)
		type(Tensor)::parameter(2)
		integer::N
		call parameter(1)%allocate([N,N],'real*8')
		call parameter(1)%random([-1d0,1d0])
		parameter(1)=(.H.parameter(1))*parameter(1)
		
		call parameter(2)%allocate([N],'real*8')
		call parameter(2)%random([-1d0,1d0])
	end subroutine

end

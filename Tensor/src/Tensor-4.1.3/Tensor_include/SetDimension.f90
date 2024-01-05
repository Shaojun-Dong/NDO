	!*********************************************************
	!
	!  reset dimension
	!
	!*********************************************************

	subroutine resetdim1(T,dimData)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dimData(:)
		if(T%getSymmetryFlag())then
			call writemess('Can not reset the Dimension in the Tensor',-1)
			call writemess('The Tensor is of symmetry but input a non-symmetry dimension',-1)
			call error_stop
		end if
		if(product(dimData).ne.T%getTotalData())then
			call writemess('ERROR in resetting dimension, dimension do not match',-1)
			call writemess(' input dimension:',-1)
			call writemess(dimData,-1)
			call T%diminfo(.true.)
			call error_stop
		end if
		T%Dimension=dimData
		return
	end subroutine

	subroutine resetdim2(T,dimData)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::dimData(:)
		type(Dimension)::NewDim
		if(.not.T%getSymmetryFlag())then
			call writemess('Can not reset the Dimension in the Tensor',-1)
			call writemess('The Tensor is not of symmetry but input a array of quantum number',-1)
			call error_stop
		end if
		NewDim=dimData
		if(product(NewDim%dim()).ne.T%getTotalData())then
			call writemess('ERROR in resetting dimension, dimension do not match',-1)
			call writemess(' input dimension:',-1)
			call NewDim%diminfo(.true.)
			call T%diminfo(.true.)
			call error_stop
		end if
		T%Dimension=dimData
		return
	end subroutine

	subroutine resetdim3(T,dimData)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dimData
		if(T%getSymmetryFlag().neqv.dimData%getSymmetryFlag())then
			call writemess('Can not reset the Dimension in the Tensor',-1)
			call writemess('T%getSymmetryFlag()='+T%getSymmetryFlag(),-1)
			call writemess('dimData%getSymmetryFlag()='+dimData%getSymmetryFlag(),-1)
			call error_stop
		end if
		T%Dimension=dimData
		return
	end subroutine
	subroutine enlargeTensor(T,newD,randomNumberScal)
		class(Tensor),intent(inout)::T
		integer,intent(in)::newD
		class(*),intent(in)::randomNumberScal(2)
		type(dimension)::dimen
		class(*),pointer::clp(:)
		integer::oldTotal,rank
		if(T%getSymmetryFlag())then
			call writemess('ERROR in calling T%enlarge(newD,scal)',-1)
			call writemess('The input T is of symmetry type',-1)
			call error_stop
		end if
		TMPenlargeTensor=T
		rank=T%getRank()
		oldTotal=T%getTotalData()
		call pasteDimension(dimen,T%Dimension,1,rank-1,[newD])
		call dimen%setName(rank,T%getName(rank))
		call T%allocate(dimen,T%getType())
		call T%random(randomNumberScal)
		call ClassPointer1DFunc(T%Data%ClassData,1,oldTotal,clp)
		call FastCopyArray(clp,TMPenlargeTensor%Data%ClassData,oldTotal)
		return
	end subroutine


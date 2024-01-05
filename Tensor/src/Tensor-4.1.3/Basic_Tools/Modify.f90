	subroutine ModifyaValue(outA,inA,ith)
		class(*),intent(in)::inA
		class(*),intent(inout)::outA(:)
		integer,intent(in)::ith
		class(*),pointer::p
		if(ith.gt.size(outA))then
			call writemess('ERROR in ModifyaValue, ith>size',-1)
			call writemess('ith='+ith,-1)
			call writemess('size='+size(outA),-1)
			call error_stop
		end if
		call ClassPointer1DFunc(outA,ith,p)
		call CopyData(p,inA)
	end subroutine

	subroutine ModifySomeValue(outA,inA,ith,jth)
		integer,intent(in)::ith,jth
		class(*),intent(in)::inA(:)
		class(*),intent(inout)::outA(:)
		integer::i,ii,leninA
		class(*),pointer::p(:)
		if(ith.le.0)then
			call writemess('ERROR in ModifySomeValue, ith<=0',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		leninA=jth-ith+1
		if(ClassSize(inA).lt.leninA)then
			call writemess('ERROR in ModifySomeValue, (jth-ith+1)>size(input)',-1)
			call writemess('(jth-ith+1)='+(jth-ith+1),-1)
			call writemess('size(input)='+ClassSize(inA),-1)
			call error_stop
		end if
		if(jth.gt.ClassSize(outA))then
			call writemess('ERROR in ModifySomeValue, ith>size',-1)
			call writemess('jth='+jth,-1)
			call writemess('size='+ClassSize(outA),-1)
			call error_stop
		end if
		call ClassPointer1DFunc(outA,ith,jth,p)
		call FastCopyArray(p,inA,leninA)
		return
	end subroutine

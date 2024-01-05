	!pasteTensor(T1,T2,.true.)
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [l,m,n,...] matrix
		!			pasteTensor(T1,T2):a [2*l,m,n,...] matrix
		!           [1--->l, m, n,...] is T1
		!			[l+1--->2*l, m, n,...] is T2
		!        /   \
		!        | T1 |
		!        |----|
		!        | T2 |
		!        \    /
		!pasteTensor(T1,T2,.false.)
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n,l] matrix
		!			pasteTensor(T1,T2):a [...,m,n,2*l] matrix
		!        [... , m, n, 1--->l] is T1
		!			[... , m, n, l+1--->2*l] is T2
		!        /         \
		!        | T1 | T2 |
		!        \         /
		!        

	function pasteTensorFunc(T1,T2,row)result(pasteTensor)
		type(Tensor)::pasteTensor
		type(Tensor),intent(in)::T1
		type(Tensor),intent(in)::T2
		logical,intent(in)::row
		integer::i,classtype,rank1,rank2,pasteDim1,pasteDim2,collen,newDim_i,chalen
		type(Dimension)::newDim
		if(.not.T2%getFlag())return
		rank1=T1%Getrank()
		rank2=T2%Getrank()
		if(rank1.ne.rank2) then
			call writemess('can not paste two Tensor,ranks are,'+rank1+','+rank2,-1)
			call error_stop()
		end if
		if(row)then
			if(rank1.eq.1)then
				!call writemess('Do not finsiehd this part,can not paste two Tensor,ranks are,'+rank1+','+rank2)
				!call error_stop()
				pasteDim1=1
				pasteDim2=1
				newDim_i=2
				collen=T1%getTotalData()
				newDim=(/2,collen/)
			else
				pasteDim1=T1%dim(1)
				pasteDim2=T2%dim(1)
				newDim_i=pasteDim1+pasteDim2
				!newDim=(/newDim_i/)
				call pasteDimension(newDim,[newDim_i],T1%Dimension,2,rank1)
				collen=1
				do i=2,rank1
				!	newDim=newDim+(T1%Dimension.subDim.i)
					collen=collen*T1%dim(i)
					if(T1%dim(i).ne.T2%dim(i))then
						call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
						call writemess('T1%dim('+i+')='+T1%dim(i),-1)
						call writemess('T2%dim('+i+')='+T2%dim(i),-1)
						call error_stop()
					end if
				end do
			end if
			if(T1%getNameFlag()) then
				call newDim%setName(1,T1%getName(1))
			end if
			classtype=T1%getType()
			if(classtype.ne.T2%getType())then
				call writemess('ERROR type in  pasting two Tensor',-1)
				call error_stop
			end if

			chalen=T1%Data%DataCharacterLen
			if(T2%Data%DataCharacterLen.ne.chalen)then
				call writemess('ERROR type in  pasting two Tensor,chalen',-1)
				call error_stop
			end if
			call pasteTensor%allocate(newDim,classtype,chalen)
			call combinationRow_TData(pasteTensor,T1,T2,newDim_i,collen,pasteDim1,pasteDim2)
		else
			pasteDim1=T1%dim(rank1)
			pasteDim2=T2%dim(rank2)
			newDim_i=pasteDim1+pasteDim2
			collen=1
			if(rank1.eq.1)then
				newDim=(/newDim_i/)
			else
				call pasteDimension(newDim,T1%Dimension,1,rank1-1,[newDim_i])
			end if
			do i=1,rank1-1
			!	if(i.eq.1)then
			!		newDim=T1%Dimension.subDim.i
			!	else
			!		newDim=newDim+(T1%Dimension.subDim.i)
			!	end if
				collen=collen*T1%dim(i)
				if(T1%dim(i).ne.T2%dim(i))then
					call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
					call writemess('T1%dim('+i+')='+T1%dim(i),-1)
					call writemess('T2%dim('+i+')='+T2%dim(i),-1)
					call error_stop()
				end if
			end do
			!newDim=newDim+(/newDim_i/)
			if(T1%getNameFlag()) then
				call newDim%setName(rank1,T1%getName(rank1))
			end if
			classtype=T1%getType()
			if(classtype.ne.T2%getType())then
				call writemess('ERROR type in  pasting two Tensor',-1)
				call error_stop
			end if
			chalen=T1%Data%DataCharacterLen
			if(T2%Data%DataCharacterLen.ne.chalen)then
				call writemess('ERROR type in  pasting two Tensor,chalen',-1)
				call error_stop
			end if
			call pasteTensor%allocate(newDim,classtype,chalen)
			call combinationCol_TData(pasteTensor,T1,T2)
		end if
		return
	end function


	subroutine combinationRow_TData(R,A,B,LDR,LDR2,LDA,LDB)
		type(Tensor),intent(inout)::R
		type(Tensor),intent(in)::A,B
		integer,intent(in)::LDR,LDR2,LDA,LDB
		class(*),pointer::Rp(:,:),Ap(:,:),Bp(:,:)
		integer::Rlength,ALength,Blength
		Rlength=R%getTotalData()
		ALength=A%getTotalData()
		BLength=B%getTotalData()
		call ClassPointer2DFunc(Rlength,R%Data%ClassData,LDR,LDR2,[1,LDA],[1,LDR2],Rp)
		call ClassPointer2DFunc(ALength,A%Data%ClassData,LDA,LDR2,Ap)
		call CopyArray(Rp,Ap,LDA,LDR2)
		call ClassPointer2DFunc(Rlength,R%Data%ClassData,LDR,LDR2,[LDA+1,LDR],[1,LDR2],Rp)
		call ClassPointer2DFunc(BLength,B%Data%ClassData,LDB,LDR2,Bp)
		call CopyArray(Rp,Bp,LDB,LDR2)
		!RData(1:LDA,:)=AData(:,:LDR2)
		!RData(LDA+1:,:)=BData(:,:LDR2)
		return
	end subroutine
	subroutine combinationCol_TData(R,A,B)
		type(Tensor),intent(inout)::R
		type(Tensor),intent(in)::A,B
		class(*),pointer::Rp(:)
		integer::lengthA,lengthB,lengthR
		lengthA=A%getTotalData()
		lengthB=B%getTotalData()
		lengthR=lengthA+lengthB
		call ClassPointer1DFunc(R%Data%ClassData,1,lengthA,Rp)
		call FastCopyArray(Rp,A%Data%ClassData,lengthA)
		call ClassPointer1DFunc(R%Data%ClassData,lengthA+1,lengthR,Rp)
		call FastCopyArray(Rp,B%Data%ClassData,lengthB)
		return
	end subroutine


	!combinationCol:
		!			T1 :a [...,l,m,n] matrix
		!			T2 :a [...,l,m,n] matrix
		!			combination(T1,T2):a [...,l,m,n,2] matrix
		!			or 
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n] matrix
		!			combination(T1,T2):a [...,m,n,l+1] matrix

	function combinationColFunc(T1,T2)result(combination)
		type(Tensor),intent(in)::T1
		type(Tensor),intent(in)::T2
		type(Tensor)::combination
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n,classtype,chalen
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%getRank().eq.T2%getRank()) then
			do i=1,T1%getRank()
				if(T1%dim(i).ne.T2%dim(i))then
					call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
					call writemess('T1%dim('+i+')='+T1%dim(i),-1)
					call writemess('T2%dim('+i+')='+T2%dim(i),-1)
					call error_stop()
				end if
			end do
			call pasteDimension(newDim,T1%Dimension,1,T1%getRank(),[2])
			classtype=T1%getType()
			if(classtype.ne.T2%getType())then
				call writemess('ERROR type in  pasting two Tensor',-1)
				call error_stop
			end if
			chalen=T1%Data%DataCharacterLen
			if(T2%Data%DataCharacterLen.ne.chalen)then
				call writemess('ERROR type in  pasting two Tensor,chalen',-1)
				call error_stop
			end if
			call combination%allocate(newDim,classtype,chalen)
			call combinationCol_TData(combination,T1,T2)
			return
		end if

		do i=1,T1%getRank()-1
			if(T1%dim(i).ne.T2%dim(i))then
				call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
				call writemess('T1%dim('+i+')='+T1%dim(i),-1)
				call writemess('T2%dim('+i+')='+T2%dim(i),-1)
				call error_stop()
			end if
		end do
		dim_n=T1%dim(T1%getRank())
		call pasteDimension(newDim,T1%Dimension,1,T1%getRank()-1,[dim_n+1])
		classtype=T1%getType()
		if(classtype.ne.T2%getType())then
			call writemess('ERROR type in  pasting two Tensor',-1)
			call error_stop
		end if
		chalen=T1%Data%DataCharacterLen
		if(T2%Data%DataCharacterLen.ne.chalen)then
			call writemess('ERROR type in  pasting two Tensor,chalen',-1)
			call error_stop
		end if
		call combination%allocate(newDim,classtype,chalen)
		call combinationCol_TData(combination,T1,T2)
		return
	end function

	function combinationrowFunc(T1,T2)result(combinationrow)
		type(Tensor),intent(in)::T1
		type(Tensor),intent(in)::T2
		type(Tensor)::combinationrow
		integer,allocatable::dim1(:),dim2(:)
		integer::total1,total2,i,dim_n,classtype,chalen
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%getRank().eq.T2%getRank()) then
			do i=1,T1%getRank()
				if(T1%dim(i).ne.T2%dim(i))then
					call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
					call writemess('T1%dim('+i+')='+T1%dim(i),-1)
					call writemess('T2%dim('+i+')='+T2%dim(i),-1)
					call error_stop()
				end if
			end do
			call pasteDimension(newDim,[2],T1%Dimension,1,T1%getRank())
			classtype=T1%getType()
			if(classtype.ne.T2%getType())then
				call writemess('ERROR type in  pasting two Tensor',-1)
				call error_stop
			end if
			chalen=T1%Data%DataCharacterLen
			if(T2%Data%DataCharacterLen.ne.chalen)then
				call writemess('ERROR type in  pasting two Tensor,chalen',-1)
				call error_stop
			end if
			call combinationrow%allocate(newDim,classtype,chalen)
			total1=T1%getTotalData()
			total2=T2%getTotalData()
			call combinationRow_TData(combinationrow,T1,T2,2,total1,1,1)
			return
		end if

		do i=2,T1%getRank()
			if(T1%dim(i).ne.T2%dim(i-1))then
				call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
				call writemess('T1%dim('+i+')='+T1%dim(i),-1)
				call writemess('T2%dim('+i+')='+T2%dim(i),-1)
				call error_stop()
			end if
		end do
		dim_n=T1%dim(1)
		call pasteDimension(newDim,[dim_n+1],T1%Dimension,2,T1%getRank())
		classtype=T1%getType()
		if(classtype.ne.T2%getType())then
			call writemess('ERROR type in  pasting two Tensor',-1)
			call error_stop
		end if
		chalen=T1%Data%DataCharacterLen
		if(T2%Data%DataCharacterLen.ne.chalen)then
			call writemess('ERROR type in  pasting two Tensor,chalen',-1)
			call error_stop
		end if
		call combinationrow%allocate(newDim,classtype,chalen)
		total1=T1%getTotalData()
			total2=T2%getTotalData()
		call combinationRow_TData(combinationrow,T1,T2,dim_n+1,total2,dim_n,1)
		return
	end function	
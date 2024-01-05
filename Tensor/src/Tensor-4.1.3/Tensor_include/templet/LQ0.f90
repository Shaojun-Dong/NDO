	subroutine LQmatrixFuncName(A,L,M,N,classtype,info)!output Q matrix in A
		DATATYPE,intent(inout)::A(:,:)
		DATATYPE,intent(inout)::L(:,:)
		integer,intent(inout)::info
		integer,intent(in)::M,N,classtype
		integer :: i,j,min_MN
		DATATYPE,pointer::WQ(:),v(:),tempv(:),tau(:),dwork(:)
		DATATYPE,pointer::WorkQ(:,:)
		real*4,pointer::QN1(:),QN2(:)
		integer,pointer::dim(:)
		integer::LWORK

		min_MN=min(M,N)
		

		LWORK=m+m+n+n
		call WorkingMemory%allocate(classType,N+N+N+LWORK)
		call WorkingMemory%get_memory(v,N)
		call WorkingMemory%get_memory(tempv,N)
		call WorkingMemory%get_memory(tau,N)
		call WorkingMemory%get_memory(dwork,LWORK)


		INFO=999

		select case(classType)
			case(2)
				call SGELQF( M, N, A, M, Tau, dwork, LWORK, INFO )
			case(3)
				call DGELQF( M, N, A, M, Tau, dwork, LWORK, INFO )
			case(4)
				call CGELQF( M, N, A, M, Tau, dwork, LWORK, INFO )
			case(5)
				call ZGELQF( M, N, A, M, Tau, dwork, LWORK, INFO )
		end select
		
		

		if(M.le.N)then
			do i=1,min_MN
				L(i:M,i)=A(i:M,i)
				if(i.gt.1)L(1:i-1,i)=0
			end do
		else
			!L=A
			select case(classType)
				case(2)
					call scopy(size(A),A,1,L,1)
				case(3)
					call dcopy(size(A),A,1,L,1)
				case(4)
					call ccopy(size(A),A,1,L,1)
				case(5)
					call zcopy(size(A),A,1,L,1)
			end select
			do i=2,min_MN
				L(1:i-1,i)=0
			end do
		end if

		select case(classType)
			case(2)
				call SORGLQ( min_MN, N, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(3)
				call DORGLQ( min_MN, N, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(4)
				call CUNGLQ( min_MN, N, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(5)
				call ZUNGLQ( min_MN, N, min_MN, A, M, tau, dwork, LWORK, INFO )
		end select
		call WorkingMemory%free()
		return
	end subroutine

	subroutine LQSubroutineFUNCNAME(A,L,Q)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::L,Q
		type(QuanNum)::NewQN(2)
		integer,pointer::dim(:)
		integer::i,j,testi,testj,indexl,indexq,ii,jj
		DATATYPE,pointer::Ap(:,:),Qp(:,:),Lp(:,:)
		real*4,pointer::QN1(:),QN2(:),LQLQN(:),LQQQN(:)
		integer::LWORK,m,n,min_MN,classType,INFO,TotalData
		type(Dimension)::LDim,QDim
		if(.not.A%getSymmetryFlag())then
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			classType=A%getType()
			TotalData=A%getTotalData()
			call pasteDimension(LDim,A%Dimension,1,1,[min_MN])
			call pasteDimension(QDim,[min_MN],A%Dimension,2,2)
			call L%allocate(LDim,classtype)
			call Q%Allocate(QDim,classtype)
			call A%pointer(Ap)
			call L%pointer(Lp)
			call Q%pointer(Qp)
			call LQmatrixFuncName(Ap,Lp,M,N,classtype,info)
			Qp=Ap(1:min_MN,:)

			if(info.ne.0) then
				call writemess('Error in LQ decomposition ,info='+info,-1)
				call writemess('output The data in ./_LQ_ERROR_LOG.err',-1)
				open(unit=9991,file='./_LQ_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call L%write(9991)
				call Q%write(9991)
				close(9991)
				call error_stop()
			end if
			return
		end if	


		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(L%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(Q%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call L%pointDim(dim)
		call L%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call Q%pointDim(dim)
		call Q%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		call L%pointQN(LQLQN,2)
		call Q%pointQN(LQQQN,1)

		classType=A%getType()


		testi=0
		testj=0
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([i,j]).and.ifNonZeroBlockDATATYPE(A,i,j))then
					if(i.lt.testi)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in QRSubroutine on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					if(j.lt.testj)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in QRSubroutine on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					testi=i
					testj=j
					call get_index(LQLQN,QN2(j),indexl)
					call get_index(LQQQN,QN1(i),indexq)
					if((indexq.eq.0).or.(indexl.eq.0))then
						call writemess('ERRRO in LQ',-1)
						call NewQN(1)%print
						call writemess('j='+j)
						call writemess('QN2(j)='+QN2(j))
						call A%diminfo(.true.)
						call L%diminfo(.true.)
						call Q%diminfo(.true.)
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call L%setBlockMomery([i,indexl])
					call L%pointer(Lp,[i,indexl])
					call Q%setBlockMomery([indexq,j])
					call Q%pointer(Qp,[indexq,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(M,N)
					call LQmatrixFuncName(Ap,Lp,M,N,classtype,info)
					!Qp=Ap(1:min_MN,:)
					do jj=1,n
						do ii=1,min_MN
							Qp(ii,jj)=Ap(ii,jj)
						end do
					end do
					if(info.ne.0) then
						call writemess('Error in LQ decomposition ,info='+info,-1)
						call writemess('output The data in ./_LQ_ERROR_LOG.err',-1)
						open(unit=9991,file='./_LQ_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
						call A%write(9991)
						call L%write(9991)
						call Q%write(9991)
						close(9991)
						call error_stop()
					end if
				end if
			end do
		end do
		return
	end subroutine
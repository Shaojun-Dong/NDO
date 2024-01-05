	subroutine QRmatrixFuncName(A,R,M,N,classtype,info)
		integer,intent(in)::M,N,classType
		integer,intent(inout)::INFO
		DATATYPE,intent(inout)::A(:,:),R(:,:)
		DATATYPE,pointer::v(:),Tau(:),tempv(:),dwork(:)
		integer :: i,j,min_MN
		real*4,pointer::QN1(:),QN2(:)
		integer,pointer::dim(:)
		integer::LWORK

		min_MN=min(M,N)
		

		LWORK=m+m+n+n
		call WorkingMemory%allocate(classType,M+M+M+LWORK)
		call WorkingMemory%get_memory(v,M)
		call WorkingMemory%get_memory(tempv,M)
		call WorkingMemory%get_memory(tau,M)
		call WorkingMemory%get_memory(dwork,LWORK)

		select case(classType)
			case(2)
				call SGEQRF( M, N, A , M, tau , dWORK, LWORK, INFO )
			case(3)
				call DGEQRF( M, N, A , M, tau , dWORK, LWORK, INFO )
			case(4)
				call CGEQRF( M, N, A , M, tau , dWORK, LWORK, INFO )
			case(5)
				call ZGEQRF( M, N, A , M, tau , dWORK, LWORK, INFO )
		end select

		if(M.lt.N)then
			
			select case(classType)
				case(2)
					call scopy(size(R),A,1,R,1)
				case(3)
					call dcopy(size(R),A,1,R,1)
				case(4)
					call ccopy(size(R),A,1,R,1)
				case(5)
					call zcopy(size(R),A,1,R,1)
			end select
			do i=1,min_MN-1
				R(i+1:min_MN,i)=0
			end do
		else
			do i=1,min_MN
				R(:i,i)=A(:i,i)
				if(i.lt.min_MN)R(i+1:min_MN,i)=0
			end do
		end if

		select case(classType)
			case(2)
				call SORGQR( M, min_MN, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(3)
				call DORGQR( M, min_MN, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(4)
				call CUNGQR( M, min_MN, min_MN, A, M, tau, dwork, LWORK, INFO )
			case(5)
				call ZUNGQR( M, min_MN, min_MN, A, M, tau, dwork, LWORK, INFO )
		end select
		call WorkingMemory%free()
		return
	end subroutine


	subroutine QRSubroutineName(A,Q,R)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::Q,R
		type(QuanNum)::NewQN(2)
		integer,pointer::dim(:)
		integer::i,j,testi,testj,indexq,indexr,ii,jj
		DATATYPE,pointer::Ap(:,:),Qp(:,:),Rp(:,:)
		real*4,pointer::QN1(:),QN2(:),QRQQN(:),QRRQN(:)
		integer::LWORK,m,n,min_MN,classType,INFO,TotalData
		type(Dimension)::QDim,RDim
		if(.not.A%getSymmetryFlag())then
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			classType=A%getType()
			TotalData=A%getTotalData()
			call pasteDimension(QDim,A%Dimension,1,1,[min_MN])
			call pasteDimension(RDim,[min_MN],A%Dimension,2,2)
			call Q%allocate(QDim,classtype)
			call R%Allocate(RDim,classtype)
			call A%pointer(Ap)
			call Q%pointer(Qp)
			call R%pointer(Rp)
			call QRmatrixFuncName(Ap,Rp,M,N,classtype,info)
			Qp=Ap(:,1:min_MN)

			if(info.ne.0) then
				call writemess('Error in QR decomposition ,info='+info,-1)
				call writemess('output The data in ./_QR_ERROR_LOG.err',-1)
				open(unit=9991,file='./_QR_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call Q%write(9991)
				call R%write(9991)
				close(9991)
				call error_stop()
			end if
			return
		end if	


		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(Q%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(R%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call Q%pointDim(dim)
		call Q%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call R%pointDim(dim)
		call R%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		call Q%pointQN(QRQQN,2)
		call R%pointQN(QRRQN,1)
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
					call get_index(QRQQN,QN2(j),indexq)
					call get_index(QRRQN,QN1(i),indexr)
					if((indexq.eq.0).or.(indexr.eq.0))then
						call writemess('ERRRO in QR',-1)
						call NewQN(1)%print
						call writemess('j='+j)
						call writemess('QN2(j)='+QN2(j))
						call A%diminfo(.true.)
						call Q%diminfo(.true.)
						call R%diminfo(.true.)
						call A%symmetryCheck
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call Q%setBlockMomery([i,indexq])
					call Q%pointer(Qp,[i,indexq])
					call R%setBlockMomery([indexr,j])
					call R%pointer(Rp,[indexr,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(M,N)
					call QRmatrixFuncName(Ap,Rp,M,N,classtype,info)
					!Qp=Ap(:,1:min_MN)
					do jj=1,min_MN
						do ii=1,m
							Qp(ii,jj)=Ap(ii,jj)
						end do
					end do
					
					if(info.ne.0) then
						call writemess('Error in QR decomposition ,info='+info,-1)
						call writemess('output The data in ./_QR_ERROR_LOG.err',-1)
						open(unit=9991,file='./_QR_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
						call A%write(9991)
						call Q%write(9991)
						call R%write(9991)
						close(9991)
						call error_stop()
					end if
				end if
			end do
		end do
		return
	end subroutine
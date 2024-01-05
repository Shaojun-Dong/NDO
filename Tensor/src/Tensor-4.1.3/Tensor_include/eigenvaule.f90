	function eigTensor(T,outvex)	result(res)
		type(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::T
		logical,optional,intent(in)::outvex
		if(present(outvex).and.outvex)then
			allocate(res(2))
			call eigvalue(T,res(1),res(2))
		else
			allocate(res(1))
			call eigvalue(T,res(1))
		end if
		return
	end function

	subroutine eigvalue(H,val,vec)
		class(Tensor),intent(in) ::H
		type(Tensor),intent(inout) ::val
		type(Tensor),optional,intent(inout)::vec
		integer::hdim(2)
		if(H%getSymmetryFlag())then
			call writemess('ERROR in eigvalue',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		hdim(1)=H%dim(1)
		hdim(2)=H%dim(2)
		if(H%getRank().ne.2) then
			write(*,*)"ERROR in eng"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in eng"
			call error_stop()
		end if
		select case(H%getType())
			case (5)
				call val%allocate((/hdim(1)/),5)
				if(present(vec))call vec%allocate(H%Dimension,5)
			case (4)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))call vec%allocate(H%Dimension,4)
			case (3)
				call val%allocate((/hdim(1)/),5)
				if(present(vec))call vec%allocate(H%Dimension,3)
			case (2)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))call vec%allocate(H%Dimension,2)
			case (1)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))call vec%allocate(H%Dimension,2)
		end select
		call eigenvalue_suboutine(H,hdim(1),val,vec)
		return
	end subroutine	

	subroutine eigenvalue_suboutine(T,N,eigvalue,eigvector)
		integer,intent(in)::N
		type(Tensor) ,intent(inout):: eigvalue
		type(Tensor) ,optional,intent(inout)::eigvector
		type(Tensor) ,intent(in):: T
		character*1::eivFlag
		complex(kind=8),pointer::zHdata(:)
		complex(kind=4),pointer::cHdata(:)
		real(kind=8),pointer::dHdata(:)
		real(kind=4),pointer::sHdata(:)
		
		complex(kind=8),pointer::zeigVal(:)
		complex(kind=4),pointer::ceigVal(:)
		real(kind=8),pointer::deigValR(:),deigValI(:)
		real(kind=4),pointer::seigValR(:),seigValI(:)
		
		complex(kind=8),pointer::zeigVec(:)
		complex(kind=4),pointer::ceigVec(:)
		real(kind=8),pointer::deigVec(:)
		real(kind=4),pointer::seigVec(:)
		
		complex(kind=8),pointer::zwork(:)
		complex(kind=4),pointer::cwork(:)
		real(kind=8),pointer::dwork(:)
		real(kind=4),pointer::swork(:)
		
		real(kind=8),pointer::dRWORK(:)
		real(kind=4),pointer::sRWORK(:)
		
		logical::bwork(N)
		integer::INFO,lwork,SDIM

		if(present(eigvector))then
			eivFlag='V'
		else
			eivFlag='N'
		end if
		select case (T%getType())
			case(5)
				lwork=2*N
				call WorkingMemory%check()
				call WorkingMemory%allocate(5,T%gettotalData()+N+lwork+T%gettotalData())
				call WorkingMemory%allocate(3,N)
				call WorkingMemory%get_memory(zHdata,T%gettotalData())
				!call WorkingMemory%get_memory(zeigVal,N)
				call WorkingMemory%get_memory(dRWORK,N)
				call WorkingMemory%get_memory(zwork,lwork)
				if(present(eigvector))then
					call eigvector%pointer(zeigVec)
				else
					call WorkingMemory%get_memory(zeigVec,T%gettotalData())
				end if
				zHdata=T
				call eigvalue%pointer(zeigVal)
				call ZGEES(eivFlag,'N',1,N,zHdata,N,SDIM,zeigVal,zeigVec,N,zwork,LWORK,dRWORK,BWORK,INFO)
			case(4)

				lwork=2*N
				call WorkingMemory%check()
				call WorkingMemory%allocate(4,T%gettotalData()+N+lwork+T%gettotalData())
				call WorkingMemory%allocate(2,N)
				call WorkingMemory%get_memory(cHdata,T%gettotalData())
				!call WorkingMemory%get_memory(ceigVal,N)
				call WorkingMemory%get_memory(sRWORK,N)
				call WorkingMemory%get_memory(cwork,lwork)
				if(present(eigvector))then
					call eigvector%pointer(ceigVec)
				else
					call WorkingMemory%get_memory(ceigVec,T%gettotalData())
				end if
				cHdata=T
				call eigvalue%pointer(ceigVal)
				call CGEES(eivFlag,'N',1,N,cHdata,N,SDIM,ceigVal,ceigVec,N,cwork,LWORK,sRWORK,BWORK,INFO)
			case(3)
				lwork=3*N
				call WorkingMemory%check()
				call WorkingMemory%allocate(3,T%gettotalData()+N+N+N+lwork+T%gettotalData())
				call WorkingMemory%get_memory(dHdata,T%gettotalData())
				call WorkingMemory%get_memory(deigValR,N)
				call WorkingMemory%get_memory(deigValI,N)
				call WorkingMemory%get_memory(dRWORK,N)
				call WorkingMemory%get_memory(dwork,lwork)
				if(present(eigvector))then
					call eigvector%pointer(deigVec)
				else
					call WorkingMemory%get_memory(deigVec,T%gettotalData())
				end if
				dHdata=T

				call DGEES(eivFlag,'N',1,N,dHdata,N,SDIM,deigValR,deigValI,deigVec,N,dwork,LWORK,BWORK,INFO)
				call eigvalue%pointer(zeigVal)
				call zcopy(N,dcmplx(deigValR,deigValI),1,zeigVal,1)
			case(2)
				lwork=3*N
				call WorkingMemory%check()
				call WorkingMemory%allocate(2,T%gettotalData()+N+N+N+lwork+T%gettotalData())
				call WorkingMemory%get_memory(sHdata,T%gettotalData())
				call WorkingMemory%get_memory(seigValR,N)
				call WorkingMemory%get_memory(seigValI,N)
				call WorkingMemory%get_memory(sRWORK,N)
				call WorkingMemory%get_memory(swork,lwork)
				if(present(eigvector))then
					call eigvector%pointer(deigVec)
				else
					call WorkingMemory%get_memory(deigVec,T%gettotalData())
				end if
				dHdata=T


				call SGEES(eivFlag,'N',1,N,sHdata,N,SDIM,seigValR,seigValI,seigVec,N,swork,LWORK,BWORK,INFO)
				call eigvalue%pointer(ceigVal)
				call ccopy(N,cmplx(seigValR,seigValI,kind=4),1,ceigVal,1)
			case(1)
				lwork=3*N
				call WorkingMemory%check()
				call WorkingMemory%allocate(2,T%gettotalData()+N+N+N+lwork+T%gettotalData())
				call WorkingMemory%get_memory(sHdata,T%gettotalData())
				call WorkingMemory%get_memory(seigValR,N)
				call WorkingMemory%get_memory(seigValI,N)
				call WorkingMemory%get_memory(sRWORK,N)
				call WorkingMemory%get_memory(swork,lwork)
				if(present(eigvector))then
					call eigvector%pointer(deigVec)
				else
					call WorkingMemory%get_memory(deigVec,T%gettotalData())
				end if
				dHdata=T


				call SGEES(eivFlag,'N',1,N,sHdata,N,SDIM,seigValR,seigValI,seigVec,N,swork,LWORK,BWORK,INFO)
				call eigvalue%pointer(ceigVal)
				call ccopy(N,cmplx(seigValR,seigValI,kind=4),1,ceigVal,1)
			case default
				call writemess('ERROR in eig, error type='+T%getType(),-1)
				call error_stop
		end select

		if(info.ne.0) then
			call writemess("Error in eig ,info="+info,-1)
			call error_stop()
		end if
		call WorkingMemory%free()
		return
	end subroutine
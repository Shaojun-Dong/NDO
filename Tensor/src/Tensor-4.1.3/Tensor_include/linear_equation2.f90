
	!A*X=B

	function  LLS(A,B,RCOND)
		type(Tensor)::LLS
		class(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		class(*),intent(in),optional::RCOND
		integer::Na,Nb,Na2,classtype
		type(dimension)::Xdim
		if(A%getSymmetryFlag())then
			call writemess('ERROR in LLS',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		if(B%getSymmetryFlag())then
			call writemess('ERROR in LLS',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		Na=A%dim(1)
		if(A%getRank().ne.2)then
			call writemess("error in LLS,A should be a matrix",-1)
			call error_stop()
		end if
		Na2=A%dim(2)
		if(B%getRank().eq.1)then
			Nb=1
		else
			Nb=B%dim(2)
		end if
		if(Na.ne.(B%dim(1))) then
			call writemess("error in LLS,dimension of A and B",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu_class_type(A%getType(),B%getType())
		if(classtype.eq.1)classtype=2

		if(B%getRank().ne.1)then
			call pasteDimension(XDim,A%Dimension,1,1,B%Dimension,2,2)
		else
			call getsubDimension(A%Dimension,[1,1],Xdim)
		end if
		call LLS%allocate(Xdim,classtype)
		call LLS2_subroute(LLS,A,B,Na,Na2,Nb,RCOND)
		return
	end function	

	subroutine  SolveLLS(X,A,B,RCOND)
		class(Tensor),intent(inout)::X
		Type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		class(*),intent(in),optional::RCOND
		integer::Na,Nb,Na2,classtype
		type(dimension)::Xdim
		if(A%getSymmetryFlag())then
			call writemess('ERROR in LLS',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		if(B%getSymmetryFlag())then
			call writemess('ERROR in LLS',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		Na=A%dim(1)
		if(A%getRank().ne.2)then
			call writemess("error in LLS,A should be a matrix",-1)
			call error_stop()
		end if
		Na2=A%dim(2)
		if(B%getRank().eq.1)then
			Nb=1
		else
			Nb=B%dim(2)
		end if
		if(Na.ne.(B%dim(1))) then
			call writemess("error in LLS,dimension of A and B",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu_class_type(A%getType(),B%getType())
		if(classtype.eq.1)classtype=2

		if(B%getRank().ne.1)then
			call pasteDimension(XDim,A%Dimension,1,1,B%Dimension,2,2)
		else
			call getsubDimension(A%Dimension,[1,1],Xdim)
		end if
		call X%allocate(Xdim,classtype)
		call LLS2_subroute(X,A,B,Na,Na2,Nb,RCOND)
		return
	end subroutine

	subroutine LLS2_subroute(R,A,B,M,N,NRHS,RCOND_class)
		type(Tensor),intent(inout)::R
		type(Tensor),intent(in)::A,B
		integer,intent(in)::M,N,NRHS
		class(*),optional,intent(in)::RCOND_class
		real*8::dRCOND
		real*4::sRCOND
		real(kind=4),pointer::sdata(:),sWORK(:),Rsp(:),sRWORK(:)
		real(kind=8),pointer::ddata(:),dWORK(:),Rdp(:),dRWORK(:)
		complex(kind=4),pointer::cdata(:),cWORK(:),Rcp(:)
		complex(kind=8),pointer::zdata(:),zWORK(:),Rzp(:)
		integer,pointer::JPVT(:)
		integer::INFO,RANK,LWORK,MN
		if(present(RCOND_class))then
			select type(RCOND_class)
				type is (real(kind=4))
				 dRCOND=RCOND_class
				 sRCOND=RCOND_class
				type is (real(kind=8))
				 dRCOND=RCOND_class
				 sRCOND=RCOND_class
				type is (integer)
				 dRCOND=RCOND_class
				 sRCOND=RCOND_class
			end select
		else
			dRCOND=1d-16
			sRCOND=1e-8
		end if
		MN=min(M,N)
		LWORK=max(MN+3*N+1,2*MN+NRHS)
		select case(R%getType())
			case(2)
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,N)
				call WorkingMemory%allocate(2,LWORK+A%gettotalData())
				call WorkingMemory%get_memory(sWORK,LWORK)
				call WorkingMemory%get_memory(sdata,A%gettotalData())
				call WorkingMemory%get_memory(JPVT,N)
				sdata=A
				call R%pointer(Rsp)
				Rsp=B
				if(dRCOND.le.0)sRCOND=1e-8
				call SGELSY( M, N, NRHS, sdata, M, Rsp, N, JPVT, sRCOND, RANK, sWORK, LWORK, INFO)
				call WorkingMemory%free()
			case(3)
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,N)
				call WorkingMemory%allocate(3,LWORK+A%gettotalData())
				call WorkingMemory%get_memory(dWORK,LWORK)
				call WorkingMemory%get_memory(ddata,A%gettotalData())
				call WorkingMemory%get_memory(JPVT,N)
				ddata=A
				call R%pointer(Rdp)
				Rdp=B
				if(dRCOND.le.0)dRCOND=1d-16
				call DGELSY( M, N, NRHS, ddata, M, Rdp, N, JPVT, dRCOND, RANK, dWORK, LWORK, INFO )
				call WorkingMemory%free()
			case(4)
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,N)
				call WorkingMemory%allocate(2,N+N)
				call WorkingMemory%allocate(4,LWORK+A%gettotalData())
				call WorkingMemory%get_memory(cWORK,LWORK)
				call WorkingMemory%get_memory(sRWORK,N+N)
				call WorkingMemory%get_memory(cdata,A%gettotalData())
				call WorkingMemory%get_memory(JPVT,N)
				cdata=A
				call R%pointer(Rcp)
				Rcp=B
				if(sRCOND.le.0)sRCOND=1e-8
				call CGELSY( M, N, NRHS, cdata, M, Rcp, N, JPVT, sRCOND, RANK,  cWORK, LWORK, sRWORK, INFO)
				call WorkingMemory%free()
			case(5)
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,N)
				call WorkingMemory%allocate(3,N+N)
				call WorkingMemory%allocate(5,LWORK+A%gettotalData())
				call WorkingMemory%get_memory(zWORK,LWORK)
				call WorkingMemory%get_memory(dRWORK,N+N)
				call WorkingMemory%get_memory(zdata,A%gettotalData())
				call WorkingMemory%get_memory(JPVT,N)
				zdata=A
				call R%pointer(Rzp)
				Rzp=B
				if(dRCOND.le.0)dRCOND=1d-16
				call ZGELSY( M, N, NRHS, zdata, M, Rzp, N, JPVT, dRCOND, RANK,  zWORK, LWORK, dRWORK, INFO)
				call WorkingMemory%free()
		end select
		if(INFO.ne.0)then
			call writemess("?GELSS is not successful",-1)
			call writemess("INFO="+INFO,-1)
			call error_stop()
		end if
		return
	end subroutine


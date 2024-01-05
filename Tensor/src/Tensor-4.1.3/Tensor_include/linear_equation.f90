
	!A*X=B
		!input A  and B
		!output X
		!if A^{-1} may not exit,input RCOND,
		!perform SVD on A, Only keep  singular values S(i) <= RCOND*S(1) 
		!if RCONDM<0 keep the S(i)>0

	function  linequ(A,B,RCOND)
		type(Tensor)::linequ
		class(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		class(*),intent(in),optional::RCOND
		integer::Na,Nb,Na2,classtype
		type(dimension)::Xdim
		if(A%getSymmetryFlag())then
			call writemess('ERROR in linequ',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		if(B%getSymmetryFlag())then
			call writemess('ERROR in linequ',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		Na=A%dim(1)
		if(A%getRank().ne.2)then
			call writemess("error in linequ,A should be a matrix",-1)
			call error_stop()
		end if
		Na2=A%dim(2)
		if(B%getRank().eq.1)then
			Nb=1
		else
			Nb=B%dim(2)
		end if
		if(Na.ne.(B%dim(1))) then
			call writemess("error in linequ,dimension of A and B",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu_class_type(A%getType(),B%getType())
		if(classtype.eq.1)classtype=2

		if(B%getRank().ne.1)then
			call pasteDimension(XDim,A%Dimension,1,1,B%Dimension,2,2)
		else
			call getsubDimension(A%Dimension,[1,1],Xdim)
		end if
		call linequ%allocate(Xdim,classtype)
		if(present(RCOND))then
			call linequ2_subroute(linequ,A,B,Na,Na2,Nb,RCOND)
		else
			if(Na.ne.Na2) then
				write(*,*)"error in linequ,dimension of A"
				call error_stop()
			end if
			call linequ_routine_TData(linequ,A,B,Na,Nb)
		end if
		return
	end function	

	subroutine  Solvelinequ(X,A,B,RCOND)
		class(Tensor),intent(inout)::X
		Type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		class(*),intent(in),optional::RCOND
		integer::Na,Nb,Na2,classtype
		type(dimension)::Xdim
		if(A%getSymmetryFlag())then
			call writemess('ERROR in linequ',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		if(B%getSymmetryFlag())then
			call writemess('ERROR in linequ',-1)
			call writemess('DO NOT finshed the symmetry case',-1)
			call error_stop
		end if
		Na=A%dim(1)
		if(A%getRank().ne.2)then
			call writemess("error in linequ,A should be a matrix",-1)
			call error_stop()
		end if
		Na2=A%dim(2)
		if(B%getRank().eq.1)then
			Nb=1
		else
			Nb=B%dim(2)
		end if
		if(Na.ne.(B%dim(1))) then
			call writemess("error in linequ,dimension of A and B",-1)
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
		if(present(RCOND))then
			call linequ2_subroute(X,A,B,Na,Na2,Nb,RCOND)
		else
			if(Na.ne.Na2) then
				write(*,*)"error in linequ,dimension of A"
				call error_stop()
			end if
			call linequ_routine_TData(X,A,B,Na,Nb)
		end if
		return
	end subroutine

	subroutine linequ2_subroute(R,A,B,Na,Na2,Nb,RCOND_class)
		type(Tensor),intent(inout)::R
		type(Tensor),intent(in)::A,B
		integer,intent(in)::Na,Na2,Nb
		class(*),intent(in)::RCOND_class
		real*8::dRCOND
		real*4::sRCOND
		real(kind=4),pointer::sdata(:),sWORK(:),sS(:),Rsp(:)
		real(kind=8),pointer::ddata(:),dWORK(:),dS(:),Rdp(:)
		complex(kind=4),pointer::cdata(:),cWORK(:),cRWORK(:),Rcp(:)
		complex(kind=8),pointer::zdata(:),zWORK(:),zRWORK(:),Rzp(:)
		integer::INFO,RANK,LWORK,maxMn
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
		maxMn=max(Na,Na2)
		LWORK=3*maxMn+MAX(2*maxMn,Nb,maxMn)+1
		select case(R%getType())
			case(2)
				call WorkingMemory%check()
				call WorkingMemory%allocate(2,A%gettotalData()+min(Na,Na2)+LWORK)
				call WorkingMemory%get_memory(sdata,A%gettotalData())
				call WorkingMemory%get_memory(ss,min(Na,Na2))
				call WorkingMemory%get_memory(sWORK,LWORK)
				sdata=A
				call R%pointer(Rsp)
				Rsp=B
				call SGELSS(Na,Na2, Nb, sdata, Na, Rsp, Na, sS, sRCOND,RANK, sWORK, LWORK, INFO )
				call WorkingMemory%free()
			case(3)
				call WorkingMemory%check()
				call WorkingMemory%allocate(3,A%gettotalData()+min(Na,Na2)+LWORK)
				call WorkingMemory%get_memory(ddata,A%gettotalData())
				call WorkingMemory%get_memory(ds,min(Na,Na2))
				call WorkingMemory%get_memory(dWORK,LWORK)
				ddata=A
				call R%pointer(Rdp)
				Rdp=B
				call DGELSS(Na,Na2, Nb, ddata, Na, Rdp, Na, dS, dRCOND,RANK, dWORK, LWORK, INFO )
				call WorkingMemory%free()
			case(4)
				call WorkingMemory%check()
				call WorkingMemory%allocate(4,A%gettotalData()+LWORK+5*min(Na,Na2))
				call WorkingMemory%allocate(2,min(Na,Na2))
				call WorkingMemory%get_memory(cdata,A%gettotalData())
				call WorkingMemory%get_memory(ss,min(Na,Na2))
				call WorkingMemory%get_memory(cWORK,LWORK)
				call WorkingMemory%get_memory(cRWORK,5*min(Na,Na2))
				cdata=A
				call R%pointer(Rcp)
				Rcp=B
				call CGELSS(Na,Na2, Nb, cdata, Na, Rcp, Na, sS, sRCOND,RANK,cWORK, LWORK,cRWORK, INFO )
				call WorkingMemory%free()
			case(5)
				call WorkingMemory%check()
				call WorkingMemory%allocate(5,A%gettotalData()+LWORK+max(5*min(Na,Na2)-4,1))
				call WorkingMemory%allocate(3,min(Na,Na2))
				call WorkingMemory%get_memory(zdata,A%gettotalData())
				call WorkingMemory%get_memory(ds,min(Na,Na2))
				call WorkingMemory%get_memory(zWORK,LWORK)
				call WorkingMemory%get_memory(zRWORK,max(5*min(Na,Na2)-4,1))
				zdata=A
				call R%pointer(Rzp)
				Rzp=B
				call ZGELSS(Na,Na2, Nb, zdata, Na, Rzp, Na, dS, dRCOND,RANK,zWORK, LWORK,zRWORK, INFO )
				call WorkingMemory%free()
		end select
		if(INFO.ne.0)then
			call writemess("?GELSS is not successful",-1)
			call writemess("INFO="+INFO,-1)
			call error_stop()
		end if
		return
	end subroutine

	!Solves a general system of linear equations
		!A*X=B,find X
		!X is a vector or a matrix,the dimension of which is the same as B
		!on output,A will change and X is in B	
		!the subroutine ZGESV will change A	

	subroutine linequ_routine_TData(R,A,B,Na,Nb)
		type(Tensor),intent(inout)::R
		type(Tensor),intent(in)::A,B
		integer,intent(in)::Na,Nb
		integer,pointer::IPIV(:)
		real(kind=4),pointer::sdata(:),Rsp(:)
		real(kind=8),pointer::ddata(:),Rdp(:)
		complex(kind=4),pointer::cdata(:),Rcp(:)
		complex(kind=8),pointer::zdata(:),Rzp(:)
		integer::INFO
		
		select case(R%getType())
			case(2)
				call WorkingMemory%check()
				call WorkingMemory%allocate(2,A%gettotalData())
				call WorkingMemory%allocate(1,Na)
				call WorkingMemory%get_memory(sdata,A%gettotalData())
				call WorkingMemory%get_memory(IPIV,Na)
				sdata=A
				call R%pointer(Rsp)
				Rsp=B
				call SGESV(Na,Nb,sdata,Na,IPIV,Rsp,Na,INFO)
			case(3)
				call WorkingMemory%check()
				call WorkingMemory%allocate(3,A%gettotalData())
				call WorkingMemory%allocate(1,Na)
				call WorkingMemory%get_memory(ddata,A%gettotalData())
				call WorkingMemory%get_memory(IPIV,Na)
				ddata=A
				call R%pointer(Rdp)
				Rdp=B
				call DGESV(Na,Nb,ddata,Na,IPIV,Rdp,Na,INFO)
			case(4)
				call WorkingMemory%check()
				call WorkingMemory%allocate(4,A%gettotalData())
				call WorkingMemory%allocate(1,Na)
				call WorkingMemory%get_memory(cdata,A%gettotalData())
				call WorkingMemory%get_memory(IPIV,Na)
				cdata=A
				call R%pointer(Rcp)
				Rcp=B
				call CGESV(Na,Nb,cdata,Na,IPIV,Rcp,Na,INFO)	
			case(5)
				call WorkingMemory%check()
				call WorkingMemory%allocate(5,A%gettotalData())
				call WorkingMemory%allocate(1,Na)
				call WorkingMemory%get_memory(zdata,A%gettotalData())
				call WorkingMemory%get_memory(IPIV,Na)
				zdata=A
				call R%pointer(Rzp)
				Rzp=B
				call ZGESV(Na,Nb,zdata,Na,IPIV,Rzp,Na,INFO)	
		end select
		call WorkingMemory%free()
		if(INFO.ne.0)then
			call writemess("ZGESV is not successful",-1)
			call writemess("INFO="+INFO,-1)
			call error_stop()
		end if
		return
	end subroutine
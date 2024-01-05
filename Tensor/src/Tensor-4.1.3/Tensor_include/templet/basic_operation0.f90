	function minFuncName(A,w)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		character(len=*),optional,intent(in)::w
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: The input tensor has no data',-1)
			call error_stop
		end if
		select case(A%getType())
			case (1)
				call A%Data%pointAllData(ip)
				if(present(w))then
					if(w.equ.'abs')then
						Res=minval(abs(ip))
					else if(w.equ.'mina')then
						Res=minval(abs(ip))
					else if(w.equ.'real')then
						Res=minval(sp)
					else
						call writemess('ERROR in put in min(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=minval(ip)
				end if
				
			case (2)
				call A%Data%pointAllData(sp)
				if(present(w))then
					if(w.equ.'abs')then
						Res=minval(abs(sp))
					else if(w.equ.'mina')then
						Res=minval(abs(sp))
					else if(w.equ.'real')then
						Res=minval(sp)
					else
						call writemess('ERROR in put in min(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=minval(sp)
				end if
				
			case (3)
				call A%Data%pointAllData(dp)
				if(present(w))then
					if(w.equ.'abs')then
						Res=minval(abs(dp))
					else if(w.equ.'mina')then
						Res=minval(abs(dp))
					else if(w.equ.'real')then
						Res=minval(dp)
					else
						call writemess('ERROR in put in min(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=minval(dp)
				end if
			case (4)
				call A%Data%pointAllData(cp)
				if(present(w))then
					if(w.equ.'imag')then
						Res=minval(aimag(cp))
					else if(w.equ.'real')then
						Res=minval(real(cp))
					else if(w.equ.'abs')then
						Res=minval(abs(cp))
					else
						call writemess('ERROR in put in min(w)',-1)
						call writemess('the valid w are:imag,real,abs',-1)
						call error_stop
					end if
				else
					Res=minval(abs(cp))
				end if
			case (5)
				call A%Data%pointAllData(zp)
				if(present(w))then
					if(w.equ.'imag')then
						Res=minval(aimag(zp))
					else if(w.equ.'real')then
						Res=minval(real(zp))
					else if(w.equ.'abs')then
						Res=minval(abs(zp))
					else
						call writemess('ERROR in put in min(w)',-1)
						call writemess('the valid w are:imag,real,abs',-1)
						call error_stop
					end if
				else
					Res=minval(abs(zp))
				end if
			case default
				call writemess('ERROR in min, datatype='+A%getType(),-1)
				call errOr_stop
		end select
		return
	end function
	function maxFuncName(A,w)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		character(len=*),optional,intent(in)::w
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: The input tensor has no data',-1)
			call error_stop
		end if
		select case(A%getType())
			case (1)
				call A%Data%pointAllData(ip)
				if(present(w))then
					if(w.equ.'abs')then
						Res=maxval(abs(ip))
					else if(w.equ.'maxa')then
						Res=maxval(abs(ip))
					else if(w.equ.'real')then
						Res=maxval(ip)
					else
						call writemess('ERROR in put in max(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=maxval(ip)
				end if
				
			case (2)
				call A%Data%pointAllData(sp)
				if(present(w))then
					if(w.equ.'abs')then
						Res=maxval(abs(sp))
					else if(w.equ.'maxa')then
						Res=maxval(abs(sp))
					else if(w.equ.'real')then
						Res=maxval(sp)
					else
						call writemess('ERROR in put in max(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=maxval(sp)
				end if
				
			case (3)
				call A%Data%pointAllData(dp)
				if(present(w))then
					if(w.equ.'abs')then
						Res=maxval(abs(dp))
					else if(w.equ.'maxa')then
						Res=maxval(abs(dp))
					else if(w.equ.'real')then
						Res=maxval(dp)
					else
						call writemess('ERROR in put in max(w)',-1)
						call writemess('the valid w are:imag,real,abs for complex number',-1)
						call writemess('the valid w are:abs for real number',-1)
						call error_stop
					end if
				else
					Res=maxval(dp)
				end if
				
			case (4)
				call A%Data%pointAllData(cp)
				if(present(w))then
					if(w.equ.'imag')then
						Res=maxval(aimag(cp))
					else if(w.equ.'real')then
						Res=maxval(real(cp))
					else if(w.equ.'abs')then
						Res=maxval(abs(cp))
					else if(w.equ.'maxa')then
						Res=maxval(abs(cp))
					else
						call writemess('ERROR in put in max(w)',-1)
						call writemess('the valid w are:imag,real,abs',-1)
						call error_stop
					end if
				else
					Res=maxval(abs(cp))
				end if
			case (5)
				call A%Data%pointAllData(zp)
				if(present(w))then
					if(w.equ.'imag')then
						Res=maxval(aimag(zp))
					else if(w.equ.'real')then
						Res=maxval(real(zp))
					else if(w.equ.'abs')then
						Res=maxval(abs(zp))
					else if(w.equ.'maxa')then
						Res=maxval(abs(zp))
					else
						call writemess('ERROR in put in max(w)',-1)
						call writemess('the valid w are:imag,real,abs',-1)
						call error_stop
					end if
				else
					Res=maxval(abs(zp))
				end if
			case default
				call writemess('ERROR in max, datatype='+A%getType(),-1)
				call errOr_stop
		end select
		return
	end function

	function sumFuncName(A)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: The input tensor has no data',-1)
			call error_stop
		end if
		select case(A%getType())
			case (1)
				call A%Data%pointAllData(ip)
				Res=sum(ip)
			case (2)
				call A%Data%pointAllData(sp)
				Res=sum(sp)
			case (3)
				call A%Data%pointAllData(dp)
				Res=sum(dp)
			case (4)
				call A%Data%pointAllData(cp)
				Res=sum(cp)
			case (5)
				call A%Data%pointAllData(zp)
				Res=sum(zp)
			case default
				call writemess('ERROR in max, datatype='+A%getType(),-1)
				call errOr_stop
		end select
		return
	end function

	function traceFuncName(A)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		integer,pointer::ip(:,:)
		real*4,pointer::sp(:,:)
		real*8,pointer::dp(:,:)
		complex*8,pointer::cp(:,:)
		complex*16,pointer::zp(:,:)
		integer::i,blocki
		integer,pointer::dim(:)
		if(.not.A%getFlag())then
			call writemess('ERROR: The input tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getDimFlag())then
			call writemess('ERROR: The dimension of the input tensor is empty',-1)
			call error_stop
		end if
		if(A%getRank().ne.2)then
			call writemess('ERROR in trace tensor, input tensor should be a matrix',-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call A%pointDim(dim)
		select case(A%getType())
			case (1)
				Res=0
				if(A%getSymmetryFlag())then
					do blocki=1,dim(1)
						call A%pointer(ip,[blocki,blocki])
						if(associated(ip))then
							do i=1,size(ip,1)
								Res=Res+ip(i,i)
							end do
						end if
					end do
				else
					call A%pointer(ip)
					do i=1,size(ip,1)
						Res=Res+ip(i,i)
					end do
				end if
			case (2)
				Res=0
				if(A%getSymmetryFlag())then
					do blocki=1,dim(1)
						call A%pointer(sp,[blocki,blocki])
						if(associated(sp))then
							do i=1,size(sp,1)
								Res=Res+sp(i,i)
							end do
						end if
					end do
				else
					call A%pointer(sp)
					do i=1,size(sp,1)
						Res=Res+sp(i,i)
					end do
				end if
			case (3)
				Res=0
				if(A%getSymmetryFlag())then
					do blocki=1,dim(1)
						call A%pointer(dp,[blocki,blocki])
						if(associated(dp))then
							do i=1,size(dp,1)
								Res=Res+dp(i,i)
							end do
						end if
					end do
				else
					call A%pointer(dp)
					do i=1,size(dp,1)
						Res=Res+dp(i,i)
					end do
				end if
			case (4)
				Res=0
				if(A%getSymmetryFlag())then
					do blocki=1,dim(1)
						call A%pointer(cp,[blocki,blocki])
						if(associated(cp))then
							do i=1,size(cp,1)
								Res=Res+cp(i,i)
							end do
						end if
					end do
				else
					call A%pointer(cp)
					do i=1,size(cp,1)
						Res=Res+cp(i,i)
					end do
				end if
			case (5)
				Res=0
				if(A%getSymmetryFlag())then
					do blocki=1,dim(1)
						call A%pointer(zp,[blocki,blocki])
						if(associated(zp))then
							do i=1,size(zp,1)
								Res=Res+zp(i,i)
							end do
						end if
					end do
				else
					call A%pointer(zp)
					do i=1,size(zp,1)
						Res=Res+zp(i,i)
					end do
				end if
			case default
				call writemess('ERROR in trace, datatype='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function normFuncName(A)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		real(kind=4),external::snrm2,scnrm2
		real(kind=8),external::dznrm2,dnrm2
		integer::i
		if(.not.A%Data%getFlag())then
			if(A%getSymmetryFlag())then
				Res=0
				return
			end if
			call writemess('ERROR: The input tensor has no data',-1)
			call error_stop
		end if
		select case(A%getType())
			case (1)
				call A%Data%pointAllData(ip)
				Res=sqrt(real(DOT_PRODUCT(ip,ip)))
			case (2)
				call A%Data%pointAllData(sp)
				Res=snrm2(size(sp),sp,1)
			case (3)
				call A%Data%pointAllData(dp)
				Res=dnrm2(size(dp),dp,1)
			case (4)
				call A%Data%pointAllData(cp)
				Res=scnrm2(size(cp),cp,1)
			case (5)
				call A%Data%pointAllData(zp)
				Res=dznrm2(size(zp),zp,1)
			case default
				call writemess('ERROR in norm2, datatype='+A%getType(),-1)
				call errOr_stop
		end select
		return
	end function

	function norm2FuncName(A)result(Res)
		DATATYPE::Res
		class(Tensor),intent(in)::A
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		real(kind=4),external::snrm2,scnrm2
		real(kind=8),external::dznrm2,dnrm2
		integer::i
		if(.not.A%Data%getFlag())then
			if(A%getSymmetryFlag())then
				Res=0
				return
			end if
			call writemess('ERROR: The input tensor has no data',-1)
			call error_stop
		end if
		select case(A%getType())
			case (1)
				call A%Data%pointAllData(ip)
				Res=sqrt(real(DOT_PRODUCT(ip,ip)))
			case (2)
				call A%Data%pointAllData(sp)
				Res=snrm2(size(sp),sp,1)
			case (3)
				call A%Data%pointAllData(dp)
				Res=dnrm2(size(dp),dp,1)
			case (4)
				call A%Data%pointAllData(cp)
				Res=scnrm2(size(cp),cp,1)
			case (5)
				call A%Data%pointAllData(zp)
				Res=dznrm2(size(zp),zp,1)
			case default
				call writemess('ERROR in norm2, datatype='+A%getType(),-1)
				call errOr_stop
		end select
		Res=Res*Res
		return
	end function

	function product_dotFuncName(A,B)Result(R)
		DATATYPE::R
		type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		integer::i,length
		real*4,External::	sdot
		real*8,External::	ddot
		complex(kind=4),External::	cdotu
		complex(kind=8),External::	zdotu
		integer,pointer::Aip(:),Bip(:)
		real*4,pointer::Asp(:),Bsp(:)
		real*8,pointer::Adp(:),Bdp(:)
		complex*8,pointer::Acp(:),Bcp(:)
		complex*16,pointer::Azp(:),Bzp(:)
		if(A%getType().ne.B%getType())then
			call writemess('ERROR in dot product',-1)
			call writemess('A%getType()='+A%getType(),-1)
			call writemess('B%getType()='+B%getType(),-1)
			call error_stop
		end if
		select case(A%GetType())
			case(1)
				R=0
				do i=1,A%getTotalBlock()
					call A%pointer(Aip,i)
					call B%pointer(Bip,i)
					if(associated(Aip).and.associated(Bip))then
						length=size(Aip)
						if(length.ne.size(Bip))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+DOT_PRODUCT(Aip,Bip)
					end if
				end do
				
			case(2)
				R=0.
				do i=1,A%getTotalBlock()
					call A%pointer(Asp,i)
					call B%pointer(Bsp,i)
					if(associated(Asp).and.associated(Bsp))then
						length=size(Asp)
						if(length.ne.size(Bsp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+sdot(length, Asp, 1, Bsp, 1)
					end if
				end do
			case(3)
				R=0d0
				do i=1,A%getTotalBlock()
					call A%pointer(Adp,i)
					call B%pointer(Bdp,i)
					if(associated(Adp).and.associated(Bdp))then
						length=size(Adp)
						if(length.ne.size(Bdp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+ddot(length, Adp, 1, Bdp, 1)
					end if
				end do
			case(4)
				R=cmplx(0.,kind=4)
				do i=1,A%getTotalBlock()
					call A%pointer(Acp,i)
					call B%pointer(Bcp,i)
					if(associated(Acp).and.associated(Bcp))then
						length=size(Acp)
						if(length.ne.size(Bcp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+cdotu(length, Acp, 1, Bcp, 1)
					end if
				end do
			case(5)
				R=dcmplx(0d0)
				do i=1,A%getTotalBlock()
					call A%pointer(Azp,i)
					call B%pointer(Bzp,i)
					if(associated(Azp).and.associated(Bzp))then
						length=size(Azp)
						if(length.ne.size(Bzp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+zdotu(length, Azp, 1, Bzp, 1)
					end if
				end do
		end select
		return
	end function

	function product_cdotFuncName(A,B)Result(R)
		DATATYPE::R
		type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		integer::i,length
		real*4,External::	sdot
		real*8,External::	ddot
		complex(kind=4),External::	cdotc
		complex(kind=8),External::	zdotc
		integer,pointer::Aip(:),Bip(:)
		real*4,pointer::Asp(:),Bsp(:)
		real*8,pointer::Adp(:),Bdp(:)
		complex*8,pointer::Acp(:),Bcp(:)
		complex*16,pointer::Azp(:),Bzp(:)
		if(A%getType().ne.B%getType())then
			call writemess('ERROR in dot product',-1)
			call writemess('A%getType()='+A%getType(),-1)
			call writemess('B%getType()='+B%getType(),-1)
			call error_stop
		end if
		select case(A%GetType())
			case(1)
				R=0
				do i=1,A%getTotalBlock()
					call A%pointer(Aip,i)
					call B%pointer(Bip,i)
					if(associated(Aip).and.associated(Bip))then
						length=size(Aip)
						if(length.ne.size(Bip))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+DOT_PRODUCT(Aip,Bip)
					end if
				end do
				
			case(2)
				R=0.
				do i=1,A%getTotalBlock()
					call A%pointer(Asp,i)
					call B%pointer(Bsp,i)
					if(associated(Asp).and.associated(Bsp))then
						length=size(Asp)
						if(length.ne.size(Bsp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+sdot(length, Asp, 1, Bsp, 1)
					end if
				end do
			case(3)
				R=0d0
				do i=1,A%getTotalBlock()
					call A%pointer(Adp,i)
					call B%pointer(Bdp,i)
					if(associated(Adp).and.associated(Bdp))then
						length=size(Adp)
						if(length.ne.size(Bdp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+ddot(length, Adp, 1, Bdp, 1)
					end if
				end do
			case(4)
				R=cmplx(0.,kind=4)
				do i=1,A%getTotalBlock()
					call A%pointer(Acp,i)
					call B%pointer(Bcp,i)
					if(associated(Acp).and.associated(Bcp))then
						length=size(Acp)
						if(length.ne.size(Bcp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+cdotc(length, Acp, 1, Bcp, 1)
					end if
				end do
			case(5)
				R=dcmplx(0d0)
				do i=1,A%getTotalBlock()
					call A%pointer(Azp,i)
					call B%pointer(Bzp,i)
					if(associated(Azp).and.associated(Bzp))then
						length=size(Azp)
						if(length.ne.size(Bzp))then
							call writemess('ERROR, the block has different length',-1)
							call error_stop
						end if
						R=R+zdotc(length, Azp, 1, Bzp, 1)
					end if
				end do
		end select
		return
	end function



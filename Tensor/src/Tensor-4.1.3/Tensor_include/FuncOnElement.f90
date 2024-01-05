	function elementSquare(A)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Rip(:)
		real*4,pointer::sp(:),Rsp(:)
		real*8,pointer::dp(:),Rdp(:)
		complex*8,pointer::cp(:),Rcp(:)
		complex*16,pointer::zp(:),Rzp(:)
		integer::i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%allocate(A)
		select case(A%getType())
			case(1)
				call Res%Data%pointAllData(Rip)
				call A%Data%pointAllData(ip)
				do i=1,A%getTotalData()
					Rip(i)=ip(i)*ip(i)
				end do
			case(2)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(sp)
				do i=1,A%getTotalData()
					Rsp(i)=sp(i)*sp(i)
				end do
			case(3)
				call Res%Data%pointAllData(Rdp)
				call A%Data%pointAllData(dp)
				do i=1,A%getTotalData()
					Rdp(i)=dp(i)*dp(i)
				end do
			case(4)
				call Res%Data%pointAllData(Rcp)
				call A%Data%pointAllData(cp)
				do i=1,A%getTotalData()
					Rcp(i)=cp(i)*cp(i)
				end do
			case(5)
				call Res%Data%pointAllData(Rzp)
				call A%Data%pointAllData(zp)
				do i=1,A%getTotalData()
					Rzp(i)=zp(i)*zp(i)
				end do
			case default
				call writemess('ERROR case in square',-1)
				call error_stop
		end select
		return
	end function
	function elementSquare2(A)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Rip(:)
		real*4,pointer::sp(:),Rsp(:)
		real*8,pointer::dp(:),Rdp(:)
		complex*8,pointer::cp(:),Rcp(:)
		complex*16,pointer::zp(:),Rzp(:)
		integer::i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%allocate(A)
		select case(A%getType())
			case(1)
				call Res%Data%pointAllData(Rip)
				call A%Data%pointAllData(ip)
				do i=1,A%getTotalData()
					Rip(i)=ip(i)*ip(i)
				end do
			case(2)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(sp)
				do i=1,A%getTotalData()
					Rsp(i)=sp(i)*sp(i)
				end do
			case(3)
				call Res%Data%pointAllData(Rdp)
				call A%Data%pointAllData(dp)
				do i=1,A%getTotalData()
					Rdp(i)=dp(i)*dp(i)
				end do
			case(4)
				call Res%Data%pointAllData(Rcp)
				call A%Data%pointAllData(cp)
				do i=1,A%getTotalData()
					Rcp(i)=cmplx(real(cp(i))*real(cp(i)),aimag(cp(i))*aimag(cp(i)))
				end do
			case(5)
				call Res%Data%pointAllData(Rzp)
				call A%Data%pointAllData(zp)
				do i=1,A%getTotalData()
					Rzp(i)=cmplx(real(zp(i))*real(zp(i)),aimag(zp(i))*aimag(zp(i)))
				end do
			case default
				call writemess('ERROR case in square',-1)
				call error_stop
		end select
		return
	end function
	function elementnorm2(A)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Rip(:)
		real*4,pointer::sp(:),Rsp(:)
		real*8,pointer::dp(:),Rdp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		integer::i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		select case(A%getType())
			case(1)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rip)
				call A%Data%pointAllData(ip)
				do i=1,A%getTotalData()
					Rip(i)=abs(ip(i))*abs(ip(i))
				end do
			case(2)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(sp)
				do i=1,A%getTotalData()
					Rsp(i)=abs(sp(i))*abs(sp(i))
				end do
			case(3)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rdp)
				call A%Data%pointAllData(dp)
				do i=1,A%getTotalData()
					Rdp(i)=dabs(dp(i))*dabs(dp(i))
				end do
			case(4)
				call Res%allocate(A,2)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(cp)
				do i=1,A%getTotalData()
					Rsp(i)=cabs(cp(i))*cabs(cp(i))
				end do
			case(5)
				call Res%allocate(A,3)
				call Res%Data%pointAllData(Rdp)
				call A%Data%pointAllData(zp)
				do i=1,A%getTotalData()
					Rdp(i)=cdabs(zp(i))*cdabs(zp(i))
				end do
			case default
				call writemess('ERROR case in square',-1)
				call error_stop
		end select
		return
	end function
	function elementSqrt(A)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:)
		real*4,pointer::sp(:),Rsp(:)
		real*8,pointer::dp(:),Rdp(:)
		complex*8,pointer::cp(:),Rcp(:)
		complex*16,pointer::zp(:),Rzp(:)
		integer::i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		select case(A%getType())
			case(1)
				call Res%allocate(A,2)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(ip)
				do i=1,A%getTotalData()
					Rsp(i)=sqrt(real(ip(i)))
				end do
			case(2)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rsp)
				call A%Data%pointAllData(sp)
				do i=1,A%getTotalData()
					Rsp(i)=sqrt(sp(i))
				end do
			case(3)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rdp)
				call A%Data%pointAllData(dp)
				do i=1,A%getTotalData()
					Rdp(i)=dsqrt(dp(i))
				end do
			case(4)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rcp)
				call A%Data%pointAllData(cp)
				do i=1,A%getTotalData()
					Rcp(i)=sqrt(cp(i))
				end do
			case(5)
				call Res%allocate(A)
				call Res%Data%pointAllData(Rzp)
				call A%Data%pointAllData(zp)
				do i=1,A%getTotalData()
					Rzp(i)=cdsqrt(zp(i))
				end do
			case default
				call writemess('ERROR case in square',-1)
				call error_stop
		end select
		return
	end function
	function elementProduct(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A,B
		integer::i,NewClassType
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		if(A%getTotalData().ne.B%getTotalData())then
			call writemess('ERROR in elementProduct,totalData',-1)
			call writemess('A%getTotalData()='+A%getTotalData(),-1)
			call writemess('B%getTotalData()='+B%getTotalData(),-1)
			call error_stop
		end if
		NewClassType=select_type_in_add_minu_class_type(A%getType(),B%getType())
		call Res%allocate(A,NewClassType)
		call ArrayTimeArray(Res%Data%ClassData,A%Data%ClassData,B%Data%ClassData,A%getTotalData())
		return
	end function
	function elementdivide(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A,B
		integer::i,NewClassType
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		if(A%getTotalData().ne.B%getTotalData())then
			call writemess('ERROR in elementProduct,totalData',-1)
			call error_stop
		end if
		NewClassType=select_type_in_add_minu_class_type(A%getType(),B%getType())
		call Res%allocate(A,NewClassType)
		call ArrayDivideArray(Res%Data%ClassData,A%Data%ClassData,B%Data%ClassData,A%getTotalData())
		return
	end function

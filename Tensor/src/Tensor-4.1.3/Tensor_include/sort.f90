

	subroutine sortTensor1(T,indices,increase,realpart)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::indices
		logical,intent(in)::increase,realpart
		integer,pointer::ip(:),indexp(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		call indices%empty()
		call indices%allocate([T%getTotalData()],'integer')
		call indices%pointer(indexp)
		select case(T%getType())
			case(1)
				call T%pointer(ip)
				call sortData(ip,indexp,increase)
			case(2)
				call T%pointer(sp)
				call sortData(sp,indexp,increase)
			case(3)
				call T%pointer(dp)
				call sortData(dp,indexp,increase)
			case(4)
				call T%pointer(cp)
				call sortData(cp,indexp,realpart,increase)
			case(5)
				call T%pointer(zp)
				call sortData(zp,indexp,realpart,increase)
			case default
				call writemess('ERROR type in sort Tensor',-1)
				call error_stop
		end 	select
		return
	end subroutine

	subroutine sortTensor2(T,increase,realpart)
		class(Tensor),intent(inout)::T
		logical,intent(in)::increase,realpart
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		select case(T%getType())
			case(1)
				call T%pointer(ip)
				call sortData(ip,increase)
			case(2)
				call T%pointer(sp)
				call sortData(sp,increase)
			case(3)
				call T%pointer(dp)
				call sortData(dp,increase)
			case(4)
				call T%pointer(cp)
				call sortData(cp,realpart,increase)
			case(5)
				call T%pointer(zp)
				call sortData(zp,realpart,increase)
			case default
				call writemess('ERROR type in sort Tensor',-1)
				call error_stop
		end 	select
		return
	end subroutine

	subroutine sortTensor3(T,indices,increase)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::indices
		logical,intent(in)::increase
		logical::realpart
		realpart=.true.
		call sortTensor1(T,indices,increase,realpart)
		return
	end subroutine

	subroutine sortTensor4(T,increase)
		class(Tensor),intent(inout)::T
		logical,intent(in)::increase
		logical::realpart
		realpart=.true.
		call sortTensor2(T,increase,realpart)
		return
	end subroutine

	subroutine sortTensor5(T)
		class(Tensor),intent(inout)::T
		logical::increase,realpart
		increase=.true.
		realpart=.true.
		call sortTensor2(T,increase,realpart)
		return
	end subroutine
	integer function iwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		integer,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (1)
				do whichindex=1,T%getTotalData()
					if(element.eq.T%ii(whichindex))return
				end do
			case default
				call writemess('The Tensor is integer, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	integer function swhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		real*4,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (2)
				do whichindex=1,T%getTotalData()
					if(abs(element-T%si(whichindex)).le. default_zero_real_number)return
				end do
			case default
				call writemess('The Tensor is real, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function dwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		real*8,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (3)
				do whichindex=1,T%getTotalData()
					if(dabs(element-T%di(whichindex)).le. default_zero_double_number)return
				end do
			case default
				call writemess('The Tensor is real*8, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function cwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		complex*8,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (4)
				do whichindex=1,T%getTotalData()
					if(cabs(element-T%ci(whichindex)).le. default_zero_real_number)return
				end do
			case default
				call writemess('The Tensor is complex(kind=4), input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function zwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		complex*16,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (5)
				do whichindex=1,T%getTotalData()
					if(cdabs(element-T%zi(whichindex)).le. default_zero_double_number)return
				end do
			case default
				call writemess('The Tensor is complex(kind=4), input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	integer function awhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (7)
				do whichindex=1,T%getTotalData()
					if(element.equ.T%ai(whichindex))return
				end do
			case default
				call writemess('The Tensor is character, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function awhichindex2(T,element,maxlen)result(whichindex)
		class(Tensor),intent(in)::T
		integer,intent(in)::maxlen
		character(len=*),intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (7)
				do whichindex=1,maxlen
					if(element.equ.T%ai(whichindex))return
				end do
			case default
				call writemess('The Tensor is character, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	logical function T_gt_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).gt.num
			case (2)
				T_eq_int=T%si(1).gt.num
			case (3)
				T_eq_int=T%di(1).gt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function DATANAME_gt_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.gt.T%ii(1)
			case (2)
				T_eq_int=num.gt.T%si(1)
			case (3)
				T_eq_int=num.gt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_ge_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ge.num
			case (2)
				T_eq_int=T%si(1).ge.num
			case (3)
				T_eq_int=T%di(1).ge.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function DATANAME_ge_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.ge.T%ii(1)
			case (2)
				T_eq_int=num.ge.T%si(1)
			case (3)
				T_eq_int=num.ge.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_lt_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).lt.num
			case (2)
				T_eq_int=T%si(1).lt.num
			case (3)
				T_eq_int=T%di(1).lt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function DATANAME_lt_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.lt.T%ii(1)
			case (2)
				T_eq_int=num.lt.T%si(1)
			case (3)
				T_eq_int=num.lt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_le_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).le.num
			case (2)
				T_eq_int=T%si(1).le.num
			case (3)
				T_eq_int=T%di(1).le.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function DATANAME_le_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.le.T%ii(1)
			case (2)
				T_eq_int=num.le.T%si(1)
			case (3)
				T_eq_int=num.le.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_eq_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function DATANAME_eq_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_equ_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).equ.real(num)
			case (3)
				T_eq_int=T%di(1).equ.dble(num)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function DATANAME_equ_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).equ.real(num)
			case (3)
				T_eq_int=T%di(1).equ.dble(num)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_ne_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.ne.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ne.num
			case (2)
				T_eq_int=T%si(1).ne.num
			case (3)
				T_eq_int=T%di(1).ne.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ne.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function DATANAME_ne_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.ne.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ne.num
			case (2)
				T_eq_int=T%si(1).ne.num
			case (3)
				T_eq_int=T%di(1).ne.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ne.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_nequ_DATANAME(T,num)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.nequ.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ne.num
			case (2)
				T_eq_int=T%si(1).nequ.real(num)
			case (3)
				T_eq_int=T%di(1).nequ.dble(num)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.nequ.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function DATANAME_nequ_T(num,T)result(T_eq_int)
		type(Tensor),intent(in) :: T
		DATATYPE,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ne.num
			case (2)
				T_eq_int=T%si(1).nequ.real(num)
			case (3)
				T_eq_int=T%di(1).nequ.dble(num)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.nequ.)!",-1)
				call error_stop()
		end select
		return
	end function

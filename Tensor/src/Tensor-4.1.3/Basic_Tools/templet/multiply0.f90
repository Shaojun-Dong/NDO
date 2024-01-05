	subroutine ClassTimenumFUNCName(Res,A,B,length)
		DATATYPE,target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in multiply, length<0',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				call ClassTimenumName(Res,A,B,length)
			type is (real(kind=4))
				call ClassTimenumName(Res,A,B,length)
			type is (real(kind=8))
				call ClassTimenumName(Res,A,B,length)
			type is (complex(kind=4))
				call ClassTimenumName(Res,A,B,length)
			type is (complex(kind=8))
				call ClassTimenumName(Res,A,B,length)
			class default
				call Class_Time_num(Res,A,B,length)
		end select
	end subroutine

	subroutine ClassTimenum2FUNCName(Res,B,A,length)
		DATATYPE,target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in multiply, length<0',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				call ClassTimenumName(Res,A,B,length)
			type is (real(kind=4))
				call ClassTimenumName(Res,A,B,length)
			type is (real(kind=8))
				call ClassTimenumName(Res,A,B,length)
			type is (complex(kind=4))
				call ClassTimenumName(Res,A,B,length)
			type is (complex(kind=8))
				call ClassTimenumName(Res,A,B,length)
			class default
				call num_Time_class(Res,B,A,length)
		end select
	end subroutine

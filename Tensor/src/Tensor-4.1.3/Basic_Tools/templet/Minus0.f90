	subroutine ClassMinusClassFUNCName(Res,A,B,length)
		DATATYPE,target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassMinusClass, length<0',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				call ClassMinusClassName(Res,A,B,length)
			type is (real(kind=4))
				call ClassMinusClassName(Res,A,B,length)
			type is (real(kind=8))
				call ClassMinusClassName(Res,A,B,length)
			type is (complex(kind=4))
				call ClassMinusClassName(Res,A,B,length)
			type is (complex(kind=8))
				call ClassMinusClassName(Res,A,B,length)
			class default
				call Class_Minus_Class(Res,A,B,length)
		end select
	end subroutine

	subroutine ClassPlusClassFUNCName(Res,A,B,length)
		DATATYPE,intent(inout)::Res(:)
		class(*),intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassPlusClass, length<0',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				call ClassPlusClassName(Res,A,B,length)
			type is (real(kind=4))
				call ClassPlusClassName(Res,A,B,length)
			type is (real(kind=8))
				call ClassPlusClassName(Res,A,B,length)
			type is (complex(kind=4))
				call ClassPlusClassName(Res,A,B,length)
			type is (complex(kind=8))
				call ClassPlusClassName(Res,A,B,length)
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine

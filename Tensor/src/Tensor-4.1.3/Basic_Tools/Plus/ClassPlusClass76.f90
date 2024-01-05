	subroutine ClassPlusClass76(Res,A,B,length)
		character(len=*),intent(inout)::Res(:)
		logical,intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(B).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		select type(B)
			type is (character(len=*))
				Res(1:length)=A(1:length)+B(1:length)
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine

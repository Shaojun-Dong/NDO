	subroutine ClassPlusClass12(Res,A,B,length)
		integer,intent(inout)::Res(:)
		real*4,intent(in)::A(:)
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
			type is (integer)
				Res(1:length)=A(1:length)+B(1:length)
			type is (real(kind=4))
				Res(1:length)=A(1:length)+B(1:length)
			type is (real(kind=8))
				Res(1:length)=A(1:length)+B(1:length)
			type is (complex(kind=4))
				Res(1:length)=A(1:length)+B(1:length)
			type is (complex(kind=8))
				Res(1:length)=A(1:length)+B(1:length)
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine

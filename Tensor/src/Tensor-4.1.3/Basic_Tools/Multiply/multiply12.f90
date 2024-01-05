	subroutine ClassTimenum12(Res,A,B,length)
		integer,target,intent(inout)::Res(:)
		real*4,target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		select type(B)
			type is (integer)
				Res(1:length)=A(1:length)*B
			type is (real(kind=4))
				Res(1:length)=A(1:length)*B
			type is (real(kind=8))
				Res(1:length)=A(1:length)*B
			type is (complex(kind=4))
				Res(1:length)=A(1:length)*B
			type is (complex(kind=8))
				Res(1:length)=A(1:length)*B
			class default
				call Class_Time_num(Res,A,B,length)
		end select
	end subroutine

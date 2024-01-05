	subroutine CopyArraya(outA,inA,length)
		integer,intent(in)::length
		character(len=*),intent(inout)::outA(:)
		class(*),intent(in)::inA(:)
		integer::i
		select type(inA)
			type is (integer)
				outA(1:length)=inA(1:length)
			type is (real(kind=4))
				outA(1:length)=inA(1:length)
			type is (real(kind=8))
				outA(1:length)=inA(1:length)
			type is (complex(kind=4))
				outA(1:length)=inA(1:length)
			type is (complex(kind=8))
				outA(1:length)=inA(1:length)
			type is (logical)
				outA(1:length)=inA(1:length)
			type is (character(len=*))
				outA(1:length)=inA(1:length)
			class default
				call copy_Class_Array(outA,inA,length)
		end select
	end subroutine

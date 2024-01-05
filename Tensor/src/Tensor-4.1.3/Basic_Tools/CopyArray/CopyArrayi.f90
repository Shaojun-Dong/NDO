	subroutine CopyArrayi(outA,inA,length)
		integer,intent(in)::length
		class(*),intent(in)::inA(:)
		integer,intent(inout)::outA(:)
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
				do i=1,length
					if(inA(i))then
						outA(i)=1
					else 
						outA(i)=0
					end if
				end do
			type is (character(len=*))
				outA(1:length)=inA(1:length)
			class default
				call copy_Class_Array(outA,inA,length)
		end select
	end subroutine

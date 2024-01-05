	subroutine CopyArrayl(outA,inA,length)
		integer,intent(in)::length
		logical,intent(inout)::outA(:)
		class(*),intent(in)::inA(:)
		integer::i
		select type(inA)
			type is (integer)
				do i=1,length
					if(inA(i).eq.1)then
						outA(i)=.true.
					else if(inA(i).eq.0)then
						outA(i)=.false.
					else
						call writemess('ERROR in logical=integer, the input should be 0 or 1',-1)
						call error_stop
					end if
				end do
			type is (real(kind=4))
				call writemess('ERROR in logical=real*4',-1)
				call error_stop
			type is (real(kind=8))
				call writemess('ERROR in logical=real*8',-1)
				call error_stop
			type is (complex(kind=4))
				call writemess('ERROR in logical=complex*8',-1)
				call error_stop
			type is (complex(kind=8))
				call writemess('ERROR in logical=complex*16',-1)
				call error_stop
			type is (logical)
				outA(1:length)=inA(1:length)
			type is (character(len=*))
				outA(1:length)=inA(1:length)
			class default
				call copy_Class_Array(outA,inA,length)
		end select
	end subroutine

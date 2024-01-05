	subroutine CopyArrayc(outA,inA,length)
		integer,intent(in)::length
		complex*8,intent(inout)::outA(:)
		class(*),intent(in)::inA(:)
		integer::i
		select type(inA)
			type is (integer)
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(inA(1:length),kind=4), 1, outA, 1)
				else
					outA(1:length)=inA(1:length)
				end if
			type is (real(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(inA(1:length),kind=4), 1, outA, 1)
				else
					outA(1:length)=inA(1:length)
				end if
			type is (real(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(inA(1:length),kind=4), 1, outA, 1)
				else
					outA(1:length)=inA(1:length)
				end if
			type is (complex(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, inA(1:length), 1, outA, 1)
				else
					outA(1:length)=inA(1:length)
				end if
			type is (complex(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(inA(1:length),kind=4), 1, outA, 1)
				else
					outA(1:length)=inA(1:length)
				end if
			type is (logical)
				call writemess('ERROR in logical=complex, the input should be 0 or 1',-1)
				call error_stop
			type is (character(len=*))
				outA(1:length)=inA(1:length)
			class default
				call copy_Class_Array(outA,inA,length)
		end select
	end subroutine

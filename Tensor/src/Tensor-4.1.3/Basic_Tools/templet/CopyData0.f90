
	subroutine copyDataFUNCNAME(outA,inA)
		DATATYPE,intent(inout)::outA
		class(*),intent(in)::inA
		select type(inA)
			type is (integer)
				outA=inA
			type is (real(kind=4))
				outA=inA
			type is (real(kind=8))
				outA=inA
			type is (complex(kind=4))
				outA=inA
			type is (complex(kind=8))
				outA=inA
			type is (logical)
				outA=inA
			type is (character(len=*))
				outA=inA
			class default
				call copy_Class_Data(outA,inA)
		end select
		return
	end subroutine

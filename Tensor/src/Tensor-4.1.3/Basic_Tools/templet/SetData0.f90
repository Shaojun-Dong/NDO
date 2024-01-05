
	subroutine SetDataFUNCNAME(outA,inA,length)
		DATATYPE,target,intent(inout)::outA(:)
		class(*),target,intent(in)::inA
		integer,intent(in)::length
		select type(inA)
			type is (integer)
				outA(1:length)=inA
			type is (real(kind=4))
				outA(1:length)=inA
			type is (real(kind=8))
				outA(1:length)=inA
			type is (complex(kind=4))
				outA(1:length)=inA
			type is (complex(kind=8))
				outA(1:length)=inA
			type is (logical)
				outA(1:length)=inA
			type is (character(len=*))
				outA(1:length)=inA
			class default
				call set_Class_Data(outA,inA,length)
		end select
		return
	end subroutine

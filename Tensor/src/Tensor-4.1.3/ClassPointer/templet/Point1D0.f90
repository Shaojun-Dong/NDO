	subroutine Point1D_All_ValueFUNCNAME(length,inA,outp)
		class(*),target::inA(:)
		DATATYPE,pointer::outp(:)
		integer,intent(in)::length
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				outp(1:length)=>inA(1:length)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine
	subroutine Point1D_Some_ValueFUNCNAME(inA,ith,jth,outp)
		integer::ith,jth
		class(*),target::inA(:)
		DATATYPE,pointer::outp(:)
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.jth)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				outp(1:(jth-ith+1))=>inA(ith:jth)
			class default
				call writemess('ERROR type in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point1D_a_ValueFUNCNAME(inA,ith,outp)
		integer::ith
		class(*),target::inA(:)
		DATATYPE,pointer::outp
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.ith)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				outp=>inA(ith)
			class default
				call writemess('ERROR type in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

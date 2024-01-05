	subroutine numMinusClassFUNCNAME(Res,A,B,length)
		DATATYPE,intent(inout)::Res(:)
		class(*),intent(in)::A,B(:)
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (-), length',-1)
			call error_stop
		end if
		if(size(B).lt.length)then
			call writemess('ERROR in (-), length',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
						Res(1:length)=A-B(1:length)
					type is (real(kind=4))
						Res(1:length)=A-B(1:length)
					type is (real(kind=8))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A-B(1:length)
					class default
						call num_Minus_Class(Res,A,B,length)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A-B(1:length)
					type is (real(kind=4))
						Res(1:length)=A-B(1:length)
					type is (real(kind=8))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A-B(1:length)
					class default
						call num_Minus_Class(Res,A,B,length)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A-B(1:length)
					type is (real(kind=4))
						Res(1:length)=A-B(1:length)
					type is (real(kind=8))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A-B(1:length)
					class default
						call num_Minus_Class(Res,A,B,length)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A-B(1:length)
					type is (real(kind=4))
						Res(1:length)=A-B(1:length)
					type is (real(kind=8))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A-B(1:length)
					class default
						call num_Minus_Class(Res,A,B,length)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A-B(1:length)
					type is (real(kind=4))
						Res(1:length)=A-B(1:length)
					type is (real(kind=8))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A-B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A-B(1:length)
					class default
						call num_Minus_Class(Res,A,B,length)
				end select

			class default
				call num_Minus_Class(Res,A,B,length)
		end select
	end subroutine


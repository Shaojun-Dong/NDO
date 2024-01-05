	subroutine ClassMinusnumFUNCNAME(Res,A,B,length)
		DATATYPE,intent(inout)::Res(:)
		class(*),intent(in)::A(:),B
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (-), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (-), length',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)-B
					type is (real(kind=4))
						Res(1:length)=A(1:length)-B
					type is (real(kind=8))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)-B
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)-B
					type is (real(kind=4))
						Res(1:length)=A(1:length)-B
					type is (real(kind=8))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)-B
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)-B
					type is (real(kind=4))
						Res(1:length)=A(1:length)-B
					type is (real(kind=8))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)-B
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)-B
					type is (real(kind=4))
						Res(1:length)=A(1:length)-B
					type is (real(kind=8))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)-B
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)-B
					type is (real(kind=4))
						Res(1:length)=A(1:length)-B
					type is (real(kind=8))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)-B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)-B
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			class default
				call Class_Plus_num(Res,A,B,length)
		end select
	end subroutine


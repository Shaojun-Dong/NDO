	subroutine NumTimenumFUNCNAME(Res,A,B)
		DATATYPE,target,intent(inout)::Res
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
							Res=A*B
					type is (real(kind=4))
							Res=A*B
					type is (real(kind=8))
							Res=A*B
					type is (complex(kind=4))
							Res=A*B
					type is (complex(kind=8))
							Res=A*B
					class default
						call Num_Time_num(Res,A,B)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
							Res=A*B
					type is (real(kind=4))
							Res=A*B
					type is (real(kind=8))
							Res=A*B
					type is (complex(kind=4))
							Res=A*B
					type is (complex(kind=8))
							Res=A*B
					class default
						call Num_Time_num(Res,A,B)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
							Res=A*B
					type is (real(kind=4))
							Res=A*B
					type is (real(kind=8))
							Res=A*B
					type is (complex(kind=4))
							Res=A*B
					type is (complex(kind=8))
							Res=A*B
					class default
						call Num_Time_num(Res,A,B)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
							Res=A*B
					type is (real(kind=4))
							Res=A*B
					type is (real(kind=8))
							Res=A*B
					type is (complex(kind=4))
							Res=A*B
					type is (complex(kind=8))
							Res=A*B
					class default
						call Num_Time_num(Res,A,B)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
							Res=A*B
					type is (real(kind=4))
							Res=A*B
					type is (real(kind=8))
							Res=A*B
					type is (complex(kind=4))
							Res=A*B
					type is (complex(kind=8))
							Res=A*B
					class default
						call Num_Time_num(Res,A,B)
				end select

			class default
				call Num_Time_num(Res,A,B)
		end select
	end subroutine


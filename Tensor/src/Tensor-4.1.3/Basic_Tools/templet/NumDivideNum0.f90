	subroutine NDivideNFUNCNAME(Res,A,B)
		DATATYPE,intent(inout)::Res
		class(*),intent(in)::A
		class(*),intent(in)::B
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
							Res=A/real(B,kind=4)
					type is (real(kind=4))
							Res=A/B
					type is (real(kind=8))
							Res=A/B
					type is (complex(kind=4))
							Res=A/B
					type is (complex(kind=8))
							Res=A/B
					class default
						call num_divide_num(Res,A,B)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
							Res=A/real(B,kind=4)
					type is (real(kind=4))
							Res=A/B
					type is (real(kind=8))
							Res=A/B
					type is (complex(kind=4))
							Res=A/B
					type is (complex(kind=8))
							Res=A/B
					class default
						call num_divide_num(Res,A,B)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
							Res=A/real(B,kind=4)
					type is (real(kind=4))
							Res=A/B
					type is (real(kind=8))
							Res=A/B
					type is (complex(kind=4))
							Res=A/B
					type is (complex(kind=8))
							Res=A/B
					class default
						call num_divide_num(Res,A,B)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
							Res=A/real(B,kind=4)
					type is (real(kind=4))
							Res=A/B
					type is (real(kind=8))
							Res=A/B
					type is (complex(kind=4))
							Res=A/B
					type is (complex(kind=8))
							Res=A/B
					class default
						call num_divide_num(Res,A,B)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
							Res=A/real(B,kind=4)
					type is (real(kind=4))
							Res=A/B
					type is (real(kind=8))
							Res=A/B
					type is (complex(kind=4))
							Res=A/B
					type is (complex(kind=8))
							Res=A/B
					class default
						call num_divide_num(Res,A,B)
				end select

			class default
				call num_divide_num(Res,A,B)
		end select
	end subroutine


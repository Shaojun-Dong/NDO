	subroutine ClassMinusClass41(Res,A,B,length)
		complex*8,intent(inout)::Res(:)
		integer,intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		complex*8,parameter::one=cmplx(-1.,kind=4)
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(B).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		select type(B)
			type is (integer)
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(A(1:length),kind=4), 1, Res, 1)
					call caxpy(length, one, cmplx(B(1:length),kind=4), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (real(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(A(1:length),kind=4), 1, Res, 1)
					call caxpy(length, one, cmplx(B(1:length),kind=4), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (real(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(A(1:length),kind=4), 1, Res, 1)
					call caxpy(length, one, cmplx(B(1:length),kind=4), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (complex(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(A(1:length),kind=4), 1, Res, 1)
					call caxpy(length, one, B, 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (complex(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call ccopy(length, cmplx(A(1:length),kind=4), 1, Res, 1)
					call caxpy(length, one, cmplx(B(1:length),kind=4), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			class default
				call Class_Minus_Class(Res,A,B,length)
		end select
	end subroutine

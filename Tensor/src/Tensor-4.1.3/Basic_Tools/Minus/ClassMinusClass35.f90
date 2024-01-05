	subroutine ClassMinusClass35(Res,A,B,length)
		real*8,intent(inout)::Res(:)
		complex*16,intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		integer::i
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
					call dcopy(length, dble(A(1:length)), 1, Res, 1)
					call daxpy(length, -1d0, dble(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (real(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call dcopy(length, dble(A(1:length)), 1, Res, 1)
					call daxpy(length, -1d0, dble(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (real(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call dcopy(length, dble(A(1:length)), 1, Res, 1)
					call daxpy(length, -1d0, B, 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (complex(kind=4))
				 if(length.gt.LAPACK_LENGTH)then
					call dcopy(length, dble(A(1:length)), 1, Res, 1)
					call daxpy(length, -1d0, dble(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			type is (complex(kind=8))
				 if(length.gt.LAPACK_LENGTH)then
					call dcopy(length, dble(A(1:length)), 1, Res, 1)
					call daxpy(length, -1d0, dble(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)-B(1:length)
				end if
			class default
				call Class_Minus_Class(Res,A,B,length)
		end select
	end subroutine

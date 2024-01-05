	subroutine ClassPlusClass53(Res,A,B,length)
		complex*16,intent(inout)::Res(:)
		real*8,intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		complex*16,parameter::one=dcmplx(1d0)
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
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zaxpy(length, one, dcmplx(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)+B(1:length)
				end if
			type is (real(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zaxpy(length, one, dcmplx(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)+B(1:length)
				end if
			type is (real(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zaxpy(length, one, dcmplx(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)+B(1:length)
				end if
			type is (complex(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zaxpy(length, one, dcmplx(B(1:length)), 1 , Res, 1)
				else
					Res(1:length)=A(1:length)+B(1:length)
				end if
			type is (complex(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zaxpy(length, one, B, 1 , Res, 1)
				else
					Res(1:length)=A(1:length)+B(1:length)
				end if
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine

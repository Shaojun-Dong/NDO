	subroutine ClassTimenum54(Res,A,B,length)
		complex*16,target,intent(inout)::Res(:)
		complex*8,target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		select type(B)
			type is (integer)
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zscal (length, dcmplx(B), Res, 1)
				else
					Res(1:length)=A(1:length)*B
				end if
			type is (real(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zscal (length, dcmplx(B), Res, 1)
				else
					Res(1:length)=A(1:length)*B
				end if
			type is (real(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zscal (length, dcmplx(B), Res, 1)
				else
					Res(1:length)=A(1:length)*B
				end if
			type is (complex(kind=4))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zscal (length, dcmplx(B), Res, 1)
				else
					Res(1:length)=A(1:length)*B
				end if
			type is (complex(kind=8))
				if(length.gt.LAPACK_LENGTH)then
					call zcopy(length, dcmplx(A(1:length)), 1, Res, 1)
					call zscal (length, B, Res, 1)
				else
					Res(1:length)=A(1:length)*B
				end if
			class default
				call Class_Time_num(Res,A,B,length)
		end select
	end subroutine

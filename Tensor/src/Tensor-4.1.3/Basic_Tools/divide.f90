	subroutine classdivideNum(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in divide, length<0',-1)
			call error_stop
		end if
		select type(B)
			type is (integer)
				call classTimeNum(Res,A,1./real(B,kind=4),length)
			type is (real(kind=4))
				call classTimeNum(Res,A,1./B,length)
			type is (real(kind=8))
				call classTimeNum(Res,A,1d0/B,length)
			type is (complex(kind=4))
				call classTimeNum(Res,A,cmplx(1.,kind=4)/B,length)
			type is (complex(kind=8))
				call classTimeNum(Res,A,dcmplx(1d0)/B,length)
			class default
				call Class_divide_num(Res,A,B,length)
		end select
	end subroutine

	subroutine default_Class_divide_num(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Class_divide_num is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_divide_num(Res(i),A(i),B)
		end do
	end subroutine

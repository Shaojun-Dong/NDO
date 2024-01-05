	type(Tensor) function expmMatrix1(H)result(expmTensor)
		class(Tensor),intent(in) ::H
		type(Tensor)::temp
		integer::i
		if(H%getRank().ne.2) then
			write(*,*)"ERROR in expm"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if(H%dim(1).ne.H%dim(2)) then
			write(*,*)"ERROR in expm"
			call error_stop()
		end if
		temp=H
		expmTensor=H
		do i=2,99999
			temp=temp*H
			if(temp%isZero())exit
			temp=temp/i
			expmTensor=expmTensor+temp
		end do
		temp=identityMatrix(H%dimension,H%getType())
		expmTensor=expmTensor+temp
		return
	end function

	function identityMatrix(dimen,classtype)result(Res)
		type(Tensor)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::classtype
		integer,pointer::ip(:,:)
		real*4,pointer::sp(:,:)
		real*8,pointer::dp(:,:)
		complex*8,pointer::cp(:,:)
		complex*16,pointer::zp(:,:)
		integer::i,ii
		if(dimen%getRank().ne.2)then
			call writemess('ERROR in identityMatrix',-1)
			call error_stop
		end if
		call Res%allocate(dimen,classtype)
		call Res%zero()
		if(Res%getSymmetryFlag())then
			select case(classtype)
				case(1)
					do i=1,Res%dim(1)
						call Res%pointer(ip,[i,i])
						if(associated(ip))then
							do ii=1,size(ip,1)
								ip(ii,ii)=1
							end do
						end if
					end do
				case(2)
					do i=1,Res%dim(1)
						call Res%pointer(sp,[i,i])
						if(associated(sp))then
							do ii=1,size(sp,1)
								sp(ii,ii)=1
							end do
						end if
					end do
				case(3)
					do i=1,Res%dim(1)
						call Res%pointer(dp,[i,i])
						if(associated(dp))then
							do ii=1,size(dp,1)
								dp(ii,ii)=1
							end do
						end if
					end do
				case(4)
					do i=1,Res%dim(1)
						call Res%pointer(cp,[i,i])
						if(associated(cp))then
							do ii=1,size(cp,1)
								cp(ii,ii)=1
							end do
						end if
					end do
				case(5)
					do i=1,Res%dim(1)
						call Res%pointer(zp,[i,i])
						if(associated(zp))then
							do ii=1,size(zp,1)
								zp(ii,ii)=1
							end do
						end if
					end do
				case default
					call writemess('ERROR in identityMatrix',-1)
			end select
		else
			select case(classtype)
				case(1)
					call Res%pointer(ip)
					do i=1,size(ip,1)
						ip(i,i)=1
					end do
				case(2)
					call Res%pointer(sp)
					do i=1,size(sp,1)
						sp(i,i)=1
					end do
				case(3)
					call Res%pointer(dp)
					do i=1,size(dp,1)
						dp(i,i)=1
					end do
				case(4)
					call Res%pointer(cp)
					do i=1,size(cp,1)
						cp(i,i)=1
					end do
				case(5)
					call Res%pointer(zp)
					do i=1,size(zp,1)
						zp(i,i)=1
					end do
				case default
					call writemess('ERROR in identityMatrix',-1)
			end select
		end if
		return
	end function

	function realTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,2)
		call Res%Data%pointAllData(sp)

		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				sp=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				sp=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				sp=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				sp=real(Acp)
			case(5)
				CALL A%Data%pointAllData(Azp)
				sp=dble(Azp)
			case default
				call writemess('ERROR in real(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function dbleTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,3)
		call Res%Data%pointAllData(dp)

		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				dp=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				dp=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				dp=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				dp=dble(Acp)
			case(5)
				CALL A%Data%pointAllData(Azp)
				dp=dble(Azp)
			case default
				call writemess('ERROR in dble(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function imagTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,2)
		call Res%Data%pointAllData(sp)

		select case(A%getType())
			case(1)
				sp=0
			case(2)
				sp=0
			case(3)
				sp=0
			case(4)
				CALL A%Data%pointAllData(Acp)
				sp=aimag(Acp)
			case(5)
				CALL A%Data%pointAllData(Azp)
				sp=aimag(Azp)
			case default
				call writemess('ERROR in imag(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function dimagTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,3)
		call Res%Data%pointAllData(dp)

		select case(A%getType())
			case(1)
				sp=0
			case(2)
				sp=0
			case(3)
				sp=0
			case(4)
				CALL A%Data%pointAllData(Acp)
				dp=aimag(Acp)
			case(5)
				CALL A%Data%pointAllData(Azp)
				dp=aimag(Azp)
			case default
				call writemess('ERROR in dimag(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function cmplxTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,4)
		call Res%Data%pointAllData(cp)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				cp=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				cp=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				cp=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				cp=Acp
			case(5)
				CALL A%Data%pointAllData(Azp)
				cp=Azp
			case default
				call writemess('ERROR in cmplx(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function

	function cmplxTensor2(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A,B
		integer,pointer::ip(:),Aip(:),Bip(:)
		real*4,pointer::sp(:),Asp(:),Bsp(:)
		real*8,pointer::dp(:),Adp(:),Bdp(:)
		complex*8,pointer::cp(:),Acp(:),Bcp(:)
		complex*16,pointer::zp(:),Azp(:),Bzp(:)
		if(A%getType().ne.B%getType())then
			call writemess('ERROR in cmplx(A,B), the data type of A and B should be the same',-1)
			call error_stop
		end if

		call Res%allocate(A,4)
		call Res%Data%pointAllData(cp)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				CALL B%Data%pointAllData(Bip)
				cp=cmplx(Aip,Bip)
			case(2)
				CALL A%Data%pointAllData(Asp)
				CALL B%Data%pointAllData(Bsp)
				cp=cmplx(Asp,Bsp)
			case(3)
				CALL A%Data%pointAllData(Adp)
				CALL B%Data%pointAllData(Bdp)
				cp=cmplx(Adp,Bdp)
			case default
				call writemess('ERROR in cmplx(A,B),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function
	function dcmplxTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,5)
		call Res%Data%pointAllData(zp)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				zp=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				zp=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				zp=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				zp=Acp
			case(5)
				CALL A%Data%pointAllData(Azp)
				zp=Azp
			case default
				call writemess('ERROR in dcmplx(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function
	function dcmplxTensor2(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A,B
		integer,pointer::ip(:),Aip(:),Bip(:)
		real*4,pointer::sp(:),Asp(:),Bsp(:)
		real*8,pointer::dp(:),Adp(:),Bdp(:)
		complex*8,pointer::cp(:),Acp(:),Bcp(:)
		complex*16,pointer::zp(:),Azp(:),Bzp(:)
		if(A%getType().ne.B%getType())then
			call writemess('ERROR in dcmplx(A,B), the data type of A and B should be the same',-1)
			call error_stop
		end if

		call Res%allocate(A,5)
		call Res%Data%pointAllData(zp)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				CALL B%Data%pointAllData(Bip)
				zp=dcmplx(Aip,Bip)
			case(2)
				CALL A%Data%pointAllData(Asp)
				CALL B%Data%pointAllData(Bsp)
				zp=dcmplx(Asp,Bsp)
			case(3)
				CALL A%Data%pointAllData(Adp)
				CALL B%Data%pointAllData(Bdp)
				zp=dcmplx(Adp,Bdp)
			case default
				call writemess('ERROR in cmplx(A,B),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function
	function intTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)

		call Res%allocate(A,1)
		call Res%Data%pointAllData(ip)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				ip=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				ip=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				ip=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				ip=Acp
			case(5)
				CALL A%Data%pointAllData(Azp)
				ip=Azp
			case default
				call writemess('ERROR in int(Tensor),classType='+A%getType(),-1)
				call error_stop
		end select
		return
	end function
	function charTensor(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,pointer::ip(:),Aip(:)
		real*4,pointer::sp(:),Asp(:)
		real*8,pointer::dp(:),Adp(:)
		complex*8,pointer::cp(:),Acp(:)
		complex*16,pointer::zp(:),Azp(:)
		logical,pointer::Alp(:)
		character(len=characterlen),pointer::ap(:),Aap(:)

		call Res%allocate(A,7)
		call Res%Data%pointAllData(Ap)
		select case(A%getType())
			case(1)
				CALL A%Data%pointAllData(Aip)
				ap=Aip
			case(2)
				CALL A%Data%pointAllData(Asp)
				ap=Asp
			case(3)
				CALL A%Data%pointAllData(Adp)
				ap=Adp
			case(4)
				CALL A%Data%pointAllData(Acp)
				ap=Acp
			case(5)
				CALL A%Data%pointAllData(Azp)
				ap=Azp
			case(6)
				CALL A%Data%pointAllData(Alp)
				ap=Alp
			case(7)
				CALL A%Data%pointAllData(Aap)
				ap=Aap
		end select
		return
	end function
	
	subroutine pauli_matrix(Sx,Sy,Sz,num)
		type(Tensor) :: Sx,Sy,Sz
		class(*),optional,intent(in)::num
		complex*16::II=(0d0,1d0)
		complex*8::sII=(0.,1.)
		complex*16::one=(1d0,0d0)
		complex*8::sone=(1.,0.)
		if(present(num))then
			select type(num)
				type is (real(kind=4))
					call Sx%allocate([2,2],2)
					call Sy%allocate([2,2],4)
					call Sz%allocate([2,2],2)
					call Sx%zero()
					call Sy%zero()
					call Sz%zero()
					call Sx%setValue([1,2],num)
					call Sx%setValue([2,1],num)
					call Sy%setValue([1,2],-num*sII)
					call Sy%setValue([2,1],num*sII)
					call Sz%setValue([1,1],num)
					call Sz%setValue([2,2],-num)
				type is (real(kind=8))
					call Sx%allocate([2,2],3)
					call Sy%allocate([2,2],5)
					call Sz%allocate([2,2],3)
					call Sx%zero()
					call Sy%zero()
					call Sz%zero()
					call Sx%setValue([1,2],num)
					call Sx%setValue([2,1],num)
					call Sy%setValue([1,2],-num*II)
					call Sy%setValue([2,1],num*II)
					call Sz%setValue([1,1],num)
					call Sz%setValue([2,2],-num)
				type is (complex(kind=4))
					call Sx%allocate([2,2],4)
					call Sy%allocate([2,2],4)
					call Sz%allocate([2,2],4)
					call Sx%zero()
					call Sy%zero()
					call Sz%zero()
					call Sx%setValue([1,2],num*sone)
					call Sx%setValue([2,1],num*sone)
					call Sy%setValue([1,2],-num*sII)
					call Sy%setValue([2,1],num*sII)
					call Sz%setValue([1,1],num*sone)
					call Sz%setValue([2,2],-num*sone)
				type is (complex(kind=8))
					call Sx%allocate([2,2],5)
					call Sy%allocate([2,2],5)
					call Sz%allocate([2,2],5)
					call Sx%zero()
					call Sy%zero()
					call Sz%zero()
					call Sx%setValue([1,2],num*one)
					call Sx%setValue([2,1],num*one)
					call Sy%setValue([1,2],-num*II)
					call Sy%setValue([2,1],num*II)
					call Sz%setValue([1,1],num*one)
					call Sz%setValue([2,2],-num*one)
				class default
					call writemess('ERROR in pauli_matrix',-1)
					call error_stop
			end select
		else
			call Sx%allocate([2,2],5)
			call Sy%allocate([2,2],5)
			call Sz%allocate([2,2],5)
			call Sx%zero()
			call Sy%zero()
			call Sz%zero()
			call Sx%setValue([1,2],one)
			call Sx%setValue([2,1],one)
			call Sy%setValue([1,2],-II)
			call Sy%setValue([2,1],II)
			call Sz%setValue([1,1],one)
			call Sz%setValue([2,2],-one)
		end if
		

		return
	end subroutine	
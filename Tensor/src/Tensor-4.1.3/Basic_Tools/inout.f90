	subroutine write_out_data(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		integer::i
		select type (A)
			type is (integer)
				write(uni,*)A(1:length)

			type is (real(kind=4))
				write(uni,*)A(1:length)

			type is (real(kind=8))
				write(uni,*)A(1:length)

			type is (complex(kind=4))
				write(uni,*)real(A(1:length),kind=4)
				write(uni,*)aimag(A(1:length))

			type is (complex(kind=8))
				write(uni,*)dble(A(1:length))
				write(uni,*)aimag(A(1:length))

			type is (logical)
				write(uni,*)A(1:length)

			type is (character(len=*))
				write(uni,*)(trim(adjustl(A(i)))//"  ",i=1,length)

			class default
				call write_class_data(A,length,uni)
		end select
	end subroutine

	subroutine read_externl_data(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		integer::i
		real*4,allocatable::rsp(:),isp(:)
		real*8,allocatable::rdp(:),idp(:)
		if(ClassSize(A).lt.length)then
			call writemess('ERROR in read_externl_data, ClassSize(A)='+ClassSize(A),-1)
			call writemess('length='+length,-1)
			call error_stop
		end if
		select type (A)
			type is (integer)
				read(uni,*)(A(i),i=1,length)

			type is (real(kind=4))
				read(uni,*)(A(i),i=1,length)

			type is (real(kind=8))
				read(uni,*)(A(i),i=1,length)

			type is (complex(kind=4))
				allocate(rsp(length))
				allocate(isp(length))
				read(uni,*)(rsp(i),i=1,length)
				read(uni,*)(isp(i),i=1,length)
				A=cmplx(rsp,isp,kind=4)

			type is (complex(kind=8))
				allocate(rdp(length))
				allocate(idp(length))
				read(uni,*)(rdp(i),i=1,length)
				read(uni,*)(idp(i),i=1,length)
				A=dcmplx(rdp,idp)

			type is (logical)
				read(uni,*)(A(i),i=1,length)

			type is (character(len=*))
				read(uni,*)A(1:length)

			class default
				call read_class_data(A,length,uni)
		end select
	end subroutine

	subroutine default_write_class_data(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		call writemess('ERROR write class data',-1)
		call writemess('DO NOT setting the default subroutine of write_class_data',-1)
		call error_stop
	end subroutine
	subroutine default_read_class_data(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		call writemess('ERROR read class data',-1)
		call writemess('DO NOT setting the default subroutine of read_class_data',-1)
		call error_stop
	end subroutine

#include "CopyArray/CopyArrayi.f90"
#include "CopyArray/CopyArrays.f90"
#include "CopyArray/CopyArrayd.f90"
#include "CopyArray/CopyArrayc.f90"
#include "CopyArray/CopyArrayz.f90"
#include "CopyArray/CopyArrayl.f90"
#include "CopyArray/CopyArraya.f90"

	subroutine FastcopyARRAY(outA,inA,length)
		integer,intent(in)::length
		class(*),target,intent(in)::inA(:)
		class(*),target,intent(inout)::outA(:)
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in FastcopyARRAY, length<0',-1)
			call error_stop
		end if
		if(ClassSize(inA).lt.length)then
			call writemess('ERROR in FastcopyARRAY, length>size',-1)
			call writemess('length='+length,-1)
			call writemess('size='+size(inA),-1)
			call error_stop
		end if
		if(ClassSize(outA).lt.length)then
			call writemess('ERROR in FastcopyARRAY, length>size',-1)
			call writemess('length='+length,-1)
			call writemess('size='+size(outA),-1)
			call error_stop
		end if
		select type(outA)
			type is (integer)
				call CopyArrayi(outA,inA,length)
			type is (real(kind=4))
				call CopyArrays(outA,inA,length)
			type is (real(kind=8))
				call CopyArrayd(outA,inA,length)
			type is (complex(kind=4))
				call CopyArrayc(outA,inA,length)
			type is (complex(kind=8))
				call CopyArrayz(outA,inA,length)
			type is (logical)
				call CopyArrayl(outA,inA,length)
			type is (character(len=*))
				call CopyArraya(outA,inA,length)
			class default
				call copy_Class_Array(outA,inA,length)
		end select
	end subroutine

	subroutine copyARRAYSameType2D(outA,inA,dim1,dim2)
		integer,intent(in)::dim1,dim2
		class(*),intent(in)::inA(:,:)
		class(*),intent(inout)::outA(:,:)
		logical::NoMatchFlag
		if(dim1.eq.0)return
		if(dim2.eq.0)return
		if(dim1.lt.0)then
			call writemess('ERROR in copyARRAYSameType, length<0',-1)
			call error_stop
		end if
		if(dim2.lt.0)then
			call writemess('ERROR in copyARRAYSameType, length<0',-1)
			call error_stop
		end if
		if(ClassSize(inA,1).lt.dim1)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim1='+dim1,-1)
			call writemess('ClassSize(inA,1)='+ClassSize(inA,1),-1)
			call error_stop
		end if
		if(ClassSize(inA,2).lt.dim2)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim2='+dim2,-1)
			call writemess('ClassSize(inA,2)='+ClassSize(inA,2),-1)
			call error_stop
		end if
		if(ClassSize(outA,1).lt.dim1)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim1='+dim1,-1)
			call writemess('ClassSize(outA,1)='+ClassSize(outA,1),-1)
			call error_stop
		end if
		if(ClassSize(outA,2).lt.dim2)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim2='+dim2,-1)
			call writemess('ClassSize(outA,2)='+ClassSize(outA,2),-1)
			call error_stop
		end if
		NoMatchFlag=.true.
		select type(inA)
			type is (integer)
				select type (outA)
					type is (integer)
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select

			type is (real(kind=4))
				select type (outA)
					type is (real(kind=4))
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select

			type is (real(kind=8))
				select type (outA)
					type is (real(kind=8))
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select

			type is (complex(kind=4))
				select type (outA)
					type is (complex(kind=4))
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select

			type is (complex(kind=8))
				select type (outA)
					type is (complex(kind=8))
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select

			type is (logical)
				select type (outA)
					type is (logical)
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select


			type is (character(len=*))
				select type (outA)
					type is (character(len=*))
						outA(1:dim1,1:dim2)=inA(1:dim1,1:Dim2)
						NoMatchFlag=.false.
				end select
		end select
		if(NoMatchFlag)then
			call copyARRAY_Same_Type_2D(outA,inA,dim1,dim2)
		end if
		return
	end subroutine


	subroutine copyARRAYSameType3D(outA,inA,dim1,dim2,dim3)
		integer,intent(in)::dim1,dim2,dim3
		class(*),intent(in)::inA(:,:,:)
		class(*),intent(inout)::outA(:,:,:)
		logical::NoMatchFlag
		if(dim1.eq.0)return
		if(dim2.eq.0)return
		if(dim3.eq.0)return
		if(dim1.lt.0)then
			call writemess('ERROR in copyARRAYSameType, length<0',-1)
			call error_stop
		end if
		if(dim2.lt.0)then
			call writemess('ERROR in copyARRAYSameType, length<0',-1)
			call error_stop
		end if
		if(dim3.lt.0)then
			call writemess('ERROR in copyARRAYSameType, length<0',-1)
			call error_stop
		end if
		if(ClassSize(inA,1).lt.dim1)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim1='+dim1,-1)
			call writemess('ClassSize(inA,1)='+ClassSize(inA,1),-1)
			call error_stop
		end if
		if(ClassSize(inA,2).lt.dim2)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim2='+dim2,-1)
			call writemess('ClassSize(inA,2)='+ClassSize(inA,2),-1)
			call error_stop
		end if
		if(ClassSize(inA,3).lt.dim3)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim3='+dim3,-1)
			call writemess('ClassSize(inA,3)='+ClassSize(inA,3),-1)
			call error_stop
		end if
		if(ClassSize(outA,1).lt.dim1)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim1='+dim1,-1)
			call writemess('ClassSize(outA,1)='+ClassSize(outA,1),-1)
			call error_stop
		end if
		if(ClassSize(outA,2).lt.dim2)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim2='+dim2,-1)
			call writemess('ClassSize(outA,2)='+ClassSize(outA,2),-1)
			call error_stop
		end if
		if(ClassSize(outA,3).lt.dim3)then
			call writemess('ERROR in copyARRAYSameType, length>size',-1)
			call writemess('dim3='+dim3,-1)
			call writemess('ClassSize(outA,3)='+ClassSize(outA,3),-1)
			call error_stop
		end if
		NoMatchFlag=.true.
		select type(inA)
			type is (integer)
				select type (outA)
					type is (integer)
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select

			type is (real(kind=4))
				select type (outA)
					type is (real(kind=4))
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select

			type is (real(kind=8))
				select type (outA)
					type is (real(kind=8))
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select

			type is (complex(kind=4))
				select type (outA)
					type is (complex(kind=4))
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select

			type is (complex(kind=8))
				select type (outA)
					type is (complex(kind=8))
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select

			type is (logical)
				select type (outA)
					type is (logical)
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select


			type is (character(len=*))
				select type (outA)
					type is (character(len=*))
						outA(1:dim1,1:dim2,1:dim3)=inA(1:dim1,1:Dim2,1:dim3)
						NoMatchFlag=.false.
				end select
		end select
		if(NoMatchFlag)then
			call copyARRAY_Same_Type_3D(outA,inA,dim1,dim2,dim3)
		end if
		return
	end subroutine

	subroutine default_copy_Class_Array(outA,inA,length)
		integer,intent(in)::length
		class(*),target,intent(in)::inA(:)
		class(*),target,intent(inout)::outA(:)
		class(*),pointer::inp,outp
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default copy_Class_Array is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call ClassPointer1DFunc(outA,i,outp)
			call ClassPointer1DFunc(inA,i,inp)
			call copy_Class_Data(outp,inp)
		end do
		return
	end subroutine
	subroutine default_copyARRAY_Same_Type_2D(outA,inA,dim1,dim2)
		integer,intent(in)::dim1,dim2
		class(*),target,intent(in)::inA(:,:)
		class(*),target,intent(inout)::outA(:,:)
		class(*),pointer::inp,outp
		integer::i,j
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default copy_Class_Data is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do j=1,dim2
			do i=1,dim1
				call ClassPointer2DFunc(outA,i,j,outp)
				call ClassPointer2DFunc(inA,i,j,inp)
				call copy_Class_Data(outp,inp)
			end do
		end do
		return
	end subroutine
	subroutine default_copyARRAY_Same_Type_3D(outA,inA,dim1,dim2,dim3)
		integer,intent(in)::dim1,dim2,dim3
		class(*),target,intent(in)::inA(:,:,:)
		class(*),target,intent(inout)::outA(:,:,:)
		class(*),pointer::inp,outp
		integer::i,j,k
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default copy_Class_Data is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do k=1,dim3
			do j=1,dim2
				do i=1,dim1
					call ClassPointer3DFunc(outA,i,j,k,outp)
					call ClassPointer3DFunc(inA,i,j,k,inp)
					call copy_Class_Data(outp,inp)
				end do
			end do
		end do
		return
	end subroutine

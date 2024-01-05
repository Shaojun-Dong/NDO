#define constructor0FuncName constructor0i
#define constructor1FuncName constructor1i
#define constructor2FuncName constructor2i
#define constructor3FuncName constructor3i
#define constructor4FuncName constructor4i
#define DATATYPE integer
#define DATATYPE2 integer
#define DATATYPENumber 1
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber

#define constructor0FuncName constructor0s
#define constructor1FuncName constructor1s
#define constructor2FuncName constructor2s
#define constructor3FuncName constructor3s
#define constructor4FuncName constructor4s
#define DATATYPE real*4
#define DATATYPE2 real*4
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber

#define constructor0FuncName constructor0d
#define constructor1FuncName constructor1d
#define constructor2FuncName constructor2d
#define constructor3FuncName constructor3d
#define constructor4FuncName constructor4d
#define DATATYPE real*8
#define DATATYPE2 real*8
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber

#define constructor0FuncName constructor0c
#define constructor1FuncName constructor1c
#define constructor2FuncName constructor2c
#define constructor3FuncName constructor3c
#define constructor4FuncName constructor4c
#define DATATYPE complex*8
#define DATATYPE2 complex*8
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber


#define constructor0FuncName constructor0z
#define constructor1FuncName constructor1z
#define constructor2FuncName constructor2z
#define constructor3FuncName constructor3z
#define constructor4FuncName constructor4z
#define DATATYPE complex*16
#define DATATYPE2 complex*16
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber


#define constructor0FuncName constructor0l
#define constructor1FuncName constructor1l
#define constructor2FuncName constructor2l
#define constructor3FuncName constructor3l
#define constructor4FuncName constructor4l
#define DATATYPE logical
#define DATATYPE2 logical
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber

#define constructor0FuncName constructor0a
#define constructor1FuncName constructor1a
#define constructor2FuncName constructor2a
#define constructor3FuncName constructor3a
#define constructor4FuncName constructor4a
#define DATATYPE character(len=*)
#define DATATYPE2 character(len=characterlen)
#define DATATYPENumber 2
#include "templet/constructor0.f90"
#undef constructor0FuncName
#undef constructor1FuncName
#undef constructor2FuncName
#undef constructor3FuncName
#undef constructor4FuncName
#undef DATATYPE
#undef DATATYPE2
#undef DATATYPENumber
	function constructor_char_scal(cha_,dimen)result(res)
		type(Tensor)::res
		character(len=*),intent(in)::cha_
		integer,optional,intent(in)::dimen(:)
		character(len=len(trim(adjustl(cha_))))::cha
		integer::TotalData,lenCha,i,j,ith,jth
		character(len=1)::w,divider,TensorType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex(kind=4),pointer::cp(:)
		complex(kind=8),pointer::zp(:)
		logical,pointer::lp(:)
		integer,allocatable::location(:)
		cha=trim(adjustl(cha_))
		lenCha=len(cha)
 		TotalData=1
 		divider=array_character_divider
 		TensorType=cha(1:1)
 		w=cha(2:2)
 		if(w.nequ.'=')then
 			Res=cha
			if(present(dimen))call Res%resetdim(dimen)
 			return
 		end if
 		do i=3,lenCha
 			w=cha(i:i)
 			if(w.equ.divider)TotalData=TotalData+1
 		end do
		allocate(location(TotalData+1))
		j=1
		location(1)=2
		do i=3,lenCha
			w=cha(i:i)
			if(w.equ.divider)then
				if(j.ge.(TotalData+1))then
					call writemess('ERROR in converting character for Tensor=character')
					call error_stop
				end if
				j=j+1
				location(j)=i
			end if
		end do
		if(j.ge.(TotalData+1))then
			call writemess('ERROR in converting character for Tensor=character')
			call error_stop
		end if
		j=j+1
		location(j)=lenCha+1
		select case(TensorType)
			case ('i')
				call Res%allocate([TotalData],'integer')
				call Res%pointer(ip)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)ip(i)
				end do
			case ('s')
				call Res%allocate([TotalData],'real*4')
				call Res%pointer(sp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)sp(i)
				end do
			case ('d')
				call Res%allocate([TotalData],'real*8')
				call Res%pointer(dp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)dp(i)
				end do
			case ('c')
				call writemess('DO NO finished this type, in Tensor=character')
				call error_stop
			case ('z')
				call writemess('DO NO finished this type, in Tensor=character')
				call error_stop
			case ('l')
				call Res%allocate([TotalData],'logical')
				call Res%pointer(lp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)lp(i)
				end do
			case ('a')
				call Res%allocate([TotalData],'character')
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					call Res%setValue(i,cha(ith:jth))
				end do
			case default
				call writemess('ERROR in converting character for Tensor=character')
				call error_stop
		end select
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
#define DataArrayithPointerFuncName DataArrayithPointeri
#define DataArrayAllPointerFuncName DataArrayAllPointeri
#define DATATYPE integer
#define DATATYPENumber 1
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber

#define DataArrayithPointerFuncName DataArrayithPointers
#define DataArrayAllPointerFuncName DataArrayAllPointers
#define DATATYPE real*4
#define DATATYPENumber 2
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber

#define DataArrayithPointerFuncName DataArrayithPointerd
#define DataArrayAllPointerFuncName DataArrayAllPointerd
#define DATATYPE real*8
#define DATATYPENumber 3
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber

#define DataArrayithPointerFuncName DataArrayithPointerc
#define DataArrayAllPointerFuncName DataArrayAllPointerc
#define DATATYPE complex*8
#define DATATYPENumber 4
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber

#define DataArrayithPointerFuncName DataArrayithPointerz
#define DataArrayAllPointerFuncName DataArrayAllPointerz
#define DATATYPE complex*16
#define DATATYPENumber 5
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber

#define DataArrayithPointerFuncName DataArrayithPointerl
#define DataArrayAllPointerFuncName DataArrayAllPointerl
#define DATATYPE logical
#define DATATYPENumber 6
#include "templet/pointer0.f90"
#undef DataArrayithPointerFuncName
#undef DataArrayAllPointerFuncName
#undef DATATYPE
#undef DATATYPENumber



	subroutine DataArrayithPointera(Da,p,ith)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::ith
		character(len=*),pointer,intent(inout)::p(:)
		integer::chalen
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		if(ith.gt.Da%getTotalBlock())then
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('ith='+ith,-1)
			call writemess('Da%getTotalBlock()='+Da%getTotalBlock(),-1)
			call error_stop
		end if
		if((Da%getType().eq.7).or.(Da%getType().eq.8))then
			if(Da%getFlag(ith))then
				chalen=Da%getCharacterlen()
				if(chalen.ne.len(p))then
					call writemess('ERROR in pointting to the character',-1)
					call error_stop
				end if
				call Pointer1DFunc(Da%ClassData,Da%starti(ith),Da%endi(ith),p,chalen)
			else
				call writemess('ERROR in pointting to the data of DataArray',-1)
				call writemess('the Block ith is empty',-1)
				call writemess('ith='+ith,-1)
				call error_stop
			end if
		else
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('Da%getType()='+Da%getType(),-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine DataArrayAllPointera(Da,p)
		class(DataArray),target,intent(in)::Da
		character(len=*),pointer,intent(inout)::p(:)
		integer::chalen
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		if((Da%getType().eq.7).or.(Da%getType().eq.8))then
			chalen=Da%getCharacterlen()
			if(chalen.ne.len(p))then
				call writemess('ERROR in pointting to the character',-1)
				call writemess('chalen='+chalen,-1)
				call writemess('len(p)='+len(p),-1)
				call error_stop
			end if
			call Pointer1DFunc(Da%totalData,Da%ClassData,p,chalen)
		else
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('Da%getType()='+Da%getType(),-1)
			call error_stop
		end if
		return
	end subroutine

	!************************************************


	subroutine pointStarti1D(Da,p)
		class(DataArray),target,intent(in)::Da
		integer,pointer,intent(inout)::p(:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer1DFunc(Da%getTotalBlock(),Da%Starti,p)
		return
	end subroutine
	subroutine pointStarti2D(Da,dim1,dim2,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2
		integer,pointer,intent(inout)::p(:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer2DFunc(Da%getTotalBlock(),Da%Starti,dim1,dim2,p)
		return
	end subroutine
	subroutine pointStarti3D(Da,dim1,dim2,dim3,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2,dim3
		integer,pointer,intent(inout)::p(:,:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer3DFunc(Da%getTotalBlock(),Da%Starti,dim1,dim2,dim3,p)
		return
	end subroutine
	subroutine pointStarti4D(Da,dim1,dim2,dim3,dim4,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2,dim3,dim4
		integer,pointer,intent(inout)::p(:,:,:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer4DFunc(Da%getTotalBlock(),Da%Starti,dim1,dim2,dim3,dim4,p)
		return
	end subroutine
	subroutine pointEndi1D(Da,p)
		class(DataArray),target,intent(in)::Da
		integer,pointer,intent(inout)::p(:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer1DFunc(Da%getTotalBlock(),Da%Endi,p)
		return
	end subroutine
	subroutine pointEndi2D(Da,dim1,dim2,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2
		integer,pointer,intent(inout)::p(:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer2DFunc(Da%getTotalBlock(),Da%Endi,dim1,dim2,p)
		return
	end subroutine
	subroutine pointEndi3D(Da,dim1,dim2,dim3,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2,dim3
		integer,pointer,intent(inout)::p(:,:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer3DFunc(Da%getTotalBlock(),Da%Endi,dim1,dim2,dim3,p)
		return
	end subroutine
	subroutine pointEndi4D(Da,dim1,dim2,dim3,dim4,p)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::dim1,dim2,dim3,dim4
		integer,pointer,intent(inout)::p(:,:,:,:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		call Pointer4DFunc(Da%getTotalBlock(),Da%Endi,dim1,dim2,dim3,dim4,p)
		return
	end subroutine
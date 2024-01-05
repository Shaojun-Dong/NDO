	subroutine Point3D_All_ValueFUNCNAME(length,inA,dim1,dim2,dim3,outp)
		integer::dim1,dim2,dim3,length
		class(*),target::inA(:)
		DATATYPE,pointer::outp(:,:,:)
		if((dim1*dim2*dim3).ne.length)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				outp(1:dim1,1:dim2,1:dim3)=>inA(1:length)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point3D_Some_ValueFUNCNAME(length,inA,dim1,dim2,dim3,ith,jth,kth,outp)
		integer::dim1,dim2,dim3,ith(2),jth(2),kth(2),length
		class(*),target::inA(:)
		DATATYPE,pointer::outp(:,:,:)
		class(*),pointer::inA2(:,:,:)
		if((dim1*dim2*dim3).ne.length)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(ith(2).gt.dim1)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(jth(2).gt.dim2)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(kth(2).gt.dim3)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(ith(1).le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(jth(1).le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(kth(1).le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				inA2(1:dim1,1:dim2,1:dim3)=>inA(1:length)
				select type(inA2)
					type is (DATATYPE2)
						outp=>inA2(ith(1):ith(2),jth(1):jth(2),kth(1):kth(2))
					class default
						call writemess('ERROR type  in pointting',-1)
						call error_stop
				end select
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point3D_a_ValueFUNCNAME(length,inA,dim1,dim2,dim3,ith,jth,kth,outp)
		integer::dim1,dim2,dim3,ith,jth,kth,length
		class(*),target::inA(:)
		DATATYPE,pointer::outp
		class(*),pointer::inA2(:,:,:)
		if((dim1*dim2*dim3).ne.length)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(ith.gt.dim1)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(jth.gt.dim2)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(kth.gt.dim3)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(jth.le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		if(kth.le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (DATATYPE2)
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				inA2(1:dim1,1:dim2,1:dim3)=>inA(1:length)

				select type(inA2)
					type is (DATATYPE2)
						outp=>inA2(ith,jth,kth)
					class default
						call writemess('ERROR type  in pointting',-1)
						call error_stop
				end select
				
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point3D_All_Value2FUNCNAME(inA,outp)
		class(*),target::inA(:,:,:)
		DATATYPE,pointer::outp(:,:,:)
		select type(inA)
			type is (DATATYPE2)
				outp=>inA
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point3D_Some_Value2FUNCNAME(inA,ith,jth,kth,outp)
		integer::ith(2),jth(2),kth(2)
		class(*),target::inA(:,:,:)
		DATATYPE,pointer::outp(:,:,:)
		select type(inA)
			type is (DATATYPE2)
				outp=>inA(ith(1):ith(2),jth(1):jth(2),kth(1):kth(2))
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point3D_a_Value2FUNCNAME(inA,ith,jth,kth,outp)
		integer::ith,jth,kth
		class(*),target::inA(:,:,:)
		DATATYPE,pointer::outp
		select type(inA)
			type is (DATATYPE2)
				if(ith.gt.size(inA,1))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(jth.gt.size(inA,2))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(kth.gt.size(inA,3))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(ith.le.0)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(jth.le.0)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(kth.le.0)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				outp=>inA(ith,jth,kth)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine
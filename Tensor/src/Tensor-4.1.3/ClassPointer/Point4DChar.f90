	subroutine Point4D_All_Valuea2(length,inA,dim1,dim2,dim3,dim4,outp,chalen)
		integer::dim1,dim2,dim3,dim4,length,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp(:,:,:,:)
		if((dim1*dim2*dim3*dim4).ne.length)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp(1:dim1,1:dim2,1:dim3,1:dim4)=>inA(1:length)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point4D_Some_Valuea2(length,inA,dim1,dim2,dim3,dim4,ith,jth,kth,lth,outp,chalen)
		integer::dim1,dim2,dim3,dim4,ith(2),jth(2),kth(2),lth(2),length,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp(:,:,:,:)
		character(len=chalen),pointer::inA2(:,:,:,:)
		if((dim1*dim2*dim3*dim4).ne.length)then
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
		if(lth(2).gt.dim4)then
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
		if(lth(1).le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				inA2(1:dim1,1:dim2,1:dim3,1:dim4)=>inA(1:length)
				outp=>inA2(ith(1):ith(2),jth(1):jth(2),kth(1):kth(2),lth(1):lth(2))
				
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point4D_a_Valuea2(length,inA,dim1,dim2,dim3,dim4,ith,jth,kth,lth,outp,chalen)
		integer::dim1,dim2,dim3,dim4,ith,jth,kth,lth,length,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp
		character(len=chalen),pointer::inA2(:,:,:,:)
		if((dim1*dim2*dim3*dim4).ne.length)then
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
		if(lth.gt.dim4)then
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
		if(lth.le.0)then
			call writemess('ERROR in pointing ',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				inA2(1:dim1,1:dim2,1:dim3,1:dim4)=>inA(1:length)
				outp=>inA2(ith,jth,kth,lth)
				
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point4D_All_Value2a2(inA,outp,chalen)
		integer,intent(in)::chalen
		class(*),target::inA(:,:,:,:)
		character(len=chalen),pointer::outp(:,:,:,:)
		select type(inA)
			type is (character(len=*))
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp=>inA
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point4D_Some_Value2a2(inA,ith,jth,kth,lth,outp,chalen)
		integer::ith(2),jth(2),kth(2),lth(2),chalen
		class(*),target::inA(:,:,:,:)
		character(len=chalen),pointer::outp(:,:,:,:)
		select type(inA)
			type is (character(len=*))
				if(ith(2).gt.size(inA,1))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(jth(2).gt.size(inA,2))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(kth(2).gt.size(inA,3))then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(lth(2).gt.size(inA,4))then
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
				if(lth(1).le.0)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp=>inA(ith(1):ith(2),jth(1):jth(2),kth(1):kth(2),lth(1):lth(2))
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point4D_a_Value2a2(inA,ith,jth,kth,lth,outp,chalen)
		integer::ith,jth,kth,lth,chalen
		class(*),target::inA(:,:,:,:)
		character(len=chalen),pointer::outp
		select type(inA)
			type is (character(len=*))
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
				if(lth.gt.size(inA,4))then
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
				if(lth.le.0)then
					call writemess('ERROR in pointing ',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp=>inA(ith,jth,kth,lth)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine
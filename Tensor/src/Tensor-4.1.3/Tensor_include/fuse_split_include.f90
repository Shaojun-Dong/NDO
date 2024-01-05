	subroutine ReshapeClassAssignment1D(outA,outith,inA,inith)
		class(*)::outA(:)
		class(*)::inA(:)
		class(*),pointer::pout(:),pin(:)
		integer::outith(2),inith(2)
		call ClassPointer1DFunc(outA,outith(1),outith(2),pout)
		call ClassPointer1DFunc(inA,inith(1),inith(2),pin)
		call FastCopyArray(pout,pin,classSize(inA))
		return
	end subroutine
	subroutine ReshapeClassAssignment2D(outA,outith,outjth,inA,inith,injth)
		class(*)::outA(:,:)
		class(*)::inA(:,:)
		integer::outith(2),outjth(2),inith(2),injth(2)
		class(*),pointer::pin(:,:),pout(:,:)
		integer::dim1,dim2
		call ClassPointer2DFunc(outA,outith,outjth,pout)
		call ClassPointer2DFunc(inA,inith,injth,pin)
		dim1=classSize(pin,1)
		dim2=classSize(pin,2)
		call CopyArray(pout,pin,dim1,dim2)
		return
	end subroutine
	subroutine ReshapeClassAssignment3D(outA,outith,outjth,outkth,inA,inith,injth,inkth)
		class(*)::outA(:,:,:)
		class(*)::inA(:,:,:)
		integer::outith(2),outjth(2),outkth(2),inith(2),injth(2),inkth(2)
		class(*),pointer::pin(:,:,:),pout(:,:,:)
		integer::dim1,dim2,dim3
		call ClassPointer3DFunc(outA,outith,outjth,outkth,pout)
		call ClassPointer3DFunc(inA,inith,injth,inkth,pin)
		dim1=classSize(pin,1)
		dim2=classSize(pin,2)
		dim3=classSize(pin,3)
		call CopyArray(pout,pin,dim1,dim2,dim3)
		return
	end subroutine

	subroutine QuantumFuseClassData(A,ith,jth,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith,jth
		type(Tensor),intent(inout)::Res
		type(Dimension),optional,intent(inout)::SubDim
		type(Tensor),optional,intent(inout)::ResOrder
		class(*),pointer::Rp(:),Ap(:),Rp2(:,:),Rp3(:,:,:),Ap2(:,:),Ap3(:,:,:)
		type(QuanNum)::NewQuanNum
		integer::rankA,caseFlag,NewRank,i,Rsi,Rei,Orderrow,lenindex,NewQNlength,TotalBlock,inforow
		integer::Rdim1,RDim2,RDim3,RI,RJ,RK,RDegI,RDegJ,RDegK
		integer::ADegI,ADegJ,ADegK,Aindex
		integer,pointer::Dim(:)
		integer,pointer::indices(:),NewDeg(:),order(:,:),fuseinfo(:,:),BlockN(:)
		real*4,pointer::TMPQN(:),NewQN(:)
		class(*),pointer::mold

		
		rankA=A%getRank()
		call A%pointDim(dim)
		TotalBlock=A%getTotalBlock()
		lenindex=jth-ith+1
		NewQNlength=product(dim(ith:jth))
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rankA+NewQNlength+TotalBlock*6+TotalBlock*10+TotalBlock*2+rankA)
		call WorkingMemory%allocate(2,NewQNlength+lenindex)
		call WorkingMemory%get_memory(indices,lenindex)
		call WorkingMemory%get_memory(NewDeg,NewQNlength)
		call WorkingMemory%get_memory(NewQN,NewQNlength)
		call WorkingMemory%get_memory(TMPQN,lenindex)
		call WorkingMemory%get_memory(order,TotalBlock,6)
		call WorkingMemory%get_memory(fuseinfo,TotalBlock,10)
		call FuseDimension(NewQuanNum,Order,A%Dimension,ith,jth,Orderrow,indices,NewDeg,NewQN,TMPQN)

		if((ith.eq.1).and.(jth.eq.rankA))then
			Res%Dimension=[NewQuanNum]
			caseFlag=1
		else if(ith.eq.1)then
			call pasteDimension(Res%Dimension,NewQuanNum,A%dimension,jth+1,rankA)
			caseFlag=2
		else if(jth.eq.rankA)then
			call pasteDimension(Res%Dimension,A%dimension,1,ith-1,NewQuanNum)
			caseFlag=3
		else
			call pasteDimension(Res%Dimension,A%dimension,1,ith-1,NewQuanNum,A%dimension,jth+1,rankA)
			caseFlag=4
		end if
		call WorkingMemory%get_memory(indices,rankA)
		call ClassPointer1DFunc(A%Data%ClassData,1,mold)
		
		
		select case(caseFlag)
			case(1)![1:rank]
				call Res%pointDim(dim)
				call WorkingMemory%get_memory(BlockN,dim(1))
				call fuseinfo1(fuseinfo,BlockN,Order,A,ith,jth,inforow)
				NewRank=Res%getRank()
				call Res%Data%allocateClassType(BlockN,mold,A%getType())
				call Res%zero()
				Rdim1=product(dim)
				!outinfo: old_index,Deg_i1i2i3,I,DegI,si,ei
				!              1         2     3   4  5 ,6
				do i=1,inforow
					Aindex=fuseinfo(i,1)
					ADegI=fuseinfo(i,2)
					RI=fuseinfo(i,3)
					RDegI=fuseinfo(i,4)
					Rsi=fuseinfo(i,5)
					Rei=fuseinfo(i,6)
					if(.not.Res%getFlag(RI))then
						call writemess('ERROR in fuse tensor',-1)
						call error_stop
					end if
					call Res%ClassPointer(Rp,RI)
					call A%ClassPointer(Ap,Aindex)
					if((Rei-Rsi+1).ne.ADegI)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call error_stop
					end if
					if(size(Rp).lt.Rei)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call writemess('size(Rp)='+size(Rp),-1)
						call writemess('Rei='+Rei,-1)
						call writemess('RI='+RI,-1)
						call writemess('Aindex='+Aindex,-1)
						call writemess('RDegI='+RDegI,-1)
						call writemess('ADegI='+ADegI,-1)
					end if
					call ReshapeClassAssignment1D(Rp,[Rsi,Rei],Ap,[1,ADegI])
				end do
			case(2)![1:jth,k,l,n....]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=dim(1)
				Rdim2=product(dim(2:NewRank))
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2)
				call fuseinfo2(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,inforow,indices)
				call Res%Data%allocateClassType(BlockN,mold,A%getType())
				call Res%zero()
				!outinfo: old_index,Deg_i1i2i3,I,j,DegI,DegJ,si,ei
				!              1     2         3 4   5  6    7  8
				do i=1,inforow
					Aindex=fuseinfo(i,1)
					ADegI=fuseinfo(i,2)
					RI=fuseinfo(i,3)
					RJ=fuseinfo(i,4)
					RDegI=fuseinfo(i,5)
					ADegJ=fuseinfo(i,6)
					RDegJ=ADegJ
					Rsi=fuseinfo(i,7)
					Rei=fuseinfo(i,8)
					
					if(.not.Res%getFlag([RDim1,RDim2],[RI,RJ]))then
						call writemess('ERROR in fuse tensor',-1)
						call error_stop
					end if
					
					call A%ClassPointer(Ap,Aindex)
					call Res%ClassPointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
					if((Rei-Rsi+1).ne.ADegI)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call error_stop
					end if
					if(size(Rp2,1).lt.Rei)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call writemess('size(Rp2,1)='+size(Rp2,1),-1)
						call writemess('Rei='+Rei,-1)
						call writemess('RDegI='+RDegI,-1)
						call writemess('RDegJ='+RDegJ,-1)
						call writemess('RI='+RI,-1)
						call writemess('RJ='+RJ,-1)
						call error_stop
					end if
					!Ap2(1:ADegI,1:ADegJ)=>Ap
					call ClassPointer2DFunc(ClassSize(Ap),Ap,ADegI,ADegJ,Ap2)
					call ReshapeClassAssignment2D(Rp2,[Rsi,Rei],[1,ADegJ],Ap2,[1,ADegI],[1,ADegJ])
				end do
			case(3)![1,2,3,...,ith:rank]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=product(dim(1:NewRank-1))
				Rdim2=dim(NewRank)
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2)
				call fuseinfo2Rank(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,inforow,indices)
				call Res%Data%allocateClassType(BlockN,mold,A%getType())
				call Res%zero()
				!outinfo: old_index,Deg_j1j2j3,i,J,DegI,DegJ,si,ei
				!              1         2     3 4  5    6    7  8
 				do i=1,inforow
 					Aindex=fuseinfo(i,1)
 					ADegJ=fuseinfo(i,2)
					RI=fuseinfo(i,3)
					RJ=fuseinfo(i,4)
					RDegI=fuseinfo(i,5)
					ADegI=RDegI
					RDegJ=fuseinfo(i,6)
					Rsi=fuseinfo(i,7)
					Rei=fuseinfo(i,8)
					
					if(.not.Res%getFlag([RDim1,RDim2],[RI,RJ]))then
						call writemess('ERROR in fuse tensor',-1)
						call error_stop
					end if
					
					
					call A%ClassPointer(Ap,Aindex)
					call Res%ClassPointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
					if((Rei-Rsi+1).ne.ADegJ)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call error_stop
					end if
					if(size(Rp2,2).lt.Rei)then
						call writemess('ERROR in QuantumFuse, case 1',-1)
						call writemess('size(Rp2,2)='+size(Rp2,2),-1)
						call writemess('Rei='+Rei,-1)
						call writemess('RDegI='+RDegI,-1)
						call writemess('RDegJ='+RDegJ,-1)
						call writemess('RI='+RI,-1)
						call writemess('RJ='+RJ,-1)
						call error_stop
					end if
					!Ap2(1:ADegI,1:ADegJ)=>Ap
					call ClassPointer2DFunc(ClassSize(Ap),Ap,ADegI,ADegJ,Ap2)
					call ReshapeClassAssignment2D(Rp2,[1,ADegI],[Rsi,Rei],Ap2,[1,ADegI],[1,ADegJ])
				end do
				
			case(4)![1,2,3,...,ith:jth,k,l,m,n]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=product(dim(1:ith-1))
				Rdim2=dim(ith)
				Rdim3=product(dim(ith+1:NewRank))
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2*Rdim3)
				call fuseinfo3(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,Rdim3,inforow,indices)
				call Res%Data%allocateClassType(BlockN,mold,A%getType())
				call Res%zero()
				!outinfo: old_index,Deg_j1j2j3,i,J,k,DegI,DegJ,DegK,si,ei
				!              1         2     3 4 5  6    7    8   9  10
				do i=1,inforow
					Aindex=fuseinfo(i,1)
					ADegJ=fuseinfo(i,2)
					RI=fuseinfo(i,3)
					RJ=fuseinfo(i,4)
					RK=fuseinfo(i,5)
					RDegI=fuseinfo(i,6)
					ADegI=RDegI
					RDegJ=fuseinfo(i,7)
					RDegK=fuseinfo(i,8)
					ADegK=RDegK
					Rsi=fuseinfo(i,9)
					Rei=fuseinfo(i,10)
					
					if(.not.Res%getFlag([RDim1,RDim2,Rdim3],[RI,RJ,RK]))then
						call writemess('ERROR in fuse tensor',-1)
						call error_stop
					end if
					
					call A%ClassPointer(Ap,Aindex)
					call Res%ClassPointer([RDim1,RDim2,Rdim3],Rp3,[RDegI,RDegJ,RDegK],[RI,RJ,RK])
					if((Rei-Rsi+1).ne.ADegJ)then
						call writemess('ERROR in QuantumFuse, case 4',-1)
						call error_stop
					end if
					if(size(Rp3,2).lt.Rei)then
						call writemess('ERROR in QuantumFuse, case 4',-1)
						call writemess('size(Rp3,2)='+size(Rp3,2),-1)
						call writemess('Rei='+Rei,-1)
						call writemess('RDegI='+RDegI,-1)
						call writemess('RDegJ='+RDegJ,-1)
						call writemess('RDegK='+RDegK,-1)
						call writemess('RI='+RI,-1)
						call writemess('RJ='+RJ,-1)
						call writemess('RK='+RK,-1)
						call error_stop
					end if
					!Ap3(1:ADegI,1:ADegJ,1:ADegK)=>Ap
					call ClassPointer3DFunc(ClassSize(Ap),Ap,ADegI,ADegJ,ADegK,Ap3)
					call ReshapeClassAssignment3D(Rp3,[1,AdegI],[Rsi,Rei],[1,ADegK],&
													Ap3,[1,ADegI],[1,ADegJ],[1,ADegK])
				end do
		end select

		if(present(ResOrder))then
			ResOrder=Order(1:Orderrow,:)
			if(Orderrow.eq.1)then
				call ResOrder%resetDim([1,ResOrder%getTotalData()])
			end if
		end if
		
		if(present(SubDim))then
			call getsubDimension(A%dimension,[ith,jth],SubDim)
		end if
		call WorkingMemory%free()
		return
	end subroutine

	subroutine QuantumSplitClassData(A,legi,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::legi
		type(Tensor),intent(inout)::Res
		type(Dimension),intent(in)::SubDim
		type(Tensor),intent(in)::ResOrder
		class(*),pointer::Rp(:),Ap(:),Ap2(:,:),Ap3(:,:,:),Rp2(:,:),Rp3(:,:,:)
		integer::NewRank,RankA,len_index_in_used
		integer,pointer::Order(:,:),dim(:),outindex(:),Deg(:),indices(:)
		integer::Si,Ei,caseFlag,len_of_new_leg,i,j,k,ii,iii
		integer::Rdim1,RDim2,RDim3,RI,RJ,RK,RDegI,RDegJ,RDegK
		integer::Adim1,ADim2,ADim3,AI,AJ,AK,ADegI,ADegJ,ADegK
		logical::goon
		class(*),pointer::mold

		call ClassPointer1DFunc(A%Data%ClassData,1,mold)

		call insertDimension(A%Dimension,SubDim,legi,Res%Dimension)
		len_of_new_leg=SubDim%getRank()
		call Res%pointDim(dim)
		call Res%Data%allocateMomeryClassType(product(dim),A%getTotalData(),mold,A%getType())
		call ResOrder%pointer(Order)
		RankA=A%getRank()
		NewRank=Res%getRank()
		if(RankA.eq.1)then
			caseFlag=1
		else if(legi.eq.1) then
			caseFlag=2
		else if(legi.eq.RankA) then
			caseFlag=3
		else
			caseFlag=4
		end if

		len_index_in_used=ResOrder%dim(1)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,len_index_in_used+rankA)
		call WorkingMemory%get_memory(outindex,len_index_in_used)
		call WorkingMemory%get_memory(indices,rankA)

		
		select case(caseFlag)
			case(1)![1:rank],AI=0,AK=0
				call Res%pointDim(Dim)
				Rdim1=product(Dim)
				call A%pointDim(Dim)
				Adim1=Dim(1)
				do i=1,ADim1
					!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
					!             1        2        3     4        5     6
					if(A%getFlag(i))then
						call Find_Data_in_order(outindex,order,3,i,len_index_in_used)
						do ii=1,len_index_in_used
							AI=order(outindex(ii),3)
							RI=order(outindex(ii),1)
							RDegI=order(outindex(ii),2)
							ADegI=order(outindex(ii),4)
							Si=order(outindex(ii),5)
							Ei=order(outindex(ii),6)
							if(.not.Res%getFlag(RI))then
								call Res%setBlockMomery([Rdim1],[RI],RDegI)
								call Res%ClassPointer(Rp,RI)
								call A%ClassPointer(Ap,AI)
								if((Ei-Si+1).ne.RDegI)then
									call writemess('ERROR in QuantumFuse, case 1',-1)
									call error_stop
								end if
								if(size(Ap).lt.Ei)then
									call writemess('ERROR in QuantumFuse, case 1',-1)
									call error_stop
								end if
								!Rp=Ap(Si:Ei)
								call ReshapeClassAssignment1D(Rp,[1,RDegI],Ap,[Si,Ei])
							end if
						end do
					end if
				end do
			case(2)![1:jth,k,l,n....] AI=0, AK!=0
				call Res%pointDim(Dim)
				Rdim1=product(Dim(1:len_of_new_leg))
				Rdim2=product(Dim(len_of_new_leg+1:NewRank))

				call A%pointDim(Dim)
				Adim1=Dim(1)
				Adim2=product(Dim(2:rankA))
				indices=1
				!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
				!             1        2        3     4        5     6
				do j=1,Adim2
					do i=1,Adim1
						if(.not.goon)then
							call writemess('ERROR in QuantumSplitDataTYPE',-1)
							call error_stop
						end if
						if(A%getFlag([Adim1,Adim2],[i,j]))then
							call Find_Data_in_order(outindex,order,3,i,len_index_in_used)
							do ii=1,len_index_in_used
								AI=order(outindex(ii),3)
								AJ=j
								RI=order(outindex(ii),1)
								RJ=j
								RDegI=order(outindex(ii),2)
								ADegI=order(outindex(ii),4)
								Si=order(outindex(ii),5)
								Ei=order(outindex(ii),6)
								call A%pointDeg(Deg,2)
								ADegJ=Deg(indices(2))
								do iii=3,rankA
									call A%pointDeg(Deg,iii)
									ADegJ=ADegJ*Deg(indices(iii))
								end do
								RDegJ=ADegJ
								if(.not.Res%getFlag([Rdim1,Rdim2],[RI,RJ]))then
									call Res%setBlockMomery([Rdim1,Rdim2],[RI,RJ],RDegI*RDegJ)
									call Res%ClassPointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
									call A%ClassPointer([ADim1,ADim2],Ap2,[ADegI,ADegJ],[AI,AJ])
									if((Ei-Si+1).ne.RDegI)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									if(size(Ap2,1).lt.Ei)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									!Rp2=Ap2(Si:Ei,:)
									call ReshapeClassAssignment2D(Rp2,[1,RDegI],[1,RDegJ],&
																Ap2,[Si,Ei],[1,ADegJ])
								end if
							end do
						end if
						goon=index_counter(indices,dim)
					end do
				end do
			case(3)![1,2,3,...,ith:rank] AI!=0,AK=0
				call Res%pointDim(Dim)
				Rdim1=product(Dim(1:legi-1))
				Rdim2=product(Dim(legi:NewRank))

				call A%pointDim(Dim)
				Adim1=product(Dim(1:rankA-1))
				Adim2=Dim(rankA)
				indices=1
				!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
				!             1        2        3     4        5     6
				do j=1,Adim2
					do i=1,Adim1
						if(.not.goon)then
							call writemess('ERROR in QuantumSplitDataTYPE',-1)
							call error_stop
						end if
						if(A%getFlag([Adim1,Adim2],[i,j]))then
							call Find_Data_in_order(outindex,order,3,j,len_index_in_used)
							do ii=1,len_index_in_used
								AI=i
								AJ=order(outindex(ii),3)
								RI=i
								RJ=order(outindex(ii),1)
								RDegJ=order(outindex(ii),2)
								ADegJ=order(outindex(ii),4)
								Si=order(outindex(ii),5)
								Ei=order(outindex(ii),6)
								call A%pointDeg(Deg,1)
								ADegI=Deg(indices(1))
								do iii=2,rankA-1
									call A%pointDeg(Deg,iii)
									ADegI=ADegI*Deg(indices(iii))
								end do
								RDegI=ADegI
								if(.not.Res%getFlag([Rdim1,Rdim2],[RI,RJ]))then
									call Res%setBlockMomery([Rdim1,Rdim2],[RI,RJ],RDegI*RDegJ)
									call Res%ClassPointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
									call A%ClassPointer([ADim1,ADim2],Ap2,[ADegI,ADegJ],[AI,AJ])
									if((Ei-Si+1).ne.RDegJ)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									if(size(Ap2,2).lt.Ei)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									!Rp2=Ap2(:,Si:Ei)
									call ReshapeClassAssignment2D(Rp2,[1,RDegI],[1,RDegJ],&
																Ap2,[1,ADegI],[Si,Ei])
								end if
							end do
						end if
						goon=index_counter(indices,dim)
					end do
				end do
			case(4)![1,2,3,...,ith:jth,k,l,m,n]
				call Res%pointDim(Dim)
				Rdim1=product(Dim(1:legi-1))
				Rdim2=product(Dim(legi:legi+len_of_new_leg))
				Rdim3=product(Dim(legi+len_of_new_leg+1:NewRank))

				call A%pointDim(Dim)
				Adim1=product(Dim(1:legi-1))
				Adim2=Dim(legi)
				Adim3=product(Dim(legi+1:rankA))
				indices=1
				!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
				!             1        2        3     4        5     6
				do k=1,Adim3
					do j=1,Adim2
						do i=1,Adim1
							if(.not.goon)then
								call writemess('ERROR in QuantumSplitDataTYPE',-1)
								call error_stop
							end if
							if(A%getFlag([Adim1,Adim2,Adim3],[i,j,k]))then
								call Find_Data_in_order(outindex,order,3,j,len_index_in_used)
								do ii=1,len_index_in_used
									AI=i
									AJ=order(outindex(ii),3)
									AK=k
									RI=i
									RJ=order(outindex(ii),1)
									RK=k
									RDegJ=order(outindex(ii),2)
									ADegJ=order(outindex(ii),4)
									Si=order(outindex(ii),5)
									Ei=order(outindex(ii),6)
									call A%pointDeg(Deg,1)
									ADegI=Deg(indices(1))
									do iii=2,legi-1
										call A%pointDeg(Deg,iii)
										ADegI=ADegI*Deg(indices(iii))
									end do
									RDegI=ADegI

									call A%pointDeg(Deg,legi+1)
									ADegK=Deg(indices(legi+1))
									do iii=legi+2,rankA
										call A%pointDeg(Deg,iii)
										ADegK=ADegK*Deg(indices(iii))
									end do
									RDegK=ADegK

									if(.not.Res%getFlag([Rdim1,Rdim2,RDim3],[RI,RJ,Rk]))then
										call Res%setBlockMomery([Rdim1,Rdim2,RDim3],[RI,RJ,RK],RDegI*RDegJ*RDegK)
										call Res%ClassPointer([RDim1,RDim2,RDim3],Rp3,[RDegI,RDegJ,RDegK],[RI,RJ,RK])
										call A%ClassPointer([ADim1,ADim2,ADim3],Ap3,[ADegI,ADegJ,ADegK],[AI,AJ,AK])
										if((Ei-Si+1).ne.RDegJ)then
											call writemess('ERROR in QuantumFuse, case 1',-1)
											call error_stop
										end if
										if(size(Ap2,2).lt.Ei)then
											call writemess('ERROR in QuantumFuse, case 1',-1)
											call error_stop
										end if
										!Rp3=Ap3(:,Si:Ei,:)
										call ReshapeClassAssignment3D(Rp3,[1,RDegI],[1,RDegJ],[1,RDegK],&
																	Ap3,[1,ADegI],[Si,Ei],[2,ADegK])
									end if
								end do
							end if
							goon=index_counter(indices,dim)
						end do
					end do
				end do
		end select
		call WorkingMemory%free()
		if(Res%getTotalData().ne.A%getTotalData())then
			call writemess('ERROR in split',-1)
			call writemess('Res%getTotalData()='+Res%getTotalData(),-1)
			call writemess('A%getTotalData()='+A%getTotalData(),-1)
			call Res%print()
			call writemess('--------------------')
			call A%print()
			call error_stop
		end if
		return
	end subroutine
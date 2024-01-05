	!QuanNum1+QuanNum2=NewQuanNum
		!
		!order of codes for do
		! outQ=Q1+Q2 
		!  first order Q1 and then Q2 at last outQ
		!Parity example:
		!order of codes for do
		! outQ      Q1  Q2 
		!  -1       -1  -1   false
		!  -1       +1  -1
		!  -1       -1  +1
		!  -1       +1  +1   false
		!  +1       -1  -1
		!  +1       +1  -1   false
		!  +1       -1  +1   false
		!  +1       +1  +1
		!
		!The U(1) symmetry
		!the rule for outQ Q1 Q2 are -1 1 1
		! outQ      Q1   Q2 
		! -1.5      -1  -0.5
		! -1.5       0  -0.5 false
		! -1.5       1  -0.5 false
		! -1.5      -1   0.5 false
		! -1.5       0   0.5 false
		! -1.5       1   0.5 false
		!
		! -0.5      -1  -0.5 false
		! -0.5       0  -0.5 
		! -0.5       1  -0.5 false
		! -0.5      -1   0.5 
		! -0.5       0   0.5 false
		! -0.5       1   0.5 false
		!
		!  0.5      -1  -0.5 false
		!  0.5       0  -0.5 false
		!  0.5       1  -0.5 
		!  0.5      -1   0.5 false
		!  0.5       0   0.5 
		!  0.5       1   0.5 false
		!
		!  1.5      -1  -0.5 false
		!  1.5       0  -0.5 false
		!  1.5       1  -0.5 false
		!  1.5      -1   0.5 false
		!  1.5       0   0.5 false
		!  1.5       1   0.5 


	subroutine QuantumFuseTYPE(A,ith,jth,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith,jth
		type(Tensor),intent(inout)::Res
		type(Dimension),optional,intent(inout)::SubDim
		type(Tensor),optional,intent(inout)::ResOrder
		DataTYPE,pointer::Rp(:),Ap(:),Rp2(:,:),Rp3(:,:,:),Ap2(:,:),Ap3(:,:,:)
		type(QuanNum)::NewQuanNum
		integer::rankA,caseFlag,NewRank,i,Rsi,Rei,Orderrow,lenindex,NewQNlength,TotalBlock,inforow
		integer::Rdim1,RDim2,RDim3,RI,RJ,RK,RDegI,RDegJ,RDegK
		integer::ADegI,ADegJ,ADegK,Aindex,ii,jj,kk,ll
		integer,pointer::Dim(:)
		integer,pointer::indices(:),NewDeg(:),order(:,:),fuseinfo(:,:),BlockN(:)
		real*4,pointer::TMPQN(:),NewQN(:)

		
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
		
		
		
		select case(caseFlag)
			case(1)![1:rank]
				call Res%pointDim(dim)
				call WorkingMemory%get_memory(BlockN,dim(1))
				call fuseinfo1(fuseinfo,BlockN,Order,A,ith,jth,inforow)
				NewRank=Res%getRank()
				call Res%Data%allocate(BlockN,A%getType())
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
					call Res%pointer(Rp,RI)
					call A%pointer(Ap,Aindex)
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
					!Rp(Rsi:Rei)=Ap

					kk=0
					do ii=Rsi,Rei
						kk=kk+1
						Rp(ii)=Ap(kk)
					end do


				end do
			case(2)![1:jth,k,l,n....]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=dim(1)
				Rdim2=product(dim(2:NewRank))
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2)
				call fuseinfo2(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,inforow,indices)
				call Res%Data%allocate(BlockN,A%getType())
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
					
					call A%pointer(Ap,Aindex)
					call Res%pointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
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
					Ap2(1:ADegI,1:ADegJ)=>Ap
					!Rp2(Rsi:Rei,:)=Ap2
					do jj=1,ADegJ
						kk=0
						do ii=Rsi,Rei
							kk=kk+1
							Rp2(ii,jj)=Ap2(kk,jj)
						end do
					end do

				end do
			case(3)![1,2,3,...,ith:rank]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=product(dim(1:NewRank-1))
				Rdim2=dim(NewRank)
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2)
				call fuseinfo2Rank(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,inforow,indices)
				call Res%Data%allocate(BlockN,A%getType())
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
					
					
					call A%pointer(Ap,Aindex)
					call Res%pointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
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
					Ap2(1:ADegI,1:ADegJ)=>Ap
					!Rp2(:,Rsi:Rei)=Ap2


					kk=0
					do jj=Rsi,Rei
						kk=kk+1
						do ii=1,ADegI
							Rp2(ii,jj)=Ap2(ii,kk)
						end do
					end do

				end do
				
			case(4)![1,2,3,...,ith:jth,k,l,m,n]
				
				NewRank=Res%getRank()
				call Res%pointDim(dim)
				Rdim1=product(dim(1:ith-1))
				Rdim2=dim(ith)
				Rdim3=product(dim(ith+1:NewRank))
				call WorkingMemory%get_memory(BlockN,Rdim1*Rdim2*Rdim3)
				call fuseinfo3(fuseinfo,BlockN,Order,A,ith,jth,Rdim1,Rdim2,Rdim3,inforow,indices)
				call Res%Data%allocate(BlockN,A%getType())
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
					
					call A%pointer(Ap,Aindex)
					call Res%pointer([RDim1,RDim2,Rdim3],Rp3,[RDegI,RDegJ,RDegK],[RI,RJ,RK])
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
					Ap3(1:ADegI,1:ADegJ,1:ADegK)=>Ap
					!Rp3(:,Rsi:Rei,:)=Ap3
					
					do kk=1,ADegK
						ll=0
						do jj=Rsi,Rei
							ll=ll+1
							do ii=1,ADegI
								Rp3(ii,jj,kk)=Ap3(ii,ll,kk)
							end do
						end do
					end do
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


	!case 1 example
		!    old dim     --->   new dim  
		!    [1 1 1]        /  [1 1 1]  \     dim(1)
		!    [2 1 1]        |-----------|     
		!    [1 2 1]        |  [2 1 1]  |
		!    [2 2 1]        |  [1 2 1]  |     dim(2)
		!    [1 1 2]        |  [1 1 2]  |
		!    [2 1 2]        |-----------|  
		!    [1 2 2]        |  [2 2 1]  |
		!    [2 2 2]        |  [2 1 2]  |     dim(3)
		!                   |  [2 2 1]  |
		!                   |-----------|
		!                   \  [2 2 2]  /     dim(4)

	!case 2 example
		!    old dim       --->   new dim  
		!    [1 1 1 1]  
		!    [2 1 1 1]           /  [1 1 1 1]   :   [1 1 1 2]  \
		!    [1 2 1 1]           |--------------:--------------|  
		!    [2 2 1 1]           |  [2 1 1 1]   :   [2 1 1 2]  |
		!    [1 1 2 1]           |  [1 2 1 1]   :   [1 2 1 1]  |
		!    [2 1 2 1]           |  [1 1 2 1]   :   [1 1 2 1]  |
		!    [1 2 2 1]           |--------------:--------------|    
		!    [2 2 2 1]           |  [2 2 1 1]   :   [2 2 1 2]  |
		!    [1 1 1 2]           |  [2 1 2 1]   :   [2 1 2 1]  | 
		!    [2 1 1 2]           |  [2 2 1 1]   :   [2 2 1 1]  |
		!    [1 2 1 2]           |--------------:--------------|
		!    [2 2 1 2]           \  [2 2 2 1]   :   [2 2 2 1]  /
		!    [1 1 2 2]   
		!    [2 1 2 2]  
		!    [1 2 2 2]   
		!    [2 2 2 2]   
		!          

	subroutine QuantumSplitTYPE(A,legi,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::legi
		type(Tensor),intent(inout)::Res
		type(Dimension),intent(in)::SubDim
		type(Tensor),intent(in)::ResOrder
		DataTYPE,pointer::Rp(:),Ap(:),Ap2(:,:),Ap3(:,:,:),Rp2(:,:),Rp3(:,:,:)
		integer::NewRank,RankA,len_index_in_used
		integer,pointer::Order(:,:),dim(:),outindex(:),Deg(:),indices(:)
		integer::Si,Ei,caseFlag,len_of_new_leg,i,j,k,ii,iii
		integer::m,n,o,p
		integer::Rdim1,RDim2,RDim3,RI,RJ,RK,RDegI,RDegJ,RDegK
		integer::Adim1,ADim2,ADim3,AI,AJ,AK,ADegI,ADegJ,ADegK
		logical::goon

		call insertDimension(A%Dimension,SubDim,legi,Res%Dimension)
		len_of_new_leg=SubDim%getRank()
		call Res%pointDim(dim)
		call Res%Data%allocateDataArrayMomery(product(dim),A%getTotalData(),A%getType())
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
								call Res%pointer(Rp,RI)
								call A%pointer(Ap,AI)
								if((Ei-Si+1).ne.RDegI)then
									call writemess('ERROR in QuantumFuse, case 1',-1)
									call error_stop
								end if
								if(size(Ap).lt.Ei)then
									call writemess('ERROR in QuantumFuse, case 1',-1)
									call error_stop
								end if
								!Rp=Ap(Si:Ei)
								n=0
								do m=Si,Ei
									n=n+1
									Rp(n)=Ap(m)
								end do
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
							call writemess('ERROR in QuantumSplitTYPE',-1)
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
									call Res%pointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
									call A%pointer([ADim1,ADim2],Ap2,[ADegI,ADegJ],[AI,AJ])
									if((Ei-Si+1).ne.RDegI)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									if(size(Ap2,1).lt.Ei)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									!Rp2=Ap2(Si:Ei,:)
									o=0
									do m=Si,Ei
										o=o+1
										do n=1,ADegJ
											Rp2(o,n)=Ap2(m,n)
										end do
									end do
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
							call writemess('ERROR in QuantumSplitTYPE',-1)
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
									call Res%pointer([RDim1,RDim2],Rp2,[RDegI,RDegJ],[RI,RJ])
									call A%pointer([ADim1,ADim2],Ap2,[ADegI,ADegJ],[AI,AJ])
									if((Ei-Si+1).ne.RDegJ)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									if(size(Ap2,2).lt.Ei)then
										call writemess('ERROR in QuantumFuse, case 1',-1)
										call error_stop
									end if
									!Rp2=Ap2(:,Si:Ei)
									do m=1,ADegI
										o=0
										do n=Si,Ei
											o=o+1
											Rp2(m,o)=Ap2(m,n)
										end do
									end do
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
								call writemess('ERROR in QuantumSplitTYPE',-1)
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
										call Res%pointer([RDim1,RDim2,RDim3],Rp3,[RDegI,RDegJ,RDegK],[RI,RJ,RK])
										call A%pointer([ADim1,ADim2,ADim3],Ap3,[ADegI,ADegJ,ADegK],[AI,AJ,AK])
										if((Ei-Si+1).ne.RDegJ)then
											call writemess('ERROR in QuantumFuse, case 1',-1)
											call error_stop
										end if
										if(size(Ap2,2).lt.Ei)then
											call writemess('ERROR in QuantumFuse, case 1',-1)
											call error_stop
										end if
										!Rp3=Ap3(:,Si:Ei,:)
										do m=1,ADegI
											p=0
											do n=Si,Ei
												p=p+1
												do o=1,ADegK
													Rp3(m,p,o)=Ap3(m,n,o)
												end do
											end do
										end do
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
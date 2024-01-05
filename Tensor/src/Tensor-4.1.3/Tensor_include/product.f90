	! Res = alpha*TRANSA( A) * TRANSB( B)  + beta*Res

#include "SymProduct.f90"

#include "nonSymProduct.f90"
	subroutine pointToTensorWithRightType(Ap,Bp,A,B,TMPA,TMPB)
		type(Tensor),target,intent(in) :: A,B
		type(Tensor),pointer,intent(inout)::Ap,Bp
		type(Tensor),target,intent(inout) :: TMPA,TMPB
		if((.not.A%data%getFlag()).or.(.not.B%data%getFlag()))then
			Ap=>A
			Bp=>B
			return
		end if
		if(A%getType().eq.B%getType())then
			Ap=>A
			Bp=>B
			return
		end if
		call TMPA%empty()
		call TMPB%empty()
		if(A%getType().eq.1)then
			call TMPA%setType(B%getType())
			TMPA=A
			Ap=>TMPA
			Bp=>B
			call TMPA%setDynamic()
			return
		end if
		if(B%getType().eq.1)then
			Ap=>A
			call TMPB%setType(A%getType())
			TMPB=B
			Bp=>TMPB
			call TMPB%setDynamic()
			return
		end if
		if(A%getType().eq.2)then
			call TMPA%setType(B%getType())
			TMPA=A
			Ap=>TMPA
			Bp=>B
			call TMPA%setDynamic()
			return
		end if
		if(B%getType().eq.2)then
			Ap=>A
			call TMPB%setType(A%getType())
			TMPB=B
			Bp=>TMPB
			call TMPB%setDynamic()
			return
		end if
		if(A%getType().eq.3)then
			if(B%getType().eq.4)then
				call TMPA%setType(5)
				call TMPB%setType(5)
				TMPA=A
				TMPB=B
				Ap=>TMPA
				Bp=>TMPB
				call TMPA%setDynamic()
				call TMPB%setDynamic()
				return
			end if
			if(B%getType().eq.5)then
				call TMPA%setType(5)
				TMPA=A
				Ap=>TMPA
				Bp=>B
				call TMPA%setDynamic()
				return
			end if
		end if
		if(B%getType().eq.3)then
			if(A%getType().eq.4)then
				call TMPA%setType(5)
				call TMPB%setType(5)
				TMPA=A
				TMPB=B
				Ap=>TMPA
				Bp=>TMPB
				call TMPA%setDynamic()
				call TMPB%setDynamic()
				return
			end if
			if(A%getType().eq.5)then
				Ap=>A
				call TMPB%setType(5)
				TMPB=B
				Bp=>TMPB
				call TMPB%setDynamic()
				return
			end if
		end if
		if(A%getType().eq.4)then
			if(B%getType().eq.5)then
				call TMPA%setType(5)
				TMPA=A
				Ap=>TMPA
				Bp=>B
				call TMPA%setDynamic()
				return
			end if
		end if
		if(B%getType().eq.4)then
			if(A%getType().eq.5)then
				Ap=>A
				call TMPB%setType(5)
				TMPB=B
				Bp=>TMPB
				call TMPB%setDynamic()
				return
			end if
		end if
		call writemess('ERROR case in product',-1)
		call error_stop
	end subroutine
	subroutine ProductTensorRoutine(Res,A_,B_,alpha,beta,contractLegLen_)
		type(Tensor),target,intent(inout)::Res
		type(Tensor),target,intent(in) :: A_,B_
		class(*),intent(in)::alpha,beta
		integer,intent(in),optional::contractLegLen_
		integer::contractLegLen
		type(Tensor),pointer::Resp,A,B
		integer::i,rankA,rankB,iType,TotalDataA,TotalDataB,TotalBlockA,TotalBlockB,TotalBlockR,rankR
		integer::Flag,TotalMomeryLength
		integer,pointer::BlockDimA1(:),BlockDimA2(:),BlockDimB1(:),BlockDimB2(:),dimA(:),dimB(:)
		integer,pointer::BlockDimR1(:),BlockDimR2(:),dimR(:),ResBlockNum(:)
		integer,pointer::index(:)
		type(Dimension)::NewDim
		logical::NewResFlag
		if(present(contractLegLen_))then
			contractLegLen=contractLegLen_
		else
			contractLegLen=1
		end if
		if(.not.A_%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		if(.not.B_%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			call error_stop()
		end if	
		
		Resp=>Res
		A=>A_
		B=>B_
		if(associated(Resp,A).or.associated(Resp,B))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%ProductTensorRoutine(Res,T1,T2,alpha,beta)')
			call writemess('Res and T1, or Res and T2, can not be a same variable')
			call error_stop
		end if
		Resp=>null()
		call pointToTensorWithRightType(A,B,A_,B_,TMPproductA,TMPproductB)

		iType=A%getType()
		if(iType.ne.B%getType())then
			if((A%Data%getFlag()).and.(B%Data%getFlag()))then
				call writemess('ERROR in ProductTensorRoutine,Type',-1)
				call writemess('A%getType()='+A%getType(),-1)
				call writemess('B%getType()='+B%getType(),-1)
				call error_stop
			end if
		end if
		if(iType.ge.6)then
			call writemess('ERROR in ProductTensorRoutine,Type',-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if



		rankA=A%getRank()
		rankB=B%getRank()
		TotalDataA=A%getTotalData()
		TotalDataB=B%getTotalData()
		TotalBlockA=A%getTotalBlock()
		TotalBlockB=B%getTotalBlock()


		if(A%getSymmetryFlag().or.B%getSymmetryFlag())then
			if(.not.A%getSymmetryFlag())then
				if(A%getTotalData().eq.1)then
					select case(A%getType())
						case(1)
							Res=TmultiplyNumi(B,A%ii(1))
						case(2)
							Res=TmultiplyNums(B,A%si(1))
						case(3)
							Res=TmultiplyNumd(B,A%di(1))
						case(4)
							Res=TmultiplyNumc(B,A%ci(1))
						case(5)
							Res=TmultiplyNumz(B,A%zi(1))
						case default
							call writemess('ERROR data type in * ',-1)
							call error_stop
					end select
					return
				end if
				call writemess('ERROR in ProductTensorRoutine,SymmetryFlag',-1)
				call writemess('A%getSymmetryFlag()='+A%getSymmetryFlag(),-1)
				call writemess('B%getSymmetryFlag()='+B%getSymmetryFlag(),-1)
				call error_stop
			end if
			if(.not.B%getSymmetryFlag())then
				if(B%getTotalData().eq.1)then
					select case(B%getType())
						case(1)
							Res=TmultiplyNumi(A,B%ii(1))
						case(2)
							Res=TmultiplyNums(A,B%si(1))
						case(3)
							Res=TmultiplyNumd(A,B%di(1))
						case(4)
							Res=TmultiplyNumc(A,B%ci(1))
						case(5)
							Res=TmultiplyNumz(A,B%zi(1))
						case default
							call writemess('ERROR data type in * ',-1)
							call error_stop
					end select
					return
				end if
				call writemess('ERROR in ProductTensorRoutine,SymmetryFlag',-1)
				call writemess('A%getSymmetryFlag()='+A%getSymmetryFlag(),-1)
				call writemess('B%getSymmetryFlag()='+B%getSymmetryFlag(),-1)
				call error_stop
			end if
		end if
		if(A%getSymmetryFlag())then
			flag=SelectCaseForSym(rankA-contractLegLen+1,TotalDataA,rankB-contractLegLen+1,TotalDataB)
		else
			flag=SelectCase(rankA-contractLegLen+1,TotalDataA,rankB-contractLegLen+1,TotalDataB)
		end if

		NewResFlag=.not.Res%getFlag()

		select case(Flag)
			case(1) !Number*number,(1) *(1)
				if(.not.Res%getFlag())then
					call Res%allocate([1],iType)
				end if
				call ProductTensorCase1234(Res,A,B,alpha,beta,iType,NewResFlag)
			case(2)!Number*Tensor,(1) * (1,1,1)
				if(.not.Res%getFlag())then
					!NewDim=B%Dimension.subdim.[contractLegLen+1,rankB]
					call Res%allocate(B%Dimension,iType)
				end if
				call ProductTensorCase1234(Res,A,B,alpha,beta,iType,NewResFlag)
			case(3)!Tensor*Number, (1,1,1) * (1)
				if(.not.Res%getFlag())then
					!NewDim=A%Dimension.subdim.[1,rankA-contractLegLen]
					call Res%allocate(A%Dimension,iType)
				end if
				call ProductTensorCase1234(Res,A,B,alpha,beta,iType,NewResFlag)
			case(4)!Tensor*Tensor, (1,1,1) * (1,1)
				if(.not.Res%getFlag())then
					call pasteTwoSubDim(NewDim,A%Dimension,1,rankA-contractLegLen,&
							B%Dimension,contractLegLen+1,rankB)
					call Res%allocate(NewDim,iType)
				end if
				call ProductTensorCase1234(Res,A,B,alpha,beta,iType,NewResFlag)
			case(5)!Number*Tensor,(1) *(3) or (1) *(2,3)  or (1) * (1,2,2)
				if(.not.Res%getFlag())then
					!if(B%dim(1).eq.1)then
					!	call Res%allocate(B%Dimension.subdim.[contractLegLen+1,rankB],iType)
					!else
						call Res%allocate(B%Dimension,iType)
					!end if
				end if
				call ProductTensorCase56(Res,A,B,alpha,beta,iType,NewResFlag)
			case(6)!Tensor*Tensor,(1,1) *(1,2,1,2)
				if(B%dim(1).ne.1)then
					call writemess("ERROR in ProductTensor,case -2,stop",-1)
					call A%diminfo()
					call B%diminfo()
					call error_stop()
				end if
				if(rankB.eq.1)then!Tensor*number,(1,1) *(3)
					call writemess("ERROR in ProductTensor,case -1,stop",-1)
					call error_stop()
				else!Tensor*Tensor,(1,1) *(1,2,1,2)
					if(B%dim(1).ne.1)then
						call writemess("ERROR in ProductTensor,case -2,stop",-1)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop()
					end if
					if(.not.Res%getFlag())then
						call pasteTwoSubDim(NewDim,A%Dimension,1,rankA-contractLegLen,&
								B%Dimension,contractLegLen+1,rankB)
						call Res%allocate(NewDim,iType)
					end if
					call ProductTensorCase56(Res,A,B,alpha,beta,iType,NewResFlag)
				end if
			case (7)!Tensor*number,(3) *(1) or (2,3) *(1)  or (2,3,1) * (1)
				if(.not.Res%getFlag())then
					!if(A%dim(rankA).eq.1)then
					!	call Res%allocate(A%Dimension.subdim.[1,rankA-contractLegLen],iType)
					!else
						call Res%allocate(A%Dimension,iType)
					!end if
				end if
				call ProductTensorCase78(Res,A,B,alpha,beta,iType,NewResFlag)
			case (8)!Tensor*Tensor,(1,2,2,1) *(1,1)
				if(rankA.eq.1)then!Tensor*number,(3)*(1,1)
					call writemess("ERROR in ProductTensor,case -3,stop",-1)
					call error_stop()
				else!Tensor*Tensor,(1,2,2,1) *(1,1)
					if(A%dim(rankA).ne.1)then
						call writemess("ERROR in ProductTensor,case -4,stop",-1)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop()
					end if
					if(.not.Res%getFlag())then
						call pasteTwoSubDim(NewDim,A%Dimension,1,rankA-contractLegLen,&
									B%Dimension,contractLegLen+1,rankB)
						call Res%allocate(NewDim,iType)
					end if
					call ProductTensorCase78(Res,A,B,alpha,beta,iType,NewResFlag)
				end if
			case(9)! vec * vec
				if(.not.Res%getFlag())then
					call Res%allocate([1],iType)
				end if
				if(A%getSymmetryFlag())then
					call ProductTensorBlockCase9(Res,A,B,alpha,beta,iType,NewResFlag)
				else
					if(TotalBlockA.ne.TotalBlockB)then
						call writemess('ERROR in ProductTensorRoutine, TotalData',-1)
						call writemess('TotalBlockA='+TotalBlockA)
						call writemess('TotalBlockB='+TotalBlockB)
						call error_stop
					end if
					call ProductTensorCase9(Res,A,B,alpha,beta,iType,NewResFlag)
				end if
			case(10)! vec * matrix [M] [M,N]
				if(A%getSymmetryFlag())then
					if((.not.A%Data%getFlag()).or.(.not.B%Data%getFlag()))then
						Res%Dimension=B%Dimension.subdim.[contractLegLen+1,rankB]
						call Res%Data%allocate([0],iType)
						return
					end if
					call WorkingMemory%check()
					call WorkingMemory%allocate(1,2+TotalBlockB+TotalBlockB+rankB)
					call WorkingMemory%get_memory(dimB,2)
					call WorkingMemory%get_memory(BlockDimB1,TotalBlockB)
					call WorkingMemory%get_memory(BlockDimB2,TotalBlockB)
					call WorkingMemory%get_memory(index,rankB)
					call MatrixDimBlock(B,contractLegLen,dimB,BlockDimB1,BlockDimB2,index)

					TotalBlockR=dimB(2)
					call allocateCheck(product_working_ResBlockNum,TotalBlockR)
					ResBlockNum=>product_working_ResBlockNum(1:TotalBlockR)
					call MatrixDimBlockBToRes(dimB,BlockDimB1,BlockDimB2,A,ResBlockNum)

					if(.not.Res%getFlag())then
						Res%Dimension=B%Dimension.subdim.[contractLegLen+1,rankB]
						call Res%Data%allocate(ResBlockNum,iType)
					end if
					if(.not.Res%Data%getFlag())then
						call WorkingMemory%free()
						return
					end if

					call ProductTensorBlockCase10(Res,A,B,alpha,beta,DimB,BlockDimB1,BlockDimB2,iType,NewResFlag)
					call WorkingMemory%free()
				else
					if(.not.Res%getFlag())then
						NewDim=B%Dimension.subdim.[contractLegLen+1,rankB]
						call Res%allocate(NewDim,iType)
					end if
					call WorkingMemory%check()
					call WorkingMemory%allocate(1,2)
					call WorkingMemory%get_memory(dimB,2)
					call MatrixDim(B,contractLegLen,dimB)
					call ProductTensorCase10(Res,A,B,alpha,beta,dimB(1),dimB(2),iType,NewResFlag)
					call WorkingMemory%free()
				end if
			case(11)!  matrix* vec  [M,N] [N]
				if(A%getSymmetryFlag())then
					if((.not.A%Data%getFlag()).or.(.not.B%Data%getFlag()))then
						Res%Dimension=A%Dimension.subdim.[1,rankA-contractLegLen]
						call Res%Data%allocate([0],iType)
						return
					end if
					call WorkingMemory%check()
					call WorkingMemory%allocate(1,2+TotalBlockA+TotalBlockA+rankA)
					call WorkingMemory%get_memory(dimA,2)
					call WorkingMemory%get_memory(BlockDimA1,TotalBlockA)
					call WorkingMemory%get_memory(BlockDimA2,TotalBlockA)
					call WorkingMemory%get_memory(index,rankA)
					call MatrixDimBlock(A,rankA-contractLegLen,dimA,BlockDimA1,BlockDimA2,index)


					TotalBlockR=DimA(1)
					call allocateCheck(product_working_ResBlockNum,TotalBlockR)
					ResBlockNum=>product_working_ResBlockNum(1:TotalBlockR)
					call MatrixDimBlockAToRes(DimA,BlockDimA1,BlockDimA2,B,ResBlockNum)

					if(.not.Res%getFlag())then
						Res%Dimension=A%Dimension.subdim.[1,rankA-contractLegLen]
						call Res%Data%allocate(ResBlockNum,iType)
					end if
					
					if(.not.Res%Data%getFlag())then
						call WorkingMemory%free()
						return
					end if

					call ProductTensorBlockCase11(Res,A,B,alpha,beta,DimA,BlockDimA1,BlockDimA2,iType,NewResFlag)
					call WorkingMemory%free()
				else
					if(.not.Res%getFlag())then
						NewDim=A%Dimension.subdim.[1,rankA-contractLegLen]
						call Res%allocate(NewDim,iType)
					end if
					call WorkingMemory%check()
					call WorkingMemory%allocate(1,2)
					call WorkingMemory%get_memory(dimA,2)
					call MatrixDim(A,rankA-contractLegLen,dimA)
					call ProductTensorCase11(Res,A,B,alpha,beta,dimA(1),dimA(2),iType,NewResFlag)
					call WorkingMemory%free()
				end if
			case(12)!  matrix* matrix  [M,K] [K,N]
				if(A%getSymmetryFlag())then
					if((.not.A%Data%getFlag()).or.(.not.B%Data%getFlag()))then
						call pasteTwoSubDim(Res%Dimension,A%Dimension,1,rankA-contractLegLen,&
									B%Dimension,contractLegLen+1,rankB)
						call Res%Data%allocate([0],iType)
						return
					end if
					rankR=rankA+rankB-contractLegLen-contractLegLen
					TotalMomeryLength=6+TotalBlockA+TotalBlockA+rankA+TotalBlockB+TotalBlockB+rankB
					TotalMomeryLength=TotalMomeryLength

					call WorkingMemory%check()
					call WorkingMemory%allocate(1,TotalMomeryLength)
					call WorkingMemory%get_memory(dimA,2)
					call WorkingMemory%get_memory(BlockDimA1,TotalBlockA)
					call WorkingMemory%get_memory(BlockDimA2,TotalBlockA)
					call WorkingMemory%get_memory(index,rankA)
					call MatrixDimBlock(A,rankA-contractLegLen,dimA,BlockDimA1,BlockDimA2,index)

					call WorkingMemory%get_memory(dimB,2)
					call WorkingMemory%get_memory(BlockDimB1,TotalBlockB)
					call WorkingMemory%get_memory(BlockDimB2,TotalBlockB)
					call WorkingMemory%get_memory(index,rankB)
					call MatrixDimBlock(B,contractLegLen,dimB,BlockDimB1,BlockDimB2,index)

					call WorkingMemory%get_memory(dimR,2)

					TotalBlockR=DimA(1)*DimB(2)

					call allocateCheck(product_working_BlockM,TotalBlockR)
					BlockDimR1=>product_working_BlockM(1:TotalBlockR)

					call allocateCheck(product_working_BlockN,TotalBlockR)
					BlockDimR2=>product_working_BlockN(1:TotalBlockR)

					call allocateCheck(product_working_ResBlockNum,TotalBlockR)
					ResBlockNum=>product_working_ResBlockNum(1:TotalBlockR)



					call MatrixDimBlockABToRes(dimA,BlockDimA1,BlockDimA2,&
											   dimB,BlockDimB1,BlockDimB2,&
											   dimR,BlockDimR1,BlockDimR2,ResBlockNum)
					
					if(.not.Res%getFlag())then
						call pasteTwoSubDim(Res%Dimension,A%Dimension,1,rankA-contractLegLen,&
									B%Dimension,contractLegLen+1,rankB)
						call Res%Data%allocate(ResBlockNum,iType)
					end if
					
					if(.not.Res%Data%getFlag())then
						call WorkingMemory%free()
						return
					end if

					call ProductTensorBlockCase12(Res,A,B,alpha,beta,DimA,BlockDimA1,&
							BlockDimA2,DimB,BlockDimB1,BlockDimB2,DimR,BlockDimR1,&
							BlockDimR2,iType,NewResFlag)
					call WorkingMemory%free()
				else
					if(.not.Res%getFlag())then
						call pasteTwoSubDim(NewDim,A%Dimension,1,rankA-contractLegLen,&
									B%Dimension,contractLegLen+1,rankB)
						call Res%allocate(NewDim,iType)
					end if
					call WorkingMemory%check()
					call WorkingMemory%allocate(1,4)
					call WorkingMemory%get_memory(dimA,2)
					call WorkingMemory%get_memory(dimB,2)
					call MatrixDim(A,rankA-contractLegLen,dimA)
					call MatrixDim(B,contractLegLen,dimB)
					call ProductTensorCase12(Res,A,B,alpha,beta,dimA(1),dimA(2),dimB(2),iType,NewResFlag)
					call WorkingMemory%free()
				end if
		end select
		return
	end subroutine

	!MatrixDim: 1 to Legith will fuse as one legs and the other as the other one

	subroutine MatrixDimBlock(A,Legith,Dim,BlockDim1,BlockDim2,index)
		type(Tensor),intent(in)::A
		integer,intent(in)::Legith
		integer,intent(inout)::dim(:),BlockDim1(:),BlockDim2(:),index(:)
		integer,pointer::Deg(:),maxindex(:)
		integer::bi,i,rank
		logical::goon
		call MatrixDim(A,Legith,Dim)
		index=1
		rank=A%getRank()
		goon=.true.
		call A%pointDim(maxindex)
		do bi=1,A%getTotalBlock()
			if(.not.goon)then
				call writemess('ERROR in MatrixDimBlock',-1)
				call error_stop
			end if
			if(A%getFlag(bi))then
				if(Legith.eq.0)then
					BlockDim1(bi)=0
				else
					call A%pointDeg(Deg,1)
					BlockDim1(bi)=Deg(index(1))
					do i=2,Legith
						call A%pointDeg(Deg,i)
						BlockDim1(bi)=BlockDim1(bi)*Deg(index(i))
					end do
				end if
				if(Legith.eq.rank)then
					BlockDim2(bi)=0
				else
					call A%pointDeg(Deg,Legith+1)
					BlockDim2(bi)=Deg(index(Legith+1))
					do i=Legith+2,rank
						call A%pointDeg(Deg,i)
						BlockDim2(bi)=BlockDim2(bi)*Deg(index(i))
					end do
				end if
			else
				BlockDim1(bi)=0
				BlockDim2(bi)=0
			end if
			goon=index_counter(index,maxindex)
		end do
		return
	end subroutine

	subroutine MatrixDim(A,Legith,Dim)
		type(Tensor),intent(in)::A
		integer,intent(in)::Legith
		integer,intent(inout)::dim(:)
		integer::bi,i,rank
		integer,pointer::dimData(:)
		rank=A%getRank()
		if(Legith.lt.0)then
			call writemess('ERROR ',-1)
			call error_stop
		end if
		call A%pointDim(dimData)
		if(Legith.eq.0)then
			dim(1)=0
		else
			!dim(1)=dimData(1)
			!do i=2,Legith
			!	dim(1)=dim(1)*dimData(i)
			!end do
			dim(1)=product(dimData(1:Legith))
		end if
		if(Legith.eq.rank)then
			dim(2)=0
		else
			dim(2)=product(dimData((Legith+1):rank))
			!dim(2)=dimData(Legith+1)
			!do i=Legith+2,rank
			!	dim(2)=dim(2)*dimData(i)
			!end do
		end if
	end subroutine

	subroutine MatrixDimBlockABToRes(DimA,MA,NA,DimB,MB,NB,DimR,MRes,NRes,BlockNum)
		integer,intent(in)::DimA(2),DimB(2)
		integer,intent(in)::MA(DimA(1),DimA(2))
		integer,intent(in)::NA(DimA(1),DimA(2))
		integer,intent(in)::MB(DimB(1),DimB(2))
		integer,intent(in)::NB(DimB(1),DimB(2))
		integer,intent(inout)::DimR(2)
		integer,intent(inout)::MRes(DimA(1),DimB(2))
		integer,intent(inout)::NRes(DimA(1),DimB(2))
		integer,intent(inout)::BlockNum(DimA(1),DimB(2))
		logical::Flag
		integer::i,j,k,M,N
		DimR=[DimA(1),DimB(2)]
		do j=1,DimR(2)
			do i=1,DimR(1)
				Flag=.false.
				do k=1,DimA(2)
					if(.not.((MA(i,k).eq.0).or.(NA(i,k).eq.0).or.(MB(k,j).eq.0).or.(NB(k,j).eq.0)))then
						M=MA(i,k)
						N=NB(k,j)
						Flag=.true.
						exit
					end if
				end do
				if(Flag)then
					MRes(i,j)=M
					NRes(i,j)=N
					BlockNum(i,j)=M*N
				else
					MRes(i,j)=0
					NRes(i,j)=0
					BlockNum(i,j)=0
				end if
			end do
		end do
		return
	end subroutine


	subroutine MatrixDimBlockAToRes(DimA,MA,NA,B,BlockNum)
		type(Tensor),intent(in)::B
		integer,intent(in)::DimA(2)
		integer,intent(in)::MA(DimA(1),DimA(2))
		integer,intent(in)::NA(DimA(1),DimA(2))
		integer,intent(inout)::BlockNum(DimA(1))
		logical::Flag
		integer::i,k,M
		do i=1,DimA(1)
			Flag=.false.
			do k=1,DimA(2)
				if(.not.((MA(i,k).eq.0).or.(NA(i,k).eq.0).or.(.not.B%getFlag(k))))then
					M=MA(i,k)
					Flag=.true.
					exit
				end if
			end do
			if(Flag)then
				BlockNum(i)=M
			else
				BlockNum(i)=0
			end if
		end do
		return
	end subroutine

	subroutine MatrixDimBlockBToRes(DimB,MB,NB,A,BlockNum)
		type(Tensor),intent(in)::A
		integer,intent(in)::DimB(2)
		integer,intent(in)::MB(DimB(1),DimB(2))
		integer,intent(in)::NB(DimB(1),DimB(2))
		integer,intent(inout)::BlockNum(DimB(2))
		logical::Flag
		integer::j,k,N
		do j=1,DimB(2)
			Flag=.false.
			do k=1,DimB(1)
				if(.not.((.not.A%getFlag(k)).or.(MB(k,j).eq.0).or.(NB(k,j).eq.0)))then
					N=NB(k,j)
					Flag=.true.
					exit
				end if
			end do
			if(Flag)then
				BlockNum(j)=N
			else
				BlockNum(j)=0
			end if
		end do
		return
	end subroutine


	function ProductTensor(A,B)result(C)
		type(Tensor)::C
		type(Tensor),intent(in)::A,B
		call ProductTensorRoutine(C,A,B,1,0)
		return
	end function

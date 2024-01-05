	function SelectCaseForSym(rank1,total1,rank2,total2)result(Flag)
		integer::Flag
		integer,intent(in)::rank1,total1,rank2,total2
		if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=9
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=10
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=11	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=12
		else
			call writemess('ERROR in ProductTensor',-1)
			call writemess('A%getRank()='+rank1,-1)
			call writemess('B%getRank()='+rank2,-1)
			call error_stop()
		end if
	end function

	subroutine ProductTensorBlockCase9(Res,A,B,alpha,beta,classType,NewResFlag)!vec*vec
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:,:),iAp(:,:),iBp(:,:),iRp0(:),iAp0(:),iBp0(:)
		real*4,pointer::sRp(:),sAp(:),sBp(:)
		real*8,pointer::dRp(:),dAp(:),dBp(:)
		complex*8,pointer::cRp(:),cAp(:),cBp(:)
		complex*16,pointer::zRp(:),zAp(:),zBp(:)
		real*4,External::	sdot
		real*8,External::	ddot
		complex(kind=4),External::	cdotu
		complex(kind=8),External::	zdotu
		real*4::salpha,sbeta
		real*4::dalpha,dbeta
		complex*8::calpha,cbeta
		complex*16::zalpha,zbeta
		integer::ialpha,ibeta
		integer::i
		select case(classType)
			case(1)
				ialpha=iselect(alpha)
				ibeta=iselect(beta)
				call Res%pointer(iRp0)
				iRp(1:1,1:1)=>iRp0(1:1)
				call A%pointer(iAp0,1)
				call B%pointer(iBp0,1)
				if(associated(iAp0).and.associated(iBp0))then
					iAp(1:1,1:size(iAp0))=>iAp0(1:size(iAp0))
					iBp(1:size(iBp0),1:1)=>iBp0(1:size(iBp0))
					if(NewResFlag)then
						iRp=ialpha*matmul(iAp,iBp)
					else
						iRp=ialpha*matmul(iAp,iBp)+ibeta*iRp
					end if
				else
					if(NewResFlag)then
						iRp=0
					else
						iRp= ibeta*iRp
					end if
				end if
				do i=2,A%getTotalBlock()
					call A%pointer(iAp0,i)
					call B%pointer(iBp0,i)
					if(associated(iAp0).and.associated(iBp0))then
						iAp(1:1,1:size(iAp0))=>iAp0(1:size(iAp0))
						iBp(1:size(iBp0),1:1)=>iBp0(1:size(iBp0))
						iRp=ialpha*matmul(iAp,iBp)+iRp
					end if
				end do
				
			case(2)
				salpha=sselect(alpha)
				sbeta=sselect(beta)
				call Res%pointer(sRp)
				call A%pointer(sAp,1)
				call B%pointer(sBp,1)
				if(associated(sAp).and.associated(sBp))then
					if(NewResFlag)then
						sRp(1)=( salpha * sdot(size(sAp), sAp, 1, sBp, 1) ) 
					else
						sRp(1)=( salpha * sdot(size(sAp), sAp, 1, sBp, 1) ) + ( sbeta*sRp(1) )
					end if
				else
					if(NewResFlag)then
						sRp(1)=0.
					else
						sRp(1)= sbeta*sRp(1) 
					end if
				end if
				do i=2,A%getTotalBlock()
					call A%pointer(sAp,i)
					call B%pointer(sBp,i)
					if(associated(sAp).and.associated(sBp))then
						sRp(1)=( salpha * sdot(size(sAp), sAp, 1, sBp, 1) ) + sRp(1) 
					end if
				end do
			case(3)
				dalpha=dselect(alpha)
				dbeta=dselect(beta)
				call Res%pointer(dRp)
				call A%pointer(dAp,1)
				call B%pointer(dBp,1)
				if(associated(dAp).and.associated(dBp))then
					if(NewResFlag)then
						dRp(1)=( dalpha * ddot(size(dAp), dAp, 1, dBp, 1) )
					else
						dRp(1)=( dalpha * ddot(size(dAp), dAp, 1, dBp, 1) ) + ( dbeta*dRp(1) )
					end if
				else
					if(NewResFlag)then
						dRp(1)=0d0
					else
						dRp(1)=dbeta*dRp(1) 
					end if
				end if
				do i=2,A%getTotalBlock()
					call A%pointer(dAp,i)
					call B%pointer(dBp,i)
					if(associated(dAp).and.associated(dBp))then
						dRp(1)=( dalpha * ddot(size(dAp), dAp, 1, dBp, 1) ) + dRp(1) 
					end if
				end do
			case(4)
				calpha=cselect(alpha)
				cbeta=cselect(beta)
				call Res%pointer(cRp)
				call A%pointer(cAp,1)
				call B%pointer(cBp,1)
				if(associated(cAp).and.associated(cBp))then
					if(NewResFlag)then
						cRp(1)=( calpha * cdotu(size(cAp), cAp, 1, cBp, 1) )
					else
						cRp(1)=( calpha * cdotu(size(cAp), cAp, 1, cBp, 1) ) + ( cbeta*cRp(1) )
					end if
				else
					if(NewResFlag)then
						cRp(1)=0.
					else
						cRp(1)=cbeta*cRp(1) 
					end if
				end if
				do i=2,A%getTotalBlock()
					call A%pointer(cAp,i)
					call B%pointer(cBp,i)
					if(associated(cAp).and.associated(cBp))then
						cRp(1)=( calpha * cdotu(size(cAp), cAp, 1, cBp, 1) ) + cRp(1) 
					end if
				end do
			case(5)
				zalpha=zselect(alpha)
				zbeta=zselect(beta)
				call Res%pointer(zRp)
				call A%pointer(zAp,1)
				call B%pointer(zBp,1)
				if(associated(zAp).and.associated(zBp))then
					if(NewResFlag)then
						zRp(1)=( zalpha * zdotu(size(zAp), zAp, 1, zBp, 1) )
					else
						zRp(1)=( zalpha * zdotu(size(zAp), zAp, 1, zBp, 1) ) + ( zbeta*zRp(1) )
					end if
				else
					if(NewResFlag)then
						zRp(1)=0d0
					else
						zRp(1)=zbeta*zRp(1)
					end if
				end if
				do i=2,A%getTotalBlock()
					call A%pointer(zAp,i)
					call B%pointer(zBp,i)
					if(associated(zAp).and.associated(zBp))then
						zRp(1)=( zalpha * zdotu(size(zAp), zAp, 1, zBp, 1) ) + zRp(1) 
					end if
				end do
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorBlockCase10(Res,A,B,alpha,beta,DimB,BlockDimB1,BlockDimB2,classType,NewResFlag)!vec*matrix
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,DimB(2)
		integer,intent(in)::BlockDimB1(DimB(1),*)
		integer,intent(in)::BlockDimB2(DimB(1),*)
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:,:),iAp(:,:),iBp(:,:),iRp0(:),iAp0(:)
		real*4,pointer::sRp(:),sAp(:),sBp(:,:)
		real*8,pointer::dRp(:),dAp(:),dBp(:,:)
		complex*8,pointer::cRp(:),cAp(:),cBp(:,:)
		complex*16,pointer::zRp(:),zAp(:),zBp(:,:)
		integer::i,j
		integer::inum,ialpha
		real*4::snum,salpha
		real*8::dnum,dalpha
		complex*8::cnum,calpha
		complex*16::znum,zalpha
		logical::Flag
		Flag=NewResFlag
		select case(classType)
			case(1)
				ialpha=iselect(alpha)
				do j=1,DimB(2)
					call Res%pointer(iRp0,j)
					if(associated(iRp0))then
						inum=iselect(beta)
						do i=1,DimB(1)
							call A%pointer(iAp0,i)
							call B%pointer(DimB,iBp,[BlockDimB1(i,j),BlockDimB2(i,j)],[i,j])
							if(associated(iAp0).and.associated(iBp))then
								iRp(1:1,1:BlockDimB2(i,j))=>iRp0
								iAp(1:1,1:BlockDimB1(i,j))=>iAp0
								if(Flag)then
									iRp=ialpha*matmul(iAp,iBp)
									Flag=.false.
								else
									iRp=ialpha*matmul(iAp,iBp)+inum*iRp
								end if
								inum=1
							end if
						end do
					end if
				end do
			case(2)
				salpha=sselect(alpha)
				do j=1,DimB(2)
					call Res%pointer(sRp,j)
					if(associated(sRp))then
						snum=sselect(beta)
						do i=1,DimB(1)
							call A%pointer(sAp,i)
							call B%pointer(DimB,sBp,[BlockDimB1(i,j),BlockDimB2(i,j)],[i,j])
							if(associated(sAp).and.associated(sBp))then
								call SGEMV('T',BlockDimB1(i,j),BlockDimB2(i,j),salpha,sBp,BlockDimB1(i,j),&
											 sAp, 1, snum, sRp, 1)
								snum=1.
							end if
						end do
					end if
				end do
			case(3)
				dalpha=dselect(alpha)
				do j=1,DimB(2)
					call Res%pointer(dRp,j)
					if(associated(dRp))then
						dnum=dselect(beta)
						do i=1,DimB(1)
							call A%pointer(dAp,i)
							call B%pointer(DimB,dBp,[BlockDimB1(i,j),BlockDimB2(i,j)],[i,j])
							if(associated(dAp).and.associated(dBp))then
								call DGEMV('T',BlockDimB1(i,j),BlockDimB2(i,j),dalpha,dBp,BlockDimB1(i,j),&
											 dAp, 1, dnum, dRp, 1)
								dnum=1d0
							end if
						end do
					end if
				end do
			case(4)
				calpha=cselect(alpha)
				do j=1,DimB(2)
					call Res%pointer(cRp,j)
					if(associated(cRp))then
						cnum=cselect(beta)
						do i=1,DimB(1)
							call A%pointer(cAp,i)
							call B%pointer(DimB,cBp,[BlockDimB1(i,j),BlockDimB2(i,j)],[i,j])
							if(associated(cAp).and.associated(cBp))then
								call CGEMV('T',BlockDimB1(i,j),BlockDimB2(i,j),calpha,cBp,BlockDimB1(i,j),&
											 cAp, 1, cnum, cRp, 1)
								cnum=cmplx(1.,0.)
							end if
						end do
					end if
				end do
			case(5)
				zalpha=zselect(alpha)
				do j=1,DimB(2)
					call Res%pointer(zRp,j)
					if(associated(zRp))then
						znum=zselect(beta)
						do i=1,DimB(1)
							call A%pointer(zAp,i)
							call B%pointer(DimB,zBp,[BlockDimB1(i,j),BlockDimB2(i,j)],[i,j])
							if(associated(zAp).and.associated(zBp))then
								call ZGEMV('T',BlockDimB1(i,j),BlockDimB2(i,j),zalpha,zBp,BlockDimB1(i,j),&
											 zAp, 1, znum, zRp, 1)
								znum=dcmplx(1d0,0d0)
							end if
						end do
					end if
				end do
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine


	subroutine ProductTensorBlockCase11(Res,A,B,alpha,beta,DimA,BlockDimA1,BlockDimA2,classType,NewResFlag)!matrix*vec
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,DimA(:)
		integer,intent(in)::BlockDimA1(DimA(1),*),BlockDimA2(DimA(1),*)
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:,:),iAp(:,:),iBp(:,:),iRp0(:),iBp0(:)
		real*4,pointer::sRp(:),sAp(:,:),sBp(:)
		real*8,pointer::dRp(:),dAp(:,:),dBp(:)
		complex*8,pointer::cRp(:),cAp(:,:),cBp(:)
		complex*16,pointer::zRp(:),zAp(:,:),zBp(:)
		integer::i,j
		integer::inum,ialpha
		real*4::snum,salpha
		real*8::dnum,dalpha
		complex*8::cnum,calpha
		complex*16::znum,zalpha
		logical::Flag
		Flag=NewResFlag
		select case(classType)
			case(1)
				ialpha=iselect(alpha)
				do i=1,DimA(1)
					call Res%pointer(iRp0,i)
					if(associated(iRp0))then
						inum=iselect(beta)
						do j=1,DimA(2)
							call A%pointer(DimA,iAp,[BlockDimA1(i,j),BlockDimA2(i,j)],[i,j])
							call B%pointer(iBp0,j)
							if(associated(iAp).and.associated(iBp0))then
								iRp(1:BlockDimA1(i,j),1:1)=>iRp0
								iBp(1:BlockDimA2(i,j),1:1)=>iBp0
								if(Flag)then
									iRp=ialpha*matmul(iAp,iBp)
									Flag=.false.
								else
									iRp=ialpha*matmul(iAp,iBp)+inum*iRp
								end if
								inum=1
							end if
						end do
					end if
				end do
			case(2)
				salpha=sselect(alpha)
				do i=1,DimA(1)
					call Res%pointer(sRp,i)
					if(associated(sRp))then
						snum=sselect(beta)
						do j=1,DimA(2)
							call A%pointer(DimA,sAp,[BlockDimA1(i,j),BlockDimA2(i,j)],[i,j])
							call B%pointer(sBp,j)
							if(associated(sAp).and.associated(sBp))then
								call SGEMV('N',BlockDimA1(i,j),BlockDimA2(i,j),salpha,sAp,BlockDimA1(i,j),&
											 sBp, 1, snum, sRp, 1)
								snum=1.
							end if
						end do
					end if
				end do
			case(3)
				dalpha=dselect(alpha)
				do i=1,DimA(1)
					call Res%pointer(dRp,i)
					if(associated(dRp))then
						dnum=dselect(beta)
						do j=1,DimA(2)
							call A%pointer(DimA,dAp,[BlockDimA1(i,j),BlockDimA2(i,j)],[i,j])
							call B%pointer(dBp,j)
							if(associated(dAp).and.associated(dBp))then
								call DGEMV('N',BlockDimA1(i,j),BlockDimA2(i,j),dalpha,dAp,BlockDimA1(i,j),&
											 dBp, 1, dnum, dRp, 1)
								dnum=1d0
							end if
						end do
					end if
				end do
			case(4)
				calpha=cselect(alpha)
				do i=1,DimA(1)
					call Res%pointer(cRp,i)
					if(associated(cRp))then
						cnum=cselect(beta)
						do j=1,DimA(2)
							call A%pointer(DimA,cAp,[BlockDimA1(i,j),BlockDimA2(i,j)],[i,j])
							call B%pointer(cBp,j)
							if(associated(cAp).and.associated(cBp))then
								call CGEMV('N',BlockDimA1(i,j),BlockDimA2(i,j),calpha,cAp,BlockDimA1(i,j),&
											 cBp, 1, cnum, cRp, 1)
								cnum=cmplx(1.,0.)
							end if
						end do
					end if
				end do
			case(5)
				zalpha=zselect(alpha)
				do i=1,DimA(1)
					call Res%pointer(zRp,i)
					if(associated(zRp))then
						znum=zselect(beta)
						do j=1,DimA(2)
							call A%pointer(DimA,zAp,[BlockDimA1(i,j),BlockDimA2(i,j)],[i,j])
							call B%pointer(zBp,j)
							if(associated(zAp).and.associated(zBp))then
								call ZGEMV('N',BlockDimA1(i,j),BlockDimA2(i,j),zalpha,zAp,BlockDimA1(i,j),&
											 zBp, 1, znum, zRp, 1)
								znum=dcmplx(1d0,0d0)
							end if
						end do
					end if
				end do
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorBlockCase12(Res,A,B,alpha,beta,DimA,BlockDimA1,&
		BlockDimA2,DimB,BlockDimB1,BlockDimB2,DimR,BlockDimR1&
		,BlockDimR2,classType,NewResFlag)!matrix*vec
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,DimA(:),DimB(:),DimR(:)
		integer,intent(in)::BlockDimA1(DimA(1),*),BlockDimA2(DimA(1),*)
		integer,intent(in)::BlockDimB1(DimB(1),*),BlockDimB2(DimB(1),*)
		integer,intent(in)::BlockDimR1(DimR(1),*),BlockDimR2(DimR(1),*)
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:,:),iAp(:,:),iBp(:,:)
		real*4,pointer::sRp(:,:),sAp(:,:),sBp(:,:)
		real*8,pointer::dRp(:,:),dAp(:,:),dBp(:,:)
		complex*8,pointer::cRp(:,:),cAp(:,:),cBp(:,:)
		complex*16,pointer::zRp(:,:),zAp(:,:),zBp(:,:)
		integer::i,j,k,LD1,LD2,LD3
		logical::First
		integer::inum,ialpha
		real*4::snum,salpha
		real*8::dnum,dalpha
		complex*8::cnum,calpha
		complex*16::znum,zalpha
		logical::Flag
		Flag=NewResFlag
		select case(classType)
			case(1)
				ialpha=iselect(alpha)
				do j=1,DimR(2)
					do i=1,DimR(1)
						call Res%pointer(DimR,iRp,[BlockDimR1(i,j),BlockDimR2(i,j)],[i,j])
						if(associated(iRp))then
							inum=iselect(beta)
							do k=1,DimA(2)
								call A%pointer(DimA,iAp,[BlockDimA1(i,k),BlockDimA2(i,k)],[i,k])
								call B%pointer(DimB,iBp,[BlockDimB1(k,j),BlockDimB2(k,j)],[k,j])
								if(associated(iAp).and.associated(iBp))then
									if(Flag)then
										iRp=ialpha*matmul(iAp,iBp)
										Flag=.false.
									else
										iRp=ialpha*matmul(iAp,iBp)+inum*iRp
									end if
									inum=1
								end if
							end do
						end if
					end do
				end do
			case(2)
				salpha=sselect(alpha)
				do j=1,DimR(2)
					do i=1,DimR(1)
						call Res%pointer(DimR,sRp,[BlockDimR1(i,j),BlockDimR2(i,j)],[i,j])
						if(associated(sRp))then
							snum=sselect(beta)
							do k=1,DimA(2)
								call A%pointer(DimA,sAp,[BlockDimA1(i,k),BlockDimA2(i,k)],[i,k])
								call B%pointer(DimB,sBp,[BlockDimB1(k,j),BlockDimB2(k,j)],[k,j])
								if(associated(sAp).and.associated(sBp))then
									call SGEMM('N','N',BlockDimA1(i,k),BlockDimB2(k,j),BlockDimB1(k,j),&
										salpha,sAp,BlockDimA1(i,k),sBp,BlockDimB1(k,j),snum,sRp,&
										BlockDimA1(i,k))
									snum=1.
								end if
							end do
						end if
					end do
				end do
			case(3)
				dalpha=dselect(alpha)
				do j=1,DimR(2)
					do i=1,DimR(1)
						call Res%pointer(DimR,dRp,[BlockDimR1(i,j),BlockDimR2(i,j)],[i,j])
						if(associated(dRp))then
							dnum=dselect(beta)
							do k=1,DimA(2)
								call A%pointer(DimA,dAp,[BlockDimA1(i,k),BlockDimA2(i,k)],[i,k])
								call B%pointer(DimB,dBp,[BlockDimB1(k,j),BlockDimB2(k,j)],[k,j])
								if(associated(dAp).and.associated(dBp))then
									call DGEMM('N','N',BlockDimA1(i,k),BlockDimB2(k,j),BlockDimB1(k,j),&
										dalpha,dAp,BlockDimA1(i,k),dBp,BlockDimB1(k,j),dnum,dRp,&
										BlockDimA1(i,k))
									dnum=1d0
								end if
							end do
						end if
					end do
				end do
			case(4)
				calpha=cselect(alpha)
				do j=1,DimR(2)
					do i=1,DimR(1)
						call Res%pointer(DimR,cRp,[BlockDimR1(i,j),BlockDimR2(i,j)],[i,j])
						if(associated(cRp))then
							cnum=cselect(beta)
							do k=1,DimA(2)
								call A%pointer(DimA,cAp,[BlockDimA1(i,k),BlockDimA2(i,k)],[i,k])
								call B%pointer(DimB,cBp,[BlockDimB1(k,j),BlockDimB2(k,j)],[k,j])
								if(associated(cAp).and.associated(cBp))then
									call CGEMM('N','N',BlockDimA1(i,k),BlockDimB2(k,j),BlockDimB1(k,j),&
										calpha,cAp,BlockDimA1(i,k),cBp,BlockDimB1(k,j),cnum,cRp,&
										BlockDimA1(i,k))
									cnum=cmplx(1.,0.)
								end if
							end do
						end if
					end do
				end do
			case(5)
				zalpha=zselect(alpha)
				do j=1,DimR(2)
					do i=1,DimR(1)
						call Res%pointer(DimR,zRp,[BlockDimR1(i,j),BlockDimR2(i,j)],[i,j])
						if(associated(zRp))then
							znum=zselect(beta)
							do k=1,DimA(2)
								call A%pointer(DimA,zAp,[BlockDimA1(i,k),BlockDimA2(i,k)],[i,k])
								call B%pointer(DimB,zBp,[BlockDimB1(k,j),BlockDimB2(k,j)],[k,j])
								if(associated(zAp).and.associated(zBp))then
									call ZGEMM('N','N',BlockDimA1(i,k),BlockDimB2(k,j),BlockDimB1(k,j),&
										zalpha,zAp,BlockDimA1(i,k),zBp,BlockDimB1(k,j),znum,zRp,&
										BlockDimA1(i,k))
									znum=dcmplx(1d0,0d0)
								end if
							end do
						end if
					end do
				end do
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine
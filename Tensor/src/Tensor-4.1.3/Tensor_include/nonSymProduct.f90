	function SelectCase(rank1,total1,rank2,total2)result(Flag)
		integer::Flag
		integer,intent(in)::rank1,total1,rank2,total2
		if((total1.eq.1).or.(total2.eq.1)) then
			if((total1.eq.1).and.(total2.eq.1)) then!there is only one element in both T1 and T2
				if((rank1.eq.1).and.(rank2.eq.1))then!Number*number,(1) *(1)
					flag=1
				else if((rank1.eq.1).and.(rank2.ne.1))then!Number*Tensor,(1) * (1,1,1)
					flag=2
				else if((rank1.ne.1).and.(rank2.eq.1))then!Tensor*Number, (1,1,1) * (1)
					flag=3
				else if((rank1.ne.1).and.(rank2.ne.1))then!Tensor*Tensor, (1,1,1) * (1,1)
					flag=4
				else 
					call writemess("ERROR in ProductTensor,case -1,stop",-1)
					call error_stop()
				end if
			else if((total1.eq.1).and.(total2.ne.1)) then!there is only one element in both T1, but not T2
					if(rank1.eq.1)then!Number*Tensor,(1) *(3) or (1) *(2,3)  or (1) * (1,2,2)
						flag=5
					else 
						if(rank2.eq.1)then!Tensor*number,(1,1) *(3)
							call writemess("ERROR in ProductTensor,case -1,stop",-1)
							call error_stop()
						else!Tensor*Tensor,(1,1) *(1,2,1,2)
							!if(T2%dim(1).ne.1)then
							!	call writemess("ERROR in ProductTensor,case -2,stop",-1)
							!	call T1%diminfo('dimension of T1')
							!	call T2%diminfo('dimension of T2')
							!	call error_stop()
							!end if
							flag=6
						end if
					end if
			else if((total1.ne.1).and.(total2.eq.1)) then!there is only one element in  T2, but not T1
					if(rank2.eq.1)then!Tensor*number,(3) *(1) or (2,3) *(1)  or (2,3,1) * (1)
						flag=7
					else
						if(rank1.eq.1)then!Tensor*number,(3)*(1,1)
							call writemess("ERROR in ProductTensor,case -3,stop",-1)
							call error_stop()
						else!Tensor*Tensor,(1,2,2,1) *(1,1)
							!if(T1%dim(rank1).ne.1)then
							!	call writemess("ERROR in ProductTensor,case -4,stop",-1)
							!	call T1%diminfo('dimension of T1')
							!	call T2%diminfo('dimension of T2')
							!	call error_stop()
							!end if
							flag=8
						end if
					end if
			else
					call writemess("ERROR in ProductTensor,case -5,stop",-1)
					call error_stop()
			end if
		else if((rank1.eq.1).and.(rank2.eq.1)) then
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
	subroutine ProductTensorCase1234(Res,A,B,alpha,beta,classType,NewResFlag)!Number*number,(1) *(1)
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:),iAp(:),iBp(:)
		real*4,pointer::sRp(:),sAp(:),sBp(:)
		real*8,pointer::dRp(:),dAp(:),dBp(:)
		complex*8,pointer::cRp(:),cAp(:),cBp(:)
		complex*16,pointer::zRp(:),zAp(:),zBp(:)
		select case(classType)
			case(1)
				call A%pointer(iAp)
				call B%pointer(iBp)
				call Res%pointer(iRp)
				if(NewResFlag)then
					iRp(1)=( iselect(alpha)*iAp(1)*iBp(1) )
				else
					iRp(1)=( iselect(alpha)*iAp(1)*iBp(1) )+( iselect(beta)*iRp(1) )
				end if
			case(2)
				call A%pointer(sAp)
				call B%pointer(sBp)
				call Res%pointer(sRp)
				if(NewResFlag)then
					sRp(1)=( sselect(alpha)*sAp(1)*sBp(1) )
				else
					sRp(1)=( sselect(alpha)*sAp(1)*sBp(1) )+( sselect(beta)*sRp(1) )
				end if
			case(3)
				call A%pointer(dAp)
				call B%pointer(dBp)
				call Res%pointer(dRp)
				if(NewResFlag)then
					dRp(1)=( dselect(alpha)*dAp(1)*dBp(1) )
				else
					dRp(1)=( dselect(alpha)*dAp(1)*dBp(1) )+( dselect(beta)*dRp(1) )
				end if
			case(4)
				call A%pointer(cAp)
				call B%pointer(cBp)
				call Res%pointer(cRp)
				if(NewResFlag)then
					cRp(1)=( cselect(alpha)*cAp(1)*cBp(1) )
				else
					cRp(1)=( cselect(alpha)*cAp(1)*cBp(1) )+( cselect(beta)*cRp(1) )
				end if
			case(5)
				call A%pointer(zAp)
				call B%pointer(zBp)
				call Res%pointer(zRp)
				if(NewResFlag)then
					zRp(1)=( zselect(alpha)*zAp(1)*zBp(1) )
				else
					zRp(1)=( zselect(alpha)*zAp(1)*zBp(1) )+( zselect(beta)*zRp(1) )
				end if
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorCase56(Res,A,B,alpha,beta,classType,NewResFlag)!Number*Tensor,(1) *(3) or (1) *(2,3) 
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:),iAp(:),iBp(:)
		real*4,pointer::sRp(:),sAp(:),sBp(:)
		real*8,pointer::dRp(:),dAp(:),dBp(:)
		complex*8,pointer::cRp(:),cAp(:),cBp(:)
		complex*16,pointer::zRp(:),zAp(:),zBp(:)
		select case(classType)
			case(1)
				call A%pointer(iAp)
				call B%pointer(iBp)
				call Res%pointer(iRp)
				if(NewResFlag)then
					iRp=( iselect(alpha)*iAp(1)*iBp )
				else
					iRp=( iselect(alpha)*iAp(1)*iBp )+( iselect(beta)*iRp )
				end if
			case(2)
				call A%pointer(sAp)
				call B%pointer(sBp)
				call Res%pointer(sRp)
				if(NewResFlag)then
					sRp=( sselect(alpha)*sAp(1)*sBp )
				else
					sRp=( sselect(alpha)*sAp(1)*sBp )+( sselect(beta)*sRp )
				end if
			case(3)
				call A%pointer(dAp)
				call B%pointer(dBp)
				call Res%pointer(dRp)
				if(NewResFlag)then
					dRp=( dselect(alpha)*dAp(1)*dBp )
				else
					dRp=( dselect(alpha)*dAp(1)*dBp )+( dselect(beta)*dRp )
				end if
			case(4)
				call A%pointer(cAp)
				call B%pointer(cBp)
				call Res%pointer(cRp)
				if(NewResFlag)then
					cRp=( cselect(alpha)*cAp(1)*cBp )
				else
					cRp=( cselect(alpha)*cAp(1)*cBp )+( cselect(beta)*cRp )
				end if

			case(5)
				call A%pointer(zAp)
				call B%pointer(zBp)
				call Res%pointer(zRp)
				if(NewResFlag)then
					zRp=( zselect(alpha)*zAp(1)*zBp )
				else
					zRp=( zselect(alpha)*zAp(1)*zBp )+( zselect(beta)*zRp )
				end if
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorCase78(Res,A,B,alpha,beta,classType,NewResFlag)!Tensor*number,(3) *(1) or (2,3) *(1) 
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType
		logical,intent(in)::NewResFlag
		integer,pointer::iRp(:),iAp(:),iBp(:)
		real*4,pointer::sRp(:),sAp(:),sBp(:)
		real*8,pointer::dRp(:),dAp(:),dBp(:)
		complex*8,pointer::cRp(:),cAp(:),cBp(:)
		complex*16,pointer::zRp(:),zAp(:),zBp(:)
		select case(classType)
			case(1)
				call B%pointer(iBp)
				call A%pointer(iAp)
				call Res%pointer(iRp)
				if(NewResFlag)then
					iRp=( iselect(alpha)*iAp*iBp(1) )
				else
					iRp=( iselect(alpha)*iAp*iBp(1) )+( iselect(beta)*iRp )
				end if
			case(2)
				call B%pointer(sBp)
				call A%pointer(sAp)
				call Res%pointer(sRp)
				if(NewResFlag)then
					sRp=( sselect(alpha)*sAp*sBp(1) )
				else
					sRp=( sselect(alpha)*sAp*sBp(1) )+( sselect(beta)*sRp )
				end if
			case(3)
				call B%pointer(dBp)
				call A%pointer(dAp)
				call Res%pointer(dRp)
				if(NewResFlag)then
					dRp=( dselect(alpha)*dAp*dBp(1) )
				else
					dRp=( dselect(alpha)*dAp*dBp(1) )+( dselect(beta)*dRp )
				end if
			case(4)
				call B%pointer(cBp)
				call A%pointer(cAp)
				call Res%pointer(cRp)
				if(NewResFlag)then
					cRp=( cselect(alpha)*cAp*cBp(1) )
				else
					cRp=( cselect(alpha)*cAp*cBp(1) )+( cselect(beta)*cRp )
				end if
			case(5)
				call B%pointer(zBp)
				call A%pointer(zAp)
				call Res%pointer(zRp)
				if(NewResFlag)then
					zRp=( zselect(alpha)*zAp*zBp(1) )
				else
					zRp=( zselect(alpha)*zAp*zBp(1) )+( zselect(beta)*zRp )
				end if
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorCase9(Res,A,B,alpha,beta,classType,NewResFlag)!vec*vec
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
		select case(classType)
			case(1)
				call Res%pointer(iRp0)
				call A%pointer(iAp0)
				call B%pointer(iBp0)
				iRp(1:1,1:1)=>iRp0(1:1)
				iAp(1:1,1:size(iAp0))=>iAp0(1:size(iAp0))
				iBp(1:size(iBp0),1:1)=>iBp0(1:size(iBp0))
				if(NewResFlag)then
					iRp=iselect(alpha)*matmul(iAp,iBp)
				else
					iRp=iselect(alpha)*matmul(iAp,iBp)+iselect(beta)*iRp
				end if
			case(2)
				call Res%pointer(sRp)
				call A%pointer(sAp)
				call B%pointer(sBp)
				if(NewResFlag)then
					sRp(1)=( sselect(alpha) * sdot(size(sAp), sAp, 1, sBp, 1) )
				else
					sRp(1)=( sselect(alpha) * sdot(size(sAp), sAp, 1, sBp, 1) ) + ( sselect(beta)*sRp(1) )
				end if
			case(3)
				call Res%pointer(dRp)
				call A%pointer(dAp)
				call B%pointer(dBp)
				if(NewResFlag)then
					dRp(1)=( dselect(alpha) * ddot(size(dAp), dAp, 1, dBp, 1) )
				else
					dRp(1)=( dselect(alpha) * ddot(size(dAp), dAp, 1, dBp, 1) ) + ( dselect(beta)*dRp(1) )
				end if
			case(4)
				call Res%pointer(cRp)
				call A%pointer(cAp)
				call B%pointer(cBp)
				if(NewResFlag)then
					cRp(1)=( cselect(alpha) * cdotu(size(cAp), cAp, 1, cBp, 1) )
				else
					cRp(1)=( cselect(alpha) * cdotu(size(cAp), cAp, 1, cBp, 1) ) + ( cselect(beta)*cRp(1) )
				end if
			case(5)
				call Res%pointer(zRp)
				call A%pointer(zAp)
				call B%pointer(zBp)
				if(NewResFlag)then
					zRp(1)=( zselect(alpha) * zdotu(size(zAp), zAp, 1, zBp, 1) )
				else
					zRp(1)=( zselect(alpha) * zdotu(size(zAp), zAp, 1, zBp, 1) ) + ( zselect(beta)*zRp(1) )
				end if
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorCase10(Res,A,B,alpha,beta,LD1,LD2,classType,NewResFlag)!vec*matrix
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,LD1,LD2
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
		select case(classType)
			case(1)
				call A%pointer(iAp0)
				call B%pointer(iBp0)
				call Res%pointer(iRp0)
				iRp(1:1,1:LD2)=>iRp0(1:LD1)
				iAp(1:1,1:LD1)=>iAp0(1:LD1)
				iBp(1:LD1,1:LD2)=>iBp0(1:size(iBp0))
				if(NewResFlag)then
					iRp=iselect(alpha)*matmul(iAp,iBp)
				else
					iRp=iselect(alpha)*matmul(iAp,iBp)+iselect(beta)*iRp
				end if
			case(2)
				call A%pointer(sAp)
				call B%pointer(sBp)
				call Res%pointer(sRp)
				call SGEMV('T', LD1, LD2, sselect(alpha), sBp, LD1, sAp, 1, sselect(beta), sRp, 1)
			case(3)
				call A%pointer(dAp)
				call B%pointer(dBp)
				call Res%pointer(dRp)
				call DGEMV('T', LD1, LD2, dselect(alpha), dBp, LD1, dAp, 1, dselect(beta), dRp, 1)
			case(4)
				call A%pointer(cAp)
				call B%pointer(cBp)
				call Res%pointer(cRp)
				call CGEMV('T', LD1, LD2, cselect(alpha), cBp, LD1, cAp, 1, cselect(beta), cRp, 1)
			case(5)
				call A%pointer(zAp)
				call B%pointer(zBp)
				call Res%pointer(zRp)
				call ZGEMV('T', LD1, LD2, zselect(alpha), zBp, LD1, zAp, 1, zselect(beta), zRp, 1)
			case default
		end select
		return
	end subroutine

	subroutine ProductTensorCase11(Res,A,B,alpha,beta,LD1,LD2,classType,NewResFlag)!matrix*vec
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,LD1,LD2
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
		select case(classType)
			case(1)
				call A%pointer(iAp0)
				call B%pointer(iBp0)
				call Res%pointer(iRp0)
				iRp(1:1,1:LD1)=>iRp0
				iAp(1:LD1,1:LD2)=>iAp0
				iBp(1:LD2,1:1)=>iBp0
				if(NewResFlag)then
					iRp=iselect(alpha)*matmul(iAp,iBp)
				else
					iRp=iselect(alpha)*matmul(iAp,iBp)+iselect(beta)*iRp
				end if
			case(2)
				call A%pointer(sAp)
				call B%pointer(sBp)
				call Res%pointer(sRp)
				call SGEMV('N', LD1, LD2, sselect(alpha), sAp, LD1, sBp, 1, sselect(beta), sRp, 1)
			case(3)
				call A%pointer(dAp)
				call B%pointer(dBp)
				call Res%pointer(dRp)
				call DGEMV('N', LD1, LD2, dselect(alpha), dAp, LD1, dBp, 1, dselect(beta), dRp, 1)
			case(4)
				call A%pointer(cAp)
				call B%pointer(cBp)
				call Res%pointer(cRp)
				call CGEMV('N', LD1, LD2, cselect(alpha), cAp, LD1, cBp, 1, cselect(beta), cRp, 1)
			case(5)
				call A%pointer(zAp)
				call B%pointer(zBp)
				call Res%pointer(zRp)
				call ZGEMV('N', LD1, LD2, zselect(alpha), zAp, LD1, zBp, 1, zselect(beta), zRp, 1)
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine ProductTensorCase12(Res,A,B,alpha,beta,LD1,LD2,LD3,classType,NewResFlag)!matrix*matrix
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: A,B
		class(*),intent(in)::alpha,beta
		integer,intent(in)::classType,LD1,LD2,LD3
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
		select case(classType)
			case(1)
				call A%pointer(iAp0)
				call B%pointer(iBp0)
				call Res%pointer(iRp0)
				iRp(1:LD1,1:LD3)=>iRp0
				iAp(1:LD1,1:LD2)=>iAp0
				iBp(1:LD2,1:LD3)=>iBp0
				if(NewResFlag)then
					iRp=iselect(alpha)*matmul(iAp,iBp)
				else
					iRp=iselect(alpha)*matmul(iAp,iBp)+iselect(beta)*iRp
				end if
			case(2)
				call A%pointer(sAp)
				call B%pointer(sBp)
				call Res%pointer(sRp)
				call SGEMM('N','N',LD1,LD3,LD2,sselect(alpha),sAp,LD1,sBp,LD2,sselect(beta),sRp,LD1)
			case(3)
				call A%pointer(dAp)
				call B%pointer(dBp)
				call Res%pointer(dRp)
				call DGEMM('N','N',LD1,LD3,LD2,dselect(alpha),dAp,LD1,dBp,LD2,dselect(beta),dRp,LD1)
			case(4)
				call A%pointer(cAp)
				call B%pointer(cBp)
				call Res%pointer(cRp)
				call CGEMM('N','N',LD1,LD3,LD2,cselect(alpha),cAp,LD1,cBp,LD2,cselect(beta),cRp,LD1)
			case(5)
				call A%pointer(zAp)
				call B%pointer(zBp)
				call Res%pointer(zRp)
				call ZGEMM('N','N',LD1,LD3,LD2,zselect(alpha),zAp,LD1,zBp,LD2,zselect(beta),zRp,LD1)
			case default
				call writemess('ERROR case in product',-1)
				call error_stop
		end select
		return
	end subroutine

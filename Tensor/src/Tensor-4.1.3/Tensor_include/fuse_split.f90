


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
	!FuseDimension,
		!    When fusing tensors: i1i2i3  fused as I
		!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
		!    regard A as dim1*Dim2 Tensor, fuse it into Dim1*Dim2 T tensor
		!
		!    1.read the data in [old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi], Ap
		!  point to the block of [old_i1i2i3,J] by call A%pointer([dim1,Dim2],Ap,[old_i1i2i3,J])
		!
		!    2. allocate momery for Result tensor, the block of [New_I,J] with total length
		!  of New_degI*DegJ, then Rp point to [New_I,J] by 
		!  call Res%pointer([Dim1,Dim2],Rp,[New_degI,degJ],[New_I,J]) 
		!
		!    3. Rp( Starti:endi , : ) = Ap
		!

	subroutine FuseDimension(outQN,Order,TDim,ith,jth,Orderrow,indices,NewDeg,NewQN,TMPQN)
		type(QuanNum),intent(inout)::outQN
		integer,intent(inout)::Order(:,:),Orderrow,indices(:),NewDeg(:)
		type(dimension),intent(in)::TDim
		real*4,intent(inout)::NewQN(:),TMPQN(:)
		integer,intent(in)::ith,jth
		integer::i,j,dimi,ii,NewI,degi,length_in_used
		integer,pointer::dim(:),deg(:)
		real*4::QNi
		real*4,pointer::QN(:)
		logical::Goon
		if(size(indices).ne.(jth-ith+1))then
			call writemess('ERROR in FuseDimension',-1)
			call error_stop
		end if
		indices=1
		call TDim%pointDim(dim)
		dimi=product(dim(ith:jth))
		Orderrow=0
		length_in_used=0
		do i=1,dimi
			Orderrow=Orderrow+1

			ii=ith-1
			do j=1,size(indices)
				ii=ii+1
				call TDim%pointQN(QN,ii)
				TMPQN(j)=QN(indices(j))
			end do
			call SymmetryNewQaunNum(QNi,TMPQN)

			Degi=TDim%getBlockDim(ith,jth,indices)

			call check_Fusedindex(QNi,NewQN,NewDeg,length_in_used,NewI)
			call Check_push_back(i,Degi,NewI,Order,Orderrow,NewDeg)
			goon=index_counter(indices,dim(ith:jth))
		end do

		call outQN%setQN(NewQN(1:length_in_used))
		call outQN%setDeg(NewDeg(1:length_in_used))
		do i=1,Orderrow
			Order(i,4)=NewDeg(Order(i,3))
		end do
		if(TDim%getFermiFlag())then
			call outQN%setFermiArrow(TDim%getArrow(ith))
		end if
		return
	end subroutine

	!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
	!             1        2        3     4        5     6

	subroutine Check_push_back(ith,Degi,NewI,Order,Orderrow,NewDeg)
		integer,intent(in)::ith,Degi,NewI,Orderrow
		integer,intent(inout)::Order(:,:),NewDeg(:)
		integer::Si,Ei,i
		do i=Orderrow-1,1,-1
			if(NewI.eq.Order(i,3))then
				Si=Order(i,6)+1
				Ei=Order(i,6)+Degi
				if(NewDeg(NewI).lt.Ei)NewDeg(NewI)=Ei
				Order(Orderrow,:)=[ith,Degi,NewI,0,Si,Ei]
				return
			end if
		end do
		Si=1
		Ei=Degi
		if(NewDeg(NewI).lt.Ei)NewDeg(NewI)=Ei
		Order(Orderrow,:)=[ith,Degi,NewI,0,Si,Ei]
		return
	end subroutine

	subroutine check_Fusedindex(QN,NewQN,NewDeg,length_in_used,indexJ)
		real*4,intent(in)::QN
		real*4,intent(inout)::NewQN(:)
		integer,intent(inout)::length_in_used,NewDeg(:),indexJ
		integer::i
		if(length_in_used.eq.0)then
			length_in_used=length_in_used+1
			NewQN(length_in_used)=QN
			NewDeg(length_in_used)=-1
			indexJ=1
			return
		end if
		do i=1,length_in_used
			if(QN.equ.NewQN(i))then
				indexJ=i
				return
			end if
		end do
		length_in_used=length_in_used+1
		indexJ=length_in_used
		NewQN(length_in_used)=QN
		NewDeg(length_in_used)=-1
		return
	end subroutine


	!fuseinfo2 [i1,i2,i3,j,k,l]-->[I,j,k,l]
	!outinfo: old_index,Deg_i1i2i3,I,j,DegI,DegJ,si,ei
	!              1     2         3 4   5  6    7  8
	!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
	!             1        2        3     4        5     6

	subroutine fuseinfo2(outinfo,BlockN,Order,T,ith,jth,NewDim1,NewDim2,inforow,indices)
		integer,intent(inout)::inforow,indices(:),outinfo(:,:)
		type(Tensor),intent(in)::T
		integer,target,intent(inout)::BlockN(:)
		integer,intent(in)::ith,jth,Order(:,:),NewDim1,NewDim2
		integer::i,dim1,dim2,dimi,dimj
		integer,pointer::dim(:),deg(:),Nij(:,:)
		integer::ii,Degj
		integer::Si,Ei,rank
		logical::goon
		if(ith.ne.1)then
			call writemess('ERROR in fuseinfo2',-1)
			call error_Stop
		end if
		call T%pointDim(dim)
		rank=T%getRank()
		if(size(indices).ne.rank)then
			call writemess('ERROR in fuseinfo2',-1)
			call error_Stop
		end if
		if(NewDim1*NewDim2.ne.size(BlockN))then
			call writemess('ERROR in fuseinfo2',-1)
			call error_Stop
		end if
		dim1=product(dim(1:jth))
		dim2=product(dim(jth+1:rank))
		Nij(1:NewDim1,1:NewDim2)=>BlockN

		ii=0
		inforow=0
		BlockN=0
		indices=1
		do dimj=1,dim2
			do dimi=1,dim1
				ii=ii+1
				if(T%getFlag(ii))then
					inforow=inforow+1
					Degj=T%getBlockDim(jth+1,rank,indices(jth+1:rank))
					!outinfo: old_index,Deg_i1i2i3,I,j,DegI,DegJ,si,ei
					!              1     2         3 4   5  6    7  8
					!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
					!             1        2        3     4        5     6
					do i=1,size(Order,1)
						if(dimi.eq.Order(i,1))then
							outinfo(inforow,1:8)=[ii,Order(i,2:3),dimj,Order(i,4),Degj,Order(i,5:6)]
							Nij(Order(i,3),dimj)=Order(i,4)*Degj
							exit
						end if
					end do
				end if
				goon=index_counter(indices,dim)
			end do
		end do
		return
	end subroutine

	!fuseinfo1 [i1,i2,i3]-->[I]
	!outinfo: old_index,Deg_i1i2i3,I,DegI,si,ei
	!              1         2     3   4  5 ,6
	!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
	!             1        2        3     4        5     6


	subroutine fuseinfo1(outinfo,BlockN,Order,T,ith,jth,inforow)
		integer,intent(inout)::inforow,BlockN(:),outinfo(:,:)
		type(Tensor),intent(in)::T
		integer,intent(in)::ith,jth,Order(:,:)
		integer::i,dim1,dimi
		integer,pointer::dim(:)
		integer::ii,rank
		rank=T%getRank()
		if(ith.ne.1)then
			call writemess('ERROR in fuseinfo1',-1)
			call error_Stop
		end if
		if(jth.ne.rank)then
			call writemess('ERROR in fuseinfo1',-1)
			call error_Stop
		end if
		call T%pointDim(dim)
		
		dim1=product(dim)

		ii=0
		inforow=0
		BlockN=0
		do dimi=1,dim1
			ii=ii+1
			if(T%getFlag(ii))then
				inforow=inforow+1
				!outinfo: old_index,Deg_i1i2i3,I,DegI,si,ei
				!              1         2     3   4  5 ,6
				!Order=[old_i1i2i3,Deg_i1i2i3,New_I,New_degI,Starti,endi]
				!             1        2        3     4        5     6
				do i=1,size(Order,1)
					if(dimi.eq.Order(i,1))then
						outinfo(inforow,1:6)=[ii,Order(i,2:3),Order(i,4),Order(i,5:6)]
						BlockN(Order(i,3))=Order(i,4)
						exit
					end if
				end do
			end if
		end do
		
		return
	end subroutine


	!fuseinfo2Rank [i,j1,j2,j3]-->[j,J]
	!outinfo: old_index,Deg_j1j2j3,i,J,DegI,DegJ,si,ei
	!              1         2     3 4  5    6    7  8
	!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
	!             1        2        3     4        5     6

	subroutine fuseinfo2Rank(outinfo,BlockN,Order,T,ith,jth,NewDim1,NewDim2,inforow,indices)
		integer,intent(inout)::inforow,indices(:),outinfo(:,:)
		type(Tensor),intent(in)::T
		integer,target,intent(inout)::BlockN(:)
		integer,intent(in)::ith,jth,Order(:,:),NewDim1,NewDim2
		integer::i,dim1,dim2,dimi,dimj
		integer,pointer::dim(:),Nij(:,:)
		integer::ii,Degi,Si,Ei,rank
		logical::goon
		rank=T%getRank()
		if(jth.ne.rank)then
			call writemess('ERROR in fuseinfo2Rank',-1)
			call error_Stop
		end if
		if(size(indices).ne.rank)then
			call writemess('ERROR in fuseinfo2Rank',-1)
			call error_Stop
		end if
		if(NewDim1*NewDim2.ne.size(BlockN))then
			call writemess('ERROR in fuseinfo2Rank',-1)
			call writemess('NewDim1='+NewDim1,-1)
			call writemess('NewDim2='+NewDim2,-1)
			call writemess('size(BlockN)='+size(BlockN),-1)
			call error_Stop
		end if
		call T%pointDim(dim)
		Nij(1:NewDim1,1:NewDim2)=>BlockN

		dim1=product(dim(1:ith-1))
		dim2=product(dim(ith:rank))

		ii=0
		inforow=0
		indices=1
		BlockN=0
		do dimj=1,dim2
			do dimi=1,dim1
				ii=ii+1
				if(T%getFlag(ii))then
					inforow=inforow+1
					Degi=T%getBlockDim(1,ith-1,indices(1:ith-1))
					!outinfo: old_index,Deg_j1j2j3,i,J,DegI,DegJ,si,ei
					!              1         2     3 4  5    6    7  8
					!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
					!             1        2        3     4        5     6
					do i=1,size(Order,1)
						if(dimj.eq.Order(i,1))then
							outinfo(inforow,1:8)=[ii,Order(i,2),dimi,Order(i,3),degi,Order(i,4),Order(i,5:6)]
							Nij(dimi,Order(i,3))=degi*Order(i,4)
							exit
						end if
					end do
				end if
				goon=index_counter(indices,dim)
			end do
		end do
		
		return
	end subroutine

	!fuseinfo3 [i,j1,j2,j3,k,l]-->[j,J,k,l]
	!outinfo: old_index,Deg_j1j2j3,i,J,k,DegI,DegJ,DegK,si,ei
	!              1         2     3 4 5  6    7    8   9  10
	!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
	!             1        2        3     4        5     6

	subroutine fuseinfo3(outinfo,BlockN,Order,T,ith,jth,NewDim1,NewDim2,NewDim3,inforow,indices)
		integer,intent(inout)::inforow,indices(:),outinfo(:,:)
		type(Tensor),intent(in)::T
		integer,target,intent(inout)::BlockN(:)
		integer,intent(in)::ith,jth,Order(:,:),NewDim1,NewDim2,NewDim3
		integer::i,dim1,dim2,dim3,dimi,dimj,dimk
		integer,pointer::dim(:),Nij(:,:,:)
		integer::ii,iii,Degi,Degk
		integer::Si,Ei,rank
		logical::goon
		rank=T%getRank()
		if(ith.eq.1)then
			call writemess('ERROR in fuseinfo3',-1)
			call error_Stop
		end if
		if(jth.eq.rank)then
			call writemess('ERROR in fuseinfo3',-1)
			call error_Stop
		end if
		if(size(indices).ne.rank)then
			call writemess('ERROR in fuseinfo3',-1)
			call error_Stop
		end if
		if(NewDim1*NewDim2*NewDim3.ne.size(BlockN))then
			call writemess('ERROR in fuseinfo3',-1)
			call error_Stop
		end if
		call T%pointDim(dim)
		Nij(1:NewDim1,1:NewDim2,1:NewDim3)=>BlockN
		
		dim1=product(dim(1:ith-1))
		dim2=product(dim(ith:jth))
		dim3=product(dim(jth+1:rank))

		ii=0
		indices=1
		inforow=0
		BlockN=0
		do dimk=1,dim3
			do dimj=1,dim2
				do dimi=1,dim1
					ii=ii+1
					if(T%getFlag(ii))then
						inforow=inforow+1
						Degi=T%getBlockDim(1,ith-1,indices(1:ith-1))
						Degk=T%getBlockDim(jth+1,rank,indices(jth+1:rank))
						!outinfo: old_index,Deg_j1j2j3,i,J,k,DegI,DegJ,DegK,si,ei
						!              1         2     3 4 5  6    7    8   9  10
						!Order=[old_j1j2j3,Deg_j1j2j3,New_J,New_degJ,Starti,endi]
						!             1        2        3     4        5     6
						do i=1,size(Order,1)
							if(dimj.eq.Order(i,1))then
								outinfo(inforow,1:10)=[ii,Order(i,2),dimi,Order(i,3),dimk,degi,Order(i,4),degk,Order(i,5:6)]
								Nij(dimi,Order(i,3),dimk)=degi*Order(i,4)*degk
								exit
							end if
						end do
					end if
					goon=index_counter(indices,dim)
				end do
			end do
		end do
		
		return
	end subroutine

	subroutine Find_Data_in_order(outindex,order,ith,targetObj,leninuse)
		integer,intent(inout)::outindex(:),leninuse
		integer,intent(in)::order(:,:),ith,targetObj
		integer::i
		leninuse=0
		do i=1,size(order,1)
			if(order(i,ith).eq.targetObj)then
				leninuse=leninuse+1
				outindex(leninuse)=i
			end if
		end do
		return
	end subroutine








#define QuantumFuseTYPE QuantumFusei
#define cpydataDEDATANAME cpydatai
#define QuantumSplitTYPE QuantumSpliti
#define DataTYPE integer
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE


#define QuantumFuseTYPE QuantumFuses
#define cpydataDEDATANAME cpydatas
#define QuantumSplitTYPE QuantumSplits
#define DataTYPE real*4
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#define QuantumFuseTYPE QuantumFused
#define cpydataDEDATANAME cpydatad
#define QuantumSplitTYPE QuantumSplitd
#define DataTYPE real*8
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#define QuantumFuseTYPE QuantumFusec
#define cpydataDEDATANAME cpydatac
#define QuantumSplitTYPE QuantumSplitc
#define DataTYPE complex*8
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#define QuantumFuseTYPE QuantumFusez
#define cpydataDEDATANAME cpydataz
#define QuantumSplitTYPE QuantumSplitz
#define DataTYPE complex*16
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#define QuantumFuseTYPE QuantumFusel
#define cpydataDEDATANAME cpydatal
#define QuantumSplitTYPE QuantumSplitl
#define DataTYPE logical
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#define QuantumFuseTYPE QuantumFusea
#define cpydataDEDATANAME cpydataa
#define QuantumSplitTYPE QuantumSplita
#define DataTYPE character(len=characterlen)
#include "templet/fuse_split0.f90"
#undef QuantumFuseTYPE
#undef cpydataDEDATANAME
#undef QuantumSplitTYPE
#undef DataTYPE

#include "fuse_split_include.f90"

	subroutine QuantumFuse(A,ith,jth,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith,jth
		type(Tensor),intent(inout)::Res
		type(Dimension),optional,intent(inout)::SubDim
		type(Tensor),optional,intent(inout)::ResOrder
		select case(A%getType())
			case(1)
				call QuantumFusei(A,ith,jth,Res,SubDim,ResOrder)
			case(2)
				call QuantumFuses(A,ith,jth,Res,SubDim,ResOrder)
			case(3)
				call QuantumFused(A,ith,jth,Res,SubDim,ResOrder)
			case(4)
				call QuantumFusec(A,ith,jth,Res,SubDim,ResOrder)
			case(5)
				call QuantumFusez(A,ith,jth,Res,SubDim,ResOrder)
			case(6)
				call QuantumFusel(A,ith,jth,Res,SubDim,ResOrder)
			case(7)
				call QuantumFusea(A,ith,jth,Res,SubDim,ResOrder)
			case default
				call QuantumFuseClassData(A,ith,jth,Res,SubDim,ResOrder)
		end select
	end subroutine
	subroutine QuantumSplit(A,legi,Res,SubDim,ResOrder)
		class(Tensor),intent(in)::A
		integer,intent(in)::legi
		type(Tensor),intent(inout)::Res
		type(Dimension),intent(in)::SubDim
		type(Tensor),intent(in)::ResOrder
		select case(A%getType())
			case(1)
				call QuantumSpliti(A,legi,Res,SubDim,ResOrder)
			case(2)
				call QuantumSplits(A,legi,Res,SubDim,ResOrder)
			case(3)
				call QuantumSplitd(A,legi,Res,SubDim,ResOrder)
			case(4)
				call QuantumSplitc(A,legi,Res,SubDim,ResOrder)
			case(5)
				call QuantumSplitz(A,legi,Res,SubDim,ResOrder)
			case(6)
				call QuantumSplitl(A,legi,Res,SubDim,ResOrder)
			case(7)
				call QuantumSplita(A,legi,Res,SubDim,ResOrder)
			case default
				call QuantumSplitClassData(A,legi,Res,SubDim,ResOrder)
		end select
	end subroutine

	subroutine Fuselegs_sign(Block,dimen,legstart,legendi,legCheckith,indices,maxinde)
		integer,intent(inout)::Block(:),indices(:),legCheckith(:)
		integer,intent(in)::maxinde(:),legstart,legendi
		type(Dimension),intent(in)::dimen
		integer::i,j,ii,ii_in_use,rank
		rank=dimen%getRank()
		ii=1
		legCheckith(ii)=legstart
		do i=legstart+1,legendi
			ii=ii+1
			legCheckith(ii)=i
			ii_in_use=ii
			call Fuse_sign(Block,dimen,i,legCheckith(1:ii_in_use),indices,maxinde)
		end do
		return
	end subroutine
	
	subroutine Fuse_sign(Block,dimen,ithleg,legCheckith,indices,maxinde)
		integer,intent(inout)::Block(:),indices(:)
		integer,intent(in)::maxinde(:),legCheckith(:),ithleg
		type(Dimension),intent(in)::dimen
		integer::iQN,i
		logical::goon,ifp,testrule
		indices=1
		goon=.true.
		i=0
		do while (goon)
			i=i+1
			if(block(i).ne.0)then
				if(check_same_name_flag)then
					testrule = if_symmetry_Rule(dimen,indices)
					if(.not.testrule)then
						call writemess('ERRRO symmetry rule',-1)
						call writemess('indices=')
						call writemess(indices,'I4')
						call dimen%diminfo(.true.)
						call error_stop
					end if
				end if
				call QaunNumParity(iQN,dimen,ithleg,indices(ithleg) )
				if(iQN .eq.(-1))then
					call ifParity(ifp,dimen,indices,legCheckith)
					if(ifp) then
						block(i)=-block(i)
					end if
				end if
			end if
			goon=index_counter(indices,maxinde)
		end do
		return
	end subroutine
	subroutine Fuse_for_fermi_subroutine(T,ith,jth,Check)
		type(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		logical,intent(in)::Check
		integer::i,length,rank
		integer,pointer::Arrow(:),ei(:),legCheckith(:),indices(:),maxinde(:)
		if(.not.Check) return
		if(.not.T%getFermiFlag())return

		call T%pointArrow(Arrow)
		do i=ith+1,jth
			if(arrow(i).ne.arrow(ith))then
				call writemess('ERROR in fusing fermi-Tensor',-1)
				call writemess(' the fusing legs should have the same fermi arrow',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call T%info(10)
				call error_stop
			end if
		end do
		if(arrow(ith).gt.0)return

		length=jth-ith+1
		rank=T%getRank()
		if(length.le.0)then
			call writemess('ERROR in Fuse_for_fermi_subroutine',-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,length+rank)
		call WorkingMemory%get_memory(legCheckith,length)
		call WorkingMemory%get_memory(indices,rank)
		call T%pointDim(maxinde)
		call T%Data%pointEndi(ei)

		call Fuselegs_sign(ei,T%Dimension,ith,jth,legCheckith,indices,maxinde)
		call WorkingMemory%free()
		call fix_fermi_sign(T)
		return
	end subroutine


	subroutine FuseLegs0(T,ith,jth,fermiCheck_)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		logical,optional,intent(in)::fermiCheck_
		integer,pointer::dim(:)
		character(len=len_of_name),pointer::Names(:)
		integer::rank,i,ii
		logical::NameFlag,FermiCheck
		integer,allocatable::TMPDim(:)
		character(len=len_of_name),allocatable::TMPNames(:)
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if

		if(.not.T%getSymmetryFlag())then
			if(ith.eq.jth)return
			if(deallocate_memory_flag)then
				allocate(TMPDim(rank))
				allocate(TMPNames(rank))
				TMPDim=T%Dimension
				NameFlag=T%getNameFlag()
				if(NameFlag)then
					do i=1,rank
						TMPNames(i)=T%getName(i)
					end do
				end if
				if((ith.gt.1).and.(jth.lt.rank))then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(ith.eq.1)then
					T%Dimension=[product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(jth.eq.rank)then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth))]
				else
					call writemess('ERROR in fusing legs',-1)
					call writemess('ith='+ith,-1)
					call writemess('jth='+jth,-1)
					call error_stop
				end if
				if(.not.NameFlag)return
				ii=ith
				call T%setName(ii,'Fuse.leg'+ii)
				do i=jth+1,rank
					ii=ii+1
					call T%setName(ii,TMPNames(i))
				end do
				return
			end if
			call T%pointDim(Dim)
			NameFlag=T%getNameFlag()
			if(NameFlag)call T%pointName(Names)
			if((ith.gt.1).and.(jth.lt.rank))then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth)),dim(jth+1:rank)]
			else if(ith.eq.1)then
				T%Dimension=[product(dim(ith:jth)),dim(jth+1:rank)]
			else if(jth.eq.rank)then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth))]
			else
				call writemess('ERROR in fusing legs',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call error_stop
			end if
			if(.not.NameFlag)return
			ii=ith
			call T%setName(ii,'Fuse.leg'+ii)
			do i=jth+1,rank
				ii=ii+1
				call T%setName(ii,Names(i))
			end do
			return
		end if

		if(ith.eq.jth) return
		if(present(fermiCheck_))then
			fermiCheck=fermiCheck_
		else
			fermiCheck=.true.
		end if

		TMPFuseTensor=T	
		call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
		call QuantumFuse(TMPFuseTensor,ith,jth,T)
		return
	end subroutine

	subroutine FuseLegs1(T,ith,jth,SubDim,ResOrder,fermiCheck_)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		type(Dimension),intent(inout)::SubDim
		type(Tensor),intent(inout)::ResOrder
		logical,optional,intent(in)::fermiCheck_
		integer,pointer::dim(:)
		character(len=len_of_name),pointer::Names(:)
		integer::rank,i,ii
		logical::NameFlag,fermiCheck
		integer,allocatable::TMPDim(:)
		character(len=len_of_name),allocatable::TMPNames(:)
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if

		if(.not.T%getSymmetryFlag())then
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			CALL ResOrder%empty
			if(ith.eq.jth)return
			if(deallocate_memory_flag)then
				allocate(TMPDim(rank))
				allocate(TMPNames(rank))
				TMPDim=T%Dimension
				NameFlag=T%getNameFlag()
				if(NameFlag)then
					do i=1,rank
						TMPNames(i)=T%getName(i)
					end do
				end if
				if((ith.gt.1).and.(jth.lt.rank))then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(ith.eq.1)then
					T%Dimension=[product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(jth.eq.rank)then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth))]
				else
					call writemess('ERROR in fusing legs',-1)
					call writemess('ith='+ith,-1)
					call writemess('jth='+jth,-1)
					call error_stop
				end if
				if(.not.NameFlag)return
				ii=ith
				call T%setName(ii,'Fuse.leg'+ii)
				do i=jth+1,rank
					ii=ii+1
					call T%setName(ii,TMPNames(i))
				end do
				return
			end if
			call T%pointDim(Dim)
			NameFlag=T%getNameFlag()
			if(NameFlag)call T%pointName(Names)
			if((ith.gt.1).and.(jth.lt.rank))then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth)),dim(jth+1:rank)]
			else if(ith.eq.1)then
				T%Dimension=[product(dim(ith:jth)),dim(jth+1:rank)]
			else if(jth.eq.rank)then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth))]
			else
				call writemess('ERROR in fusing legs',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call error_stop
			end if
			if(.not.NameFlag)return
			ii=ith
			call T%setName(ii,'Fuse.leg'+ii)
			do i=jth+1,rank
				ii=ii+1
				call T%setName(ii,Names(i))
			end do
			return
		end if

		if(ith.eq.jth) return
		if(present(fermiCheck_))then
			fermiCheck=fermiCheck_
		else
			fermiCheck=.true.
		end if

		TMPFuseTensor=T	
		call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
		call QuantumFuse(TMPFuseTensor,ith,jth,T,SubDim,ResOrder)
		return
	end subroutine

	subroutine FuseLegs2(T,ith,jth,SubDim,fermiCheck_)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		type(Dimension),intent(inout)::SubDim
		logical,optional,intent(in)::fermiCheck_
		integer,pointer::dim(:)
		character(len=len_of_name),pointer::Names(:)
		integer::rank,i,ii
		logical::NameFlag,fermiCheck
		integer,allocatable::TMPDim(:)
		character(len=len_of_name),allocatable::TMPNames(:)
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if

		if(.not.T%getSymmetryFlag())then
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			if(ith.eq.jth)return
			if(deallocate_memory_flag)then
				allocate(TMPDim(rank))
				allocate(TMPNames(rank))
				TMPDim=T%Dimension
				NameFlag=T%getNameFlag()
				if(NameFlag)then
					do i=1,rank
						TMPNames(i)=T%getName(i)
					end do
				end if
				if((ith.gt.1).and.(jth.lt.rank))then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(ith.eq.1)then
					T%Dimension=[product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(jth.eq.rank)then
					T%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth))]
				else
					call writemess('ERROR in fusing legs',-1)
					call writemess('ith='+ith,-1)
					call writemess('jth='+jth,-1)
					call error_stop
				end if
				if(.not.NameFlag)return
				ii=ith
				call T%setName(ii,'Fuse.leg'+ii)
				do i=jth+1,rank
					ii=ii+1
					call T%setName(ii,TMPNames(i))
				end do
				return
			end if
			call T%pointDim(Dim)
			NameFlag=T%getNameFlag()
			if(NameFlag)call T%pointName(Names)
			if((ith.gt.1).and.(jth.lt.rank))then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth)),dim(jth+1:rank)]
			else if(ith.eq.1)then
				T%Dimension=[product(dim(ith:jth)),dim(jth+1:rank)]
			else if(jth.eq.rank)then
				T%Dimension=[dim(1:ith-1),product(dim(ith:jth))]
			else
				call writemess('ERROR in fusing legs',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call error_stop
			end if
			if(.not.NameFlag)return
			ii=ith
			call T%setName(ii,'Fuse.leg'+ii)
			do i=jth+1,rank
				ii=ii+1
				call T%setName(ii,Names(i))
			end do
			return
		end if

		if(ith.eq.jth) return
		if(present(fermiCheck_))then
			fermiCheck=fermiCheck_
		else
			fermiCheck=.true.
		end if
		TMPFuseTensor=T	
		call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
		call QuantumFuse(TMPFuseTensor,ith,jth,T,SubDim)
		return
	end subroutine

	subroutine FuseLegs1Subroutine(T,ith,jth,fermiCheck,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith,jth
		logical,intent(in)::fermiCheck
		type(Tensor),intent(inout)::Res
		integer,pointer::dim(:)
		character(len=len_of_name),pointer::Names(:)
		integer::rank,i,ii
		logical::NameFlag
		integer,allocatable::TMPDim(:)
		character(len=len_of_name),allocatable::TMPNames(:)
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if

		if(.not.T%getSymmetryFlag())then
			if(ith.eq.jth)then
				Res=T
				return
			end if
			if(deallocate_memory_flag)then
				allocate(TMPDim(rank))
				allocate(TMPNames(rank))
				TMPDim=T%Dimension
				NameFlag=T%getNameFlag()
				if(NameFlag)then
					do i=1,rank
						TMPNames(i)=T%getName(i)
					end do
				end if
				if((ith.gt.1).and.(jth.lt.rank))then
					Res%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(ith.eq.1)then
					Res%Dimension=[product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(jth.eq.rank)then
					Res%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth))]
				else
					call writemess('ERROR in fusing legs',-1)
					call writemess('ith='+ith,-1)
					call writemess('jth='+jth,-1)
					call error_stop
				end if
				Res%Data=T%Data
				if(.not.NameFlag)return
				ii=ith
				call Res%setName(ii,'Fuse.leg'+ii)
				do i=jth+1,rank
					ii=ii+1
					call Res%setName(ii,TMPNames(i))
				end do
				return
			end if
			call T%pointDim(Dim)
			NameFlag=T%getNameFlag()
			if(NameFlag)call T%pointName(Names)
			if((ith.gt.1).and.(jth.lt.rank))then
				Res%Dimension=[dim(1:ith-1),product(dim(ith:jth)),dim(jth+1:rank)]
			else if(ith.eq.1)then
				Res%Dimension=[product(dim(ith:jth)),dim(jth+1:rank)]
			else if(jth.eq.rank)then
				Res%Dimension=[dim(1:ith-1),product(dim(ith:jth))]
			else
				call writemess('ERROR in fusing legs',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call error_stop
			end if
			Res%Data=T%Data
			if(.not.NameFlag)return

			ii=0
			do i=1,ith-1
				ii=ii+1
				call Res%setName(ii,Names(i))
			end do

			ii=ith
			call Res%setName(ii,'Fuse.leg'+ii)
			do i=jth+1,rank
				ii=ii+1
				call Res%setName(ii,Names(i))
			end do
			return
		end if

		if(ith.eq.jth)then
			Res=T
			return
		end if
		if(fermiCheck)then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res)
		else
			call QuantumFuse(T,ith,jth,Res)
		end if
	end subroutine
	subroutine FuseLegs2Subroutine(T,ith,jth,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith,jth
		type(Tensor),intent(inout)::Res
		integer,pointer::dim(:)
		character(len=len_of_name),pointer::Names(:)
		integer::rank,i,ii
		logical::NameFlag
		integer,allocatable::TMPDim(:)
		character(len=len_of_name),allocatable::TMPNames(:)
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if

		if(.not.T%getSymmetryFlag())then
			if(ith.eq.jth)then
				Res=T
				return
			end if
			if(deallocate_memory_flag)then
				allocate(TMPDim(rank))
				allocate(TMPNames(rank))
				TMPDim=T%Dimension
				NameFlag=T%getNameFlag()
				if(NameFlag)then
					do i=1,rank
						TMPNames(i)=T%getName(i)
					end do
				end if
				if((ith.gt.1).and.(jth.lt.rank))then
					Res%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(ith.eq.1)then
					Res%Dimension=[product(TMPDim(ith:jth)),TMPDim(jth+1:rank)]
				else if(jth.eq.rank)then
					Res%Dimension=[TMPDim(1:ith-1),product(TMPDim(ith:jth))]
				else
					call writemess('ERROR in fusing legs',-1)
					call writemess('ith='+ith,-1)
					call writemess('jth='+jth,-1)
					call error_stop
				end if
				Res%Data=T%Data
				if(.not.NameFlag)return
				ii=ith
				call Res%setName(ii,'Fuse.leg'+ii)
				do i=jth+1,rank
					ii=ii+1
					call Res%setName(ii,TMPNames(i))
				end do
				return
			end if
			call T%pointDim(Dim)
			NameFlag=T%getNameFlag()
			if(NameFlag)call T%pointName(Names)
			if((ith.gt.1).and.(jth.lt.rank))then
				Res%Dimension=[dim(1:ith-1),product(dim(ith:jth)),dim(jth+1:rank)]
			else if(ith.eq.1)then
				Res%Dimension=[product(dim(ith:jth)),dim(jth+1:rank)]
			else if(jth.eq.rank)then
				Res%Dimension=[dim(1:ith-1),product(dim(ith:jth))]
			else
				call writemess('ERROR in fusing legs',-1)
				call writemess('ith='+ith,-1)
				call writemess('jth='+jth,-1)
				call error_stop
			end if
			Res%Data=T%Data
			if(.not.NameFlag)return
			ii=0
			do i=1,ith-1
				ii=ii+1
				call Res%setName(ii,Names(i))
			end do
			
			ii=ith
			call Res%setName(ii,'Fuse.leg'+ii)
			do i=jth+1,rank
				ii=ii+1
				call Res%setName(ii,Names(i))
			end do
			return
		end if

		if(ith.eq.jth)then
			Res=T
			return
		end if
		if(T%getFermiFlag())then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,.true.)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res)
		else
			call QuantumFuse(T,ith,jth,Res)
		end if
	end subroutine

	subroutine FuseLegs3Subroutine(T,ith,jth,SubDim,fermiCheck,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith,jth
		logical,intent(in)::fermiCheck
		type(Dimension),intent(inout)::SubDim
		type(Tensor),intent(inout)::Res
		if(.not.T%getSymmetryFlag())then
			call FuseLegs1Subroutine(T,ith,jth,.false.,Res)
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			return
		end if
		if(ith.eq.jth)then
			Res=T
			return
		end if
		if(fermiCheck)then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res,SubDim)
		else
			call QuantumFuse(T,ith,jth,Res,SubDim)
		end if
	end subroutine

	subroutine FuseLegs4Subroutine(T,ith,jth,SubDim,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith,jth
		type(Dimension),intent(inout)::SubDim
		type(Tensor),intent(inout)::Res
		if(.not.T%getSymmetryFlag())then
			call FuseLegs1Subroutine(T,ith,jth,.false.,Res)
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			return
		end if
		if(ith.eq.jth)then
			Res=T
			return
		end if
		if(T%getFermiFlag())then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,.true.)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res,SubDim)
		else
			call QuantumFuse(T,ith,jth,Res,SubDim)
		end if
	end subroutine


	subroutine FuseLegs5Subroutine(T,ith,jth,SubDim,ResOrder,fermiCheck,Res)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		type(Dimension),intent(inout)::SubDim
		logical,intent(in)::fermiCheck
		type(Tensor),intent(inout)::Res,ResOrder
		integer::rank

		if(.not.T%getSymmetryFlag())then
			call FuseLegs1Subroutine(T,ith,jth,.false.,Res)
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			call ResOrder%empty()
			return
		end if
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if
		if(ith.eq.jth)then
			Res=T
			call SubDim%empty()
			call ResOrder%empty()
			return
		end if

		if(fermiCheck)then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,fermiCheck)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res,SubDim,ResOrder)
		else
			call QuantumFuse(T,ith,jth,Res,SubDim,ResOrder)
		end if
		return
	end subroutine

	subroutine FuseLegs6Subroutine(T,ith,jth,SubDim,ResOrder,Res)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith,jth
		type(Dimension),intent(inout)::SubDim
		type(Tensor),intent(inout)::Res,ResOrder
		integer::rank

		if(.not.T%getSymmetryFlag())then
			call FuseLegs1Subroutine(T,ith,jth,.false.,Res)
			call getsubDimension(T%Dimension,[ith,jth],SubDim)
			call ResOrder%empty()
			return
		end if
		rank=T%getRank()
		if(jth.gt.rank)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('jth='+jth,-1)
			call writemess('rank='+rank,-1)
			call error_stop
		end if
		if(ith.le.0)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		if(ith.gt.jth)then
			call writemess('ERROR in fusing legs',-1)
			call writemess('ith='+ith,-1)
			call writemess('jth='+jth,-1)
			call error_stop
		end if
		if(ith.eq.jth)then
			Res=T
			call SubDim%empty()
			call ResOrder%empty()
			return
		end if

		if(T%getFermiFlag())then
			TMPFuseTensor=T
			call Fuse_for_fermi_subroutine(TMPFuseTensor,ith,jth,.true.)
			call QuantumFuse(TMPFuseTensor,ith,jth,Res,SubDim,ResOrder)
		else
			call QuantumFuse(T,ith,jth,Res,SubDim,ResOrder)
		end if
		return
	end subroutine


	subroutine SplitTensor1Subroutine(T,ith,SubDim,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith
		type(Dimension)::subDim
		type(Tensor),intent(inout)::Res
		if(T%getSymmetryFlag())then
			call writemess('ERROR in split tensor, DO NOT input the FuseOrder for the symmetry tensor to split',-1)
			call error_stop
		end if
		call insertDimension(T%Dimension,SubDim,ith,Res%Dimension)
		Res%Data=T%Data
		return
	end subroutine

	subroutine SplitTensor1(T,ith,SubDim)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith
		type(Dimension)::subDim
		if(T%getSymmetryFlag())then
			call writemess('ERROR in split tensor, DO NOT input the FuseOrder for the symmetry tensor to split',-1)
			call error_stop
		end if
		if(.not.SubDim%getDimFlag()) return
		TMPSplitDimension=T%Dimension
		call insertDimension(TMPSplitDimension,SubDim,ith,T%Dimension)
		return
	end subroutine

	subroutine SplitTensor2Subroutine(T,ith,SubDim,ResOrder,fermiCheck,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith
		type(Dimension),intent(in)::SubDim
		logical,intent(in)::fermiCheck
		type(Tensor),intent(in)::ResOrder
		type(Tensor),intent(inout)::Res
		if(.not.T%getSymmetryFlag())then
			call SplitTensor1Subroutine(T,ith,SubDim,Res)
			return
		end if
		if(.not.SubDim%getDimFlag())then
			Res=T
			return
		end if
		call QuantumSplit(T,ith,Res,SubDim,ResOrder)
		if(fermiCheck)then
			call Fuse_for_fermi_subroutine(Res,ith,ith+SubDim%getRank()-1,fermiCheck)
		end if
		return
	end subroutine

	subroutine SplitTensor3Subroutine(T,ith,SubDim,ResOrder,Res)
		class(Tensor),intent(in)::T
		integer,intent(in)::ith
		type(Dimension),intent(in)::SubDim
		type(Tensor),intent(in)::ResOrder
		type(Tensor),intent(inout)::Res
		if(.not.T%getSymmetryFlag())then
			call SplitTensor1Subroutine(T,ith,SubDim,Res)
			return
		end if
		if(.not.SubDim%getDimFlag())then
			Res=T
			return
		end if
		call QuantumSplit(T,ith,Res,SubDim,ResOrder)
		if(T%getFermiFlag())then
			call Fuse_for_fermi_subroutine(Res,ith,ith+SubDim%getRank()-1,.true.)
		end if
		return
	end subroutine

	subroutine SplitTensor2(T,ith,SubDim,ResOrder,fermiCheck)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith
		type(Dimension),intent(in)::SubDim
		logical,intent(in)::fermiCheck
		type(Tensor),intent(in)::ResOrder
		if(.not.T%getSymmetryFlag())then
			call SplitTensor1(T,ith,SubDim)
			return
		end if
		if(.not.SubDim%getDimFlag()) return
			
		TMPSplitTensor=T
		call QuantumSplit(TMPSplitTensor,ith,T,SubDim,ResOrder)
		if(fermiCheck)then
			call Fuse_for_fermi_subroutine(T,ith,ith+SubDim%getRank()-1,fermiCheck)
		end if
		return
	end subroutine

	subroutine SplitTensor3(T,ith,SubDim,ResOrder)
		class(Tensor),intent(inout)::T
		integer,intent(in)::ith
		type(Dimension),intent(in)::SubDim
		type(Tensor),intent(in)::ResOrder
		if(.not.T%getSymmetryFlag())then
			call SplitTensor1(T,ith,SubDim)
			return
		end if
		if(.not.SubDim%getDimFlag()) return
			
		TMPSplitTensor=T
		call QuantumSplit(TMPSplitTensor,ith,T,SubDim,ResOrder)
		if(T%getFermiFlag())then
			call Fuse_for_fermi_subroutine(T,ith,ith+SubDim%getRank()-1,.true.)
		end if
		return
	end subroutine

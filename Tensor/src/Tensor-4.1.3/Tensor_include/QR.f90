#define QRmatrixFuncName QRmatrixs
#define QRSubroutineName QRSubroutines
#define ifNonZeroBlockDATATYPE ifNonZeroBlocks
#define DATATYPE real*4
#include "templet/QR0.f90"
#undef QRmatrixFuncName
#undef QRSubroutineName
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define QRmatrixFuncName QRmatrixd
#define QRSubroutineName QRSubroutined
#define ifNonZeroBlockDATATYPE ifNonZeroBlockd
#define DATATYPE real*8
#include "templet/QR0.f90"
#undef QRmatrixFuncName
#undef QRSubroutineName
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define QRmatrixFuncName QRmatrixc
#define QRSubroutineName QRSubroutinec
#define ifNonZeroBlockDATATYPE ifNonZeroBlockc
#define DATATYPE complex*8
#include "templet/QR0.f90"
#undef QRmatrixFuncName
#undef QRSubroutineName
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define QRmatrixFuncName QRmatrixz
#define QRSubroutineName QRSubroutinez
#define ifNonZeroBlockDATATYPE ifNonZeroBlockz
#define DATATYPE complex*16
#include "templet/QR0.f90"
#undef QRmatrixFuncName
#undef QRSubroutineName
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE


	subroutine QRSubroutine(A,Q,R)
		class(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::Q,R
		select case(A%getType())
			case(2)
				call QRSubroutines(A,Q,R)
			case(3)
				call QRSubroutined(A,Q,R)
			case(4)
				call QRSubroutinec(A,Q,R)
			case(5)
				call QRSubroutinez(A,Q,R)
			case default
				call writemess('ERROR type in QR',-1)
				call writemess('A%getType()='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine


	subroutine QRSubroutineName(A,Q,R,nameU,nameV)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::Q,R
		character(len=*),intent(in)::nameU,nameV
		integer::rankU,rankV,i,j,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		rank=A%getRank()
		call allocateCheck(SVDPermute,rank+rank)
		forwardi=>SVDPermute(1:rank)
		backwarki=>SVDPermute(rank+1:rank+rank)
		call A%pointName(AName)
		rankU=0
		rankV=0
		do i=1,rank
			if((AName(i).subl.indexsymbol).equ.nameU)then
				rankU=rankU+1
				forwardi(rankU)=i
			end if
			if((AName(i).subl.indexsymbol).equ.nameV)then
				rankV=rankV+1
				backwarki(rankV)=i
			end if
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in QRSubroutineName",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in QRSubroutineName,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in QRSubroutineName,no such name",-1)
			call writemess(nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call A%forward(forwardi(1:rankU),TMPQR1)
		else
			call A%backward(backwarki(1:rankV),TMPQR1)
		end if

		call TMPQR1%FuseTensor(rank-rankV+1,rank,dimen(2),order(2),.false.,TMPQR2)
		call TMPQR2%FuseTensor(1,rankU,dimen(1),order(1),.false.,TMPQR1)
		if(A%getSymmetryFlag())then
			len_dim2=TMPQR1%dim(2)
			call allocateCheck(reorder_in_QR,len_dim2)
			call TMPQR1%reOrderToDiag(reorder_in_QR(1:len_dim2),ifDiagFlag)


			call QRSubroutine(TMPQR1,TMPQR2,TMPQR3)

			call TMPQR3%reOrderTensor(reorder_in_QR(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPQR2%setName(2,QR_Q_leg)
				call TMPQR3%setName(1,QR_R_leg)
			end if
			call TMPQR2%SplitTensor(1,dimen(1),order(1),.false.,Q)
			call TMPQR3%SplitTensor(2,dimen(2),order(2),.false.,R)
			if(A%getFermiFlag())then
				call Q%setFermiArrow(Q%getRank(),-1)
				call R%setFermiArrow(1,1)
			end if
		else
			call QRSubroutine(TMPQR1,Q,R)
			if(A%getNameFlag())then
				call Q%setName(2,QR_Q_leg)
				call R%setName(1,QR_R_leg)
			end if
			call Q%Split(1,dimen(1),order(1),.false.)
			call R%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine

	subroutine QRSubroutineNumLegs(A,Q,R,legs,LeftFlag_)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::Q,R
		character(len=*),intent(in)::legs(:)
		character(len=*),optional,intent(in)::LeftFlag_
		character(len=15)::LeftFlag
		integer::lenlegs,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		if(present(LeftFlag_))then
			LeftFlag=LeftFlag_
		else
			LeftFlag='left'
		end if
		rank=A%getRank()
		lenlegs=size(legs)
		if((LeftFlag.eq.'left').or.(LeftFlag.eq.'row').or.(LeftFlag.eq.'r'))then
			call A%forward(legs,TMPQR1)
			call TMPQR1%FuseTensor(lenlegs+1,rank,dimen(2),order(2),.false.,TMPQR2)
			call TMPQR2%FuseTensor(1,lenlegs,dimen(1),order(1),.false.,TMPQR1)
		else if((LeftFlag.eq.'right').or.(LeftFlag.eq.'col').or.(LeftFlag.eq.'c'))then
			call A%backward(legs,TMPQR1)
			call TMPQR1%FuseTensor(rank-lenlegs+1,rank,dimen(2),order(2),.false.,TMPQR2)
			call TMPQR2%FuseTensor(1,rank-lenlegs,dimen(1),order(1),.false.,TMPQR1)
		else
			call writemess('ERROR in QRSubroutineNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then
			len_dim2=TMPQR1%dim(2)
			call allocateCheck(reorder_in_QR,len_dim2)
			call TMPQR1%reOrderToDiag(reorder_in_QR(1:len_dim2),ifDiagFlag)

			call QRSubroutine(TMPQR1,TMPQR2,TMPQR3)

			call TMPQR3%reOrderTensor(reorder_in_QR(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPQR2%setName(2,QR_Q_leg)
				call TMPQR3%setName(1,QR_R_leg)
			end if
			call TMPQR2%SplitTensor(1,dimen(1),order(1),.false.,Q)
			call TMPQR3%SplitTensor(2,dimen(2),order(2),.false.,R)
			if(A%getFermiFlag())then
				call Q%setFermiArrow(Q%getRank(),-1)
				call R%setFermiArrow(1,1)
			end if
		else
			call QRSubroutine(TMPQR1,Q,R)
			if(A%getNameFlag())then
				call Q%setName(2,QR_Q_leg)
				call R%setName(1,QR_R_leg)
			end if
			call Q%Split(1,dimen(1),order(1),.false.)
			call R%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine

	subroutine QRMatrix(A,Q,R)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Q,R
		TMPQR1=A
		select case(A%getType())
			case(2)
				call QRSubroutines(TMPQR1,Q,R)
			case(3)
				call QRSubroutined(TMPQR1,Q,R)
			case(4)
				call QRSubroutinec(TMPQR1,Q,R)
			case(5)
				call QRSubroutinez(TMPQR1,Q,R)
			case default
				call writemess('ERROR type in QR',-1)
				call writemess('A%getType()='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine
#define LQmatrixFuncName LQmatrixs
#define LQSubroutineFUNCNAME LQSubroutines
#define ifNonZeroBlockDATATYPE ifNonZeroBlocks
#define DATATYPE real*4
#include "templet/LQ0.f90"
#undef LQmatrixFuncName
#undef LQSubroutineFUNCNAME
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define LQmatrixFuncName LQmatrixd
#define LQSubroutineFUNCNAME LQSubroutined
#define ifNonZeroBlockDATATYPE ifNonZeroBlockd
#define DATATYPE real*8
#include "templet/LQ0.f90"
#undef LQmatrixFuncName
#undef LQSubroutineFUNCNAME
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define LQmatrixFuncName LQmatrixc
#define LQSubroutineFUNCNAME LQSubroutinec
#define ifNonZeroBlockDATATYPE ifNonZeroBlockc
#define DATATYPE complex*8
#include "templet/LQ0.f90"
#undef LQmatrixFuncName
#undef LQSubroutineFUNCNAME
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE

#define LQmatrixFuncName LQmatrixz
#define LQSubroutineFUNCNAME LQSubroutinez
#define ifNonZeroBlockDATATYPE ifNonZeroBlockz
#define DATATYPE complex*16
#include "templet/LQ0.f90"
#undef LQmatrixFuncName
#undef LQSubroutineFUNCNAME
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE


	subroutine LQSubroutine(A,L,Q)
		class(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::L,Q
		select case(A%getType())
			case(2)
				call LQSubroutines(A,L,Q)
			case(3)
				call LQSubroutined(A,L,Q)
			case(4)
				call LQSubroutinec(A,L,Q)
			case(5)
				call LQSubroutinez(A,L,Q)
			case default
				call writemess('ERROR type in LQ',-1)
				call writemess('A%getType()='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine LQSubroutineName(A,L,Q,nameU,nameV)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::L,Q
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
			call writemess("ERROR in LQSubroutineName",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call writemess('nameU='+nameU,-1)
			call writemess('nameV='+nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in LQSubroutineName,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in LQSubroutineName,no such name",-1)
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
			call allocateCheck(reorder_in_LQ,len_dim2)
			call TMPQR1%reOrderToDiag(reorder_in_LQ(1:len_dim2),ifDiagFlag)


			call LQSubroutine(TMPQR1,TMPQR2,TMPQR3)

			call TMPQR3%reOrderTensor(reorder_in_LQ(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPQR2%setName(2,LQ_L_leg)
				call TMPQR3%setName(1,LQ_Q_leg)
			end if
			call TMPQR2%SplitTensor(1,dimen(1),order(1),.false.,L)
			call TMPQR3%SplitTensor(2,dimen(2),order(2),.false.,Q)
			if(A%getFermiFlag())then
				call L%setFermiArrow(L%getRank(),-1)
				call Q%setFermiArrow(1,1)
			end if
		else
			call LQSubroutine(TMPQR1,L,Q)
			if(A%getNameFlag())then
				call L%setName(2,LQ_L_leg)
				call Q%setName(1,LQ_Q_leg)
			end if
			call L%Split(1,dimen(1),order(1),.false.)
			call Q%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine

	subroutine LQSubroutineNumLegs(A,L,Q,legs,LeftFlag_)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::L,Q
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
			call writemess('ERROR in LQSubroutineNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then

			len_dim2=TMPQR1%dim(2)
			call allocateCheck(reorder_in_LQ,len_dim2)
			call TMPQR1%reOrderToDiag(reorder_in_LQ(1:len_dim2),ifDiagFlag)



			call LQSubroutine(TMPQR1,TMPQR2,TMPQR3)

			call TMPQR3%reOrderTensor(reorder_in_LQ(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPQR2%setName(2,LQ_L_leg)
				call TMPQR3%setName(1,LQ_Q_leg)
			end if
			call TMPQR2%SplitTensor(1,dimen(1),order(1),.false.,L)
			call TMPQR3%SplitTensor(2,dimen(2),order(2),.false.,Q)
			if(A%getFermiFlag())then
				call L%setFermiArrow(L%getRank(),-1)
				call Q%setFermiArrow(1,1)
			end if
		else
			call LQSubroutine(TMPQR1,L,Q)
			if(A%getNameFlag())then
				call L%setName(2,LQ_L_leg)
				call Q%setName(1,LQ_Q_leg)
			end if
			call L%Split(1,dimen(1),order(1),.false.)
			call Q%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine

	subroutine LQmatrix(A,L,Q)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::L,Q
		TMPQR1=A
		select case(A%getType())
			case(2)
				call LQSubroutines(TMPQR1,L,Q)
			case(3)
				call LQSubroutined(TMPQR1,L,Q)
			case(4)
				call LQSubroutinec(TMPQR1,L,Q)
			case(5)
				call LQSubroutinez(TMPQR1,L,Q)
			case default
				call writemess('ERROR type in LQ',-1)
				call writemess('A%getType()='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine
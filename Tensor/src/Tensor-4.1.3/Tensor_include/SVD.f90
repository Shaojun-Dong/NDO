	subroutine choosSameQN(outQN,Dimen)
		type(Dimension),intent(in)::Dimen
		type(QuanNum),intent(inout)::outQN(2)
		real*4::Qi
		real*4,pointer::Q1(:),QN1(:),QN2(:),Q2(:)
		integer,pointer::Alldeg(:),Deg1(:),Deg2(:),dim(:)
		integer::maxlen,lenQ,i,j
		logical::Flag
		call Dimen%pointDim(Dim)
		maxlen=maxval(Dim)
		call WorkingMemory%allocate(2,maxlen+maxlen)
		call WorkingMemory%allocate(1,maxlen)
		call WorkingMemory%get_memory(Q1,maxlen)
		call WorkingMemory%get_memory(Q2,maxlen)
		call WorkingMemory%get_memory(Alldeg,maxlen)
		lenQ=0
		call Dimen%pointQN(QN1,1)
		call Dimen%pointQN(QN2,2)
		call Dimen%pointDeg(Deg1,1)
		call Dimen%pointDeg(Deg2,2)
		do i=1,Dim(1)
			do j=1,Dim(2)
				Flag = if_symmetry_Rule(dimen,[i,j])
				if(Flag)then
					lenQ=lenQ+1
					Q1(lenQ)=QN2(j)
					Q2(lenQ)=QN1(i)
					Alldeg(lenQ)=min(Deg1(i),Deg2(j))
					exit
				end if
			end do
		end do
		if(lenQ.eq.0)then
			call writemess(' ERROR in choosSameQN, there is no Quantum numbers that are the same')
			call Dimen%diminfo(.true.)
			call error_stop
		end if
		call outQN(1)%setQN(Q1(1:lenQ))
		call outQN(1)%setDeg(Alldeg(1:lenQ))
		call outQN(2)%setQN(Q2(1:lenQ))
		call outQN(2)%setDeg(Alldeg(1:lenQ))
		call WorkingMemory%free()
		return
	end subroutine

#define SVDDataSubroutineDATATYPE SVDDataSubroutines
#define SVDSavingSingleValueDegDATATYPE SVDSavingSingleValueDegs
#define SVDSavingSingleValueDegValueDATATYPE SVDSavingSingleValueDegValues
#define SVDSavingSingleValueNonSymmetryDATATYPE SVDSavingSingleValueNonSymmetrys
#define reorderDataDATATYPE reorderDatas
#define SVDMatrixKillNumDATATYPE SVDMatrixKillNums
#define SVDMatrixKillNum2DATATYPE SVDMatrixKillNum2s
#define SVDMatrixKillValueDATATYPE SVDMatrixKillValues
#define SVDMatrixKillValue2DATATYPE SVDMatrixKillValue2s
#define ifNonZeroBlockDATATYPE ifNonZeroBlocks
#define DATATYPE real*4
#define DATATYPE2 real*4
#include "templet/SVD0.f90"
#undef SVDDataSubroutineDATATYPE
#undef SVDSavingSingleValueDegDATATYPE
#undef SVDSavingSingleValueDegValueDATATYPE
#undef SVDSavingSingleValueNonSymmetryDATATYPE
#undef reorderDataDATATYPE
#undef SVDMatrixKillNumDATATYPE
#undef SVDMatrixKillNum2DATATYPE
#undef SVDMatrixKillValueDATATYPE
#undef SVDMatrixKillValue2DATATYPE
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE
#undef DATATYPE2

#define SVDDataSubroutineDATATYPE SVDDataSubroutined
#define SVDSavingSingleValueDegDATATYPE SVDSavingSingleValueDegd
#define SVDSavingSingleValueDegValueDATATYPE SVDSavingSingleValueDegValued
#define SVDSavingSingleValueNonSymmetryDATATYPE SVDSavingSingleValueNonSymmetryd
#define reorderDataDATATYPE reorderDatad
#define SVDMatrixKillNumDATATYPE SVDMatrixKillNumd
#define SVDMatrixKillNum2DATATYPE SVDMatrixKillNum2d
#define SVDMatrixKillValueDATATYPE SVDMatrixKillValued
#define SVDMatrixKillValue2DATATYPE SVDMatrixKillValue2d
#define ifNonZeroBlockDATATYPE ifNonZeroBlockd
#define DATATYPE real*8
#define DATATYPE2 real*8
#include "templet/SVD0.f90"
#undef SVDDataSubroutineDATATYPE
#undef SVDSavingSingleValueDegDATATYPE
#undef SVDSavingSingleValueDegValueDATATYPE
#undef SVDSavingSingleValueNonSymmetryDATATYPE
#undef reorderDataDATATYPE
#undef SVDMatrixKillNumDATATYPE
#undef SVDMatrixKillNum2DATATYPE
#undef SVDMatrixKillValueDATATYPE
#undef SVDMatrixKillValue2DATATYPE
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE
#undef DATATYPE2

#define SVDDataSubroutineDATATYPE SVDDataSubroutinec
#define SVDSavingSingleValueDegDATATYPE SVDSavingSingleValueDegc
#define SVDSavingSingleValueDegValueDATATYPE SVDSavingSingleValueDegValuec
#define SVDSavingSingleValueNonSymmetryDATATYPE SVDSavingSingleValueNonSymmetryc
#define reorderDataDATATYPE reorderDatac
#define SVDMatrixKillNumDATATYPE SVDMatrixKillNumc
#define SVDMatrixKillNum2DATATYPE SVDMatrixKillNum2c
#define SVDMatrixKillValueDATATYPE SVDMatrixKillValuec
#define SVDMatrixKillValue2DATATYPE SVDMatrixKillValue2c
#define ifNonZeroBlockDATATYPE ifNonZeroBlockc
#define DATATYPE complex*8
#define DATATYPE2 real*4
#include "templet/SVD0.f90"
#undef SVDDataSubroutineDATATYPE
#undef SVDSavingSingleValueDegDATATYPE
#undef SVDSavingSingleValueDegValueDATATYPE
#undef SVDSavingSingleValueNonSymmetryDATATYPE
#undef reorderDataDATATYPE
#undef SVDMatrixKillNumDATATYPE
#undef SVDMatrixKillNum2DATATYPE
#undef SVDMatrixKillValueDATATYPE
#undef SVDMatrixKillValue2DATATYPE
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE
#undef DATATYPE2

#define SVDDataSubroutineDATATYPE SVDDataSubroutinez
#define SVDSavingSingleValueDegDATATYPE SVDSavingSingleValueDegz
#define SVDSavingSingleValueDegValueDATATYPE SVDSavingSingleValueDegValuez
#define SVDSavingSingleValueNonSymmetryDATATYPE SVDSavingSingleValueNonSymmetryz
#define reorderDataDATATYPE reorderDataz
#define SVDMatrixKillNumDATATYPE SVDMatrixKillNumz
#define SVDMatrixKillNum2DATATYPE SVDMatrixKillNum2z
#define SVDMatrixKillValueDATATYPE SVDMatrixKillValuez
#define SVDMatrixKillValue2DATATYPE SVDMatrixKillValue2z
#define ifNonZeroBlockDATATYPE ifNonZeroBlockz
#define DATATYPE complex*16
#define DATATYPE2 real*8
#include "templet/SVD0.f90"
#undef SVDDataSubroutineDATATYPE
#undef SVDSavingSingleValueDegDATATYPE
#undef SVDSavingSingleValueDegValueDATATYPE
#undef SVDSavingSingleValueNonSymmetryDATATYPE
#undef reorderDataDATATYPE
#undef SVDMatrixKillNumDATATYPE
#undef SVDMatrixKillNum2DATATYPE
#undef SVDMatrixKillValueDATATYPE
#undef SVDMatrixKillValue2DATATYPE
#undef ifNonZeroBlockDATATYPE
#undef DATATYPE
#undef DATATYPE2


	subroutine SVDMatrixKillNum(A,inoutU,inoutS,inoutV,CUTOFF)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV
		integer,optional,intent(in)::CUTOFF
		select case(A%getType())
			case(2)
				call SVDMatrixKillNums(A,inoutU,inoutS,inoutV,CUTOFF)
			case(3)
				call SVDMatrixKillNumd(A,inoutU,inoutS,inoutV,CUTOFF)
			case(4)
				call SVDMatrixKillNumc(A,inoutU,inoutS,inoutV,CUTOFF)
			case(5)
				call SVDMatrixKillNumz(A,inoutU,inoutS,inoutV,CUTOFF)
			case default
				call writemess('ERROR in SVD, data type',-1)
				call writemess('A%getType()='+A%getType())
				call error_stop
		end select
		return
	end subroutine

	subroutine SVDMatrixKillNum2(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV,WT1,WT2,WT3
		integer,optional,intent(in)::CUTOFF
		select case(A%getType())
			case(2)
				call SVDMatrixKillNum2s(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
			case(3)
				call SVDMatrixKillNum2d(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
			case(4)
				call SVDMatrixKillNum2c(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
			case(5)
				call SVDMatrixKillNum2z(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
			case default
				call writemess('ERROR in SVD, data type',-1)
				call writemess('A%getType()='+A%getType())
				call error_stop
		end select
		return
	end subroutine

	subroutine SVDMatrixKillValue(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		select case(A%getType())
			case(2)
				call SVDMatrixKillValues(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(3)
				call SVDMatrixKillValued(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(4)
				call SVDMatrixKillValuec(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(5)
				call SVDMatrixKillValuez(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case default
				call writemess('ERROR in SVD, data type',-1)
				call writemess('A%getType()='+A%getType())
				call error_stop
		end select
		return
	end subroutine

	subroutine SVDMatrixKillValue2(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV,WT1,WT2,WT3
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		select case(A%getType())
			case(2)
				call SVDMatrixKillValue2s(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(3)
				call SVDMatrixKillValue2d(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(4)
				call SVDMatrixKillValue2c(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case(5)
				call SVDMatrixKillValue2z(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			case default
				call writemess('ERROR in SVD, data type',-1)
				call writemess('A%getType()='+A%getType())
				call error_stop
		end select
		return
	end subroutine


	subroutine SVDNumName(A,U,S,V,nameU,nameV,CUTOFF)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,optional,intent(in)::CUTOFF
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
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call A%forward(forwardi(1:rankU),TMPSVD1)
		else
			call A%backward(backwarki(1:rankV),TMPSVD1)
		end if

		call TMPSVD1%FuseTensor(rank-rankV+1,rank,dimen(2),order(2),.false.,TMPSVD2)
		call TMPSVD2%FuseTensor(1,rankU,dimen(1),order(1),.false.,TMPSVD1)
		if(A%getSymmetryFlag())then

			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillNum2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,CUTOFF)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillNum(TMPSVD1,U,S,V,CUTOFF)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine
	subroutine SVDNumNameMatrixS(A,U,S,V,nameU,nameV,MatrixS,CUTOFF)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		logical,intent(in)::MatrixS
		integer,optional,intent(in)::CUTOFF
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
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call A%forward(forwardi(1:rankU),TMPSVD1)
		else
			call A%backward(backwarki(1:rankV),TMPSVD1)
		end if

		call TMPSVD1%FuseTensor(rank-rankV+1,rank,dimen(2),order(2),.false.,TMPSVD2)
		call TMPSVD2%FuseTensor(1,rankU,dimen(1),order(1),.false.,TMPSVD1)
		if(A%getSymmetryFlag())then

			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillNum2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,CUTOFF)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillNum(TMPSVD1,U,S,V,CUTOFF)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		if(MatrixS)then
			if(A%getSymmetryFlag())then
				call S%eye(V%getQuantumNumber(1),.true.)
			else
				call S%eye(V%dim(1),V%dim(1))
			end if
			call S%setName(1,SVD_S_leg1)
			call S%setName(2,SVD_S_leg2)
		end if
		return
	end subroutine

	subroutine SVDValueName(A,U,S,V,nameU,nameV,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
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
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call A%forward(forwardi(1:rankU),TMPSVD1)
		else
			call A%backward(backwarki(1:rankV),TMPSVD1)
		end if

		call TMPSVD1%FuseTensor(rank-rankV+1,rank,dimen(2),order(2),.false.,TMPSVD2)
		call TMPSVD2%FuseTensor(1,rankU,dimen(1),order(1),.false.,TMPSVD1)
		if(A%getSymmetryFlag())then


			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillValue2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillValue(TMPSVD1,U,S,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine
	subroutine SVDValueNameMatrixS(A,U,S,V,nameU,nameV,minNumSave,maxNumSave,inmaxValue,VType,&
		MatrixS,outNumSave,outmaxValue)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		logical,intent(in)::MatrixS
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
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
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call A%diminfo(.true.)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call A%forward(forwardi(1:rankU),TMPSVD1)
		else
			call A%backward(backwarki(1:rankV),TMPSVD1)
		end if

		call TMPSVD1%FuseTensor(rank-rankV+1,rank,dimen(2),order(2),.false.,TMPSVD2)
		call TMPSVD2%FuseTensor(1,rankU,dimen(1),order(1),.false.,TMPSVD1)
		if(A%getSymmetryFlag())then


			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillValue2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillValue(TMPSVD1,U,S,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		if(MatrixS)then
			if(A%getSymmetryFlag())then
				call S%eye(V%getQuantumNumber(1),.true.)
			else
				call S%eye(V%dim(1),V%dim(1))
			end if
			call S%setName(1,SVD_S_leg1)
			call S%setName(2,SVD_S_leg2)
		end if
		return
	end subroutine


	subroutine SVDNumLegs(A,U,S,V,legs,LeftFlag,CUTOFF)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,optional,intent(in)::CUTOFF
		character(len=*),intent(in)::legs(:)
		character(len=*),intent(in)::LeftFlag
		integer::lenlegs,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		rank=A%getRank()
		lenlegs=size(legs)
		if((LeftFlag.eq.'left').or.(LeftFlag.eq.'row').or.(LeftFlag.eq.'r'))then
			call A%forward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else if((LeftFlag.eq.'right').or.(LeftFlag.eq.'col').or.(LeftFlag.eq.'c'))then
			call A%backward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(rank-lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,rank-lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else
			call writemess('ERROR in SVDMatrixKillNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then
			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)



			call SVDMatrixKillNum2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,CUTOFF)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillNum(TMPSVD1,U,S,V,CUTOFF)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine
	subroutine SVDNumLegsMatrixS(A,U,S,V,legs,LeftFlag,MatrixS,CUTOFF)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		logical,intent(in)::MatrixS
		integer,optional,intent(in)::CUTOFF
		character(len=*),intent(in)::legs(:)
		character(len=*),intent(in)::LeftFlag
		integer::lenlegs,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		rank=A%getRank()
		lenlegs=size(legs)
		if((LeftFlag.eq.'left').or.(LeftFlag.eq.'row').or.(LeftFlag.eq.'r'))then
			call A%forward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else if((LeftFlag.eq.'right').or.(LeftFlag.eq.'col').or.(LeftFlag.eq.'c'))then
			call A%backward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(rank-lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,rank-lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else
			call writemess('ERROR in SVDMatrixKillNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then
			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)



			call SVDMatrixKillNum2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,CUTOFF)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillNum(TMPSVD1,U,S,V,CUTOFF)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		if(MatrixS)then
			if(A%getSymmetryFlag())then
				call S%eye(V%getQuantumNumber(1),.true.)
			else
				call S%eye(V%dim(1),V%dim(1))
			end if
			call S%setName(1,SVD_S_leg1)
			call S%setName(2,SVD_S_leg2)
		end if
		return
	end subroutine

	subroutine SVDValueLegs(A,U,S,V,legs,LeftFlag,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		character(len=*),intent(in)::legs(:)
		character(len=*),intent(in)::LeftFlag
		integer::lenlegs,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		rank=A%getRank()
		lenlegs=size(legs)
		if((LeftFlag.eq.'left').or.(LeftFlag.eq.'row').or.(LeftFlag.eq.'r'))then
			call A%forward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else if((LeftFlag.eq.'right').or.(LeftFlag.eq.'col').or.(LeftFlag.eq.'c'))then
			call A%backward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(rank-lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,rank-lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else
			call writemess('ERROR in SVDMatrixKillNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then
			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillValue2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillValue(TMPSVD1,U,S,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		return
	end subroutine
	subroutine SVDValueLegsMatrixS(A,U,S,V,legs,LeftFlag,minNumSave,maxNumSave,inmaxValue,VType&
		,MatrixS,outNumSave,outmaxValue)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::U,S,V
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		logical,intent(in)::MatrixS
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		character(len=*),intent(in)::legs(:)
		character(len=*),intent(in)::LeftFlag
		integer::lenlegs,rank,len_dim2
		type(Tensor)::order(2)
		type(Dimension)::dimen(2)
		logical::ifDiagFlag
		character(len=len_of_name),pointer::AName(:)
		integer,pointer::forwardi(:),backwarki(:)
		rank=A%getRank()
		lenlegs=size(legs)
		if((LeftFlag.eq.'left').or.(LeftFlag.eq.'row').or.(LeftFlag.eq.'r'))then
			call A%forward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else if((LeftFlag.eq.'right').or.(LeftFlag.eq.'col').or.(LeftFlag.eq.'c'))then
			call A%backward(legs,TMPSVD1)
			call TMPSVD1%FuseTensor(rank-lenlegs+1,rank,dimen(2),order(2),.false.,TMPSVD2)
			call TMPSVD2%FuseTensor(1,rank-lenlegs,dimen(1),order(1),.false.,TMPSVD1)
		else
			call writemess('ERROR in SVDMatrixKillNumLegs',-1)
			call writemess('LeftFlag='+LeftFlag)
			call error_stop
		end if

		
		if(A%getSymmetryFlag())then
			len_dim2=TMPSVD1%dim(2)
			call allocateCheck(reorder_in_SVD,len_dim2)
			call TMPSVD1%reOrderToDiag(reorder_in_SVD(1:len_dim2),ifDiagFlag)


			call SVDMatrixKillValue2(TMPSVD1,TMPSVD2,S,TMPSVD3,U,TMPSVD4,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)

			call TMPSVD3%reOrderTensor(reorder_in_SVD(1:len_dim2),ifDiagFlag)

			if(A%getNameFlag())then
				call TMPSVD2%setName(2,SVD_U_leg)
				call TMPSVD3%setName(1,SVD_V_leg)
			end if
			call TMPSVD2%SplitTensor(1,dimen(1),order(1),.false.,U)
			call TMPSVD3%SplitTensor(2,dimen(2),order(2),.false.,V)
			if(A%getFermiFlag())then
				call U%setFermiArrow(U%getRank(),-1)
				call V%setFermiArrow(1,1)
			end if
		else
			call SVDMatrixKillValue(TMPSVD1,U,S,V,minNumSave,maxNumSave,inmaxValue,VType,outNumSave,outmaxValue)
			if(A%getNameFlag())then
				call U%setName(2,SVD_U_leg)
				call V%setName(1,SVD_V_leg)
			end if
			call U%Split(1,dimen(1),order(1),.false.)
			call V%Split(2,dimen(2),order(2),.false.)
		end if
		if(MatrixS)then
			if(A%getSymmetryFlag())then
				call S%eye(V%getQuantumNumber(1),.true.)
			else
				call S%eye(V%dim(1),V%dim(1))
			end if
			call S%setName(1,SVD_S_leg1)
			call S%setName(2,SVD_S_leg2)
		end if
		return
	end subroutine

	subroutine SVDMatrixNum(A,inoutU,inoutS,inoutV,CUTOFF)
		class(Tensor),intent(in)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV
		integer,optional,intent(in)::CUTOFF
		TMPSVD1=A
		select case(A%getType())
			case(2)
				call SVDMatrixKillNums(TMPSVD1,inoutU,inoutS,inoutV,CUTOFF)
			case(3)
				call SVDMatrixKillNumd(TMPSVD1,inoutU,inoutS,inoutV,CUTOFF)
			case(4)
				call SVDMatrixKillNumc(TMPSVD1,inoutU,inoutS,inoutV,CUTOFF)
			case(5)
				call SVDMatrixKillNumz(TMPSVD1,inoutU,inoutS,inoutV,CUTOFF)
			case default
				call writemess('ERROR in SVD, data type',-1)
				call writemess('A%getType()='+A%getType())
				call error_stop
		end select
		return
	end subroutine
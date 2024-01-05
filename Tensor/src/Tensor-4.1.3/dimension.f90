!
!                   _ooOoo_
!                  o8888888o
!                  88" . "88
!                  (| -_- |)
!                  O\  =  /O
!               ____/`---'\____
!             .'  \\|     |//  `.
!            /  \\|||  :  |||//  \
!           /  _||||| -:- |||||-  \
!           |   | \\\  -  /// |   |
!           | \_|  ''\---/''  |   |
!           \  .-\__  `-`  ___/-. /
!         ___`. .'  /--.--\  `. . __
!      ."" '<  `.___\_<|>_/___.'  >'"".
!     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
!     \  \ `-.   \_ __\ /__ _/   .-` /  /
!======`-.____`-.___\_____/___.-`____.-'======
!                   `=---='
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Buddha blessed , no BUG 
! Report bugs of the package to sj.dong@outlook.com
!

module dimension_tools
	use data_base_Tools
	use QuantumNumber_Tools
	use Tools
	use mpi
	use memory_type
	implicit none
	private
	type(memory),private::WorkingMemory
	integer,private::characterlen_in_print=500
	character(len=len_of_Name),private,allocatable::WorkingMemory_Name(:)
	logical,private::old_version_read=.false.

	public::Dimension
	type Dimension
		private
		type(DataArray)::QuanNum
		type(DataArray)::Deg
		integer,allocatable::dimData(:)
		integer,allocatable::FermiArrow(:)
		character(len=len_of_Name),allocatable::DimName(:)
		logical::Flag(4)=.false.
		!Flag(1):dimension 
		!Flag(2):symmetry dimension: dimension,QuanNum,Deg
		!Flag(3):fermi dimension:    dimension,QuanNum,Deg,fermi-arrow
		!Flag(4):DimName
		integer::fixType=0
		!fixType=0, the dimension can store any data
		!fixType=1, the dimension can only be the non-symmetry and non-fermionic dimension
		!fixType=2, the dimension can only be the symmetry dimension
		!fixType=3, the dimension can only be the fermionic dimension
		integer::rank=0
	contains
		generic,public::deallocate=>deallocatableDimension
		procedure::emptyDimension,deallocatableDimension
		generic,public::empty=>emptyDimension
		procedure,public::GetNameFlag
		procedure::GetDimensionFlag
		generic,public::GetDimFlag=>GetDimensionFlag
		procedure,public::GetSymmetryFlag
		procedure,public::GetFermiFlag
		procedure,public::getRank
		procedure::setNameAll,setNameith,setNameAll2,setNameChaith
		generic,public::setName=>setNameAll,setNameith,setNameAll2,setNameChaith
		procedure::getNameIth,getNameCha,getNameAll
		generic,public::getName=>getNameIth,getNameCha,getNameAll
		procedure,public::outAllNameChar
		procedure::setQNith,setQNChar
		procedure::setDegith,setDegchar,setDegithAll,setDegAll,setDegCharAll
		generic,public::setQN=>setQNith,setQNChar
		generic,public::setDeg=>setDegith,setDegchar,setDegithAll,setDegAll,setDegCharAll
		procedure::setFermiArrowith,setFermiArrowAll,setFermiArrowChar
		generic,public::setFermiArrow=>setFermiArrowith,setFermiArrowAll,setFermiArrowChar
		procedure::FindOrderith,FindOrderAll
		procedure::NameOrderith,NameOrderAll
		generic,public::FindOrder=>FindOrderith,FindOrderAll
		generic,public::NameOrder=>NameOrderith,NameOrderAll
		procedure,public::pointDim
		procedure,public::pointArrow
		procedure,public::pointFermiArrow=>pointArrow
		procedure::pointQNith,pointQNChar
		generic,public::pointQN=>pointQNith,pointQNChar
		procedure::pointDegith,pointDegchar
		generic,public::pointDeg=>pointDegith,pointDegchar
		procedure,public::pointName

		procedure::writeExternlData,readExternalData
		procedure::printDataDetail,printDataSimple
		generic,public::diminfo=>printDataDetail,printDataSimple
		generic,public::write=>writeExternlData
		generic,public::read=>readExternalData


		procedure::getDimi,getDimchar,getDimAll
		generic,public::dim=>getDimi,getDimChar,getDimAll

		procedure::getBlockDimAll,getBlockDimith,getBlockDim_legi,getBlockDim_leg2
		generic,public::getBlockDim=>getBlockDimAll,getBlockDimith,getBlockDim_legi,getBlockDim_leg2

		procedure::getQNAll,getQNOne,getQNChaOne
		procedure::getDegAll,getDegOne,getDegChaOne,getDegAll_name
		procedure::getArrowAll,getArrowith,getArrowCha
		generic,public::getQN=>getQNAll,getQNOne,getQNChaOne
		generic,public::getDeg=>getDegAll,getDegOne,getDegChaOne,getDegAll_name
		generic,public::getArrow=>getArrowAll,getArrowith,getArrowCha
		generic,public::getFermiArrow=>getArrowAll,getArrowith,getArrowCha

		procedure::QN2nonSyminde_one,QN2nonSyminde_one2,QN2nonSyminde,QN2nonSyminde2
		generic,public::NonSymIndex=>QN2nonSyminde_one,QN2nonSyminde_one2,QN2nonSyminde,QN2nonSyminde2

		procedure,public::unsetFermiFlag
		generic,public::unsetFermiArrow=>unsetFermiFlag
		procedure,public::killZeroDeg
		procedure::getQuantumNumber1,getQuantumNumber2
		generic,public::getQuantumNumber=>getQuantumNumber1,getQuantumNumber2
		generic,public::QuantumNumber=>getQuantumNumber1,getQuantumNumber2
		procedure,public::ifName

		procedure::index2QNinfo1
		generic,public::index2QNinfo=>index2QNinfo1
		procedure::QNinfo2index1
		generic,public::QNinfo2index=>QNinfo2index1
		procedure::getRule1,getRule2,getAllRule
		generic,public::getRule=>getRule1,getRule2,getAllRule
		procedure,public::nonSymDimension
		procedure,public::getTotalVolume
		procedure,public::FixnonSymmetryType,FixNonFermionicType,FixSymmetryType,FixFermionicType
		procedure,public::unFixSymmetryType,unFixfermionicType,unFixnonSymmetryType,unFixNonFermionicType
	end type dimension

	public::assignment(=)
	interface assignment(=)
		module procedure ArrayToDim
		module procedure dimToArray
		module procedure DimToDim
		module procedure QuanNumToDim
		module procedure DimToQuanNum
	end interface

	public::paste
	interface paste
		module procedure pasteTwoDimension
	end interface

	public::operator(.subdim.)
	interface operator(.subdim.)
		module procedure subDimension
		module procedure subDimension2
		module procedure subDimension4
	end interface
	

	public::permutationDimension,pasteTwoSubDim,insertDimension,send_Dimension,BCAST_Dimension

	public::pasteDimension
	interface pasteDimension
		module procedure pasteTwoSubDim
		module procedure pasteDimQNDim
		module procedure pasteQNDim
		module procedure pasteDimQN	
		module procedure pastevecSubDim
		module procedure pasteSubDimvec
	end interface

	public::getsubDimension
	interface getsubDimension
		module procedure subDimensionSubroutine
	end interface

	INTERFACE
	  SUBROUTINE SetDimensionRuleRoutine(dimen,legi)
	  	import :: dimension
		Type(dimension),intent(inout)::dimen
			integer,optional,intent(in)::legi
	  END SUBROUTINE SetDimensionRuleRoutine
	END INTERFACE
	procedure(SetDimensionRuleRoutine),pointer::ReverDimensionRule=>defaultSetRule
	
	public::set_SetDimensionRule_function,set_old_version_read
contains
	subroutine set_old_version_read()
		old_version_read=.true.
	end subroutine

	subroutine error_mess()
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@                 ERROR                   @@@  ')
		call writemess('@@@ You have not set the symmetry group yet @@@  ')
		call writemess('@@@                                         @@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call error_stop
	end subroutine

	subroutine defaultSetRule(dimen,legi)
		Type(dimension),intent(inout)::dimen
			integer,optional,intent(in)::legi
		call error_mess()
	end subroutine
	subroutine set_SetDimensionRule_function(Func)
		procedure(SetDimensionRuleRoutine)::Func
		ReverDimensionRule=>Func
	end subroutine
	!*********************************************************
	!
	! assignment
	!
	!*********************************************************

	subroutine ArrayToDim(Dimen,DimData)
		type(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		call Dimen%empty()
		call check_type(dimen,1)
		Dimen%rank=size(DimData)
		call allocateCheck(Dimen%DimData,Dimen%rank)
		Dimen%DimData(1:Dimen%rank)=DimData
		Dimen%Flag(1)=.true.
		Dimen%Flag(2)=.false.
		Dimen%Flag(3)=.false.
		Dimen%Flag(4)=.false.
		return
	end subroutine
	subroutine dimToArray(intDimen,Dimen)
		integer,intent(inout) :: intDimen(:)
		type(Dimension),intent(in) ::Dimen
		if(.not.Dimen%getDimFlag())then
			call writemess('ERROR the dimension is empty',-1)
			call error_stop
		end if
		intDimen=Dimen%DimData(1:Dimen%rank)
		return
	end subroutine
	subroutine DimToDim(outDimen,inDimen)
		type(Dimension),intent(inout) ::outDimen
		type(Dimension),intent(in) ::inDimen
		integer,pointer::ip(:)
		integer::rank
		call outDimen%empty()
		if(.not.inDimen%GetDimFlag())return

		if(inDimen%getFermiFlag())then
			call check_type(outDimen,3)
		else if(inDimen%getSymmetryFlag())then
			call check_type(outDimen,2)
		else
			call check_type(outDimen,1)
		end if

		outDimen%Flag=inDimen%Flag
		outDimen%rank=inDimen%rank

		call allocateCheck(outDimen%DimData,outDimen%rank)
		outDimen%DimData(1:outDimen%rank)=inDimen%DimData(1:inDimen%rank)

		if(outDimen%Flag(2))then
			outDimen%QuanNum=inDimen%QuanNum
			outDimen%Deg=inDimen%Deg
		else
			call outDimen%QuanNum%empty()
			call outDimen%Deg%empty()
		end if
		if(outDimen%Flag(3))then
			call allocateCheck(outDimen%FermiArrow,outDimen%rank)
			outDimen%FermiArrow(1:outDimen%rank)=inDimen%FermiArrow(1:inDimen%rank)
		end if
		outDimen%fixType=inDimen%fixType
		if(inDimen%getNameFlag())then
			rank=inDimen%getRank()
			call allocateCheck(outDimen%DimName,rank)
			outDimen%DimName(1:rank)=inDimen%DimName(1:rank)
		end if
		return
	end subroutine
	subroutine QuanNumToDim(dimen,QN)
		type(Dimension),intent(inout)::dimen
		type(QuanNum),intent(in)::QN(:)
		integer,pointer::dim(:),arrow(:),ip(:)
		real*4,pointer::sp(:)
		integer::rank,i,FixFlag
		call Dimen%empty()
		FixFlag=2
		rank=size(QN)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank)
		call WorkingMemory%get_memory(dim,rank)
		call WorkingMemory%get_memory(arrow,rank)
		do i=1,rank
			dim(i)=QN(i)%getQNlength()
			arrow(i)=QN(i)%GetFermiArrow()
		end do

		Dimen%rank=size(dim)
		call allocateCheck(Dimen%DimData,Dimen%rank)
		Dimen%DimData(1:Dimen%rank)=dim

		call dimen%QuanNum%allocate(dim,'real*4')
		call dimen%Deg%allocate(dim,'integer')
		dimen%Flag(1)=.true.
		dimen%Flag(2)=.true.
		do i=1,rank
			call dimen%QuanNum%pointer(sp,i)
			call dimen%Deg%pointer(ip,i)
			sp=QN(i)%getQN()
			ip=QN(i)%getDeg()
			if(QN(i)%getRule().lt.0)then
				call ReverDimensionRule(Dimen,i)
			end if
		end do
		if(arrow(1).ne.0)then
			dimen%Flag(3)=.true.
			call allocateCheck(Dimen%FermiArrow,Dimen%rank)
			Dimen%FermiArrow(1:Dimen%rank)=arrow
			FixFlag=3
		else
			dimen%Flag(3)=.false.
		end if
		dimen%Flag(4)=.false.
		call WorkingMemory%free()
		call check_type(dimen,FixFlag)
		return
	end subroutine

	subroutine DimToQuanNum(QN,Dimen)
		type(QuanNum),intent(inout)::QN(:)
		type(Dimension),intent(in)::Dimen
		integer::i,rank
		integer,pointer::arrow(:),deg(:)
		real*4,pointer::sp(:)
		if(.not.Dimen%GetDimFlag())then
			do i=1,size(QN)
				call QN(i)%empty()
			end do
			return
		end if
		rank=Dimen%getRank()
		if(size(QN).lt.rank)then
			call writemess('ERROR in QN(:)=dimension',-1)
			call writemess('size(QN)='+size(QN),-1)
			call writemess('Dimen%getRank()='+Dimen%getRank(),-1)
			call error_Stop
		end if
		call Dimen%pointArrow(arrow)
		do i=1,rank
			call Dimen%pointQN(sp,i)
			call Dimen%pointDeg(deg,i)
			call QN(i)%setQN(sp)
			call QN(i)%setDeg(deg)
			call QN(i)%setRule(1)
			if(dimen%getFermiFlag())then
				call QN(i)%setFermiArrow(arrow(i)) 
			else
				call QN(i)%setFermiArrow(0)  
			end if
		end do
		return
	end subroutine

	subroutine allocateDimension(dimen,dimenData)
		class(Dimension),intent(inout)::dimen
		integer,intent(in)::dimenData(:)
		integer,pointer::iTMP(:)
		integer,pointer::ip(:)
		call dimen%QuanNum%allocate(dimenData,'real*4')
		call dimen%Deg%allocate(dimenData,'integer')
		dimen%rank=size(dimenData)
		call allocateCheck(dimen%dimData,dimen%rank)
		call allocateCheck(dimen%FermiArrow,dimen%rank)
		dimen%dimData(1:dimen%rank)=dimenData
		dimen%FermiArrow=0
		dimen%Flag(1)=.true.
		dimen%Flag(2:4)=.false.
		return
	end subroutine
	subroutine deallocatableDimension(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%Flag=.false.
		dimen%rank=0
		if(allocated(dimen%DimName))deallocate(dimen%DimName)
		if(allocated(dimen%dimData))deallocate(dimen%dimData)
		if(allocated(dimen%FermiArrow))deallocate(dimen%FermiArrow)
		call dimen%QuanNum%deallocate()
		call dimen%Deg%deallocate()
		return
	end subroutine
	subroutine emptyDimension(dimen)
		class(Dimension),intent(inout)::dimen
		if(deallocate_memory_flag)then
			call deallocatableDimension(dimen)
			return
		end if
		dimen%rank=0
		dimen%Flag=.false.
		call dimen%QuanNum%empty()
		call dimen%Deg%empty()
		return
	end subroutine

	!*********************************************************
	!
	!   print data
	!
	!*********************************************************

	subroutine printDataDetail(A,printdelta)
		class(Dimension),intent(in)::A
		logical,intent(in)::printdelta
		integer::i,j
		character(len=characterlen_in_print)::w
		integer,pointer::ip(:),ip2(:),arrowp(:),deg(:)
		character(len=len_of_Name),pointer::wp(:)
		real*4,pointer::sp(:)
		call writemess('==================================')
		if(.not.A%GetDimFlag())then
			call writemess(' The dimension is empty')
			return
		end if
		call A%pointDim(ip)
		w='dimension:'+(' '+ip(1))
		do i=2,A%getRank()
			w=w+' ,'+(' '+ip(i))
		end do
		call writemess(w)

		if(A%getSymmetryFlag())then
			call A%pointDeg(deg,1)
			w='non-Sym dimension:'+(' '+sum(deg))
			do i=2,A%getRank()
				call A%pointDeg(deg,i)
				w=w+' ,'+(' '+sum(deg))
			end do
			call writemess(w)
		end if

		if(A%getNameFlag())then
			call A%pointName(wp)
			w=' Name:'+(' '+wp(1))
			do i=2,A%getRank()
				w=w+' ,'+(' '+wp(i))
			end do
			call writemess(w)
		end if


		if(A%getFermiFlag())then
			call A%pointArrow(Arrowp)
			w='Arrow:'+(' '+Arrowp(1))
			do i=2,A%getRank()
				w=w+' ,'+(' '+Arrowp(i))
			end do
			call writemess(w)
		end if
		call writemess('Rank='+(' '+A%getRank()))

		if(printdelta.and.(A%getSymmetryFlag().or.A%getNameFlag()))then
			do i=1,A%getRank()
				call writemess(i)
				if(A%getNameFlag())then
					w='    Name:'+(' '+wp(i))
					call writemess(w) 
				end if

				w='     dim:'+(' '+ip(i))
				call writemess(w)

				if(A%getSymmetryFlag())then
					call A%pointQN(sp,i)
					call A%pointDeg(ip2,i)
					w='      QN:'+(' '+sp(1))
					do j=2,ip(i)
						w=w+','+(' '+sp(j))
					end do
					call writemess(w)
					w='     Deg:'+(' '+ip2(1))
					do j=2,ip(i)
						w=w+','+(' '+ip2(j))
					end do
					call writemess(w)
				end if

				if(A%getFermiFlag())then
					call A%pointArrow(Arrowp)
					call writemess('   arrow:'+(' '+Arrowp(i)))
				end if
				call writemess(' ')
			end do
		end if
		call writemess('--------------------------------')
		return
	end subroutine
	subroutine printDataSimple(A)
		class(Dimension),intent(in)::A
		integer::i,j
		character(len=characterlen_in_print)::w
		integer,pointer::ip(:),ip2(:),arrowp(:),deg(:)
		character(len=len_of_Name),pointer::wp(:)
		real*4,pointer::sp(:)
		call writemess('==================================')
		if(.not.A%GetDimFlag())then
			call writemess(' The dimension is empty')
			return
		end if
		call A%pointDim(ip)
		w='dimension:'+(' '+ip(1))
		do i=2,A%getRank()
			w=w+' ,'+(' '+ip(i))
		end do
		call writemess(w)

		if(A%getSymmetryFlag())then
			call A%pointDeg(deg,1)
			w='non-Sym dimension:'+(' '+sum(deg))
			do i=2,A%getRank()
				call A%pointDeg(deg,i)
				w=w+' ,'+(' '+sum(deg))
			end do
			call writemess(w)
		end if

		if(A%getNameFlag())then
			call A%pointName(wp)
			w=' Name:'+(' '+wp(1))
			do i=2,A%getRank()
				w=w+' ,'+(' '+wp(i))
			end do
			call writemess(w)
		end if


		if(A%getFermiFlag())then
			call A%pointArrow(Arrowp)
			w='Arrow:'+(' '+Arrowp(1))
			do i=2,A%getRank()
				w=w+' ,'+(' '+Arrowp(i))
			end do
			call writemess(w)
		end if
		call writemess('Rank='+(' '+A%getRank()))

		call writemess('--------------------------------')
		return
	end subroutine

	subroutine writeExternlData(A,uni)
		class(Dimension),intent(in)::A
		integer,intent(in)::uni
		integer::i
		if(old_version_read)then
			call writeExternlData_old(A,uni)
			return
		end if
		write(uni,*)A%Flag,A%rank
		if(.not.A%GetDimFlag())then
			write(uni,*)'Empty dimension'
			return
		end if
		call A%QuanNum%write(uni)
		call A%Deg%write(uni)
		write(uni,*)A%dimData(1:A%rank)
		if(A%Flag(3))then
			write(uni,*)A%FermiArrow(1:A%rank)
		end if
		if(A%getNameFlag())then
			write(uni,*)(trim(adjustl(A%DimName(i)))//"  ",i=1,A%GetRank())
		end if
		return
	end subroutine
	subroutine readExternalData(A,uni)
		class(Dimension),intent(inout)::A
		integer,intent(in)::uni
		integer::i
		character(len=50)::NoUsed
		if(old_version_read)then
			call readExternalData_old(A,uni)
			return
		end if
		call A%empty()
		read(uni,*)A%Flag,A%rank
		if(.not.A%GetDimFlag())then
			read(uni,*)NoUsed
			return
		end if
		call A%QuanNum%read(uni)
		call A%Deg%read(uni)
		call allocateCheck(A%dimData,A%rank)
		read(uni,*)(A%dimData(i),i=1,A%rank)
		if(A%Flag(3))then
			call allocateCheck(A%FermiArrow,A%rank)
			read(uni,*)(A%FermiArrow(i),i=1,A%rank)
		end if
		if(A%getNameFlag())then
			call allocateCheck(A%DimName,A%GetRank())
			read(uni,*)(A%DimName(i),i=1,A%GetRank())
		end if
		return
	end subroutine
	subroutine readExternalData_old(A,uni)
		class(Dimension),intent(inout)::A
		integer,intent(in)::uni
		integer::i
		character(len=50)::NoUsed
		type(DataArray)::TMP
		integer,pointer::ip(:)
		call A%empty()
		read(uni,*)A%Flag
		if(.not.A%GetDimFlag())then
			read(uni,*)NoUsed
			return
		end if
		call A%QuanNum%read(uni)
		call A%Deg%read(uni)
		call TMP%read(uni)

		call TMP%pointer(ip,1)
		call allocateCheck(A%dimData,size(ip))
		A%dimData(1:size(ip))=ip
		A%rank=size(ip)
		if(TMP%getTotalBlock().gt.1)then
			call TMP%pointer(ip,2)
			call allocateCheck(A%FermiArrow,size(ip))
			A%FermiArrow(1:size(ip))=ip
		end if
		if(A%getNameFlag())then
			call allocateCheck(A%DimName,A%GetRank())
			read(uni,*)(A%DimName(i),i=1,A%GetRank())
		end if
		return
	end subroutine

	subroutine writeExternlData_old(A,uni)
		class(Dimension),intent(in)::A
		integer,intent(in)::uni
		integer::i,rank
		type(DataArray)::TMP
		integer,pointer::ip(:)
		write(uni,*)A%Flag
		if(.not.A%GetDimFlag())then
			write(uni,*)'Empty dimension'
			return
		end if
		call A%QuanNum%write(uni)
		call A%Deg%write(uni)
		rank=A%rank
		if(A%Flag(3))then
			call TMP%allocate([rank,rank],'integer')
			call TMP%pointer(ip,1)
			ip=A%DimData(1:rank)
			call TMP%pointer(ip,2)
			ip=A%FermiArrow(1:rank)
		else
			call TMP%allocate([rank],'integer')
			call TMP%pointer(ip,1)
			ip=A%DimData(1:rank)
		end if
		call TMP%write(uni)
		if(A%getNameFlag())then
			write(uni,*)(trim(adjustl(A%DimName(i)))//"  ",i=1,A%GetRank())
		end if
		return
	end subroutine



	!*********************************************************
	!
	!  
	!
	!*********************************************************

	function GetNameFlag(dimen)
		logical::GetNameFlag
		class(Dimension),intent(in)::dimen
		GetNameFlag=dimen%Flag(4)
		return
	end function
	function GetDimensionFlag(dimen)
		logical::GetDimensionFlag
		class(Dimension),intent(in)::dimen
		GetDimensionFlag=dimen%Flag(1)
		return
	end function
	function GetSymmetryFlag(dimen)
		logical::GetSymmetryFlag
		class(Dimension),intent(in)::dimen
		GetSymmetryFlag=dimen%Flag(2)
		return
	end function
	function GetFermiFlag(dimen)
		logical::GetFermiFlag
		class(Dimension),intent(in)::dimen
		GetFermiFlag=dimen%Flag(3)
		return
	end function
	function getRank(Dimen)
		integer::getRank
		class(Dimension),intent(in)::dimen
		getRank=Dimen%Rank
		return
	end function

	function getDimi(Dimen,ith)
		integer::getDimi
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::ith
		integer,pointer::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('There is no data in the dimension',-1)
			call error_stop
		end if
		if(ith.gt.dimen%getRank())then
			call writemess('ERROR in getting dimension data,ith>rank',-1)
			call writemess('ith='+ith)
			call writemess('dimen%getRank()='+dimen%getRank())
			call error_stop
		end if
		call Dimen%pointDim(ip)
		getDimi=ip(ith)
		return
	end function

	function getDimchar(Dimen,cha)
		integer::getDimchar
		class(Dimension),intent(in)::Dimen
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('There is no data in the dimension',-1)
			call error_stop
		end if
		ith=dimen%FindOrder(cha)
		call Dimen%pointDim(ip)
		getDimchar=ip(ith)
		return
	end function

	function getDimAll(Dimen)
		integer,allocatable::getDimAll(:)
		class(Dimension),intent(in)::Dimen
		integer,pointer::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('There is no data in the dimension',-1)
			call error_stop
		end if
		allocate(getDimAll(Dimen%getRank()))
		call Dimen%pointDim(ip)
		getDimAll=ip
		return
	end function

	subroutine outAllNameChar(dimen,outchar,lenofdata,w,ch_)
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		character(len=*),intent(in),optional::ch_
		integer,intent(inout)::lenofdata
		character(len=*),intent(inout)::outchar(:)
		character(len=10)::ch
		character(len=len_of_name),pointer::Names(:)
		integer::i,rank
		if(.not.dimen%getDimFlag())then
			call writemess("There is no data in the dimension",-1)
			call error_stop
		end if
		if(.not.Dimen%getnameflag())then
			call writemess("There is no CHARACTER name in the dimension",-1)
			call error_stop()
		end if
		if(present(ch_))then
			ch=ch_
		else
			ch='dim'
		end if
		rank=dimen%getRank()
 		call Dimen%pointName(Names)
		lenofdata=0
		do i=1,rank
			if(ch.equ.'Tensor')then
				if(w.equ.(Names(i).subl.indexsymbol))then
					if(lenofdata.ge.size(outchar))then
						call writemess('ERROR in outAllNameChar',-1)
						call writemess('lenofdata='+lenofdata,-1)
						call writemess('size(outchar)='+size(outchar),-1)
						call error_stop
					end if
					lenofdata=lenofdata+1
					outchar(lenofdata)=Names(i)
				end if
			else if(ch.equ.'dim')then
				if(w.equ.(Names(i).subr.indexsymbol))then
					if(lenofdata.ge.size(outchar))then
						call writemess('ERROR in outAllNameChar',-1)
						call writemess('lenofdata='+lenofdata,-1)
						call writemess('size(outchar)='+size(outchar),-1)
						call error_stop
					end if
					lenofdata=lenofdata+1
					outchar(lenofdata)=Names(i)
				end if
			else if(ch.equ.'fullname') then
				if(w.equ.Names(i))then
					if(lenofdata.ge.size(outchar))then
						call writemess('ERROR in outAllNameChar',-1)
						call writemess('lenofdata='+lenofdata,-1)
						call writemess('size(outchar)='+size(outchar),-1)
						call error_stop
					end if
					lenofdata=lenofdata+1
					outchar(lenofdata)=Names(i)
				end if
			else
				call writemess('ERROR in getAllName(),A.B',-1)
				call writemess('input type can be:',-1)
				call writemess('  Tensor : Find the leg accoding to A',-1)
				call writemess('    dim  : Find the leg accoding to B',-1)
				call writemess('fullname : Find the leg accoding to A.B',-1)
				call error_stop
			end if
		end do
		return
	end subroutine

	function getRule1(Dimen,ith)result(Res)
		integer::Res
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		Res=1
		return
	end function
	function getRule2(Dimen,ith)result(Res)
		integer::Res
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::ith
		Res=1
		return
	end function
	function getAllRule(Dimen)Result(Res)
		integer,allocatable::Res(:)
		class(Dimension),intent(in)::dimen
		allocate(Res(Dimen%getRank()))
		Res=1
		return
	end function
	!*********************************************************
	!
	!   function for block dimension
	!
	!*********************************************************

	function getBlockDimAll(Dimen,dimi)
		integer,allocatable::getBlockDimAll(:)
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::dimi(:)
		integer::rank,i
		integer,pointer::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('There is no data in the dimension',-1)
			call error_stop
		end if
		if(.not.Dimen%getSymmetryFlag())then
			call writemess('There is NOT a symmetry dimension',-1)
			call error_stop
		end if
		rank=Dimen%getRank()
		if(size(dimi).ne.rank)then
			call writemess('ERROR in getting block dimenion, length of input',-1)
			call writemess('size(dimi)='+size(dimi),-1)
			call writemess('Dimen%getRank()='+Dimen%getRank(),-1)
			call error_stop
		end if
		allocate(getBlockDimAll(rank))
		do i=1,rank
			call Dimen%pointDeg(ip,i)
			if(dimi(i).gt.size(ip))then
				call writemess('ERROR in getting block dimenion, length of input',-1)
				call writemess('dimi(i)='+dimi(i),-1)
				call writemess('size(Deg)='+size(ip),-1)
				call error_stop
			end if
			getBlockDimAll(i)=ip(dimi(i))
		end do
		return
	end function
	function getBlockDimith(Dimen,dimi,degi)
		integer::getBlockDimith
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::dimi,degi
		integer,pointer::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('There is no data in the dimension',-1)
			call error_stop
		end if
		if(.not.Dimen%getSymmetryFlag())then
			call writemess('There is NOT a symmetry dimension',-1)
			call error_stop
		end if
		call Dimen%pointDeg(ip,dimi)
		if(degi.gt.size(ip))then
			call writemess('ERROR in getting block dimenion, length of input',-1)
			call writemess('degi='+degi,-1)
			call writemess('size(Deg)='+size(ip),-1)
			call error_stop
		end if
		getBlockDimith=ip(degi)
		return
	end function
	function getBlockDim_legi(Dimen,legith,indices)
		integer::getBlockDim_legi
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::legith(:),indices(:)
		integer,pointer::Deg(:)
		integer::i
		call Dimen%pointDeg(Deg,legith(1))
		getBlockDim_legi=Deg(indices(1))
		do i=2,size(indices)
			call Dimen%pointDeg(Deg,legith(i))
			getBlockDim_legi=getBlockDim_legi*Deg(indices(i))
		end do
		return
	end function

	function getBlockDim_leg2(Dimen,startleg,endleg,indices)
		integer::getBlockDim_leg2
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::startleg,endleg,indices(:)
		integer,pointer::Deg(:)
		integer::i,ii
		if((endleg-startleg+1).ne.size(indices))then
			call writemess('ERROR in getBlockDim')
			call error_stop
		end if
		if(endleg.gt.Dimen%getRank())then
			call writemess('ERROR in getBlockDim,endleg>rank',-1)
			call error_stop
		end if
		call Dimen%pointDeg(Deg,startleg)
		getBlockDim_leg2=Deg(indices(1))
		ii=1
		do i=startleg+1,endleg
			ii=ii+1
			call Dimen%pointDeg(Deg,i)
			if(indices(ii).gt.size(Deg))then
				call writemess('ERROR in getBlockDim,endleg>rank',-1)
				call writemess("i="+i)
				call writemess('indices(i)='+indices(i),-1)
				call writemess('size(Deg)='+size(Deg),-1)
				call error_stop
			end if
			getBlockDim_leg2=getBlockDim_leg2*Deg(indices(ii))
		end do
		return
	end function

	function nonSymDimension(Dimen)Result(REs)
		integer,allocatable::Res(:)
		class(Dimension),intent(in)::Dimen
		integer::i,rank
		integer,pointer::deg(:)
		rank=Dimen%getRank()
		allocate(Res(Rank))
		do i=1,rank
			call Dimen%pointDeg(deg,i)
			Res(i)=sum(deg)
		end do
		return
	end function
	function getTotalVolume(Dimen)Result(Res)
		integer::Res
		class(Dimension),intent(in)::Dimen
		integer::i
		integer,pointer::deg(:)
		if(Dimen%getSymmetryFlag())then
			Res=1
			do i=1,Dimen%getrank()
				call Dimen%pointDeg(deg,i)
				Res=Res*sum(deg)
			end do
		else
			call Dimen%pointDim(deg)
			Res=product(deg)
		end if
	end function

	!*********************************************************
	!
	!   function for Name
	!
	!*********************************************************
	
	subroutine FixnonSymmetryType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=1
		if(Dimen%getFermiFlag())then
			call writemess('Can NOT fix the dimension to nonSymmetryType, the dimenison is of fermionic',-1)
			call error_stop
		end if
		if(Dimen%getSymmetryFlag())then
			call writemess('Can NOT fix the dimension to nonSymmetryType, the dimenison is of symmetry',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine FixNonFermionicType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=1
		if(Dimen%getFermiFlag())then
			call writemess('Can NOT fix the dimension to nonSymmetryType, the dimenison is of fermionic',-1)
			call error_stop
		end if
		if(Dimen%getSymmetryFlag())then
			call writemess('Can NOT fix the dimension to nonSymmetryType, the dimenison is of symmetry',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine FixSymmetryType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=2
		if(Dimen%getFermiFlag())then
			call writemess('Can NOT fix the dimension to nonSymmetryType, the dimenison is of fermionic',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine FixFermionicType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=3
		return
	end subroutine
	subroutine unFixSymmetryType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=0
		return
	end subroutine
	subroutine unFixfermionicType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=0
		return
	end subroutine
	subroutine unFixnonSymmetryType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=0
		return
	end subroutine
	subroutine unFixNonFermionicType(dimen)
		class(Dimension),intent(inout)::dimen
		dimen%fixType=0
		return
	end subroutine
	subroutine check_type(dimen,SymCase)
		class(Dimension),intent(in)::dimen
		integer,intent(in)::SymCase
		if(dimen%fixType.eq.0)return
		if(dimen%fixType.eq.1)then
			if(SymCase.ne.1)then
				call writemess('ERROR: the tensor of dimension is fix to a non-symmetry one',-1)
				call writemess('SymCase='+SymCase,-1)
				call error_stop
			end if
			return
		end if
		if(dimen%fixType.eq.2)then
			if(SymCase.ne.2)then
				call writemess('ERROR: the tensor of dimension is fix to a symmetry one',-1)
				call writemess('SymCase='+SymCase,-1)
				call error_stop
			end if
			return
		end if
		if(dimen%fixType.eq.3)then
			if(SymCase.ne.3)then
				call writemess('ERROR: the tensor of dimension is fix to a fermionic one',-1)
				call writemess('SymCase='+SymCase,-1)
				call error_stop
			end if
			return
		end if
		if(dimen%fixType.eq.4)then
			if(SymCase.eq.1)then
				call writemess('ERROR: the tensor of dimension is fix to a fermionic or a symmery one',-1)
				call writemess('SymCase='+SymCase,-1)
				call error_stop
			end if
			return
		end if
		call writemess('ERROR case in check_type',-1)
		call writemess('SymCase='+SymCase,-1)
		call error_stop
	end subroutine

	subroutine setNameith(dimen,ith,cha)
		class(Dimension),intent(inout)::dimen
		integer,intent(in)::ith
		character(len=*),intent(in)::Cha
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call allocateCheck(dimen%DimName,dimen%getRank())
			dimen%Flag(4)=.true.
		end if
		if(ith.gt.dimen%getRank())then
			call writemess('ERROR in setName,ith>rank',-1)
			call writemess('ith='+ith)
			call writemess('dimen%getRank()='+dimen%getRank())
			call error_stop
		end if
		if(index(cha,indexsymbol).ne.0)then
			dimen%DimName(ith)=cha
		else
			if(index(dimen%DimName(ith),indexsymbol).ne.0)then
				dimen%DimName(ith)=cha+indexsymbol+(dimen%DimName(ith).subr.indexsymbol)
			else
				dimen%DimName(ith)=cha+indexsymbol+ith
			end if
		end if
		return
	end subroutine
	subroutine setNameAll(dimen,cha)
		class(Dimension),intent(inout)::dimen
		character(len=*),intent(in)::Cha
		logical::NewFlag
		integer::i
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		NewFlag=.false.
		if(.not.Dimen%getNameFlag())then
			call allocateCheck(dimen%DimName,dimen%getRank())
			NewFlag=.true.
			dimen%Flag(4)=.true.
		end if
		if(index(cha,indexsymbol).ne.0)then
			call writemess('ERROR in setting name to all the leg, you could in put call dim%setName("A")')
			call error_stop
		else
			if(NewFlag)then
				do i=1,Dimen%getRank()
					dimen%DimName(i)=cha+indexsymbol+i
				end do
			else
				do i=1,Dimen%getRank()
					if(index(dimen%DimName(i),indexsymbol).ne.0)then
						dimen%DimName(i)=cha+indexsymbol+(dimen%DimName(i).subr.indexsymbol)
					else
						dimen%DimName(i)=cha+indexsymbol+i
					end if
				end do
			end if
		end if
		return
	end subroutine
	subroutine setNameAll2(dimen,cha)
		class(Dimension),intent(inout)::dimen
		character(len=*),intent(in)::Cha(:)
		integer::i
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call allocateCheck(dimen%DimName,dimen%getRank())
			dimen%Flag(4)=.true.
		end if
		do i=1,Dimen%getRank()
			call setNameith(dimen,i,cha(i))
		end do
		return
	end subroutine
	subroutine setNameChaith(dimen,chaith,cha)
		class(Dimension),intent(inout)::dimen
		character(len=*),intent(in) :: chaith
		integer::ith
		character(len=*),intent(in)::Cha
		character(len=len_of_name),pointer::names(:)
		character(len=len_of_name)::TMPNAME
		logical::Flag
		if(index(cha,indexsymbol).ne.0)then
			ith=dimen%FindOrder(chaith)
			call setNameith(dimen,ith,cha)
		else
			call dimen%pointName(Names)
			Flag=.true.
			do ith=1,dimen%getRank()
				TMPNAME=Names(ith).subl.indexsymbol
				if(TMPNAME.equ.chaith)then
					Names(ith)=cha+indexsymbol+(Names(ith).subr.indexsymbol)
					Flag=.false.
				end if
			end do
			if(Flag)then
				call writemess('ERROR in setName, Can not find the name='+chaith)
				call dimen%diminfo(.true.)
				call error_stop
			end if
		end if
		return
	end subroutine

	function FindOrderith(dimen,cha)
		integer::FindOrderith
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha
		integer::i
		if(index(cha,indexsymbol).eq.0)then
			call writemess('ERROR in input name, one should input looks like A.B',-1)
			call writemess('input='+cha,-1)
			call error_stop
		end if
		if(.not.dimen%getNameFlag())then
			if(.not.dimen%GetDimFlag())then
				call writemess('The Dimension is empty',-1)
				call error_stop
			end if
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		do i=1,dimen%getRank()
			if(cha.eq.dimen%DimName(i))then
				FindOrderith=i
				return
			end if
		end do
		call writemess('Can Not find the input name',-1)
		call writemess('The name of the input leg is:'+cha,-1)
		call Dimen%diminfo(.true.)
		call error_stop
	end function
	function FindOrderAll(dimen,cha)
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha(:)
		integer::FindOrderAll(size(cha))
		integer::i
		do i=1,size(cha)
			FindOrderAll(i)=FindOrderith(dimen,cha(i))
		end do
	end function

	function Nameorderith(dimen,cha)
		integer::Nameorderith
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha
		integer::i
		if(index(cha,indexsymbol).eq.0)then
			call writemess('ERROR in input name, one should input looks like A.B',-1)
			call writemess('input='+cha,-1)
			call error_stop
		end if
		if(.not.dimen%getNameFlag())then
			if(.not.dimen%GetDimFlag())then
				call writemess('The Dimension is empty',-1)
				call error_stop
			end if
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		do i=1,dimen%getRank()
			if(cha.eq.dimen%DimName(i))then
				Nameorderith=i
				return
			end if
		end do
		Nameorderith=0
		return
	end function

	function NameorderAll(dimen,cha)
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha(:)
		integer::NameorderAll(size(cha))
		integer::i
		do i=1,size(cha)
			NameorderAll(i)=Nameorderith(dimen,cha(i))
		end do
	end function

	function getNameIth(dimen,ith)
		character(len_of_name)::getNameIth
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		if(ith.gt.dimen%getRank())then
			call writemess('ERROR in getName,ith>rank',-1)
			call writemess('ith='+ith)
			call writemess('dimen%getRank()='+dimen%getRank())
			call error_stop
		end if
		getNameIth=dimen%DimName(ith)
		return
	end function

	function getNameAll(dimen)
		character(len_of_name),allocatable::getNameAll(:)
		class(Dimension),intent(in)::dimen
		integer::i
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		allocate(getNameAll(dimen%getRank()))
		do i=1,dimen%getRank()
			getNameAll(i)=dimen%DimName(i)
		end do
		return
	end function

	function getNameCha(dimen,cha,TensorName)
		character(len_of_name)::getNameCha
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha
		logical,optional,intent(in)::TensorName
		integer::i,checkflag
		logical::Flag
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		if(present(TensorName))then
			flag=TensorName
		else
			flag=.false.
		end if
		checkflag=0
		do i=1,dimen%getRank()
			if(Flag)then
				if((Dimen%DimName(i).subl.indexsymbol).equ.cha)then
					getNameCha=Dimen%DimName(i)
					checkflag=checkflag+1
				end if
			else
				if((Dimen%DimName(i).subr.indexsymbol).equ.cha)then
					getNameCha=Dimen%DimName(i)
					checkflag=checkflag+1
				end if
			end if
		end do
		if(checkflag.eq.0)then
			call writemess('Can Not Find the name:'+cha)
			call dimen%diminfo(.true.)
			call error_stop
		end if
		if(checkflag.ne.1)then
			call writemess('There are more than 1 legd which name:'+cha)
			call dimen%diminfo(.true.)
			call error_stop
		end if
		return
	end function

	function ifName(dimen,w,ch_)
		logical::ifName
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		character(len=*),intent(in),optional::ch_
		character(len=10)::ch
		integer::i
		character(len=len_of_name),pointer::names(:)
		if(.not.dimen%getDimFlag())then
			call writemess("There is data in the dimension",-1)
			call error_stop()
		end if
		if(.not.dimen%getNameFlag())then
			call writemess("There is no name in the dimension",-1)
			call error_stop()
		end if
		call dimen%pointName(names)
		if(present(ch_))then
			ch=ch_
		else
			ch='dimension'
		end if
		ifName=.true.
		if(index(w,indexsymbol).ne.0)then
			do i=1,dimen%getRank()
				if(w.equ.names(i))return
			end do
		else
			if(ch.equ.'Tensor') then
				do i=1,dimen%getrank()
					if(w.equ.(names(i).subl.indexsymbol))return
				end do
			else
				do i=1,dimen%getrank()
					if(w.equ.(names(i).subr.indexsymbol))return
				end do
			end if
		end if
		ifName=.false.
		return
	end function

	!*********************************************************
	!
	!   setting symmetry infomation
	!
	!*********************************************************

	subroutine setQNith(Dimen,ith,QuantumNumber)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::ith
		real*4,intent(in)::QuantumNumber(:)
		real*4,pointer::sp(:)
		integer,pointer::dim(:)
		call check_type(Dimen,4)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(2))then
			Dimen%Flag(2)=.true.
			call Dimen%pointDim(Dim)
			call dimen%QuanNum%allocate(Dim,'real*4')
			call dimen%Deg%allocate(Dim,'integer')
		end if
		call Dimen%QuanNum%pointer(sp,ith)
		if(size(sp).ne.size(QuantumNumber))then
			call writemess('ERROR in setQN',-1)
			call writemess('size(QuantumNumber)='+size(QuantumNumber))
			call writemess('size(Dimen%QuanNum)='+size(sp))
			call error_stop
		end if
		sp=QuantumNumber
		
		return
	end subroutine
	subroutine setQNchar(Dimen,chaith,QuantumNumber)
		class(Dimension),intent(inout) ::Dimen
		character(len=*),intent(in)::chaith
		real*4,intent(in)::QuantumNumber(:)
		integer::ith
		ith=Dimen%FindOrder(chaith)
		call setQNith(Dimen,ith,QuantumNumber)
		return
	end subroutine
	subroutine setDegith(Dimen,ith,Deg)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::ith
		integer,intent(in)::Deg(:)
		integer,pointer::ip(:)
		integer,pointer::dim(:)
		call check_type(Dimen,4)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(2))then
			Dimen%Flag(2)=.true.
			call Dimen%pointDim(Dim)
			call dimen%QuanNum%allocate(Dim,'real*4')
			call dimen%Deg%allocate(Dim,'integer')
		end if
		call Dimen%Deg%pointer(ip,ith)
		if(size(ip).ne.size(Deg))then
			call writemess('ERROR in setQN',-1)
			call writemess('size(Deg)='+size(Deg))
			call writemess('size(Dimen%Deg)='+size(ip))
			call error_stop
		end if
		ip=Deg
		return
	end subroutine
	subroutine setDegChar(Dimen,chaith,Deg)
		class(Dimension),intent(inout) ::Dimen
		character(len=*),intent(in)::chaith
		integer,intent(in)::Deg(:)
		integer::ith
		ith=Dimen%FindOrder(chaith)
		call setDegith(Dimen,ith,Deg)
		return
	end subroutine
	subroutine setDegCharAll(Dimen,chaith,Deg)
		class(Dimension),intent(inout) ::Dimen
		character(len=*),intent(in)::chaith
		integer,intent(in)::Deg
		integer::ith
		ith=Dimen%FindOrder(chaith)
		call setDegithAll(Dimen,ith,Deg)
		return
	end subroutine
	subroutine setDegithAll(Dimen,ith,Deg)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::ith
		integer,intent(in)::Deg
		integer,pointer::ip(:)
		integer,pointer::dim(:)
		call check_type(Dimen,4)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(2))then
			Dimen%Flag(2)=.true.
			call Dimen%pointDim(Dim)
			call dimen%QuanNum%allocate(Dim,'real*4')
			call dimen%Deg%allocate(Dim,'integer')
		end if
		call Dimen%Deg%pointer(ip,ith)
		ip=Deg
		return
	end subroutine
	subroutine setDegAll(Dimen,Deg)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::Deg
		integer,pointer::ip(:)
		integer,pointer::dim(:)
		integer::i
		call check_type(Dimen,4)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(2))then
			Dimen%Flag(2)=.true.
			call Dimen%pointDim(Dim)
			call dimen%QuanNum%allocate(Dim,'real*4')
			call dimen%Deg%allocate(Dim,'integer')
		end if
		do i=1,Dimen%getRank()
			call Dimen%Deg%pointer(ip,i)
			ip=Deg
		end do
		return
	end subroutine

	subroutine setFermiArrowith(Dimen,ith,Arrow)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::ith
		integer,intent(in)::Arrow
		integer,pointer::ip(:)
		call check_type(Dimen,3)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(3))then
			Dimen%Flag(3)=.true.
			call allocateCheck(Dimen%FermiArrow,Dimen%rank)
		end if
		call Dimen%pointArrow(ip)
		if(size(ip).lt.ith)then
			call writemess('ERROR in setFermiArrowith',-1)
			call writemess('size(Arrow)='+size(ip))
			call writemess('ith='+ith)
			call error_stop
		end if
		ip(ith)=Arrow
		return
	end subroutine
	subroutine setFermiArrowChar(Dimen,chaith,Arrow)
		class(Dimension),intent(inout) ::Dimen
		character(len=*),intent(in)::chaith
		integer,intent(in)::Arrow
		integer::ith
		ith=Dimen%FindOrder(chaith)
		call setFermiArrowith(Dimen,ith,Arrow)
		return
	end subroutine
	subroutine setFermiArrowAll(Dimen,Arrow)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in)::Arrow(:)
		integer,pointer::ip(:)
		call check_type(Dimen,3)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is not empty',-1)
			call error_stop
		end if
		if(.not.Dimen%Flag(3))then
			Dimen%Flag(3)=.true.
			call allocateCheck(Dimen%FermiArrow,Dimen%rank)
		end if
		call Dimen%pointArrow(ip)
		if(size(ip).ne.size(Arrow))then
			call writemess('ERROR in setRuleith',-1)
			call writemess('size(dim%Arrow)='+size(ip))
			call writemess('size(Arrow)='+size(Arrow))
			call error_stop
		end if
		ip=Arrow
		return
	end subroutine

	subroutine unsetFermiFlag(Dimen)
		class(Dimension),intent(inout)::Dimen
		integer,pointer::arrow(:)
		if(.not.Dimen%getDimFlag())return
		if(Dimen%getFermiFlag())then
			call Dimen%pointArrow(arrow)
			arrow=0
			Dimen%Flag(3)=.false.
		end if
		return
	end subroutine

	subroutine killZeroDeg(Dimen)
		class(Dimension),intent(inout)::Dimen
		integer,pointer::dim(:),deg(:),Arrow(:),Newdim(:),NewArrow(:)
		real*4,pointer::QN(:)
		integer::i,rank,j,lenDeg,ii,Newrank,iii
		type(DataArray)::NewQuanNum
		type(DataArray)::NewDeg
		type(DataArray)::NewDimData
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		integer::ZeroDimNum
		character(len=len_of_name),pointer::Nams(:)
		if(.not.ifZeroDeg(Dimen)) return

		rank=Dimen%getRank()
		call Dimen%pointDim(dim)
		ZeroDimNum=0
		do i=1,rank
			call Dimen%pointDeg(deg,i)
			lenDeg=dim(i)
			do j=1,lenDeg
				if(deg(j).eq.0)then
					dim(i)=dim(i)-1
				end if
			end do
			if(dim(i).le.0)ZeroDimNum=ZeroDimNum+1
		end do
		if(ZeroDimNum.eq.0)then
			call NewQuanNum%allocate(dim,'real*4')
			call NewDeg%allocate(dim,'integer')
			do i=1,rank
				call Dimen%pointDeg(deg,i)
				call Dimen%pointQN(QN,i)
				call NewQuanNum%pointer(sp,i)
				call NewDeg%pointer(ip,i)
				ii=0
				do j=1,Dimen%QuanNum%getTotalData(i)
					if(deg(j).ne.0)then
						ii=ii+1
						sp(ii)=QN(j)
						ip(ii)=deg(j)
					end if
				end do
			end do
			Dimen%QuanNum=NewQuanNum
			Dimen%Deg=NewDeg
			return
		end if

		Newrank=rank-ZeroDimNum
		call NewDimData%allocate([Newrank,Newrank],'integer')
		call NewDimData%pointer(NewDim,1)
		call NewDimData%pointer(NewArrow,2)
		if(Dimen%getFermiFlag())call Dimen%pointArrow(Arrow)
		ii=0
		do i=1,rank
			if(dim(i).gt.0)then
				ii=ii+1
				NewDim(ii)=dim(i)
				if(Dimen%getFermiFlag())NewArrow(ii)=Arrow(i)
			end if
		end do

		call NewQuanNum%allocate(NewDim,'real*4')
		call NewDeg%allocate(NewDim,'integer')
		if(Dimen%Flag(4))then
			call allocateCheck(WorkingMemory_Name,newRank)
			call Dimen%pointName(Nams)
		end if

		ii=0
		do i=1,rank
			if(dim(i).gt.0)then
				ii=ii+1
				call Dimen%pointDeg(deg,i)
				call Dimen%pointQN(QN,i)
				call NewQuanNum%pointer(sp,ii)
				call NewDeg%pointer(ip,ii)
				iii=0
				do j=1,Dimen%QuanNum%getTotalData(i)
					if(deg(j).ne.0)then
						iii=iii+1
						sp(iii)=QN(j)
						ip(iii)=deg(j)
					end if
				end do
				if(Dimen%Flag(4))then
					WorkingMemory_Name(ii)=Nams(i)
				end if
			end if
		end do
		Dimen%QuanNum=NewQuanNum
		Dimen%Deg=NewDeg
		call NewDimData%pointer(ip,1)
		call allocateCheck(Dimen%dimData,size(ip))
		Dimen%dimData(1:size(ip))=ip
		if(Dimen%getFermiFlag())then
			call NewDimData%pointer(ip,2)
			call allocateCheck(Dimen%FermiArrow,size(ip))
			Dimen%FermiArrow(1:size(ip))=ip
		end if
		if(Dimen%Flag(4))Dimen%DimName(1:newRank)=WorkingMemory_Name(1:newRank)
		return
	end subroutine

	function ifZeroDeg(Dimen)
		logical::ifZeroDeg
		type(Dimension),intent(in)::dimen
		integer::i,j
		integer,pointer::deg(:),dim(:)
		ifZeroDeg=.false.
		call dimen%pointDim(dim)
		do i=1,dimen%getRank()
			call dimen%pointDeg(Deg,i)
			do j=1,dim(i)
				if(deg(j).eq.0)then
					ifZeroDeg=.true.
					return
				end if
			end do
		end do
		return
	end function

	!*********************************************************
	!
	!   getting symmetry infomation
	!
	!*********************************************************

	function getQNAll(dimen,ith)
		real*4,allocatable::getQNAll(:)
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		real*4,pointer::sp(:)
		integer::rank
		rank=dimen%getRank()
		if(ith.gt.rank)then
			call writemess('ERROR in getQN, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
		end if
		allocate(getQNAll(dimen%dim(ith)))
		call dimen%pointQN(sp,ith)
		getQNAll=sp
		return
	end function

	function getQNOne(dimen,ith,degi)
		real*4::getQNOne
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith,degi
		real*4,pointer::sp(:)
		integer::rank,dimi
		rank=dimen%getRank()
		if(ith.gt.rank)then
			call writemess('ERROR in getQN, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
		end if
		dimi=dimen%dim(ith)
		if(degi.gt.dimi)then
			call writemess('ERROR in getQN, degeneracy',-1)
			call writemess('dimen%dim(ith)='+dimen%dim(ith),-1)
			call writemess('degi='+degi,-1)
		end if
		call dimen%pointQN(sp,ith)
		getQNOne=sp(degi)
		return
	end function
	function getQNChaOne(dimen,w,degi)
		real*4::getQNChaOne
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::w
		integer,intent(in)::degi
		real*4,pointer::sp(:)
		integer::rank,dimi,ith
		ith=dimen%FindOrder(w)
		dimi=dimen%dim(ith)
		if(degi.gt.dimi)then
			call writemess('ERROR in getQN, degeneracy',-1)
			call writemess('dimen%dim(ith)='+dimen%dim(ith),-1)
			call writemess('degi='+degi,-1)
		end if
		call dimen%pointQN(sp,ith)
		getQNChaOne=sp(degi)
		return
	end function

	function getDegAll(dimen,ith)
		integer,allocatable::getDegAll(:)
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		integer,pointer::ip(:)
		integer::rank
		rank=dimen%getRank()
		if(ith.gt.rank)then
			call writemess('ERROR in getDeg, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
		end if
		allocate(getDegAll(dimen%dim(ith)))
		call dimen%pointDeg(ip,ith)
		getDegAll=ip
		return
	end function
	function getDegAll_name(dimen,aith)result(getDegAll)
		integer,allocatable::getDegAll(:)
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::aith
		integer,pointer::ip(:)
		integer::rank,ith
		rank=dimen%getRank()
		ith=dimen%FindOrder(aith)
		if(ith.gt.rank)then
			call writemess('ERROR in getDeg, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
		end if
		allocate(getDegAll(dimen%dim(ith)))
		call dimen%pointDeg(ip,ith)
		getDegAll=ip
		return
	end function

	function getDegOne(dimen,ith,degi)
		integer::getDegOne
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith,degi
		integer,pointer::ip(:)
		integer::rank,dimi
		rank=dimen%getRank()
		if(ith.gt.rank)then
			call writemess('ERROR in getDeg, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		dimi=dimen%dim(ith)
		if(degi.gt.dimi)then
			call writemess('ERROR in getDeg, degeneracy',-1)
			call writemess('dimen%dim(ith)='+dimen%dim(ith),-1)
			call writemess('degi='+degi,-1)
			call error_stop
		end if
		call dimen%pointDeg(ip,ith)
		getDegOne=ip(degi)
		return
	end function
	function getDegChaOne(dimen,w,degi)
		integer::getDegChaOne
		class(Dimension),intent(in)::dimen
		integer,intent(in)::degi
		character(len=*),intent(in)::w
		integer,pointer::ip(:)
		integer::ith,dimi
		ith=dimen%FindOrder(w)
		call dimen%pointDeg(ip,ith)
		dimi=dimen%dim(ith)
		if(degi.gt.dimi)then
			call writemess('ERROR in getDeg, degeneracy',-1)
			call writemess('dimen%dim(ith)='+dimen%dim(ith),-1)
			call writemess('degi='+degi,-1)
			call error_stop
		end if
		getDegChaOne=ip(degi)
		return
	end function

	function getQuantumNumber1(dimen,ith)result(Res)
		type(QuanNum)::Res
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		real*4,pointer::QN(:)
		integer,pointer::Deg(:),arrow(:)
		call dimen%pointQN(QN,ith)
		call dimen%pointDeg(Deg,ith)
		call Res%setQN(QN)
		call Res%setDeg(Deg)
		call Res%setRule(1)
		if(dimen%getFermiFlag())then
			call dimen%pointArrow(arrow)
			call Res%setFermiArrow(arrow(ith))
		end if
		return
	end function
	function getQuantumNumber2(dimen,cha)result(Res)
		type(QuanNum)::Res
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::cha
		integer::ith
		real*4,pointer::QN(:)
		integer,pointer::Deg(:),arrow(:)
		ith=dimen%FindOrder(cha)
		call dimen%pointQN(QN,ith)
		call dimen%pointDeg(Deg,ith)
		call Res%setQN(QN)
		call Res%setDeg(Deg)
		call Res%setRule(1)
		if(dimen%getFermiFlag())then
			call dimen%pointArrow(arrow)
			call Res%setFermiArrow(arrow(ith))
		end if
		return
	end function



	function getArrowAll(dimen)
		integer,allocatable::getArrowAll(:)
		class(Dimension),intent(in)::dimen
		integer,pointer::ip(:)
		allocate(getArrowAll(dimen%getRank()))
		call dimen%pointArrow(ip)
		getArrowAll=ip
		return
	end function

	function getArrowith(dimen,ith)
		integer::getArrowith
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		integer,pointer::ip(:)
		if(ith.gt.dimen%getRank())then
			call writemess('ERROR in getArrow, rank',-1)
			call writemess('dimen%getRank()='+dimen%getRank(),-1)
			call writemess('ith='+ith,-1)
		end if
		call dimen%pointArrow(ip)
		getArrowith=ip(ith)
		return
	end function
	function getArrowCha(dimen,w)
		integer::getArrowCha
		class(Dimension),intent(in)::dimen
		character(len=*),intent(in)::w
		integer::ith
		integer,pointer::ip(:)
		ith=dimen%FindOrder(w)
		call dimen%pointArrow(ip)
		getArrowCha=ip(ith)
		return
	end function

	!Note:
		!U(1) symmetry
		!|S,rs>,S=-1,0,1,degeneracy: 1,2,1
		!in the symmetry base
		! |-1,1>,|0,1>,|0,2>,|1,1>
		!transfrom to non-symmety index will be
		! |1>  , |2>  ,|3>  , |4>
		!input s, output imin and imax
		! example,input s=0,output imin=2,imax=3, that is output a array of (2,3)
		! example,input s=1,output imin=4,imax=4, that is output a array of (4,4)
		!
		! in Parity
		! |p,d>, p=+1,-1. degeneracy are 3 for -1 and 4 for 1
		! in the symmetry base
		! |-1,1>,|-1,2>,|-1,3>,|1,1>,|1,2>,|1,3>,|1,4>
		!transfrom to non-symmety index will be
		! |1>  , |2>  ,|3>  , |4> ,  |5>,  |6>,   |7>
		! input p=1, output imin=4,imax=7, that is output a array of (4,7)

	function QN2nonSyminde_one(dimen,ith,inQN)result(vec)
		integer,allocatable::vec(:)
		class(Dimension),intent(in)::dimen
		real*4,intent(in)::inQN
		integer,intent(in)::ith
		integer::imin,imax
		integer::i
		integer,pointer::Deg(:)
		allocate(vec(2))
		imin=1
		call Dimen%pointDeg(Deg,ith)
		do i=1,dimen%getRank()
			if(abs(Deg(i)-inQN).gt.1d-15)then
				imin=imin+Deg(i)
			else
				imax=Deg(i)
				exit
			end if
		end do
		imax=imax+imin-1
		vec=[imin,imax]
		return
	end function

	!Note:
		!U(1) symmetry
		!|S,rs>,S=-1,0,1,degeneracy: 1,2,1
		!in the symmetry base
		! |-1,1>,|0,1>,|0,2>,|1,1>
		!transfrom to non-symmety index will be
		! |1>  , |2>  ,|3>  , |4>
		!input i,where s are store in order, i is the index of s, output imin and imax
		! example,input i=2, which mean s(2)=0,output imin=2,imax=3
		! example,input i=3, which mean s(3)=1,output imin=4,imax=4

	function QN2nonSyminde_one2(dimen,ith,inde)result(vec)
		integer,allocatable::vec(:)
		class(Dimension),intent(in)::dimen
		integer,intent(in)::inde
		integer,intent(in)::ith
		integer::imin,imax
		integer::i
		integer,pointer::Deg(:)
		allocate(vec(2))
		call Dimen%pointDeg(Deg,ith)
		imin=1
		do i=1,inde-1
			imin=imin+Deg(i)
		end do
		imax=Deg(inde)+imin-1
		vec=[imin,imax]
		return
	end function

	!Note:
		!|S,rs>,S=-1,0,1,degeneracy: 1,2,1
		!in the symmetry base
		! |-1,1>,|0,1>,|0,2>,|1,1>
		!transfrom to non-symmety index will be
		! |1>  , |2>  ,|3>  , |4>
		!input s, output imin and imax,thay are array
		! example,dimen=[2,2,3],QN=[0.5,0.5,1],degeneracy=[(1,1),(1,1),(1,2,1)]
		! input [0.5,-0.5,0] output imin:[1,2,2] imax:[1,2,3] there store in a array
		!    / 1 2 2 \
		!    \ 1 2 3 /

	function QN2nonSyminde(dimen,QN) result(vec)
		integer,allocatable::vec(:,:)
		class(Dimension),intent(in)::dimen
		real,intent(in)::QN(:)
		integer::i,lendim,temp(2)
		lendim=dimen%getRank()
		allocate(vec(2,lendim))
		do i=1,dimen%getRank()
			temp=QN2nonSyminde_one(dimen,i,QN(i))
			vec(1,i)=temp(1)
			vec(2,i)=temp(2)
		end do
		return
	end function
	function QN2nonSyminde2(dimen,ith) result(vec)
		integer,allocatable::vec(:,:)
		class(Dimension),intent(in)::dimen
		integer,intent(in)::ith(:)
		integer::i,lendim,temp(2)
		lendim=dimen%getRank()
		allocate(vec(2,lendim))
		do i=1,dimen%getRank()
			temp=QN2nonSyminde_one2(dimen,i,ith(i))
			vec(1,i)=temp(1)
			vec(2,i)=temp(2)
		end do
		return
	end function


	subroutine index2QNinfo_ith(maxDeg,outQN,outDegi,inindex)
		integer,intent(in)::maxDeg(:)
		integer,intent(inout)::outQN,outDegi
		integer,intent(in)::inindex
		integer::lenQN,i,Qni,degi
		lenQN=size(maxDeg)
		outQN=1
		degi=0
		do i=1,inindex
			degi=degi+1
			if(degi.gt.maxDeg(outQN))then
				outQN=outQN+1
				degi=1
			end if
			if(Qni.gt.lenQN)then
				call writemess('ERROR, the index is larger the max index of the dimension')
				call error_stop
			end if
			outDegi=degi
		end do
		return
	end subroutine


	subroutine index2QNinfo1(Dimen,outQN,outDegi,inindex)
		class(Dimension),intent(in)::Dimen
		integer,intent(inout)::outQN(:),outDegi(:)
		integer,intent(in)::inindex(:)
		integer::i,rank
		rank=Dimen%getRank()
		if(size(outQN).lt.rank)then
			call writemess(' The length of array can not store the output data')
			call writemess('size(outQN)='+size(outQN))
			call writemess('the length output data='+rank)
			call error_stop
		end if
		if(size(outDegi).lt.rank)then
			call writemess(' The length of array can not store the output data')
			call writemess('size(outDegi)='+size(outDegi))
			call writemess('the length output data='+rank)
			call error_stop
		end if
		do i=1,rank
			call index2QNinfo_ith(Dimen%getDeg(i),outQN(i),outDegi(i),inindex(i))
		end do
		return
	end subroutine


	subroutine QNinfo2index_ith(maxDeg,outindex,inQN,inDeg)
		integer,intent(in)::maxDeg(:)
		integer,intent(in)::inQN,inDeg
		integer,intent(inout)::outindex
		integer::i,lenQN
		lenQN=size(maxDeg)
		if(inQN.gt.lenQN)then
			call writemess('ERROR in QNinfo2index_ith',-1)
			call error_stop
		end if
		outindex=0
		do i=1,inQN-1
			outindex=outindex+maxDeg(i)
		end do
		outindex=outindex+inDeg
		return
	end subroutine

	function QNinfo2index1(Dimen,QN,Deg)result(outindex)
		integer,allocatable::outindex(:)
		class(Dimension),intent(in)::Dimen
		integer,intent(in)::QN(:),Deg(:)
		integer::i,rank
		rank=Dimen%getRank()
		if(size(QN).lt.rank)then
			call writemess('ERROR in the length of input QN array')
			call writemess('size(QN)='+size(QN))
			call writemess('the rank of the dimension='+rank)
			call error_stop
		end if
		if(size(Deg).lt.rank)then
			call writemess('ERROR in the length of input Deg array')
			call writemess('size(Deg)='+size(Deg))
			call writemess('the rank of the dimension='+rank)
			call error_stop
		end if
		allocate(outindex(rank))
		do i=1,rank
			call QNinfo2index_ith(Dimen%getDeg(i),outindex(i),QN(i),Deg(i))
		end do
		return
	end function


	!*********************************************************
	!
	!   pointer
	!
	!*********************************************************

	subroutine pointDim(Dimen,ip)
		class(Dimension),target,intent(in) ::Dimen
		integer,pointer,intent(inout)::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		ip=>Dimen%dimData(1:Dimen%rank)
		return
	end subroutine
	subroutine pointArrow(Dimen,ip)
		class(Dimension),target,intent(in) ::Dimen
		integer,pointer,intent(inout)::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getFermiFlag())then
			call writemess('Do Not set arrow to the dimension yet',-1)
			call error_stop
		end if
		ip=>Dimen%FermiArrow(1:Dimen%rank)
		return
	end subroutine
	subroutine pointQNith(Dimen,ip,ith)
		class(Dimension),target,intent(in) ::Dimen
		real*4,pointer,intent(inout)::ip(:)
		integer,intent(in)::ith
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getSymmetryFlag())then
			call writemess('Do Not set quantum number to the dimension yet',-1)
			call error_stop
		end if
		call Dimen%QuanNum%pointer(ip,ith)
		return
	end subroutine
	subroutine pointQNChar(Dimen,ip,charith)
		class(Dimension),target,intent(in) ::Dimen
		real*4,pointer,intent(inout)::ip(:)
		character(len=*),intent(in)::charith
		integer::ith
		ith=Dimen%FindOrder(charith)
		call pointQNith(Dimen,ip,ith)
		return
	end subroutine
	subroutine pointDegith(Dimen,ip,ith)
		class(Dimension),target,intent(in) ::Dimen
		integer,pointer,intent(inout)::ip(:)
		integer,intent(in)::ith
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getSymmetryFlag())then
			call writemess('Do Not set degeneracies to the dimension yet',-1)
			call error_stop
		end if
		call Dimen%Deg%pointer(ip,ith)
		return
	end subroutine
	subroutine pointDegChar(Dimen,ip,charith)
		class(Dimension),target,intent(in) ::Dimen
		integer,pointer,intent(inout)::ip(:)
		character(len=*),intent(in)::charith
		integer::ith
		ith=Dimen%FindOrder(charith)
		call pointDegith(Dimen,ip,ith)
		return
	end subroutine
	subroutine pointName(Dimen,ip)
		class(Dimension),target,intent(in) ::Dimen
		character(len=len_of_Name),pointer,intent(inout)::ip(:)
		if(.not.Dimen%GetDimFlag())then
			call writemess('The dimension is empty',-1)
			call error_stop
		end if
		if(.not.Dimen%getNameFlag())then
			call writemess('There is no name in the dimension',-1)
			call error_stop
		end if
		ip=>Dimen%DimName(1:Dimen%getRank())
		return
	end subroutine


	!*********************************************************
	!
	!     paste for two dimension
	!
	!*********************************************************

	function pasteTwoDimension(A,B)result(C)
		type(Dimension)::C
		type(Dimension),intent(in)::A,B
		logical::SymmetryFlag,nameFlag
		integer::lenA,lenB,lenC,i
		integer::si,ei,AMaxi,starti_si,starti_ei
		integer,pointer::dimA(:),dimB(:),dimC(:),Cstarti(:),Cendi(:),Astarti(:),Aendi(:),Bstarti(:),Bendi(:)
		real*4,pointer::AQN(:),BQN(:),CQN(:)
		integer,pointer::Adeg(:),BDeg(:),CDeg(:)
		integer,pointer::Aip(:),Bip(:),Cip(:)
		character(len=len_of_Name),pointer::newName(:),AName(:),BName(:)

		if(.not.A%GetDimFlag())then
			call writemess('There is no data in the dimension A',-1)
			call error_stop
		end if
		if(.not.B%GetDimFlag())then
			call writemess('There is no data in the dimension B',-1)
			call error_stop
		end if
		SymmetryFlag=A%getSymmetryFlag().or.B%getSymmetryFlag()
		lenA=A%getRank()
		lenB=B%getRank()
		lenC=lenA+lenB
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,lenc+lenC+3)
		call WorkingMemory%get_memory(dimA,lenA)
		call WorkingMemory%get_memory(dimB,lenB)
		call WorkingMemory%get_memory(dimC,lenC)
		dimA=A
		dimB=B
		dimC=[dimA,dimB]
		if(SymmetryFlag)then
			call allocateDimension(C,dimC)
			C%Flag(1)=.true.
			C%Flag(2)=A%Flag(2).or.B%Flag(2)
			C%Flag(3)=A%Flag(3).or.B%Flag(3)
			C%Flag(4)=.false.
			call C%QuanNum%pointAllData(CQN)
			call C%Deg%pointAllData(CDeg)
			CQN=0
			CDeg=0
			AMaxi=0
			if(A%getSymmetryFlag())then
				call A%QuanNum%pointAllData(AQN)
				si=1
				ei=size(AQN)
				CQN(si:ei)=AQN
				call A%Deg%pointAllData(ADeg)
				CDeg(si:ei)=ADeg


				call A%QuanNum%pointStarti(Astarti)
				call A%QuanNum%pointEndi(Aendi)
				call C%QuanNum%pointStarti(Cstarti)
				call C%QuanNum%pointEndi(Cendi)
				starti_si=1
				starti_ei=size(Aendi)
				Cstarti(starti_si:starti_ei)=Astarti
				Cendi(starti_si:starti_ei)=Aendi
				call A%Deg%pointStarti(Astarti)
				call A%Deg%pointEndi(Aendi)
				call C%Deg%pointStarti(Cstarti)
				call C%Deg%pointEndi(Cendi)
				Cstarti(starti_si:starti_ei)=Astarti
				Cendi(starti_si:starti_ei)=Aendi
				AMaxi=AMaxi+maxval(Aendi)

			end if
			if(B%getSymmetryFlag())then
				call B%QuanNum%pointAllData(BQN)
				ei=size(CQN)
				si=ei-size(BQN)+1
				CQN(si:ei)=BQN
				call B%Deg%pointAllData(BDeg)
				CDeg(si:ei)=BDeg

				call B%QuanNum%pointStarti(Bstarti)
				call B%QuanNum%pointEndi(Bendi)
				call C%QuanNum%pointStarti(Cstarti)
				call C%QuanNum%pointEndi(Cendi)
				starti_ei=size(Cstarti)
				starti_si=starti_ei-size(Bstarti)+1
				Cstarti(starti_si:starti_ei)=Bstarti+AMaxi
				Cendi(starti_si:starti_ei)=Bendi+AMaxi
				call B%Deg%pointStarti(Bstarti)
				call B%Deg%pointEndi(Bendi)
				call C%Deg%pointStarti(Cstarti)
				call C%Deg%pointEndi(Cendi)
				Cstarti(starti_si:starti_ei)=Bstarti+AMaxi
				Cendi(starti_si:starti_ei)=Bendi+AMaxi


			end if
			if(C%Flag(3))then
				C%FermiArrow=0
			end if
			if(A%getFermiFlag())then
				call A%pointArrow(Aip)
				C%FermiArrow(1:lenA)=Aip
			end if
			if(B%getFermiFlag())then
				call B%pointArrow(Bip)
				C%FermiArrow(lenA+1:)=Bip
			end if


		else
			C=dimC
		end if
		call WorkingMemory%free()
		nameFlag=A%getNameFlag().or.B%getNameFlag()
		if(.not.nameFlag)return
		C%Flag(4)=.true.
		allocate(C%DimName(lenC))
		call C%pointName(NewName)
		if(A%getNameFlag())then
			call A%pointName(AName)
			newName(1:lenA)=AName
		else
			do i=1,lenA
				newName(i)='NewName.'+i
			end do
		end if
		if(B%getNameFlag())then
			call B%pointName(BName)
			newName(lenA+1:)=BName
		else
			do i=lenA+1,lenC
				newName(i)='NewName.'+i
			end do
		end if
		return
	end function


	!*********************************************************
	!
	! permutation
	!
	!*********************************************************

	subroutine permutationDimension(inoutD,newOrder)
		type(Dimension),intent(inout)::inoutD
		integer,intent(in)::newOrder(:)
		integer, pointer :: oip(:),ip(:)
		character(len=len_of_name),pointer::ap(:)
		integer::rank,i
		if(.not.inoutD%GetDimFlag())then
			call writemess(' The dimension is empty',-1)
			call error_stop
		end if
		rank=inoutD%getRank()
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(oip,rank)
		oip=inoutD%dim()
		call inoutD%pointDim(ip)
		do i=1,rank
			ip(i)=oip(newOrder(i))
		end do
		if(inoutD%getSymmetryFlag())then
			call reorderBlock(inoutD%QuanNum,newOrder)
			call reorderBlock(inoutD%Deg,newOrder)
		end if
		if(inoutD%getFermiFlag())then
			call inoutD%pointArrow(ip)
			oip=ip
			do i=1,rank
				ip(i)=oip(newOrder(i))
			end do
		end if
		if(inoutD%getNameFlag())then
			call inoutD%pointName(ap) 
			call allocateCheck(WorkingMemory_Name,rank)
			WorkingMemory_Name(1:rank)=ap
			do i=1,rank
				ap(i)=WorkingMemory_Name(newOrder(i))
			end do
		end if
		call WorkingMemory%free()
		if(check_same_name_flag)call check_input_order(newOrder)
		return
	end subroutine


	!*********************************************************
	!
	!  sub-dimension
	!
	!*********************************************************

	function subDimension(Dimen,vec)Result(Res)
		type(Dimension)::Res
		type(Dimension),intent(in)::Dimen
		integer,intent(in)::vec(:)
		integer::ith,jth
		integer,pointer::dim(:),deg(:),NewDeg(:),NewDim(:)
		real*4,pointer::QN(:),NewQN(:)
		integer,pointer::Arrow(:),NewArrow(:)
		character(len=len_of_Name),pointer::Names(:)
		integer,pointer::QNBlockLen(:)
		integer::rank,i,ii,NewRank
		if(.not.Dimen%Flag(1))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		ith=vec(1)
		jth=vec(2)
		rank=Dimen%getRank()
		if((ith.eq.1).and.(jth.eq.rank))then
			Res=Dimen
			return
		end if
		if(jth.gt.rank)then
			call writemess('ERROR in subDimension',-1)
			call writemess('Dimen%getRank()='+Dimen%getRank(),-1)
			call writemess('jth='+jth,-1)
		end if
		if(Dimen%Flag(2))then
			call Dimen%pointDim(Dim)
			NewRank=jth-ith+1
			call allocateDimension(Res,Dim(ith:jth))
			Res%Flag(2)=.true.
			call WorkingMemory%check()
			call WorkingMemory%get_memory(QNBlockLen,NewRank)
			ii=0
			do i=ith,jth
				ii=ii+1
				QNBlockLen(ii)=Dim(i)
			end do
			call Res%QuanNum%allocate(QNBlockLen,2)
			call Res%Deg%allocate(QNBlockLen,1)
			ii=0
			do i=ith,jth
				ii=ii+1
				call Dimen%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call Dimen%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			call Dimen%pointDim(Dim)
			call Res%pointDim(NewDim)
			NewDim=Dim(ith:jth)
			if(Dimen%Flag(3))then
				Res%Flag(3)=.true.
				call Dimen%pointArrow(Arrow)
				call Res%pointArrow(NewArrow)
				NewArrow=Arrow(ith:jth)
			else
				Res%Flag(3)=.false.
			end if
			call WorkingMemory%free()
		else
			call Dimen%pointDim(Dim)
			Res=Dim(ith:jth)
		end if

		if(Dimen%Flag(4))then
			call allocateCheck(Res%DimName,rank)
			call Dimen%pointName(Names)
			ii=0
			do i=ith,jth
				ii=ii+1
				Res%DimName(ii)=Names(i)
			end do
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
		return
	end function
	function subDimension2(Dimen,ith)Result(Res)
		type(Dimension)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::ith
		call subDimensionSubroutine(Dimen,[ith,ith],Res)
	end function
	function subDimension4(Dimen,aith)Result(Res)
		type(Dimension)::Res
		type(Dimension),intent(in)::dimen
		character(len=*),intent(in)::aith
		integer::ith
		ith=dimen%FindOrder(aith)
		call subDimensionSubroutine(Dimen,[ith,ith],Res)
	end function


	subroutine subDimensionSubroutine(Dimen,vec,Res)
		type(Dimension),intent(inout)::Res
		type(Dimension),intent(in)::Dimen
		integer,intent(in)::vec(:)
		integer::ith,jth
		integer,pointer::dim(:),deg(:),NewDeg(:),NewDim(:)
		real*4,pointer::QN(:),NewQN(:)
		integer,pointer::Arrow(:),NewArrow(:)
		character(len=len_of_Name),pointer::Names(:)
		integer,pointer::QNBlockLen(:)
		integer::rank,i,ii,NewRank
		if(.not.Dimen%Flag(1))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		ith=vec(1)
		jth=vec(2)
		rank=Dimen%getRank()
		if((ith.eq.1).and.(jth.eq.rank))then
			Res=Dimen
			return
		end if
		if(jth.gt.rank)then
			call writemess('ERROR in subDimension',-1)
			call writemess('Dimen%getRank()='+Dimen%getRank(),-1)
			call writemess('jth='+jth,-1)
		end if
		if(Dimen%Flag(2))then
			call Dimen%pointDim(Dim)
			NewRank=jth-ith+1
			call allocateDimension(Res,Dim(ith:jth))
			Res%Flag(2)=.true.
			call WorkingMemory%check()
			call WorkingMemory%get_memory(QNBlockLen,NewRank)
			ii=0
			do i=ith,jth
				ii=ii+1
				QNBlockLen(ii)=Dim(i)
			end do
			call Res%QuanNum%allocate(QNBlockLen,2)
			call Res%Deg%allocate(QNBlockLen,1)
			ii=0
			do i=ith,jth
				ii=ii+1
				call Dimen%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call Dimen%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			call Dimen%pointDim(Dim)
			call Res%pointDim(NewDim)
			NewDim=Dim(ith:jth)
			if(Dimen%Flag(3))then
				Res%Flag(3)=.true.
				call Dimen%pointArrow(Arrow)
				call Res%pointArrow(NewArrow)
				NewArrow=Arrow(ith:jth)
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
			call WorkingMemory%free()
		else
			call Dimen%pointDim(Dim)
			Res=Dim(ith:jth)
			call check_type(Res,1)
		end if

		if(Dimen%Flag(4))then
			call allocateCheck(Res%DimName,rank)
			call Dimen%pointName(Names)
			ii=0
			do i=ith,jth
				ii=ii+1
				Res%DimName(ii)=Names(i)
			end do
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
		return
	end subroutine

	subroutine pasteTwoSubDim(Res,A,Aith,Ajth,B,Bith,Bjth)
		type(Dimension),intent(inout)::Res
		type(Dimension),intent(in)::A,B
		integer,intent(in)::Aith,Ajth,Bith,Bjth
		integer,pointer::dimA(:),dimB(:),NewDim(:),Arrow(:)
		integer,pointer::NewArrow(:),Deg(:),NewDeg(:)
		character(len=len_of_Name),pointer::NameA(:),NameB(:)
		integer::rankA,rankB,i,ii,NewRank,sublenA,sublenB
		real*4,pointer::QN(:),NewQN(:)
		if((.not.A%Flag(1)).or.(.not.B%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		rankA=A%getRank()
		rankB=B%getRank()
		sublenA=Ajth-Aith+1
		sublenB=Bjth-Bith+1
		NewRank=sublenA+sublenB
		if(Ajth.gt.rankA)then
			call writemess('ERROR in subDimension',-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call writemess('Ajth='+Ajth,-1)
		end if
		if(Bjth.gt.rankB)then
			call writemess('ERROR in subDimension',-1)
			call writemess('B%getRank()='+B%getRank(),-1)
			call writemess('Bjth='+Bjth,-1)
		end if
		!call Res%empty()
		if(A%Flag(2).or.B%Flag(2))then
			if(.not.(A%Flag(2).and.B%Flag(2)))then
				call writemess('ERROR in pasteTwoSubDim',-1)
				call A%diminfo(.true.)
				call B%diminfo(.true.)
				call error_stop
			end if
			call A%pointDim(dimA)
			call B%pointDim(dimB)
			call allocateDimension(Res,[dimA(Aith:Ajth),dimB(Bith:Bjth)])
			Res%Flag(2)=.true.
			call Res%pointDim(NewDim)
			!ii=0
			!do i=Aith,Ajth
			!	ii=ii+1
			!	QNBlockLen(ii)=DimA(i)
			!end do
			!do i=Bith,Bjth
			!	ii=ii+1
			!	QNBlockLen(ii)=DimB(i)
			!end do
			call Res%pointDim(NewDim)
			call Res%QuanNum%allocate(NewDim,2)
			call Res%Deg%allocate(NewDim,1)
			ii=0
			do i=Aith,Ajth
				ii=ii+1
				call A%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call A%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			do i=Bith,Bjth
				ii=ii+1
				call B%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call B%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			!call Res%pointDim(NewDim)
			!NewDim=[dimA(Aith:Ajth),dimB(Bith:Bjth)]
			if(A%Flag(3))then
				Res%Flag(3)=.true.
				call A%pointArrow(Arrow)
				call Res%pointArrow(NewArrow)
				NewArrow(1:sublenA)=Arrow(Aith:Ajth)
				call B%pointArrow(Arrow)
				NewArrow(sublenA+1:)=Arrow(Bith:Bjth)
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
		else
			call A%pointDim(dimA)
			call B%pointDim(dimB)
			Res=[dimA(Aith:Ajth),dimB(Bith:Bjth)]
			call check_type(Res,1)
		end if

		if(A%Flag(4).or.B%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			if(A%Flag(4))then
				call A%pointName(NameA)
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)=NameA(i)
				end do
			else
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			if(B%Flag(4))then
				call B%pointName(NameB)
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)=NameB(i)
				end do
			else
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine

	subroutine pasteDimQNDim(Res,A,Aith,Ajth,QuanNumber,B,Bith,Bjth)
		type(Dimension),target,intent(inout)::Res
		type(Dimension),target,intent(in)::A,B
		type(QuanNum),intent(in)::QuanNumber
		integer,intent(in)::Aith,Ajth,Bith,Bjth
		integer,pointer::dimA(:),dimB(:),NewDim(:),Arrow(:)
		integer,pointer::NewArrow(:),Deg(:),NewDeg(:)
		character(len=len_of_Name),pointer::NameA(:),NameB(:)
		integer::rankA,rankB,i,ii,NewRank,sublenA,sublenB
		real*4,pointer::QN(:),NewQN(:)
		type(Dimension),pointer::DP1,DP2
		if((.not.A%Flag(1)).or.(.not.B%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		DP1=>Res
		DP2=>A
		if(associated(DP1,DP2))then
			call writemess('input dimension can not be the same variable',-1)
			call error_stop
		end if
		DP1=>Res
		DP2=>B
		if(associated(DP1,DP2))then
			call writemess('input dimension can not be the same variable',-1)
			call error_stop
		end if
		rankA=A%getRank()
		rankB=B%getRank()
		sublenA=Ajth-Aith+1
		sublenB=Bjth-Bith+1
		NewRank=sublenA+sublenB+1
		if(Ajth.gt.rankA)then
			call writemess('ERROR in subDimension',-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call writemess('Ajth='+Ajth,-1)
		end if
		if(Bjth.gt.rankB)then
			call writemess('ERROR in subDimension',-1)
			call writemess('B%getRank()='+B%getRank(),-1)
			call writemess('Bjth='+Bjth,-1)
		end if
		if(A%Flag(2).or.B%Flag(2))then
			if(.not.(A%Flag(2).and.B%Flag(2)))then
				call writemess('ERROR in pasteDimQNDim',-1)
				call A%diminfo(.true.)
				call B%diminfo(.true.)
				call error_stop
			end if
			call A%pointDim(dimA)
			call B%pointDim(dimB)
			call allocatedimension(Res,[dimA(Aith:Ajth),QuanNumber%getQNlength(),dimB(Bith:Bjth)])
			Res%Flag(2)=.true.
			call Res%pointDim(NewDim)
			call Res%QuanNum%allocate(NewDim,2)
			call Res%Deg%allocate(NewDim,1)
			ii=0
			do i=Aith,Ajth
				ii=ii+1
				call A%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call A%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			ii=ii+1
			call Res%pointQN(NewQN,ii)
			NewQN=QuanNumber%getQN()
			call Res%pointDeg(NewDeg,ii)
			NewDeg=QuanNumber%getDeg()
			do i=Bith,Bjth
				ii=ii+1
				call B%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call B%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			if(A%Flag(3))then
				Res%Flag(3)=.true.
				call A%pointArrow(Arrow)
				call Res%pointArrow(NewArrow)
				NewArrow(1:sublenA)=Arrow(Aith:Ajth)
				NewArrow(sublenA+1)=QuanNumber%getFermiArrow()
				call B%pointArrow(Arrow)
				NewArrow(sublenA+2:)=Arrow(Bith:Bjth)
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
		else
			call writemess('ERROR in pasteDimQNDim',-1)
			call A%diminfo(.true.)
			call B%diminfo(.true.)
			call error_stop
		end if

		if(A%Flag(4).or.B%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			if(A%Flag(4))then
				call A%pointName(NameA)
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)=NameA(i)
				end do
			else
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			ii=ii+1
			Res%DimName(ii)='NoName.'+ii
			if(B%Flag(4))then
				call B%pointName(NameB)
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)=NameB(i)
				end do
			else
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine

	subroutine pasteQNDim(Res,QuanNumber,B,Bith,Bjth)
		type(Dimension),target,intent(inout)::Res
		type(Dimension),target,intent(in)::B
		type(QuanNum),intent(in)::QuanNumber
		integer,intent(in)::Bith,Bjth
		integer,pointer::dimB(:),NewDim(:),Arrow(:)
		integer,pointer::NewArrow(:),Deg(:),NewDeg(:)
		character(len=len_of_Name),pointer::NameB(:)
		integer::rankB,i,ii,NewRank,sublenB
		real*4,pointer::QN(:),NewQN(:)
		type(Dimension),pointer::DP1,DP2
		if((.not.B%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		DP1=>Res
		DP2=>B
		if(associated(DP1,DP2))then
			call writemess('input dimension can not be the same variable',-1)
			call error_stop
		end if
		rankB=B%getRank()
		sublenB=Bjth-Bith+1
		NewRank=sublenB+1
		if(Bjth.gt.rankB)then
			call writemess('ERROR in subDimension',-1)
			call writemess('B%getRank()='+B%getRank(),-1)
			call writemess('Bjth='+Bjth,-1)
		end if
		if(B%Flag(2))then
			call B%pointDim(dimB)
			call allocatedimension(Res,[QuanNumber%getQNlength(),dimB(Bith:Bjth)])
			Res%Flag(2)=.true.
			call Res%pointDim(NewDim)
			call Res%QuanNum%allocate(NewDim,2)
			call Res%Deg%allocate(NewDim,1)
			ii=0
			ii=ii+1
			call Res%pointQN(NewQN,ii)
			NewQN=QuanNumber%getQN()
			call Res%pointDeg(NewDeg,ii)
			NewDeg=QuanNumber%getDeg()
			do i=Bith,Bjth
				ii=ii+1
				call B%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call B%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			if(B%Flag(3))then
				Res%Flag(3)=.true.
				call Res%pointArrow(NewArrow)
				NewArrow(1)=QuanNumber%getFermiArrow()
				call B%pointArrow(Arrow)
				NewArrow(2:)=Arrow(Bith:Bjth)
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
		else
			call writemess('ERROR in pasteDimQNDim',-1)
			call B%diminfo(.true.)
			call error_stop
		end if

		if(B%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			ii=0
			ii=ii+1
			Res%DimName(ii)='NoName.'+ii
			if(B%Flag(4))then
				call B%pointName(NameB)
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)=NameB(i)
				end do
			else
				do i=Bith,Bjth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine

	subroutine pasteDimQN(Res,A,Aith,Ajth,QuanNumber)
		type(Dimension),target,intent(inout)::Res
		type(Dimension),target,intent(in)::A
		type(QuanNum),intent(in)::QuanNumber
		integer,intent(in)::Aith,Ajth
		integer,pointer::dimA(:),NewDim(:),Arrow(:)
		integer,pointer::NewArrow(:),Deg(:),NewDeg(:)
		character(len=len_of_Name),pointer::NameA(:)
		integer::rankA,i,ii,NewRank,sublenA
		real*4,pointer::QN(:),NewQN(:)
		type(Dimension),pointer::DP1,DP2
		if((.not.A%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		DP1=>Res
		DP2=>A
		if(associated(DP1,DP2))then
			call writemess('input dimension can not be the same variable',-1)
			call error_stop
		end if
		rankA=A%getRank()
		sublenA=Ajth-Aith+1
		NewRank=sublenA+1
		if(Ajth.gt.rankA)then
			call writemess('ERROR in subDimension',-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call writemess('Ajth='+Ajth,-1)
		end if
		if(A%Flag(2))then
			call A%pointDim(dimA)
			call allocatedimension(Res,[dimA(Aith:Ajth),QuanNumber%getQNlength()])
			Res%Flag(2)=.true.
			call Res%pointDim(NewDim)
			call Res%QuanNum%allocate(NewDim,2)
			call Res%Deg%allocate(NewDim,1)
			ii=0
			do i=Aith,Ajth
				ii=ii+1
				call A%pointQN(QN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=QN
				call A%pointDeg(Deg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=Deg
			end do
			ii=ii+1
			call Res%pointQN(NewQN,ii)
			NewQN=QuanNumber%getQN()
			call Res%pointDeg(NewDeg,ii)
			NewDeg=QuanNumber%getDeg()

			if(A%Flag(3))then
				Res%Flag(3)=.true.
				call A%pointArrow(Arrow)
				call Res%pointArrow(NewArrow)
				NewArrow(1:sublenA)=Arrow(Aith:Ajth)
				NewArrow(sublenA+1)=QuanNumber%getFermiArrow()
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
		else
			call writemess('ERROR in pasteDimQNDim',-1)
			call A%diminfo(.true.)
			call error_stop
		end if

		if(A%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			if(A%Flag(4))then
				call A%pointName(NameA)
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)=NameA(i)
				end do
			else
				ii=0
				do i=Aith,Ajth
					ii=ii+1
					Res%DimName(ii)='NoName.'+ii
				end do
			end if
			ii=ii+1
			Res%DimName(ii)='NoName.'+ii
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine

	subroutine pastevecSubDim(Res,dimA,B,Bith,Bjth)
		type(Dimension),intent(inout)::Res
		type(Dimension),intent(in)::B
		integer,intent(in)::dimA(:)
		integer,intent(in)::Bith,Bjth
		integer,pointer::dimB(:),NewDim(:)
		character(len=len_of_Name),pointer::NameB(:)
		integer::rankA,rankB,i,ii,NewRank,sublenA,sublenB
		if((.not.B%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		rankA=size(dimA)
		rankB=B%getRank()
		sublenA=rankA
		sublenB=Bjth-Bith+1
		NewRank=sublenA+sublenB
		if(Bjth.gt.rankB)then
			call writemess('ERROR in subDimension',-1)
			call writemess('B%getRank()='+B%getRank(),-1)
			call writemess('Bjth='+Bjth,-1)
		end if
		if(B%Flag(2))then
			call writemess('ERROR for paste [vec,dimension], dimension cannot be of symmetry',-1)
			call error_Stop
		else
			call B%pointDim(dimB)
			Res=[dimA,dimB(Bith:Bjth)]
			call check_type(Res,1)
		end if

		if(B%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			ii=0
			do i=1,rankA
				ii=ii+1
				Res%DimName(ii)='NoName.'+ii
			end do
			call B%pointName(NameB)
			do i=Bith,Bjth
				ii=ii+1
				Res%DimName(ii)=NameB(i)
			end do
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine
	subroutine pasteSubDimVec(Res,A,Aith,Ajth,DimB)
		type(Dimension),intent(inout)::Res
		type(Dimension),intent(in)::A
		integer,intent(in)::DimB(:)
		integer,intent(in)::Aith,Ajth
		integer,pointer::dimA(:),NewDim(:)
		character(len=len_of_Name),pointer::NameA(:)
		integer::rankA,rankB,i,ii,NewRank,sublenA,sublenB
		real*4,pointer::QN(:),NewQN(:)
		if((.not.A%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		rankA=A%getRank()
		rankB=size(DimB)
		sublenA=Ajth-Aith+1
		sublenB=rankB
		NewRank=sublenA+sublenB
		if(Ajth.gt.rankA)then
			call writemess('ERROR in subDimension',-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call writemess('Ajth='+Ajth,-1)
		end if
		if(A%Flag(2))then
			call writemess('ERROR for paste [dimension,vec], dimension cannot be of symmetry',-1)
			call error_Stop
		else
			call A%pointDim(dimA)
			Res=[dimA(Aith:Ajth),dimB]
			call check_type(Res,1)
		end if

		if(A%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			call A%pointName(NameA)
			ii=0
			do i=Aith,Ajth
				ii=ii+1
				Res%DimName(ii)=NameA(i)
			end do
			do i=1,rankB
				ii=ii+1
				Res%DimName(ii)='NoName.'+ii
			end do
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
	end subroutine




	!insertDimension: use in Quantumsplit
		! the legi of the Dimen replace by subDim

	subroutine insertDimension(A,subDim,legi,Res)
		type(Dimension),target,intent(in)::A,subDim
		type(Dimension),target,intent(inout)::Res
		integer,intent(in)::legi
		type(Dimension),pointer::DP1,DP2
		integer::NewRank,ARank,SubRank,ii,i
		integer,pointer::ADim(:),SDim(:),NewDim(:)
		integer,pointer::AArrow(:),SArrow(:),NewArrow(:)
		integer,pointer::ADeg(:),SDeg(:),NewDeg(:)
		real*4,pointer::NewQN(:),AQN(:),SQN(:)
		character(len=len_of_Name),pointer::NameA(:),NameS(:)
		if((.not.A%Flag(1)))then
			call writemess('The dimension is empty',-1)
			call error_Stop
		end if
		DP1=>Res
		DP2=>A
		if(associated(DP1,DP2))then
			call writemess('input dimension can not be the same variable',-1)
			call error_stop
		end if
		ARank=A%getRank()
		if(Arank.eq.1)then
			if(legi.ne.1)then
				call writemess('ERROR in insertDimension',-1)
				call error_stop
			end if
			Res=subDim
			return
		end if
		SubRank=subDim%getRank()
		NewRank=ARank+SubRank-1
		if(A%Flag(2).or.SubDim%Flag(2))then
			call A%pointDim(ADim)
			call subDim%pointDim(SDim)
			if(legi.eq.1)then
				call allocatedimension(Res,[SDim,ADim(2:Arank)])
			else if(legi.eq.Arank)then
				call allocatedimension(Res,[ADim(1:Arank-1),SDim])
			else
				call allocatedimension(Res,[ADim(1:legi-1),SDim,ADim(legi+1:Arank)])
			end if
			Res%Flag(2)=.true.
			call Res%pointDim(NewDim)
			call Res%QuanNum%allocate(NewDim,2)
			call Res%Deg%allocate(NewDim,1)
			ii=0
			if(legi.gt.1)then
				do i=1,legi-1
					ii=ii+1
					call A%pointQN(AQN,i)
					call Res%pointQN(NewQN,ii)
					NewQN=AQN
					call A%pointDeg(ADeg,i)
					call Res%pointDeg(NewDeg,ii)
					NewDeg=ADeg
				end do
			end if
			do i=1,SubRank
				ii=ii+1
				call subDim%pointQN(SQN,i)
				call Res%pointQN(NewQN,ii)
				NewQN=SQN
				call subDim%pointDeg(SDeg,i)
				call Res%pointDeg(NewDeg,ii)
				NewDeg=SDeg
			end do
			if(legi.lt.Arank)then
				do i=legi+1,Arank
					ii=ii+1
					call A%pointQN(AQN,i)
					call Res%pointQN(NewQN,ii)
					NewQN=AQN
					call A%pointDeg(ADeg,i)
					call Res%pointDeg(NewDeg,ii)
					NewDeg=ADeg
				end do
			end if

			if(A%Flag(3).or.SubDim%Flag(3))then
				Res%Flag(3)=.true.
				call A%pointArrow(AArrow)
				call Res%pointArrow(NewArrow)
				call subDim%pointArrow(SArrow)
				if(legi.eq.1)then
					NewArrow=[SArrow,AArrow(2:ARank)]
				else if(legi.eq.Arank)then
					NewArrow=[AArrow(1:ARank-1),SArrow]
				else
					NewArrow=[AArrow(1:legi-1),SArrow,AArrow(legi+1:Arank)]
				end if
				call check_type(Res,3)
			else
				Res%Flag(3)=.false.
				call check_type(Res,2)
			end if
		else
			call A%pointDim(ADim)
			call subDim%pointDim(SDim)
			if(legi.eq.1)then
				Res=[SDim,ADim(2:Arank)]
			else if(legi.eq.Arank)then
				Res=[ADim(1:Arank-1),SDim]
			else
				Res=[ADim(1:legi-1),SDim,ADim(legi+1:Arank)]
			end if
			call check_type(Res,1)
		end if

		if(A%Flag(4).or.SubDim%Flag(4))then
			call allocateCheck(Res%DimName,NewRank)
			call A%pointName(NameA)
			call subDim%pointName(NameS)
			ii=0
			if(legi.gt.1)then
				do i=1,legi-1
					ii=ii+1
					Res%DimName(ii)=NameA(i)
				end do
			end if
			do i=1,SubRank
				ii=ii+1
				Res%DimName(ii)=NameS(i)
			end do
			if(legi.lt.Arank)then
				do i=legi+1,Arank
					ii=ii+1
					Res%DimName(ii)=NameA(i)
				end do
			end if
			Res%Flag(4)=.true.
		else
			Res%Flag(4)=.false.
		end if
		return
	end subroutine

	!*********************************************************
	!
	!  code for MPI
	!
	!*********************************************************


	subroutine send_Dimension(Dimen1,Dimen2,ID1,ID2,MPIcommon)
		type(Dimension),intent(in)::Dimen1
		type(Dimension),intent(inout)::Dimen2
		integer,intent(in)::ID1,ID2
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,lenDimData,istatus(MPI_STATUS_SIZE),i,nameflag,mpi_comm,rank
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		lenDimData=0
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(present(MPIcommon))then
			if((ID1.ge.proNum).or.(ID2.ge.proNum))return
		end if
		
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
		!******************   Flag  ***********************************************************
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%Flag,4,MPI_logical,ID2,tag,mpi_comm,ierr)
			if(.not.Dimen1%Flag(1)) return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%Flag,4,MPI_logical,ID1,tag,mpi_comm,istatus,ierr)
			if(.not.Dimen2%Flag(1)) then
				call Dimen2%empty()
				return
			end if
		end if

		call send_DataArray(Dimen1%QuanNum,Dimen2%QuanNum,ID1,ID2,MPIcommon)
		call send_DataArray(Dimen1%Deg,Dimen2%Deg,ID1,ID2,MPIcommon) 
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%rank,1,MPI_integer,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%rank,1,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		if(proID.eq.ID1) then
			call mpi_send(Dimen1%dimdata,Dimen1%rank,MPI_integer,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%dimdata,Dimen2%rank,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		if(proID.eq.ID1) then
			if(Dimen1%Flag(3)) then
				call mpi_send(Dimen1%FermiArrow,Dimen1%rank,MPI_integer,ID2,tag,MPI_Comm,ierr)
			end if
		end if

		if(proID.eq.ID2) then
			if(Dimen2%Flag(3)) then
				call allocateCheck(Dimen2%FermiArrow,Dimen2%rank)
				call mpi_recv(Dimen2%FermiArrow,Dimen2%rank,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
			end if
		end if

		
		if(proID.eq.ID1) then
			if(.not.Dimen1%Flag(4)) return
			rank=Dimen1%getRank()
			call mpi_send(Dimen1%DimName(1:rank),len_of_name*rank,&
										MPI_character,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			if(.not.Dimen2%Flag(4)) return

			rank=Dimen2%getRank()
			call allocateCheck(Dimen2%DimName,rank)
			call mpi_recv(Dimen2%DimName(1:rank),len_of_name*rank,&
										MPI_character,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		return
	end subroutine
	
	subroutine BCAST_Dimension(Dimen1,ID,MPIcommon)
		type(Dimension),intent(inout)::Dimen1
		integer,intent(in)::ID
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE),i,mpi_comm,rank
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(present(MPIcommon))then
			if(ID.ge.proNum)return
		end if
		
		!******************   sample_dimension_flag  *******************************
		call MPI_BCAST(Dimen1%Flag,4,MPI_logical,ID,mpi_comm,ierr)
		call BCAST_DataArra(Dimen1%QuanNum,ID,MPIcommon)
		call BCAST_DataArra(Dimen1%Deg,ID,MPIcommon)

		call mpi_BCAST(Dimen1%rank,1,MPI_integer,ID,MPI_Comm,ierr)
		call mpi_BCAST(Dimen1%dimdata,Dimen1%rank,MPI_integer,ID,MPI_Comm,ierr)
		if(Dimen1%Flag(3))then
			if(proID.ne.ID) then
				call allocateCheck(Dimen1%FermiArrow,Dimen1%rank)
			end if
			call mpi_BCAST(Dimen1%FermiArrow,Dimen1%rank,MPI_integer,ID,MPI_Comm,ierr)
		end if

		if(.not.Dimen1%Flag(4))return

		rank=Dimen1%getRank()
		if(proID.ne.ID) then
			call allocateCheck(Dimen1%DimName,rank)
		end if
		call mpi_BCAST(Dimen1%DimName(1:rank),len_of_name*rank,&
												MPI_character,ID,MPI_Comm,ierr)
		return
	end subroutine


end module

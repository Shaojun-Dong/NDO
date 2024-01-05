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

module Parity_Tools
	use Dimension_Tools
	use data_base_Tools
	use QuantumNumber_Tools 
	use Tools
	use memory_type
	implicit none
	private
	type(memory),private::WorkingMemory
	public::checkParityRule
	public::ParityRule
	public::ParityifParity
	public::ParityParity
	public::ParityNewQN
	public::hermitian_conjugate_parity_dimension
	public::ParityQN
	
	
contains

	function ParityQN(deg,arrow)
		type(QuanNum)::ParityQN
		integer,intent(in)::deg(2)
		integer,optional,intent(in)::arrow
		call ParityQN%setQN([-1.,1.])
		call ParityQN%setdeg(deg)
		call ParityQN%setRule(1)
		if(present(arrow))then
			call ParityQN%setFermiArrow(arrow)
		end if
		return
	end function
	
	subroutine checkParityRule(D1,i1,D2,i2)!use in contract
		type(Dimension),intent(in)::D1,D2
		integer,intent(in)::i1(:),i2(:)
		integer::i,j
		real*4,pointer::QN1(:),QN2(:)
		do i=1,size(i1)
			call D1%pointQN(QN1,i1(i))
			call D2%pointQN(QN2,i2(i))
			do j=1,size(QN1)
				if(QN1(j).ne.QN2(j))then
					call writemess('ERROR in contract for the parity Symmetry Rule',-1)
					call D1%diminfo(.true.)
					call D2%diminfo(.true.)
					call writemess(QN1,-1)
					call writemess(QN2,-1)
					call error_stop
				end if
			end do
		end do
		return
	end subroutine
	
	
	function ParityRule(dimen,indices)
		logical::ParityRule
		type(Dimension),intent(in)::dimen
		integer,intent(in)::indices(:)
		real*4,pointer::Q(:)
		integer::i,QNcounter
		call dimen%pointQN(Q,1)
		QNcounter=0
		if(Q(indices(1)).lt.0)QNcounter=QNcounter+1

		do i=2,dimen%getRank()
			call dimen%pointQN(Q,i)
			if(Q(indices(i)).lt.0)QNcounter=QNcounter+1
		end do
		if(mod(QNcounter,2).eq.0)then
			ParityRule=.true.
		else
			ParityRule=.false.
		end if
		return
	end function

	subroutine ParityifParity(REs,dimen,vec,legi)
		logical,intent(inout)::REs
		type(Dimension),intent(in)::dimen
		integer,intent(in)::vec(:),legi(:)
		integer::i,rank,sizelegi,QNcounter
		real*4,pointer::Q(:)
		rank=dimen%getRank()
		sizelegi=size(legi)
		if(rank.ne.size(vec))then
			call writemess("ERROR in calculate the parity",-1)
			call error_stop
		end if
		QNcounter=0
		do i=1,sizelegi
			call Dimen%pointQN(Q,legi(i))
			if(Q(vec(legi(i))).lt.0)QNcounter=QNcounter+1
		end do
		if(mod(QNcounter,2).eq.0)then
			REs=.true.
		else
			REs=.false.
		end if
		return
	end subroutine

	subroutine ParityParity(Res,dimen,ith,jth)
		integer,intent(inout)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		real*4,pointer::sQN(:)
		call Dimen%pointQN(sQN,ith)
		if(sQN(jth).gt.0)then
			Res=1
		else
			Res=-1
		end if
		return
	end subroutine
	subroutine ParityNewQN(NewQN,QN)
		real*4,intent(inout)::NewQN
		real*4,intent(in)::QN(:)
		integer::i,QNcounter
		QNcounter=0
		do i=1,size(QN)
			if(QN(i).lt.0)QNcounter=QNcounter+1
		end do
		if(mod(QNcounter,2).eq.0)then
			NewQN=1.
		else
			NewQN=-1.
		end if
		return
	end subroutine

	subroutine hermitian_conjugate_parity_dimension(dimen,legi)
		Type(dimension),intent(inout)::dimen
		integer,optional,intent(in)::legi
		return
	end subroutine

end module


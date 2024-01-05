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

module U1_Tools
	use Dimension_Tools
	use data_base_Tools
	use QuantumNumber_Tools 
	use Tools
	use memory_type
	implicit none
	private
	type(memory),private::WorkingMemory
	public::checkU1Rule
	public::U1RuleFunc
	public::U1ifParity
	public::U1Parity
	public::U1NewQN
	public::hermitian_conjugate_U1_dimension
	public::U1QN



contains
	function U1QN(QN,deg,arrow)
		type(QuanNum)::U1QN
		integer,intent(in)::deg(:)
		real*4,intent(in)::QN(:)
		integer,optional,intent(in)::arrow
		call U1QN%setQN(QN)
		call U1QN%setdeg(deg)
		if(present(arrow))then
			call U1QN%setFermiArrow(arrow)
		end if
		return
	end function

	subroutine checkU1Rule(D1,i1,D2,i2)
		type(Dimension),intent(in)::D1,D2
		integer,intent(in)::i1(:),i2(:)
		integer::i,j
		real*4,pointer::QN1(:),QN2(:)
		do i=1,size(i1)
			call D1%pointQN(QN1,i1(i))
			call D2%pointQN(QN2,i2(i))
			do j=1,size(QN1)
				if((QN1(j)+QN2(j)).nequ.0e0)then
					call writemess('ERROR in contract for the U(1) Symmetry Rule',-1)
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






	function U1RuleFunc(dimen,indices)
		logical::U1RuleFunc
		type(Dimension),intent(in)::dimen
		integer,intent(in)::indices(:)
		real*4,pointer::Q(:)
		integer::i
		real*4::temp
		call dimen%pointQN(Q,1)
		temp=Q(indices(1))
		do i=2,dimen%getRank()
			call dimen%pointQN(Q,i)
			temp=temp+Q(indices(i))
		end do
		U1RuleFunc=temp.equ.0e0
		return
	end function


	subroutine U1ifParity(REs,dimen,vec,legi)
		logical,intent(inout)::REs
		type(Dimension),intent(in)::dimen
		integer,intent(in)::vec(:),legi(:)
		real*4,pointer::sQN(:)
		integer::i,rank,sizelegi,intSumQN
		rank=dimen%getRank()
		sizelegi=size(legi)
		if(rank.ne.size(vec))then
			call writemess("ERROR in calculate the parity",-1)
			call error_stop
		end if
		intSumQN=0
		do i=1,sizelegi
			call Dimen%pointQN(sQN,legi(i))
			intSumQN=intSumQN+int(anint(sQN(vec(legi(i)))))
		end do
		if(mod(intSumQN,2).eq.0)then
			REs=.true.
		else
			REs=.false.
		end if
		return
	end subroutine

	subroutine U1Parity(Res,dimen,ith,jth)
		integer,intent(inout)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		integer::QN
		real*4,pointer::sQN(:)
		call Dimen%pointQN(sQN,ith)
		QN=int(anint(sQN(jth)))
		if(mod(QN,2).eq.0)then
			Res=1
		else
			Res=-1
		end if
		return
	end subroutine

	subroutine U1NewQN(NewQN,QN)
		real*4,intent(inout)::NewQN
		real*4,intent(in)::QN(:)
		integer::i
		NewQN=QN(1)
		do i=2,size(QN)
			NewQN=NewQN+QN(i)
		end do
		return
	end subroutine



	subroutine hermitian_conjugate_U1_dimension(dimen,legi)
		Type(dimension),intent(inout)::dimen
		integer,optional,intent(in)::legi
		integer::i,rank
		real*4,pointer::QN(:)
		if(present(legi))then
			call dimen%pointQN(QN,legi)
			QN=-QN
		else
			rank=dimen%getRank()
			do i=1,rank
				call dimen%pointQN(QN,i)
				QN=-QN
			end do
		end if
		return
	end subroutine




end module
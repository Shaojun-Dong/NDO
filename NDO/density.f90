

module Density_type
	use Tensor_Tools
	use Tools
	implicit none


	real*8,private::exp_factor=1d0
	real*8,private::max_norm_scal=1d300
	real*8,private::min_norm_scal=1d-300
	real*8,private::max_norm_scal_warning=1d300
	logical,private::NO_rescal_Flag=.true.

	interface GammaFuncPlus
		procedure GammaFuncPlus1
		procedure GammaFuncPlus2
	end interface

	interface GammaFuncMinus
		procedure GammaFuncMinus1
		procedure GammaFuncMinus2
	end interface

	interface PiFunc
		procedure PiFunc1
		procedure PiFunc2
	end interface

	interface rhoFunc0
		procedure rhoFunc0_
		procedure rhoFunc0_diag
	end interface
	interface rescal_subroutine
		procedure rescal_subroutine1
		procedure rescal_subroutine2
		procedure rescal_subroutine3
		procedure rescal_subroutine4
	end interface

contains

	subroutine rescal_subroutine1(a)
		real*8::a
		if(NO_rescal_Flag) return
		a=a*exp_factor
	end subroutine
	subroutine rescal_subroutine4(a)
		complex*16::a
		if(NO_rescal_Flag) return
		a=a*exp_factor
	end subroutine
	subroutine rescal_subroutine2(a)
		type(Tensor)::a
		if(NO_rescal_Flag) return
		a=a*exp_factor
	end subroutine
	subroutine rescal_subroutine3(a)
		type(Tensor)::a(:)
		integer::i
		if(NO_rescal_Flag) return
		do i=1,size(a)
			if(a(i)%getFlag()) a(i)=a(i)*exp_factor
		end do
		return
	end subroutine

	subroutine find_exp_factor(W1,W2,U1,U2,C1,C2,B1,B2,D)
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D
		logical::goon
		real*8::norm
		integer::counter
		norm=normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
		goon=norm.gt.max_norm_scal
		counter=0
		do while(goon)
			exp_factor=exp_factor*0.5
			norm=normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
			goon=norm.gt.max_norm_scal
			counter=counter+1
			if(counter.gt.100)then
				call writemess('ERROR in find_exp_factor',-1)
				call error_stop
			end if
			if(.not.goon)then
				call writemess('*************************************')
				call writemess('  Rescall the desity matrix,exp(N*x)')
				call writemess('  the rescal factor N='+exp_factor)
				call writemess('  norm='+norm)
				call writemess('*************************************')
			end if
		end do
		if(norm.le.min_norm_scal)then
			call writemess('ERROR, too small norm',-1)
			call error_stop
		end if
		return
	end subroutine

	function logWSigam(W,sigma,C)
		real*8::logWSigam
		type(Tensor),intent(in)::W,sigma,C
		integer::M,N,i
		type(Tensor)::temp
		if(W%getRank().ne.2)then
			call writemess('ERROR in logWSigam')
			call error_stop
		end if
		M=W%dim(1)
		N=W%dim(2)
		if(sigma%getTotalData().ne.N)then
			call writemess('ERROR in logWSigam,2')
			call error_stop
		end if
		if(C%getTotalData().ne.M)then
			call writemess('ERROR in logWSigam,3')
			call error_stop
		end if
		temp=W*sigma+C
		logWSigam=0d0
		do i=1,M
			logWSigam=logWSigam+log(1d0+exp(temp%di(i)))
		end do
		return
	end function

	function GammaFuncPlus2(W,C,B,S1,S2)result(GammaFunc)
		real*8::GammaFunc
		type(Tensor),intent(in)::W,C,B,S1,S2
		type(Tensor)::temp
		GammaFunc=logWSigam(W,S1,C)+logWSigam(W,S2,C)
		temp=B.dot.(S1+S2)
		GammaFunc=GammaFunc+temp
		GammaFunc=GammaFunc*0.5d0
		call rescal_subroutine(GammaFunc)
		return
	end function
	function GammaFuncPlus1(W,C,B,S1)result(GammaFunc)
		real*8::GammaFunc
		type(Tensor),intent(in)::W,C,B,S1
		type(Tensor)::temp
		GammaFunc=logWSigam(W,S1,C)
		temp=B.dot.S1
		GammaFunc=GammaFunc+temp
		call rescal_subroutine(GammaFunc)
		return
	end function
	function GammaFuncMinus2(W,C,B,S1,S2)result(GammaFunc)
		real*8::GammaFunc
		type(Tensor),intent(in)::W,C,B,S1,S2
		type(Tensor)::temp
		GammaFunc=logWSigam(W,S1,C)-logWSigam(W,S2,C)
		temp=B.dot.(S1-S2)
		GammaFunc=GammaFunc+temp
		GammaFunc=GammaFunc*0.5d0
		call rescal_subroutine(GammaFunc)
		return
	end function
	function GammaFuncMinus1(W,C,B,S1)result(GammaFunc)
		real*8::GammaFunc
		type(Tensor),intent(in)::W,C,B,S1
		type(Tensor)::temp
		GammaFunc=0d0
		return
	end function


	!PartialGammaFunc:
		!GammaFunc(1):  d GammaFunc^+
		!              -------------
		!                 d W 
		!
		!GammaFunc(2):  d GammaFunc^-
		!              -------------
		!                 d W 
		!
		!GammaFunc(3):  d GammaFunc^+
		!              -------------
		!                 d c 
		!
		!GammaFunc(4):  d GammaFunc^-
		!              -------------
		!                 d c 
		!
		!GammaFunc(5):  d GammaFunc^+
		!              -------------
		!                 d b 
		!
		!GammaFunc(6):  d GammaFunc^-
		!              -------------
		!                 d b 
		!

	function PartialGammaFunc(W,C,S1,S2)result(GammaFunc)
		type(Tensor),allocatable::GammaFunc(:)
		type(Tensor),intent(in)::W,C,S1,S2
		type(Tensor)::temp1,temp2,WsigmaC1,WsigmaC2
		integer::LD1,LD2,i,j,lenC,LenS
		complex*16,pointer::dp1(:,:),dp2(:,:),dp3(:),dp4(:),dp5(:),dp6(:)
		allocate(GammaFunc(6))
		LD1=W%dim(1)
		LD2=W%dim(2)
		lenC=C%getTotalData()
		LenS=S1%getTotalData()
		call GammaFunc(1)%allocate([LD1,LD2],'complex*16')
		call GammaFunc(1)%pointer(dp1)
		call GammaFunc(2)%allocate([LD1,LD2],'complex*16')
		call GammaFunc(2)%pointer(dp2)
		call GammaFunc(3)%allocate([lenC],'complex*16')
		call GammaFunc(3)%pointer(dp3)
		call GammaFunc(4)%allocate([lenC],'complex*16')
		call GammaFunc(4)%pointer(dp4)
		call GammaFunc(5)%allocate([LenS],'complex*16')
		call GammaFunc(5)%pointer(dp5)
		call GammaFunc(6)%allocate([LenS],'complex*16')
		call GammaFunc(6)%pointer(dp6)
		WsigmaC1=W*S1+C
		WsigmaC2=W*S2+C
		do i=1,LD1
			temp1=exp(WsigmaC1%di(i))
			temp1=temp1/(1d0+temp1)
			temp2=exp(WsigmaC2%di(i))
			temp2=temp2/(1d0+temp2)
			do j=1,LD2
				dp1(i,j)=( temp1*S1%di(j) ) + ( temp2*S2%di(j) )
				dp1(i,j)=dp1(i,j)*0.5d0
				dp2(i,j)=( temp1*S1%di(j) ) - ( temp2*S2%di(j) )
				dp2(i,j)=dp2(i,j)*0.5d0
			end do
			dp3(i)=(temp1+temp2)*0.5d0
			dp4(i)=(temp1-temp2)*0.5d0
		end do
		do i=1,LenS
			dp5(i)=(S1%di(i)+S2%di(i))*0.5d0
			dp6(i)=(S1%di(i)-S2%di(i))*0.5d0
		end do
		call rescal_subroutine(GammaFunc)
		return
	end function

	!PartialGammaFunc:
		!GammaFunc(1):  d GammaFunc^+
		!              -------------
		!                 d W 
		!
		!GammaFunc(2):  empty()
		!
		!GammaFunc(3):  d GammaFunc^+
		!              -------------
		!                 d c 
		!
		!GammaFunc(4):  empty()
		!
		!GammaFunc(5):  d GammaFunc^+
		!              -------------
		!                 d b 
		!
		!GammaFunc(6):  empty()
		!

	function PartialGammaPlusFunc(W,C,S1,S2)result(GammaFunc)
		type(Tensor),allocatable::GammaFunc(:)
		type(Tensor),intent(in)::W,C,S1,S2
		type(Tensor)::temp1,temp2,WsigmaC1,WsigmaC2
		integer::LD1,LD2,i,j,lenC,LenS
		complex*16,pointer::dp1(:,:),dp2(:,:),dp3(:),dp4(:),dp5(:),dp6(:)
		allocate(GammaFunc(6))
		LD1=W%dim(1)
		LD2=W%dim(2)
		lenC=C%getTotalData()
		LenS=S1%getTotalData()
		call GammaFunc(1)%allocate([LD1,LD2],'complex*16')
		call GammaFunc(1)%pointer(dp1)
		call GammaFunc(2)%empty()
		call GammaFunc(3)%allocate([lenC],'complex*16')
		call GammaFunc(3)%pointer(dp3)
		call GammaFunc(4)%empty()
		call GammaFunc(5)%allocate([LenS],'complex*16')
		call GammaFunc(5)%pointer(dp5)
		call GammaFunc(6)%empty()
		WsigmaC1=W*S1+C
		WsigmaC2=W*S2+C
		do i=1,LD1
			temp1=exp(WsigmaC1%di(i))
			temp1=(temp1/(1d0+temp1))*0.5d0
			temp2=exp(WsigmaC2%di(i))
			temp2=(temp2/(1d0+temp2))*0.5d0
			do j=1,LD2
				dp1(i,j)=( temp1*S1%di(j) ) + ( temp2*S2%di(j) )
			end do
			dp3(i)=(temp1+temp2)
		end do
		do i=1,LenS
			dp5(i)=(S1%di(i)+S2%di(i))*0.5d0
		end do
		call rescal_subroutine(GammaFunc)
		return
	end function


	!PartialGammaFunc:
		!GammaFunc(1):  empty()
		!
		!GammaFunc(2):  d GammaFunc^-
		!              -------------
		!                 d W 
		!
		!GammaFunc(3):  empty()
		!
		!GammaFunc(4):  d GammaFunc^-
		!              -------------
		!                 d c 
		!
		!GammaFunc(5): empty()
		!
		!GammaFunc(6):  d GammaFunc^-
		!              -------------
		!                 d b 
		!

	function PartialGammaMinusFunc(W,C,S1,S2)result(GammaFunc)
		type(Tensor),allocatable::GammaFunc(:)
		type(Tensor),intent(in)::W,C,S1,S2
		type(Tensor)::temp1,temp2,WsigmaC1,WsigmaC2
		integer::LD1,LD2,i,j,lenC,LenS
		complex*16,pointer::dp1(:,:),dp2(:,:),dp3(:),dp4(:),dp5(:),dp6(:)
		allocate(GammaFunc(6))
		LD1=W%dim(1)
		LD2=W%dim(2)
		lenC=C%getTotalData()
		LenS=S1%getTotalData()
		call GammaFunc(1)%empty()
		call GammaFunc(2)%allocate([LD1,LD2],'complex*16')
		call GammaFunc(2)%pointer(dp2)
		call GammaFunc(3)%empty()
		call GammaFunc(4)%allocate([lenC],'complex*16')
		call GammaFunc(4)%pointer(dp4)
		call GammaFunc(5)%empty()
		call GammaFunc(6)%allocate([LenS],'complex*16')
		call GammaFunc(6)%pointer(dp6)
		WsigmaC1=W*S1+C
		WsigmaC2=W*S2+C
		do i=1,LD1
			temp1=exp(WsigmaC1%di(i))
			temp1=(temp1/(1d0+temp1))*0.5d0
			temp2=exp(WsigmaC2%di(i))
			temp2=(temp2/(1d0+temp2))*0.5d0
			do j=1,LD2
				dp2(i,j)=( temp1*S1%di(j) ) - ( temp2*S2%di(j) )
			end do
			dp4(i)=(temp1-temp2)
		end do
		do i=1,LenS
			dp6(i)=(S1%di(i)-S2%di(i))*0.5d0
		end do
		call rescal_subroutine(GammaFunc)
		return
	end function


	function PiFunc2(U1,U2,D,S1,S2)result(PiFunc)
		complex*16::PiFunc
		type(Tensor),intent(in)::U1,U2,S1,S2,D
		type(Tensor)::US1,US2,S_plus_S,S_S,temp
		integer::M,N,i
		if(U1%getRank().ne.2)then
			call writemess('ERROR in PiFunc')
			call error_stop
		end if
		if(U2%getRank().ne.2)then
			call writemess('ERROR in PiFunc')
			call error_stop
		end if
		M=U1%dim(1)
		N=U1%dim(2)
		if(S1%getTotalData().ne.N)then
			call writemess('ERROR in PiFunc,2')
			call error_stop
		end if
		if(D%getTotalData().ne.M)then
			call writemess('ERROR in PiFunc,3')
			call error_stop
		end if


		S_plus_S=S1+S2
		S_S=S1-S2
		US1=(U1*S_plus_S)*0.5d0
		US2=(U2*S_S)*dcmplx(0d0,0.5d0)
		temp=US1+US2+D

		PiFunc=0d0
		do i=1,M
			PiFunc=PiFunc+log(1d0+exp(temp%zi(i)))
		end do
		call rescal_subroutine(PiFunc)
	end function

	function PiFunc1(U1,U2,D,S1)result(PiFunc)
		real*8::PiFunc
		type(Tensor),intent(in)::U1,U2,S1,D
		type(Tensor)::US1,temp
		integer::M,N,i
		if(U1%getRank().ne.2)then
			call writemess('ERROR in PiFunc')
			call error_stop
		end if
		if(U2%getRank().ne.2)then
			call writemess('ERROR in PiFunc')
			call error_stop
		end if
		M=U1%dim(1)
		N=U1%dim(2)
		if(S1%getTotalData().ne.N)then
			call writemess('ERROR in PiFunc,2')
			call error_stop
		end if
		if(D%getTotalData().ne.M)then
			call writemess('ERROR in PiFunc,3')
			call error_stop
		end if


		US1=U1*S1
		temp=US1+D

		PiFunc=0d0
		do i=1,M
			PiFunc=PiFunc+log(1d0+exp(temp%di(i)))
		end do
		call rescal_subroutine(PiFunc)
	end function

	!PartialPiFunc:
		!PiFunc(1):  d PiFunc
		!              -------------
		!                 d U_\lambeda 
		!
		!PiFunc(2):  d PiFunc
		!              -------------
		!                 d U_\mu  
		!
		!PiFunc(3):  d PiFunc
		!              -------------
		!                 d d
		!

	function PartialPiFunc(U1,U2,D,S1,S2)result(PiFunc)
		type(Tensor),allocatable::PiFunc(:)
		type(Tensor),intent(in)::U1,U2,S1,S2,D
		type(Tensor):: US1,US2,S_plus_S,S_S,temp
		integer::LD1,LD2,i,j,lenD
		complex*16,pointer::dp1(:,:),dp2(:,:),dp3(:)
		complex*16::znum,znum2
		allocate(PiFunc(3))
		LD1=U1%dim(1)
		LD2=U1%dim(2)
		lenD=D%getTotalData()
		call PiFunc(1)%allocate([LD1,LD2],'complex*16')
		call PiFunc(2)%allocate([LD1,LD2],'complex*16')
		call PiFunc(3)%allocate([lenD],'complex*16')
		call PiFunc(1)%pointer(dp1)
		call PiFunc(2)%pointer(dp2)
		call PiFunc(3)%pointer(dp3)
		S_plus_S=S1+S2
		S_S=S1-S2
		US1=(U1*S_plus_S)*0.5d0
		US2=(U2*S_S)*dcmplx(0d0,0.5d0)
		temp=US1+US2+D
		do i=1,LD1
			znum=exp(temp%zi(i))
			dp3(i)=znum/(1d0+znum)
			znum=dp3(i)*0.5d0
			znum2=znum*dcmplx(0d0,1d0)
			do j=1,LD2
				dp1(i,j)=znum*S_plus_S%di(j)  !(S1%di(j)+S2%di(j))
				dp2(i,j)=znum2*S_S%di(j)  !(S1%di(j)-S2%di(j))
			end do
			
		end do
		call rescal_subroutine(PiFunc)
		return
	end function


	function rhoFunc0_(W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2)
		complex*16::rhoFunc0_
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2
		rhoFunc0_=GammaFuncPlus(W1,C1,B1,S1,S2)+Dcmplx(0d0,1d0)*GammaFuncMinus(W2,C2,B2,S1,S2)+PiFunc(U1,U2,D,S1,S2)
		rhoFunc0_=exp(rhoFunc0_)
		return
	end function

	function rhoFunc0_diag(W1,W2,U1,U2,C1,C2,B1,B2,D,S1)
		real*8::rhoFunc0_diag
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D,S1
		rhoFunc0_diag=GammaFuncPlus(W1,C1,B1,S1)+PiFunc(U1,U2,D,S1)
		rhoFunc0_diag=exp(rhoFunc0_diag)
		return
	end function

	function normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
		real*8::normRho
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D
		integer::i,lenS
		type(Tensor)::S
		normRho=0d0
		lenS=W1%dim(2)
		call S%allocate([lenS],'real*8')
		call S%zero()
		do i=1,lenS
			call S%zero()
			call S%setValue(i,i)
			normRho=normRho+rhoFunc0(W1,W2,U1,U2,C1,C2,B1,B2,D,S)
		end do
		return
	end function


	function RhoFunc(W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2)
		complex*16::RhoFunc
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2
		real*8::norm
		norm=normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
		RhoFunc=rhoFunc0(W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2)/norm
		return
	end function

	function DensityMatrix(W1,W2,U1,U2,C1,C2,B1,B2,D)
		type(Tensor)::DensityMatrix
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D
		complex*16::rhoij
		integer::i,j,lenS
		type(Tensor)::S1,S2
		real*8::norm,OneOverNorm
		lenS=W1%dim(2)
		call DensityMatrix%empty()
		call DensityMatrix%allocate([lenS,lenS],'complex*16')
		norm=normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
		OneOverNorm=1d0/norm
		call S1%allocate([lenS],'real*8')
		call S2%allocate([lenS],'real*8')
		do i=1,lenS
			call S1%zero()
			call S1%setValue(i,i)
			do j=1,lenS
				call S2%zero()
				call S2%setValue(j,j)
				rhoij=rhoFunc0(W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2)
				call DensityMatrix%setValue([i,j],rhoij)
			end do
		end do
		DensityMatrix=DensityMatrix*OneOverNorm
		return
	end function


	!PartialGammaFunc:
		!Res(1): Rho(S1,S2)
		!
		!Res(2):  d Gamma^+(S1,S2)
		!         ---------------  * Rho(S1,S2)
		!           d W(i,j)
		!
		!Res(3):  d Gamma^-(S1,S2)
		!         ---------------  * Rho(S1,S2) *i
		!           d W(i,j)
		!
		!Res(4):  d  Pi(S1,S2)
		!         ---------------  * Rho(S1,S2)
		!         d U_lambda(i,j)
		!
		!Res(5):  d  Pi(S1,S2)
		!         ---------------  * Rho(S1,S2) *i
		!           d U_mu(i,j)
		!
		!Res(6):  d Gamma^+(S1,S2)
		!         ---------------  * Rho(S1,S2)
		!           d C(i)
		!
		!Res(7):  d Gamma^-(S1,S2)
		!         ---------------  * Rho(S1,S2) *i
		!           d C(i)
		!
		!Res(8):  d Gamma^+(S1,S2)
		!         ---------------  * Rho(S1,S2)
		!           d b(i)
		!
		!Res(6):  d Gamma^-(S1,S2)
		!         ---------------  * Rho(S1,S2) *i
		!           d b(i)
		!
		!Res(10):  d Pi(S1,S2)
		!         ---------------  * Rho(S1,S2)
		!           d d(i)
		!

	function DensityMatrixAndGradient(W1,W2,U1,U2,C1,C2,B1,B2,D)result(Res)
		type(Tensor),allocatable::Res(:)
		type(Tensor),intent(in)::W1,W2,U1,U2,C1,C2,B1,B2,D
		type(Tensor)::PartialGamma(6),PartialPi(3)
		complex*16::rhoij
		integer::i,j,k,lenS
		type(Tensor)::S1,S2
		real*8::norm,OneOverNorm
		allocate(Res(10))
		lenS=W1%dim(2)
		call Res(1)%empty()
		call Res(1)%allocate([lenS,lenS],'complex*16')
		call Res(1)%setName(1,'rho.S1')
		call Res(1)%setName(2,'rho.S2')
		call Res(2)%allocate([lenS,lenS,W1%dim(1),W1%dim(2)],'complex*16')!W1
		call Res(2)%setName(1,'A.S1')
		call Res(2)%setName(2,'A.S2')
		call Res(2)%setName(3,'A.i')
		call Res(2)%setName(4,'A.j')
		call Res(3)%allocate([lenS,lenS,W2%dim(1),W2%dim(2)],'complex*16')!W2
		call Res(3)%setName(1,'A.S1')
		call Res(3)%setName(2,'A.S2')
		call Res(3)%setName(3,'A.i')
		call Res(3)%setName(4,'A.j')
		call Res(4)%allocate([lenS,lenS,U1%dim(1),U1%dim(2)],'complex*16')!U1
		call Res(4)%setName(1,'A.S1')
		call Res(4)%setName(2,'A.S2')
		call Res(4)%setName(3,'A.i')
		call Res(4)%setName(4,'A.j')
		call Res(5)%allocate([lenS,lenS,U2%dim(1),U2%dim(2)],'complex*16')!U2
		call Res(5)%setName(1,'A.S1')
		call Res(5)%setName(2,'A.S2')
		call Res(5)%setName(3,'A.i')
		call Res(5)%setName(4,'A.j')
		call Res(6)%allocate([lenS,lenS,C1%getTotalData()],'complex*16')!C1
		call Res(6)%setName(1,'A.S1')
		call Res(6)%setName(2,'A.S2')
		call Res(6)%setName(3,'A.j')
		call Res(7)%allocate([lenS,lenS,C2%getTotalData()],'complex*16')!C2
		call Res(7)%setName(1,'A.S1')
		call Res(7)%setName(2,'A.S2')
		call Res(7)%setName(3,'A.i')
		call Res(8)%allocate([lenS,lenS,B1%getTotalData()],'complex*16')!B1
		call Res(8)%setName(1,'A.S1')
		call Res(8)%setName(2,'A.S2')
		call Res(8)%setName(3,'A.i')
		call Res(9)%allocate([lenS,lenS,B2%getTotalData()],'complex*16')!B2
		call Res(9)%setName(1,'A.S1')
		call Res(9)%setName(2,'A.S2')
		call Res(9)%setName(3,'A.i')
		call Res(10)%allocate([lenS,lenS,D%getTotalData()],'complex*16')!D
		call Res(10)%setName(1,'A.S1')
		call Res(10)%setName(2,'A.S2')
		call Res(10)%setName(3,'A.i')
		norm=normRho(W1,W2,U1,U2,C1,C2,B1,B2,D)
		if(norm.gt.max_norm_scal_warning)then
			call writemess('the trace of the density matrix is too large, trace=')
			call writemess(norm)
			call find_exp_factor(W1,W2,U1,U2,C1,C2,B1,B2,D)
		end if
		OneOverNorm=1d0/norm
		call S1%allocate([lenS],'real*8')
		call S2%allocate([lenS],'real*8')
		do i=1,lenS
			call S1%zero()
			call S1%setValue(i,i)
			do j=1,lenS
				call S2%zero()
				call S2%setValue(j,j)
				rhoij=rhoFunc0(W1,W2,U1,U2,C1,C2,B1,B2,D,S1,S2)
				rhoij=rhoij*OneOverNorm
				call Res(1)%setValue([i,j],rhoij)

				PartialGamma=PartialGammaPlusFunc(W1,C1,S1,S2)
				call addvaluedim4(Res(2),i,j,PartialGamma(1)*rhoij)!W1
				call addvaluedim3(Res(6),i,j,PartialGamma(3)*rhoij)!C1
				call addvaluedim3(Res(8),i,j,PartialGamma(5)*rhoij)!B1
				PartialGamma=PartialGammaMinusFunc(W2,C2,S1,S2)
				call addvaluedim4(Res(3),i,j,PartialGamma(2)*dcmplx(0d0,1d0)*rhoij)!W2
				call addvaluedim3(Res(7),i,j,PartialGamma(4)*dcmplx(0d0,1d0)*rhoij)!C2
				call addvaluedim3(Res(9),i,j,PartialGamma(6)*dcmplx(0d0,1d0)*rhoij)!B2
				PartialPi=PartialPiFunc(U1,U2,D,S1,S2)
				call addvaluedim4(Res(4),i,j,PartialPi(1)*rhoij)!U1
				call addvaluedim4(Res(5),i,j,PartialPi(2)*rhoij)!U2
				call addvaluedim3(Res(10),i,j,PartialPi(3)*rhoij)!D
			end do
		end do
		return
	end function

	subroutine addvaluedim4(A,S1,S2,B)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(in)::B
		integer::S1,S2
		complex*16,pointer::zp4(:,:,:,:),zp2(:,:)
		call A%pointer(zp4)
		call B%pointer(zp2)
		zp4(S1,S2,:,:)=zp2(:,:)
		return
	end subroutine
	subroutine addvaluedim3(A,S1,S2,B)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(in)::B
		integer::S1,S2
		complex*16,pointer::zp3(:,:,:),zp1(:)
		call A%pointer(zp3)
		call B%pointer(zp1)
		zp3(S1,S2,:)=zp1(:)
		return
	end subroutine

	!Test for PartialGammaFunc


	function GraGammaPlusMatrix(inW,inC,inB,S1,S2)
		type(Tensor),allocatable::GraGammaPlusMatrix(:)
		type(Tensor),intent(in)::inW,inC,inB,S1,S2
		real*8::delta
		type(Tensor)::W,C,B
		complex*16::Gamma0,Gamma
		integer::i,j
		complex*16,pointer::zp(:,:),zp1(:)
		real*8,pointer::dp(:,:),indp(:,:),dp1(:),indp1(:)
		allocate(GraGammaPlusMatrix(3))
		delta=0.0005
		W=inW
		Gamma0=GammaFuncPlus2(inW,inC,inB,S1,S2)
		call GraGammaPlusMatrix(1)%allocate([inW%dim(1),inW%dim(2)],'complex*16')
		call GraGammaPlusMatrix(1)%pointer(zp)
		W=inW
		C=inC
		B=inB
		call W%pointer(dp)
		call inW%pointer(indp)
		do i=1,inW%dim(1)
			do j=1,inW%dim(2)
				dp=indp
				dp(i,j)=dp(i,j)+delta
				zp(i,j)=(GammaFuncPlus2(W,C,B,S1,S2)-Gamma0)/delta
			end do
		end do

		call GraGammaPlusMatrix(2)%allocate([inC%getTotalData()],'complex*16')
		call GraGammaPlusMatrix(2)%pointer(zp1)
		call C%pointer(dp1)
		call inC%pointer(indp1)
		W=inW
		C=inC
		B=inB
		do i=1,C%getTotalData()
			dp1=indp1
			dp1(i)=dp1(i)+delta
			zp1(i)=(GammaFuncPlus2(W,C,B,S1,S2)-Gamma0)/delta
		end do

		call GraGammaPlusMatrix(3)%allocate([inB%getTotalData()],'complex*16')
		call GraGammaPlusMatrix(3)%pointer(zp1)
		call B%pointer(dp1)
		call inB%pointer(indp1)
		W=inW
		C=inC
		B=inB
		do i=1,B%getTotalData()
			dp1=indp1
			dp1(i)=dp1(i)+delta
			zp1(i)=(GammaFuncPlus2(W,C,B,S1,S2)-Gamma0)/delta
		end do
		return
	end function

	function GraGammaMunisMatrix(inW,inC,inB,S1,S2)
		type(Tensor),allocatable::GraGammaMunisMatrix(:)
		type(Tensor),intent(in)::inW,inC,inB,S1,S2
		real*8::delta
		type(Tensor)::W,C,B
		complex*16::Gamma0,Gamma
		integer::i,j
		complex*16,pointer::zp(:,:),zp1(:)
		real*8,pointer::dp(:,:),indp(:,:),dp1(:),indp1(:)
		allocate(GraGammaMunisMatrix(3))
		delta=0.0005
		W=inW
		Gamma0=GammaFuncMinus2(inW,inC,inB,S1,S2)
		call GraGammaMunisMatrix(1)%allocate([inW%dim(1),inW%dim(2)],'complex*16')
		call GraGammaMunisMatrix(1)%pointer(zp)
		W=inW
		C=inC
		B=inB
		call W%pointer(dp)
		call inW%pointer(indp)
		do i=1,inW%dim(1)
			do j=1,inW%dim(2)
				dp=indp
				dp(i,j)=dp(i,j)+delta
				zp(i,j)=(GammaFuncMinus2(W,C,B,S1,S2)-Gamma0)/delta
			end do
		end do

		call GraGammaMunisMatrix(2)%allocate([inC%getTotalData()],'complex*16')
		call GraGammaMunisMatrix(2)%pointer(zp1)
		call C%pointer(dp1)
		call inC%pointer(indp1)
		W=inW
		C=inC
		B=inB
		do i=1,C%getTotalData()
			dp1=indp1
			dp1(i)=dp1(i)+delta
			zp1(i)=(GammaFuncMinus2(W,C,B,S1,S2)-Gamma0)/delta
		end do

		call GraGammaMunisMatrix(3)%allocate([inB%getTotalData()],'complex*16')
		call GraGammaMunisMatrix(3)%pointer(zp1)
		call B%pointer(dp1)
		call inB%pointer(indp1)
		W=inW
		C=inC
		B=inB
		do i=1,B%getTotalData()
			dp1=indp1
			dp1(i)=dp1(i)+delta
			zp1(i)=(GammaFuncMinus2(W,C,B,S1,S2)-Gamma0)/delta
		end do
		return
	end function

	function GraPiFunc(inU1,inU2,inD,S1,S2)
		type(Tensor),allocatable::GraPiFunc(:)
		type(Tensor),intent(in)::inU1,inU2,inD,S1,S2
		real*8::delta
		type(Tensor)::U1,U2,D
		complex*16::Gamma0,Gamma
		integer::i,j
		complex*16,pointer::zp(:,:),zp1(:)
		real*8,pointer::dp(:,:),indp(:,:),dp1(:),indp1(:)
		allocate(GraPiFunc(3))
		delta=0.0005
		Gamma0=PiFunc2(inU1,inU2,inD,S1,S2)
		call GraPiFunc(1)%allocate([inU1%dim(1),inU1%dim(2)],'complex*16')
		call GraPiFunc(1)%pointer(zp)
		U1=inU1
		U2=inU2
		D=inD
		call U1%pointer(dp)
		call inU1%pointer(indp)
		do i=1,inU1%dim(1)
			do j=1,inU1%dim(2)
				dp=indp
				dp(i,j)=dp(i,j)+delta
				zp(i,j)=(PiFunc2(U1,U2,D,S1,S2)-Gamma0)/delta
			end do
		end do

		call GraPiFunc(2)%allocate([inU2%dim(1),inU2%dim(2)],'complex*16')
		call GraPiFunc(2)%pointer(zp)
		U1=inU1
		U2=inU2
		D=inD
		call U2%pointer(dp)
		call inU2%pointer(indp)
		do i=1,inU2%dim(1)
			do j=1,inU2%dim(2)
				dp=indp
				dp(i,j)=dp(i,j)+delta
				zp(i,j)=(PiFunc2(U1,U2,D,S1,S2)-Gamma0)/delta
			end do
		end do

		call GraPiFunc(3)%allocate([inD%dim(1)],'complex*16')
		call GraPiFunc(3)%pointer(zp1)
		U1=inU1
		U2=inU2
		D=inD
		call D%pointer(dp1)
		call inD%pointer(indp1)
		do i=1,inD%dim(1)
				dp1=indp1
				dp1(i)=dp1(i)+delta
				zp1(i)=(PiFunc2(U1,U2,D,S1,S2)-Gamma0)/delta
		end do

	end function

end MODULE 
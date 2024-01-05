module CG_For_LLS_tools
	use Tensor_Tools
	implicit none
	private
	!use for the solve LLS with CG method
	!
	! AX=b, A is a n times n positive definite matrix
	!
	! Use for the two situations below
	!
	! 1.
	!
	!   A=Y' * G * Y,  Y' is the transposition of Y
	!
	!    Y is N time n matrix, N may not be equation to n
	!
	! 2. 
	!   A=a1 * ( Z1 *  Y1 ) + a2 * (Z2 * Y2) + ...
	!    Yi are N time n matrix, Zi are n times N matrix, ai are number
	!
	!      N may not be equation to n
	!
	!   1) set_CGLLS_ZYp set the subroutine to return
	!        (a1 * ( Z1 *  Y1 ) + a2 * (Z2 * Y2) )*x
	!       Yi are N time n matrix, Zi are n times N matrix, ai are number
	!       and x is a vector of length=n
	!   2) set_CGLLS_py set the subroutine to return
	!             p.dot.y, both p and y are vectors with length=n
	!   3) set_CGLLS_norm2 set the subroutine to return
	!             |x|^2, x is a vector with length=n
	!   4) set_CGLLS_xalphap set the subroutine to return
	!       if input output_in_x=.true.
	!             x=x+alpha*p
	!       else output_in_x =.false.   
	!             p=x+alpha*p
	!      x and p are vectors with length=n and alpha is real*8

	interface
		subroutine ZY_times_p_interface(ai,Z,Y,p,res)
			import :: Tensor
			type(Tensor),target,intent(inout)::Z(:),Y(:),p,res
			real*8,intent(in)::ai(:)
		end subroutine ZY_times_p_interface
	end interface
	procedure(ZY_times_p_interface),pointer,private::ZY_times_p=>default_ZY_times_p

	interface
		subroutine p_dot_y_interface(p1,p2,res)
			import :: Tensor
			type(Tensor),target,intent(inout)::p1,p2
		real*8::res
		end subroutine p_dot_y_interface
	end interface
	procedure(p_dot_y_interface),pointer,private::p_dot_y=>default_p_dot_y

	interface
		subroutine p_norm2_interface(p,res)
			import :: Tensor
			type(Tensor),target,intent(inout)::p
			real*8::res
		end subroutine p_norm2_interface
	end interface
	procedure(p_norm2_interface),pointer,private::p_norm2=>default_p_norm2


	interface
		subroutine x_plus_alpha_p_interface(inoutx,alpha,p,output_in_x)
			import :: Tensor
			type(Tensor),target,intent(inout)::inoutx,p
			logical,intent(in)::output_in_x
			real*8,intent(in)::alpha
		end subroutine x_plus_alpha_p_interface
	end interface
	procedure(x_plus_alpha_p_interface),pointer,private::x_plus_alpha_p=>default_x_plus_alpha_p


	interface
		subroutine matrix_time_vector_interface(outAp,inp)!outAp=A*inp
			import :: Tensor
			type(Tensor),target,intent(inout)::outAp
			type(Tensor),target,intent(in)::inp
		end subroutine matrix_time_vector_interface
	end interface


	public::runCGLLS
	interface runCGLLS
		module procedure runCGLLS1
		module procedure runCGLLS2
		module procedure runCGLLS3
	end interface

	public::set_CGLLS_ZYp,set_CGLLS_py,set_CGLLS_norm2,set_CGLLS_xalphap,matrix_time_vector_interface
contains
	subroutine set_CGLLS_ZYp(Func)
		procedure(ZY_times_p_interface)::Func
		call writemess('Set the subroutine in CG_For_LLS_tools')
		call writemess(' return Z*Y*p')
		ZY_times_p=>Func
	end subroutine
	subroutine set_CGLLS_py(Func)
		procedure(p_dot_y_interface)::Func
		call writemess('Set the subroutine in CG_For_LLS_tools')
		call writemess(' return py')
		p_dot_y=>Func
	end subroutine
	subroutine set_CGLLS_norm2(Func)
		procedure(p_norm2_interface)::Func
		call writemess('Set the subroutine in CG_For_LLS_tools')
		call writemess(' return |p|^2')
		p_norm2=>Func
	end subroutine
	subroutine set_CGLLS_xalphap(Func)
		procedure(x_plus_alpha_p_interface)::Func
		call writemess('Set the subroutine in CG_For_LLS_tools')
		call writemess(' return x=x+alpha*p or p=x+alpha*p')
		x_plus_alpha_p=>Func
	end subroutine

	subroutine default_ZY_times_p(ai,Z,Y,p,res)
		type(Tensor),target,intent(inout)::Z(:),Y(:),p,res
		real*8,intent(in)::ai(:)
		integer::i
		type(Tensor)::TMP
		call res%empty
		do i=1,size(Z)
			TMP=Y(i)*p
			TMP=Z(i)*TMP
			res=res+(TMP*ai(i))
		end do
		
		return
	end subroutine

	subroutine default_p_dot_y(p1,p2,res)
		type(Tensor),target,intent(inout)::p1,p2
		real*8::res
		res=p1.dot.p2
		return
	end subroutine

	subroutine default_p_norm2(p,res)
		type(Tensor),target,intent(inout)::p
		real*8::res
		res=p%norm2()
		return
	end subroutine

	subroutine default_x_plus_alpha_p(inoutx,alpha,inputp,output_in_x)
		type(Tensor),target,intent(inout)::inoutx,inputp
		real*8,intent(in)::alpha
		logical,intent(in)::output_in_x
		if(output_in_x)then
			inoutx=inoutx+(alpha*inputp)
		else
			inputp=inoutx+(alpha*inputp)
		end if
		return
	end subroutine


	subroutine runOneStep(inoutx,inoutp,inoutr,ai,Z,Y,memoryAp)
		type(Tensor)::inoutx,inoutp,inoutr,Z(:),Y(:),memoryAp
		real*8::ai(:)
		real*8::rr,alpha,beta
		call p_norm2(inoutr,rr)!rr= r*r
		call ZY_times_p(ai,Z,Y,inoutp,memoryAp)!Ap=A*p
		call p_dot_y(inoutp,memoryAp,alpha)!alpha=p*A*p
		alpha=rr/alpha                         !alpha=r*r/p*A*p
		call x_plus_alpha_p(inoutx,alpha,inoutp,.true.)!x=x+alpha*p
		call x_plus_alpha_p(inoutr,-alpha,memoryAp,.true.)!r=r-alpha*A*P
		call p_norm2(inoutr,beta)             
		beta=beta/rr                                
		call x_plus_alpha_p(inoutr,beta,inoutp,.false.)!p=r+beta*P
		return
	end subroutine

	subroutine inital_r0(r0,x0,b,ai,Z,Y,memoryAp)
		type(Tensor)::r0,x0,Z(:),Y(:),memoryAp,b
		real*8::ai(:)
		r0=b
		call ZY_times_p(ai,Z,Y,x0,memoryAp)
		call x_plus_alpha_p(r0,-1d0,memoryAp,.true.)
		return
	end subroutine

	subroutine runCGLLS1(num,err,inoutx,b,ai,Z,Y,memoryAp,memory_r,memory_p)
		type(Tensor)::inoutx,b,Z(:),Y(:),memoryAp,memory_r,memory_p
		integer,intent(in)::num
		real*8,intent(in)::ai(:),err
		integer::i
		real*8::norm2
		if(.not.inoutx%getFlag())then
			call inoutx%allocate([b%getTotalData()],'real*8')
			call inoutx%random([-1d0,1d0])
		end if
		call inital_r0(memory_r,inoutx,b,ai,Z,Y,memoryAp)
		memory_p=memory_r
		do i=1,num
			call runOneStep(inoutx,memory_p,memory_r,ai,Z,Y,memoryAp)
			call p_norm2(memory_p,norm2)
			if(norm2.le.err)exit
		end do
		return
	end subroutine


	subroutine YGYp(Y,G,p,res)
		type(Tensor),target,intent(inout)::Y,G,p,res
		res=(G*(Y*p))*Y
		return

		res=Y*p
		res=G*res
		res=res*Y
		return
	end subroutine

	subroutine runCGLLS2(num,err,inoutx,b,Y,G,memoryAp,memory_r,memory_p)
		type(Tensor)::inoutx,b,Y,G,memoryAp,memory_r,memory_p
		integer,intent(in)::num
		real*8,intent(in)::err
		integer::i
		real*8::norm2,rr,alpha,beta
		if(.not.inoutx%getFlag())then
			call inoutx%allocate([b%getTotalData()],'real*8')
			call inoutx%random([-1d0,1d0])
		end if
		call YGYp(Y,G,inoutx,memory_r)
		memory_r=b-memory_r
		memory_p=memory_r
		do i=1,num
			rr=memory_r%dnorm2()
			call YGYp(Y,G,memory_p,memoryAp)
			alpha=rr/(memory_p.dot.memoryAp)
			inoutx=inoutx+(alpha*memory_p)
			memory_r=memory_r-(alpha*memoryAp)
			beta=memory_r%dnorm2()/rr
			memory_p=memory_r+(beta*memory_p)
			norm2=memory_p%dnorm2()
			if(norm2.le.err)exit
		end do
		return
	end subroutine

	subroutine runCGLLS3(num,err,inoutx,b,FuncAp,memoryAp,memory_r,memory_p)
		type(Tensor)::inoutx,b,memoryAp,memory_r,memory_p
		procedure(matrix_time_vector_interface)::FuncAp
		integer,intent(in)::num
		real*8,intent(in)::err
		integer::i
		real*8::norm2,rr,alpha,beta
		if(.not.inoutx%getFlag())then
			call inoutx%allocate([b%getTotalData()],'real*8')
			call inoutx%random([-1d0,1d0])
		end if
		call FuncAp(memory_r,inoutx)
		memory_r=b-memory_r
		memory_p=memory_r
		do i=1,num
			rr=memory_r%dnorm2()
			call FuncAp(memoryAp,memory_p)
			alpha=rr/(memory_p.dot.memoryAp)
			inoutx=inoutx+(alpha*memory_p)
			memory_r=memory_r-(alpha*memoryAp)
			beta=memory_r%dnorm2()/rr
			memory_p=memory_r+(beta*memory_p)
			norm2=memory_p%dnorm2()
			if(norm2.le.err)exit
		end do
		return
	end subroutine
end module

module LineSearchEnergy
	use Tensor_Tools
	use Tools
	implicit none

	type(Tensor),private::M,abc,y
	real*8::zero_a=1d-2
	real*8::zero_x=1d-5

	!input fi are the energies, do not use the gradients
	!
	!

contains

	!input x1>0
	!if f1 > f0, return  (x1+x0)/2 
	!
	!    f 
	!  /|
	!   |         !      x1
	!   |         ! 
	!   |         !
	!   |  x0     !
	!   |_________!_____________\ x
	!             x2
	!
	!if f1 < f0 , return  x1+(x1-x0)
	!
	!    f 
	!  /|
	!   |  x0        !     
	!   |       x1   !
	!   |            !
	!   |            !  
	!   |____________!__________\ x
	!                x2
	!

	function find_the_second_point_(inx,inf)result(Res)
		real*8::Res
		real*8::inx(3),inf(3)
		real*8::x0,x1,f0,f1
		x0=inx(1)
		x1=inx(2)
		f0=inf(1)
		f1=inf(2)
		if(f1.ge.f0)then
			Res=(x1+x0)/2d0
		else
			Res=x1+(x1-x0)
		end if
		return
	end function

	!
		!input x1,f1=f(x1), f0=f(0),g0=g(0)
		!
		! set    f(x)=a*x*x + b x +c
		!   then g(x)=2*a*x + b
		!
		!   f0,g0 ==>   c=f0, b=g0
		!
		!   x1,f1 ==>  a = (f1-f0-g0*x1)/(x1*x1)
		!
		!  if a>zero_a  , x2 = -b/(a+a)
		!  if a<=zero_a , the point of x1 is too large, x2=x1/2  
		!

	function find_the_second_point(x1,f1,f0,g0)result(Res)
		real*8::Res
		real*8::x1,f1,f0,g0
		real*8::a
		a=(f1-f0-g0*x1)/(x1*x1)
		if(a.gt.zero_a)then
			Res=-g0/(a+a)
		else
			Res=x1/2d0
		end if
		if(Res.le.zero_x)then
			Res=x1*0.5d0
		end if
		return
	end function

	!
		!input x0<x1<x2
		!there are 6 situations
		! 1.
		!      f0 < f1 <f2
		!
		!    f 
		!  /|
		!   |                x2
		!   |         x1
		!   |
		!   |  x0
		!   |_______________________\ x
		!
		!   if x0>0
		!    x3=x0/2
		!   else
		!    x3=x1/2
		!
		!---------------------------------
		!
		!  x3= x2/2
		!
		! 2.
		!      f0 < f2 <f1
		!
		!    f 
		!  /|
		!   |         x1       
		!   |                x2
		!   |
		!   |  x0
		!   |_______________________\ x
		!
		!   if x0>0
		!    x3=x0/2
		!   else
		!    x3=x1/2
		!
		!---------------------------------
		!
		!
		! 3.
		!      f2 < f0 <f1
		!
		!    f 
		!  /|         x1
		!   |                  
		!   |  x0              
		!   |     
		!   |                x2
		!   |_______________________\ x
		!
		!
		!   if x0>0
		!    x3=x0/2
		!   else
		!    x3=x1/2
		!
		!
		!---------------------------------
		!
		!
		! 4.
		!      f1 < f0 <f2
		!
		!    f 
		!  /|
		!   |               x2 
		!   |                  
		!   |  x0
		!   |         x1
		!   |_______________________\ x
		!
		!
		!   quadratic function  to find the min point 
		!
		!
		!---------------------------------
		!
		!
		!
		! 5.
		!      f1 < f2 <f0
		!
		!    f 
		!  /|                
		!   |  x0              
		!   |                x2
		!   |     
		!   |         x1
		!   |_______________________\ x
		!
		!  
		!   quadratic function  to find the min point 
		!
		!---------------------------------
		!
		! 6.
		!      f2 < f1 <f0
		!
		!    f 
		!  /|
		!   |  x0              
		!   |        x1        
		!   |       
		!   |               x2
		!   |_______________________\ x
		!
		!  
		!  x3=x2+(x2-x0)
		!
	
	function find_the_next_point_(inx,inf)result(Res)
		real*8::Res
		real*8::inx(3),inf(3)
		real*8::x0,f0,x1,f1,x2,f2
		x0=inx(1)
		x1=inx(2)
		x2=inx(3)
		f0=inf(1)
		f1=inf(2)
		f2=inf(3)
		if((f0.le.f1).and.(f1.le.f2))then   !(f0 < f1 <f2) case 1
			if(x0.gt.zero_x)then
				Res=x0/2d0
			else
				Res=x1/2d0
			end if
		else if((f0.le.f2).and.(f2.le.f1))then   !f0 < f2 <f1 case 2
			if(x0.gt.zero_x)then
				Res=x0/2d0
			else
				Res=x1/2d0
			end if
		else if((f2.le.f0).and.(f0.le.f1))then   !f2 < f0 <f1 case 3
			if(x0.gt.zero_x)then
				Res=x0/2d0
			else
				Res=x1/2d0
			end if
		else if((f1.le.f0).and.(f0.le.f2))then   !f1 < f0 <f2 case 4
			call quadratic_function_([x0,x1,x2],[f0,f1,f2],Res)
		else if((f1.le.f2).and.(f2.le.f0))then   !f1 < f2 <f0 case 5
			call quadratic_function_([x0,x1,x2],[f0,f1,f2],Res)
		else if((f2.le.f1).and.(f1.le.f0))then   !f2 < f1 <f0 case 6
			Res=x2+(x2-x0)
		else
			call writemess('ERROPR in line search')
			call writemess('x0,x1,x2 are')
			call writemess(x0+','+x1+','+x2)
			call writemess('f0,f1,f2 are')
			call writemess(f0+','+f1+','+f2)
			call error_stop
		end if
		return
	end function


	subroutine quadratic_function_(xin,fin,newX)
		real*8::xin(3),fin(3),newX,newf
		real*8,pointer::Ap(:,:),yp(:)
		real*8::a,b,c,delta
		if(.not.M%getFlag())then
			call M%allocate([3,3],'real*8')
			call y%allocate([3],'real*8')
		end if
		call M%pointer(Ap)
		call y%pointer(yp)
		Ap(1,:)=[xin(1)*xin(1),xin(1),1d0]
		Ap(2,:)=[xin(2)*xin(2),xin(2),1d0]
		Ap(3,:)=[xin(3)*xin(3),xin(3),1d0]
		yp=fin
		call abc%SolveLLS(M,y,0d0)
		a=abc%di(1)
		b=abc%di(2)
		c=abc%di(3)
		if(a.gt.0)then
			newX=-b/(a+a)
		else
			call writemess('WRONNING in quadratic_function of line search',-1)
			call writemess('x0,x1,x2 are',-1)
			call writemess(xin(1)+','+xin(2)+','+xin(3),-1)
			call writemess('f0,f1,f2 are',-1)
			call writemess(fin(1)+','+fin(2)+','+fin(3),-1)
			call writemess('a in quadraticfunction is'+a,-1)
			newX=(xin(1)+xin(2))/2d0
		end if
		return
	end subroutine


	function find_the_next_point(inx,inf)result(Res)
		real*8::Res
		real*8::inx(3),inf(3)
		real*8::x0,f0,x1,f1,x2,f2
		logical::resetNewX
		x0=inx(1)
		x1=inx(2)
		x2=inx(3)
		f0=inf(1)
		f1=inf(2)
		f2=inf(3)
		if((f0.le.f1).and.(f1.le.f2))then   !(f0 < f1 <f2) case 1
			resetNewX=quadratic_function([x0,x1,x2],[f0,f1,f2],Res)
			if(resetNewX)then
				if(x0.gt.zero_x)then
					Res=x0/2d0
				else
					Res=x1/2d0
				end if
			end if
			if(Res.le.zero_x)then
				Res=x1*0.8
			end if
		else if((f0.le.f2).and.(f2.le.f1))then   !f0 < f2 <f1 case 2
			if(x0.gt.zero_x)then
				Res=x0/2d0
			else
				Res=x1/2d0
			end if
			if(Res.le.zero_x)then
				Res=x1*0.8
			end if
		else if((f2.le.f0).and.(f0.le.f1))then   !f2 < f0 <f1 case 3
			if(x0.gt.zero_x)then
				Res=x0/2d0
			else
				Res=x1/2d0
			end if
			if(Res.le.zero_x)then
				Res=x1*0.8
			end if
		else if((f1.le.f0).and.(f0.le.f2))then   !f1 < f0 <f2 case 4
			resetNewX=quadratic_function([x0,x1,x2],[f0,f1,f2],Res)
			if(resetNewX)then
				call writemess('WRONNING in quadratic_function of line search',-1)
				call writemess('x0,x1,x2 are',-1)
				call writemess(x0+','+x1+','+x2,-1)
				call writemess('f0,f1,f2 are',-1)
				call writemess(f0+','+f1+','+f2,-1)
				Res=(x1+x2)/2d0
			end if
		else if((f1.le.f2).and.(f2.le.f0))then   !f1 < f2 <f0 case 5
			resetNewX=quadratic_function([x0,x1,x2],[f0,f1,f2],Res)
			if(resetNewX)then
				call writemess('WRONNING in quadratic_function of line search',-1)
				call writemess('x0,x1,x2 are',-1)
				call writemess(x0+','+x1+','+x2,-1)
				call writemess('f0,f1,f2 are',-1)
				call writemess(f0+','+f1+','+f2,-1)
				Res=(x1+x2)/2d0
			end if
		else if((f2.le.f1).and.(f1.le.f0))then   !f2 < f1 <f0 case 6
			resetNewX=quadratic_function([x0,x1,x2],[f0,f1,f2],Res)
			if(resetNewX)then
				Res=x2+(x2-x0)
			end if
			if(Res.le.x2)then
				Res=x2+(x2-x0)
			end if
		else
			call writemess('ERROPR in line search')
			call writemess('x0,x1,x2 are')
			call writemess(x0+','+x1+','+x2)
			call writemess('f0,f1,f2 are')
			call writemess(f0+','+f1+','+f2)
			call error_stop
		end if
		return
	end function

	! if quadratic_function= .true. can not find newX, one should resetX after the calling

	logical function quadratic_function(xin,fin,newX)result(Res)
		real*8::xin(3),fin(3),newX,newf
		real*8,pointer::Ap(:,:),yp(:)
		real*8::a,b,c,delta
		if(.not.M%getFlag())then
			call M%allocate([3,3],'real*8')
			call y%allocate([3],'real*8')
		end if
		call M%pointer(Ap)
		call y%pointer(yp)
		Ap(1,:)=[xin(1)*xin(1),xin(1),1d0]
		Ap(2,:)=[xin(2)*xin(2),xin(2),1d0]
		Ap(3,:)=[xin(3)*xin(3),xin(3),1d0]
		yp=fin
		call abc%SolveLLS(M,y,0d0)
		a=abc%di(1)
		b=abc%di(2)
		c=abc%di(3)
		if(a.gt.zero_a)then
			newX=-b/(a+a)
			Res=.false.
		else
			Res=.true.
		end if
		return
	end function

	subroutine select_new_point(inoutx,inoutf,newx,newf)
		real*8,intent(inout)::inoutx(3),inoutf(3),newx,newf
		integer::maxi
		maxi=Find_max_out_index(inoutf)
		inoutx(maxi)=newx
		inoutf(maxi)=newf
		return
	end subroutine
	function Find_max_out_index(inf)result(maxi)
		real*8,intent(in)::inf(3)
		integer::maxi
		real*8::maxf
		maxi=1
		maxf=inf(1)
		if(inf(2).ge.maxf)then
			maxi=2
			maxf=inf(2)
		end if
		if(inf(3).ge.maxf)then
			maxi=3
			maxf=inf(3)
		end if
		return
	end function

	function Find_min_out_index(inf)result(mini)
		real*8,intent(in)::inf(3)
		integer::mini
		real*8::minf
		mini=1
		minf=inf(1)
		if(inf(2).le.minf)then
			mini=2
			minf=inf(2)
		end if
		if(inf(3).le.minf)then
			mini=3
			minf=inf(3)
		end if
		return
	end function
	subroutine reorderdata(allx,allf)
		real*8,intent(inout)::allx(3),allf(3)
		integer::indices(3)
		real*8::tmpf(3)
		call sortdata(allx,indices,.true.)
		tmpf=allf
		allf(1)=tmpf(indices(1))
		allf(2)=tmpf(indices(2))
		allf(3)=tmpf(indices(3))
		return
	end subroutine

end module

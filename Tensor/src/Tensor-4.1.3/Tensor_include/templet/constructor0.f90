	function constructor0FuncName(val,dimen)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::val
		integer,optional,intent(in)::dimen(:)
		Res%Data=[val]
		if(present(dimen))then
			call Res%resetdim(dimen)
		else
			Res%Dimension=[1]
		end if
		return
	end function
	function constructor1FuncName(val,dimen)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::val(:)
		integer,optional,intent(in)::dimen(:)
		Res%Data=val
		if(present(dimen))then
			call Res%resetdim(dimen)
		else
			Res%Dimension=[size(val)]
		end if
		return
	end function
	function constructor2FuncName(B,dimen)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::B(:,:)
		integer,optional,intent(in)::dimen(:)
		DATATYPE2,pointer::Bp(:),BP2(:,:)
		integer::m,n
		m=size(B,1)
		n=size(B,2)
		call Res%Data%allocate(size(B),DATATYPENumber)
		call Res%Data%pointer(Bp)
		BP2(1:m,1:n)=>BP
		BP2=B
		if(present(dimen))then
			call Res%resetdim(dimen)
		else
			Res%Dimension=[m,n]
		end if
		return
	end function

	function constructor3FuncName(B,dimen)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::B(:,:,:)
		integer,optional,intent(in)::dimen(:)
		DATATYPE2,pointer::Bp(:),BP2(:,:,:)
		integer::m,n,o
		m=size(B,1)
		n=size(B,2)
		o=size(B,3)
		call Res%Data%allocate(size(B),DATATYPENumber)
		call Res%Data%pointer(Bp)
		BP2(1:m,1:n,1:o)=>BP
		BP2=B
		if(present(dimen))then
			call Res%resetdim(dimen)
		else
			Res%Dimension=[m,n,o]
		end if
		return
	end function

	function constructor4FuncName(B,dimen)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::B(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		DATATYPE2,pointer::Bp(:),BP2(:,:,:,:)
		integer::m,n,o,p
		m=size(B,1)
		n=size(B,2)
		o=size(B,3)
		p=size(B,4)
		call Res%Data%allocate(size(B),DATATYPENumber)
		call Res%Data%pointer(Bp)
		BP2(1:m,1:n,1:o,1:p)=>BP
		BP2=B
		if(present(dimen))then
			call Res%resetdim(dimen)
		else
			Res%Dimension=[m,n,o,p]
		end if
		return
	end function

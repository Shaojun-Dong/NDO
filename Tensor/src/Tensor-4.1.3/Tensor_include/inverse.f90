
	!! The inverse of a matrix: the input tensor should be a square matrix 
		!A*X=E ==>X=A^-1

	function inverse(T)
		type(Tensor)::inverse
		type(Tensor),intent(in) :: T
		type(Tensor):: E,TMPT,TMPinv
		integer :: M,N,i,j
		real*4,pointer::sp(:,:)
		real*8,pointer::dp(:,:)
		complex*8,pointer::cp(:,:)
		complex*16,pointer::zp(:,:)
		if(T%getRank().ne.2) then
			call writemess("ERROR in calculating the inverse of a Tensor",-1)
			call writemess("input Tensor should be a square matrix",-1)
			call error_stop()
		endif
		if(.not.T%getSymmetryFlag())then
			M = T%dim(1)
			N = T%dim(2)
			if(M.ne.N) then
				call writemess("ERROR in calculating the inverse of a Tensor",-1)
				call writemess("input Tensor should be a square matrix",-1)
				call error_stop()
			endif
			E=eye(M,N,T%getType())
			inverse=linequ(T,E)
			call inverse%resetdim(T%Dimension)
			return
		end if
		call inverse%allocate(T,T%getType())
		select case(T%getType())
			case(2)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(sp,[i,j])
							TMPT=sp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(sp,[i,j])
							sp=TMPinv
						end if
					end do
				end do
			case(3)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(dp,[i,j])
							TMPT=dp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(dp,[i,j])
							dp=TMPinv
						end if
					end do
				end do
			case(4)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(cp,[i,j])
							TMPT=cp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(cp,[i,j])
							cp=TMPinv
						end if
					end do
				end do
			case(5)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(zp,[i,j])
							TMPT=zp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(zp,[i,j])
							zp=TMPinv
						end if
					end do
				end do
			case default
				call writemess('ERROR data type in inverse',-1)
				call error_stop
		end select
		return
	end function	
	function inverseTen(T,RCOND)result(inverse)
		type(Tensor)::inverse
		type(Tensor),intent(in) :: T
		class(*),intent(in)::RCOND
		type(Tensor):: E,TMPT,TMPinv
		integer :: M,N,i,j
		real*4,pointer::sp(:,:)
		real*8,pointer::dp(:,:)
		complex*8,pointer::cp(:,:)
		complex*16,pointer::zp(:,:)
		if(T%getRank().ne.2) then
			call writemess("ERROR in calculating the inverse of a Tensor",-1)
			call writemess("input Tensor should be a square matrix",-1)
			call error_stop()
		endif
		if(.not.T%getSymmetryFlag())then
			M = T%dim(1)
			N = T%dim(2)
			if(M.ne.N) then
				call writemess("ERROR in calculating the inverse of a Tensor",-1)
				call writemess("input Tensor should be a square matrix",-1)
				call error_stop()
			endif
			E=eye(M,N,T%getType())
			inverse=linequ(T,E,RCOND)
			call inverse%resetdim(T%dimension)
			return
		end if
		call inverse%allocate(T,T%getType())
		select case(T%getType())
			case(2)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(sp,[i,j])
							TMPT=sp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(sp,[i,j])
							sp=TMPinv
						end if
					end do
				end do
			case(3)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(dp,[i,j])
							TMPT=dp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(dp,[i,j])
							dp=TMPinv
						end if
					end do
				end do
			case(4)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(cp,[i,j])
							TMPT=cp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(cp,[i,j])
							cp=TMPinv
						end if
					end do
				end do
			case(5)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(zp,[i,j])
							TMPT=zp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(zp,[i,j])
							zp=TMPinv
						end if
					end do
				end do
			case default
				call writemess('ERROR data type in inverse',-1)
				call error_stop
		end select
		return
	end function	
	
	function inverseTensor(T,RCOND)result(inverse)
		type(Tensor)::inverse
		class(Tensor),intent(in) :: T
		class(*),intent(in),optional::RCOND
		type(Tensor):: E,TMPT,TMPinv
		integer :: M,N,i,j
		real*4,pointer::sp(:,:)
		real*8,pointer::dp(:,:)
		complex*8,pointer::cp(:,:)
		complex*16,pointer::zp(:,:)
		if(T%getRank().ne.2) then
			call writemess("ERROR in calculating the inverse of a Tensor",-1)
			call writemess("input Tensor should be a square matrix",-1)
			call error_stop()
		endif
		if(.not.T%getSymmetryFlag())then
			M = T%dim(1)
			N = T%dim(2)
			if(M.ne.N) then
				call writemess("ERROR in calculating the inverse of a Tensor",-1)
				call writemess("input Tensor should be a square matrix",-1)
				call error_stop()
			endif
			E=eye(M,N,T%getType())
			inverse=linequ(T,E,RCOND)
			call inverse%resetdim(T%Dimension)
			return
		end if
		call inverse%allocate(T,T%getType())
		select case(T%getType())
			case(2)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(sp,[i,j])
							TMPT=sp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(sp,[i,j])
							sp=TMPinv
						end if
					end do
				end do
			case(3)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(dp,[i,j])
							TMPT=dp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(dp,[i,j])
							dp=TMPinv
						end if
					end do
				end do
			case(4)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(cp,[i,j])
							TMPT=cp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(cp,[i,j])
							cp=TMPinv
						end if
					end do
				end do
			case(5)
				do i=1,T%dim(1)
					do j=1,T%dim(2)
						if(T%GetFlag([i,j]))then
							call T%pointer(zp,[i,j])
							TMPT=zp
							M = TMPT%dim(1)
							N = TMPT%dim(2)
							if(M.ne.N) then
								call writemess("ERROR in calculating the inverse of a Tensor",-1)
								call writemess("input Tensor block should be a square matrix",-1)
								call error_stop()
							endif
							E=eye(M,N,T%getType())
							TMPinv=linequ(TMPT,E,RCOND)
							call TMPinv%resetdim(TMPT%Dimension)
							call inverse%pointer(zp,[i,j])
							zp=TMPinv
						end if
					end do
				end do
			case default
				call writemess('ERROR data type in inverse',-1)
				call error_stop
		end select
		return
	end function	
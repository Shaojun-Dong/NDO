program aaa
	use Density_type
	use targetFunction_Tools
	use Tensor_Tools
	use Parameter_Tools
	use Tools
	implicit none
	type(Tensor)::temp,AllTau
	integer::parameterLen,Nstep,LBFGSLen,delta_step_for_save,seed,lineSearchNum
	type(Tensor),allocatable::DataValue(:),Umatrix(:)
	real*8::E,stop_error,first_step,max_step,elementScal,tau
	character(len=200)::data_address,search_method,read_address
	type(List)::Parameter(3)
	type(OptimEngine)::OE
	type(OptData),target::RBMpoint
	type(OptimTensor),target::MLEPoint
	integer::lenBasic,i,j,k
	logical::replace
	class(OptimElement_structure),pointer::point
	call set_error_pointer
	call initial_log('output/runninglog')
	call open_file(111,'input/input.dat','old')
	call Parameter(1)%read(111)
	call Parameter(2)%read(111)
	call Parameter(3)%read(111)
	close(111)
	call writemess('Input Parameters are:')
	call writemess(Parameter(1))
	call writemess(Parameter(2))
	call writemess(Parameter(3))
	seed=Parameter(1)%ii('random_Seed')
	if(seed.ne.0)then
		call writemess('Set the random seed as seed='+seed)
		call set_seed(seed)
	end if
	data_address=Parameter(3)%ai('input_address')
	read_address=Parameter(3)%ai('read_address')
	call readData('input/'+data_address,DataValue,Umatrix)
	call initialOptim(OE,DataValue,Umatrix,Parameter)

	if(Parameter(3)%ai('model').equ.'RBM')then
		if(read_address.nequ.'null')then
			call writemess('Read the point from save/'+read_address)
			call open_file(112,'save/'+read_address,'old')
			call RBMpoint%read(112)
			close(112)
		else
			call RBMpoint%initialData(OE%ParameterLen,OE%phy_length,Parameter(2)%di('elementScal'))
			call inputcheck(RBMpoint)
		end if
		point=>RBMpoint
	else if(Parameter(3)%ai('model').equ.'MLE')then
		if(read_address.nequ.'null')then
			call writemess('Read the point from save/'+read_address)
			call open_file(112,'save/'+read_address,'old')
			call MLEPoint%Data%read(112)
			close(112)
		else
			call initialOptimTensor(MLEPoint,OE%phy_length,Parameter(2)%di('elementScal'))
		end if
		point=>MLEPoint
	else
		call writemess('ERROR model')
		call error_stop
	end if
	

	AllTau=Parameter(3)%Ti('AllTau')
	call temp%setType('real*8')
	Nstep=parameter(1)%ii('Nstep')
	lineSearchNum=Parameter(1)%ii('lineSearchNum')
	do k=1,AllTau%getTotalData()
		
		
		if(lineSearchNum.gt.0)then
			call OE%set_first_step_in_Line_search(AllTau%di(k))
			call writemess('first_step_in_Line_search='+AllTau%di(k))
			call reset_time_calculator(Nstep,30)
			call writemess('Running for i='+k+',Total i='+AllTau%getTotalData()+', lineSearchNum='+lineSearchNum)
			call OE%Optim(lineSearchNum,Nstep,point,E)
		else
			tau=AllTau%di(k)
			call reset_time_calculator(Nstep,30)
			call writemess('Running for i='+k+',Total i='+AllTau%getTotalData()+', tau='+tau)
			call OE%Optim(tau,Nstep,point,E)
		end if

		call writemess('divergences='+E)

		if(Parameter(3)%ai('model').equ.'RBM')then
			call open_file(115,'output/point.dat','replace')
			call RBMpoint%write(115)
			close(115)
			call open_file(115,'save/point.dat','replace')
			call RBMpoint%write(115)
			close(115)
		else if(Parameter(3)%ai('model').equ.'MLE')then
			call open_file(115,'output/point.dat','replace')
			call MLEPoint%data%write(115)
			close(115)
			call open_file(115,'save/point.dat','replace')
			call MLEPoint%data%write(115)
			close(115)
		else
			call writemess('ERROR model')
			call error_stop
		end if

		call open_file(116,'output/density.dat','replace')
		call OE%density(1)%write(116)
		close(116)


		temp=real(OE%density(1))
		call open_file(116,'output/density_real.dat','replace')
		call temp%write(116)
		close(116)

		temp=imag(OE%density(1))
		call open_file(116,'output/density_imag.dat','replace')
		call temp%write(116)
		close(116)

		do i=1,size(OE%density)
			call open_file(117,'output/probabity'+i+'.dat','replace')
			lenBasic=OE%density(i)%dim(1)
			do j=1,lenBasic
				write(117,*)j,OE%density(i)%di([j,j]),DataValue(i)%di(j)
			end do
			close(117)
		end do

		call open_file(117,'output/error.dat','replace')
		do i=1,size(OE%density)
			lenBasic=OE%density(i)%dim(1)
			do j=1,lenBasic
				write(117,*)abs(OE%density(i)%di([j,j])-DataValue(i)%di(j))
			end do
		end do
		close(117)

	end do
	


	stop

contains
	subroutine readData(address,probabity,BasicMatrix)
		character(len=*),intent(in)::address
		type(Tensor),intent(inout),allocatable::probabity(:),BasicMatrix(:)
		integer::TotalData,lenBasic,i,j,k
		real*8,allocatable::tmp1(:,:),tmp2(:,:)
		real*8,pointer::dp(:)
		complex*16,pointer::zp(:,:)
		call open_file(1234,address,'old')
		read(1234,*)lenBasic
		read(1234,*)TotalData
		if(allocated(probabity))then
			deallocate(probabity)
		end if
		if(allocated(BasicMatrix))then
			deallocate(BasicMatrix)
		end if
		allocate(probabity(TotalData))
		allocate(BasicMatrix(TotalData))
		allocate(tmp1(lenBasic,lenBasic))
		allocate(tmp2(lenBasic,lenBasic))
		do i=1,TotalData
			call probabity(i)%allocate([lenBasic],'real*8')
			call BasicMatrix(i)%allocate([lenBasic,lenBasic],'complex*16')
		end do
		do i=1,TotalData
			call probabity(i)%pointer(dp)
			read(1234,*)(dp(j),j=1,lenBasic)
			do j=1,lenBasic
				if(dp(j).lt.0)then
					dp(j)=-1d0*dp(j)
					if(dp(j).gt.1d-16)then
						call writemess('ERROR in input probabity')
						call error_stop
					end if
				end if
			end do
			do j=1,lenBasic
				read(1234,*)(tmp1(j,k),k=1,lenBasic)
			end do
			do j=1,lenBasic
				read(1234,*)(tmp2(j,k),k=1,lenBasic)
			end do
			call BasicMatrix(i)%pointer(zp)
			zp=dcmplx(tmp1,tmp2)
		end do
		return
	end subroutine

	
end program

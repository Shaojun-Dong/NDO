	interface
		subroutine copyData_interface(outA,inA)
			class(*),target,intent(in)::inA
			class(*),target,intent(inout)::outA
		end subroutine copyData_interface
	end interface
	procedure(copyData_interface),public,pointer::copy_Class_Data=>default_copy_Class_Data

	interface
		subroutine num_plus_num_interface(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		end subroutine num_plus_num_interface
	end interface
	procedure(num_plus_num_interface),public,pointer::num_plus_num=>default_num_plus_num

	interface
		subroutine num_minus_num_interface(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		end subroutine num_minus_num_interface
	end interface
	procedure(num_minus_num_interface),public,pointer::num_minus_num=>default_num_minus_num

	interface
		subroutine Num_Time_num_interface(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		end subroutine Num_Time_num_interface
	end interface
	procedure(Num_Time_num_interface),public,pointer::Num_Time_num=>default_Num_Time_num

	interface
		subroutine num_divide_num_interface(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		end subroutine num_divide_num_interface
	end interface
	procedure(num_divide_num_interface),public,pointer::num_divide_num=>default_num_divide_num

	interface
		subroutine write_class_data_interface(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		end subroutine write_class_data_interface
	end interface
	procedure(write_class_data_interface),public,pointer::write_class_data=>default_write_class_data

	interface
		subroutine read_class_data_interface(A,length,uni)
		class(*),target::A(:)
		integer::uni,length
		end subroutine read_class_data_interface
	end interface
	procedure(read_class_data_interface),public,pointer::read_class_data=>default_read_class_data

	interface
		subroutine ClassMPISend_interface(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		end subroutine ClassMPISend_interface
	end interface
	procedure(ClassMPISend_interface),public,pointer::mpi_send_class=>default_mpi_send_class

	interface
		subroutine ClassMpiRecv_interface(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		end subroutine ClassMpiRecv_interface
	end interface
	procedure(ClassMpiRecv_interface),public,pointer::mpi_recv_class=>default_mpi_recv_class

	interface
		subroutine ClassMpiBcast_interface(mess,length,ID,mpi_comm,ierr)
		class(*)::mess(:)
		integer::length,ID,MPI_Comm
		integer::ierr
		end subroutine ClassMpiBcast_interface
	end interface
	procedure(ClassMpiBcast_interface),public,pointer::mpi_bcast_class=>default_mpi_bcast_class

	interface
		subroutine ClassMpiAllreduce_interface(inmess,outmess,length,op,mpi_comm,ierr)
		class(*)::inmess(:),outmess(:)
		integer::length,op,mpi_comm,ierr
		end subroutine ClassMpiAllreduce_interface
	end interface
	procedure(ClassMpiAllreduce_interface),public,pointer::mpi_Allreduce_class=>default_mpi_Allreduce_class

!************************************

	interface
		subroutine Array_Divide_Arrayinterface(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		end subroutine Array_Divide_Arrayinterface
	end interface
	procedure(Array_Divide_Arrayinterface),public,pointer::Array_Divide_Array=>default_Array_Divide_Array
	
	interface
		subroutine Array_Time_Arrayinterface(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		end subroutine Array_Time_Arrayinterface
	end interface
	procedure(Array_Time_Arrayinterface),public,pointer::Array_Time_Array=>default_Array_Time_Array

	interface
		subroutine Set_Class_Datainterface(outA,inA,length)
			class(*),target,intent(inout)::outA(:)
			class(*),target,intent(in)::inA
			integer,intent(in)::length
		end subroutine Set_Class_Datainterface
	end interface
	procedure(Set_Class_Datainterface),public,pointer::Set_Class_Data=>default_Set_Class_Data
	
	interface
		subroutine Class_divide_numinterface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A(:)
			class(*),target,intent(in)::B
			integer,intent(in)::length
		end subroutine Class_divide_numinterface
	end interface
	procedure(Class_divide_numinterface),public,pointer::Class_divide_num=>default_Class_divide_num

	interface
		subroutine num_time_Class_interface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B(:)
			integer,intent(in)::length
		end subroutine num_time_Class_interface
	end interface
	procedure(num_time_Class_interface),public,pointer::num_time_Class=>default_num_time_Class

	
	interface
		subroutine Class_Time_numinterface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A(:)
			class(*),target,intent(in)::B
			integer,intent(in)::length
		end subroutine Class_Time_numinterface
	end interface
	procedure(Class_Time_numinterface),public,pointer::Class_Time_num=>default_Class_Time_num
	
	interface
		subroutine num_Minus_Class_interface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B(:)
			integer,intent(in)::length
		end subroutine num_Minus_Class_interface
	end interface
	procedure(num_Minus_Class_interface),public,pointer::num_Minus_Class=>default_num_Minus_Class

	
	interface
		subroutine Class_Minus_numinterface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A(:)
			class(*),target,intent(in)::B
			integer,intent(in)::length
		end subroutine Class_Minus_numinterface
	end interface
	procedure(Class_Minus_numinterface),public,pointer::Class_Minus_num=>default_Class_Minus_num

	interface
		subroutine num_Plus_Class_interface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B(:)
			integer,intent(in)::length
		end subroutine num_Plus_Class_interface
	end interface
	procedure(num_Plus_Class_interface),public,pointer::num_Plus_Class=>default_num_Plus_Class
	
	interface
		subroutine Class_Plus_num_interface(R,A,B,length)
			class(*),target,intent(inout)::R(:)
			class(*),target,intent(in)::A(:)
			class(*),target,intent(in)::B
			integer,intent(in)::length
		end subroutine Class_Plus_num_interface
	end interface
	procedure(Class_Plus_num_interface),public,pointer::Class_Plus_num=>default_Class_Plus_num

	interface
		subroutine Class_Plus_Class_interface(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		end subroutine Class_Plus_Class_interface
	end interface
	procedure(Class_Plus_Class_interface),public,pointer::Class_Plus_Class=>default_ClassPlusClass

	interface
		subroutine Class_Minus_Class_interface(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		end subroutine Class_Minus_Class_interface
	end interface
	procedure(Class_Minus_Class_interface),public,pointer::Class_Minus_Class=>default_ClassMinusClass


	interface
		subroutine copyArray_interface(outA,inA,length)
			integer,intent(in)::length
			class(*),target,intent(in)::inA(:)
			class(*),target,intent(inout)::outA(:)
		end subroutine copyArray_interface
	end interface
	procedure(copyArray_interface),public,pointer::copy_Class_Array=>default_copy_Class_Array

	interface
		subroutine copyARRAYSameType2D_interface(outA,inA,dim1,dim2)
		integer,intent(in)::dim1,dim2
		class(*),target,intent(in)::inA(:,:)
		class(*),target,intent(inout)::outA(:,:)
		end subroutine copyARRAYSameType2D_interface
	end interface
	procedure(copyARRAYSameType2D_interface),public,pointer::copyARRAY_Same_Type_2D=>default_copyARRAY_Same_Type_2D

	interface
		subroutine copyARRAYSameType3D_interface(outA,inA,dim1,dim2,dim3)
		integer,intent(in)::dim1,dim2,dim3
		class(*),target,intent(in)::inA(:,:,:)
		class(*),target,intent(inout)::outA(:,:,:)
		end subroutine copyARRAYSameType3D_interface
	end interface
	procedure(copyARRAYSameType3D_interface),public,pointer::copyARRAY_Same_Type_3D=>default_copyARRAY_Same_Type_3D

	interface
		subroutine Class_Writemess0_interface(mess,cpu_number)
		class(*),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		end subroutine Class_Writemess0_interface
	end interface
	procedure(Class_Writemess0_interface),public,pointer::Class_Writemess0=>default_Class_Writemess0

	interface
		subroutine Class_Writemess_from_interface(mess,form_,cpu_number)
		class(*),intent(in)::mess
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		end subroutine Class_Writemess_from_interface
	end interface
	procedure(Class_Writemess_from_interface),public,pointer::Class_Writemess_from=>default_Class_Writemess_from

	interface
		subroutine Class_Writemess_array_form_interface(mess,form_,cpu_number)
		class(*),intent(in)::mess(:)
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		end subroutine Class_Writemess_array_form_interface
	end interface
	procedure(Class_Writemess_array_form_interface),public,pointer::Class_Writemess_array_form=>default_Class_Writemess_array_form

	interface
		subroutine Class_Writemess_array_interface(mess,cpu_number)
		class(*),intent(in)::mess(:)
		integer,optional,intent(in)::cpu_number
		end subroutine Class_Writemess_array_interface
	end interface
	procedure(Class_Writemess_array_interface),public,pointer::Class_Writemess_array=>default_Class_Writemess_array
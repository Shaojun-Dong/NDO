	subroutine ClassMPISend(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		if(Length.eq.0)return
		select type(mess)
			type is (integer)
				call mpi_send(mess,length,MPI_integer,ID,tag,MPI_Comm,ierr)
			type is (real(kind=4))
				call mpi_send(mess,length,MPI_real,ID,tag,MPI_Comm,ierr)
			type is (real(kind=8))
				call mpi_send(mess,length,MPI_double_precision,ID,tag,MPI_Comm,ierr)
			type is (complex(kind=4))
				call mpi_send(mess,length,MPI_complex,ID,tag,MPI_Comm,ierr)
			type is (complex(kind=8))
				call mpi_send(mess,length,MPI_double_complex,ID,tag,MPI_Comm,ierr)
			type is (logical)
				call mpi_send(mess,length,MPI_logical,ID,tag,MPI_Comm,ierr)
			type is (character(len=*))
				call mpi_send(mess,len(mess)*length,MPI_character,ID,tag,MPI_Comm,ierr)
			class default
				call mpi_send_class(mess,Length,ID,tag,MPI_Comm,ierr)
		end select
		return
	end subroutine
	subroutine ClassMpiRecv(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		integer::istatus(MPI_STATUS_SIZE)
		if(Length.eq.0)return
		select type(mess)
			type is (integer)
				call mpi_recv(mess,length,MPI_integer,ID,tag,MPI_Comm,istatus,ierr)
			type is (real(kind=4))
				call mpi_recv(mess,length,MPI_real,ID,tag,MPI_Comm,istatus,ierr)
			type is (real(kind=8))
				call mpi_recv(mess,length,MPI_double_precision,ID,tag,MPI_Comm,istatus,ierr)
			type is (complex(kind=4))
				call mpi_recv(mess,length,MPI_complex,ID,tag,MPI_Comm,istatus,ierr)
			type is (complex(kind=8))
				call mpi_recv(mess,length,MPI_double_complex,ID,tag,MPI_Comm,istatus,ierr)
			type is (logical)
				call mpi_recv(mess,length,MPI_logical,ID,tag,MPI_Comm,istatus,ierr)
			type is (character(len=*))
				call mpi_recv(mess,len(mess)*length,MPI_character,ID,tag,MPI_Comm,istatus,ierr)
			class default
				call mpi_recv_class(mess,Length,ID,tag,MPI_Comm,ierr)
		end select
		return
	end subroutine

	subroutine ClassMpiBcast(mess,length,ID,mpi_comm,ierr)
		class(*)::mess(:)
		integer::length,ID,MPI_Comm
		integer::ierr
		if(Length.eq.0)return
		select type(mess)
			type is (integer)
				call MPI_BCAST(mess,length,MPI_integer,ID,mpi_comm,ierr)
			type is (real(kind=4))
				call MPI_BCAST(mess,length,MPI_real,ID,mpi_comm,ierr)
			type is (real(kind=8))
				call MPI_BCAST(mess,length,MPI_double_precision,ID,mpi_comm,ierr)
			type is (complex(kind=4))
				call MPI_BCAST(mess,length,MPI_complex,ID,mpi_comm,ierr)
			type is (complex(kind=8))
				call MPI_BCAST(mess,length,MPI_double_complex,ID,mpi_comm,ierr)
			type is (logical)
				call MPI_BCAST(mess,length,MPI_logical,ID,mpi_comm,ierr)
			type is (character(len=*))
				call MPI_BCAST(mess,len(mess)*length,MPI_character,ID,mpi_comm,ierr)
			class default
				call mpi_bcast_class(mess,length,ID,mpi_comm,ierr)
		end select
		return
	end subroutine

	subroutine ClassMpiAllreduce(inmess,outmess,length,op,mpi_comm,ierr)
		class(*)::inmess(:),outmess(:)
		integer::length,op,mpi_comm,ierr
		if(Length.eq.0)return
		select type(inmess)
			type is (integer)
				select type(outmess)
					type is (integer)
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_integer,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (real(kind=4))
				select type(outmess)
					type is (real(kind=4))
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_real,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (real(kind=8))
				select type(outmess)
					type is (real(kind=8))
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_double_precision,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (complex(kind=4))
				select type(outmess)
					type is (complex(kind=4))
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_complex,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (complex(kind=8))
				select type(outmess)
					type is (complex(kind=8))
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_double_complex,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (logical)
				select type(outmess)
					type is (logical)
						call MPI_ALLREDUCE(inmess,outmess,length,MPI_logical,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			type is (character(len=*))
				select type(outmess)
					type is (character(len=*))
						if(len(inmess).ne.len(outmess))then
							call writemess('ERROR in MPI_ALLREDUCE',-1)
							call error_stop
						end if
						call MPI_ALLREDUCE(inmess,outmess,len(inmess)*length,MPI_character,OP,mpi_comm,ierr)
					class default
						call writemess('ERROR in mpi_allreduce',-1)
						call error_stop
				end select
			class default
				call mpi_Allreduce_class(inmess,outmess,length,op,mpi_comm,ierr)
		end select
		return
	end subroutine

	subroutine default_mpi_send_class(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		call writemess('ERROR in mpi_send, error class type',-1)
		call error_stop
	end subroutine 
	subroutine default_mpi_recv_class(mess,Length,ID,tag,MPI_Comm,ierr)
		class(*)::mess(:)
		integer::length,ID,tag,MPI_Comm
		integer::ierr
		call writemess('ERROR in mpi_recv, error class type',-1)
		call error_stop
	end subroutine 
	subroutine default_mpi_bcast_class(mess,length,ID,mpi_comm,ierr)
		class(*)::mess(:)
		integer::length,ID,MPI_Comm
		integer::ierr
		call writemess('ERROR in mpi_bcast, error class type',-1)
		call error_stop
	end subroutine 
	subroutine default_mpi_Allreduce_class(inmess,outmess,length,op,mpi_comm,ierr)
		class(*)::inmess(:),outmess(:)
		integer::length,op,mpi_comm,ierr
		call writemess('ERROR in mpi_Allreduce, error class type',-1)
		call error_stop
	end subroutine 

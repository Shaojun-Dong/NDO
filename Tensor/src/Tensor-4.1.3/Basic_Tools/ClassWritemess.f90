	subroutine ClassWritemess0(mess,cpu_number)
		class(*),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		select type(mess)
			type is (integer)
				call writemess(mess,cpu_number)
			type is (real(kind=4))
				call writemess(mess,cpu_number)
			type is (real(kind=8))
				call writemess(mess,cpu_number)
			type is (complex(kind=4))
				call writemess(mess,cpu_number)
			type is (complex(kind=8))
				call writemess(mess,cpu_number)
			type is (logical)
				call writemess(mess,cpu_number)
			type is (character(len=*))
				call writemess(mess,cpu_number)
			class default
				call Class_Writemess0(mess,cpu_number)
		end select
		return
	end subroutine
	subroutine ClassWritemess_form(mess,form_,cpu_number)
		class(*),intent(in)::mess
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		select type(mess)
			type is (integer)
				call writemess(mess,cpu_number)
			type is (real(kind=4))
				call writemess(mess,cpu_number)
			type is (real(kind=8))
				call writemess(mess,cpu_number)
			type is (complex(kind=4))
				call writemess(mess,cpu_number)
			type is (complex(kind=8))
				call writemess(mess,cpu_number)
			type is (logical)
				call writemess(mess,cpu_number)
			type is (character(len=*))
				call writemess(mess,cpu_number)
			class default
				call Class_Writemess_from(mess,form_,cpu_number)
		end select
		return
	end subroutine
	subroutine ClassWritemess_array_form(mess,form_,cpu_number)
		class(*),intent(in)::mess(:)
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		select type(mess)
			type is (integer)
				call writemess(mess,cpu_number)
			type is (real(kind=4))
				call writemess(mess,cpu_number)
			type is (real(kind=8))
				call writemess(mess,cpu_number)
			type is (complex(kind=4))
				call writemess(mess,cpu_number)
			type is (complex(kind=8))
				call writemess(mess,cpu_number)
			type is (logical)
				call writemess(mess,cpu_number)
			type is (character(len=*))
				call writemess(mess,cpu_number)
			class default
				call Class_Writemess_array_form(mess,form_,cpu_number)
		end select
		return
	end subroutine
	subroutine ClassWritemess_array(mess,cpu_number)
		class(*),intent(in)::mess(:)
		integer,optional,intent(in)::cpu_number
		select type(mess)
			type is (integer)
				call writemess(mess,cpu_number)
			type is (real(kind=4))
				call writemess(mess,cpu_number)
			type is (real(kind=8))
				call writemess(mess,cpu_number)
			type is (complex(kind=4))
				call writemess(mess,cpu_number)
			type is (complex(kind=8))
				call writemess(mess,cpu_number)
			type is (logical)
				call writemess(mess,cpu_number)
			type is (character(len=*))
				call writemess(mess,cpu_number)
			class default
				call Class_Writemess_array(mess,cpu_number)
		end select		
		return
	end subroutine

	subroutine default_Class_Writemess0(mess,cpu_number)
		class(*),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		character(len=characterlen)::TMPa
		call copyData(TMPa,mess)
		call writemess(TMPa,cpu_number)
	end subroutine 

	subroutine default_Class_Writemess_from(mess,form_,cpu_number)
		class(*),intent(in)::mess
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		character(len=characterlen)::TMPa
		call copyData(TMPa,mess)
		call writemess(TMPa,form_,cpu_number)
	end subroutine 

	subroutine default_Class_Writemess_array_form(mess,form_,cpu_number)
		class(*),intent(in)::mess(:)
		character(len=*),intent(in)::form_
		integer,optional,intent(in)::cpu_number
		character(len=characterlen),allocatable::TMPa(:)
		integer::leng_mess
		leng_mess=ClassSize(mess)
		allocate(TMPa(leng_mess))
		call FastcopyARRAY(TMPa,mess,leng_mess)
		call writemess(TMPa,form_,cpu_number)
	end subroutine 

	subroutine default_Class_Writemess_array(mess,cpu_number)
		class(*),intent(in)::mess(:)
		integer,optional,intent(in)::cpu_number
		character(len=characterlen),allocatable::TMPa(:)
		integer::leng_mess
		leng_mess=ClassSize(mess)
		allocate(TMPa(leng_mess))
		call FastcopyARRAY(TMPa,mess,leng_mess)
		call writemess(TMPa,cpu_number)
	end subroutine 
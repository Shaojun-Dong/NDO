#define T_gt_DATANAME T_gt_i
#define T_ge_DATANAME T_ge_i
#define T_lt_DATANAME T_lt_i
#define T_le_DATANAME T_le_i
#define T_eq_DATANAME T_eq_i
#define T_equ_DATANAME T_equ_i
#define T_ne_DATANAME T_ne_i
#define T_nequ_DATANAME T_nequ_i
#define DATANAME_gt_T i_gt_T
#define DATANAME_ge_T i_ge_T
#define DATANAME_lt_T i_lt_T
#define DATANAME_le_T i_le_T
#define DATANAME_eq_T i_eq_T
#define DATANAME_equ_T i_equ_T
#define DATANAME_ne_T i_ne_T
#define DATANAME_nequ_T i_nequ_T
#define DATATYPE integer
#include "templet/logicalFunc0.f90"
#undef T_gt_DATANAME
#undef T_ge_DATANAME
#undef T_lt_DATANAME
#undef T_le_DATANAME
#undef T_eq_DATANAME
#undef T_equ_DATANAME
#undef T_ne_DATANAME 
#undef T_nequ_DATANAME 
#undef DATANAME_gt_T
#undef DATANAME_ge_T
#undef DATANAME_lt_T
#undef DATANAME_le_T
#undef DATANAME_eq_T
#undef DATANAME_equ_T
#undef DATANAME_ne_T
#undef DATANAME_nequ_T
#undef DATATYPE

#define T_gt_DATANAME T_gt_s
#define T_ge_DATANAME T_ge_s
#define T_lt_DATANAME T_lt_s
#define T_le_DATANAME T_le_s
#define T_eq_DATANAME T_eq_s
#define T_equ_DATANAME T_equ_s
#define T_ne_DATANAME T_ne_s
#define T_nequ_DATANAME T_nequ_s
#define DATANAME_gt_T s_gt_T
#define DATANAME_ge_T s_ge_T
#define DATANAME_lt_T s_lt_T
#define DATANAME_le_T s_le_T
#define DATANAME_eq_T s_eq_T
#define DATANAME_equ_T s_equ_T
#define DATANAME_ne_T s_ne_T
#define DATANAME_nequ_T s_nequ_T
#define DATATYPE real*4
#include "templet/logicalFunc0.f90"
#undef T_gt_DATANAME
#undef T_ge_DATANAME
#undef T_lt_DATANAME
#undef T_le_DATANAME
#undef T_eq_DATANAME
#undef T_equ_DATANAME
#undef T_ne_DATANAME 
#undef T_nequ_DATANAME 
#undef DATANAME_gt_T
#undef DATANAME_ge_T
#undef DATANAME_lt_T
#undef DATANAME_le_T
#undef DATANAME_eq_T
#undef DATANAME_equ_T
#undef DATANAME_ne_T
#undef DATANAME_nequ_T
#undef DATATYPE

#define T_gt_DATANAME T_gt_d
#define T_ge_DATANAME T_ge_d
#define T_lt_DATANAME T_lt_d
#define T_le_DATANAME T_le_d
#define T_eq_DATANAME T_eq_d
#define T_equ_DATANAME T_equ_d
#define T_ne_DATANAME T_ne_d
#define T_nequ_DATANAME T_nequ_d
#define DATANAME_gt_T d_gt_T
#define DATANAME_ge_T d_ge_T
#define DATANAME_lt_T d_lt_T
#define DATANAME_le_T d_le_T
#define DATANAME_eq_T d_eq_T
#define DATANAME_equ_T d_equ_T
#define DATANAME_ne_T d_ne_T
#define DATANAME_nequ_T d_nequ_T
#define DATATYPE real*8
#include "templet/logicalFunc0.f90"
#undef T_gt_DATANAME
#undef T_ge_DATANAME
#undef T_lt_DATANAME
#undef T_le_DATANAME
#undef T_eq_DATANAME
#undef T_equ_DATANAME
#undef T_ne_DATANAME 
#undef T_nequ_DATANAME 
#undef DATANAME_gt_T
#undef DATANAME_ge_T
#undef DATANAME_lt_T
#undef DATANAME_le_T
#undef DATANAME_eq_T
#undef DATANAME_equ_T
#undef DATANAME_ne_T
#undef DATANAME_nequ_T
#undef DATATYPE



	logical function gt_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.gt.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.2)then
			call writemess("ONLY FOR Tensor with one element,(.gt.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
			call error_stop()
		end if
		gt_of_Tensor=T1%di(1).gt.T2%di(1)
		return
	end function
	logical function le_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.le.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.le.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
			call error_stop()
		end if
		le_of_Tensor=T1%di(1).le.T2%di(1)
		return
	end function
	logical function lt_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.lt.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.lt.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
			call error_stop()
		end if
		lt_of_Tensor=T1%di(1).lt.T2%di(1)
		return
	end function
	logical function ge_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.ge.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.ge.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
			call error_stop()
		end if
		ge_of_Tensor=T1%di(1).ge.T2%di(1)
		return
	end function
	logical function equal_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Tensor) :: T
		integer :: l,flag
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.eq.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.eq.)!!",-1)
			call error_stop()
		end if
		select case(T1%getType())
			case(1)
				equal_of_Tensor=T1%ii(1).eq.T2%ii(1)
			case(2)
				equal_of_Tensor=T1%si(1).eq.T2%si(1)
			case(3)
				equal_of_Tensor=T1%di(1).eq.T2%di(1)
			case(4)
				equal_of_Tensor=T1%ci(1).eq.T2%ci(1)
			case(5)
				equal_of_Tensor=T1%zi(1).eq.T2%zi(1)
			case(6)
				equal_of_Tensor=T1%li(1).eqv.T2%li(1)
			case(7)
				equal_of_Tensor=T1%ai(1).equ.T2%ai(1)
		end select
		return
	end function

	logical function equ_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Tensor) :: T
		integer :: l,flag
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.equ.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.equ.)',-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.equ.)!!",-1)
			call error_stop()
		end if
		if(T1%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.equ.)!!",-1)
			call error_stop()
		end if
		select case(T1%getType())
			case(1)
				equ_of_Tensor=T1%ii(1).eq.T2%ii(1)
			case(2)
				equ_of_Tensor=T1%si(1).equ.T2%si(1)
			case(3)
				equ_of_Tensor=T1%di(1).equ.T2%di(1)
			case(4)
				equ_of_Tensor=T1%ci(1).equ.T2%ci(1)
			case(5)
				equ_of_Tensor=T1%zi(1).equ.T2%zi(1)
			case(6)
				equ_of_Tensor=T1%li(1).eqv.T2%li(1)
			case(7)
				equ_of_Tensor=T1%ai(1).equ.T2%ai(1)
		end select
		return
	end function

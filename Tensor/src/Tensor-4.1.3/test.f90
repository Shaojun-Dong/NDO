program main
	use Basic_Tools
	use Pointer_Tools
	use memory_type
	implicit none
	integer,target::idata(3)
	real*8,target::sdata(3),rdata(3)
	class(*),pointer::clp1(:),clp2(:),clp3(:)
	clp1=>idata
	clp2=>sdata
	clp3=>rdata
	call FastcopyARRAY(clp1,[1.,2.,3.],3)
	call FastcopyARRAY(clp2,[1,1,1],3)
	call ClassMinusClass(clp3,clp1,clp2,3)
	write(*,*)rdata
	call ClasstimeNum(clp3,clp1,0.5,3)
	write(*,*)rdata,sdata*[1,2,3]
contains


end

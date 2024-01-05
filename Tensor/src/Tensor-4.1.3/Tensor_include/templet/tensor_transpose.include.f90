subroutine TAT_FUNCNAME(idata, odata, dim_i, plan, rank)
   TAT_TYPENAME, intent(in) :: idata(:)
   TAT_TYPENAME, intent(out) :: odata(:)
   integer, intent(in) :: dim_i(:), plan(:), rank
   integer, allocatable :: ldi(:), ldo(:), dim_o(:), ldi_i(:)
   integer :: i, ii, io, active
   integer, allocatable :: index(:)

   allocate(ldi(rank))
   allocate(ldo(rank))
   allocate(dim_o(rank))
   allocate(ldi_i(rank))
   allocate(index(rank))

   do i=1, rank
      dim_o(i) = dim_i(plan(i))
   end do
   ldo(1) = 1
   do i=1, rank- 1
      ldo(i+1) = ldo(i)*dim_o(i)
      end do
   ldi_i(1) = 1
   do i=1, rank - 1
      ldi_i(i+1) = ldi_i(i)*dim_i(i)
   end do
   do i=1, rank
      ldi(i) = ldi_i(plan(i))
   end do

   do i=1, rank
      index(i) = 1
   end do
   ii = 1
   io = 1
   do
      odata(io) = idata(ii)
      active = 1
      index(active) = index(active) + 1
      ii = ii + ldi(active)
      io = io + ldo(active)
      do while (index(active) > dim_o(active))
         index(active) = 1
         ii = ii - dim_o(active) * ldi(active)
         io = io - dim_o(active) * ldo(active)
         if (active == rank) then
            return
         end if
         active = active + 1
         index(active) = index(active) + 1
         ii = ii + ldi(active)
         io = io + ldo(active)
      end do
   end do
end subroutine

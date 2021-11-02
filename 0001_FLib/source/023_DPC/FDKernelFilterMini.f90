module FDKernelFilterMini
	!****************************************************************************
	!  Kernel filter module
	!****************************************************************************
	!============================================================================
	!　Includes
	!============================================================================
	implicit none
	
	!============================================================================
	!  Module Procedures
	!============================================================================
	contains
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  apply kernel filter
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine apply_filter(hpc, vpc, pc, ph, pv, pw, org, edge, image, mask)
		!**************************************************************************
		!  apply kernel filter to a image
		!**************************************************************************
		!==========================================================================
		!  Includes
		!==========================================================================
		use FImgInterpolationMini
		
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  image pixel count
		integer, intent(in) :: hpc, vpc
		!  2.  kernel point count
		integer, intent(in) :: pc
		!  3.  kernel point coordinates & weight
		real(8), intent(in) :: ph(pc), pv(pc), pw(pc)
		!  4.  original image data
		real(8), intent(in) :: org(hpc, vpc)
		!  5.  excluding edge region
		real(8), intent(in) :: edge
		!  6.  Filtered image data
		real(8), intent(out) :: image(hpc, vpc)
		!  7.  Valid region data
		integer, intent(out) :: mask(hpc, vpc)
		
		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.  loop variable
		integer :: i, j, k
		!  2.  interpolation buffer
		real(8) :: val(1)
		logical :: is_in
		
		!==========================================================================
		!　Operations
		!==========================================================================
		do i = 1, vpc
			do j = 1, hpc
				mask(j, i) = 1
				image(j, i) = 0.0d0
				do k = 1, pc
					call bicubic_intp(dble(j) + ph(k), dble(i) + pv(k), &
						& 1, hpc, vpc, org, edge, is_in, val)
					if(is_in) then
						image(j, i) = image(j, i) + pw(k) * val(1)
					else
						mask(j, i) = 0
						image(j, i) = 0.0d0
						exit
					end if
				end do				
			end do
		end do
		
	end subroutine apply_filter

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  kernel filter of various filter type
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine antisymm_filter_data(pave, have, vave, wave, rp, rq, pc, ph, pv, pw)
		!**************************************************************************
		!  calculate antisymmetric filter data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  Averaging point count
		integer, intent(in) :: pave
		!  2.  Averaging point coordinate and weight
		real(8), intent(in) :: have(pave), vave(pave), wave(pave)
		!  3.  Lattice vector
		real(8), intent(in) :: rp(2), rq(2)
		!  4.  filter point count
		integer, intent(out) :: pc
		!  5.  filter point coordinate and weight
		real(8), allocatable, intent(out) :: ph(:), pv(:), pw(:)
		
		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.  loop variable
		integer :: i
		!  2.  kernel filter point count
		integer, parameter :: kfpc = 5
		
		!==========================================================================
		!　Operations
		!==========================================================================
		pc = pave * kfpc
		allocate(ph(pc))
		allocate(pv(pc))
		allocate(pw(pc))
		do i = 1, pave
			ph((i - 1) * kfpc + 1) = have(i)
			pv((i - 1) * kfpc + 1) = vave(i)
			pw((i - 1) * kfpc + 1) = wave(i) * 0.5d0
			
			ph((i - 1) * kfpc + 2) = have(i) + rp(1)
			pv((i - 1) * kfpc + 2) = vave(i) + rp(2)
			pw((i - 1) * kfpc + 2) = wave(i) * (- 0.125d0)
			
			ph((i - 1) * kfpc + 3) = have(i) - rp(1)
			pv((i - 1) * kfpc + 3) = vave(i) - rp(2)
			pw((i - 1) * kfpc + 3) = wave(i) * (- 0.125d0)
			
			ph((i - 1) * kfpc + 4) = have(i) + rq(1)
			pv((i - 1) * kfpc + 4) = vave(i) + rq(2)
			pw((i - 1) * kfpc + 4) = wave(i) * (- 0.125d0)
			
			ph((i - 1) * kfpc + 5) = have(i) - rq(1)
			pv((i - 1) * kfpc + 5) = vave(i) - rq(2)
			pw((i - 1) * kfpc + 5) = wave(i) * (- 0.125d0)			
		end do
		
	end subroutine antisymm_filter_data

end module FDKernelFilterMini
	
!******************************************************************************
!  Library Interface 
!******************************************************************************
subroutine fd_apply_kernel_filter_lif(filter_type, hpc, vpc, rp, rq, &
	& org, edge, pave, have, vave, wave, image, mask) &
	& bind(c, name = "FDApplyKernelFilterMini")
	!dec$ attributes dllexport :: fd_apply_kernel_filter_lif
	!****************************************************************************
	!  Remove electric potential component
	!****************************************************************************
	!============================================================================
	!　Includes
	!============================================================================
	use FDKernelFilterMini
	implicit none
	
	!============================================================================
	!  Arguments
	!============================================================================
	!  1.  Filter type
	integer, value :: filter_type
	!  2.  Pixel count of the image
	integer, value :: hpc, vpc
	!  3.  Lattice vector
	real(8) :: rp(2), rq(2)
	!  4.  Original image data
	real(8) :: org(hpc, vpc)
	
	real(8), value :: edge
	!
	integer, value :: pave
	!
	real(8) :: have(pave), vave(pave), wave(pave)
	!  5.  Filter applied image data
	real(8) :: image(hpc, vpc)
	!  6.  Valid region mask
	!      1 -> valid pixel
	!      0 -> invalid pixel
	integer :: mask(hpc, vpc)
	
	!============================================================================
	!  Locals
	!============================================================================
	integer :: pc
	real(8), allocatable :: ph(:), pv(:), pw(:)
	
	!============================================================================
	!　Operations
	!============================================================================
	call antisymm_filter_data(pave, have, vave, wave, rp, rq, pc, ph, pv, pw)
	call apply_filter(hpc, vpc, pc, ph, pv, pw, org, edge, image, mask)
	
end subroutine fd_apply_kernel_filter_lif	
	


	
	
	
	
	
	
	
	
	
	
	
	
	
	
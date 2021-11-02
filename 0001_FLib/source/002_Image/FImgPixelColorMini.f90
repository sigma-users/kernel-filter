module FImgPixelColorMini
	!****************************************************************************
	!  Plot color module
	!****************************************************************************
	!============================================================================
	!  Includes
	!============================================================================
	implicit none
	
	!============================================================================
	!  Module Procedures
	!============================================================================
	contains
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  Gray scale image
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	function gray_u8(val, plot_vmin, plot_vmax) result(ret)
		!**************************************************************************
		!  Convert double to 1 byte gray scale data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  input data
		real(8), intent(in) :: val
		!  2.  plot range
		real(8), intent(in) :: plot_vmin, plot_vmax
		!  3.  output
		integer :: ret
		
		!==========================================================================
		!  Locals
		!==========================================================================
		real(8) :: buf
		
		!==========================================================================
		!  Operations
		!==========================================================================
		if(plot_vmax <= plot_vmin) then
			ret = 255
			return
		end if
		
		buf = 256.0d0 * (val - plot_vmin) / (plot_vmax - plot_vmin)
		if(buf > 255.0d0) then
			ret = 255
		else if(buf < 0.0d0) then
			ret = 0
		else
			ret = int(buf)
		end if
	end function gray_u8
	
	function gray_u32(val, plot_vmin, plot_vmax) result(ret)
		!************************************************************************
		!  Convert double to 4 byte RGB gray scale data
		!************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  input data
		real(8), intent(in) :: val
		!  2.  plot range
		real(8), intent(in) :: plot_vmin, plot_vmax
		!  3.  output
		integer :: ret
		
		!==========================================================================
		!  Locals
		!==========================================================================
		integer :: buf
		
		!==========================================================================
		!  Operations
		!==========================================================================
		buf = gray_u8(val, plot_vmin, plot_vmax)
		
		
		ret = buf
		ret = ishft(ret, 8)
		ret = or(ret, buf)
		ret = ishft(ret, 8)
		ret = or(ret, buf)
		
	end function gray_u32
		
	subroutine calc_scalar_plot_data(color_index, width, height, image, &
		& plot_vmin, plot_vmax, gray)
		!**************************************************************************
		!  Convert scalar image data to 4 byte RGB gray scale data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  Color tone index
		integer, intent(in) :: color_index
		!  2.  image pixel size
		integer, intent(in) :: width, height
		!  3.  image intensity
		real(8), intent(in) :: image(width, height)
		!  4.  plot range
		real(8), intent(in) :: plot_vmin, plot_vmax
		!  5.  gray scale data
		integer, intent(out) :: gray(width, height)

		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.  loop index
		integer :: i, j
		
		!==========================================================================
		!　Operations
		!==========================================================================
		!!$omp parallel do private(i, j) 
		do i = 1, height
			do j = 1, width
				gray(j, i) = gray_u32(image(j, i), plot_vmin, plot_vmax)
			end do
		end do
		!!$omp end parallel do
		
	end subroutine calc_scalar_plot_data
		
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  Color image
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	function rgb_vs_u32(vx, vy, vs, plot_vsmin, plot_vsmax) result(ret)
		!**************************************************************************
		!  Convert vector data with absolute value to 4 byte RGB color data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  vector data(vx, vy) and absolute value(vs)
		real(8), intent(in) :: vx, vy, vs
		!  2.  plot range
		real(8), intent(in) :: plot_vsmin, plot_vsmax
		!  3.  output
		integer(4) :: ret
		
		!==========================================================================
		!  Locals
		!==========================================================================
		real(8), parameter :: kc_pi = 3.14159265358979323846264338328d0
		real(8) :: buf
		integer :: li
		real(8) :: r, g, b, h, s, v, f
		
		!==========================================================================
		!  Operations
		!==========================================================================
		!  hue
		h = atan2(vy, vx) * 3.0d0 / kc_pi
		if(h < 0.0d0) then
			h = h + 6.0d0
		end if
		
		!  saturation
		s = 1.0d0
		
		!  value
		if(plot_vsmax <= plot_vsmin) then
		  v = 1.0d0
		else
			v = (vs - plot_vsmin) / (plot_vsmax - plot_vsmin)
			if(v > 1.0d0) then
				v = 1.0d0
			else if(v < 0.0d0) then
				v = 0.0d0
			end if
		end if
		
		!  convert hsv to rgb
		li = int(h)
		f = h - dble(li)
			
		r = v
		g = v
		b = v
		if(li == 0) then
			g = g * (1.0d0 - s * (1.0d0 - f))
			b = b * (1.0d0 - s)
		else if(li == 1) then
			r = r * (1.0d0 - s * f)
			b = b * (1.0d0 - s)
		else if(li == 2) then
			r = r * (1.0d0 - s)
			b = b * (1.0d0 - s * (1.0d0 - f))
		else if(li == 3) then
			r = r * (1.0d0 - s)
			g = g * (1.0d0 - s * f)
		else if(li == 4) then
			r = r * (1.0d0 - s * (1.0d0 - f))
			g = g * (1.0d0 - s)
		else
			g = g * (1.0d0 - s)
			b = b * (1.0d0 - s * f)
		end if
		
		ret = int(255.0d0 * r + 0.5d0)
		ret = ishft(ret, 8)
		ret = or(ret, int(255.0d0 * g + 0.5d0))
		ret = ishft(ret, 8)
		ret = or(ret, int(255.0d0 * b + 0.5d0))
			
	end function rgb_vs_u32
	
	function rgb_u32(vx, vy, plot_vsmin, plot_vsmax) result(ret)
		!**************************************************************************
		!  Convert vector data to 4 byte RGB color data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  vector data
		real(8), intent(in) :: vx, vy
		!  2.  plot range
		real(8), intent(in) :: plot_vsmin, plot_vsmax
		!  3.  output
		integer(4) :: ret
		
		!==========================================================================
		!  Operations
		!==========================================================================
		ret = rgb_vs_u32(vx, vy, dsqrt(vx * vx + vy * vy), plot_vsmin, plot_vsmax)
		
	end function rgb_u32
	
	subroutine calc_vector_plot_data(width, height, &
		& vx, vy, plot_vsmin, plot_vsmax, rgb)
		!**************************************************************************
		!  convert vector image data to 4 byte RGB color data
		!**************************************************************************
		!==========================================================================
		!  Arguments
		!==========================================================================
		!  1.  image pixel size
		integer, intent(in) :: width, height
		!  2.  vector data
		real(8), intent(in) :: vx(width, height), vy(width, height)
		!  3.  plot range
		real(8), intent(in) :: plot_vsmin, plot_vsmax
		!  4.  RGB color data
		integer, intent(out) :: rgb(width, height)

		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.  loop index
		integer :: i, j
		
		!==========================================================================
		!　Operations
		!==========================================================================
		!$omp parallel do private(i, j) 
		do i = 1, height
			do j = 1, width
				rgb(j, i) = rgb_u32(vx(j, i), vy(j, i), plot_vsmin, plot_vsmax)
			end do
		end do
		!$omp end parallel do
		
	end subroutine calc_vector_plot_data
		
end module FImgPixelColorMini
	
!==============================================================================
!　Library interface
!==============================================================================
subroutine fimg_scalar_plot_data_lif(color_index, width, height, image, smin, smax, gray) &
	& bind(c, name = "FImgScalarPlotData")
	!dec$ attributes dllexport :: fimg_scalar_plot_data_lif
	!****************************************************************************
	!　Convert scalar image data to 4 byte RGB gray scale data
	!****************************************************************************
	!============================================================================
	!　Includes
	!============================================================================
	use FImgPixelColorMini
	implicit none
	
	!============================================================================
	!  Arguments
	!============================================================================
	!  1.  Color tone index
	integer, value :: color_index
	!  2.  image pixel size
	integer, value :: width, height
	!  3.  image intensity
	real(8) :: image(width, height)
	!  4.  plot range
	real(8), value :: smin, smax
	!  5.  gray scale data
	integer :: gray(width, height)

	!============================================================================
	!　Operations
	!============================================================================
	call calc_scalar_plot_data(color_index, width, height, image, smin, smax, gray)
		
end subroutine fimg_scalar_plot_data_lif		
	
subroutine fimg_vector_plot_data_lif(width, height, fxa, fya, fmin, fmax, rgb) &
	& bind(c, name = "FImgVectorPlotData")
	!dec$ attributes dllexport :: fimg_vector_plot_data_lif
	!****************************************************************************
	!　 convert vector image data to 4 byte RGB color data
	!****************************************************************************
	!============================================================================
	!　Includes
	!============================================================================
	use FImgPixelColorMini
	implicit none
	
	!============================================================================
	!  Arguments
	!============================================================================
	!  1.  image pixel size
	integer, value :: width, height
	!  2.  vector data
	real(8) :: fxa(width, height), fya(width, height)
	!  3.  plot range
	real(8), value :: fmin, fmax
	!  4.  RGB color data
	integer :: rgb(width, height)

	!============================================================================
	!　Operations
	!============================================================================
	call calc_vector_plot_data(width, height, fxa, fya, fmin, fmax, rgb)
		
end subroutine fimg_vector_plot_data_lif		
	


	
module FImgInterpolationMini
	!============================================================================
	!  Includes
	!============================================================================
	implicit none
	
	!============================================================================
	!  Module Procedure
	!============================================================================
	contains
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  Bicubic interpolation
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	function bicubic_check_point(h, v, hpc, vpc, edge) result(is_in)
		!**************************************************************************
		!  Check if the point resides in valid region.
		!**************************************************************************
		!==========================================================================
		!　Arguments
		!==========================================================================
		!  1.  Poisition in the image
		real(8), intent(in) :: h, v
		!  2.　Pixel count of the image
		integer, intent(in) :: hpc, vpc
		!  3.  Exclusion area at the edge
		real(8), intent(in) :: edge		
		!  4.  flag for the point resides in the valid region
		logical :: is_in
		
		!==========================================================================
		!　Operations
		!==========================================================================
		if(h < 1.0d0 + edge .or. h > dble(hpc) - edge &
			& .or. v < 1.0d0 + edge .or. v > dble(vpc) - edge) then
			is_in = .false.
		else
			is_in = .true.
		end if
		
	end function bicubic_check_point
		
	subroutine bicubic_weight(h, v, hpc, vpc, hs, vs, weight)
		!**************************************************************************
		!  Calculating interpolation weight nearby the point
		!**************************************************************************
		!==========================================================================
		!　Arguments
		!==========================================================================
		!  1.  Poisition in the image
		real(8), intent(in) :: h, v
		!  2.　Pixel count of the image
		integer, intent(in) :: hpc, vpc
		!  3.  Position of the weight matrix
		integer, intent(out) :: hs, vs
		!  4.  The interpolation weight value
		real(8), intent(out) :: weight(4, 4)
		
		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.	 loop index
		integer :: i
		!  2.  weight
		real(8) :: wh1, wh2, wh3, wh4, wv1, wv2, wv3, wv4
		real(8) :: w11, w12, w13, w14, w21, w22, w23, w24
		real(8) :: w31, w32, w33, w34, w41, w42, w43, w44
		!  3.  pixel position
		integer :: hi, vi, h1, h2, h3, h4, v1, v2, v3, v4
		
		!==========================================================================
		!　Operations
		!==========================================================================
		hi = int(h)
		vi = int(v)
		
		if(hi == 1) then
			hs = 1
			h1 = 0
			h2 = 1
			h3 = 2
			h4 = 3
			wh1 = weight_1d(h - h1)
			wh2 = weight_1d(h - h2)
			wh3 = weight_1d(h3 - h)
			wh4 = weight_1d(h4 - h)
				
			if(vi == 1) then
				vs = 1
				v1 = 0
				v2 = 1
				v3 = 2
				v4 = 3
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w22 =	4.0d0 * wh1 * wv1 + 2.0d0 * wh1 * wv2 + 2.0d0 * wh2 * wv1 + wh2 * wv2
				w23 = -2.0d0 * wh1 * wv1 + 2.0d0 * wh1 * wv3 - wh2 * wv1 + wh2 * wv3
				w24 =	2.0d0 * wh1 * wv4 + wh2 * wv4
				w32 = -2.0d0 * wh1 * wv1 - wh1 * wv2 + 2.0d0 * wh3 * wv1 + wh3 * wv2
				w33 =	wh1 * wv1 - wh1 * wv3 - wh3 * wv1 + wh3 * wv3
				w34 = - wh1 * wv4 + wh3 * wv4
				w42 =	2.0d0 * wh4 * wv1 + wh4 * wv2
				w43 =	- wh4 * wv1 + wh4 * wv3
				w44 =	wh4 * wv4				

				weight(1, 1) = w22
				weight(1, 2) = w23
				weight(1, 3) = w24
				weight(1, 4) = 0.0d0
				weight(2, 1) = w32
				weight(2, 2) = w33
				weight(2, 3) = w34
				weight(2, 4) = 0.0d0
				weight(3, 1) = w42
				weight(3, 2) = w43
				weight(3, 3) = w44
				weight(3, 4) = 0.0d0
				weight(4, 1) = 0.0d0
				weight(4, 2) = 0.0d0
				weight(4, 3) = 0.0d0
				weight(4, 4) = 0.0d0
				
			else if(vi >= vpc - 1) then
				vs = vpc - 3
				v1 = vpc - 2
				v2 = vpc - 1
				v3 = vpc
				v4 = vpc + 1				
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w21 = 2.0d0 * wh1 * wv1 + wh2 * wv1
				w22 = 2.0d0 * wh1 * wv2 - 2.0d0 * wh1 * wv4 + wh2 * wv2 - wh2 * wv4
				w23 = 2.0d0 * wh1 * wv3 + 4.0d0 * wh1 * wv4 + wh2 * wv3 + 2.0d0 * wh2 * wv4
				w31 = - wh1 * wv1 + wh3 * wv1
				w32 = - wh1 * wv2 + wh1 * wv4 + wh3 * wv2 - wh3 * wv4
				w33 = - wh1 * wv3 - 2.0d0 * wh1 * wv4 + wh3 * wv3 + 2.0d0 * wh3 * wv4
				w41 = wh4 * wv1
				w42 = wh4 * wv2 - wh4 * wv4
				w43 = wh4 * wv3 + 2.0d0 * wh4 * wv4
				
				weight(1, 1) = 0.0d0
				weight(1, 2) = w21
				weight(1, 3) = w22
				weight(1, 4) = w23
				weight(2, 1) = 0.0d0
				weight(2, 2) = w31
				weight(2, 3) = w32
				weight(2, 4) = w33
				weight(3, 1) = 0.0d0
				weight(3, 2) = w41
				weight(3, 3) = w42
				weight(3, 4) = w43
				weight(4, 1) = 0.0d0
				weight(4, 2) = 0.0d0
				weight(4, 3) = 0.0d0
				weight(4, 4) = 0.0d0
				
			else
				vs = vi - 1
				v1 = vi - 1
				v2 = vi
				v3 = vi + 1
				v4 = vi + 2
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w21 = 2.0d0 * wh1 * wv1 + wh2 * wv1
				w22 = 2.0d0 * wh1 * wv2 + wh2 * wv2
				w23 = 2.0d0 * wh1 * wv3 + wh2 * wv3
				w24 = 2.0d0 * wh1 * wv4 + wh2 * wv4
				w31 = - wh1 * wv1 + wh3 * wv1
				w32 = - wh1 * wv2 + wh3 * wv2
				w33 = - wh1 * wv3 + wh3 * wv3
				w34 = - wh1 * wv4 + wh3 * wv4
				w41 = wh4 * wv1
				w42 = wh4 * wv2
				w43 = wh4 * wh3
				w44 = wh4 * wv4
				
				weight(1, 1) = w21
				weight(1, 2) = w22
				weight(1, 3) = w23
				weight(1, 4) = w24
				weight(2, 1) = w31
				weight(2, 2) = w32
				weight(2, 3) = w33
				weight(2, 4) = w34
				weight(3, 1) = w41
				weight(3, 2) = w42
				weight(3, 3) = w43
				weight(3, 4) = w44
				weight(4, 1) = 0.0d0
				weight(4, 2) = 0.0d0
				weight(4, 3) = 0.0d0
				weight(4, 4) = 0.0d0
				
			end if
		else if(hi >= hpc - 1) then
			hs = hpc - 3
			h1 = hpc - 2
			h2 = hpc - 1
			h3 = hpc
			h4 = hpc + 1
			wh1 = weight_1d(h - h1)
			wh2 = weight_1d(h - h2)
			wh3 = weight_1d(h3 - h)
			wh4 = weight_1d(h4 - h)
			
			if(vi == 1) then				
				vs = 1
				v1 = 0
				v2 = 1
				v3 = 2
				v4 = 3				
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w12 = 2.0d0 * wh1 * wv1 + wh1 * wv2
				w13 = - wh1 * wv1 + wh1 * wv3
				w14 = wh1 * wv4
				w22 = 2.0d0 * wh2 * wv1 + wh2 * wv2 - 2.0d0 * wh4 * wv1 - wh4 * wv2
				w23 = - wh2 * wv1 + wh2 * wv3 + wh4 * wv1 - wh4 * wv3
				w24 = wh2 * wv4 - wh4 * wv4
				w32 = 2.0d0 * wh3 * wv1 + wh3 * wv2 + 4.0d0 * wh4 * wv1 + 2.0d0 * wh4 * wv2
				w33 = - wh3 * wv1 + wh3 * wv3 - 2.0d0 * wh4 * wv1 + 2.0d0 * wh4 * wv3
				w34 = wh3 * wv4 + 2.0d0 * wh4 * wv4
				
				weight(1, 1) = 0.0d0
				weight(1, 2) = 0.0d0
				weight(1, 3) = 0.0d0
				weight(1, 4) = 0.0d0
				weight(2, 1) = w12
				weight(2, 2) = w13
				weight(2, 3) = w14
				weight(2, 4) = 0.0d0
				weight(3, 1) = w22
				weight(3, 2) = w23
				weight(3, 3) = w24
				weight(3, 4) = 0.0d0
				weight(4, 1) = w32
				weight(4, 2) = w33
				weight(4, 3) = w34
				weight(4, 4) = 0.0d0
				
			else if(vi >= vpc - 1) then
				vs = vpc - 3
				v1 = vpc - 2
				v2 = vpc - 1
				v3 = vpc
				v4 = vpc + 1
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w11 = wh1 * wv1
				w12 = wh1 * wv2 - wh1 * wv4
				w13 = wh1 * wv3 + 2.0d0 * wh1 * wv4
				w21 = wh2 * wv1 - wh4 * wv1
				w22 = wh2 * wv2 - wh2 * wv4 - wh4 * wv2 + wh4 * wv4
				w23 = wh2 * wv3 + 2.0d0 * wh2 * wv4 - wh4 * wv3 - 2.0d0 * wh4 * wv4
				w31 = wh3 * wv1 + 2.0d0 * wh4 * wv1
				w32 = wh3 * wv2 - wh3 * wv4 + 2.0d0 * wh4 * wv2 - 2.0d0 * wh4 * wv4
				w33 = wh3 * wv3 + 2.0d0 * wh3 * wv4 + 2.0d0 * wh4 * wv3 + 4.0d0 * wh4 * wv4
				
				weight(1, 1) = 0.0d0
				weight(1, 2) = 0.0d0
				weight(1, 3) = 0.0d0
				weight(1, 4) = 0.0d0
				weight(2, 1) = 0.0d0
				weight(2, 2) = w11
				weight(2, 3) = w12
				weight(2, 4) = w13
				weight(3, 1) = 0.0d0
				weight(3, 2) = w21
				weight(3, 3) = w22
				weight(3, 4) = w23
				weight(4, 1) = 0.0d0
				weight(4, 2) = w31
				weight(4, 3) = w32
				weight(4, 4) = w33
				
			else
				vs = vi - 1
				v1 = vi - 1
				v2 = vi
				v3 = vi + 1
				v4 = vi + 2
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w11 = wh1 * wv1
				w12 = wh1 * wv2
				w13 = wh1 * wv3
				w14 = wh1 * wv4
				w21 = wh2 * wv1 - wh4 * wv1
				w22 = wh2 * wv2 - wh4 * wv2
				w23 = wh2 * wv3 - wh4 * wv3
				w24 = wh2 * wv4 - wh4 * wv4
				w31 = wh3 * wv1 + 2.0d0 * wh4 * wv1
				w32 = wh3 * wv2 + 2.0d0 * wh4 * wv2
				w33 = wh3 * wv3 + 2.0d0 * wh4 * wv3
				w34 = wh3 * wv4 + 2.0d0 * wh4 * wv4
								
				weight(1, 1) = 0.0d0
				weight(1, 2) = 0.0d0
				weight(1, 3) = 0.0d0
				weight(1, 4) = 0.0d0
				weight(2, 1) = w11
				weight(2, 2) = w12
				weight(2, 3) = w13
				weight(2, 4) = w14
				weight(3, 1) = w21
				weight(3, 2) = w22
				weight(3, 3) = w23
				weight(3, 4) = w24
				weight(4, 1) = w31
				weight(4, 2) = w32
				weight(4, 3) = w33
				weight(4, 4) = w34
				
			end if
			
		else
			hs = hi - 1			
			h1 = hi - 1
			h2 = hi
			h3 = hi + 1
			h4 = hi + 2
			wh1 = weight_1d(h - h1)
			wh2 = weight_1d(h - h2)
			wh3 = weight_1d(h3 - h)
			wh4 = weight_1d(h4 - h)
			
			if(vi == 1) then
				vs = 1
				v1 = 0
				v2 = 1
				v3 = 2
				v4 = 3
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w12 = 2.0d0 * wh1 * wv1 + wh1 * wv2
				w13 = - wh1 * wv1 + wh1 * wv3
				w14 = wh1 * wv4
				w22 = 2.0d0 * wh2 * wv1 + wh2 * wv2
				w23 = - wh2 * wv1 + wh2 * wv3
				w24 = wh2 * wv4
				w32 = 2.0d0 * wh3 * wv1 + wh3 * wv2
				w33 = - wh3 * wv1 + wh3 * wv3
				w34 = wh3 * wv4
				w42 = 2.0d0 * wh4 * wv1 + wh4 * wv2
				w43 = - wh4 * wv1 + wh4 * wv3
				w44 = wh4 * wv4
								
				weight(1, 1) = w12
				weight(1, 2) = w13
				weight(1, 3) = w14
				weight(1, 4) = 0.0d0
				weight(2, 1) = w22
				weight(2, 2) = w23
				weight(2, 3) = w24
				weight(2, 4) = 0.0d0
				weight(3, 1) = w32
				weight(3, 2) = w33
				weight(3, 3) = w34
				weight(3, 4) = 0.0d0
				weight(4, 1) = w42
				weight(4, 2) = w43
				weight(4, 3) = w44
				weight(4, 4) = 0.0d0
				
			else if(vi >= vpc - 1) then
				vs = vpc - 3
				v1 = vpc - 2
				v2 = vpc - 1
				v3 = vpc
				v4 = vpc + 1				
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)
				
				w11 = wh1 * wv1
				w12 = wh1 * wv2 - wh1 * wv4
				w13 = wh1 * wv3 + 2.0d0 * wh1 * wv4
				w21 = wh2 * wv1
				w22 = wh2 * wv2 - wh2 * wv4
				w23 = wh2 * wv3 + 2.0d0 * wh2 * wv4
				w31 = wh3 * wv1
				w32 = wh3 * wv2 - wh3 * wv4
				w33 = wh3 * wv3 + 2.0d0 * wh3 * wv4
				w41 = wh4 * wv1
				w42 = wh4 * wv2 - wh4 * wv4
				w43 = wh4 * wv3 + 2.0d0 * wh4 * wv4				
								
				weight(1, 1) = 0.0d0
				weight(1, 2) = w11
				weight(1, 3) = w12
				weight(1, 4) = w13
				weight(2, 1) = 0.0d0
				weight(2, 2) = w21
				weight(2, 3) = w22
				weight(2, 4) = w23
				weight(3, 1) = 0.0d0
				weight(3, 2) = w31
				weight(3, 3) = w32
				weight(3, 4) = w33
				weight(4, 1) = 0.0d0
				weight(4, 2) = w41
				weight(4, 3) = w42
				weight(4, 4) = w43
				
			else
				vs = vi - 1
				v1 = vi - 1
				v2 = vi
				v3 = vi + 1
				v4 = vi + 2				
				wv1 = weight_1d(v - v1)
				wv2 = weight_1d(v - v2)
				wv3 = weight_1d(v3 - v)
				wv4 = weight_1d(v4 - v)

				w11 = wh1 * wv1
				w12 = wh1 * wv2
				w13 = wh1 * wv3
				w14 = wh1 * wv4
				w21 = wh2 * wv1
				w22 = wh2 * wv2
				w23 = wh2 * wv3
				w24 = wh2 * wv4
				w31 = wh3 * wv1
				w32 = wh3 * wv2
				w33 = wh3 * wv3
				w34 = wh3 * wv4
				w41 = wh4 * wv1
				w42 = wh4 * wv2
				w43 = wh4 * wv3
				w44 = wh4 * wv4
								
				weight(1, 1) = w11
				weight(1, 2) = w12
				weight(1, 3) = w13
				weight(1, 4) = w14
				weight(2, 1) = w21
				weight(2, 2) = w22
				weight(2, 3) = w23
				weight(2, 4) = w24
				weight(3, 1) = w31
				weight(3, 2) = w32
				weight(3, 3) = w33
				weight(3, 4) = w34
				weight(4, 1) = w41
				weight(4, 2) = w42
				weight(4, 3) = w43
				weight(4, 4) = w44

			end if
		end if
		
		!==========================================================================
		!  Inner functions
		!==========================================================================
		contains
		function weight_1d(d) result(w)
			!************************************************************************
			!  calculate the weight for the bicubic interpolation
			!************************************************************************
			!========================================================================
			!  Arguments
			!========================================================================
			!  1.  length to the nearby point
			real(8) :: d
			!  2.  weight
			real(8) :: w
			
			!========================================================================
			!  Locals
			!========================================================================
			!  1.  sharpness of bicubic interpolation.
			!      -0.5 and -1.0 are often used.
			real(8), parameter :: a = - 1.0d0

			!========================================================================
			!  Operations
			!========================================================================
			d = dabs(d)
			if(d < 1) then
				w = (a + 2.0d0) * d * d * d - (a + 3.0d0) * d * d + 1.0d0
			else if(dabs(d) < 2) then
				w = a * d * d * d - 5.0d0 * a * d * d + 8.0d0 * a * d - 4.0d0 * a			
			else
				w = 0.0d0
			end if		
		end function weight_1d
		
	end subroutine bicubic_weight
		
	subroutine bicubic_intp(h, v, count, hpc, vpc, image, edge, is_in, val)
		!**************************************************************************
		!  Calculate bicubic interpolation value of the image
		!**************************************************************************
		!==========================================================================
		!　Arguments
		!==========================================================================
		!  1.  Poisition in the image
		real(8), intent(in) :: h, v
		!  2.  Image count
		integer, intent(in) :: count
		!  3.　Pixel count of the image
		integer, intent(in) :: hpc, vpc
		!  4.  Image data
		real(8), intent(in) :: image(hpc, vpc, count)
		!  5.  Exclusion area at the edge
		real(8), intent(in) :: edge
		!  6.  flag for the point resides in the valid region
		logical, intent(out) :: is_in
		!  7.  The interpolated value
		real(8), intent(out) :: val(count)
		
		!==========================================================================
		!　Locals
		!==========================================================================
		!  1.  loop index
		integer :: i
		!  2.  weight matrix
		real(8) :: cw(4, 4)
		!  3.  pixel position
		integer :: hs, vs
		
		!==========================================================================
		!　Operations
		!==========================================================================
		is_in = bicubic_check_point(h, v, hpc, vpc, edge)
		if(.not. is_in) then
			val = 0.0d0
			return
		end if
		
		call bicubic_weight(h, v, hpc, vpc, hs, vs, cw)
		
		do i = 1, count
			val(i) = &
				&   cw(1, 1) * image(hs, vs, i) &
				& + cw(1, 2) * image(hs, vs + 1, i) &
				& + cw(1, 3) * image(hs, vs + 2, i) &
				& + cw(1, 4) * image(hs, vs + 3, i) &
				& + cw(2, 1) * image(hs + 1, vs, i) &
				& + cw(2, 2) * image(hs + 1, vs + 1, i) &
				& + cw(2, 3) * image(hs + 1, vs + 2, i) &
				& + cw(2, 4) * image(hs + 1, vs + 3, i) &
				& + cw(3, 1) * image(hs + 2, vs, i) &
				& + cw(3, 2) * image(hs + 2, vs + 1, i) &
				& + cw(3, 3) * image(hs + 2, vs + 2, i) &
				& + cw(3, 4) * image(hs + 2, vs + 3, i) &
				& + cw(4, 1) * image(hs + 3, vs, i) &
				& + cw(4, 2) * image(hs + 3, vs + 1, i) &
				& + cw(4, 3) * image(hs + 3, vs + 2, i) &
				& + cw(4, 4) * image(hs + 3, vs + 3, i)
		end do
		
	end subroutine bicubic_intp
		
end module FImgInterpolationMini
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
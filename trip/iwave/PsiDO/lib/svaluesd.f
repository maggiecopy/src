       	subroutine svaluesd(x1i,hi,vi,ni,s,iend,work,wl,ierr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Second Derivative Estimator
c     by Ann Campbell
c
c     function:  This program generates estimates for the
c                values of the second derivative of a function
c                at given points.
c
c	the  parameters are:
c 	x1i:	initial x value for input data 
c	hi:	sample rate for input data
c	vi:	function values for input data
c	ni:	number of input data points
c	s:      output vector, contains second derivative estimates 
c	iend:	type of end condition to be used
c		iend=1	linear ends, s(1)=s(n)=0
c		iend=2	parabolic ends, s(1)=s(2),s(n)=s(n-1)
c		iend=3  cubic ends, s(1),s(n) are extrapolated
c			end point
c	work:	workspace vector
c	wl:	length of work vector
c	ierr:	error flag, will equal 0 if no errors
c       i:      do loop counter
c       first:  first element used in gaussian elimination
c       last:   last element used in gaussian elimination
c       j:      do loop counter, used to hold intermediate values
c       dv1,dv2:stores intermediate difference calculations
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       


	integer ni,iend,i,first,last,j,wl,ierr
	double precision x1i,hi,vi(ni),work(wl)
	double precision dv1,dv2,s(ni)


c
c	first test to see if error flag equals 0 (no error).
c       if it does not, it should return to the calling
c       program.
c

	if (ierr.eq.0) then
c
c	then test to see if the length of the workspace vector
c	is sufficient to calculate the spline.
c       if it is not, set an error code and return to the calling
c       program.
c
	if (wl.lt.4*ni-8) then
		ierr=10
c		write(6,*)'** Error in second derivative estimating program **'
c		write(6,*)'       not enough workspace'
c		write(6,*)
c		write(6,*)wl,' was length provided'
c		write(6,*)4*ni-8,' was length needed'
		return
	end if

	do 5 i=1,wl
	work(i)=0.0D0
 5	continue

	 
c
c     to determine the coefficients of the spline, we need
c     estimates of the second derivatives at the different
c     points of the input vector.
c     these estimates are determined by multiplying the inverse
c     of a tridiagonal matrix with the product of another matrix
c     with the input vector.
c     because we know the tridiagonal matrix, rather than its
c     inverse, and the vector produced by the other multiplication,
c     we can use gaussian elimination to solve for the value of
c     the second derivatives.

c     the components of the tridiagonal matrix are stored in
c     work(1..ni-2), work(ni-1..2*ni-4), and work(2*ni-3...3*ni-6).
c     the other vector is stored in work(3*ni-5...4*ni-8).
c
	dv1= (vi(2)-vi(1))/hi*6.0D0
	do 10 i=1,ni-2
	dv2= (vi(i+2)-vi(i+1))/hi*6.0D0
	work(i)= hi
	work(i+ni-2)= 4.0D0*hi
	work(i+2*ni-4)= hi
	work(i+3*ni-6)= dv2-dv1
	dv1=dv2
 10	continue
	first= 2
	last= ni-2

c	adjust for different end conditions
c	iend = 1:no changes
c	iend = 2:

	if (iend.eq.2) then
 50		work(ni-1)= work(ni-1) + hi
		work(2*ni-4)= work(2*ni-4)+hi
c	iend = 3:

	else if (iend.eq.3) then
		work(ni-1)= 6D0*hi
		work(2*ni-3)= 0.0D0
		work(ni-2)= 0.0D0
		work(2*ni-4)=6D0*hi
	end if
c	
c	this is where the gaussian elimination is performed
c

	do 110 i=first,last
	work(i)= work(i)/work(i+ni-3)
	work(i+ni-2)= work(i+ni-2)-work(i)*work(i+2*ni-5)
	work(i+3*ni-6)=work(i+3*ni-6)-work(i)*work(i+3*ni-7)
 110	continue

	work(last+3*ni-6)=work(last+3*ni-6)/work(last+ni-2)
	do 120 j=last-1,first-1,-1
	work(j+3*ni-6)=(work(j+3*ni-6)-work(j+2*ni-4)*work(j+3*ni-5))
     1	    /work(j+ni-2)
 120	continue

c
c     now that the second derivative values have been calculated,
c     they will be stored in work(2..ni-1)
c
	do 130 i=first-1,last
	work(i+1)= work(i+3*ni-6)
 130	continue

c
c	the second derivative values must be adjusted for 
c       different end conditions, making the second derivative
c	estimates now occupy work(1...ni).
c
	if (iend.eq.1) then
		work(1)= 0.0D0
		work(ni)= 0.0D0
	else if (iend.eq.2) then
		work(1)= work(2)
		work(ni)= work(ni-1)
	else if (iend.eq.3) then
		work(1)=2D0*work(2)-work(3)
		work(ni)=2D0*work(ni-1)-work(ni-2)
        end if

	do 140 i=1,ni
	s(i)=work(i)
 140	continue

	end if
	return
	end





















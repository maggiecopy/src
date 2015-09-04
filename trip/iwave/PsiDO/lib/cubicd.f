       	subroutine cubicd(x1i,hi,vi,ni,x1o,ho,vo,no,iend,work,wl,ierr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Cubic Spline
c     by Ann Campbell
c
c     function:  This program generates a cubic spline using an
c                input vector and uses it to create an output
c                vector with a specific initial value, sample rate,
c                and number of data points.
c
c	the  parameters are:
c 	x1i:	initial x value for input data 
c	hi:	sample rate for input data
c	vi:	function values for input data
c	ni:	number of input data points
c	x1o:	initial x value for output data
c	ho:	sample rate for output data
c	vo:	function values as calculated by spline
c	no:	number of output data points 
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
c       k:      used to hold intermediate values
c       int:    interval of the input values containing a specific
c               output value
c       dv1,dv2:stores intermediate difference calculations
c       dx:     the difference between the x value of a component
c               of the output vector and the x value of the input
c               vector that is closest to it
c       u:      current output x value
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       


	integer ni,no,iend,i,first,last,j,k,int,wl,ierr
	double precision x1i,hi,vi(ni),work(wl)
	double precision x1o,ho,vo(no)
	double precision dv1,dv2,dx,u


c
c	first test to see if error flag equals 0 (no error).
c       if it does not, it should return to the calling
c       program.
c

	if (ierr.eq.0) then

c
c	next, check if an appropriate end condition has been supplied
c
	if (iend.lt.1.or.iend.gt.3) then
		ierr=11
c		write(6,*)'** Error in cubic spline program **'
c		write(6,*)'   incorrect end condition given'
c		write(6,*)'  choose end condition between 1 and 3'
		return
	end if

c
c	now check if there are enough input data points for the
c	program to work correctly
c
	if (ni.lt.3) then
		ierr=12
c		write(6,*)'** Error in cubic spline program **'
c		write(6,*)'  not enough input data points'
		return
	end if
c
c	then test to see if the length of the workspace vector
c	is sufficient to calculate the spline.
c       if it is not, set an error code and return to the calling
c       program.
c
	if (wl.lt.4*ni-3) then
		ierr=10
c		write(6,*)'** Error in cubic spline program **'
c		write(6,*)'       not enough workspace'
c		write(6,*)
c		write(6,*)wl,' was length provided'
c		write(6,*)4*ni-3,' was length needed'
		return
	end if
	 
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
	dv1= (vi(2)-vi(1))/hi*6.0
	do 10 i=1,ni-2
	dv2= (vi(i+2)-vi(i+1))/hi*6.0
	work(i)= hi
	work(i+ni-2)= 4.0*hi
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
		work(ni-1)= 6*hi
		work(2*ni-3)= 0.0
		work(ni-2)= 0.0
		work(2*ni-4)=6*hi
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
		work(1)= 0.0
                work(ni)= 0.0
	else if (iend.eq.2) then
		work(1)= work(2)
		work(ni)= work(ni-1)
	else if (iend.eq.3) then
		work(1)=2*work(2)-work(3)
		work(ni)=2*work(ni-1)-work(ni-2)
	end if


c
c     the second derivative values are now used to calculate
c     the coefficients of the spline.
c     the first coefficient(a) is stored in work(ni+1..2*ni-1),
c     the second coefficient(b) is stored in work(2*ni...3*ni-2),
c     and the third coefficient(c) is stored in work(3*ni-1,4*ni-3).
c	
	do 200 i=1,ni-1
	work(i+ni)= (work(i+1)-work(i))/(6*hi)
	work(i+2*ni-1)= work(i)/2
	work(i+3*ni-2)=((vi(i+1)-vi(i))/hi)-((2*hi*work(i)+hi*
     1	      work(i+1))/6)
 200	continue


c
c     the spline is now used to produce the output vector.
c     
c     for each component of the output vector, the interval of
c     the input vector which contains this value must first 
c     be determined.
c     this interval determines which coefficients of the spline
c     are used, as well as which input data value is added.
c

	int= 1
	u=x1o 
	do 300 i=1,no
	if (u.ge.x1i+(ni-1)*hi) then
		dx= u-(x1i+(ni-2)*hi)
		vo(i)=vi(ni-1)+dx*(work(4*ni-3)+dx*(work(3*ni-2)
     1		      +dx*work(2*ni-1)))
		u=u+ho
		goto 300
	end if
	if (u.ge.(x1i+(int-1)*hi)) then
		if (u.le.x1i+int*hi) then
		dx= u-(x1i+(int-1)*hi)
		vo(i)=vi(int)+dx*(work(int+3*ni-2)+dx*(work(int+2*ni-1)+
     1		      dx*work(int+ni)))
		u=u+ho
		goto 300
		end if
	end if
	int=1
	j= ni+1
 250	k=(int+j)/2
	if (u.lt.x1i+(k-1)*hi) then
		j=k
	else
		int=k
	end if
	if (j.gt.int+1) then
		goto 250
	else
		dx= u- (x1i+(int-1)*hi)
		vo(i)=vi(int)+dx*(work(int+3*ni-2)+dx*(work(int+2*ni-1)+
     1		      dx*work(int+ni)))
		u=u+ho
		goto 300
	end if
 300	continue	
	end if
	return
	end


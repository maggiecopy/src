	subroutine bicubicadj(x1in,dxin,nxin,y1in,dyin,nyin,
     1      zin,x1out,dxout,nxout,y1out,dyout,nyout,zout,iend,
     1	    work,wl,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Bicubic Spline Adjoint
c	by Ann Campbell
c
c	function:  This program takes data in two-dimensions, and
c		   by using transposes to cubic splines generates data
c		   in two-dimensions with defined initial values,     
c		   sample rates, and number of data points
c
c	the parameters are:
c	x1in:	initial x value(first dimension) for input data
c	dxin:	sample rate for x values for input data
c	nxin:	number of x values for input data
c	y1in:	initial y value(second dimension) for input data
c	dyin:	sample rate for y values for input data
c	nyin:	number of y values for input data
c	zin:	function values for input data, a matrix of size
c		nxin by nyin
c	x1out:	initial x value(first dimension) for output data
c	dxout:	sample rate for x values for output data
c	nxout:	number of x values for output data
c	y1out:	initial y value(second dimension) for output data
c	dyout:	sample rate for y values for output data
c	nyout:	number of y values for output data
c	zout:	output vector, of size nxout by nyout
c	iend:	type of end condition to be used
c		iend=1 linear ends
c		iend=2	parabolic ends
c		iend=3	cubic ends
c	work:   workspace vector
c	wl:	length of workspace vector
c	ierr:	error code, will equal 0 if no errors
c	gin/gout/gwork: used to determine allotment of workspace
c	ptr1,ptr2,ptr3: used as pointers in workspace, for use in
c		1-D programs
c	wl2:	workspace allowed in 1-D program
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	integer nxin,nyin,nxout,nyout,iend
	real x1in,dxin,y1in,dyin,zin(nxin,nyin)  
	real x1out,dxout,y1out,dyout,zout(nxout,nyout)

	integer gin,gout,gwork,ptr1,ptr2,ptr3,wl2,wl
	real work(wl)

        integer ierr,i,j

c
c	first, check if error flag equals 0 (no error).
c	if it does not, return to the calling program.


	if (ierr.eq.0) then

c
c	check if an appropriate end condition has been given
c
	if (iend.lt.1.or.iend.gt.3) then
		ierr=11
c		write(6,*)'** Error in bicubic adjoint program **'
c		write(6,*)'   incorrect end condition given'
c		write(6,*)'choose end condition between 1 and 3'
		return
	end if

c
c	check if there are a sufficient number of input and output
c	data points for the program to work properly.
c
	if (nxin.lt.3.or.nyin.lt.3) then
		ierr=12
c		write(6,*)'** Error in bicubic adjoint program **'
c		write(6,*)'   not enough input data points'
		return
	end if

	if (nxout.lt.4.or.nyout.lt.4) then
		ierr=12
c		write(6,*)'** Error in bicubic adjoint program **'
c		write(6,*)'   not enough output data points'
		return
	end if

	if (nxin.ge.nyin) then
		gin=nxin
	else
		gin=nyin
	end if
	if (nxout.ge.nyout) then
		gout=nxout
	else
		gout=nyout
	end if 

c
c	now check if a sufficient amount of workspace has
c	been provided
c
	if (7*nxout-8+nxin.ge.7*nyout-8+nyin) then
		gwork=7*nxout-8+nxin
	else
		gwork=7*nyout-8+nyin
	end if

	if (wl.lt.gin+gout+gwork+nxout*nyin) then
		ierr=10
c		write(6,*)'**Error in bicubic spline adjoint program**'
c		write(6,*)'       not enough workspace'
c		write(6,*)wl,' was length provided'
c		write(6,*)gin+gout+gwork+nxout*nyin,' was length needed'
		return
	end if

	do 5 i=1,wl
	work(i)=0.0
 5	continue

c
c	each column of the input matrix(work(1..nxin)) is converted to a vector
c	(work(1+gin..gin+nxout) of length nxout by applying the transpose 
c	of the cubic spline.
c	these vectors are stored in work(gin+gout+1..gin+gout+nxout*nyin).
c
	do 100 i=1,nyin
	do 50 j=1,nxin
	work(j)=zin(j,i)
 50	continue
	ptr1=1
	ptr2=1+gin
	ptr3=gin+gout+nxout*nyin+1
	wl2=wl-ptr3+1
	call cubictrans(x1in,dxin,work(ptr1),nxin,x1out,
     1           dxout,work(ptr2),nxout,iend,work(ptr3),wl2,ierr)
	do 60 j=1,nxout
	work(j+gin+gout+(i-1)*nxout)=work(gin+j)
 60	continue
 100	continue


c
c	each 'row' of the 'matrix' stored in work(gin+gout+1..
c	gin+gout+nxout*nyin) is converted to a vector of
c	length nyout by applying the transpose of the cubic spline.
c	these vectors are stored in the output matrix.
c

	do 200 i=1,nxout
	do 150 j=1,nyin
	work(j)=work(i+gin+gout+(j-1)*nxout)
 150	continue
	call cubictrans(y1in,dyin,work(ptr1),nyin,y1out,
     1       dyout,work(ptr2),nyout,iend,work(ptr3),wl2,ierr)

c
c	the next multiplication finally makes the matrix the adjoint
c
	do 175 j=1,nyout
	zout(i,j)=(dxin*dyin/(dxout*dyout))*work(j+gin)
 175	continue
 200	continue
	end if
	return
	end

	

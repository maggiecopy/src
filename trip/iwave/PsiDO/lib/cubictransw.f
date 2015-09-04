	subroutine cubictrans(x1i,hi,vi,ni,x1o,
     1     ho,vo,no,iend,work,wl,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Cubic Spline Transpose
c     by Ann Campbell
c
c     function:  This program generates the transpose of the cubic
c                spline and uses it to create an output vector
c                with a specific initial value, sample rate, and
c                number of data points.
c
c
c	the parameters are:
c	x1i:	initial x value for input data
c	hi:	sample rate for input data
c	vi:	function values for input data
c	ni:	number of input data points
c	x1o:	initial x value for output data
c	ho:	sample rate for output data
c	vo:	output vector
c	no:	number of output data points
c	iend:	type of end condition to be used
c		=1	linear ends, s(1)=s(n)=0
c		=2	parabolic ends,s(1)=s(2),s(n)=s(n-1)
c		=3	cubic ends,s(1),s(n) are extrapolated
c	work:	vector of workspace
c	wl:	length of work vector
c	ierr:	error flag, should equal 0 if no error
c       i:      do loop counter
c       first:  first component used in gaussian
c               elimination
c       last:   last component used in gaussian
c               elimination
c       j,k:    used in intermediate calculations
c       interval:the interval of the output vector that contains
c               a specific component of the input vector
c       trial:  do loop counter for the gaussian elimination
c       u:      the current input x value
c       dx:     the difference between the x value of a component
c               of the input vector and the x value of the output
c               vector that is closest to it
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	integer ni,no,iend,i,first,last,j,k,interval,wl,ierr
	integer trial
	real x1i,hi,vi(ni),x1o,ho,vo(no),work(wl),u,dx
        real temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
c
c	first test if error flag equals 0 (no error).
c       if it does not, it should return to the calling
c       program.
c

	if (ierr.eq.0) then
c
c	then test to see if the length of the workspace vector
c	is sufficient to calculate the adjoint.
c       if it is not, set an error code and return to the calling
c       program.	
c	*note: for this program to work correctly, the input vector
c	must have length at least 4 and the output vector must have
c	length at least 3
c
	if (wl.lt.ni+7*no-8) then
		ierr=10   
c		write(6,*)'**Error in cubic transpose program**'
c		write(6,*)'      not enough workspace'
c		write(6,*)
c		write(6,*)wl,' was length provided'
c		write(6,*)ni+7*no-8,'was length needed'
		return
	end if 

c	
c 	now, set everything to zero
c
	do 1 i=1,wl
	work(i)=0.0
 1	continue

	first=0
	last=0
	interval=0
	trial=0
	u=0
	dx=0

c
c     to calculate the adjoint, the first step is to determine
c     what interval of the output values that each of the input
c     data points lies in.
c     the intervals are stored in work(1 . .ni)
c

	interval=1
	u=x1i
	do 10 i=1,ni
	if (u.ge.x1o+(no-1)*ho) then
		work(i)=no-1
		u=u+hi
		goto 10
	end if
	if (u.ge.x1o+(interval-1)*ho) then
	  if (u.le.x1o+interval*ho) then
		work(i)=interval
		u=u+hi
		goto 10
          end if
	end if
	interval=1
	j=no+1
 5	k=(interval+j)/2
	if (u.lt.x1o+(k-1)*ho) then
		j=k
	else 
		interval=k
	end if
	if (j.gt.interval+1) then
		goto 5
	else
		work(i)=interval
		u=u+hi
		goto 10
	end if
 10	continue


c
c     each component of the final output vector is the sum of
c     five numbers, all stored in the work vector.
c     three of the five sets of numbers that sum to the final 
c     vector are calculated in the next part of the program, using
c     the values calculated above.
c     the first set (a) is stored in work(ni+1..ni+no), the second
c     set (b) is stored in (1+ni+no..ni+2*no), the third set(c) is 
c     stored in (1+ni+2*no..ni+3*no).  
c     the components of each of these sets initially equals the sum of
c     the elements in the input vector in a specific interval of
c     the output vector multiplied by dx itself(c), dx squared(b),
c     or dx cubed(a).
c
	do 20 i=1,ni
	j=work(i)
	dx=x1i+(i-1)*hi-(x1o+(j-1)*ho)
	work(j+ni)=work(j+ni)+vi(i)*dx**3
	work(j+ni+no)=work(j+ni+no)+vi(i)*dx**2
	work(j+ni+2*no)=work(j+ni+2*no)+vi(i)*dx
 20	continue

	temp1=work(ni+no-1)
	temp2=work(ni+3*no-1)

	do 30 i=no-1,2,-1
	work(i+ni)=work(i+ni-1)/(6*ho)-work(i+ni)/(6*ho)
	work(i+ni+no)=work(i+ni+no)/2
	work(i+ni+2*no)=(ho/6)*work(ni+2*no+i-1)+(2*ho)/6*work(ni+2*no+i)
 30	continue

	work(ni+1)=work(ni+1)*(-1/(6*ho))
	work(ni+no)=temp1/(6*ho)
	work(ni+no+1)=work(ni+no+1)/2
	work(ni+2*no)=0.0
	work(ni+2*no+1)=(2*ho)/6*work(ni+2*no+1)
	work(ni+3*no)=(ho/6)*temp2



c
c     now I make changes based on the end conditions which
c     transform the length of the sets of numbers to no-2
c     [(start+2...start+no-1) instead of (start+1...start+no)]
c
c	for iend=1:  no changes
c	for iend=2,
c
	if (iend.eq.2) then
		work(ni+2)=work(ni+1)+work(ni+2)
		work(ni+no-1)=work(ni+no-1)+work(ni+no)
		work(ni+no+2)=work(ni+no+1)+work(ni+no+2)
		work(ni+2*no-1)=work(ni+2*no-1)+work(ni+2*no)
		work(ni+2*no+2)=work(ni+2*no+1)+work(ni+2*no+2)
		work(ni+3*no-1)=work(ni+3*no-1)+work(ni+3*no)
	end if


	if (iend.eq.3) then
		work(ni+2)=work(ni+2)+2*work(ni+1)
		work(ni+3)=work(ni+3)-work(ni+1)
		work(ni+no-2)=work(ni+no-2)-work(ni+no)
		work(ni+no-1)=work(ni+no-1)+2*work(ni+no)
		work(ni+no+2)=work(ni+no+2)+2*work(ni+no+1)
		work(ni+no+3)=work(ni+no+3)-work(ni+no+1)
		work(ni+2*no-2)=work(ni+2*no-2)-work(ni+2*no)
		work(ni+2*no-1)=work(ni+2*no-1)+2*work(ni+2*no)
	 	work(ni+2*no+2)=work(ni+2*no+2)+2*work(ni+2*no+1)
		work(ni+2*no+3)=work(ni+2*no+3)-work(ni+2*no+1)
		work(ni+3*no-2)=work(ni+3*no-2)-work(ni+3*no)
		work(ni+3*no-1)=work(ni+3*no-1)+2*work(ni+3*no)
	end if

c
c     for the three sets of numbers that I have been dealing with,
c     [work(ni+1..ni+no),work(ni+no+1..ni+2*no) and work(ni+2*no+1..
c     ni+3*no)], the next step is for each of them to be multiplied
c     by the inverse of a certain matrix.  Since I know the matrix,
c     rather than its inverse, gaussian elimination is used.  
c
c     the components of the tridiagonal matrix are stored in 
c     work(ni+3*no+1..ni+4*no-2),work(ni+4*no-1..ni+5*no-4),
c     and work(ni+5*no-3..ni+6*no-6).
c
c     the 3 trials represent the multiplication of the inverse
c     on the matrix on the three sets of numbers.  In each trial, 
c     the set of numbers is stored in work(ni+6*no-5...ni+7*no-8).
c

	do 100 trial=1,3	
	do 40 i=1,no-2
	work(ni+3*no+i)=ho
	work(ni+4*no-2+i)=4*ho
	work(ni+5*no-4+i)=ho
	if (trial.eq.1) then
	   work(ni+6*no-6+i)=work(i+1+ni)
        else if (trial.eq.2) then
	   work(ni+6*no-6+i)=work(i+1+ni+no)
	else if (trial.eq.3) then
	   work(ni+6*no-6+i)=work(i+1+ni+2*no)
	end if
 40	continue
	first=2
	last=no-2

c
c     again, I need to make adjustments for end conditions
c
c     iend=1  no changes
c

	if (iend.eq.2) then
		work(1+ni+4*no-2)=work(1+ni+4*no-2)+ho 
		work(ni+5*no-4)=work(ni+5*no-4)+ho
	else if (iend.eq.3) then
		work(ni+4*no-1)=6*ho
		work(ni+3*no+2)=0.0
		work(ni+6*no-7)=0
		work(ni+5*no-4)=6*ho
	end if

c
c     this is the gaussian elimination
c
	do 60 i=first,last
	work(ni+3*no+i)=work(ni+3*no+i)/work(ni+4*no+i-3)
	work(ni+4*no-2+i)=work(ni+4*no-2+i)-work(ni+3*no+i)*
     1			  work(ni+5*no+i-5)
	work(ni+6*no-6+i)=work(ni+6*no-6+i)-work(ni+3*no+i)*
     1			  work(ni+6*no+i-7)
 60	continue
	work(last+ni+6*no-6)=work(last+ni+6*no-6)/work(last+ni+4*no-2)


	do 70 j=last-1,first-1,-1
	work(j+ni+6*no-6)=(work(j+ni+6*no-6)-work(j+ni+5*no-4)*
     1			   work(j+ni+6*no-5))/work(j+ni+4*no-2)
 70	continue


	do 80 i=first-1,last
	if (trial.eq.1) then
	work(ni+i)=work(ni+6*no-6+i)
	else if (trial.eq.2) then
	work(ni+no+i)=work(ni+6*no-6+i)
	else if (trial.eq.3) then
	work(ni+2*no+i)=work(ni+6*no-6+i)
	end if
 80	continue
 100	continue

c
c     to get the final values for the three sets of numbers,
c     each of the sets must be multiplied by the same matrix.
c     this is accomplished in the following steps.
c


	temp3=work(ni+no-2)
	temp4=work(ni+2*no-2)
	temp5=work(ni+3*no-2)
	temp6=work(ni+no-3)
	temp7=work(ni+2*no-3)
	temp8=work(ni+3*no-3)

	do 110 i=no-2,3,-1
	work(ni+i)=6/ho*work(ni+i-2)-12/ho*work(ni+i-1)+6/ho*work(ni+i)
	work(ni+no+i)=6/ho*work(ni+no+i-2)-12/ho*work(ni+no+i-1)+
     1			6/ho*work(ni+no+i)
	work(ni+2*no+i)=6/ho*work(ni+2*no+i-2)-12/ho*work(ni+2*no+i-1)+
     1			6/ho*work(ni+2*no+i)
 110	continue


	work(ni+2)=-12/ho*work(ni+1)+6/ho*work(ni+2)
	work(ni+1)=6/ho*work(ni+1)
	work(ni+no-1)=6/ho*temp6-12/ho*temp3
	work(ni+no)=6/ho*temp3
	work(ni+no+2)=-12/ho*work(ni+no+1)+6/ho*work(ni+no+2)
	work(ni+no+1)=6/ho*work(ni+no+1)
	work(ni+2*no-1)=6/ho*temp7-12/ho*temp4
	work(ni+2*no)=6/ho*temp4
	work(ni+2*no+2)=-12/ho*work(ni+2*no+1)+6/ho*work(ni+2*no+2)
	work(ni+2*no+1)=6/ho*work(ni+2*no+1)
	work(ni+3*no-1)=6/ho*temp8-12/ho*temp5
	work(ni+3*no)=6/ho*temp5

c
c     now, the object is to determine the remaining two sets of
c     numbers that sum with the three already calculated to 
c     produce the output vector.
c     for efficiency, I will 'recycle' the part of the workspace 
c     used in calculating the previous three sets that is
c     no longer needed.  
c     Thus, I need to zero out this part of the work vector.
c

	do 115 i=ni+3*no+1,ni+7*no-8
	work(i)=0.0
 115	continue

c
c     the two sets of numbers to be calculated are stored in
c     work(ni+3*no+1..ni+4*no)(d) and work(ni+4*no..ni+5*no)(e).
c     the components of each set should equal the sum of the 
c     elements in the input vector in a specific interval of the
c     output vector (d) or this amount multiplied by dx(e).
c     

	do 120 i=1,ni
	j=work(i)
	work(ni+3*no+j)=work(ni+3*no+j)+vi(i)
	dx=x1i+(i-1)*hi-(x1o+(j-1)*ho)
	work(ni+4*no+j)=work(ni+4*no+j)+vi(i)*dx
 120	continue

c
c     the last set (e) is also multiplied by a matrix, as reflected
c     in the following calculations
c
	do 130 i=no-1,2,-1
	work(ni+4*no+i)=work(ni+4*no+i-1)/ho-work(ni+4*no+i)/ho
 130	continue
	work(ni+4*no+1)=-1/ho*work(ni+4*no+1)
	work(ni+5*no)=temp2/ho

c
c     now all five sets of numbers can be summed to produce my
c     output vector.
c
	do 140 i=1,no
	vo(i)=work(ni+3*no+i)+work(ni+i)+work(ni+no+i)-work(ni+2*no+i)
     1		+ work(ni+4*no+i)

 140	continue
	end if
	return
	end

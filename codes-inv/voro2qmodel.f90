subroutine voro2qmodel(voro,npt,npt_max,d_min,d_max,beta,h,vpvs,qa,qb,alpha,rho)

	real voro(npt_max,3),d_min,d_max
	real beta(npt_max),alpha(npt_max),vpvs(npt_max),rho(npt_max),qa(npt_max),qb(npt_max)
	real maxx,minn,summ,h(npt_max) 
	integer npt,ind,i,j,order(npt),npt_max

	beta=0
	vpvs=0
	qa=0
	qb=0
	h=0

!	do i=1,npt
!	write(*,*) voro(i,:)
!	enddo

	ind=1


	do i=1,npt
		maxx=d_max
		if (i==1) then
			minn=d_min
		else
			minn=voro(order(i-1),1)
		endif
!		write(*,*) "****",minn,maxx
	do j=1,npt
		if ((minn<voro(j,1)).and.(voro(j,1)<=maxx)) then
!			write(*,*)minn,voro(j,1),maxx
			ind=j
			maxx=voro(j,1)
!			write(*,*)ind,maxx
		endif
	enddo
	order(i)=ind

	enddo


	summ=0
	do i=1,npt-1

	h(i)= (voro(order(i),1)+voro(order(i+1),1))/2 - ssumm
	ssumm = ssumm + h(i)
	beta(i) = voro(order(i),2)
        vpvs(i) = voro(order(i),3)
	if (vpvs(i)<=0) then
        vpvs(i)=-999.
	alpha(i) = 1.50
	rho(i) = 1.00
        else
	vpvs(i) = voro(order(i),3)
	alpha(i)=vpvs(i)*beta(i)
	rho(i)= 1.6612*alpha(i)-0.4721*alpha(i)*alpha(i)+&
		&0.0671*alpha(i)*alpha(i)*alpha(i)-&
		&0.0043*alpha(i)*alpha(i)*alpha(i)*alpha(i)+&
		&0.000106*alpha(i)*alpha(i)*alpha(i)*alpha(i)*alpha(i)
        endif
	qa(i)= 1450
	qb(i)= 600 
	enddo

	h(npt)= 0.0
	beta(npt) = voro(order(npt),2)
	!!!jiahui 
	if (beta(npt)<=4.3) then
		beta(npt)=4.305
	endif
	
	vpvs(npt) = 1.85
	alpha(npt) = beta(npt)*vpvs(npt)
	rho(npt)= 1.6612*alpha(npt)-0.4721*alpha(npt)*alpha(npt)+&
		&0.0671*alpha(npt)*alpha(npt)*alpha(npt)-&
		&0.0043*alpha(npt)*alpha(npt)*alpha(npt)*alpha(npt)+&
		&0.000106*alpha(npt)*alpha(npt)*alpha(npt)*alpha(npt)*alpha(npt)

	qa(npt)= 1450
	qb(npt)= 600 
	
return
end

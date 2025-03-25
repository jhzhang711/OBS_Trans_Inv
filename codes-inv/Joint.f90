program Joint_inversion
  use cls_param
  use cls_vmodel
  use cls_disper
  use cls_recv_func
  implicit none
	include 'mpif.h'

	real , EXTERNAL    ::    gasdev,ran3
	real log, sqrt
	integer, parameter :: inv_type = 1 ! -1 for RF; -2 for SWD; 1 for RF+SWD
	integer, parameter :: ndatar = 300  ! sampling points for RF
	integer, parameter :: ndatad = 46   ! sampling points for SWD_phase
	integer, parameter :: ndataHoV = 14   ! sampling points for SWD_group
	real, parameter :: water_depth = 5.0 ! water depth for OBS station
	real, parameter :: w1 = 0 ! weights for RF1 0.2
	real, parameter :: w2 = 0.6 ! weights for RF2 0.55
	real, parameter :: w3 = 0.4 ! weights for RF3 0.25
	real, parameter :: gauss_a1 = 1.0 ! RF1 with gaussian parameter of 1.0
	real, parameter :: gauss_a2 = 2.0 ! RF2 with gaussian parameter of 2.0
	real, parameter :: gauss_a3 = 3.0 ! RF3 with gaussian parameter of 3.0
	double precision, parameter ::    Ar_max=0.15 !Upper bound for noise parameter for RF
	double precision, parameter ::    Ad_max=0.02 !Upper bound for noise parameter for SWD_phase
	double precision, parameter ::    AHoV_max=0.02 !Upper bound for noise parameter for SWD_group
	real, parameter :: pd1 = 0.5    ! proposal on change in position  
	real, parameter :: pv1 = 0.1    ! proposal on velocity
	real, parameter :: pd2 = 2.0    ! proposal on change in position    
	real, parameter :: pv2 = 0.1    ! proposal on velocity  
	real, parameter :: pvpvs = 0.02

! Each chain is run for 'burn_in + nsample' steps in total. The first burn-in samples are discarded as burn-in steps, only after which the sampling algorithm is assumed to have converged. To eliminate dependent samples in the ensemble solution, every thin model visited in the second part is selected for the ensemble. The convergence of the algorithm is monitored with a number of indicators such as acceptance rates, and sampling efficiency is optimized by adjusting the variance of the Gaussian proposal functions 
	integer, parameter :: burn_in = 5500 !Burn-in period
	integer, parameter :: nsample = 5500 !Post burn-in
	integer, parameter :: thin = 2 !Thinning of the chain 
	integer, parameter :: display = 500

!Proposals
! sigmav is the most tricky to "tune". Birth and Death jumps are difficult to be accepted (often < 10%)
! Try to find a sigmav that maximizes Birth and Death acceptances rates.
	real, parameter :: sigmav = 0.05  ! proposal on velocity when Birth move  
	real, parameter :: pAr = 0.0005    ! proposal for change in noise parameter       
	real, parameter :: pAd = 0.0005    ! proposal for change in noise parameter    
	real, parameter :: pAHoV = 0.0005  ! proposal for change in HoV parameter      
	real, parameter :: ped = 0.0005
	double precision, parameter ::    sig=2!fs/a


! The algorithm is Hierarchical Bayes, and the level of data errors 
! (the required level of data fit) is an unknown paramete rto be inverted for. 
! Here the noise is supposed indepemdant and Cd is diagonal.
! A Uniform prior distribution is defined for two parameter (RF and SWD). 
! The noise parameter for SWD is a % of observed value. That is for Vphase(i), the error is Vphase(i)*Ad/100
	double precision, parameter ::    Ar_min=0.02
	double precision, parameter ::    Ad_min=0.01
	double precision, parameter ::    AHoV_min=0.01
	double precision, parameter ::    epsd_max=0.0025
	double precision, parameter ::    epsd_min=0.001

!------------------------------------------------
! Prior distribution (Bounds odf the model space)
!------------------------------------------------
	real, parameter :: d_min = 0      ! depth 
	real, parameter :: d_max = 120      
	real, parameter :: width = 0.0    ! Lower and upper bound of the prior
	real, parameter :: beta_min = 0.0
	real, parameter :: beta_max = 5.0
	real, parameter :: vpvs_min = 1.5
	real, parameter :: vpvs_max = 2.2
	integer,parameter  :: npt_min = 15 ! min number of the layers in inversion
	integer, parameter :: npt_max = 30

!discretezation for the posterior distribution.
!specifies the number of velocity pixels in (Vs, depth)
!This is because the solution is visualized as an histogram, 
!and hence we need to define the number of bins

	integer, parameter :: disd = 201 !nb bins in depth
	integer, parameter :: disv =  91 !nb bins in velocity
	integer, parameter :: disA = 150 !nb bins for noise parameter  
	integer, parameter :: dise = 100
	integer, parameter :: disc = 150
	integer, parameter :: disvpvs = 50
!depth of model for display
	real, parameter :: prof = 100

!input parameters (geomtry)
	real, parameter :: fs = 10       !sampling frequency
	real, parameter :: water_c = 0.001
	real, parameter ::    pi = 3.14159265
	real, parameter ::    rad = 0.017453292 
	real, parameter :: time_shift = 5

!***************************************************************

! DECLARATIONS OF VARIABLES

!****************************************************************
	integer i,npt,sample,ind,accept,l,th,ount,npt_prop,nb,k,tes,npt_best,npt_pre
	integer birth, death,move,value,noisr,noisd,noisHoV,ind2,v,Ivpvs,j,chwidd,chwidr
	integer histo(npt_max),histos(npt_max),histoch(disd),histochs(disd)
	real beta(npt_max),h(npt_max),vpvs(npt_max),qa(npt_max),noiser(ndatar),peri(ndatad)
	real alpha(npt_max),rho(npt_max)
	real mod_iter(201,3),mod_tmp
	real periHoV(ndataHoV)
        real depth(npt_max)
	real ppara,pvpvs_min,pvpvs_max,ppara1,ppara2,ppara3 !jiahui
	real qb(npt_max),av(disd),post(disd,disv),PrP(2),PrV(2),AcP(2),AcV(2),mo(disd)
	real avpvs(disd),avvpvs(disd),postvpvs(disd,disvpvs),postsvpvs(disd,disvpvs)
	real PrB,AcB,PrD,AcD,Prnr,Prnd,PrnHoV,Acnr,Acnd,AcnHoV,avs(disd),posts(disd,disv),out,&
		&Prwd,Acwd,Prwr,Acwr,LSr,LSr_prop,LSr_min, LSd,LSd_prop,LSd_min,&
		&LSHoV,LSHoV_prop,LSHoV_min,AcwHoV,PrwHoV
	real LSr_prop1,LSr_prop2,LSr_prop3
	real  din,ht,d,logprob,dd,convAr(nsample+burn_in),convArs(nsample+burn_in),&
		&conver(nsample+burn_in),convers(nsample+burn_in),&
		&conveHoV(nsample+burn_in),conveHoVs(nsample+burn_in)
	real convAHoV(nsample+burn_in),convAHoVs(nsample+burn_in)
	real convAd(nsample+burn_in),convAds(nsample+burn_in),conved(nsample+burn_in),conveds(nsample+burn_in)
	real best_datar1(ndatar),d_obsd(ndatad),voro(npt_max,3),d_obsHoV(ndataHoV)
	real best_datar2(ndatar),best_datar3(ndatar)
	real voro_best(npt_max,3)
	real d_obsr1(ndatar),d_obsr2(ndatar),d_obsr3(ndatar)
	real wdata1(ndatar),wdata2(ndatar),wdata3(ndatar)
	real voro_prop(npt_max,3),conv(nsample+burn_in),ncell(nsample+burn_in)
	real convs(nsample+burn_in),ncells(nsample+burn_in),t1,t2
	real cvs(nsample/thin,disd),cvs1(disd),cvs2(disd),sort,cvss1(disd),cvss2(disd)
	integer q1,q2

	real like,like_prop,liked,liked_prop,liker,liker_prop,u,pv_min,pv_max,&
		&like_min,xi,period(ndatad+1),Vphase(ndatad),best_datad(ndatad),HoV(ndataHoV),&
		likeHoV,likeHoV_prop,best_dataHoV(ndataHoV),periodHoV(ndataHoV+1)
	real tmpHoV(ndatad),tmpVphase(ndataHoV)
	real liker1,liker2,liker3,liker1_prop,liker2_prop,liker3_prop
	real ML_Ar(disA),ML_Ars(disA),ML_ers(dise),ML_AHoVs(disA)
	real ML_Ad(disA),ML_Ads(disA),ML_ed(dise),ML_eds(dise)
	real ML_AHoV(disA),ML_AdHoV(disA)
	double precision bd(ndatad),bHoV(ndataHoV),Ar,Ar_prop,logrsig,logreps,&
		&AHoV,AHoV_prop,Ad,Ad_prop,epsd,epsd_prop
	double precision br1(ndatar),br2(ndatar),br3(ndatar)
!       Parameters for rj_mcmc
	double precision AId(ndatad,ndatad),besd(ndatad),besr(ndatar),AId_prop(ndatad,ndatad),&
		&AIr(ndatar,ndatar),AIr_prop(ndatar,ndatar),Cd(ndatar,ndatar),UU(ndatar,ndatar),&
		&W(ndatar),x(ndatar),vv(ndatar,ndatar),bes(ndatar),b(ndatar),wmin,besHoV(ndataHoV),&
		&AIHoV(ndataHoV,ndataHoV),AIHoV_prop(ndataHoV,ndataHoV)
	real Moy(ndatad,ndatad),Moys(ndatad,ndatad),ress(ndatad)

!	for MPI
	integer ra,ran,rank, nbproc, ierror, tag, status(MPI_STATUS_SIZE)
	integer cor(disc,disv),cors(disc,disv)
	character filename*11, number*3
	type(vmodel) :: vm,vm2
  	type(disper) :: disp1,disp2
  	type(recv_func) :: rf1,rf2,rf3
 	logical :: is_ok

	CALL cpu_time(t1) 
! Start Parralelization of the code. From now on, the code is run on each processor independently,
! with ra = the number of the proc.
	call MPI_INIT(ierror)
 	call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror)
 	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

 	nb=nbproc
        ra=rank
        ran=rank
        conv=0
        ncell=0

!*********************************************************************
! READ THE DATA
	open(55,file="RF1.obs",status='old')
        read(55,*) ppara1 !jiahui
	do i=1,ndatar
		read(55,*) u,d_obsr1(i)
	end do
	close(55)! close the file

	open(65,file="RF2.obs",status='old')
        read(65,*) ppara2 !jiahui
	do i=1,ndatar
		read(65,*) u,d_obsr2(i)
	end do
	close(65)! close the file

	open(75,file="RF3.obs",status='old')
        read(75,*) ppara3 !jiahui
	do i=1,ndatar
		read(75,*) u,d_obsr3(i)
	end do
	close(75)! close the file

	open(45,file="SWD_phase.obs",status='old')
	do i=1,ndatad
		read(45,*) peri(i),d_obsd(i)
	end do
	close(45)! close the file

	open(35,file="SWD_group.obs",status='old')
	do i=1,ndataHoV
		read(35,*) periHoV(i),d_obsHoV(i)
	end do
	close(35)! close the file

        open(85,file="REF_in.mod",status='old')
        read(85,*) npt_pre
        do i=1,npt_pre
        read(85,*) mod_iter(i,1),mod_iter(i,2),mod_iter(i,3)
        end do
        close(85)

!*********************************************************************
	Ar = Ar_min+ran3(ra)*(Ar_max-Ar_min)
	Ad = Ad_min+ran3(ra)*(Ad_max-Ad_min)
	AHoV = AHoV_min+ran3(ra)*(AHoV_max-AHoV_min)
! Initial number of cells
!------------------------------------
	j=0
	tes=0
	do while(tes==0)
		npt = npt_min+0.5*ran3(ra)*(npt_max-npt_min)
		j=j+1
		tes=1

		! the first layer would be water
                voro(1,1) = water_depth-0.01 ! water layer  
                voro(1,2) = -3.0 !
                voro(1,3) = -999 ! 

                voro(2,1) = water_depth+0.01 ! water layer  
		call priorvalue(voro(2,1),water_depth,pv_min,pv_max)
                voro(2,2) = pv_min+(pv_max-pv_min)*ran3(ra)
                voro(2,3) = 1.78

!		do i=3,npt-1
!			voro(i,1)=water_depth+0.01+ran3(ra)*(d_max-water_depth)
!			call priorvalue(voro(i,1),water_depth,pv_min,pv_max)
!			voro(i,2)=pv_min+(pv_max-pv_min)*ran3(ra)
!                       voro(i,3) = 1.78
!			! voro(:,1) depth of upper disconotinuity
!			! voro(:,2) Vs
!			! voro(:,3) Vp/Vs ratio
!		enddo
                do i=3,npt
                       k = 10*(i-3)+floor(10*ran3(ra))
                       voro(i,1)=mod_iter(k,1)+2.0*ran3(ra)-1.0
                        if(voro(i,1) < water_depth+0.5) then
                                voro(i,1) = water_depth+0.5
                        endif
                       voro(i,2)=mod_iter(k,2)+0.5*ran3(ra)-0.25
                       voro(i,3) = 1.78
                end do
!------------------------------------
		! Set the last half-space layer
		do i=npt+1,npt_max
			voro(i,1)=0.0
			voro(i,2)=voro(i-1,2)
			voro(i,3)=voro(i-1,3)
		enddo
!------------------------------------
		call voro2qmodel(voro,npt,npt_max,d_min,d_max,beta,h,vpvs,qa,qb,alpha,rho)

  		vm = init_vmodel()
        	call vm%setvmodel(npt,alpha,beta,rho,h,vm2)

!------------------------------------
!		din=asin(ppara*4.25*1.8)/rad
!		call theo(&
!			npt-1, beta(2:npt), h(2:npt), vpvs(2:npt), qa(2:npt), qb(2:npt), fs, din,&
!			gauss_a1, water_c, time_shift, ndatar, wdata1 )
!		call theo(&
!			npt-1, beta(2:npt), h(2:npt), vpvs(2:npt), qa(2:npt), qb(2:npt), fs, din,&
!			gauss_a2, water_c, time_shift, ndatar, wdata2 )
!		call theo(&
!			npt-1, beta(2:npt), h(2:npt), vpvs(2:npt), qa(2:npt), qb(2:npt), fs, din,&
!			gauss_a3, water_c, time_shift, ndatar, wdata3 )

 		rf1 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara1,8), &
       			& a_gauss = real(gauss_a1,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf1%compute()

  		rf2 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara2,8), &
       			& a_gauss = real(gauss_a2,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf2%compute()

 		rf3 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara3,8), &
       			& a_gauss = real(gauss_a3,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf3%compute()

		do i = 1,ndatar
			wdata1(i)=rf1%get_rf_data(i)
			wdata2(i)=rf2%get_rf_data(i)
			wdata3(i)=rf3%get_rf_data(i)
		end do

		! phase 
		disp1 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(peri(1),8), &
       			& xmax   = real(peri(ndatad),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(1.0,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp1%dispersion(is_ok=is_ok)
		do i = 1,ndatad
			Vphase(i)=disp1%get_c(i)
		end do

		! group 
		disp2 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(periHoV(1),8), &
       			& xmax   = real(periHoV(ndataHoV),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(0.5,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp2%dispersion(is_ok=is_ok)
		do i = 1,ndataHoV
			HoV(i)=disp2%get_u(i)
		end do

		
                if (disp1%get_c(1) <=0.0 ) tes=0
                if (disp2%get_c(1) <=0.0 ) tes=0
		!do i=1,ndatar
		!	if (isnan(wdata3(i))) tes=0
		!enddo

	! Get the data covariance Matrix for dispersion
		epsd =  0.0
		do i=1,ndatad
			do j=1,ndatad
				AId(i,j)=0
				if (i==j) AId(i,j)=1+epsd**2
				if (abs(i-j)==1) AId(i,j)=-epsd
			enddo		
		enddo
		AId(1,1)=1
		AId(ndatad,ndatad)=1
		AId=AId/(1-epsd**2)
		do i=1,ndatad
 			bd(i)=(d_obsd(i)-Vphase(i))
		enddo
		liked=0
		do i=1,ndatad
 			besd(i)=0
 			do j=1,ndatad
				besd(i)=besd(i)+AId(i,j)*bd(j)
			enddo
			liked=liked+besd(i)*bd(i)
		enddo
		liked=liked/(2*Ad**2)

		if (liked>=100000.0) tes=0


		if (j>1000) write(*,*) ran,j,npt
		if (j>1001) then
			write(*,*)'STOP'	
			stop
		endif
	enddo

!***************************************************
! Get the data covariance Matrix for RF
	do i=1,ndatar
		do j=1,ndatar
			Cd(i,j)=exp((-(i-j)**2)/(2*sig**2))
			if (Cd(i,j)<0.0001) Cd(i,j)=0
		enddo		
	enddo		
	UU=Cd
 	call dsvdcmp(UU,ndatar,ndatar,ndatar,ndatar,w,vv)
! Remove small eigenvalues
	wmin=0.000001
	do i=1,ndatar
		if (w(i)<wmin) w(i)=0
	enddo
! test the svd decomposition ability to solve the equation Cd *x = b
! construct a b 
	do i=1,ndatar
 		b(i)=gasdev(ra)
	enddo
	call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,b,x)
	do i=1,ndatar
 		bes(i)=0
         	do j=1,ndatar
         		bes(i)=bes(i)+x(j)*Cd(i,j)
        	enddo
	enddo
	like=0
	do i=1,ndatar
		like=like+(b(i)-bes(i))**2
	enddo
	like=like/ndatar
	write(*,*)'... Test SVD decomposition',sqrt(like),'should be << 1'

!Get initial likelihood----------------------------------
	LSr=0
	do i=1,ndatar
 		br1(i)=(d_obsr1(i)-wdata1(i))
 		LSr=LSr+br1(i)**2
 	enddo
	LSr_min=LSr
 	br1=br1/(Ar**2)
! MIsfit INIT with SVD
	liker1=0
	do i=1,ndatar
		liker1=liker1+br1(i)**2
	enddo
!//////
	LSr=0
	do i=1,ndatar
 		br2(i)=(d_obsr2(i)-wdata2(i))
 		LSr=LSr+br2(i)**2
 	enddo
	LSr_min=LSr
 	br2=br2/(Ar**2)
! MIsfit INIT with SVD
	call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,br2,x)
	liker2=0
	do i=1,ndatar
		liker2=liker2+br2(i)**2
	enddo
!//////
	LSr=0
	do i=1,ndatar
 		br3(i)=(d_obsr3(i)-wdata3(i))
 		LSr=LSr+br3(i)**2
 	enddo
	LSr_min=LSr
 	br3=br3/(Ar**2)
! MIsfit INIT with SVD
	liker3=0
	do i=1,ndatar
		liker3=liker3+br3(i)**2
	enddo
!//////

	liker=(w1*liker1+w2*liker2+w3*liker3)/(2)
	write(*,*)'RF  likelihood',liker,"(",ran,")"

!***************************************************
! Get the data covariance Matrix for dispersion
	epsd =  epsd_min+ran3(ra)*(epsd_max-epsd_min)
	do i=1,ndatad
		do j=1,ndatad
			AId(i,j)=0
			if (i==j) AId(i,j)=1+epsd**2
			if (abs(i-j)==1) AId(i,j)=-epsd
		enddo		
	enddo
	AId(1,1)=1
	AId(ndatad,ndatad)=1
	AId=AId/(1-epsd**2)
	Lsd=0
	do i=1,ndatad
 		bd(i)=(d_obsd(i)-Vphase(i))
 		LSd=LSd+bd(i)**2
	enddo
	LSd_min=LSd
	liked=0
	do i=1,ndatad
 		besd(i)=0
 		do j=1,ndatad
			besd(i)=besd(i)+AId(i,j)*bd(j)
		enddo
		liked=liked+besd(i)*bd(i)
	enddo
	liked=liked/(2*Ad**2)
	write(*,*)'SWD likelihood',liked,"(",ran,")"
!***************************************************
! Get the data covariance Matrix for HoV
	epsd =  epsd_min+ran3(ra)*(epsd_max-epsd_min)
	do i=1,ndataHoV
		do j=1,ndataHoV
			AIHoV(i,j)=0
			if (i==j) AIHoV(i,j)=1+epsd**2
			if (abs(i-j)==1) AIHoV(i,j)=-epsd
		enddo		
	enddo
	AIHoV(1,1)=1
	AIHoV(ndataHoV,ndataHoV)=1
	AIHoV=AIHoV/(1-epsd**2)
	LsHoV=0
	do i=1,ndataHoV
 		bHoV(i)=(d_obsHoV(i)-HoV(i))
 		LSHoV=LSHoV+bHoV(i)**2
	enddo
	LSHoV_min=LsHoV
	likeHoV=0
	do i=1,ndataHoV
 		besHoV(i)=0
 		do j=1,ndataHoV
			besHoV(i)=besHoV(i)+AIHoV(i,j)*bHoV(j)
		enddo
		likeHoV=likeHoV+besHoV(i)*bHoV(i)
	enddo
	likeHoV=likeHoV/(2*AHoV**2)
	write(*,*)'HoV likelihood',likeHoV,"(",ran,")"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (inv_type == -1) then
		like=liker
	else if (inv_type == -2) then
		like=liked + likeHoV
	else if (inv_type == 1) then
		like = 0.9*liker + 1.1*liked + likeHoV !jiahui
	else
		stop
	endif 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*******************************************************************************
!\\\\\\\\\\//////////\\\\\\\\\\//////////\\\\\\\\\\//////////\\\\\\\\\/////////*
!*******************************************************************************
!!!!!! START THE SAMPLING

	av=0
	avpvs = 0
	sample=0
	th=0
	ount=0
	PrP=0
	PrV=0
	PrB=0
	PrD=0
	AcP=0
	AcV=0
	AcB=0
	AcD=0
	post=0
	postvpvs = 0

	Acnr=0
	Prnr=0
	Acnd=0
	Prnd=0
	AcnHoV=0
	PrnHoV=0

	Prwd=0
	Acwd=0
	Prwr=0
	Acwr=0
	PrwHoV=0
	AcwHoV=0

	histo=0
	histoch=0
 	cor = 0
	ML_Ad=0
	ML_Ar=0
	ML_AHoV=0
	Moy=0

	best_datar1 = wdata1
	best_datar2 = wdata2
	best_datar3 = wdata3
	best_datad = Vphase
	best_dataHoV = HoV
	like_min = like

!	write(number,1000)ran
!1000  format(I3)
!	filename=number//'model.dat'
!	open(ran,file='OUT/'//filename,status='replace')


	do while (sample<nsample)
		ount=ount+1
		voro_prop=voro
		like_prop=like
		npt_prop = npt
		Ar_prop = Ar
		AIr_prop = AIr
		Ad_prop = Ad
		AId_prop = AId
		AHoV_prop = AHoV
		AIHoV_prop = AIHoV

		u=ran3(ra)

		out=1
		move=0
		value=0
		birth=0
		death=0
		noisr=0
		noisd=0
		noisHoV=0
		chwidd=0
		chwidr=0
		logprob=0

!*************************************************************

!                   Propose a new model

!*************************************************************

	if (u<0.3) then !change position------------------
                move=1
                ind=1+ceiling(ran3(ra)*(npt-1))
                do while (voro(ind,1)>=60.0 .AND. i<=100)
                        ind=1+ceiling(ran3(ra)*(npt-1))
                enddo
                if (ind<=3) ind=3
                if (ind>=npt) ind=npt
                if (ount.GT.burn_in) then
                        if (voro(ind,1)<(d_max/2)) then
                                PrP(1)=PrP(1)+1
                        else
                                PrP(2)=PrP(2)+1
                        endif
                endif
                if (voro(ind,1)<(d_max/2)) then
                         voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd1
                else
                         voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd2
                endif
                !Check if oustide bounds of prior
                if ((voro_prop(ind,1)<=voro_prop(2,1)).or.(voro_prop(ind,1)>d_max)) then
                        out=0
                endif
                call priorvalue(voro_prop(ind,1),water_depth,pv_min,pv_max)
                if ((voro_prop(ind,2)<=pv_min).or.(voro_prop(ind,2)>=pv_max)) then
                        out=0
                endif
	
        elseif (u<0.4) then
		move=1
		ind=1+ceiling(ran3(ra)*(npt-1))
                if (ind<=3) ind=3
                if (ind>=npt) ind=npt
		if (ount.GT.burn_in) then 
			if (voro(ind,1)<(d_max/2)) then
				PrP(1)=PrP(1)+1
			else
				PrP(2)=PrP(2)+1
			endif
		endif
		if (voro(ind,1)<(d_max/2)) then
			 voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd1
		else
			 voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd2
		endif
		!Check if oustide bounds of prior
		if ((voro_prop(ind,1)<=voro_prop(2,1)).or.(voro_prop(ind,1)>d_max)) then
			out=0
		endif
		call priorvalue(voro_prop(ind,1),water_depth,pv_min,pv_max)
		if ((voro_prop(ind,2)<pv_min).or.(voro_prop(ind,2)>pv_max)) then
			out=0
		endif
        elseif (u<0.475) then
		value=1
		ind=ceiling(ran3(ra)*(npt))
                if (ind<=2) ind=2
                if (ind>=npt) ind=npt
		if (voro(ind,1)<(d_max/2)) then
			PrV(1)=PrV(1)+1
		else
			PrV(2)=PrV(2)+1
		endif
		voro_prop(ind,3)=voro(ind,3)+gasdev(ra)*pvpvs
		!Check if oustide bounds of prior
                call priorvpvs(voro_prop(ind,1),voro_prop(ind,2),water_depth,pvpvs_min,pvpvs_max)
		if ((voro_prop(ind,3)<pvpvs_min).or.(voro_prop(ind,3)>pvpvs_max)) then
			out=0
		endif
	elseif (u<0.50) then ! Change noise parameter for receiver function
		noisr=1
		Prnr = Prnr + 1
		Ar_prop = Ar+gasdev(ra)*pAr
		!Check if oustide bounds of prior
		if ((Ar_prop<Ar_min).or.(Ar_prop>Ar_max)) then
			out=0
		endif
 	elseif (u<0.525) then ! Change noise parameter for dispersion
 		noisd=1
 		Prnd = Prnd + 1
 		Ad_prop = Ad+gasdev(ra)*pAd
 		!Check if oustide bounds of prior
 		if ((Ad_prop<Ad_min).or.(Ad_prop>Ad_max)) then
 			out=0
 		endif
 	elseif (u<0.55) then ! Change noise parameter for HoV
 		noisHoV=1
 		PrnHoV = PrnHoV + 1
 		AHoV_prop = AHoV+gasdev(ra)*pAd
 		!Check if oustide bounds of prior
 		if ((AHoV_prop<AHoV_min).or.(AHoV_prop>AHoV_max)) then
 			out=0
 		endif
	elseif (u<0.9) then ! change value---------------------------------------------------------
		value=1
		ind=1+ceiling(ran3(ra)*(npt-1))
                if (ind<=2) ind=2
                if (ind>=npt) ind=npt
		if (voro(ind,1)<(d_max/2)) then
			PrV(1)=PrV(1)+1
		else
			PrV(2)=PrV(2)+1
		endif

		if (voro(ind,1)<(d_max/2)) then
			voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv1
		else
			voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv2
		endif

		!Check if oustide bounds of prior
		call priorvalue(voro_prop(ind,1),water_depth,pv_min,pv_max)
		if ((voro_prop(ind,2)<pv_min).or.(voro_prop(ind,2)>pv_max)) then
			out=0
		endif
	elseif (u<0.95) then !Birth----------------------------------------
		birth = 1
		PrB = PrB + 1
		npt_prop = npt + 1		
		if (npt_prop>npt_max) then
				out=0
		else
			voro_prop(npt_prop,1) = water_depth+0.5+ran3(ra)*(d_max-water_depth)
			call whichcell(voro_prop(npt_prop,1),voro,npt,npt_max,ind)
			voro_prop(npt_prop,2) = voro(ind,2)+gasdev(ra)*sigmav
                        voro_prop(npt_prop,3) = voro(ind,3)+gasdev(ra)*0.1
			logprob=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(npt_prop,2))**2)/(2*sigmav**2)
			!Check bounds					
			call priorvalue(voro_prop(npt_prop,1),water_depth,pv_min,pv_max)
			if ((voro_prop(npt_prop,2)<pv_min).or.(voro_prop(npt_prop,2)>pv_max)) out=0
                        call priorvpvs(voro_prop(npt_prop,1),voro_prop(npt_prop,2),water_depth,pvpvs_min,pvpvs_max)
		        if ((voro_prop(npt_prop,3)<pvpvs_min).or.(voro_prop(npt_prop,3)>pvpvs_max)) then
		          out=0
		        endif
			
		endif

	else !death!---------------------------------------	
		death = 1
		PrD = PrD + 1
		ind=ceiling(ran3(ra)*npt)
                if (ind<=3) ind=3
                if (ind>=npt) ind=npt
		npt_prop=npt-1
		if (npt_prop<npt_min) then
		 out=0
		else
			voro_prop(ind,:)=voro(npt,:)
			call whichcell(voro(ind,1),voro_prop,npt_prop,npt_max,ind2)
			logprob=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(ind2,2))**2)/(2*sigmav**2)
		endif
	endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!********Compute like_prop
	if (out==1) then

		! Call Forward MOdel ------------------------------
		call voro2qmodel(voro_prop,npt_prop,npt_max,d_min,d_max,beta,h,vpvs,qa,qb,alpha,rho)
  		vm = init_vmodel()
        	call vm%setvmodel(npt_prop,alpha,beta,rho,h,vm2)


!		din=asin(ppara*4.25*1.8)/rad
!		call theo(&
!		npt_prop-1, beta(2:npt_prop), h(2:npt_prop), vpvs(2:npt_prop), qa(2:npt_prop),&
!               qb(2:npt_prop), fs, din,gauss_a1, water_c, time_shift, ndatar, wdata1 )
!                call theo(&
!                npt_prop-1, beta(2:npt_prop), h(2:npt_prop), vpvs(2:npt_prop), qa(2:npt_prop),&
!                qb(2:npt_prop), fs, din,gauss_a2, water_c, time_shift, ndatar, wdata2 )
!                call theo(&
!                npt_prop-1, beta(2:npt_prop), h(2:npt_prop), vpvs(2:npt_prop), qa(2:npt_prop),&
!                qb(2:npt_prop), fs, din,gauss_a3, water_c, time_shift, ndatar, wdata3 )
		!call dispersion(npt_prop,beta,vpvs,h,period,Vphase,tmpHoV,peri,ndatad)
		!call dispersion(npt_prop,beta,vpvs,h,periodHoV,tmpVphase,HoV,periHoV,ndataHoV)

 		rf1 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara1,8), &
       			& a_gauss = real(gauss_a1,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf1%compute()

  		rf2 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara2,8), &
       			& a_gauss = real(gauss_a2,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf2%compute()

 		rf3 = recv_func(&
       			& vm    = vm2, &
       			& n     = ndatar, &
       			& delta = real(1/fs,8), &
       			& rayp  = real(ppara3,8), &
       			& a_gauss = real(gauss_a3,8), &
       			& rf_phase = "P", &
      			& deconv_flag = .true., &
       			& t_pre = real(time_shift,8), &
       			& correct_amp = .true., &
       			& noise_added = 0.d0, &
       			& damp = real(water_c,8) &
       			& )
  		call rf3%compute()

		do i = 1,ndatar
			wdata1(i)=rf1%get_rf_data(i)
			wdata2(i)=rf2%get_rf_data(i)
			wdata3(i)=rf3%get_rf_data(i)
		end do
		! phase 
		disp1 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(peri(1),8), &
       			& xmax   = real(peri(ndatad),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(1.0,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp1%dispersion(is_ok=is_ok)
		do i = 1,ndatad
			Vphase(i)=disp1%get_c(i)
		end do

		! group 
		disp2 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(periHoV(1),8), &
       			& xmax   = real(periHoV(ndataHoV),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(0.5,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp2%dispersion(is_ok=is_ok)
		do i = 1,ndataHoV
			HoV(i)=disp2%get_u(i)
		end do

                if (disp1%get_c(1) <=0.0 ) out=0
                if (disp2%get_c(1) <=0.0 ) out=0
		!do i=1,ndatar
		!	if (isnan(wdata3(i))) out=0
		!enddo
		!--------------------------------------------------
		!--------compute liker_prop -----------------------
		LSr_prop1=0
		do i=1,ndatar
			br1(i)=(d_obsr1(i)-wdata1(i))
			LSr_prop1=LSr_prop1+br1(i)**2
		enddo
			br1=br1/(Ar_prop**2)
		liker1_prop=0
		do i=1,ndatar
			liker1_prop=liker1_prop+br1(i)**2
		enddo
		!/////////////
		LSr_prop2=0
		do i=1,ndatar
			br2(i)=(d_obsr2(i)-wdata2(i))
			LSr_prop2=LSr_prop2+br2(i)**2
		enddo
			br2=br2/(Ar_prop**2)
		call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,br2,x)
		liker2_prop=0
		do i=1,ndatar
			liker2_prop=liker2_prop+br2(i)**2
		enddo
		!/////////////
		LSr_prop3=0
		do i=1,ndatar
			br3(i)=(d_obsr3(i)-wdata3(i))
			LSr_prop3=LSr_prop3+br3(i)**2
		enddo
			br3=br3/(Ar_prop**2)
		call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,br3,x)
		liker3_prop=0
		do i=1,ndatar
			liker3_prop=liker3_prop+br3(i)**2
		enddo
		!/////////////	
		LSr_prop = 0
		LSr_prop = w1*LSr_prop1+w2*LSr_prop2+w3*LSr_prop3
		liker_prop=(w1*liker1_prop+w2*liker2_prop+w3*liker3_prop)/(2)
		!----------------------------------------------------
		LSd_prop=0
		do i=1,ndatad
			bd(i)=(d_obsd(i)-Vphase(i))
			LSd_prop=LSd_prop+bd(i)**2
		enddo
		liked_prop=0
		do i=1,ndatad
		besd(i)=0
			do j=1,ndatad
				besd(i)=besd(i)+AId_prop(i,j)*bd(j)
			enddo
		liked_prop=liked_prop+besd(i)*bd(i)
		enddo
		liked_prop=liked_prop/(2*Ad_prop**2)
		!----------------------------------------------------
		LSHoV_prop=0
		do i=1,ndataHoV
			bHoV(i)=(d_obsHoV(i)-HoV(i))
			LSHoV_prop=LSHoV_prop+bHoV(i)**2
		enddo

		likeHoV_prop=0
		do i=1,ndataHoV
		besHoV(i)=0
			do j=1,ndataHoV
				besHoV(i)=besHoV(i)+AIHoV_prop(i,j)*bHoV(j)
			enddo
		likeHoV_prop=likeHoV_prop+besHoV(i)*bHoV(i)
		enddo

		likeHoV_prop=likeHoV_prop/(2*AHoV_prop**2)
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if (inv_type == -1) then
			like_prop= liker_prop
		else if (inv_type == -2) then
			like_prop= liked_prop + likeHoV_prop
		else if (inv_type == 1) then
			like_prop= liker_prop + liked_prop + likeHoV_prop
		else
			stop
		endif 
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	endif

!        like_his(ount,1) = liker_prop + liked_prop + likeHoV_prop
!        like_his(ount,2) = liker_prop
!        like_his(ount,3) = liked_prop
!        like_his(ount,4) = likeHoV_prop
!********************************************************************************
!  Now, depending on the type of move, compute acceptance term, and see if we accept or rejetc the model
!  See if we accept the proposed change to the model using acceptance ratio.
!   And update the Acceptance ratios, for each type of move.
 
!  Here the likelihood need to be normalised by the determinent of the matrix of datat errors ! lgsigma is the log of the ratio of the determinents of the matrix of data errors.
!********************************************************************************
	Accept = 0
	if (birth==1)then
        	if (log(ran3(ra))<log(out)-like_prop+like) then
            		accept=1
            		AcB=AcB+1
		endif
        elseif (death==1) then
		if (log(ran3(ra))<log(out)-like_prop+like) then
            		accept=1
            		AcD=AcD+1
        	endif
	elseif (noisr==1) then
 		logrsig=ndatar*log(Ar/Ar_prop)
        	! FOR SAMPLING PRIOR if (log(ran3(ra))<log(out)-like_prop+like) then
        	if (log(ran3(ra))<log(out)-like_prop+like) then
            		accept=1
			Acnr=Acnr+1
        	endif
	elseif (noisd==1) then
 		logrsig=ndatad*log(Ad/Ad_prop)
        	! FOR SAMPLING PRIOR if (log(ran3(ra))<log(out)-like_prop+like) then
        	if (log(ran3(ra))<log(out)-like_prop+like) then
            		accept=1
			Acnd=Acnd+1
        	endif
	elseif (noisHoV==1) then
 		logrsig=ndataHoV*log(AHoV/AHoV_prop)
        	! FOR SAMPLING PRIOR if (log(ran3(ra))<log(out)-like_prop+like) then
        	if (log(ran3(ra))<log(out)-like_prop+like) then
            		accept=1
			AcnHoV=AcnHoV+1
        	endif	
	else !NO JUMP
		if (log(ran3(ra))<log(out)-like_prop+like)then
			accept=1
             		if((value==1).and.(ount.GT.burn_in)) then
				if (voro(ind,1)<(d_max/2)) then
					AcV(1)=AcV(1)+1
				else
					AcV(2)=AcV(2)+1
				endif
             		elseif((move==1).and.(ount.GT.burn_in))then
				if (voro(ind,1)<(d_max/2)) then
					AcP(1)=AcP(1)+1
				else
					AcP(2)=AcP(2)+1
				endif
			endif
	     	endif!accept
        endif !KIND OF JUMP

!***********************************************************************************
!   If we accept the proposed model, update the status of the Markov Chain
!***********************************************************************************
	if (accept==1) then
        voro=voro_prop
        like=like_prop
	liker=liker_prop
	liked=liked_prop
	likeHoV=likeHoV_prop
	LSr=LSr_prop
	LSd=LSd_prop
	LSHoV=LSHoV_prop

	npt=npt_prop
	Ar=Ar_prop
	AIr=AIr_prop	
	Ad=Ad_prop
	AId=AId_prop
	AHoV=AHoV_prop
	AIHoV=AIHoV_prop

		if (LSr<LSr_min)then
			 LSr_min = LSr			
			 best_datar1=wdata1
 			 best_datar2=wdata2
 			 best_datar3=wdata3
			 voro_best=voro
			 npt_best=npt
		endif
		if (LSd<LSd_min)then
			 LSd_min = LSd			
			 best_datad=Vphase	
		endif
		if (LSHoV<LSHoV_min)then
			 LSHoV_min = LSHoV			
			 best_dataHoV=HoV	
		endif
!		if (LSr*w1+LSHoV*w2<LSr_min*w1+LSHoV_min*w2) then
!			 LSr_min = LSr			
!			 best_datar=wdata
!			 LSd_min = LSd			
!			 best_datad=Vphase
!			 LSHoV_min = LSHoV			
!			 best_dataHoV=HoV	
!		endif

	endif
!	write(*,*) "Like----",liker_prop,likeHoV_prop
!****************************************************************
!                  Store models for ensemble solution
!****************************************************************
	IF (ount.GT.burn_in) THEN
	sample=sample+1
	IF (mod(ount,thin)==0) THEN
!        IF (accept==1) THEN
		th = th + 1
		histo(npt)=histo(npt)+1
		call voro2qmodel(voro,npt,npt_max,d_min,d_max,beta,h,vpvs,qa,qb,alpha,rho)
  		vm = init_vmodel()
        	call vm%setvmodel(npt,alpha,beta,rho,h,vm2)

		!call dispersion(npt,beta,vpvs,h,period,Vphase,tmpHoV,peri,ndatad)
		!call dispersion(npt,beta,vpvs,h,periodHoV,tmpVphase,HoV,periHoV,ndataHoV)
		! phase 
		disp1 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(peri(1),8), &
       			& xmax   = real(peri(ndatad),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(1.0,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp1%dispersion(is_ok=is_ok)
		do i = 1,ndatad
			Vphase(i)=disp1%get_c(i)
		end do

		! group 
		disp2 = disper(&
       			& vm     = vm2, &
      			& freq_or_period = "period", &
       			& xmin   = real(periHoV(1),8), &
       			& xmax   = real(periHoV(ndataHoV),8), &
      			& dx     = real(1.0,8), &
      			& cmin   = real(0.5,8), &
       			& cmax   = real(4.3,8), &
       			& dc     = real(0.02,8), &
     	  		& disper_phase = "R", &
       			& n_mode = 0, &
       			& noise_added = 0.d0 &
       			& )
        	call disp2%dispersion(is_ok=is_ok)
		do i = 1,ndataHoV
			HoV(i)=disp2%get_u(i)
		end do

	ress=Vphase-d_obsd
	do i=1,ndatad
    		do j=1,ndatad
			!write(*,*)Moy(i,j)
        	        Moy(i,j)=Moy(i,j)+ress(i)*ress(j);

        	enddo
	enddo


	l=1
	ht=h(l)
	do i=1,disd
		d=(i-1)*prof/int(disd-1)
		if (d<ht)then
			av(i)=av(i)+beta(l)
                        if(vpvs(l)<0) vpvs(l)=vpvs_min
			avpvs(i) = avpvs(i) + vpvs(l)
			mo(i)=beta(l)
			cvs(th,i)=beta(l)
		else	
			l=l+1
			av(i)=av(i)+beta(l)
                        if(vpvs(l)<0) vpvs(l)=vpvs_min
			avpvs(i) = avpvs(i) + vpvs(l)
			mo(i)=beta(l)
			cvs(th,i)=beta(l)
			if (l<npt) then
				ht=ht+h(l)
			else
				ht=1000
			endif
		endif
		
	enddo


	ht=0
	l=1
 	ht=h(l)
 	do i=1,disd
 		d=(i-1)*prof/int(disd-1)
 		if (d<ht)then
			v=ceiling((beta(l)-beta_min+width)*&
			disv/(beta_max+2*width-beta_min))
			post(i,v)=post(i,v)+1

                        if(vpvs(l)<0) vpvs(l)=vpvs_min
			Ivpvs = ceiling((vpvs(l)-vpvs_min)*&
			disvpvs/(vpvs_max-vpvs_min))
			postvpvs(i,Ivpvs) = postvpvs(i,Ivpvs)+1
 		else	
 			l=l+1
 			v=ceiling((beta(l)-beta_min+width)*&
			disv/(beta_max+2*width-beta_min))
			post(i,v)=post(i,v)+1

                        if(vpvs(l)<0) vpvs(l)=vpvs_min
			Ivpvs = ceiling((vpvs(l)-vpvs_min)*&
			disvpvs/(vpvs_max-vpvs_min))
			postvpvs(i,Ivpvs) = postvpvs(i,Ivpvs)+1
 			if (l<npt) then
 				ht=ht+h(l)
 			else
 				ht=1000
 			endif
 		endif
 	enddo

        i=ceiling((Ar-Ar_min)*disA/(Ar_max-Ar_min))
	ML_Ar(i) = ML_Ar(i)+1
   	i=ceiling((Ad-Ad_min)*disA/(Ad_max-Ad_min))
	ML_Ad(i) = ML_Ad(i)+1
   	i=ceiling((AHoV-AHoV_min)*disA/(AHoV_max-AHoV_min))
	ML_AHoV(i) = ML_AHoV(i)+1

	!Get distribution on changepoint locations.
	ht=0
	do i=1,npt-1
		ht=ht+h(i)
		j=ceiling((ht)*disd/(prof))
		histoch(j)=histoch(j)+1
	enddo	
!        endif
	endif
	endif



!******************
! get convergence
	conv(ount)=LSr
	ncell(ount)=npt
	convAr(ount)=Ar	
	convAd(ount)=Ad	
	convAHoV(ount)=AHoV	

!*********************************************************
!Write stat
	IF (mod(ount,display).EQ.0) THEN
!		write(*,*)'processor number',ran
		write(*,*)'sample:',ount,'/',burn_in+nsample
		write(*,*)'like:',like,'liker:',liker,'SWD:',liked,'likeHoV:',likeHoV
                write(*,*)'noise',Ar,Ad,AHoV
		write(*,*)
!		write(*,*)'Acceptance rates'
!		if (ount.GT.burn_in) write(*,*)'AR_move',100*AcP(1)/PrP(1),100*AcP(2)/PrP(2)
!		if (ount.GT.burn_in) write(*,*)'AR_value',100*AcV(1)/PrV(1),100*AcV(2)/PrV(2)
	END IF

	enddo !end Markov chain

!***************************************************************************

! Collect information from all the chains and average everything

!***************************************************************************
	Moy=Moy/th
	av=av/th
	avpvs = avpvs/th
	q1=ANINT(nsample/thin*0.5*0.05)
	q2=nsample/thin-q1+1
	do i=1,disd
	do k=1,nsample/thin-1
	do l=k+1,nsample/thin
		if(cvs(k,i) .gt. cvs(l,i)) then
		sort=cvs(k,i)
		cvs(k,i)=cvs(l,i)
		cvs(l,i)=sort
		endif
	enddo
	enddo
		cvs1(i)=cvs(q1,i)
		cvs2(i)=cvs(q2,i)
	enddo
	call MPI_REDUCE(Moy,Moys,ndatad*ndatad,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(histoch,histochs,disd,MPI_Integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(ML_Ar,ML_Ars,disA,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(ML_ed,ML_eds,dise,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(ML_Ad,ML_Ads,disA,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(ML_AHoV,ML_AHoVs,disA,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(av,avs,disd,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(avpvs,avvpvs,disd,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(cvs1,cvss1,disd,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(cvs2,cvss2,disd,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(post,posts,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(postvpvs,postsvpvs,disd*disvpvs,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(ncell,ncells,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(conv,convs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(histo,histos,npt_max,MPI_Integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(convAd,convAds,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(conved,conveds,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(convAr,convArs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(conver,convers,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(convAHoV,convAHoVs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
 	call MPI_REDUCE(conveHoV,conveHoVs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(cor,cors,disv*disc,MPI_integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)


	Moys=Moys/nb
	avs=avs/nb
	avvpvs = avvpvs/nb
 	convs=convs/nb
	cvss1=cvss1/nb
	cvss2=cvss2/nb
 	convArs=convArs/nb
 	convers=convers/nb
 	convAds=convAds/nb
 	conveds=conveds/nb
 	convAHoVs=convAHoVs/nb
 	conveHoVs=conveHoVs/nb
 	ncells=ncells/nb


!***********************************************************************

!                      Write the results

!***********************************************************************
	IF (ran==0) THEN

	open(56,file='OUT/Averageg.out',status='replace')
	do i=1,disd
		d=(i-1)*prof/int(disd-1)
		write(56,*)d,avs(i),cvss1(i),cvss2(i),avvpvs(i)
	enddo
	close(56)

	open(71,file='OUT/posteriorg.out',status='replace')
	write(71,*) beta_min,beta_max,0.0,prof

	do j=1,disv
	do i=1,disd
		d=(i-1)*prof/int(disd-1)
		write(71,*) beta_min-width+&
		&(j-1)*(beta_max-beta_min+2*width)/(disv-1),d,posts(i,j)
	enddo
	enddo
	close(71)! close the file with Evidence

	open(61,file='OUT/posteriorgvpvs.out',status='replace')
	write(61,*) vpvs_min,vpvs_max,0.0,prof

	do j=1,disvpvs
	do i=1,disd
		d=(i-1)*prof/int(disd-1)
		write(61,*) vpvs_min+&
		&(j-1)*(vpvs_max-vpvs_min)/(disvpvs-1),d,postsvpvs(i,j)
	enddo
	enddo
	close(61)! close the file with Evidence


	open(72,file='OUT/data_bestrg1.out',status='replace')
	do i=1,ndatar
		xi=-time_shift+(i-1)/fs
		write(72,*)xi,best_datar1(i)
	enddo
	close(72)

	open(73,file='OUT/data_bestrg2.out',status='replace')
	do i=1,ndatar
		xi=-time_shift+(i-1)/fs
		write(73,*)xi,best_datar2(i)
	enddo
	close(73)

	open(74,file='OUT/data_bestrg3.out',status='replace')
	do i=1,ndatar
		xi=-time_shift+(i-1)/fs
		write(74,*)xi,best_datar3(i)
	enddo
	close(74)

	open(27,file='OUT/data_best_phase.out',status='replace')
	do i=1,ndatad
		write(27,*)peri(i),best_datad(i)
	enddo
	close(27)

	open(28,file='OUT/data_best_group.out',status='replace')
	do i=1,ndataHoV
		write(28,*) periHoV(i),best_dataHoV(i)
	enddo
	close(28)

	open(75,file='OUT/model_best.out',status='replace')
	do i=1,npt_best
		write(75,*) voro_best(i,1),voro_best(i,2),voro_best(i,3)
	enddo
	close(75)

	endif

!        open(104,file='OUT/mcmcHis.out',status='replace')
!        do i = 1,nb
!                write(104,*) "> cpu"
!                do j = 1,nsample+burn_in,5
!                        write(104,*) j,alllike_his(j+(i-1)*(nsample+burn_in),1),&
!			&alllike_his(j+(i-1)*(nsample+burn_in),2),&
!			&alllike_his(j+(i-1)*(nsample+burn_in),3),&
!			&alllike_his(j+(i-1)*(nsample+burn_in),4)
!                end do
!        end do
!        close(104)

	close(ran)
	CALL cpu_time(t2)
	if (ran==0) write(*,*)'time taken by the code was',t2-t1,'seconds'

	call MPI_FINALIZE(ierror)

end program Joint_inversion













!----------------------------------------------------------------------

!               FUNCTIONS USED BY THE MAIN CODE  

!---------------------------------------------------------------------


!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------


FUNCTION GASDEV(idum)

!     ..Arguments..
integer          idum
real GASDEV

!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10   v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END


!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------

FUNCTION ran3(idum)
INTEGER idum
INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803478,MZ=0,FAC=1./MBIG)
!PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
	iff=1
	mj=MSEED-iabs(idum)
	mj=mod(mj,MBIG)
	ma(55)=mj
	mk=1
	do 11 i=1,54
	ii=mod(21*i,55)
	ma(ii)=mk
	mk=mj-mk
	if(mk.lt.MZ)mk=mk+MBIG
	mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
	do 13 k=1,4
	do 12 i=1,55
	ma(i)=ma(i)-ma(1+mod(i+30,55))
	if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
	inext=0
	inextp=31
	idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
	
return
END

!the code is finished. 

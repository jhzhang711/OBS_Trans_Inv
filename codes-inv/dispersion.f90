subroutine dispersion(nlyrs,beta,vpvs,thick,t,c,u,peri,NTMAX)
! C ----------------------------------------------------------------------
! C Hrvoje Tkalcic, February 25, 2005, LLNL
! C The only input file is model.0              
! C Format - more or less self-explanatory, first line being number of layers
! C Uses the same logics as programs by Julia and Ammon but modified so that
! C it is useful for a grid search. This program can be executed from a 
! C shell script, which will form models in a grid search manner before each run.
! C The output will be 2 files, for Love and Rayleigh wave dispersion with
! C first column representing period, second - group and third - phase velocity.
! C Subroutines: lovdsp     -> DISPER80 (Saito)
! C              lovmrx     -> DISPER80 (Saito)
! C              raydsp     -> DISPER80 (Saito)
! C              raymrx     -> DISPER80 (Saito)
! C Libraries:   sac.a      -> SAC10.6
! C Sample run:  program 'dispersion' is invoked as follows
! C              dispersion > disp.out
! C  c-phase velocity; u group velocity
! C ----------------------------------------------------------------------
      PARAMETER (NLMAX=50,NCUMAX=2)
!C
      INTEGER   nlyrs, ilay, idum1, itr, iper, NTMAX
      REAL      alpha(nlyrs), beta(nlyrs),vpvs(nlyrs),rho(nlyrs),&
               thick(nlyrs), smooth(nlyrs), beta_ap(nlyrs),&
               weight(nlyrs), cmin, cmax,& 
               dc, tol, t(NTMAX+1),&
               pi2, w, ax(nlyrs), c(NTMAX), u(NTMAX), ekd, y0l(6),&
               vp(2*nlyrs), vs(2*nlyrs), z(2*nlyrs), za,&
               y0r(3), yij(15), ap(nlyrs), ae(nlyrs), peri (NTMAX),&
	       HoV(NTMAX)
      CHARACTER*20 title
!C
      EXTERNAL  lovmrx, raymrx
!C
      DATA      pi2/6.283185/, ia/0/
      cmin=0.5
      cmax=6.5
      dc=0.05 
      tol=0.0001
      itr=10
! ! ! ! !       OPEN (8,file='disp_love.output')
! ! ! ! !       OPEN (9,file='disp_rayl.output')
! CCCCCCCCCCCCCCCCCCCCCCCC
! C Reads in Earth model C
! CCCCCCCCCCCCCCCCCCCCCCCC

    za=0.
  !    OPEN (1,file='model.0',status='Old')
  !    READ (1,*) nlyrs, title

  DO ilay=1,nlyrs
  !        READ (1,*) idum1,alpha(ilay),beta(ilay),rho(ilay),thick(ilay),
   !  &         smooth(ilay),beta_ap(ilay),weight(ilay)
	
	if (vpvs(ilay) < 0.0 ) then
	alpha(ilay) = 1.50
	rho(ilay) = 1.00

	else
	alpha(ilay)=vpvs(ilay)*beta(ilay)
	rho(ilay)= 1.6612*alpha(ilay)-0.4721*alpha(ilay)*alpha(ilay)+&
		&0.0671*alpha(ilay)*alpha(ilay)*alpha(ilay)-&
		&0.0043*alpha(ilay)*alpha(ilay)*alpha(ilay)*alpha(ilay)+&
		&0.000106*alpha(ilay)*alpha(ilay)*alpha(ilay)*alpha(ilay)*alpha(ilay)
	endif

  END DO

DO ilay=1,nlyrs
	  z(2*ilay-1)=za
	  vs(2*ilay-1)=beta(ilay)
	  vp(2*ilay-1)=alpha(ilay)
	  z(2*ilay)=z(2*ilay-1)-thick(ilay)
	  vs(2*ilay)=beta(ilay)
	  vp(2*ilay)=alpha(ilay)
	  za=z(2*ilay)
       END DO

      !write (*,*)alpha,beta,rho,thick
!       !CLOSE(1)
! 
! CCCCCCCCCCCCCCCCCCCCC
! C Loop over periods C
! CCCCCCCCCCCCCCCCCCCCC
	    !t(1)=7.0	
	    DO iper=1,NTMAX
		t(iper)= peri(iper)
		!write(*,*)t(iper)
	       w=pi2/t(iper)
! C	       write (*,*)t(iper),w	
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C Call DISPER80 for Love dispersion C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ! ! ! ! 	          CALL lovdsp (lovmrx,thick,rho,beta,ax,nlyrs,w,
! ! ! ! !      &                         cmin,cmax,dc,
! ! ! ! !      &                         tol,itr,ia,c(iper),
! ! ! ! !      &                         u(iper),ekd,y0l,ier)
! ! ! ! ! 	          write (8,77)t(iper),u(iper),c(iper)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C Call DISPER80 for Rayleigh dispersion C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			!write(*,*)'------------------------'
			c(iper)=100
			u(iper)=100
                  CALL raydsp (raymrx,thick,rho,alpha,beta,ap,ae,nlyrs,&
                              w,cmin,cmax,dc,&
                              tol,itr,ia,c(iper),u(iper),&
                              ekd,y0r,yij,ier)
			HoV(iper) = -1*y0r(3)
			!write (*,*)t(iper),u(iper),c(iper)
			!write(*,*)'------------------------'
	        !!!!!!  write (9,77)t(iper),u(iper),c(iper)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C End of loop over periods C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCC
		!t(iper+1)=5.0+iper*5.0
               !t(iper+1)=t(1)+iper*5
            END DO
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C Writes predicted dispersion C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ! ! ! ! 77    format(f6.2,f7.4,f7.4)
! ! ! ! !       close(8)
! ! ! ! !       close(9)
return
END

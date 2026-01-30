SUBROUTINE VG_INV(H_PRESS,THETA,THETAS,THETAR,ALPHA,M2)  ! INVERSE VG (FROM SOIL MOISTURE CONTENT TO PRESSURE HEAD)
 IMPLICIT NONE 
 REAL (KIND = 8) :: H_PRESS
 REAL (KIND = 8) :: THETA,S_THETA,THETAS,THETAR
 REAL (KIND = 8) :: ALPHA,M1,M2
 
 M1 = 1.0-1.0/M2
 
 IF (THETA.GT.THETAS) THEN 
 	S_THETA = 1.0
 ELSE 
 	S_THETA = (THETA-THETAR)/(THETAS-THETAR)
 END IF
 
 H_PRESS = -((S_THETA**M1-1.0)**(-M2))/ALPHA

END SUBROUTINE VG_INV

SUBROUTINE gardner_inv(H_PRESS,THETA,THETAS,THETAR,ALPHA,M2)  ! Inverse Gardner model 
 IMPLICIT NONE 
 REAL (KIND = 8) :: H_PRESS
REAL (KIND = 8) :: THETA,S_THETA,THETAS,THETAR
 REAL (KIND = 8) :: ALPHA,M1,M2
 
 IF (THETA.GT.THETAS) THEN 
 	S_THETA = 1.0
 ELSE 
 	S_THETA = (THETA-THETAR)/(THETAS-THETAR)
 END IF
 
 H_PRESS = LOG(S_THETA)/M2
 
END SUBROUTINE 

subroutine vg(theta,ktheta,c_theta,h_press,alpha,thetas,thetar,k_s,m2)  ! van Genuchten
    implicit none
    real(kind=8) :: h_press
    real(kind=8) :: theta,ktheta,c_theta,s_theta
    real(kind=8) :: alpha,thetas,thetar,k_s,m1,m2
    real(kind=8) :: kd,ktheta1,ktheta2

    m1 = 1.0-1.0/m2

    if (h_press.lt.0.0) then
        s_theta = 1.0/((1.0+(alpha*abs(h_press))**m2)**m1)
    else
        s_theta = 1.0
    end if

    if (h_press.lt.0.0) then
        theta = thetar+(thetas-thetar)*s_theta
    else
        theta = thetas
    end if

    if (h_press.lt.0.0) then
        kd = (alpha*abs(h_press))**m2
      !  c_theta = (thetas-thetar)*(m2*m1*alpha)*((1.0+kd)**(-m1-1.0))*kd/(-alpha*h_press)
        c_theta = alpha*m1*(thetas-thetar)/(1.0-m1)*(s_theta**(1.0/m1))*((1.0-s_theta**(1.0/m1))**m1)
    else
        c_theta = 0.0
    end if

    if (h_press.lt.0.0) then
        ktheta = sqrt(s_theta)*(1.0-(1.0-s_theta**(1.0/m1))**m1)**2.0
    else
        ktheta = 1.0
    end if

 end subroutine

SUBROUTINE MOD_SOIL(THETA,KTHETA,C_THETA,H_PRESS,H_BOT,H_TOP,THETAS,THETAR)  ! Design for the modflow soil 
 IMPLICIT NONE 
 REAL (KIND = 8) :: THETA,KTHETA,C_THETA
 REAL (KIND = 8) :: H_PRESS,H_BOT,H_TOP,THETAS,THETAR
 
 IF (H_PRESS.GT.H_TOP) THEN 
 	KTHETA = 1.0; THETA = THETAS; C_THETA = 0.0
 ELSE 
 	C_THETA = (THETAS-THETAR)/(H_TOP-H_BOT)
 	IF (H_PRESS.LT.H_BOT) THEN 
 		KTHETA = 0.0; THETA = THETAR
 	ELSE
 		KTHETA = (H_PRESS-H_BOT)/(H_TOP-H_BOT)*1.0
 		THETA = THETAR+(THETAS-THETAR)*KTHETA
 	END IF
 END IF
 

END SUBROUTINE 


 subroutine gardner(theta,ktheta,c_theta,h_press,alpha,thetas,thetar,k_s,m2)  ! For Gardner's relationship
  implicit none 
  real (kind = 8) :: h_press
  real (kind = 8) :: theta,ktheta,c_theta
  real (kind = 8) :: alpha,thetas,thetar,k_s,m2

  if (h_press.lt.0.0) then 
		ktheta = exp(alpha*h_press)
		theta = thetar + (thetas-thetar)*exp(m2*h_press)
		c_theta = (thetas-thetar)*m2*exp(m2*h_press)
  else 
		ktheta = 1.0
		theta = thetas
		c_theta = 0.0
  end if

 end subroutine gardner

 subroutine Haverkamp(theta,ktheta,c_theta,h_press,alpha,thetas,thetar,k_s,a,m2,b)   ! This relationship is proposed by Haverkamp et al. 1977
    implicit none
    real(kind=8) :: h_press
    real(kind=8) :: theta,ktheta,c_theta
    real(kind=8) :: alpha,thetas,thetar,k_s,a,m2,b

    h_press = h_press*100.0

    if (h_press.lt.-1.0) then
        theta = alpha*(thetas-thetar)/(alpha+(log(abs(h_press)))**m2)+thetar
    else
        theta = thetas
    end if

    if (h_press.lt.-1.0) then
        ktheta = a/(a+(abs(h_press))**b)
    else
        ktheta = 1.0
    end if

    if (h_press.lt.-1.0) then
        c_theta = alpha*(thetas-thetar)/(alpha+(log(abs(h_press)))**m2)/(alpha+(log(abs(h_press)))**m2)
        c_theta = -c_theta*m2*((log(abs(h_press)))**(m2-1.0))/h_press*100.0
    else
        c_theta = 0.0
    end if

end subroutine

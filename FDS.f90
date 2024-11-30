!This code solves the 1D Euler equations with the flux difference splitting method.
!Equation numbers are from Fujii "Numerical Methods for Computational Fluid Dynamics".
!See Section 6.1 in Fujii's textbook for details.
program FDS
implicit none

integer::imax,i,j,n,nmax
real(8),allocatable::rho(:),u(:),etot(:),p(:),Q(:,:),E(:,:),E_flux(:,:),c(:),A_abs(:,:,:),x(:)
real(8),allocatable::R(:,:,:),Rinv(:,:,:),Lambda(:,:,:),H(:),rho_ave(:),u_ave(:),H_ave(:),c_ave(:)
real(8)::gamma,b1,b2,dx,dt,t,CFL,c_temp,etot_temp,rho_temp,u_temp,p_temp

!Initial condition
!-------------------------------------------------
imax=500
nmax=1000
CFL=0.2
gamma=1.4

allocate(rho   (0:imax+1))
allocate(u     (0:imax+1))
allocate(etot  (0:imax+1))
allocate(p     (0:imax+1))
allocate(c     (0:imax+1))
allocate(H     (0:imax+1))
allocate(Q     (0:imax+1,3))
allocate(E     (0:imax+1,3))
allocate(E_flux(0:imax,3))
allocate(R     (0:imax,3,3))
allocate(Rinv  (0:imax,3,3))
allocate(Lambda(0:imax,3,3))
allocate(A_abs (0:imax,3,3))
allocate(rho_ave(0:imax))
allocate(u_ave  (0:imax))
allocate(H_ave  (0:imax))
allocate(c_ave  (0:imax))
allocate(x      (1:imax))

!Test: Shock tube
do i=0,imax/2
	rho(i) = 1d0
	p  (i) = 1d0
	u  (i) = 0d0
end do

do i=imax/2+1,imax+1
	rho(i) = 0.125
	p  (i) = 0.1
	u  (i) = 0d0
end do
!-------------------------------------------------

!Grid generation
!-------------------------------------------------
dx=10d0/imax
do i=1,imax
	x(i) = dx*i
end do
!-------------------------------------------------

!Calculating conserved quantities
!-------------------------------------------------
do i=0,imax+1
	rho_temp = rho(i)
	p_temp   = p(i)
	u_temp   = u(i)
	call EoS(gamma,rho_temp,u_temp,p_temp,etot_temp,c_temp)
	etot(i) = etot_temp
	c(i) = c_temp
		
	H(i) = (etot(i)+p(i))/rho(i)
	!Eq. (6.2)
	!-------------------------------------------------
	Q(i,1) = rho(i)
	Q(i,2) = rho(i)*u(i)
	Q(i,3) = etot(i)
	E(i,1) = rho(i)*u(i)
	E(i,2) = p(i)+rho(i)*u(i)**2d0
	E(i,3) = (etot(i)+p(i))*u(i)
	!-------------------------------------------------
!		write(*,*) c_temp,p(i),gamma
end do
!-------------------------------------------------

t=0d0
do n=1,nmax

!write(*,*) Q,E
	dt = CFL*dx/maxval(c(1:imax))
	t=t+dt

	Lambda=0d0
	do i=0,imax
		!Roe's average
		!-------------------------------------------------
		rho_ave(i) = sqrt(rho(i)*rho(i+1)) !Eq. (6.20)
		u_ave  (i) = (sqrt(rho(i))*u(i)+sqrt(rho(i+1))*u(i+1))/(sqrt(rho(i))+sqrt(rho(i+1))) !Eq. (6.21)
		H_ave  (i) = (sqrt(rho(i))*H(i)+sqrt(rho(i+1))*H(i+1))/(sqrt(rho(i))+sqrt(rho(i+1))) !Eq. (6.22)
		c_ave  (i) = sqrt((gamma-1d0)*(H_ave(i)-0.5*u_ave(i)**2d0)) !Eq. (6.23)
		!-------------------------------------------------

		b1 = 0.5*u_ave(i)**2d0*(gamma-1d0)/c_ave(i)**2d0
		b2 = (gamma-1d0)/c_ave(i)**2d0
		
		!Eq. (5.76)
		!-------------------------------------------------
		R(i,1,1) = 1d0
		R(i,1,2) = 1d0
		R(i,1,3) = 1d0
		R(i,2,1) = u_ave(i)-c_ave(i)
		R(i,2,2) = u_ave(i)
		R(i,2,3) = u_ave(i)+c_ave(i)
		R(i,3,1) = H_ave(i)-u_ave(i)*c_ave(i)
		R(i,3,2) = 0.5*u_ave(i)**2d0
		R(i,3,3) = H_ave(i)+u_ave(i)*c_ave(i)
		!-------------------------------------------------

		!Eq. (5.77)
		!-------------------------------------------------
		Rinv(i,1,1) = 0.5*(b1+u_ave(i)/c_ave(i))
		Rinv(i,1,2) = -0.5*(c_ave(i)**(-1d0)+b2*u_ave(i))
		Rinv(i,1,3) = 0.5*b2
		Rinv(i,2,1) = 1d0-b1
		Rinv(i,2,2) = b2*u_ave(i)
		Rinv(i,2,3) = -b2
		Rinv(i,3,1) = 0.5*(b1-u_ave(i)/c_ave(i))
		Rinv(i,3,2) = 0.5*(c_ave(i)**(-1d0)-b2*u_ave(i))
		Rinv(i,3,3) = 0.5*b2
		!-------------------------------------------------

		Lambda(i,1,1) = abs(u_ave(i)-c_ave(i))
		Lambda(i,2,2) = abs(u_ave(i))
		Lambda(i,3,3) = abs(u_ave(i)+c_ave(i))
		A_abs(i,:,:) = matmul(R(i,:,:),matmul(Lambda(i,:,:),Rinv(i,:,:)))
!		write(*,*) i,rho_ave(i),u_ave(i),H_ave(i),c_ave(i)
	end do

	!Numerical flux
	!-------------------------------------------------
	do i=0,imax
		E_flux(i,:) = 0.5*(E(i+1,:)+E(i,:)-matmul(A_abs(i,:,:),(Q(i+1,:)-Q(i,:)))) !Eq. (6.25)
	end do
	!-------------------------------------------------
	
!	write(*,*) Rinv
	
	do i=1,imax
		Q(i,:) = Q(i,:)-dt/dx*(E_flux(i,:)-E_flux(i-1,:))
		!Eq. (6.5)
		!-------------------------------------------------
		E(i,1) = Q(i,2)
		E(i,2) = (gamma-1d0)*Q(i,3)+0.5*(3d0-gamma)*Q(i,2)**2d0/Q(i,1)
		E(i,3) = gamma*Q(i,3)*Q(i,2)/Q(i,1)-0.5*(gamma-1d0)*Q(i,2)**3d0/Q(i,1)**2d0
		!-------------------------------------------------
	end do
	
	!Boundary condition
	!-------------------------------------------------
	Q(0     ,:) = Q(1   ,:)
	Q(imax+1,:) = Q(imax,:)
	E(0     ,:) = E(1   ,:)
	E(imax+1,:) = E(imax,:)
	!-------------------------------------------------
	
	do i=0,imax+1
		rho(i) = Q(i,1)
		u  (i) = Q(i,2)/rho(i)
		p  (i) = E(i,2)-rho(i)*u(i)**2d0
	end do
	
	write(*,*) "#Time=",t
	do i=1,imax
		write(*,*) x(i),rho(i),u(i),p(i)
	end do
	write(*,*)
	write(*,*)
	
end do

stop

contains

!Equation of state for ideal gas
subroutine EoS(gamma,rho,u,p,etot,c)
	real(8),intent(in)::gamma,rho,u,p
	real(8),intent(out)::etot,c
	etot  = p/(gamma-1d0)+0.5*rho*u**2d0
	c     = sqrt(gamma*p/rho)
end subroutine

end program
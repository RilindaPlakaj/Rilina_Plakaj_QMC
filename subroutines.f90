!**********************Fisrt subroutine*************************

subroutine pdmc(a,N,nn,Z,Nuc,dt,nmax,energy,accep,tau,E_ref)

  implicit none
  real*8, intent(in)             :: a,dt,tau,Nuc(3*nn)
  integer, intent(in)            :: N,nn,Z
  integer*8, intent(in)          :: nmax
  real*8, intent(out)            :: energy, accep
  real*8, intent(in)             :: E_ref
  integer*8                      :: istep
  integer*8                      :: n_accep
  real*8, external               :: e_loc, psi
  real*8                         :: sq_dt, chi(3*N), B2_old, B2_new,prod, u, q
  real*8                         :: psi_old, psi_new,argexpo
  real*8                         :: R_old(3*N), R_new(3*N)
  real*8                         :: e, w, norm, tau_current, B_old(3*N),B_new(3*N)
  integer                        :: i,j,k



  sq_dt = dsqrt(dt)   !dt is the variance of the transition prob. gaussi initialization

  energy = 0.d0
  n_accep = 0.d0
  norm = 0.d0
  w = 1.d0
  tau_current = 0.d0
  call random_gauss(R_old,3*N)

 call drift(a,R_old,N,nn,Nuc,B_old)

  do i = 1,3*N
    B2_old = B2_old + B_old(i)*B_old(i)   !here we calculate the square of the drift

 end do

  psi_old = psi(a,R_old,N,nn,Nuc)

 do istep = 1,nmax
    e =  e_loc(a,R_old,N,nn,Nuc,Z)   !the local energy in the old position with
    w = w * dexp(-dt*(e - E_ref))   !it's weight

    norm = norm + w
    energy = energy + w*e

    tau_current = tau_current + dt

    ! Reset when tau is reached
    if (tau_current > tau) then
      w = 1.d0
      tau_current = 0.d0
    end if

    call random_gauss(chi,3*N)
    R_new(:) = R_old(:) + dt*B_old(:) + chi(:)*sq_dt

    call drift(a,R_new,N,nn,Nuc,B_new)

    do i = 1,3*N
      B2_new = B2_new + B_new(i)*B_new(i)   !here we calculate the square of the drift
    end do

    psi_new = psi(a,R_new,N,nn,Nuc)

    !Metropolis, here we calculate the ratio of the acceptances as
    ! T(r_old|r_new)P(r_new)/T(r_new|r_old)P(r_old)
    do i = 1,3*N
       prod = prod + (B_new(i) + B_old(i))*(R_new(i) - R_old(i))
    end do

    argexpo = 0.5d0 * (B2_new - B2_old)*dt + prod

    q = psi_new / psi_old
    q = dexp(-argexpo) * q*q

    call random_number(u)

        if (u <= q) then                !the new move is accepted with a
                                         !probability that match with
                                         !out
                                         !probability density

            n_accep = n_accep + 1

            R_old(:) = R_new(:)
            B_old(:) = B_new(:)
            B2_old = B2_new
            psi_old = psi_new

    end if
  end do

energy = energy / norm

  accep = dble(n_accep) / dble(nmax)

end subroutine pdmc


!*****************************2nd subroutine**************************************

subroutine random_gauss(z,k)
  implicit none
  integer, intent(in)            :: k
  real*8, intent(out)            :: z(k)
  real*8                         :: u(k+1)
  real*8, parameter              :: two_pi = 2.d0*dacos(-1.d0)
  integer                        :: i

  call random_number(u)

  if (iand(k,1) == 0) then
    ! k is even
    do i=1,k,2
      z(i) = dsqrt(-2.d0*dlog(u(i)))
      z(i+1) = z(i) * dsin( two_pi*u(i+1))
      z(i) = z(i) * dcos( two_pi*u(i+1))
    end do

  else
    !k is odd
    do i=1,k-1,2
      z(i) = dsqrt(-2.d0*dlog(u(i)))
      z(i+1) = z(i) * dsin( two_pi*u(i+1))
      z(i) = z(i) * dcos( two_pi*u(k+1))
    end do

    z(k) = dsqrt( -2.d0*dlog(u(k)))
    z(k) = z(k) * dcos( two_pi*u(k+1))

  end if

end subroutine random_gauss




!****************************3d Subroutine****************************************


subroutine drift(a,Rtot,N,nn,Nuc,Bi)
  implicit none !this subroutine evaluate the drift vector as
                ! (gradient of psi)/psi for everi electron
  real*8, intent(in)             :: a
  integer, intent(in)            :: N,nn
  real*8, intent(in)             :: Rtot(3*N), Nuc(3*nn)
  real*8, intent(out)            :: Bi(3*N)                    !the drift vector
  real*8                         :: r(3),b(3), Numerator(3)
  integer                        :: i,j,k
  real*8                         :: distance_e_nuc,Den

  do i = 1,3*N,3
    do j = 1,3*nn,3
      r(1) = Rtot(i) - Nuc(j)
      r(2) = Rtot(i+1) - Nuc(j+1)
      r(3) = Rtot(i+2) - Nuc(j+2)
      distance_e_nuc = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
      b(:) = r(:)*(-a * (1/distance_e_nuc))*dexp(-a * distance_e_nuc)
       Den = Den + dexp(-a * distance_e_nuc)
       Numerator(:) = Numerator(:) + b(:)
       r(:) = 0.d0
       distance_e_nuc = 0.d0
       b(:) = 0.d0
    end do
    Numerator(:) = Numerator(:) / Den
    Bi(i) = Numerator(1)
    Bi(i+1) = Numerator(2)
    Bi(1+2) = Numerator(3)
    Numerator(:) = 0.d0
    Den = 0.d0
  end do

end subroutine drift


!*********************************4th subroutine****************************************

subroutine ave_error(x,n,ave,err)
  implicit none ! the subroutine calculate the average and the statistical err.
                ! for the elements of a given array x(n)
  integer, intent(in)            :: n
  real*8, intent(in)             :: x(n)
  real*8, intent(out)            :: ave, err
  real*8                         :: variance

  if (n < 1) then
    stop 'n<1 in ave_error' !it has no sense having an array of negative dimen.

  else if (n == 1) then
    ave = x(1)
    err = 0.d0

  else
    ave = sum(x(:)) / dble(n)
    variance = sum((x(:) - ave)**2)/ dble(n-1)
    err = dsqrt(variance/dble(n))

  end if

end subroutine ave_error





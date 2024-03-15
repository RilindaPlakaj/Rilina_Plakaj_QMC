!*********************Function1***************************************
real*8 function potential(Rtot,N,nn,Z,Nuc)

  implicit none   !this function evaluate the potential at R(3N)
  real*8, intent(in)            :: Rtot(3*N) 
  real*8, intent(in)            :: Nuc(3*nn)
  integer, intent(in)           :: Z      !the nuc. charge
  integer, intent(in)           :: N, nn  !number of electrons and nucleus
  real*8                        :: distance_ee, distance_nn
  real*8                        :: distance_e_nuc
  real*8                        :: r(3)
  integer                       :: i,k,j,p,f,pn,fn

 potential = 0.d0

 do i = 1,3*N,3
    do j = 1,3*nn,3
      r(1) = Rtot(i) - Nuc(j)
      r(2) = Rtot(i+1) - Nuc(j+1)
      r(3) = Rtot(i+2) - Nuc(j+2)
      distance_e_nuc = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
        if (distance_e_nuc > 0.d0) then
          potential = potential - (Z / distance_e_nuc)
          r(:) = 0.d0
          distance_e_nuc = 0.d0
        else
           stop 'potential at R diverges'
        end if
    end do
  end do

   if (N.eq.2) then
    r(1) = Rtot(1) - Rtot(4)
    r(2) = Rtot(2) - Rtot(5)
    r(3) = Rtot(3) - Rtot(6)
     distance_ee = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
     if (distance_ee > 0) then
              potential = potential + 1 / distance_ee
            r(:) = 0.d0
            distance_ee = 0.d0
            else
              stop 'potential al R ee  diverges'
            end if
          end if
 

  if (nn.gt.1) then
  do pn=1,3*nn,3
     do fn=1,3*nn,3
         if (pn.lt.fn) then
              r(1) = Rtot(pn) - Rtot(fn)
              r(2) = Rtot(pn+1) - Rtot(fn+1)
              r(3) = Rtot(pn+2) - Rtot(fn+2)
              distance_nn = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
            if (distance_nn > 0) then
              potential = potential + Z*Z / distance_nn
            r(:) = 0.d0
            distance_nn = 0.d0
            else
              stop 'potential al R nn  diverges'
            end if
       else 
         potential = potential + 0
          end if
         end do
        end do
   end if    
 
end function potential


!***************************Function2****************************************

real*8 function psi(a,Rtot,N,nn,Nuc)

  implicit none  !this function compute the wave func. at R
  real*8, intent(in)             :: a
  integer, intent(in)            :: N,nn
  real*8, intent(in)             :: Rtot(3*N), Nuc(3*nn)
  integer                        :: i,j
  real*8                         :: r(3),phi,Norm,distance_e_nuc
  
  psi = 1.d0
  phi = 0.d0
  r(:) = 0.d0
  distance_e_nuc = 0.d0
  Norm = 0.d0

   do i = 1,3*N,3
    do j = 1,3*nn,3
      r(1) = Rtot(i) - Nuc(j)
      r(2) = Rtot(i+1) - Nuc(j+1)
      r(3) = Rtot(i+2) - Nuc(j+2)
      distance_e_nuc = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
        phi = phi + dexp(-a * distance_e_nuc)
          r(:) = 0.d0
          distance_e_nuc = 0.d0
    end do
    psi = psi * phi
  end do

  Norm = (1/(N**(1/2)))**N
  psi = psi * Norm

 end function psi


!*****************************Function3********************************


real*8 function kinetic(a,Rtot,N,nn,Nuc)

  implicit none  !this function computes the local kin. energy at R
  real*8, intent(in)            :: a
  integer, intent(in)           :: N,nn
  real*8, intent(in)            :: Rtot(3*N), Nuc(3*nn)
  real*8                        :: distance_e_nuc
  real*8                        :: r(3),Num,Den,Arg
  integer                       :: i,j

  kinetic = 0
  Num = 0
  Den = 0
  Arg = 0

  do i = 1,3*N,3
    do j = 1,3*nn,3
      r(1) = Rtot(i) - Nuc(j)
      r(2) = Rtot(i+1) - Nuc(j+1)
      r(3) = Rtot(i+2) - Nuc(j+2)
      distance_e_nuc = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
      if (distance_e_nuc > 0.d0) then
       Num = Num + (a**2 - (2*a / distance_e_nuc))*dexp(-a * distance_e_nuc)
       Den = Den + dexp(-a * distance_e_nuc)
       r(:) = 0.d0
       distance_e_nuc = 0.d0
        else
          stop 'kinetic energy diverges'
      end if
    end do
    Arg = Num / Den
    kinetic = kinetic + Arg
    Num=0.d0
    Den=0.d0
    Arg=0.d0
  end do

    kinetic = -0.5d0 * kinetic

end function kinetic

!*******************************Function4************************************

real*8 function e_loc(a,Rtot,N,nn,Nuc,Z)

  implicit none !this function compute the local energy in R

  real*8, intent(in)            :: a
  integer, intent(in)           :: N,nn,Z
  real*8, intent(in)            :: Rtot(3*N),Nuc(3*nn)
  real*8, external              :: kinetic   ! we use the kinetic local energy
  real*8, external              :: potential ! and the potential calculated
                                             ! before
  e_loc = kinetic(a,Rtot,N,nn,Nuc)+potential(Rtot,N,nn,Z,Nuc)

end function e_loc

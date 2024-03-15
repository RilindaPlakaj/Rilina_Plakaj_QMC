program main_qmc
  implicit none

  integer             :: molecule
  real*8              :: a
  integer             :: N,i
  integer             :: nn
  integer             :: Z
  integer*8           :: nmax     !number of MC iteration for run
  integer             :: nruns    !number of MC runs
  real*8              :: dt
  real*8              :: tau
  real*8              :: E_ref
  integer             :: irun
  real*8, allocatable                :: X(:), accep(:)
  real*8                             :: ave, err
  real*8, allocatable                :: Nuc(:)


  write(*,*) '*********************************************************************'
  write(*,*) 'This program computes the ground state energy for:'
  write(*,*) '   H       H2       H2+       H3+       He'
  write(*,*) '*********************************************************************'
  write(*,*) 'For system:                    press:'
  write(*,*) '  H                              1'
  write(*,*) '  H2                             2'
  write(*,*) '  H2+                            3'
  write(*,*) '  H3+                            4'
  write(*,*) '  He                             5'
  write(*,*) '*********************************************************************'
  read(*,*) molecule

!  if (molecule.eq.1)then
! open(unit = 1, file = 'H_atom.dat', status = 'old', action = 'read')
 ! else if(molecule.eq.5) then
! open(unit = 1, file = 'He_atom.dat', status = 'old', action = 'read')
 ! else if(molecule.eq.3) then
 !open(unit = 1, file = 'H2+_ion.dat', status = 'old', action = 'read')
 ! else if(molecule.eq.2) then
 !open(unit = 1, file = 'H2_molecule.dat', status = 'old', action = 'read')
 ! elseif(molecule.eq.4) then
! open(unit = 1, file = 'H3+_ion.dat', status = 'old', action = 'read')
 ! end if

select case (molecule)
case(1)
open(unit=1, file='H_atom.dat', status = 'old', action = 'read')
case(2)
open(unit = 1, file = 'H2_molecule.dat', status = 'old', action = 'read')
case(3)
open(unit = 1, file = 'H2+_ion.dat', status = 'old', action = 'read')
case(4)
open(unit = 1, file = 'H3+_ion.dat', status = 'old', action = 'read')
case(5)
open(unit = 1, file = 'He_atom.dat', status = 'old', action = 'read')
end select


  read(1,*) a
  read(1,*) N
  read(1,*) nn
  read(1,*) Z
  read(1,*) nmax
  read(1,*) nruns
  read(1,*) dt
  read(1,*) tau
  read(1,*) E_ref
  allocate(Nuc(3*nn))
  do i=1,3*nn
  read(1,*) Nuc(i)
  enddo
 close(1)

  allocate(X(nruns))
  allocate(accep(nruns))

  do irun = 1, nruns
    call pdmc(a,N,nn,Z,Nuc,dt, nmax,X(irun), accep(irun), tau, E_ref)
  end do

  call ave_error(X, nruns, ave, err)

write(*,*) '************************** Result ***********************************'
write(*,*) 'The ground state energy is:'

  print *, 'E = ', ave, '+/-' , err

  call ave_error(accep, nruns, ave, err)

write(*,*) '                                     '
write(*,*) 'With an average error of:'
 print *, 'A = ', ave, '+/-' , err




end program main_qmc

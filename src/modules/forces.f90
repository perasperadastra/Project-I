!==============================================================================!
!                             MODULE FORCES
! This module contains all the subroutines regarding the forces that drive the
! trajectories of the molecular dynamics.
!==============================================================================!
! The module contains:
!                        -> force()
!            Calculates the forces between neighboor particles (using a cutoff)
!             and returns the forces, the potential energy, the pressure
!                                         Input:
!                                                        - natom (number of atoms) (in): integer scalar
!                                                        - r(coordinater)(inout): double precision array
!                                                        - boxlength (in) : double precision scalar
!                                                        - rc(cutoff radius)(in):  double precision scalar
!                                                        - F (forces) (inout): double precision array
!                                                        - deltag (in): double precision scalar
!                                                        - gr (g(r)) (out): double precision array
!                                         Output:
!                                                        - epot (E potential) (out): double precision scalar
!                                                        - press (pressure) (out): double precision scalar
!                                         Depencency:
!                                                        - lj() : Tool in thi module
!
!                        ->lj()
!                                         Input variables:
!                                                        -r(coordinater)(inout): double precision array
!                                                        -boxlength (in) : double precision scalar
!                                              -rc(cutoff radius)(in):  double precision scalar
!                                                        -ii (in): integer scalar
!                                                        -jj (in): integer scalar
!                                                        -F (forces) (inout): double precision array
!                                         Output:
!                                                        -pot (potential energy) (out): double precision array
!                                                         -piter (potential pressure) (out):
!                                             -d  (distance) (in)
!                                         Depencency:
!                                             -pbc(): Tool in boundary module
!
! Dependency:
!                        -> boundary.f90 module
!                        -> constants.h file
!==============================================================================!

module forces
    use boundary
    use mpi
    implicit none

contains
!==============================================================================!
!                                                       FORCES SUBROUTINE
!==============================================================================!
! Input:
!                - natom (number of atoms) (in): integer scalar
!                                        Total numbe of atoms in the system.
!                - r(coordinater)(inout): double precision array
!                                        Array which contains the positions of the atoms in the lattice
!
!                - boxlength (in) : double precision scalar
!
!                - rc(cutoff radius)(in):  double precision scalar
!                                        Cutoff radius from which interactions are neglected
!
!                - F (forces) (inout): double precision array
!                                        Array of the interacting forces on each atom of the simulation box.
!
!                - deltag (in): double precision scalar
!                                        Width of the bin of the rdf
!
!                - gr (g(r)) (out): double precision array
!                                        Radial pair distribution of particles
! Output:
!                - epot (E potential) (out): double precision scalar
!                                        Potential energy term  calculated at each time step
!
!                - press (pressure) (out): double precision scalar
!                                        Potential pressure term  calculated at each time step
! Depencency:
!                - lj(): Tool in thi module
!
!==============================================================================!
    subroutine force(natoms, r, boxlength, rc, F, particle_range)
        integer, intent(in)::natoms, particle_range(2)
        double precision, allocatable, intent(in) :: r(:, :)
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, intent(in) :: boxlength, rc
        integer :: ii, is, js, kk, M

        F = 0.d0
        do is = particle_range(1), particle_range(2)
            do js = 1, natoms
                if (js .ne. is) then
                    call lj(r, boxlength, rc, is, js, F, particle_range)
                    ! calling function that computes the Lennard-Jones interaction between
                    ! pair of particles i and j
                end if
            end do
        end do
    end subroutine force

!==============================================================================!
!                        POTENTIAL SUBROUTINE
!==============================================================================!
    subroutine potential(natoms, r, boxlength, rc, Upot, pressp, gr, deltag, particle_range, taskid)
        include "../declaration_variables/parallel_variables.h"
        include "constants.h"
        integer, intent(in)::natoms, particle_range(2)
        double precision, allocatable, intent(in) :: r(:, :)
        double precision, allocatable, intent(inout) :: gr(:)
        double precision, intent(in) :: boxlength, rc, deltag
        double precision, intent(out) :: Upot, pressp
        double precision :: vol, rho, factp, cutoff_press, cutoff_pot
        double precision :: dx, dy, dz, d, dU
        double precision :: pot, piter, epot, press
        integer :: ig
        integer :: is, js, kk, M

        vol = boxlength**3.; rho = dble(natoms)/vol
        factp = (16.d0/3.d0)*PI*(rho**2)

        cutoff_pot = 4.d0*(1.d0/(rc**12.d0) - 1.d0/(rc**6.d0))

        press = 0.d0; epot = 0.d0; Upot = 0.d0; pressp = 0.d0

        do is = particle_range(1), particle_range(2)
            do js = 1, natoms
                if (js .ne. is) then
                    pot = 0.d0; piter = 0.d0
                    dx = r(is, 1) - r(js, 1); dy = r(is, 2) - r(js, 2); dz = r(is, 3) - r(js, 3)
                    ! Apply the boundary conditions to the particles distance
                    call pbc(dx, boxlength, 0.d0)
                    call pbc(dy, boxlength, 0.d0)
                    call pbc(dz, boxlength, 0.d0)
                    d = (dx**2.+dy**2.+dz**2.)**0.5 ! Distance between particles i and j

                    if (d .lt. rc) then
                        dU = (48.d0/(d**14.0) - 24.d0/(d**8.0))
                        pot = pot + 4.d0*(1.d0/(d**12.d0) - 1.d0/(d**6.d0)) - cutoff_pot
                        piter = piter + dU*dx; piter = piter + dU*dy; piter = piter + dU*dz
                    end if
                    press = press + piter; epot = epot + pot

                    if (d .lt. rc) then ! computation of the radial distribution function
                        ! adding the each pair of interaction into the corresponding bin
                        ig = int(d/deltag)
                        gr(ig) = gr(ig) + 1 ! in this case the addition to the histogram is 1
                        ! since the loop is for all particles instead of interactions.
                    end if
                end if
            end do
        end do

        call MPI_REDUCE(epot, Upot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        call MPI_REDUCE(press, pressp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        if (taskid .eq. 0) then
            Upot = Upot/2.d0; pressp = pressp/2
            cutoff_press = factp*((2.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
            pressp = pressp/(3.d0*vol) + cutoff_press
        end if

    end subroutine potential

!==============================================================================!
!                  LENNARD-JONES INTERACTION COMPUTATION
!==============================================================================!
! Input variables:
!                -r(coordinater)(inout): double precision array
!                              Array which contains the positions of the atoms in the lattice
!                -boxlength (in) : double precision scalar
!                                        Length of the box of simulation
!                -rc(cutoff radius)(in):  double precision scalar
!                                        Cutoff radius from which interactions are neglected
!                -ii (in): integer scalar
!                                        Index of the first particle which is looped around
!                -jj (in): integer scalar
!                                        Index of pair particle of i at the actual interaction computation
!                -F (forces) (inout): double precision array
!                                        Array of the interacting forces on each atom of the simulation box.
!
! Output:
!                -pot (potential energy) (out): double precision array
!                                        Potential energy interaction between i and j elements
!                -piter (potential pressure) (out):
!                                        Potential pressure term from interaction between particles i and j
!                -d  (distance) (in)
!                                        Distance between particles i and j
!
! Depencency:
!                -pbc(): Tool in boundary module
!==============================================================================!
    subroutine lj(r, boxlength, rc, ii, jj, F, particle_range)
        integer, intent(in) :: particle_range(2)
        double precision, allocatable, intent(in) :: r(:, :)
        double precision, allocatable, intent(inout) :: F(:, :)
        double precision, intent(in) :: boxlength, rc
        integer, intent(in) :: ii, jj
        double precision :: dx, dy, dz, d, dU
        integer :: iv

        iv = ii - particle_range(1) + 1
        dx = r(ii, 1) - r(jj, 1); dy = r(ii, 2) - r(jj, 2); dz = r(ii, 3) - r(jj, 3)
        ! Apply the boundary conditions to the particles distance
        call pbc(dx, boxlength, 0.d0)
        call pbc(dy, boxlength, 0.d0)
        call pbc(dz, boxlength, 0.d0)
        d = (dx**2.+dy**2.+dz**2.)**0.5 ! Distance between particles i and j

        if (d .lt. rc) then
            dU = (48.d0/(d**14.0) - 24.d0/(d**8.0))
            F(iv, 1) = F(iv, 1) + dU*dx
            F(iv, 2) = F(iv, 2) + dU*dy
            F(iv, 3) = F(iv, 3) + dU*dz
        end if
    end subroutine

end module forces

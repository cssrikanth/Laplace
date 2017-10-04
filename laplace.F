!
!
!  Solves a nonlinear system in parallel with a user-defined
!  Newton method that uses KSP to solve the linearized Newton sytems.  This solver
!  is a very simplistic inexact Newton method.  The intent of this code is to
!  demonstrate the repeated solution of linear sytems with the same nonzero pattern.

!
!  The command line options include:
!  -par <parameter>, where <parameter> indicates the problem's nonlinearity
!     problem SFI:  <parameter> = Bratu parameter (0 <= par <= 6.81)
!  -mx <xg>, where <xg> = number of grid points in the x-direction
!  -my <yg>, where <yg> = number of grid points in the y-direction
!  -Nx <npx>, where <npx> = number of processors in the x-direction
!  -Ny <npy>, where <npy> = number of processors in the y-direction
!  -mf use matrix free for matrix vector product
!
!/*T
!   Concepts: KSP^writing a user-defined nonlinear solver
!   Concepts: DMDA^using distributed arrays
!   Processors: n
!T*/
!  ------------------------------------------------------------------------
!
!    Solid Fuel Ignition (SFI) problem.  This problem is modeled by
!    the partial differential equation
!
!            -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1,
!
!    with boundary conditions
!
!             u = 0  for  x = 0, x = 1, y = 0, y = 1.
!
!    A finite difference approximation with the usual 5-point stencil
!    is used to discretize the boundary value problem to obtain a nonlinear
!    system of equations.
!
!    The SNES version of this problem is:  snes/examples/tutorials/ex5f.F
!
!  -------------------------------------------------------------------------

      program main
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                    Include files
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     petscsys.h       - base PETSc routines   petscvec.h - vectors
!     petscmat.h - matrices
!     petscis.h     - index sets            petscksp.h - Krylov subspace methods
!     petscviewer.h - viewers               petscpc.h  - preconditioners

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

      MPI_Comm comm
      Vec      X,Y,F,localX,Fold,Fnew
      Mat      J,B
      DM       da
      KSP      ksp

      PetscInt  Nx,Ny,N,mx,my,ifive,ithree
      PetscBool  flg,nooutput,usemf
      common   /mycommon/ mx,my,B,localX,da
!
!
!      This is the routine to use for matrix-free approach
!


!     --------------- Data to define nonlinear solver --------------
      PetscReal   rtol,ttol,error
      PetscReal   fnorm,ynorm,xnorm,fnorm1
      PetscInt            max_nonlin_its,one
      PetscInt            lin_its,ncount
      PetscInt           i,m
      PetscScalar        mone
      PetscErrorCode ierr

      mone           = -1.0
      rtol           = 1.e-8
      max_nonlin_its = 100
      one            = 1
      ifive          = 5
      ithree         = 3
      ncount	     = 1
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      comm = PETSC_COMM_WORLD

!  Initialize problem parameters

!
      mx = 4
      my = 4
      call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,    &
     &                        '-mx',mx,flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,    &
     &                        '-my',my,flg,ierr)
      N = mx*my

      nooutput = .false.
      call PetscOptionsHasName(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,    &
     &                         '-no_output',nooutput,ierr)

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create linear solver context
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call KSPCreate(comm,ksp,ierr)

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create vector data structures
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!
!  Create distributed array (DMDA) to manage parallel grid and vectors
!
      Nx = PETSC_DECIDE
      Ny = PETSC_DECIDE
      call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,   &
     &                        '-Nx',Nx,flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,   &
     &                         '-Ny',Ny,flg,ierr)
      call DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,         &
     &     DMDA_STENCIL_STAR,mx,my,Nx,Ny,one,one,                        &
     &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)

!
!  Extract global and local vectors from DMDA then duplicate for remaining
!  vectors that are the same types
!
       call DMCreateGlobalVector(da,X,ierr)
       call DMCreateLocalVector(da,localX,ierr)
       call VecDuplicate(X,F,ierr)
       call VecDuplicate(X,Y,ierr)
       call VecDuplicate(X,Fold,ierr)
       call VecDuplicate(X,Fnew,ierr)




!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Customize linear solver set runtime options
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     Set runtime options (e.g., -ksp_monitor -ksp_rtol <rtol> -ksp_type <type>)
!

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Evaluate initial guess
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       call FormInitialGuess(X,ierr)
       call ComputeFunction(X,F,ierr)
       call VecNorm(F,NORM_2,fnorm,ierr)
       ttol = fnorm*rtol
       if (.not. nooutput) then
         print*, 'Initial function norm ',fnorm
       endif
!       call VecView(F,PETSC_VIEWER_STDOUT_WORLD,ierr)

!

       do 10 


         
!  Compute updated iterate
!         call VecNorm(Y,NORM_2,ynorm,ierr)
!         call VecAYPX(Y,mone,X,ierr)
!         call VecCopy(Y,X,ierr)
!         call VecNorm(X,NORM_2,xnorm,ierr)
         


         call VecCopy(F,Fold,ierr)
!         call VecView(Fold,PETSC_VIEWER_STDOUT_WORLD,ierr)
         call ComputeFunction(Fold,F,ierr)
!         call VecView(F,PETSC_VIEWER_STDOUT_WORLD,ierr)
!         call VecAXPY(Fold,-1,F,ierr)
!         call VecView(Fold,PETSC_VIEWER_STDOUT_WORLD,ierr)
         call VecNorm(F,NORM_1,fnorm1,ierr)
!         call VecNorm(Fnew,NORM_1,fnorm,ierr)
         error = fnorm1
         
!         call VecView(Fold,PETSC_VIEWER_STDOUT_WORLD,ierr)
         print *,error

!         if (.not. nooutput) then
!           print*, 'Iteration ',ncount,' function norm',fnorm
!         endif
!         ncount=ncount+1
!         print *,ncount,fnorm
!  Test for convergence

!       if (fnorm .le. ttol) then
!         if (.not. nooutput) then
!           print*,'Converged: function norm ',fnorm,' tolerance ',ttol
!         endif
!         exit
!       endif
 10   continue
 20   continue

      write(6,100) i+1
 100  format('Number of SNES iterations =',I2)
!       call VecView(X,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Free work space.  All PETSc objects should be destroyed when they
!     are no longer needed.
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       call VecDestroy(localX,ierr)
       call VecDestroy(X,ierr)
       call VecDestroy(Y,ierr)
       call VecDestroy(F,ierr)
       call VecDestroy(Fold,ierr)
       call VecDestroy(Fnew,ierr)
       call KSPDestroy(ksp,ierr)
       call DMDestroy(da,ierr)
       call PetscFinalize(ierr)
       end

! -------------------------------------------------------------------
!
!   FormInitialGuess - Forms initial approximation.
!
!   Input Parameters:
!   X - vector
!
!   Output Parameter:
!   X - vector
!
      subroutine FormInitialGuess(X,ierr)
      implicit none

!     petscsys.h       - base PETSc routines   petscvec.h - vectors
!     petscmat.h - matrices
!     petscis.h     - index sets            petscksp.h - Krylov subspace methods
!     petscviewer.h - viewers               petscpc.h  - preconditioners

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
      PetscErrorCode    ierr
      PetscOffset      idx
      Vec       X,localX
      PetscInt  i,j,row,mx
      PetscInt  my, xs,ys,xm
      PetscInt  ym
      PetscReal one,lambda,temp1,temp,hx,hy
      PetscScalar      xx(2)
      DM               da
      Mat              B
      common   /mycommon/ mx,my,B,localX,da
      one    = 1.0
      lambda = 0.0
      hx     = one/(mx-1)
      hy     = one/(my-1)
      temp1  = lambda/(lambda + one)

!  Get a pointer to vector data.
!    - VecGetArray() returns a pointer to the data array.
!    - You MUST call VecRestoreArray() when you no longer need access to
!      the array.
       call VecGetArray(X,xx,idx,ierr)
       

!  Get local grid boundaries (for 2-dimensional DMDA):
!    xs, ys   - starting grid indices (no ghost points)
!    xm, ym   - widths of local grid (no ghost points)

       call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,             &
     &      PETSC_NULL_INTEGER,ierr)
       

!  Compute initial guess over the locally owned part of the grid
     
      do 30 j=ys,ys+ym-1
        temp = (min(j,my-j-1))*hy
        !print *,'temp',temp
        do 40 i=xs,xs+xm-1
          row = i - xs + (j - ys)*xm + 1
         ! print*,'row',row
          if (i .eq. 0 .or. j .eq. 0) then
          
            xx(idx+row) = 50
	  elseif (i .eq. mx-1 .or. j .eq. my-1) then
            xx(idx+row) = 30
          else
            xx(idx+row) = 10.0
          endif
         
          
 40     continue
 30   continue
!       call VecView(X,PETSC_VIEWER_STDOUT_WORLD,ierr)


!     Restore vector

       call VecRestoreArray(X,xx,idx,ierr)
       
       return
       end

! -------------------------------------------------------------------
!
!   ComputeFunction - Evaluates nonlinear function, F(x).
!
!   Input Parameters:
!.  X - input vector
!
!   Output Parameter:
!.  F - function vector
!
      subroutine  ComputeFunction(X,F,ierr)
      implicit none

!     petscsys.h       - base PETSc routines   petscvec.h - vectors
!     petscmat.h - matrices
!     petscis.h     - index sets            petscksp.h - Krylov subspace methods
!     petscviewer.h - viewers               petscpc.h  - preconditioners

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

      Vec              X,F,localX
      PetscInt         gys,gxm,gym
      PetscOffset      idx,idf
      PetscErrorCode ierr
      PetscInt i,j,row,mx,my,xs,ys,xm,ym,gxs
      PetscInt rowf
      PetscReal two,one,lambda,hx
      PetscReal hy,hxdhy,hydhx,sc
      PetscScalar      u,uxx,uyy,xx(2),ff(2)
      DM               da
      Mat              B
      common   /mycommon/ mx,my,B,localX,da

      two    = 2.0
      one    = 1.0
      lambda = 6.0

      hx     = one/(mx-1)
      hy     = one/(my-1)
      sc     = hx*hy*lambda
      hxdhy  = hx/hy
      hydhx  = hy/hx

!  Scatter ghost points to local vector, using the 2-step process
!     DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
!  By placing code between these two statements, computations can be
!  done while messages are in transition.
!
      call DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX,ierr)
      call DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX,ierr)

!  Get pointers to vector data

      call VecGetArray(localX,xx,idx,ierr)
      call VecGetArray(F,ff,idf,ierr)
      
!  Get local grid boundaries

      call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,              &
     &     PETSC_NULL_INTEGER,ierr)
      call DMDAGetGhostCorners(da,gxs,gys,PETSC_NULL_INTEGER,gxm,gym,     &
     &     PETSC_NULL_INTEGER,ierr)
!      print *,xs,ys,xm,ym,gxs,gys,gxm,gym

!  Compute function over the locally owned part of the grid
      rowf = 0
      do 50 j=ys,ys+ym-1

        row  = (j - gys)*gxm + xs - gxs
!        print *,row
        do 60 i=xs,xs+xm-1
          row  = row + 1
          rowf = rowf + 1

          if (i .eq. 0 .or. j .eq. 0 .or. i .eq. mx-1 .or.              &
     &        j .eq. my-1) then
            ff(idf+rowf) = xx(idx+row)
            goto 60
          endif
          u   = xx(idx+row)
          uxx = ( xx(idx+row-1) + xx(idx+row+1))*hydhx
          uyy = ( xx(idx+row-gxm) + xx(idx+row+gxm))*hxdhy
!check here
          ff(idf+rowf) = 0.25*(uxx + uyy) 
!imp
!          print *,ff(idf+rowf),xx(idx+row)
 60     continue
 50   continue
!  Restore vectors
!       call VecView(F,PETSC_VIEWER_STDOUT_WORLD,ierr)
       call VecRestoreArray(localX,xx,idx,ierr)
       call VecRestoreArray(F,ff,idf,ierr)
       return
       end



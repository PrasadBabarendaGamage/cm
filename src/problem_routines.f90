!> \file
!> \author Chris Bradley
!> \brief This module handles all problem routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all problem routines.
MODULE PROBLEM_ROUTINES

  USE BASE_ROUTINES
  USE BIOELECTRIC_ROUTINES
  USE CLASSICAL_FIELD_ROUTINES
  USE CMISS_PETSC
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ELASTICITY_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE FITTING_ROUTINES
  USE FLUID_MECHANICS_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_CONDITIONS_ROUTINES
  USE INTERFACE_OPERATORS_ROUTINES
  USE INTERFACE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE MULTI_PHYSICS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE RIGID_BODY_ROUTINES
  USE SOLVER_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(PROBLEMS_TYPE), TARGET :: PROBLEMS
  
  !Interfaces

  INTERFACE PROBLEM_CELLML_EQUATIONS_GET
    MODULE PROCEDURE PROBLEM_CELLML_EQUATIONS_GET_0
    MODULE PROCEDURE PROBLEM_CELLML_EQUATIONS_GET_1
  END INTERFACE !PROBLEM_CELLML_EQUATIONS_GET

  INTERFACE PROBLEM_CONTROL_LOOP_GET
    MODULE PROCEDURE PROBLEM_CONTROL_LOOP_GET_0
    MODULE PROCEDURE PROBLEM_CONTROL_LOOP_GET_1
  END INTERFACE !PROBLEM_CONTROL_LOOP_GET

  INTERFACE PROBLEM_SOLVER_EQUATIONS_GET
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_GET_0
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_GET_1
  END INTERFACE !PROBLEM_SOLVER_EQUATIONS_GET

  INTERFACE PROBLEM_SOLVER_GET
    MODULE PROCEDURE PROBLEM_SOLVER_GET_0
    MODULE PROCEDURE PROBLEM_SOLVER_GET_1
  END INTERFACE !PROBLEM_SOLVER_GET
  
  PUBLIC PROBLEMS_INITIALISE,PROBLEMS_FINALISE
  
  PUBLIC PROBLEM_CELLML_EQUATIONS_CREATE_START,PROBLEM_CELLML_EQUATIONS_CREATE_FINISH
  
  PUBLIC PROBLEM_CELLML_EQUATIONS_GET
  
  PUBLIC PROBLEM_CREATE_START,PROBLEM_CREATE_FINISH,PROBLEM_DESTROY
  
  PUBLIC PROBLEM_SPECIFICATION_GET,PROBLEM_SPECIFICATION_SET
  
  PUBLIC PROBLEM_CONTROL_LOOP_CREATE_START,PROBLEM_CONTROL_LOOP_CREATE_FINISH
  
  PUBLIC PROBLEM_CONTROL_LOOP_DESTROY
  
  PUBLIC PROBLEM_CONTROL_LOOP_GET
  
  PUBLIC PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC

  PUBLIC PROBLEM_SOLVER_EQUATIONS_CREATE_START,PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH
  
  PUBLIC PROBLEM_SOLVER_EQUATIONS_DESTROY
  
  PUBLIC PROBLEM_SOLVER_EQUATIONS_GET
  
  PUBLIC PROBLEM_SOLVER_JACOBIAN_EVALUATE,PROBLEM_SOLVER_RESIDUAL_EVALUATE
  
  PUBLIC PROBLEM_SOLVER_GET
  
  PUBLIC Problem_SolverNonlinearMonitor
  
  PUBLIC PROBLEM_SOLVE
  
  PUBLIC PROBLEM_SOLVERS_CREATE_START,PROBLEM_SOLVERS_CREATE_FINISH
  
  PUBLIC PROBLEM_SOLVERS_DESTROY
  
  PUBLIC PROBLEM_USER_NUMBER_FIND
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of the CellML equations for the problem solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsCreateFinish
  SUBROUTINE PROBLEM_CELLML_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the CellML equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_CELLML_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN      
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_FINISH_ACTION
      !Finish problem specific startup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CELLML_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_CELLML_EQUATIONS_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of CellML equations for a problem solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsCreateStart
  SUBROUTINE PROBLEM_CELLML_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variablesg
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of the CellML equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_CELLML_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
      !Start the problem specific control setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CELLML_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_CELLML_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations defined with a solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsGet
  SUBROUTINE PROBLEM_CELLML_EQUATIONS_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,CELLML_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier to get the solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the solver CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS !<On exit, a pointer to the specified solver CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_CELLML_EQUATIONS_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_CELLML_EQUATIONS_GET_1(PROBLEM,[CONTROL_LOOP_IDENTIFIER],SOLVER_INDEX,CELLML_EQUATIONS,ERR,ERROR,*999)
    
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_CELLML_EQUATIONS_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_CELLML_EQUATIONS_GET_0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver CellML equations defined with a solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsGet
  SUBROUTINE PROBLEM_CELLML_EQUATIONS_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,CELLML_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_CELLML_EQUATIONS_GET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(CELLML_EQUATIONS)) THEN
        CALL FLAG_ERROR("The CellML equations is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(CELLML_EQUATIONS)
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN            
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(ASSOCIATED(SOLVER)) THEN
                CELLML_EQUATIONS=>SOLVER%CELLML_EQUATIONS
                IF(.NOT.ASSOCIATED(CELLML_EQUATIONS)) CALL FLAG_ERROR("CellML equations is not associated.",ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)          
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_GET_1")
    RETURN
999 CALL ERRORS("PROBLEM_CELLML_EQUATIONS_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_CELLML_EQUATIONS_GET_1
  
  !
  !================================================================================================================================
  !

  !>Solves CellML equations for a problem.
  SUBROUTINE PROBLEM_CELLML_EQUATIONS_SOLVE(CELLML_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS !<A pointer to the CellML equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    
    CALL ENTERS("PROBLEM_CELLML_EQUATIONS_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CELLML_EQUATIONS)) THEN
      IF(CELLML_EQUATIONS%CELLML_EQUATIONS_FINISHED) THEN
        SOLVER=>CELLML_EQUATIONS%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          
          CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
          
        ELSE
          CALL FLAG_ERROR("CellML equations solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_CELLML_EQUATIONS_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_CELLML_EQUATIONS_SOLVE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_CELLML_EQUATIONS_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves a problem control loop.
  RECURSIVE SUBROUTINE PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: iteration_idx,loop_idx,solver_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER :: SIMPLE_LOOP
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER :: LOAD_INCREMENT_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_CONTROL_LOOP_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        !Solve this control loop
        IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Control loop: ",CONTROL_LOOP%LABEL,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Control loop level = ",CONTROL_LOOP%CONTROL_LOOP_LEVEL,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Sub loop index     = ",CONTROL_LOOP%SUB_LOOP_INDEX,ERR,ERROR,*999)
        ENDIF
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          SIMPLE_LOOP=>CONTROL_LOOP%SIMPLE_LOOP
          IF(ASSOCIATED(SIMPLE_LOOP)) THEN
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Simple control loop: ",ERR,ERROR,*999)
            ENDIF
            CALL PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
              !If there are no sub loops then solve.
              SOLVERS=>CONTROL_LOOP%SOLVERS
              IF(ASSOCIATED(SOLVERS)) THEN
                DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                  SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR

                  CALL PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)

                ENDDO !solver_idx
              ELSE
                CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              !If there are sub loops the recursively solve those control loops
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
              ENDDO !loop_idx
            ENDIF
            CALL PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Control loop simple loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
          FIXED_LOOP=>CONTROL_LOOP%FIXED_LOOP
          IF(ASSOCIATED(FIXED_LOOP)) THEN
            DO iteration_idx=FIXED_LOOP%START_ITERATION,FIXED_LOOP%STOP_ITERATION,FIXED_LOOP%ITERATION_INCREMENT
              IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Fixed control loop iteration: ",iteration_idx,ERR,ERROR,*999)
              ENDIF
              FIXED_LOOP%ITERATION_NUMBER=iteration_idx
              CALL PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    
                    CALL PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)

                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
              CALL PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            ENDDO !iteration_idx
          ELSE
            CALL FLAG_ERROR("Control loop fixed loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            !Set the current time to be the start time. Solvers should use the first time step to do any initialisation.
            TIME_LOOP%CURRENT_TIME=TIME_LOOP%START_TIME
            TIME_LOOP%ITERATION_NUMBER=0
            DO WHILE(TIME_LOOP%CURRENT_TIME<=TIME_LOOP%STOP_TIME)
              IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Time control loop iteration: ",TIME_LOOP%ITERATION_NUMBER, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Current time   = ",TIME_LOOP%CURRENT_TIME, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Stop time      = ",TIME_LOOP%STOP_TIME, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Time increment = ",TIME_LOOP%TIME_INCREMENT, &
                  & ERR,ERROR,*999)
              ENDIF
              !Perform any pre-loop actions.
              CALL PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve.
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    
                    CALL PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
                    
                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
              !Perform any post loop actions.
              CALL PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
              !Increment loop counter and time
              TIME_LOOP%ITERATION_NUMBER=TIME_LOOP%ITERATION_NUMBER+1
              TIME_LOOP%GLOBAL_ITERATION_NUMBER=TIME_LOOP%GLOBAL_ITERATION_NUMBER+1
              TIME_LOOP%CURRENT_TIME=TIME_LOOP%CURRENT_TIME+TIME_LOOP%TIME_INCREMENT
            ENDDO !time loop
          ELSE
            CALL FLAG_ERROR("Control loop time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          WHILE_LOOP=>CONTROL_LOOP%WHILE_LOOP
          IF(ASSOCIATED(WHILE_LOOP)) THEN
            WHILE_LOOP%ITERATION_NUMBER=0
            WHILE_LOOP%CONTINUE_LOOP=.TRUE.
            DO WHILE(WHILE_LOOP%CONTINUE_LOOP.AND.WHILE_LOOP%ITERATION_NUMBER &
              & <WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS)
              WHILE_LOOP%ITERATION_NUMBER=WHILE_LOOP%ITERATION_NUMBER+1
              IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"While control loop iteration: ",WHILE_LOOP%ITERATION_NUMBER, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of iterations = ", &
                  & WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
              ENDIF
              CALL PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    
                    CALL PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
                    
                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
              CALL PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            ENDDO !while loop
          ELSE
            CALL FLAG_ERROR("Control loop while loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          LOAD_INCREMENT_LOOP=>CONTROL_LOOP%LOAD_INCREMENT_LOOP
          IF(ASSOCIATED(LOAD_INCREMENT_LOOP)) THEN
            LOAD_INCREMENT_LOOP%ITERATION_NUMBER=0
            IF (LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS<1) THEN
              ! automatic stepping
              CALL FLAG_ERROR("Automatic load incrementing is not implemented yet.",ERR,ERROR,*999)
            ELSE
              ! fixed number of steps
              DO WHILE(LOAD_INCREMENT_LOOP%ITERATION_NUMBER<LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS)
                LOAD_INCREMENT_LOOP%ITERATION_NUMBER=LOAD_INCREMENT_LOOP%ITERATION_NUMBER+1
                IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Load increment control loop iteration: ", &
                    & LOAD_INCREMENT_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of iterations = ", &
                    & LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
                ENDIF
                CALL PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
                IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                  !If there are no sub loops then solve
                  SOLVERS=>CONTROL_LOOP%SOLVERS
                  IF(ASSOCIATED(SOLVERS)) THEN
                    DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                      SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                      IF(ASSOCIATED(SOLVER)) THEN
                        IF(ASSOCIATED(SOLVER%SOLVER_EQUATIONS)) THEN
                          !Apply incremented boundary conditions here => 
                          CALL PROBLEM_SOLVER_LOAD_INCREMENT_APPLY(SOLVER%SOLVER_EQUATIONS,LOAD_INCREMENT_LOOP%ITERATION_NUMBER, &
                            & LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
                        ENDIF
                        CALL PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO !solver_idx
                  ELSE
                    CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  !If there are sub loops the recursively solve those control loops
                  DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                    CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                    CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                  ENDDO !loop_idx
                ENDIF
                CALL PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
              ENDDO !while loop
            ENDIF
          ELSE
            CALL FLAG_ERROR("Control loop while loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The control loop loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_SOLVE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a problem. \see OPENCMISS::CMISSProblemCreateFinish
  SUBROUTINE PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_INITIAL_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_FINISH_ACTION
      !Finish the problem specific setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finish the problem creation
      PROBLEM%PROBLEM_FINISHED=.TRUE.
    ELSE        
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of problems = ",PROBLEMS%NUMBER_OF_PROBLEMS,ERR,ERROR,*999)
      DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Problem number  = ",problem_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  User number     = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global number   = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem class   = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%CLASS, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem type    = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%TYPE, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem subtype = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%SUBTYPE, &
          & ERR,ERROR,*999)
      ENDDO !problem_idx    
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE PROBLEM_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a problem defined by USER_NUMBER. \see OPENCMISS::CMISSProblemCreateStart
  !>The default values of the PROBLEM attributes are:
  !>- CLASS: 4 (PROBLEM_CLASSICAL_FIELD_CLASS)
  !>- TYPE: 1 (PROBLEM_LAPLACE_EQUATION_TYPE)
  !>- SUBTYPE: 1 (PROBLEM_STANDARD_LAPLACE_SUBTYPE)
  SUBROUTINE PROBLEM_CREATE_START(USER_NUMBER,PROBLEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to create
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<On return, a pointer to the created problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(PROBLEM_TYPE), POINTER :: NEW_PROBLEM
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    NULLIFY(NEW_PROBLEM)
    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL FLAG_ERROR("Problem is already associated.",ERR,ERROR,*999)
    ELSE
      NULLIFY(PROBLEM)
      CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*999)
      IF(ASSOCIATED(PROBLEM)) THEN
        LOCAL_ERROR="Problem number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))//" has already been created."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        !Allocate the new problem
        ALLOCATE(NEW_PROBLEM,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problem.",ERR,ERROR,*999)
        !Initalise problem
        CALL PROBLEM_INITIALISE(NEW_PROBLEM,ERR,ERROR,*999)
        !Set default problem values
        NEW_PROBLEM%USER_NUMBER=USER_NUMBER
        NEW_PROBLEM%GLOBAL_NUMBER=PROBLEMS%NUMBER_OF_PROBLEMS+1
        NEW_PROBLEM%PROBLEMS=>PROBLEMS
        !Default to a standardised Laplace.
        NEW_PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        NEW_PROBLEM%TYPE=PROBLEM_LAPLACE_EQUATION_TYPE
        NEW_PROBLEM%SUBTYPE=PROBLEM_STANDARD_LAPLACE_SUBTYPE
        NEW_PROBLEM%PROBLEM_FINISHED=.FALSE.
        !Initialise the problem setup information
        CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_INITIAL_TYPE
        PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
        !Start problem specific setup
        CALL PROBLEM_SETUP(NEW_PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the problem setup information
        CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        !Add new problem into list of problems
        ALLOCATE(NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems.",ERR,ERROR,*999)
        DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
          NEW_PROBLEMS(problem_idx)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
        ENDDO !problem_idx
        NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS+1)%PTR=>NEW_PROBLEM
        IF(ASSOCIATED(PROBLEMS%PROBLEMS)) DEALLOCATE(PROBLEMS%PROBLEMS)
        PROBLEMS%PROBLEMS=>NEW_PROBLEMS
        PROBLEMS%NUMBER_OF_PROBLEMS=PROBLEMS%NUMBER_OF_PROBLEMS+1
        PROBLEM=>NEW_PROBLEM
      ENDIF
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN 1   
  END SUBROUTINE PROBLEM_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys a problem. \see OPENCMISS::CMISSProblemDestroy
  SUBROUTINE PROBLEM_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx,problem_position
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)

    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEMS%PROBLEMS)) THEN
        
        problem_position=PROBLEM%GLOBAL_NUMBER
      
        !Destroy all the problem components
        CALL PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*999)
        
        !Remove the problem from the list of problems
        IF(PROBLEMS%NUMBER_OF_PROBLEMS>1) THEN
          ALLOCATE(NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems.",ERR,ERROR,*999)
          DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
            IF(problem_idx<problem_position) THEN
              NEW_PROBLEMS(problem_idx)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
            ELSE IF(problem_idx>problem_position) THEN
              PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER=PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER-1
              NEW_PROBLEMS(problem_idx-1)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
            ENDIF
          ENDDO !problem_idx
          DEALLOCATE(PROBLEMS%PROBLEMS)
          PROBLEMS%PROBLEMS=>NEW_PROBLEMS
          PROBLEMS%NUMBER_OF_PROBLEMS=PROBLEMS%NUMBER_OF_PROBLEMS-1
        ELSE
          DEALLOCATE(PROBLEMS%PROBLEMS)
          PROBLEMS%NUMBER_OF_PROBLEMS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Problem problems is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("PROBLEM_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_PROBLEMS)) DEALLOCATE(NEW_PROBLEMS)
998 CALL ERRORS("PROBLEM_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_DESTROY")
    RETURN 1   
  END SUBROUTINE PROBLEM_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalise the problem setup and deallocate all memory.
  SUBROUTINE PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SETUP_TYPE), INTENT(OUT) :: PROBLEM_SETUP_INFO !<The problem setup to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SETUP_FINALISE",ERR,ERROR,*999)

    PROBLEM_SETUP_INFO%SETUP_TYPE=0
    PROBLEM_SETUP_INFO%ACTION_TYPE=0
       
    CALL EXITS("PROBLEM_SETUP_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SETUP_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SETUP_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SETUP_FINALISE

 !
  !================================================================================================================================
  !

  !>Initialise the problem setup.
  SUBROUTINE PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SETUP_TYPE), INTENT(OUT) :: PROBLEM_SETUP_INFO !<The problem setup to intialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SETUP_INITIALISE",ERR,ERROR,*999)

    PROBLEM_SETUP_INFO%SETUP_TYPE=0
    PROBLEM_SETUP_INFO%ACTION_TYPE=0
        
    CALL EXITS("PROBLEM_SETUP_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SETUP_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SETUP_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SETUP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the problem and deallocate all memory.
  SUBROUTINE PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) CALL CONTROL_LOOP_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      DEALLOCATE(PROBLEM)
    ENDIF
       
    CALL EXITS("PROBLEM_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a problem.
  SUBROUTINE PROBLEM_INITIALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<The pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The errror code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      PROBLEM%USER_NUMBER=0
      PROBLEM%GLOBAL_NUMBER=0
      PROBLEM%PROBLEM_FINISHED=.FALSE.
      NULLIFY(PROBLEM%PROBLEMS)
      PROBLEM%CLASS=PROBLEM_NO_CLASS
      PROBLEM%TYPE=PROBLEM_NO_TYPE
      PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE
      NULLIFY(PROBLEM%CONTROL_LOOP)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finish the creation of the control for the problem. \see OPENCMISS::CMISSProblemControlLoopCreateFinish
  SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the control for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_CONTROL_LOOP_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        IF(PROBLEM%CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
          CALL FLAG_ERROR("Problem control loop has already been finished.",ERR,ERROR,*999)
        ELSE
          !Initialise the problem setup information
          CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
          PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_CONTROL_TYPE
          PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_FINISH_ACTION
          !Finish problem specific startup
          CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the problem setup information
          CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
          !Finish problem control creation
          PROBLEM%CONTROL_LOOP%CONTROL_LOOP_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a control loop for a problem. \see OPENCMISS::CMISSProblemControlLoopCreateStart
  !>The default values of the PROBLEM CONTROL LOOP attributes are:
  !>- LOOP_TYPE: PROBLEM_CONTROL_SIMPLE_TYPE
  !>- CONTROL_LOOP_LEVEL: 1
  !>- NUMBER_OF_SUB_LOOPS: 0
  SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of a control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_CONTROL_LOOP_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("The problem control loop is already associated.",ERR,ERROR,*999)        
      ELSE
        !Initialise the problem setup information
        CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_CONTROL_TYPE
        PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
        !Start the problem specific control setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the problem setup information
        CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the control loop for a problem. \see OPENCMISS::CMISSProblemControlLoopDestroy
  SUBROUTINE PROBLEM_CONTROL_LOOP_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_LOOP_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN        
        CALL CONTROL_LOOP_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_DESTROY

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a problem. \see OPENCMISS::CMISSProblemControlLoopGet
  SUBROUTINE PROBLEM_CONTROL_LOOP_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_CONTROL_LOOP_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_CONTROL_LOOP_GET_1(PROBLEM,[CONTROL_LOOP_IDENTIFIER],CONTROL_LOOP,ERR,ERROR,*999) 
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_GET_0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control_loop for a problem. \see OPENCMISS::CMISSProblemControlLoopGet
  SUBROUTINE PROBLEM_CONTROL_LOOP_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT
 
    CALL ENTERS("PROBLEM_CONTROL_LOOP_GET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("Control loop is already associated.",ERR,ERROR,*999)
      ELSE
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_1")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_GET_1
  
  !
  !================================================================================================================================
  !

  !>Sets up the specifices for a problem.
  SUBROUTINE PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP_INFO !<The problem setup information.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        CALL ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE(PROBLEM_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE(PROBLEM_BIOELECTRICS_CLASS)
        CALL BIOELECTRIC_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE(PROBLEM_FITTING_CLASS)
        CALL FITTING_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE(PROBLEM_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OPENCMISS::CMISSProblemSolverEquationsGet
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_SOLVER_EQUATIONS_GET_1(PROBLEM,[CONTROL_LOOP_IDENTIFIER],SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*999)
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OPENCMISS::CMISSProblemSolverEquationsGet
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_GET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        CALL FLAG_ERROR("The solver equations is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(SOLVER_EQUATIONS)
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN            
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(ASSOCIATED(SOLVER)) THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(.NOT.ASSOCIATED(SOLVER_EQUATIONS)) CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)          
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_1")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_1
  
  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear problem solver.
  SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE(SOLVER,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, LINKING_SOLVER !<A pointer to the solver to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,solver_matrix_idx,noComp
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVER_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            IF(SOLVER%OUTPUT_TYPE>=SOLVER_MATRIX_OUTPUT) THEN
              SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Solver vector values:",ERR,ERROR,*999)
                DO solver_matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                  SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                  IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Solver matrix is not associated for solver matrix index "// &
                      & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !solver_matrix_idx
              ELSE
                CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
              !Check if the nonlinear solver is linked to a dynamic solver 
              LINKING_SOLVER=>SOLVER%LINKING_SOLVER
              IF(ASSOCIATED(LINKING_SOLVER)) THEN
                IF(LINKING_SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  !Update the field values from the dynamic factor * current solver values AND add in mean predicted displacements/
                  CALL SOLVER_VARIABLES_DYNAMIC_NONLINEAR_UPDATE(SOLVER,ERR,ERROR,*999)
                  !Calculate the Jacobian
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    !Assemble the equations for dynamic problems
                    CALL EQUATIONS_SET_JACOBIAN_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                  ENDDO !equations_set_idx
                  !Assemble the dynamic nonlinear solver matrices
                  CALL SOLVER_MATRICES_DYNAMIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_JACOBIAN_ONLY,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Solver equations linking solver mapping is not dynamic.",ERR,ERROR,*999)
                END IF
              ELSE
                !Otherwise perform as steady nonlinear
                !Copy the current solution vector to the dependent field
                CALL SOLVER_VARIABLES_FIELD_UPDATE(SOLVER,ERR,ERROR,*999)
                !Calculate the Jacobian
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Jacobian evaluation******************",ERR,ERROR,*999)
                  !\todo: XY rigid-deformable contact, temporarily change the number of components to exclude rigid body dofs
!                  noComp=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%NUMBER_OF_COMPONENTS
!                  EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%NUMBER_OF_COMPONENTS= &
!                    & noComp-6
                  !Assemble the equations for linear problems
                  CALL EQUATIONS_SET_JACOBIAN_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                  !\todo: XY rigid-deformable contact, restore number of components
!                  EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%NUMBER_OF_COMPONENTS=noComp
                  !\todo: XY -Modify Jacobian for contact
                  IF(SOLVER%SOLVERS%CONTROL_LOOP%PROBLEM%TYPE==PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)THEN
                    ! rigid body-deformable body contact
                    IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_MULTI_PHYSICS_CLASS .AND.  &
                        & EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_RIGID_BODY_TYPE) THEN
                      CALL EquationsSet_JacobianRigidBodyContactUpdateStaticFEM(EQUATIONS_SET,ERR,ERROR,*999)
                    ! deformable-deformable body contact
                    ELSEIF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                      CALL EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*999)
                    ENDIF !Rigid or deformable contact
                  ENDIF !contact problem
                ENDDO !equations_set_idx
                !Update interface matrices
!                DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
!                  interfaceCondition=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
!                  !Assemble the interface condition for the Jacobian LHS
!                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Jacobian evaluation******************",ERR,ERROR,*999)
!                  CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
!                ENDDO
                !Assemble the static nonlinear solver matrices
                CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_JACOBIAN_ONLY,ERR,ERROR,*999)
              END IF       
            ELSE
              CALL FLAG_ERROR("Solver equations solver type is not associated.",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a static equations set with contact using the finite element method
  SUBROUTINE EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    !TYPE(EQUATIONS_SET_TYPE), POINTER :: coupledRegionEquationsSet(2) !<A pointer to the equations set to evaluate the Jacobian for
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition  !<A pointer to the equations set to evaluate the element Jacobian for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    TYPE(InterfaceContactPointMetricsType), POINTER :: contactPointMetrics
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(BASIS_TYPE), POINTER :: rowDependentBasis,colDependentBasis,rowDomainFaceBasis,colDomainFaceBasis
    TYPE(DOMAIN_FACE_TYPE), POINTER :: rowDomainFace,colDomainFace
    TYPE(EQUATIONS_SET_TYPE), POINTER :: multipleRegionEquationsSet
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: jacobian
    !TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    !TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: residualParameterSet
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: jacobianNumber,bodyIdx,equationSetNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber
    INTEGER(INTG) :: globalDataPointNum,rowElementNum,colElementNum, &
      & rowConnectedFace,colConnectedFace,rowFieldComp, colFieldComp, &
      & rowMeshComp,colMeshComp,rowDecompositionFaceNumber,colDecompositionFaceNumber, &
      & rowLocalFaceNodeIdx,colLocalFaceNodeIdx,rowFaceLocalElemNode,colFaceLocalElemNode,rowGlobalNode,colGlobalNode, &
      & rowFaceDerivative,colFaceDerivative,rowDerivative,colDerivative,rowVersion,colVersion,rowIdx,colIdx, &
      & subMatrix,rowBodyIdx,colBodyIdx
    INTEGER(INTG) :: xiIdxAlpha,xiIdxBeta,xiIdxGamma,rowElemParameterNo,colElemParameterNo,rowPreviousFaceNo,colPreviousFaceNo
    REAL(DP) :: matrixValue,rowPhi,colPhi
    REAL(DP) :: coefficient,forceTerm,geometricTerm,tempA,tempB,rowDofScaleFactor,colDofScaleFactor
    REAL(DP) :: rowXi(2),colXi(2) !\todo generalise xi allocations for 1D,2D and 3D points connectivity
    REAL(DP) :: kappa(2,2),phiDeriRow(2),phiDeriCol(2),TRow(2),TCol(2),NRow(2),NCol(2),DRow(2),DCol(2)
  
    CALL ENTERS("EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(dependentField)) THEN
        equations=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          equationsMatrices=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(equationsMatrices)) THEN
            nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
            nonlinearMapping=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING

            interfaceGlobalNumber=1
            interfaceConditionGlobalNumber=1
            interface=>EQUATIONS_SET%REGION%PARENT_REGION%INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR
            interfaceCondition=>interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
            pointsConnectivity=>interface%pointsConnectivity
            contactMetrics=>interfaceCondition%interfaceContactMetrics

            dependentVariable=>dependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)
            jacobianNumber=1
            jacobian=>nonlinearMatrices%JACOBIANS(jacobianNumber)%PTR%JACOBIAN
            dependentVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(jacobianNumber)%VARIABLE
            !Setup pointer to the equation set of the coupled bodies which are setup in thier own separate regions
            !(note that these regions has not been added to the solver equations and are merely here for convinence if needed)
            !equationSetNumber=1
            !DO bodyIdx=1,2
            !  coupledRegionEquationsSet(bodyIdx)=>EQUATIONS_SET%REGION%PARENT_REGION%SUB_REGIONS(bodyIdx)%PTR% &
            !    & EQUATIONS_SETS%EQUATIONS_SETS(equationSetNumber)%PTR
            !ENDDO
            !Since we are computing the contact term in a single region, we do not need to determine the dependent field
            !through the interface condition. We can use the dependentField pointer defined for this single region
            !dependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(interfaceMatrixIdx)%PTR% &
            ! & DEPENDENT%DEPENDENT_FIELD

            DO subMatrix=1,4
              SELECT CASE(subMatrix)
              CASE(1) !Contact subMatrix11
                rowBodyIdx=1
                colBodyIdx=1
                coefficient=1;
              CASE(2) !Contact subMatrix12
                rowBodyIdx=1
                colBodyIdx=2
                coefficient=-1;
              CASE(3) !Contact subMatrix21
                rowBodyIdx=2
                colBodyIdx=1
                coefficient=-1;
              CASE(4) !Contact subMatrix22
                rowBodyIdx=2
                colBodyIdx=2
                coefficient=1;
              END SELECT
              !Loop over each data point and find the connected element and their dofs
              DO globalDataPointNum=1,SIZE(pointsConnectivity%pointsConnectivity,1)
                IF(contactMetrics%inContact(globalDataPointNum)) THEN
                  ! Get the metric structure for this contact point
                  contactPointMetrics=>contactMetrics%contactPointMetrics(globalDataPointNum)
                  rowElementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%coupledMeshElementNumber
                  rowConnectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%elementLineFaceNumber
                  rowXi=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%reducedXi
                  colElementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%coupledMeshElementNumber
                  colConnectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%elementLineFaceNumber
                  colXi=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%reducedXi
                  
                  !################################################################################################################
                  IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                    !Calculate quantities that do not vary w.r.t xyz or node
                    kappa=0.0_DP
                    !Calculate kappa (see Jae's cm implementation)
                    DO xiIdxAlpha=1,2
                      DO xiIdxBeta=1,2
                        kappa(xiIdxAlpha,xiIdxBeta)= & 
                          & DOT_PRODUCT(contactPointMetrics%tangentDerivatives(xiIdxAlpha,xiIdxBeta,:), &
                          & contactPointMetrics%normal(:))
                      ENDDO !xiIdxBeta
                    ENDDO !xiIdxAlpha
                  ENDIF !addGeometricTerm
                  !################################################################################################################
                  rowPreviousFaceNo=0
                  !Find the row dof 
                  DO rowFieldComp=1,3
                    rowMeshComp=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                      & COMPONENTS(rowFieldComp)%MESH_COMPONENT_NUMBER
                    rowDependentBasis=>dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR% &
                      & TOPOLOGY%ELEMENTS%ELEMENTS(rowElementNum)%BASIS
                    rowDecompositionFaceNumber=dependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(rowElementNum)%ELEMENT_FACES(rowConnectedFace)
                    rowDomainFace=>dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                      & FACES%FACES(rowDecompositionFaceNumber)
                    rowDomainFaceBasis=>rowDomainFace%BASIS
                    !Only interpolate for the first field component and when face number changes
                    IF((rowFieldComp==1) .AND. (rowDecompositionFaceNumber/=rowPreviousFaceNo)) THEN
                      CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(rowDecompositionFaceNumber, &
                        & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                      rowPreviousFaceNo=rowDecompositionFaceNumber
                    ENDIF
                    DO rowLocalFaceNodeIdx=1,rowDependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(rowConnectedFace)
                      rowFaceLocalElemNode=rowDependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(rowLocalFaceNodeIdx,rowConnectedFace)
                      rowGlobalNode=dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(rowElementNum)%ELEMENT_NODES(rowFaceLocalElemNode)
                      DO rowFaceDerivative=1,rowDomainFace%BASIS%NUMBER_OF_DERIVATIVES(rowLocalFaceNodeIdx)
                        rowDerivative=rowDependentBasis% &
                          & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(rowFaceDerivative,rowLocalFaceNodeIdx,rowConnectedFace)
                        rowVersion=dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                          & ELEMENTS%ELEMENTS(rowElementNum)%elementVersions(rowDerivative,rowFaceLocalElemNode)
                        !Find the face parameter's element parameter index 
                        rowElemParameterNo=rowDomainFaceBasis%ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx)
                        rowDofScaleFactor=equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% PTR% &
                          & SCALE_FACTORS(rowElemParameterNo,rowFieldComp)
                        !Find dof associated with this particular field, component, node, derivative and version.
                        rowIdx=dependentVariable%components(rowFieldComp)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES( &
                          & rowGlobalNode)%DERIVATIVES(rowDerivative)%VERSIONS(rowVersion)
                        !Evaluate the basis at the projected/connected xi
                        rowPhi=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                          & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),NO_PART_DERIV,rowXi,err,error)
                        
                        !###########################################################################################################
                        IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                          !Calculate row metrics
                          !\todo: generalise contact xi direction, at the moment assume in xi1 and xi2
                          phiDeriRow(1)=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                            & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),PART_DERIV_S1,rowXi,err,error)
                          phiDeriRow(2)=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                            & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),PART_DERIV_S2,rowXi,err,error)  
                            
                          !NRow TRow changes at every row dof
                          DO xiIdxAlpha=1,2
                            IF (rowBodyIdx==1) THEN
                              NRow(xiIdxAlpha)=0.0_DP
                              TRow(xiIdxAlpha)=rowPhi*contactPointMetrics%tangents(xiIdxAlpha,rowFieldComp)
                            ELSE
                              NRow(xiIdxAlpha)=-phiDeriRow(xiIdxAlpha)*contactPointMetrics%normal(rowFieldComp) 
                              TRow(xiIdxAlpha)=-rowPhi*contactPointMetrics%tangents(xiIdxAlpha,rowFieldComp)
                            ENDIF !rowBodyIdx
                          ENDDO !xiIdxAlpha    
                          
                          !DRow varies at every row dof
                          DRow=0.0_DP
                          DO xiIdxAlpha=1,2
                            DO xiIdxBeta=1,2
                              DRow(xiIdxAlpha)=DRow(xiIdxAlpha)+contactPointMetrics%inverseA(xiIdxAlpha,xiIdxBeta)* &
                                & (TRow(xiIdxBeta)+contactPointMetrics%signedGapNormal*NRow(xiIdxBeta))
                            ENDDO !xiIdxBeta
                          ENDDO !xiIdxAlpha   
                        ENDIF !addGeometricTerm
                        !###########################################################################################################

                        colPreviousFaceNo=0
                        DO colFieldComp=1,3
                          !Find the col 
                          colMeshComp=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & COMPONENTS(colFieldComp)%MESH_COMPONENT_NUMBER
                          colDependentBasis=>dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR% &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(colElementNum)%BASIS
                          colDecompositionFaceNumber=dependentField%DECOMPOSITION%TOPOLOGY% &
                            & ELEMENTS%ELEMENTS(colElementNum)%ELEMENT_FACES(colConnectedFace)
                          colDomainFace=>dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                            & FACES%FACES(colDecompositionFaceNumber)
                          colDomainFaceBasis=>colDomainFace%BASIS
                          !Only interpolate for the first field component and when face number changes
                          IF((colFieldComp==1) .AND. (colDecompositionFaceNumber/=colPreviousFaceNo)) THEN
                            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(colDecompositionFaceNumber, &
                              & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            colPreviousFaceNo=colDecompositionFaceNumber
                          ENDIF
                          DO colLocalFaceNodeIdx=1,colDependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(colConnectedFace)
                            colFaceLocalElemNode=colDependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(colLocalFaceNodeIdx,colConnectedFace)
                            colGlobalNode=dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                              & ELEMENTS%ELEMENTS(colElementNum)%ELEMENT_NODES(colFaceLocalElemNode)
                            DO colFaceDerivative=1,colDomainFace%BASIS%NUMBER_OF_DERIVATIVES(colLocalFaceNodeIdx)
                              colDerivative=colDependentBasis% &
                                & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(colFaceDerivative,colLocalFaceNodeIdx,colConnectedFace)
                              colVersion=dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                                & ELEMENTS%ELEMENTS(colElementNum)%elementVersions(colDerivative,colFaceLocalElemNode)
                              !Find the face parameter's element parameter index 
                              colElemParameterNo=colDomainFaceBasis%ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx)
                              colDofScaleFactor=equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% PTR% &
                                & SCALE_FACTORS(colElemParameterNo,colFieldComp)
                              !Find dof associated with this particular field, component, node, derivative and version.
                              colIdx=dependentVariable%components(colFieldComp)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES( &
                                & colGlobalNode)%DERIVATIVES(colDerivative)%VERSIONS(colVersion)
                              !Evaluate the basis at the projected/connected xi
                              colPhi=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),NO_PART_DERIV,colXi,err,error)
                              
                              !Calculate the force term --idx 1 for frictionless, normal direction
                              forceTerm=coefficient*rowPhi*contactPointMetrics%normal(rowFieldComp)* & 
                                & colPhi*contactPointMetrics%normal(colFieldComp)*contactPointMetrics%contactStiffness(1) 
                              geometricTerm=0.0_DP
                              
                              !#####################################################################################################
                              IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                                !Calculate col metrics
                                phiDeriCol(1)=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                  & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),PART_DERIV_S1,colXi,err,error)
                                phiDeriCol(2)=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                  & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),PART_DERIV_S1,colXi,err,error) 
                                
                                !NCol and TCol changes at every col dof
                                DO xiIdxAlpha=1,2
                                  IF (colBodyIdx==1) THEN
                                    NCol(xiIdxAlpha)=0.0_DP
                                    TCol(xiIdxAlpha)=colPhi*contactPointMetrics%tangents(xiIdxAlpha,colFieldComp)
                                  ELSE
                                    NCol(xiIdxAlpha)=-phiDeriCol(xiIdxAlpha)*contactPointMetrics%normal(colFieldComp) 
                                    TCol(xiIdxAlpha)=-colPhi*contactPointMetrics%tangents(xiIdxAlpha,colFieldComp)
                                  ENDIF
                                ENDDO      
                                
                                !DCol varies at every col dof
                                DCol=0.0_DP
                                DO xiIdxAlpha=1,2
                                  DO xiIdxBeta=1,2
                                    DCol(xiIdxAlpha)=DCol(xiIdxAlpha)+contactPointMetrics%inverseA(xiIdxAlpha,xiIdxBeta)* &
                                      & (TCol(xiIdxBeta)+contactPointMetrics%signedGapNormal*NCol(xiIdxBeta))
                                  ENDDO !xiIdxBeta
                                ENDDO !xiIdxAlpha  
                                
                                !****************************************************************************************************
                                
                                !Calculate geometric term, see Jae's thesis equation 4.40
                                geometricTerm=0.0_DP
                                DO xiIdxGamma=1,2
                                  DO xiIdxBeta=1,2
                                    tempA=0.0_DP
                                    tempB=0.0_DP
                                    DO xiIdxAlpha=1,2
                                      tempA=tempA+kappa(xiIdxAlpha,xiIdxGamma)*DRow(xiIdxAlpha) !For row variable
                                      tempB=tempB+kappa(xiIdxAlpha,xiIdxBeta)*DCol(xiIdxAlpha) !For col variable
                                    ENDDO !xiIdxGamma
                                    geometricTerm=geometricTerm+contactPointMetrics%signedGapNormal* &
                                      & contactPointMetrics%covariantMetricTensor(xiIdxGamma,xiIdxBeta)* &
                                      & (NRow(xiIdxGamma)-tempA)*(NCol(xiIdxBeta)-tempB) + &
                                      & kappa(xiIdxBeta,xiIdxGamma)*DRow(xiIdxGamma)*DCol(xiIdxBeta)
                                  ENDDO !xiIdxBeta
                                  geometricTerm=geometricTerm-DRow(xiIdxGamma)*NCol(xiIdxGamma)-NRow(xiIdxGamma)*DCol(xiIdxGamma)
                                ENDDO !xiIdxAlpha  
                                
                                geometricTerm=geometricTerm*contactPointMetrics%contactForce
!                                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Geometric term : ",geometricTerm,ERR,ERROR,*999)
                              ENDIF !addGeometricTerm
                              !#####################################################################################################  
                              
                              !Multiply by the scale factor
  !                            matrixValue=(forceTerm+geometricTerm)*rowDofScaleFactor*colDofScaleFactor
                              matrixValue=(forceTerm+geometricTerm)

                              !Multiply by scale factors
                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(jacobian,rowIdx,colIdx,matrixValue,err,error,*999)

                            ENDDO !colFaceDerivative
                          ENDDO !colLocalFaceNodeIdx
                        ENDDO !colFieldComp
                      ENDDO !rowFaceDerivative
                    ENDDO !rowLocalFaceNodeIdx
                  ENDDO !rowFieldComp
                ENDIF !inContact
              ENDDO !globalDataPointNum

            ENDDO !subMatrix

            !Set all jacobian values to 0.0. Only for testing
            !CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(jacobian,0.0_DP,err,error,*999)

            CALL DISTRIBUTED_MATRIX_UPDATE_START(jacobian,err,error,*999)
            CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(jacobian,err,error,*999)

            !Output equations matrices and RHS vector if required
            !\todo Uncomment below after EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM is moved to equations_set_routines.
            !IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
            ! CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            !ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_JACOBIAN_CONTACT_UPDATE_STATIC_FEM
  
  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a static equations set with contact using the finite element method
  SUBROUTINE EquationsSet_JacobianRigidBodyContactUpdateStaticFEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    !TYPE(EQUATIONS_SET_TYPE), POINTER :: coupledRegionEquationsSet(2) !<A pointer to the equations set to evaluate the Jacobian for
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition  !<A pointer to the equations set to evaluate the element Jacobian for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    TYPE(InterfaceContactPointMetricsType), POINTER :: contactPointMetrics
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(BASIS_TYPE), POINTER :: rowDependentBasis,colDependentBasis,rowDomainFaceBasis,colDomainFaceBasis
    TYPE(DOMAIN_FACE_TYPE), POINTER :: rowDomainFace,colDomainFace
    TYPE(EQUATIONS_SET_TYPE), POINTER :: multipleRegionEquationsSet
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: jacobian
    !TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    !TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: residualParameterSet
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: jacobianNumber,bodyIdx,equationSetNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber
    INTEGER(INTG) :: globalDataPointNum,rowElementNum,colElementNum, &
      & rowConnectedFace,colConnectedFace,rowFieldComp, colFieldComp, &
      & rowMeshComp,colMeshComp,rowDecompositionFaceNumber,colDecompositionFaceNumber, &
      & rowLocalFaceNodeIdx,colLocalFaceNodeIdx,rowFaceLocalElemNode,colFaceLocalElemNode,rowGlobalNode,colGlobalNode, &
      & rowFaceDerivative,colFaceDerivative,rowDerivative,colDerivative,rowVersion,colVersion,rowIdx,colIdx, &
      & subMatrix,rowBodyIdx,colBodyIdx
    INTEGER(INTG) :: xiIdxAlpha,xiIdxBeta,xiIdxGamma,rowElemParameterNo,colElemParameterNo,rowPreviousFaceNo,colPreviousFaceNo
    REAL(DP) :: matrixValue,rowPhi,colPhi
    REAL(DP) :: coefficient,forceTerm,geometricTerm,tempA,tempB,rowDofScaleFactor,colDofScaleFactor
    REAL(DP) :: rowXi(2),colXi(2) !\todo generalise xi allocations for 1D,2D and 3D points connectivity
    REAL(DP) :: kappa(2,2),phiDeriRow(2),phiDeriCol(2),TRow(2),TCol(2),NRow(2),NCol(2),DRow(2),DCol(2)
  
    CALL ENTERS("EquationsSet_JacobianRigidBodyContactUpdateStaticFEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(dependentField)) THEN
        equations=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          equationsMatrices=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(equationsMatrices)) THEN
            nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
            nonlinearMapping=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING

            interfaceGlobalNumber=1
            interfaceConditionGlobalNumber=1
            interface=>EQUATIONS_SET%REGION%PARENT_REGION%INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR
            interfaceCondition=>interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
            pointsConnectivity=>interface%pointsConnectivity
            contactMetrics=>interfaceCondition%interfaceContactMetrics

            dependentVariable=>dependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)
            jacobianNumber=1
            jacobian=>nonlinearMatrices%JACOBIANS(jacobianNumber)%PTR%JACOBIAN
            dependentVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(jacobianNumber)%VARIABLE
            !Setup pointer to the equation set of the coupled bodies which are setup in thier own separate regions
            !(note that these regions has not been added to the solver equations and are merely here for convinence if needed)
            !equationSetNumber=1
            !DO bodyIdx=1,2
            !  coupledRegionEquationsSet(bodyIdx)=>EQUATIONS_SET%REGION%PARENT_REGION%SUB_REGIONS(bodyIdx)%PTR% &
            !    & EQUATIONS_SETS%EQUATIONS_SETS(equationSetNumber)%PTR
            !ENDDO
            !Since we are computing the contact term in a single region, we do not need to determine the dependent field
            !through the interface condition. We can use the dependentField pointer defined for this single region
            !dependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(interfaceMatrixIdx)%PTR% &
            ! & DEPENDENT%DEPENDENT_FIELD

            DO subMatrix=1,4
              SELECT CASE(subMatrix)
              CASE(1) !Contact subMatrix11
                rowBodyIdx=1
                colBodyIdx=1
                coefficient=1;
              CASE(2) !Contact subMatrix12
                rowBodyIdx=1
                colBodyIdx=2
                coefficient=-1;
              CASE(3) !Contact subMatrix21
                rowBodyIdx=2
                colBodyIdx=1
                coefficient=-1;
              CASE(4) !Contact subMatrix22
                rowBodyIdx=2
                colBodyIdx=2
                coefficient=1;
              END SELECT
              !Loop over each data point and find the connected element and their dofs
              DO globalDataPointNum=1,SIZE(pointsConnectivity%pointsConnectivity,1)
                IF(contactMetrics%inContact(globalDataPointNum)) THEN
                  ! Get the metric structure for this contact point
                  contactPointMetrics=>contactMetrics%contactPointMetrics(globalDataPointNum)
                  rowElementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%coupledMeshElementNumber
                  rowConnectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%elementLineFaceNumber
                  rowXi=pointsConnectivity%pointsConnectivity(globalDataPointNum,rowBodyIdx)%reducedXi
                  colElementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%coupledMeshElementNumber
                  colConnectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%elementLineFaceNumber
                  colXi=pointsConnectivity%pointsConnectivity(globalDataPointNum,colBodyIdx)%reducedXi
                  
                  !################################################################################################################
                  IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                    !Calculate quantities that do not vary w.r.t xyz or node
                    kappa=0.0_DP
                    !Calculate kappa (see Jae's cm implementation)
                    DO xiIdxAlpha=1,2
                      DO xiIdxBeta=1,2
                        kappa(xiIdxAlpha,xiIdxBeta)= & 
                          & DOT_PRODUCT(contactPointMetrics%tangentDerivatives(xiIdxAlpha,xiIdxBeta,:), &
                          & contactPointMetrics%normal(:))
                      ENDDO !xiIdxBeta
                    ENDDO !xiIdxAlpha
                  ENDIF !addGeometricTerm
                  !################################################################################################################
                  rowPreviousFaceNo=0
                  !Find the row dof 
                  DO rowFieldComp=1,3
                    rowMeshComp=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                      & COMPONENTS(rowFieldComp)%MESH_COMPONENT_NUMBER
                    rowDependentBasis=>dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR% &
                      & TOPOLOGY%ELEMENTS%ELEMENTS(rowElementNum)%BASIS
                    rowDecompositionFaceNumber=dependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(rowElementNum)%ELEMENT_FACES(rowConnectedFace)
                    rowDomainFace=>dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                      & FACES%FACES(rowDecompositionFaceNumber)
                    rowDomainFaceBasis=>rowDomainFace%BASIS
                    !Only interpolate for the first field component and when face number changes
                    IF((rowFieldComp==1) .AND. (rowDecompositionFaceNumber/=rowPreviousFaceNo)) THEN
                      CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(rowDecompositionFaceNumber, &
                        & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                      rowPreviousFaceNo=rowDecompositionFaceNumber
                    ENDIF
                    DO rowLocalFaceNodeIdx=1,rowDependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(rowConnectedFace)
                      rowFaceLocalElemNode=rowDependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(rowLocalFaceNodeIdx,rowConnectedFace)
                      rowGlobalNode=dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(rowElementNum)%ELEMENT_NODES(rowFaceLocalElemNode)
                      DO rowFaceDerivative=1,rowDomainFace%BASIS%NUMBER_OF_DERIVATIVES(rowLocalFaceNodeIdx)
                        rowDerivative=rowDependentBasis% &
                          & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(rowFaceDerivative,rowLocalFaceNodeIdx,rowConnectedFace)
                        rowVersion=dependentField%DECOMPOSITION%DOMAIN(rowMeshComp)%PTR%TOPOLOGY% &
                          & ELEMENTS%ELEMENTS(rowElementNum)%elementVersions(rowDerivative,rowFaceLocalElemNode)
                        !Find the face parameter's element parameter index 
                        rowElemParameterNo=rowDomainFaceBasis%ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx)
                        rowDofScaleFactor=equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% PTR% &
                          & SCALE_FACTORS(rowElemParameterNo,rowFieldComp)
                        !Find dof associated with this particular field, component, node, derivative and version.
                        rowIdx=dependentVariable%components(rowFieldComp)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES( &
                          & rowGlobalNode)%DERIVATIVES(rowDerivative)%VERSIONS(rowVersion)
                        !Evaluate the basis at the projected/connected xi
                        rowPhi=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                          & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),NO_PART_DERIV,rowXi,err,error)
                        
                        !###########################################################################################################
                        IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                          !Calculate row metrics
                          !\todo: generalise contact xi direction, at the moment assume in xi1 and xi2
                          phiDeriRow(1)=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                            & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),PART_DERIV_S1,rowXi,err,error)
                          phiDeriRow(2)=BASIS_EVALUATE_XI(rowDomainFaceBasis,rowDomainFaceBasis% &
                            & ELEMENT_PARAMETER_INDEX(rowFaceDerivative,rowLocalFaceNodeIdx),PART_DERIV_S2,rowXi,err,error)  
                            
                          !NRow TRow changes at every row dof
                          DO xiIdxAlpha=1,2
                            IF (rowBodyIdx==1) THEN
                              NRow(xiIdxAlpha)=0.0_DP
                              TRow(xiIdxAlpha)=rowPhi*contactPointMetrics%tangents(xiIdxAlpha,rowFieldComp)
                            ELSE
                              NRow(xiIdxAlpha)=-phiDeriRow(xiIdxAlpha)*contactPointMetrics%normal(rowFieldComp) 
                              TRow(xiIdxAlpha)=-rowPhi*contactPointMetrics%tangents(xiIdxAlpha,rowFieldComp)
                            ENDIF !rowBodyIdx
                          ENDDO !xiIdxAlpha    
                          
                          !DRow varies at every row dof
                          DRow=0.0_DP
                          DO xiIdxAlpha=1,2
                            DO xiIdxBeta=1,2
                              DRow(xiIdxAlpha)=DRow(xiIdxAlpha)+contactPointMetrics%inverseA(xiIdxAlpha,xiIdxBeta)* &
                                & (TRow(xiIdxBeta)+contactPointMetrics%signedGapNormal*NRow(xiIdxBeta))
                            ENDDO !xiIdxBeta
                          ENDDO !xiIdxAlpha   
                        ENDIF !addGeometricTerm
                        !###########################################################################################################

                        colPreviousFaceNo=0
                        DO colFieldComp=1,3
                          !Find the col 
                          colMeshComp=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & COMPONENTS(colFieldComp)%MESH_COMPONENT_NUMBER
                          colDependentBasis=>dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR% &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(colElementNum)%BASIS
                          colDecompositionFaceNumber=dependentField%DECOMPOSITION%TOPOLOGY% &
                            & ELEMENTS%ELEMENTS(colElementNum)%ELEMENT_FACES(colConnectedFace)
                          colDomainFace=>dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                            & FACES%FACES(colDecompositionFaceNumber)
                          colDomainFaceBasis=>colDomainFace%BASIS
                          !Only interpolate for the first field component and when face number changes
                          IF((colFieldComp==1) .AND. (colDecompositionFaceNumber/=colPreviousFaceNo)) THEN
                            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(colDecompositionFaceNumber, &
                              & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            colPreviousFaceNo=colDecompositionFaceNumber
                          ENDIF
                          DO colLocalFaceNodeIdx=1,colDependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(colConnectedFace)
                            colFaceLocalElemNode=colDependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(colLocalFaceNodeIdx,colConnectedFace)
                            colGlobalNode=dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                              & ELEMENTS%ELEMENTS(colElementNum)%ELEMENT_NODES(colFaceLocalElemNode)
                            DO colFaceDerivative=1,colDomainFace%BASIS%NUMBER_OF_DERIVATIVES(colLocalFaceNodeIdx)
                              colDerivative=colDependentBasis% &
                                & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(colFaceDerivative,colLocalFaceNodeIdx,colConnectedFace)
                              colVersion=dependentField%DECOMPOSITION%DOMAIN(colMeshComp)%PTR%TOPOLOGY% &
                                & ELEMENTS%ELEMENTS(colElementNum)%elementVersions(colDerivative,colFaceLocalElemNode)
                              !Find the face parameter's element parameter index 
                              colElemParameterNo=colDomainFaceBasis%ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx)
                              colDofScaleFactor=equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% PTR% &
                                & SCALE_FACTORS(colElemParameterNo,colFieldComp)
                              !Find dof associated with this particular field, component, node, derivative and version.
                              colIdx=dependentVariable%components(colFieldComp)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES( &
                                & colGlobalNode)%DERIVATIVES(colDerivative)%VERSIONS(colVersion)
                              !Evaluate the basis at the projected/connected xi
                              colPhi=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),NO_PART_DERIV,colXi,err,error)
                              
                              !Calculate the force term --idx 1 for frictionless, normal direction
                              forceTerm=coefficient*rowPhi*contactPointMetrics%normal(rowFieldComp)* & 
                                & colPhi*contactPointMetrics%normal(colFieldComp)*contactPointMetrics%contactStiffness(1) 
                              geometricTerm=0.0_DP
                              
                              !#####################################################################################################
                              IF(contactMetrics%addGeometricTerm) THEN !Only calculate if geometric term is included
                                !Calculate col metrics
                                phiDeriCol(1)=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                  & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),PART_DERIV_S1,colXi,err,error)
                                phiDeriCol(2)=BASIS_EVALUATE_XI(colDomainFaceBasis,colDomainFaceBasis% &
                                  & ELEMENT_PARAMETER_INDEX(colFaceDerivative,colLocalFaceNodeIdx),PART_DERIV_S1,colXi,err,error) 
                                
                                !NCol and TCol changes at every col dof
                                DO xiIdxAlpha=1,2
                                  IF (colBodyIdx==1) THEN
                                    NCol(xiIdxAlpha)=0.0_DP
                                    TCol(xiIdxAlpha)=colPhi*contactPointMetrics%tangents(xiIdxAlpha,colFieldComp)
                                  ELSE
                                    NCol(xiIdxAlpha)=-phiDeriCol(xiIdxAlpha)*contactPointMetrics%normal(colFieldComp) 
                                    TCol(xiIdxAlpha)=-colPhi*contactPointMetrics%tangents(xiIdxAlpha,colFieldComp)
                                  ENDIF
                                ENDDO      
                                
                                !DCol varies at every col dof
                                DCol=0.0_DP
                                DO xiIdxAlpha=1,2
                                  DO xiIdxBeta=1,2
                                    DCol(xiIdxAlpha)=DCol(xiIdxAlpha)+contactPointMetrics%inverseA(xiIdxAlpha,xiIdxBeta)* &
                                      & (TCol(xiIdxBeta)+contactPointMetrics%signedGapNormal*NCol(xiIdxBeta))
                                  ENDDO !xiIdxBeta
                                ENDDO !xiIdxAlpha  
                                
                                !****************************************************************************************************
                                
                                !Calculate geometric term, see Jae's thesis equation 4.40
                                geometricTerm=0.0_DP
                                DO xiIdxGamma=1,2
                                  DO xiIdxBeta=1,2
                                    tempA=0.0_DP
                                    tempB=0.0_DP
                                    DO xiIdxAlpha=1,2
                                      tempA=tempA+kappa(xiIdxAlpha,xiIdxGamma)*DRow(xiIdxAlpha) !For row variable
                                      tempB=tempB+kappa(xiIdxAlpha,xiIdxBeta)*DCol(xiIdxAlpha) !For col variable
                                    ENDDO !xiIdxGamma
                                    geometricTerm=geometricTerm+contactPointMetrics%signedGapNormal* &
                                      & contactPointMetrics%covariantMetricTensor(xiIdxGamma,xiIdxBeta)* &
                                      & (NRow(xiIdxGamma)-tempA)*(NCol(xiIdxBeta)-tempB) + &
                                      & kappa(xiIdxBeta,xiIdxGamma)*DRow(xiIdxGamma)*DCol(xiIdxBeta)
                                  ENDDO !xiIdxBeta
                                  geometricTerm=geometricTerm-DRow(xiIdxGamma)*NCol(xiIdxGamma)-NRow(xiIdxGamma)*DCol(xiIdxGamma)
                                ENDDO !xiIdxAlpha  
                                
                                geometricTerm=geometricTerm*contactPointMetrics%contactForce
!                                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Geometric term : ",geometricTerm,ERR,ERROR,*999)
                              ENDIF !addGeometricTerm
                              !#####################################################################################################  
                              
                              !Multiply by the scale factor
  !                            matrixValue=(forceTerm+geometricTerm)*rowDofScaleFactor*colDofScaleFactor
                              matrixValue=(forceTerm+geometricTerm)

                              !Multiply by scale factors
                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(jacobian,rowIdx,colIdx,matrixValue,err,error,*999)

                            ENDDO !colFaceDerivative
                          ENDDO !colLocalFaceNodeIdx
                        ENDDO !colFieldComp
                      ENDDO !rowFaceDerivative
                    ENDDO !rowLocalFaceNodeIdx
                  ENDDO !rowFieldComp
                ENDIF !inContact
              ENDDO !globalDataPointNum

            ENDDO !subMatrix

            !Set all jacobian values to 0.0. Only for testing
            !CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(jacobian,0.0_DP,err,error,*999)

            CALL DISTRIBUTED_MATRIX_UPDATE_START(jacobian,err,error,*999)
            CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(jacobian,err,error,*999)

            !Output equations matrices and RHS vector if required
            !\todo Uncomment below after EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM is moved to equations_set_routines.
            !IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
            ! CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            !ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("EquationsSet_JacobianRigidBodyContactUpdateStaticFEM")
    RETURN
999 CALL ERRORS("EquationsSet_JacobianRigidBodyContactUpdateStaticFEM",ERR,ERROR)
    CALL EXITS("EquationsSet_JacobianRigidBodyContactUpdateStaticFEM")
    RETURN 1
  END SUBROUTINE EquationsSet_JacobianRigidBodyContactUpdateStaticFEM

  !
  !================================================================================================================================
  ! 

  !>Evaluates the residual for a nonlinear problem solver.
  SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,solver_matrix_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,equationSetRigidNodal
    TYPE(SOLVER_TYPE), POINTER :: CELLML_SOLVER,LINKING_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: residualVectorDis
    REAL(DP), POINTER :: residualVector(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    !\todo Temporarily added the variables below to allow the interface condition to be used in the single region contact problem
    !to be manually specified. Need to Generalise.
    INTEGER(INTG) :: equationsSetGlobalNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber,iterationNumber, &
      & rigidBodyRegionNumber
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_TYPE), POINTER :: interface
    
    NULLIFY(CELLML_SOLVER)
    NULLIFY(LINKING_SOLVER)

    CALL ENTERS("PROBLEM_SOLVER_RESIDUAL_EVALUATE",ERR,ERROR,*999)

!    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Residual evaluation****************",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            IF(SOLVER%OUTPUT_TYPE>=SOLVER_MATRIX_OUTPUT) THEN
              SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Solver vector values:",ERR,ERROR,*999)
                DO solver_matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                  SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                  IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Solver matrix is not associated for solver matrix index "// &
                      & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !solver_matrix_idx
              ELSE
                CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
              !Check if the nonlinear solver is linked to a dynamic solver 
              LINKING_SOLVER=>SOLVER%LINKING_SOLVER
              IF(ASSOCIATED(LINKING_SOLVER)) THEN
                IF(LINKING_SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  !Update the field values from the dynamic factor*current solver values AND add in predicted displacements
                  CALL SOLVER_VARIABLES_DYNAMIC_NONLINEAR_UPDATE(SOLVER,ERR,ERROR,*999)
                  !Caculate the strain field for an CellML evaluator solver
                  CALL PROBLEM_PRE_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*999)
                  !check for a linked CellML solver 
                  CELLML_SOLVER=>SOLVER%NONLINEAR_SOLVER%NEWTON_SOLVER%CELLML_EVALUATOR_SOLVER
                  IF(ASSOCIATED(CELLML_SOLVER)) THEN
                    CALL SOLVER_SOLVE(CELLML_SOLVER,ERR,ERROR,*999)
                  ENDIF
                  !Calculate the residual for each element (M, C, K and g)
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    SELECT CASE(EQUATIONS_SET%EQUATIONS%LINEARITY)
                    CASE(EQUATIONS_LINEAR)
                      !Assemble the equations for linear equations
                      CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
                    CASE(EQUATIONS_NONLINEAR)
                      !Evaluate the residual for nonlinear equations
                      CALL EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                    END SELECT
                  ENDDO !equations_set_idx
                  !Assemble the final solver residual.
                  CALL SOLVER_MATRICES_DYNAMIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Solver equations linking solver mapping is not dynamic.",ERR,ERROR,*999)
                END IF
              ELSE
                !Perform as normal nonlinear solver
                !Copy the current solution vector to the dependent field
                CALL SOLVER_VARIABLES_FIELD_UPDATE(SOLVER,ERR,ERROR,*999)
                !Caculate the strain field for an CellML evaluator solver
                CALL PROBLEM_PRE_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*999)
                !check for a linked CellML solver 
                CELLML_SOLVER=>SOLVER%NONLINEAR_SOLVER%NEWTON_SOLVER%CELLML_EVALUATOR_SOLVER
                IF(ASSOCIATED(CELLML_SOLVER)) THEN
                  CALL SOLVER_SOLVE(CELLML_SOLVER,ERR,ERROR,*999)
                ENDIF
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Residual evaluation******************",ERR,ERROR,*999)
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  SELECT CASE(EQUATIONS_SET%EQUATIONS%LINEARITY)
                  CASE(EQUATIONS_LINEAR)
                    !Assemble the equations for linear equations
                    CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
                  CASE(EQUATIONS_NONLINEAR)
                    !Evaluate the residual for nonlinear equations
                    CALL EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                  END SELECT
                ENDDO !equations_set_idx

                !\todo Temporarily comment out the looping through of interface conditions added to the solver as these are not
                !present in the single region contact problem. Needs to generalised.
                !DO interfaceConditionIdx=1,solverMapping%NUMBER_OF_INTERFACE_CONDITIONS
                  !interfaceCondition=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
                IF(SOLVER%SOLVERS%CONTROL_LOOP%PROBLEM%TYPE==PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)THEN
                  equationsSetGlobalNumber=1
                  rigidBodyRegionNumber=2
                  interfaceGlobalNumber=1
                  interfaceConditionGlobalNumber=1
                  interfaceCondition=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetGlobalNumber)%PTR%REGION%PARENT_REGION% &
                    & INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR%INTERFACE_CONDITIONS% &
                    & INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
                    
!                    SOLVER_MAPPING%equations_sets(1)%ptr%region%parent_region%sub_regions(2)%ptr%equations_sets%equations_sets(1)%ptr
                    
                  equationSetRigidNodal=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetGlobalNumber)%PTR%REGION% &
                    & PARENT_REGION%SUB_REGIONS(rigidBodyRegionNumber)%PTR%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(interfaceCondition)) THEN
                    IF(interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR .OR. &
                        & interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_OPERATOR) THEN !Only reproject for contact operator
                      IF(interfaceCondition%integrationType==INTERFACE_CONDITION_DATA_POINTS_INTEGRATION) THEN !Only reproject for data point interpolated field
                        interface=>interfaceCondition%INTERFACE
                        IF(ASSOCIATED(interface)) THEN
                          CALL PETSC_SNESGETITERATIONNUMBER(SOLVER%NONLINEAR_SOLVER%NEWTON_SOLVER%LINESEARCH_SOLVER%SNES, &
                            & iterationNumber,ERR,ERROR,*999)
                            IF(SOLVER%SOLVERS%CONTROL_LOOP%PROBLEM%SUBTYPE==PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE .OR. &
                                & SOLVER%SOLVERS%CONTROL_LOOP%PROBLEM%SUBTYPE==PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE .OR.  &
                                & iterationNumber==0) THEN
!                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************  Reproject! ***************",ERR,ERROR,*999)
                              CALL RigidBody_ApplyTransformation(EQUATIONS_SET,equationSetRigidNodal,err,error,*999)
                              CALL InterfacePointsConnectivity_DataReprojection(interface,interfaceCondition,err,error,*999)
                              
                              
                            ENDIF
                          ! iteration+1 since iterationNumber is counting the iterations completed
                          CALL FrictionlessContact_contactMetricsCalculate(interfaceCondition,iterationNumber+1,err,error,*999)
  !                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************  Contact residual! ***********",ERR,ERROR,*999)
                          !\todo: generalise when LHS mapping is in place
                          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_MULTI_PHYSICS_CLASS .AND.  &
                             & EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_RIGID_BODY_TYPE) THEN !Rigid-deformable contact
                            !Modify residual for rigid-deformable contact
                            CALL EquationsSet_ResidualRigidBodyContactUpdateStaticFEM(SOLVER_MAPPING% &
                              & EQUATIONS_SETS(equationsSetGlobalNumber)%PTR,ERR,ERROR,*999)
                          ELSEIF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN !Deformable-deformable contact
                            !Modify residual for deformable bodies contact
                            CALL EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM(SOLVER_MAPPING% &
                              & EQUATIONS_SETS(equationsSetGlobalNumber)%PTR,ERR,ERROR,*999)
                          ENDIF
  !                        CALL SolverEquations_ResidualVectorGet(SOLVER_EQUATIONS,residualVectorDis,err,error,*999)
                          CALL DISTRIBUTED_VECTOR_DATA_GET(SOLVER_EQUATIONS%solver_matrices%residual,residualVector,err,error,*999)
                          
                          !\todo Temporarily commented out INTERFACE_CONDITION_ASSEMBLE as the interface matrices are not
                          ! required for the single region contact problem. Needs to generalised.
                          !CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
                        ELSE
                          CALL FLAG_ERROR("Interface is not associated for nonlinear solver equations mapping.", &
                            & err,error,*999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface condition is not associated for nonlinear solver equations mapping.", &
                      & err,error,*999)
                  ENDIF 
                ENDIF !problemType==contact
              !ENDDO !interfaceConditionIdx

                !\todo Temporarily comment out the looping through of interface conditions added to the solver as these are not
                !present in the single region contact problem. Needs to generalised.
                !DO interfaceConditionIdx=1,solverMapping%NUMBER_OF_INTERFACE_CONDITIONS
                  !interfaceCondition=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
                !Note that the linear interface matrices are not required to be updated since these matrices do not change
                !Update interface matrices
!                DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
!                  interfaceCondition=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
!                  !Assemble the interface condition for the Jacobian LHS
!                  CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
!                ENDDO
                !Assemble the solver matrices
                CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,ERR,ERROR,*999)
              END IF
            ELSE
               CALL FLAG_ERROR("Solver equations solver type is not associated.",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations mapping is not associated.",ERR,ERROR,*999)
        ENDIF
        CALL PROBLEM_POST_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Updates the equation set residual for a static equations set which includes contact using the finite element method
  SUBROUTINE EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition  !<A pointer to the equations set to evaluate the element Jacobian for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(FIELD_TYPE), POINTER :: dependentField,penaltyField
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,domainFaceBasis
    TYPE(DOMAIN_FACE_TYPE), POINTER :: domainFace
!    TYPE(EQUATIONS_SET_TYPE), POINTER :: multipleRegionEquationsSet
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: residualParameterSet
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    TYPE(InterfaceContactPointMetricsType), POINTER :: contactPointMetrics
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: bodyIdx,equationSetNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber
    INTEGER(INTG) :: globalDataPointNum,elementNum,connectedFace,fieldComponent,meshComp, &
      & decompositionFaceNumber,localFaceNodeIdx,faceLocalElemNode,globalNode,faceDerivative,derivative,versionNumber, &
      & residualVariableIdx,dofIdx
    INTEGER(INTG) :: coefficient,previousFaceNo,elemParameterNo
    REAL(DP) :: residualValue,phi,contactForce
    REAL(DP) :: xi(2) !\todo generalise xi allocations for 1D,2D and 3D points connectivity

    CALL ENTERS("EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(dependentField)) THEN
        equations=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          equationsMatrices=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(equationsMatrices)) THEN
            nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
            !nonlinearResidual=>nonlinearMatrices%RESIDUAL
            nonlinearMapping=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING

            interfaceGlobalNumber=1
            interfaceConditionGlobalNumber=1
            interface=>EQUATIONS_SET%REGION%PARENT_REGION%INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR
            interfaceCondition=>interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
            pointsConnectivity=>interface%pointsConnectivity
            
            contactMetrics=>interfaceCondition%interfaceContactMetrics

            residualVariableIdx=1
            residualVariable=>nonlinearMapping%RESIDUAL_VARIABLES(residualVariableIdx)%PTR
            IF(ASSOCIATED(residualVariable)) THEN
              residualValue=0.0_DP
              !Loop over each coupled body and add the contact contribution associated with each contact point
              DO bodyIdx=1,2
                !Setup pointer to the equation set of the coupled bodies which are setup in thier own separate regions
                !(note that these regions has not been added to the solver equations and are merely here for convinence if needed)
                equationSetNumber=1
!                multipleRegionEquationsSet=>EQUATIONS_SET%REGION%PARENT_REGION%SUB_REGIONS(bodyIdx)%PTR%EQUATIONS_SETS% &
!                  & EQUATIONS_SETS(equationSetNumber)%PTR
                
                ! Residual is +ve for body 1 and -ve for body 2  
                SELECT CASE(bodyIdx)
                CASE(1)
                  coefficient=1;
                CASE(2)
                  coefficient=-1;
                CASE DEFAULT
                  CALL FLAG_ERROR("Contact for 3 or more bodies is not implemented",err,error,*999)
                END SELECT 

                !Since we are computing the contact term in a single region, we do not need to determine the dependent field
                !through the interface condition. We can simply use the dependentField pointer defined above
                !dependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(interfaceMatrixIdx)%PTR% &
                ! & DEPENDENT%DEPENDENT_FIELD

                !Loop over each data point and find the connected element and their dofs
                previousFaceNo=0
                DO globalDataPointNum=1,SIZE(pointsConnectivity%pointsConnectivity,1)
                  IF(contactMetrics%inContact(globalDataPointNum)) THEN
                    elementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%coupledMeshElementNumber
                    connectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%elementLineFaceNumber
                    xi=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%reducedXi
                    contactPointMetrics=>contactMetrics%contactPointMetrics(globalDataPointNum)
                    DO fieldComponent=1,3
                      meshComp=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                        & COMPONENTS(fieldComponent)%MESH_COMPONENT_NUMBER
                      dependentBasis=>dependentField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNum)%BASIS
                      decompositionFaceNumber=dependentField%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(elementNum)%ELEMENT_FACES(connectedFace)
                      domainFace=>dependentField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY%FACES%FACES(decompositionFaceNumber)
                      domainFaceBasis=>domainFace%BASIS
                      !Only interpolate for the first field component and when face number changes
                      IF((fieldComponent==1) .AND. (decompositionFaceNumber/=previousFaceNo)) THEN
                        CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(decompositionFaceNumber, &
                          & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        previousFaceNo=decompositionFaceNumber
                      ENDIF
                      DO localFaceNodeIdx=1,dependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        faceLocalElemNode=dependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        globalNode=dependentField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY% &
                          & ELEMENTS%ELEMENTS(elementNum)%ELEMENT_NODES(faceLocalElemNode)
                        DO faceDerivative=1,domainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                          derivative=dependentBasis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(faceDerivative,localFaceNodeIdx,connectedFace)
                          versionNumber=dependentField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY% &
                            & ELEMENTS%ELEMENTS(elementNum)%elementVersions(derivative,faceLocalElemNode)

                          !Evaluate the basis at the projected/connected xi
                          phi=BASIS_EVALUATE_XI(domainFaceBasis,domainFaceBasis% &
                            & ELEMENT_PARAMETER_INDEX(faceDerivative,localFaceNodeIdx),NO_PART_DERIV,xi,err,error)

                          !Find dof associated with this particular field, component, node, derivative and version.
                          dofIdx=residualVariable%components(fieldComponent)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(globalNode)%DERIVATIVES(derivative)%VERSIONS(versionNumber)
                          ! See Jae's thesis equation 4.34
                          residualValue=-coefficient*phi*contactPointMetrics%normal(fieldComponent)*contactPointMetrics%contactForce
                          !Get the face parameter index in the element
                          elemParameterNo=domainFace%BASIS%ELEMENT_PARAMETER_INDEX(faceDerivative,localFaceNodeIdx)
                          !Multiply the contribution by scale factor
                          
  !                        residualValue=residualValue*equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
  !                          & PTR%SCALE_FACTORS(elemParameterNo,fieldComponent)
                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(nonlinearMatrices%RESIDUAL,dofIdx,residualValue,err,error,*999)
                        ENDDO !faceDerivative
                      ENDDO !localFaceNodeIdx
                    ENDDO !fieldComponent
                  ENDIF !inContact
                ENDDO !globalDataPointNum
              ENDDO !bodyIdx

              !Update the residual parameter set
              residualParameterSet=>residualVariable%PARAMETER_SETS%SET_TYPE(FIELD_RESIDUAL_SET_TYPE)%PTR
              IF(ASSOCIATED(residualParameterSet)) THEN
                !Residual parameter set exists
                !Copy the residual vector to the residuals parameter set.
                CALL DISTRIBUTED_VECTOR_COPY(nonlinearMatrices%RESIDUAL,residualParameterSet%PARAMETERS,1.0_DP, &
                  & err,error,*999)
              ENDIF
            ELSE
              localError="Nonlinear mapping residual variable for residual variable index "// &
                & TRIM(NUMBER_TO_VSTRING(residualVariableIdx,"*",err,error))//" is not associated."
              CALL FLAG_ERROR(localError,err,error,*999)
            ENDIF

            !Output equations matrices and RHS vector if required
            !\todo Uncomment below after EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM is moved to equations_set_routines.
            !IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
            ! CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            !ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM",err,error)
    CALL EXITS("EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM
  
  !
  !================================================================================================================================
  !

  !>Updates the equation set residual for a static equations set which includes rigid body contact using the finite element method
  SUBROUTINE EquationsSet_ResidualRigidBodyContactUpdateStaticFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: defDepField,rigidGeoField,LagrangeField
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition  !<A pointer to the equations set to evaluate the element Jacobian for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: residualParameterSet
    TYPE(InterfaceContactPointMetricsType), POINTER :: contactPointMetrics
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,domainFaceBasis
    TYPE(DOMAIN_FACE_TYPE), POINTER :: domainFace
    INTEGER(INTG) :: interfaceGlobalNumber,interfaceConditionGlobalNumber,rigidBodyRegionNumber,residualVariableIdx
    INTEGER(INTG) :: coefficient,elementNum,connectedFace,previousFaceNo,meshComp,decompositionFaceNumber,faceLocalElemNode, &
      & globalNode,versionNumber,elemParameterNo
    INTEGER(INTG) :: bodyIdx,globalDataPointNum,fieldComponent,dofIdx,localFaceNodeIdx,faceDerivative,derivative,rigidBodyDofCompIdx
    REAL(DP) :: residualValue,phi
    REAL(DP) :: xi(2), rigidBodyMatrix(3,6),contactPtPosition(3) !\todo generalise xi allocations for 1D,2D and 3D points connectivity
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("EquationsSet_ResidualRigidBodyContactUpdateStaticFEM",ERR,ERROR,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      defDepField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(defDepField)) THEN
        equations=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          equationsMatrices=>equations%EQUATIONS_MATRICES
          IF(ASSOCIATED(equationsMatrices)) THEN
            nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
            !nonlinearResidual=>nonlinearMatrices%RESIDUAL
            nonlinearMapping=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING
            
            interfaceGlobalNumber=1
            interfaceConditionGlobalNumber=1
            interface=>equationsSet%REGION%PARENT_REGION%INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR
            interfaceCondition=>interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
            pointsConnectivity=>interface%pointsConnectivity
            contactMetrics=>interfaceCondition%interfaceContactMetrics
            
            rigidBodyRegionNumber=2
            rigidGeoField=>equationsSet%REGION%PARENT_REGION%SUB_REGIONS(rigidBodyRegionNumber)%PTR%EQUATIONS_SETS% &
              & EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
            LagrangeField=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD
            
            residualVariableIdx=1
            residualVariable=>nonlinearMapping%RESIDUAL_VARIABLES(residualVariableIdx)%PTR
            IF(ASSOCIATED(residualVariable)) THEN
              previousFaceNo=0
              DO globalDataPointNum=1,SIZE(pointsConnectivity%pointsConnectivity,1)
                contactPointMetrics=>contactMetrics%contactPointMetrics(globalDataPointNum)
                IF(contactMetrics%inContact(globalDataPointNum)) THEN
                
                  !###########################################################################################################
                  !                                             body 1 - deformable
                  coefficient=1 !Residual is +ve for body 1
                  bodyIdx=1
                  residualValue=0.0_DP
                  elementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%coupledMeshElementNumber
                  connectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%elementLineFaceNumber
                  xi=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%reducedXi
                  DO fieldComponent=1,3
                    meshComp=defDepField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                      & COMPONENTS(fieldComponent)%MESH_COMPONENT_NUMBER
                    dependentBasis=>defDepField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNum)%BASIS
                    decompositionFaceNumber=defDepField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(elementNum)%ELEMENT_FACES(connectedFace)
                    domainFace=>defDepField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY%FACES%FACES(decompositionFaceNumber)
                    domainFaceBasis=>domainFace%BASIS
                    
                    !Only interpolate for the first field component and when face number changes
                    IF((fieldComponent==1) .AND. (decompositionFaceNumber/=previousFaceNo)) THEN
                      CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_FACE_GET(decompositionFaceNumber, &
                        & equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                      previousFaceNo=decompositionFaceNumber
                    ENDIF
                    !\todo: connectedFace is the local face no the global face, need to check if face is on the current domain?
                    DO localFaceNodeIdx=1,dependentBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                      faceLocalElemNode=dependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                      globalNode=defDepField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(elementNum)%ELEMENT_NODES(faceLocalElemNode)
                      DO faceDerivative=1,domainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                        derivative=dependentBasis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(faceDerivative,localFaceNodeIdx,connectedFace)
                        versionNumber=defDepField%DECOMPOSITION%DOMAIN(meshComp)%PTR%TOPOLOGY% &
                          & ELEMENTS%ELEMENTS(elementNum)%elementVersions(derivative,faceLocalElemNode)
                        !Evaluate the basis at the projected/connected xi
                        phi=BASIS_EVALUATE_XI(domainFaceBasis,domainFaceBasis% &
                          & ELEMENT_PARAMETER_INDEX(faceDerivative,localFaceNodeIdx),NO_PART_DERIV,xi,err,error)
                        !Find dof associated with this particular field, component, node, derivative and version.
                        dofIdx=residualVariable%components(fieldComponent)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(globalNode)%DERIVATIVES(derivative)%VERSIONS(versionNumber)
                        ! See Jae's thesis equation 4.34
                        residualValue=-coefficient*phi*contactPointMetrics%normal(fieldComponent)*contactPointMetrics%contactForce
                        !Get the face parameter index in the element
                        elemParameterNo=domainFace%BASIS%ELEMENT_PARAMETER_INDEX(faceDerivative,localFaceNodeIdx)
                        !Multiply the contribution by scale factor
!                        residualValue=residualValue*equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
!                          & PTR%SCALE_FACTORS(elemParameterNo,fieldComponent)
                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(nonlinearMatrices%RESIDUAL,dofIdx,residualValue,err,error,*999)
                      ENDDO !faceDerivative
                    ENDDO !localFaceNodeIdx  
                  ENDDO !fieldComponent
                  
                  !###########################################################################################################
                  !                                             body 2 - rigid body
                  coefficient=-1 !Residual is -ve for body 2
                  residualValue=0.0_DP
                  DO fieldComponent=1,3
                    CALL Field_ParameterSetGetDataPoint(LagrangeField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & globalDataPointNum,fieldComponent,contactPtPosition(fieldComponent),err,error,*999)
                  ENDDO !fieldComponent
                  rigidBodyMatrix=0.0_DP
                  rigidBodyMatrix(1,1)=1.0_DP
                  rigidBodyMatrix(2,2)=1.0_DP
                  rigidBodyMatrix(3,3)=1.0_DP
                  rigidBodyMatrix(1,5)=-contactPtPosition(3)
                  rigidBodyMatrix(1,6)=contactPtPosition(2)
                  rigidBodyMatrix(2,4)=contactPtPosition(3)
                  rigidBodyMatrix(2,6)=-contactPtPosition(1)
                  rigidBodyMatrix(3,4)=-contactPtPosition(2)
                  rigidBodyMatrix(3,5)=contactPtPosition(1)
                  !\todo: generalise the offset for deformable body components, i.e. 4
                  DO fieldComponent=1,3
                    DO rigidBodyDofCompIdx=1,6
                      dofIdx=residualVariable%components(rigidBodyDofCompIdx+4)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP 
                      residualValue=-coefficient*rigidBodyMatrix(fieldComponent,rigidBodyDofCompIdx)* &
                        & contactPointMetrics%normal(fieldComponent)*contactPointMetrics%contactForce
                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(nonlinearMatrices%RESIDUAL,dofIdx, &
                        & residualValue,err,error,*999)
                    ENDDO !dofIdx
                  ENDDO !fieldComponent
                ENDIF !inContact
              ENDDO !globalDataPointNum
              
              !Update the residual parameter set
              residualParameterSet=>residualVariable%PARAMETER_SETS%SET_TYPE(FIELD_RESIDUAL_SET_TYPE)%PTR
              IF(ASSOCIATED(residualParameterSet)) THEN
                !Residual parameter set exists
                !Copy the residual vector to the residuals parameter set.
                CALL DISTRIBUTED_VECTOR_COPY(nonlinearMatrices%RESIDUAL,residualParameterSet%PARAMETERS,1.0_DP, &
                  & err,error,*999)
              ENDIF
            ELSE
              localError="Nonlinear mapping residual variable for residual variable index "// &
                & TRIM(NUMBER_TO_VSTRING(residualVariableIdx,"*",err,error))//" is not associated."
              CALL FLAG_ERROR(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Deformable dependent field is not associated",err,error,*999)
      ENDIF
    ELSE
      
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("EquationsSet_ResidualRigidBodyContactUpdateStaticFEM")
    RETURN
999 CALL ERRORS("EquationsSet_ResidualRigidBodyContactUpdateStaticFEM",err,error)
    CALL EXITS("EquationsSet_ResidualRigidBodyContactUpdateStaticFEM")
    RETURN 1
  END SUBROUTINE EquationsSet_ResidualRigidBodyContactUpdateStaticFEM

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for the solver
  SUBROUTINE PROBLEM_PRE_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to pre-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_PRE_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                    SELECT CASE(EQUATIONS%LINEARITY)
                    CASE(EQUATIONS_LINEAR)            
                      CALL FLAG_ERROR("Can not pre-evaluate a residual for linear equations.",ERR,ERROR,*999)
                    CASE(EQUATIONS_NONLINEAR)
                      SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC,EQUATIONS_FIRST_ORDER_DYNAMIC) ! quasistatic handled like static
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                          SELECT CASE(EQUATIONS_SET%CLASS)
                          CASE(EQUATIONS_SET_ELASTICITY_CLASS)
                            CALL ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
                            !Pre residual evaluate not used
                          CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
                            !Pre residual evaluate not used
                          CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
                            !Pre residual evaluate not used
                          CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
                            !Pre residual evaluate not used
                          CASE(EQUATIONS_SET_MODAL_CLASS)
                            !Pre residual evaluate not used
                          CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
                            !Pre residual evaluate not used
                          CASE DEFAULT
                            LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))// &
                              & " is not valid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT !EQUATIONS_SET%CLASS
                        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                          SELECT CASE(EQUATIONS_SET%CLASS)
                          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
                            !Pre residual evaluate not used
                          CASE DEFAULT
                            LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))// &
                              & " is not valid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT !EQUATIONS_SET%CLASS
                        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The equations set solution method  of "// &
                            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                            & " is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT !EQUATIONS_SET%SOLUTION_METHOD
                      CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(EQUATIONS_TIME_STEPPING)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The equations set time dependence type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE(EQUATIONS_NONLINEAR_BCS)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The equations linearity of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Equations have not been finished.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                ENDIF      
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_set_idx
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
       
    CALL EXITS("PROBLEM_PRE_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_PRE_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_PRE_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_PRE_RESIDUAL_EVALUATE
     
  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for the solver
  SUBROUTINE PROBLEM_POST_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to post-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_POST_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                    SELECT CASE(EQUATIONS%LINEARITY)
                    CASE(EQUATIONS_LINEAR)            
                      CALL FLAG_ERROR("Can not post-evaluate a residual for linear equations.",ERR,ERROR,*999)
                    CASE(EQUATIONS_NONLINEAR)
                      SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC,EQUATIONS_FIRST_ORDER_DYNAMIC) ! quasistatic handled like static
                        SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                          SELECT CASE(EQUATIONS_SET%CLASS)
                          CASE(EQUATIONS_SET_ELASTICITY_CLASS)
                            CALL ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
                            !Post residual evaluate not used
                          CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
                            !Post residual evaluate not used
                          CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
                            !Post residual evaluate not used
                          CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
                            !Post residual evaluate not used
                          CASE(EQUATIONS_SET_MODAL_CLASS)
                            !Post residual evaluate not used
                          CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
                            !Post residual evaluate not used
                          CASE DEFAULT
                            LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))// &
                              & " is not valid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT !EQUATIONS_SET%CLASS
                        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                          SELECT CASE(EQUATIONS_SET%CLASS)
                          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
                            !Post residual evaluate not used
                          CASE DEFAULT
                            LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))// &
                              & " is not valid with the nodal solution method."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT !EQUATIONS_SET%CLASS
                        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The equations set solution method  of "// &
                            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                            & " is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT !EQUATIONS_SET%SOLUTION_METHOD
                      CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(EQUATIONS_TIME_STEPPING)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The equations set time dependence type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE(EQUATIONS_NONLINEAR_BCS)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The equations linearity of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Equations have not been finished.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                ENDIF      
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_set_idx
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
       
    CALL EXITS("PROBLEM_POST_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_POST_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_POST_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_POST_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Finish the creation of solvers for a problem. \see OPENCMISS::CMISSProblemSolversCreateFinish
  SUBROUTINE PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the creation of the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO
     
    CALL ENTERS("PROBLEM_SOLVERS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN              
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_SOLVERS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_FINISH_ACTION
      !Finish the problem specific solvers setup.
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVERS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a solvers for the problem. \see OPENCMISS::CMISSProblemSolversCreateStart
  SUBROUTINE PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to create the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_SOLVERS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN    
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_SOLVERS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
      !Start the problem specific solvers setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVERS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Solves a problem. \see OPENCMISS::CMISSProblemSolve
  SUBROUTINE PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    
    CALL ENTERS("PROBLEM_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVE

  !
  !================================================================================================================================
  !

  !> Apply the load increment for each equations_set associated with solver.
  SUBROUTINE PROBLEM_SOLVER_LOAD_INCREMENT_APPLY(SOLVER_EQUATIONS,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG) :: equations_set_idx

    CALL ENTERS("PROBLEM_SOLVER_LOAD_INCREMENT_APPLY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
      IF(ASSOCIATED(SOLVER_MAPPING)) THEN
        !Make sure the equations sets are up to date
        DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
          EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
          CALL EQUATIONS_SET_LOAD_INCREMENT_APPLY(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ITERATION_NUMBER, &
            & MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
        ENDDO !equations_set_idx
      ELSE
        CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_LOAD_INCREMENT_APPLY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_LOAD_INCREMENT_APPLY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_LOAD_INCREMENT_APPLY")
    RETURN 1

  END SUBROUTINE PROBLEM_SOLVER_LOAD_INCREMENT_APPLY

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE PROBLEM_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        !For all time loops, update the previous values from the current values
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE(CONTROL_LOOP,ERR,ERROR,*999)
        ENDIF
        SELECT CASE(CONTROL_LOOP%PROBLEM%CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL ELASTICITY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_BIOELECTRICS_CLASS)
          !do nothing
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          !do nothing
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          !do nothing
        CASE(PROBLEM_FITTING_CLASS)
          !do nothing
        CASE(PROBLEM_MODAL_CLASS)
          !do nothing
        CASE(PROBLEM_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%CLASS,"*",ERR,ERROR))//" &
            & is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("PROBLEM_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE PROBLEM_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        SELECT CASE(CONTROL_LOOP%PROBLEM%CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL Elasticity_ControlLoopPostLoop(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_BIOELECTRICS_CLASS)
          CALL BIOELECTRIC_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          !Do nothing
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          !Do nothing
        CASE(PROBLEM_FITTING_CLASS)
          !Do nothing
        CASE(PROBLEM_MODAL_CLASS)
          !Do nothing
        CASE(PROBLEM_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%CLASS,"*",ERR,ERROR))//" &
            & is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("PROBLEM_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Executes pre solver routines for a problem.
  SUBROUTINE PROBLEM_SOLVER_PRE_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SOLVER_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            SELECT CASE(PROBLEM%CLASS)
            CASE(PROBLEM_ELASTICITY_CLASS)
              CALL ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_BIOELECTRICS_CLASS)
              CALL BIOELECTRIC_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_FLUID_MECHANICS_CLASS)
              CALL FLUID_MECHANICS_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
              !Do nothing???
            CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
              CALL CLASSICAL_FIELD_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_FITTING_CLASS)
              CALL FITTING_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MODAL_CLASS)
              !Do nothing???
            CASE(PROBLEM_MULTI_PHYSICS_CLASS)
              CALL MULTI_PHYSICS_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem class of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))//" &
                & is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_PRE_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Executes post solver routines for a problem.
  SUBROUTINE PROBLEM_SOLVER_POST_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVER_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            SELECT CASE(PROBLEM%CLASS)
            CASE(PROBLEM_ELASTICITY_CLASS)
              CALL ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_BIOELECTRICS_CLASS)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_FLUID_MECHANICS_CLASS)
              CALL FLUID_MECHANICS_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
              !Do nothing???
            CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
              CALL CLASSICAL_FIELD_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_FITTING_CLASS)
              CALL FITTING_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MODAL_CLASS)
              !Do nothing???
            CASE(PROBLEM_MULTI_PHYSICS_CLASS)
              CALL MULTI_PHYSICS_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem class of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%CLASS,"*",ERR,ERROR))//" &
                & is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_SOLVER_POST_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_POST_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_POST_SOLVE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_SOLVER_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves solver equations for a problem.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        SELECT CASE(SOLVER_EQUATIONS%TIME_DEPENDENCE)
        CASE(SOLVER_EQUATIONS_STATIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(SOLVER_EQUATIONS_QUASISTATIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The solver equations time dependence type of "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Solver equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves dynamic linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,loop_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_TIME_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
              !Get current control loop times. The control loop may be a sub loop below a time loop, so iterate up
              !through loops checking for the time loop
              CONTROL_TIME_LOOP=>CONTROL_LOOP
              DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
                IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                  CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
                  EXIT
                ENDIF
                IF(ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                  CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
                ELSE
                  CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
                ENDIF
              ENDDO
              !Set the solver time
              CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
              !Back-substitute to find flux values for linear problems
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ENDDO !equations_set_idx
            ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves dynamic nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,loop_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_TIME_LOOP
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
        IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
          SOLVERS=>SOLVER%SOLVERS
          IF(ASSOCIATED(SOLVER)) THEN
            CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
            IF(ASSOCIATED(CONTROL_LOOP)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  IF(DYNAMIC_SOLVER%RESTART.OR..NOT.DYNAMIC_SOLVER%SOLVER_INITIALISED) THEN
                    !If we need to restart or we haven't initialised yet, make sure the equations sets are up to date
                    EQUATIONS=>EQUATIONS_SET%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      SELECT CASE(EQUATIONS%LINEARITY)
                      CASE(EQUATIONS_LINEAR)
                        !Assemble the equations
                        CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
                      CASE(EQUATIONS_NONLINEAR)
                        !Evaluate the residuals
                        CALL EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
                      CASE(EQUATIONS_NONLINEAR_BCS)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The equations linearity type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
                          & " is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ENDDO !equations_set_idx
                !Get current control loop times. The control loop may be a sub loop below a time loop, so iterate up
                !through loops checking for the time loop
                CONTROL_TIME_LOOP=>CONTROL_LOOP
                DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
                  IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                    CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
                    EXIT
                  ENDIF
                  IF(ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                    CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
                  ELSE
                    CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
                  ENDIF
                ENDDO
                !Set the solver time
                CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
                !Solve for the next time i.e., current time + time increment
                CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver dynamic solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves quasistatic linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
!     REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
     
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
!               !Get current control loop times
!               CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the current time
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
              !Back-substitute to find flux values for linear problems
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ENDDO !equations_set_idx
            ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves quasistatic nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
   
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
              ! sander - this gives an error, and current time seems to be updated without it
              !Get current control loop times
              !CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Set the solver time
              !CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
             ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves static linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,interface_condition_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
#ifdef TAUPROF
    CHARACTER(12) :: CVAR
    INTEGER :: PHASE(2) = [ 0, 0 ]
    SAVE PHASE
#endif

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Make sure the equations sets are up to date
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
#ifdef TAUPROF
            WRITE (CVAR,'(a8,i2)') 'Assemble',equations_set_idx
            CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
            CALL TAU_PHASE_START(PHASE)
#endif
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)
            !Assemble the equations for linear problems
            CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(PHASE)
#endif
          ENDDO !equations_set_idx
          !Make sure the interface matrices are up to date
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
#ifdef TAUPROF
            WRITE (CVAR,'(a8,i2)') 'Interface',interface_condition_idx
            CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
            CALL TAU_PHASE_START(PHASE)
#endif
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            CALL INTERFACE_CONDITION_ASSEMBLE(INTERFACE_CONDITION,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(PHASE)
#endif
          ENDDO !interface_condition_idx

          !Solve
          CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)

#ifdef TAUPROF
          CALL TAU_STATIC_PHASE_START('EQUATIONS_SET_BACKSUBSTITUTE()')
#endif
          !Back-substitute to find flux values for linear problems
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          ENDDO !equations_set_idx
#ifdef TAUPROF
          CALL TAU_STATIC_PHASE_STOP('EQUATIONS_SET_BACKSUBSTITUTE()')
#endif
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE
  
  !
  !================================================================================================================================
  !

  !>Solves static nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,interface_condition_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
#ifdef TAUPROF
    CHARACTER(12) :: CVAR
    INTEGER :: PHASE(2) = [ 0, 0 ]
    SAVE PHASE
#endif
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Apply boundary conditition
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            !Assemble the equations set
            CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
          ENDDO !equations_set_idx
          !Make sure the interface matrices are up to date
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
#ifdef TAUPROF
            WRITE (CVAR,'(a8,i2)') 'Interface',interface_condition_idx
            CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
            CALL TAU_PHASE_START(PHASE)
#endif
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
!            CALL FrictionlessContact_contactMetricsCalculate(INTERFACE_CONDITION,err,error,*999)
            CALL INTERFACE_CONDITION_ASSEMBLE(INTERFACE_CONDITION,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(PHASE)
#endif
          ENDDO !interface_condition_idx
          !Solve
          CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
          !Update the rhs field variable with residuals or backsubstitute for any linear
          !equations sets
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            EQUATIONS=>EQUATIONS_SET%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              SELECT CASE(EQUATIONS%LINEARITY)
              CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              CASE(EQUATIONS_NONLINEAR)
                CALL EQUATIONS_SET_NONLINEAR_RHS_UPDATE(EQUATIONS_SET,SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              CASE DEFAULT
                CALL FLAG_ERROR("Invalid linearity for equations set equations",ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !equations_set_idx
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !


  !>Solves a solver for a problem.
  SUBROUTINE PROBLEM_SOLVER_SOLVE(SOLVER,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("PROBLEM_SOLVER_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER)) THEN

      IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Solver: ",SOLVER%LABEL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Solver index = ",SOLVER%GLOBAL_NUMBER,ERR,ERROR,*999)
      ENDIF
      
#ifdef TAUPROF
      CALL TAU_STATIC_PHASE_START('Pre solve')
#endif
     CALL PROBLEM_SOLVER_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
#ifdef TAUPROF
      CALL TAU_STATIC_PHASE_STOP('Pre solve')
      
      CALL TAU_STATIC_PHASE_START('Solve')
#endif
      
      IF(ASSOCIATED(SOLVER%SOLVER_EQUATIONS)) THEN
        !A solver with solver equations.
        CALL PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
      ELSE
        !Check for other equations.
        IF(ASSOCIATED(SOLVER%CELLML_EQUATIONS)) THEN
          !A solver with CellML equations.
          CALL PROBLEM_CELLML_EQUATIONS_SOLVE(SOLVER%CELLML_EQUATIONS,ERR,ERROR,*999)
        ELSEIF(SOLVER%SOLVE_TYPE==SOLVER_GEOMETRIC_TRANSFORMATION_TYPE) THEN
!          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Geometric transformation******************",ERR,ERROR,*999)
          CALL Problem_SolverGeometricTransformationSolve(SOLVER%geometricTransformationSolver,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Solver does not have any equations associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF

#ifdef TAUPROF
      CALL TAU_STATIC_PHASE_STOP('Solve')
      
      CALL TAU_STATIC_PHASE_START('Post solve')
#endif
      CALL PROBLEM_SOLVER_POST_SOLVE(SOLVER,ERR,ERROR,*999)
#ifdef TAUPROF
      CALL TAU_STATIC_PHASE_STOP('Post solve')
#endif
      
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_SOLVE")
    RETURN 1
    
  END SUBROUTINE PROBLEM_SOLVER_SOLVE

  !
  !================================================================================================================================
  !

  !>Destroy the solvers for a problem. \see OPENCMISS::CMISSProblemSolversDestroy
  SUBROUTINE PROBLEM_SOLVERS_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the solvers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVERS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN        
        CALL CONTROL_LOOP_SOLVERS_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVERS_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_DESTROY

  !
  !================================================================================================================================
  !

  !>Set boundary conditions for solver equations according to the analytic equations. \see OPENCMISS_CMISSProblemSolverEquationsBoundaryConditionsAnalytic
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to get the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set is not associated for index "//TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*", &
                  & ERR,ERROR))//".",ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver equations boundary conditions are not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN 1

  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC

  !
  !================================================================================================================================
  !

  !>Finish the creation of the solver equations for the problem. \see OPENCMISS::CMISSProblemSolverEquationsCreateFinish
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the solver equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN      
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_FINISH_ACTION
      !Finish problem specific startup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of solver equations for a problem. \see OPENCMISS::CMISSProblemSolverEquationsCreateStart
  !>The default values of the SOLVER attributes are:
  !>- SOLVE_TYPE: 1 (SOLVER_LINEAR_TYPE)
  !>- OUTPUT_TYPE: 0 (SOLVER_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (SOLVER_SPARSE_MATRICES)
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variablesg
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Initialise the problem setup information
      CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE
      PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
      !Start the problem specific control setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the problem setup information
      CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !!TODO: this should be removed - just call the solver equations destroy directly???
  
  !>Destroy the solver equations for a problem. \see OPENCMISS::CMISSProblemSolverEquationsDestroy
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY(CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Solves geometric transformation for a field 
  SUBROUTINE Problem_SolverGeometricTransformationSolve(geometricTransformationSolver,err,error,*) !\todo: Add rotation operations.
    
   !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: GeometricTransformationSolver !<A pointer to the geometric transformation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER :: loadIncrementLoop
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER :: simpleLoop
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: fixedLoop
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: whileLoop
    INTEGER(INTG) :: componentIdx,versionIdx,derivativeIdx,nodeIdx,domainNodeIdx,noGeomComp
    INTEGER(INTG) :: localNodeNumber,userNodeNumber,incrementIdx,iterationNumber,noNodes
    REAL(DP) :: nodalParameters(3),nodalParametersTrans(3),transformationMatrix(4,4)
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    LOGICAL :: transformBC=.FALSE.,sameBases=.TRUE.,transformAllNodes=.TRUE.,nodeExist,ghostNode
    
    CALL ENTERS("Problem_SolverGeometricTransformationSolve",err,error,*999) 
    
    IF(ASSOCIATED(geometricTransformationSolver)) THEN
      IF(ASSOCIATED(geometricTransformationSolver%field)) THEN
        fieldVariable=>geometricTransformationSolver%field%VARIABLE_TYPE_MAP(geometricTransformationSolver%fieldVariableType)%PTR
        IF(ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)%PTR)) transformBC=.TRUE. !if the BC is defined on the field variable to be transformed
        noGeomComp=SIZE(geometricTransformationSolver%transformationMatrices,1)-1 ! Number of geometric components
        !**********************************************************************************************************************
        !Determine iteration/load increment number 
        IF(geometricTransformationSolver%numberOfIncrements>1) THEN
          solver=>geometricTransformationSolver%solver
          IF(ASSOCIATED(solver)) THEN
            solvers=>solver%SOLVERS
            IF(ASSOCIATED(solvers)) THEN
              controlLoop=>solvers%CONTROL_LOOP
              IF(ASSOCIATED(controlLoop)) THEN
                SELECT CASE(controlLoop%LOOP_TYPE)
                CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
                  simpleLoop=>controlLoop%SIMPLE_LOOP
                  IF(ASSOCIATED(simpleLoop)) THEN
                    iterationNumber=1
                  ELSE
                    CALL FLAG_ERROR("Simple loop is not associated.",err,error,*999)
                  ENDIF
                CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
                  fixedLoop=>controlLoop%FIXED_LOOP
                  IF(ASSOCIATED(fixedLoop)) THEN
                    iterationNumber=fixedLoop%ITERATION_NUMBER
                  ELSE
                    CALL FLAG_ERROR("Fixed loop is not associated.",err,error,*999)
                  ENDIF
                CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
                  CALL FLAG_ERROR("Geometric transformation for time loop is not implemented.",err,error,*999)
                CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
                  whileLoop=>controlLoop%WHILE_LOOP
                  IF(ASSOCIATED(whileLoop)) THEN
                    iterationNumber=whileLoop%ITERATION_NUMBER
                  ELSE
                    CALL FLAG_ERROR("Simple loop is not associated.",err,error,*999)
                  ENDIF
                CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
                  loadIncrementLoop=>controlLoop%LOAD_INCREMENT_LOOP
                  IF(ASSOCIATED(loadIncrementLoop)) THEN
                    iterationNumber=loadIncrementLoop%ITERATION_NUMBER
                  ELSE
                    CALL FLAG_ERROR("Load increment loop is not associated.",err,error,*999)
                  ENDIF
                END SELECT
                IF(iterationNumber>geometricTransformationSolver%numberOfIncrements) THEN
                  !If load increment is not specified for that iteration, loop around
                  incrementIdx=MOD(iterationNumber-1,geometricTransformationSolver%numberOfIncrements)+1
                ELSE
                  incrementIdx=iterationNumber !If load increment is specified for that iteration, use that load increment
                ENDIF
              ELSE
                CALL FLAG_ERROR("Control loop is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solvers is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver is not associated.",err,error,*999)
          ENDIF
        ELSE
          incrementIdx=1
        ENDIF
        !Determine the transformation matrix to use
        IF(geometricTransformationSolver%arbitraryPath .OR. geometricTransformationSolver%numberOfIncrements==1) THEN
          transformationMatrix(1:noGeomComp+1,1:noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
            & (1:noGeomComp+1,1:noGeomComp+1,incrementIdx)
        ELSE !If need to scale transformation matrix (i.e. transformation applied through several load increment.)
          IF(incrementIdx==1) THEN ! 1st load increment, rotation is applied
            transformationMatrix(1:noGeomComp,1:noGeomComp)=geometricTransformationSolver%transformationMatrices &
              & (1:noGeomComp,1:noGeomComp,1)
          ELSE !No rotation operation in any other load increments
            DO componentIdx=1,noGeomComp
              transformationMatrix(componentIdx,componentIdx)=1.0_DP
            ENDDO !componentIdx
          ENDIF
          !Translation is scaled for every load increment 
          IF(ALLOCATED(geometricTransformationSolver%scalings)) THEN
            transformationMatrix(1:noGeomComp,noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
              & (1:noGeomComp,noGeomComp+1,1)*geometricTransformationSolver%scalings(incrementIdx)
          ELSE !if no scaling just take 1/numberOfIncrements as scaling
            transformationMatrix(1:noGeomComp,noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
              & (1:noGeomComp,noGeomComp+1,1)/geometricTransformationSolver%numberOfIncrements
          ENDIF
        ENDIF
        !**********************************************************************************************************************
        ! Transform the field
        ! Determine if the all components have the same mesh components/ bases
        DO componentIdx=1,noGeomComp-1
          IF(fieldVariable%COMPONENTS(componentIdx)%MESH_COMPONENT_NUMBER/= &
            & fieldVariable%COMPONENTS(componentIdx+1)%MESH_COMPONENT_NUMBER) sameBases=.FALSE.
        ENDDO
        IF(sameBases) THEN
          domain=>fieldVariable%COMPONENTS(1)%DOMAIN !Use the 1st component domain since they are the same for all components
          domainNodes=>domain%TOPOLOGY%NODES
          IF(ASSOCIATED(domain)) THEN
            !Determine how many nodes are to be transformed
            IF(ALLOCATED(geometricTransformationSolver%nodeUserNumbers)) THEN
              noNodes=SIZE(geometricTransformationSolver%nodeUserNumbers)
              transformAllNodes=.FALSE.
            ELSE
              noNodes=domainNodes%NUMBER_OF_NODES !Not including ghost nodes
            ENDIF
            DO nodeIdx=1,noNodes
              !Get user number for the nodes to be transformed
              IF(transformAllNodes) THEN
                domainNodeIdx=nodeIdx
                localNodeNumber=domainNodes%NODES(domainNodeIdx)%LOCAL_NUMBER !localNodeNumber is the same as domainNodeIdx
                userNodeNumber=domainNodes%NODES(domainNodeIdx)%USER_NUMBER
                nodeExist=.TRUE. ! Transform all local nodes on the domain, ghost nodes are already excluded
                ghostNode=.FALSE.
              ELSE
                userNodeNumber=geometricTransformationSolver%nodeUserNumbers(nodeIdx)
                CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(domain%TOPOLOGY,userNodeNumber,nodeExist,domainNodeIdx, &
                  & ghostNode,err,error,*999)
              ENDIF
              IF((nodeExist) .AND. (.NOT.ghostNode)) THEN !only transform local nodes, not ghost nodes
                DO derivativeIdx=1,domainNodes%NODES(domainNodeIdx)%NUMBER_OF_DERIVATIVES
                  DO versionIdx=1,domainNodes%NODES(domainNodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                    DO componentIdx=1,noGeomComp !Get all component for a nodal derivative
                      CALL FIELD_PARAMETER_SET_GET_NODE(geometricTransformationSolver%field,geometricTransformationSolver% &
                        & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                        & nodalParameters(componentIdx),err,error,*999)
                    ENDDO !componentIdx
                    !Rotate the nodal parameters
                    nodalParametersTrans(1:noGeomComp)=MATMUL(transformationMatrix(1:noGeomComp,1:noGeomComp), &
                      & nodalParameters(1:noGeomComp))
                    DO componentIdx=1,noGeomComp !Update all component for a nodal derivative
                      CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricTransformationSolver%field,geometricTransformationSolver% &
                        & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                        & nodalParametersTrans(componentIdx),err,error,*999)
                      IF(derivativeIdx==1) THEN ! Translate nodal coordinate
                        CALL FIELD_PARAMETER_SET_ADD_NODE(geometricTransformationSolver%field,geometricTransformationSolver% &
                          & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                          & transformationMatrix(componentIdx,1+noGeomComp),err,error,*999)
                      ENDIF !derivativeIdx==1
                      IF(transformBC) THEN
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricTransformationSolver%field,geometricTransformationSolver% &
                          & fieldVariableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber, &
                          & componentIdx,nodalParametersTrans(componentIdx),err,error,*999)
                        IF(derivativeIdx==1) THEN ! Translate nodal coordinate for BC
                          CALL FIELD_PARAMETER_SET_ADD_NODE(geometricTransformationSolver%field,geometricTransformationSolver% &
                            & fieldVariableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber, &
                            & componentIdx,transformationMatrix(componentIdx,1+noGeomComp),err,error,*999)
                        ENDIF !derivativeIdx==1
                      ENDIF !transformBC
                    ENDDO !componentIdx
                  ENDDO !versionIdx
                ENDDO !derivativeIdx
              ENDIF !only transform local nodes, not ghost nodes
            ENDDO !nodeIdx
          ELSE
            CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Transformation for different component bases not implemented.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The field of geometric transformation solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Geometric transformation solver is not associated.",err,error,*999)
    ENDIF
      
    CALL EXITS("Problem_SolverGeometricTransformationSolve")
    RETURN
999 CALL ERRORS("Problem_SolverGeometricTransformationSolve",err,error)
    CALL EXITS("Problem_SolverGeometricTransformationSolve")
    RETURN 1
  END SUBROUTINE Problem_SolverGeometricTransformationSolve

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OPENCMISS::CMISSProblemSolverGet
  SUBROUTINE PROBLEM_SOLVER_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLVER_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_SOLVER_GET_1(PROBLEM,[CONTROL_LOOP_IDENTIFIER],SOLVER_INDEX,SOLVER,ERR,ERROR,*999) 
       
    CALL EXITS("PROBLEM_SOLVER_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_GET_0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OPENCMISS::CMISSProblemSolverGet
  SUBROUTINE PROBLEM_SOLVER_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier to get the solver for.
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLVER_GET_1",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*998)
      ELSE
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(.NOT.ASSOCIATED(SOLVER)) CALL FLAG_ERROR("Solvers solver is not associated.",ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_GET_1")
    RETURN
999 NULLIFY(SOLVER)
998 CALL ERRORS("PROBLEM_SOLVER_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_GET_1
  
  !
  !================================================================================================================================
  !

  !>Monitors the problem nonlinear solve
  SUBROUTINE Problem_SolverNonlinearMonitor(solver,iterationNumber,residualNorm,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to monitor
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The number of iterations
    REAL(DP), INTENT(IN) :: residualNorm !<The residual norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    LOGICAL :: reproject
    TYPE(VARYING_STRING) :: localError
    !\todo Temporarily added the variables below to allow the interface condition to be used in the single region contact problem
    !to be manually specified. Need to Generalise.
    INTEGER(INTG) :: equationsSetGlobalNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_TYPE), POINTER :: interface

    CALL ENTERS("Problem_SolverNonlinearMonitor",err,error,*998)
    
    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%SOLVERS
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop)) THEN
          problem=>controlLoop%PROBLEM
          IF(ASSOCIATED(problem)) THEN
            SELECT CASE(problem%CLASS)
            CASE(PROBLEM_ELASTICITY_CLASS)
              SELECT CASE(problem%TYPE)
              CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE)
                !Output meshes at iterations
                IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
                  nonlinearSolver=>solver%NONLINEAR_SOLVER
                  IF(ASSOCIATED(nonlinearSolver)) THEN
                    CALL Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*999)
                  ENDIF
                ENDIF
              CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
                SELECT CASE(problem%SUBTYPE)
                CASE(PROBLEM_LE_CONTACT_TRANSFORM_SUBTYPE,PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE) !Reproject at iteration 0 before the nonlinear solve to update xi location since the field is transformed.
                  IF(iterationNumber==0) THEN
                    reproject=.TRUE.
                  ELSE
                    reproject=.FALSE.
                  ENDIF
                CASE(PROBLEM_LE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_LE_CONTACT_REPROJECT_SUBTYPE, &
                    & PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)
                  reproject=.TRUE.
                CASE DEFAULT
                  localError="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(problem%SUBTYPE,"*",err,error))//" &
                    & is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
                IF(Reproject) THEN
                  solverEquations=>solver%SOLVER_EQUATIONS
                  IF(ASSOCIATED(solverEquations)) THEN
                    solverMapping=>solverEquations%SOLVER_MAPPING
                    IF(ASSOCIATED(solverMapping)) THEN
                      !\todo Temporarily comment out the looping through of interface conditions added to the solver as these are not
                      !present in the single region contact problem. Needs to generalised.
                      !DO interfaceConditionIdx=1,solverMapping%NUMBER_OF_INTERFACE_CONDITIONS
                        !interfaceCondition=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
                        equationsSetGlobalNumber=1
                        interfaceGlobalNumber=1
                        interfaceConditionGlobalNumber=1
                        interfaceCondition=>solverMapping%EQUATIONS_SETS(equationsSetGlobalNumber)%PTR%REGION%PARENT_REGION% &
                          & INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR%INTERFACE_CONDITIONS% &
                          & INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
                        IF(ASSOCIATED(interfaceCondition)) THEN
                          IF(interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR .OR. &
                              & interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_OPERATOR) THEN !Only reproject for contact operator
                            IF(interfaceCondition%integrationType==INTERFACE_CONDITION_DATA_POINTS_INTEGRATION) THEN !Only reproject for data point interpolated field
                              interface=>interfaceCondition%INTERFACE
                              IF(ASSOCIATED(interface)) THEN
                                !CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"**************** Reproject! ****************",ERR,ERROR,*999)
                                !CALL InterfacePointsConnectivity_DataReprojection(interface,interfaceCondition,err,error,*999)
                                !CALL FrictionlessContact_contactMetricsCalculate(interfaceCondition,err,error,*999)
                                !CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"************** Contact residual! ***********",ERR,ERROR,*999)
                                !CALL EQUATIONS_SET_RESIDUAL_CONTACT_UPDATE_STATIC_FEM(solverMapping%EQUATIONS_SETS(1)%PTR,&
                                !  & ERR,ERROR,*999)
                                !\todo Temporarily commented out INTERFACE_CONDITION_ASSEMBLE as the interface matrices are not
                                ! required for the single region contact problem. Needs to generalised.
                                !CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
                              ELSE
                                CALL FLAG_ERROR("Interface is not associated for nonlinear solver equations mapping.", &
                                  & err,error,*999)
                              ENDIF
                            ENDIF
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Interface condition is not associated for nonlinear solver equations mapping.", &
                            & err,error,*999)
                        ENDIF
                      !ENDDO !interfaceConditionIdx
                    ELSE
                      CALL FLAG_ERROR("Nonlinear solver equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Nonlinear solver equations is not associated.",err,error,*999)
                  ENDIF
                ENDIF !Reproject
                !Output meshes at iterations
                IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
                  nonlinearSolver=>solver%NONLINEAR_SOLVER
                  IF(ASSOCIATED(nonlinearSolver)) THEN
                    CALL Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*999)
                  ENDIF
                ENDIF
              CASE DEFAULT
                localError="The problem type of "//TRIM(NUMBER_TO_VSTRING(problem%TYPE,"*",err,error))//" &
                  & is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
                & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
              !Do nothing???
            CASE DEFAULT
              localError="The problem class of "//TRIM(NUMBER_TO_VSTRING(problem%CLASS,"*",err,error))//" &
                & is invalid."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",err,error,*999)
        ENDIF
      ENDIF
      !Nonlinear solve monitor--progress output if required
      IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
        nonlinearSolver=>solver%NONLINEAR_SOLVER
        IF(ASSOCIATED(nonlinearSolver)) THEN
          CALL SOLVER_NONLINEAR_MONITOR(nonlinearSolver,iterationNumber,residualNorm,err,error,*999)
        ELSE
          CALL FLAG_ERROR("Nonlinear solver is not associated.",err,error,*999)
        ENDIF
      ELSE
        localError="Invalid solve type. The solve type of "//TRIM(NUMBER_TO_VSTRING(solver%SOLVE_TYPE,"*",err,error))// &
          & " does not correspond to a nonlinear solver."
        CALL FLAG_ERROR(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("Problem_SolverNonlinearMonitor")
    RETURN
999 NULLIFY(SOLVER)
998 CALL ERRORS("Problem_SolverNonlinearMonitor",err,error)
    CALL EXITS("Problem_SolverNonlinearMonitor")
    RETURN 1
  END SUBROUTINE Problem_SolverNonlinearMonitor
  
  !
  !================================================================================================================================
  !

  !> Output fields at Newton iterations. This is in temporarily for debug output. It may be removed at a later date.
  SUBROUTINE Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to solver to output the fields for
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<Iteration number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,load_step
    LOGICAL :: dirExists
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to region to output the fields for
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping 
    TYPE(FIELDS_TYPE), POINTER :: fields
    TYPE(VARYING_STRING) :: fileName,method,directory
    
    INTEGER(INTG) :: interfaceConditionIdx, interfaceElementNumber, dataPointIdx, globalDataPointNum, elementNum, &
      & coupledMeshFaceLineNumber, coupledMeshIdx,component
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(DATA_POINTS_TYPE), POINTER :: interfaceDatapoints
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection

    TYPE(PROBLEM_TYPE), POINTER :: problem

    INTEGER(INTG) :: IUNIT
    CHARACTER(LEN=100) :: filenameOutput,groupname

    TYPE(VARYING_STRING) :: fileToCheck,localError
    LOGICAL :: fileExists
    INTEGER(INTG) :: firstIterationNumber, solve_call, max_solve_calls

    !\todo Temporarily added the variables below to allow the interface condition to be used in the single region contact problem
    !to be manually specified. Need to Generalise.
    INTEGER(INTG) :: equationsSetGlobalNumber,interfaceGlobalNumber,interfaceConditionGlobalNumber,connectedFace,bodyidx
    REAL(DP) :: xi(2)

    CALL ENTERS("Problem_SolverNewtonFieldsOutput",err,error,*999)
    
    IF(ASSOCIATED(solver%SOLVER_EQUATIONS))THEN
      solverMapping=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING
      problem=>solver%SOLVERS%CONTROL_LOOP%PROBLEM

      SELECT CASE(problem%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        SELECT CASE(problem%TYPE)
        CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE, &
          & PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)

          !\todo Temporarily commenting out diagnostics for single region contact problem. Needs to be generalised.
          !IF(DIAGNOSTICS1) THEN
            directory="results_iter/"
            INQUIRE(FILE=CHAR(directory),EXIST=dirExists)
            IF(.NOT.dirExists) THEN
              CALL SYSTEM(CHAR("mkdir "//directory))
            ENDIF

            ! Find how many times the problem solve command has been issued.
            max_solve_calls=100
            coupledMeshIdx=1
            load_step=1
            firstIterationNumber=0
            DO solve_call=1,max_solve_calls
              fileToCheck=directory// &
                & "mesh"//TRIM(NUMBER_TO_VSTRING(coupledMeshIdx,"*",err,error))// &
                & "_solveCall"//TRIM(NUMBER_TO_VSTRING(solve_call,"*",err,error))// &
                & "_load"//TRIM(NUMBER_TO_VSTRING(load_step,"*",err,error))// &
                & "_iter"//TRIM(NUMBER_TO_VSTRING(firstIterationNumber,"*",err,error))//".part0.exnode"
              INQUIRE(FILE=CHAR(fileToCheck),EXIST=fileExists)
              IF(.NOT.fileExists) THEN
                EXIT
              ENDIF
            ENDDO

            load_step=solver%SOLVERS%CONTROL_LOOP%LOAD_INCREMENT_LOOP%ITERATION_NUMBER

            IF((iterationNumber > 0).OR.(load_step > 1))THEN
              solve_call = solve_call - 1
            ENDIF

            WRITE(*,'(1X,''SolveCall: '',I4)') solve_call
            WRITE(*,'(1X,''  LoadStep: '',I4)') load_step
            WRITE(*,'(1X,''    Iteration: '',I4)') iterationNumber

            DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
              region=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%PTR%REGION
              IF(ASSOCIATED(region))THEN
                NULLIFY(fields)
                fields=>region%FIELDS
                fileName=directory//"mesh"//TRIM(NUMBER_TO_VSTRING(equationsSetIdx,"*",err,error))// &
                  & "_solveCall"//TRIM(NUMBER_TO_VSTRING(solve_call,"*",err,error))// &
                  & "_load"//TRIM(NUMBER_TO_VSTRING(load_step,"*",err,error))// &
                  & "_iter"//TRIM(NUMBER_TO_VSTRING(iterationNumber,"*",err,error))
                method="FORTRAN"
                CALL FIELD_IO_ELEMENTS_EXPORT(fields,fileName,method,err,error,*999)
                CALL FIELD_IO_NODES_EXPORT(fields,fileName,method,err,error,*999)
              ELSE
                CALL FLAG_ERROR("Region is not associated.",err,error,*999)
              ENDIF
            ENDDO
          !ENDIF

        CASE DEFAULT
          localError="The problem type of "//TRIM(NUMBER_TO_VSTRING(problem%TYPE,"*",err,error))//" &
            & is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
          & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
        !Do nothing???
      CASE DEFAULT
        localError="The problem class of "//TRIM(NUMBER_TO_VSTRING(problem%CLASS,"*",err,error))//" &
          & is invalid."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT

      SELECT CASE(problem%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        SELECT CASE(problem%TYPE)
        CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE)
          ! Pass
        CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)

!\todo Temporarily commenting out the looping through of interface conditions added to the solver as these are not
!present in the single region contact problem. Needs to be generalised.
!          IF(DIAGNOSTICS1) THEN
!            IUNIT = 300
!            DO interfaceConditionIdx=1,solverMapping%NUMBER_OF_INTERFACE_CONDITIONS
!              interfaceCondition=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
!              interface=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR%interface
!              pointsConnectivity=>interface%pointsConnectivity
!              interfaceDatapoints=>interface%DATA_POINTS
!              IF(ASSOCIATED(pointsConnectivity)) THEN
!                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
!                  filenameOutput=directory//"PointsConnectivity"//TRIM(NUMBER_TO_VSTRING(coupledMeshIdx,"*",err,error))// &
!                    & "_solveCall"//TRIM(NUMBER_TO_VSTRING(solve_call,"*",err,error))// &
!                    & "_load"//TRIM(NUMBER_TO_VSTRING(load_step,"*",err,error))// &
!                    & "_iter"//TRIM(NUMBER_TO_VSTRING(iterationNumber,"*",err,error))//".exdata"
!                  OPEN(UNIT=IUNIT,FILE=filenameOutput,STATUS="UNKNOWN",ACTION="WRITE",IOSTAT=ERR)
!                  groupname="PointsConnectivity"//TRIM(NUMBER_TO_VSTRING(coupledMeshIdx,"*",err,error))
!                  WRITE(IUNIT,'( '' Group name: '',A)') groupname
!                  WRITE(IUNIT,'(1X,''#Fields=4'')')
!                  WRITE(IUNIT,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
!                  WRITE(IUNIT,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''2) error, field, rectangular cartesian, #Components=3'')')
!                  WRITE(IUNIT,'(1X,''  x.  Value index= 4, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  y.  Value index= 5, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  z.  Value index= 6, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''3) projectedCoordinate, field, rectangular cartesian, #Components=3'')')
!                  WRITE(IUNIT,'(1X,''  x.  Value index= 7, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  y.  Value index= 8, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''  z.  Value index= 9, #Derivatives=0'')')
!                  WRITE(IUNIT,'(1X,''4) exitTag, field, rectangular cartesian, #Components=1'')')
!                  WRITE(IUNIT,'(1X,''  tag.  Value index= 10, #Derivatives=0'')')
!                  dependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
!                    & DEPENDENT%DEPENDENT_FIELD
!                  NULLIFY(interpolationParameters)
!                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(dependentField,interpolationParameters,err,error, &
!                    & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
!                  NULLIFY(interpolatedPoints)
!                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
!                    & FIELD_GEOMETRIC_COMPONENTS_TYPE)
!                  interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
!                  dataProjection=>interfaceDatapoints%DATA_PROJECTIONS(coupledMeshIdx+1)%PTR
!                  DO interfaceElementNumber=1,SIZE(pointsConnectivity%coupledElements,1)
!                    decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
!                      & elementDataPoint(interfaceElementNumber)
!                    DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
!                      globalDataPointNum=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
!                      WRITE(IUNIT,'(1X,''Node:'',I4)') globalDataPointNum
!                      DO component=1,3
!                        WRITE(IUNIT,'(1X,3E25.15)') interfaceDatapoints%DATA_POINTS(globalDataPointNum)%position(component)
!                      ENDDO !component
!                      elementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,coupledMeshIdx)% &
!                        & coupledMeshElementNumber
!                      coupledMeshFaceLineNumber=dependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS% &
!                        & ELEMENTS(elementNum)% &
!                        & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(globalDataPointNum,coupledMeshIdx)% &
!                        & elementLineFaceNumber)
!                      CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
!                        & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
!                      CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,pointsConnectivity%pointsConnectivity(globalDataPointNum, &
!                        & coupledMeshIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
!                      DO component=1,3
!                        WRITE(IUNIT,'(1X,3E25.15)') interpolatedPoint%VALUES(component,NO_PART_DERIV) - &
!                          & interfaceDatapoints%DATA_POINTS(globalDataPointNum)%position(component)
!                      ENDDO !component
!                      DO component=1,3
!                        WRITE(IUNIT,'(1X,3E25.15)') interpolatedPoint%VALUES(component,NO_PART_DERIV)
!                      ENDDO !component
!                      WRITE(IUNIT,'(1X,I2)') dataProjection%DATA_PROJECTION_RESULTS(globalDataPointNum)%EXIT_TAG
!                    ENDDO !dataPointIdx
!                  ENDDO !interfaceElementNumber
!                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
!                  CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoints,err,error,*999)
!                  OPEN(UNIT=IUNIT)
!                ENDDO !coupledMeshIdx
!              ENDIF
!            ENDDO !interfaceConditionIdx
!          ENDIF

          !Output fields for single region contact problem
          !\todo Needs to be generalised.
          IUNIT = 300
          equationsSetGlobalNumber=1
          interfaceGlobalNumber=1
          interfaceConditionGlobalNumber=1
          interfaceCondition=>solverMapping%EQUATIONS_SETS(equationsSetGlobalNumber)%PTR%REGION%PARENT_REGION% &
            & INTERFACES%INTERFACES(interfaceGlobalNumber)%PTR%INTERFACE_CONDITIONS% &
            & INTERFACE_CONDITIONS(interfaceConditionGlobalNumber)%PTR
          interface=>interfaceCondition%interface
          pointsConnectivity=>interface%pointsConnectivity
          interfaceDatapoints=>interface%DATA_POINTS
          IF(ASSOCIATED(pointsConnectivity)) THEN
            DO bodyidx=1,2
              filenameOutput=directory//"PointsConnectivity"//TRIM(NUMBER_TO_VSTRING(bodyidx,"*",err,error))// &
                & "_solveCall"//TRIM(NUMBER_TO_VSTRING(solve_call,"*",err,error))// &
                & "_load"//TRIM(NUMBER_TO_VSTRING(load_step,"*",err,error))// &
                & "_iter"//TRIM(NUMBER_TO_VSTRING(iterationNumber,"*",err,error))//".exdata"
              OPEN(UNIT=IUNIT,FILE=filenameOutput,STATUS="UNKNOWN",ACTION="WRITE",IOSTAT=ERR)
              groupname="PointsConnectivity"//TRIM(NUMBER_TO_VSTRING(bodyidx,"*",err,error))
              WRITE(IUNIT,'( '' Group name: '',A)') groupname
              WRITE(IUNIT,'(1X,''#Fields=2'')')
              WRITE(IUNIT,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
              WRITE(IUNIT,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
              WRITE(IUNIT,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
              WRITE(IUNIT,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
              WRITE(IUNIT,'(1X,''2) exitTag, field, rectangular cartesian, #Components=1'')')
              WRITE(IUNIT,'(1X,''  tag.  Value index= 4, #Derivatives=0'')')
              dependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(bodyidx)%PTR% &
                & DEPENDENT%DEPENDENT_FIELD
              NULLIFY(interpolationParameters)
              CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(dependentField,interpolationParameters,err,error, &
                & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              NULLIFY(interpolatedPoints)
              CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
              dataProjection=>interfaceDatapoints%DATA_PROJECTIONS(bodyidx+1)%PTR
              DO globalDataPointNum=1,SIZE(pointsConnectivity%pointsConnectivity,1)
                WRITE(IUNIT,'(1X,''Node:'',I4)') globalDataPointNum
                elementNum=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%coupledMeshElementNumber
                connectedFace=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%elementLineFaceNumber
                xi=pointsConnectivity%pointsConnectivity(globalDataPointNum,bodyIdx)%reducedXi
                coupledMeshFaceLineNumber=dependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNum)% &
                  & ELEMENT_FACES(connectedFace)
                CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                DO component=1,3
                  WRITE(IUNIT,'(1X,3E25.15)') interpolatedPoint%VALUES(component,NO_PART_DERIV)
                ENDDO !component
                WRITE(IUNIT,'(1X,I2)') dataProjection%DATA_PROJECTION_RESULTS(globalDataPointNum)%EXIT_TAG
              ENDDO !dataPointIdx
              CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
              CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoints,err,error,*999)
              OPEN(UNIT=IUNIT)
            ENDDO !bodyidx
          ENDIF

        CASE DEFAULT
          localError="The problem type of "//TRIM(NUMBER_TO_VSTRING(problem%TYPE,"*",err,error))//" &
            & is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
          & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
        !Do nothing???
      CASE DEFAULT
        localError="The problem class of "//TRIM(NUMBER_TO_VSTRING(problem%CLASS,"*",err,error))//" &
          & is invalid."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT

    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("Problem_SolverNewtonFieldsOutput")
    RETURN
999 CALL ERRORS("Problem_SolverNewtonFieldsOutput",err,error)
    CALL EXITS("Problem_SolverNewtonFieldsOutput")
    RETURN 1
  END SUBROUTINE Problem_SolverNewtonFieldsOutput
  
  !
  !================================================================================================================================
  !

  !>Gets the problem specification i.e., problem class, type and subtype for a problem identified by a pointer. \see OPENCMISS::CMISSProblemSpecificationGet
  SUBROUTINE PROBLEM_SPECIFICATION_GET(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_CLASS !<On return, The problem class to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem equation type to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the problem subtype to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SPECIFICATION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        PROBLEM_CLASS=PROBLEM%CLASS
        SELECT CASE(PROBLEM_CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FITTING_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_OPTIMISATION_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SPECIFICATION_GET")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_GET",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_GET")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem specification i.e., problem class, type and subtype for a problem identified by a pointer. \see OPENCMISS::CMISSProblemSpecificationSet
  SUBROUTINE PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_CLASS !<The problem class to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem equation type to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SETUP_TYPE) :: PROBLEM_SETUP_INFO
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SPECIFICATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CALL FLAG_ERROR("Problem has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(PROBLEM_CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL ELASTICITY_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_BIOELECTRICS_CLASS)
          CALL BIOELECTRIC_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FITTING_CLASS)
          CALL FITTING_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_OPTIMISATION_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        !Initialise the problem setup information
        CALL PROBLEM_SETUP_INITIALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        PROBLEM_SETUP_INFO%SETUP_TYPE=PROBLEM_SETUP_INITIAL_TYPE
        PROBLEM_SETUP_INFO%ACTION_TYPE=PROBLEM_SETUP_START_ACTION
        !Finish the problem specific setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the problem setup information
        CALL PROBLEM_SETUP_FINALISE(PROBLEM_SETUP_INFO,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SPECIFICATION_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in PROBLEM a pointer to the problem identified by USER_NUMBER. If no problem with that USER_NUMBER exists PROBLEM is left nullified.
  SUBROUTINE PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<On return a pointer to the problem with the given user number. If no problem with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx

    CALL ENTERS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL FLAG_ERROR("Problem is already associated.",ERR,ERROR,*999)
    ELSE
      problem_idx=1
      DO WHILE(problem_idx<=PROBLEMS%NUMBER_OF_PROBLEMS.AND..NOT.ASSOCIATED(PROBLEM))
        IF(PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          PROBLEM=>PROBLEMS%PROBLEMS(problem_idx)%PTR
        ELSE
          problem_idx=problem_idx+1
        ENDIF
      ENDDO
    ENDIF
    
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE PROBLEM_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all problems and deallocates all memory.
  SUBROUTINE PROBLEMS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEMS_FINALISE",ERR,ERROR,*999)

    DO WHILE(PROBLEMS%NUMBER_OF_PROBLEMS>0)
      CALL PROBLEM_DESTROY(PROBLEMS%PROBLEMS(1)%PTR,ERR,ERROR,*999)
    ENDDO !problem_idx
    
    CALL EXITS("PROBLEMS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEMS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_FINALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all problems.
  SUBROUTINE PROBLEMS_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEMS_INITIALISE",ERR,ERROR,*999)

    PROBLEMS%NUMBER_OF_PROBLEMS=0
    NULLIFY(PROBLEMS%PROBLEMS)
    
    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEMS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Updates the dependent variables from the solver solution for all dynamic solvers under the time control loop
  RECURSIVE SUBROUTINE PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer the time control loop to update the variables from
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer the solvers to update the variables from
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to update the variables from
    INTEGER(INTG) :: solver_idx

    NULLIFY(SOLVER)

    CALL ENTERS("PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        !If there are no sub loops then get the solvers for this loop.
        SOLVERS=>CONTROL_LOOP%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
            SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_DYNAMIC_TYPE)
              CALL SOLVER_VARIABLES_DYNAMIC_FIELD_PREVIOUS_VALUES_UPDATE(SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              !Do nothing
            END SELECT
          ENDDO !solver_idx
        ELSE
          CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE")
    RETURN 1

  END SUBROUTINE PROBLEM_CONTROL_LOOP_PREVIOUS_VALUES_UPDATE

  !
  !================================================================================================================================
  !

  
END MODULE PROBLEM_ROUTINES

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the Jacobian for a Newton like nonlinear solver
SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC(SNES,X,A,B,FLAG,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
  INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    NONLINEAR_SOLVER=>CTX%NONLINEAR_SOLVER
    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        SOLVER_EQUATIONS=>CTX%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
              SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
              IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
                IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                  CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(SOLVER_VECTOR,X,ERR,ERROR,*999)
                  
                  CALL PROBLEM_SOLVER_JACOBIAN_EVALUATE(CTX,ERR,ERROR,*999)
                  
                  CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*998)              
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              LOCAL_ERROR="The number of solver matrices of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                & " is invalid. There should be 1 solver matrix."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*998)
        ENDIF
!!TODO: move this to PROBLEM_SOLVER_JACOBIAN_EVALUATE or elsewhere?
        NEWTON_SOLVER%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS=NEWTON_SOLVER%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS+1
      ELSE
        CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver nonlinear solver is not associated.",ERR,ERROR,*998)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*998)
  ENDIF
  
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL WRITE_ERROR(ERR,ERROR,*997)
997 CALL FLAG_WARNING("Error evaluating nonlinear Jacobian.",ERR,ERROR,*996)
996 RETURN 
END SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the Jacobian for a Newton like nonlinear solver using PETSc's FD Jacobian calculation
SUBROUTINE PROBLEM_SOLVER_JACOBIAN_FD_CALCULATE_PETSC(SNES,X,A,B,FLAG,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
  INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
  TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(PETSC_MATFDCOLORING_TYPE), POINTER :: JACOBIAN_FDCOLORING
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    NONLINEAR_SOLVER=>CTX%NONLINEAR_SOLVER
    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
        IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
          SOLVER_EQUATIONS=>CTX%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
            IF(ASSOCIATED(SOLVER_MATRICES)) THEN
              IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  SELECT CASE(SOLVER_EQUATIONS%SPARSITY_TYPE)
                  CASE(SOLVER_SPARSE_MATRICES)
                    JACOBIAN_FDCOLORING=>LINESEARCH_SOLVER%JACOBIAN_FDCOLORING
                    IF(ASSOCIATED(JACOBIAN_FDCOLORING)) THEN
                      CALL PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR(SNES,X,A,B,FLAG,JACOBIAN_FDCOLORING,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Linesearch solver FD colouring is not associated.",ERR,ERROR,*998)
                    ENDIF
                  CASE(SOLVER_FULL_MATRICES)
                    CALL PETSC_SNESDEFAULTCOMPUTEJACOBIAN(SNES,X,A,B,FLAG,CTX,ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The specified solver equations sparsity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  IF(CTX%OUTPUT_TYPE>=SOLVER_MATRIX_OUTPUT) THEN
                    CALL DISTRIBUTED_MATRIX_OVERRIDE_SET_ON(SOLVER_MATRICES%MATRICES(1)%PTR%MATRIX,A,ERR,ERROR,*999)
                    CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES_JACOBIAN_ONLY,SOLVER_MATRICES,ERR,ERROR,*998)
                    CALL DISTRIBUTED_MATRIX_OVERRIDE_SET_OFF(SOLVER_MATRICES%MATRICES(1)%PTR%MATRIX,ERR,ERROR,*999)
                  ENDIF

                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*998)
                ENDIF
              ELSE
                LOCAL_ERROR="The number of solver matrices of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                  & " is invalid. There should be 1 solver matrix."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Nonlinear solver is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Newton solver is not associated.",ERR,ERROR,*998)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Newton linesearch solver context is not associated.",ERR,ERROR,*998)
  ENDIF

  RETURN
999 CALL DISTRIBUTED_MATRIX_OVERRIDE_SET_OFF(SOLVER_MATRIX%MATRIX,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL WRITE_ERROR(ERR,ERROR,*997)
997 CALL FLAG_WARNING("Error evaluating nonlinear Jacobian.",ERR,ERROR,*996)
996 RETURN
END SUBROUTINE PROBLEM_SOLVER_JACOBIAN_FD_CALCULATE_PETSC

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the residual for a Newton like nonlinear solver
SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC(SNES,X,F,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: F !<The PETSc F Vec type
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL_VECTOR,SOLVER_VECTOR
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    NONLINEAR_SOLVER=>CTX%NONLINEAR_SOLVER
    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        SOLVER_EQUATIONS=>CTX%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
              SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
              IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
                IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                  RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                  IF(ASSOCIATED(RESIDUAL_VECTOR)) THEN
                    CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(SOLVER_VECTOR,X,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(RESIDUAL_VECTOR,F,ERR,ERROR,*999)                
                    
                    CALL PROBLEM_SOLVER_RESIDUAL_EVALUATE(CTX,ERR,ERROR,*999)
                    
                    CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(RESIDUAL_VECTOR,ERR,ERROR,*999)                
                  ELSE
                    CALL FLAG_ERROR("Residual vector is not associated.",ERR,ERROR,*997)                
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*997)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*997)
              ENDIF
            ELSE
              LOCAL_ERROR="The number of solver matrices of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                & " is invalid. There should be 1 solver matrix."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)          
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*997)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*997)
        ENDIF
!!TODO: move this to PROBLEM_SOLVER_RESIDUAL_EVALUATE or elsewhere?
        NEWTON_SOLVER%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS=NEWTON_SOLVER%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS+1
      ELSE
        CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver nonlinear solver is not associated.",ERR,ERROR,*997)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*997)
  ENDIF
  
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)  
998 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(RESIDUAL_VECTOR,DUMMY_ERR,DUMMY_ERROR,*997)
997 CALL WRITE_ERROR(ERR,ERROR,*996)
996 CALL FLAG_WARNING("Error evaluating nonlinear residual.",ERR,ERROR,*995)
995 RETURN

END SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to test convergence for a Newton like nonlinear solver
SUBROUTINE ProblemSolver_ConvergenceTestPetsc(snes,iterationNumber,xnorm,gnorm,fnorm,reason,ctx,err)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE INPUT_OUTPUT
  USE KINDS
  USE PROBLEM_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES
  USE CMISS_PETSC

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: snes !<The PETSc SNES type
  INTEGER(INTG), INTENT(INOUT) :: iterationNumber !< The current iteration (1 is the first and is before any Newton step)
  REAL(DP), INTENT(INOUT) :: xnorm !<The 2-norm of current iterate
  REAL(DP), INTENT(INOUT) :: gnorm !<The 2-norm of current step
  REAL(DP), INTENT(INOUT) :: fnorm !<The 2-norm of function
  INTEGER(INTG), INTENT(INOUT) :: reason !<The reason for convergence/divergence
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(PETSC_VEC_TYPE) :: x,f,y,w,g
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(PetscSnesLinesearchType) :: lineSearch
  REAL(DP) :: energy,normalisedEnergy
  TYPE(VARYING_STRING) :: error,localError

  IF(ASSOCIATED(ctx)) THEN
    nonlinearSolver=>CTX%NONLINEAR_SOLVER
    IF(ASSOCIATED(nonlinearSolver)) THEN
      newtonSolver=>nonlinearSolver%NEWTON_SOLVER
      IF(ASSOCIATED(newtonSolver)) THEN 
        reason=PETSC_SNES_CONVERGED_ITERATING
        SELECT CASE(newtonSolver%convergenceTestType)
        CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM) 
          IF(iterationNumber>0) THEN
            CALL Petsc_SnesLineSearchInitialise(lineSearch,err,error,*999)
            CALL Petsc_SnesGetSnesLineSearch(snes,lineSearch,err,error,*999)
            CALL PETSC_VECINITIALISE(x,err,error,*999)
            CALL PETSC_VECINITIALISE(f,err,error,*999)
            CALL PETSC_VECINITIALISE(y,err,error,*999)
            CALL PETSC_VECINITIALISE(w,err,error,*999)
            CALL PETSC_VECINITIALISE(g,err,error,*999)
            CALL Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*999)
            CALL Petsc_VecDot(y,g,energy,err,error,*999)
            IF(iterationNumber==1) THEN
              IF(ABS(energy)<ZERO_TOLERANCE) THEN
                reason=PETSC_SNES_CONVERGED_FNORM_ABS
              ELSE
                newtonSolver%convergenceTest%energyFirstIter=energy
              ENDIF
            ELSE
              normalisedEnergy=energy/newtonSolver%convergenceTest%energyFirstIter
              newtonSolver%convergenceTest%normalisedEnergy=normalisedEnergy
              IF(ABS(normalisedEnergy)<newtonSolver%ABSOLUTE_TOLERANCE*newtonSolver%SOLUTION_TOLERANCE) THEN
                reason=PETSC_SNES_CONVERGED_FNORM_ABS
              ENDIF
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Normalised energy = ",normalisedEnergy,err,error,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
            ENDIF
            CALL Petsc_SnesLineSearchFinalise(lineSearch,err,error,*999)
          ELSE
            newtonSolver%convergenceTest%energyFirstIter=0.0_DP
            newtonSolver%convergenceTest%normalisedEnergy=0.0_DP
          ENDIF
        CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
          CALL FLAG_ERROR("Differentiated ratio convergence test not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified convergence test type of "//TRIM(NUMBER_TO_VSTRING( &
            & newtonSolver%convergenceTestType,"*",err,error))//" is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver nonlinear solver is not associated.",err,error,*999)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",err,error,*999)
  ENDIF
  
  RETURN
999 CALL WRITE_ERROR(err,error,*998)
998 CALL FLAG_WARNING("Error in convergence test.",err,error,*997)
997 RETURN    

END SUBROUTINE ProblemSolver_ConvergenceTestPetsc

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to monitor a nonlinear solver
SUBROUTINE Problem_SolverNonlinearMonitorPETSC(SNES,iterationNumber,residualNorm,context,err)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES type
  INTEGER(INTG), INTENT(INOUT) :: iterationNumber !<The iteration number
  REAL(DP), INTENT(INOUT) :: residualNorm !<The residual norm
  TYPE(SOLVER_TYPE), POINTER :: context !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(SOLVER_TYPE), POINTER :: solver
  TYPE(VARYING_STRING) :: error

  IF(ASSOCIATED(context)) THEN
    nonlinearSolver=>context%NONLINEAR_SOLVER
    IF(ASSOCIATED(nonlinearSolver)) THEN
      solver=>nonlinearSolver%SOLVER
      IF(ASSOCIATED(solver)) THEN
        CALL Problem_SolverNonlinearMonitor(solver,iterationNumber,residualNorm,err,error,*999)
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver nonlinear solver is not associated.",err,error,*999)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",err,error,*999)
  ENDIF
  
  RETURN

999 CALL WRITE_ERROR(err,error,*998)
998 CALL FLAG_WARNING("Error evaluating nonlinear residual.",err,error,*997)
997 RETURN    

END SUBROUTINE Problem_SolverNonlinearMonitorPETSC

!> \file
!> $Id: analytic_analysis_routines.f90 28 2008-09-02 15:35:14Z cpb $
!> \author Ting Yu
!> \brief This module handles all Analytic solutions.
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
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
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

!>This module handles all Laplace equations routines.
MODULE ANALYTIC_ANALYSIS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ANALYTIC_ANALYSIS_EXPORT
  
  PUBLIC ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET,ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET, &
    & ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET,ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET
 
  PUBLIC ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET,ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET
  
  PUBLIC ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET,ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET
    
CONTAINS  

  !
  !================================================================================================================================
  !  

  !> Calculate the analytic analysis data.
  SUBROUTINE ANALYTIC_ANALYSIS_CALCULATE(FIELD,FILE_ID,ERR,ERROR,*)
  
    !Argument variables 
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: FIELD  
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: var_idx,comp_idx,node_idx,dev_idx,NUM_OF_NODAL_DEV, pow_idx
    TYPE(VARYING_STRING) :: STRING_DATA
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    REAL(DP), ALLOCATABLE :: VALUE_BUFFER(:)
    REAL(DP) :: RMS_PERCENT, RMS_ABSOLUTE, RMS_RELATIVE, INTEGRAL_NUM, INTEGRAL_ANA
    
    CALL ENTERS("ANALYTIC_ANALYSIS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      ALLOCATE(VALUE_BUFFER(5),STAT=ERR)
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        IF(var_idx==1) STRING_DATA="Dependent variable"
        IF(var_idx==2) STRING_DATA="Normal Derivative"
        RMS_PERCENT=0.0_DP
        RMS_ABSOLUTE=0.0_DP
        RMS_RELATIVE=0.0_DP
        INTEGRAL_NUM=0.0_DP
        INTEGRAL_ANA=0.0_DP
        CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_TRIM(STRING_DATA), ERR,ERROR,*999)
        STRING_DATA="Node #              Numerical      Analytic      % error      Absolute error Relative error"
        CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_TRIM(STRING_DATA), ERR,ERROR,*999)
        DO comp_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
          DOMAIN_NODES=>FIELD%VARIABLES(var_idx)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
          DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
            NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            DO dev_idx=1,NUM_OF_NODAL_DEV 
              CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,var_idx,comp_idx,node_idx,dev_idx,VALUE_BUFFER(1),ERR,ERROR,*999)
              CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,var_idx,comp_idx,node_idx,dev_idx,VALUE_BUFFER(2),ERR,ERROR,*999)
              CALL ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,var_idx,comp_idx,node_idx,dev_idx,VALUE_BUFFER(3),ERR,ERROR,*999)
              CALL ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,var_idx,comp_idx,node_idx,dev_idx,VALUE_BUFFER(4),ERR,ERROR,*999)
              CALL ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,var_idx,comp_idx,node_idx,dev_idx,VALUE_BUFFER(5),ERR,ERROR,*999)
              
              INTEGRAL_ANA=VALUE_BUFFER(2)**3/3+INTEGRAL_ANA
              CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER, &
                CHAR('("     '//NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR)//'",5(X,D13.4))'),'(20X,5(X,D13.4))', &
                & ERR,ERROR,*999)
            ENDDO
          ENDDO 
        ENDDO
        
        CALL ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET(FIELD,var_idx,RMS_PERCENT,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Percent) = ",RMS_PERCENT,ERR,ERROR,*999)
        CALL ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET(FIELD,var_idx,RMS_ABSOLUTE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Absolute) = ",RMS_ABSOLUTE,ERR,ERROR,*999)
        CALL ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET(FIELD,var_idx,RMS_RELATIVE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Relative) = ",RMS_RELATIVE,ERR,ERROR,*999)
        
      ENDDO
 
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
      ! Integral error
        IF(var_idx==1) THEN
          STRING_DATA="Dependent variable integral error"
        ELSE
          STRING_DATA="Flux integral error"
        END IF
        CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_TRIM(STRING_DATA), ERR,ERROR,*999)
        STRING_DATA="Node #              Numerical      Analytic      % error      Absolute error Relative error"
        CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_TRIM(STRING_DATA), ERR,ERROR,*999)
        DO pow_idx=1,2 
          CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,var_idx,pow_idx,VALUE_BUFFER(1),ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,var_idx,pow_idx,VALUE_BUFFER(2),ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET(FIELD,var_idx,pow_idx,VALUE_BUFFER(3),ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET(FIELD,var_idx,pow_idx,VALUE_BUFFER(4),ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET(FIELD,var_idx,pow_idx,VALUE_BUFFER(5),ERR,ERROR,*999)
          SELECT CASE(pow_idx)
          CASE(1)
            CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER,'("  Intgl         ",5(X,D13.4))','(20X,5(X,D13.4))', &
              & ERR,ERROR,*999)
          CASE(2)
            CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER,'("  Int^2         ",5(X,D13.4))','(20X,5(X,D13.4))', &
              & ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid power value!",ERR,ERROR,*999)     
          END SELECT
        ENDDO
      ENDDO
      IF(ALLOCATED(VALUE_BUFFER)) DEALLOCATE(VALUE_BUFFER)
    ELSE
       CALL FLAG_ERROR("The field is not associated!",ERR,ERROR,*999)     
    ENDIF
          
    CALL EXITS("ANALYTIC_ANALYSIS_CALCULATE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_CALCULATE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_CALCULATE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Export analytic information \see{ANALYTIC_ANALYSIS_ROUTINES::ANALYTIC_ANALYSIS_EXPORT}.                 
  SUBROUTINE ANALYTIC_ANALYSIS_EXPORT(FIELD,FILE_NAME, METHOD, ERR,ERROR,*)
    !Argument variables       
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field object
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
    TYPE(VARYING_STRING), INTENT(IN):: METHOD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: FILE_STATUS 
    INTEGER(INTG) :: FILE_ID

    CALL ENTERS("ANALYTIC_ANALYSIS_EXPORT", ERR,ERROR,*999)    

    IF(METHOD=="FORTRAN") THEN
       FILE_STATUS="REPLACE"
       FILE_ID=1245+FIELD%GLOBAL_NUMBER
       CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
       CALL ANALYTIC_ANALYSIS_CALCULATE(FIELD,FILE_ID,ERR,ERROR,*999)
    ELSE IF(METHOD=="MPIIO") THEN
       CALL FLAG_ERROR("MPI IO has not been implemented yet!",ERR,ERROR,*999)
    ELSE 
       CALL FLAG_ERROR("Unknown method!",ERR,ERROR,*999)   
    ENDIF   
    
    CALL EXITS("ANALYTIC_ANALYSIS_EXPORT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_EXPORT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_EXPORT")
    RETURN 1  
  END SUBROUTINE ANALYTIC_ANALYSIS_EXPORT
  
  !
  !================================================================================================================================
  !

  !>Open a file using Fortran. TODO should we use method in FIELD_IO??     
  SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<the name of file.
    TYPE(VARYING_STRING), INTENT(IN) :: FILE_STATUS !<status for opening a file
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
        
    CALL ENTERS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN",ERR,ERROR,*999)       

    !CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"OPEN FILE",ERR,ERROR,*999)

    OPEN(UNIT=FILE_ID, FILE=CHAR(FILE_NAME), STATUS=CHAR(FILE_STATUS), FORM="FORMATTED", ERR=999)   
        
    
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN
  
  !
  !================================================================================================================================
  !  

  !>Write a string using FORTRAN IO. TODO should we use FIELD_IO_FORTRAN_FILE_WRITE_STRING??    
  SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: STRING_DATA !<the string data.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING",ERR,ERROR,*999)
        
    IF(LEN_OF_DATA==0) THEN
       CALL FLAG_ERROR("leng of string is zero",ERR,ERROR,*999) 
    ENDIF
   
    WRITE(FILE_ID, "(A)") CHAR(STRING_DATA) 
    
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_WRITE_STRING
  
  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    REAL(DP), INTENT(OUT) :: VALUE !<the integral absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ABS(ANALYTIC_VALUE-NUMERICAL_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get integral analytic value for the field TODO should we use analytical formula to calculate the integration?
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    REAL(DP), INTENT(OUT) :: VALUE !<the integral analytic value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ANALYTIC_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
	        CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,comp_idx,node_idx,dev_idx,ANALYTIC_VALUE,ERR,ERROR,*999)
	        SELECT CASE(POWER)
	        CASE(1)
	          VALUE=VALUE+ANALYTIC_VALUE/2
	        CASE(2)
	          VALUE=VALUE+ANALYTIC_VALUE**2/2
	        CASE DEFAULT
	          CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
	        END SELECT
	      ENDDO
        ENDDO 
      ENDDO 
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral numerical value for the field, TODO check integral calculation
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    REAL(DP), INTENT(OUT) :: VALUE !<the integral numerical value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,VARIABLE_NUMBER,comp_idx,node_idx,dev_idx,NUMERICAL_VALUE,ERR,ERROR,*999)
            SELECT CASE(POWER)
            CASE(1)
              VALUE=VALUE+NUMERICAL_VALUE/2
            CASE(2)
              VALUE=VALUE+NUMERICAL_VALUE**2/2
            CASE DEFAULT
              CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
            END SELECT
          ENDDO
        ENDDO 
      ENDDO 
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral percentage error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    REAL(DP), INTENT(OUT) :: VALUE !<the integral percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ANALYTIC_VALUE/=0.0_DP) THEN
        VALUE=(ANALYTIC_VALUE-NUMERICAL_VALUE)/ANALYTIC_VALUE*100
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get integral relative error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    REAL(DP), INTENT(OUT) :: VALUE !<the integral relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ANALYTIC_VALUE/=-1) THEN
        VALUE=ABS((NUMERICAL_VALUE-ANALYTIC_VALUE)/(1+ANALYTIC_VALUE))
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ABS(NUMERICAL_VALUE-ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get Analytic value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the analytic value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:)
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET(FIELD,FIELD_ANALYTIC_SET_TYPE,ANALYTIC_PARAMETERS,ERR,ERROR,*999)
      IF(ASSOCIATED(ANALYTIC_PARAMETERS)) THEN      
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES
        IF(ASSOCIATED(DOMAIN_NODES)) THEN
          VALUE=ANALYTIC_PARAMETERS((VARIABLE_NUMBER-1)*DOMAIN_NODES%NUMBER_OF_NODES+NODE_NUMBER+DERIVATIVE_NUMBER-1)  
        ELSE
          CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("ANALYTIC_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF 
    
    IF(ASSOCIATED(ANALYTIC_PARAMETERS)) THEN
      NULLIFY(ANALYTIC_PARAMETERS)
    ELSE
      CALL FLAG_ERROR("ANALYTIC_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get Numerical value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the numerical value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: NUMERICAL_PARAMETERS(:)
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET(FIELD,FIELD_VALUES_SET_TYPE,NUMERICAL_PARAMETERS,ERR,ERROR,*999)
            
      IF(ASSOCIATED(NUMERICAL_PARAMETERS)) THEN      
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES
        IF(ASSOCIATED(DOMAIN_NODES)) THEN
          VALUE=NUMERICAL_PARAMETERS((VARIABLE_NUMBER-1)*DOMAIN_NODES%NUMBER_OF_NODES+NODE_NUMBER+DERIVATIVE_NUMBER-1)  
        ELSE
          CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("NUMERICAL_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF 
    
    IF(ASSOCIATED(NUMERICAL_PARAMETERS)) THEN
      NULLIFY(NUMERICAL_PARAMETERS)
    ELSE
      CALL FLAG_ERROR("NUMERICAL_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET
  
  
  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !< derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ANALYTIC_VALUE/=0.0_DP) THEN
        VALUE=(ANALYTIC_VALUE-NUMERICAL_VALUE)/ANALYTIC_VALUE*100
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET
  
  
 
  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !< derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,DERIVATIVE_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ANALYTIC_VALUE/=-1) THEN
        VALUE=ABS((NUMERICAL_VALUE-ANALYTIC_VALUE)/(1+ANALYTIC_VALUE))
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get rms percentage error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the rms percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: PERCENT_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            CALL ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,comp_idx,node_idx,dev_idx,PERCENT_ERROR_VALUE,ERR,ERROR,*999)
            VALUE=VALUE+PERCENT_ERROR_VALUE**2
          ENDDO
        ENDDO 
        VALUE=SQRT(VALUE/DOMAIN_NODES%NUMBER_OF_NODES)
      ENDDO 
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET
  
  
  !
  !================================================================================================================================
  !

  !>Get rms absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the rms absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ABSOLUTE_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx, dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            CALL ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,comp_idx,node_idx,dev_idx,ABSOLUTE_ERROR_VALUE,ERR,ERROR,*999)
            VALUE=VALUE+ABSOLUTE_ERROR_VALUE**2
          ENDDO
        ENDDO 
        VALUE=SQRT(VALUE/DOMAIN_NODES%NUMBER_OF_NODES)
      ENDDO 
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get rms relative error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the rms relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: RELATIVE_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx, dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            CALL ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,comp_idx,node_idx,dev_idx,RELATIVE_ERROR_VALUE,ERR,ERROR,*999)
            VALUE=VALUE+RELATIVE_ERROR_VALUE**2
          ENDDO
        ENDDO 
        VALUE=SQRT(VALUE/DOMAIN_NODES%NUMBER_OF_NODES)
      ENDDO 
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET

  !
  !================================================================================================================================
  !  
  
 
END MODULE ANALYTIC_ANALYSIS_ROUTINES

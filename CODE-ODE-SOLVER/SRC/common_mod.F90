!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMOUN MODULES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  in this file, there are the main and most important parameters and variables used in this !!!!!!!
!                                            program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE global_variables


 real    :: t !  temporal  variable
 real    :: dt
 parameter (dt  = 0.0001)
 real    ::  time_max
 integer :: number_iteration 


 real, allocatable,  dimension(:)   :: phi_k 
 real,  allocatable, dimension(:) :: k 
 real, allocatable, dimension(:)  :: X
 real, allocatable, dimension(:) :: F
 real, allocatable, dimension(:) :: RKf1
 real, allocatable, dimension(:) :: RKf2
 real, allocatable, dimension(:) :: RKf3
 real, allocatable, dimension(:) :: RKf4
 real, allocatable, dimension(:):: X_intermdiaire
 real, allocatable, dimension(:,:) ::  Matrix_time_evol_state
 real, allocatable, dimension(:,:) ::  Matrix_time_evol_SD

 ! nombre total  de variables;  dimension du  système incluant les variables constantes
 integer :: Nbre_total_x
 parameter ( Nbre_total_x = 20)
 ! nombre de réactions élémentaires, nbre_réact du  fichier csv, réaction_max
 integer :: nbre_react
 parameter (nbre_react = 36)
 
END MODULE global_variables



program  odesolver

use global_variables
    
implicit none
real     :: start_time, stop_time, CteW
integer ::  start_yn , purpose, jk,freq_saving
integer :: iter , i,j  ,nrun,save_step
 logical :: a,b,c,d
  character(len=70)  :: filename 
  character(len=14)  :: rnum
 
  
  ! ici 2*X+1 o X est le nombre de variables, complet (c'est Nbre_total_x= variables + perturbation externe)
   character(len=14) ,dimension(2*20+1) :: names
   a = .TRUE.
   b = .TRUE.
   c = .FALSE.
   d = .TRUE.
      write(*,*)' ***********************************************************'
      write(*,*)' ***********************************************************'
      write(*,*)' *****             WELCOME IN ODE SOLVER            ********'
      write(*,*)' **** A PROGRAM TO PERFORM DETERMINISTIC SIMULATIONS *******'
      write(*,*)' *****       OF  BIOLOGICAL  MOLECULAR NETWORKS     ********'
      write(*,*)' ***********************************************************'
      write(*,*)' ***********************************************************'
      write(*,*)''
     
     
      write(*,*)''
      write(*,*)' ****************************************************************************'
      write(*,*)' Purpose of the simulation?(1 = full simulation, 0 = Generation of SRES data)'
      write(*,*)' ****************************************************************************'
      read(*,*), purpose
     
      write(*,*)''
      write(*,*)' ***********************************************************'
      write(*,*)' Enter the duration of the simulation (in minutes):'
      write(*,*)' ***********************************************************'
      read(*,*),  time_max
      number_iteration  = floor( time_max/dt)
      if (purpose.eq. 0) freq_saving  = 200000
       if (purpose .eq. 1) freq_saving = 2000
      write(*,*)''
      write(*,*)' ********************************************************************'
      write(*,*)' With the current dt =', dt, 'the number of iterations will be ',number_iteration
       write(*,*)' the number of save point is', floor(number_iteration/(freq_saving*1.0))
      write(*,*)' ********************************************************************'
   
      write(*,*)''
      write(*,*)' ********************************************************************'
      write(*,*)' Do  you  want to  proceed?(1 = yes, 0 = no)'
      write(*,*)' ********************************************************************'
      read(*,*), start_yn
  if (start_yn == 1) then    
 
  allocate(X( Nbre_total_x),X_intermdiaire( Nbre_total_x))
  allocate(F( Nbre_total_x ))
  allocate(phi_k(nbre_react), k(nbre_react))
  allocate(RKf1( Nbre_total_x ),RKf2(Nbre_total_x),RKf3( Nbre_total_x),RKf4(Nbre_total_x))
  allocate(Matrix_time_evol_state(floor(number_iteration/(freq_saving*1.0)), Nbre_total_x+1))
  allocate(Matrix_time_evol_SD( floor(number_iteration/(freq_saving*1.0)), 2*Nbre_total_x+1))
 
  t = 0.0
  F(:) = 0.0
  phi_k(:)  =0.0
  X(:) = 0.0
  X_intermdiaire(:) = 0.0

  CALL initial_conditions()

  X_intermdiaire(:) = X(:)
  
  CALL CPU_TIME(start_time) 
  save_step = 0
  do iter =1,number_iteration
  t = t +dt
  ! Ici, on peut initialiser les constantes et perturbations externes au systme
  if (a) then
  if (t >= 350.0) then 
  X_intermdiaire(20) =10
  a = .FALSE.
  end if 
  end if
  if (b) then
  if (t >= 600.0) then 
  X_intermdiaire(20) =0
  b = .FALSE.
  end if 
  end if
  if (c) then
  if (t >= 1000.0) then 
  X_intermdiaire(20) =1000
  c = .FALSE.
  end if 
  end if
  if (d) then
  if (t >= 1500.0) then 
  X_intermdiaire(20) =1000
  d = .FALSE.
  end if 
  end if
  
  !  initialization RK4
  RKf1(:) = 0.0
  RKf2(:) = 0.0
  RKf3(:) = 0.0
  RKf4(:) = 0.0
  X(:) =  X_intermdiaire(:)

  ! step  1 RK4
   CALL phi_F_functions()
   
   RKf1(:) = dt*F(:)
   
  ! step  2 RK4
  X(:) =  X_intermdiaire(:) + 0.5*RKf1(:)
  
   CALL phi_F_functions()
   RKf2(:) = dt*F(:)
    X(:) =  X_intermdiaire(:) + 0.5*RKf2(:)
    
  ! step  3 RK4
   CALL phi_F_functions()
   RKf3(:) = dt*F(:)
    X(:) =  X_intermdiaire(:) + RKf3(:)
   
  ! step  4 RK4
   CALL  phi_F_functions()
   
   RKf4(:) = dt*F(:)
   
   X_intermdiaire(:) =  X_intermdiaire(:) +(1./6.)*( RKf1(:)+ 2.0* RKf2(:)+2.0* RKf3(:) +RKf4(:)) 
    if (purpose == 1) then
        if (mod(iter, freq_saving ).eq.0) then
         save_step = save_step+1
         Matrix_time_evol_state( save_step,1:Nbre_total_x) =  X_intermdiaire(:) 
         Matrix_time_evol_state( save_step,Nbre_total_x+1) = t
        endif 
     else
        if (mod(iter, freq_saving ).eq.0) then
         save_step = save_step+1
         Matrix_time_evol_state( save_step,1) = t
         Matrix_time_evol_state( save_step,2:Nbre_total_x+1) =  X_intermdiaire(:) 
        end if
     end if
  end do
 
  
  CALL CPU_TIME(stop_time)
  if (purpose == 1) then
   nrun  = 1
   write (rnum,'(i2.2)') nrun 
   filename = '../CODE-ODE-SOLVER-LIF/LIF-2022-ODE'
   open(1,file=filename, position="append")  
   do  i =1,save_step
    write(1,1008)( Matrix_time_evol_state(i,j) , j = 1,Nbre_total_x+1)
   end do
   !  format Nbre_total_x+1
   1008  format(21(2x, e14.8))
   close(1)
  else
  
  ! SRES generation of data

   do i=1,floor(number_iteration/(freq_saving*1.0))
     Matrix_time_evol_SD(i,1) = Matrix_time_evol_state(i,1)
   end do
   jk = 2
   do j=2,Nbre_total_x+1
    do i=1,floor(number_iteration/(freq_saving*1.0))
      Matrix_time_evol_SD(i,jk) = Matrix_time_evol_state(i,j)
      if (Matrix_time_evol_SD(i,jk) < 1.0e-15)   Matrix_time_evol_SD(i,jk) = 0.0
    end do
    jk  = jk+1
    do i=1,floor(number_iteration/(freq_saving*1.0))
      Matrix_time_evol_SD(i,jk) = 0.0
    end do
    jk  = jk+1
  end do
  write (rnum,'(i1.1)') nrun 
  filename = '../ODE_SOLVER_TF/SRES_TF'
  names(1) = 'Time'
  do i=1,2*Nbre_total_x+1
   write (rnum,'(i2.2)') i
    names(2*i) = 'species'//rnum
    names(2*i+1) = 'SD'
  end do
  open(1,file=filename, position="append")  
  write(1,*)(names(i), i=1,2*Nbre_total_x+1)
  do  i =1,save_step
    write(1,1009)( Matrix_time_evol_SD(i,j) , j = 1,2*Nbre_total_x+1)
  end do
   write(1,*)(1, i=1,Nbre_total_x)
   !  format 2*Nbre_total_x+1
   1009 format(113(2x, e14.8))
  
   close(1)
  end if
      write(*,*)''
      write(*,*)'*****************************************************************************************************************************************'
      write(*,*)'                                         The simulation is finished.'
      write(*,*)'                           Physical duration of the simulation is t=',t, 'minutes'
      write(*,*)'                         The Elapsed CPU time is ', stop_time-start_time, 'seconds' 
      write(*,*)'                          or approximately',(stop_time-start_time)/60.,'minutes'
      write(*,*)'*****************************************************************************************************************************************'
    
 deallocate (X,X_intermdiaire)
 deallocate (F)
 deallocate (phi_k,k)
 deallocate (RKf1,RKf2,RKf3,RKf4 )
 deallocate (Matrix_time_evol_state)
 deallocate (Matrix_time_evol_SD)
 else 
      write(*,*)' ********************************************************************'
      write(*,*)' Please choose another setting for the parameters'
      write(*,*)' ********************************************************************'
 end if
end


!===================================================================================================

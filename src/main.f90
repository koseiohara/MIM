program main

    use namelist     , only : namelist_init
    use timestep     , only : get_nt
    use com_var      , only : com_var_ini, com_var_end                                                        , &
                            & get_pin, get_pout, get_zd, get_rho, get_alat, get_sintbl, get_costbl, get_tantbl
    use mim_var      , only : mim_var_ini, mim_var_end
    use status_output, only : warn_open, warn_close
    use io_main      , only : files_open, files_close
    use mim          , only : mim_exec

    implicit none

    real(4) :: clock_begin
    real(4) :: clock_end
    integer :: elapse_day
    integer :: elapse_hr
    integer :: elapse_min
    integer :: elapse_sec

    integer :: nt

    write(*,*)
    write(*,'(A)') 'MIM varsion 24'
    write(*,*)

    call cpu_time(clock_begin)

    ! Read Namelist
    call namelist_init()

    ! Get Number of Time Steps
    call get_nt(nt)

    ! Allocate Constant Variables
    call com_var_ini()

    ! Set Constant Variables
    call get_pin()
    call get_pout()
    call get_zd()
    call get_rho()
    call get_alat()
    call get_sintbl()
    call get_costbl()
    call get_tantbl()

    ! Allocate The Other Variables
    call mim_var_ini()

    ! Open Alert File
    call warn_open()

    ! Open All Input and Output Files
    call files_open()

    ! Main Routine Execution
    call mim_exec(nt)
    
    ! Close All Input and Output Fiels
    call files_close()

    ! Close Alert File
    call warn_close()

    ! Deallocate The MIM-Variables
    call mim_var_end()

    ! Deallocate Constant Variables
    call com_var_end()

    call cpu_time(clock_end)

    !! ELAPSE
    elapse_sec = int(clock_end - clock_begin)
    elapse_min = elapse_sec / 60
    elapse_sec = mod(elapse_sec, 60)
    elapse_hr  = elapse_min / 60
    elapse_min = mod(elapse_min, 60)
    elapse_day = elapse_hr  / 24
    elapse_hr  = mod(elapse_hr , 24)

    write(*,*)
    write(*,'(A)') 'PROCESS COMPLETE'
    write(*,'(A,I0,"d ",I0,"hr ",I0,"min ",I0,"s")') 'EXECUTION TIME : ', elapse_day, elapse_hr, elapse_min, elapse_sec
    write(*,*)

end program main


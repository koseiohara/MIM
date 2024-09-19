module BinIO

    use params       , only : rkp, wkp, kp
    use status_output, only : nan_write

    implicit none

    private
    public :: finfo


    type finfo
        private

        integer       :: unit       ! Unit number of the file
        character(64) :: varname    ! Name of variable
        integer       :: record     ! File record
        integer       :: nx         ! Number of grids in x-direction
        integer       :: ny         ! Number of grids in y-direction
        integer       :: nz         ! Number of grids in z-direction
        logical       :: xrev       ! X Reverse
        logical       :: yrev       ! Y Reverse
        logical       :: zrev       ! Z REverse

        contains

        procedure     :: fclose
        procedure     :: fread
        procedure     :: fwrite

    end type finfo


    interface finfo
        module procedure init
    end interface finfo


    contains


    !!!!! CONSTRUCTOR
    function init(fname, action, nx, ny, nz, xrev, yrev, zrev, varname, endian) result(self)
        type(finfo) :: self

        character(*), intent(in)  :: fname
        character(*), intent(in)  :: action
        integer     , intent(in)  :: nx
        integer     , intent(in)  :: ny
        integer     , intent(in)  :: nz
        logical     , intent(in)  :: xrev
        logical     , intent(in)  :: yrev
        logical     , intent(in)  :: zrev
        character(*), intent(in)  :: varname
        character(*), intent(in)  :: endian

        integer        :: recl
        integer        :: iostat
        character(128) :: iomsg

        if (trim(action) == 'READ') then
            recl = rkp*nx*ny*nz
        else if (trim(action) == 'WRITE') then
            recl = wkp*nx*ny*nz
        else
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid action name in fopen()'
            write(0,'(A)') 'Variable : ' // trim(varname)
            write(0,'(A)') 'Input action name : ' // trim(action)
            ERROR STOP
        endif

        open(NEWUNIT=self%unit    , &
           & FILE   =fname        , &
           & ACTION =action       , &
           & FORM   ='UNFORMATTED', &
           & ACCESS ='DIRECT'     , &
           & RECL   =recl         , &
           & IOSTAT =iostat       , &
           & IOMSG  =iomsg        , &
           & CONVERT=endian         )

        if (iostat /= 0) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Failed to open ' // trim(fname)
            write(0,'(A)') trim(iomsg)
            ERROR sTOP
        endif

        self%varname = varname
        self%record  = 1
        self%nx      = nx
        self%ny      = ny
        self%nz      = nz
        self%xrev    = xrev
        self%yrev    = yrev
        self%zrev    = zrev

    end function init


    subroutine fclose(self)
        class(finfo) :: self
        
        logical :: opened

        INQUIRE(self%unit, OPENED=opened)

        if (opened) then
            close(self%unit)

            self%unit = -1
            self%varname = 'ERROR (closed file)'
            self%record  = 0
            self%nx      = 0
            self%ny      = 0
            self%nz      = 0
        endif

    end subroutine fclose


    subroutine fread(self, rdata)
        class(finfo), intent(inout) :: self
        real(rkp)   , intent(out)   :: rdata(self%nx,self%ny,self%nz)

        !real(rkp) :: reader(self%nx,self%ny,self%nz)

        integer :: xstart
        integer :: xend
        integer :: xstep
        integer :: ystart
        integer :: yend
        integer :: ystep
        integer :: zstart
        integer :: zend
        integer :: zstep

        call get_order(self%xrev, &  !! IN
                     & self%nx  , &  !! IN
                     & xstart   , &  !! OUT
                     & xend     , &  !! OUT
                     & xstep      )  !! OUT

        call get_order(.NOT. self%yrev, &  !! IN
                     & self%ny        , &  !! IN
                     & ystart         , &  !! OUT
                     & yend           , &  !! OUT
                     & ystep            )  !! OUT

        call get_order(.NOT. self%zrev, &  !! IN
                     & self%nz        , &  !! IN
                     & zstart         , &  !! OUT
                     & zend           , &  !! OUT
                     & zstep            )  !! OUT

        !write(*,*) xstart, xend, xstep, ystart, yend, ystep, zstart, zend, zstep
        read(self%unit,rec=self%record) rdata(xstart:xend:xstep,ystart:yend:ystep,zstart:zend:zstep)
        !rdata(1:self%nx,1:self%ny,1:self%nz) = real(reader(1:self%nx,1:self%ny,1:self%nz), kind=kp)

        self%record = self%record + 1

    end subroutine fread


    subroutine fwrite(self, tt, var, wdata)
        class(finfo), intent(inout) :: self
        integer     , intent(in)    :: tt
        character(*), intent(in)    :: var
        real(kp)    , intent(in)    :: wdata(self%nx,self%ny,self%nz)

        real(wkp) :: writer(self%nx,self%ny,self%nz)

        integer :: xstart
        integer :: xend
        integer :: xstep
        integer :: ystart
        integer :: yend
        integer :: ystep
        integer :: zstart
        integer :: zend
        integer :: zstep

        call get_order(self%xrev, &  !! IN
                     & self%nx  , &  !! IN
                     & xstart   , &  !! OUT
                     & xend     , &  !! OUT
                     & xstep      )  !! OUT

        call get_order(self%yrev, &  !! IN
                     & self%ny  , &  !! IN
                     & ystart   , &  !! OUT
                     & yend     , &  !! OUT
                     & ystep      )  !! OUT

        call get_order(self%zrev, &  !! IN
                     & self%nz  , &  !! IN
                     & zstart   , &  !! OUT
                     & zend     , &  !! OUT
                     & zstep      )  !! OUT

        call nan_write(tt                                  , &  !! IN
                     & self%nx                             , &  !! IN
                     & self%ny                             , &  !! IN
                     & self%nz                             , &  !! IN
                     & wdata(1:self%nx,1:self%ny,1:self%nz), &  !! IN
                     & var                                   )  !! IN

        writer(1:self%nx,1:self%ny,1:self%nz) = real(wdata(1:self%nx,1:self%ny,1:self%nz), kind=kp)
        write(self%unit,rec=self%record) writer(xstart:xend:xstep,ystart:yend:ystep,zstart:zend:zstep)

        self%record = self%record + 1

    end subroutine fwrite


    subroutine get_order(rev, n, istart, iend, istep)
        logical, intent(in)  :: rev         ! True (reverse) or False (normal)
        integer, intent(in)  :: n           ! Number of grids
        integer, intent(out) :: istart      ! Start index
        integer, intent(out) :: iend        ! End index
        integer, intent(out) :: istep       ! Index step

        if (rev) then
            istart = n
            iend   = 1
            istep  = -1
        else
            istart = 1
            iend   = n
            istep  = 1
        endif

    end subroutine get_order


end module BinIO


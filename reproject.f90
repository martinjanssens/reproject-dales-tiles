program reproject

  implicit none

  ! TODO:
  ! - Will not work if there are hetero fields
  ! - Will not work if isurf==1

  ! INPUTS SET BY USER
  character(256) :: inpath = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/code/dalesruns/bomex'
  character(256) :: outpath = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/code/dalesruns/bomex_test'
  ! integer :: ierr

  ! Could read this from namoptions.001, but then you'll have to parse all the inputs (and we don't want that)
  character(50) :: name = 'initd000h10mx000y000.001'
  integer :: nsv = 0
  integer :: nprocx = 2
  integer :: nprocy = 2
  integer :: itot = 64
  integer :: jtot = 64
  integer :: kmax = 80
  integer :: iadv_mom = 2
  integer :: iadv_tke = 2
  integer :: iadv_thl = 2
  integer :: iadv_qt = 2
  integer :: iadv_sv(100) = -1

  ! DECLARATIONS
  integer, parameter :: ifinput = 1
  integer, parameter :: ifoutput = 2
  integer, parameter :: longint = 8
  integer :: procx, procy
  integer :: imax, jmax, i1, j1, k1, i2, j2, k2, ih, jh, kh, i, j, k, n
  integer :: advarr(4)
  character(8) :: cmyid

  real, allocatable :: u0(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0h(:,:,:)
  real, allocatable :: e120(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdz(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekm(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0(:,:,:) 		!<   temperature at full level
  real, allocatable :: esl(:,:,:)
  real, allocatable :: qvsl(:,:,:)
  real, allocatable :: qvsi(:,:,:)
  real, allocatable :: ustar (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlflux (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtflux  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldz(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdz (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: presf(:)         !<   hydrostatic pressure at full level
  real, allocatable :: presh(:)         !<   hydrostatic pressure at half level
  real, allocatable :: initial_presf(:) !<   initial hydrostatic pressure at full level
  real, allocatable :: initial_presh(:) !<   initial hydrostatic pressure at half level
  real              :: ps               !<  Surface pressure [Pa]
  real              :: thls             !<  Surface liquid water potential temperature [K]
  real              :: qts              !<  Surface specific humidity [kg/kg]
  real              :: thvs             !<  Surface virtual temperature [K]
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real              :: dtheta           !<     * applied gradient of theta at top of model
  real              :: dqt              !<     * applied gradient of qt at top of model
  integer(kind=longint) :: timee        !<     * elapsed time since the "cold" start
  integer(kind=longint) :: dt           !<     * time integration interval
  real              :: tres
  real, allocatable :: obl(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskin(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskin(:,:)       !<  Skin specific humidity [kg/kg]

  ! Radiation
  integer (kind=longint) :: tnext_radiation  !<  time of the first upcoming call of the radiation scheme
  real, allocatable :: thlprad(:,:,:)!<   the radiative tendencies
  real, allocatable :: swd(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdir(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdif(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwc(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swu(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwd(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwu(:,:,:)    !<   longwave upward radiative flux
!
  real, allocatable :: swdca(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swuca(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdca(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwuca(:,:,:)  !<  clear air longwave upward radiative flux

  real, allocatable :: SW_up_TOA(:,:), SW_dn_TOA(:,:), LW_up_TOA(:,:), LW_dn_TOA(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOA(:,:), SW_dn_ca_TOA(:,:), LW_up_ca_TOA(:,:), LW_dn_ca_TOA(:,:)

  ! Scalars
  real, allocatable :: sv0(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svflux  (:,:,:)  !<  Kinematic scalar flux [- m/s]
  real, allocatable :: dsv(:)           !<     * applied gradient of sv(n) at top of model

  ! Get ghost cells for input simulation
  imax = itot/nprocx
  jmax = jtot/nprocy
  i1=imax+1
  j1=jmax+1
  k1=kmax+1
  k2=kmax+2
  i2=imax+2
  j2=jmax+2
  advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
  if     (any(advarr==6).or.any(iadv_sv(1:nsv)==6)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==62).or.any(iadv_sv(1:nsv)==62)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==5).or.any(iadv_sv(1:nsv)==5)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==52).or.any(iadv_sv(1:nsv)==52)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==55).or.any(iadv_sv(1:nsv)==55)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==555).or.any(iadv_sv(1:nsv)==555)) then
    ih = 3
    jh = 3
    kh = 1
  elseif (any(advarr==7).or.any(iadv_sv(1:nsv)==7)) then
    ih = 2
    jh = 2
    kh = 1
  elseif (any(advarr==2).or.any(iadv_sv(1:nsv)==2)) then
    ih = 1
    jh = 1
    kh = 1
  end if

  ! DALES has already initialised all fields with their correct sizes before reading from files and assigning them -> need to do this too
  allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(thl0 (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qt0  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0h  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(e120 (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(dthvdz(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(tmp0  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(esl (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsl(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsi(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ustar   (i2,j2))
  allocate(thlflux (i2,j2))
  allocate(qtflux  (i2,j2))
  allocate(dqtdz   (i2,j2))
  allocate(dthldz  (i2,j2))
  allocate(presf        (k1))
  allocate(presh        (k1))
  allocate(initial_presf(k1))
  allocate(initial_presh(k1))
  allocate(obl(i2,j2))
  allocate(tskin(i2,j2))
  allocate(qskin(i2,j2))
  allocate(thlprad   (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swd       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swu       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwd       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwu       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swuca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwuca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdir     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdif     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwc       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(SW_up_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_dn_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_up_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_dn_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_up_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(SW_dn_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_up_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_dn_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(sv0  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  allocate(svflux  (i2,j2,nsv))
  allocate(dsv(nsv))

  ! Loop over processors
  do procx = 0, nprocx-1
  	do procy = 0, nprocy-1

      ! Get current processor
      write(cmyid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', procy

	  ! Read restartfiles
	  name(5:5) = 'd'
	  name(13:20) = cmyid
	  write(6,*) 'loading ',name
	  open(unit=ifinput,file=trim(inpath)//'/'//name,form='unformatted', status='old')

	  read(ifinput)  (((u0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((v0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((w0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((thl0  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((qt0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((ql0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((ql0h  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((e120  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((dthvdz(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((ekm   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((ekh   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((tmp0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((esl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((qvsl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((qvsi   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)   ((ustar (i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((thlflux (i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((qtflux  (i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((dthldz(i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((dqtdz (i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)  (  presf (    k)                            ,k=1,k1)
	  read(ifinput)  (  presh (    k)                            ,k=1,k1)
	  read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
	  read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
	  read(ifinput)  ps,thls,qts,thvs,oblav
	  read(ifinput)  dtheta,dqt,timee,dt,tres
	  read(ifinput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
	  read(ifinput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )

	!!!!! radiation quantities
	  read(ifinput)  tnext_radiation
	  read(ifinput)  (((thlprad (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((lwd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((lwu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((lwdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((lwuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swdir   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((swdif   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  read(ifinput)  (((lwc     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

	  read(ifinput)  ((SW_up_TOA    (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((SW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((LW_up_TOA    (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((LW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((SW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((SW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((LW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  read(ifinput)  ((LW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)

	  close(ifinput)

	  if (nsv>0) then
	    name(5:5) = 's'
	    write(6,*) 'loading ',name
	    open(unit=ifinput,file=trim(outpath)//'/'//name,form='unformatted')
	    read(ifinput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
	    read(ifinput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
	    read(ifinput) (dsv(n),n=1,nsv)
	    read(ifinput)  timee
	    close(ifinput)
	  end if

	  ! And write
	  name(13:20)= cmyid
	  write(6,*) 'writing ',name
	  open(ifoutput,file=trim(outpath)//'/'//name,form='unformatted',status='replace')

	  write(ifoutput)  (((u0 (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((v0 (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((w0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((thl0  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((qt0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((ql0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((ql0h  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((e120  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((dthvdz(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((ekm   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((ekh   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((tmp0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((esl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((qvsl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((qvsi   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)   ((ustar (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((thlflux (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((qtflux  (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((dthldz(i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((dqtdz (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  (  presf (    k)                            ,k=1,k1)
	  write(ifoutput)  (  presh (    k)                            ,k=1,k1)
	  write(ifoutput)  (  initial_presf (    k)                            ,k=1,k1)
	  write(ifoutput)  (  initial_presh (    k)                            ,k=1,k1)
	  write(ifoutput)  ps,thls,qts,thvs,oblav
	  write(ifoutput)  dtheta,dqt,timee,  dt,tres
	  write(ifoutput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )

	!!!!! radiation quantities
	  write(ifoutput)  tnext_radiation
	  write(ifoutput)  (((thlprad (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((lwd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((lwu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((lwdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((lwuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swdir   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((swdif   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	  write(ifoutput)  (((lwc     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

	  write(ifoutput)  ((SW_up_TOA    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_up_TOA    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)

	  close (ifoutput)

	  if (nsv>0) then
	    name(5:5)='s'
	    open  (ifoutput,file=trim(outpath)//'/'//name,form='unformatted')
	    write(ifoutput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
	    write(ifoutput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
	    write(ifoutput) (dsv(n),n=1,nsv)
	    write(ifoutput)  timee
	    close (ifoutput)
	  end if

	end do ! nprocy
  end do ! nprocx

  deallocate(u0, v0, w0, thl0, qt0, ql0, ql0h, e120, dthvdz, ekm, ekh, tmp0, esl, qvsl, qvsi)
  deallocate(ustar, thlflux, qtflux, dthldz, dqtdz, presf, presh, initial_presf, initial_presh)
  deallocate(obl, tskin, qskin, thlprad, swd, swu, lwd, lwu, swdca, swuca)
  deallocate(lwdca, lwuca, swdir, swdif, lwc, SW_up_TOA, SW_dn_TOA, LW_up_TOA, LW_dn_TOA)
  deallocate(SW_up_ca_TOA, SW_dn_ca_TOA, LW_up_ca_TOA, LW_dn_ca_TOA)
  deallocate(sv0, svflux, dsv)

end program reproject
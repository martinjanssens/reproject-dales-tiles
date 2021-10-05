program reproject

  implicit none

  ! TODO:
  ! - Will not work if there are hetero fields
  ! - Will not work if isurf==1
  ! - Assumes periodic BCs

  ! INPUTS SET BY USER
  character(256) :: inpath = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/code/dalesruns/bomex'
  character(256) :: outpath = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/code/dalesruns/bomex_test'
  ! integer :: ierr

  ! SETTINGS FROM INPUT SIMULATION
  ! Could read this from namoptions.001, but then you'll have to parse all the inputs (and we don't want that)
  character(50) :: name = 'initd000h10mx000y000.001'
  integer :: nsv = 0
  integer :: nprocx = 4
  integer :: nprocy = 2
  integer :: itot = 64
  integer :: jtot = 64
  integer :: kmax = 80
  integer :: iadv_mom = 2
  integer :: iadv_tke = 2
  integer :: iadv_thl = 2
  integer :: iadv_qt = 2
  integer :: iadv_sv(100) = -1

  ! SETTINGS FOR DESIRED OUTPUT SIMULATION
  integer :: itoto = 32
  integer :: jtoto = 32
  integer :: kmaxo = 80
  integer :: iadv_momo = 2
  integer :: iadv_tkeo = 2
  integer :: iadv_thlo = 2
  integer :: iadv_qto = 2
  integer :: iadv_svo(100) = -1  

  !!
  !! DECLARATIONS
  !!

  !! INTEGERS, CHARACTERS, REALS, SWITCHES
  integer, parameter :: ifinput = 1
  integer, parameter :: ifoutput = 2
  integer, parameter :: longint = 8
  integer :: procx, procy
  integer :: i, j, k, n, ii, ij
  integer :: imax, jmax, i1, j1, k1, i2, j2, k2, ih, jh, kh
  integer :: imaxo, jmaxo, i1o, j1o, i2o, j2o, iho, jho
  integer :: ncoarse
  integer :: advarr(4)
  character(8) :: cmyid
  character(8) :: ceid
  character(8) :: cwid
  character(8) :: cnid
  character(8) :: csid
  real    :: ncoarse2i
  logical :: ladv
  logical :: lcoarse

  !! INPUT ALLOCATABLE FIELDS

  !! CENTRAL PROCESSOR
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
  integer (kind=longint) :: tnext_radiation  !<  time of the first upcoming call of the radiation scheme
  real, allocatable :: thlprad(:,:,:)!<   the radiative tendencies
  real, allocatable :: swd(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdir(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdif(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwc(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swu(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwd(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwu(:,:,:)    !<   longwave upward radiative flux
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

  !! EAST PROCESSOR
  real, allocatable :: u0e(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0e(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0e(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0e(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0e(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0e(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0he(:,:,:)
  real, allocatable :: e120e(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdze(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekme(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekhe(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0e(:,:,:) 		!<   temperature at full level
  real, allocatable :: esle(:,:,:)
  real, allocatable :: qvsle(:,:,:)
  real, allocatable :: qvsie(:,:,:)
  real, allocatable :: ustare (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlfluxe (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtfluxe  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldze(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdze (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: oble(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskine(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskine(:,:)       !<  Skin specific humidity [kg/kg]
  real, allocatable :: thlprade(:,:,:)!<   the radiative tendencies
  real, allocatable :: swde(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdire(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdife(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwce(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swue(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwde(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwue(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swdcae(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swucae(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdcae(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwucae(:,:,:)  !<  clear air longwave upward radiative flux
  real, allocatable :: SW_up_TOAe(:,:), SW_dn_TOAe(:,:), LW_up_TOAe(:,:), LW_dn_TOAe(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOAe(:,:), SW_dn_ca_TOAe(:,:), LW_up_ca_TOAe(:,:), LW_dn_ca_TOAe(:,:)
  ! Scalars
  real, allocatable :: sv0e(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svfluxe  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  !! WEST PROCESSOR
  real, allocatable :: u0w(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0w(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0w(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0w(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0w(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0w(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0hw(:,:,:)
  real, allocatable :: e120w(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdzw(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekmw(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekhw(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0w(:,:,:) 		!<   temperature at full level
  real, allocatable :: eslw(:,:,:)
  real, allocatable :: qvslw(:,:,:)
  real, allocatable :: qvsiw(:,:,:)
  real, allocatable :: ustarw (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlfluxw (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtfluxw  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldzw(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdzw (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: oblw(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskinw(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskinw(:,:)       !<  Skin specific humidity [kg/kg]
  real, allocatable :: thlpradw(:,:,:)!<   the radiative tendencies
  real, allocatable :: swdw(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdirw(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdifw(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwcw(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swuw(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwdw(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwuw(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swdcaw(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swucaw(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdcaw(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwucaw(:,:,:)  !<  clear air longwave upward radiative flux
  real, allocatable :: SW_up_TOAw(:,:), SW_dn_TOAw(:,:), LW_up_TOAw(:,:), LW_dn_TOAw(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOAw(:,:), SW_dn_ca_TOAw(:,:), LW_up_ca_TOAw(:,:), LW_dn_ca_TOAw(:,:)
  ! Scalars
  real, allocatable :: sv0w(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svfluxw  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  !! NORTH PROCESSOR
  real, allocatable :: u0n(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0n(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0n(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0n(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0n(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0n(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0hn(:,:,:)
  real, allocatable :: e120n(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdzn(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekmn(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekhn(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0n(:,:,:) 		!<   temperature at full level
  real, allocatable :: esln(:,:,:)
  real, allocatable :: qvsln(:,:,:)
  real, allocatable :: qvsin(:,:,:)
  real, allocatable :: ustarn (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlfluxn (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtfluxn  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldzn(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdzn (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: obln(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskinn(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskinn(:,:)       !<  Skin specific humidity [kg/kg]
  real, allocatable :: thlpradn(:,:,:)!<   the radiative tendencies
  real, allocatable :: swdn(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdirn(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdifn(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwcn(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swun(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwdn(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwun(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swdcan(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swucan(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdcan(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwucan(:,:,:)  !<  clear air longwave upward radiative flux
  real, allocatable :: SW_up_TOAn(:,:), SW_dn_TOAn(:,:), LW_up_TOAn(:,:), LW_dn_TOAn(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOAn(:,:), SW_dn_ca_TOAn(:,:), LW_up_ca_TOAn(:,:), LW_dn_ca_TOAn(:,:)
  ! Scalars
  real, allocatable :: sv0n(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svfluxn  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  !! SOUTH PROCESSOR
  real, allocatable :: u0s(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0s(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0s(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0s(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0s(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0s(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0hs(:,:,:)
  real, allocatable :: e120s(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdzs(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekms(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekhs(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0s(:,:,:) 		!<   temperature at full level
  real, allocatable :: esls(:,:,:)
  real, allocatable :: qvsls(:,:,:)
  real, allocatable :: qvsis(:,:,:)
  real, allocatable :: ustars (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlfluxs (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtfluxs  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldzs(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdzs (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: obls(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskins(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskins(:,:)       !<  Skin specific humidity [kg/kg]
  real, allocatable :: thlprads(:,:,:)!<   the radiative tendencies
  real, allocatable :: swds(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdirs(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdifs(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwcs(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swus(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwds(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwus(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swdcas(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swucas(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdcas(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwucas(:,:,:)  !<  clear air longwave upward radiative flux
  real, allocatable :: SW_up_TOAs(:,:), SW_dn_TOAs(:,:), LW_up_TOAs(:,:), LW_dn_TOAs(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOAs(:,:), SW_dn_ca_TOAs(:,:), LW_up_ca_TOAs(:,:), LW_dn_ca_TOAs(:,:)
  ! Scalars
  real, allocatable :: sv0s(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svfluxs (:,:,:)  !<  Kinematic scalar flux [- m/s]

  !! OUTPUT
  real, allocatable :: u0o(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0o(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0o(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0o(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0o(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0o(:,:,:)  		!<   liquid water content
  real, allocatable :: ql0ho(:,:,:)
  real, allocatable :: e120o(:,:,:)      !<   square root of turb. kin. energy at time step t
  real, allocatable :: dthvdzo(:,:,:)	!<   theta_v at half level
  real, allocatable :: ekmo(:,:,:)   	!< k-coefficient for momentum
  real, allocatable :: ekho(:,:,:)  		!< k-coefficient for heat and q_tot
  real, allocatable :: tmp0o(:,:,:) 		!<   temperature at full level
  real, allocatable :: eslo(:,:,:)
  real, allocatable :: qvslo(:,:,:)
  real, allocatable :: qvsio(:,:,:)
  real, allocatable :: ustaro (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlfluxo (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtfluxo  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: dthldzo(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]
  real, allocatable :: dqtdzo (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: oblo(:,:)         !<  Obukhov length [m]
  real, allocatable :: tskino(:,:)       !<  Skin temperature [K]
  real, allocatable :: qskino(:,:)       !<  Skin specific humidity [kg/kg]
  real, allocatable :: thlprado(:,:,:)!<   the radiative tendencies
  real, allocatable :: swdo(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swdiro(:,:,:)  !<   Direct shortwave downward radiative flux
  real, allocatable :: swdifo(:,:,:)  !<   Difuse shortwave downward radiative flux
  real, allocatable :: lwco(:,:,:)    !<   Liquid water content calculated in rrtmg
  real, allocatable :: swuo(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwdo(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwuo(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swdcao(:,:,:)  !<  clear air shortwave downward radiative flux
  real, allocatable :: swucao(:,:,:)  !<  clear air shortwave upward radiative flux
  real, allocatable :: lwdcao(:,:,:)  !<  clear air longwave downward radiative flux
  real, allocatable :: lwucao(:,:,:)  !<  clear air longwave upward radiative flux
  real, allocatable :: SW_up_TOAo(:,:), SW_dn_TOAo(:,:), LW_up_TOAo(:,:), LW_dn_TOAo(:,:) !< Top of the atmosphere radiative fluxes
  real, allocatable :: SW_up_ca_TOAo(:,:), SW_dn_ca_TOAo(:,:), LW_up_ca_TOAo(:,:), LW_dn_ca_TOAo(:,:)
  ! Scalars
  real, allocatable :: sv0o(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svfluxo  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  !!
  !! GET GHOST CELLS
  !!

  !! INPUT SIMULATION
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

  ! OUTPUT SIMULATION
  imaxo = itoto/nprocx
  jmaxo = jtoto/nprocy
  i1o=imaxo+1
  j1o=jmaxo+1
  i2o=imaxo+2
  j2o=jmaxo+2
  advarr = (/iadv_momo,iadv_tkeo,iadv_thlo,iadv_qto/)
  if     (any(advarr==6).or.any(iadv_svo(1:nsv)==6)) then
    iho = 3
    jho = 3
  elseif (any(advarr==62).or.any(iadv_svo(1:nsv)==62)) then
    iho = 3
    jho = 3
  elseif (any(advarr==5).or.any(iadv_svo(1:nsv)==5)) then
    iho = 3
    jho = 3
  elseif (any(advarr==52).or.any(iadv_svo(1:nsv)==52)) then
    iho = 3
    jho = 3
  elseif (any(advarr==55).or.any(iadv_svo(1:nsv)==55)) then
    iho = 3
    jho = 3
  elseif (any(advarr==555).or.any(iadv_svo(1:nsv)==555)) then
    iho = 3
    jho = 3
  elseif (any(advarr==7).or.any(iadv_svo(1:nsv)==7)) then
    iho = 2
    jho = 2
  elseif (any(advarr==2).or.any(iadv_svo(1:nsv)==2)) then
    iho = 1
    jho = 1
  end if

  ! SET ADVECTION SWITCH
  if (iho.ne.ih .or. jho.ne.jh) then
  	ladv=.true.
  	print *, 'Number of output ghost cells:', iho
  	print *, 'Number of input ghost cells:', ih
  	print *, 'Running advection scheme switch routine'
  else
  	ladv=.false.
  end if

  ! SET COARSE-GRAINING SWITCH
  ! FIXME not tested whether you could do both advection and mesh switch
  ! FIXME what happens to integer division?
  ncoarse = itot / itoto
  ncoarse2i = 1./(ncoarse*ncoarse)
  if (ncoarse > 1) then
  	lcoarse = .true.
  	print *, 'Number of cells per input processor is ', ncoarse, 'times larger than output, running coarse graining'
  else
  	lcoarse = .false.
  end if

  !!
  !! INITIALISE FIELDS
  !!

  !! CENTRAL PROCESSOR
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

  !! EAST PROCESSOR
  allocate(u0e   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(v0e   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(w0e   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(thl0e (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qt0e  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0e   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0he  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(e120e (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(dthvdze(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekme(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekhe(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(tmp0e  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(esle (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsle(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsie(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ustare   (i2,j2))
  allocate(thlfluxe (i2,j2))
  allocate(qtfluxe  (i2,j2))
  allocate(dqtdze   (i2,j2))
  allocate(dthldze  (i2,j2))
  allocate(oble(i2,j2))
  allocate(tskine(i2,j2))
  allocate(qskine(i2,j2))
  allocate(thlprade   (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swde       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swue       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwde       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwue       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdcae     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swucae     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdcae     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwucae     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdire     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdife     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwce       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(SW_up_TOAe (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_dn_TOAe (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_up_TOAe (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_dn_TOAe (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_up_ca_TOAe(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(SW_dn_ca_TOAe(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_up_ca_TOAe(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_dn_ca_TOAe(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(sv0e  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  allocate(svfluxe  (i2,j2,nsv))

  !! WEST PROCESSOR
  allocate(u0w   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(v0w   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(w0w   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(thl0w (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qt0w  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0w   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0hw  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(e120w (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(dthvdzw(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekmw(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekhw(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(tmp0w  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(eslw (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvslw(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsiw(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ustarw   (i2,j2))
  allocate(thlfluxw (i2,j2))
  allocate(qtfluxw  (i2,j2))
  allocate(dqtdzw   (i2,j2))
  allocate(dthldzw  (i2,j2))
  allocate(oblw(i2,j2))
  allocate(tskinw(i2,j2))
  allocate(qskinw(i2,j2))
  allocate(thlpradw   (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdw       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swuw       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdw       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwuw       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdcaw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swucaw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdcaw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwucaw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdirw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdifw     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwcw       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(SW_up_TOAw (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_dn_TOAw (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_up_TOAw (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_dn_TOAw (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_up_ca_TOAw(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(SW_dn_ca_TOAw(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_up_ca_TOAw(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_dn_ca_TOAw(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(sv0w  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  allocate(svfluxw  (i2,j2,nsv))

  !! NORTH PROCESSOR
  allocate(u0n   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(v0n   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(w0n   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(thl0n (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qt0n  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0n   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0hn  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(e120n (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(dthvdzn(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekmn(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekhn(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(tmp0n  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(esln (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsln(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsin(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ustarn   (i2,j2))
  allocate(thlfluxn (i2,j2))
  allocate(qtfluxn  (i2,j2))
  allocate(dqtdzn   (i2,j2))
  allocate(dthldzn  (i2,j2))
  allocate(obln(i2,j2))
  allocate(tskinn(i2,j2))
  allocate(qskinn(i2,j2))
  allocate(thlpradn   (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdn       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swun       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdn       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwun       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdcan     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swucan     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdcan     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwucan     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdirn     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdifn     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwcn       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(SW_up_TOAn (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_dn_TOAn (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_up_TOAn (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_dn_TOAn (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_up_ca_TOAn(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(SW_dn_ca_TOAn(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_up_ca_TOAn(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_dn_ca_TOAn(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(sv0n  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  allocate(svfluxn  (i2,j2,nsv))

  !! SOUTH PROCESSOR
  allocate(u0s   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(v0s   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(w0s   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(thl0s (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qt0s  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0s   (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ql0hs  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(e120s (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(dthvdzs(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekms(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ekhs(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(tmp0s  (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(esls (2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsls(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(qvsis(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(ustars   (i2,j2))
  allocate(thlfluxs (i2,j2))
  allocate(qtfluxs  (i2,j2))
  allocate(dqtdzs   (i2,j2))
  allocate(dthldzs  (i2,j2))
  allocate(obls(i2,j2))
  allocate(tskins(i2,j2))
  allocate(qskins(i2,j2))
  allocate(thlprads   (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swds       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swus       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwds       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwus       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdcas     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swucas     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwdcas     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwucas     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdirs     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(swdifs     (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(lwcs       (2-ih:i1+ih,2-jh:j1+jh,k1) )
  allocate(SW_up_TOAs (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_dn_TOAs (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_up_TOAs (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(LW_dn_TOAs (2-ih:i1+ih,2-jh:j1+jh)    )
  allocate(SW_up_ca_TOAs(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(SW_dn_ca_TOAs(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_up_ca_TOAs(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(LW_dn_ca_TOAs(2-ih:i1+ih,2-jh:j1+jh)  )
  allocate(sv0s  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  allocate(svfluxs  (i2,j2,nsv))

  !! OUTPUT PROCESSOR
  allocate(u0o    (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(v0o    (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(w0o    (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(thl0o  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(qt0o   (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(ql0o   (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(ql0ho  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(e120o  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(dthvdzo(2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(ekmo   (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(ekho   (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(tmp0o  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(eslo   (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(qvslo  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(qvsio  (2-iho:i1o+iho,2-jho:j1o+jho,k1))
  allocate(ustaro   (i2o,j2o))
  allocate(thlfluxo (i2o,j2o))
  allocate(qtfluxo  (i2o,j2o))
  allocate(dqtdzo   (i2o,j2o))
  allocate(dthldzo  (i2o,j2o))
  allocate(oblo(i2o,j2o))
  allocate(tskino(i2o,j2o))
  allocate(qskino(i2o,j2o))
  allocate(thlprado   (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swdo       (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swuo       (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(lwdo       (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(lwuo       (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swdcao     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swucao     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(lwdcao     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(lwucao     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swdiro     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(swdifo     (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(lwco       (2-iho:i1o+iho,2-jho:j1o+jho,k1) )
  allocate(SW_up_TOAo (2-iho:i1o+iho,2-jho:j1o+jho)    )
  allocate(SW_dn_TOAo (2-iho:i1o+iho,2-jho:j1o+jho)    )
  allocate(LW_up_TOAo (2-iho:i1o+iho,2-jho:j1o+jho)    )
  allocate(LW_dn_TOAo (2-iho:i1o+iho,2-jho:j1o+jho)    )
  allocate(SW_up_ca_TOAo(2-iho:i1o+iho,2-jho:j1o+jho)  )
  allocate(SW_dn_ca_TOAo(2-iho:i1o+iho,2-jho:j1o+jho)  )
  allocate(LW_up_ca_TOAo(2-iho:i1o+iho,2-jho:j1o+jho)  )
  allocate(LW_dn_ca_TOAo(2-iho:i1o+iho,2-jho:j1o+jho)  )
  allocate(sv0o  (2-iho:i1o+iho,2-jho:j1o+jho,k1,nsv))
  allocate(svfluxo  (i2o,j2o,nsv))

  !!
  !! LOOP OVER TILES
  !!
  do procx = 0, nprocx-1
  	do procy = 0, nprocy-1

      !! Get processor ids -> FIXME assumes periodic BCs

      ! Central processor
      write(cmyid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', procy ! Current
      
      ! East
      if (procx < nprocx-1) then
        write(ceid,'(a,i3.3,a,i3.3)') 'x', procx+1, 'y', procy
      else
      	write(ceid,'(a,i3.3,a,i3.3)') 'x', 0, 'y', procy
      end if

      ! West
      if (procx > 0) then
        write(cwid,'(a,i3.3,a,i3.3)') 'x', procx-1, 'y', procy
      else
      	write(cwid,'(a,i3.3,a,i3.3)') 'x', nprocx-1, 'y', procy
      end if

      ! North
      if (procy > 0) then
        write(cnid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', procy-1
      else
      	write(cnid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', nprocy-1
      end if

      ! South
      if (procy < nprocy-1) then
        write(csid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', procy+1
      else
      	write(csid,'(a,i3.3,a,i3.3)') 'x', procx, 'y', 0
      end if

      !!
	  !! READ RESTARTFILES
	  !!

	  !! CURRENT PROCESSOR
	  name(5:5) = 'd'
	  name(13:20) = cmyid
	  print *, 'Reading restartfiles for ', name
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

      ! If the input fields are bigger than the output fields, you don't need the adjacent procs
      ! Assumes ih and jh behave consistently
      if (iho > ih .or. lcoarse) then

		!! EAST PROCESSOR
		name(5:5) = 'd'
		name(13:20) = ceid
		write(6,*) 'East processor: ',name
		open(unit=ifinput,file=trim(inpath)//'/'//name,form='unformatted', status='old')

		read(ifinput)  (((u0e    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((v0e    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((w0e    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((thl0e  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qt0e   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0e   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0he  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((e120e  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((dthvdze(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekme   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekhe   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((tmp0e   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((esle   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsle   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsie   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)   ((ustare (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((thlfluxe (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qtfluxe  (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dthldze(i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dqtdze (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)  (  presf (    k)                            ,k=1,k1)
		read(ifinput)  (  presh (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
		read(ifinput)  ps,thls,qts,thvs,oblav
		read(ifinput)  dtheta,dqt,timee,dt,tres
		read(ifinput)   ((oble (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((tskine(i,j ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qskine(i,j ),i=1,i2      ),j=1,j2      )

		!!!!! radiation quantities
		read(ifinput)  tnext_radiation
		read(ifinput)  (((thlprade (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
	    read(ifinput)  (((swde     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swue     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwde     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwue     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdcae   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swucae   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdcae   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwucae   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdire   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdife   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwce     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

		read(ifinput)  ((SW_up_TOAe    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_TOAe    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_TOAe    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_TOAe    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_up_ca_TOAe (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_ca_TOAe (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_ca_TOAe (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_ca_TOAe (i,j ),i=1,i2),j=1,j2)

		close(ifinput)

		if (nsv>0) then
		  name(5:5) = 's'
		  write(6,*) 'East processor s: ',name
		  open(unit=ifinput,file=trim(outpath)//'/'//name,form='unformatted')
		  read(ifinput) ((((sv0e(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
		  read(ifinput) (((svfluxe(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
		  read(ifinput) (dsv(n),n=1,nsv)
		  read(ifinput)  timee
		  close(ifinput)
		end if

		!! WEST PROCESSOR
		name(5:5) = 'd'
		name(13:20) = cwid
		write(6,*) 'West processor: ',name
		open(unit=ifinput,file=trim(inpath)//'/'//name,form='unformatted', status='old')

		read(ifinput)  (((u0w    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((v0w    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((w0w    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((thl0w  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qt0w   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0w   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0hw  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((e120w  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((dthvdzw(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekmw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekhw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((tmp0w   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((eslw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvslw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsiw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)   ((ustarw (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((thlfluxw (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qtfluxw  (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dthldzw(i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dqtdzw (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)  (  presf (    k)                            ,k=1,k1)
		read(ifinput)  (  presh (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
		read(ifinput)  ps,thls,qts,thvs,oblav
		read(ifinput)  dtheta,dqt,timee,dt,tres
		read(ifinput)   ((oblw (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((tskinw(i,j ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qskinw(i,j ),i=1,i2      ),j=1,j2      )

		!!!!! radiation quantities
		read(ifinput)  tnext_radiation
		read(ifinput)  (((thlpradw (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdw     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swuw     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdw     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwuw     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdcaw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swucaw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdcaw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwucaw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdirw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdifw   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwcw     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

		read(ifinput)  ((SW_up_TOAw    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_TOAw    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_TOAw    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_TOAw    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_up_ca_TOAw (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_ca_TOAw (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_ca_TOAw (i,j ),i=1,i2),j=1,j2)
	    read(ifinput)  ((LW_dn_ca_TOAw (i,j ),i=1,i2),j=1,j2)

		close(ifinput)

		if (nsv>0) then
		  name(5:5) = 's'
		  write(6,*) 'West processor s: ',name
		  open(unit=ifinput,file=trim(outpath)//'/'//name,form='unformatted')
		  read(ifinput) ((((sv0w(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
		  read(ifinput) (((svfluxw(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
		  read(ifinput) (dsv(n),n=1,nsv)
		  read(ifinput)  timee
		  close(ifinput)
		end if

		!! NORTH PROCESSOR
		name(5:5) = 'd'
		name(13:20) = cnid
		write(6,*) 'North processor: ',name
		open(unit=ifinput,file=trim(inpath)//'/'//name,form='unformatted', status='old')

		read(ifinput)  (((u0n    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((v0n    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((w0n    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((thl0n  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qt0n   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0n   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0hn  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((e120n  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((dthvdzn(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekmn   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekhn   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((tmp0n   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((esln   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsln   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsin   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)   ((ustarn (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((thlfluxn (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qtfluxn  (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dthldzn(i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dqtdzn (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)  (  presf (    k)                            ,k=1,k1)
		read(ifinput)  (  presh (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
		read(ifinput)  ps,thls,qts,thvs,oblav
		read(ifinput)  dtheta,dqt,timee,dt,tres
		read(ifinput)   ((obln (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((tskinn(i,j ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qskinn(i,j ),i=1,i2      ),j=1,j2      )

		read(ifinput)  tnext_radiation
		read(ifinput)  (((thlpradn (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdn     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swun     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdn     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwun     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdcan   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swucan   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdcan   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwucan   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdirn   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdifn   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwcn     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

		read(ifinput)  ((SW_up_TOAn    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_TOAn    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_TOAn    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_TOAn    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_up_ca_TOAn (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_ca_TOAn (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_ca_TOAn (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_ca_TOAn (i,j ),i=1,i2),j=1,j2)

		close(ifinput)

		if (nsv>0) then
		  name(5:5) = 's'
		  write(6,*) 'North processor n: ',name
		  open(unit=ifinput,file=trim(outpath)//'/'//name,form='unformatted')
		  read(ifinput) ((((sv0n(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
		  read(ifinput) (((svfluxn(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
		  read(ifinput) (dsv(n),n=1,nsv)
		  read(ifinput)  timee
		  close(ifinput)
		end if

		!! SOUTH PROCESSOR
		name(5:5) = 'd'
		name(13:20) = csid
		write(6,*) 'South processor: ',name
		open(unit=ifinput,file=trim(inpath)//'/'//name,form='unformatted', status='old')

		read(ifinput)  (((u0s    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((v0s    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((w0s    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((thl0s  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qt0s   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0s   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ql0hs  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((e120s  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((dthvdzs(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekms   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((ekhs   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((tmp0s   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((esls   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsls   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((qvsis   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)   ((ustars (i,j  ),i=1,i2      ),j=1,j2      )
	    read(ifinput)   ((thlfluxs (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qtfluxs  (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dthldzs(i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((dqtdzs (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)  (  presf (    k)                            ,k=1,k1)
		read(ifinput)  (  presh (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
		read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
		read(ifinput)  ps,thls,qts,thvs,oblav
		read(ifinput)  dtheta,dqt,timee,dt,tres
		read(ifinput)   ((obls (i,j  ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((tskins(i,j ),i=1,i2      ),j=1,j2      )
		read(ifinput)   ((qskins(i,j ),i=1,i2      ),j=1,j2      )

		read(ifinput)  tnext_radiation
		read(ifinput)  (((thlprads (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swds     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swus     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwds     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwus     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdcas   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swucas   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwdcas   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwucas   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdirs   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((swdifs   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
		read(ifinput)  (((lwcs     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

		read(ifinput)  ((SW_up_TOAs    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_TOAs    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_TOAs    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_TOAs    (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_up_ca_TOAs (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((SW_dn_ca_TOAs (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_up_ca_TOAs (i,j ),i=1,i2),j=1,j2)
		read(ifinput)  ((LW_dn_ca_TOAs (i,j ),i=1,i2),j=1,j2)

		close(ifinput)

		if (nsv>0) then
		  name(5:5) = 's'
		  write(6,*) 'East processor s: ',name
		  open(unit=ifinput,file=trim(outpath)//'/'//name,form='unformatted')
		  read(ifinput) ((((sv0s(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
		  read(ifinput) (((svfluxs(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
		  read(ifinput) (dsv(n),n=1,nsv)
		  read(ifinput)  timee
		  close(ifinput)
		end if

      end if ! iho > ih .or. lcoarse

      !!
      !! ADVECTION SCHEME SWITCH
      !!
      if (ladv) then
        ! fields that are unaffected by ghost cells
        ustaro = ustar
        thlfluxo = thlflux
        qtfluxo = qtflux
        dthldzo = dthldz
        dqtdzo = dqtdz
        oblo = obl
        tskino = tskin
        qskino = qskin
        SW_up_TOAo = SW_up_TOA
        SW_dn_TOAo = SW_dn_TOA
        LW_up_TOAo = LW_up_TOA
        LW_dn_TOAo = LW_dn_TOA
        SW_up_ca_TOAo = SW_up_ca_TOA
        SW_dn_ca_TOAo = SW_dn_ca_TOA
        LW_up_ca_TOAo = LW_up_ca_TOA
        LW_dn_ca_TOAo = LW_dn_ca_TOA
      
        if (iho < ih) then
          ! Your output is smaller than your input -> You just need to sample the central processor 

          ! Indexing always starts at 2, i.e. x(2,2,:) always gives the first value in the processor,
          ! regardless of the number of ghost cells. Hence, you just need to sample iho, jho values to
          ! the left and right of this common indexing: xo = xi(2-iho:i1+iho, 2-jho:j1+jho, k1)
          u0o = u0(2-iho:i1+iho, 2-jho:j1+jho,:)
          v0o = v0(2-iho:i1+iho, 2-jho:j1+jho,:)
          w0o = w0(2-iho:i1+iho, 2-jho:j1+jho,:)
          thl0o = thl0(2-iho:i1+iho, 2-jho:j1+jho,:)
          qt0o = qt0(2-iho:i1+iho, 2-jho:j1+jho,:)
          ql0o = ql0(2-iho:i1+iho, 2-jho:j1+jho,:)
          ql0ho = ql0h(2-iho:i1+iho, 2-jho:j1+jho,:)
          e120o = e120(2-iho:i1+iho, 2-jho:j1+jho,:)
          dthvdzo = dthvdz(2-iho:i1+iho, 2-jho:j1+jho,:)
          ekmo = ekm(2-iho:i1+iho, 2-jho:j1+jho,:)
          ekho = ekh(2-iho:i1+iho, 2-jho:j1+jho,:)
          tmp0o = tmp0(2-iho:i1+iho, 2-jho:j1+jho,:)
          eslo = esl(2-iho:i1+iho, 2-jho:j1+jho,:)
          qvslo = qvsl(2-iho:i1+iho, 2-jho:j1+jho,:)
          qvsio = qvsi(2-iho:i1+iho, 2-jho:j1+jho,:)
          thlprado = thlprad(2-iho:i1+iho, 2-jho:j1+jho,:)
          swdo = swd(2-iho:i1+iho, 2-jho:j1+jho,:)
          swuo = swu(2-iho:i1+iho, 2-jho:j1+jho,:)
          lwdo = lwd(2-iho:i1+iho, 2-jho:j1+jho,:)
          lwuo = lwu(2-iho:i1+iho, 2-jho:j1+jho,:)
          swdcao = swdca(2-iho:i1+iho, 2-jho:j1+jho,:)
          swucao = swuca(2-iho:i1+iho, 2-jho:j1+jho,:)
          lwdcao = lwdca(2-iho:i1+iho, 2-jho:j1+jho,:)
          lwucao = lwuca(2-iho:i1+iho, 2-jho:j1+jho,:)
          swdiro = swdir(2-iho:i1+iho, 2-jho:j1+jho,:)
          swdifo = swdif(2-iho:i1+iho, 2-jho:j1+jho,:)
          lwco = lwc(2-iho:i1+iho, 2-jho:j1+jho,:)

          ! TODO ADD A CLAUSE FOR JUST RETURNING THE INPUT IF BOTH LCOARSE AND LADV ARE FALSE
!         else if (iho == ih) then
!           u0o = u0
!           v0o = v0
!           w0o = w0
!           thl0o = thl0
!           qt0o = qt0
!           ql0o = ql0
!           ql0ho = ql0h
!           e120o = e120
!           dthvdzo = dthvdz
!           ekmo = ekm
!           ekho = ekh
!           tmp0o = tmp0
!           eslo = esl
!           qvslo = qvsl
!           qvsio = qvsi
!           thlprado = thlprad
!           swdo = swd
!           swuo = swu
!           lwdo = lwd
!           lwuo = lwu
!           swdcao = swdca
!           swucao = swuca
!           lwdcao = lwdca
!           lwucao = lwuca
!           swdiro = swdir
!           swdifo = swdif
!           lwco = lwc

        else
      	  ! Your output fields are larger than your input fields -> Need to sample adjacent processors

      	  ! First insert the central processor
      	  do i=2-ih, i1+ih
      	    do j=2-jh, j1+jh
      	      do k=1, k1
      	        u0o(i,j,k) = u0(i,j,k)
      	        v0o(i,j,k) = v0(i,j,k)
      	        w0o(i,j,k) = w0(i,j,k)
      	        thl0o(i,j,k) = thl0(i,j,k)
      	        qt0o(i,j,k) = qt0(i,j,k)
      	        ql0o(i,j,k) = ql0(i,j,k)
      	        ql0ho(i,j,k) = ql0h(i,j,k)
      	        e120o(i,j,k) = e120(i,j,k)
      	        dthvdzo(i,j,k) = dthvdz(i,j,k)
      	        ekmo(i,j,k) = ekm(i,j,k)
      	        ekho(i,j,k) = ekh(i,j,k)
      	        tmp0o(i,j,k) = tmp0(i,j,k)
      	        eslo(i,j,k) = esl(i,j,k)
      	        qvslo(i,j,k) = qvsl(i,j,k)
      	        qvsio(i,j,k) = qvsi(i,j,k)
      	        thlprado(i,j,k) = thlprad(i,j,k)
      	        swdo(i,j,k) = swd(i,j,k)
      	        swuo(i,j,k) = swu(i,j,k)
      	        lwdo(i,j,k) = lwd(i,j,k)
      	        lwuo(i,j,k) = lwu(i,j,k)
      	        swdcao(i,j,k) = swdca(i,j,k)
      	        swucao(i,j,k) = swuca(i,j,k)
      	        lwdcao(i,j,k) = lwdca(i,j,k)
      	        lwucao(i,j,k) = lwuca(i,j,k)
      	        swdiro(i,j,k) = swdir(i,j,k)
      	        swdifo(i,j,k) = swdif(i,j,k)
      	        lwco(i,j,k) = lwc(i,j,k)
      	      end do
            end do
          end do

          ! Insert ghost cells. Account for:
          ! - counter i runs from 2-iho, not 0
          ! - need first/last/last/first points in e/w/n/s-proc, minus ghost cells from the central proc

          ! Insert missing ghost cells on the east boundary
          do i=i1+ih+1, i1+iho
!             print *, 'EAST PROC: output i: ', i, 'input i: ', 1+i-i1
            do j=2, j1
              do k=1, k1
                u0o(i,j,k) = u0e(1+i-i1,j,k)
        	    v0o(i,j,k) = v0e(1+i-i1,j,k)
        	    w0o(i,j,k) = w0e(1+i-i1,j,k)
        	    thl0o(i,j,k) = thl0e(1+i-i1,j,k)
        	    qt0o(i,j,k) = qt0e(1+i-i1,j,k)
        	    ql0o(i,j,k) = ql0e(1+i-i1,j,k)
        	    ql0ho(i,j,k) = ql0he(1+i-i1,j,k)
        	    e120o(i,j,k) = e120e(1+i-i1,j,k)
        	    dthvdzo(i,j,k) = dthvdze(1+i-i1,j,k)
        	    ekmo(i,j,k) = ekme(1+i-i1,j,k)
        	    ekho(i,j,k) = ekhe(1+i-i1,j,k)
        	    tmp0o(i,j,k) = tmp0e(1+i-i1,j,k)
        	    eslo(i,j,k) = esle(1+i-i1,j,k)
        	    qvslo(i,j,k) = qvsle(1+i-i1,j,k)
        	    qvsio(i,j,k) = qvsie(1+i-i1,j,k)
        	    thlprado(i,j,k) = thlprade(1+i-i1,j,k)
        	    swdo(i,j,k) = swde(1+i-i1,j,k)
        	    swuo(i,j,k) = swue(1+i-i1,j,k)
        	    lwdo(i,j,k) = lwde(1+i-i1,j,k)
        	    lwuo(i,j,k) = lwue(1+i-i1,j,k)
        	    swdcao(i,j,k) = swdcae(1+i-i1,j,k)
        	    swucao(i,j,k) = swucae(1+i-i1,j,k)
        	    lwdcao(i,j,k) = lwdcae(1+i-i1,j,k)
        	    lwucao(i,j,k) = lwucae(1+i-i1,j,k)
        	    swdiro(i,j,k) = swdire(1+i-i1,j,k)
        	    swdifo(i,j,k) = swdife(1+i-i1,j,k)
        	    lwco(i,j,k) = lwce(1+i-i1,j,k)
              end do
		    end do
          end do

      	  ! Insert missing ghost cells on the west boundary
      	  do i=2-iho, 1-ih
!             print *, 'WEST PROC: output i: ', i, 'input i: ', i1-(2-ih)+i
      	    do j=2, j1
      	  	  do k=1, k1
                u0o(i,j,k) = u0w(i1-(2-ih)+i, j, k)
                v0o(i,j,k) = v0w(i1-(2-ih)+i, j, k)
                w0o(i,j,k) = w0w(i1-(2-ih)+i, j, k)
                thl0o(i,j,k) = thl0w(i1-(2-ih)+i, j, k)
                qt0o(i,j,k) = qt0w(i1-(2-ih)+i, j, k)
                ql0o(i,j,k) = ql0w(i1-(2-ih)+i, j, k)
                ql0ho(i,j,k) = ql0hw(i1-(2-ih)+i, j, k)
                e120o(i,j,k) = e120w(i1-(2-ih)+i, j, k)
                dthvdzo(i,j,k) = dthvdzw(i1-(2-ih)+i, j, k)
                ekmo(i,j,k) = ekmw(i1-(2-ih)+i, j, k)
                ekho(i,j,k) = ekhw(i1-(2-ih)+i, j, k)
                tmp0o(i,j,k) = tmp0w(i1-(2-ih)+i, j, k)
                eslo(i,j,k) = eslw(i1-(2-ih)+i, j, k)
                qvslo(i,j,k) = qvsle(i1-(2-ih)+i, j, k)
                qvsio(i,j,k) = qvsiw(i1-(2-ih)+i, j, k)
                thlprado(i,j,k) = thlpradw(i1-(2-ih)+i, j, k)
                swdo(i,j,k) = swdw(i1-(2-ih)+i, j, k)
                swuo(i,j,k) = swuw(i1-(2-ih)+i, j, k)
                lwdo(i,j,k) = lwdw(i1-(2-ih)+i, j, k)
                lwuo(i,j,k) = lwuw(i1-(2-ih)+i, j, k)
                swdcao(i,j,k) = swdcaw(i1-(2-ih)+i, j, k)
                swucao(i,j,k) = swucaw(i1-(2-ih)+i, j, k)
                lwdcao(i,j,k) = lwdcaw(i1-(2-ih)+i, j, k)
                lwucao(i,j,k) = lwucaw(i1-(2-ih)+i, j, k)
                swdiro(i,j,k) = swdirw(i1-(2-ih)+i, j, k)
                swdifo(i,j,k) = swdifw(i1-(2-ih)+i, j, k)
                lwco(i,j,k) = lwcw(i1-(2-ih)+i, j, k)
              end do
            end do
          end do

          ! Insert missing ghost cells on the north boundary
          do j=j1+jh+1, j1+jho
!             print *, 'NORTH PROC: output j: ', j, 'input j: ', 1+j-j1
            do i=2, i1
              do k=1, k1
        	    u0o(i,j,k) = u0n(i,1+j-j1,k)
        	    v0o(i,j,k) = v0n(i,1+j-j1,k)
        	    w0o(i,j,k) = w0n(i,1+j-j1,k)
        	    thl0o(i,j,k) = thl0n(i,1+j-j1,k)
        	    qt0o(i,j,k) = qt0n(i,1+j-j1,k)
        	    ql0o(i,j,k) = ql0n(i,1+j-j1,k)
        	    ql0ho(i,j,k) = ql0hn(i,1+j-j1,k)
        	    e120o(i,j,k) = e120n(i,1+j-j1,k)
        	    dthvdzo(i,j,k) = dthvdzn(i,1+j-j1,k)
        	    ekmo(i,j,k) = ekmn(i,1+j-j1,k)
        	    ekho(i,j,k) = ekhn(i,1+j-j1,k)
        	    tmp0o(i,j,k) = tmp0n(i,1+j-j1,k)
        	    eslo(i,j,k) = esln(i,1+j-j1,k)
        	    qvslo(i,j,k) = qvsln(i,1+j-j1,k)
        	    qvsio(i,j,k) = qvsin(i,1+j-j1,k)
        	    thlprado(i,j,k) = thlpradn(i,1+j-j1,k)
        	    swdo(i,j,k) = swdn(i,1+j-j1,k)
        	    swuo(i,j,k) = swun(i,1+j-j1,k)
        	    lwdo(i,j,k) = lwdn(i,1+j-j1,k)
        	    lwuo(i,j,k) = lwun(i,1+j-j1,k)
        	    swdcao(i,j,k) = swdcan(i,1+j-j1,k)
        	    swucao(i,j,k) = swucan(i,1+j-j1,k)
        	    lwdcao(i,j,k) = lwdcan(i,1+j-j1,k)
        	    lwucao(i,j,k) = lwucan(i,1+j-j1,k)
        	    swdiro(i,j,k) = swdirn(i,1+j-j1,k)
        	    swdifo(i,j,k) = swdifn(i,1+j-j1,k)
        	    lwco(i,j,k) = lwcn(i,1+j-j1,k)
    		  end do
		    end do
          end do

          ! Insert missing ghost cells on the south boundary
      	  do j=2-jho, 1-jh
!             print *, 'SOUTH PROC: output j: ', j, 'input i: ', j1-(2-jh)+j
            do i=2, i1
              do k=1, k1
                u0o(i,j,k) = u0s(i, j1-(2-jh)+j, k)
                v0o(i,j,k) = v0s(i, j1-(2-jh)+j, k)
                w0o(i,j,k) = w0s(i, j1-(2-jh)+j, k)
                thl0o(i,j,k) = thl0s(i, j1-(2-jh)+j, k)
                qt0o(i,j,k) = qt0s(i, j1-(2-jh)+j, k)
                ql0o(i,j,k) = ql0s(i, j1-(2-jh)+j, k)
                ql0ho(i,j,k) = ql0hs(i, j1-(2-jh)+j, k)
                e120o(i,j,k) = e120s(i, j1-(2-jh)+j, k)
                dthvdzo(i,j,k) = dthvdzs(i, j1-(2-jh)+j, k)
                ekmo(i,j,k) = ekms(i, j1-(2-jh)+j, k)
                ekho(i,j,k) = ekhs(i, j1-(2-jh)+j, k)
                tmp0o(i,j,k) = tmp0s(i, j1-(2-jh)+j, k)
                eslo(i,j,k) = esls(i, j1-(2-jh)+j, k)
                qvslo(i,j,k) = qvsls(i, j1-(2-jh)+j, k)
                qvsio(i,j,k) = qvsis(i, j1-(2-jh)+j, k)
                thlprado(i,j,k) = thlprads(i, j1-(2-jh)+j, k)
                swdo(i,j,k) = swds(i, j1-(2-jh)+j, k)
                swuo(i,j,k) = swus(i, j1-(2-jh)+j, k)
                lwdo(i,j,k) = lwds(i, j1-(2-jh)+j, k)
                lwuo(i,j,k) = lwus(i, j1-(2-jh)+j, k)
                swdcao(i,j,k) = swdcas(i, j1-(2-jh)+j, k)
                swucao(i,j,k) = swucas(i, j1-(2-jh)+j, k)
                lwdcao(i,j,k) = lwdcas(i, j1-(2-jh)+j, k)
                lwucao(i,j,k) = lwucas(i, j1-(2-jh)+j, k)
                swdiro(i,j,k) = swdirs(i, j1-(2-jh)+j, k)
                swdifo(i,j,k) = swdifs(i, j1-(2-jh)+j, k)
                lwco(i,j,k) = lwcs(i, j1-(2-jh)+j, k)
              end do
            end do
          end do

	    end if ! iho > ih

	  !!
	  !! COARSE-GRAINING SWITCH
	  !!
	  elseif (lcoarse) then
	  	
	    ! Central processor
! 	  	u0 = 1.
	  	ii = 2-ncoarse
	  	ij = 2-ncoarse
	  	do i=2,i1o
	  	  ii = ii+ncoarse
	  	  do j=2,j1o
	  	  	ij = ij+ncoarse
! 	  	  	if (procx == 0 .and. procy == 0) then
! 	  	  	  print *, 'i:',i,'j:',j,'ii:',ii,'ij:',ij
! 	  	  	  print *, sum(u0(ii:ii+ncoarse-1,ij:ij+ncoarse-1,1))*ncoarse2i
! 	  	  	end if
	  	  	do k=1,k1
	  	  	  u0o(i,j,k) = sum(u0(ii:ii+ncoarse-1,ij:ij+ncoarse-1,k))*ncoarse2i
	  	  	end do
	  	  end do
	  	  ij = 0
	  	end do

	  	! East boundary
! 	  	u0e = 1.
	  	ii = 2-ncoarse ! Start at west boundary of eastern proc
	  	ij = 2-ncoarse
	  	do i=i2o,i1o+iho
	  	  ii = ii+ncoarse
	  	  do j=2,j1o
	  	  	ij = ij+ncoarse
! 	  	  	if (procx == 0 .and. procy == 0) then
! 	  	  	  print *, 'i:',i,'j:',j,'ii:',ii,'ij:',ij
! 	  	  	  print *, sum(u0e(ii:ii+ncoarse-1,ij:ij+ncoarse-1,1))*ncoarse2i
! 	  	  	end if
	  	  	do k=1,k1
	  	  	  u0o(i,j,k) = sum(u0e(ii:ii+ncoarse-1,ij:ij+ncoarse-1,k))*ncoarse2i
	  	  	end do
	  	  end do
	  	  ij = 0
	  	end do

	  end if ! ladv / lcoarse
	  

	  !!
	  !! WRITE
	  !!
	  name(13:20)= cmyid
	  write(6,*) 'writing ',name
	  open(ifoutput,file=trim(outpath)//'/'//name,form='unformatted',status='replace')

	  write(ifoutput)  (((u0o    (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((v0o    (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((w0o    (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((thl0o  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((qt0o   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((ql0o   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((ql0ho  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((e120o  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((dthvdzo(i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((ekmo   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((ekho   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((tmp0o  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((eslo   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((qvslo  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((qvsio  (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  ((ustaro  (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  ((thlfluxo(i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  ((qtfluxo (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  ((dthldzo (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  ((dqtdzo  (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)  (  presf (    k)                            ,k=1,k1)
	  write(ifoutput)  (  presh (    k)                            ,k=1,k1)
	  write(ifoutput)  (  initial_presf (    k)                            ,k=1,k1)
	  write(ifoutput)  (  initial_presh (    k)                            ,k=1,k1)
	  write(ifoutput)  ps,thls,qts,thvs,oblav
	  write(ifoutput)  dtheta,dqt,timee,  dt,tres
	  write(ifoutput)   ((oblo (i,j  ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((tskino(i,j ),i=1,i2      ),j=1,j2      )
	  write(ifoutput)   ((qskino(i,j ),i=1,i2      ),j=1,j2      )

	!!!!! radiation quantities
	  write(ifoutput)  tnext_radiation
	  write(ifoutput)  (((thlprado (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swdo     (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swuo     (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((lwdo     (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((lwuo     (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swdcao   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swucao   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((lwdcao   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((lwucao   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swdiro   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((swdifo   (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)
	  write(ifoutput)  (((lwco     (i,j,k),i=2-iho,i1+iho),j=2-jho,j1+jho),k=1,k1)

	  write(ifoutput)  ((SW_up_TOAo    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_dn_TOAo    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_up_TOAo    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_dn_TOAo    (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_up_ca_TOAo (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((SW_dn_ca_TOAo (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_up_ca_TOAo (i,j ),i=1,i2),j=1,j2)
	  write(ifoutput)  ((LW_dn_ca_TOAo (i,j ),i=1,i2),j=1,j2)

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

  deallocate(u0e, v0e, w0e, thl0e, qt0e, ql0e, ql0he, e120e, dthvdze, ekme, ekhe, tmp0e, esle, qvsle, qvsie)
  deallocate(ustare, thlfluxe, qtfluxe, dthldze, dqtdze)
  deallocate(oble, tskine, qskine, thlprade, swde, swue, lwde, lwue, swdcae, swucae)
  deallocate(lwdcae, lwucae, swdire, swdife, lwce, SW_up_TOAe, SW_dn_TOAe, LW_up_TOAe, LW_dn_TOAe)
  deallocate(SW_up_ca_TOAe, SW_dn_ca_TOAe, LW_up_ca_TOAe, LW_dn_ca_TOAe)
  deallocate(sv0e, svfluxe)

  deallocate(u0w, v0w, w0w, thl0w, qt0w, ql0w, ql0hw, e120w, dthvdzw, ekmw, ekhw, tmp0w, eslw, qvslw, qvsiw)
  deallocate(ustarw, thlfluxw, qtfluxw, dthldzw, dqtdzw)
  deallocate(oblw, tskinw, qskinw, thlpradw, swdw, swuw, lwdw, lwuw, swdcaw, swucaw)
  deallocate(lwdcaw, lwucaw, swdirw, swdifw, lwcw, SW_up_TOAw, SW_dn_TOAw, LW_up_TOAw, LW_dn_TOAw)
  deallocate(SW_up_ca_TOAw, SW_dn_ca_TOAw, LW_up_ca_TOAw, LW_dn_ca_TOAw)
  deallocate(sv0w, svfluxw)

  deallocate(u0n, v0n, w0n, thl0n, qt0n, ql0n, ql0hn, e120n, dthvdzn, ekmn, ekhn, tmp0n, esln, qvsln, qvsin)
  deallocate(ustarn, thlfluxn, qtfluxn, dthldzn, dqtdzn)
  deallocate(obln, tskinn, qskinn, thlpradn, swdn, swun, lwdn, lwun, swdcan, swucan)
  deallocate(lwdcan, lwucan, swdirn, swdifn, lwcn, SW_up_TOAn, SW_dn_TOAn, LW_up_TOAn, LW_dn_TOAn)
  deallocate(SW_up_ca_TOAn, SW_dn_ca_TOAn, LW_up_ca_TOAn, LW_dn_ca_TOAn)
  deallocate(sv0n, svfluxn)

  deallocate(u0s, v0s, w0s, thl0s, qt0s, ql0s, ql0hs, e120s, dthvdzs, ekms, ekhs, tmp0s, esls, qvsls, qvsis)
  deallocate(ustars, thlfluxs, qtfluxs, dthldzs, dqtdzs)
  deallocate(obls, tskins, qskins, thlprads, swds, swus, lwds, lwus, swdcas, swucas)
  deallocate(lwdcas, lwucas, swdirs, swdifs, lwcs, SW_up_TOAs, SW_dn_TOAs, LW_up_TOAs, LW_dn_TOAs)
  deallocate(SW_up_ca_TOAs, SW_dn_ca_TOAs, LW_up_ca_TOAs, LW_dn_ca_TOAs)
  deallocate(sv0s, svfluxs)

  deallocate(u0o, v0o, w0o, thl0o, qt0o, ql0o, ql0ho, e120o, dthvdzo, ekmo, ekho, tmp0o, eslo, qvslo, qvsio)
  deallocate(ustaro, thlfluxo, qtfluxo, dthldzo, dqtdzo)
  deallocate(oblo, tskino, qskino, thlprado, swdo, swuo, lwdo, lwuo, swdcao, swucao)
  deallocate(lwdcao, lwucao, swdiro, swdifo, lwco, SW_up_TOAo, SW_dn_TOAo, LW_up_TOAo, LW_dn_TOAo)
  deallocate(SW_up_ca_TOAo, SW_dn_ca_TOAo, LW_up_ca_TOAo, LW_dn_ca_TOAo)
  deallocate(sv0o, svfluxo)

end program reproject
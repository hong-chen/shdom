MODULE molec_absorption_interface

! This module should contain the single column forward, tangent linear,
! and adjoint routines for the calculation of monochromatic or band average
! layer molecular optical depth for a single satellite channel.
!   The layers input to this module go from the top of the atmosphere down.
! This is a dummy module to show the interface.  It returns zero 
! optical depth and gradients.


IMPLICIT NONE

PRIVATE  ! All variables are private to this module
PUBLIC :: initialize_molec_abs, dealloc_molec_abs_arrays, &
          forward_molec_absorption, tanlin_molec_absorption, &
          adjoint_molec_absorption


CONTAINS


SUBROUTINE initialize_molec_abs ()
 ! Initializes the molecular absorption calculation for a
 ! particular set of satellite channels.
 !
  IMPLICIT NONE

END SUBROUTINE initialize_molec_abs




SUBROUTINE forward_molec_absorption (channel_index, nlay, &
                                     preslev, pres, temp, vapmixr, o3mixr, &
                                     molectau)
 ! Calculates the layer optical depths for molecular absorption for one 
 ! satellite channel for the input atmospheric profile.  The profile layer 
 ! arrays go from TOA down.
 !
 ! Inputs:
 !  channel_index     satellite channel index
 !  nlay              number of layers
 !  preslev           layer interface pressures (from top of top layer downwards) (mb)
 !  pres              layer midpoint pressures (mb)
 !  temp              layer midpoint temperatures (K)
 !  vapmixr           layer midpoint water vapor mixing ratios (g/kg)
 !  o3mixr            layer midpoint ozone mixing ratios (ppmv)
 !
 ! Outputs:
 !  molectau          layer optical depth from molecular absorption
 !                      (equivalent vertical optical path)
 !
  INTEGER, INTENT(IN) :: channel_index, nlay
  REAL,    INTENT(IN) :: preslev(0:nlay)
  REAL,    INTENT(IN) :: pres(nlay), temp(nlay), vapmixr(nlay), o3mixr(nlay)
  REAL,    INTENT(OUT) :: molectau(nlay)

  molectau(:) = 0.0
END SUBROUTINE forward_molec_absorption




SUBROUTINE tanlin_molec_absorption (channel_index, nlay, &
                                     preslev, pres, temp, vapmixr, o3mixr, &
                                     g_temp, g_vapmixr, g_o3mixr, &
                                     g_molectau)
 ! Calculates the tangent linear of the layer optical depths for 
 ! molecular absorption for one satellite channel for the input 
 ! atmospheric profile.  The profile layer 
 ! arrays go from TOA down.
 !
 ! Inputs:
 !  channel_index     satellite channel index
 !  nlay              number of layers
 !  preslev           layer interface pressures (from top of top layer downwards) (mb)
 !  pres              layer midpoint pressures (mb)
 !  temp              layer midpoint temperatures (K)
 !  vapmixr           layer midpoint water vapor mixing ratios (g/kg)
 !  o3mixr            layer midpoint ozone mixing ratios (ppmv)
 !  g_temp            perturbations in layer temperatures (K)
 !  g_vapmixr         perturbations in layer water vapor mixing ratios (g/kg)
 !  g_o3mixr          perturbations in layer ozone mixing ratios (ppmv)
 !
 ! Outputs:
 !  g_molectau        output change in layer optical depth from molecular 
 !                      absorption
 !
  INTEGER, INTENT(IN) :: channel_index, nlay
  REAL,    INTENT(IN) :: preslev(0:nlay)
  REAL,    INTENT(IN) :: pres(nlay), temp(nlay), vapmixr(nlay), o3mixr(nlay)
  REAL,    INTENT(IN) :: g_temp(nlay), g_vapmixr(nlay), g_o3mixr(nlay)
  REAL,    INTENT(OUT) :: g_molectau(nlay)

  g_molectau(:) = 0.0
END SUBROUTINE tanlin_molec_absorption








SUBROUTINE adjoint_molec_absorption (channel_index, nlay, &
                                      preslev, pres, temp, vapmixr, o3mixr, &
                                      admolectau, &
                                      adtemp, advapmixr, ado3mixr)
 ! Calculates the adjoint atmospheric profile (temperature, water vapor, 
 ! and ozone) from the adjoint optical depths for molecular absorption for 
 ! one satellite channel for the input atmospheric profile.  The profile 
 ! layer arrays go from TOA down.  Must be called after a call to 
 ! forward_molec_absorption with all the same parameters.  The three 
 ! output arrays are added to.
 !
 ! Inputs:
 !  channel_index     satellite channel index
 !  nlay              number of layers
 !  preslev           layer interface pressures (from top of top layer downwards) (mb)
 !  pres              layer midpoint pressures (mb)
 !  temp              layer midpoint temperatures (K)
 !  vapmixr           layer midpoint water vapor mixing ratios (g/kg)
 !  o3mixr            layer midpoint ozone mixing ratios (ppmv)
 !  admolectau        adjoint molecular absorption layer optical depths
 !                      (equivalent vertical optical path)
 !
 ! Outputs:
 !  adtemp           adjoint temperature profile
 !  advapmixr        adjoint water vapor mixing ratio profile
 !  ado3mixr         adjoint ozone mixing ratio profile
 !
  INTEGER, INTENT(IN) :: channel_index, nlay
  REAL,    INTENT(IN) :: preslev(0:nlay)
  REAL,    INTENT(IN) :: pres(nlay), temp(nlay), vapmixr(nlay), o3mixr(nlay)
  REAL,    INTENT(IN) :: admolectau(nlay)
  REAL,    INTENT(INOUT) :: adtemp(nlay), advapmixr(nlay), ado3mixr(nlay)

END SUBROUTINE adjoint_molec_absorption



SUBROUTINE dealloc_molec_abs_arrays
 ! Deallocates arrays allocated in forward_molec_absorption

END SUBROUTINE dealloc_molec_abs_arrays


END MODULE molec_absorption_interface


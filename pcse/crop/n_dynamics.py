#!/usr/bin/env python

from .. import exceptions as exc
from ..traitlets import Float, Int, Instance
from ..decorators import prepare_rates, prepare_states
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from ..util import AfgenTrait
from .nutrients import N_Demand_Uptake


class N_Crop_Dynamics(SimulationObject):
    """Implementation of overall N crop dynamics.

    NPK_Crop_Dynamics implements the overall logic of N book-keeping within the
    crop.

    **Simulation parameters**
    
    =============  ================================================= =======================
     Name           Description                                        Unit
    =============  ================================================= =======================
    NMAXLV_TB      Maximum N concentration in leaves as               kg N kg-1 dry biomass
                   function of dvs
    NMAXRT_FR      Maximum N concentration in roots as fraction       -
                   of maximum N concentration in leaves
    NMAXST_FR      Maximum N concentration in stems as fraction       -
                   of maximum N concentration in leaves
    NRESIDLV       Residual N fraction in leaves                      kg N kg-1 dry biomass
    NRESIDRT       Residual N fraction in roots                       kg N kg-1 dry biomass
    NRESIDST       Residual N fraction in stems                       kg N kg-1 dry biomass
    =============  ================================================= =======================

    **State variables**

    ==========  ================================================== ============
     Name        Description                                          Unit
    ==========  ================================================== ============
    NamountLV     Actual N amount in living leaves                  |kg N ha-1|
    NamountST     Actual N amount in living stems                   |kg N ha-1|
    NamountSO     Actual N amount in living storage organs          |kg N ha-1|
    NamountRT     Actual N amount in living roots                   |kg N ha-1| 
    Nuptake_T    total absorbed N amount                            |kg N ha-1|
    Nfix_T       total biological fixated N amount                  |kg N ha-1|
    ==========  ================================================== ============

    **Rate variables**

    ===========  =================================================  ================
     Name         Description                                           Unit
    ===========  =================================================  ================
    RNamountLV     Weight increase (N) in leaves                    |kg N ha-1 d-1|    
    RNamountST     Weight increase (N) in stems                     |kg N ha-1 d-1|
    RNamountRT     Weight increase (N) in roots                     |kg N ha-1 d-1|    
    RNamountSO     Weight increase (N) in storage organs            |kg N ha-1 d-1|
    RNdeathLV      Rate of N loss in leaves                         |kg N ha-1 d-1|
    RNdeathST      Rate of N loss in roots                          |kg N ha-1 d-1|
    RNdeathRT      Rate of N loss in stems                          |kg N ha-1 d-1|
    RNloss         N loss due to senescence                         |kg N ha-1 d-1|
    ===========  =================================================  ================
    
    **Signals send or handled**
    
    None
    
    **External dependencies**
    
    =======  =================================== ====================  ==============
     Name     Description                         Provided by            Unit
    =======  =================================== ====================  ==============
    DVS      Crop development stage              DVS_Phenology           -
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    WRT      Dry weight of living roots          WOFOST_Root_Dynamics  |kg ha-1|
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    DRLV     Death rate of leaves                WOFOST_Leaf_Dynamics  |kg ha-1 d-1|
    DRRT     Death rate of roots                 WOFOST_Root_Dynamics  |kg ha-1 d-1|
    DRST     Death rate of stems                 WOFOST_Stem_Dynamics  |kg ha-1 d-1|
    =======  =================================== ====================  ==============
    """

    demand_uptake = Instance(SimulationObject)

    NamountLVI = Float(-99.)  # initial soil N amount in leaves
    NamountSTI = Float(-99.)  # initial soil N amount in stems
    NamountRTI = Float(-99.)  # initial soil N amount in roots
    NamountSOI = Float(-99.)  # initial soil N amount in storage organs
    
    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()
        NMAXST_FR = Float(-99.)
        NMAXRT_FR = Float(-99.)
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NRESIDRT = Float(-99.)  # residual N fraction in roots [kg N kg-1 dry biomass]

    class StateVariables(StatesTemplate):
        NamountLV = Float(-99.) # N amount in leaves [kg N ha-1]
        NamountST = Float(-99.) # N amount in stems [kg N ]      
        NamountSO = Float(-99.) # N amount in storage organs [kg N ]
        NamountRT = Float(-99.) # N amount in roots [kg N ]        
        NuptakeTotal = Float(-99.) # total absorbed N amount [kg N ]
        NfixTotal = Float(-99.) # total biological fixated N amount [kg N ]
    
        NlossesTotal = Float(-99.)

    class RateVariables(RatesTemplate):
        RNamountLV = Float(-99.)  # Net rates of N in different plant organs 
        RNamountST = Float(-99.)
        RNamountRT = Float(-99.)
        RNdeathLV = Float(-99.)  # N loss rate leaves [kg ha-1 d-1]
        RNdeathST = Float(-99.)  # N loss rate stems  [kg ha-1 d-1]
        RNdeathRT = Float(-99.)  # N loss rate roots  [kg ha-1 d-1]        
        RNloss = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameters as key/value pairs
        """  
        
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # Initialize components of the npk_crop_dynamics
        self.demand_uptake = N_Demand_Uptake(day, kiosk, parvalues)

        p = self.params
        k = kiosk
        self.NamountLVI = k.WLV * p.NMAXLV_TB(k.DVS)
        self.NamountSTI = k.WST * p.NMAXLV_TB(k.DVS) * p.NMAXST_FR
        self.NamountRTI = k.WRT * p.NMAXLV_TB(k.DVS) * p.NMAXRT_FR
        self.NamountSOI = 0.

        s = dict(
            NamountLV=self.NamountLVI,
            NamountST=self.NamountSTI,
            NamountRT=self.NamountRTI,
            NamountSO=self.NamountSOI,
            NuptakeTotal=0,
            NfixTotal=0.,
            NlossesTotal=0
        )
        
        self.states = self.StateVariables(kiosk, publish=["NamountLV", "NamountST", "NamountRT", "NamountSO"], **s)

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        p = self.params
        k = self.kiosk
        
        self.demand_uptake.calc_rates(day, drv)

        # Compute loss of NPK due to death of plant material
        r.RNdeathLV = p.NRESIDLV * k.DRLV
        r.RNdeathST = p.NRESIDST * k.DRST
        r.RNdeathRT = p.NRESIDRT * k.DRRT

        # N rates in leaves, stems, root and storage organs computed as
        # uptake - translocation - death.
        # except for storage organs which only take up as a result of translocation.
        r.RNamountLV = k.RNuptakeLV - k.RNtranslocationLV - r.RNdeathLV
        r.RNamountST = k.RNuptakeST - k.RNtranslocationST - r.RNdeathST
        r.RNamountRT = k.RNuptakeRT - k.RNtranslocationRT - r.RNdeathRT
        r.RNamountSO = k.RNuptakeSO + k.RNtranslocation
        r.RNloss = r.RNdeathLV + r.RNdeathST + r.RNdeathRT

        self._check_N_balance(day)
        
    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        s = self.states
        k = self.kiosk

        # N amount in leaves, stems, root and storage organs
        s.NamountLV += r.RNamountLV * delt
        s.NamountST += r.RNamountST * delt
        s.NamountRT += r.RNamountRT * delt
        s.NamountSO += r.RNamountSO * delt
                
        self.demand_uptake.integrate(day, delt)

        # total NPK uptake from soil
        s.NuptakeTotal += k.RNuptake * delt
        s.NfixTotal += k.RNfixation * delt
        s.NlossesTotal += r.RNloss * delt

    def _check_N_balance(self, day):
        s = self.states
        checksum = abs(s.NuptakeTotal + s.NfixTotal +
                       (self.NamountLVI + self.NamountSTI + self.NamountRTI + self.NamountSOI) -
                       (s.NamountLV + s.NamountST + s.NamountRT + s.NamountSO + s.NlossesTotal))

        if abs(checksum) >= 1.0:
            msg = "N flows not balanced on day %s\n" % day
            msg += "Checksum: %f, Nuptake_T: %f, Nfix_T: %f\n" % (checksum, s.NuptakeTotal, s.NfixTotal)
            msg += "NamountLVI: %f, NamountSTI: %f, NamountRTI: %f, NamountSOI: %f\n"  % \
                   (self.NamountLVI, self.NamountSTI, self.NamountRTI, self.NamountSOI)
            msg += "NamountLV: %f, NamountST: %f, NamountRT: %f, NamountSO: %f\n" % \
                   (s.NamountLV, s.NamountST, s.NamountRT, s.NamountSO)
            msg += "NLOSST: %f\n" % (s.NlossesTotal)
            raise exc.NutrientBalanceError(msg)
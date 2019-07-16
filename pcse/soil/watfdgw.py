from math import sqrt
import numpy as np
from dotmap import DotMap


from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, HasTraits, List
from ..decorators import prepare_rates, prepare_states
from ..util import limit, Afgen, merge_dict, doy
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import exceptions as exc
from .. import signals

from .soil_profile import SoilProfile


class WaterBalanceLayered(SimulationObject):

    _default_RD = 10.  # default rooting depth at 10 cm
    _RDold = None
    _RINold = Float(-99.)
    _DSLR = Int(-99)
    _RDM = None
    _RIRR = Float(-99)
    _RAIN = Float(-99)
    _WCI = Float(-99)

    # Max number of flow iterations
    MaxFlowIter = 50
    TinyFlow = 0.001

    # Maximum upward flow is 50% of amount needed to reach equilibrium between layers
    # see documentation Kees Rappoldt - page 80
    UpwardFlowLimit = 0.50

    crop_start = Bool(False)

    soil_profile = None

    class Parameters(ParamTemplate):
        IFUNRN = Int(-99)
        NOTINF = Float(-99.)
        SSI = Float(-99.)
        SSMAX = Float(-99.)
        SMLIM = Float(-99.)
        WAV = Float(-99.)

    class StateVariables(StatesTemplate):
        WTRAT = Float(-99)
        EVST = Float(-99)
        EVWT = Float(-99)
        TSR = Float(-99)
        RAINT = Float(-99)
        WDRT = Float(-99)
        TOTINF = Float(-99)
        TOTIRR = Float(-99)
        # PERCT = Float(-99)
        # LOSST = Float(-99)
        CRT = Float(-99)
        SM = Instance(np.ndarray)
        SM_MEAN = Float(-99.)
        WC = Instance(np.ndarray)
        W = Float(-99)
        WLOW = Float(-99)
        WWLOW = Float(-99)
        WBOT = Float(-99)
        WAVUPP = Float(-99)
        WAVLOW = Float(-99)
        WAVBOT = Float(-99)
        SS = Float(-99)
        BOTTOMFLOWT = Float(-99.)


    class RateVariables(RatesTemplate):
        RIN = Float(-99)
        WTRALY = Instance(np.ndarray)
        WTRA = Float(-99)
        EVS = Float(-99)
        EVW = Float(-99)
        RIRR = Float(-99)
        DWC = Instance(np.ndarray)
        # PERC = Float(-99)
        # LOSS = Float(-99)
        DRAINT = Float(-99)
        DSS = Float(-99)
        DTSR = Float(-99)
        BOTTOMFLOW = Float(-99.)

    def initialize(self, day, kiosk, parvalues):

        self.soil_profile = SoilProfile(parvalues)
        parvalues._soildata["soil_profile"] = self.soil_profile

        # Maximum rootable depth, this has to change later because RDMCR
        # is often not known at this point.
        RDMsoil = self.soil_profile.get_max_rootable_depth()
        self._RDM = min(parvalues["RDMCR"], RDMsoil)
        self.soil_profile.validate_max_rooting_depth(self._RDM)

        self.params = self.Parameters(parvalues)
        p = self.params

        self.soil_profile.determine_rooting_status(self._default_RD, self._RDM)
        self.soil_profile.compute_layer_weights(self._default_RD, self._RDM)

        if self.soil_profile.GroundWater:
            raise NotImplementedError("Groundwater influence not yet implemented.")
        else:
            # AVMAX -  maximum available content of layer(s)
            # This is calculated first to achieve an even distribution of water in the rooted top
            # if WAV is small. Note the separate limit for initial SM in the rooted zone.
            TOPLIM = 0.0
            LOWLIM = 0.0
            AVMAX = []
            # for il in range(self.ILM+1):
            for il, layer in enumerate(self.soil_profile):
                if layer.rooting_status in ["rooted", "partially rooted"]:
                    # Check whether SMLIM is within boundaries
                    SML = limit(layer.SMW, layer.SM0, p.SMLIM)
                    # the IAIRDU mess is disabled here.
                    # if (IAIRDU.EQ.1) SML = SM0(il)
                    AVMAX.append((SML - layer.SMW) * layer.Thickness)   # available in cm
                    # also if partly rooted, the total layer capacity counts in TOPLIM
                    # this means the water content of layer ILR is set as if it would be
                    # completely rooted. This water will become available after a little
                    # root growth and through numerical mixing each time step.
                    TOPLIM += AVMAX[il]
                elif layer.rooting_status == "potentially rooted":
                    # below the rooted zone the maximum is saturation (see code for WLOW in one-layer model)
                    # again the full layer capacity adds to LOWLIM.
                    SML = layer.SM0
                    AVMAX.append((SML - layer.SMW) * layer.Thickness)   # available in cm
                    LOWLIM += AVMAX[il]
                else:  # Below the potentially rooted zone
                    break

        if p.WAV <= 0.0:
            # no available water
            TOPRED = 0.0
            LOWRED = 0.0
        elif p.WAV <= TOPLIM:
            # available water fits in layer(s) 1..ILR, these layers are rooted or almost rooted
            # reduce amounts with ratio WAV / TOPLIM
            TOPRED = p.WAV / TOPLIM
            LOWRED = 0.0
        elif p.WAV < TOPLIM + LOWLIM:
            # available water fits in potentially rooted layer
            # rooted zone is filled at capacity ; the rest reduced
            TOPRED = 1.0
            LOWRED = (p.WAV-TOPLIM) / LOWLIM
        else:
            # water does not fit ; all layers "full"
            TOPRED = 1.0
            LOWRED = 1.0

        W = 0.0    ; WAVUPP = 0.0
        WLOW = 0.0 ; WAVLOW = 0.0
        SM = np.zeros(len(self.soil_profile))
        WC = np.zeros_like(SM)
        for il, layer in enumerate(self.soil_profile):
            if layer.rooting_status in ["rooted", "partially rooted"]:
                # Part of the water assigned to ILR may not actually be in the rooted zone, but it will
                # be available shortly through root growth (and through numerical mixing).
                SM[il] = layer.SMW + AVMAX[il] * TOPRED / layer.Thickness
                W      += SM[il] * layer.Thickness * layer.Wtop
                WLOW   += SM[il] * layer.Thickness * layer.Wpot
                # available water
                WAVUPP += (SM[il] - layer.SMW) * layer.Thickness * layer.Wtop
                WAVLOW += (SM[il] - layer.SMW) * layer.Thickness * layer.Wpot
            elif layer.rooting_status == "potentially rooted":
                SM[il] = layer.SMW + AVMAX[il] * LOWRED / layer.Thickness
                WLOW += SM[il] * layer.Thickness * layer.Wpot
                # available water
                WAVLOW += (SM[il] - layer.SMW) * layer.Thickness * layer.Wpot
            else:
                # below the maximum rooting depth, set SM content to wiltint point
                SM[il] = layer.SMW
            WC[il] = SM[il] * layer.Thickness

        # set groundwater depth far away for clarity ; this prevents also
        # the root routine to stop root growth when they reach the groundwater
        ZT = 999.0

        # Initial values for rootzone and subsoil water
        self._WCI = WC.sum()

        # soil evaporation, days since last rain
        layer1_half_wet = self.soil_profile[0].SMW + 0.5 * \
                                   (self.soil_profile[0].SMFCF - self.soil_profile[0].SMW)
        self._DSLR = 5 if SM[0] <= layer1_half_wet else 1

        # all summation variables of the water balance are set at zero.
        states = {
            "WTRAT": 0., "EVST": 0., "EVWT": 0., "TSR": 0.,
            "RAINT": 0., "WDRT": 0., "TOTINF": 0., "TOTIRR": 0.,
            "PERCT": 0., "LOSST": 0., "SS":0., "BOTTOMFLOWT": 0.,
            "CRT": 0., "RAINT": 0., "WLOW": WLOW, "W": W, "WC": WC, "SM":SM,
            "SS": p.SSI, "WWLOW": W+WLOW, "WBOT":0., "SM_MEAN": W/self._default_RD,
            "WAVUPP":WAVUPP, "WAVLOW": WAVLOW, "WAVBOT":0.
        }
        self.states = self.StateVariables(kiosk, publish=["WC", "SM"], **states)

        self._RINold = 0.
        self._RIRR = 0.
        self._RDold = self._default_RD

        # rate variables
        self.rates = self.RateVariables(kiosk)

        # Connect to CROP_START/CROP_FINISH signals for water balance to
        # search for crop transpiration values
        self._connect_signal(self._on_CROP_START, signals.crop_start)
        self._connect_signal(self._on_CROP_FINISH, signals.crop_finish)
        # signal for irrigation
        self._connect_signal(self._on_IRRIGATE, signals.irrigate)


    @prepare_rates
    def calc_rates(self, day, drv):
        p = self.params
        s = self.states
        k = self.kiosk
        r = self.rates

        delt = 1.0

        # Rate of irrigation (RIRR)
        r.RIRR = self._RIRR
        self._RIRR = 0.

        self._RAIN = drv.RAIN

        # Transpiration and maximum soil and surface water evaporation rates
        # are calculated by the crop evapotranspiration module.
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if "TRALY" not in self.kiosk:
            WTRALY = np.zeros_like(s.SM)
            r.WTRA = 0.
            EVWMX = drv.E0
            EVSMX = drv.ES0
        else:
            WTRALY = k.TRALY
            r.WTRA = k.TRA
            EVWMX = k.EVWMX
            EVSMX = k.EVSMX

        # Actual evaporation rates
        r.EVW = 0.
        r.EVS = 0.
        if s.SS > 1.:
            # If surface storage > 1cm then evaporate from water layer on
            # soil surface
            r.EVW = EVWMX
        else:
            # else assume evaporation from soil surface
            if self._RINold >= 1:
                # If infiltration >= 1cm on previous day assume maximum soil
                # evaporation
                r.EVS = EVSMX
                self._DSLR = 1
            else:
                # Else soil evaporation is a function days-since-last-rain (DSLR)
                EVSMXT = EVSMX * (sqrt(self._DSLR + 1) - sqrt(self._DSLR))
                r.EVS = min(EVSMX, EVSMXT + self._RINold)
                self._DSLR += 1

        # conductivities and Matric Flux Potentials for all layers
        pF = np.zeros_like(s.SM)
        conductivity = np.zeros_like(s.SM)
        matricfluxpot = np.zeros_like(s.SM)
        for i, layer in enumerate(self.soil_profile):
            pF[i] = layer.PFfromSM(s.SM[i])
            conductivity[i] = 10**layer.CONDfromPF(pF[i])
            matricfluxpot[i] = layer.MFPfromPF(pF[i])
            if self.soil_profile.GroundWater:
                raise NotImplementedError("Groundwater influence not yet implemented.")

        # Potentially infiltrating rainfall
        if p.IFUNRN == 0:
            RINPRE = (1. - p.NOTINF) * drv.RAIN
        else:
            # infiltration is function of storm size (NINFTB)
            RINPRE = (1. - p.NOTINF * self.NINFTB(drv.RAIN)) * drv.RAIN


        # Second stage preliminary infiltration rate (RINPRE)
        # including surface storage and irrigation
        RINPRE = RINPRE + r.RIRR + s.SS
        if s.SS > 0.1:
            # with surface storage, infiltration limited by SOPE
            AVAIL = RINPRE + r.RIRR - r.EVW
            RINPRE = min(self.soil_profile.SurfaceConductivity, AVAIL)

        # maximum flow at Top Boundary of each layer
        # ------------------------------------------
        # DOWNWARD flows are calculated in two ways,
        # (1) a "dry flow" from the matric flux potentials
        # (2) a "wet flow" under the current layer conductivities and downward gravity.
        # Clearly, only the dry flow may be negative (=upward). The dry flow accounts for the large
        # gradient in potential under dry conditions (but neglects gravity). The wet flow takes into
        # account gravity only and will dominate under wet conditions. The maximum of the dry and wet
        # flow is taken as the downward flow, which is then further limited in order the prevent
        # (a) oversaturation and (b) water content to decrease below field capacity.
        #
        # UPWARD flow is just the dry flow when it is negative. In this case the flow is limited
        # to a certain fraction of what is required to get the layers at equal potential, taking
        # into account, however, the contribution of an upward flow from further down. Hence, in
        # case of upward flow from the groundwater, this upward flow in propagated upward if the
        # suction gradient is sufficiently large.

        FlowMX = np.zeros(len(s.SM) + 1)
        # first get flow through lower boundary of bottom layer
        if self.soil_profile.GroundWater:
            raise NotImplementedError("Groundwater influence not yet implemented.")
        #    the old capillairy rise routine is used to estimate flow to/from the groundwater
        #    note that this routine returns a positive value for capillairy rise and a negative
        #    value for downward flow, which is the reverse from the convention in WATFDGW.

        # is = SubSoilType
        # if (ZT >= LBSL(NSL)) then
        #     # groundwater below the layered system ; call the old capillairty rise routine
        #     # the layer PF is allocated at 1/3 * TSL above the lower boundary ; this leeds
        #     # to a reasonable result for groundwater approaching the bottom layer
        #     call SUBSOL (PF(NSL), ZT-LBSL(NSL)+TSL(NSL)/3.0, SubFlow, Soil(is)%CONDfromPF, Soil(is)%ilCOND)
        #     # write (*,*) 'call SUBSOL ', PF(NSL), ZT-LBSL(NSL)+TSL(NSL)/3.0, SubFlow
        #     if (SubFlow >= 0.0) then
        #         # capillairy rise is limited by the amount required to reach equilibrium:
        #         # step 1. calculate equilibrium ZT for all air between ZT and top of layer
        #         EqAir   = WSUB0 - WSUB + (WC0(NSL)-WC(NSL))
        #         # step 2. the groundwater level belonging to this amount of air in equilibrium
        #         ZTeq1   = (LBSL(NSL)-TSL(NSL)) + AFGEN(Soil(is)%HeightFromAir, EquilTableLEN, EqAir)
        #         # step 3. this level should normally lie below the current level (otherwise there should
        #         # not be capillairy rise). In rare cases however, due to the use of a mid-layer height
        #         # in subroutine SUBSOL, a deviation could occur
        #         ZTeq2   = MAX(ZT, ZTeq1)
        #         # step 4. calculate for this ZTeq2 the equilibrium amount of water in the layer
        #         WCequil = AFGEN(Soil(is)%WaterFromHeight, EquilTableLEN, ZTeq2-LBSL(NSL)+TSL(NSL)) - &
        #                   AFGEN(Soil(is)%WaterFromHeight, EquilTableLEN, ZTeq2-LBSL(NSL))
        #         # step5. use this equilibrium amount to limit the upward flow
        #         FlowMX(NSL+1) = -1.0 * MIN (SubFlow, MAX(WCequil-WC(NSL),0.0)/DELT)
        #     else:
        #         # downward flow ; air-filled pore space of subsoil limits downward flow
        #         AirSub = (ZT-LBSL(NSL))*SubSM0 - AFGEN(Soil(is)%WaterFromHeight, EquilTableLEN, ZT-LBSL(NSL))
        #         FlowMX(NSL+1) = MIN (ABS(SubFlow), MAX(AirSub,0.0)/DELT)
        #         # write (*,*) 'Limiting downward flow: AirSub, FlowMX(NSL+1) = ', AirSub, FlowMX(NSL+1)
        # else:
        #     # groundwater is in the layered system ; no further downward flow
        #     FlowMX(NSL+1) = 0.0
        else:
            # Bottom layer conductivity limits the flow. Below field capacity there is no
            # downward flow, so downward flow through lower boundary can be guessed as
            FlowMX[-1] = max(self.soil_profile[-1].CondFC, conductivity[-1])

        # drainage
        DMAX = 0.0

        LIMDRY = np.zeros_like(s.SM)
        LIMWET = np.zeros_like(s.SM)
        # LIMWET = self._compute_wet_flow_limit(conductivity)
        TSL = [l.Thickness for l in self.soil_profile]
        for il in reversed(range(len(s.SM))):
            # if this layers contains maximum rooting depth and if rice, downward water loss is limited
            # THIS IS HACK FOR RICE!!! -> REMOVE IT
            # if (IAIRDU==1 .and. il==ILM) FlowMX(il+1) = 0.05 * CondK0(il)

            # limiting DOWNWARD flow rate
            # == wet conditions: the soil conductivity is larger
            #    the soil conductivity is the flow rate for gravity only
            #    this limit is DOWNWARD only
            # == dry conditions: the MFP gradient
            #    the MFP gradient is larger for dry conditions
            #    allows SOME upward flow
            if il == 0:
                LIMWET[il] = self.soil_profile.SurfaceConductivity
                LIMDRY[il] = 0.0
            else:
                if self.soil_profile[il-1] == self.soil_profile[il]:
                    # Layers il-1 and il have same properties: flow rates are estimated from
                    # the gradient in Matric Flux Potential
                    LIMDRY[il] = 2.0 * (matricfluxpot[il-1]-matricfluxpot[il]) / (TSL[il-1]+TSL[il])
                    if LIMDRY[il] < 0.0:
                        # upward flow rate ; amount required for equal water content is required below
                        MeanSM = (s.WC[il-1] + s.WC[il]) / (TSL[il-1]+TSL[il])
                        EqualPotAmount = s.WC[il-1] - TSL[il-1] * MeanSM  # should be negative like the flow
                else:
                    # iterative search to PF at layer boundary (by bisection)
                    il1  = il-1               ; il2 = il
                    PF1  = pF[il1]            ; PF2 = pF[il2]
                    MFP1 = matricfluxpot[il1] ; MFP2 = matricfluxpot[il2]
                    for z in range(self.MaxFlowIter):  # Loop counter not used here
                        PFx = (PF1 + PF2) / 2.0
                        Flow1 = 2.0 * (+ MFP1 - self.soil_profile[il1].MFPfromPF(PFx)) / TSL[il1]
                        Flow2 = 2.0 * (- MFP2 + self.soil_profile[il2].MFPfromPF(PFx)) / TSL[il2]
                        if abs(Flow1-Flow2) < self.TinyFlow:
                            # sufficient accuracy
                            break
                        elif abs(Flow1) > abs(Flow2):
                            # flow in layer 1 is larger ; PFx must shift in the direction of PF1
                            PF2 = PFx
                        elif abs(Flow1) < abs(Flow2):
                            # flow in layer 2 is larger ; PFx must shift in the direction of PF2
                            PF1 = PFx
                    else:  # No break
                        msg = 'WATFDGW: LIMDRY flow iteration failed. Are your soil moisture and ' + \
                              'conductivity curves decreasing with increase pF?'
                        raise exc.PCSEError(msg)
                    LIMDRY[il] = (Flow1 + Flow2) / 2.0

                    if LIMDRY[il] < 0.0:
                        # upward flow rate ; amount required for equal potential is required below
                        Eq1 = -s.WC[il2]; Eq2 = 0.0
                        for z in range(self.MaxFlowIter):
                            EqualPotAmount = (Eq1 + Eq2) / 2.0
                            SM1 = (s.WC[il1] - EqualPotAmount) / TSL[il1]
                            SM2 = (s.WC[il2] + EqualPotAmount) / TSL[il2]
                            PF1 = self.soil_profile[il1].SMfromPF(SM1)
                            PF2 = self.soil_profile[il2].SMfromPF(SM2)
                            if abs(Eq1-Eq2) < self.TinyFlow:
                                # sufficient accuracy
                                break
                            elif PF1 > PF2:
                                # suction in top layer 1 is larger ; absolute amount should be larger
                                Eq2 = EqualPotAmount
                            else:
                                # suction in bottom layer 1 is larger ; absolute amount should be reduced
                                Eq1 = EqualPotAmount
                        else:
                            msg = "WATFDGW: Limiting amount iteration in dry flow failed. Are your soil moisture " \
                                  "and conductivity curves decreasing with increase pF?"
                            raise exc.PCSEError(msg)

                # the limit under wet conditions is a unit gradient
                LIMWET[il] = (TSL[il-1]+TSL[il]) / (TSL[il-1]/conductivity[il-1] + TSL[il]/conductivity[il])

            FlowDown = True  # default

            if LIMDRY[il] < 0.0:
                # upward flow (negative !) is limited by fraction of amount required for equilibrium
                FlowMax = max(LIMDRY[il], EqualPotAmount * self.UpwardFlowLimit)
                if il > 0:
                    # upward flow is limited by amount required to bring target layer at equilibrium/field capacity
                    # if (il==2) write (*,*) '2: ',FlowMax, LIMDRY(il), EqualPotAmount * UpwardFlowLimit
                    if self.soil_profile.GroundWater:
                        # soil does not drain below equilibrium with groundwater
                        # FCequil = MAX(WCFC(il-1), EquilWater(il-1))
                        raise NotImplementedError("Groundwater influence not implemented yet.")
                    else:
                        # free drainage
                        FCequil = self.soil_profile[il-1].WCFC

                    TargetLimit = WTRALY[il-1] + FCequil - s.WC[il-1]/delt
                    if TargetLimit > 0.0:
                        # target layer is "dry": below field capacity ; limit upward flow
                        FlowMax = max(FlowMax, -1.0 * TargetLimit)
                        # there is no saturation prevention since upward flow leads to a decrease of WC[il]
                        # instead flow is limited in order to prevent a negative water content
                        FlowMX[il] = max(FlowMax, FlowMX[il+1] + WTRALY[il] - s.WC[il]/delt)
                        FlowDown = False
                    elif self.soil_profile.GroundWater:
                        # target layer is "wet", above field capacity. Since gravity is neglected
                        # in the matrix potential model, this "wet" upward flow is neglected.
                        FlowMX[il] = 0.0
                        FlowDown   = True
                    else:
                        # target layer is "wet", above field capacity, without groundwater
                        # The free drainage model implies that upward flow is rejected here.
                        # Downward flow is enabled and the free drainage model applies.
                        FlowDown = True

            if FlowDown:
                # maximum downward flow rate (LIMWET is always a positive number)
                FlowMax = max(LIMDRY[il], LIMWET[il])
                # this prevents saturation of layer il
                # maximum top boundary flow is bottom boundary flow plus saturation deficit plus sink
                FlowMX[il] = min(FlowMax, FlowMX[il+1] + (self.soil_profile[il].WC0 - s.WC[il])/delt + WTRALY[il])
        # end for

        r.RIN = min(RINPRE, FlowMX[0])

        # contribution of layers to soil evaporation in case of drought upward flow is allowed
        EVSL = np.zeros_like(s.SM)
        for il, layer in enumerate(self.soil_profile):
            if il == 0:
                EVSL[il] = min(r.EVS, (s.WC[il] - layer.WCW) / delt + r.RIN - WTRALY[il])
                EVrest = r.EVS - EVSL[il]
            else:
                Available = max(0.0, (s.WC[il] - layer.WCW)/delt - WTRALY[il])
                if Available >= EVrest:
                    EVSL[il] = EVrest
                    EVrest   = 0.0
                    break
                else:
                    EVSL[il] = Available
                    EVrest   = EVrest - Available
        # reduce evaporation if entire profile becomes airdry
        # there is no evaporative flow through lower boundary of layer NSL
        r.EVS = r.EVS - EVrest

        # Convert contribution of soil layers to EVS as an upward flux
        # evaporative flow (taken positive !!!!) at layer boundaries
        NSL = len(s.SM)
        EVflow = np.zeros_like(FlowMX)
        EVflow[0] = r.EVS
        for il in range(1, NSL):
           EVflow[il] = EVflow[il-1] - EVSL[il-1]
        EVflow[NSL] = 0.0  # see comment above

        # limit downward flows as to not get below field capacity / equilibrium content
        Flow = np.zeros_like(FlowMX)
        r.DWC = np.zeros_like(s.SM)
        Flow[0] = r.RIN - EVflow[0]
        for il, layer in enumerate(self.soil_profile):
            if self.soil_profile.GroundWater:
                # soil does not drain below equilibrium with groundwater
                #WaterLeft = max(self.WCFC[il], EquilWater[il])
                raise NotImplementedError("Groundwater influence not implemented yet.")
            else:
                # free drainage
                WaterLeft = layer.WCFC
            MXLOSS = (s.WC[il] - WaterLeft)/delt               # maximum loss
            Excess = max(0.0, MXLOSS + Flow[il] - WTRALY[il])  # excess of water (positive)
            Flow[il+1] = min(FlowMX[il+1], Excess - EVflow[il+1])  # note that a negative (upward) flow is not affected
            # rate of change
            r.DWC[il] = Flow[il] - Flow[il+1] - WTRALY[il]

        # Flow at the bottom of the profile
        r.BOTTOMFLOW = Flow[-1]

        # # Percolation and Loss.
        # # Equations were derived from the requirement that in the same layer, above and below
        # # depth RD (or RDM), the water content is uniform. Note that transpiration above depth
        # # RD (or RDM) may require a negative percolation (or loss) for keeping the layer uniform.
        # # This is in fact a numerical dispersion. After reaching RDM, this negative (LOSS) can be
        # # prevented by choosing a layer boundary at RDM.
        # RD = self._determine_rooting_depth()
        # ILR, ILM = self.soil_profile.find_layer_indices_for_rooting(RD, self._RDM)
        # if ILR < ILM:
        #     # layer ILR is divided into rooted part (where the sink is) and a below-roots part
        #     # The flow in between is PERC
        #     f1 = self.soil_profile[ILR].Wtop
        #     r.PERC = (1.0-f1) * (Flow[ILR] - WTRALY[ILR]) + f1 * Flow[ILR+1]
        #
        #     # layer ILM is divided as well ; the flow in between is LOSS
        #     f1 = self.soil_profile[ILR].Wpot
        #     r.LOSS = (1.0-f1) * Flow[ILM] + f1 * Flow[ILM+1]
        # elif ILR == ILM:
        #     # depths RD and RDM in the same soil layer: there are three "sublayers":
        #     # - the rooted sublayer with fraction f1
        #     # - between RD and RDM with fraction f2
        #     # - below RDM with fraction f3
        #     # PERC goes from 1->2, LOSS from 2->3
        #     # PERC and LOSS are calculated in such a way that the three sublayers have equal SM
        #     f1 = self.soil_profile[ILR].Wtop
        #     f2 = self.soil_profile[ILM].Wpot
        #     f3 = 1.0 - f1 - f2
        #     r.LOSS = f3 * (Flow[ILR] - WTRALY[ILR]) + (1.0-f3) * Flow[ILR+1]
        #     r.PERC = (1.0-f1) * (Flow[ILR] - WTRALY[ILR]) + f1 * Flow[ILR+1]
        # else:
        #     raise exc.PCSEError("Failure calculating LOSS/PERC")
        #
        # # rates of change in amounts of moisture W and WLOWI
        # r.DW = -sum(WTRALY) - r.EVS - r.PERC + r.RIN
        # r.DWLOW = r.PERC - r.LOSS



        if self.soil_profile.GroundWater:
            # groundwater influence
            # DWBOT = LOSS - Flow[self.NSL+1]
            # DWSUB = Flow[self.NSL+1]
            raise NotImplementedError("Groundwater influence not implemented yet.")

        # Computation of rate of change in surface storage and surface runoff
        # SStmp is the layer of water that cannot infiltrate and that can potentially
        # be stored on the surface. Here we assume that RAIN_NOTINF automatically
        # ends up in the surface storage (and finally runoff).
        SStmp = drv.RAIN + r.RIRR - r.EVW - r.RIN
        # rate of change in surface storage is limited by SSMAX - SS
        r.DSS = min(SStmp, (p.SSMAX - s.SS))
        # Remaining part of SStmp is send to surface runoff
        r.DTSR = SStmp - r.DSS
        # incoming rainfall rate
        r.DRAINT = drv.RAIN

        self._RINold = r.RIN

    @prepare_states
    def integrate(self, day, delt):
        p = self.params
        s = self.states
        k = self.kiosk
        r = self.rates

        # amount of water in soil layers ; soil moisture content
        SM = np.zeros_like(s.SM)
        WC = np.zeros_like(s.WC)
        for il, layer in enumerate(self.soil_profile):
            WC[il] = s.WC[il] + r.DWC[il] * delt
            SM[il] = WC[il] / layer.Thickness
        # NOTE: We cannot replace WC[il] with s.WC[il] above because the kiosk will not
        # be updated since traitlets cannot monitor changes within lists/arrays.
        # So we have to assign:
        s.SM = SM
        s.WC = WC

        # total transpiration
        s.WTRAT += r.WTRA * delt

        # total evaporation from surface water layer and/or soil
        s.EVWT += r.EVW * delt
        s.EVST += r.EVS * delt

        # totals for rainfall, irrigation and infiltration
        s.TOTINF += r.RIN * delt
        s.TOTIRR += r.RIRR * delt

        # surface storage and runoff
        s.SS += r.DSS * delt
        s.TSR += r.DTSR * delt

        # loss of water by outflow through bottom of profile
        s.BOTTOMFLOWT += r.BOTTOMFLOW * delt

        # percolation from rootzone ; interpretation depends on mode
        if self.soil_profile.GroundWater:
            # with groundwater this flow is either percolation or capillary rise
            if r.PERC > 0.0:
                s.PERCT = s.PERCT + r.PERC * delt
            else:
                s.CRT = s.CRT - r.PERC * delt
        else:
            # without groundwater this flow is always called percolation
            # s.PERCT += r.PERC * delt
            s.CRT   = 0.0

        # loss of water by flow from the potential rootzone
        # s.LOSST += r.LOSS * delt

        s.RAINT += self._RAIN

        # change of rootzone subsystem boundary
        RD = self._determine_rooting_depth()
        RDchange = RD - self._RDold
        if abs(RDchange) > 0.001:
            self.soil_profile.determine_rooting_status(RD, self._RDM)
            self.soil_profile.compute_layer_weights(RD, self._RDM)

        # compute summary values for rooted, potentially rooted and unrooted soil compartments
        W = 0.0 ; WAVUPP = 0.0
        WLOW = 0.0 ; WAVLOW = 0.0
        WBOT = 0.0 ; WAVBOT = 0.0
        # get W and WLOW and available water amounts
        for il, layer in enumerate(self.soil_profile):
            W += s.WC[il] * layer.Wtop
            WLOW += s.WC[il] * layer.Wpot
            WBOT += s.WC[il] * layer.Wund
            WAVUPP += (s.WC[il] - layer.WCW) * layer.Wtop
            WAVLOW += (s.WC[il] - layer.WCW) * layer.Wpot
            WAVBOT += (s.WC[il] - layer.WCW) * layer.Wund

        # Update states
        s.W = W
        s.WLOW = WLOW
        s.WWLOW = s.W + s.WLOW
        s.WBOT = WBOT
        s.WAVUPP = WAVUPP
        s.WAVLOW = WAVLOW
        s.WAVBOT = WAVBOT

        # save rooting depth for which layer contents have been determined
        self._RDold = RD

        s.SM_MEAN = s.W/RD

    @prepare_states
    def finalize(self, day):
        s = self.states
        p = self.params
        if self.soil_profile.GroundWater:
            # checksums waterbalance for system Groundwater version
            # WBALRT_GW = TOTINF + CRT + WI - W + WDRT - EVST - TRAT - PERCT
            # WBALTT_GW = SSI + RAINT + TOTIRR + WI - W + WZI - WZ - TRAT - EVWT - EVST - TSR - DRAINT - SS
            pass
        else:
            # checksums waterbalance for system Free Drainage version
            checksum = (p.SSI + s.RAINT + s.TOTIRR + self._WCI - s.WC.sum() -
                        s.WTRAT - s.EVWT - s.EVST - s.TSR - s.BOTTOMFLOWT - s.SS)
            if abs(checksum) > 0.0001:
                msg = "Waterbalance not closing on %s with checksum: %f" % (day, checksum)
                raise exc.WaterBalanceError(msg)

    def _determine_rooting_depth(self):
        """Determines appropriate use of the rooting depth (RD)

        This function includes the logic to determine the depth of the upper (rooted)
        layer of the water balance. See the comment in the code for a detailed description.
        """
        if "RD" in self.kiosk:
            return self.kiosk["RD"]
        else:
            # Hold RD at default value
            return self._default_RD

    def _on_CROP_START(self):
        pass
        self.crop_start = True
        # self.rooted_layer_needs_reset = True

    def _on_CROP_FINISH(self):
        pass
        # self.in_crop_cycle = False
        # self.rooted_layer_needs_reset = True

    def _on_IRRIGATE(self, amount, efficiency):
        self._RIRR = amount * efficiency

    def _setup_new_crop(self, day):
        """Retrieves the crop maximum rootable depth, validates it and updates the layer weights
        in order to have a correct calculation of summary waterbalance states.

        """
        self._RDM = self.parameter_provider["RDMCR"]
        self.soil_profile.validate_max_rooting_depth(self._RDM)
        self.soil_profile.determine_rooting_status(self._default_RD, self._RDM)
        self.soil_profile.compute_layer_weights(self._default_RD, self._RDM)

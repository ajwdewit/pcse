from math import sqrt
import numpy as np
from dotmap import DotMap


from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, HasTraits, List
from ..decorators import prepare_rates, prepare_states
from ..util import limit, Afgen, merge_dict
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject

from .soil_profile import SoilProfile


class WaterBalanceLayered(SimulationObject):

    _default_RD = 15.  # default rooting depth at 10 cm
    RDold = None
    RINold = Float(-99.)
    DSLR = Int(-99)
    RDM = None
    _RIRR = Float(-99)

    # output from layer weights function
    Wtop = None
    Wpot = None
    Wund = None
    ILR = None
    ILM = None

    # Max number of flow iterations
    MaxFlowIter = 50
    TinyFlow = 0.001

    # Maximum upward flow is 50% of amount needed to reach equilibrium between layers
    # see documentatin Kees Rappoldt - page 80
    UpwardFlowLimit = 0.50

    soil_profile = None

    class Parameters(ParamTemplate):
        # SMFCF = List()
        # SM0 = List()
        # SMW = List()
        # CRAIRC = List()
        # SOIL_LAYERS = List()
        # SOIL_PROFILE = Instance(DotMap)
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
        SUMSM = Float(-99)
        PERCT = Float(-99)
        LOSST = Float(-99)
        CRT = Float(-99)
        RAINT = Float(-99)
        # DSLR = Int(-1)
        SM = Instance(np.ndarray)
        WC = Instance(np.ndarray)
        WWLOW = Float(-99)
        WLOW = Float(-99)
        WBOT = Float(-99)
        W = Float(-99)
        SS = Float(-99)

    class RateVariables(RatesTemplate):
        RIN = Float(-99)
        WTRA = Instance(np.ndarray)
        EVS = Float(-99)
        EVW = Float(-99)
        RIRR = Float(-99)
        DWC = Instance(np.ndarray)
        PERC = Float(-99)
        LOSS = Float(-99)
        DRAINT = Float(-99)
        DSS = Float(-99)
        DTSR = Float(-99)

    def initialize(self, day, kiosk, parvalues):

        self.soil_profile = SoilProfile(parvalues)
        self._LayerThickness = [l.Thickness for l in self.soil_profile]
        LayerLowerBoundary = list(np.cumsum(self._LayerThickness))


        # TEMPORARY INPUTS FROM RERUN FILES.
        # Maximum rootable depth, this has to change layer because RDMCR
        # is often not known at this point.
        RDMCR = 125.
        self.RDM = min(RDMCR, max(LayerLowerBoundary))
        self.soil_profile.validate_max_rooting_depth(self.RDM)

        self.params = self.Parameters(parvalues)
        p = self.params

        r = self._LayerWeights(RD=self._default_RD, RDM=self.RDM, LayerTHK=self._LayerThickness)
        self.Wtop, self.Wpot, self.Wund, self.ILR, self.ILM = r
        self.soil_profile.determine_rooting_status(self._default_RD, self.RDM)
        self.soil_profile.compute_layer_weights(self._default_RD, self.RDM)

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

        # within rootzone
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
                WAVUPP += (SM[il] - layer.SMW) * layer.Thickness * layer.Wtot
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
        self._WI    = W
        self._WLOWI = WLOW
        self._WWLOWI = W + WLOW

        # soil evaporation, days since last rain
        self.DSLR = 5 if SM[1] <= (SMW[1] + 0.5 * (SMFCF[1]-SMW[1])) else 1

        # all summation variables of the water balance are set at zero.
        states = {
            "WTRAT": 0., "EVST": 0., "EVWT": 0., "TSR": 0.,
            "RAINT": 0., "WDRT": 0., "TOTINF": 0., "TOTIRR": 0.,
            "SUMSM": 0., "PERCT": 0., "LOSST": 0., "SS":0.,
            "CRT": 0., "RAINT": 0., "WLOW": WLOW, "W": W, "WC": WC, "SM":SM,
            "SS": p.SSI, "WWLOW": W+WLOW, "WBOT":0.
        }
        self.states = self.StateVariables(kiosk, publish=["WC"], **states)

        self.RINold = 0.
        self._RIRR = 0.
        self.RDold = self._default_RD

        # rate variables
        self.rates = self.RateVariables(kiosk)

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

        # Transpiration and maximum soil and surface water evaporation rates
        # are calculated by the crop evapotranspiration module.
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if "TRA" not in self.kiosk:
            r.WTRA = np.zeros_like(s.SM)
            EVWMX = drv.E0
            EVSMX = drv.ES0
        else:
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
            if self.RINold >= 1:
                # If infiltration >= 1cm on previous day assume maximum soil
                # evaporation
                r.EVS = EVSMX
                self._DSLR = 1.
            else:
                # Else soil evaporation is a function days-since-last-rain (DSLR)
                EVSMXT = EVSMX * (sqrt(self.DSLR + 1) - sqrt(self.DSLR))
                r.EVS = min(EVSMX, EVSMXT + self.RINold)
                self.DSLR += 1

        # conductivities and Matric Flux Potentials for all layers
        pF = np.zeros_like(s.SM)
        conductivity = np.zeros_like(s.SM)
        matricfluxpot = np.zeros_like(s.SM)
        for i, layer in enumerate(p.SOIL_LAYERS):
            pF[i] = layer.PFfromSM(s.SM[i])
            conductivity[i] = 10**layer.CONDfromPF(pF[i])
            matricfluxpot[i] = layer.MFPfromPF(pF[i])
            if p.SOIL_PROFILE.GroundWater:
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
            RINPRE = min(p.SOPE, AVAIL)

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
        if p.SOIL_PROFILE.GroundWater:
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
            FlowMX[-1] = max(self.CondFC[-1], conductivity[-1])

        # drainage
        DMAX = 0.0

        LIMDRY =np.zeros_like(s.SM)
        LIMWET =np.zeros_like(s.SM)
        TSL = self._LayerThickness
        for il in reversed(range(len(p.SOIL_LAYERS))):
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
                LIMWET[il] = p.SOIL_PROFILE.SurfaceConductivity
                LIMDRY[il] = 0.0
            else:
                if p.SOIL_LAYERS[il-1] == p.SOIL_LAYERS[il]:
                    # flow rate estimate from gradient in Matric Flux Potential
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
                        Flow1 = 2.0 * (+ MFP1 - p.SOIL_LAYERS[il1].MFPfromPF(PFx)) / TSL[il1]
                        Flow2 = 2.0 * (- MFP2 + p.SOIL_LAYERS[il2].MFPfromPF(PFx)) / TSL[il2]
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
                        raise RuntimeError('WATFDGW: LIMDRY flow iteration failed')
                    LIMDRY[il] = (Flow1 + Flow2) / 2.0

                    if LIMDRY[il] < 0.0:
                        # upward flow rate ; amount required for equal potential is required below
                        Eq1 = -s.WC[il2]; Eq2 = 0.0
                        for z in range(self.MaxFlowIter):
                            EqualPotAmount = (Eq1 + Eq2) / 2.0
                            SM1 = (s.WC[il1] - EqualPotAmount) / TSL[il1]
                            SM2 = (s.WC[il2] + EqualPotAmount) / TSL[il2]
                            PF1 = p.SOIL_LAYERS[il1].SMfromPF(SM1)
                            PF2 = p.SOIL_LAYERS[il2].SMfromPF(SM2)
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
                            raise RuntimeError('WATFDGW Limiting amount iteration failed')

                # the limit under wet conditions is a unit gradient
                LIMWET[il] = (TSL[il-1]+TSL[il]) / (TSL[il-1]/conductivity[il-1] + TSL[il]/conductivity[il])

            FlowDown = True  # default

            if LIMDRY[il] < 0.0:
                # upward flow (negative !) is limited by fraction of amount required for equilibrium
                FlowMax = max(LIMDRY[il], EqualPotAmount * self.UpwardFlowLimit)
                if il > 0:
                    # upward flow is limited by amount required to bring target layer at equilibrium/field capacity
                    # if (il==2) write (*,*) '2: ',FlowMax, LIMDRY(il), EqualPotAmount * UpwardFlowLimit
                    if p.SOIL_PROFILE.GroundWater:
                        # soil does not drain below equilibrium with groundwater
                        # FCequil = MAX(WCFC(il-1), EquilWater(il-1))
                        raise NotImplementedError("Groundwater influence not implemented yet.")
                    else:
                        # free drainage
                        FCequil = self.WCFC[il-1]

                    TargetLimit = r.WTRA[il-1] + FCequil - s.WC[il-1]/delt
                    if TargetLimit > 0.0:
                        # target layer is "dry": below field capacity ; limit upward flow
                        FlowMax = max(FlowMax, -1.0 * TargetLimit)
                        # there is no saturation prevention since upward flow leads to a decrease of WC[il]
                        # instead flow is limited in order to prevent a negative water content
                        FlowMX[il] = max(FlowMax, FlowMX[il+1] + r.WTRA[il] - s.WC[il]/delt)
                        FlowDown = False
                    elif p.SOIL_PROFILE.GroundWater:
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
                FlowMX[il] = min(FlowMax, FlowMX[il+1] + (self.WC0[il] - s.WC[il])/delt + r.WTRA[il])
        # end for
        print(1)

        r.RIN = min(RINPRE, FlowMX[0])

        # contribution of layers to soil evaporation in case of drought upward flow is allowed
        EVSL    = np.zeros_like(s.SM)
        EVSL[0] = min(r.EVS, (s.WC[0] - self.WCW[0])/delt + r.RIN - r.WTRA[0])
        EVrest  = r.EVS - EVSL[0]
        for il in range(1, len(p.SOIL_LAYERS)):
            Available = max(0.0, (s.WC[il] - self.WCW[il])/delt - r.WTRA[il])
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
        EVflow = np.zeros(NSL + 1)
        EVflow[0] = r.EVS
        for il in range(1, NSL):
           EVflow[il] = EVflow[il-1] - EVSL[il-1]
        EVflow[NSL] = 0.0  # see comment above

        # limit downward flows as to not get below field capacity / equilibrium content
        Flow = np.zeros_like(FlowMX)
        r.DWC = np.zeros_like(s.SM)
        Flow[0] = r.RIN - EVflow[0]
        for il in range(NSL):
            if p.SOIL_PROFILE.GroundWater:
                # soil does not drain below equilibrium with groundwater
                #WaterLeft = max(self.WCFC[il], EquilWater[il])
                raise NotImplementedError("Groundwater influence not implemented yet.")
            else:
                # free drainage
                WaterLeft = self.WCFC[il]
            MXLOSS = (s.WC[il] - WaterLeft)/delt               # maximum loss
            Excess = max(0.0, MXLOSS + Flow[il] - r.WTRA[il])  # excess of water (positive)
            Flow[il+1] = min(FlowMX[il+1], Excess - EVflow[il+1])  # note that a negative (upward) flow is not affected
            # rate of change
            r.DWC[il] = Flow[il] - Flow[il+1] - r.WTRA[il]

        # Percolation and Loss.
        # Equations were derived from the requirement that in the same layer, above and below
        # depth RD (or RDM), the water content is uniform. Note that transpiration above depth
        # RD (or RDM) may require a negative percolation (or loss) for keeping the layer uniform.
        # This is in fact a numerical dispersion. After reaching RDM, this negative (LOSS) can be
        # prevented by choosing a layer boundary at RDM.
        if self.ILR < self.ILM:
            # layer ILR is divided into rooted part (where the sink is) and a below-roots part
            # The flow in between is PERC
            f1 = self.Wtop[self.ILR]
            r.PERC = (1.0-f1) * (Flow[self.ILR] - r.WTRA[self.ILR]) + f1 * Flow[self.ILR+1]

            # layer ILM is divided as well ; the flow in between is LOSS
            f1 = self.Wpot[self.ILM]
            r.LOSS = (1.0-f1) * Flow[self.ILM] + f1 * Flow[self.ILM+1]
        elif self.ILR == self.ILM:
            # depths RD and RDM in the same soil layer: there are three "sublayers":
            # - the rooted sublayer with fraction f1
            # - between RD and RDM with fraction f2
            # - below RDM with fraction f3
            # PERC goes from 1->2, LOSS from 2->3
            # PERC and LOSS are calculated in such a way that the three sublayers have equal SM
            f1 = self.Wtop[self.ILR]
            f2 = self.Wpot[self.ILM]
            f3 = 1.0 - f1 - f2
            r.LOSS = f3 * (Flow[self.ILR] - r.WTRA[self.ILR]) + (1.0-f3) * Flow[self.ILR+1]
            r.PERC = (1.0-f1) * (Flow[self.ILR] - r.WTRA[self.ILR]) + f1 * Flow[self.ILR+1]
        else:
            raise RuntimeError("Failure calculating LOSS/PERC")

        # rates of change in amounts of moisture W and WLOWI
        r.DW = -sum(r.WTRA) - r.EVS - r.PERC + r.RIN
        r.DWLOW = r.PERC - r.LOSS

        if p.SOIL_PROFILE.GroundWater:
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

    @prepare_states
    def integrate(self, day, delt):
        p = self.params
        s = self.states
        k = self.kiosk
        r = self.rates

        # amount of water in soil layers ; soil moisture content
        for il, layer in enumerate(p.SOIL_LAYERS):
            s.WC[il] = s.WC[il] + r.DWC[il] * delt
            s.SM[il] = s.WC[il] / layer.Thickness

        # total transpiration
        s.WTRAT = s.WTRAT + sum(r.WTRA) * delt

        # total evaporation from surface water layer and/or soil
        s.EVWT = s.EVWT + r.EVW * delt
        s.EVST = s.EVST + r.EVS * delt

        # totals for rainfall, irrigation and infiltration
        s.TOTINF = s.TOTINF + r.RIN * delt
        s.TOTIRR = s.TOTIRR + r.RIRR * delt

        # surface storage and runoff
        s.SS    += r.DSS * delt
        s.TSR   += r.DTSR * delt

        # amount of water in rooted zone
        s.W      = s.W + r.DW * delt
        # WAVUPP = s.WAVUPP + r.DW * delt

        # amount of water in unrooted, lower part of rootable zone
        s.WLOW   = s.WLOW   + r.DWLOW*delt
        #s.WAVLOW = s.WAVLOW + r.DWLOW*delt

        # total amount of water in the whole rootable zone
        s.WWLOW = s.W + s.WLOW

        # and in layered soil below RDM
        # s.WBOT = s.WBOT + r.DWBOT * delt

        # percolation from rootzone ; interpretation depends on mode
        if p.SOIL_PROFILE.GroundWater:
            # with groundwater this flow is either percolation or capillairy rise
            if r.PERC > 0.0:
                s.PERCT = s.PERCT + r.PERC * delt
            else:
                s.CRT = s.CRT - r.PERC * delt
        else:
            # without groundwater this flow is always called percolation
            s.PERCT = s.PERCT + r.PERC * delt
            s.CRT   = 0.0

        # loss of water by flow from the potential rootzone
        s.LOSST = s.LOSST + r.LOSS * delt


        #----------------------------------------------
        # change of rootzone subsystem boundary
        #----------------------------------------------
        # calculation of amount of soil moisture in new rootzone
        RD = self._determine_rooting_depth()
        RDchange = RD - self.RDold

        if abs(RDchange) > 0.001:
            # roots have grown find new values ; overwrite W, WLOW, WAVUPP, WAVLOW, WBOT, WAVBOT
            r = self._LayerWeights(RD=RD, RDM=self.RDM, LayerTHK=self._LayerThickness)
            self.Wtop, self.Wpot, self.Wund, self.ILR, self.ILM = r

            WOLD = s.W
            W    = 0.0 #; WAVUPP = 0.0
            WLOW = 0.0 #; WAVLOW = 0.0
            # WBOT = 0.0 ; WAVBOT = 0.0
            # get W and WLOW and available water amounts
            for il in range(len(p.SOIL_LAYERS)):
                W    += s.WC[il] * self.Wtop[il]
                WLOW += s.WC[il] * self.Wpot[il]
                # WBOT += s.WC[il] * self.Wund[il]
                #
                # WAVUPP += (s.WC[il] - self.WCW[il]) * self.Wtop[il]
                # WAVLOW += (s.WC[il] - self.WCW[il]) * self.Wpot[il]

            # Update states
            s.W = W
            s.WLOW = WLOW
            s.WWLOW = s.W + s.WLOW

            # save rooting depth for which layer contents have been determined
            self.RDold = RD

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

    def _LayerWeights(self,RD, RDM, LayerTHK):
        # !-----------------------------------------------------------------------
        # !
        # ! weight factors for rooted- and sub-layer calculations
        # ! RD       - the rooting depth
        # ! RDM      - max rooting depth
        # ! LayerTHK - the layerthickness
        # ! LayerLB  - the Lower Boundaries of the NL soil layers
        # ! NL       - number of layers
        # ! LRD      - deepest layer containing roots
        # ! LRM      - deepest layer that will contain roots
        Wtop = []  # weights for contribution to rootzone
        Wpot = []  # weights for contribution to potentially rooted zone
        Wund = []  # weights for contribution to never rooted layers
        # find deepest layer with roots
        NL = len(LayerTHK)
        LayerLB = list(np.cumsum(LayerTHK))
        LRD = NL
        LRM = NL
        for il in reversed(range(NL)):
            if LayerLB[il] >= RD: LRD = il
            if LayerLB[il] >= RDM: LRM = il

        for il in range(NL):
            if il < LRD:
                # rooted layer
                Wtop.append(1.0)
                Wpot.append(0.0)
                Wund.append(0.0)
            elif il == LRD and il < LRM:
                # partly rooted, at the end fully rooted
                Wtop.append(1.0 - (LayerLB[il]-RD) / LayerTHK[il])
                Wpot.append(1.0 - Wtop[il])
                Wund.append(0.0)
            elif il == LRD and il == LRM:
                # partly rooted, at the end partly rooted
                Wtop.append(1.0 - (LayerLB[il]-RD ) / LayerTHK[il])
                Wund.append((LayerLB[il]-RDM) / LayerTHK[il])
                Wpot.append(1.0 - Wund[il] - Wtop[il])
            elif il < LRM:
                # not rooted, at the end fully rooted
                Wtop.append(0.0)
                Wpot.append(1.0)
                Wund.append(0.0)
            elif il == LRM:
                # not rooted, at the end partly rooted
                Wtop.append(0.0)
                Wund.append((LayerLB[il]-RDM) / LayerTHK[il])
                Wpot.append(1.0 - Wund[il])
            else:
                # never rooted
                Wtop.append(0.0)
                Wpot.append(0.0)
                Wund.append(1.0)

        return Wtop, Wpot, Wund, LRD, LRM

from math import sqrt
import numpy as np
from dotmap import DotMap


from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, HasTraits, List
from ..decorators import prepare_rates, prepare_states
from ..util import limit, Afgen, merge_dict
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject


class pFCurve(Afgen):
    """Pf curve should check that:
    - length of array is even
    - has at least 3 xy pairs
    - pF should start at pF = -1 and end at pF=6
    - SM/cond values should decrease with increase pF

    """
    pass

class MFPCurve(Afgen):
    """Computes Matrix Flow Potential using a Gaussian integration over the pfCurve

    """
    elog10 = 2.302585092994
    Pgauss = (0.0469100770, 0.2307653449, 0.5000000000, 0.7692346551, 0.9530899230)
    Wgauss = (0.1184634425, 0.2393143352, 0.2844444444, 0.2393143352, 0.1184634425)

    def __init__(self, SMfromPF, CONDfromPF):
        SMfromPF = np.array(SMfromPF)
        CONDfromPF = pFCurve(CONDfromPF)
        MFPfromPF = np.zeros_like(SMfromPF)

        # zero MFP at highest pF
        MFPfromPF[-1] = 0.
        MFPfromPF[-2] = SMfromPF[-2]

        for ip in range(len(SMfromPF)-3, 0, -2):
            # Copy corresponding pF value
            MFPfromPF[ip-1] = SMfromPF[ip-1]
            # get integral over PF range
            add = 0.0
            DeltaPF = SMfromPF[ip + 1] - SMfromPF[ip - 1]
            for i in range(len(self.Pgauss)):
                PFg = SMfromPF[ip - 1] + self.Pgauss[i] * DeltaPF
                CON = 10.0 ** CONDfromPF(PFg)
                add += CON * 10.0 ** PFg * self.elog10 * self.Wgauss[i]
            MFPfromPF[ip] = add * DeltaPF + MFPfromPF[ip + 2]
        Afgen.__init__(self, MFPfromPF)


class SoilLayer(HasTraits):
    """Contains the intrinsic and derived properties for each soil layers
    """
    SMfromPF = Instance(pFCurve)  # soil moisture content as function of pF
    CONDfromPF = Instance(pFCurve)  # conductivity as function of pF
    PFfromSM = Instance(Afgen)  # Inverse from SMfromPF
    MFPfromPF = Instance(MFPCurve)  # Matrix Flux Potential as function of pF
    CRAIRC = Float()  # Critical air content
    SMsat = Float()  # Saturation soil moisture content
    Thickness = Float()

    def __init__(self, layer):
        self.SMfromPF = pFCurve(layer.SMfromPF)
        self.CONDfromPF = pFCurve(layer.CONDfromPF)
        self.PFfromSM = self._invert_pF(layer.SMfromPF)
        self.CRAIRC = layer.CRAIRC
        self.SMsat = self.SMfromPF(-1.0)
        self.MFPfromPF = MFPCurve(layer.SMfromPF, layer.CONDfromPF)
        self.Thickness = layer.Thickness
        self.SoilID = layer.SoilID
        
        # compute hash value of this layer
        #self._hash = hash((tuple(layer.SMfromPF), tuple(layer.CONDfromPF)))
        self._hash = hash(self.SoilID)

    def _invert_pF(self, SMfromPF):
        """Inverts the SMfromPF table to get pF from SM
        """
        l = []
        for t in zip(reversed(SMfromPF[1::2]), reversed(SMfromPF[0::2])):
            l.extend(t)
        return Afgen(l)
    
    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and
            self._hash == other._hash
        )


class WaterBalanceLayered(SimulationObject):

    _default_RD = 10.  # default rooting depth at 10 cm
    RINold = Float(-99.)
    DSLR = Int(-99)
    _RIRR = Float(-99)

    WC0 = None  # TODO: these 5 lines are soil properties: move to soil layer class
    WCW = None
    WCFC = None
    CondFC = None
    CondK0 = None

    # Max number of flow iterations
    MaxFlowIter = 50
    TinyFlow = 0.001

    # Maximum upward flow is 50% of amount needed to reach equilibrium between layers
    # see documentatin Kees Rappoldt - page 80
    UpwardFlowLimit = 0.50

    class Parameters(ParamTemplate):
        SMFCF = List()
        SM0 = List()
        SMW = List()
        CRAIRC = List()
        SOIL_LAYERS = List()
        SOIL_PROFILE = Instance(DotMap)
        IFUNRN = Int(-99)
        NOTINF = Float(-99.)

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
        DRAINT = Float(-99)
        # DSLR = Int(-1)
        SM = Instance(np.ndarray)
        WC = Instance(np.ndarray)
        # WWLOW = Float(-99)
        WLOW = Float(-99)
        W = Float(-99)
        SS = Float(-99)

    class RateVariables(RatesTemplate):
        RIN = Float(-99)
        WTRA = Instance(np.ndarray)
        EVS = Float(-99)
        EVW = Float(-99)
        RIRR = Float(-99)


    def initialize(self, day, kiosk, parvalues):

        soil_profile = DotMap(parvalues["SoilProfileDescription"])
        SOIL_LAYERS = []
        for layer_properties in soil_profile.SoilLayers:
            soil_layer = SoilLayer(layer_properties)
            SOIL_LAYERS.append(soil_layer)

        SMFCF = [l.SMfromPF(soil_profile.PFFieldCapacity) for l in SOIL_LAYERS]
        SM0 = [l.SMsat for l in SOIL_LAYERS]
        SMW = [l.SMfromPF(soil_profile.PFWiltingPoint) for l in SOIL_LAYERS]
        CRAIRC = [l.CRAIRC for l in SOIL_LAYERS]

        parvalues.set_derived("SOIL_LAYERS", SOIL_LAYERS)
        parvalues.set_derived("SOIL_PROFILE", soil_profile)
        parvalues.set_derived("SM0", SM0)
        parvalues.set_derived("SMW", SMW)
        parvalues.set_derived("SMFCF", SMFCF)
        parvalues.set_derived("CRAIRC", CRAIRC)
        self.params = self.Parameters(parvalues)

        self._LayerThickness = [l.Thickness for l in SOIL_LAYERS]
        LayerLowerBoundary = list(np.cumsum(self._LayerThickness))

        # TEMPORARY INPUTS FROM RERUN FILES.
        # Maximum rootable depth, this has to change layer because RDMCR
        # is often not known at this point.
        RDMCR = 125.
        RDM = min(RDMCR, max(LayerLowerBoundary))
        # SMLIM limit on moisture content in top layer
        SMLIM = 0.36
        WAV = 22.

        r = self._LayerWeights(RD=self._default_RD, RDM=RDM, LayerTHK=self._LayerThickness)
        Wtop, Wpot, Wund, ILR, ILM = r

        if soil_profile.GroundWater:
            raise NotImplementedError("Groundwater influence not yet implemented.")
        else:
            # AVMAX -  maximum available content of layer(s)
            # This is calculated first to achieve an even distribution of water in the rooted top
            # if WAV is small. Note the separate limit for initial SM in the rooted zone.
            TOPLIM = 0.0
            LOWLIM = 0.0
            AVMAX = []
            for il in range(ILM+1):
                # determine maximum content for this layer
                if il <= ILR:
                    # Check whether SMLIM is within boundaries (TvdW, 24-jul-97)
                    SML = limit(SMW[il], SM0[il], SMLIM)
                    # the IAIRDU mess is disabled here.
                    # if (IAIRDU.EQ.1) SML = SM0(il)
                    AVMAX.append( (SML - SMW[il]) * self._LayerThickness[il])   # available in cm
                    # also if partly rooted, the total layer capacity counts in TOPLIM
                    # this means the water content of layer ILR is set as if it would be
                    # completely rooted. This water will become available after a little
                    # root growth and through numerical mixing each time step.
                    TOPLIM = TOPLIM + AVMAX[il]
                else:
                    # below the rooted zone the maximum is saturation (see code for WLOW in one-layer model)
                    # again the full layer capacity adds to LOWLIM.
                    SML = SM0[il]
                    AVMAX.append((SML-SMW[il]) * self._LayerThickness[il])   # available in cm
                    LOWLIM = LOWLIM + AVMAX[il]


        if WAV <= 0.0:
            # no available water
            TOPRED = 0.0
            LOWRED = 0.0
        elif WAV <= TOPLIM:
            # available water fits in layer(s) 1..ILR, these layers are rooted or almost rooted
            # reduce amounts with ratio WAV / TOPLIM
            TOPRED = WAV / TOPLIM
            LOWRED = 0.0
        elif WAV < TOPLIM + LOWLIM:
            # available water fits in potentially rooted layer
            # rooted zone is filled at capacity ; the rest reduced
            TOPRED = 1.0
            LOWRED = (WAV-TOPLIM) / LOWLIM
        else:
            # water does not fit ; all layers "full"
            TOPRED = 1.0
            LOWRED = 1.0

        # within rootzone
        W = 0.0    ; WAVUPP = 0.0
        WLOW = 0.0 ; WAVLOW = 0.0
        SM = np.zeros(len(self._LayerThickness))
        for il in range(ILR+1):
            # Part of the water assigned to ILR may not actually be in the rooted zone, but it will
            # be available shortly through root growth (and through numerical mixing).
            SM[il] = SMW[il] + AVMAX[il] * TOPRED / self._LayerThickness[il]
            W      = W    + SM[il] * self._LayerThickness[il] * Wtop[il]
            WLOW   = WLOW + SM[il] * self._LayerThickness[il] * Wpot[il]
            # available water
            WAVUPP = WAVUPP + (SM[il]-SMW[il]) * self._LayerThickness[il] * Wtop[il]
            WAVLOW = WAVLOW + (SM[il]-SMW[il]) * self._LayerThickness[il] * Wpot[il]
        # between initial and maximum rooting depth. In case RDM is not a layer boundary (it should be!!)
        # layer ILM contains additional water in unrooted part. Only rooted part contributes to WAV.
        for il in range(ILR+1, ILM+1):
            SM[il] = SMW[il] + AVMAX[il] * LOWRED / self._LayerThickness[il]
            WLOW = WLOW + SM[il] * self._LayerThickness[il] * Wpot[il]
            # available water
            WAVLOW = WAVLOW + (SM[il]-SMW[il]) * self._LayerThickness[il] * Wpot[il]
        # below the maximum rooting depth
        for il in range(ILM+1, len(self._LayerThickness)):
            SM[il] = SMW[il]
   
        # set groundwater depth far away for clarity ; this prevents also
        # the root routine to stop root growth when they reach the groundwater
        ZT = 999.0

        WC = np.zeros_like(SM)
        self.WC0 = np.zeros_like(SM)  #TODO: these 5 lines are soil properties: move to soil layer class
        self.WCW = np.zeros_like(SM)
        self.WCFC = np.zeros_like(SM)
        self.CondFC = np.zeros_like(SM)
        self.CondK0 = np.zeros_like(SM)
        for il, layer in enumerate(SOIL_LAYERS):
            WC[il]     = SM[il] * self._LayerThickness[il]  # State variable!
            self.WC0[il]    = SM0[il] * self._LayerThickness[il]
            self.WCW[il]    = SMW[il] * self._LayerThickness[il]
            self.WCFC[il]   = SMFCF[il] * self._LayerThickness[il]
            self.CondFC[il] = 10.0**layer.CONDfromPF(soil_profile.PFFieldCapacity)
            self.CondK0[il] = 10.0**layer.CONDfromPF(-1.0)

        # rootzone and subsoil water
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
            "CRT": 0., "DRAINT": 0., "WLOW": WLOW, "W": W, "WC": WC, "SM":SM
        }
        self.states = self.StateVariables(kiosk, publish=["WC"], **states)

        # rate variables
        self.RINold = 0.
        self._RIRR = 0.
        self.rates = self.RateVariables(kiosk)

    def calc_rates(self, day, drv):
        p = self.params
        s = self.states
        k = self.kiosk
        r = self.rates

        DELT = 1.0

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

                    TargetLimit = r.WTRA[il-1] + FCequil - s.WC[il-1]/DELT
                    if TargetLimit > 0.0:
                        # target layer is "dry": below field capacity ; limit upward flow
                        FlowMax = max(FlowMax, -1.0 * TargetLimit)
                        # there is no saturation prevention since upward flow leads to a decrease of WC[il]
                        # instead flow is limited in order to prevent a negative water content
                        FlowMX[il] = max(FlowMax, FlowMX(il+1) + k.TRA[il] - s.WC[il]/DELT)
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
                FlowMX[il] = min(FlowMax, FlowMX[il+1] + (self.WC0[il] - s.WC[il])/DELT + r.WTRA[il])
        # end for
        print(1)







    def integrate(self, day, delt):
        pass

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

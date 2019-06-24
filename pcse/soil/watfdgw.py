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
    - pF should start at pF = -1 and end at pH=6
    - SM/cond values should decrease with increase pF

    """
    pass

class MFPCurve(Afgen):
    """Computes Matrix Flow Potential using the P

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
    PFfromSM = Instance(pFCurve)  # Inverse from SMfromPF
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

    def _invert_pF(self, SMfromPF):
        """Inverts the SMfromPF table to get pF from SM
        """
        l = []
        for t in zip(reversed(SMfromPF[1::2]), reversed(SMfromPF[0::2])):
            l.extend(t)
        return pFCurve(l)


class WaterBalanceLayered(SimulationObject):

    _default_RD = 10.  # default rooting depth at 10 cm

    class Parameters(ParamTemplate):
        SMFCF = List()
        SM0 = List()
        SMW = List()
        CRAIRC = List()
        SOIL_LAYERS = List()


    def initialize(self, day, kiosk, parvalues):

        groundwater = False
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
        parvalues.set_derived("SM0", SM0)
        parvalues.set_derived("SMW", SMW)
        parvalues.set_derived("SMFCF", SMFCF)
        parvalues.set_derived("CRAIRC", CRAIRC)
        self.params = self.Parameters(parvalues)

        LayerThickness = [l.Thickness for l in SOIL_LAYERS]
        LayerLowerBoundary = list(np.cumsum(LayerThickness))

        # TEMPORARY INPUTS FROM RERUN FILES.
        # Maximum rootable depth, this has to change layer because RDMCR
        # is often not known at this point.
        RDMCR = 125.
        RDM = min(RDMCR, max(LayerLowerBoundary))
        # SMLIM limit on moisture content in top layer
        SMLIM = 0.36
        WAV = 22.

        r = self._LayerWeights(RD=self._default_RD, RDM=RDM, LayerTHK=LayerThickness)
        Wtop, Wpot, Wund, ILR, ILM = r

        if groundwater:
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
                    AVMAX.append( (SML - SMW[il]) * LayerThickness[il])   # available in cm
                    # also if partly rooted, the total layer capacity counts in TOPLIM
                    # this means the water content of layer ILR is set as if it would be
                    # completely rooted. This water will become available after a little
                    # root growth and through numerical mixing each time step.
                    TOPLIM = TOPLIM + AVMAX[il]
                else:
                    # below the rooted zone the maximum is saturation (see code for WLOW in one-layer model)
                    # again the full layer capacity adds to LOWLIM.
                    SML = SM0[il]
                    AVMAX.append((SML-SMW[il]) * LayerThickness[il])   # available in cm
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
        SM = np.zeros(len(LayerThickness))
        for il in range(ILR+1):
            # Part of the water assigned to ILR may not actually be in the rooted zone, but it will
            # be available shortly through root growth (and through numerical mixing).
            SM[il] = SMW[il] + AVMAX[il] * TOPRED / LayerThickness[il]
            W      = W    + SM[il] * LayerThickness[il] * Wtop[il]
            WLOW   = WLOW + SM[il] * LayerThickness[il] * Wpot[il]
            # available water
            WAVUPP = WAVUPP + (SM[il]-SMW[il]) * LayerThickness[il] * Wtop[il]
            WAVLOW = WAVLOW + (SM[il]-SMW[il]) * LayerThickness[il] * Wpot[il]
        # between initial and maximum rooting depth. In case RDM is not a layer boundary (it should be!!)
        # layer ILM contains additional water in unrooted part. Only rooted part contributes to WAV.
        for il in range(ILR+1, ILM+1):
            SM[il] = SMW[il] + AVMAX[il] * LOWRED / LayerThickness[il]
            WLOW = WLOW + SM[il] * LayerThickness[il] * Wpot[il]
            # available water
            WAVLOW = WAVLOW + (SM[il]-SMW[il]) * LayerThickness[il] * Wpot[il]
        # below the maximum rooting depth
        for il in range(ILM+1, len(LayerThickness)):
            SM[il] = SMW[il]
   
        # set groundwater depth far away for clarity ; this prevents also
        # the root routine to stop root growth when they reach the groundwater
        ZT = 999.0

        WC = np.zeros_like(SM)
        WC0 = np.zeros_like(SM)
        WCW = np.zeros_like(SM)
        WCFC = np.zeros_like(SM)
        CondFC = np.zeros_like(SM)
        CondK0 = np.zeros_like(SM)
        for il, layer in enumerate(SOIL_LAYERS):
            WC[il]     = SM[il] * LayerThickness[il]  # State variable!
            WC0[il]    = SM0[il] * LayerThickness[il]
            WCW[il]    = SMW[il] * LayerThickness[il]
            WCFC[il]   = SMFCF[il] * LayerThickness[il]
            CondFC[il] = 10.0**layer.CONDfromPF(soil_profile.PFFieldCapacity)
            CondK0[il] = 10.0**layer.CONDfromPF(-1.0)

        # rootzone and subsoil water
        WI    = W
        WLOWI = WLOW
        WWLOW = W + WLOW

        # soil evaporation, days since last rain
        self.DLSR = 5 if SM[1] <= (SMW[1] + 0.5 * (SMFCF[1]-SMW[1])) else 1

        # all summation variables of the water balance are set at zero.
        TRAT  = 0.  ;  EVST  = 0.  ;  EVWT   = 0.  ; TSR    = 0.
        RAINT = 0.  ;  WDRT  = 0.  ;  TOTINF = 0.  ; TOTIRR = 0.
        SUMSM = 0.  ;  PERCT = 0.  ;  LOSST  = 0.
        CRT   = 0.  ; DRAINT = 0.  #  added for groundwater version

        for k, v in locals():
            setattr(self, k, v)

    def calc_rates(self, day, drv):
        pass

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

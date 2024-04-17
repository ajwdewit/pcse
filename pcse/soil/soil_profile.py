from math import sqrt
import numpy as np
from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, HasTraits, List
from ..util import limit, Afgen, merge_dict, DotMap

from .. import exceptions as exc


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

        for ip in range(len(SMfromPF) - 3, 0, -2):
            # Copy corresponding pF value
            MFPfromPF[ip - 1] = SMfromPF[ip - 1]
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
    """Contains the intrinsic and derived properties for each soil layers as required by the multilayer
    waterbalance and SNOMIN.

    :param layer: a soil layer definition providing parameter values for this layer, see table below.
    :param PFFieldcapacity: the pF value for defining Field Capacity
    :param PFWiltingPoint: the pF value for defining the Wilting Point

    The following properties have to be defined for each layer.

    =============== ================================================================  =======================
    Name             Description                                                      Unit
    =============== ================================================================  =======================
    CONDfromPF      Table function of the 10-base logarithm of the unsaturated
                    hydraulic conductivity as a function of pF.                       log10(cm water d-1), -
    SMfromPF        Table function that describes the soil moisture
                    content as a function of pF                                       m3 water m-3 soil, -
    Thickness       Layer thickness                                                   cm
    FSOMI           Initial fraction of soil organic matter in soil                   kg OM kg-1 soil
    CNRatioSOMI     Initial C:N ratio of soil organic matter                          kg C kg-1 N
    RHOD            Bulk density of the soil                                          g soil cm-3 soil
    Soil_pH         pH of the soil layer                                              -
    CRAIRC          Critical air content for aeration for root system                 m3 air m-3 soil
    =============== ================================================================  =======================


    Based on the soil layer definition the following properties are derived from the parameters in the
    table above.

    =============== ================================================================  =======================
     Name             Description                                                      Unit
    =============== ================================================================  =======================
    PFfromSM         Afgen table providing the inverted SMfromPF curve                 m3 water m-3 soil, -
    MFPfromPF        AfTen table describing the Matric Flux Potential as a             cm2 d-1
                     function of the hydraulic head (pF).
    SM0              The volumetric soil moisture content at saturation (pF = -1)      m3 water m-3 soil
    SMW              The volumetric soil moisture content at wilting point             m3 water m-3 soil
    SMFCF            The volumetric soil moisture content at field capacity            m3 water m-3 soil
    WC0              The soil moisture amount (cm) at saturation (pF = -1)             cm water
    WCW              The soil moisture amount (cm) at wilting point                    cm water
    WCFC             The soil moisture amount (cm) at field capacity                   cm water
    CondFC           Soil hydraulic conductivity at field capacity                     cm water d-1
    CondK0           Soil hydraulic conductivity at saturation                         cm water d-1
    =============== ================================================================  =======================

    Finally `rooting_status` is initialized to None (not yet known at initialization).
    """
    SMfromPF = Instance(pFCurve)  # soil moisture content as function of pF
    CONDfromPF = Instance(pFCurve)  # conductivity as function of pF
    PFfromSM = Instance(Afgen)  # Inverse from SMfromPF
    MFPfromPF = Instance(MFPCurve)  # Matrix Flux Potential as function of pF
    CNRatioSOMI = Float()  # Initial C:N ratio of soil organic matter
    FSOMI = Float()  # Initial fraction of soil organic matter in soil
    RHOD = Float()  # Bulk density of the soil
    Soil_pH = Float()  # pH of the soil layer
    CRAIRC = Float()  # Critical air content
    Thickness = Float()
    rooting_status = Enum(["rooted","partially rooted","potentially rooted","never rooted",])

    def __init__(self, layer, PFFieldCapacity, PFWiltingPoint):
        self.SMfromPF = pFCurve(layer.SMfromPF)
        self.CONDfromPF = pFCurve(layer.CONDfromPF)
        self.PFfromSM = self._invert_pF(layer.SMfromPF)
        self.MFPfromPF = MFPCurve(layer.SMfromPF, layer.CONDfromPF)
        self.CNRatioSOMI = layer.CNRatioSOMI
        self.FSOMI = layer.FSOMI
        self.RHOD = layer.RHOD
        self.CRAIRC = layer.CRAIRC
        self.Soil_pH = layer.Soil_pH

        if 5 <= layer.Thickness <= 250:
            self.Thickness = layer.Thickness
        else:
            msg = "Soil layer should have thickness between 5 and 250 cm. Current value: %f" % layer.Thickness
            raise exc.PCSEError(msg)

        self.SM0 = self.SMfromPF(-1.0)
        self.SMFCF = self.SMfromPF(PFFieldCapacity)
        self.SMW = self.SMfromPF(PFWiltingPoint)
        self.WC0 = self.SM0 * self.Thickness
        self.WCW = self.SMW * self.Thickness
        self.WCFC = self.SMFCF * self.Thickness
        self.CondFC = 10.0 ** self.CONDfromPF(PFFieldCapacity)
        self.CondK0 = 10.0 ** self.CONDfromPF(-1.0)
        self.rooting_status = None

        # compute hash value of this layer based on pF curves for conductivity and SM
        self._hash = hash((tuple(layer.SMfromPF), tuple(layer.CONDfromPF)))

    @property
    def Thickness_m(self):
        return self.Thickness * 1e-2

    @property
    def RHOD_kg_per_m3(self):
        return self.RHOD * 1e-3 * 1e06

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


class SoilProfile(list):
    """A component that represents the soil column as required by the multilayer waterbalance and SNOMIN.

    :param parvalues: a ParameterProvider instance that will be used to retrieve the description of the
        soil profile which is assumed to be available under the key `SoilProfileDescription`.

    This class is basically a container that stores `SoilLayer` instances with some additional logic mainly related
    to root growth, root status and root water extraction. For detailed information on the properties of the soil
    layers have a look at the description of the `SoilLayer` class.

    The description of the soil profile can be defined in YAML format with an example of the structure
    given below. In this case first three soil physical types are defined under the section `SoilLayerTypes`.
    Next, these SoilLayerTypes are used to define an actual soil profile using two upper layers of 10 cm of
    type `TopSoil`, three layers of 10, 20 and 30 cm of type `MidSoil`, a lower layer of 45 cm of type `SubSoil`
    and finally a SubSoilType with a thickness of 200 cm. Only the top 3 layers contain a certain amaount of
    organic carbon (FSOMI).

    An example of a data structure for the soil profile::

        SoilLayerTypes:
            TopSoil: &TopSoil
                SMfromPF: [-1.0,     0.366,
                            1.0,     0.338,
                            1.3,     0.304,
                            1.7,     0.233,
                            2.0,     0.179,
                            2.3,     0.135,
                            2.4,     0.123,
                            2.7,     0.094,
                            3.0,     0.073,
                            3.3,     0.059,
                            3.7,     0.046,
                            4.0,     0.039,
                            4.17,    0.037,
                            4.2,     0.036,
                            6.0,     0.02]
                CONDfromPF: [-1.0,     1.8451,
                              1.0,     1.02119,
                              1.3,     0.51055,
                              1.7,    -0.52288,
                              2.0,    -1.50864,
                              2.3,    -2.56864,
                              2.4,    -2.92082,
                              2.7,    -4.01773,
                              3.0,    -5.11919,
                              3.3,    -6.22185,
                              3.7,    -7.69897,
                              4.0,    -8.79588,
                              4.17,   -9.4318,
                              4.2,    -9.5376,
                              6.0,   -11.5376]
                CRAIRC:  0.090
                CNRatioSOMI: 9.0
                RHOD: 1.406
                Soil_pH: 7.4
                SoilID: TopSoil
            MidSoil: &MidSoil
                SMfromPF: [-1.0,     0.366,
                            1.0,     0.338,
                            1.3,     0.304,
                            1.7,     0.233,
                            2.0,     0.179,
                            2.3,     0.135,
                            2.4,     0.123,
                            2.7,     0.094,
                            3.0,     0.073,
                            3.3,     0.059,
                            3.7,     0.046,
                            4.0,     0.039,
                            4.17,    0.037,
                            4.2,     0.036,
                            6.0,     0.02]
                CONDfromPF: [-1.0,     1.8451,
                              1.0,     1.02119,
                              1.3,     0.51055,
                              1.7,    -0.52288,
                              2.0,    -1.50864,
                              2.3,    -2.56864,
                              2.4,    -2.92082,
                              2.7,    -4.01773,
                              3.0,    -5.11919,
                              3.3,    -6.22185,
                              3.7,    -7.69897,
                              4.0,    -8.79588,
                              4.17,   -9.4318,
                              4.2,    -9.5376,
                              6.0,   -11.5376]
                CRAIRC:  0.090
                CNRatioSOMI: 9.0
                RHOD: 1.406
                Soil_pH: 7.4
                SoilID: MidSoil_10
            SubSoil: &SubSoil
                SMfromPF: [-1.0,     0.366,
                            1.0,     0.338,
                            1.3,     0.304,
                            1.7,     0.233,
                            2.0,     0.179,
                            2.3,     0.135,
                            2.4,     0.123,
                            2.7,     0.094,
                            3.0,     0.073,
                            3.3,     0.059,
                            3.7,     0.046,
                            4.0,     0.039,
                            4.17,    0.037,
                            4.2,     0.036,
                            6.0,     0.02]
                CONDfromPF: [-1.0,     1.8451,
                              1.0,     1.02119,
                              1.3,     0.51055,
                              1.7,    -0.52288,
                              2.0,    -1.50864,
                              2.3,    -2.56864,
                              2.4,    -2.92082,
                              2.7,    -4.01773,
                              3.0,    -5.11919,
                              3.3,    -6.22185,
                              3.7,    -7.69897,
                              4.0,    -8.79588,
                              4.17,   -9.4318,
                              4.2,    -9.5376,
                              6.0,   -11.5376]
                CRAIRC:  0.090
                CNRatioSOMI: 9.0
                RHOD: 1.406
                Soil_pH: 7.4
                SoilID: SubSoil_10
        SoilProfileDescription:
            PFWiltingPoint: 4.2
            PFFieldCapacity: 2.0
            SurfaceConductivity: 70.0 # surface conductivity cm / day
            SoilLayers:
            -   <<: *TopSoil
                Thickness: 10
                FSOMI: 0.02
            -   <<: *TopSoil
                Thickness: 10
                FSOMI: 0.02
            -   <<: *MidSoil
                Thickness: 10
                FSOMI: 0.01
            -   <<: *MidSoil
                Thickness: 20
                FSOMI: 0.00
            -   <<: *MidSoil
                Thickness: 30
                FSOMI: 0.00
            -   <<: *SubSoil
                Thickness: 45
                FSOMI: 0.00
            SubSoilType:
                <<: *SubSoil
                Thickness: 200
            GroundWater: null
    """
    
    def __init__(self, parvalues):
        list.__init__(self)

        sp = DotMap(parvalues["SoilProfileDescription"])
        for layer_properties in sp.SoilLayers:
            layer = SoilLayer(layer_properties, sp.PFFieldCapacity, sp.PFWiltingPoint)
            self.append(layer)
        for attr, value in sp.items():
            if attr == "SoilLayers":
                continue
            if attr == "SubSoilType":
                value = SoilLayer(value, sp.PFFieldCapacity, sp.PFWiltingPoint)
            setattr(self, attr, value)

    def determine_rooting_status(self, RD, RDM):
        """Determines the rooting status of the soil layers and update layer weights.

        Soil layers can be rooted, partially rooted, potentially rooted or never rooted.
        This is stored in layer.rooting_status

        Note that this routine expected that the maximum rooting depth coincides
        with a layer boundary.

        :param RD:the current rooting depth
        :param RDM: the maximum rooting depth
        """
        upper_layer_boundary = 0
        lower_layer_boundary = 0
        for layer in self:
            lower_layer_boundary += layer.Thickness
            if lower_layer_boundary <= RD:
                layer.rooting_status = "rooted"
            elif upper_layer_boundary < RD < lower_layer_boundary:
                layer.rooting_status = "partially rooted"
            elif RD < lower_layer_boundary <= RDM:
                layer.rooting_status = "potentially rooted"
            else:
                layer.rooting_status = "never rooted"
            upper_layer_boundary = lower_layer_boundary

        self._compute_layer_weights(RD)

    def _compute_layer_weights(self, RD):
        """computes the layer weights given the current rooting depth.

        :param RD: The current rooting depth
        """
        lower_layer_boundary = 0
        for layer in self:
            lower_layer_boundary += layer.Thickness
            if layer.rooting_status == "rooted":
                layer.Wtop = 1.0
                layer.Wpot = 0.0
                layer.Wund = 0.0
            elif layer.rooting_status == "partially rooted":
                layer.Wtop = 1.0 - (lower_layer_boundary - RD) / layer.Thickness
                layer.Wpot = 1.0 - layer.Wtop
                layer.Wund = 0.0
            elif layer.rooting_status == "potentially rooted":
                layer.Wtop = 0.0
                layer.Wpot = 1.0
                layer.Wund = 0.0
            elif layer.rooting_status == "never rooted":
                layer.Wtop = 0.0
                layer.Wpot = 0.0
                layer.Wund = 1.0
            else:
                msg = "Unknown rooting status: %s" % layer.rooting_status
                raise exc.PCSEError(msg)

    def validate_max_rooting_depth(self, RDM):
        """Validate that the maximum rooting depth coincides with a layer boundary.

        :param RDM: The maximum rootable depth
        :return: True or False
        """
        tiny = 0.01
        lower_layer_boundary = 0
        for layer in self:
            lower_layer_boundary += layer.Thickness
            if abs(RDM - lower_layer_boundary) < tiny:
                break
        else:  # no break
            msg = "Current maximum rooting depth (%f) does not coincide with a layer boundary!" % RDM
            raise exc.PCSEError(msg)

    def get_max_rootable_depth(self):
        """Returns the maximum soil rootable depth.

        here we assume that the max rootable depth is equal to the lower boundary of the last layer.

        :return: the max rootable depth in cm
        """
        LayerThickness = [l.Thickness for l in self]
        LayerLowerBoundary = list(np.cumsum(LayerThickness))
        return max(LayerLowerBoundary)
# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Herman Berghuijs (herman.berghuijs@wur.nl) and Allard de Wit (allard.dewit@wur.nl), January 2024

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

class SoilNLayer(HasTraits):
    """
    Contains the intrinsic and derived properties for each soil layer that are used by
    SNOMIN. The following properties are read from the *.soil input file or derived from
    variables that are read from this file and defined for each layer:

    =============== ================================================ =======================
    Name             Description                                     Unit
    =============== ================================================ =======================
    CNRatioSOMI     Initial C:N ratio of soil organic matter         kg C kg-1 N
    CONDfromPF      Table function of the 10-base logarithm of the 
                    unsaturated hydraulic conductivity as a function 
                    of pF                                            log10(cm water d-1), -
    FSOMI           Initial fraction of soil organic matter in soil  kg OM kg-1 soil
    PFFieldCapacity pF value for which the soil moisture content is 
                    at field capacity                                -
    PFWiltingPoint  pF value for which the soil moisture content is 
                    at its permanent wilting point                   -
    RHOD            Bulk density of the soil                         g soil cm-3 soil
    SMfromPF        Table function that describes the soil moisture 
                    content as a function of pF                      m3 water m-3 soil, -
    SMsat           Soil moisture content at saturation              m3 water m-3 soil
    Soil_pH         pH of the soil layer                             -
    Thickness       Layer thickness                                  cm
    =============== ================================================ =======================
    """

    SMfromPF = Instance(pFCurve)  # soil moisture content as function of pF
    PFfromSM = Instance(Afgen)  # Inverse from SMfromPF
    SMsat = Float()  # Saturation soil moisture content
    Thickness = Float()
    rooting_status = Enum(["rooted","partially rooted","potentially rooted","never rooted",])

    # Parameters soil N model

    def __init__(self, layer, PFFieldCapacity, PFWiltingPoint):
        self.SMfromPF = pFCurve(layer.SMfromPF)
        self.PFfromSM = self._invert_pF(layer.SMfromPF)

        if 5 <= layer.Thickness <= 250:
            self.Thickness = layer.Thickness
        else:
            msg = "Soil layer should have thickness between 5 and 250 cm. Current value: %f" % layer.Thickness
            raise exc.PCSEError(msg)

        self.SM0 = self.SMfromPF(-1.0)
        self.SMFCF = self.SMfromPF(PFFieldCapacity)
        self.SMW = self.SMfromPF(PFWiltingPoint)
        self.rooting_status = None
        self.Soil_pH = layer.Soil_pH

        self.CNRatioSOMI = layer.CNRatioSOMI
        self.FSOMI = layer.FSOMI
        self.RHOD = layer.RHOD

        # compute hash value of this layer
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

class SoilNProfile(list):
    
    def __init__(self, parvalues):
        list.__init__(self)

        sp = DotMap(parvalues["SoilProfileDescription"])
        for layer_properties in sp.SoilLayers:
            layer = SoilNLayer(layer_properties, sp.PFFieldCapacity, sp.PFWiltingPoint)
            self.append(layer)
        for attr, value in sp.items():
            if attr == "SoilLayers":
                continue
            if attr == "SubSoilType":
                value = SoilNLayer(value, sp.PFFieldCapacity, sp.PFWiltingPoint)
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
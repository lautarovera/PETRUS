#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Corrections.py:
# This is the Corrections Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           Corrections.py
#  Date(YY/MM/DD): 16/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
from math import sin, sqrt, exp
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Iono import computeIonoMappingFunction
from InputOutput import RcvrIdx, SatIdx, LosIdx
import numpy as np


def runCorrectMeas(Conf, Rcvr, PreproObsInfo, SatInfo, LosInfo):

    # Purpose: correct GNSS preprocessed measurements and compute
    #          pseudo range residuals

    #          More in detail, this function handles the following:
    #          tasks:

    #             *  Correct the satellite navigation position and clock using EGNOS Fast-Long-Term (FLT) corrections: FC and LTC.
    #             *  Estimate the Slant Ionospheric Delay (UISD) using MOPS guidelines interpolation criteria for IGP Selection
    #             *  Estimate the Slant Troposphere delay (STD) using MOPS model (ZTD) and its mapping function. 
    #             *  Correct the Pre-processed measurements from Geometrical Range, Satellite clock, ionosphere and troposphere. 
    #             *  Build the Corrected Measurements and Measurement Residuals
    #             *  Estimate all Range level Sigmas to build the Sigma UERE:
    #                   -  Estimate the SigmaUIRE from MT26 information
    #                   -  Estimate the SigmaFLT from UDRE and MT28 
    #                   -  Estimate the SigmaTRO budget in line with MOPS.
    #                   -  Estimate the SigmaAirborne budget in line with MOPS 

    #             *  Estimate the Sigma UERE budget in line with MOPS


    # Parameters
    # ==========
    # Conf: dict
    #         Configuration dictionary
    # Rcvr: list
    #         Receiver information: position, masking angle...
    # PreproObsInfo: dict
    #         Preprocessed observations for current epoch per sat
    #         PreproObsInfo["G01"]["C1"]
    # SatInfo: dict
    #         dictionary containing the split lines of the SAT file
    #         SatInfo["G01"][1] is the second field of the line
    #         containing G01 info
    # LosInfo: dict
    #         dictionary containing the split lines of the LOS file
    #         SatInfo["G01"][1] is the second field of the line
    #         containing G01 info

    # Returns
    # =======
    # CorrInfo: dict
    #         Corrected measurements for current epoch per sat
    #         CorrInfo["G01"]["CorrectedPsr"]

    # Initialize output
    CorrInfo = OrderedDict({})

    # Initialize some values
    ResSum = 0.0
    ResN = 0
    EntGpsSum = 0.0
    EntGpsN = 0

    # Loop over satellites
    for prn, SatPrepro in PreproObsInfo.items():
        # If satellite is in convergence
        if(SatPrepro["Status"] == 1):
            # Initialize output info
            SatCorrInfo = {
                "Sod": 0.0,             # Second of day
                "Doy": 0,               # Day of year
                "Elevation": 0.0,       # Elevation
                "Azimuth": 0.0,         # Azimuth
                "IppLon": 0.0,          # IPP Longitude
                "IppLat": 0.0,          # IPP Latitude
                "Flag": 0,              # 0: Not Used 1: Used for PA 2: Used for NPA
                "SatX": 0.0,            # X-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatY": 0.0,            # Y-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatZ": 0.0,            # Z-Component of the Satellite Position 
                                        # corrected with SBAS LTC
                "SatClk": 0.0,          # Satellite Clock corrected with SBAS FLT
                "Uisd": 0.0,            # User Ionospheric Slant Delay
                "Std": 0.0,             # Slant Tropospheric Delay
                "CorrPsr": 0.0,         # Pseudo Range corrected from delays
                "GeomRange": 0.0,       # Geometrical Range (distance between Satellite 
                                        # Position and Receiver Reference Position)
                "PsrResidual": 0.0,     # Pseudo Range Residual
                "RcvrClk": 0.0,         # Receiver Clock estimation
                "SigmaFlt": 0,          # Sigma of the residual error associated to the 
                                        # fast and long-term correction (FLT)
                "SigmaUire": 0,         # User Ionospheric Range Error Sigma
                "SigmaTropo": 0,        # Sigma of the Tropospheric error 
                "SigmaAirborne": 0.0,   # Sigma Airborne Error
                "SigmaNoiseDiv": 0.0,   # Sigma of the receiver noise + divergence
                "SigmaMultipath": 0.0,  # Sigma of the receiver multipath
                "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of 
                                        # the total residual error associated to the 
                                        # satellite)
                "EntGps": 0.0,          # ENT to GPS Offset

            } # End of SatCorrInfo

            # Prepare outputs
            # Get SoD
            SatCorrInfo["Sod"] = SatPrepro["Sod"]
            # Get DoY
            SatCorrInfo["Doy"] = SatPrepro["Doy"]
            # Get Elevation
            SatCorrInfo["Elevation"] = SatPrepro["Elevation"]
            # Get Azimuth
            SatCorrInfo["Azimuth"] = SatPrepro["Azimuth"]

            # If SBAS information is available for current satellite
            if (prn in SatInfo) and (prn in LosInfo):
                # Get IPP Longitude
                SatCorrInfo["IppLon"] = float(LosInfo[prn][LosIdx["IPPLON"]])
                # Get IPP Latitude
                SatCorrInfo["IppLat"] = float(LosInfo[prn][LosIdx["IPPLAT"]])
                # Get STD
                SatCorrInfo["Std"] = float(LosInfo[prn][LosIdx["STD"]])

                # If satellite is monitored
                if int(SatInfo[prn][SatIdx["UDREI"]]) < 14:
                    if int(SatInfo[prn][SatIdx["UDREI"]]) < 12:
                        # Satellite Monitored and Usable for PA
                        SatCorrInfo["Flag"] = 1                         
                    else:
                        # Satellite Monitored and Usable for NPA
                        SatCorrInfo["Flag"] = 2
                # Satellite Non-Monitored (UDREI = 14) or Don't Use (UDREI = 15)
                else:          
                    SatCorrInfo["Flag"] = 0                           

                # If the satellite can be used in PA service
                if SatCorrInfo["Flag"] == 1:                                     
                    # -------------------
                    # PETRUS-CORR-REQ-010
                    # -------------------

                    # Correct Satellite X Position
                    SatCorrInfo["SatX"] = float(SatInfo[prn][SatIdx["SAT-X"]]) + float(SatInfo[prn][SatIdx["LTC-X"]])
                    # Correct Satellite Y Position
                    SatCorrInfo["SatY"] = float(SatInfo[prn][SatIdx["SAT-Y"]]) + float(SatInfo[prn][SatIdx["LTC-Y"]])
                    # Correct Satellite Z Position
                    SatCorrInfo["SatZ"] = float(SatInfo[prn][SatIdx["SAT-Z"]]) + float(SatInfo[prn][SatIdx["LTC-Z"]])
                    # Calculate DTR
                    SatOrb = np.array([float(SatInfo[prn][SatIdx["SAT-X"]]), float(SatInfo[prn][SatIdx["SAT-Y"]]), float(SatInfo[prn][SatIdx["SAT-Z"]])])
                    SatVel = np.array([float(SatInfo[prn][SatIdx["VEL-X"]]), float(SatInfo[prn][SatIdx["VEL-Y"]]), float(SatInfo[prn][SatIdx["VEL-Z"]])])
                    DTR = - 2 * (np.dot(SatOrb, SatVel)/Const.SPEED_OF_LIGHT)
                    # Correct Satellite CLK
                    SatCorrInfo["SatClk"] = float(SatInfo[prn][SatIdx["SAT-CLK"]]) + DTR \
                                                - float(SatInfo[prn][SatIdx["TGD"]])     \
                                                + float(SatInfo[prn][SatIdx["FC"]])      \
                                                + float(SatInfo[prn][SatIdx["LTC-B"]])

                    # -------------------
                    # PETRUS-CORR-REQ-030
                    # -------------------

                    # Error bound estimation
                    SatCorrInfo["SigmaFlt"] = ((float(SatInfo[prn][SatIdx["SIGMAUDRE"]]) * float(SatInfo[prn][SatIdx["DELTAUDRE"]])) \
                                                    + float(SatInfo[prn][SatIdx["EPS-FC"]]) + float(SatInfo[prn][SatIdx["EPS-RRC"]]) \
                                                    + float(SatInfo[prn][SatIdx["EPS-LTC"]]) + float(SatInfo[prn][SatIdx["EPS-ER"]])) ** 2 \
                                                    if float(SatInfo[prn][SatIdx["RSS"]]) == 0.00 \
                                                    else (float(SatInfo[prn][SatIdx["SIGMAUDRE"]]) * float(SatInfo[prn][SatIdx["DELTAUDRE"]])) ** 2 \
                                                    + float(SatInfo[prn][SatIdx["EPS-FC"]]) ** 2 + float(SatInfo[prn][SatIdx["EPS-RRC"]]) ** 2 \
                                                    + float(SatInfo[prn][SatIdx["EPS-LTC"]]) ** 2 + float(SatInfo[prn][SatIdx["EPS-ER"]]) ** 2
                
                    # -----------------------------------------------------------------------
                    # Computation of UISD and Sigma UIRE to comply: 
                    #                                                   - PETRUS-CORR-REQ-050
                    #                                                   - PETRUS-CORR-REQ-060
                    #                                                   - PETRUS-CORR-REQ-070
                    #                                                   - PETRUS-CORR-REQ-080 
                    #                  NORTH
                    #           v2..............v1
                    #           .                .
                    #           .                .
                    #           .                .
                    #   WEST    .       IPP      .      EAST
                    #           .        *       .
                    #           .                .
                    #           .                .
                    #           v3..............v4
                    #                  SOUTH
                    # -----------------------------------------------------------------------

                    LonV1 = float(LosInfo[prn][LosIdx["IGP_NE_LON"]])
                    LatV1 = float(LosInfo[prn][LosIdx["IGP_NE_LAT"]])
                    TauV1 = float(LosInfo[prn][LosIdx["GIVD_NE"]])
                    ErrTauV1 = float(LosInfo[prn][LosIdx["GIVE_NE"]])

                    LonV2 = float(LosInfo[prn][LosIdx["IGP_NW_LON"]])
                    LatV2 = float(LosInfo[prn][LosIdx["IGP_NW_LAT"]])
                    TauV2 = float(LosInfo[prn][LosIdx["GIVD_NW"]])
                    ErrTauV2 = float(LosInfo[prn][LosIdx["GIVE_NW"]])

                    LonV3 = float(LosInfo[prn][LosIdx["IGP_SW_LON"]])
                    LatV3 = float(LosInfo[prn][LosIdx["IGP_SW_LAT"]])
                    TauV3 = float(LosInfo[prn][LosIdx["GIVD_SW"]])
                    ErrTauV3 = float(LosInfo[prn][LosIdx["GIVE_SW"]])

                    LonV4 = float(LosInfo[prn][LosIdx["IGP_SE_LON"]])
                    LatV4 = float(LosInfo[prn][LosIdx["IGP_SE_LAT"]])
                    TauV4 = float(LosInfo[prn][LosIdx["GIVD_SE"]])
                    ErrTauV4 = float(LosInfo[prn][LosIdx["GIVE_SE"]])

                    LonWest = LonV2
                    LonEast = LonV4

                    LatNorth = LatV1
                    LatSouth = LatV3

                    # Square interpolation
                    if int(LosInfo[prn][LosIdx["INTERP"]]) == 0:
                        if abs(SatCorrInfo["IppLat"]) < 85: 
                            Xpp = (SatCorrInfo["IppLon"] - LonWest) / (LonEast - LonWest)
                            Ypp = (SatCorrInfo["IppLat"] - LatSouth) / (LatNorth - LatSouth)
                        else:
                            Ypp = (abs(SatCorrInfo["IppLat"]) - 85.0) / 10.0
                            Xpp = ((SatCorrInfo["IppLon"] - LonV3) / 90.0) * (1 - 2 * Ypp) + Ypp

                        W1 = Xpp * Ypp
                        W2 = (1 - Xpp) * Ypp
                        W3 = (1 - Xpp) * (1 - Ypp)
                        W4 = Xpp * (1 - Ypp)

                    # Triangular interpolation, v2 is always the vertex opposite the hypotenuse
                    else:
                        if abs(SatCorrInfo["IppLat"]) < 75: 
                            Xpp = (SatCorrInfo["IppLon"] - LonWest) / (LonEast - LonWest)
                            Ypp = (SatCorrInfo["IppLat"] - LatSouth) / (LatNorth - LatSouth)
                        else:
                            Ypp = (abs(SatCorrInfo["IppLat"]) - 85.0) / 10.0
                            Xpp = ((SatCorrInfo["IppLon"] - LonV3) / 90.0) * (1 - 2 * Ypp) + Ypp

                        W1 = Ypp
                        W2 = 1 - Xpp - Ypp
                        W3 = Xpp
                        W4 = 0

                        # Triangular interpolation without v1
                        if int(LosInfo[prn][LosIdx["INTERP"]]) == 1:
                            TauV1 = TauV2
                            ErrTauV1 = ErrTauV2

                            TauV2 = TauV3
                            ErrTauV2 = ErrTauV3

                            TauV3 = TauV4
                            ErrTauV3 = ErrTauV4

                        # Triangular interpolation without v2
                        elif int(LosInfo[prn][LosIdx["INTERP"]]) == 2:
                            TauV1 = TauV3
                            ErrTauV1 = ErrTauV3

                            TauV2 = TauV4
                            ErrTauV2 = ErrTauV4

                            TauV3 = TauV1
                            ErrTauV3 = ErrTauV1

                        # Triangular interpolation without v3
                        elif int(LosInfo[prn][LosIdx["INTERP"]]) == 3:
                            TauV1 = TauV4
                            ErrTauV1 = ErrTauV4

                            TauV2 = TauV1
                            ErrTauV2 = ErrTauV1

                            TauV3 = TauV2
                            ErrTauV3 = ErrTauV2

                        # Triangular interpolation without v4
                        elif int(LosInfo[prn][LosIdx["INTERP"]]) == 4:
                            pass

                        # No interpolate
                        else:
                            pass
                    
                    # Compute UIVD
                    TauVpp = (W1 * TauV1) + (W2 * TauV2) + (W3 * TauV3) + (W4 * TauV4)
                    # Compute UIVE
                    ErrTauVpp = sqrt((W1 * (ErrTauV1 ** 2)) + (W2 * (ErrTauV2 ** 2)) + (W3 * (ErrTauV3 ** 2)) + (W4 * (ErrTauV4 ** 2)))
                    # Compute mapping function
                    Fpp = computeIonoMappingFunction(SatCorrInfo["Elevation"])
                    # Compute UISD
                    SatCorrInfo["Uisd"] = round(Fpp * TauVpp, 4)
                    # Compute Sigma UIRE
                    SatCorrInfo["SigmaUire"] = round(Fpp * ErrTauVpp, 4)

                    # -------------------
                    # PETRUS-CORR-REQ-130
                    # -------------------

                    SigmaTVE = 0.12
                    Mpp = 1.001 / sqrt(0.002001 + sin(SatCorrInfo["Elevation"]) ** 2)
                    # Compute Sigma Tropo
                    SatCorrInfo["SigmaTropo"] = SigmaTVE * Mpp

                    # -------------------
                    # PETRUS-CORR-REQ-110
                    # -------------------

                    # Compute Sigma Multi Path
                    SatCorrInfo["SigmaMultipath"] = 0.13 + 0.53 * exp(-SatCorrInfo["Elevation"] / 10.0) if SatCorrInfo["Elevation"] > 2 else "ERROR"
                    # Compute Sigma Noise
                    SatCorrInfo["SigmaNoiseDiv"] = 0.36 if SatCorrInfo["Elevation"] < Conf["ELEV_NOISE_TH"] else 0.15
                    # Compute Sigma Airborne
                    SatCorrInfo["SigmaAirborne"] = sqrt(SatCorrInfo["SigmaMultipath"] ** 2 + SatCorrInfo["SigmaNoiseDiv"] ** 2)

                    # -------------------
                    # PETRUS-CORR-REQ-140
                    # -------------------
                    SatCorrInfo["SigmaUere"] = sqrt(SatCorrInfo["SigmaFlt"] ** 2 \
                                                + SatCorrInfo["SigmaUire"] ** 2 \
                                                + SatCorrInfo["SigmaTropo"] ** 2 \
                                                + SatCorrInfo["SigmaAirborne"] ** 2)

                    # Corrected Measurements from previous information
                    SatCorrInfo["CorrPsr"] = SatPrepro["SmoothC1"] + SatCorrInfo["SatClk"] - SatCorrInfo["Uisd"] - SatCorrInfo["Std"]
                    # Compute Geometrical Range
                    SATxyz = np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])
                    RCVRxyz = np.array(Rcvr[RcvrIdx["XYZ"]])
                    GEOMRANGExyz = np.subtract(SATxyz, RCVRxyz)
                    SatCorrInfo["GeomRange"] = np.linalg.norm(GEOMRANGExyz)
                    # Compute Geometrical Range
                    SatCorrInfo["PsrResidual"] = SatCorrInfo["CorrPsr"] - SatCorrInfo["GeomRange"]

                    # Compute residuals sum
                    Weight = 1 / (SatCorrInfo["SigmaUere"] ** 2)
                    ResSum += Weight * SatCorrInfo["PsrResidual"]
                    ResN += Weight

                    # Compute ENT-GPS offset
                    LTCxyz = np.array([float(SatInfo[prn][SatIdx["LTC-X"]]), (float(SatInfo[prn][SatIdx["LTC-Y"]])), (float(SatInfo[prn][SatIdx["LTC-Z"]]))])
                    Ulos = GEOMRANGExyz / SatCorrInfo["GeomRange"]
                    EntGpsSum += np.dot(LTCxyz, Ulos) - (float(SatInfo[prn][SatIdx["FC"]]) + float(SatInfo[prn][SatIdx["LTC-B"]]))
                    EntGpsN += 1

            # Prepare output for the satellite
            CorrInfo[prn] = SatCorrInfo

        # End of if(SatPrepro["Status"] == 1):

    # End of for prn, SatPrepro in PreproObsInfo.items():

    # Loop over satellites
    for prn in CorrInfo:
        # If SBAS information is available for current satellite
        if (prn in SatInfo) and (prn in LosInfo):
            # Receiver clock first guess
            CorrInfo[prn]["RcvrClk"] = ResSum / ResN
            CorrInfo[prn]["PsrResidual"] -= CorrInfo[prn]["RcvrClk"]
            CorrInfo[prn]["EntGps"] = EntGpsSum / EntGpsN

    return CorrInfo

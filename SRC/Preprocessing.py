#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Preprocessing.py:
# This is the Preprocessing Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           Preprocessing.py
#  Date(YY/MM/DD): 01/02/21
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
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from InputOutput import RcvrIdx, ObsIdx, REJECTION_CAUSE
from InputOutput import FLAG, VALUE, TH, CSNEPOCHS
import numpy as np
from COMMON.Iono import computeIonoMappingFunction
from operator import itemgetter

# Preprocessing internal functions
#-----------------------------------------------------------------------
def getElevationCut(Conf, PreproObs):
    ElevationCut = 0.0

    SatsNumberToReject = len(PreproObs) - int(Conf["NCHANNELS_GPS"])
    
    if SatsNumberToReject > 0:
        SatsElevation = [Sat["Elevation"] for Sat in PreproObs]
        print(SatsElevation)
        SatsElevation = sorted(SatsElevation)
        ElevationCut = SatsElevation[SatsNumberToReject]
    
    return ElevationCut


def detectCycleSlips(PreproObs, PrevPreproObs, Threshold):
    # Compute the Third order difference (TOD)
    CsFlag = False
    # Current and previous phase measurements
    CP_n = PreproObs["L1"]
    CP_n_1 = PrevPreproObs["L1_n_1"]
    CP_n_2 = PrevPreproObs["L1_n_2"]
    CP_n_3 = PrevPreproObs["L1_n_3"]
    # Previous measurements instants
    t1 = PreproObs["Sod"] - PrevPreproObs["t_n_1"]
    t2 = PrevPreproObs["t_n_1"] - PrevPreproObs["t_n_2"]
    t3 = PrevPreproObs["t_n_2"] - PrevPreproObs["t_n_3"]
    # Residuals equation factors
    if PrevPreproObs["t_n_3"] > 0:
        R1 = float((t1 + t2) * (t1 + t2 + t3)) / (t2 * (t2 + t3))
        R2 = float(-t1 * (t1 + t2 + t3)) / (t2 * t3)
        R3 = float(t1 * (t1 + t2)) / ((t2 + t3) * t3)
        # Compute TOD residuals
        # TOD = L1 -3*L1_n_1 + 3*L1_n_2 - L1_n_3
        CsResidual = abs(CP_n - R1 * CP_n_1 - R2 * CP_n_2 - R3 * CP_n_3)
        # Compute CS Flag
        CsFlag = CsResidual > Threshold

    return CsFlag


def updateCsBuffer(PrevPreproObs):
    PrevPreproObs["CsBuff"][PrevPreproObs["CsIdx"]] = 1
    PrevPreproObs["CsIdx"] += 1
    PrevPreproObs["CsIdx"] %= len(PrevPreproObs["CsBuff"])


def computeRate(Measurement, PrevMeasurement, dT):
    Rate = 0

    if dT > 0:
        Rate = abs(int(Measurement) - int(PrevMeasurement)) / int(dT)

    return Rate


def computeStep(Rate, PrevRate, dT):
    Step = 0

    if dT > 0:
        Step = abs(int(Rate) - int(PrevRate)) / int(dT)

    return Step


def resetHatchFilter(PreproObs, PrevPreproObs):
    PrevPreproObs["PrevEpoch"] = PreproObs["Sod"]
    PrevPreproObs["GapCounter"] = 0
    PreproObs["SmoothC1"] = PreproObs["C1"]
    PrevPreproObs["PrevSmoothC1"] = PreproObs["C1"]
    PrevPreproObs["Ksmooth"] = 1
    PrevPreproObs["PrevL1"] = PreproObs["L1Meters"]
    PrevPreproObs["PrevPhaseRateL1"] = -9999.9
    PrevPreproObs["PrevRangeRateL1"] = -9999.9
    PrevPreproObs["L1_n_1"] = 0.0
    PrevPreproObs["L1_n_2"] = 0.0
    PrevPreproObs["L1_n_3"] = 0.0
    PrevPreproObs["t_n_1"] = 0.0
    PrevPreproObs["t_n_2"] = 0.0
    PrevPreproObs["t_n_3"] = 0.0
    PrevPreproObs["CsBuff"] = [0] * len(PrevPreproObs["CsBuff"])
    PreproObs["Status"] = 0
    PrevPreproObs["ResetHatchFilter"] = 0
    

def runPreProcMeas(Conf, Rcvr, ObsInfo, PrevPreproObsInfo):
    
    # Purpose: preprocess GNSS raw measurements from OBS file
    #          and generate PREPRO OBS file with the cleaned,
    #          smoothed measurements

    #          More in detail, this function handles:
             
    #          * Measurements cleaning and validation and exclusion due to different 
    #          criteria as follows:
    #             - Minimum Masking angle
    #             - Maximum Number of channels
    #             - Minimum Carrier-To-Noise Ratio (CN0)
    #             - Pseudo-Range Output of Range 
    #             - Maximum Pseudo-Range Step
    #             - Maximum Pseudo-Range Rate
    #             - Maximum Carrier Phase Increase
    #             - Maximum Carrier Phase Increase Rate
    #             - Data Gaps checks and handling 
    #             - Cycle Slips detection

    #         * Filtering/Smoothing of Code-Phase Measurements with a Hatch filter 

    # Parameters
    # ==========
    # Conf: dict
    #         Configuration dictionary
    # Rcvr: list
    #         Receiver information: position, masking angle...
    # ObsInfo: list
    #         OBS info for current epoch
    #         ObsInfo[1][1] is the second field of the 
    #         second satellite
    # PrevPreproObsInfo: dict
    #         Preprocessed observations for previous epoch per sat
    #         PrevPreproObsInfo["G01"]["C1"]

    # Returns
    # =======
    # PreproObsInfo: dict
    #         Preprocessed observations for current epoch per sat
    #         PreproObsInfo["G01"]["C1"]
    
    # Initialize output
    PreproObsInfo = OrderedDict({})
    
    # Loop over satellites
    for SatObs in ObsInfo:
        # Initialize output info
        SatPreproObsInfo = {
            "Sod": 0.0,             # Second of day
            "Doy": 0,               # Day of year
            "Elevation": 0.0,       # Elevation
            "Azimuth": 0.0,         # Azimuth
            "C1": 0.0,              # GPS L1C/A pseudorange
            "P1": 0.0,              # GPS L1P pseudorange
            "L1": 0.0,              # GPS L1 carrier phase (in cycles)
            "L1Meters": 0.0,        # GPS L1 carrier phase (in m)
            "S1": 0.0,              # GPS L1C/A C/No
            "P2": 0.0,              # GPS L2P pseudorange
            "L2": 0.0,              # GPS L2 carrier phase 
            "S2": 0.0,              # GPS L2 C/No
            "SmoothC1": 0.0,        # Smoothed L1CA 
            "GeomFree": 0.0,        # Geom-free in Phases
            "GeomFreePrev": 0.0,    # t-1 Geom-free in Phases
            "ValidL1": 1,           # L1 Measurement Status
            "RejectionCause": 0,    # Cause of rejection flag
            "StatusL2": 0,          # L2 Measurement Status
            "Status": 0,            # L1 Smoothing status
            "RangeRateL1": 0.0,     # L1 Code Rate
            "RangeRateStepL1": 0.0, # L1 Code Rate Step
            "PhaseRateL1": 0.0,     # L1 Phase Rate
            "PhaseRateStepL1": 0.0, # L1 Phase Rate Step
            "VtecRate": 0.0,        # VTEC Rate
            "iAATR": 0.0,           # Instantaneous AATR
            "Mpp": 0.0,             # Iono Mapping
            "NSats": 0

        } # End of SatPreproObsInfo

        # Get satellite label
        SatLabel = SatObs[ObsIdx["CONST"]] + "%02d" % int(SatObs[ObsIdx["PRN"]])
        # Prepare outputs
        # Get SoD
        SatPreproObsInfo["Sod"] = float(SatObs[ObsIdx["SOD"]])
        # Get DoY
        SatPreproObsInfo["Doy"] = int(SatObs[ObsIdx["DOY"]])
        # Get Elevation
        SatPreproObsInfo["Elevation"] = float(SatObs[ObsIdx["ELEV"]])
        # Get Azimuth
        SatPreproObsInfo["Azimuth"] = float(SatObs[ObsIdx["AZIM"]])
        # Get C1
        SatPreproObsInfo["C1"] = float(SatObs[ObsIdx["C1"]])
        # Get L1
        SatPreproObsInfo["L1"] = float(SatObs[ObsIdx["L1"]])
        # Get L1 in meters
        SatPreproObsInfo["L1Meters"] = float(SatObs[ObsIdx["L1"]]) * Const.GPS_L1_WAVE
        # Get S1
        SatPreproObsInfo["S1"] = float(SatObs[ObsIdx["S1"]])
        # Get L2
        SatPreproObsInfo["L2"] = float(SatObs[ObsIdx["L2"]])
        # Compute Iono mapping function
        SatPreproObsInfo["Mpp"] = computeIonoMappingFunction(SatPreproObsInfo["Elevation"])

        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

    # Limit the satellites to the Number of Channels
    # ---------------------------------------------------------------------------
    # [PETRUS-PPVE-REQ-010]

    # Get elevation cut if the number of satellites in view exceeds the maximum allowed channels
    ChannelElevationCut = getElevationCut(Conf, PreproObsInfo)

    # Loop over all GPS PRNs
    for SatLabel, PreproObs in PreproObsInfo.items():
        # Check if the satellite was rejected by number of channels
        # ------------------------------------------------------------------------
        if PreproObs["Elevation"] < ChannelElevationCut:
            PreproObs["ValidL1"] = 0
            PreproObs["RejectionCause"] = REJECTION_CAUSE["NCHANNELS_GPS"]
            PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["NCHANNELS_GPS"]

            continue

        # Check satellite Elevation angle in front of the minimum by configuration
        # ------------------------------------------------------------------------
        if PreproObs["Elevation"] < Rcvr[RcvrIdx["MASK"]]:
            PreproObs["ValidL1"] = 0
            PreproObs["RejectionCause"] = REJECTION_CAUSE["MASKANGLE"]
            PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["MASKANGLE"]

            continue

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated) 
        # ------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-020]

        if Conf["MIN_CNR"][0] == 1 and PreproObs["S1"] < float(Conf["MIN_CNR"][1]):
            PreproObs["ValidL1"] = 0
            PreproObs["RejectionCause"] = REJECTION_CAUSE["MIN_CNR"]
            PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["MIN_CNR"]

            continue
            # print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Noise Rejection')

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration 
        # ------------------------------------------------------------------------------
        if (Conf["MAX_PSR_OUTRNG"][0] == 1 and PreproObs["C1"] > Conf["MAX_PSR_OUTRNG"][1]):
            PreproObs["ValidL1"] = 0
            PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PSR_OUTRNG"]
            PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["MAX_PSR_OUTRNG"]

            continue
            
        
        # Check Measurement Data gaps 
        # ------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-090]

        DeltaT = PreproObs["Sod"] - PrevPreproObsInfo[SatLabel]["PrevEpoch"]

        if (DeltaT > Conf["SAMPLING_RATE"]):
            # Increment gap counter
            PrevPreproObsInfo[SatLabel]["GapCounter"] = DeltaT
            # Check if gap is longer than the allowed maximum
            if PrevPreproObsInfo[SatLabel]["GapCounter"] > Conf["HATCH_GAP_TH"]:
                PrevPreproObsInfo[SatLabel]["GapCounter"] = 0
                # Reset filter
                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                if PrevPreproObsInfo[SatLabel]["PrevRej"] != REJECTION_CAUSE["MASKANGLE"]:
                    # PreproObsInfo[SatLabel]["ValidL1"] = 0
                    PreproObs["RejectionCause"] = REJECTION_CAUSE["DATA_GAP"]
                    PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["DATA_GAP"]
                    print('[TESTING][runPreProcMeas]' + ' SoD ' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Hatch filter reset (gap=' + "%.2f" % DeltaT + ')')
        else:
            PrevPreproObsInfo[SatLabel]["GapCounter"] = 0

        # Check Cycle Slips, if activated, with Third Order Difference algorithm 
        # -------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-080]

        if (not PrevPreproObsInfo[SatLabel]["ResetHatchFilter"]) and (Conf["MIN_NCS_TH"][0] == 1):
            # Check the existence of a Cycle Slip with 3rd order difference
            # IMPORTANT: The Phase shall not be propagated with Not Valid measurement
            # If some of the 3 last measurements were Not Valid, the most recent
            # valid measurements shall be used to propagate the Phase
            CsFlag = detectCycleSlips(PreproObs, PrevPreproObsInfo[SatLabel], float(Conf["MIN_NCS_TH"][1]))

            # Cumulate the Number of Consecutive Cycle Slips
            PrevPreproObsInfo[SatLabel]["CsBuff"][PrevPreproObsInfo[SatLabel]["CsIdx"]] = CsFlag

            # If L1 measurement diverged from the propagated value
            if CsFlag == True:
                # Set measurement as Not Valid
                PreproObsInfo[SatLabel]["ValidL1"] = 0
                # print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' CS')

            # If three last consecutive flags were raised
            if np.sum(PrevPreproObsInfo[SatLabel]["CsBuff"]) == Conf["MIN_NCS_TH"][CSNEPOCHS]:
                PreproObs["RejectionCause"] = REJECTION_CAUSE["CYCLE_SLIP"]
                PrevPreproObsInfo[SatLabel]["PrevRej"] = REJECTION_CAUSE["CYCLE_SLIP"]
                # Reset filter
                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                # print('[TESTING][runPreProcMeas]' + ' SoD ' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Hatch filter reset (CS)')
            else:
                # Update CsBuffer index
                PrevPreproObsInfo[SatLabel]["CsIdx"] += 1
                PrevPreproObsInfo[SatLabel]["CsIdx"] %= Conf["MIN_NCS_TH"][CSNEPOCHS]

                continue
            
            # Update CsBuffer index
            PrevPreproObsInfo[SatLabel]["CsIdx"] += 1
            PrevPreproObsInfo[SatLabel]["CsIdx"] %= Conf["MIN_NCS_TH"][CSNEPOCHS]

            if np.sum(PrevPreproObsInfo[SatLabel]["CsBuff"]) == 0:
                PrevPreproObsInfo[SatLabel]["L1_n_3"] = PrevPreproObsInfo[SatLabel]["L1_n_2"]
                PrevPreproObsInfo[SatLabel]["L1_n_2"] = PrevPreproObsInfo[SatLabel]["L1_n_1"]
                PrevPreproObsInfo[SatLabel]["L1_n_1"] = PreproObs["L1"]
                PrevPreproObsInfo[SatLabel]["t_n_3"] = PrevPreproObsInfo[SatLabel]["t_n_2"]
                PrevPreproObsInfo[SatLabel]["t_n_2"] = PrevPreproObsInfo[SatLabel]["t_n_1"]
                PrevPreproObsInfo[SatLabel]["t_n_1"] = PreproObs["Sod"]

        # Hatch filter (re)initialization
        # -------------------------------------------------------------------------------
        if PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] == 1:
            # Reset filter
            resetHatchFilter(PreproObsInfo[SatLabel], PrevPreproObsInfo[SatLabel])

            continue

        # Perform the Code Carrier Smoothing with a Hatch Filter of 100 seconds
        # -------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-100]

        # Count the number of seconds from Smoothing Start
        PrevPreproObsInfo[SatLabel]["Ksmooth"] = PrevPreproObsInfo[SatLabel]["Ksmooth"] + DeltaT

        # Update the Smoothing time if below the 100 seconds
        SmoothingTime = PrevPreproObsInfo[SatLabel]["Ksmooth"]

        # Update the Smoothing time if above the 100 seconds
        if PrevPreproObsInfo[SatLabel]["Ksmooth"] >= Conf["HATCH_TIME"]:
            SmoothingTime = Conf["HATCH_TIME"]

        # Weighting factor of the smoothing filter
        Alpha = float(DeltaT) / SmoothingTime
        # Compute the new Smoothed Code
        PreproObs["SmoothC1"] = (Alpha * PreproObs["C1"]) + (1 - Alpha) * PrevPreproObsInfo[SatLabel]["PrevSmoothC1"] + \
                                (PreproObs["L1Meters"] - PrevPreproObsInfo[SatLabel]["PrevL1"])

        # Check Phase Rate (if activated)
        # ------------------------------------------------------------------------------- 
        # [PETRUS-PPVE-REQ-040]

        # Compute the Phase Rate in m/s
        PreproObs["PhaseRateL1"] = computeRate(PreproObs["L1Meters"], PrevPreproObsInfo[SatLabel]["PrevL1"], DeltaT)

        if Conf["MAX_PHASE_RATE"][0] == 1:
            # Check Phase Jump
            if PreproObsInfo[SatLabel]["PhaseRateL1"] > Conf["MAX_PHASE_RATE"][1]:
                # Reset smoothing filter and raise not valid flag
                PreproObsInfo["ValidL1"] = 0
                PreproObsInfo["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE"]
                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                continue
        
        # Check Phase Rate Step (if activated) 
        # -------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-050]

        # Compute the Phase Rate Step in m/s2
        if (PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] != -9999.99):
            PreproObs["PhaseRateStepL1"] = computeStep(PreproObs["PhaseRateL1"], PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"], DeltaT)

            if Conf["MAX_PHASE_RATE_STEP"][0] == 1:
                # Check Phase Rate Jump
                if PreproObs["PhaseRateStepL1"] > Conf["MAX_PHASE_RATE_STEP"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PreproObs["ValidL1"] = 0
                    PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE_STEP"]
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                    continue
        
        # Check Code Rate detector (if activated) 
        # -------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-070]
        
        # Compute the Code Rate in m/s as the first derivative of Smoothed Codes
        PreproObs["RangeRateL1"] = computeRate(PreproObsInfo[SatLabel]["SmoothC1"], PrevPreproObsInfo[SatLabel]["PrevSmoothC1"], DeltaT)

        if Conf["MAX_CODE_RATE"][0]:
            # Check Code Jump
            if abs(PreproObs["RangeRateL1"]) > Conf["MAX_CODE_RATE"][1]:
                # Reset smoothing filter and raise not valid flag
                PreproObs["ValidL1"] = 0
                PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE"]
                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                    
                continue

        # Check Code Rate Step detector (if activated)
        # -------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-060]

        if (PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] != -9999.99):
            # Compute the Code Rate step in m/s2 as the second derivative of Smoothed Codes
            PreproObs["RangeRateStepL1"] = computeStep(PreproObs["RangeRateL1"], PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"], DeltaT)

            if Conf["MAX_CODE_RATE_STEP"][0] == 1:
                # Check Code Rate Step
                if abs(PreproObs["RangeRateStepL1"]) > Conf["MAX_CODE_RATE_STEP"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PreproObs["ValidL1"] = 0
                    PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE_STEP"]
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                    continue

        # Update Measurement Smoothing Status and Hatch filter Convergence 
        # -------------------------------------------------------------------------------
        if (PrevPreproObsInfo[SatLabel]["Ksmooth"] > (Conf["HATCH_STATE_F"] * Conf["HATCH_TIME"])) and (PreproObs["ValidL1"] == 1):
            PreproObs["Status"] = 1
        else:
            PreproObs["Status"] = 0

        # Update previous values
        PrevPreproObsInfo[SatLabel]["PrevEpoch"] = PreproObs["Sod"]
        PrevPreproObsInfo[SatLabel]["PrevL1"] = PreproObs["L1"]
        PrevPreproObsInfo[SatLabel]["PrevSmoothC1"] = PreproObs["SmoothC1"]
        PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] = PreproObs["RangeRateL1"]
        PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] = PreproObs["PhaseRateL1"]
        PrevPreproObsInfo[SatLabel]["PrevRej"] = PreproObs["RejectionCause"]

        # Build Geometry-Free combination of Phases
        # Check if L1 and L2 are OK
        if PreproObs["ValidL1"] == 1 and PreproObs["L2"] > 0:
            # Compute the Geometry-Free Observable
            PreproObs["GeomFree"] = (Const.GPS_L1_WAVE * PreproObs["L1"]) - (Const.GPS_L2_WAVE * PreproObs["L2"])
            # Obtain the final Geom-free dividing by 1 - GAMMA
            PreproObs["GeomFree"] = PreproObs["GeomFree"] / (1 - Const.GPS_GAMMA_L1L2)

        # Compute the VTEC Rate
        if PrevPreproObsInfo[SatLabel]["PrevGeomFree"] > 0:
            # Compute STEC Gradient
            DeltaSTEC = computeRate(PreproObs["GeomFree"], PrevPreproObsInfo[SatLabel]["PrevGeomFree"], 
                                    PreproObs["Sod"] - PrevPreproObsInfo[SatLabel]["PrevGeomFreeEpoch"])
            # Compute VTEC Gradient
            DeltaVTEC = DeltaSTEC / PreproObs["Mpp"]
            # Store Delta VTEC in mm/s
            PreproObs["VtecRate"] = DeltaVTEC * 1000
            # Compute AATR
            PreproObs["iAATR"] = PreproObs["VtecRate"] / PreproObs["Mpp"]

        # Update previous geometry-free values
        PrevPreproObsInfo[SatLabel]["PrevGeomFreeEpoch"] = PreproObs["Sod"]
        PrevPreproObsInfo[SatLabel]["PrevGeomFree"] = PreproObs["GeomFree"]

    return PreproObsInfo

# End of function runPreProcMeas()

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################

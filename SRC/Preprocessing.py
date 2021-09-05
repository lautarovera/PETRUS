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

# Globals
#-----------------------------------------------------------------------
NSATS_GPS = NSATS_GPS = 0

# Preprocessing internal functions
#-----------------------------------------------------------------------
def rejectGpsSatsMinElevation(Conf, Obs, PreproObs):
    SatsElevation = {Sat[ObsIdx["PRN"]]:float(Sat[ObsIdx["ELEV"]]) for Sat in Obs if Sat[ObsIdx["CONST"]] == "G"}
    SatsNumberToReject = NSATS_GPS - Conf["NCHANNELS_GPS"]
    SatsRejected = dict(sorted(SatsElevation.items(), key = itemgetter(1))[:SatsNumberToReject])

    for Sat in SatsRejected:
        SatLabel = "G" + "%02d" % int(Sat)
        PreproObs[SatLabel]["ValidL1"] = 0
        PreproObs[SatLabel]["RejectionCause"] = REJECTION_CAUSE["NCHANNELS_GPS"]

def rejectGalSatsMinElevation(Conf, Obs, PreproObs):
    SatsElevation = {Sat[ObsIdx["PRN"]]:float(Sat[ObsIdx["ELEV"]]) for Sat in Obs if Sat[ObsIdx["CONST"]] == "E"}
    SatsNumberToReject = NSATS_GPS - Conf["NCHANNELS_GAL"]
    SatsRejected = dict(sorted(SatsElevation.items(), key = itemgetter(1))[:SatsNumberToReject])

    for Sat in SatsRejected:
        SatLabel = "E" + "%02d" % int(Sat)
        PreproObs[SatLabel]["ValidL1"] = 0
        PreproObs[SatLabel]["RejectionCause"] = REJECTION_CAUSE["NCHANNELS_GPS"]

def detectCycleSlips(PreproObs, PrevPreproObs, Prn, Threshold):
    # Compute the Third order difference (TOD)
    result = False
    SatLabel = "G" + "%02d" % int(Prn)
    PrevPreproObs[SatLabel]["CsBuff"][PrevPreproObs[SatLabel]["CsIdx"]] = 0
    # Discard detection if not enough measurements
    if not PrevPreproObs[SatLabel]["t_n_3"] > 0:
        return result
    # Current and previous phase measurements
    CP_n = PreproObs[SatLabel]["L1"]
    CP_n_1 = PrevPreproObs[SatLabel]["L1_n_1"]
    CP_n_2 = PrevPreproObs[SatLabel]["L1_n_2"]
    CP_n_3 = PrevPreproObs[SatLabel]["L1_n_3"]
    # Previous measurements instants
    t1 = PreproObs[SatLabel]["Sod"] - PrevPreproObs[SatLabel]["t_n_1"]
    t2 = PrevPreproObs[SatLabel]["t_n_1"] - PrevPreproObs[SatLabel]["t_n_2"]
    t3 = PrevPreproObs[SatLabel]["t_n_2"] - PrevPreproObs[SatLabel]["t_n_3"]
    # Residuals equation factors
    R1 = float((t1 + t2) * (t1 + t2 + t3)) / (t2 * (t2 + t3))
    R2 = float(-t1 * (t1 + t2 + t3)) / (t2 * t3)
    R3 = float(t1 * (t1 + t2)) / ((t2 + t3) * t3)
    # Compute TOD residuals
    # TOD = L1 -3*L1_n_1 + 3*L1_n_2 - L1_n_3
    CsResidual = abs(CP_n - R1 * CP_n_1 - R2 * CP_n_2 - R3 * CP_n_3)
    # Compute CS Flag
    if CsResidual > Threshold:
        # Cumulate the Number of Consecutive Cycle Slips
        PrevPreproObs[SatLabel]["CsBuff"][PrevPreproObs[SatLabel]["CsIdx"]] = 1
        result = True

    PrevPreproObs[SatLabel]["CsIdx"] += 1
    PrevPreproObs[SatLabel]["CsIdx"] %= len(PrevPreproObs[SatLabel]["CsBuff"])

    return result

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
    PrevPreproObs["Ksmooth"] = 0
    PrevPreproObs["PrevSmoothC1"] = 0
    PrevPreproObs["PrevPhaseRateL1"] = 0
    PrevPreproObs["PrevRangeRateL1"] = 0
    PrevPreproObs["L1_n_3"] = 0
    PrevPreproObs["L1_n_2"] = 0
    PrevPreproObs["L1_n_1"] = PreproObs["L1"]
    PrevPreproObs["t_n_3"] = 0
    PrevPreproObs["t_n_2"] = 0
    PrevPreproObs["t_n_1"] = PreproObs["Sod"]

def updatePrevPrepro(PreproObs, PrevPreproObs):
    # Propagate the current valid measurements to previous epoch
    # -------------------------------------------------------------------------------
    PrevPreproObs["PrevRej"] = PreproObs["RejectionCause"]
    PrevPreproObs["PrevEpoch"] = PreproObs["Sod"]
    
    if PreproObs["ValidL1"] == 0 and PrevPreproObs["ResetHatchFilter"] == 0:
        PrevPreproObs["L1_n_3"] = PrevPreproObs["L1_n_2"]
        PrevPreproObs["L1_n_2"] = PrevPreproObs["L1_n_1"]
        PrevPreproObs["L1_n_1"] = PreproObs["L1"]
        PrevPreproObs["t_n_3"] = PrevPreproObs["t_n_2"]
        PrevPreproObs["t_n_2"] = PrevPreproObs["t_n_1"]
        PrevPreproObs["t_n_1"] = PreproObs["Sod"]
        PrevPreproObs["PrevL1"] = PreproObs["L1Meters"]
        PrevPreproObs["PrevSmoothC1"] = PreproObs["SmoothC1"]
        PrevPreproObs["PrevPhaseRateL1"] = PreproObs["PhaseRateL1"]
        PrevPreproObs["PrevRangeRateL1"] = PreproObs["RangeRateL1"]


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
    GapCounter = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    ResetHatchFilter = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    PhaseRate = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    PhaseRateStep = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    CodeRate = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    CodeRateStep = {int(Sat[ObsIdx["PRN"]]):0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}

    # Number of satellites visibles for each constellation
    NSATS_GPS = NSATS_GAL = 0

    # Loop over satellites
    for SatObs in ObsInfo:
        # Get GPS satellites in view 
        if SatObs[ObsIdx["CONST"]] == "G":
            NSATS_GPS += 1
        # Get GAL satellites in view 
        if SatObs[ObsIdx["CONST"]] == "E":
            NSATS_GAL += 1

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

        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

    # Limit the satellites to the Number of Channels
    # ---------------------------------------------------------------------------
    # [PETRUS-PPVE-REQ-010]

    # If the number of satellites in view exceeds the maximum allowed channels
    if NSATS_GPS > Conf["NCHANNELS_GPS"]:
        # Remove those satellites with the lower Elevation
        rejectGpsSatsMinElevation(Conf, ObsInfo, PreproObsInfo)

    if NSATS_GAL > Conf["NCHANNELS_GAL"]:
        # Remove those satellites with the lower Elevation
        rejectGalSatsMinElevation(Conf, ObsInfo, PreproObsInfo)

    # Loop over all GPS PRNs
    for prn in range(1, Const.MAX_NUM_SATS_CONSTEL + 1):
        SatLabel = "G" + "%02d" % int(prn)
        # Check if the satellite is in view
        # ------------------------------------------------------------------------
        if not SatLabel in PreproObsInfo:
            continue

        # Check if the satellite was rejected by number of channels
        # ------------------------------------------------------------------------
        if PreproObsInfo[SatLabel]["ValidL1"] == 0:
            continue

        # Check satellite Elevation angle in front of the minimum by configuration
        # ------------------------------------------------------------------------
        if PreproObsInfo[SatLabel]["Elevation"] < Conf["RCVR_MASK"]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MASKANGLE"]

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated) 
        # ------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-020]

        elif Conf["MIN_CNR"][0] == 1 and PreproObsInfo[SatLabel]["S1"] < Conf["MIN_CNR"][1]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MIN_CNR"]
            # print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Noise Rejection')

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration 
        # ------------------------------------------------------------------------------
        elif Conf["MAX_PSR_OUTRNG"][0] == 1 and PreproObsInfo[SatLabel]["C1"] > Conf["MAX_PSR_OUTRNG"][1]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MAX_PSR_OUTRNG"]
            
        else:
            # Check Measurement Data gaps 
            # ------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-090]

            DeltaT = int(PreproObsInfo[SatLabel]["Sod"] - PrevPreproObsInfo[SatLabel]["PrevEpoch"])

            if DeltaT > int(Conf["SAMPLING_RATE"]):
                # Increment gap counter
                GapCounter[prn] = DeltaT
                # Check if gap is longer than the allowed maximum
                if GapCounter[prn] > Conf["HATCH_GAP_TH"]:
                    # Reset filter
                    ResetHatchFilter[prn] = 1
                    if PrevPreproObsInfo[SatLabel]["PrevRej"] != REJECTION_CAUSE["MASKANGLE"]:
                        PreproObsInfo[SatLabel]["ValidL1"] = 0
                        PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["DATA_GAP"]
                        print('[TESTING][runPreProcMeas]' + ' SoD ' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Hatch filter reset (gap=' + "%.2f" % DeltaT + ')')

            # Check Cycle Slips, if activated, with Third Order Difference algorithm 
            # -------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-080]

            if (ResetHatchFilter[prn] == 0) and (Conf["MIN_NCS_TH"][0] == 1):
                # Check the existence of a Cycle Slip with 3rd order difference
                # IMPORTANT: The Phase shall not be propagated with Not Valid measurement
                # If some of the 3 last measurements were Not Valid, the most recent
                # valid measurements shall be used to propagate the Phase
                CsFlag = detectCycleSlips(PreproObsInfo, PrevPreproObsInfo, prn, Conf["MIN_NCS_TH"][1])

                # If L1 measurement diverged from the propagated value
                if CsFlag == True:
                    # Set measurement as Not Valid
                    PreproObsInfo[SatLabel]["ValidL1"] = 0
                    # updateCsBuffer(PrevPreproObsInfo[SatLabel])
                    # print('[TESTING][runPreProcMeas]' + ' epoch' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' CS')

                # If three last consecutive flags were raised
                if sum(PrevPreproObsInfo[SatLabel]["CsBuff"]) == Conf["MIN_NCS_TH"][2]:
                    # Reset filter
                    ResetHatchFilter[prn] = 1
                    PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["CYCLE_SLIP"]
                    print('[TESTING][runPreProcMeas]' + ' SoD ' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Hatch filter reset (CS)')
            
            # Reset hatch filter requested from previous epoch
            if PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] == 1:
                ResetHatchFilter[prn] = PrevPreproObsInfo[SatLabel]["ResetHatchFilter"]
                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 0

            # Hatch filter (re)initialization
            # -------------------------------------------------------------------------------
            if ResetHatchFilter[prn] == 1:
                # Reset filter
                resetHatchFilter(PreproObsInfo[SatLabel], PrevPreproObsInfo[SatLabel])

            # Perform the Code Carrier Smoothing with a Hatch Filter of 100 seconds
            # -------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-100]

            else:
                # Count the number of seconds from Smoothing Start
                PrevPreproObsInfo[SatLabel]["Ksmooth"] = int(PrevPreproObsInfo[SatLabel]["Ksmooth"]) + DeltaT

                # Update the Smoothing time if below the 100 seconds
                SmoothingTime = int(PrevPreproObsInfo[SatLabel]["Ksmooth"])

                # Update the Smoothing time if above the 100 seconds
                if PrevPreproObsInfo[SatLabel]["Ksmooth"] >= Conf["HATCH_TIME"]:
                    SmoothingTime = int(Conf["HATCH_TIME"])

                # Weighting factor of the smoothing filter
                Alpha = DeltaT / SmoothingTime
                # Compute the new Smoothed Code
                PreproObsInfo[SatLabel]["SmoothC1"] = Alpha * PreproObsInfo[SatLabel]["C1"] + (1 - Alpha) * PrevPreproObsInfo[SatLabel]["PrevSmoothC1"]

            # Check Phase Rate (if activated)
            # ------------------------------------------------------------------------------- 
            # [PETRUS-PPVE-REQ-040]

            # Compute the Phase Rate in m/s
            PreproObsInfo[SatLabel]["PhaseRateL1"] = computeRate(PreproObsInfo[SatLabel]["L1Meters"], PrevPreproObsInfo[SatLabel]["PrevL1"], DeltaT)

            if Conf["MAX_PHASE_RATE"][0] == 1 and PreproObsInfo[SatLabel]["ValidL1"] == 1 and ResetHatchFilter[prn] == 0:
                # Check Phase Jump
                if PhaseRate[prn] > Conf["MAX_PHASE_RATE"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                    PreproObsInfo[SatLabel]["ValidL1"] = 0
                    PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE"]
            
            # Check Phase Rate Step (if activated) 
            # -------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-050]

            # Compute the Phase Rate Step in m/s2
            PreproObsInfo[SatLabel]["PhaseRateStepL1"] = computeStep(PreproObsInfo[SatLabel]["PhaseRateL1"], PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"], DeltaT)

            if Conf["MAX_PHASE_RATE_STEP"][0] == 1 and PreproObsInfo[SatLabel]["ValidL1"] == 1 and ResetHatchFilter[prn] == 0:
                # Check Phase Rate Jump
                if PreproObsInfo[SatLabel]["PhaseRateStepL1"] > Conf["MAX_PHASE_RATE_STEP"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                    PreproObsInfo[SatLabel]["ValidL1"] = 0
                    PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE_STEP"]
            
            # Check Code Rate detector (if activated) 
            # -------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-070]
            
            # Compute the Code Rate in m/s as the first derivative of Smoothed Codes
            PreproObsInfo[SatLabel]["RangeRateL1"] = computeRate(PreproObsInfo[SatLabel]["SmoothC1"], PrevPreproObsInfo[SatLabel]["PrevSmoothC1"], DeltaT)

            if Conf["MAX_CODE_RATE"][0] == 1 and PreproObsInfo[SatLabel]["ValidL1"] == 1 and ResetHatchFilter[prn] == 0:
                # Check Code Jump
                if PreproObsInfo[SatLabel]["RangeRateL1"] > Conf["MAX_CODE_RATE"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                    PreproObsInfo[SatLabel]["ValidL1"] = 0
                    PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE"]

            # Check Code Rate Step detector (if activated)
            # -------------------------------------------------------------------------------
            # [PETRUS-PPVE-REQ-060]

            # Compute the Code Rate step in m/s2 as the second derivative of Smoothed Codes
            PreproObsInfo[SatLabel]["RangeRateStepL1"] = computeStep(PreproObsInfo[SatLabel]["RangeRateL1"], PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"], DeltaT)

            if Conf["MAX_CODE_RATE_STEP"][0] == 1 and PreproObsInfo[SatLabel]["ValidL1"] == 1 and ResetHatchFilter[prn] == 0:
                # Check Code Rate Step
                if PreproObsInfo[SatLabel]["RangeRateStepL1"] > Conf["MAX_CODE_RATE_STEP"][1]:
                    # Reset smoothing filter and raise not valid flag
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1
                    PreproObsInfo[SatLabel]["ValidL1"] = 0
                    PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE_STEP"]

        updatePrevPrepro(PreproObsInfo[SatLabel], PrevPreproObsInfo[SatLabel])

    return PreproObsInfo

# End of function runPreProcMeas()

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################

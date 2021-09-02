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

def updatePrevPrepro(PreproObs, PrevPreproObs, ResetHatchFilter):
    # Propagate the current valid measurements to previous epoch
    # -------------------------------------------------------------------------------
    PrevPreproObs["PrevRej"] = PreproObs["RejectionCause"]
    
    if PreproObs["ValidL1"] == 1:
        PrevPreproObs["PrevEpoch"] = PreproObs["Sod"]
        PrevPreproObs["L1_n_3"] = PrevPreproObs["L1_n_2"]
        PrevPreproObs["L1_n_2"] = PrevPreproObs["L1_n_1"]
        PrevPreproObs["L1_n_1"] = PreproObs["L1"]
        PrevPreproObs["t_n_3"] = PrevPreproObs["t_n_2"]
        PrevPreproObs["t_n_2"] = PrevPreproObs["t_n_1"]
        PrevPreproObs["t_n_1"] = PreproObs["Sod"]

    if ResetHatchFilter == 1:
        PrevPreproObs["PrevEpoch"] = PreproObs["Sod"]
        PrevPreproObs["L1_n_3"] = 0
        PrevPreproObs["L1_n_2"] = 0
        PrevPreproObs["L1_n_1"] = PreproObs["L1"]
        PrevPreproObs["t_n_3"] = 0
        PrevPreproObs["t_n_2"] = 0
        PrevPreproObs["t_n_1"] = PreproObs["Sod"]

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
            
        # Measurement quality monitoring 

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
            DeltaT = int(PreproObsInfo[SatLabel]["Sod"] - PrevPreproObsInfo[SatLabel]["PrevEpoch"])

            if DeltaT > int(Conf["SAMPLING_RATE"]):
                # Increment gap counter
                GapCounter[prn] = DeltaT

                # Check if gap is longer than the allowed maximum
                if GapCounter[prn] > Conf["HATCH_GAP_TH"]:
                    # Reset filter
                    ResetHatchFilter[prn] = 1

                    if PrevPreproObsInfo[SatLabel]["PrevRej"] != REJECTION_CAUSE["MASKANGLE"]:
                        PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["DATA_GAP"]
                        print('[TESTING][runPreProcMeas]' + ' SoD ' + ObsInfo[0][0] + ' Satellite ' + SatLabel + ' Hatch filter reset (gap=' + "%.2f" % DeltaT + ')')

            # Check Cycle Slips, if activated, with Third Order Difference algorithm 
            # -------------------------------------------------------------------------------
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

        updatePrevPrepro(PreproObsInfo[SatLabel], PrevPreproObsInfo[SatLabel], ResetHatchFilter[prn])

    return PreproObsInfo

# End of function runPreProcMeas()

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################

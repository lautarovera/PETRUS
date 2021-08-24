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
def rejectGpsSatsMinElevation(Conf, ObsInfo, PreproObsInfo):
    SatsElevation = {Sat[ObsIdx["PRN"]]:float(Sat[ObsIdx["ELEV"]]) for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    SatsNumberToReject = NSATS_GPS - Conf["NCHANNELS_GPS"]
    SatsRejected = dict(sorted(SatsElevation.items(), key = itemgetter(1))[:SatsNumberToReject])

    for Sat in SatsRejected:
        SatLabel = "G" + "%02d" % int(Sat)
        PreproObsInfo[SatLabel]["ValidL1"] = 0
        PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["NCHANNELS_GPS"]

def rejectGalSatsMinElevation(Conf, ObsInfo, PreproObsInfo):
    SatsElevation = {Sat[ObsIdx["PRN"]]:float(Sat[ObsIdx["ELEV"]]) for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "E"}
    SatsNumberToReject = NSATS_GPS - Conf["NCHANNELS_GAL"]
    SatsRejected = dict(sorted(SatsElevation.items(), key = itemgetter(1))[:SatsNumberToReject])

    for Sat in SatsRejected:
        SatLabel = "E" + "%02d" % int(Sat)
        PreproObsInfo[SatLabel]["ValidL1"] = 0
        PreproObsInfo[SatLabel]["RejectionCause"] = REJECTION_CAUSE["NCHANNELS_GPS"]


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
    GapCounter = {Sat[ObsIdx["PRN"]]:0 for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}
    ResetHatchFilter = {Sat[ObsIdx["PRN"]]:False for Sat in ObsInfo if Sat[ObsIdx["CONST"]] == "G"}

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
            "ValidL1": 1,          # L1 Measurement Status
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

        # Check if the satellite is valid
        # ------------------------------------------------------------------------
        if PreproObsInfo[SatLabel]["ValidL1"] == 0:
            continue

        # Check satellite Elevation angle in front of the minimum by configuration
        # ------------------------------------------------------------------------
        if PreproObsInfo[SatLabel]["Elevation"] < Conf["RCVR_MASK"]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = 2

        # Measurement quality monitoring 

        # Check Signal To Noise Ratio in front of Minimum by configuration (if activated) 
        # ------------------------------------------------------------------------------
        # [PETRUS-PPVE-REQ-020]
        
        if Conf["MIN_CNR"][0] == 1 and PreproObsInfo[SatLabel]["S1"] < Conf["MIN_CNR"][1]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = 3

        # Check Pseudo-ranges Out-of-Range in front of Maximum by configuration 
        # ------------------------------------------------------------------------------
        if Conf["MAX_PSR_OUTRNG"][0] == 1 and PreproObsInfo[SatLabel]["C1"] > Conf["MAX_PSR_OUTRNG"][1]:
            PreproObsInfo[SatLabel]["ValidL1"] = 0
            PreproObsInfo[SatLabel]["RejectionCause"] = 4

        # Check Measurement Data gaps 
        # ------------------------------------------------------------------------------
        DeltaT = PreproObsInfo[SatLabel]["Sod"] - PrevPreproObsInfo[SatLabel]["PrevEpoch"]

        if DeltaT > Conf["SAMPLING_RATE"]:
            # Increment gap counter
            GapCounter[prn] = DeltaT
            # Check if gap is longer than the allowed maximum
            if GapCounter[prn] > Conf["HATCH_GAP_TH"]:
                # Reset filter
                ResetHatchFilter[prn] = True

        PrevPreproObsInfo[SatLabel]["PrevEpoch"] = PreproObsInfo[SatLabel]["Sod"]

    return PreproObsInfo

# End of function runPreProcMeas()

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################

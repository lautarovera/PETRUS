#!/usr/bin/env python

########################################################################
# PETRUS/SRC/PreprocessingPlots.py:
# This is the PreprocessingPlots Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           PreprocessingPlots.py
#  Date(YY/MM/DD): 05/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import CorrIdx, SatIdx, RcvrIdx
from InputOutput import REJECTION_CAUSE_DESC
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot
from COMMON.Coordinates import xyz2llh
import numpy as np
from collections import OrderedDict

def initPlot(CorrFile, PlotConf, Title, Label):
    CorrFileName = os.path.basename(CorrFile)
    CorrFileNameSplit = CorrFileName.split('_')
    Rcvr = CorrFileNameSplit[1]
    DatepDat = CorrFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/CORR/figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)


# Plot Satellite Visibility
def plotSatTracks(CorrFile, CorrData, RcvrInfo):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Monitored Satellite Tracks", "MONSAT_TRACKS")
    
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (16.8,15.2)

    RcvrLat = int(RcvrInfo[RcvrIdx["LAT"]])
    RcvrLon = int(RcvrInfo[RcvrIdx["LON"]])

    if RcvrLat > 55:
        PlotConf["LonMin"] = -180
        PlotConf["LonMax"] = 180
        PlotConf["LatMin"] = (RcvrLat - 70)
        PlotConf["LatMax"] = (RcvrLat + 10) if RcvrLat < 80 else 90
    else:
        PlotConf["LonMin"] = (RcvrLon - 115) 
        PlotConf["LonMax"] = (RcvrLon + 115) 
        PlotConf["LatMin"] = (RcvrLat - 70)
        PlotConf["LatMax"] = (RcvrLat + 35) if RcvrLat < 55 else 90

    PlotConf["LonStep"] = 15
    PlotConf["LatStep"] = 10

    PlotConf["zLabel"] = "Elevation [deg]"

    PlotConf["yTicks"] = range(PlotConf["LatMin"], PlotConf["LatMax"] + 1, PlotConf["LatStep"])
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"], PlotConf["LonMax"] + 1, PlotConf["LonStep"])
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    # Transform ECEF to Geodetic
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    x = CorrData[CorrIdx["SAT-X"]][FilterCond].to_numpy()
    y = CorrData[CorrIdx["SAT-Y"]][FilterCond].to_numpy()
    z = CorrData[CorrIdx["SAT-Z"]][FilterCond].to_numpy()
    DataLen = len(x)
    Longitude = np.zeros(DataLen)
    Latitude = np.zeros(DataLen)

    for index in range(DataLen):
        Longitude[index], Latitude[index], h = xyz2llh(x[index], y[index], z[index])

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = Longitude
    PlotConf["yData"][Label] = Latitude
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Number of Satellites
def plotNumSats(PreproObsFile, PreproObsData):
    pass
# Plot Satellite Polar View
def plotSatPolarView(PreproObsFile, PreproObsData):
    pass

# Plot C1 - C1Smoothed
def plotC1C1Smoothed(PreproObsFile, PreproObsData):
    pass

# Plot Rejection Flags
def plotRejectionFlags(PreproObsFile, PreproObsData):
    pass

# Plot Code Rate
def plotCodeRate(PreproObsFile, PreproObsData):
    pass

# Plot Phase Rate
def plotPhaseRate(PreproObsFile, PreproObsData):
    pass

# VTEC Gradient
def plotVtecGradient(PreproObsFile, PreproObsData):
    pass

# AATR index
def plotAatr(PreproObsFile, PreproObsData):
    pass

def generateCorrPlots(CorrFile, SatFile, RcvrInfo):
    
    # Purpose: generate output plots regarding Corrections results

    # Parameters
    # ==========
    # CorrFile: str
    #         Path to CORRECTIONS output file
    # SatFile: str
    #         Path to SATELLITE output file
    # RcvrInfo: str
    #         Path to RECEIVER info file

    # Returns
    # =======
    # Nothing
    
    # Satellite Tracks
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SAT-X"],
    CorrIdx["SAT-Y"],
    CorrIdx["SAT-Z"],
    CorrIdx["ELEV"],
    CorrIdx["FLAG"]])

    # Configure plot and call plot generation function
    plotSatTracks(CorrFile, CorrData, RcvrInfo)

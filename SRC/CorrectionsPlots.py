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
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/CORR/figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)


# Plot Satellite Tracks
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
    PlotConf["LineWidth"] = 0.5

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
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = Longitude
    PlotConf["yData"][Label] = Latitude
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Fast and Long Term Corrections
def plotLTCandENTGPS(CorrFile, CorrData, SatData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Fast and Long Term Corrections", "SAT_FLT")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "Fast and Long Term Corrections [m]"
    PlotConf["yTicks"] = range(-7, 8)
    PlotConf["yLim"] = [-7, 7]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Label"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["Label"][Label] = 'LTC-X'
    PlotConf["Color"][Label] = 'skyblue'
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SatData[SatIdx["LTC-X"]]
    Label = 1
    PlotConf["Label"][Label] = 'LTC-Y'
    PlotConf["Color"][Label] = 'black'
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SatData[SatIdx["LTC-Y"]]
    Label = 2
    PlotConf["Label"][Label] = 'LTC-Z'
    PlotConf["Color"][Label] = 'grey'
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SatData[SatIdx["LTC-Z"]]
    Label = 3
    PlotConf["Label"][Label] = 'LTC-B'
    PlotConf["Color"][Label] = 'orange'
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SatData[SatIdx["LTC-B"]]
    Label = 4
    PlotConf["Label"][Label] = 'FC'
    PlotConf["Color"][Label] = 'blueviolet'
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SatData[SatIdx["FC"]]
    Label = 5
    PlotConf["Label"][Label] = 'ENT to GPS'
    PlotConf["Color"][Label] = 'seagreen'
    PlotConf["MarkerSize"][Label] = 1.5
    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["ENTtoGPS"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot ENT-GPS Offset
def plotENTtoGPS(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "ENT-GPS", "ENT_GPS")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "ENT-GPS [m]"
    PlotConf["yTicks"] = np.arange(-0.8, 0.8, 0.2)
    PlotConf["yLim"] = [-0.8, 0.6]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["Color"][Label] = 'orange'
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["ENTtoGPS"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Sigma FLT vs Elevation
def plotSigmaFLTvsElev(PreproObsFile, PreproObsData):
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
    usecols=[CorrIdx["SAT-X"], CorrIdx["SAT-Y"], CorrIdx["SAT-Z"], CorrIdx["ELEV"], CorrIdx["FLAG"]])

    # Configure plot and call plot generation function
    plotSatTracks(CorrFile, CorrData, RcvrInfo)

    # LTC and ENT-GPS
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["FLAG"], CorrIdx["ENTtoGPS"]])

    SatData = read_csv(SatFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[SatIdx["SOD"], SatIdx["LTC-X"], SatIdx["LTC-Y"], SatIdx["LTC-Z"], SatIdx["LTC-B"], SatIdx["FC"]])
    
    # Configure plot and call plot generation function
    plotLTCandENTGPS(CorrFile, CorrData, SatData)

    # ENT-GPS Offset
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["FLAG"], CorrIdx["ENTtoGPS"]])
    
    # Configure plot and call plot generation function
    plotENTtoGPS(CorrFile, CorrData)

    # Sigma FLT vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["SFLT"]])
    
    # Configure plot and call plot generation function
    plotSigmaFLTvsElev(CorrFile, CorrData)
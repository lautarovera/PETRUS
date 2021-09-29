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
from COMMON.Iono import computeIonoMappingFunction
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
    PlotConf["FigSize"] = (10.4,6.6)

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
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 90
    PlotConf["ColorBarStep"] = 10

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
def plotSigmaFLTvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma FLT", "SFLTvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "SigmaFLT [m]"
    PlotConf["yTicks"] = np.arange(0.8, 2.4, 0.2)
    PlotConf["yLim"] = [0.8, 2.2]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SFLT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot UIVD
def plotUIVD(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "UIVD", "UIVD")
    
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["LonStep"] = 15
    PlotConf["LatStep"] = 10

    PlotConf["zLabel"] = "Elevation [deg]"

    PlotConf["yTicks"] = np.arange(0.25, 2.5, 0.25)
    PlotConf["yLim"] = [0.25, 2.25]
    PlotConf["yLabel"] = "UIVD [m]"

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 90
    PlotConf["ColorBarStep"] = 10

    # Transform ECEF to Geodetic
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["UISD"]][FilterCond] / computeIonoMappingFunction(CorrData[CorrIdx["ELEV"]][FilterCond])
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot UISD
def plotUISD(CorrFile, CorrData, RcvrInfo):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "UISD", "UISD")
    
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    RcvrLat = int(RcvrInfo[RcvrIdx["LAT"]])
    RcvrLon = int(RcvrInfo[RcvrIdx["LON"]])

    PlotConf["LonMin"] = (RcvrLon - 25) 
    PlotConf["LonMax"] = (RcvrLon + 25) 
    PlotConf["LatMin"] = (RcvrLat - 15)
    PlotConf["LatMax"] = (RcvrLat + 15) if RcvrLat < 75 else 90

    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["zLabel"] = "UISD [m]"

    PlotConf["yTicks"] = range(PlotConf["LatMin"], PlotConf["LatMax"] + 1, PlotConf["LatStep"])
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"], PlotConf["LonMax"] + 1, PlotConf["LonStep"])
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "UISD [m]"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 6
    PlotConf["ColorBarStep"] = 0.5

    # Transform ECEF to Geodetic
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = CorrData[CorrIdx["IPPLON"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["IPPLAT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["UISD"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma UIRE
def plotSigmaUIREvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma UIRE", "SUIREvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma UIRE [m]"
    PlotConf["yTicks"] = np.arange(0.5, 3.5, 0.5)
    PlotConf["yLim"] = [0.5, 3.0]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SUIRE"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot STD
def plotSTD(CorrFile, CorrData, RcvrInfo):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Slant Tropo Delay", "STD")
    
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    RcvrLat = int(RcvrInfo[RcvrIdx["LAT"]])
    RcvrLon = int(RcvrInfo[RcvrIdx["LON"]])

    PlotConf["LonMin"] = (RcvrLon - 25) 
    PlotConf["LonMax"] = (RcvrLon + 25) 
    PlotConf["LatMin"] = (RcvrLat - 15)
    PlotConf["LatMax"] = (RcvrLat + 15) if RcvrLat < 75 else 90

    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["zLabel"] = "STD [m]"

    PlotConf["yTicks"] = range(PlotConf["LatMin"], PlotConf["LatMax"] + 1, PlotConf["LatStep"])
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"], PlotConf["LonMax"] + 1, PlotConf["LonStep"])
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "STD [m]"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 14
    PlotConf["ColorBarStep"] = 2

    # Transform ECEF to Geodetic
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 0.5
    PlotConf["xData"][Label] = CorrData[CorrIdx["IPPLON"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["IPPLAT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["STD"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Sigma TROPO
def plotSigmaTROPOvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma TROPO", "STROPOvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma TROPO [m]"
    PlotConf["yTicks"] = np.arange(0.0, 0.8, 0.1)
    PlotConf["yLim"] = [0.1, 0.7]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["STROPO"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Sigma Multipath
def plotSigmaMPvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma Multipath", "SMPvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma Multipath [m]"
    PlotConf["yTicks"] = np.arange(0.100, 0.375, 0.025)
    PlotConf["yLim"] = [0.125, 0.350]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SMP"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Sigma Noise Division
def plotSigmaNOISEDIVvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma Noise + Divergence", "SNOISEDIVvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma Noise + Divergence [m]"
    PlotConf["yTicks"] = np.arange(0.10, 0.45, 0.05)
    PlotConf["yLim"] = [0.10, 0.40]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SNOISEDIV"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Sigma Airborne
def plotSigmaAIRBORNEvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma Airborne", "SAIRvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma Airborne [m]"
    PlotConf["yTicks"] = np.arange(0.15, 0.55, 0.05)
    PlotConf["yLim"] = [0.15, 0.50]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SAIR"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Sigma UERE
def plotSigmaUEREvsElev(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Sigma UERE", "SUEREvsELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["yLabel"] = "Sigma UERE [m]"
    PlotConf["yTicks"] = np.arange(1.00, 3.50, 0.25)
    PlotConf["yLim"] = [1.00, 3.25]

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 100, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarStep"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SUERE"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Receiver Clock
def plotRcvrClk(CorrFile, CorrData):
    PlotConf = {}

    initPlot(CorrFile, PlotConf, "Receiver Clock Estimation", "RCVR_CLK")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "Receiver Clock Estimation [m]"

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1

    FilterCond = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["MarkerSize"] = {}
    Label = 0
    PlotConf["Color"][Label] = 'seagreen'
    PlotConf["MarkerSize"][Label] = 1
    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["RCVR-CLK"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


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
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SFLT"]])
    
    # Configure plot and call plot generation function
    plotSigmaFLTvsElev(CorrFile, CorrData)

    # UIVD
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["UISD"]])
    
    # Configure plot and call plot generation function
    plotUIVD(CorrFile, CorrData)

    # UISD
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["IPPLON"], CorrIdx["IPPLAT"], CorrIdx["FLAG"], CorrIdx["UISD"]])
    
    # Configure plot and call plot generation function
    plotUISD(CorrFile, CorrData, RcvrInfo)

    # Sigma UIRE vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SUIRE"]])
    
    # Configure plot and call plot generation function
    plotSigmaUIREvsElev(CorrFile, CorrData)

    # STD
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["IPPLON"], CorrIdx["IPPLAT"], CorrIdx["FLAG"], CorrIdx["STD"]])
    
    # Configure plot and call plot generation function
    plotSTD(CorrFile, CorrData, RcvrInfo)

    # Sigma TROPO vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["STROPO"]])
    
    # Configure plot and call plot generation function
    plotSigmaTROPOvsElev(CorrFile, CorrData)

    # Sigma Multipath vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SMP"]])
    
    # Configure plot and call plot generation function
    plotSigmaMPvsElev(CorrFile, CorrData)

    # Sigma Noise Division vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SNOISEDIV"]])
    
    # Configure plot and call plot generation function
    plotSigmaNOISEDIVvsElev(CorrFile, CorrData)

    # Sigma Airborne vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SAIR"]])
    
    # Configure plot and call plot generation function
    plotSigmaAIRBORNEvsElev(CorrFile, CorrData)

    # Sigma UERE vs Elevation
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["PRN"], CorrIdx["ELEV"], CorrIdx["FLAG"], CorrIdx["SUERE"]])
    
    # Configure plot and call plot generation function
    plotSigmaUEREvsElev(CorrFile, CorrData)

    # Receiver Clock
    # ----------------------------------------------------------
    # Read the cols we need from CORRECTIONS file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["FLAG"], CorrIdx["RCVR-CLK"]])
    
    # Configure plot and call plot generation function
    plotRcvrClk(CorrFile, CorrData)
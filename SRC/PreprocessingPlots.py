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
from InputOutput import PreproIdx
from InputOutput import REJECTION_CAUSE_DESC
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot, generateChallengePlot
import numpy as np
from collections import OrderedDict

def initPlot(PreproObsFile, PlotConf, Title, Label):
    PreproObsFileName = os.path.basename(PreproObsFile)
    PreproObsFileNameSplit = PreproObsFileName.split('_')
    Rcvr = PreproObsFileNameSplit[2]
    DatepDat = PreproObsFileNameSplit[3]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    PlotConf["xLabel"] = "Hour of Day %s" % Doy
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PPVE/figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)


# Plot Satellite Visibility
def plotSatVisibility(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Satellite Visibility", "SAT_VISIBILITY")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["zLabel"] = "Elevation [deg]"

    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yTicks"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yTicksLabels"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yLim"] = [0, max(unique(PreproObsData[PreproIdx["PRN"]]))]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = False
    PlotConf["DoubleAxis"] = False
    PlotConf["NotConv"] = True

    PlotConf["Marker"] = 'o'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["xDataNotConv"] = {}
    PlotConf["yDataNotConv"] = {}
    PlotConf["zDataNotConv"] = {}
    for prn in unique(PreproObsData[PreproIdx["PRN"]]):
        Label = "G" + ("%02d" % prn)
        FilterCond = ((PreproObsData[PreproIdx["PRN"]] == prn) & (PreproObsData[PreproIdx["STATUS"]] == 1))
        PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = PreproObsData[PreproIdx["PRN"]][FilterCond]
        PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

        FilterCond = ((PreproObsData[PreproIdx["PRN"]] == prn) & (PreproObsData[PreproIdx["STATUS"]] == 0))
        PlotConf["xDataNotConv"][Label] = PreproObsData[PreproIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yDataNotConv"][Label] = PreproObsData[PreproIdx["PRN"]][FilterCond]
        PlotConf["zDataNotConv"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
    
# Plot Number of Satellites
def plotNumSats(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Number of Satellites", "SAT_NUM")
    
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10.4,6.6)

    PlotConf["yLabel"] = "Number of Satellites"
    PlotConf["yTicks"] = range(0, 15, 2)
    PlotConf["yLim"] = [0, 14]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = True
    PlotConf["DoubleAxis"] = False

    PlotConf["Marker"] = '-'
    PlotConf["LineWidth"] = 1

    RawSats = []
    SmoothedSats = []

    for sod in sorted(unique(PreproObsData[PreproIdx["SOD"]])):
        RawFilterCond = (PreproIdx["SOD"] == sod)
        SmoothedFilterCond = ((PreproIdx["SOD"] == sod) & (PreproObsData[PreproIdx["STATUS"]] == 1))
        RawSats.append(len(PreproObsData[PreproIdx["SOD"]][RawFilterCond]))
        SmoothedSats.append(len(PreproObsData[PreproIdx["SOD"]][SmoothedFilterCond]))

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Label"] = {}
    Label = 0
    PlotConf["Label"][Label] = 'Raw'
    PlotConf["Color"][Label] = 'orange'
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = RawSats
    Label = 1
    PlotConf["Label"][Label] = 'Smoothed'
    PlotConf["Color"][Label] = 'green'
    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = SmoothedSats

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite Polar View
def plotSatPolarView(PreproObsFile, PreproObsData):
    PlotConf = {}

    initPlot(PreproObsFile, PlotConf, "Satellite Polar View", "SAT_POLAR_VIEW")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (10,10)

    PlotConf["zLabel"] = "GPS-PRN"

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 0
    PlotConf["ColorBarMax"] = 32

    PlotConf["tData"] = {}
    PlotConf["rData"] = {}
    PlotConf["zData"] = {}
    PlotConf["tData"] = np.deg2rad(PreproObsData[PreproIdx["AZIM"]])
    PlotConf["rData"] = PreproObsData[PreproIdx["ELEV"]]
    PlotConf["zData"] = PreproObsData[PreproIdx["PRN"]]

    # Call generatePlot from Plots library
    generateChallengePlot(PlotConf)

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

def generatePreproPlots(PreproObsFile):
    
    # Purpose: generate output plots regarding Preprocessing results

    # Parameters
    # ==========
    # PreproObsFile: str
    #         Path to PREPRO OBS output file

    # Returns
    # =======
    # Nothing
    
    # Satellite Visibility
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["PRN"], PreproIdx["ELEV"], PreproIdx["STATUS"]])

    plotSatVisibility(PreproObsFile, PreproObsData)

    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"], PreproIdx["STATUS"]])

    plotNumSats(PreproObsFile, PreproObsData)

    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["ELEV"], PreproIdx["AZIM"], PreproIdx["PRN"]])

    plotSatPolarView(PreproObsFile, PreproObsData)
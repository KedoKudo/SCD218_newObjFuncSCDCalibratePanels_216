# necessary import
import numpy as np
from collections import namedtuple
from mantid.simpleapi import *


def convert(dictionary):
    return namedtuple('GenericDict', dictionary.keys())(**dictionary)


import os
directory = os.path.dirname(os.path.realpath(__file__))
lc_natrolite = {
    "a": 18.29,  # A
    "b": 18.64,  # A
    "c": 6.56,  # A
    "alpha": 90,  # deg
    "beta": 90,  # deg
    "gamma": 90,  # deg
}
natrolite = convert(lc_natrolite)
pws = Load("data/Natrolite_FindPeaks_unrounded.nxs")
LoadIsawUB('pws','data/Natrolite_300K.mat')

SCDCalibratePanels(
    PeakWorkspace='pws',
    a=natrolite.a,
    b=natrolite.b,
    c=natrolite.c,
    alpha=natrolite.alpha,
    beta=natrolite.beta,
    gamma=natrolite.gamma,
    CalibrateT0=False,
    CalibrateL1=True,
    CalibrateBanks=False,
    OutputWorkspace="testCaliTable",
    DetCalFilename=directory + "/test.DetCal",
    CSVFilename=directory + "/test.csv",
    XmlFilename=directory + "/test.xml",
    ToleranceOfTranslation=1e-16,
    ToleranceOfReorientation=1e-16,
    TranslationSearchRadius=0.05,
    RotationSearchRadius=5,
    SourceShiftSearchRadius=0.1,
    VerboseOutput=True,
)
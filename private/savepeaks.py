from mantid.simpleapi import *
from mantid.geometry import SymmetryOperationFactory
import numpy as np

IPTS = 23019
iptsfolder = "/SNS/CORELLI/IPTS-"+str(IPTS)
nxfiledir = iptsfolder+"/nexus/"
sharedir = iptsfolder+"/shared/"
ccfiledir = iptsfolder+"/shared/autoreduce/"
scriptdir = iptsfolder + "/shared/Natrolite/"
outputdir = scriptdir

start = 133752
stop = 133812
stop = 133760

UBfile="Natrolite_300K.mat"

runs = range(start, stop+1,1)

FirstPass = 1
toMerge = []
LoadCC = False

combinedPeaks = CreatePeaksWorkspace()
for r in runs:
    print('Processing run : %s' %r)
    ows = 'COR_'+str(r)
    matrixfile = outputdir+ows+'_UB.mat'
    toMerge.append(ows)

    if LoadCC :
        filename = ccfiledir+'CORELLI_'+str(r)+'_elastic.nxs'
        if not mtd.doesExist(ows):
            LoadNexus(Filename=filename, OutputWorkspace=ows)
            #MaskDetectors(Workspace=ows, MaskedWorkspace=maskfile)
    else:
        filename = nxfiledir+'CORELLI_'+str(r)+'.nxs.h5'
        if not mtd.doesExist(ows):
            #print('bank'+str(bankidx) )
            #LoadEventNexus(Filename=filename, OutputWorkspace=ows,BankName='bank'+str(bankidx) )
            LoadEventNexus(Filename=filename, OutputWorkspace=ows)
        #MaskDetectors(Workspace=ows, MaskedWorkspace=maskfile)

    # get total proton_charge from run log
    owshandle = mtd[ows]
    lrun = owshandle.getRun()
    pclog = lrun.getLogData('proton_charge')
    pc = sum(pclog.value)/1e12
    #owshandle /= pc
    print('the current proton charge :'+ str(pc))


    #CopyInstrumentParameters(calibration, OutputWorkspace=ows)
    SetGoniometer(ows,Axis0="BL9:Mot:Sample:Axis3,0,1,0,1")
    if FirstPass == True:
        LoadIsawUB(InputWorkspace=ows, Filename=scriptdir+UBfile)
        #LoadIsawUB(InputWorkspace=ows, Filename=matrixfile)
    else:
        LoadIsawUB(InputWorkspace=ows, Filename=matrixfile)

    outputmd = ows+'_MDqsample'
    #if not mtd.doesExist(outputmd):
    ConvertToMD(InputWorkspace=ows,
                OutputWorkspace=outputmd,
                QDimensions="Q3D",
                dEAnalysisMode="Elastic",
                Q3DFrames="Q_sample",
                LorentzCorrection=1,
                PreprocDetectorsWS='-',
                MinValues="-15,-15,-15",
                MaxValues="15,15,15",
                Uproj='1,0,0',
                Vproj='0,1,0',
                Wproj='0,0,1')

    PredictPeaks(InputWorkspace=ows,
                 MinDSpacing=.3,
                 MaxDSpacing=10,
                 #ReflectionCondition='Rhombohedrally centred, obverse',
                 #ReflectionCondition='Primitive',
                 ReflectionCondition='All-face centred',
                 OutputWorkspace=ows+'_predictedPeaks')

    IntegratePeaksMD(InputWorkspace=ows+'_MDqsample',
                     PeakRadius=0.10,
                     BackgroundInnerRadius=0.10,
                     BackgroundOuterRadius=0.14,
                     PeaksWorkspace=ows+'_predictedPeaks',
                     OutputWorkspace=ows+'_integratedPeaks',
                     IntegrateIfOnEdge=True)
    FilterPeaks(InputWorkspace=ows+'_integratedPeaks',
                OutputWorkspace=ows+'_filteredPeaks',
                FilterVariable='Intensity',
                FilterValue=100,
                Operator='>')

    CentroidPeaksMD(InputWorkspace=ows+'_MDqsample',
                    PeakRadius=0.14,
                    PeaksWorkspace=ows+'_filteredPeaks',
                    OutputWorkspace=ows+'_centeredPeaks')

    #IndexPeaks(PeaksWorkspace=ows+'_centeredPeaks', RoundHKLs=True)

    CombinePeaksWorkspaces(LHSWorkspace=ows+'_centeredPeaks',
                           RHSWorkspace='combinedPeaks',
                           OutputWorkspace='combinedPeaks')

OptimizeLatticeForCellType(PeaksWorkspace='combinedPeaks',
                            #CellType = 'Orthorhombic',
                            CellType = 'Hexagonal',
                            #CellType = 'Cubic',
                            Apply=True)

lattice = mtd['combinedPeaks'].sample().getOrientedLattice()
print(lattice)
SaveIsawPeaks(InputWorkspace='combinedPeaks',
              Filename=scriptdir+'Natrolite_300K.peaks')

#DeleteWorkspace(ows)
#DeleteWorkspace(ows+'_centeredPeaks')
#DeleteWorkspace(ows+'_integratedPeaks')
#DeleteWorkspace(ows+'_filteredPeaks')
#DeleteWorkspace(ows+'_predictedPeaks')
#DeleteWorkspace(ows+'_MDqsample')

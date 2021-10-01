'''
Created on 10 juin 2021

@author: Fab
'''
import os
from LMT.USV.vocImporter3.vocImporter import fileBiggerSplit, getDataFileMatch
from experimental.voc.analysis.quantif.Voc import Voc
import numpy as np

def getAllWavData( folder , limit=None ):
    print( "Extracting wav files numbers in " , folder )
    files = os.listdir( folder )
    
    wavDataList = []    
    
    for file in files:

        # file name example: 2021-03-25_17-33-43_0000048_p_099_l_3_c_2.wav
        # file is split by _ symbol.
        # we take the biggest string as the number of the string.
        if file.endswith(".wav"):
            print( "Loading data for ", file )
            number = fileBiggerSplit ( file )            
            dataFile = getDataFileMatch( files, number )
            
            if dataFile == None:
                print("getAllWavData: Can't load data. Quit()")
                quit()
                
            wavDataList.append ( WavData( folder+"/"+file, folder+"/"+dataFile ) )
            if limit!=None:
                if len( wavDataList ) == limit:
                    break
                
    
    return wavDataList

def safe( number ):
    if np.isnan( number ):
        return 0
    if np.isinf( number ):
        return 0

    return number

class WavData(object):
    
    def __init__(self, wavFile, dataFile ):
        
        self.wavFile = wavFile
        self.dataFile = dataFile
        self.vocList = []

        with open( dataFile ) as f:
            lines = f.readlines()
            
        for line in lines:
            #print( line )
            data = line.split( ";")
            if data[0].isnumeric():
                #print( "Voc number : " , data[0] )
                
                voc = Voc( )
                                                
                voc.durationMs = float( data[3] )
                
                voc.frequencyDynamicHz = float( data[4] )                
                
                voc.startFrequencyHz = float( data[5] )
                
                voc.endFrequencyHz = float( data[6] )    
                
                voc.diffStartEndFrequencyHz = voc.startFrequencyHz-voc.endFrequencyHz            
                
                voc.meanFrequencyHz= float( data[7] )
                voc.frequencyTVHz = float ( data[8] )                
                voc.meanFrequencyTVHz = float ( data[9] )
                voc.linearityIndex = float( data[10] )                
                voc.meanPower = float( data[11] )
                
                voc.nbModulation = float( data[12] )
                voc.nbPtHarmonics = int( data[13] )
                voc.nbJump = int ( data[14] )
                voc.isModulated = data[15] == "true"
                voc.isShort = data[16] == "true"
                voc.isUpward = data[17] == "true"
                voc.isDownward = data[18] == "true"

                voc.minFrequency = float( data[21] )
                voc.maxFrequency = float( data[22] )
                voc.peakPower = float( data[23] )
                voc.peakFrequency = float( data[24] )
                voc.minPower = float( data[25] )
                
                self.vocList.append ( voc )
                
        self.nbVoc = safe ( len( self.vocList ) )

        self.meanPower = safe ( np.mean( [v.meanPower for v in self.vocList] ) )
        self.stdPower = safe ( np.mean( [v.meanPower for v in self.vocList] ) ) 
        
        self.meanDuration = safe ( np.mean( [v.durationMs for v in self.vocList] ) )
        self.stdDuration = safe ( np.std( [v.durationMs for v in self.vocList] ) )
        
        self.meanFrequency = safe ( np.std( [v.meanFrequencyHz for v in self.vocList] ) )
        self.stdFrequency = safe ( np.std( [v.meanFrequencyHz for v in self.vocList] ) )
        
        self.meanPeakPower = safe ( np.std( [v.peakPower for v in self.vocList] ) )
        self.stdPeakPower = safe ( np.std( [v.peakPower for v in self.vocList] ) )

        self.meanMinPower = safe ( np.std( [v.minPower for v in self.vocList] ) )
        self.stdMinPower = safe ( np.std( [v.minPower for v in self.vocList] ) )
        
        self.meanMinFrequency = safe ( np.std( [v.minFrequency for v in self.vocList] ) )
        self.stdMinFrequency = safe ( np.std( [v.minFrequency for v in self.vocList] ) )

        self.meanMaxFrequency = safe ( np.std( [v.maxFrequency for v in self.vocList] ) )
        self.stdMaxFrequency = safe ( np.std( [v.maxFrequency for v in self.vocList] ) )
        
        self.meanFrequencyTVHz = safe ( np.std( [v.frequencyTVHz for v in self.vocList] ) )
        self.stdFrequencyTVHz = safe ( np.std( [v.frequencyTVHz for v in self.vocList] ) )


        
        
        
        
    def getAttributes(self):
        
        attributes = [self.nbVoc ,
                      self.meanPower,
                      self.stdPower,
                      self.meanDuration,
                      self.stdDuration,
                      self.meanFrequency,
                      self.stdFrequency,
                      self.meanPeakPower,
                      self.stdPeakPower,
        
                      self.meanMinPower,
                      self.stdMinPower,
                      self.meanMinFrequency,
                      self.stdMinFrequency,
                      self.meanMaxFrequency,
                      self.stdMaxFrequency,
                      self.meanFrequencyTVHz,
                      self.stdFrequencyTVHz
                      
                      ]
                      
                      
                      
            
                   
        return attributes
        
    def getAttributeLabels(self):
        return [
            "nbVoc",
            "meanPower",
            "stdPower",
            "meanDuration",
            "stdDuration",
            "meanFrequency",
            "stdFrequency",
            "meanPeakPower",
            "stdPeakPower",
            "meanMinPower",
            "stdMinPower",
            "meanMinFrequency",
            "stdMinFrequency",
            "meanMaxFrequency",
            "stdMaxFrequency",
            "meanFrequencyTVHz",
            "stdFrequencyTVHz",
            ]
        
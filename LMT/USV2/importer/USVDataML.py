'''
Created on 10 juin 2021

@author: Fab
'''
import os

from LMT.USV2.importer.Voc import Voc
import numpy as np
from LMT.USV2.importer.importUtil import fileBiggerSplit, getDataFileMatch



def getAllUSV_ML_DataForWav( folder , limit=None ):
    print( "Extracting wav files numbers in " , folder )
    files = os.listdir( folder )
    
    USVDataList = []    
    
    for file in files:

        # file name example: 2021-03-25_17-33-43_0000048_p_099_l_3_c_2.wav
        # file is split by _ symbol.
        # we take the biggest string as the number of the string.
        if file.endswith(".wav"):
            print( "Loading data for ", file )
            number = fileBiggerSplit ( file )            
            dataFile = getDataFileMatch( files, number )
            
            if dataFile == None:
                print("getAllUSV_ML_DataForWav: Can't load data. Quit()")
                quit()
                
            USVDataList.extend ( grabUSVDataML( folder+"/"+file, folder+"/"+dataFile ) )
            if limit!=None:
                if len( USVDataList ) >= limit:
                    print( "Limit reached : " , limit )
                    break                
    
    return USVDataList

def safe( number ):
    if np.isnan( number ):
        return 0
    if np.isinf( number ):
        return 0

    return number

def grabUSVDataML( wavFile, dataFile ):
    
    USVDataML_List = []
    
    with open( dataFile ) as f:
        lines = f.readlines()
        
    for line in lines:
        data = line.split( ";")
        if data[0].isnumeric():
            
            voc = Voc( )
                                            
            voc.startOffsetMs = float( data[2] )
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

            voc.minFrequency = float( data[21] )
            voc.maxFrequency = float( data[22] )
            voc.peakPower = float( data[23] )
            voc.peakFrequency = float( data[24] )
            voc.minPower = float( data[25] )
            
            USVDataML_List.append ( USVDataML( voc ) )
    
    return USVDataML_List        
            
class USVDataML(object):
    
    def __init__(self, voc ):
        
        self.voc = voc
                                      
    def getAttributes(self):
        
        attributes = [
            self.voc.durationMs,                
            self.voc.frequencyDynamicHz,                                
            self.voc.startFrequencyHz,
            self.voc.endFrequencyHz,    
            self.voc.diffStartEndFrequencyHz,            
            self.voc.meanFrequencyHz,
            self.voc.frequencyTVHz,                
            self.voc.meanFrequencyTVHz,
            self.voc.linearityIndex,                
            self.voc.meanPower,
            self.voc.nbModulation,
            self.voc.nbPtHarmonics,
            self.voc.nbJump,
    
            self.voc.minFrequency,
            self.voc.maxFrequency,
            self.voc.peakPower,
            self.voc.peakFrequency,
            self.voc.minPower
        ]
                   
        return attributes
    
    def getAttributeLabels(self):
        return [
            "durationMs",                
            "frequencyDynamicHz",                                
            "startFrequencyHz",
            "endFrequencyHz",    
            "diffStartEndFrequencyHz",            
            "meanFrequencyHz",
            "frequencyTVHz",                
            "meanFrequencyTVHz",
            "linearityIndex",                
            "meanPower",
            "nbModulation",
            "nbPtHarmonics",
            "nbJump",
    
            "minFrequency",
            "maxFrequency",
            "peakPower",
            "peakFrequency",
            "minPower"
            ]  
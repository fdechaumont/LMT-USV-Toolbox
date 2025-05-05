'''
Created on 23 janv. 2025

@author: Fab
'''



from LMT.USV2.importer.Voc import Voc

import pickle
import os

import time
from LMT.USV2.importer.USVDataML import USVDataML

from LMT.USV2.importer.Wav import Wav
from scipy.io import wavfile

import glob


class MachineLearningTester(object):
    


    def __init__(self ):
    

        print("This script will export detected vocalizations as spectrum image files.")
        print("It uses machine learning to select if USV should be exported or not")
        
        print("Loading predictor...")
        self.rf = pickle.load( open("trainingSet.bin", 'rb'))
        
        self.accuracyList = []
        self.accuracyStdList = []
        self.featureImportanceList = []
        
        self.accuracyList.append( self.rf.accuracy )
        self.accuracyStdList.append( self.rf.accuracyError )
        self.featureImportanceList.append( self.rf.feature_importances )
                
        print("******** Results")
        
        print( "Accuracy: ", self.accuracyList )
        print( "Accuracy std: " , self.accuracyStdList )
        print( "Feature importance: " , self.featureImportanceList )
        
        print("Machine learning : Predictor Loaded.")
        
        print( "--------")
        
        print("Select wav file to export vocalizations as images.")
        
        self.preTriggerFrameMs = 1000
        self.nbOfPreTriggerFrame = self.preTriggerFrameMs / 33.33
        
        print("The pre trigger is set to " , self.preTriggerFrameMs , " ms")


    def grabDataForFile(self, inputFile):
        # look for any data concerning the file, and also the crop like file.wav_crop_0.0_50.0.txt and others
        
        folder = os.path.dirname( inputFile )
        #files = os.listdir( folder )
        print ( folder )
        files = glob.glob( f"{folder}/*.txt")
        
        targetingFile = inputFile.split(".")[0]
        
        #print( files )
        lines = []
        
        for file in files:
            #print( file )
            #startFile=inputFile.split("/")[-1]
            
            #print( "start file " , targetingFile )
            
            if targetingFile in file:
                print( "--> ok: " + file )
                
            
                with open( file ) as f:
                    l = f.readlines()
                    lines.extend( l )
        
        return lines

    def process( self, audioFile ):
        
        
        
        print( "File is : ", audioFile )

        
        lines = self.grabDataForFile( audioFile )
        
        print( f"lines:{len(lines)}")
        
        vocList = []
        
        for line in lines:
            #print( line )
            data = line.split( ";")
            if data[0].isnumeric():
                #print( "Voc number : " , data[0] )
                
                voc = Voc( )
                voc.fileName = audioFile
                voc.vocNumber = int ( data[0] )
                voc.nbPoint = int( data[1] )
                startOffsetMs = float( data[2] )
                voc.startOffsetMs = startOffsetMs                
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
                #voc.power = float( cols[23] )
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
                voc.isInBadRepeat = data[26] == "true"
                
                # the avisoft trigger to consider is the first one in the database
                # we remove the number of frame used in pre-trigger.
                                
                #voc.setStartFrame( int ( startFrame - nbOfPreTriggerFrame + startOffsetMs/33.3333 ) )
        
                # machine learning USV test:
                USVmachineLearningAttributes = USVDataML( voc ).getAttributes()
                pred = self.rf.clf.predict( [ USVmachineLearningAttributes ] )[0]
                
                if pred == 1: 
                    #eventTimeLineVoc.addEvent( voc.getAsEvent(), noCheck=True )
                    vocList.append( voc )
                    #renderVoc( voc, f"{audioFile}" )                
                    print("->USV")
                else:
                    print("->not USV")
 
        return vocList
    
    
    
    
    
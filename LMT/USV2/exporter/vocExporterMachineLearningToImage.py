'''
Created on 15 avr. 2021

@author: Fab
'''
from tkinter.filedialog import askopenfile
from tkinter import Tk


import sqlite3
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import EventTimeLine
from LMT.USV2.importer.Voc import Voc

import pickle
import os

import time
from LMT.USV2.importer.USVDataML import USVDataML
from LMT.USV2.importer.importUtil import fileBiggerSplit, getDataFileMatch
from LMT.USV2.importer.Wav import Wav
from scipy.io import wavfile



    
def renderVoc( vocList , saveImageFile ):
    
           
        
        wav = Wav( vocList[0].fileName ) # TODO: cache this buffer just like in the player
            
        index = 0    
        for voc in vocList:        
        
            print(  voc.fileName )
            startMs = voc.startOffsetMs 
            durationMs = voc.durationMs 
            
            minT = startMs/1000
            maxT = (startMs+durationMs) / 1000
            
            print( "Save image : ", saveImageFile )
            text = str( startMs ) + " - " + str( durationMs ) 
            print( text )
            
            #text = time.strftime('%d-%H:%M:%S', time.gmtime( voc.startFrame*(1/33) ) )
            #text = text.replace( "-","\n" )
            #text+= "\nf:"+ str( voc.startFrame )
            #text+="\n"+str( voc.metadata["burstNumber"] )+" - " + str( voc.metadata["vocNumber"] )
            try:                
                wav.saveSpectrumImage( f"{saveImageFile}_spectrum_{index}.png", minT, maxT , text = text,textColor="black", save=True, vmin=-15, vmax=35 )
                index+=1
            
            except:
                print("Can't save spectrum: " , saveImageFile )


def grabDataForFile(inputFile):
    # look for any data concerning the file, and also the crop like file.wav_crop_0.0_50.0.txt and others
    
    folder = os.path.dirname( inputFile )
    files = os.listdir( folder )
    #print( files )
    lines = []
    
    for file in files:
        print( file )
        startFile=inputFile.split("/")[-1]
        #print( "start file " , startFile )
        
        if file.endswith(".txt") and file.startswith( startFile ):
            print( "--> ok: " + file )
        
            with open( folder+"/"+file ) as f:
                l = f.readlines()
                lines.extend( l )
    return lines
    
if __name__ == '__main__':
    
    print("This script will export detected vocalizations as spectrum image files.")
    print("It uses machine learning to select if USV should be exported or not")
    
    print("Loading predictor...")
    rf = pickle.load( open("trainingSet.bin", 'rb'))
    
    accuracyList = []
    accuracyStdList = []
    featureImportanceList = []
    
    accuracyList.append( rf.accuracy )
    accuracyStdList.append( rf.accuracyError )
    featureImportanceList.append( rf.feature_importances )
            
    print("******** Results")
    
    print( "Accuracy: ", accuracyList )
    print( "Accuracy std: " , accuracyStdList )
    print( "Feature importance: " , featureImportanceList )
    
    print("Machine learning : Predictor Loaded.")
    
    print( "--------")
    
    print("Select wav file to export vocalizations as images.")
    
    preTriggerFrameMs = 1000
    nbOfPreTriggerFrame = preTriggerFrameMs / 33.33
    
    print("The pre trigger is set to " , preTriggerFrameMs , " ms")
    
    root = Tk()
    root.lift()
    root.withdraw()
    inputFile = askopenfile( mode="r" , defaultextension=".wav" )
    
    if inputFile == None:
        print("No file selected")
        quit() 
    
    startTime = time.time()
    
    audioFile = inputFile.name
    print( "File is : ", audioFile )
    #path = os.path.dirname( audioFile )
    #print( path )
    
    #audioFilePrefix = inputFile.name.split(".wav")[0]
    #print( audioFilePrefix )
    
    lines = grabDataForFile( audioFile )
    
    
    
    fileIndex = 0
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
            pred = rf.clf.predict( [ USVmachineLearningAttributes ] )[0]
            
            if pred == 1: 
                #eventTimeLineVoc.addEvent( voc.getAsEvent(), noCheck=True )
                vocList.append( voc )
                #renderVoc( voc, f"{audioFile}" )                
                print("->USV")
            else:
                print("->not USV")
 
    renderVoc( vocList , audioFile )
    
    totalTime = time.time() - startTime
    import datetime    
    print("Execution time (h:m:s):" , str(datetime.timedelta(seconds=totalTime)) )
    print("All done.")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
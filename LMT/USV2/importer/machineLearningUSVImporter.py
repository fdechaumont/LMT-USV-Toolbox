'''
Created on 10 juin 2021

@author: Fab
'''

import pandas as pd

import pickle
import os
from shutil import copyfile
from pathlib import Path


import matplotlib.pyplot as plt
from matplotlib import patches
from LMT.USV2.importer.USVDataML import getAllUSV_ML_DataForWav, grabUSVDataML

from LMT.USV2.lib.Wav import Wav
from LMT.USV2.importer.WavData import getAllWavData
from LMT.USV2.importer.randomForestTester import RandomForestTester

def training( ):
    
    folderVoc = "D:/USV training stuff with experiment 843/voc - strict/voc"
    folderNoise = "D:/USV training stuff with experiment 843/noise"

    dataVoc = getAllUSV_ML_DataForWav( folderVoc , limit = None ) #limit = 10000
    dataNoise = getAllUSV_ML_DataForWav( folderNoise , limit = None )
    
    attributeVoc = []
    for d in dataNoise:
        a = d.getAttributes()
        a.append( 0 ) # noise targetId
        attributeVoc.append( a )
    
    for d in dataVoc:
        a = d.getAttributes()
        a.append( 1 ) # voc targetId
        attributeVoc.append( a )
    
    # machine learning
    
    rows_list = attributeVoc
    
    df = pd.DataFrame(rows_list)
    
    labels = dataVoc[0].getAttributeLabels()
    labels.append( "target" )
    df.columns = labels #["A","B","C","target"]
    
    print ( df ) 
    
    
    rf = RandomForestTester( df, testCLF=True , showAsATree = False , showConfusionMatrix=False,realClassNames=["Noise","Voc"] )    
        
    pickle.dump( rf, open("trainingSet.bin", 'wb'))
    
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

def copyFile( source, target ):
    print( source , "->" , target )
    copyfile(source, target )
        

def predictor( saveJPG = False ):

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
    
    print("- Starting prediction ")
    
    #folderToProcess = "D:/20210421_usv_2days_B6_F23_F24_Experiment 843/usv/all"
    folderToProcess = "D:/20210423_usv_2days_B6_F25_F26_Experiment 6791/usv"
    
    folderVoc = str ( os.path.dirname( folderToProcess ) ) +"/autoVoc"
    folderNoise = str ( os.path.dirname( folderToProcess ) ) +"/autoNoise"
    
    print( folderVoc )
    print( folderNoise )
    
    Path(folderVoc).mkdir(parents=True, exist_ok=True)
    Path(folderNoise).mkdir(parents=True, exist_ok=True)        

    
    '''
    folderVoc = "D:/20210421_usv_2days_B6_F23_F24_Experiment 843/usv/voc"
    folderNoise = "D:/20210421_usv_2days_B6_F23_F24_Experiment 843/usv/noise"
    '''

    dataWavFileList = getAllWavData( folderToProcess , limit = 1000 )
    #dataWavFileList = getAllWavData( folderToProcess , limit = None )
    
    for dataWavFile in dataWavFileList:
        print( "File: " , dataWavFile.wavFile )
        
        USVList = grabUSVDataML( dataWavFile.wavFile , dataWavFile.dataFile )
        redBars = [] # pred 0
        greenBars = [] # pred 1
        patcheList = []
        for usvData in USVList:
            a = usvData.getAttributes()
            pred = rf.clf.predict( [ a ] )[0]
            print( pred )
            minX = usvData.voc.startOffsetMs/1000
            maxX = usvData.voc.startOffsetMs/1000 + usvData.voc.durationMs/1000
            minY = usvData.voc.minFrequency
            maxY = usvData.voc.maxFrequency
            
            minX-=0.02
            maxX+=0.02
            minY-=10000
            maxY+=10000
            
            if pred == 0: # noise
                redBars.append( minX )
                redBars.append( maxX )
                patcheList.append( patches.Rectangle((minX, minY), maxX-minX, maxY-minY, linewidth=1, edgecolor='r', facecolor='none') )
            if pred == 1: # voc
                greenBars.append( minX )
                greenBars.append( maxX )
                patcheList.append( patches.Rectangle((minX, minY), maxX-minX, maxY-minY, linewidth=1, edgecolor='g', facecolor='none') )            
        
        if ( saveJPG == True ):
            if len( USVList ) > 0:
                wav = Wav( dataWavFile.wavFile )
                wav.saveSpectrumImage( dataWavFile.wavFile+".jpg" , save = False )        
                ax = plt.gca()
                for patch in patcheList:
                    ax.add_patch( patch )
                
                for v in redBars:
                    ax.axvline( v , color="red" , alpha=0.8  )  # , lw=1000
                for v in greenBars:
                    ax.axvline( v , color="green" , alpha=0.8  )
                
                try:
                    plt.savefig( dataWavFile.wavFile+".jpg", bbox_inches='tight')
                except:
                    print( "Figure too large !")
                    
                plt.close()
        
        
    
if __name__ == '__main__':
    
    print("Test of machine learning per USV in each file")
    #predictor( saveJPG = False )
    
    # train the system with your own USVs classification
    
    training()

    print("All Done.")
    
    
    
    
    
'''
Created on 15 avr. 2021
Modified by Elodie on 23 Jan 2026

@author: Fab
'''
from tkinter.filedialog import askopenfile
from tkinter import Tk
import os

import sqlite3
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import EventTimeLine
from LMT.USV2.importer.importUtil import fileBiggerSplit, getDataFileMatch
from LMT.USV2.importer.Voc import Voc
from copy import deepcopy

import pandas as pd
import numpy as np



        
def getAllUSVSeq( connection ):
    
    # build the dictionary containing all the USV seq recorded live in LMT.
    # extracts from the description of the USV seq the startFrame and the number contained in the file. 
    # as the USV seq are queried in startframe DESC, we only keep the first record of each voc record attempt from avisoft, which is the one to consider.
    # result is a dictionary with k= the number in the file name and the value= corresponding startFrame.
    
    c = connection.cursor()

    query = "SELECT * FROM EVENT WHERE NAME='USV seq' ORDER BY STARTFRAME DESC";
    print( query )
    c.execute( query )
    all_rows = c.fetchall()

    resultDic = {}
    
    for row in all_rows:
        
        startFrame = row[3]
        file = row[2]
        number = fileBiggerSplit ( file )
        resultDic[number] = startFrame
        
    return resultDic

        
if __name__ == '__main__':
    
    print("This script will import vocalizations analysed with Usvseg in a LMT database.")    
    print("Your USVs should be located in databasefile.sqlite/usv/ folder.")
    print("All analysed files that are exported in csv format in this folder will be imported.")
    print("!!!! check the pre trigger duration before launching the script.")
    print("Select database to import vocalizations.")
    
    extension = "_dat.csv"
    
    preTriggerFrameMs = 1000
    nbOfPreTriggerFrame = preTriggerFrameMs / 33.33
    
    print("The pre trigger is set to " , preTriggerFrameMs , " ms")
    
    root = Tk()
    root.lift()
    root.withdraw()
    inputFile = askopenfile( mode="r" , defaultextension=".sqlite" )
    
    if inputFile == None:
        print("No file selected")
        quit() 
    
    
    dataBaseFile = inputFile.name
    print( "File is : ", dataBaseFile )
    path = os.path.dirname(dataBaseFile)

    # open database
    print("Open database...")
    connection = sqlite3.connect( dataBaseFile )
    animalPool = AnimalPool()
    animalPool.loadAnimals( connection )


    print( "Path: " , path )
    #allUSVPath = path + "/usv/voc/"
    allUSVSelected = path + "/usv/"
    
    print( "Extracting wav files numbers in usv/")
    numberList = []
    files = os.listdir( allUSVSelected )
    
    number2WaveFile = {}
    
    for file in files:
        print( file )
        # file name example: 2021-03-25_17-33-43_0000048_p_099_l_3_c_2.wav
        # file is split by _ symbol.
        # we take the biggest string as the number of the string.
        if file.endswith(".wav"):
            number = fileBiggerSplit ( file )
            numberList.append ( number )
            number2WaveFile[ number ] = allUSVSelected+"/"+file
        
    print("Selected numbers: " , numberList, len(numberList) )
    
    print("Number of voc selected : " , len( numberList ))
    
    # Import
    
    # get all .xls files in the USV folder:
    dataFiles = []
    
    dataFiles = os.listdir( allUSVSelected )
    numberListShort = deepcopy(numberList)
    for number in numberList:
        
        dataFile = getDataFileMatch( dataFiles, number, extension=extension )
        
        if dataFile == None:
            #if the analysis file is not found, remove this wav file from the list to be considered
            print( f"Error: data {extension} file not found for number: " , number )
            numberListShort.remove(number)
            
        print( "Match : " , dataFile )
        
    print(f"All data files found ({len(numberListShort)}).")
    
    
    # creating VOC timeline
    print ( "creating VOC timeline" )
    
    eventTimeLineVoc = EventTimeLine( connection, "Voc US", loadEvent=False )
    
    # get USV Seq recorded with avisoft triggering
    print("Loading USV seq from database...")
    USVSeqDic = getAllUSVSeq( connection )
    for k in USVSeqDic:
        #print( "USV seq: n=" , k , "  startFrame=" , USVSeqDic[k] )
        pass
        
    # parsing DATA
    
    print("Parsing data...")
    for number in numberListShort:
        print("####", int(number))
        if int(number) < 300:
            if not number in USVSeqDic:
                print("Warning: the key #" , number , " does not exist in wav files.")
                continue
            
            # retrieve database USV seq        
            startFrame = USVSeqDic[number]
            
            dataFile = allUSVSelected + "/" + getDataFileMatch( dataFiles, number, extension=extension )
            print ( "Processing datafile: " , dataFile )
            if extension == "_dat.csv":
                automaticDf = pd.read_csv(f'{dataFile}', sep=',', decimal=".")
                #pd.read_csv(f'{path}/{fileName}_{software}.csv', sep=',', decimal=".")
                for index, row in automaticDf.iterrows():
                    voc = Voc( )
                    startTime = np.round( row["start"] *1000 ) #convert into milliseconds and round the values
                    endTime = np.round( row["end"] *1000 ) #convert into milliseconds and round the values
                    voc.startOffsetMs = startTime
                    voc.endTime = endTime
                    voc.durationMs = np.round( row["duration"] ) #convert into milliseconds and round the values
                    #voc.fileName = row["File"]
                    voc.vocNumber = row["#"]
                    voc.label = row["#"]
                    voc.maxFrequency = row["maxfreq"]
                    voc.maxAmplitude = row["maxamp"]
                    voc.meanPeakFrequency = row["meanfreq"]
                    voc.cvFrequency = row["cvfreq"]
                    
                    
                    # the avisoft trigger to consider is the first one in the database
                    # we remove the number of frame used in pre-trigger.
                    
                    
                    voc.setStartFrame( int ( startFrame - nbOfPreTriggerFrame + startTime/33.3333 ) )
     
                    eventTimeLineVoc.addEvent( voc.getAsEvent(), noCheck=True )
                    
                    print(voc.getAsEvent().metadata)
 
    # remove the vocs that are doubled because of pre-filter effect
    
    for voc1 in eventTimeLineVoc.eventList:
        for voc2 in eventTimeLineVoc.eventList:
            if voc1!=voc2:
                if not "remove" in voc2.metadata:
                    if voc1.overlapInT( voc2.startFrame-1, voc2.endFrame+1 ):
                        print( "overlap...")
                        print( "voc1 : " , voc1.metadata["fileName"] , voc1.metadata["vocNumber"] , voc1.metadata["durationMs"] )
                        print( "voc2 : " , voc2.metadata["fileName"] , voc2.metadata["vocNumber"] , voc2.metadata["durationMs"] )
                        print( "Check duration...")
                        if voc1.metadata["durationMs"] == voc2.metadata["durationMs"]:
                            print("Duration match > remove first one")
                            voc1.metadata["remove"] = True
    
    print( "------ remove list:")
    for voc in eventTimeLineVoc.eventList[:]:
        if "remove" in voc.metadata:
            print( "voc : " , voc.metadata["fileName"] , voc.metadata["vocNumber"] , "(removing)")
            eventTimeLineVoc.eventList.remove( voc )
 
    print( "Saving vocs...")
    eventTimeLineVoc.endRebuildEventTimeLine(connection, deleteExistingEvent = True )
               
    
    print("All done.")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
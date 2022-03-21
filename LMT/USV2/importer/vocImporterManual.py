'''
Created on 15 avr. 2021

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


        
def getAllUSVSeq( connection ):
    
    # build the dictionnary containing all the USV seq recorded live in LMT.
    # extracts from the description of the USV seq the startFrame and the number contained in the file. 
    # as the USV seq are queried in startframe DESC, we only keep the first record of each voc record attempt from avisoft, which is the one to consider.
    # result is a dictionnary with k= the number in the file name and the value= corresponding startFrame.
    
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
    
    print("This script will import vocalizations in an LMT database.")    
    print("Your USVs should be located in databasefile.sqlite/usv/voc/ folder")
    print("All wav file in this folder will be imported")
    print("Select database to import vocalizations.")
    
    
    
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
    allUSVSelected = path + "/usv/voc/"
    
    print( "Extracting wav files numbers in usv/voc")
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
        
    print("Selected numbers: " , numberList )
    
    print("Number of voc selected : " , len( numberList ))
    
    # Import
    
    # get all .txt files in the USV folder:
    dataFiles = []
    
    dataFiles = os.listdir( allUSVSelected )
    for number in numberList:
        
        dataFile = getDataFileMatch( dataFiles, number )
        
        if dataFile == None:
            print( "Error: data text file not found for number: " , number )
            print( "Quit." )
            quit()
            
        print( "Match : " , dataFile )
        
    print("All data files found.")
    
    # creating VOC timeline
    print ( "creating VOC timeline" )
    
    eventTimeLineVoc = EventTimeLine( connection, "Voc", loadEvent=False )
    
    # get USV Seq recorded with avisoft triggering
    print("Loading USV seq from database...")
    USVSeqDic = getAllUSVSeq( connection )
    for k in USVSeqDic:
        print( "USV seq: n=" , k , "  startFrame=" , USVSeqDic[k] )
        
    # parsing DATA
    
    print("Parsing data...")
    for number in numberList:
        
        if not number in USVSeqDic:
            print("Warning: they key #" , number , " does not exist in wav files.")
            continue
        
        # retreive database USV seq        
        startFrame = USVSeqDic[number]
        
        dataFile = allUSVSelected + "/" + getDataFileMatch( dataFiles, number )
        print ( "Processing datafile: " , dataFile )
        with open( dataFile ) as f:
            lines = f.readlines()

        for line in lines:
            #print( line )
            data = line.split( ";")
            if data[0].isnumeric():
                #print( "Voc number : " , data[0] )
                
                voc = Voc( )
                voc.fileName = number2WaveFile[number]
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
                
                
                voc.setStartFrame( int ( startFrame - nbOfPreTriggerFrame + startOffsetMs/33.3333 ) )
 
                eventTimeLineVoc.addEvent( voc.getAsEvent(), noCheck=True )
 
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
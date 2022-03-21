'''
Created on 17 janv. 2020

@author: Fab
'''

import os

from lmtanalysis.Util import convert_to_d_h_m_s
from lmtanalysis.Animal import AnimalPool


def getWavFileName( experimentFile, USVFile ):
    '''
    As the user may move the folder containing the datafile, this function provides the real location of wavfiles
    using the experimentFile and the USVFile. It provides the wav file located in the voc subfolder of the corresponding database.sqlite
    
    somefolder/voc/wavfile.wav
    somefolder/database.sqlite
    '''    
    wavFileName = os.path.dirname(os.path.abspath( experimentFile ) )
    wavFileName += "\\voc\\" + os.path.basename(USVFile)
    return wavFileName


def getStrainAgeSexPairGenoPerFile ( connection ):
    """
    This function returns the information of the animals tested in the file given as argument.
    """
    pool = AnimalPool( )
    pool.loadAnimals( connection )
    
    strainFile = pool.animalDictionnary[1].strain
    #print( "strain: ", strainFile )
    
    ageFile = pool.animalDictionnary[1].age
    #print( "age: ", ageFile )
    
    if pool.animalDictionnary[1].sex == pool.animalDictionnary[2].sex:
            sexFile = pool.animalDictionnary[1].sex
    else:
        sexFile = "mixed"
    #print( "sex: ", sexFile )
    
    rfid1 = int(pool.animalDictionnary[1].RFID[-7:])
    rfid2 = int(pool.animalDictionnary[2].RFID[-7:])
    minRFID = min ( rfid1, rfid2 )
    maxRFID = max ( rfid1, rfid2 )
     
    pairFile = "{}-{}".format( minRFID, maxRFID )
    #print( "pair ", pairFile )
    
    geno1 = pool.animalDictionnary[1].genotype
    geno2 = pool.animalDictionnary[2].genotype
    
    if geno1 == geno2:
        genotype = "{}-{}".format( geno1, geno2 )
    else:
        genoSortedList = sorted ( [geno1, geno2] )
        genotype = "{}-{}".format( genoSortedList[0], genoSortedList[1] )

    strainAgeSexPairGeno = ( strainFile, ageFile, sexFile, pairFile, genotype )
    print( "file info: ", *strainAgeSexPairGeno)
    return strainAgeSexPairGeno
    
def cleanVocLegacy( eventTimeLineVoc ):
    
    toRemove = []
    for event in eventTimeLineVoc.eventList:
        
        if "excluded" in event.metadata:
        
            if event.metadata["excluded"]==True:
                toRemove.append( event )
    
                     
    print("Number of voc events removed: ", len(toRemove ))
    for event in toRemove:
        eventTimeLineVoc.eventList.remove( event )
    
def cleanVocAdvanced2( eventTimeLineVoc ):
    
    '''
    Shift voc of -30 frames for hold of 1 second
    Shift voc of -9 frames because of 300ms latency of AviSoft
    '''
    
    for event in eventTimeLineVoc.eventList:
        event.startFrame-=39
        event.endFrame-=39
            
    ''' Correct voc with doublon's duration '''
    # {'vocNumber': 1, 'nbPoint': 55, 'startOffsetMs': 1347.4133, 'durationMs': 46.08, 'frequencyDynamicHz': 15234.375, 'startFrequencyHz': 74414.06, 
    # 'endFrequencyHz': 88476.56, 'diffStartEndFrequencyHz': -14062.5, 'meanFrequencyHz': 82952.766, 'frequencyTVHz': 18164.062, 'meanFrequencyTVHz': 330.25568, 'linearityIndex': 3.124535913905807, 
    # 'meanPower': 1.3660032, 'nbModulation': 4.0, 'nbPtHarmonics': 0, 'nbJump': 0, 'isModulated': True, 'isShort': False, 'isUpward': False, 'isDownward': False, 'minFrequency': 74414.06,
    # 'maxFrequency': 89648.44, 'peakPower': 4.8890347, 'peakFrequency': 88476.56, 'minPower': 0.26254624, 'isInBadRepeat': False, 
    # 'fileName': 'D:\\all_data_usv_pairs\\20190125_usv_lmt_pair_F13-F14_2mo_1\\usv_021\\F13-F14_2mo_ch1\\voc\\T2019-01-28_10-12-17_0000497_p_084_l_1_c_2.wav', 
    # 'excluded': False, 'griffIndex': 0.0, 'slope': 23171.37369497132, 'slopeNormalized': 502.8509916443429, 'maxMeanSpeedInCage': 0.01687314391374318,
    # 'maxMeanSpeedInCageW150': 0.05838735967925878, 'maxMaxSpeedInCage': 0.03374628782748636, 'maxSpeedInCage': 0.03374628782748636, 'maxSpeedInCageW60': 0.12140807211786962}

    
        
    oldFile = "other val"
    currentFile = "some val"
    offsetMs = 0
    nbRepeat = 0
    toRemove = []
    
    for event in eventTimeLineVoc.eventList:
        
        # test if we changed file:
        currentFile = event.metadata["fileName"]
        if currentFile != oldFile:
            print("new file: " , currentFile[-30:] )
            oldFile = currentFile
            offsetMs = 0 # reset offset
            nbRepeat=0 # reset repeat number            
        
        # shift event
        frameCorrection = int ( offsetMs/33 )
        event.startFrame+= frameCorrection
        print( "Correction: " , frameCorrection , " ms:", offsetMs )
        
        if event.metadata["isInBadRepeat"]:
            n = event.metadata["vocNumber"]
            print( "repeat found at position " , n , "duration : " , event.metadata["durationMs"] )
            
            # each repeat is about 120ms. So for each repeat found we remove 60ms.
            # this is an estimation. Most double buffer will affect 2 vocs
            # la correction n'est pas idéale, il faudrait avoir l'info des zones exactes des doublons dans le fichier,
            # mais c'est pas enregistré à la detection.
            offsetMs -= 60 
            
            if nbRepeat == 1:
                print("removed")                
                toRemove.append( event )
                nbRepeat = 0
            else:
                nbRepeat+=1
                     
    print("Number of voc events removed: ", len(toRemove ))
    for event in toRemove:
        eventTimeLineVoc.eventList.remove( event )

'''    
def cleanVocAdvanced2( eventTimeLineVoc ):
    
    
    # Shift voc of -30 frames for hold of 1 second
    # Shift voc of -9 frames because of 300ms latency of AviSoft
    
    
    for event in eventTimeLineVoc.eventList:
        event.startFrame-=39
        event.endFrame-=39
            
    # Correct voc with doublon's duration
    # {'vocNumber': 1, 'nbPoint': 55, 'startOffsetMs': 1347.4133, 'durationMs': 46.08, 'frequencyDynamicHz': 15234.375, 'startFrequencyHz': 74414.06, 
    # 'endFrequencyHz': 88476.56, 'diffStartEndFrequencyHz': -14062.5, 'meanFrequencyHz': 82952.766, 'frequencyTVHz': 18164.062, 'meanFrequencyTVHz': 330.25568, 'linearityIndex': 3.124535913905807, 
    # 'meanPower': 1.3660032, 'nbModulation': 4.0, 'nbPtHarmonics': 0, 'nbJump': 0, 'isModulated': True, 'isShort': False, 'isUpward': False, 'isDownward': False, 'minFrequency': 74414.06,
    # 'maxFrequency': 89648.44, 'peakPower': 4.8890347, 'peakFrequency': 88476.56, 'minPower': 0.26254624, 'isInBadRepeat': False, 
    # 'fileName': 'D:\\all_data_usv_pairs\\20190125_usv_lmt_pair_F13-F14_2mo_1\\usv_021\\F13-F14_2mo_ch1\\voc\\T2019-01-28_10-12-17_0000497_p_084_l_1_c_2.wav', 
    # 'excluded': False, 'griffIndex': 0.0, 'slope': 23171.37369497132, 'slopeNormalized': 502.8509916443429, 'maxMeanSpeedInCage': 0.01687314391374318,
    # 'maxMeanSpeedInCageW150': 0.05838735967925878, 'maxMaxSpeedInCage': 0.03374628782748636, 'maxSpeedInCage': 0.03374628782748636, 'maxSpeedInCageW60': 0.12140807211786962}

    
        
    oldFile = "other val"
    currentFile = "some val"
    offsetMs = 0
    nbRepeat = 0
    toRemove = []
    
    for event in eventTimeLineVoc.eventList:
        
        # test if we changed file:
        currentFile = event.metadata["fileName"]
        if currentFile != oldFile:
            print("new file: " , currentFile[-30:] )
            oldFile = currentFile
            offsetMs = 0 # reset offset
            nbRepeat=0 # reset repeat number            
        
        # shift event
        frameCorrection = int ( offsetMs/33 )
        event.startFrame+= frameCorrection
        print( "Correction: " , frameCorrection , " ms:", offsetMs )
        
        if event.metadata["isInBadRepeat"]:
            n = event.metadata["vocNumber"]
            print( "repeat found at position " , n , "duration : " , event.metadata["durationMs"] )
            
            # each repeat is about 120ms. So for each repeat found we remove 60ms.
            # this is an estimation. Most double buffer will affect 2 vocs
            # la correction n'est pas idéale, il faudrait avoir l'info des zones exactes des doublons dans le fichier,
            # mais c'est pas enregistré à la detection.
            offsetMs -= 60 
            
            if nbRepeat == 1:
                print("removed")                
                toRemove.append( event )
                nbRepeat = 0
            else:
                nbRepeat+=1
                     
    print("Number of voc events removed: ", len(toRemove ))
    for event in toRemove:
        eventTimeLineVoc.eventList.remove( event )
'''

def cleanVoc( eventTimeLineVoc , advancedCorrection = True ):
    
    if advancedCorrection:
        cleanVocAdvanced2(eventTimeLineVoc)
        cleanVocLegacy(eventTimeLineVoc) # temp
    else:
        cleanVocLegacy(eventTimeLineVoc)
   

def sortVocInTime( eventTimeLineVoc ):
    #print("Sorting vocs...")
    eventTimeLineVoc.eventList.sort(key=lambda x: x.startFrame )
    
def frameToTimeTicker(x, pos):
   
    vals= convert_to_d_h_m_s( x )
    return "D{0} - {1:02d}:{2:02d}".format( int(vals[0])+1, int(vals[1]), int(vals[2]) )
    
def formatTimeAxisForTimeLine( ax ):
    
    import matplotlib.ticker as ticker
    formatter = ticker.FuncFormatter( frameToTimeTicker )
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(labelsize=6 )
    ax.xaxis.set_major_locator(ticker.MultipleLocator( 30 * 60 * 60 * 12 ))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator( 30 * 60 * 60 ))


def checkIfOverlapWith(eventBehavToCheck, vocDictionnary):
    for t in range(eventBehavToCheck.startFrame, eventBehavToCheck.endFrame + 1):
        if t in vocDictionnary:
            return True

    return False
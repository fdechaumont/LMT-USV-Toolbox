'''
Created on 17 févr. 2020

@author: Fab
'''
#from LMT.USV.lib.vocUtil import sortVocInTime
#from LMT.USV.burster.USVBurst import USVBurst
from LMT.USV.lib.vocUtil import sortVocInTime
from LMT.USV.burster.USVBurst import USVBurst



'''
The burster creates burst from voc.
It takes a vocEventList as entry, and then create bursts 
'''


def createBurstFromVoc( eventTimeLineVoc, silenceBetweenBurstMs = 750 ):
        
    sortVocInTime( eventTimeLineVoc )
    
    burstList = []
    
    currentBurst = None
    
    for vocEvent in eventTimeLineVoc.eventList:
        
        if currentBurst != None:  
            #rationale behind the division by 33: the difference between endFrame and startFrame is in frames, so the silenceBetweenBurstMs
            #should be converted into frames, and one frame lasts 33 ms
            if vocEvent.startFrame - currentBurst.getEndFrame() < silenceBetweenBurstMs/33:
                currentBurst.vocEventList.append( vocEvent )
                continue
            else:
                #burstList.append( currentBurst )
                currentBurst = None
                
        if currentBurst == None:
            
            currentBurst = USVBurst()
            burstList.append( currentBurst )
            currentBurst.vocEventList.append( vocEvent )        
            continue
    
    # indexing voc in burst:
    burstNumber = 0
    for burst in burstList:
        vocNumber=0
        for voc in burst.vocEventList:
            voc.metadata["vocNumber"]= vocNumber
            voc.metadata["burstNumber"]= burstNumber            
            vocNumber+=1
        burstNumber+=1
    
    return burstList


if __name__ == '__main__':
    
    pass
    
'''
Created on 17 févr. 2020

@author: Fab
'''
import numpy as np
from lmtanalysis.Event import Event

class USVBurst(object):
    '''
    classdocs
    '''

    def __init__(self ):
        
        '''
        Constructor
        '''
        self.vocEventList = []
        self.metadata = {}
        
    def getEndFrame(self):

        endFrame = None
        
        for vocEvent in self.vocEventList:
            if ( endFrame == None or vocEvent.endFrame > endFrame ):
                endFrame = vocEvent.endFrame
                
        return endFrame

    def getStartFrame(self):

        startFrame = None
        
        for vocEvent in self.vocEventList:
            if ( startFrame == None or vocEvent.startFrame < startFrame ):
                startFrame = vocEvent.startFrame
                
        return startFrame
    
    def getDurationMs(self):
        
        startFrame = self.vocEventList[0].startFrame
        endFrame = self.vocEventList[-1].endFrame
        
        return (endFrame-startFrame) * 33
    
    def getDuty(self):
        
        startFrame = self.vocEventList[0].startFrame
        endFrame = self.vocEventList[-1].endFrame

        totalDurationFrame = 0
        for vocEvent in self.vocEventList:
            totalDurationFrame+= vocEvent.duration()
        
        ratio = totalDurationFrame / ( endFrame-startFrame +1 )
        return ratio
        
    def getDurationList(self):
        
        durationList = []
        for vocEvent in self.vocEventList:                
            durationList.append( vocEvent.metadata["durationMs"] )
        return durationList
        
    def getIntervalList(self):
        
        intervalList = []
        
        previous=None
        
        for vocEvent in self.vocEventList:
            if previous==None:
                previous=vocEvent
                continue
            intervalList.append( vocEvent.startFrame - previous.endFrame )
            previous = vocEvent            
            
        return intervalList
        
    def getTraitList(self):
        return [ "durationMs", "nbUSV", "meanDuration", "stdDuration", "meanInterval","stdInterval","meanFrequency" ]
    
    def getValue(self, valueName ):
        
        if valueName == "durationMs":
            return self.getDurationMs()

        if valueName == "nbUSV":
            return len( self.vocEventList )
        
        if valueName == "meanDuration":
            
            return np.mean( self.getDurationList( ) )

        if valueName == "stdDuration":
            
            durationList = self.getDurationList( )
            if len(durationList ) > 0:
                return np.std( durationList )
            
            return None
        
        if valueName == "meanInterval":
            intervalList = self.getIntervalList()
            
            if len( intervalList ) > 0:
                return np.mean( intervalList )
            
            return None
        
        if valueName == "stdInterval":
            
            intervalList = self.getIntervalList()
            
            if len( intervalList ) > 0:
                return np.std( intervalList )
            
            return None

        if valueName == "meanFrequency":

            meanFrequencyList = []
            for vocEvent in self.vocEventList:
                meanFrequencyList.append( vocEvent.metadata["meanFrequencyHz"] )
                                 
            return np.mean( meanFrequencyList )
        
        return None
        
            
    def __str__(self):
        
        return "Burst: nbVoc:\t"+str( len( self.vocEventList ) )+ "\tstartFrame:\t" + str( self.getStartFrame() )+ "\tendFrame:\t" + str(  self.getEndFrame( ) )
    
    def getAsEvent( self ):
        
        usvBurstEvent = Event( self.getStartFrame(), self.getEndFrame() )
        for trait in self.getTraitList():
            usvBurstEvent.metadata[trait]= self.getValue( trait )
    
    
        return usvBurstEvent
    
    
    
    
    
    
    
    
    
    
    
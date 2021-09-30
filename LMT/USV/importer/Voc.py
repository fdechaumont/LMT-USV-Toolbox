'''
Created on 23 juil. 2019

@author: Fab
'''
from lmtanalysis.Event import Event


class Voc:
    '''
    a Voc class that store voc quantification
    '''

    def __init__(self ):
        '''
        Constructor
        '''
        self.startFrame = None
        self.endFrame = None
        
        self.vocNumber = 0
        self.nbPoint = 0
        self.startOffsetMs = 0
        self.durationMs = 0
        self.frequencyDynamicHz = 0
        self.startFrequencyHz = 0    
        self.endFrequencyHz = 0
        self.diffStartEndFrequencyHz = 0    
        self.meanFrequencyHz = 0
        self.frequencyTVHz = 0
        self.meanFrequencyTVHz = 0
        self.linearityIndex = 0
        self.meanPower = 0
        self.nbModulation = 0
        self.nbPtHarmonics = 0 
        self.nbJump = 0
        self.isModulated = False
        self.isShort = False
        self.isUpward = False
        self.isDownward = False
        self.minFrequency = 0
        self.maxFrequency = 0
        self.peakPower = 0
        self.peakFrequency = 0
        self.minPower = 0
        self.isInBadRepeat = False
        self.fileName = ""
        
    def setStartFrame(self,startFrame):
        self.startFrame = startFrame
        self.endFrame = int( startFrame + self.durationMs/33.3 )
    
    def getAsEvent(self):
        if self.startFrame == None:
            return None
        if self.endFrame == None:
            return None
        event = Event( self.startFrame , self.endFrame )
        
        event.metadata["vocNumber"] = self.vocNumber 
        event.metadata["nbPoint"] = self.nbPoint
        event.metadata["startOffsetMs"] = self.startOffsetMs
        event.metadata["durationMs"] = self.durationMs        
        event.metadata["frequencyDynamicHz"] = self.frequencyDynamicHz
        event.metadata["startFrequencyHz"] = self.startFrequencyHz
        event.metadata["endFrequencyHz"] = self.endFrequencyHz
        event.metadata["diffStartEndFrequencyHz"] = self.diffStartEndFrequencyHz
        event.metadata["meanFrequencyHz"] = self.meanFrequencyHz
        event.metadata["frequencyTVHz"] = self.frequencyTVHz        
        event.metadata["meanFrequencyTVHz"] = self.meanFrequencyTVHz
        event.metadata["linearityIndex"] = self.linearityIndex
        event.metadata["meanPower"] = self.meanPower
        event.metadata["nbModulation"] = self.nbModulation
        event.metadata["nbPtHarmonics"] = self.nbPtHarmonics
        event.metadata["nbJump"] = self.nbJump
        event.metadata["isModulated"] = self.isModulated
        event.metadata["isShort"] = self.isShort
        event.metadata["isUpward"] = self.isUpward
        event.metadata["isDownward"] = self.isDownward
        event.metadata["minFrequency"] = self.minFrequency
        event.metadata["maxFrequency"] = self.maxFrequency
        event.metadata["peakPower"] = self.peakPower
        event.metadata["peakFrequency"] = self.peakFrequency
        event.metadata["minPower"] = self.minPower
        event.metadata["isInBadRepeat"] = self.isInBadRepeat
        event.metadata["fileName"] = self.fileName
        return event
    
    def __str__(self):
        tab="\t"
        text ="Voc number:" + tab + str( self.vocNumber ) + "\n"
        text ="number of point defined:" + tab + str( self.nbPoint ) + "\n"
        text+="Start offset:" +tab + str( self.startOffsetMs ) + " ms \n"
        text+="Duration:"+ tab + str( self.durationMs ) + " ms \n" 
        
        text+="Frequency dynamic:"+ tab + str( self.frequencyDynamicHz ) + " Hz \n" 
        text+="Start Frequency:"+ tab + str( self.startFrequencyHz ) + " Hz \n" 
        text+="End Frequency:"+ tab + str( self.endFrequencyHz ) + " Hz \n" 
        text+="Mean Frequency:"+ tab + str( self.meanFrequencyHz ) + " Hz \n" 

        text+="Frequency total variation:"+ tab + str( self.frequencyTVHz ) + " Hz \n" 
        text+="Mean Frequency TV:"+ tab + str( self.meanFrequencyTVHz ) + " Hz \n" 
        
        text+="Linearity index:"+ tab + str( self.linearityIndex ) + " \n" 
        text+="mean power:"+ tab + str( self.meanPower ) + " \n" 

        text+="Number of modulation:"+ tab + str( self.nbModulation ) + " \n" 
        text+="Number of harmonics points:"+ tab + str( self.nbPtHarmonics ) + " \n" 
        text+="Number of jumps:"+ tab + str( self.nbJump ) + " \n" 
        text+="isModulated:"+ tab + str( self.isModulated ) + " \n" 
        text+="isShort:"+ tab + str( self.isShort ) + " \n" 
        text+="isUpward:"+ tab + str( self.isUpward ) + " \n" 
        text+="isDownward:"+ tab + str( self.isDownward ) + " \n" 
                    
        text+="minFrequency(Hz):"+ tab + str( self.minFrequency ) + " \n" 
        text+="maxFrequency(Hz):"+ tab + str( self.maxFrequency ) + " \n" 
        text+="peak power:"+ tab + str( self.peakPower ) + " \n" 
        text+="peak frequency(Hz):"+ tab + str( self.peakFrequency ) + " \n" 
    
        text+="min power:"+ tab + str( self.minPower ) + " \n" 
        text+="isInBadRepeat:"+ tab + str( self.isInBadRepeat ) + " \n" 
    
        text+="startFrame in experiment:"+ tab + str( self.startFrame ) + " \n" 
        text+="endFrame in experiment:"+ tab + str( self.endFrame ) + " \n"
        
        text+="filename:"+ tab + str( self.fileName ) + " \n"  
    
                    
    
        return text
        
    
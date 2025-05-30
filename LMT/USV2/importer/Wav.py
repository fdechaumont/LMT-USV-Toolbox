'''
Created on 23 mars 2020

@author: Fab
'''


from experimental.voc.analysis.quantif.VocPool import VocPool
from experimental.voc.analysis.quantif.Voc import Voc
from experimental.voc.analysis.quantif.Burst import Burst
from experimental.voc.analysis.quantif.ExperimentList import getDataBaseList
import os
import sqlite3
from experimental.voc.analysis.quantif.VocPool import VocPool
from experimental.voc.analysis.quantif.Voc import Voc
from experimental.voc.analysis.quantif.Burst import Burst
from tkinter.filedialog import askopenfilename
import sqlite3
import os
import numpy as np
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import Event, EventTimeLine, deleteEventTimeLineInBase
from lmtanalysis.Measure import oneDay, oneHour

from lmtanalysis.Util import getMinTMaxTAndFileNameInput
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
import pandas
import scipy.stats as stats

from numpy import sum

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from lmtanalysis.EventTimeLineCache import EventTimeLineCached


from scipy import signal
from scipy.io import wavfile
import os
import wave
import pylab

import numpy as np
from random import random
from experimental.voc.analysis.quantif.playAudio import VocTools
from experimental.voc.analysis.quantif.playAudio.VocTools import removeDoublon,\
    get_wav_info, overlap

#from LMT.USV.doublon.tagDoublon import process
#from LMT.USV.lib.vocUtil import cleanVoc, getWavFileName
#from LMT.USV.experimentList.experimentList import getExperimentList
import time
import matplotlib


# TODO Constant to put in another file

Fs = 300000    
'''
NFFT = int(Fs*0.005)    
noverlap = int( Fs*0.004 )
'''

Fs = 300000
NFFT = 1024
noverlap = int ( 1024 * 0.75 )

class Wav():
    
    def __init__ (self , wavFileName ):
        
        print( "Opening wav: " + str( wavFileName ) )
        '''
        self.__wavFile = wavFileName
                
        wav = wave.open( wavFileName, 'r')
                
        frames = wav.readframes(-1)
        #frames = wav.readframes(300000*3)
        self.__wavBuffer = pylab.frombuffer(frames, 'int16')
        
        self.__frameRate = wav.getframerate()
        self.__durationS = len( self.__wavBuffer ) / self.__frameRate      
        wav.close()
        '''
        self.__wavFile = wavFileName
        print( "wav: reading " , wavFileName )
        rate, audio = wavfile.read( wavFileName, mmap=True )
        print( rate )
        #print( rate, audio )
        print("quit")
        #quit()
        try:            
            audio = np.mean( audio, axis=1 )
        except:
            # file is not multi-channel
            print( "File is not multi channel")            
        
        self.__wavBuffer= audio # 
        
        self.__frameRate = rate
        self.__durationS = len( self.__wavBuffer ) / self.__frameRate
        

    def __str__ (self):
        durationString = str( round ( self.__durationS , 3 ) ) 
        return "Wav file: " + str( self.__wavFile ) + " rate: " + str( self.__frameRate ) + "Hz / Duration (s): " + str ( durationString ) 
    
    def getWavBuffer(self):
        return self.__wavBuffer
        
    def getFrameRate(self):
        return self.__frameRate

    def getWaveFile(self):
        return self.__waveFile

    def getDurationS(self):
        return self.__durationS
    
    def showWaveform( self , ax=None, minT=None, maxT = None ):
        data = self.getWavBuffer()
        
        if minT!=None:
            data = data[int(minT*300000):int(maxT*300000)]
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        timePoint = []
        dataPercent = []
        t = minT
        for d in data:
            dataPercent.append( d / 320 )
            timePoint.append( t )
            t+= 1/300000
                    
        ax.plot( timePoint, dataPercent , linewidth=0.5, color="black" )
                
        #ax.set_ylim( -32000,32000 )
        ax.set_ylim( -100,100 )
        ax.set_ylabel("Dynamic (%)" )    
        ax.set_xlabel("time (s)")
    
    def getMaxWaveformValue(self , minT , maxT ):
        # return the max percentage of dynamic
        
        data = self.getWavBuffer()
        
        if minT!=None:
            data = data[int(minT*300000):int(maxT*300000)]
        
        absData = []
        for d in data:
            absData.append( abs( d ) )
        
        return max(absData) / 320
    
    def show( self , display=True, ax=None, minT=None, maxT = None, fastMode = False , minLUT = -5, maxLUT = 10 ):
        
        '''
        wavBufferRenderedSize = len(self.__wavBuffer)
        if minT != None:
            wavBufferRenderedSize = ( maxT - minT ) * 300000
        '''
        
        if ax!=None:
            display=False
        
        if ax==None:    
            pylab.figure(num=None, figsize=( int ( self.getDurationS() *2 ) , 4))
            pylab.subplot(111)
            pylab.set_cmap( 'Greys')
            ax = pylab.gca()            
                
        if fastMode:
            #ax.specgram( self.getWavBuffer()[int(minT*300000):int(maxT*300000)], Fs=Fs , NFFT=NFFT, noverlap=noverlap ,vmin=-5, vmax=10 , cmap= 'Greys' )
            ax.specgram( self.getWavBuffer()[int(minT*300000):int(maxT*300000)], Fs=Fs , NFFT=NFFT, noverlap=noverlap ,vmin= minLUT, vmax=maxLUT , cmap= "Greys" )                
            ax.set_xlim( 0, maxT-minT )
        else:
            #pylab.specgram( self.getWavBuffer(), Fs=Fs , NFFT=NFFT, noverlap=noverlap ,  vmin=-5, vmax = 10 )
            ax.specgram( self.getWavBuffer(), Fs=Fs , NFFT=NFFT, noverlap=noverlap ,vmin=-5, vmax=10 , cmap= 'Greys' )        
            #ax.specgram( self.getWavBuffer(), Fs=Fs , NFFT=NFFT, noverlap=noverlap ,vmin=-5, vmax=10  )            
            if minT != None:            
                ax.set_xlim( minT, maxT )

        ax.set_xlabel("Seconds")
        # Only show ticks on the left and bottom spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
            
        '''
        x= range(0,10)
        y = range( 0, 150000 , int( 150000/10 ) )
        plt.plot( x , y )
        '''
        
        if display:
            pylab.show()

    def saveSpectrumImageSlowDown(self, fileName , minT=None, maxT = None ):

        
        slowDownFactor = 16    

        wavBufferRenderedSize = len(self.__wavBuffer)
        
        if minT != None:
            wavBufferRenderedSize = ( maxT - minT ) * 300000 / 16
        
            
        pylab.figure(num=None, figsize=( 4*16*10*wavBufferRenderedSize / 19041900 , 4))
        pylab.subplot( 111 )
        
        if minT != None:
            # convert to second 
            minT /= slowDownFactor
            maxT /= slowDownFactor
        
        pylab.set_cmap( 'Greys')
        pylab.specgram( self.__wavBuffer, Fs=Fs , NFFT=NFFT, noverlap=noverlap ,  vmin=-5, vmax = 10 ) 

        ax = plt.gca()
        ax.xaxis.set_ticks(np.arange(0, len(self.__wavBuffer)/300000, 1 ))
        
        #ax.xaxis.set_ticks(np.arange(0, len(self.__wavBuffer)/300000, 1/ ( (300/slowDownFactor) / 5) )) # 1 tick each 5 seconds.

        #formatter = matplotlib.ticker.FuncFormatter(lambda s, x: time.strftime('%M:%S', time.gmtime( s*(300/slowDownFactor) )))
        
        
        formatter = matplotlib.ticker.FuncFormatter(lambda s, x: time.strftime('%M:%S', time.gmtime( s * slowDownFactor ) ) )
        ax.xaxis.set_major_formatter(formatter)
        
        # *( 19041900 *25*60 ) ))
        
        
        '''
        from matplotlib.ticker import MultipleLocator #, NullLocator
        minorLocator   = MultipleLocator(1/((300/slowDownFactor)))
        ax.xaxis.set_minor_locator(minorLocator)
        '''

        if minT != None:            
            ax.set_xlim( minT, maxT )

        pylab.savefig( fileName, bbox_inches='tight')
        pylab.close()
            
        
    def saveSpectrumImage(self, fileName , minT=None, maxT = None , text= None , textColor= "black", save=True , vmin= -5 , vmax = 10 ):
        # default vmin and vmax are high contrast.
        # more subtle view with vmin=-15,vmax=35: 
                
        import time
        import matplotlib
        
        marginS = 0.05 # margin in second
        
        if minT != None:
            minT-=marginS
        
        if maxT != None:
            maxT+=marginS
        
        slowDownFactor = 1

        wavBufferRenderedSize = len(self.__wavBuffer)
        if minT != None:
            wavBufferRenderedSize = ( maxT - minT ) * 300000
            
        pylab.figure(num=None, figsize=( 8*16*10*wavBufferRenderedSize / 19041900 , 4))
        #pylab.figure(num=None, figsize=( 4*16*10*wavBufferRenderedSize / 19041900 , 4))
        pylab.subplot( 111 )
        
        if minT != None:
            # convert to second 
            minT /= slowDownFactor
            maxT /= slowDownFactor
        
        pylab.set_cmap( 'Greys')
        pylab.specgram( self.__wavBuffer, Fs=Fs , NFFT=NFFT, noverlap=noverlap ,  vmin=vmin, vmax = vmax ) 

        ax = plt.gca()
        ax.xaxis.set_ticks(np.arange(0, len(self.__wavBuffer)/300000, 1 ))
        
        #ax.xaxis.set_ticks(np.arange(0, len(self.__wavBuffer)/300000, 1/ ( (300/slowDownFactor) / 5) )) # 1 tick each 5 seconds.

        #formatter = matplotlib.ticker.FuncFormatter(lambda s, x: time.strftime('%M:%S', time.gmtime( s*(300/slowDownFactor) )))
        
        
        formatter = matplotlib.ticker.FuncFormatter(lambda s, x: time.strftime('%M:%S', time.gmtime( s * slowDownFactor ) ) )
        ax.xaxis.set_major_formatter(formatter)
        
        # *( 19041900 *25*60 ) ))
        
        
        '''
        from matplotlib.ticker import MultipleLocator #, NullLocator
        minorLocator   = MultipleLocator(1/((300/slowDownFactor)))
        ax.xaxis.set_minor_locator(minorLocator)
        '''

        # Only show ticks on the left and bottom spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_yaxis().set_visible(False)

        if minT != None:            
            ax.set_xlim( minT, maxT )
            
        if text != None:
            ax.text( minT+marginS, 145000 , text , color=textColor, va="top" )            

        if minT != None and maxT!=None:
            ax.axvline( minT+marginS , color="lightblue" , alpha=0.5 )
            ax.axvline( maxT-marginS , color="lightblue" , alpha=0.5 )

        if save:
            pylab.savefig( fileName, bbox_inches='tight')
            pylab.close()

        
            
def getWavFileName( experimentFile, USVFile ):
    '''
    As the user may move the folder containing the datafile, this function provides the real location of wavfiles
    using the experimentFile and the USVFile. It provides the wav file located in the voc subfolder of the corresponding database.sqlite
    
    somefolder/voc/wavfile.wav
    somefolder/database.sqlite
    '''    
    wavFileName = os.path.dirname(os.path.abspath( experimentFile ) )
    wavFileName += "\\usv\\" + os.path.basename(USVFile)
    return wavFileName
    
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
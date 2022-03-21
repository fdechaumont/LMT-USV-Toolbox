'''
This code produces the test of voc trait usage differences per event context.
'''






import os
import sqlite3
from tkinter.filedialog import askopenfilename
import sqlite3
import os
import numpy as np
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import Event, EventTimeLine
from lmtanalysis.Measure import oneDay, oneHour, oneMinute

from lmtanalysis.Util import getMinTMaxTAndFileNameInput, convert_to_d_h_m_s
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
import pandas
import scipy.stats as stats
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
from numpy import sum

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from lmtanalysis.EventTimeLineCache import EventTimeLineCached
from lmtanalysis.Animal import *
import matplotlib.pyplot as plt
from lmtanalysis.Event import *
from lmtanalysis.Measure import *

import lmtanalysis
from tkinter.filedialog import askopenfilename
from tabulate import tabulate
from collections import Counter
import collections
import xlsxwriter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from scipy import signal
from scipy.io import wavfile
import os
import wave
import pylab

import pandas as pd


from LMT.USV2.figure.figParameter import getFigureVocTraits,\
    getFigureBehaviouralEvents, getFigureLabelTrait,\
    getFigureBehaviouralEventsLabels, getColorKO, getColorWT
from anaconda_project.internal import slugify
from scripts.Compute_Bad_Orientation_Estimation import isSameSign
from LMT.USV2.experimentList.experimentList import getAllExperimentList,\
    getExperimentList
from LMT.USV2.figure.featureHeatMap import featureHeatMap


def getVocCommonToEventTimeLine( eventTimeLineVoc, eventTimeLineEvent ):
    
    vocEventList = []
    eventDic = eventTimeLineEvent.getDictionary()
    
    for vocEvent in eventTimeLineVoc.eventList:
        if eventTimeLineEvent.overlap( vocEvent , eventDic ):
            vocEventList.append( vocEvent )
    '''
    # ELODIE TEST:
    for vocEvent in eventTimeLineVoc.eventList:
        for t in range( vocEvent.startFrame, vocEvent.endFrame ): 
            if t in eventDic:
                vocEventList.append( vocEvent )
    '''
    return vocEventList
    
def getTraitValueList( vocEventList , vocTrait ):

    valueList = []
    for vocEvent in vocEventList:
        valueList.append( vocEvent.metadata[vocTrait] )
        
    return valueList
        
def getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait ):
    vocList = getVocCommonToEventTimeLine(eventTimeLineVoc, eventTimeLine )
    traitData = getTraitValueList( vocList , vocTrait )
    return traitData



def getTraitData( eventTimeLineVoc , mainEventTimeLine , secondaryEventTimeLine, vocTrait ):
    
    
    # select voc common with mainTraitData:
    mainVocList = getVocCommonToEventTimeLine(eventTimeLineVoc, mainEventTimeLine )
    secondaryVocList = getVocCommonToEventTimeLine(eventTimeLineVoc, secondaryEventTimeLine )
    
    # extract data:
    mainTraitData = getTraitValueList( mainVocList , vocTrait )
    secondaryTraitData = getTraitValueList( secondaryVocList , vocTrait )
        
    return mainTraitData, secondaryTraitData

def getStarValue( p , correction ):
    
    stars = 0
    if p<0.05/correction:
        stars=1
    if p<0.01/correction:
        stars=2
    if p<0.001/correction:
        stars=3
    return stars

def getSign( mainTraitData, secondaryTraitData ):
    
    if np.mean( mainTraitData ) < np.mean( secondaryTraitData ):
        return 1
    return -1


def getJsonFile():
    return __file__[:-2]+"json"

def loadResult():
    
    print("Loading json...")
    try:
        with open( getJsonFile() ) as json_data:
            result = json.load(json_data)
        print("json file loaded.")
        return result
    except:
        print("Can't load result")
        result = {}
        return result


def loadResultFromJson(jsonFile):
    print("Loading json...")
    try:
        with open(jsonFile) as json_data:
            result = json.load(json_data)
        print("json file loaded.")
        return result
    except:
        print("Can't load result")
        result = {}
        return result


def saveResult( result ):

    with open( getJsonFile() , 'w') as f:
        json.dump(result, f, indent=4 )

'''
def getVocTraitsAndEvents():
    
    vocTraits = [ "durationMs","frequencyDynamicHz","startFrequencyHz",
                    "endFrequencyHz","diffStartEndFrequencyHz", 
                    "minFrequency","maxFrequency",
                   "meanFrequencyHz","frequencyTVHz", "linearityIndex",
                   "nbModulation","nbJump","meanPower","peakPower",
                   "slope","griffIndex"                  
                   ]
    
    eventListToTest = [ 
        "Stop isolated",
        "Break contact",
        "Approach rear",
        "Approach contact",
        "Oral-oral Contact",
        "Oral-genital Contact",
        "Side by side Contact",
        "Side by side Contact, opposite way", 
        "FollowZone Isolated",
        "Train2",
        "longChase",
        "Urinate USV" ]
    
    #vocTraits = sorted( vocTraits )
    #eventListToTest= sorted ( eventListToTest )
    
    return vocTraits, eventListToTest
'''
         
def setResult( result, expName, vocTrait, eventName, value ):
        
    if expName not in result:
        result[expName] = {}
    
    if vocTrait not in result[expName]:
        result[expName][vocTrait] = {}
            
    result[expName][vocTrait][eventName] = value
    
    
def computeVocTraitUsagePerEventContext( experiments ):
    
    results = loadResult()
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    vocTraits = getFigureVocTraits()
    #vocTraits = ["minFrequency"]
    eventListToTest = getFigureBehaviouralEvents( longList=True, withUrinate = False )
    
    for experiment in experiments:
        
        file = experiment.file
        print(file)
        
        expName = experiment.getFullName()
        
        print ("Current file: " , file )
        chrono = Chronometer("File time")
        connection = sqlite3.connect( file )        
        print("Load all vocs...")
        eventTimeLineVoc = EventTimeLine( connection, "Voc", loadEventIndependently = True )

        print( "Prepare and preload data...")
        allData = {}
        for vocTrait in vocTraits:
            allData[vocTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + vocTrait + " / " + eventName )
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[vocTrait][eventName] = getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait )
                setResult(results, expName, vocTrait, eventName , allData[vocTrait][eventName] )
        
        
        '''
        # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.
        
        axIndex = 0
        nbrows=4
        fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
        
        # flatten axes to ease ax indexing
        print( axes )
        axes2 = []
        for value in axes:
            axes2.extend( value )
        axes = axes2    
        print( "axes2", axes2 )
        
        # let's go
        plt.title( file )
        
        # bonfferoni correction 
        correction = (len(eventListToTest)*(len(eventListToTest)-1)/2) # number of eventTested
        
        
        for vocTrait in vocTraits:
            
            print ("Current trait: " , vocTrait )
                            
            textsToPlot = []
            
            data = {}
            
            plt.sca( axes[axIndex ])
            plt.gca().set_title( vocTrait )
    
            # kruskal test
            kruskalList = []
            for eventName in eventListToTest:
            
                traitData = allData[vocTrait][eventName]
                print("Event: ", eventName)
                #print(list(traitData))
                print("length: ", len(list(traitData)))
                kruskalList.append( list( traitData ) )
                print("---------------")
            
            
            kruskalTest = stats.kruskal( *kruskalList )
            kruskalTestOk = True
            if kruskalTest.pvalue > 0.05/ (len(vocTraits)):
            
                kruskalTestOk = False
                print("KRUSKAL FAILED FOR VOCTRAIT ", vocTrait )
               
            if not kruskalTestOk:
                for mainEventName in eventListToTest:
                    data[ mainEventName ] = []                    
                    for secondaryEventName in eventListToTest:
                        data[ mainEventName ].append( 0 )
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "K" ) )
                
            if kruskalTestOk:
                for mainEventName in eventListToTest:
                    
                    
                    #mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
                    data[ mainEventName ] = []
                    
                    kruskalList = []
                    
                    print("Computing ", mainEventName )
                    
                    # pour un contexte:
                    # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
                    
                    
                    for secondaryEventName in eventListToTest:
                                        
                        mainTraitData = allData[vocTrait][mainEventName]
                        secondaryTraitData = allData[vocTrait][secondaryEventName] 
                                                
                        stat = 0
                        p = 0
                        error = False
                        try:
                            #debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                            stat, p = mannwhitneyu( mainTraitData, secondaryTraitData )
                        except:
                            error= True
                            print("Mann Whitney Error.")
                            x = eventListToTest.index( mainEventName )
                            y = eventListToTest.index( secondaryEventName )
                            textsToPlot.append( ( x , y , "MWE" ) )
                        
                        sign = getSign( mainTraitData, secondaryTraitData )
                        
                        starValue = getStarValue( p , correction )*sign
                        data[ mainEventName ].append( starValue )
                        
                        
            df = pd.DataFrame( data=data )
            
            ax = sns.heatmap( df , linewidths=.5, yticklabels=eventListToTest , cmap= sns.color_palette("coolwarm", 7) )
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            
            for text in textsToPlot:
                plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=6,horizontalalignment='center', verticalalignment='center' )
                   
            axIndex += 1 
            
        #plt.tight_layout()
    
        plt.savefig( "vocTraitUsagePerEventContext.pdf" , dpi=100 )
        #plt.show()
        #plt.text( 0 , 0 , file )
        #debugOut.close()
        print( "Time for file (s): " , chrono.getTimeInS() )
        '''
        connection.close()
    
    saveResult( results )



def plotVocTraitUsagePerEventContextPerFile( experiments ):

    results = loadResult()
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    vocTraits = getFigureVocTraits()
    eventListToTest = getFigureBehaviouralEvents( longList=True, withUrinate = False )

    sns.set(font_scale=1.5)
    
    for experiment in experiments:
    
        file = experiment.file
        print(file)
        
        expName = experiment.getFullName()
        
        '''
        print ("Current file: " , file )
        chrono = Chronometer("File time")
        connection = sqlite3.connect( file )        
        print("Load all vocs...")
        eventTimeLineVoc = EventTimeLine( connection, "Voc", loadEventIndependently = True )
        cleanVoc( eventTimeLineVoc )
        '''
        
        print( "Prepare and preload data...")
        allData = {}
        for vocTrait in vocTraits:
            allData[vocTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + vocTrait + " / " + eventName )
                '''
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[vocTrait][eventName] = getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait )
                setResult(results, expName, vocTrait, eventName , allData[vocTrait][eventName] )
                '''
        
                allData[vocTrait][eventName] = results[expName][vocTrait][eventName] 
        
        
        # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.
        
        axIndex = 0
        nbrows=4
        fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
        
        # flatten axes to ease ax indexing
        print( axes )
        axes2 = []
        for value in axes:
            axes2.extend( value )
        axes = axes2    
        print( "axes2", axes2 )
        
        # let's go
        plt.title( file )
        
        # bonfferoni correction 
        # correction = (len(eventListToTest)*(len(eventListToTest)-1)/2) # number of eventTested
        correction = getBonfferoniCorrection( len( eventListToTest ) , 3, len( vocTraits ) ) # 3 for 5we, 3mo, 7mo
        
        for vocTrait in vocTraits:
            
            print ("Current trait: " , vocTrait )
                            
            textsToPlot = []
            
            data = {}
            
            plt.sca( axes[axIndex ])
            plt.gca().set_title( vocTrait )
    
            # kruskal test
            kruskalList = []
            for eventName in eventListToTest:
            
                traitData = allData[vocTrait][eventName]
                print("Event: ", eventName)
                #print(list(traitData))
                print("length: ", len(list(traitData)))
                kruskalList.append( list( traitData ) )
                print("---------------")
            
            
            kruskalTest = stats.kruskal( *kruskalList )
            kruskalTestOk = True
            if kruskalTest.pvalue > 0.05/ (len(vocTraits)):
            
                kruskalTestOk = False
                print("KRUSKAL FAILED FOR VOCTRAIT ", vocTrait )
               
            if not kruskalTestOk:
                for mainEventName in eventListToTest:
                    data[ mainEventName ] = []                    
                    for secondaryEventName in eventListToTest:
                        data[ mainEventName ].append( 0 )
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "K" ) )
                
            if kruskalTestOk:
                for mainEventName in eventListToTest:
                    
                    
                    #mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
                    data[ mainEventName ] = []
                    
                    kruskalList = []
                    
                    print("Computing ", mainEventName )
                    
                    # pour un contexte:
                    # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
                    
                    
                    for secondaryEventName in eventListToTest:
                                        
                        mainTraitData = allData[vocTrait][mainEventName]
                        secondaryTraitData = allData[vocTrait][secondaryEventName] 
                                                
                        stat = 0
                        p = 0
                        error = False
                        try:
                            #debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                            stat, p = mannwhitneyu( mainTraitData, secondaryTraitData )
                        except:
                            error= True
                            print("Mann Whitney Error.")
                            x = eventListToTest.index( mainEventName )
                            y = eventListToTest.index( secondaryEventName )
                            textsToPlot.append( ( x , y , "MWE" ) )
                        
                        sign = getSign( mainTraitData, secondaryTraitData )
                        
                        starValue = getStarValue( p , correction )*sign
                        data[ mainEventName ].append( starValue )
                        
                        
            df = pd.DataFrame( data=data )
            
            showColorBar = False
            if (axIndex-3) % nbrows == 0:    
                showColorBar = True  
            
            sns.heatmap( df , linewidths=.5, yticklabels=eventListToTest , xticklabels=eventListToTest, cmap= sns.color_palette("coolwarm", 7) ,
                               vmin = -3, vmax= 3, cbar = showColorBar )
            
            
            plt.xticks(rotation=90, fontsize=18 )
            plt.yticks(rotation=0, fontsize=18 )
            
            for text in textsToPlot:
                plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=6,horizontalalignment='center', verticalalignment='center' )
                   
            axIndex += 1 
            
        plt.tight_layout()
                
        plt.savefig( "vocTraitUsagePerEventContext_"+slugify( expName) +".pdf" , dpi=100 )
        #plt.show()
        #plt.text( 0 , 0 , file )
        #debugOut.close()
        #print( "Time for file (s): " , chrono.getTimeInS() )

def getBonfferoniCorrection( numberOfEventTested , numberOfAgeClassTested, numberOfVocTraitTested ):
    
    correction = ( numberOfEventTested* (numberOfEventTested-1) ) / 2 # number of eventTested
    correction*= numberOfAgeClassTested
    correction*= numberOfVocTraitTested
    return correction

def plotVocTraitUsagePerEventContextPerSet( experiments, genotype ):

    results = loadResult()

    if "allSet" not in results:
        results["allSet"] = {}
    results["allSet"][genotype] = {}

    #vocTraits, eventListToTest = getVocTraitsAndEvents()
    vocTraits = getFigureVocTraits()
    vocTraitsLabels = []
    for trait in vocTraits:
        vocTraitsLabels.append(getFigureLabelTrait(trait))

    eventListToTest = getFigureBehaviouralEvents( longList=False, withUrinate = False )
    eventLabels = []
    for event in eventListToTest:
        eventLabels.append(getFigureBehaviouralEventsLabels(event))


    print( "Prepare and preload data...")
    allData = {}

    sns.set(font_scale=1.5)

    for experiment in experiments:

        file = experiment.file
        print(file)

        expName = experiment.getFullName()
        '''
        print ("Current file: " , file )
        chrono = Chronometer("File time")
        connection = sqlite3.connect( file )        
        print("Load all vocs...")
        eventTimeLineVoc = EventTimeLine( connection, "Voc", loadEventIndependently = True )
        cleanVoc( eventTimeLineVoc )
        '''
        for vocTrait in vocTraits:
            if vocTrait not in allData:
                allData[vocTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + vocTrait + " / " + eventName )
                '''
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[vocTrait][eventName] = getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait )
                setResult(results, expName, vocTrait, eventName , allData[vocTrait][eventName] )
                '''
                if eventName not in allData[vocTrait]:
                    allData[vocTrait][eventName] = []
                allData[vocTrait][eventName].extend( results[expName][vocTrait][eventName] )


    # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.

    axIndex = 0
    nbrows=4
    fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)

    # flatten axes to ease ax indexing
    print( axes )
    axes2 = []
    for value in axes:
        axes2.extend( value )
    axes = axes2
    print( "axes2", axes2 )

    # let's go
    plt.title( file )

    # bonfferoni correction
    ''' 
    correction = (len(eventListToTest)*(len(eventListToTest)-1)/2) # number of eventTested
    correction*= 3 # number of age class
    correction*= len( vocTraits )
    '''
    correction = getBonfferoniCorrection( len( eventListToTest ) , 3, len( vocTraits ) ) # 3 for 5we, 3mo, 7mo

    for vocTrait in vocTraits:

        print ("Current trait: " , vocTrait )

        textsToPlot = []

        data = {}

        plt.sca( axes[axIndex ])
        plt.gca().set_title( getFigureLabelTrait(vocTrait) )

        results["allSet"][genotype][vocTrait] = {}

        # kruskal test
        kruskalList = []
        for eventName in eventListToTest:

            traitData = allData[vocTrait][eventName]
            print("Event: ", eventName)
            #print(list(traitData))
            print("length: ", len(list(traitData)))
            kruskalList.append( list( traitData ) )
            print("---------------")


        kruskalTest = stats.kruskal( *kruskalList )
        kruskalTestOk = True
        if kruskalTest.pvalue > 0.05/ (len(vocTraits)):

            kruskalTestOk = False
            print("KRUSKAL FAILED FOR VOCTRAIT ", vocTrait )

        if not kruskalTestOk:
            for mainEventName in eventListToTest:
                data[ mainEventName ] = []
                results["allSet"][genotype][vocTrait][ mainEventName ] = []
                for secondaryEventName in eventListToTest:
                    data[ mainEventName ].append( 0 )
                    results["allSet"][genotype][vocTrait][ mainEventName ].append( 0 )
                    x = eventListToTest.index( mainEventName )
                    y = eventListToTest.index( secondaryEventName )
                    textsToPlot.append( ( x , y , "K" ) )

        if kruskalTestOk:
            for mainEventName in eventListToTest:


                #mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
                data[ mainEventName ] = []
                results["allSet"][genotype][vocTrait][ mainEventName ] = []

                kruskalList = []

                print("Computing ", mainEventName )

                # pour un contexte:
                # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes


                for secondaryEventName in eventListToTest:

                    mainTraitData = allData[vocTrait][mainEventName]
                    secondaryTraitData = allData[vocTrait][secondaryEventName]

                    stat = 0
                    p = 0
                    error = False
                    try:
                        #debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                        stat, p = mannwhitneyu( mainTraitData, secondaryTraitData )
                    except:
                        error= True
                        print("Mann Whitney Error.")
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "MWE" ) )

                    sign = getSign( mainTraitData, secondaryTraitData )

                    starValue = getStarValue( p , correction )*sign
                    data[ mainEventName ].append( starValue )
                    results["allSet"][genotype][vocTrait][ mainEventName ].append( starValue )


        df = pd.DataFrame( data=data )

        showColorBar = False
        if (axIndex-3) % nbrows == 0:
            showColorBar = True

        g = sns.heatmap( df , linewidths=.5, yticklabels=eventLabels , xticklabels=eventLabels, cmap= sns.color_palette("coolwarm", 7) ,
                           vmin = -3, vmax= 3, cbar = showColorBar )




        plt.xticks(rotation=90, fontsize=18 )
        plt.yticks(rotation=0, fontsize=18 )

        for text in textsToPlot:
            plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=6,horizontalalignment='center', verticalalignment='center' )

        axIndex += 1

    plt.tight_layout()

    name = ""
    for experiment in experiments:
        name+= slugify ( experiment.getFullName() ) + "__"
    plt.savefig( "vocTraitUsagePerEventContext_AllSet"+name[:20]+".pdf" , dpi=100 )
    plt.savefig( "vocTraitUsagePerEventContext_AllSet"+name[:20]+".svg" , dpi=100 )
    #plt.show()
    #plt.text( 0 , 0 , file )
    #debugOut.close()
    #print( "Time for file (s): " , chrono.getTimeInS() )

    saveResult(results)
    print( results["allSet"] )


def plotVocTraitUsagePerEventContextPerSetShort(jsonFile, experiments, genotype, eventListToTest, variable, ax):
    results = loadResultFromJson(jsonFile=jsonFile)

    if "allSet" not in results:
        results["allSet"] = {}
    results["allSet"][genotype] = {}

    vocTrait = variable

    eventLabels = []
    for event in eventListToTest:
        eventLabels.append(getFigureBehaviouralEventsLabels(event))

    print("Prepare and preload data...")
    allData = {}

    sns.set(font_scale=1.5)

    for experiment in experiments:

        file = experiment.file
        print(file)
        expName = experiment.getFullName()

        if vocTrait not in allData:
            allData[vocTrait] = {}
        for eventName in eventListToTest:
            print("PRELOAD " + vocTrait + " / " + eventName)

            if eventName not in allData[vocTrait]:
                allData[vocTrait][eventName] = []
            allData[vocTrait][eventName].extend(results[expName][vocTrait][eventName])

            # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.

    # bonferroni correction
    correction = getBonfferoniCorrection(len(getFigureBehaviouralEvents(longList=False, withUrinate=False)), 3, len(getFigureVocTraits()))  # 3 for 5we, 3mo, 7mo

    print("Current trait: ", vocTrait)

    textsToPlot = []

    data = {}

    #ax.set_title(getFigureLabelTrait(vocTrait), fontsize=12)

    results["allSet"][genotype][vocTrait] = {}

    # kruskal test
    kruskalList = []
    for eventName in eventListToTest:
        traitData = allData[vocTrait][eventName]
        print("Event: ", eventName)
        # print(list(traitData))
        print("length: ", len(list(traitData)))
        kruskalList.append(list(traitData))
        print("---------------")

    kruskalTest = stats.kruskal(*kruskalList)
    kruskalTestOk = True
    if kruskalTest.pvalue > 0.05 / (len(getFigureVocTraits())):
        kruskalTestOk = False
        print("KRUSKAL FAILED FOR VOCTRAIT ", vocTrait)

    if not kruskalTestOk:
        for mainEventName in eventListToTest:
            data[mainEventName] = []
            results["allSet"][genotype][vocTrait][mainEventName] = []
            for secondaryEventName in eventListToTest:
                data[mainEventName].append(0)
                results["allSet"][genotype][vocTrait][mainEventName].append(0)
                x = eventListToTest.index(mainEventName)
                y = eventListToTest.index(secondaryEventName)
                textsToPlot.append((x, y, "K"))

    if kruskalTestOk:
        for mainEventName in eventListToTest:

            # mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
            data[mainEventName] = []
            results["allSet"][genotype][vocTrait][mainEventName] = []

            kruskalList = []

            print("Computing ", mainEventName)

            # pour un contexte:
            # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
            for secondaryEventName in eventListToTest:

                mainTraitData = allData[vocTrait][mainEventName]
                secondaryTraitData = allData[vocTrait][secondaryEventName]

                stat = 0
                p = 0
                error = False
                try:
                    # debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                    stat, p = mannwhitneyu(mainTraitData, secondaryTraitData)

                except:
                    error = True
                    print("Mann Whitney Error.")
                    x = eventListToTest.index(mainEventName)
                    y = eventListToTest.index(secondaryEventName)
                    textsToPlot.append((x, y, "MWE"))

                sign = getSign(mainTraitData, secondaryTraitData)

                starValue = getStarValue(p, correction) * sign
                data[mainEventName].append(starValue)
                results["allSet"][genotype][vocTrait][mainEventName].append(starValue)

    df = pd.DataFrame(data=data)

    showColorBar = False

    g = sns.heatmap(df, linewidths=.5, yticklabels=False, xticklabels=False,
                    cmap=sns.color_palette("coolwarm", 7),
                    vmin=-3, vmax=3, cbar=showColorBar, ax=ax)
    g.set_ylim(9,-1.5)
    ax.set_facecolor('white')

    #g.set_xticklabels(g.get_xticklabels(), rotation=45, fontsize=10, horizontalalignment='right')
    #g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=10)

    #plt.xticks(rotation=45, fontsize=12)
    #plt.yticks(rotation=0, fontsize=12)

    for text in textsToPlot:
        ax.text(text[0] + 0.5, text[1] + 0.5, text[2], fontsize=6, horizontalalignment='center', verticalalignment='center')

    saveResult(results)
    print(results["allSet"])
    return g


def plotVocTraitUsagePerEventContextPerSetShort2(jsonFile, experiments, genotype, eventListToTest, variable, ax):
    results = loadResultFromJson(jsonFile=jsonFile)

    if "allSet" not in results:
        results["allSet"] = {}
    results["allSet"][genotype] = {}

    vocTrait = variable

    eventLabels = []
    for event in eventListToTest:
        eventLabels.append(getFigureBehaviouralEventsLabels(event))

    print("Prepare and preload data...")
    allData = {}

    for experiment in experiments:

        file = experiment.file
        print(file)
        expName = experiment.getFullName()

        if vocTrait not in allData:
            allData[vocTrait] = {}
        for eventName in eventListToTest:
            print("PRELOAD " + vocTrait + " / " + eventName)

            if eventName not in allData[vocTrait]:
                allData[vocTrait][eventName] = []
            allData[vocTrait][eventName].extend(results[expName][vocTrait][eventName])

            # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.

    # bonferroni correction
    correction = getBonfferoniCorrection(len(getFigureBehaviouralEvents(longList=False, withUrinate=False)), 3, len(getFigureVocTraits()))  # 3 for 5we, 3mo, 7mo

    print("Current trait: ", vocTrait)

    textsToPlot = []

    dataPValue = {}
    dataSizeEffect = {}

    #ax.set_title(getFigureLabelTrait(vocTrait), fontsize=12)

    results["allSet"][genotype][vocTrait] = {}

    # kruskal test
    kruskalList = []
    for eventName in eventListToTest:
        traitData = allData[vocTrait][eventName]
        print("Event: ", eventName)
        # print(list(traitData))
        print("length: ", len(list(traitData)))
        kruskalList.append(list(traitData))
        print("---------------")

    kruskalTest = stats.kruskal(*kruskalList)
    kruskalTestOk = True
    if kruskalTest.pvalue > 0.05 / (len(getFigureVocTraits())):
        kruskalTestOk = False
        print("KRUSKAL FAILED FOR VOCTRAIT ", vocTrait)

    if not kruskalTestOk:
        for mainEventName in eventListToTest:
            dataPValue[mainEventName] = []
            dataSizeEffect[mainEventName] = []
            
            results["allSet"][genotype][vocTrait][mainEventName] = []
            for secondaryEventName in eventListToTest:
                dataPValue[mainEventName].append(10000)
                dataSizeEffect[mainEventName].append(0)
                
                results["allSet"][genotype][vocTrait][mainEventName].append(10000)
                '''
                x = eventListToTest.index(mainEventName)
                y = eventListToTest.index(secondaryEventName)
                textsToPlot.append((x, y, "K"))
                '''

    if kruskalTestOk:
        for mainEventName in eventListToTest:

            # mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
            dataPValue[mainEventName] = []
            dataSizeEffect[mainEventName] = []
            results["allSet"][genotype][vocTrait][mainEventName] = []

            kruskalList = []

            print("Computing ", mainEventName)

            # pour un contexte:
            # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
            for secondaryEventName in eventListToTest:

                mainTraitData = allData[vocTrait][mainEventName]
                secondaryTraitData = allData[vocTrait][secondaryEventName]

                stat = 0
                p = 0
                error = False
                try:
                    # debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                    stat, p = mannwhitneyu(mainTraitData, secondaryTraitData)
                except:
                    error = True
                    print("Mann Whitney Error.")
                    
                    x = eventListToTest.index(mainEventName)
                    y = eventListToTest.index(secondaryEventName)
                    textsToPlot.append((x, y, "MWE"))
                    

                sign = getSign(mainTraitData, secondaryTraitData)
                
                # compute sizeEffectNormalized
                meanMain = np.mean( mainTraitData )
                meanSecondary = np.mean( secondaryTraitData )
                valList= []
                for v in mainTraitData:
                    dif = (v-meanMain)**2
                    valList.append( dif ) 
                for v in secondaryTraitData:
                    dif = (v-meanSecondary)**2
                    valList.append( dif ) 
                stdCombined = math.sqrt( np.mean( valList ) )
                sizeEffect = ( meanSecondary-meanMain ) / stdCombined
                 
                pValueCorrected = p*correction
                dataPValue[mainEventName].append( pValueCorrected )
                dataSizeEffect[mainEventName].append( sizeEffect )
                results["allSet"][genotype][vocTrait][mainEventName].append( pValueCorrected )

    dfPValue = pd.DataFrame(data=dataPValue , index = eventListToTest )
    dfSizeEffect = pd.DataFrame(data=dataSizeEffect , index = eventListToTest  )
    g = featureHeatMap( dfSizeEffect, dfPValue , ax =ax , title = variable )
    
    '''
    showColorBar = False

    g = sns.heatmap(df, linewidths=.5, yticklabels=False, xticklabels=False,
                    cmap=sns.color_palette("coolwarm", 7),
                    vmin=-3, vmax=3, cbar=showColorBar, ax=ax)
    '''
    #g.set_xticklabels(g.get_xticklabels(), rotation=45, fontsize=10, horizontalalignment='right')
    #g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=10)

    #plt.xticks(rotation=45, fontsize=12)
    #plt.yticks(rotation=0, fontsize=12)
    
    
    for text in textsToPlot:
        ax.text(text[0] + 0.5, text[1] + 0.5, text[2], fontsize=6, horizontalalignment='center', verticalalignment='center')
    
    #saveResult(results)
    #print(results["allSet"])
    return ax
   

def plotVocTraitUsagePerEventContextPerSetCommon():

    results = loadResult()
    
    #results["allSet"][genotype][vocTrait][ mainEventName ].append( starValue )

    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    vocTraits = getFigureVocTraits()
    eventListToTest = getFigureBehaviouralEvents( withUrinate = True )
    
    print( "Prepare and preload data...")
    allData = {}
    
    sns.set(font_scale=1.5)

    axIndex = 0
    nbrows=4
    fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
    
    # flatten axes to ease ax indexing
    print( axes )
    axes2 = []
    for value in axes:
        axes2.extend( value )
    axes = axes2    
    print( "axes2", axes2 )
    
    
    for vocTrait in vocTraits:
        
        print ("Current trait: " , vocTrait )
                        
        textsToPlot = []
        
        data = {}
        
        plt.sca( axes[axIndex ])
        plt.gca().set_title( vocTrait )
        
        if True:
            
            for mainEventName in eventListToTest:
                
                data[ mainEventName ] = []
                #results["allSet"][genotype][vocTrait][ mainEventName ] = []
                
                kruskalList = []
                
                print("Computing ", mainEventName )
                                
                for secondaryEventName in eventListToTest:
                                    
                    resultWT = results["allSet"]["WT"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]
                    resultKO = results["allSet"]["KO"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]

                    starValue = 0
                    
                    if resultWT!=0 and resultKO!=0 and isSameSign( resultWT , resultKO ):
                        if resultWT > 0:
                            starValue = 1
                        else:
                            starValue = -1
                    
                    data[ mainEventName ].append( starValue )

                    
                    
        df = pd.DataFrame( data=data )

        showColorBar = False
        if (axIndex-3) % nbrows == 0:    
            showColorBar = True  
        
        sns.heatmap( df , linewidths=.5, yticklabels=eventListToTest , xticklabels=eventListToTest, cmap= sns.color_palette("coolwarm", 7) , 
                           vmin = -1, vmax= 1, cbar = showColorBar )
        
        plt.xticks(rotation=90, fontsize=18 )
        plt.yticks(rotation=0, fontsize=18 )
        
        for text in textsToPlot:
            plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=6,horizontalalignment='center', verticalalignment='center' )
               
        axIndex += 1 
        
    plt.tight_layout()
    
    plt.savefig( "vocTraitUsagePerEventContext_AllSet_Common.pdf" , dpi=100 )


def plotVocTraitUsagePerEventContextPerSetDifference():

    results = loadResult()

    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    vocTraits = getFigureVocTraits()
    eventListToTest = getFigureBehaviouralEvents( longList=False, withUrinate = False )
    
    print( "Prepare and preload data...")
    allData = {}
    
    sns.set(font_scale=1.5)

    axIndex = 0
    nbrows=4
    fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
    
    # flatten axes to ease ax indexing
    print( axes )
    axes2 = []
    for value in axes:
        axes2.extend( value )
    axes = axes2    
    print( "axes2", axes2 )
    
    
    for vocTrait in vocTraits:
        
        print ("Current trait: " , vocTrait )
                        
        textsToPlot = []
        
        data = {}
        
        plt.sca( axes[axIndex ])
        plt.gca().set_title( vocTrait )
        
        if True:
            
            for mainEventName in eventListToTest:
                
                data[ mainEventName ] = []
                
                print("Computing ", mainEventName )
                                
                for secondaryEventName in eventListToTest:
                                    
                    resultWT = results["allSet"]["WT"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]
                    resultKO = results["allSet"]["KO"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]

                    starValue = 0
                    
                    '''
                    color code:
                        0 - gray = no difference
                        1 - black = opposed
                    
                    '''
                    x = eventListToTest.index( mainEventName )
                    y = eventListToTest.index( secondaryEventName )
                    
                
                    #print( "-- ")
                    #print( resultWT , resultKO )
                    if resultWT!=0 and resultKO!=0:                        
                        if isSameSign ( resultWT , resultKO ):
                            starValue = 0 # same 
                        else:
                            starValue = 1 # opposed
                                        
                    if bool( resultWT ) ^ bool( resultKO ):
                        #print( "xor passed" )
                        if resultWT == 0:
                            if resultKO > 0:
                                starValue = 2 # KO Positive
                                textsToPlot.append( ( x , y , "+" ) )
                            else:
                                starValue = 2 # KO Negative
                                textsToPlot.append( ( x , y , "-" ) )
                        else:
                            if resultWT > 0:
                                starValue = 3 # WT Positive
                                textsToPlot.append( ( x , y , "+" ) )
                            else:
                                starValue = 3 # WT Negative
                                textsToPlot.append( ( x , y , "-" ) )
                    
                    #print( starValue )        
                     
                    
                    data[ mainEventName ].append( starValue )
                    #input("wait")
                    
                    
        df = pd.DataFrame( data=data )

        showColorBar = False
        '''
        if (axIndex-3) % nbrows == 0:    
            showColorBar = True
        '''  
        
        #flatui = ["#dddddd" , "#34495e" , "#ffff55", "#dddd55", "#55ff55", "#55dd55" ]
        flatui = ["#dddddd" , "#34495e" , getColorKO(), getColorWT() ]
        sns.heatmap( df , linewidths=.5, yticklabels=eventListToTest , xticklabels=eventListToTest, cmap= sns.color_palette(flatui) , 
                           vmin = 0, vmax= len( flatui ), cbar = showColorBar ) # , annot = True )
        
        plt.xticks(rotation=90, fontsize=18 )
        plt.yticks(rotation=0, fontsize=18 )
        
        for text in textsToPlot:
            plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=30,horizontalalignment='center', verticalalignment='center' )
               
        axIndex += 1 
        
    plt.tight_layout()
    
    plt.savefig( "vocTraitUsagePerEventContext_AllSet_Difference.pdf" , dpi=100 )


def plotVocTraitUsagePerEventContextCommonAndDifference(  ):

    results = loadResult()
    
    vocTraits = getFigureVocTraits()
    eventListToTest = getFigureBehaviouralEvents( longList=False, withUrinate = False )
    
    print( "Prepare and preload data...")
    allData = {}
    
    sns.set(font_scale=1.5)

    axIndex = 0
    nbrows=4
    fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(vocTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
    
    # flatten axes to ease ax indexing
    print( axes )
    axes2 = []
    for value in axes:
        axes2.extend( value )
    axes = axes2    
    print( "axes2", axes2 )
    
    
    for vocTrait in vocTraits:
        
        print ("Current trait: " , vocTrait )
                        
        textsToPlot = []
        textsToPlotOpposite = []
        
        data = {}
        
        plt.sca( axes[axIndex ])
        plt.gca().set_title( vocTrait )
        
        if True:
            
            for mainEventName in eventListToTest:
                
                data[ mainEventName ] = []
                #results["allSet"][genotype][vocTrait][ mainEventName ] = []
                
                print("Computing ", mainEventName )
                                
                for secondaryEventName in eventListToTest:
                                    
                    resultWT = results["allSet"]["WT"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]
                    resultKO = results["allSet"]["KO"][vocTrait][ mainEventName ][ eventListToTest.index( secondaryEventName )]

                    starValue = 0
                    
                    if resultWT!=0 and resultKO!=0 and isSameSign ( resultWT , resultKO ):
                        if resultWT > 0:
                            starValue = 4
                        else:
                            starValue = 5
                    else:
                        
                        '''
                        color code:
                            0 - gray = no difference
                            1 - black = opposed
                        
                        '''
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        
                    
                        #print( "-- ")
                        #print( resultWT , resultKO )
                        if resultWT!=0 and resultKO!=0:                        
                            if isSameSign ( resultWT , resultKO ):
                                starValue = 0 # same 
                            else:
                                
                                starValue = 1 # opposed
                                if resultWT > resultKO:
                                    textsToPlotOpposite.append( ( x , y , "WT" ) )
                                else:
                                    textsToPlotOpposite.append( ( x , y , "KO" ) )
                                
                                            
                        if bool( resultWT ) ^ bool( resultKO ):
                            #print( "xor passed" )
                            if resultWT == 0:
                                if resultKO > 0:
                                    starValue = 2 # KO Positive
                                    textsToPlot.append( ( x , y , "+" ) )
                                else:
                                    starValue = 2 # KO Negative
                                    textsToPlot.append( ( x , y , "-" ) )
                            else:
                                if resultWT > 0:
                                    starValue = 3 # WT Positive
                                    textsToPlot.append( ( x , y , "+" ) )
                                else:
                                    starValue = 3 # WT Negative
                                    textsToPlot.append( ( x , y , "-" ) )
                        
                        #print( starValue )        
                         
                    
                    data[ mainEventName ].append( starValue )
                    #input("wait")
                    
                    
        df = pd.DataFrame( data=data )

        showColorBar = False
        '''
        if (axIndex-3) % nbrows == 0:    
            showColorBar = True
        '''  
        
        #flatui = ["#dddddd" , "#34495e" , "#ffff55", "#dddd55", "#55ff55", "#55dd55" ]
        #flatui = ["#dddddd" , "#34495e" , getColorKO(), getColorWT() , "#DC5D4A" , "#617BBC" ]
        flatui = ["#dddddd" , "#34495e" , getColorKO(), getColorWT() , "#E29991" , "#839FDB" ]
        
        sns.heatmap( df , linewidths=.5, yticklabels=eventListToTest , xticklabels=eventListToTest, cmap= sns.color_palette(flatui) , 
                           vmin = 0, vmax= len( flatui ), cbar = showColorBar ) # , annot = True )
        
        plt.xticks(rotation=90, fontsize=18 )
        plt.yticks(rotation=0, fontsize=18 )
        
        for text in textsToPlot:
            plt.text( text[0]+0.5, text[1]+0.5, text[2], fontsize=30,horizontalalignment='center', verticalalignment='center' )

        for text in textsToPlotOpposite:
            plt.text( text[0]+0.5, text[1]+0.5, text[2], color="white", fontsize=13,horizontalalignment='center', verticalalignment='center' )

               
        axIndex += 1
        plt.gca().set_ylim(12,0)
        
    plt.tight_layout()
    
    plt.savefig( "vocTraitUsagePerEventContext_AllSet_CommonAndDif.pdf" , dpi=100 )
    plt.savefig( "vocTraitUsagePerEventContext_AllSet_CommonAndDif.svg" , dpi=100 )


if __name__ == '__main__':
    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})

    while True:
        question = "Do you want to:"
        question += "\n\t [c]ompute all?"
        question += "\n\t [p]lot heatmaps of acoustic variations between contexts?"
        question += "\n\t [pd]lot heatmaps of differences between two sets of data?"
        question += "\n\t [bp]plot boxplots of a subset of acoustic traits between a subset of contexts?"
        question += "\n"
        answer = input(question)

        if answer == "c":
            # recompute all
            #experiments = getExperimentList( age="3mo", sex="female" , genotype="WT" )
            experiments = getAllExperimentList()
            computeVocTraitUsagePerEventContext(experiments)
            print('Computation done.')
            break

        if answer == "p":
            experiments = getExperimentList(age="3mo", sex="female", genotype="WT")
            plotVocTraitUsagePerEventContextPerSet(experiments, "WT")
            print("All done")
            break

        if answer == "pd":
            '''Check if only B6 females of 3 mo are taken into account!!!'''
            plotVocTraitUsagePerEventContextCommonAndDifference()

            print("All done")
            break









    
    '''
    plotVocTraitUsagePerEventContextPerSetDifference()    
    '''
    

    
    print("All done")


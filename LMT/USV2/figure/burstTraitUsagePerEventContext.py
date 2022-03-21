'''
This code produces the test of voc trait usage differences per event context.
'''



#from experimental.voc.analysis.migrateVoc import vocUtil



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
from scripts.ComputeMeasuresIdentityProfileOneMouseAutomatic import *

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



from scripts.Compute_Bad_Orientation_Estimation import isSameSign



from LMT.USV2.figure.figParameter import getFigureBurstTraits,\
    getFigureBehaviouralEvents, getFigureLabelBurstTraits,\
    getFigureBehaviouralEventsLabels, getColorWT, getColorKO, getColorAge
from LMT.USV2.lib.burster import createBurstFromVoc
from anaconda_project.internal import slugify
from LMT.USV2.figure.featureHeatMap import featureHeatMap
from LMT.USV2.experimentList.experimentList import getExperimentList




def getBurstCommonToEventTimeLine( burstList, eventTimeLine ):
    
    burstListCommon = []
    eventDic = eventTimeLine.getDictionary()
    
    for burst in burstList:
        if eventTimeLine.overlap( burst.getAsEvent() , eventDic ):
            burstListCommon.append( burst )
    
    return burstListCommon
    
def getBurstTraitValueList( burstList , burstTrait ):

    valueList = []
    for burstEvent in burstList:
        valueList.append( burstEvent.getValue(burstTrait) )
        
    return valueList
        
def getBurstTraitDataSolo( burstList , eventTimeLine , burstTrait ):
    burstListCommon = getBurstCommonToEventTimeLine( burstList, eventTimeLine )
    traitData = getBurstTraitValueList( burstListCommon , burstTrait )
    return traitData


'''
def getTraitData( eventTimeLineVoc , mainEventTimeLine , secondaryEventTimeLine, vocTrait ):
    
    
    # select voc common with mainTraitData:
    mainVocList = getVocCommonToEventTimeLine(eventTimeLineVoc, mainEventTimeLine )
    secondaryVocList = getVocCommonToEventTimeLine(eventTimeLineVoc, secondaryEventTimeLine )
    
    # extract data:
    mainTraitData = getTraitValueList( mainVocList , vocTrait )
    secondaryTraitData = getTraitValueList( secondaryVocList , vocTrait )
        
    return mainTraitData, secondaryTraitData
'''

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
    
    if np.nanmean( np.array(mainTraitData, dtype=np.float64) ) < np.nanmean( np.array(secondaryTraitData, dtype=np.float64) ):
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
        print(result[list(result.keys())[0]].keys())
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
         
def setResult( result, expName, burstTrait, eventName, value ):
        
    if expName not in result:
        result[expName] = {}
    
    if burstTrait not in result[expName]:
        result[expName][burstTrait] = {}
            
    result[expName][burstTrait][eventName] = value
    
    
def computeBurstTraitUsagePerEventContext( experiments ):
    
    results = loadResult()
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    #vocTraits = getFigureVocTraits()
    burstTraits = getFigureBurstTraits()
    #vocTraits = ["minFrequency"]
    eventListToTest = getFigureBehaviouralEvents( withUrinate = True )
    
    for experiment in experiments:
        
        file = experiment.file
        print(file)
        
        expName = experiment.getFullName()
        
        print ("Current file: " , file )
        chrono = Chronometer("File time")
        connection = sqlite3.connect( file )        
        print("Load all vocs...")
        eventTimeLineVoc = EventTimeLine( connection, "Voc", loadEventIndependently = True )
        
        burstList = createBurstFromVoc( eventTimeLineVoc )

        '''
        durationList = []
        for burst in burstList:
            duration =burst.getDurationMs()
            if duration == 0:
                print ("-- zero duration")
                print ( len( burst.vocEventList )  )
                print( burst.vocEventList[0].metadata["durationMs"])        
        quit()
        '''
        
        print( "Prepare and preload data...")
        allData = {}
        for burstTrait in burstTraits:
            allData[burstTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + burstTrait + " / " + eventName )
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[burstTrait][eventName] = getBurstTraitDataSolo( burstList , eventTimeLine , burstTrait )
                setResult(results, expName, burstTrait, eventName , allData[burstTrait][eventName] )
        connection.close()
        
    saveResult( results )



def plotBurstTraitUsagePerEventContextPerFile( experiments ):

    results = loadResult()
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    #vocTraits = getFigureVocTraits()
    burstTraits = getFigureBurstTraits()
    
    eventListToTest = getFigureBehaviouralEvents( withUrinate = True )

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
        for burstTrait in burstTraits:
            allData[burstTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + burstTrait + " / " + eventName )
                '''
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[vocTrait][eventName] = getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait )
                setResult(results, expName, vocTrait, eventName , allData[vocTrait][eventName] )
                '''
        
                allData[burstTrait][eventName] = results[expName][burstTrait][eventName] 
        
        
        # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.
        
        axIndex = 0
        nbrows=4
        fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(burstTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
        
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
        correction = getBonferroniCorrection( len( eventListToTest ) , 3, len( burstTraits ) ) # 3 for 5we, 3mo, 7mo
        
        for burstTrait in burstTraits:
            
            print ("Current burst trait: " , burstTrait )
                            
            textsToPlot = []
            
            data = {}
            
            plt.sca( axes[axIndex ])
            plt.gca().set_title( burstTrait )
    
            # kruskal test
            kruskalList = []
            for eventName in eventListToTest:
            
                traitData = allData[burstTrait][eventName]
                print("Event: ", eventName)
                #print(list(traitData))
                print("length: ", len(list(traitData)))
                kruskalList.append( list( traitData ) )
                print("---------------")
            
            
            kruskalTest = stats.kruskal( *kruskalList )
            kruskalTestOk = True
            if kruskalTest.pvalue > 0.05/ (len(burstTraits)):
            
                kruskalTestOk = False
                print("KRUSKAL FAILED FOR VOCTRAIT ", burstTrait )
               
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
                                        
                        mainTraitData = allData[burstTrait][mainEventName]
                        secondaryTraitData = allData[burstTrait][secondaryEventName] 
                                                
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
                
        plt.savefig( "burstTraitUsagePerEventContext_"+slugify( expName) +".pdf" , dpi=100 )
        #plt.show()
        #plt.text( 0 , 0 , file )
        #debugOut.close()
        #print( "Time for file (s): " , chrono.getTimeInS() )

def getBonferroniCorrection( numberOfEventTested , numberOfAgeClassTested, numberOfVocTraitTested ):
    
    correction = ( numberOfEventTested* (numberOfEventTested-1) ) / 2 # number of eventTested
    correction*= numberOfAgeClassTested
    correction*= numberOfVocTraitTested
    return correction

def extractPValueFromResult( result ):
    r = result.summary().as_text()
    for l in r.split("\n"):
        if "main_secondary" in l:
            print (l)
            lineWithoutSpace = ' '.join(l.split())            
            pValue = float( lineWithoutSpace.split(" ")[4] )
            sign = 1
            print( "test: " , lineWithoutSpace.split(" ")[1] )
            if float( lineWithoutSpace.split(" ")[1] ) < 0:
                sign=-1
            print ( "P VALUE :" , pValue )
            print ( "SIGN :" , sign )
            return pValue, sign

def plotBurstTraitUsagePerEventContextPerSet( experiments, genotype, burstTrait, ax, letter ):

    results = loadResult()
    
    if "allSet" not in results:
        results["allSet"] = {}
    results["allSet"][genotype] = {}
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    burstTraits = getFigureBurstTraits()

    burstTraits = [burstTrait]

    #burstTraits.remove( "meanInterval" )
    #burstTraits.remove( "stdInterval" )
    burstTraitLabels = []
    for trait in burstTraits:
        burstTraitLabels.append(getFigureLabelBurstTraits(trait))
    
    eventListToTest = getFigureBehaviouralEvents( longList = False , withUrinate = False )
    #eventListToTest.reverse()
    #print( eventListToTest )
    eventLabels = []
    for event in eventListToTest:
        eventLabels.append(getFigureBehaviouralEventsLabels(event))
    
    print( "Prepare and preload data...")
    allData = {}
    
    sns.set(font_scale=1.5)
    file = ""
    
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
        for burstTrait in burstTraits:
            if burstTrait not in allData:
                allData[burstTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + burstTrait + " / " + eventName )
                '''
                eventTimeLine = EventTimeLineCached( connection, file, eventName )
                print("Nb of events {}: {}".format(eventName, len(eventTimeLine.getEventList())))
                allData[vocTrait][eventName] = getTraitDataSolo( eventTimeLineVoc , eventTimeLine , vocTrait )
                setResult(results, expName, vocTrait, eventName , allData[vocTrait][eventName] )
                '''
                if eventName not in allData[burstTrait]:
                    allData[burstTrait][eventName] = {}
                    for experiment in experiments: ## Add experiment name for pairs.
                        allData[burstTrait][eventName][experiment.name] = []
                    
                allData[burstTrait][eventName][experiment.name].extend( results[expName][burstTrait][eventName] ) 
    
    
        
    # pour chaque trait de voc, on regarde si par events ils sont utilisés de la même facon.

    '''
    axIndex = 0
    nbrows=4
    #fig, axes = plt.subplots( nrows = nbrows, ncols = int(len(burstTraits)/nbrows) , sharex=True, sharey=True, figsize = (20,20) ) # , figsize = (12,12)
    fig, axes = plt.subplots( nrows = 1, ncols = len(burstTraitLabels) , sharex=True, sharey=True, figsize = (40,10) ) # , figsize = (12,12)
    '''
    
    
    # flatten axes to ease ax indexing
    

    '''
    print( axes )
    axes2 = []
    for value in axes:
        axes2.extend( value )
    axes = axes2    
    print( "axes2", axes2 )
    '''

    
    
    
    # let's go
    plt.title( file )
    
    # bonferroni correction
    ''' 
    correction = (len(eventListToTest)*(len(eventListToTest)-1)/2) # number of eventTested
    correction*= 3 # number of age class
    correction*= len( vocTraits )
    '''
    #correction = getBonfferoniCorrection( len( eventListToTest ) , 3, len( burstTraits ) ) # 3 for 5we, 3mo, 7mo

    correction = getBonferroniCorrection( len( eventListToTest ) , 3, len( getFigureBurstTraits() ) ) # 3 for 5we, 3mo, 7mo

    for burstTrait in burstTraits:
        
        print ("Current burst trait: " , burstTrait )
                        
        textsToPlot = []
        
        data = {}
        
        plt.sca( ax )
        plt.gca().set_title( getFigureLabelBurstTraits(burstTrait), fontsize = 24 )
        
        results["allSet"][genotype][burstTrait] = {}

        # kruskal test
        kruskalList = []
        for eventName in eventListToTest:
            for experiment in experiments: ## Add experiment name for pairs.                
                traitData = allData[burstTrait][eventName][experiment.name]
                
                print("Event: ", eventName)
                #print(list(traitData))
                print("length: ", len(list(traitData)))
                kruskalList.append( list( traitData ) )
                print("---------------")
        
        print (" filter kruskal" )
        kruskalList2 = []
        print ( len ( kruskalList ) )
        for val in kruskalList:
            res = []
            for v in val:
                if v != None:            
                    res.append(v)
                    
            kruskalList2.append( res )
        kruskalList = kruskalList2
              
        #print( "KruskalList: " , kruskalList )
        kruskalTest = stats.kruskal( *kruskalList )
        # remove None in kruskalList:
        
        kruskalTestOk = True
        if kruskalTest.pvalue > 0.05/ (len(burstTraits)):
        
            kruskalTestOk = False
            print("KRUSKAL FAILED FOR BURST TRAIT ", burstTrait )
           
        if not kruskalTestOk:
            for mainEventName in eventListToTest:
                data[ mainEventName ] = []
                results["allSet"][genotype][burstTrait][ mainEventName ] = []                    
                for secondaryEventName in eventListToTest:
                    data[ mainEventName ].append( 0 )
                    results["allSet"][genotype][burstTrait][ mainEventName ].append( 0 )                    
                    x = eventListToTest.index( mainEventName )
                    y = eventListToTest.index( secondaryEventName )
                    textsToPlot.append( ( x , y , "K" ) )
                    
        
        
        if kruskalTestOk:
            for mainEventName in eventListToTest:
                
                
                #mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
                data[ mainEventName ] = []
                results["allSet"][genotype][burstTrait][ mainEventName ] = []
                
                kruskalList = []
                
                print("Computing ", mainEventName )
                
                # pour un contexte:
                # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
                
                
                for secondaryEventName in eventListToTest:

                    error = False
                    stat = 0
                    p = 1000000000
                    sign = 1
                    
                    '''
                    mainTraitData = allData[burstTrait][mainEventName][experiment.name]
                    secondaryTraitData = allData[burstTrait][secondaryEventName][experiment.name]
                    '''
                    
                    pairs = []
                    values = []
                    main_secondary = []
                    
                    # mainEventName
                    for experiment in experiments:
                        for value in allData[burstTrait][mainEventName][experiment.name]:
                            values.append( value )
                            pairs.append( experiment.name )
                            main_secondary.append("MainEvent")

                    # secondaryEventName
                    for experiment in experiments:
                        for value in allData[burstTrait][secondaryEventName][experiment.name]:
                            values.append( value )
                            pairs.append( experiment.name )
                            main_secondary.append("SecondaryEvent")
                                                
                    dfData = pandas.DataFrame({'pairs': pairs,
                               'main_secondary': main_secondary,
                               'value': values})
                    
                    try:
                        # create model: value as a function of main/secondary, with the factor of the pairs:
                        model = smf.mixedlm("value ~ main_secondary", dfData, groups=dfData["pairs"])
                        # run model:
                        result = model.fit()    
                        # get_methods( result )
                        
                        # print summary
                        print(result.summary())
                        
                        p, sign = extractPValueFromResult( result )
                        
                        #quit()
                    
                    except:
                        error= True
                        print("Stat Error.")
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "LMME" ) )
                    
                    '''                                    
                    mainTraitData = allData[burstTrait][mainEventName]
                    secondaryTraitData = allData[burstTrait][secondaryEventName] 
                                            
                    stat = 0
                    p = 0
                    error = False
                    try:
                        #debugOut.write( "{} - {} : N: {} {}\n".format( mainEventName , secondaryEventName, len( mainTraitData ) , len( secondaryTraitData ) ) )
                        mainTraitDataToTest = [x for x in mainTraitData if x is not None]
                        secondaryTraitDataToTest = [x for x in secondaryTraitData if x is not None]
                        #stat, p = mannwhitneyu( mainTraitDataToTest, secondaryTraitDataToTest )
                        
                        dfData = pandas.DataFrame({'group': groups,
                               'genotype': genotypes,
                               'value': values})


                        

                        
                        
                    except:
                        error= True
                        print("Mann Whitney Error.")
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "MWE" ) )
                    '''
                    '''
                    if burstTrait == 'meanInterval':
                        print('main trait data: ', mainTraitDataToTest[0:40])
                        print('secondary trait data: ', secondaryTraitDataToTest[0:40])
                    '''
                    #sign = getSign( mainTraitData, secondaryTraitData )
                    
                    starValue = getStarValue( p , correction )*sign
                    data[ mainEventName ].append( starValue )
                    results["allSet"][genotype][burstTrait][ mainEventName ].append( starValue )
                    
                    
        df = pd.DataFrame( data=data )

        showColorBar = False
        '''
        if (axIndex-4) % nbrows == 0:    
            showColorBar = True
        '''

        g = sns.heatmap(df, linewidths=.5, yticklabels=False, xticklabels=False,
                        cmap=sns.color_palette("coolwarm", 7),
                        vmin=-3, vmax=3, cbar=showColorBar, ax=ax)
        ax.set_ylim(9, -1.5)
        tickPos = []
        for i in list(range(len(eventLabels))):
            tickPos.append(i + 0.5)
        ax.set_xticks(tickPos)
        ax.set_yticks(tickPos)
        ax.set_xticklabels(eventLabels, rotation=45, fontsize=12, horizontalalignment='right')
        ax.set_yticklabels(eventLabels, rotation=0, fontsize=12, verticalalignment='center')

        for text in textsToPlot:
            ax.text(text[0] + 0.5, text[1] + 0.5, text[2], fontsize=6, horizontalalignment='center',
                    verticalalignment='center')

        ax.text(-7, -1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')

        #axIndex += 1
        
    #plt.tight_layout()
    
    name = ""
    for experiment in experiments:
        name+= slugify ( experiment.getFullName() ) + "__"

    #plt.savefig( "750ms-burstTraitUsagePerEventContext_AllSet_LMM2"+name[:20]+".pdf" , dpi=100 )
    #plt.savefig( "750ms-burstTraitUsagePerEventContext_AllSet_LMM2"+name[:20]+".svg" , dpi=100 )
    #plt.show()
    #plt.text( 0 , 0 , file )
    #debugOut.close()
    #print( "Time for file (s): " , chrono.getTimeInS() )

    saveResult(results)
    print( results["allSet"] )

def plotBurstTraitUsagePerEventContextPerSet2( experiments, genotype, burstTrait, ax, letter ):

    results = loadResult()
    
    if "allSet" not in results:
        results["allSet"] = {}
    results["allSet"][genotype] = {}
    
    #vocTraits, eventListToTest = getVocTraitsAndEvents() 
    burstTraits = getFigureBurstTraits()

    burstTraits = [burstTrait]

    #burstTraits.remove( "meanInterval" )
    #burstTraits.remove( "stdInterval" )
    burstTraitLabels = []
    for trait in burstTraits:
        burstTraitLabels.append(getFigureLabelBurstTraits(trait))
    
    eventListToTest = getFigureBehaviouralEvents( longList = False , withUrinate = False )
    #eventListToTest.reverse()
    #print( eventListToTest )
    eventLabels = []
    for event in eventListToTest:
        eventLabels.append(getFigureBehaviouralEventsLabels(event))
    
    print( "Prepare and preload data...")
    allData = {}
    
    sns.set(font_scale=1.5)
    file = ""
    
    for experiment in experiments:
    
        file = experiment.file
        print(file)
        
        expName = experiment.getFullName()
        
        for burstTrait in burstTraits:
            if burstTrait not in allData:
                allData[burstTrait] = {}
            for eventName in eventListToTest:
                print("PRELOAD " + burstTrait + " / " + eventName )
                
                if eventName not in allData[burstTrait]:
                    allData[burstTrait][eventName] = {}
                    for experiment in experiments: ## Add experiment name for pairs.
                        allData[burstTrait][eventName][experiment.name] = []
                    
                allData[burstTrait][eventName][experiment.name].extend( results[expName][burstTrait][eventName] ) 
    
    
    # let's go
    #plt.title( file )
    
    # bonferroni correction

    correction = getBonferroniCorrection( len( eventListToTest ) , 3, len( getFigureBurstTraits() ) ) # 3 for 5we, 3mo, 7mo

    for burstTrait in burstTraits:
        
        print ("Current burst trait: " , burstTrait )
                        
        textsToPlot = []
        
        dataPValue = {}
        dataSizeEffect = {}
        
        plt.sca( ax )
        plt.gca().set_title( getFigureLabelBurstTraits(burstTrait), fontsize = 24 )
        
        results["allSet"][genotype][burstTrait] = {}

        # kruskal test
        kruskalList = []
        for eventName in eventListToTest:
            for experiment in experiments: ## Add experiment name for pairs.                
                traitData = allData[burstTrait][eventName][experiment.name]
                
                print("Event: ", eventName)
                #print(list(traitData))
                print("length: ", len(list(traitData)))
                kruskalList.append( list( traitData ) )
                print("---------------")
        
        print (" filter kruskal" )
        kruskalList2 = []
        print ( len ( kruskalList ) )
        for val in kruskalList:
            res = []
            for v in val:
                if v != None:            
                    res.append(v)
                    
            kruskalList2.append( res )
        kruskalList = kruskalList2
              
        #print( "KruskalList: " , kruskalList )
        kruskalTest = stats.kruskal( *kruskalList )
        # remove None in kruskalList:
        
        kruskalTestOk = True
        if kruskalTest.pvalue > 0.05/ (len(burstTraits)):
        
            kruskalTestOk = False
            print("KRUSKAL FAILED FOR BURST TRAIT ", burstTrait )
           
        if not kruskalTestOk:
            for mainEventName in eventListToTest:
                dataPValue[ mainEventName ] = []
                dataSizeEffect[ mainEventName ] = []
                #data[ mainEventName ] = []
                results["allSet"][genotype][burstTrait][ mainEventName ] = []                    
                for secondaryEventName in eventListToTest:
                    dataPValue[ mainEventName ].append( 10000 )
                    dataSizeEffect[ mainEventName ].append( 0 )
                    #data[ mainEventName ].append( 0 )
                    results["allSet"][genotype][burstTrait][ mainEventName ].append( 0 )                    
                    x = eventListToTest.index( mainEventName )
                    y = eventListToTest.index( secondaryEventName )
                    textsToPlot.append( ( x , y , "K" ) )
                    
        
        
        if kruskalTestOk:
            for mainEventName in eventListToTest:
                
                
                #mainEventTimeLine = EventTimeLineCached( connection, file, mainEventName )
                #data[ mainEventName ] = []
                dataPValue[ mainEventName ] = []
                dataSizeEffect[ mainEventName ] = []
                    
                results["allSet"][genotype][burstTrait][ mainEventName ] = []
                
                kruskalList = []
                
                print("Computing ", mainEventName )
                
                # pour un contexte:
                # on fait la liste des traits de toutes les vocs qui sont dans les sous-contextes
                
                
                for secondaryEventName in eventListToTest:

                    error = False
                    stat = 0
                    p = 1000000000
                    sign = 1
                    
                    '''
                    mainTraitData = allData[burstTrait][mainEventName][experiment.name]
                    secondaryTraitData = allData[burstTrait][secondaryEventName][experiment.name]
                    '''
                    
                    pairs = []
                    values = []
                    main_secondary = []
                    
                    # compute sizeEffectNormalized
                    mainEventValues = []
                    secondaryEventValues = []
                    # mainEventName
                    for experiment in experiments:
                        for value in allData[burstTrait][mainEventName][experiment.name]:
                            values.append( value )
                            mainEventValues.append( value )
                            pairs.append( experiment.name )
                            main_secondary.append("MainEvent")

                    # secondaryEventName
                    for experiment in experiments:
                        for value in allData[burstTrait][secondaryEventName][experiment.name]:
                            values.append( value )
                            secondaryEventValues.append( value )
                            pairs.append( experiment.name )
                            main_secondary.append("SecondaryEvent")
                                                
                    dfData = pandas.DataFrame({'pairs': pairs,
                               'main_secondary': main_secondary,
                               'value': values})
                    
                    # sizeEffectNormalized
                    meanMain = np.mean( mainEventValues )
                    meanSecondary = np.mean( secondaryEventValues )
                    valList= []
                    for v in mainEventValues:
                        dif = (v-meanMain)**2
                        valList.append( dif ) 
                    for v in secondaryEventValues:
                        dif = (v-meanSecondary)**2
                        valList.append( dif ) 
                    stdCombined = math.sqrt( np.mean( valList ) )
                    sizeEffect = ( meanSecondary-meanMain ) / stdCombined
                    
                    try:
                        # create model: value as a function of main/secondary, with the factor of the pairs:
                        model = smf.mixedlm("value ~ main_secondary", dfData, groups=dfData["pairs"])
                        # run model:
                        result = model.fit()    
                        # get_methods( result )
                        
                        # print summary
                        print(result.summary())
                        
                        p, sign = extractPValueFromResult( result )
                        
                        
                        #quit()
                    
                    except:
                        error= True
                        print("Stat Error.")
                        x = eventListToTest.index( mainEventName )
                        y = eventListToTest.index( secondaryEventName )
                        textsToPlot.append( ( x , y , "LMME" ) )
                    
                    
                    pValueCorrected = p*correction
                    #data[ mainEventName ].append( starValue )
                    
                    '''
                    pValueCorrected = p*correction
                    dataPValue[mainEventName].append( pValueCorrected )
                    dataSizeEffect[mainEventName].append( sizeEffect )
                    results["allSet"][genotype][vocTrait][mainEventName].append( pValueCorrected )
    
                    dfPValue = pd.DataFrame(data=dataPValue , index = eventListToTest )
                    dfSizeEffect = pd.DataFrame(data=dataSizeEffect , index = eventListToTest  )
                    g = featureHeatMap( dfSizeEffect, dfPValue , ax =ax , title = variable )
                    '''
                    #starValue = getStarValue( p , correction )*sign
                    dataPValue[ mainEventName ].append( pValueCorrected )
                    dataSizeEffect[ mainEventName ].append( sizeEffect )
                    #results["allSet"][genotype][burstTrait][ mainEventName ].append( starValue )
                    
                    
        dfPValue = pd.DataFrame(data=dataPValue , index = eventListToTest )
        dfSizeEffect = pd.DataFrame(data=dataSizeEffect , index = eventListToTest  )
        g = featureHeatMap( dfSizeEffect, dfPValue , ax =ax , title = "nb USVs / burst" )
        
        '''
        df = pd.DataFrame( data=data )
        showColorBar = False
        g = sns.heatmap(df, linewidths=.5, yticklabels=False, xticklabels=False,
                        cmap=sns.color_palette("coolwarm", 7),
                        vmin=-3, vmax=3, cbar=showColorBar, ax=ax)
        '''
        
        #ax.set_ylim(10, 0)
        '''
        tickPos = []
        
        for i in list(range(len(eventLabels))):
            tickPos.append( i )
        
        ax.set_xticks(tickPos)
        ax.set_yticks(tickPos)
        '''
        ax.set_xticklabels(eventLabels, rotation=45, fontsize=12, horizontalalignment='right')
        ax.set_yticklabels(eventLabels, rotation=0, fontsize=12, verticalalignment='center')


        for text in textsToPlot:
            ax.text(text[0] + 0.5, text[1] + 0.5, text[2], fontsize=6, horizontalalignment='center',
                    verticalalignment='center')

        ax.text(-7, -1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')

        #axIndex += 1
        
    #plt.tight_layout()
    
    name = ""
    for experiment in experiments:
        name+= slugify ( experiment.getFullName() ) + "__"

    #plt.savefig( "750ms-burstTraitUsagePerEventContext_AllSet_LMM2"+name[:20]+".pdf" , dpi=100 )
    #plt.savefig( "750ms-burstTraitUsagePerEventContext_AllSet_LMM2"+name[:20]+".svg" , dpi=100 )
    #plt.show()
    #plt.text( 0 , 0 , file )
    #debugOut.close()
    #print( "Time for file (s): " , chrono.getTimeInS() )

    saveResult(results)
    print( results["allSet"] )




def plotBurstTraitsBoxplotsPerAgePerContexts(letter, data, experiments, ax, behavEventList, burstTrait, color=None):
    expList = []
    for experiment in experiments:
        expList.append(experiment.getFullName())
    print(expList)

    expCol = []
    eventCol = []
    valuesCol = []
    genoCol = []

    for exp in expList:
        for event in getFigureBehaviouralEvents(longList=False, withUrinate=False):
            valuesList = data[exp][burstTrait][event]
            expName = [exp] * len(valuesList)
            eventName = [event] * len(valuesList)
            genotypeList = ['WT'] * len(valuesList)

            expCol.extend(expName)
            eventCol.extend(eventName)
            valuesCol.extend(valuesList)
            genoCol.extend(genotypeList)


    burstDic = {'exp': expCol, 'event': eventCol, 'values': valuesCol, 'genotype': genoCol}
    burstDf = pd.DataFrame.from_dict(burstDic)

    yLabel = "nb USV / bursts"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel(yLabel, fontsize=15)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(behavEventList)+1)
    ax.set_ylim(0, 80)
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1.5, 80 + 0.06 * (80 - 0), letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='event', y='values', data=burstDic, ax=ax,
                     showfliers=False, showmeans=True,
                     notch=True, meanprops=meanprops, width=0.4, dodge=True)

    ax.set_xticklabels([getFigureBehaviouralEventsLabels(ev) for ev in behavEventList], rotation=45, fontsize=16, horizontalalignment='right', verticalalignment='top')
    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)

    colorList = color
    edgeList = 'black'
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList)
        box.set_edgecolor(edgeList)
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.legend().set_visible(False)



    print('Done')


def plotBurstTraitsBoxplotsPerAgePerContextsPerGeno(letter, data, experimentsWT, experimentsKO, ax, behavEventList, burstTrait, color=None):
    expListWT = []
    for experiment in experimentsWT:
        expListWT.append(experiment.getFullName())
    print(expListWT)

    expListKO = []
    for experiment in experimentsKO:
        expListKO.append(experiment.getFullName())
    print(expListKO)

    expCol = []
    eventCol = []
    valuesCol = []
    genoCol = []

    for exp in expListWT:
        for event in getFigureBehaviouralEvents(longList=False, withUrinate=False):
            valuesList = data[exp][burstTrait][event]
            expName = [exp] * len(valuesList)
            eventName = [event] * len(valuesList)
            genotypeList = ['WT'] * len(valuesList)

            expCol.extend(expName)
            eventCol.extend(eventName)
            valuesCol.extend(valuesList)
            genoCol.extend(genotypeList)

    for exp in expListKO:
        for event in getFigureBehaviouralEvents(longList=False, withUrinate=False):
            valuesList = data[exp][burstTrait][event]
            expName = [exp] * len(valuesList)
            eventName = [event] * len(valuesList)
            genotypeList = ['KO'] * len(valuesList)

            expCol.extend(expName)
            eventCol.extend(eventName)
            valuesCol.extend(valuesList)
            genoCol.extend(genotypeList)

    burstDic = {'exp': expCol, 'event': eventCol, 'values': valuesCol, 'genotype': genoCol}
    burstDf = pd.DataFrame.from_dict(burstDic)

    yLabel = getFigureLabelBurstTraits(burstTrait)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel(yLabel, fontsize=15)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(behavEventList))
    ax.set_ylim(0, 80)
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, 80 + 0.06 * (80 - 0), letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='event', y='values', hue='genotype', hue_order=['WT', 'KO'], data=burstDf, ax=ax,
                     showfliers=False, showmeans=True,
                     notch=True, meanprops=meanprops, width=0.4, dodge=True)

    # Mixed model: variable to explain: value; fixed factor = genotype; random effect: group
    n = 0
    for event in behavEventList:
        testDf = burstDf[burstDf['event']==event]
        # create model:
        model = smf.mixedlm("values ~ genotype", testDf, groups=testDf['exp'])
        # run model:
        result = model.fit()
        # print summary
        print(burstTrait)
        print(result.summary())
        p, sign = extractPValueFromLMMResult(result=result, keyword='WT')
        # add p-values on the plot
        correction = len(behavEventList)
        if p * correction >= 0.05:
            stars = getStarsFromPvalues(p, 1)
        elif p * correction < 0.05:
            stars = getStarsFromPvalues(p, 1) + '°'
        ax.text(n, 80 - 0.06 * (80 - 0), stars,
                fontsize=16, horizontalalignment='center', color='black', weight='bold')
        n += 1

    ax.set_xticklabels([getFigureBehaviouralEventsLabels(ev) for ev in behavEventList], rotation=45, fontsize=16, horizontalalignment='right', verticalalignment='top')
    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)

    colorList = [getColorWT(), getColorKO()] * len(behavEventList)
    edgeList = 'black'
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList)
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.legend().set_visible(False)



    print('Done')



if __name__ == '__main__':
    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})

    experiments = getExperimentList()
    
    
    while True:
        question = "Do you want to:"
        question += "\n\t [c]ompute burst characteristics?"
        question += "\n\t [bp] plot boxplots of burst trait variations between contexts?"
        question += "\n"
        answer = input(question)

        if answer == "c":
            experiments = getExperimentList()
            computeBurstTraitUsagePerEventContext( experiments )
            break


        if answer == 'bp':
            jsonFile = "burstTraitUsagePerEventContext.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file re-imported.")

            fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(8, 10), sharey='row', sharex=True)
            row = 0
            ax = axes[row]
            experiments = getExperimentList(age='5we', sex='female', genotype='WT')
            plotBurstTraitsBoxplotsPerAgePerContexts(letter='A', data=data, experiments=experiments, ax=ax, behavEventList=getFigureBehaviouralEvents(longList=False, withUrinate=False), burstTrait='nbUSV', color=getColorAge('5we'))
            row += 1

            ax = axes[row]
            experiments = getExperimentList(age='3mo', sex='female', genotype='WT')
            plotBurstTraitsBoxplotsPerAgePerContexts(letter='B', data=data, experiments=experiments, ax=ax,
                                                     behavEventList=getFigureBehaviouralEvents(longList=False, withUrinate=False), burstTrait='nbUSV', color=getColorAge('3mo'))

            row += 1

            ax = axes[row]
            experiments = getExperimentList(age='7mo', sex='female', genotype='WT')
            plotBurstTraitsBoxplotsPerAgePerContexts(letter='C', data=data, experiments=experiments, ax=ax,
                                                     behavEventList=getFigureBehaviouralEvents(longList=False, withUrinate=False), burstTrait='nbUSV', color=getColorAge('7mo'))

            row += 1

            ax = axes[row]
            experimentsWT = getExperimentList(age='3mo', sex='female', genotype='WT')
            experimentsKO = getExperimentList(age='3mo', sex='female', strain='Shank3')
            plotBurstTraitsBoxplotsPerAgePerContextsPerGeno(letter='D', data=data, experimentsWT=experimentsWT, experimentsKO=experimentsKO, ax=ax,
                                                     behavEventList=getFigureBehaviouralEvents(longList=False, withUrinate=False), burstTrait='nbUSV', color=getColorKO())

            fig.tight_layout()
            fig.savefig('Fig_burstTrait.pdf', dpi=300)
            plt.close(fig)

            break

    print("All done")

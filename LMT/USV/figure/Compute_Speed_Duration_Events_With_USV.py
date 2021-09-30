'''
Created on 18 avr. 2020

@author: Elodie
'''

from matplotlib.offsetbox import OffsetImage, AnnotationBbox


from lmtanalysis.Util import getMinTMaxTAndFileNameInput, getMinTMaxTInput

from scipy.stats import mannwhitneyu

from LMT.USV.experimentList.experimentList import getExperimentList    

from scipy.spatial.transform import rotation
from LMT.USV.figure.figParameter import *

from LMT.USV.figure.figUtil import getStarsFromPvalues, addJitter

from scipy.stats.morestats import wilcoxon
import string
import matplotlib.image as mpimg
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import EventTimeLine
import json
import numpy as np
import sqlite3

def checkIfOverlapWith( eventBehavToCheck , vocDictionnary ):
    
    for t in range ( eventBehavToCheck.startFrame, eventBehavToCheck.endFrame+1 ) :
        if t in vocDictionnary:
            return True
            
    return False


def getDurationSpeed(event, animal):
                    
    duration = event.duration()
    sum = 0
    for t in range ( event.startFrame, event.endFrame+1 ) :
        speed = animal.getSpeed(t)
        if ( speed != None ):
            sum+= speed
    
    meanSpeed = sum / event.duration()
    
    return ( duration, meanSpeed )


def getDurationMeanSpeed(event, animal):
                    
    duration = event.duration()
    distance = animal.getDistance(tmin=event.startFrame, tmax=event.endFrame)
    
    meanSpeed = distance / duration
    
    return ( duration, meanSpeed )


def computeEventsOverlappingWithUsv( strain = None, age = None, sex = None, eventListToTest = None):
    
    print("Compute the number of events overlapping with USV sequences ")
    experiments = getExperimentList( strain = strain, age = age, sex= sex )
    tmin, tmax = getMinTMaxTInput()
    
    expNameList = []
    for experiment in experiments:

        expName = experiment.getFullName()
        expNameList.append(expName)
        
    print("Experiments: ", expNameList)
    ###########################################################################################
    durationSpeedListWithOverlap = {}
    durationSpeedListWithoutOverlap = {}
        
    for eventToTest in eventListToTest:
        durationSpeedListWithOverlap[eventToTest] = {}
        durationSpeedListWithoutOverlap[eventToTest] = {}  
        
        for experiment in experiments:
        
            file = experiment.file
            print(file)
            expName = experiment.getFullName()
            expNameList.append(expName)
            
            connection = sqlite3.connect( file )
            pool = AnimalPool( )
            pool.loadAnimals( connection )
            pool.loadDetection(tmin, tmax, lightLoad=True )
            
            #load USV timeline
            print("loading complete usv timeline dictionary")
            # all USVs
            usvTimeLine = EventTimeLine( connection, "Voc", idA=None, minFrame = tmin, maxFrame = tmax, loadEventIndependently=True )
            cleanVoc(usvTimeLine)
            USVTimeLineCompleteDictionnary = usvTimeLine.getDictionnary()
            
            #load all event timelines
            print("Loading all timelines")
            
            
            eventToTestTimeLine = {}
            
            for animal in pool.animalDictionnary.keys():
    
                print ( "loading  for animal: ", pool.animalDictionnary[animal].RFID)
                eventToTestTimeLine[animal] = EventTimeLine( connection, eventToTest, idA=animal, minFrame=tmin, maxFrame=tmax )
            
            # dic per animal and per file        
            eventToTestTimeLineWithoutOverlap = {}
            eventToTestTimeLineWithOverlap = {}
            
            for animal in pool.animalDictionnary.keys():
    
                eventToTestTimeLineWithoutOverlap[animal] = EventTimeLine( connection, "EvToTestWithoutOverlap", idA=animal, loadEvent=False )
                eventToTestTimeLineWithOverlap[animal] = EventTimeLine( connection, "EvToTestWithOverlap", idA=animal, loadEvent=False )
                
            for animal in pool.animalDictionnary.keys():
                eventToBuildOverlapDic = {}
                eventToBuildNoOverlapDic = {}
                print("Processing animal ",animal )
                
                for eventBehavToCheck in eventToTestTimeLine[animal].getEventList():
                                        
                    over = checkIfOverlapWith( eventBehavToCheck , USVTimeLineCompleteDictionnary )    
    
                    if (over == True ):
                        for t in range ( eventBehavToCheck.startFrame, eventBehavToCheck.endFrame+1 ) :
                            eventToBuildOverlapDic[t] = True;
                    else:
                        for t in range ( eventBehavToCheck.startFrame, eventBehavToCheck.endFrame+1 ) :
                            eventToBuildNoOverlapDic[t] = True;
    
        
                print( "build" )
        
                eventToTestTimeLineWithOverlap[animal].reBuildWithDictionnary( eventToBuildOverlapDic )
                eventToTestTimeLineWithoutOverlap[animal].reBuildWithDictionnary(eventToBuildNoOverlapDic)
                

            print("Event time lines done")
            
            
            ''' Compute the results '''
            durationSpeedListWithOverlap[eventToTest][expName] = {}
            durationSpeedListWithoutOverlap[eventToTest][expName] = {}
            
            for id in pool.animalDictionnary.keys():
                animal = pool.animalDictionnary[id].RFID
                print("Computing (duration,speed) with overlap")
                durationSpeedListWithOverlap[eventToTest][expName][animal] = []
                
                for event in eventToTestTimeLineWithOverlap[id].getEventList():
                    tuplePerEventWithOverlap = getDurationSpeed( event, pool.animalDictionnary[id] )
                    
                    durationSpeedListWithOverlap[eventToTest][expName][animal].append( tuplePerEventWithOverlap )
                
                print("Computing (duration,speed) no overlap")
                durationSpeedListWithoutOverlap[eventToTest][expName][animal] = []    
                for event in eventToTestTimeLineWithoutOverlap[id].getEventList():
                    tuplePerEventWithoutOverlap = getDurationSpeed( event, pool.animalDictionnary[id] )
                    
                    durationSpeedListWithoutOverlap[eventToTest][expName][animal].append( tuplePerEventWithoutOverlap ) 
                
            print("(duration, speed) lists calculated for ", eventToTest)
            connection.close()

    print("Creating json file")                    
    durationSpeedData = {}
    durationSpeedData['with overlap'] = durationSpeedListWithOverlap
    durationSpeedData['no overlap'] = durationSpeedListWithoutOverlap
    
    with open("durationSpeedData_all.json", 'w') as jFile:
        json.dump(durationSpeedData, jFile, indent=4)
        
    print("json file created for all data")


def computeEventsOverlappingWithUsvFast(strain=None, age=None, sex=None, tmin=0, tmax=None, eventListToTest=None, strainList=None, sexList=None, ageList=None):
    print("Compute the number of events overlapping with USV sequences ")
    experiments = getExperimentList(strain=strain, age=age, sex=sex)

    expNameList = []
    for experiment in experiments:
        expName = experiment.getFullName()
        expNameList.append(expName)

    print("Experiments: ", expNameList)
    ###########################################################################################
    durationSpeedListWithOverlap = {}
    durationSpeedListWithoutOverlap = {}

    for eventToTest in eventListToTest:
        durationSpeedListWithOverlap[eventToTest] = {}
        durationSpeedListWithoutOverlap[eventToTest] = {}
        for strain in strainList:
            durationSpeedListWithOverlap[eventToTest][strain] = {}
            durationSpeedListWithoutOverlap[eventToTest][strain] = {}
            for sex in sexList:
                durationSpeedListWithOverlap[eventToTest][strain][sex] = {}
                durationSpeedListWithoutOverlap[eventToTest][strain][sex] = {}
                for age in ageList:
                    durationSpeedListWithOverlap[eventToTest][strain][sex][age] = {}
                    durationSpeedListWithoutOverlap[eventToTest][strain][sex][age] = {}

    for experiment in experiments:
        file = experiment.file
        strain = experiment.strain
        sex = experiment.sex
        age = experiment.age
        expName = experiment.getFullName()

        durationSpeedListWithOverlap[eventToTest][strain][sex][age][expName] = {}
        durationSpeedListWithoutOverlap[eventToTest][strain][sex][age][expName] = {}

        connection = sqlite3.connect(file)

        pool = AnimalPool()
        pool.loadAnimals(connection)
        pool.loadDetection(tmin, tmax, lightLoad=True)

        # load USV timeline
        print("loading complete usv timeline dictionary")
        # all USVs
        usvTimeLine = EventTimeLine(connection, "Voc", idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        cleanVoc(usvTimeLine)
        USVTimeLineCompleteDictionnary = usvTimeLine.getDictionnary()

        # load all event timelines
        print("Loading all timelines")
        for eventToTest in eventListToTest:
            durationSpeedListWithOverlap[eventToTest][strain][sex][age][expName] = {}
            durationSpeedListWithoutOverlap[eventToTest][strain][sex][age][expName] = {}

            eventToTestTimeLine = {}
            # dic per animal and per file
            eventToTestTimeLineWithoutOverlap = {}
            eventToTestTimeLineWithOverlap = {}

            for animal in pool.animalDictionnary.keys():
                print("loading  for animal: ", pool.animalDictionnary[animal].RFID)
                eventToTestTimeLine[animal] = EventTimeLine(connection, eventToTest, idA=animal, minFrame=tmin, maxFrame=tmax)

                eventToTestTimeLineWithoutOverlap[animal] = EventTimeLine(connection, "EvToTestWithoutOverlap",
                                                                          idA=animal, loadEvent=False)
                eventToTestTimeLineWithOverlap[animal] = EventTimeLine(connection, "EvToTestWithOverlap", idA=animal,
                                                                       loadEvent=False)


                eventToBuildOverlapDic = {}
                eventToBuildNoOverlapDic = {}
                print("Processing animal ", animal)

                for eventBehavToCheck in eventToTestTimeLine[animal].getEventList():

                    over = checkIfOverlapWith(eventBehavToCheck, USVTimeLineCompleteDictionnary)

                    if (over == True):
                        for t in range(eventBehavToCheck.startFrame, eventBehavToCheck.endFrame + 1):
                            eventToBuildOverlapDic[t] = True
                    else:
                        for t in range(eventBehavToCheck.startFrame, eventBehavToCheck.endFrame + 1):
                            eventToBuildNoOverlapDic[t] = True

                print("build")

                eventToTestTimeLineWithOverlap[animal].reBuildWithDictionnary(eventToBuildOverlapDic)
                eventToTestTimeLineWithoutOverlap[animal].reBuildWithDictionnary(eventToBuildNoOverlapDic)

            print("Event time lines done")

            ''' Compute the results '''
            for id in pool.animalDictionnary.keys():
                animal = pool.animalDictionnary[id].RFID
                print("Computing (duration,speed) with overlap")
                durationSpeedListWithOverlap[eventToTest][strain][sex][age][expName][animal] = []

                for event in eventToTestTimeLineWithOverlap[id].getEventList():
                    tuplePerEventWithOverlap = getDurationSpeed(event, pool.animalDictionnary[id])

                    durationSpeedListWithOverlap[eventToTest][strain][sex][age][expName][animal].append(tuplePerEventWithOverlap)

                print("Computing (duration,speed) no overlap")
                durationSpeedListWithoutOverlap[eventToTest][strain][sex][age][expName][animal] = []
                for event in eventToTestTimeLineWithoutOverlap[id].getEventList():
                    tuplePerEventWithoutOverlap = getDurationSpeed(event, pool.animalDictionnary[id])

                    durationSpeedListWithoutOverlap[eventToTest][strain][sex][age][expName][animal].append(tuplePerEventWithoutOverlap)

            print("(duration, speed) lists calculated for ", eventToTest)

    print("Creating json file")
    durationSpeedData = {}
    durationSpeedData['with overlap'] = durationSpeedListWithOverlap
    durationSpeedData['no overlap'] = durationSpeedListWithoutOverlap

    with open("durationSpeedData_all.json", 'w') as jFile:
        json.dump(durationSpeedData, jFile, indent=3)

    print("json file created for all data")


def plotDataPerAnimal(axes, row, col, xPos, ylim, yLabelText, dataPerAnimalUsv, dataPerAnimalNoUsv, stars):
    ax = axes[row][col]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)  
    ax.set_xticks( [ 1-0.3,1+0.3,2-0.3, 2+0.3  ] )
    ax.set_xticklabels( [ "USV" , "noUSV", "USV" , "noUSV"] )
    ax.set_xlim( 0 , 3 )
    ax.set_ylim( 0 , ylim )
    
    if col==0:
        ax.set_ylabel( yLabelText )
    ax.scatter(addJitter([xPos-0.3]*len(dataPerAnimalUsv), 0.05), dataPerAnimalUsv, s = 5, c = "orange")
    ax.scatter(xPos-0.3, np.mean(dataPerAnimalUsv), marker= '+', s=8, c= "black")
    ax.scatter(addJitter([xPos+0.3]*len(dataPerAnimalNoUsv), 0.05), dataPerAnimalNoUsv, s = 5, c = "blue")
    ax.text(xPos, np.mean(dataPerAnimalNoUsv), stars, horizontalalignment='center')
    ax.scatter(xPos+0.3, np.mean(dataPerAnimalNoUsv), marker= '+', s=8, c= "black")
                 
   
def plotDataPerAnimalViolin(axes, row, col, xPos, ylim, yLabelText, dataPerAnimalUsv, dataPerAnimalNoUsv, experiment, stars): 
    ax = axes[row][col]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)  
    pos = [xPos-0.3, xPos+0.3]
    ax.set_xticks( [ 1-0.3,1+0.3,3-0.3, 3+0.3 ] )
    ax.set_xticklabels( [ "USV" , "no USV", "USV" , "no USV"] , FontSize=12, rotation=45)
    ax.xaxis.set_tick_params( direction = "in")
    ax.set_xlim( 0 , 4 )
    ax.set_ylim( 0 , ylim )
    
    data = [ dataPerAnimalUsv, dataPerAnimalNoUsv ]
    
    if col==0:
        ax.set_ylabel( yLabelText, fontsize=15 )
        
    if row==0:
        ax.set_title( experiment[-7:] , fontsize=15)
    
           
    parts = ax.violinplot(data, pos, points=100, widths=0.3, showmeans=False,
                  showextrema=False, showmedians=False)
    
    colorAgeListBis = ["#efB1ff", "#7f78d2", "#481380"]
    colorText = ["white", "white", "white"]
    ageClass = ["5 weeks", "3 months", "7 months"]
    alphas = [0.8, 0.4, 0.8, 0.4]
    
    if '5we' in experiment:
        i = 0
    if '3mo' in experiment:
        i = 1
    if '7mo' in experiment:
        i = 2
    j=0    
    for pc in parts['bodies']:
        pc.set_facecolor(colorAgeList[i])
        pc.set_edgecolor(colorAgeList[i])
        pc.set_alpha(alphas[j])
        j+=1
    
    if col == 0:   
        ax.vlines( 0 , 0, ylim, color= colorAgeList[i], lw=30 )
        ax.text(0.12, ylim/3, ageClass[i], horizontalalignment='center', color=colorText[i], fontsize=12, fontweight='bold', rotation=90)
    '''
    quartile1, medians, quartile3 = np.percentile(dataPerAnimalUsv, [25, 50, 75], axis=1)
    
    ax.scatter(pos, medians, marker='o', color='white', s=30, zorder=3)
    ax.vlines(pos, quartile1, quartile3, color='k', linestyle='-', lw=5)
    '''
            
    ax.scatter(pos, [np.mean(dataPerAnimalUsv), np.mean(dataPerAnimalNoUsv)], c='black', marker= '_', s=20)
    ax.text(xPos, ylim*0.9, stars, fontsize=12, horizontalalignment='center')
    
    
def plotViolinDuration( strain=None, age=None, sex=None, eventListToTest = None ):
    #open the json file to have the data
    with open('durationSpeedData.json') as json_data:
        durationSpeedData = json.load(json_data)
    print("json file re-imported.")
    
    durationSpeedListWithOverlap = durationSpeedData['with overlap']
    durationSpeedListWithoutOverlap = durationSpeedData['no overlap']
    
    #define the list of experiments to be computed
    
    expNameList = []
    experiments = getExperimentList( strain = strain, age = age, sex = sex )
    for experiment in experiments:
        expName = experiment.getFullName()
        expNameList.append(expName)
        
    print(expNameList)
    
    
    for event in eventListToTest:
        print('event: ', event)
        
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(15, 10), sharey = 'row', sharex = True)
        
        ylimDuration = {'longChase': 1000, 'FollowZone Isolated': 50, 'Approach contact': 30, 
                        'Break contact': 100, 'Contact': 1000, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 50}

        col = 0
        
        for experiment in expNameList:
            if "5we" in experiment:
                row = 0
            if "3mo" in experiment:
                row = 1
            if "7mo" in experiment:
                row = 2
            if col > 3:
                col = 0
            
            speedUsv = {}
            speedNoUsv = {}
            durUsv = {}
            durNoUsv = {}
            
            xPos = 1
            for animal in sorted(durationSpeedListWithOverlap[event][experiment].keys()) :
                #print("Event ", event, experiment, animal)
                
                if len( durationSpeedListWithOverlap[event][experiment][animal] ) > 0:
                    
                    durUsv[animal] = []
                    durNoUsv[animal] = []
                    
                    for i in range( len(durationSpeedListWithOverlap[event][experiment][animal])):
                        durUsv[animal].append(durationSpeedListWithOverlap[event][experiment][animal][i][0])
                    
                    for j in range( len(durationSpeedListWithoutOverlap[event][experiment][animal])):
                        durNoUsv[animal].append(durationSpeedListWithoutOverlap[event][experiment][animal][j][0])
                    
                    Wdur, Pdur = mannwhitneyu(x=durUsv[animal], y=durNoUsv[animal], alternative='greater')
                    print("Duration with USV: ", len(durUsv[animal]), "no USV: ", len(durNoUsv[animal]), "U=", Wdur, " ", "p=", Pdur)
                    starsDur = getStarsFromPvalues(Pdur, 24)
                    
                    plotDataPerAnimalViolin(axes=axes, row=row, col=col, xPos=xPos, ylim = ylimDuration[event], yLabelText="duration (frames)", dataPerAnimalUsv=durUsv[animal], dataPerAnimalNoUsv=durNoUsv[animal], experiment=experiment, stars=starsDur )
                    
                    if row == 0:
                        axes[row,col].text(xPos, -ylimDuration[event]*0.08, animal[-7:], horizontalalignment='center')
                    xPos+=2
                                  
                else:
                    print("Values: no data")
                    
            
            col +=1
            
        plt.suptitle(event, y=1, fontsize=18)
        
        plt.tight_layout()
        fig.savefig( "Fig_Duration_Violin_{}.pdf".format(event) , dpi=100)             
        #plt.show()   
        
    print("job done")


def plotViolinSpeed( strain = None, age = None, sex = None, eventListToTest = None ):
    #open the json file to have the data
    with open('durationSpeedData.json') as json_data:
        durationSpeedData = json.load(json_data)
    print("json file re-imported.")
    
    durationSpeedListWithOverlap = durationSpeedData['with overlap']
    durationSpeedListWithoutOverlap = durationSpeedData['no overlap']
    
    #define the list of experiments to be computed
    expNameList = []
    experiments = getExperimentList( strain = strain, age = age, sex = sex )
    for experiment in experiments:
        expName = experiment.getFullName()
        expNameList.append(expName)
        
    print(expNameList)
    
    
    for event in eventListToTest:
        print('event: ', event)
        
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(15, 10), sharey = 'row', sharex = True)
        
        ylimSpeed = {'longChase': 30, 'FollowZone Isolated': 50, 'Approach contact': 50, 
                     'Break contact': 100, 'Contact': 50, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 60}
        
        col = 0
        
        for experiment in expNameList:
            if "5we" in experiment:
                row = 0
            if "3mo" in experiment:
                row = 1
            if "7mo" in experiment:
                row = 2
            if col > 3:
                col = 0
            
            speedUsv = {}
            speedNoUsv = {}
            
            xPos = 1
            for animal in sorted(durationSpeedListWithOverlap[event][experiment].keys()) :
                #print("Event ", event, experiment, animal)
                
                if len( durationSpeedListWithOverlap[event][experiment][animal] ) > 0:
                    
                    speedUsv[animal] = []
                    speedNoUsv[animal] = []
                    
                    for i in range( len(durationSpeedListWithOverlap[event][experiment][animal])):
                        
                        if durationSpeedListWithOverlap[event][experiment][animal][i][1] <= 70:
                            speedUsv[animal].append(durationSpeedListWithOverlap[event][experiment][animal][i][1])
                        
                    for j in range( len(durationSpeedListWithoutOverlap[event][experiment][animal])):
                        
                        if durationSpeedListWithoutOverlap[event][experiment][animal][j][1] <=70:
                            speedNoUsv[animal].append(durationSpeedListWithoutOverlap[event][experiment][animal][j][1])
                        
                    Wspeed, Pspeed = mannwhitneyu(x=speedUsv[animal], y=speedNoUsv[animal], alternative='greater')
                    print("Speed with USV: ", len(speedUsv[animal]), "no USV: ", len(speedNoUsv[animal]), "U=", Wspeed, " ", "p=", Pspeed)
                    starsSpeed = getStarsFromPvalues(Pspeed, 24)
                    
                    
                    #plotDataPerAnimal(axes=axes, row=row, col=col, xPos=xPos, ylim = 50, yLabelText="mean speed (cm/s)", dataPerAnimalUsv=speedUsv[animal], dataPerAnimalNoUsv=speedNoUsv[animal], stars=starsSpeed )
                    
                    plotDataPerAnimalViolin(axes=axes, row=row, col=col, xPos=xPos, ylim = ylimSpeed[event], yLabelText="mean speed (cm/s)", dataPerAnimalUsv=speedUsv[animal], dataPerAnimalNoUsv=speedNoUsv[animal], experiment=experiment, stars=starsSpeed )
                    
                    if row == 0:
                        axes[row,col].text(xPos, -ylimSpeed[event]*0.08, animal[-7:], horizontalalignment='center')
                        
                    xPos+=2
                                  
                else:
                    print("Values: no data")
                    
            
            col +=1
            
        plt.suptitle(event, y=1, fontsize=18)
        
        plt.tight_layout()
        fig.savefig( "Fig_Speed_Violin_{}.pdf".format(event) , dpi=100) 
                  
        #plt.show()
    
    print( "Job done.")
     

def plotMeanSpeedDuration( strain= None, age = None, sex = None, eventListToTest = None ):
    
    #open the json file to have the data
    with open('durationSpeedData.json') as json_data:
        durationSpeedData = json.load(json_data)
    print("json file re-imported.")
    
    durationSpeedListWithOverlap = durationSpeedData['with overlap']
    durationSpeedListWithoutOverlap = durationSpeedData['no overlap']
    
    #define the list of experiments to be computed
    expNameList = []
    experiments = getExperimentList( strain = strain, age = age, sex = sex )
    for experiment in experiments:
        expName = experiment.getFullName()
        expNameList.append(expName)
        
    print(expNameList)
    
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15, 10), sharey = False, sharex = True)
        
    ylimSpeed = {'longChase': 15, 'FollowZone Isolated': 25, 'Approach contact': 15, 
                     'Break contact': 40, 'Contact': 50, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 60}
    ylimDuration = {'longChase': 800, 'FollowZone Isolated': 8, 'Approach contact': 8, 
                        'Break contact': 30, 'Contact': 1000, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 50}
    
    colNb = [0, 2, 0, 2]
    rowNb = [0, 0, 1, 1]
    
    dataToTest = {}
    for var in ['speed', 'duration']:
        dataToTest[var] = {}
        for age in ['5we', '3mo', '7mo']:
            dataToTest[var][age] = {}
            for condition in ['usv', 'no usv']:
                dataToTest[var][age][condition] = []
                    
    k = 0
    
    for event in eventListToTest:
        print('event: ', event)
        print('k: ', k)
        col = colNb[k]
        row = rowNb[k]
        
        dataToTest = {}
        for var in ['speed', 'duration']:
            dataToTest[var] = {}
            for age in ['5we', '3mo', '7mo']:
                dataToTest[var][age] = {}
                for condition in ['usv', 'no usv']:
                    dataToTest[var][age][condition] = []
                
        print('col: ', col, 'row: ', row)
           
        for experiment in expNameList:
            if "5we" in experiment:
                x1 = 1
                x2 = 2
                colorPoints = getColorAge('5we')
                listSpeedToComplement = dataToTest['speed']['5we']
                listDurationToComplement = dataToTest['duration']['5we']
                
            if "3mo" in experiment:
                x1 = 4
                x2 = 5
                colorPoints = getColorAge('3mo')
                listSpeedToComplement = dataToTest['speed']['3mo']
                listDurationToComplement = dataToTest['duration']['3mo']
                
            if "7mo" in experiment:
                x1 = 7
                x2 = 8
                colorPoints = getColorAge('7mo')
                listSpeedToComplement = dataToTest['speed']['7mo']
                listDurationToComplement = dataToTest['duration']['7mo']
            
            speedUsv = {}
            speedNoUsv = {}
            durUsv = {}
            durNoUsv = {}
            
            for animal in sorted(durationSpeedListWithOverlap[event][experiment].keys()) :
                
                if len( durationSpeedListWithOverlap[event][experiment][animal] ) > 0:
                    
                    speedUsv[animal] = []
                    speedNoUsv[animal] = []
                    durUsv[animal] = []
                    durNoUsv[animal] = []
                    
                    for i in range( len(durationSpeedListWithOverlap[event][experiment][animal])):
                        
                        if durationSpeedListWithOverlap[event][experiment][animal][i][1] <= 70:
                            speedUsv[animal].append(durationSpeedListWithOverlap[event][experiment][animal][i][1])
                        
                        durUsv[animal].append(durationSpeedListWithOverlap[event][experiment][animal][i][0])
                        
                    for j in range( len(durationSpeedListWithoutOverlap[event][experiment][animal])):
                        
                        if durationSpeedListWithoutOverlap[event][experiment][animal][j][1] <=70:
                            speedNoUsv[animal].append(durationSpeedListWithoutOverlap[event][experiment][animal][j][1])
                            
                        durNoUsv[animal].append(durationSpeedListWithoutOverlap[event][experiment][animal][j][0])
                        
                    Wspeed, Pspeed = mannwhitneyu(x=speedUsv[animal], y=speedNoUsv[animal], alternative='greater')
                    starsSpeed = getStarsFromPvalues(Pspeed, 24)
                    print("Speed with USV: ", len(speedUsv[animal]), "no USV: ", len(speedNoUsv[animal]), "U=", Wspeed, " ", "p=", Pspeed, " ", starsSpeed)
                    
                    colorLineSpeed = getColorLine( starsSpeed )
                    
                    
                    Wdur, Pdur = mannwhitneyu(x=durUsv[animal], y=durNoUsv[animal], alternative='greater')
                    starsDur = getStarsFromPvalues(Pdur, 24)
                    print("Duration with USV: ", len(durUsv[animal]), "no USV: ", len(durNoUsv[animal]), "U=", Wdur, " ", "p=", Pdur, " ", starsDur)
                    
                    colorLineDur = getColorLine( starsDur )
                        
                    
                    #speed
                    ax = axes[row][col]
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)  
                    ax.set_xticks( [ 1,2, 4,5, 7,8  ] )
                    ax.set_xticklabels( [ "USV" , "no USV", "USV" , "no USV", "USV" , "no USV" ], rotation = 45, FontSize=15 )
                    ax.xaxis.set_tick_params( direction = "in")
                    ax.set_xlim( 0 , 9 )
                    ax.set_ylim( 0 , ylimSpeed[event] )
                    ax.set_ylabel( "mean speed (cm/s)", FontSize = 15 )
                    ax.tick_params(axis='y', labelsize=14)
        
                    ax.scatter( x1, np.mean(speedUsv[animal]), marker= 'o', s=16, c= colorPoints )
                    ax.scatter( x2, np.mean(speedNoUsv[animal]), marker= 'o', s=16, c= colorPoints )
                    ax.plot( [x1 , x2], [np.mean(speedUsv[animal]), np.mean(speedNoUsv[animal])] , zorder=0, color=colorLineSpeed, linewidth=0.6 )
                    
                    
                    semUsvSpeed = np.std(speedUsv[animal])/np.sqrt(len(speedUsv[animal]))
                    semNoUsvSpeed = np.std(speedNoUsv[animal])/np.sqrt(len(speedNoUsv[animal]))
                    #sem with USV
                    ax.vlines( x1 , np.mean(speedUsv[animal]) , np.mean(speedUsv[animal])+semUsvSpeed , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(speedUsv[animal])+semUsvSpeed  , x1-0.1 , x1+0.1 , linewidth = 0.2, color = colorPoints)
                    ax.vlines( x1 , np.mean(speedUsv[animal]) , np.mean(speedUsv[animal])-semUsvSpeed , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(speedUsv[animal])-semUsvSpeed  , x1-0.1 , x1+0.1 , linewidth = 0.2, color = colorPoints)
                    #sem no USV
                    ax.vlines( x2 , np.mean(speedNoUsv[animal]) , np.mean(speedNoUsv[animal])+semNoUsvSpeed , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(speedNoUsv[animal])+semNoUsvSpeed  , x2-0.1 , x2+0.1 , linewidth = 0.2, color = colorPoints)
                    ax.vlines( x2 , np.mean(speedNoUsv[animal]) , np.mean(speedNoUsv[animal])-semNoUsvSpeed , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(speedNoUsv[animal])-semNoUsvSpeed  , x2-0.1 , x2+0.1 , linewidth = 0.2, color = colorPoints)
                    
                    listSpeedToComplement['usv'].append( np.mean(speedUsv[animal]) )
                    listSpeedToComplement['no usv'].append( np.mean(speedNoUsv[animal]) )
                    
                    #duration
                    ax = axes[row][col+1]
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)  
                    ax.set_xticks( [ 1,2, 4,5, 7,8  ] )
                    ax.set_xticklabels( [ "USV" , "no USV", "USV" , "no USV", "USV" , "no USV" ], rotation = 45, FontSize=15 )
                    ax.xaxis.set_tick_params( direction = "in")
                    ax.set_xlim( 0 , 9 )
                    ax.set_ylim( 0 , ylimDuration[event] )
                    ax.set_ylabel( "duration (frames)", FontSize = 15 )
                    ax.tick_params(axis='y', labelsize=14)
        
                    ax.scatter( x1, np.mean(durUsv[animal]), marker= 'o', s=16, c= colorPoints )
                    ax.scatter( x2, np.mean(durNoUsv[animal]), marker= 'o', s=16, c= colorPoints )
                    ax.plot( [x1 , x2], [np.mean(durUsv[animal]), np.mean(durNoUsv[animal])] , zorder=0, color=colorLineDur, linewidth=1 )
                    
                    
                    semUsvDur = np.std(durUsv[animal])/np.sqrt(len(durUsv[animal]))
                    semNoUsvDur = np.std(durNoUsv[animal])/np.sqrt(len(durNoUsv[animal]))
                    #sem with USV
                    ax.vlines( x1 , np.mean(durUsv[animal]) , np.mean(durUsv[animal])+semUsvDur , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(durUsv[animal])+semUsvDur  , x1-0.1 , x1+0.1 , linewidth = 0.2, color = colorPoints)
                    ax.vlines( x1 , np.mean(durUsv[animal]) , np.mean(durUsv[animal])-semUsvDur , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(durUsv[animal])-semUsvDur  , x1-0.1 , x1+0.1 , linewidth = 0.2, color = colorPoints)
                    #sem no USV
                    ax.vlines( x2 , np.mean(durNoUsv[animal]) , np.mean(durNoUsv[animal])+semNoUsvDur , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(durNoUsv[animal])+semNoUsvDur  , x2-0.1 , x2+0.1 , linewidth = 0.2, color = colorPoints)
                    ax.vlines( x2 , np.mean(durNoUsv[animal]) , np.mean(durNoUsv[animal])-semNoUsvDur , linestyle='solid', linewidth = 0.2, color = colorPoints)
                    ax.hlines( np.mean(durNoUsv[animal])-semNoUsvDur  , x2-0.1 , x2+0.1 , linewidth = 0.2, color = colorPoints)
                    
                    listDurationToComplement['usv'].append( np.mean(durUsv[animal]) )
                    listDurationToComplement['no usv'].append( np.mean(durNoUsv[animal]) )
                    
                else:
                    print("Values: no data")
                    
        stars = {}
        
        for var in ['speed', 'duration']:
            for age in ['5we', '3mo', '7mo']:
                print ( "*************** " , event , " -var: " , var , " -geno: " , age )
                T, P = wilcoxon(x=dataToTest[var][age]['usv'], y=dataToTest[var][age]['no usv'], alternative="greater" )
                stars[var,age] = getStarsFromPvalues(P, 1)
                print("Wilcoxon ", var, " ", age, " : ", "T=", T, " ", "p=", P, " ", stars[var,age])    
        
        ax = axes[row][col]
        ax.text(1.5, ylimSpeed[event]*0.95, stars['speed', '5we'], FontSize = 12, horizontalalignment='center', color = 'black')
        ax.text(4.5, ylimSpeed[event]*0.95, stars['speed', '3mo'], FontSize = 12, horizontalalignment='center', color = 'black')
        ax.text(7.5, ylimSpeed[event]*0.95, stars['speed', '7mo'], FontSize = 12, horizontalalignment='center', color = 'black')
        
        ax = axes[row][col+1]
        ax.text(1.5, ylimDuration[event]*0.95, stars['duration', '5we'], FontSize = 12, horizontalalignment='center', color = 'black')
        ax.text(4.5, ylimDuration[event]*0.95, stars['duration', '3mo'], FontSize = 12, horizontalalignment='center', color = 'black')
        ax.text(7.5, ylimDuration[event]*0.95, stars['duration', '7mo'], FontSize = 12, horizontalalignment='center', color = 'black')
        
        k +=1
    
    
    
    m = 0        
    for event in eventListToTest:
        print('event: ', event)
        print('m: ', m)
        col = colNb[m]
        row = rowNb[m]
        
        print('col: ', col, 'row: ', row)        
        
        ax = axes[row, col]
        ax.set_title( event, loc = 'right', FontSize = 18 )
        
        if row == 0:
            ax = axes[row, col]
            ax.hlines( 0 , 0.8, 2.2, color= getColorAge('5we'), lw=30 )
            ax.text(1.5, ylimSpeed[event]*0.01, "5 we", FontSize = 12, horizontalalignment='center', color = 'white')
            ax.hlines( 0 , 3.8, 5.2, color= getColorAge('3mo'), lw=30 )
            ax.text(4.5, ylimSpeed[event]*0.01, "3 mo", FontSize = 12, horizontalalignment='center', color = 'white')
            ax.hlines( 0 , 6.8, 8.2, color= getColorAge('7mo'), lw=30 )
            ax.text(7.5, ylimSpeed[event]*0.01, "7 mo", FontSize = 12, horizontalalignment='center', color = 'white')
                    
            ax = axes[row, col+1]
            ax.hlines( 0 , 0.8, 2.2, color= getColorAge('5we'), lw=30 )
            ax.text(1.5, ylimDuration[event]*0.01, "5 we", FontSize = 12, horizontalalignment='center', color = 'white')
            ax.hlines( 0 , 3.8, 5.2, color= getColorAge('3mo'), lw=30 )
            ax.text(4.5, ylimDuration[event]*0.01, "3 mo", FontSize = 12, horizontalalignment='center', color = 'white')
            ax.hlines( 0 , 6.8, 8.2, color= getColorAge('7mo'), lw=30 )
            ax.text(7.5, ylimDuration[event]*0.01, "7 mo", FontSize = 12, horizontalalignment='center', color = 'white')
        
        m +=1                
        
    plt.tight_layout()
    fig.savefig( "Fig_Speed_Duration_Mean.pdf" , dpi=100) 
    fig.savefig( "Fig_Speed_Duration_Mean.svg" , dpi=100) 
              
    #plt.show()
    
    print( "Job done.")


def plotMeanSpeedDurationFast(jsonFile='durationSpeedData_all_classes.json', strain='C57BL/6J', sex=None, eventListToTest=None):
    # open the json file to have the data
    with open(jsonFile) as json_data:
        durationSpeedData = json.load(json_data)
    print("json file re-imported.")

    durationSpeedListWithOverlap = durationSpeedData['with overlap']
    durationSpeedListWithoutOverlap = durationSpeedData['no overlap']

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15, 10), sharey=False, sharex=True)

    ylimSpeed = {'longChase': 15, 'FollowZone Isolated': 25, 'Approach contact': 15,
                 'Break contact': 40, 'Contact': 50, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 60}
    ylimDuration = {'longChase': 800, 'FollowZone Isolated': 8, 'Approach contact': 8,
                    'Break contact': 30, 'Contact': 1000, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 50}

    colNb = [0, 2, 0, 2]
    rowNb = [0, 0, 1, 1]

    imgPosDic = {'Train2': (4.5, 12), 'FollowZone Isolated': (4.5, 5), 'Approach contact': (4.5, 2), 'Break contact': (4.5, 6)}
    imgZoomDic = {'Train2': 0.25, 'FollowZone Isolated': 0.27, 'Approach contact': 0.22, 'Break contact': 0.22}

    k = 0

    for event in eventListToTest:
        print('event: ', event)
        print('k: ', k)
        imgPos = imgPosDic[event]
        imgZoom = imgZoomDic[event]
        col = colNb[k]
        row = rowNb[k]

        dataToTest = {}
        for var in ['speed', 'duration']:
            dataToTest[var] = {}
            dataToTest[var][sex] = {}
            for ageClass in ['5we', '3mo', '7mo']:
                dataToTest[var][sex][ageClass] = {}
                for condition in ['usv', 'no usv']:
                    dataToTest[var][sex][ageClass][condition] = []

        print('col: ', col, 'row: ', row)

        x1Dic = {'5we': 1, '3mo': 4, '7mo': 7}
        x2Dic = {'5we': 2, '3mo': 5, '7mo': 8}

        speedUsv = {}
        speedNoUsv = {}
        durUsv = {}
        durNoUsv = {}

        for age in ['5we', '3mo', '7mo']:
            colorPoints = getColorAge(age)
            listSpeedToComplement = dataToTest['speed'][sex][age]
            listDurationToComplement = dataToTest['duration'][sex][age]
            x1 = x1Dic[age]
            x2 = x2Dic[age]
            for experiment in durationSpeedListWithOverlap[event][strain][sex][age].keys():
                for animal in sorted(durationSpeedListWithOverlap[event][strain][sex][age][experiment].keys()):

                    if len(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal]) > 0:

                        speedUsv[animal] = []
                        speedNoUsv[animal] = []
                        durUsv[animal] = []
                        durNoUsv[animal] = []

                        for i in range(len(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal])):

                            if durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][1] <= 70:
                                speedUsv[animal].append(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][1])

                            durUsv[animal].append(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][0])

                        for j in range(len(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal])):

                            if durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][1] <= 70:
                                speedNoUsv[animal].append(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][1])

                            durNoUsv[animal].append(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][0])

                        Wspeed, Pspeed = mannwhitneyu(x=speedUsv[animal], y=speedNoUsv[animal], alternative='greater')
                        starsSpeed = getStarsFromPvalues(Pspeed, 24)
                        print("Speed with USV: ", len(speedUsv[animal]), "no USV: ", len(speedNoUsv[animal]), "U=", Wspeed,
                              " ", "p=", Pspeed, " ", starsSpeed)

                        colorLineSpeed = getColorLine(starsSpeed)

                        Wdur, Pdur = mannwhitneyu(x=durUsv[animal], y=durNoUsv[animal], alternative='greater')
                        starsDur = getStarsFromPvalues(Pdur, 24)
                        print("Duration with USV: ", len(durUsv[animal]), "no USV: ", len(durNoUsv[animal]), "U=", Wdur,
                              " ", "p=", Pdur, " ", starsDur)

                        colorLineDur = getColorLine(starsDur)

                        # speed
                        ax = axes[row][col]
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.set_xticks([1, 2, 4, 5, 7, 8])
                        ax.set_xticklabels(["USV", "no USV", "USV", "no USV", "USV", "no USV"], rotation=45, FontSize=15)
                        ax.xaxis.set_tick_params(direction="in")
                        ax.set_xlim(0, 9)
                        ax.set_ylim(0, ylimSpeed[event])
                        ax.set_ylabel("mean speed (cm/s)", FontSize=15)
                        ax.tick_params(axis='y', labelsize=14)

                        ax.scatter(x1, np.mean(speedUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.scatter(x2, np.mean(speedNoUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.plot([x1, x2], [np.mean(speedUsv[animal]), np.mean(speedNoUsv[animal])], zorder=0,
                                color=colorLineSpeed, linewidth=0.6)

                        semUsvSpeed = np.std(speedUsv[animal]) / np.sqrt(len(speedUsv[animal]))
                        semNoUsvSpeed = np.std(speedNoUsv[animal]) / np.sqrt(len(speedNoUsv[animal]))
                        # sem with USV
                        ax.vlines(x1, np.mean(speedUsv[animal]), np.mean(speedUsv[animal]) + semUsvSpeed, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedUsv[animal]) + semUsvSpeed, x1 - 0.1, x1 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x1, np.mean(speedUsv[animal]), np.mean(speedUsv[animal]) - semUsvSpeed, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedUsv[animal]) - semUsvSpeed, x1 - 0.1, x1 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        # sem no USV
                        ax.vlines(x2, np.mean(speedNoUsv[animal]), np.mean(speedNoUsv[animal]) + semNoUsvSpeed,
                                  linestyle='solid', linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedNoUsv[animal]) + semNoUsvSpeed, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x2, np.mean(speedNoUsv[animal]), np.mean(speedNoUsv[animal]) - semNoUsvSpeed,
                                  linestyle='solid', linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedNoUsv[animal]) - semNoUsvSpeed, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)

                        listSpeedToComplement['usv'].append(np.mean(speedUsv[animal]))
                        listSpeedToComplement['no usv'].append(np.mean(speedNoUsv[animal]))

                        # duration
                        ax = axes[row][col + 1]
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.set_xticks([1, 2, 4, 5, 7, 8])
                        ax.set_xticklabels(["USV", "no USV", "USV", "no USV", "USV", "no USV"], rotation=45, FontSize=15)
                        ax.xaxis.set_tick_params(direction="in")
                        ax.set_xlim(0, 9)
                        ax.set_ylim(0, ylimDuration[event])
                        ax.set_ylabel("duration (frames)", FontSize=15)
                        ax.tick_params(axis='y', labelsize=14)

                        ax.scatter(x1, np.mean(durUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.scatter(x2, np.mean(durNoUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.plot([x1, x2], [np.mean(durUsv[animal]), np.mean(durNoUsv[animal])], zorder=0,
                                color=colorLineDur, linewidth=1)

                        semUsvDur = np.std(durUsv[animal]) / np.sqrt(len(durUsv[animal]))
                        semNoUsvDur = np.std(durNoUsv[animal]) / np.sqrt(len(durNoUsv[animal]))
                        # sem with USV
                        ax.vlines(x1, np.mean(durUsv[animal]), np.mean(durUsv[animal]) + semUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durUsv[animal]) + semUsvDur, x1 - 0.1, x1 + 0.1, linewidth=0.2, color=colorPoints)
                        ax.vlines(x1, np.mean(durUsv[animal]), np.mean(durUsv[animal]) - semUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durUsv[animal]) - semUsvDur, x1 - 0.1, x1 + 0.1, linewidth=0.2, color=colorPoints)
                        # sem no USV
                        ax.vlines(x2, np.mean(durNoUsv[animal]), np.mean(durNoUsv[animal]) + semNoUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durNoUsv[animal]) + semNoUsvDur, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x2, np.mean(durNoUsv[animal]), np.mean(durNoUsv[animal]) - semNoUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durNoUsv[animal]) - semNoUsvDur, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)

                        listDurationToComplement['usv'].append(np.mean(durUsv[animal]))
                        listDurationToComplement['no usv'].append(np.mean(durNoUsv[animal]))

                    else:
                        print("Values: no data")

        stars = {}

        for var in ['speed', 'duration']:
            correction = 3
            for age in ['5we', '3mo', '7mo']:
                print("*************** ", event, " -var: ", var, " -geno: ", age)
                T, P = wilcoxon(x=dataToTest[var][sex][age]['usv'], y=dataToTest[var][sex][age]['no usv'], alternative="greater")
                if P * correction >= 0.05:
                    stars[var, age] = getStarsFromPvalues(P, 1)
                elif P * correction < 0.05:
                    stars[var, age] = getStarsFromPvalues(P, 1) + '°'

                print("Wilcoxon ", var, " ", age, " : ", "T=", T, " ", "p=", P, " ", stars[var, age])

        ax = axes[row][col]
        ax.text(1.5, ylimSpeed[event] * 0.95, stars['speed', '5we'], FontSize=12, horizontalalignment='center',
                color='black')
        ax.text(4.5, ylimSpeed[event] * 0.95, stars['speed', '3mo'], FontSize=12, horizontalalignment='center',
                color='black')
        ax.text(7.5, ylimSpeed[event] * 0.95, stars['speed', '7mo'], FontSize=12, horizontalalignment='center',
                color='black')
        image = '{}.png'.format(event)
        behavSchema = mpimg.imread(image)
        imgBox = OffsetImage(behavSchema, zoom=imgZoom)
        imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
        ax.add_artist(imageBox)

        ax = axes[row][col + 1]
        ax.text(1.5, ylimDuration[event] * 0.95, stars['duration', '5we'], FontSize=12, horizontalalignment='center',
                color='black')
        ax.text(4.5, ylimDuration[event] * 0.95, stars['duration', '3mo'], FontSize=12, horizontalalignment='center',
                color='black')
        ax.text(7.5, ylimDuration[event] * 0.95, stars['duration', '7mo'], FontSize=12, horizontalalignment='center',
                color='black')

        k += 1

    m = 0
    letterList = ['A', 'B', 'C', 'D']
    for event in eventListToTest:
        print('event: ', event)
        print('m: ', m)
        col = colNb[m]
        row = rowNb[m]

        print('col: ', col, 'row: ', row)

        ax = axes[row, col]
        ax.set_title(getFigureBehaviouralEventsLabels(event), loc='right', FontSize=18)
        ax.text(-2, ylimSpeed[event], letterList[m], fontsize=20, horizontalalignment='center', color='black', weight='bold')

        if row == 0:
            ax = axes[row, col]
            ax.hlines(0, 0.8, 2.2, color=getColorAge('5we'), lw=30)
            ax.text(1.5, ylimSpeed[event] * 0.01, "5 we", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 3.8, 5.2, color=getColorAge('3mo'), lw=30)
            ax.text(4.5, ylimSpeed[event] * 0.01, "3 mo", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 6.8, 8.2, color=getColorAge('7mo'), lw=30)
            ax.text(7.5, ylimSpeed[event] * 0.01, "7 mo", FontSize=12, horizontalalignment='center', color='white')

            ax = axes[row, col + 1]
            ax.hlines(0, 0.8, 2.2, color=getColorAge('5we'), lw=30)
            ax.text(1.5, ylimDuration[event] * 0.01, "5 we", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 3.8, 5.2, color=getColorAge('3mo'), lw=30)
            ax.text(4.5, ylimDuration[event] * 0.01, "3 mo", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 6.8, 8.2, color=getColorAge('7mo'), lw=30)
            ax.text(7.5, ylimDuration[event] * 0.01, "7 mo", FontSize=12, horizontalalignment='center', color='white')

        m += 1

    plt.tight_layout()
    plt.show()
    fig.savefig("Fig_Speed_Duration_Mean_{}.pdf".format(sex), dpi=300)
    fig.savefig("Fig_Speed_Duration_Mean_{}.jpg".format(sex), dpi=300)
    #fig.savefig("Fig_Speed_Duration_Mean.svg", dpi=100)

    # plt.show()

    print("Job done.")


def getColorLine( stars ):
    colorLine = 'whitesmoke'
    
    if stars == '*':
        colorLine = 'lightgrey'
    
    if stars == '**':
        colorLine = 'darkgrey'
        
    if stars == '***':
        colorLine = 'black'
    
    return colorLine

def plotMeanSpeedDurationGenoFast(jsonFile='durationSpeedData_all_classes.json', age='3mo', sex='female', eventListToTest=None):
    # open the json file to have the data
    with open(jsonFile) as json_data:
        durationSpeedData = json.load(json_data)
    print("json file re-imported.")
    letterList = list(string.ascii_uppercase)
    durationSpeedListWithOverlap = durationSpeedData['with overlap']
    durationSpeedListWithoutOverlap = durationSpeedData['no overlap']

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15, 10), sharey=False, sharex=True)

    ylimSpeed = {'longChase': 15, 'FollowZone Isolated': 25, 'Approach contact': 15,
                 'Break contact': 40, 'Contact': 50, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 60}
    ylimDuration = {'longChase': 800, 'FollowZone Isolated': 8, 'Approach contact': 8,
                    'Break contact': 30, 'Contact': 1000, 'Get away': 50, 'Oral-genital Contact': 50, 'Train2': 50}

    colNb = [0, 2, 0, 2]
    rowNb = [0, 0, 1, 1]

    k = 0

    for event in eventListToTest:
        print('event: ', event)
        print('k: ', k)
        col = colNb[k]
        row = rowNb[k]

        dataToTest = {}
        for var in ['speed', 'duration']:
            dataToTest[var] = {}
            dataToTest[var][sex] = {}
            for strain in ['C57BL/6J', 'Shank3']:
                dataToTest[var][sex][strain] = {}
                for condition in ['usv', 'no usv']:
                    dataToTest[var][sex][strain][condition] = []

        print('col: ', col, 'row: ', row)

        x1Dic = {'C57BL/6J': 1, 'Shank3': 4}
        x2Dic = {'C57BL/6J': 2, 'Shank3': 5}

        speedUsv = {}
        speedNoUsv = {}
        durUsv = {}
        durNoUsv = {}

        for strain in ['C57BL/6J', 'Shank3']:
            if strain == 'C57BL/6J':
                colorPoints = getColorWT()
            elif strain == 'Shank3':
                colorPoints = getColorKO()

            listSpeedToComplement = dataToTest['speed'][sex][strain]
            listDurationToComplement = dataToTest['duration'][sex][strain]
            x1 = x1Dic[strain]
            x2 = x2Dic[strain]
            for experiment in durationSpeedListWithOverlap[event][strain][sex][age].keys():
                for animal in sorted(durationSpeedListWithOverlap[event][strain][sex][age][experiment].keys()):

                    if len(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal]) > 0:

                        speedUsv[animal] = []
                        speedNoUsv[animal] = []
                        durUsv[animal] = []
                        durNoUsv[animal] = []

                        for i in range(len(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal])):

                            if durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][1] <= 70:
                                speedUsv[animal].append(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][1])

                            durUsv[animal].append(durationSpeedListWithOverlap[event][strain][sex][age][experiment][animal][i][0])

                        for j in range(len(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal])):

                            if durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][1] <= 70:
                                speedNoUsv[animal].append(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][1])

                            durNoUsv[animal].append(durationSpeedListWithoutOverlap[event][strain][sex][age][experiment][animal][j][0])

                        Wspeed, Pspeed = mannwhitneyu(x=speedUsv[animal], y=speedNoUsv[animal], alternative='greater')
                        starsSpeed = getStarsFromPvalues(Pspeed, 24)
                        print("Speed with USV: ", len(speedUsv[animal]), "no USV: ", len(speedNoUsv[animal]), "U=", Wspeed,
                              " ", "p=", Pspeed, " ", starsSpeed)

                        colorLineSpeed = getColorLine(starsSpeed)

                        Wdur, Pdur = mannwhitneyu(x=durUsv[animal], y=durNoUsv[animal], alternative='greater')
                        starsDur = getStarsFromPvalues(Pdur, 24)
                        print("Duration with USV: ", len(durUsv[animal]), "no USV: ", len(durNoUsv[animal]), "U=", Wdur,
                              " ", "p=", Pdur, " ", starsDur)

                        colorLineDur = getColorLine(starsDur)

                        # speed
                        ax = axes[row][col]
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.set_xticks([1, 2, 4, 5, 7, 8])
                        ax.set_xticklabels(["USV", "no USV", "USV", "no USV", "USV", "no USV"], rotation=45, FontSize=15)
                        ax.xaxis.set_tick_params(direction="in")
                        ax.set_xlim(0, 9)
                        ax.set_ylim(0, ylimSpeed[event])
                        ax.set_ylabel("mean speed (cm/s)", FontSize=15)
                        ax.tick_params(axis='y', labelsize=14)


                        ax.scatter(x1, np.mean(speedUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.scatter(x2, np.mean(speedNoUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.plot([x1, x2], [np.mean(speedUsv[animal]), np.mean(speedNoUsv[animal])], zorder=0,
                                color=colorLineSpeed, linewidth=0.6)

                        semUsvSpeed = np.std(speedUsv[animal]) / np.sqrt(len(speedUsv[animal]))
                        semNoUsvSpeed = np.std(speedNoUsv[animal]) / np.sqrt(len(speedNoUsv[animal]))
                        # sem with USV
                        ax.vlines(x1, np.mean(speedUsv[animal]), np.mean(speedUsv[animal]) + semUsvSpeed, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedUsv[animal]) + semUsvSpeed, x1 - 0.1, x1 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x1, np.mean(speedUsv[animal]), np.mean(speedUsv[animal]) - semUsvSpeed, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedUsv[animal]) - semUsvSpeed, x1 - 0.1, x1 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        # sem no USV
                        ax.vlines(x2, np.mean(speedNoUsv[animal]), np.mean(speedNoUsv[animal]) + semNoUsvSpeed,
                                  linestyle='solid', linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedNoUsv[animal]) + semNoUsvSpeed, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x2, np.mean(speedNoUsv[animal]), np.mean(speedNoUsv[animal]) - semNoUsvSpeed,
                                  linestyle='solid', linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(speedNoUsv[animal]) - semNoUsvSpeed, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)

                        listSpeedToComplement['usv'].append(np.mean(speedUsv[animal]))
                        listSpeedToComplement['no usv'].append(np.mean(speedNoUsv[animal]))

                        # duration
                        ax = axes[row][col + 1]
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.set_xticks([1, 2, 4, 5])
                        ax.set_xticklabels(["USV", "no USV", "USV", "no USV"], rotation=45, FontSize=15)
                        ax.xaxis.set_tick_params(direction="in")
                        ax.set_xlim(0, 6)
                        ax.set_ylim(0, ylimDuration[event])
                        ax.set_ylabel("duration (frames)", FontSize=15)
                        ax.tick_params(axis='y', labelsize=14)

                        ax.scatter(x1, np.mean(durUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.scatter(x2, np.mean(durNoUsv[animal]), marker='o', s=16, c=colorPoints)
                        ax.plot([x1, x2], [np.mean(durUsv[animal]), np.mean(durNoUsv[animal])], zorder=0,
                                color=colorLineDur, linewidth=1)

                        semUsvDur = np.std(durUsv[animal]) / np.sqrt(len(durUsv[animal]))
                        semNoUsvDur = np.std(durNoUsv[animal]) / np.sqrt(len(durNoUsv[animal]))
                        # sem with USV
                        ax.vlines(x1, np.mean(durUsv[animal]), np.mean(durUsv[animal]) + semUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durUsv[animal]) + semUsvDur, x1 - 0.1, x1 + 0.1, linewidth=0.2, color=colorPoints)
                        ax.vlines(x1, np.mean(durUsv[animal]), np.mean(durUsv[animal]) - semUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durUsv[animal]) - semUsvDur, x1 - 0.1, x1 + 0.1, linewidth=0.2, color=colorPoints)
                        # sem no USV
                        ax.vlines(x2, np.mean(durNoUsv[animal]), np.mean(durNoUsv[animal]) + semNoUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durNoUsv[animal]) + semNoUsvDur, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)
                        ax.vlines(x2, np.mean(durNoUsv[animal]), np.mean(durNoUsv[animal]) - semNoUsvDur, linestyle='solid',
                                  linewidth=0.2, color=colorPoints)
                        ax.hlines(np.mean(durNoUsv[animal]) - semNoUsvDur, x2 - 0.1, x2 + 0.1, linewidth=0.2,
                                  color=colorPoints)

                        listDurationToComplement['usv'].append(np.mean(durUsv[animal]))
                        listDurationToComplement['no usv'].append(np.mean(durNoUsv[animal]))

                    else:
                        print("Values: no data")

        stars = {}

        for var in ['speed', 'duration']:
            for strain in ['C57BL/6J', 'Shank3']:
                print("*************** ", event, " -var: ", var, " -strain: ", strain)
                T, P = wilcoxon(x=dataToTest[var][sex][strain]['usv'], y=dataToTest[var][sex][strain]['no usv'], alternative="greater")
                stars[var, strain] = getStarsFromPvalues(P, 1)
                print("Wilcoxon ", var, " ", strain, " : ", "T=", T, " ", "p=", P, " ", stars[var, strain])

        ax = axes[row][col]
        ax.text(1.5, ylimSpeed[event] * 0.95, stars['speed', 'C57BL/6J'], fontsize=12, horizontalalignment='center',
                color='black')
        ax.text(4.5, ylimSpeed[event] * 0.95, stars['speed', 'Shank3'], fontsize=12, horizontalalignment='center',
                color='black')
        ax.text(-2, ylimSpeed[event], letterList[k], fontsize=20, horizontalalignment='center', color='black',
                weight='bold')


        ax = axes[row][col + 1]
        ax.text(1.5, ylimDuration[event] * 0.95, stars['duration', 'C57BL/6J'], fontsize=12, horizontalalignment='center',
                color='black')
        ax.text(4.5, ylimDuration[event] * 0.95, stars['duration', 'Shank3'], fontsize=12, horizontalalignment='center',
                color='black')

        k += 1

    m = 0
    for event in eventListToTest:
        print('event: ', event)
        print('m: ', m)
        col = colNb[m]
        row = rowNb[m]

        print('col: ', col, 'row: ', row)

        ax = axes[row, col]
        ax.set_title(getFigureBehaviouralEventsLabels(event), loc='right', FontSize=18)

        if row == 0:
            ax = axes[row, col]
            ax.hlines(0, 0.8, 2.2, color=getColorWT(), lw=30)
            ax.text(1.5, ylimSpeed[event] * 0.01, "B6", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 3.8, 5.2, color=getColorKO(), lw=30)
            ax.text(4.5, ylimSpeed[event] * 0.01, "KO", FontSize=12, horizontalalignment='center', color='white')

            ax = axes[row, col + 1]
            ax.hlines(0, 0.8, 2.2, color=getColorWT(), lw=30)
            ax.text(1.5, ylimDuration[event] * 0.01, "B6", FontSize=12, horizontalalignment='center', color='white')
            ax.hlines(0, 3.8, 5.2, color=getColorKO(), lw=30)
            ax.text(4.5, ylimDuration[event] * 0.01, "KO", FontSize=12, horizontalalignment='center', color='white')

        m += 1

    plt.tight_layout()
    plt.show()
    fig.savefig("Fig_Speed_Duration_Mean_Shank3.pdf", dpi=300)
    fig.savefig("Fig_Speed_Duration_Mean_Shank3.jpg", dpi=300)

    # plt.show()

    print("Job done.")


if __name__ == '__main__':
    '''
    This codes allows to test whether events are different in speed and duration if they occur with USVs or if they occur without USVs.
    '''
    print("Code launched.")
    
    # set font
    from matplotlib import rc    
    rc('font',**{'family':'serif','serif':['Arial']})
    
    eventListToTest = ["Contact", "Approach contact", "Break contact", "Get away", "FollowZone Isolated", "Oral-genital Contact", "Train2", "longChase"]
    behavEventListShort = ["Contact", "Approach contact", "Break contact", "FollowZone Isolated", "Oral-genital Contact", "Train2"]
    behavEventListShortPoints = ["Train2", "FollowZone Isolated", "Approach contact", "Break contact"]
    
    while True:
        
        question = "Do you want to:"
        question +="\n\t [c] compute the number of events overlapping with USVs"
        question +="\n\t [pd] violinplot figure for duration of events per ind"
        question +="\n\t [ps] violinplot figure for speed of events per ind"
        question +="\n\t [pp] plot scatterpoint figure for WT 3 age classes for speed and duration (mean per ind)"
        question +="\n\t [ppko] plot scatterpoint figure for Shank3 and WT for speed and duration (mean per ind)"
        question +="\n"
        answer = input("Action:")
        
        if answer=="c":
            #computeEventsOverlappingWithUsv( strain='C57BL/6J', age = '3mo', sex = 'male', eventListToTest=behavEventListShort )
            tmin, tmax = getMinTMaxTInput()
            strainList = ['C57BL/6J', 'Shank3']
            sexList = ['male', 'female']
            ageList = ['5we', '3mo', '7mo']

            computeEventsOverlappingWithUsvFast(strain=None, age=None, sex=None, tmin=tmin, tmax=tmax, strainList=strainList, sexList=sexList, ageList=ageList, eventListToTest=behavEventListShort)

            break

        if answer=="pd":
            plotViolinDuration( strain='C57BL/6J', age = None, sex = 'female', eventListToTest=behavEventListShort )
            break
        
        if answer=="ps":
            plotViolinSpeed( strain='C57BL/6J', age = None, sex = 'female', eventListToTest=behavEventListShort )
            break
        
        if answer=="pp":
            #plotMeanSpeedDuration( strain='C57BL/6J', age = None, sex = 'female', eventListToTest=behavEventListShortPoints )
            plotMeanSpeedDurationFast(jsonFile='durationSpeedData_all_classes.json', strain='C57BL/6J', sex='female', eventListToTest=behavEventListShortPoints)
            break
        
        if answer=="ppko":
            #plotMeanSpeedDurationGeno( age = '3mo', sex = 'female', eventListToTest=behavEventListShortPoints )
            plotMeanSpeedDurationGenoFast(jsonFile='durationSpeedData_all_classes.json', age='3mo', sex='female',
                                      eventListToTest=behavEventListShortPoints)

            break
    
    print("All done.")
    
    
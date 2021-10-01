'''
Created on 28 avr. 2020

@author: Fab
'''
import matplotlib.pyplot as plt

# colors per age:

colorWT = "#7f78d2"
colorKO = "#33cc33"

color5we = "#efB1ff"
color3mo = colorWT
color7mo = "#481380"

colorAgeList = [ color5we, color3mo, color7mo ]

def getColorGeno( geno ):
    if geno=="WT-WT":
        return 'steelblue'
    if geno=="KO-KO":
        return 'darkorange'
    if geno=="Del/+-Del/+":
        return 'darkorange'

def getColorAge( age ):
    if age=="5we":
        return color5we
    if age=="3mo":
        return color3mo
    if age=="7mo":
        return color7mo

def getColorEvent( event ):
    if event == 'Stop isolated':
        return 'dodgerblue'
    if event == 'Break contact':
        return 'skyblue'
    if event == 'Approach contact':
        return 'darkmagenta'
    if event == 'Oral-oral Contact':
        return 'indianred'
    if event == 'Oral-genital Contact':
        return 'darkorange'
    if event == 'Side by side Contact':
        return 'goldenrod'
    if event == 'Side by side Contact, opposite way':
        return 'darkgoldenrod'
    if event == 'FollowZone Isolated':
        return 'darkolivegreen'
    if event == 'Train2':
        return 'ForestGreen'

def getColorWT():
    return colorWT

def getColorKO():
    return colorKO

def getFigureBehaviouralEvents( longList = False , withUrinate = False ):
    
    behaviouralEvents = []
    
    if not longList:
    
        behaviouralEvents = [ "Stop isolated", 
                             "Break contact",
                             #"Approach rear",
                             "Approach contact",
                             "Oral-oral Contact",
                             "Oral-genital Contact",
                             "Side by side Contact",
                             "Side by side Contact, opposite way",
                             "FollowZone Isolated",
                             "Train2"
                             #, "longChase"
                              ]
    else:
        
        behaviouralEvents = [ "Stop isolated",
                             "Move isolated", 
                             "Break contact",
                             "Get away",
                             "Social approach",
                             #"Approach rear",
                             "Approach contact",
                             "Contact",
                             "Oral-oral Contact",
                             "Oral-genital Contact",
                             "Side by side Contact",
                             "Side by side Contact, opposite way",
                             "seq oral oral - oral genital",
                             "seq oral geni - oral oral",   
                             "FollowZone Isolated",
                             "Train2"
                             #, "longChase"
                              ]

    
    if withUrinate:
        behaviouralEvents.append("Urinate USV")
    
    return behaviouralEvents


def getFigureBehaviouralEventsLabels(event):
    behaviouralEventsLabels = {"Stop isolated": 'single idle',
                         "Move isolated": 'single move',
                         "Break contact": 'break contact',
                         "Get away": 'get away',
                         "Social approach": 'approach social range',
                         # "Approach rear": 'approach reared mouse',
                         "Approach contact": 'approach contact',
                         "Contact": 'contact',
                         "Oral-oral Contact": 'nose-nose',
                         "Oral-genital Contact": 'nose-anogenital',
                         "Side by side Contact": 'side-side',
                         "Side by side Contact, opposite way": 'side-side, head-to-tail',
                         "seq oral oral - oral genital": 'nose-nose & nose-anogenital',
                         "seq oral geni - oral oral": 'nose-anogenital & nose-nose',
                         "FollowZone Isolated": 'follow',
                         "Train2": 'train2'
                         # , "longChase": 'long chase'
                         }
    return behaviouralEventsLabels[event]



def getFigureVocTraits():

    vocTraits = [ "durationMs","startFrequencyHz", "endFrequencyHz","diffStartEndFrequencyHz",
                    "minFrequency","maxFrequency", "meanFrequencyHz",
                  "frequencyDynamicHz","frequencyTVHz", "linearityIndex",
                   "nbModulation","nbJump","meanPower","peakPower",
                   "slope","griffIndex"                  
                   ]
        
    return vocTraits


def getFigureLabelTrait(vocTrait):
    labelTraits = {"durationMs": 'duration (ms)', "startFrequencyHz": 'start freq (Hz)', "endFrequencyHz": 'end freq (Hz)', "diffStartEndFrequencyHz": 'start freq - end freq (Hz)',
                 "minFrequency": 'min freq (Hz)', "maxFrequency": 'max freq (Hz)', "meanFrequencyHz": 'mean freq (Hz)',
                 "frequencyDynamicHz": 'freq range (Hz)', "frequencyTVHz": 'freq total variation (Hz)', "linearityIndex": 'linearity index',
                 "nbModulation": 'freq modulations (nb)', "nbJump": 'freq jumps (nb)', "meanPower": 'mean power (AU)', "peakPower": 'peak power (AU)',
                 "slope": 'slope', "griffIndex": 'harsh index'}

    label = labelTraits[vocTrait]
    return label


def getFigureBurstTraits():

    return [ "durationMs", "nbUSV", "meanDuration", "stdDuration", "meanInterval","stdInterval","meanFrequency" ]


def getFigureLabelBurstTraits(burstTrait):
    labelBurstTrait = {"durationMs": 'duration (ms)', "nbUSV": 'nb USVs / burst', "nbUsvBurst": 'nb USVs / burst', "meanDuration": 'mean duration (ms)', "stdDuration": 'std duration (ms)', "meanInterval": 'mean interval (ms)', "stdInterval": 'std interval (ms)', "meanFrequency": 'mean freq (Hz)'}
    label = labelBurstTrait[burstTrait]
    return label


def getPaperColor( index, totalNumberOfColor , name='Set2' ): #name='hsv'
    
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, totalNumberOfColor+1)(index)        

    

if __name__ == '__main__':
    pass
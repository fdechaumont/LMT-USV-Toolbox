'''
Created on 7 mai 2020

@author: Fabrice de Chaumont

'''
from random import random

import json


colorAge = ["#C6D4E1", "#BDB8AD", "#44749D"]


ageList = ['5we', '3mo', '7mo']
sexList = ['male', 'female']
strainList = ['C57BL/6J', 'Shank3']

acousticVariables = [ "vocNumber", "durationMs", "startFrequencyHz", "endFrequencyHz", "diffStartEndFrequencyHz",
                         "peakFrequency", "meanFrequencyHz", "minFrequency", "maxFrequency", "frequencyDynamicHz", "frequencyTVHz", "linearityIndex",
                         "nbModulation", "nbJump", "peakPower", "meanPower", "minPower", 'griffIndex', 'slope', 'slopeNormalized'
                   ]
acousticVariablesShort = [ "durationMs", "startFrequencyHz", "endFrequencyHz", "diffStartEndFrequencyHz",
                         "peakFrequency", "meanFrequencyHz", "minFrequency", "maxFrequency", "frequencyDynamicHz", "frequencyTVHz", "linearityIndex",
                         "nbModulation", "nbJump", "peakPower", "meanPower", "minPower"
                   ]

yMinVar = {"numberOfVoc": 0, "durationMs": 0, "startFrequencyHz": 0, "endFrequencyHz": 0,
           "diffStartEndFrequencyHz": -40000,
           "peakFrequency": 0, "meanFrequencyHz": 0, "minFrequency": 0, "maxFrequency": 0,
           "frequencyDynamicHz": 0, "frequencyTVHz": 0, "linearityIndex": 0,
           "nbModulation": 0, "nbJump": 0, "peakPower": 0, "meanPower": 0, "minPower": 0, "griffIndex": 0,
           "slope": -80000, "slopeNormalized": 0}

yMaxVar = {"numberOfVoc": 100, "durationMs": 400, "startFrequencyHz": 150000, "endFrequencyHz": 150000,
           "diffStartEndFrequencyHz": 50000,
           "peakFrequency": 250000, "meanFrequencyHz": 150000, "minFrequency": 150000, "maxFrequency": 150000,
           "frequencyDynamicHz": 80000, "frequencyTVHz": 200000, "linearityIndex": 50,
           "nbModulation": 12, "nbJump": 4, "peakPower": 20, "meanPower": 4, "minPower": 200, "griffIndex": 40,
           "slope": 80000, "slopeNormalized": 20}

yMinVarAge = {"numberOfVoc": 0, "durationMs": 0, "startFrequencyHz": 30000, "endFrequencyHz": 30000,
              "diffStartEndFrequencyHz": -25000,
              "peakFrequency": 30000, "meanFrequencyHz": 30000, "minFrequency": 30000, "maxFrequency": 30000,
              "frequencyDynamicHz": 0, "frequencyTVHz": 0, "linearityIndex": 0,
              "nbModulation": 0, "nbJump": 0, "peakPower": 0, "meanPower": 0, "minPower": 0, "griffIndex": 0,
              "slope": -45000, "slopeNormalized": 0}

yMaxVarAge = {"numberOfVoc": 100, "durationMs": 200, "startFrequencyHz": 100000, "endFrequencyHz": 100000,
              "diffStartEndFrequencyHz": 25000,
              "peakFrequency": 250000, "meanFrequencyHz": 100000, "minFrequency": 100000, "maxFrequency": 100000,
              "frequencyDynamicHz": 40000, "frequencyTVHz": 160000, "linearityIndex": 30,
              "nbModulation": 6, "nbJump": 2, "peakPower": 15, "meanPower": 3, "minPower": 200, "griffIndex": 16,
              "slope": 45000, "slopeNormalized": 20}


def getJsonFile( file ):
    return file[:-2]+"json"

def loadResult( file ):
    '''
    Call with __file__ as argument
    '''
    print("Loading json...")
    try:
        with open( getJsonFile( file ) ) as json_data:
            result = json.load(json_data)
        print("json file loaded.")
        return result
    except:
        print("Can't load result")
        result = {}
        return result

def saveResult( file, result ):
    '''
    Call with __file__ as argument
    '''
    with open( getJsonFile( file ) , 'w') as f:
        json.dump(result, f, indent=4 )
    
def getStarsFromPvalues(pvalue, numberOfTests):
    stars = "ns"
    
    s1 = 0.05/numberOfTests
    s2 = 0.01/numberOfTests
    s3 = 0.001/numberOfTests
    
    if pvalue < s3:
        stars = "***"
    if pvalue >= s3 and pvalue < s2:
        stars = "**"
    if pvalue >= s2 and pvalue < s1:
        stars = "*"
    
    return stars

def addJitter(x, jit):
        
    newX = []
    for item in x:
        addedJitter = (random()*2-1)*jit
        newX.append(item+addedJitter)
    
    return newX
    
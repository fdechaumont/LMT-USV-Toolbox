'''
Created on 7 mai 2020

@author: Fabrice de Chaumont

'''
from random import random

import json

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
    
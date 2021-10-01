'''
Created on 19 févr. 2020

@author: Fab
'''

from os import path


class Experiment(object):
    '''
    describe an experiment
    '''
    
    def __init__(self, genotype, strain, age, sex, name, file ):
        
        self.genotype = genotype
        self.strain = strain
        self.age = age
        self.sex = sex
        self.name = name
        
        self.file = file
        
        if not path.exists( self.file ):
            
            fileD = "d:"+ file[2:]            
            if path.exists( fileD ):
                self.file = fileD                

            fileE = "e:"+ file[2:]
            if path.exists( fileE ):
                self.file = fileE
            
            fileE = "f:"+ file[2:]
            if path.exists( fileE ):
                self.file = fileE
            
            fileE = "g:"+ file[2:]
            if path.exists( fileE ):
                self.file = fileE
            
        
        
        
    def __str__(self):
        
        s = "Experiment "
        s+= self.genotype
        s+="\t"
        s+= self.strain
        s+="\t"
        s+= self.age
        s+="\t"
        s+= self.sex
        s+="\t"
        s+= self.name
        s+="\t"
        s+= self.file
        return s        
    
    def getFullName(self):
        
        fullName = self.sex + " " + self.age + " " + self.strain + " " + self.genotype + " " + self.name 
        return fullName
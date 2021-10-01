'''
Created on 19 févr. 2020

@author: Fab
'''
from LMT.USV2.experimentList.Experiment import Experiment

from os import path


def getExperimentForFile( file ):
    experiments = getExperimentList()
    for experiment in experiments:
        if experiment.file == file:
            return experiment
    
    return None

def checkIfAllExperimentFilesExist():
    
    experiments = getExperimentList()
    
    for experiment in experiments:
        print( experiment )
        exist = path.exists( experiment.file )
        if not exist:
            print("********** ERROR DOES NO EXIST")
        else:
            print("Ok.")

def getExperimentList( genotype=None, strain=None, age=None, sex=None ):

    all = getAllExperimentList()
    
    for experiment in list(all):
        
        if strain!=None:
            
            if experiment.strain != strain:
                all.remove( experiment )
                continue

        if genotype!=None:
            
            if experiment.genotype != genotype:
                all.remove( experiment )
                continue

        if age!=None:
            
            if experiment.age != age:
                all.remove( experiment )
                continue

        if sex!=None:
            
            if experiment.sex != sex:
                all.remove( experiment )
                continue

    return all


def getExperimentWTvsKO():
    
    experimentsKO = getExperimentList( genotype="KO" )
    experimentsWT = getExperimentList( genotype="WT", age="3mo", sex="female" )

    return experimentsWT, experimentsKO

def getAllExperimentList():

    wt = "WT"
    ko = "KO"
    
    #strain
    shank3 = "Shank3"    
    bl6 = "C57BL/6J"
    
    age5we = "5we"
    age3mo = "3mo"
    age7mo = "7mo"
    
    f = "female"
    m = "male"
      
    experimentList = [

    Experiment(wt,bl6,age5we,f,"F13-F14","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181204_usv_lmt_paired_F13_F14_5we_Experiment 8959/usv_F13_F14/ch1_F13_F14/20181204_usv_lmt_pair_F13_F14_5we.sqlite"),
    Experiment(wt,bl6,age5we,f,"F15-F16","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181204_usv_lmt_pair_F15_F16_5we_Experiment 2797/usv_F15_F16/20181204_usv_lmt_pair_F15_F16_5we.sqlite"),
    Experiment(wt,bl6,age5we,f,"F17-F18","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181207_usv_lmt_pair_F17_F18_5we_Experiment 8484/usv_011_F17-F18/ch1/20181207_usv_lmt_pair_F17_F18_5we.sqlite"),
    Experiment(wt,bl6,age5we,f,"F19-F20","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181217_usv_lmt_pair_F19_F20_5we/ch1_F19_F20/20181217_usv_lmt_pair_F19_F20_5we.sqlite"),
    
    Experiment(wt,bl6,age5we,m,"M40-M41","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181211_usv_lmt_pair_M40_M41_5we/ch1/20181211_usv_lmt_pair_M40_M41_5we.sqlite"),
    Experiment(wt,bl6,age5we,m,"M42-M43","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181211_usv_lmt_pair_M42_M43_5we/usv_010/ch1/20181211_usv_lmt_pair_M42_M43_5we.sqlite"),
    Experiment(wt,bl6,age5we,m,"M44-M45","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181214_usv_lmt_pair_M44_M45_5we/ch1_M44_M45/20181214_usv_lmt_pair_M44_M45_5we.sqlite"),
    Experiment(wt,bl6,age5we,m,"M46-M47","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20181214_usv_lmt_pair_M46_M47_5we/ch1_M46_M47/20181214_usv_lmt_pair_M46_M47_5we.sqlite"),
    
    Experiment(wt,bl6,age3mo,m,"M40-M41","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190118_usv_lmt_pair_M40-M41_2mo_1/ch1/20190118_usv_lmt_pair_M40-M41_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,m,"M42-M43","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190122_usv_lmt_pair_M42_M43_2mo_Experiment 7259/ch1/20190122_usv_lmt_pair_M42_M43_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,m,"M44-M45","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190125-usv_lmt_pair_M44_M45_2mo_1/usv_lmt004/ch1/20190125-usv_lmt_pair_M44_M45_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,m,"M46-M47","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190129_usv_lmt_pair_M46_M47_2mo_2_1/ch1/20190129_usv_lmt_pair_M46_M47_2mo.sqlite"),
    
    Experiment(wt,bl6,age3mo,f,"F13-F14","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190125_usv_lmt_pair_F13-F14_2mo_1/usv_021/F13-F14_2mo_ch1/20190125_usv_lmt_pair_F13-F14_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,f,"F15-F16","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190129_usv_lmt_pair_F15_F16_1/F15-F16_2mo_ch1/20190129_usv_lmt_pair_F15_F16_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,f,"F17-F18","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190201_usv_lmt_pair_F17_F18_2mo/F17-F18_2mo_ch1/20190201_usv_lmt_pair_F17_F18_2mo.sqlite"),
    Experiment(wt,bl6,age3mo,f,"F19-F20","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/20190201_usv_lmt_pair_F19_F20_2mo/F19-F20_2mo_ch1/20190201_usv_lmt_pair_F19_F20_2mo.sqlite"),
    
    Experiment(wt,bl6,age7mo,m,"M40-M41","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190611_usv_lmt_pair_M40_M41_8mo/20190611_usv_lmt_pair_M40_M41_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,m,"M42-M43","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190611_usv_lmt_pair_M42_M43_8mo/20190611_usv_lmt_pair_M42_M43_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,m,"M44-M45","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190614_usv_lmt_pair_M44_M45_8mo/20190614_usv_lmt_pair_M44_M45_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,m,"M46-M47","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190614_usv_lmt_pair_M46_M47_8mo/20190614_usv_lmt_pair_M46_M47_8mo.sqlite"),
    
    Experiment(wt,bl6,age7mo,f,"F13-F14","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190621_usv_lmt_pair_F13_F14_8mo/20190621_usv_lmt_pair_F13_F14_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,f,"F15-F16","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190618_usv_lmt_pair_F15_F16_8mo/20190618_usv_lmt_pair_F15_F16_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,f,"F17-F18","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190621_usv_lmt_pair_F17_F18_8mo/20190621_usv_lmt_pair_F17_F18_8mo.sqlite"),
    Experiment(wt,bl6,age7mo,f,"F19-F20","D:/USV_spontaneous_rec_pairs/all_data_usv_pairs/data_7mo/20190621_usv_lmt_pair_F19_F20_8mo/20190621_usv_lmt_pair_F19_F20_8mo.sqlite"),
              
    Experiment(ko,shank3,age3mo,f,"4724356-4724113","D:/usv_shank3/20190510_usv_lmt_pair_shank3_4724356_4724113/20190510_usv_lmt_pair_shank3_4724356_4724113.sqlite" ),
    Experiment(ko,shank3,age3mo,f,"4849069-4849197","D:/usv_shank3/20191014_usv_lmt_pair_shank3_4849069_4849197/20191014_usv_lmt_pair_shank3_4849069_4849197.sqlite" ),    
    Experiment(ko,shank3,age3mo,f,"4849144-4849294","D:/usv_shank3/20191014_usv_lmt_pair_shank3_4849144_4849294/20191014_usv_lmt_pair_shank3_4849144_4849294.sqlite" ),
    Experiment(ko,shank3,age3mo,f,"4849155_4849457","D:/usv_shank3/20191121_shank3_usv_lmt_pair_4849155_4849457_Experiment 3323/20191121_shank3_usv_lmt_pair_4849155_4849457_Experiment 3323.sqlite" ), 
    Experiment(ko,shank3,age3mo,f,"4849484_4849378","D:/usv_shank3/20200225_usv_lmt_shank3_pair_4849484_4849378_Experiment 4992/20200225_usv_lmt_shank3_pair_4849484_4849378_Experiment 4992.sqlite" ),
    Experiment(ko,shank3,age3mo,f,"4417496_4419810","D:/usv_shank3/20200306_usv_lmt_shank3_pair_4417496_4419810_Experiment 2910/20200306_usv_lmt_shank3_pair_4417496_4419810_Experiment 2910.sqlite" ) 


    ]
    

    return experimentList

if __name__ == '__main__':
    pass




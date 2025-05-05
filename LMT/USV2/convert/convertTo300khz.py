'''
Created on 30 nov. 2020

@author: Fab
'''
import wave
import librosa
import soundfile as sf
from pathlib import Path
import os
        
from tkinter.filedialog import askopenfilenames

def convert( file ):
    
    print("---> Processing file: " + file )
    print( "Checking sampling rate.")
    try:
        wav = wave.open( file, 'r')        
        #frames = wav.readframes(-1)
        #frames = wav.readframes(300000*3)
        #self.__wavBuffer = pylab.frombuffer(frames, 'int16')        
        frameRate = wav.getframerate()
        print( "framerate: " , frameRate )
        wav.close()
    except:
        print("Unknown sound format.")
        return    

    
    if frameRate == 300000:
        print("File is already in 300kHz: " , file , "... (converting anyway)")
                                 
    print("Converting sampling rate to 300kHz...")
    
    try:
        duration = librosa.get_duration(filename=file)
        print("Duration: " , duration , " seconds")
        y, sr = librosa.load( file , mono=True, sr=300000)
        
        dir = os.path.dirname( file )
        fileName = os.path.basename(file)
                
        convertedFile = f"{dir}/converted/{fileName}"
        print( convertedFile )
        Path(dir+"/converted/").mkdir(parents=True, exist_ok=True)
        #quit()
                         
        sf.write( convertedFile, y, sr, subtype='PCM_16')
        print("Converting framerate done")
    
    except:
        print("Error in conversion")                
        return
        
    # end convert

if __name__ == '__main__':
    
    print("Please provide wav file to convert.")
    print("Original files will be kept. converted files will be located in the subdirectory 'converted'")
    filenames = set( askopenfilenames() )
    
    for file in filenames:
        convert( file )
    
    print("---> Done.")
    
    
from code.a_CreateClimate import MainA
from code.b_CreateGlaciers import MainB 
from code.c_CombineData import MainC

def Main():
    MainA() # Create climate data
    MainB() # Create glacier data
    MainC() # Combine climate, glacier and stream data

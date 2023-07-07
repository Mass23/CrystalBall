from code.a_CreateClimate import MainA
from code.b_CreateGlaciers import MainB 
from code.c_CombineData import MainC

def Main():
    print('Running a_CreateClimate.py ...')
    MainA() # Create climate data

    print('Running b_CreateGlaciers.py ...')
    MainB() # Create glacier data

    print('Running c_CombineData.py ...')
    MainC() # Combine climate, glacier and stream data

if __name__ == "__main__":
    Main()
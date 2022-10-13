import sys
import clr
import time
import numpy as np
import os

# Add Newport.DLS.CommandInterface.dll in References of your script
sys.path.append(r'C:\\Windows\\Microsoft.NET\\assembly\\GAC_64\\Newport.DLS.CommandInterface\\v4.0_1.0.1.0__90ac4f829985d2bf')
    
#Add library path
clr.AddReference("C:\\Windows\\Microsoft.NET\\assembly\\GAC_64\\Newport.DLS.CommandInterface\\v4.0_1.0.1.0__90ac4f829985d2bf\\Newport.DLS.CommandInterface.dll")

from CommandInterfaceDLS import *
from System import String, Double

dummy_out0 = Double(0.)
dummy_out1 = String('')
myDLS = DLS()


class DLS_cmd():
    def __init__(self,port):
        try:
            myDLS.OpenInstrument(port)
        except:
            print("could not connect to newport delay line, check connection and power on the device")

    #delay line function
    def update(self):
        current_pos=myDLS.TP(dummy_out0,dummy_out1)[1] # current position
        print("current position is:", current_pos)
        current_vel=myDLS.VA_Get(dummy_out0,dummy_out1)[1] # current velocity
        print("current velocity is:", current_vel)
        current_acc=myDLS.AC_Get(dummy_out0,dummy_out1)[1] # current acceleration
        print("current acceleration is:", current_acc)

    def current_pos(self):
        cur_pos = myDLS.TP(dummy_out0,dummy_out1)[1] # current position
        return cur_pos

    def current_vel(self):
        cur_vel = myDLS.VA_Get(dummy_out0,dummy_out1)[1] # current velocity
        return cur_vel

    def current_acc(self):
        cur_acc = myDLS.AC_Get(dummy_out0,dummy_out1)[1] # current acceleration
        return cur_acc

    def initialize_DLS(self):
        my_dummy = myDLS.IE(dummy_out1) #home DLS

    def home_DLS(self):
        my_dummy = myDLS.OR(dummy_out1) #home DLS

    def move_stage(self,pos):
        my_dummy = myDLS.PA_Set(pos, dummy_out1)[1] # move the stage
        print(my_dummy)
        
    def change_velocity(self,vel):
        my_dummy = myDLS.VA_Set(vel, dummy_out1)[1] # change the velocity
        
    def change_acceleration(self,acc):
        my_dummy = myDLS.AC_Set(acc, dummy_out1)[1] # change the acceleration

    def enable_DLS(self):
        my_dummy = myDLS.MM_Set(1,dummy_out1) #enable DLS

    def disable_DLS(self):
        my_dummy = myDLS.MM_Set(0,dummy_out1) #disable DLS

    def stop_DLS(self):
        my_dummy = myDLS.ST(dummy_out1) #stop DLS

    def close(self):
        myDLS.CloseInstrument()

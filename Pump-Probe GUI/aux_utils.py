# -*- coding: utf-8 -*-
"""
@author: VINIT
Auxilliary Utilities
"""
import os
from tkinter import filedialog


def createfolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print('Directory Created')
        else:
            print('Directory already exists')
    except OSError:
        print('Error: Creating directory. ' + directory)

def getFileName(myDict):
    filename = ''
    return filename

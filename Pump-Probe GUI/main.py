# -*- coding: utf-8 -*-
"""
@author: VINIT
main console
To DO:
    # TODO Initialization Flags
    # TODO Add unit Markers
"""
import pickle
import time
from datetime import date
from os.path import exists

import aux_utils
import gui

config_filename = "config.pkl"
if exists(config_filename):
    with open(config_filename, 'rb') as f:

        general_param_dict, exp_param_dict = pickle.load(f)
        exp_param_dict['region_non-volatile']['Init_Time'] = time.perf_counter()
        exp_param_dict['region_non-volatile']['today'] = date.today().strftime("%Y%m%d")

else:
    general_param_dict = {'filepath': "../data",
                          'file_parameters':
                              {
                                  'Repetion Rate': "1000kHz",
                                  'Material': "crosscorr_SHG",
                                  'Temperature': "RT",
                                  'Amplifier Gain/ Mode': "none",
                                  'Lock-In Amplifier': "lockinAC_100u_1s_541Hz",
                                  'Test Suffix': "8.0mmx_test"
                              }
                          }

    aux_utils.createfolder(general_param_dict['filepath'])

    exp_param_dict = {'region_non-volatile': {
        'Init_Time': time.perf_counter(),
        'c': 299792458,
        'today': date.today().strftime("%Y%m%d")
    },
        'region_volatile': {
            'DLS Port': 'COM9',
            'Initial Position': 96.3,
            'Resolution': 10,
            'Time Constant': 0.2,
            'Scanning Period / s': 2,
            'Number of Averaging Scans': 1,
            'Background Reading Light': 0,
            'Background Sampling Time': 5,
            'Pump Wavelength': 920,
            'Pump Power': 0.5,
            'Probe Wavelength': 800,
            'Probe Power': 0.05,
            'Point Sampling Time': 0.2 / 2,
            'Positions to Scan': '2,9,5',
            'Motor Wait Time': 2
        },
        'DLS': {
            'Velocity': 10,
            'Acceleration': 10,
            'Desired Position': 150
        },
        'Lock In': {
            'Device ID': "dev4402",
            'apilevel': 6,
            'Server Host': "192.168.68.202"
        }

    }

    my_list = [general_param_dict, exp_param_dict]

    with open(config_filename, 'wb') as f:
        pickle.dump(my_list, f)

gui.Gui(general_param_dict, exp_param_dict)

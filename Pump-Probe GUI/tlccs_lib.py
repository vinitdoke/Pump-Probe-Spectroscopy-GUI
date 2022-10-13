import os
import numpy as np
import pyvisa
import ctypes as ct

SPECTROMETER_MODELS = {'0x8081': 'CCS100', '0x8083': 'CCS125', '0x8085': 'CCS150', '0x8087': 'CCS175',
                       '0x8089': 'CCS200'}
DEVICE_STATUS = {'0x0002': 'IDLE', '0x0004': 'SCANNING', '0x0008': 'SCAN START', '0x0010': 'SCAN COMPLETE, DATA READY',
                 '0x0080': 'WAIT FOR EXT TRIGGER'}
NUM_PIXELS = 3648


class tlccs:

    def __init__(self, model=None, serial=None):
        # load C library containing functions for interacting with Thorlabs CCS Spectrometers
        self.lib = ct.cdll.LoadLibrary(os.getcwd() + r"\instr_libs\TL_OSA_DLLs\TLCCS_64.dll")

        # Scan for connected spectrometers
        self.visa_rsrc_str = self.find_device(model, serial)
        self.vendor_id, self.product_id, self.serial = self.parse_rsrc_string(self.visa_rsrc_str)
        print('Connected to Thorlabs %s Spectrometer (%s)' % (SPECTROMETER_MODELS[self.product_id], self.serial))

        # Connect to spectrometer (ie open session)
        self.instr_handle = ct.c_int()
        self.openSession(self.visa_rsrc_str, self.instr_handle)

        self.background = np.zeros(NUM_PIXELS)

    def openSession(self, rsrc_str, instr_handle, idQuery=1, reset=1):
        self.lib.tlccs_init(rsrc_str.encode(), idQuery, reset, ct.byref(instr_handle))

    def closeSession(self):
        self.lib.tlccs_close(self.instr_handle)

    def reset(self):
        self.lib.tlccs_reset(self.instr_handle)

    def setIntegrationTime(self, integration_time=0.01):
        c_int_time = ct.c_double(integration_time)
        self.lib.tlccs_setIntegrationTime(self.instr_handle, c_int_time)
        return self.getIntegrationTime()

    def getIntegrationTime(self):
        integration_time = ct.c_double()
        self.lib.tlccs_getIntegrationTime(self.instr_handle, ct.byref(integration_time))
        return integration_time.value

    def startScan(self):
        self.lib.tlccs_startScan(self.instr_handle)

    def startScanCont(self):
        self.lib.tlccs_startScanCont(self.instr_handle)

    def getWavelengthData(self):
        wavelengths = (ct.c_double * NUM_PIXELS)()
        self.lib.tlccs_getWavelengthData(self.instr_handle, 0, ct.byref(wavelengths), ct.c_void_p(None),
                                         ct.c_void_p(None))
        return np.ctypeslib.as_array(wavelengths)

    def getScanData(self):
        data_array = (ct.c_double * NUM_PIXELS)()
        self.lib.tlccs_getScanData(self.instr_handle, ct.byref(data_array))
        return np.ctypeslib.as_array(data_array) - self.background

    def getRawScanData(self):
        # Returns data without doing any processing
        data_array = (ct.c_double * NUM_PIXELS)()
        self.lib.tlccs_getScanData(self.instr_handle, ct.byref(data_array))
        return np.ctypeslib.as_array(data_array)

    def setBackground(self):
        self.startScan()
        self.background = self.getRawScanData()

    def getDeviceStatus(self):
        status = ct.c_int32()
        self.lib.tlccs_getDeviceStatus(self.instr_handle, ct.byref(status))
        return status.value

    def find_device(self, model, serial):
        dev_list = self.list_devices(model)
        if not dev_list:
            raise Exception(
                'No devices found, please ensure spectrometer is connected and available for new VISA session (i.e. close Thorlabs OSA software).')
        else:
            if len(dev_list) > 1:
                if serial:
                    dev_match = [dev for dev in dev_list if serial in dev]
                    return dev_match
                else:
                    print(
                        'Multiple devices found. To choose different spectrometer, specify model or serial number when initiating tlccs object.')
            return dev_list[0]

    def list_devices(self, model):
        # if model:
        # search_str = "USB?*?{VI_ATTR_MANF_ID==0x1313 && (VI_ATTR_MODEL_CODE==%s)}"%(model)
        # else:
        search_str = "USB?*?{VI_ATTR_MANF_ID==0x1313 && ((VI_ATTR_MODEL_CODE==0x8081) || (VI_ATTR_MODEL_CODE==0x8083) || (VI_ATTR_MODEL_CODE==0x8085) || (VI_ATTR_MODEL_CODE==0x8087) || (VI_ATTR_MODEL_CODE==0x8089))}"
        rm = pyvisa.ResourceManager()
        dev_list = rm.list_resources(search_str)
        return dev_list

    def parse_rsrc_string(self, rsrc_str):
        _, vid, pid, serial, _ = rsrc_str.split('::')
        return vid, pid, serial

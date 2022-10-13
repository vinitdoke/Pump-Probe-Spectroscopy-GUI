import os.path
from tkinter import ttk
import pickle
import time
import tkinter as tk
from tkinter import filedialog, messagebox

import matplotlib.pyplot as plt
import numpy as np
import zhinst.utils as zu
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from scipy.optimize import curve_fit

import aux_utils
import newportDLS_cmd as Dls

general_param_dict = {}
exp_param_dict = {}
myDLS = None
daq = None
device = None
demod_path = None
myCanvas = None
generated_params = {
    'R': 1.0,
    'Background Data': [],
}

global filetab


class DLSTab:
    def __init__(self, notebook):
        global exp_param_dict
        self.port = None
        self.frame = ttk.Frame(notebook)
        self.label = ttk.Label(self.frame, text="DLS PORT")
        self.label.grid(row=0, column=0, pady=2, padx=2)
        self.portVar = tk.StringVar()
        self.portVar.set(exp_param_dict['region_volatile']['DLS Port'])
        self.portEntry = ttk.Entry(self.frame, textvariable=self.portVar)
        self.portEntry.grid(row=0, column=1, pady=2, padx=2)
        self.initbt = ttk.Button(self.frame, text="Initialize", command=self.initDLS)
        self.initbt.grid(row=1, column=1, pady=2, padx=2)
        self.terminatebt = ttk.Button(self.frame, text="Terminate", command=lambda: myDLS.close())
        self.terminatebt.grid(row=1, column=2, pady=2, padx=2)
        self.enablebt = ttk.Button(self.frame, text="Enable", command=lambda: myDLS.enable_DLS())
        self.enablebt.grid(row=2, column=1, pady=2, padx=2)
        self.disablebt = ttk.Button(self.frame, text="Disable", command=lambda: myDLS.disable_DLS())
        self.disablebt.grid(row=2, column=2, pady=2, padx=2)
        self.label2 = ttk.Label(self.frame, text="SET VALUES : ")
        self.label2.grid(row=3, column=0, pady=2, padx=2)
        self.entryfields = []
        self.entrylabels = []
        self.stringVars = []
        i = 4
        for param, value in exp_param_dict['DLS'].items():
            self.stringVars.append(tk.StringVar(notebook, value))
            self.entryfields.append(ttk.Entry(self.frame, textvariable=self.stringVars[-1]))
            self.entryfields[-1].grid(row=i + 1, column=1, pady=2, padx=2)
            self.entrylabels.append(ttk.Label(self.frame, text=param))
            self.entrylabels[-1].grid(row=i + 1, column=0, padx=2, pady=2)
            i += 1
        self.lockButton = ttk.Button(self.frame, text="Set Both", command=self.lock_params)
        self.lockButton.grid(row=i - 1, column=2, padx=2, pady=2)

        self.moveButton = ttk.Button(self.frame, text="Move",
                                     command=lambda: myDLS.move_stage(float(self.stringVars[-1].get())))
        self.moveButton.grid(row=i, column=2, padx=2, pady=2)

        self.stopButton = ttk.Button(self.frame, text="Stop", command=lambda: myDLS.stop_DLS())
        self.stopButton.grid(row=i, column=3, padx=2, pady=2)

        self.label3 = ttk.Label(self.frame, text="CURRENT TELEMETRY")
        self.label3.grid(row=i + 2, column=0, padx=2, pady=2)
        self.vellb = ttk.Label(self.frame, text="Current velocity :")
        self.acclb = ttk.Label(self.frame, text="Current acceleration :")
        self.poslb = ttk.Label(self.frame, text="Current position :")
        self.vellb.grid(row=i + 3, column=0, padx=2, pady=2)
        self.acclb.grid(row=i + 4, column=0, padx=2, pady=2)
        self.poslb.grid(row=i + 5, column=0, padx=2, pady=2)

        self.vel = tk.StringVar(notebook)
        self.acc = tk.StringVar(notebook)
        self.pos = tk.StringVar(notebook)

        self.vellb2 = ttk.Label(self.frame, textvariable=self.vel)
        self.acclb2 = ttk.Label(self.frame, textvariable=self.acc)
        self.poslb2 = ttk.Label(self.frame, textvariable=self.pos)
        self.vellb2.grid(row=i + 3, column=1, padx=2, pady=2)
        self.acclb2.grid(row=i + 4, column=1, padx=2, pady=2)
        self.poslb2.grid(row=i + 5, column=1, padx=2, pady=2)

        self.frame.pack(fill="both", expand=1)
        self.initDLS()

    def initDLS(self):
        global exp_param_dict, myDLS

        self.port = self.portVar.get()
        myDLS = Dls.DLS_cmd(self.port)
        myDLS.initialize_DLS()
        self.update()

    def lock_params(self):
        global exp_param_dict, myDLS
        i = 0
        for param in exp_param_dict['DLS'].keys():
            exp_param_dict['DLS'][param] = self.stringVars[i].get()
            i += 1
        myDLS.change_velocity(float(exp_param_dict['DLS']['Velocity']))
        myDLS.change_acceleration(float(exp_param_dict['DLS']['Acceleration']))
        print(exp_param_dict['DLS'])

    def update(self):
        self.vel.set(myDLS.current_vel())
        self.acc.set(myDLS.current_acc())
        self.pos.set(myDLS.current_pos())
        self.frame.after(300, self.update)
        # print("updated")


class FileTab:
    def __init__(self, notebook):
        global general_param_dict
        self.frame = ttk.Frame(notebook)
        # self.label = ttk.Label(self.frame, text="Select Root Directory :")
        # self.label.grid(row=0, column=0, pady=2, padx=2)
        # self.browsebt = ttk.Button(self.frame, text="Browse", command=self.pathchange)
        # self.browsebt.grid(row=0, column=1, pady=2, padx=2)
        self.DirectoryLabel = ttk.Label(self.frame, text="Current Directory")
        self.DirectoryLabel.grid(row=0, column=0, pady=2, padx=2)

        self.directory = tk.StringVar(self.frame, value=general_param_dict['filepath'])
        self.DirectoryField = ttk.Entry(self.frame, textvariable=self.directory)
        self.DirectoryField.grid(row=0, column=1, pady=2, padx=2)
        self.SetFilePathbt = ttk.Button(self.frame, text='Set Root Directory', command=self.setPath)
        self.SetFilePathbt.grid(row=0, column=2, pady=2, padx=2)

        self.entryfields = []
        self.entrylabels = []
        self.stringVars = []
        self.IntVars = []
        self.checkBoxes = []
        self.label2 = ttk.Label(self.frame, text="File Prefix Details:")
        self.label2.grid(row=1, column=0, pady=2, padx=2)
        self.label3 = ttk.Label(self.frame, text="Override Name")
        self.label3.grid(row=2, column=0, pady=2, padx=2)
        self.overridetext = tk.StringVar()
        self.entryOverride = ttk.Entry(self.frame, textvariable=self.overridetext)
        self.entryOverride.grid(row=2, column=1, pady=2, padx=2)
        i = 2
        for param, value in general_param_dict['file_parameters'].items():
            self.stringVars.append(tk.StringVar(notebook, value))
            self.entryfields.append(ttk.Entry(self.frame, textvariable=self.stringVars[-1]))
            self.entryfields[-1].grid(row=i + 1, column=2, pady=2, padx=2)
            self.entrylabels.append(ttk.Label(self.frame, text=param))
            self.entrylabels[-1].grid(row=i + 1, column=1, padx=2, pady=2)
            self.IntVars.append(tk.IntVar())
            self.checkBoxes.append(ttk.Checkbutton(self.frame, text='', variable=self.IntVars[-1]))
            self.checkBoxes[-1].grid(row=i + 1, column=0, padx=2, pady=2)
            i += 1

        self.lockButton = ttk.Button(self.frame, text="Save Parameters", command=self.lock_params)
        self.lockButton.grid(row=i + 1, column=2, padx=2, pady=2)

        self.current_filename = tk.StringVar(self.frame)
        self.CurrentFilenameLabel = ttk.Label(self.frame, text="File Prefix :")
        self.CurrentFilenameLabel.grid(row=i + 2, column=0, padx=2, pady=2)
        self.filenamedisplay = ttk.Label(self.frame, textvariable=self.current_filename)
        self.filenamedisplay.grid(row=i + 2, column=1, padx=2, pady=2)

        self.frame.pack(fill="both", expand=1)

    def setPath(self):
        new_directory = self.directory.get()
        general_param_dict['filepath'] = new_directory
        aux_utils.createfolder(new_directory)
        print('Filepath is ', general_param_dict['filepath'])

    def pathchange(self):
        global root
        # print('here')
        # root.withdraw()
        filepath = filedialog.askopenfilename()
        filepath += '/data'
        # print('here2')
        # aux_utils.createfolder(filepath)
        # print('here3')
        global general_param_dict
        print(filepath)
        general_param_dict['filepath'] = filepath

    def lock_params(self):
        global general_param_dict
        i = 0
        for param in general_param_dict['file_parameters'].keys():
            general_param_dict['file_parameters'][param] = self.stringVars[i].get()
            i += 1

        if len(self.overridetext.get()) == 0:
            filename = []
            for i in range(len(self.stringVars)):
                if self.IntVars[i].get():
                    filename.append(self.stringVars[i].get())
            filename = '_'.join(filename)
        else:
            filename = self.overridetext.get()

        self.current_filename.set(filename)


class Scanning:
    def __init__(self, notebook):
        self.breakFlag = False
        self.points = None
        self.b = None
        self.a = None
        self.frame = ttk.Frame(notebook)
        self.label1 = ttk.Label(self.frame, text="SCANNING COMMANDS")
        self.label1.grid(row=0, column=0, pady=2, padx=2)
        self.label2 = ttk.Label(self.frame, text="Background Scan :")
        self.label2.grid(row=1, column=0, pady=2, padx=2)

        self.bgscanbt = ttk.Button(self.frame, text="Scan Background", command=self.bg_scan)
        self.bgscanbt.grid(row=1, column=1, pady=2, padx=2)

        self.Rvalue = tk.StringVar(notebook)
        self.Rvalue.set(generated_params['R'])
        self.Rfield = ttk.Entry(self.frame, textvariable=self.Rvalue)
        self.SetRButton = ttk.Button(self.frame, text="Set R Value", command=self.setR)
        self.SetRButton.grid(row=2, column=2, pady=2, padx=2)
        self.Rfield.grid(row=2, column=1, padx=2, pady=2)
        self.Rlabel = ttk.Label(self.frame, text="R (background)")
        self.Rlabel.grid(row=2, column=0, padx=2, pady=2)

        self.label3 = ttk.Label(self.frame, text="MAIN SCAN")
        self.label3.grid(row=3, column=0, pady=2, padx=2)

        self.label4 = ttk.Label(self.frame, text="Choose Scanning Function :")
        self.label4.grid(row=4, column=0, pady=2, padx=2)

        self.scan_flag = tk.IntVar()

        self.r1 = tk.Radiobutton(self.frame, text='Point Scan', variable=self.scan_flag, value=0)
        self.r1.grid(row=5, column=0, pady=2, padx=2)
        self.pointlabel = ttk.Label(self.frame, text='Points To Scan')
        self.pointlabel.grid(row=6, column=0, pady=2, padx=2)
        self.PointsToScan = tk.StringVar(self.frame, value=exp_param_dict['region_volatile']['Positions to Scan'])
        self.EntryPointsToScan = ttk.Entry(self.frame, textvariable=self.PointsToScan)
        self.EntryPointsToScan.grid(row=6, column=1, pady=2, padx=2)

        self.r2 = tk.Radiobutton(self.frame, text='Segmented Scan', variable=self.scan_flag, value=1)
        self.r2.grid(row=7, column=0, pady=2, padx=2)

        j = 8

        self.finescanlen = tk.StringVar()
        self.scanlenafter = tk.StringVar()
        self.resolution = tk.StringVar()
        self.vfactor = tk.StringVar()
        self.initpos = tk.StringVar()
        self.scanDir = tk.IntVar()

        self.finescanlen.set('10')
        self.scanlenafter.set('20')
        self.vfactor.set('5')
        self.resolution.set(str(exp_param_dict['region_volatile']['Resolution']))
        self.initpos.set(str(exp_param_dict['region_volatile']['Initial Position']))
        self.scanDir.set(1)

        self.initposlabel = ttk.Label(self.frame, text="Initial Position (mm)")
        self.initposlabel.grid(row=j, column=0, pady=2, padx=2)
        self.initposEntry = ttk.Entry(self.frame, textvariable=self.initpos)
        self.initposEntry.grid(row=j, column=1, pady=2, padx=2)

        j += 2
        self.reslabel = ttk.Label(self.frame, text="Resolution (fs)")
        self.reslabel.grid(row=j - 1, column=0, padx=2, pady=2)
        self.resEntry = ttk.Entry(self.frame, textvariable=self.resolution)
        self.resEntry.grid(row=j - 1, column=1, padx=2, pady=2)
        self.finelabel = ttk.Label(self.frame, text="Fine Scan Length (ps)")
        self.finelabel.grid(row=j, column=0, padx=2, pady=2)
        self.afterlabel = ttk.Label(self.frame, text="Coarse Scan Length (ps)")
        self.afterlabel.grid(row=j + 1, column=0, padx=2, pady=2)
        self.FineEntry = ttk.Entry(self.frame, textvariable=self.finescanlen)
        self.FineEntry.grid(row=j, column=1, padx=2, pady=2)
        self.CoarseEntry = ttk.Entry(self.frame, textvariable=self.scanlenafter)
        self.CoarseEntry.grid(row=j + 1, column=1, padx=2, pady=2)

        self.vfactorlabel = ttk.Label(self.frame, text="Velocity Factor")
        self.vfactorlabel.grid(row=j + 2, column=0, padx=2, pady=2)
        self.vFactorEntry = ttk.Entry(self.frame, textvariable=self.vfactor)
        self.vFactorEntry.grid(row=j + 2, column=1, padx=2, pady=2)

        self.scanDirlabel = ttk.Label(self.frame, text="Scan Direction")
        self.scanDirlabel.grid(row=j + 3, column=0, padx=2, pady=2)
        self.r3 = tk.Radiobutton(self.frame, text='Forward', variable=self.scanDir, value=1)
        self.r3.grid(row=j + 3, column=1, pady=2, padx=2)
        self.r4 = tk.Radiobutton(self.frame, text='Backward', variable=self.scanDir, value=0)
        self.r4.grid(row=j + 3, column=2, pady=2, padx=2)

        j += 4

        self.scanbt = ttk.Button(self.frame, text="Start Scanning", command=self.scan)
        self.scanbt.grid(row=j, column=1, pady=2, padx=2)
        self.stopbt = ttk.Button(self.frame, text="Stop", command=self.stopscan)
        self.stopbt.grid(row=j, column=2, pady=2, padx=2)
        self.frame.pack(fill="both", expand=1)

    def stopscan(self):
        myDLS.stop_DLS()
        time.sleep(0.1)
        myDLS.disable_DLS()
        time.sleep(0.1)
        self.breakFlag = True

    def setR(self):
        global generated_params
        generated_params['R'] = float(self.Rvalue.get())
        print(generated_params)

    def bg_scan(self):
        global generated_params
        t_init = time.perf_counter()
        while time.perf_counter() - t_init < float(exp_param_dict['region_volatile']['Background Sampling Time']):
            sample = daq.getSample(demod_path)
            generated_params['Background Data'].append(sample['auxin1'][0])
            time.sleep(0.1)

        generated_params['R'] = abs(np.mean(generated_params['Background Data']) - int(
            exp_param_dict['region_volatile']['Background Reading Light'] / 1000))
        self.Rvalue.set(generated_params['R'])
        messagebox.showinfo("R is", generated_params['R'])
        MsgBox = tk.messagebox.askquestion('Save Scan',
                                           'Calibration Complete! Do you want to save the R file (under parent directory)',
                                           icon='info')
        if MsgBox == 'yes':
            np.savetxt(general_param_dict['filepath'] + '/R.txt', generated_params['Background Data'], fmt='%.9e')
        else:
            tk.messagebox.showinfo('Return', 'You will now return to the application screen')

    def scan(self):
        global filetab
        # Overwrite Check
        filename = filetab.current_filename.get()
        name = general_param_dict['filepath'] + '/' + filename + '0' + '.txt'
        if os.path.exists(name):
            MsgBox = tk.messagebox.askquestion('File already exists!',
                                               'Do you want to overwrite ?',
                                               icon='info')
            if MsgBox == 'yes':
                pass
            else:
                return

        if self.scan_flag.get() == 1:
            self.breakFlag = False
            self.segmentedscan2()
        elif self.scan_flag.get() == 0:
            self.breakFlag = False
            self.pointscan()
        else:
            messagebox.showinfo('Issue', 'Select scanning function')

    def waitUntilReached(self, destination):
        while abs(myDLS.current_pos() - destination) > 0.001:
            pass

    def pointscan(self):
        global exp_param_dict, generated_params, myCanvas, general_param_dict, filetab
        exp_param_dict['region_volatile']['Positions to Scan'] = self.PointsToScan.get()
        self.points = [float(temp) for temp in
                       exp_param_dict['region_volatile']['Positions to Scan'].strip().split(',')]
        # print(self.points)
        print("Starting Pointscan....")
        for point in self.points:
            position = float(point) * float(exp_param_dict['region_non-volatile']['c']) * 1e-12 * 1e3 / 4 + \
                       float(exp_param_dict['region_volatile'][
                                 'Initial Position'])  # convert delay time to actual delay line positions
            myDLS.move_stage(position)
            time.sleep(int(exp_param_dict['region_volatile']['Motor Wait Time']))
            # self.waitUntilReached(position)
            print("Now scanning at", position)
            x_temp = []
            temp_data = []
            datax = []

            t_start = time.perf_counter()
            while time.perf_counter() - t_start < float(exp_param_dict['region_volatile']['Point Sampling Time']):
                temp_pos = np.float64(myDLS.current_pos())
                delay = np.float64((temp_pos - float(exp_param_dict['region_volatile']['Initial Position']) * 1e-3 /
                                    float(exp_param_dict['region_non-volatile']['c']) * 1e12 * 4))
                sample = daq.getSample(demod_path)
                data_x = np.float64(sample['x'][0])
                data_y = np.float64(sample['y'][0])
                data_phase = np.float64(np.arctan(data_y / data_x))
                r = np.float64(np.sqrt(sample['x'] ** 2 + sample['y'] ** 2))
                cur_time = np.float64(time.perf_counter() - float(exp_param_dict['region_non-volatile']['Init_Time']))
                temp_data.append([temp_pos, delay, data_x, data_y, r, data_phase, cur_time])
                x_temp.append(data_x)

                myCanvas.dataplot.append((temp_pos, delay, np.mean(x_temp) / generated_params['R']))
                datax = np.array(myCanvas.dataplot)

                myCanvas.line1.set_data(datax[:, 1], datax[:, 2] / generated_params['R'])

                plt.pause(0.01)
                time.sleep(0.01)
                # TODO Change Plotting for pointscan

            myCanvas.data = temp_data

            # save current scan
            filename = filetab.current_filename.get()

            file = open(general_param_dict['filepath'] + '/' + filename + str(point) + 'ps.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + str(point) + 'ps.txt', np.array(myCanvas.data),
                       fmt='%.9e')
            file.close()
            # save avg file
            file = open(general_param_dict['filepath'] + '/' + filename + '_avg.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + '_avg.txt', datax, fmt='%.9e')
            file.close()
            myCanvas.fig.canvas.flush_events()

    def segmentedscan2(self):
        global exp_param_dict, generated_params, general_param_dict, myCanvas
        res = float(self.resolution.get())  # fs
        # time_constant = 0.2  # TC of lock-in, second
        # t_wait = 10 * float(exp_param_dict['region_volatile']['Point Sampling Time'])
        sampling_time = float(exp_param_dict['region_volatile']['Point Sampling Time'])  # second
        c = float(exp_param_dict['region_non-volatile']['c'])

        vel = res * 1e-15 * c / 2 / sampling_time * 1e3
        print(f'resolution is {res}')
        acc = 5 * vel
        scan_fine_t = float(self.finescanlen.get())  # ps
        scan_len_t_after = float(self.scanlenafter.get())
        v_factor = int(self.vfactor.get())

        # scan_len_t = scan_fine_t + scan_len_t_after

        if self.scanDir.get():
            print(' Moving Forward')
            pos_0 = float(self.initpos.get())
            pos_1 = pos_0 + scan_fine_t * c * 1e-12 / 4 * 1000
            pos_2 = pos_1 + scan_len_t_after * c * 1e-12 / 4 * 1000
        else:
            print(' Moving Backward')
            pos_0 = float(self.initpos.get())
            pos_1 = pos_0 - scan_fine_t * c * 1e-12 / 4 * 1000
            pos_2 = pos_1 - scan_len_t_after * c * 1e-12 / 4 * 1000

        pos_0 = max(pos_0, 0)
        pos_2 = min(pos_2, 300)

        data = []

        for i in range(0, int(exp_param_dict['region_volatile']['Number of Averaging Scans'])):
            data.append([])
            print("scan number:", i + 1)
            j = 0
            if i == 1:
                npdata = np.array(data[0])
                j_length = len(npdata[:, 2])  # the array length all the datasets need to match

            myDLS.stop_DLS()  # in case the stage is still moving
            time.sleep(1)
            myDLS.enable_DLS()
            myDLS.change_velocity(25)
            print('changed velocity to 25')
            myDLS.change_acceleration(100)
            print(f'moving to {pos_0}')
            myDLS.move_stage(pos_0)
            print(f'moving to {pos_0}')
            # time.sleep(5)
            self.waitUntilReached(pos_0)
            print('reached pos_0')
            dataplot = []
            temp_data = []

            myDLS.update()
            time.sleep(1)
            print(f'changing velocity to {vel}')
            print(f'velocity : {vel} ; acceleration : {acc}')
            myDLS.change_velocity(vel * 1.0)
            myDLS.change_acceleration(acc * 1.0)
            time.sleep(1)
            myDLS.update()
            print(f'moving to pos_1 : {pos_1}')
            myDLS.move_stage(pos_1)

            while abs(myDLS.current_pos() - pos_1) > 0.0001:
                # print(myDLS.current_pos())
                if self.breakFlag:
                    # myDLS.stop()
                    break
                temp_pos = np.float(myDLS.current_pos())
                # print(temp_pos, f'destination {pos_1}')
                delay = np.float((temp_pos - pos_0) * 1e-3 / c * 1e12 * 4)

                # t_start = time.perf_counter()

                sample = daq.getSample(demod_path)
                data_x = np.float(sample['x'][0])
                data_y = np.float(sample['y'][0])
                data_phase = np.float(np.arctan(data_y / data_x))
                r = np.float(np.sqrt(sample['x'] ** 2 + sample['y'] ** 2))
                cur_time = np.float(time.perf_counter() - float(exp_param_dict['region_non-volatile']['Init_Time']))
                # t_start = time.perf_counter()

                temp_data.append([temp_pos, delay, data_x, data_y, r, data_phase, cur_time])
                # time.sleep(sampling_time / 100)

                # if i > 0:
                #     x_temp[j] = x_temp[j] + data_x
                # else:
                #     x_temp.append(data_x / float(generated_params['R']))

                dataplot.append((temp_pos, delay, data_x / (i + 1)))
                datax = np.array(dataplot)

                myCanvas.line1.set_data(datax[:, 1], datax[:, 2] / float(generated_params['R']))

                # if len(datax[:, 2]) > 2 and (
                #         datax[-1, 2] / float(
                #         generated_params['R']) <= myCanvas.line1.axes.get_ylim()[0] or datax[-1, 2] / float(
                #         generated_params['R']) >=
                #         myCanvas.line1.axes.get_ylim()[1]):
                #     plt.ylim([np.amin(datax[:, 2] / float(generated_params['R'])),
                #               np.amax(datax[:, 2]) / float(generated_params['R'])])
                # plt.pause(sampling_time / 20)
                # myCanvas.dataplot.append((0, x, y))

                # datax = np.array(myCanvas.dataplot)
                # print(datax)

                # myCanvas.line1.set_data(datax[:, 1], datax[:, 2])
                if np.amin(datax[:, 2]) != np.amax(datax[:, 2]):
                    myCanvas.ax.set_ylim(np.amin(datax[:, 2] / float(generated_params['R'])),
                                         np.amax(datax[:, 2] / float(generated_params['R'])))
                if np.amin(datax[:, 1]) != np.amax(datax[:, 1]):
                    myCanvas.ax.set_xlim(np.amin(datax[:, 1]), np.amax(datax[:, 1]))
                # plt.pause(0.01)
                # time.sleep(0.01)

                myCanvas.ax2.set_xlim(myCanvas.ax.get_xlim())
                x2_min = np.min(datax[:, 1])
                x2_max = np.max(datax[:, 1])
                pos_min = datax[:, 0][0]
                pos_max = datax[:, 0][-1]
                ticks = np.linspace(x2_min, x2_max, 6)
                if self.scanDir.get():
                    ticklabels = np.linspace(pos_min, pos_max, 6)
                else:
                    ticklabels = np.linspace(pos_min, pos_max, 6)[::-1]
                myCanvas.ax2.set_xticks(ticks)
                myCanvas.ax2.set_xticklabels(np.around(ticklabels, 3))

                myCanvas.fig.canvas.flush_events()
                myCanvas.canva.draw()
                j += 1

            time.sleep(1)
            myDLS.update()
            print(f'vfactor is {v_factor}')
            myDLS.change_velocity(vel * v_factor)
            myDLS.change_acceleration(acc * v_factor)
            # print(f' velocity is {vel*v_factor}')
            myDLS.update()
            time.sleep(1)
            myDLS.update()
            print('starting coarse scan')
            myDLS.move_stage(pos_2)
            while abs(myDLS.current_pos() - pos_2) > 0.0001:
                if self.breakFlag:
                    break
                temp_pos = np.float(myDLS.current_pos())
                # print(temp_pos, f'destination {pos_2}')
                delay = np.float((temp_pos - pos_0) * 1e-3 / c * 1e12 * 4)
                # t_start = time.perf_counter()
                sample = daq.getSample(demod_path)
                data_x = np.float(sample['x'][0])
                data_y = np.float(sample['y'][0])
                data_phase = np.float(np.arctan(data_y / data_x))
                r = np.float(np.sqrt(sample['x'] ** 2 + sample['y'] ** 2))
                cur_time = np.float(time.perf_counter() - float(exp_param_dict['region_non-volatile']['Init_Time']))
                # t_start = time.perf_counter()

                temp_data.append([temp_pos, delay, data_x, data_y, r, data_phase, cur_time])
                # time.sleep(sampling_time / 10)

                # if i > 0:
                #     x_temp[j] = x_temp[j] + data_x
                # else:
                #     x_temp.append(data_x / float(generated_params['R']))

                dataplot.append((temp_pos, delay, data_x / (i + 1)))
                datax = np.array(dataplot)

                myCanvas.line1.set_data(datax[:, 1], datax[:, 2] / float(generated_params['R']))

                # if len(datax[:, 2]) > 2 and (
                #         datax[-1, 2] / float(generated_params['R']) <= myCanvas.line1.axes.get_ylim()[0] or datax[
                #         -1, 2] / float(generated_params['R']) >= myCanvas.line1.axes.get_ylim()[1]):
                #     # plt.ylim([np.amin(datax[:, 2] / float(generated_params['R'])),
                #     #           np.amax(datax[:, 2]) / float(generated_params['R'])])
                #     myCanvas.ax.set_ylim(np.amin(datax[:, 2] / float(generated_params['R'])),
                #                      np.amax(datax[:, 2] / float(generated_params['R'])))
                # plt.pause(sampling_time / 20)
                # myCanvas.dataplot.append((0, x, y))

                # datax = np.array(myCanvas.dataplot)
                # print(datax)

                # myCanvas.line1.set_data(datax[:, 1], datax[:, 2])
                myCanvas.ax.set_ylim(np.amin(datax[:, 2] / float(generated_params['R'])),
                                     np.amax(datax[:, 2] / float(generated_params['R'])))
                myCanvas.ax.set_xlim(np.amin(datax[:, 1]), np.amax(datax[:, 1]))

                myCanvas.ax2.set_xlim(myCanvas.ax.get_xlim())
                x2_min = np.min(datax[:, 1])
                x2_max = np.max(datax[:, 1])
                # pos_min = np.min(datax[:, 0])
                # pos_max = np.max(datax[:, 0])
                pos_min = datax[:, 0][0]
                pos_max = datax[:, 0][-1]
                ticks = np.linspace(x2_min, x2_max, 6)
                if self.scanDir.get():
                    ticklabels = np.linspace(pos_min, pos_max, 6)
                else:
                    ticklabels = np.linspace(pos_min, pos_max, 6)[::-1]
                myCanvas.ax2.set_xticks(ticks)
                myCanvas.ax2.set_xticklabels(np.around(ticklabels, 3))

                myCanvas.fig.canvas.flush_events()
                myCanvas.canva.draw()
                j += 1
                if i > 0:
                    if j == j_length:
                        break

            if i > 0:
                if j < j_length:
                    for k in range(j, j_length):
                        temp_data.append([0, 0, 0, 0, 0, 0, 0])

            data[i] = temp_data

            # save current scan
            # if len(filetab.overridetext.get()) == 0:
            #     filename = []
            #     for i in range(len(filetab.stringVars)):
            #         if filetab.IntVars[i].get():
            #             filename.append(filetab.stringVars[i].get())
            #     filename = '_'.join(filename)
            # else:
            #     filename = filetab.overridetext.get()
            filename = filetab.current_filename.get()

            file = open(general_param_dict['filepath'] + '/' + filename + str(i) + '.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + str(i) + '.txt', np.array(data[i]), fmt='%.9e')
            file.close()

            # save avg file
            file = open(general_param_dict['filepath'] + '/' + filename + '_avg.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + '_avg.txt', datax, fmt='%.9e')
            file.close()
            myCanvas.fig.canvas.flush_events()

        myDLS.update()
        # plt.savefig(general_param_dict['filepath'] + '/' + 'plot' + '.png')
        print("Moving back to initial position.")
        myDLS.change_velocity(25)
        myDLS.change_acceleration(100)
        myDLS.move_stage(pos_0)
        myDLS.disable_DLS()
        # myDLS.di
        self.breakFlag = True
        time.sleep(2)


class Correlation:
    def __init__(self, notebook):
        self.frame = ttk.Frame(notebook)
        self.label = ttk.Label(self.frame, text="Correlation Measurements")
        self.label.grid(row=0, column=0, pady=2, padx=2)

        self.typeFlag = tk.IntVar()
        self.scanData = None
        self.resolution = tk.StringVar()

        self.startPoint = tk.StringVar()
        self.scanlen = tk.StringVar()
        self.FWHM = tk.StringVar()
        self.scanDir = tk.IntVar()

        self.typeFlag.set(1)
        self.startPoint.set('92.6')
        self.scanlen.set('20')
        self.resolution.set('10')
        self.FWHM.set("Not Measured")
        self.scanDir.set(1)

        self.l1 = ttk.Label(self.frame, text="Start Point")
        self.l1.grid(row=1, column=0, pady=2, padx=2)
        self.l3 = ttk.Label(self.frame, text="Resolution (fs)")
        self.l3.grid(row=3, column=0, pady=2, padx=2)
        self.l2 = ttk.Label(self.frame, text="Scan Length (ps)")
        self.l2.grid(row=2, column=0, pady=2, padx=2)

        self.e1 = ttk.Entry(self.frame, textvariable=self.startPoint)
        self.e1.grid(row=1, column=1, pady=2, padx=2)
        self.e2 = ttk.Entry(self.frame, textvariable=self.scanlen)
        self.e2.grid(row=2, column=1, pady=2, padx=2)
        self.e3 = ttk.Entry(self.frame, textvariable=self.resolution)
        self.e3.grid(row=3, column=1, pady=2, padx=2)

        j = 4
        self.scanDirlabel = ttk.Label(self.frame, text="Scan Direction")
        self.scanDirlabel.grid(row=j, column=0, padx=2, pady=2)
        self.r3 = tk.Radiobutton(self.frame, text='Forward', variable=self.scanDir, value=1)
        self.r3.grid(row=j, column=1, pady=2, padx=2)
        self.r4 = tk.Radiobutton(self.frame, text='Backward', variable=self.scanDir, value=0)
        self.r4.grid(row=j, column=2, pady=2, padx=2)

        j += 1
        self.r1 = tk.Radiobutton(self.frame, text='Auto Correlation', variable=self.typeFlag, value=1)
        self.r1.grid(row=j, column=0, pady=2, padx=2)

        j += 1
        self.r2 = tk.Radiobutton(self.frame, text='Cross Correlation', variable=self.typeFlag, value=0)
        self.r2.grid(row=j, column=0, pady=2, padx=2)

        j += 1
        self.startbt = ttk.Button(self.frame, text="Start Scan", command=self.scan)
        self.startbt.grid(row=j, column=1, pady=2, padx=2)

        j += 1
        self.PulseWidth = tk.StringVar()
        self.PulseWidth.set("Not Measured")

        self.label2 = ttk.Label(self.frame, text="PULSE WIDTH : ")
        self.label2.grid(row=j, column=0, padx=2, pady=2)

        self.label3 = ttk.Label(self.frame, textvariable=self.PulseWidth)
        self.label3.grid(row=j, column=1, padx=2, pady=2)

        self.label4 = ttk.Label(self.frame, text="FWHM")
        self.label4.grid(row=j + 1, column=0, padx=2, pady=2)

        self.label5 = ttk.Label(self.frame, textvariable=self.FWHM)
        self.label5.grid(row=j + 1, column=1, padx=2, pady=2)

        self.frame.pack(fill="both", expand=1)

    def waitUntilReached(self, destination):
        # t = time.perf_counter()
        while abs(myDLS.current_pos() - destination) > 0.001:
            pass

    def scan(self):

        global filetab
        # Overwrite Check
        filename = filetab.current_filename.get()
        name = general_param_dict['filepath'] + '/' + filename + '0' + '.txt'
        if os.path.exists(name):
            MsgBox = tk.messagebox.askquestion('File already exists!',
                                               'Do you want to overwrite ?',
                                               icon='info')
            if MsgBox == 'yes':
                pass
            else:
                return

        if self.typeFlag.get():
            print("Starting Auto-Correlation Scan")
            self.autoCorr()
        else:
            print("Starting Cross-Correlation Scan")
            self.crossCorr()

    def autoCorr(self):
        global exp_param_dict, generated_params, general_param_dict, myCanvas
        res = float(self.resolution.get())  # fs
        # time_constant = 0.2  # TC of lock-in, second
        # t_wait = 10 * float(exp_param_dict['region_volatile']['Point Sampling Time'])
        sampling_time = float(exp_param_dict['region_volatile']['Point Sampling Time'])  # second
        c = float(exp_param_dict['region_non-volatile']['c'])

        vel = res * 1e-15 * c / 2 / sampling_time * 1e3
        print(f'resolution is {res}')
        acc = 5 * vel
        scan_t = float(self.scanlen.get())  # ps

        if self.scanDir.get():
            print(' Moving Forward')
            pos_0 = float(self.startPoint.get())
            pos_1 = pos_0 + scan_t * c * 1e-12 / 4 * 1000
        else:
            print(' Moving Backward')
            pos_0 = float(self.startPoint.get())
            pos_1 = pos_0 - scan_t * c * 1e-12 / 4 * 1000

        pos_0 = max(pos_0, 0)
        pos_1 = min(pos_1, 300)

        data = []

        for i in range(0, int(exp_param_dict['region_volatile']['Number of Averaging Scans'])):
            data.append([])
            print("scan number:", i + 1)
            j = 0
            if i == 1:
                npdata = np.array(data[0])
                j_length = len(npdata[:, 2])  # the array length all the datasets need to match

            myDLS.stop_DLS()  # in case the stage is still moving
            time.sleep(1)
            myDLS.enable_DLS()
            myDLS.change_velocity(25)
            print('changed velocity to 25')
            myDLS.change_acceleration(100)
            print(f'moving to {pos_0}')
            myDLS.move_stage(pos_0)
            print(f'moving to {pos_0}')
            # time.sleep(5)
            self.waitUntilReached(pos_0)
            print('reached pos_0')
            dataplot = []
            temp_data = []

            myDLS.update()
            time.sleep(1)
            print(f'changing velocity to {vel}')
            print(f'velocity : {vel} ; acceleration : {acc}')
            myDLS.change_velocity(vel * 1.0)
            myDLS.change_acceleration(acc * 1.0)
            time.sleep(1)
            myDLS.update()
            print(f'moving to pos_1 : {pos_1}')
            self.scanData = []
            myDLS.move_stage(pos_1)

            while abs(myDLS.current_pos() - pos_1) > 0.0001:
                # print(myDLS.current_pos())
                temp_pos = np.float(myDLS.current_pos())
                # print(temp_pos, f'destination {pos_1}')
                delay = np.float((temp_pos - pos_0) * 1e-3 / c * 1e12 * 4)

                # t_start = time.perf_counter()

                sample = daq.getSample(demod_path)
                data_x = np.float(sample['x'][0])
                data_y = np.float(sample['y'][0])
                data_phase = np.float(np.arctan(data_y / data_x))
                r = np.float(np.sqrt(sample['x'] ** 2 + sample['y'] ** 2))
                cur_time = np.float(time.perf_counter() - float(exp_param_dict['region_non-volatile']['Init_Time']))
                # t_start = time.perf_counter()

                self.scanData.append([delay, r])
                temp_data.append([temp_pos, delay, data_x, data_y, r, data_phase, cur_time])
                # time.sleep(sampling_time / 100)

                # if i > 0:
                #     x_temp[j] = x_temp[j] + data_x
                # else:
                #     x_temp.append(data_x / float(generated_params['R']))

                dataplot.append((temp_pos, delay, data_x / (i + 1)))
                datax = np.array(dataplot)

                myCanvas.line1.set_data(datax[:, 1], datax[:, 2] / float(generated_params['R']))

                # myCanvas.line1.set_data(datax[:, 1], datax[:, 2])
                if np.amin(datax[:, 2]) != np.amax(datax[:, 2]):
                    myCanvas.ax.set_ylim(np.amin(datax[:, 2] / float(generated_params['R'])),
                                         np.amax(datax[:, 2] / float(generated_params['R'])))
                if np.amin(datax[:, 1]) != np.amax(datax[:, 1]):
                    myCanvas.ax.set_xlim(np.amin(datax[:, 1]), np.amax(datax[:, 1]))

                myCanvas.ax2.set_xlim(myCanvas.ax.get_xlim())
                x2_min = np.min(datax[:, 1])
                x2_max = np.max(datax[:, 1])
                pos_min = np.min(datax[:, 0])
                pos_max = np.max(datax[:, 0])
                ticks = np.linspace(x2_min, x2_max, 6)
                print(self.scanDir.get())
                if self.scanDir.get():
                    ticklabels = np.linspace(pos_min, pos_max, 6)
                else:
                    ticklabels = np.linspace(pos_min, pos_max, 6)[::-1]

                myCanvas.ax2.set_xticks(ticks)
                myCanvas.ax2.set_xticklabels(np.around(ticklabels, 3))

                myCanvas.fig.canvas.flush_events()
                myCanvas.canva.draw()
                j += 1

            time.sleep(1)
            myDLS.update()

            data[i] = temp_data

            # save current scan
            # if len(filetab.overridetext.get()) == 0:
            #     filename = []
            #     for i in range(len(filetab.stringVars)):
            #         if filetab.IntVars[i]:
            #             filename.append(filetab.stringVars.get())
            #     filename = '_'.join(filename)
            # else:
            #     filename = filetab.overridetext.get()
            filename = filetab.current_filename.get()
            file = open(general_param_dict['filepath'] + '/' + filename + str(i) + '_autocorr.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + str(i) + '_autocorr.txt', np.array(data[i]),
                       fmt='%.9e')
            file.close()

            # save avg file
            file = open(general_param_dict['filepath'] + '/' + filename + '_autocorr_avg.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + '_autocorr_avg.txt', datax, fmt='%.9e')
            file.close()
            myCanvas.fig.canvas.flush_events()

        myDLS.update()
        # plt.savefig(general_param_dict['filepath'] + '/' + 'plot' + '.png')
        print("Moving back to initial position.")
        myDLS.change_velocity(25)
        myDLS.change_acceleration(100)
        myDLS.move_stage(pos_0)
        myDLS.disable_DLS()
        # myDLS.di
        # time.sleep(2)
        self.getWidth()

    def getWidth(self):
        if self.typeFlag.get():
            if self.scanData is None:
                print('Scan Data not available')
                return

            data = np.array(self.scanData)
            x_data = data[:, 0]
            y_data = data[:, 1]

            def gaussian(x, amp, mean, sigma):
                return amp * np.exp(-(((x - mean) ** 2) / (2 * sigma ** 2))) / (sigma * np.sqrt(2 * np.pi))

            params, covar = curve_fit(gaussian, x_data, y_data)
            sigma = params[2]
            self.FWHM.set(str(2.355 * sigma))
            self.PulseWidth.set(str(2.355 * sigma / 1.41))
        else:
            if self.scanData is None:
                print('Scan Data not available')
                return

            data = np.array(self.scanData)
            x_data = data[:, 0]
            y_data = data[:, 1]

            def gaussian(x, amp, mean, sigma):
                return amp * np.exp(-(((x - mean) ** 2) / (2 * sigma ** 2))) / (sigma * np.sqrt(2 * np.pi))

            params, covar = curve_fit(gaussian, x_data, y_data)
            sigma = params[2]
            self.FWHM.set(str(2.355 * sigma))
            # self.PulseWidth.set(str(2.355 * sigma / 1.41))
            self.PulseWidth.set('Not Available')

    def crossCorr(self):
        global exp_param_dict, generated_params, general_param_dict, myCanvas
        res = float(self.resolution.get())  # fs
        # time_constant = 0.2  # TC of lock-in, second
        # t_wait = 10 * float(exp_param_dict['region_volatile']['Point Sampling Time'])
        sampling_time = float(exp_param_dict['region_volatile']['Point Sampling Time'])  # second
        c = float(exp_param_dict['region_non-volatile']['c'])

        vel = res * 1e-15 * c / 2 / sampling_time * 1e3
        print(f'resolution is {res}')
        acc = 5 * vel
        scan_t = float(self.scanlen.get())  # ps
        # scan_len_t_after = float(self.scanlenafter.get())
        # v_factor = int(self.vfactor.get())

        # scan_len_t = scan_fine_t + scan_len_t_after
        if self.scanDir.get():
            print(' Moving Forward')
            pos_0 = float(self.startPoint.get())
            pos_1 = pos_0 + scan_t * c * 1e-12 / 4 * 1000
        else:
            print(' Moving Backward')
            pos_0 = float(self.startPoint.get())
            pos_1 = pos_0 - scan_t * c * 1e-12 / 4 * 1000

        pos_0 = max(pos_0, 0)
        pos_1 = min(pos_1, 300)

        data = []

        for i in range(0, int(exp_param_dict['region_volatile']['Number of Averaging Scans'])):
            data.append([])
            print("scan number:", i + 1)
            j = 0
            if i == 1:
                npdata = np.array(data[0])
                j_length = len(npdata[:, 2])  # the array length all the datasets need to match

            myDLS.stop_DLS()  # in case the stage is still moving
            time.sleep(1)
            myDLS.enable_DLS()
            myDLS.change_velocity(25)
            print('changed velocity to 25')
            myDLS.change_acceleration(100)
            print(f'moving to {pos_0}')
            myDLS.move_stage(pos_0)
            print(f'moving to {pos_0}')
            # time.sleep(5)
            self.waitUntilReached(pos_0)
            print('reached pos_0')
            dataplot = []
            temp_data = []

            myDLS.update()
            time.sleep(1)
            print(f'changing velocity to {vel}')
            print(f'velocity : {vel} ; acceleration : {acc}')
            myDLS.change_velocity(vel * 1.0)
            myDLS.change_acceleration(acc * 1.0)
            time.sleep(1)
            myDLS.update()
            print(f'moving to pos_1 : {pos_1}')
            self.scanData = []
            myDLS.move_stage(pos_1)

            while abs(myDLS.current_pos() - pos_1) > 0.0001:
                # print(myDLS.current_pos())
                temp_pos = np.float(myDLS.current_pos())
                # print(temp_pos, f'destination {pos_1}')
                delay = np.float((temp_pos - pos_0) * 1e-3 / c * 1e12 * 4)

                # t_start = time.perf_counter()

                sample = daq.getSample(demod_path)
                data_x = np.float(sample['x'][0])
                data_y = np.float(sample['y'][0])
                data_phase = np.float(np.arctan(data_y / data_x))
                r = np.float(np.sqrt(sample['x'] ** 2 + sample['y'] ** 2))
                cur_time = np.float(time.perf_counter() - float(exp_param_dict['region_non-volatile']['Init_Time']))
                # t_start = time.perf_counter()

                self.scanData.append([delay, r])
                temp_data.append([temp_pos, delay, data_x, data_y, r, data_phase, cur_time])
                # time.sleep(sampling_time / 100)

                # if i > 0:
                #     x_temp[j] = x_temp[j] + data_x
                # else:
                #     x_temp.append(data_x / float(generated_params['R']))

                dataplot.append((temp_pos, delay, data_x / (i + 1)))
                datax = np.array(dataplot)

                myCanvas.line1.set_data(datax[:, 1], datax[:, 2] / float(generated_params['R']))

                if np.amin(datax[:, 2]) != np.amax(datax[:, 2]):
                    myCanvas.ax.set_ylim(np.amin(datax[:, 2] / float(generated_params['R'])),
                                         np.amax(datax[:, 2] / float(generated_params['R'])))
                if np.amin(datax[:, 1]) != np.amax(datax[:, 1]):
                    myCanvas.ax.set_xlim(np.amin(datax[:, 1]), np.amax(datax[:, 1]))

                myCanvas.ax2.set_xlim(myCanvas.ax.get_xlim())
                x2_min = np.min(datax[:, 1])
                x2_max = np.max(datax[:, 1])
                pos_min = np.min(datax[:, 0])
                pos_max = np.max(datax[:, 0])
                ticks = np.linspace(x2_min, x2_max, 6)
                direc = self.scanDir.get()

                if direc:
                    ticklabels = np.linspace(pos_min, pos_max, 6)
                else:
                    ticklabels = np.linspace(pos_max, pos_min, 6)

                myCanvas.ax2.set_xticks(ticks)
                myCanvas.ax2.set_xticklabels(np.around(ticklabels, 3))

                myCanvas.fig.canvas.flush_events()
                myCanvas.canva.draw()
                j += 1

            time.sleep(1)
            myDLS.update()

            data[i] = temp_data

            # save current scan
            # if len(filetab.overridetext.get()) == 0:
            #     filename = []
            #     for i in range(len(filetab.stringVars)):
            #         if filetab.IntVars[i]:
            #             filename.append(filetab.stringVars.get())
            #     filename = '_'.join(filename)
            # else:
            #     filename = filetab.overridetext.get()
            filename = filetab.current_filename.get()

            file = open(general_param_dict['filepath'] + '/' + filename + str(i) + '_crosscorr.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + str(i) + '_crosscorr.txt', np.array(data[i]),
                       fmt='%.9e')
            file.close()

            # save avg file
            file = open(general_param_dict['filepath'] + '/' + filename + '_autocorr_avg.txt', 'wb')
            np.savetxt(general_param_dict['filepath'] + '/' + filename + '_autocorr_avg.txt', datax, fmt='%.9e')
            file.close()
            myCanvas.fig.canvas.flush_events()

        myDLS.update()
        # plt.savefig(general_param_dict['filepath'] + '/' + 'plot' + '.png')
        print("Moving back to initial position.")
        myDLS.change_velocity(25)
        myDLS.change_acceleration(100)
        myDLS.move_stage(pos_0)
        myDLS.disable_DLS()
        # myDLS.di
        # time.sleep(2)
        self.getWidth()
        pass


class Parameters:
    def __init__(self, notebook):
        self.frame = ttk.Frame(notebook)
        self.label = ttk.Label(self.frame, text="Set Experiment Parameters")
        self.label.grid(row=0, column=0, pady=2, padx=2)
        self.entryfields = []
        self.entrylabels = []
        self.stringVars = []
        i = 0
        for param, value in exp_param_dict['region_volatile'].items():
            self.stringVars.append(tk.StringVar(notebook, value))
            self.entryfields.append(ttk.Entry(self.frame, textvariable=self.stringVars[-1]))
            self.entryfields[-1].grid(row=i + 1, column=2, pady=2, padx=2)
            self.entrylabels.append(ttk.Label(self.frame, text=param))
            self.entrylabels[-1].grid(row=i + 1, column=1, padx=2, pady=2)
            i += 1
        self.lockButton = ttk.Button(self.frame, text="Lock Parameters", command=self.lock_params)
        self.lockButton.grid(row=i + 1, column=2, padx=2, pady=2)
        self.frame.pack(fill="both", expand=1)

    def lock_params(self):
        global exp_param_dict
        i = 0
        for param in exp_param_dict['region_volatile'].keys():
            exp_param_dict['region_volatile'][param] = self.stringVars[i].get()
            i += 1
        print(exp_param_dict['region_volatile'])


class LockInTab:
    def __init__(self, notebook):
        self.frame = ttk.Frame(notebook)
        self.label = ttk.Label(self.frame, text="Lock-In Amplifier Settings :")
        self.label.grid(row=0, column=0, pady=2, padx=2)

        # self.label2 = ttk.Label(self.frame, text = "Device ID : ")
        # self.label2.grid(row = 2, column = 0, pady = 2, padx = 2)
        self.entryfields = []
        self.entrylabels = []
        self.stringVars = []
        i = 1
        for param, value in exp_param_dict['Lock In'].items():
            self.stringVars.append(tk.StringVar(notebook, value))
            self.entryfields.append(ttk.Entry(self.frame, textvariable=self.stringVars[-1]))
            self.entryfields[-1].grid(row=i + 1, column=2, pady=2, padx=2)
            self.entrylabels.append(ttk.Label(self.frame, text=param))
            self.entrylabels[-1].grid(row=i + 1, column=1, padx=2, pady=2)
            i += 1
        self.lockButton = ttk.Button(self.frame, text="Update", command=self.lock_params)
        self.lockButton.grid(row=i + 1, column=2, padx=2, pady=2)

        self.label2 = ttk.Label(self.frame, text="INITIALIZATION :")
        self.label2.grid(row=i + 2, column=0)

        self.initializebt = ttk.Button(self.frame, text="Initialize", command=self.initialize)
        self.initializebt.grid(row=i + 2, column=1)

        self.frame.pack(fill="both", expand=1)
        try:
            self.initialize()
        except:
            print('Lock-In not Connected')

    def lock_params(self):
        global exp_param_dict
        i = 0
        for param in exp_param_dict['Lock In'].keys():
            exp_param_dict['Lock In'][param] = self.stringVars[i].get()
            i += 1
        print(exp_param_dict['Lock In'])

    def initialize(self):
        global exp_param_dict, daq, device, demod_path
        (daq, device, _) = zu.create_api_session(exp_param_dict['Lock In']['Device ID'],
                                                 int(exp_param_dict['Lock In']['apilevel']))
        # server_host = exp_param_dict['Lock In']['Server Host']
        demod_path = '/{}/demods/0/sample'.format(device)


class Notebook:
    def __init__(self, root):
        global filetab
        self.notebook = ttk.Notebook(root, width=root.winfo_screenwidth(), height=root.winfo_screenheight())
        self.hardwaretab = DLSTab(self.notebook)
        filetab = FileTab(self.notebook)
        self.scanningtab = Scanning(self.notebook)
        self.parameterstab = Parameters(self.notebook)
        self.corr_tab = Correlation(self.notebook)
        self.lockintab = LockInTab(self.notebook)
        self.notebook.add(filetab.frame, text="Save File")
        self.notebook.add(self.parameterstab.frame, text="Parameters")
        self.notebook.add(self.hardwaretab.frame, text="DLS Control")
        self.notebook.add(self.lockintab.frame, text="Lock-In Amp")
        self.notebook.add(self.scanningtab.frame, text="Scan")
        self.notebook.add(self.corr_tab.frame, text="Correlation")
        self.notebook.pack(pady=10)


class CanvasFrame:
    def __init__(self, root):
        global exp_param_dict
        self.data = []
        self.dataplot = []
        self.x_temp = []
        self.x_dev = []
        self.delay_temp = []
        self.delay_dev = []
        self.pos_temp = []
        self.pos_dev = []
        plt.ion()
        self.fig = plt.figure(figsize=(9, 7))
        self.ax = self.fig.add_subplot(111)
        self.ax2 = self.ax.twiny()
        self.x = np.linspace(0, float(exp_param_dict['region_volatile']['Scanning Period / s']), 100)
        self.y = np.linspace(1e-12, 1e-12, 100)
        plt.ylabel('reflection')
        self.ax.set_xlabel('delay time reference(ps)')
        self.ax.set_ylabel('dR/R')

        self.line1, = self.ax.plot(self.x, self.y, 'b.')  # set x axis range

        self.new_tick_locations = np.linspace(0, float(exp_param_dict['region_volatile']['Scanning Period / s']), 6)

        def tick_function(X):
            global generated_params, exp_param_dict
            generated_params['Scanning Length'] = float(
                exp_param_dict['region_volatile']['Scanning Period / s']) * float(
                exp_param_dict['region_non-volatile']['c']) * 1e-12 / 4 * 1000
            V = float(exp_param_dict['region_volatile']['Initial Position']) + (
                    X / float(exp_param_dict['region_volatile']['Scanning Period / s'])) * generated_params[
                    'Scanning Length']
            return ["%.3f" % z for z in V]

        # add second x-axis for position
        self.ax2.set_xlim(self.ax.get_xlim())
        self.ax2.set_xticks(self.new_tick_locations)
        self.ax2.set_xticklabels(tick_function(self.new_tick_locations))
        self.ax2.set_xlabel('Delay stage position(mm)')

        # start live plotting
        plt.pause(0.01)  # pause the plt for it to update
        self.fig.canvas.flush_events()
        self.canvasFrame = tk.Frame(root)
        self.canvasFrame.pack(side="right")
        self.canva = FigureCanvasTkAgg(self.fig, master=self.canvasFrame)
        self.canva.draw()
        self.canva.get_tk_widget().pack(side="right")
        self.toolbar = NavigationToolbar2Tk(self.canva, self.canvasFrame)
        self.toolbar.update()
        self.canva.get_tk_widget().pack(side="bottom")

    def tick_function(self, X):
        global generated_params, exp_param_dict
        generated_params['Scanning Length'] = float(
            exp_param_dict['region_volatile']['Scanning Period / s']) * float(
            exp_param_dict['region_non-volatile']['c']) * 1e-12 / 4 * 1000
        V = float(exp_param_dict['region_volatile']['Initial Position']) + (
                X / float(exp_param_dict['region_volatile']['Scanning Period / s'])) * generated_params[
                'Scanning Length']
        return ["%.3f" % z for z in V]


global root


def Gui(new_dict1, new_dict2):
    global general_param_dict, exp_param_dict, myCanvas, root
    general_param_dict = new_dict1
    exp_param_dict = new_dict2
    root = tk.Tk()
    root.title("Pump-Probe GUI")
    myCanvas = CanvasFrame(root)
    Notebook(root)

    def on_closing():
        if messagebox.askyesno("Quit", "Save Parameters before closing?"):
            my_list = [general_param_dict, exp_param_dict]
            with open("config.pkl", 'wb') as f:
                pickle.dump(my_list, f)
            root.destroy()
        else:
            root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()

import numpy as np
from scipy.stats import linregress
from scipy.interpolate import PchipInterpolator
import elliptec
import tlccs_lib
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Spectrometer
spect_model = tlccs_lib.SPECTROMETER_MODELS['0x8087']
mySerial = 'M00301586'
# mySerial = 'M00301586'

# Linear Stage
LS_com = 'COM11'

print('Available Spectrometer Models : \n', tlccs_lib.SPECTROMETER_MODELS)
print('Using', spect_model)

abs_min = 20
abs_max = 60
integration_time = 0.01
slope = None
intercept = None
calFunc = None

global ls, tl_spec


def initSpect():
    global tl_spec
    tl_spec = tlccs_lib.tlccs(spect_model, serial=mySerial)
    tl_spec.setIntegrationTime(integration_time)
    print(f'Integration Time : {tl_spec.getIntegrationTime()}')
    # print(getCenterWavelength())
    return tl_spec


def initLS():
    global ls
    controller = elliptec.Controller(LS_com)
    ls = elliptec.Linear(controller)
    # ls.home()
    ls.set_distance(abs_min)
    # for d in [0,10,20,30,50,60]:
    #     ls.set_distance(d)
    #     time.sleep(1)
    # ls.home()
    return ls


# Linear Stage Go-To Control
def setLSdistance(dist):
    ls.set_distance(dist)
    print('Current Distance :', dist)
    pass


# Get Center Wavelength from the Spectrometer
def getCenterWavelength():
    # Scanning Function
    print('Getting Wavelength')
    tl_spec.startScan()
    time.sleep(0.1)
    intensity_data = tl_spec.getScanData()
    wavelength_data = tl_spec.getWavelengthData()
    # print(intensity_data)
    wavelengths, intensities = np.array(wavelength_data), np.array(intensity_data)
    index_maxintensity = np.argmax(intensities)

    def gaussian(x, amp, mean, sigma):
        return amp * np.exp(-(((x - mean) ** 2) / (2 * sigma ** 2))) / (sigma * np.sqrt(2 * np.pi))

    params, covar = curve_fit(gaussian, wavelengths, intensities)
    mean = params[1]
    peak_wavelength2 = mean
    peak_wavelength = wavelengths[index_maxintensity]
    # wavelength = None
    return peak_wavelength
    # return peak_wavelength2


# Fitting Function
def fitLinear(x, y):
    result = linregress(x, y)
    return result.slope, result.intercept


def fitPoly(x, y, order):
    z = np.polyfit(y, x, order)
    z2 = np.polyfit(x, y, order)
    fitFunc = np.poly1d(z)
    distToLambda = np.poly1d(z2)
    return fitFunc, distToLambda


def fitMonotoneCubicSpline(x, y):
    fitFunc = PchipInterpolator(x, y)
    invFunc = PchipInterpolator(y, x)
    return fitFunc, invFunc


# Calibration Loop
def calibrate(n_points=20):
    global slope, intercept
    pos = []
    wavelengths = []
    # ls.home()
    ls.set_distance(abs_min)
    time.sleep(1)
    order = 3
    i = 0
    for pos_temp in np.linspace(abs_min, abs_max, n_points):
        print(i)
        i += 1
        setLSdistance(pos_temp)
        time.sleep(1)
        act_dist = ls.get_distance()
        pos.append(act_dist)
        wavelength = getCenterWavelength()
        wavelengths.append(wavelength)
        # wavelengths.append(0)

    pos = np.array(pos)
    wavelengths = np.array(wavelengths)

    slope, intercept = fitLinear(pos, wavelengths)
    calFunc, invFunc = fitPoly(pos, wavelengths, order)
    # splineFunc, invSplineFunc = fitMonotoneCubicSpline(pos, wavelengths)

    x_points = np.linspace(abs_min, abs_max, 100)
    y_points = slope * x_points + intercept

    # Plotting
    plt.scatter(pos, wavelengths, s='0.2')

    plt.plot(x_points, y_points, label='Linear Fit')
    plt.plot(x_points, invFunc(x_points), label=f'Order : {order} fit')
    # plt.plot(x_points, invSplineFunc(x_points), label='Monotone Cubic Spline')

    plt.xlabel('Linear Stage Travel (mm)')
    plt.ylabel('Peak Wavelength (nm)')
    plt.title('')

    plt.legend()
    plt.show()

    return slope, intercept
    # return calFunc


# After Calibration
def setWavelength(wavelength):
    if slope is not None:
        position = (wavelength - intercept) / slope
        if abs_min <= position <= abs_max:
            setLSdistance(position)
        else:
            print('Out of Bounds')
        return position
    else:
        print('Slope not found')
        return

    # if calFunc is not None:
    # position = calFunc(wavelength)
    # if abs_min <= position <= abs_max:
    #     setLSdistance(position)
    # else:
    #     print('Out of Bounds')
    # return position


def inputLoop():
    while True:
        wavelength = input('Set Wavelength to : (nm) [Enter "x" to close ]')
        if wavelength == 'x':
            # print(f'Position : {pos} mm \n Wavelength : {wavelength} nm')
            break
        wavelength = float(wavelength)
        pos = setWavelength(wavelength)
        print(f'Position : {pos} mm \n Wavelength : {wavelength} nm ')


def test_accuracy():
    test_wavelengths = np.linspace(600, 950, 20)
    measured_wavelengths = []
    for wavelength in test_wavelengths:
        setWavelength(wavelength)
        measured_wavelengths.append(getCenterWavelength())
    measured_wavelengths = np.array(measured_wavelengths)
    deviations = np.abs(measured_wavelengths - test_wavelengths)
    mean_deviation = np.mean(deviations)
    std_deviation = np.std(deviations)
    print('Mean of Deviations :', mean_deviation)
    plt.hist(deviations)
    plt.show()
    # plt.scatter()


if __name__ == "__main__":
    global ls, tl_spec
    try:
        ls = initLS()
    except:
        print('Linear Stage Connection Failed')

    try:
        tl_spec = initSpect()
    except:
        print('Spectrometer Connection Failed')

    try:
        slope, intercept = calibrate()
        # calFunc = calibrate()
    except:
        print('Not Calibrated')

    print(f'Slope = {slope} ; Intercept = {intercept}')

    # if slope is not None:
    #     inputLoop()
    # else:
    #     print('Exiting')

    # ls.home()

print('Script Ended')

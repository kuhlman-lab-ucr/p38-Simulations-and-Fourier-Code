from tkinter import*
from tkinter import ttk
import scipy.optimize
from ttkthemes import ThemedTk
from Bookshelf import*
import ttkbootstrap as ttkb
from ttkbootstrap.constants import*

import os
import math
import time
import openpyxl
import re
import serial
import numpy as np
import scipy
from scipy import stats

################################
#     Functions for RS-232 communicaiton of Pumpit peripumps

def connect_serial(port_ID):
    ser = serial.Serial(
        port=port_ID,
        baudrate=19200,
        parity=serial.PARITY_NONE,
        stopbits=serial.STOPBITS_ONE,
        bytesize=serial.EIGHTBITS
    )
    return ser

def Debug_read(pump_list, comport):                   #read from the comport and return the raw data as well as 
    port = pump_list[comport]['Connection']      # decoded and cleaned version of the reading
    the_read = {}          
    ret = port.read_all()
    the_read.update({'raw':ret})
    if len(ret) != 0:
        bya = list(ret)
        the_read.update({'Byte Array':bya})
        first = bya.pop(0)
        the_read.update({'First':first})
        last = bya.pop()
        the_read.update({'Last':last})
        res = ret.decode('ascii')
        the_read.update({'ascii':res})
        res_clean = bytes(bya).decode('ascii')
        the_read.update({'Cleaned':res_clean})
    return the_read

def Debug_write(pump_list, comport, command):               #write to a pump in a the correct format
    ser_con = pump_list[comport]['Connection']
    pump_number = pump_list[comport]['Pump Number']
    build = str(pump_number) + str(command) + '*\r'   # needs the com port and the pump number
    mess = bytes(build ,'utf-8')
    ser_con.write(mess)
    return build

def Debug_query(pump_list, comport,command):           #write a command and read the response after a small pause
    sent = Debug_write(pump_list, comport ,command)
    time.sleep(0.1)              #delay very nessisary
    rec = Debug_read(pump_list,comport)
    rec.update({'Sent': sent})
    return rec


######################################################################

def process_responce(ret): #functon to parse peripump responses into usable strings
    response = {}
    bya = list(ret)
    bya.pop(0) # take off first and last hexidecimal
    bya.pop()
    res_clean = bytes(bya).decode('ascii')  #turn into a string

    ads = re.search(r'\d\d[a-zA-Z]', res_clean) #header, adress (00-99) and status Letter
    promts = re.search(r'[a-zA-Z]', ads.group(0))   #isolate Status letter
    promt =  promts.group(0)

    res_str = re.split(r'\d\d[a-zA-Z]', res_clean.strip())

    if type(res_str) == list:  #just cleaning up weird junk responses to something remotely useful
        usable = ''
        for l in range(len(res_str)):
            if str(res_str[l]) == '':
                pass
            elif str(res_str[l]) == '\x03\x02':
                pass
            elif '\x03\x02' in str(res_str[l]):
                cut = str(res_str[l]).replace('\x03\x02', '')
                usable = usable + cut
            else:
                usable = usable + str(res_str[l])

        
    #
    response.update({'Address': ads.group(0)})
    #
    response.update({'Status': promt})
    #
    response.update({'Cleaned':usable})


    return response



def write(test_con, command): #Write Command com ports
    comport = list(test_con.keys())[0]
    ser_con = test_con[comport]['Connection']
    pump_number = test_con[comport]['Pump Number']
    build = str(pump_number) + str(command) + '*\r'   # needs the com port and the pump number
    mess = bytes(build ,'utf-8')
    ser_con.write(mess)
    return build

def read(test_con):                   #read from the comport and return the raw data as well as 
    comport = list(test_con.keys())[0]
    port = test_con[comport]['Connection']      # decoded and cleaned version of the reading
    the_read = {}          
    ret = port.read_all()
    the_read.update({'raw':ret})
    if len(ret) != 0:
        bya = list(ret)
        the_read.update({'Byte Array':bya})
        first = bya.pop(0)
        the_read.update({'First':first})
        last = bya.pop()
        the_read.update({'Last':last})
        res_clean = bytes(bya).decode('ascii')

        useable = process_responce(ret) #function to parse response
        #
        the_read.update({'Address':useable['Address']})
        #
        the_read.update({'Status':useable['Status']})
        #
        the_read.update({'Cleaned':useable['Cleaned']})

    else:
        the_read.update({'Cleaned':"No return"})
        the_read.update({'Status':" "})
    return the_read


def query(test_con, command):           #write a command and read the response after a small pause
    sent = write(test_con, command)
    time.sleep(0.1)              #delay very nessisary
    rec = read(test_con)
    rec.update({'Sent': sent})
    return rec



def connect_com_ports(pump_list, portlist): #attempts to get adress to verify come port connects to pump
    for i in range(len(portlist)):
        try:
            connection = connect_serial(portlist[i])
            mess = bytes('*ADR\r' ,'utf-8')
            connection.write(mess)
            time.sleep(0.1)
            raw_pn = connection.read_all()
            bya = list(raw_pn)
            try:
                bya.pop(0)
                bya.pop()
                cleaned_pn = bytes(bya).decode('ascii')
            
                pn = int((cleaned_pn.split('P'))[1])
                pump_list.update({portlist[i]:{'Connection':connection, 'Pump Number': pn}})
            except IndexError:
                #connection.close()
                pass
        except serial.SerialException:
            pass
            #connection.close()
    return pump_list
    

def make_beep(pump_list, comport):
    pump_list[comport]['Connection'].write(b'00BEP*\r')
    time.sleep(0.1)
    response = Debug_read(pump_list, comport)
    return response
    
##########


#############

def com_check(pump_list, comport): #check the values of the com connection
    port = pump_list[comport]['Connection']
    checks = {}
    checks.update({'cts': str( port.cts )})
    checks.update({'dsr': str( port.dsr )})
    checks.update({'ri': str( port.ri )})
    checks.update({'cd': str( port.cd )})
    checks.update({'bytesize': str( port.bytesize )})
    checks.update({'parity': str( port.parity )})
    checks.update({'stopbits': str( port.stopbits )})
    checks.update({'timeout': str( port.timeout )})
    checks.update({'xonxoff': str( port.xonxoff )})
    checks.update({'rtscts': str( port.rtscts )})
    checks.update({'readable': str( port.readable )})
    checks.update({'in_waiting': str( port.in_waiting )})
    return checks
    

def reset(pump):
    pump.write(b'*RESET\r')
    res = str(pump.read_all())


def close_all(pumps):
    worked = {}
    for i, value in pumps.items():
        pumps[i]['Connection'].close()
        worked.update({i: pumps[i]['Connection'].is_open})
    return worked
    
# Dictionary of all info needed to run experiment
experiment_setting = ({
    'Blank Media Pump': str,
    'Solution Pump': str,
    'Experiment Flow rate': float,
    'Frequency (mHz)': float,
    'Experment Time (S)': int,
    'Time Incriments (S)': float,
    'Flow Settings Set': {}

})

def make_flow_values(experiment_setting):
    total_time = experiment_setting['Experment Time (S)']
    t_incr = experiment_setting['Time Incriments (S)']
    flow_max = experiment_setting['Experiment Flow rate']
    freq_hz = float(experiment_setting['Frequency (mHz)'])/ 1000.0
    wavelen = 1.0/freq_hz
    set_array = {}
    

    return set_array


def run_experimnet_2pump(pump_list, experiment_setting):
    experiment_setting_conferm = {}
    return experiment_setting_conferm


def vol_RfT(Rate, t1, t2 ): #return a volume from a  rate and a change of time
    delta = round((t2 - t1), 3)
    v = abs(delta* Rate)
    return v


###########################################################################################################
######### Signal Calculation and simulation

#Current defaults Volume in mL, time in Seconds, rates in mL/s
# Concentrations technically unit-less functionally a percentage (0-100)
'''
Constant_dic = ({
    "V_m": float, #Volume of mixing chamber in mL
    "V_e": float, #Volume of expermiment chamber in mL
    "V_b": float, #Volume of Blank Solution entering mixing chamber in mL
    "V_f": float, #Volume of Full Solution entering mixing chamber in mL
    "V_mr": float, #Volume of mixing chamber removed in mL
    "V_in": float, #Volume entering experminment chamber in mL
    "V_out": float, #Volume exiting Experminent Chamber in mL
    "R_b": float, #Rate of blank solution pump's flow in mL/S
    "R_f": float, #Rate of  full solution pump's flow in mL/S
    "R_in": float, #Rate of pump going from mixing chamber to experiment chamber in mL/S
    "R_out": float, #Rate of pump going from expiremtn chamber to waste in mL/S
    "C_b": float, # Concentration of blank Solution
    "C_f": float, # Concenteation of full solution
    "Phase Delay (S)": float, #Delay of concentation from the mixing chamber to the Experiment chamber
    "Phase Delay Places": int, #Delay in places back of the generated data set of the phase
    "Time Delta": float #time between each input t value
})
'''
def phase_tube(ID, length, rate): #returns the phase delay of the function in seconds
    #Inputs should me in milimeters to convert to ml 
    #since the rates are in mL per sec
    A = ( (ID/2.0)*(ID/2.0) ) * np.pi
    v_mm3 = A * length
    v_mL = v_mm3 / 1000.0 #Convert from cubic mm to mL
    t_raw = v_mL / rate
    t = round((v_mL/rate), 2)
    return t

def phase_delay(Phase_time, delta_t):   #Deterime the delay as a number of places back in data list
    delay = int
    calc = round(Phase_time / delta_t)
    delay = int(calc)
    return delay

####################################
## this is the step in the proccess (for the program) where
## the addition of blank of concentrate is determined
def M_volumes(Constant_dic, pump_binary): #function for mixing chamebr to turn rates and time delta into volumes
    R_b = Constant_dic["R_b"]
    R_f = Constant_dic["R_f"]
    R_in = Constant_dic["R_in"]
    delta_t = Constant_dic["Time Delta"]

    #Determin if blank or full being pumped into mixing chamber
    #set up  is currently either/or so if one is runnign the other is not
    if pump_binary >= 0.0:
        V_b = R_b * delta_t
        V_f = 0.0
    else:
        V_b = 0.0
        V_f = R_f * delta_t
    V_mr = R_in * delta_t

    #set the volumes that have entered and existed this delta
    Constant_dic["V_b"] = V_b
    Constant_dic["V_f"] = V_f
    Constant_dic["V_mr"] = V_mr

    return Constant_dic

def E_volumes(Constant_dic):   #deterime volume of experiment chamber
    R_in = Constant_dic["R_in"]
    R_out = Constant_dic["R_out"]
    delta_t = Constant_dic["Time Delta"]

    V_in = R_in * delta_t
    V_out = R_out * delta_t

    Constant_dic["V_in"] = V_in
    Constant_dic["V_out"] = V_out

    return Constant_dic



def mix_concentration(point_dic, Constant_dic, Current_Time):
    #generate dictionary for current data point being made
    data_point = ({'Time': Current_Time, 
                   'Mix Cons': float, 
                   'Mix Vol': float, 
                   'Exp Cons':float, 
                   'Exp Vol':float})

    C_b = Constant_dic["C_b"]
    C_f = Constant_dic["C_f"]
    
    #system starts off assuming everything is primed blank
    if point_dic == []:
        C_prev = C_b
    else:
        C_prev = point_dic[-1]["Mix Cons"]

    V_b = Constant_dic["V_b"]
    V_f = Constant_dic["V_f"]
    V_mr = Constant_dic["V_mr"]
    V_m = Constant_dic["V_m"]

    #if structure to handle if a rate imbalence emptys out a chamber
    # basically preventing a divide by zero error
    if (V_m - V_mr + V_b + V_f) < 0.0:
        V_m_now  = 0.0
        if V_b >= V_f:
            con_now = C_b
        else:
            con_now = C_f
    else:
        V_m_now  = V_m - V_mr + V_b + V_f
        con_now = ((C_b * V_b) + (C_f * V_f) + (C_prev * (V_m - V_mr))) / V_m_now

    #Update the constant for future equations in case theres a change
    Constant_dic["V_m"] = V_m_now

    #Load data points, Time and these two are the three things this fuction creates
    # along with the dictionary for this point itself
    data_point['Mix Cons'] = con_now
    data_point['Mix Vol'] =  V_m_now

    #Add data point to list
    point_dic.append(data_point)
    #return data point list and constants
    return point_dic, Constant_dic


def Experiment_Concentration(point_dic, Constant_dic):
    #current data point being established
    current_data = point_dic[-1]
    if len(point_dic) <= 1: #Assumes starting concentration is zero (blank)
        e_con_prev = 0.0
    else:
        e_con_prev = point_dic[-2]['Exp Cons']
    #because the mixing chamber function delcares the dictionary for this current time point in the list
    # List[-1] in the that function is previous time point BUT List[-2] in THIS function is previous time
    #  point and list[-1] is current time point being added to

    #Get Phase delayed concentration
    Phase_delay = Constant_dic["Phase Delay Places"]
    if len(point_dic) < (Phase_delay + 1):  #the Plus ones or becasue THIS data point is list[-1]
        mix_con_now = 0.0
    else:
        mix_con_now = point_dic[-(Phase_delay + 1)]['Mix Cons']
    
    #Volume in experiment chamber
    V_in = Constant_dic["V_in"]
    V_out = Constant_dic["V_out"]
    V_e = Constant_dic["V_e"]

    if (V_e - V_out + V_in) < 0.0:
        V_e_now = 0.0
        con_e_now = 0.0
    else:
        V_e_now = V_e - V_out + V_in
        con_e_now = ((mix_con_now * V_in) + (e_con_prev  * (V_e -V_out))) / V_e_now

    #Update Constant in case of change
    Constant_dic["V_e"] = V_e_now
    #Load final two data points for this dictionary of the lsit
    current_data["Exp Cons"] = con_e_now
    current_data["Exp Vol"] = V_e_now

    return point_dic, Constant_dic

#function to simulate the concnetration over time of a current test set up
def Generate_simulation(Constant_dic, Total_time, period):
    #First Genrate all the time (x axis) points for the simmulation
    delta_t = float(Constant_dic["Time Delta"])
    if type(Total_time) == (int or float):
        times_x = np.linspace(0, int(Total_time), int(Total_time*int(1.00/delta_t)))
        #From period input determine sine function for on/offs
        # sin((2pi/period) * t) = sine wave with period 
        theta = 2*3.14/period   #Normalize sine wave input to output target period
        x_sin = []
        for s in range(len(times_x)):
            x_sin.append(theta *  times_x[s])
        y_sin = np.sin(x_sin) #sine wave values at set period
    else:# type(Total_time) == dict:
        times_x = Total_time["Time"]
        y_sin = Total_time["Signal"]


    #Current the M_volume function uses these numbers to deterime when to switch pumps
    # it uses a simple check:
    #       sine greater than or equal to zero is blank pump
    #       sine less than zero is full pump

    #With all time values to be used and on/off setting for the pumps
    # now we generate the simulation data

    point_dic = [] #the eventual list of data points
    #list of dictionarys
    # {'Time': float, 'Mix Cons': float, 'Mix Vol': float, 'Exp Cons':float, 'Exp Vol':float}

    for i in range(len(times_x)):   #for every time generated to measure
        
        #establish volumes goin in and out of mixing chamber
        Constant_dic = M_volumes(Constant_dic, y_sin[i])

        #Deterimine concnetration of mixing chamber at current time
        # Additional sets up dictionary of current data point
        point_dic, Constant_dic = mix_concentration(point_dic, Constant_dic, times_x[i])
    
        #Determine Volumes moving in and out of experiment chamber
        Constant_dic = E_volumes(Constant_dic)

        #Deterim Concentration of Experiment chamber
        point_dic, Constant_dic = Experiment_Concentration(point_dic, Constant_dic)

    #return completeded simmulation data for graphing and analysis
    # and sine values used to switch pumps
    return point_dic, y_sin

######################################################################################
# Experimental functions

General_constants = ({
    "Volume Container":float,
    "Rate In": float,
    "Rate Out": float,
    "Concentration In": float
})

def General_concentration_change(point_dic, General_constants, Current_Time):
    #generate dictionary for current data point being made
    data_point = ({"Time": Current_Time, 
                   "Concentration": float, 
                   "Volume": float, 
                   })

    C_in = General_constants["Concentration In"]
    Rate_in = General_constants["Rate In"]
    Rate_out = General_constants["Rate Out"]
    Volume_C = General_constants["Volume Container"]

    #
    time_delta = float
    if point_dic == []:
        C_prev = 0.0
        time_delta = Current_Time
    else:
        C_prev = point_dic[-1]["Concentration"]
        time_delta = Current_Time - point_dic[-1]["Time"]

    #Calculate Volumes
    v_in = Rate_in * time_delta
    v_out = Rate_out * time_delta

    #if structure to handle if a rate imbalence emptys out a chamber
    # basically preventing a divide by zero error
    if (v_in + Volume_C - v_out) < 0.0:
        V_m_now  = 0.0
        con_now = 0.0
    else:
        V_m_now  = v_in + Volume_C - v_out
        con_now = ( (v_in * C_in) + (C_prev * (Volume_C - v_out))) / V_m_now

    data_point["Concentration"] = con_now
    data_point["Volume"] = V_m_now
    General_constants["Volume Container"] = V_m_now

    #Add data point to list
    point_dic.append(data_point)
    #return data point list and constants
    return point_dic, General_constants

#Function that gives out a new concetration
def concentration_change(current_conc, Current_vol, concentration_in, rate_in, rate_out, delta_t):
    v_in = rate_in * delta_t
    v_out = rate_out * delta_t
    v_new = Current_vol - v_out + v_in
    concentration_new = (concentration_in * v_in) + (current_conc * (Current_vol - v_out))/ v_new
    return concentration_new, v_new

def Mix_concentration_change(current_conc, Current_vol, con1,  rate1, con2, rate2, rate_out, delta_t):
    v1 = rate1 * delta_t
    v2 = rate2 * delta_t
    v_out = rate_out * delta_t
    v_new = Current_vol - v_out + v1 + v2
    concentration_new = ((con1 * v1) + (con2 * v2) + (current_conc * (Current_vol - v_out)))/ v_new
    return concentration_new, v_new

def get_linear_trend(x_set, y_set, r_lim): #
    end = {'Slope':float, 'Intercept': float, 'R': float, 'Begin x': 0 , 'Begin y': 0 ,'End x': 0, 'End y': 0}

    i_h = 1
    x_test = x_set[5:-i_h]
    y_test = y_set[5:-i_h]
    
    slope, inter, r_value, p_val, stder = stats.linregress(x_test, y_test)
    while (r_value < r_lim) and ((i_h+6) < len(x_set)):
        x_test = x_set[5:-i_h]
        y_test = y_set[5:-i_h]
        slope, inter, r_value, p_val, stder = stats.linregress(x_test, y_test)
        i_h = i_h + 1
        
    
    x_final = x_set[5:-i_h]
    y_final = y_set[5:-i_h]
    slope, inter, r_value, p_val, stder = stats.linregress(x_final, y_final)
    
    end['Slope'] = slope
    end['Intercept'] = inter
    end['R'] = r_value
    end['Begin x'] = 5
    end['Begin y'] = 5
    end['End x'] = len(x_test) - i_h
    end['End y'] = len(y_test) - i_h
    
    return end

#Generate a KL divergence graph from two fourier graphs
def KL_div(xf, yf, compair_x, compair_y):
    x_out = []
    compaired = []
    if len(yf) < len(compair_y):
        for i in range(len(yf)):
            x_out.append(xf[i])
            base = yf[i]
            comp = compair_y[i]
            compaired.append(scipy.special.kl_div(base, comp))
    else:
        for i in range(len(compair_y)):
            x_out.append(compair_x[i])
            base = yf[i]
            comp = compair_y[i]
            compaired.append(scipy.special.kl_div(base, comp))

    return compaired, x_out

def a_linear(x,m,b):
    return ((float(m)*float(x))+float(b))

def an_exponential(x,a,t,b): #Exponential increase
    return a*np.exp(-t*x)+b

def Exponential_fit(x_set, y_set): #Fit an exponential function
    # fit to y = m * exp(t*x) + b
    a_start = 1.0
    t_start = 1.0
    b_start = y_set[0]
    p0 = (a_start,t_start, b_start)
    params, cv = scipy.optimize.curve_fit(an_exponential, x_set, y_set, p0)
    perr = np.sqrt(np.diag(cv))
    Results = ({
        'A': params[0],
        't': params[1],
        'b': params[2],
        'Fit': perr
    })
    return Results

def linear_fit(x_set, y_set):  #fit a linear funtion
    # fit to y = m*x + b
    params, cv = scipy.optimize.curve_fit(a_linear, x_set, y_set)
    perr = np.sqrt(np.diag(cv))
    Results = ({
        'm': params[0],
        'b': params[1],
        'Fit': perr
    })
    return Results

###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################



def Concentration_change(con_now, con_in, vol_change, vol_total): 
    #get new concetration from a total volume of a now concentration with a change volume removed
    # and a replaced to total with an incoming concentration  
    con_new = ( (con_now * (vol_total-vol_change)) + (con_in * vol_change) ) / vol_total
    return con_new


def Pump_volume(test_con, Volume): #trigger the Pump to dispense a set volume
    if float(Volume) <= 0.0:
        rec = ""
    else:
        command = "VOL " + str(Volume)
        sent = write(test_con, command)
        time.sleep(0.1)              #delay very nessisary
        rec = read(test_con)
        time.sleep(0.1)
        rec.update({'Sent': sent})

        sent = write(test_con, "RUN")
        time.sleep(0.1)
        rec_2nd = read(test_con)
        rec.update({'Sent Run': sent})

    return rec


def Pump_volume_Dual(Pump_con_1, Volume_1, Pump_con_2, Volume_2): #Special Program to trigger two pumps to dispence an input amount of volume
    #if one of the volumes is 0, the volume is set but the run command is not
    command_1 = "VOL " + str( round(Volume_1, 3) )  #rounding to 3 places to prevent pump errors
    command_2 = "VOL " + str( round(Volume_2, 3) )

    sent = write(Pump_con_1, command_1)
    sent = write(Pump_con_2, command_2)
    
    #messages can be written or read from seperate pumps
    time.sleep(0.1)  
    rec = []

    if float(Volume_1) <= 0.0:
        rec.append(read(Pump_con_1))
    else:
        sent = write(Pump_con_1, "RUN")

    if float(Volume_2) <= 0.0:
        rec.append(read(Pump_con_2))
    else:
        sent = write(Pump_con_2, "RUN")
        
    time.sleep(0.1)

    rec.append(read(Pump_con_1))
    rec.append(read(Pump_con_2))

    return rec


def wave_set_v3(Concentration_set): #Create settings for

    amp = Concentration_set['C2'] 
    t = np.linspace(0, Concentration_set['Period (s)'], 21)
    interval = int(t[1]-t[0])
    ct = {}
    omg = (2.0 * np.pi) / Concentration_set['Period (s)']
    for i in range(len(t)): #Create the steps in concentration needed to create the wave
        raw = np.sin( ( omg * t[i] ) - ( np.pi / 2.0 ) )
        shift = (raw + 1.0) / 2.0
        scale = amp * shift
        ct.update({t[i]:scale})

    #determing the total amount of media to be added every 
    #step of the interal
    ve = round( Concentration_set['Use Rate'] * interval , 3 )
    vf = ve + Concentration_set['Vb mL']
   
   #create the list of  volumes of each solution at each interval neeed
   #to be added in order to create the wave function
    orders = []
    for i in range(len(t)): 
        if i == 0: #create the first entry
            v2 = round( ( (vf * ct[t[i]]) /  amp ) ,  3  )
            v1 = round(ve - v2, 3)
            con = round( ( (amp * v2) / vf ), 2)
            orders.append({'Time': interval, 'v1': v1, 'v2':v2, 'Con': con })
        else:   #iterations after the first use info from the previous entry on the list
            v2 = round( ( ( (vf * ct[t[i]]) - (orders[-1]['Con'] * Concentration_set['Vb mL'] ) ) /  amp ),  3  )
            if v2 >= ve:
                v2 = round(ve, 3)
            elif v2 <= 0.0:
                v2 = 0.00
            else:
                v2 = v2
            v1 = round(ve - v2, 3)
            con = round( ( ( (orders[-1]['Con'] * Concentration_set['Vb mL'] ) + (amp * v2) ) / vf ),  3)
            orders.append({'Time': interval, 'v1': v1, 'v2':v2, 'Con': con })

    return orders


###########################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################

def make_wave_set(Period, rate, set_delay): #Create a list of full and blank rates matching a sine wave
    samps = round(float(Period) / float(set_delay))
    runset = np.linspace(0, Period, samps)
    total_set = {'Full':[], 'Blank':[]}
    adjust_rate = rate * 0.5    #input is half total rate to compensate for cos wave doubling to make range above zero
    k = (2*np.pi) / Period
    for i in range(len(runset)):
        f_r = (-1.0*(np.cos(k*runset[i])) + 1.0)
        f = round((f_r * adjust_rate), 3)
        b = round((rate - f), 3)
        total_set['Full'].append(f)
        total_set['Blank'].append(b)

    return total_set

def make_complex_wave_set(x, y, max_rate, set_delay): #Create a list of full and blank rates matching an input wave
    total_set = {'Full':[], 'Blank':[]}
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    y_adjusted = []
    if y_min < 0.0: #raise the wave above zero
        adj1 = -1.0 * y_min
        adj2 = y_max + adj1
        for a in range(len(y)):
            y_adjusted.append(y[a] + adj1)
    else:
        y_adjusted = y
        adj2 = y_max
    
    y_real = []
    for r in range(len(y_adjusted)):
        if y_adjusted[r] == 0.0:    #avoid divide by zero error
            y_real.append(0.0)
        else:
            y_real.append(round(max_rate*(adj2/y_adjusted[r]), 3))  #round to 0.001, smallest mL per sec pumps can be set

    
    for i in range(len(y_real)):
        f = y_real[i]
        b = max_rate - f

        total_set['Full'].append(f)
        total_set['Blank'].append(b)
        

    return total_set

def generate_data_dic(time):
    data_point = ({'Time': time, 
                   'Mix Cons': float, 
                   'Mix Vol': float, 
                   'Exp Cons': float, 
                   'Exp Vol':float})
    return data_point


def Generate_ratiod_simulation(wave_set, e_vol, m_vol, rate_in, bcon, fcon, phase_time, sim_time, sim_delta):

    x = np.linspace(0, sim_time, int( sim_time / sim_delta) )
    b_y = wave_set["Blank"]
    f_y = wave_set["Full"]
    delta_t = float(wave_set["Time"][1] - wave_set["Time"][0])
    start_con = 0.0
    point_dic = [] 
    delay = phase_delay(phase_time, sim_delta)
    

    
    for i in range(len(x)):   #for every time generated to measure
        data_point = generate_data_dic(x[i])
        #establish volumes goin in and out of mixing chamber
        if i == 0:   
            con_m_new, v_m_new = Mix_concentration_change(start_con, m_vol, bcon, b_y[i], fcon, f_y[i], rate_in, sim_delta)
            con_e_new, v_e_new = concentration_change(start_con, e_vol, start_con, rate_in, rate_in, sim_delta)
        else:
            con_m_new, v_m_new = Mix_concentration_change(point_dic[-1]['Mix Cons'], point_dic[-1]['Mix Vol'], bcon, b_y[i], fcon, f_y[i], rate_in, sim_delta)
            if len(point_dic) <= delay:
                con_e_new, v_e_new = concentration_change(start_con, e_vol, start_con, rate_in, rate_in, sim_delta)
            else:
                con_e_new, v_e_new = concentration_change(point_dic[-1]['Exp Cons'], point_dic[-1]['Exp Vol'], point_dic[-delay]['Mix Cons'], rate_in, rate_in, sim_delta)
                                    
        data_point['Mix Cons'] = con_m_new
        data_point['Mix Vol'] = v_m_new
        data_point['Exp Cons'] = con_e_new 
        data_point['Exp Vol'] = v_e_new
                              
        point_dic.append(data_point)

        
    #return completeded simmulation data for graphing and analysis
    # and sine values used to switch pumps
    return point_dic

def just_sine(x, period): #function to just throw out a sine value with a specific period
    w = (2.00*np.pi)/period
    n  = round( np.sin(w*x), 3) #in this context 3 places is more than enough
    return n

def sine_val(amp, per, phase, time): #return a sine value
    w = (2.0*np.pi) / per
    pc = time - phase
    s = amp * np.sin(w *pc)
    return round(s, 3)

def sine_wave(amp, per, phase, time_set): #Make the x and y outputs of a sine wave from a list of time valuea and a period
    y = []
    for i in range(len(time_set)):
        y.append(sine_val(amp, per, phase, time_set[i]))
    return y

def Sine_combine(dic_set):
    sine_total = []
    for s in range(len(dic_set)):
        if s == 0:
            first = dic_set[s]
            for i in range(len(first)):
                sine_total.append(float(first[i]))
        else:
            after = dic_set[s]
            for i in range(len(after)):
                sine_total[i] = sine_total[i] + float(after[i])
    return sine_total



######################## fitting

def concentration_change_funtion(time, rate_in_out, concentration_in, current_conc, Current_vol):
    v_in = rate_in_out * time
    v_out = rate_in_out * time
    v_new = Current_vol - v_out + v_in
    concentration_new = (concentration_in * v_in) + (current_conc * (Current_vol - v_out))/ v_new

    return concentration_new


def find_ratio_time(change_needed, source_concentration, current_conc, fractions_pos):
    t = 1
    rate_in_out = 0.010
    concentration_in = source_concentration
    Current_vol = 10.0
    
    if change_needed < 0.0: #lowering concentration
        change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)
        if change < change_needed:
            while change < change_needed:
                concentration_in = concentration_in +  fractions_pos
                change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)
            con_in = concentration_in - fractions_pos
            while change >= change_needed:
                t = t + 1
                change = (concentration_change_funtion(t, rate_in_out, con_in, current_conc, Current_vol) - current_conc)
        while change >= change_needed:
            t = t + 1
            change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)

    else:   #raising concentration
        change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)
        if change < change_needed:
            while change < change_needed:
                t = t + 1
                change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)
        elif change > change_needed:
            while change > change_needed:
                concentration_in = concentration_in - (source_concentration * fractions_pos)
                change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)
            while change < change_needed:
                t = t + 1
                change = (concentration_change_funtion(t, rate_in_out, concentration_in, current_conc, Current_vol) - current_conc)


    return t, concentration_in

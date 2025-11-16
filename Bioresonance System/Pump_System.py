########
#
#       This is a UI
#
#
#
from tkinter import*
from tkinter.font import nametofont
from tkinter import ttk
from ttkthemes import ThemedTk
from Bookshelf import*
import ttkbootstrap as ttkb
from ttkbootstrap.constants import*
from Pump_functions import*



import os
import math
import time
from datetime import datetime, timedelta
import serial
from serial.tools import list_ports
import threading
from queue import Queue

main = ttkb.Window(themename="superhero")
main.title("Pump system")

default_font = nametofont("TkDefaultFont")
default_font.configure(size=10)



#    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   Ths code assumes one pump per com port (not daisy chained like manufacturer suggests)
#
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

### Globals
global pump_list
pump_list = {}

global data_save_path
data_save_path = ""

#################
#       This is the global that need to be loaded and passd into the test thread
#       All testing function need this filled and checked
##########
global test_info
test_info = ({
    "Blank Pump": {},
    "Solution Pump": {},
    "In Pump": {},
    "Out Pump": {},
    "Min Mix Volume": float,
    "Max Mix Volume": float,
    "Max Rate mLpS": 0.010,
    "Dish Volume": float,
    "Cycles": float,
    "Period":float,
    "Experiment Duration": float,
    "Blank Start?": TRUE,
    "Full End?": TRUE,
    "S1 Concentration": float,
    "S2 Concentration": float,
    "Compensation Volume": float
})

#These queus are to allow the update treads to recieve
#and act on information from the testing threads
global data_q
data_q = Queue()
global settings_q
settings_q = Queue()
global status_q
status_q = Queue()


###########################################################################################

#     Test Functions

###########################################################################################

##########        
def Main_test_protocol(test_info, settings_q, status_q): #Main test program

    #load info
    cycle = int(test_info["Cycles"])
    Rate_in_out = test_info["Max Rate mLpS"]
    
    Solution_rates = 0.500

    #Generate pump settings
    Concentration_set =(
        {'C1':  test_info["S1 Concentration"],
         'C2':  test_info["S2 Concentration"],
         'Use Rate': Rate_in_out,
         'Vb mL': test_info['Min Mix Volume'] ,
         'Period (s)': test_info["Period"]
        }
    )

    #creates the volumes and times needed to make the concentration wave
    #from the lowest input concentration (C1) to the maximum (C2)
    order_set = wave_set_v3(Concentration_set)  

    #Set rates for experiment
    
    #In pump
    volume_limit_off = "VOL 0.00"   #make sure Vol limit is off for In/Out pumps
    resp = query(test_info["In Pump"], volume_limit_off)
    Max_rate_set_command = "RAT " + str(Rate_in_out)    #Set rate
    resp = query(test_info["In Pump"], Max_rate_set_command)
    #Out Pump with adjustment
    resp = query(test_info["Out Pump"], volume_limit_off)
    rate_set_command_adj = "RAT " + str(float(Rate_in_out)+ float(0.003))
    resp = query(test_info["Out Pump"], rate_set_command_adj)
    #Solutions 1 and 2
    rate_set_command_b = "RAT " + str(Solution_rates)
    rate_set_command_f = "RAT " + str(Solution_rates)
    resp = query(test_info["Blank Pump"], rate_set_command_b)
    resp = query(test_info["Solution Pump"], rate_set_command_f)
    
    time.sleep(0.1)
    
    #Start test
    resp = query(test_info["In Pump"], "RUN")
    resp = query(test_info["Out Pump"], "RUN")
    #Once In/Out pumps are running at the correct rate, they should not be altered

    
    #Check start Condition
    if test_info["Blank Start?"] == TRUE:
        First_hour = round(Rate_in_out * 3600.0, 2)
        resp = Pump_volume_Dual(test_info["Blank Pump"], First_hour, test_info["Solution Pump"], 0.00)
        update_text = "Running Blank, 1 hour\n" + Right_now()
        status_q.put(update_text)
        time.sleep(3600.0)

    #Add + 5mL every 12 hours check
    now = datetime.now()

    #The code runs for a number of cycles
    for c in range(cycle):  
        for s in range(len(order_set)): #each cycle it stepps through the times and concentrations to create the sin wave
            orders = order_set[s]
            #Mixing chamber Drift
            if s == (int(len(order_set)/2.0)): #middle of cycle (peak concentration)
                resp = Pump_volume_Dual(test_info["Blank Pump"], orders['v1'], test_info["Solution Pump"], float(orders['v2'] + test_info["Compensation Volume"]))
            elif s == int(len(order_set) - 1): #End of cycle (lowest concentration)
                resp = Pump_volume_Dual(test_info["Blank Pump"] , float(orders['v1']+ test_info["Compensation Volume"] ), test_info["Solution Pump"], orders['v2'])
            else:
                resp = Pump_volume_Dual(test_info["Blank Pump"], orders['v1'], test_info["Solution Pump"], orders['v2'])

            update_text =  (
                "Cycle " + str(int(c+1)) + " of " + str(cycle)  +
                "\nStep " + str(s) + " of " + str(len(order_set)) +
                "\nCurrent Concentration " + str(orders['Con']) + " mM"
                "\n" + Right_now()
                )
            status_q.put(update_text)
            time.sleep(orders['Time'])
        
    
    #Bring the dish to full cocentration over an hour at the end
    if test_info["Full End?"] == TRUE:
        last_hour = round(Rate_in_out * 3600.0, 2)
        resp = Pump_volume_Dual(test_info["Blank Pump"], 0.0, test_info["Solution Pump"], last_hour)
        update_text = "Going to max, 1 hour\n" + Right_now()
        status_q.put(update_text)
        time.sleep(3600.0)
    
    #End
    #Stops all Pumps
    update_text =  ("End at " + Right_now())
    status_q.put(update_text)
    
    #Stop pumps
    write(test_info["In Pump"],"STP")
    write(test_info["Out Pump"],"STP")
    
    #Clear buffers
    time.sleep(0.1)
    read(test_info["Blank Pump"])
    time.sleep(0.1)
    read(test_info["Solution Pump"])
    time.sleep(0.1)
    read(test_info["In Pump"])
    time.sleep(0.1)
    read(test_info["Out Pump"])

    main.update()




#########################################################################################################################################################
###########   Threading
#########################################################################################################################################################


def Run_custom_test_thread(): #function that triggers the Custom main test program to run in its own thread
    global test_info
    global settings_q
    global status_q

    test_program = threading.Thread(target=Main_test_protocol, args= (test_info, settings_q, status_q,))

    test_program.start()


def Run_Preset_test_thread(): #function that triggers the Custom main test program to run in its own thread
    global test_info
    global settings_q
    global status_q
    Load_presets()

    test_program = threading.Thread(target=Main_test_protocol, args= (test_info, settings_q, status_q,))

    test_program.start()

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################



########################
## Background Updaters
###########

def updater_data_status(data_q):    #background thread function to allow status updates on the main GUI
    while TRUE:
        if data_q.empty() == True:
            pass
        else:
            update_post = data_q.get()
            Status_update(L_Data_status, update_post)
            data_q.task_done()
        time.sleep(0.1)

def updater_test_set_stat(settings_q, status_q): #background thread function to allow status updates on the main GUI
    while TRUE:
        if settings_q.empty() == True:
            pass
        else:
            update_post = settings_q.get()
            Status_update(L_Test_settings, update_post)
            settings_q.task_done()

        if status_q.empty() == True:
            pass
        else:
            update_post = status_q.get()
            Status_update(L_Test_status, update_post)
            status_q.task_done()

        time.sleep(0.11)

############

#########################################################################################################################################################
def find_com_ports():   #returns any detected com ports
    check_for = serial.tools.list_ports.comports()
    ports = []
    for i in range(len(check_for)):
        ent = str(check_for[i])
        ent2 = ent.split('-')[0]
        ports.append(ent2.strip())
    return ports
   

def find_and_connect():         #attemps to connect to any open com port,
    global pump_list            #only returns connection to pumps
    portlist = find_com_ports()
    Status_update(L_Setup_indicator, portlist)
    pump_list = connect_com_ports(pump_list, portlist)
    input_list_dic_pumps(pump_list, l1)

def close_them():       #close the com comunications
    global pump_list
    didworked = close_all(pump_list)
    Status_update(L_Setup_indicator, didworked)

def read_from():    #Read the buffer of the selected pump
    global pump_list
    comport = get_sel_pump(l1)
    response = Debug_read(pump_list, comport)
    Status_update(L_Setup_indicator,response)

def write_from():   #Write to the selected pump manualy
    global pump_list
    comport = get_sel_pump(l1)
    command = str(E_RS232_Query.get())
    Debug_write(pump_list, comport, command)


def query_pump():           #Function Querys what is written in the list box
    global pump_list
    comport = get_sel_pump(l1)      #gets the connection from the selected pump
    command = (E_RS232_Query.get()).strip()  #formats the ipput
    response = Debug_query(pump_list, comport, command)   #performs the query
    post(T_RS232_Respond, response)             #posts response

def com_check_trig():
    global pump_list
    comport = get_sel_pump(l1)
    checked = com_check(pump_list,comport)
    post(T_RS232_Respond, checked)

def ver_check():
    global pump_list
    comport = get_sel_pump(l1)
    answer = Debug_query(pump_list, comport,'VER')
    post(T_RS232_Respond, answer)

###########################   Program Fuctions  #############################

def input_list_dic_pumps(pump_list, box):	    #resets listbox and loads the keys of a dictionary
    box.delete(0, END)                      # designed to handle the comport/adress system of the pumps
    for f, value in pump_list.items() :
        ent = f + ' (ad' + str(pump_list[f]['Pump Number']) + ')'
        box.insert(END, ent)

def get_sel_pump(Lbox):	#gets corrent entry from a listbox of pumps
    sel = Lbox.curselection()
    pump_ent = Lbox.get(sel)
    ent = (pump_ent.split(' ('))[0]
    return ent

def get_d_from_sel_pumps(Lbox, dic):	#returns the dictionary by key taken from a listbox entry
    sel = Lbox.curselection()           # designed to handle the comport/adress system of the pumps
    pump_ent = Lbox.get(sel)
    dic_key = (pump_ent.split(' ('))[0]
    dic_ent = dic[dic_key]
    return dic_ent, dic_key

def dup_check_pump(new):     #checks peripump hasnt been double assigned
    global test_info
    # dictionary of pump roles and connection to thier labels
    pumps = {"Blank Pump":label_blank,"Solution Pump": label_full,"In Pump":label_in,"Out Pump":label_out}
    for d, value in pumps.items():
        p_check = d
        if p_check == new:
            pass
        else:
            if test_info[p_check] == test_info[new]:
                pumps[d].config(background = "Red")


def set_blank_trig():       #assigns the com port for the blank pump
    global pump_list
    global test_info
    pumpname = l1.get(l1.curselection())
    label_blank.config(text= pumpname, background="limegreen", foreground="black")
    connection_dic, key = get_d_from_sel_pumps(l1, pump_list)
    test_info["Blank Pump"] = {key: connection_dic}
    dup_check_pump("Blank Pump")

def set_full_trig():       #assigns the com port for the solution pump
    global pump_list
    global test_info
    pumpname = l1.get(l1.curselection())
    label_full.config(text= pumpname, background="limegreen", foreground="black")
    connection_dic, key = get_d_from_sel_pumps(l1, pump_list)
    test_info["Solution Pump"] = {key: connection_dic}
    dup_check_pump("Solution Pump")

def set_in_trig():       #assigns the com port for the pump that feeds into the experiment dish
    global pump_list
    global test_info
    pumpname = l1.get(l1.curselection())
    label_in.config(text= pumpname, background="limegreen", foreground="black")
    connection_dic, key = get_d_from_sel_pumps(l1, pump_list)
    test_info["In Pump"] = {key: connection_dic}
    dup_check_pump("In Pump")

def set_out_trig():       #assigns the com port for the pump that drains the experiment dish
    global pump_list
    global test_info
    pumpname = l1.get(l1.curselection())
    label_out.config(text= pumpname, background="limegreen", foreground="black")
    connection_dic, key = get_d_from_sel_pumps(l1, pump_list)
    test_info["Out Pump"] = {key: connection_dic}
    dup_check_pump("Out Pump")



########################## Setup and Debug Protocols #####################


def load_custom_test_info(test_info):
    errors = []
    #Conatiner sizes
    num = mix_min_ent.get()
    if str(num).isdigit() == TRUE:
        test_info["Min Mix Volume"] = float(num)
    else:
        errors.append("Max Volume of Sample Dish not loaded")

    num = mix_max_ent.get()
    if str(num).isdigit() == TRUE:
        test_info["Max Mix Volume"] = float(num)
    else:
        errors.append("Max Volume of Mixing Chamber not loaded")

    num = exp_ent.get()  
    if str(num).isdigit() == TRUE:
        test_info["Dish Volume"] = float(num)
    else:
        errors.append("Experimnet volume not loaded")

    num = Rate_set_ent.get()
    try:
        test_info["Max Rate mLpS"] = float(num)
    except ValueError:
        errors.append("Max Rate not loaded")
        
    
    num = wave_ent.get()
    if str(num).isdigit() == TRUE:
        test_info["Period"] = float(num)
    else:
        errors.append("Experiment Period not loaded")

    num = E_s1_concentration.get()
    try:
        test_info["S1 Concentration"] = float(num)
    except ValueError:
        errors.append("Solution 1 Concentraton not loaded")

    num = E_s2_concentration.get()
    try:
        test_info["S2 Concentration"] = float(num)
    except ValueError:
        errors.append("Solution 2 Concentraton not loaded")
    

    if int(start_box.get()) == 0:
        test_info["Blank Start?"] = FALSE
    else:
        test_info["Blank Start?"] = TRUE

    if int(end_box.get()) == 0:
        test_info["Full End?"] = FALSE
    else:
        test_info["Full End?"] = TRUE

    #Experiment duration 
    num = dur_ent.get()
    if str(num).isdigit() == TRUE:
        test_info["Cycles"] = int(num)

        ex_time = float(num) * test_info["Period"]

        if test_info["Blank Start?"] == TRUE:
            ex_time = ex_time + 3600

        if test_info["Full End?"] == TRUE:
            ex_time = ex_time + 3600

        test_info["Experiment Duration"] = ex_time
    else:
        errors.append("Cycle Time not loaded")

    comp_vol = E_Custom_comp_vol.get()
    if str(comp_vol).isdigit() == TRUE:
        test_info["Compensation Volume"] = round(float(comp_vol), 3)
        E_Custom_comp_vol.delete()
        E_Custom_comp_vol.insert(0, str(round(float(comp_vol), 3)) )
    else:
        test_info["Compensation Volume"] = 0.000
        E_Custom_comp_vol.delete()
        E_Custom_comp_vol.insert(0, str(0.000))

    return errors




def test_set_verify(test_info): #Check for any errors in the test perameters
    errors = []
    #Pump assignment
    if (test_info["Blank Pump"] == {}  or
        test_info["Solution Pump"] == {} or
        test_info["In Pump"] == {} or
        test_info["Out Pump"]=={} ):         #check there are pumps assigned
        errors.append("All needed Pumps have not been assigned")
    else:   #check none of the pumps are double assigned
        if test_info["Blank Pump"] == test_info["Solution Pump"]:
            errors.append("Blank pump and Solution Pump Identical")
        if test_info["Blank Pump"] == test_info["In Pump"]:
            errors.append("Blank pump and In Pump Identical")
        if test_info["Blank Pump"] == test_info["Out Pump"]:
            errors.append("Blank pump and Out Pump Identical")
        if test_info["Solution Pump"] == test_info["In Pump"]:
            errors.append("Solution pump and In Pump Identical")
        if test_info["Solution Pump"] == test_info["Out Pump"]:
            errors.append("Solution pump and Out Pump Identical")
        if test_info["In Pump"] == test_info["Out Pump"]:
            errors.append("In Pump and Out Pump Identical")
    
    #Container details
    if test_info["Min Mix Volume"] >= test_info["Max Mix Volume"]:
        errors.append("Mixing Volume input unrealistic")

    change_vol = test_info["Max Rate mLpS"] * (test_info["Period"] / 21.0)
    available_vol = test_info["Max Mix Volume"] - test_info["Min Mix Volume"]
    if change_vol > available_vol:
        errors.append("Mixing Volume insufficiant")

    return errors
    

def test_info_verify(test_info):   #run checks to make sure the test settings are reasonable (eventualy)
    global settings_q
    
    test_time_min = int( (test_info["Experiment Duration"]/60.0) ) 
    test_time_hour = round(test_time_min/60.0, 2)
    text_message = (str(test_info["Max Rate mLpS"]) + " mL per Second\n" +
                    str(test_info["Period"]) + " Cycle time\n" +
                    str(test_time_min) + " Minutes, " +
                    str(test_time_hour) + " Hours Test Time\n" +
                    str(round( (test_info["Experiment Duration"] * test_info["Max Rate mLpS"])/1000.0, 2) ) +
                    " Liters of media Total"
    ) 
    settings_q.put(text_message)



def test_con_ver_trig():    #Do all of the pre-test checks
    global test_info
    global settings_q
    global status_q
    
    errors1 = load_custom_test_info(test_info)
    if errors1 == []:
        errors2 = test_set_verify(test_info)
        if errors2 == []:
            test_info_verify(test_info)
            status_q.put( str("Test Parameters Ready" ))
            B_Custom_Test["state"] = "normal"
        else:
            Status_update(L_Test_settings,errors2)
    else:
        Status_update(L_Test_settings,errors1)       

def full_test_trig():   #trigger the test to run after all the checks are done
    global test_info
    global settings_q
    global status_q
    errors1 = load_custom_test_info(test_info)
    if errors1 == []:
        errors2 = test_set_verify(test_info)
        if errors2 == []:
            test_info_verify(test_info)
            #if no errors are found, trigger the test thread to begin
            Run_custom_test_thread()  
        else:
            Status_update(L_Test_status, errors2)
    else:
        Status_update(L_Test_status, errors1)


######################   Presets


def Load_presets():
    global test_info
    global settings_q
    global status_q
    global data_q

    
    Frequ_numbers = ({  #Periods rounded to the nearest second
    "0.09 h": 40000, "0.18 h": 20000, "0.27 h":13333, "0.36 h":10000, 
    "0.45 h":8000, "0.54 h":6667, "0.63 h":5714, "0.72 h":5000, 
    "0.81 h":4444, "0.90 h":4000
    
    })

    Duration_seconds = ({   #Numbers are hours converted to seconds
        "24 Hours":86400 , "48 hour":172800,  "36 Hours":129600, "72 Hours":259200
    })

    Rate_numbers = ({"0.010 mL Per Sec": 0.010, "0.100 mL Per Sec":0.100, "1.000 mL Per Sec": 1.000 })
    
    in_rate = Rate_numbers[V_Rate_type.get()]

    test_info["Min Mix Volume"] = float(mix_min_ent.get())
    test_info["Max Mix Volume"] = float(mix_max_ent.get())
    test_info["Max Rate mLpS"] = in_rate
    test_info["Dish Volume"] = float(exp_ent.get())
    test_info["S1 Concentration"] = float(E_s1_concentration.get())
    test_info["S2 Concentration"] = float(E_s2_concentration.get())

    #Mixing Chamber compensator
    
    compensation_Volume = {"None":0.000, "+0.5mL ":0.500, "+1.0mL":1.000, "+1.5mL":1.500,"+2.0mL":2.000,"+2.5mL":2.500}
    test_info["Compensation Volume"]  = compensation_Volume[V_Comp_type.get()]

    #experiment time calculated in seconds 

    #Start and End conditions
    if V_Begin_type.get() == "1 hour Solution 1":
        start_blank = TRUE
    else:
        start_blank = FALSE

    if V_End_type.get() == "1 hour Solution 2":
        end_full = TRUE
    else:
        end_full = FALSE

    period = Frequ_numbers[V_freq_type.get()]

    if (start_blank == TRUE and end_full==FALSE)  or (start_blank==FALSE and end_full==TRUE):
        cycles = int( (Duration_seconds[V_Duration_type.get()] - 3600) // period)
        Experiment_time =  period * cycles + 3600
    elif (start_blank==TRUE) and (end_full== TRUE):
        cycles = int( ( (Duration_seconds[V_Duration_type.get()] - ( 2 * 3600) ) // period) )
        Experiment_time =  period * cycles + (2 * 3600)
    else:
        cycles = int( (Duration_seconds[V_Duration_type.get()]) // period)
        Experiment_time =  period * cycles


    test_info["Cycles"] = cycles
    test_info["Period"] = period
    test_info["Blank Start?"] = start_blank
    test_info["Full End?"] = end_full
    test_info["Experiment Duration"] = Experiment_time

    #Check set up
    setup_error = test_set_verify(test_info)
    if setup_error == []:
        test_info_verify(test_info)
        #allow Test start if everything passes
        B_Preset_test["state"] = "normal"
        

    else:
        Status_update(L_Test_status, setup_error)
        
    



######################   UI   ###############################
label_lbox = Label(main, text='Connected Pumps')
label_lbox.grid(row=3, column=2, padx=5, pady=5, sticky=NW)
l1 = Listbox(main, selectmode = SINGLE, exportselection=False, width = 15, height = 10)
l1.grid(row=4, column=2, columnspan=2, padx=5, pady=5, rowspan=6, sticky=NW)


############## Tab setup
nbtest = ttkb.Notebook(main)
nbtest.grid(row=2, column = 4, rowspan=10, columnspan= 6, padx=5, pady=5, sticky=N)
f1 = ttkb.Frame(nbtest)
nbtest.add(f1, text= "Setup")
f2 = ttkb.Frame(nbtest)
nbtest.add(f2, text= "Debug")
f3 = ttkb.Frame(nbtest)
nbtest.add(f3, text= "Testing")

###### Setup Tab

Connect_dev = Button(f1, text="Find Devices", command=find_and_connect)
Connect_dev.grid(row=1, column=1, columnspan=2, padx=5, pady=5, sticky=W)

set_blank = Button(f1, text="Set Blank Pump", command=set_blank_trig)
set_blank.grid(row=2, column=1, columnspan=2, padx=5, pady=5, sticky=W)
label_blank = Label(f1, text='Unconfirmed')
label_blank.grid(row=2, column=3, columnspan=2, padx=5, pady=5, sticky=NW)

set_full = Button(f1, text="Set Solution Pump", command=set_full_trig)
set_full.grid(row=3, column=1,columnspan=2, padx=5, pady=5, sticky=W)
label_full = Label(f1, text='Unconfirmed')
label_full.grid(row=3, column=3, columnspan=2, padx=5, pady=5, sticky=NW)

set_in = Button(f1, text="Set Sample In Pump", command=set_in_trig)
set_in.grid(row=4, column=1, columnspan=2, padx=5, pady=5, sticky=W)
label_in = Label(f1, text='Unconfirmed')
label_in.grid(row=4, column=3, columnspan=2, padx=5, pady=5, sticky=NW)

set_out = Button(f1, text="Set Sample Out Pump", command=set_out_trig)
set_out.grid(row=5, column=1, columnspan=2, padx=5, pady=5, sticky=W)
label_out = Label(f1, text='Unconfirmed')
label_out.grid(row=5, column=3, columnspan=2, padx=5, pady=5, sticky=NW)

label_dish_size = Label(f1, text='Min Volume of Mixing Chamber:')
label_dish_size.grid(row=6, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
mix_min_ent = Entry(f1, width = 6)
mix_min_ent.grid(row=6, column=3, padx=5, pady=5, sticky=NW)
label_dish_size_units = Label(f1, text='mL')
label_dish_size_units.grid(row=6, column=4, padx=5, pady=5, sticky=W)

label_Mix_size = Label(f1, text='Max Volume of Mixing Chamber:')
label_Mix_size.grid(row=7, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
mix_max_ent = Entry(f1, width = 6)
mix_max_ent.grid(row=7, column=3, padx=5, pady=5, sticky=NW)
label_mix_size_units = Label(f1, text='mL')
label_mix_size_units.grid(row=7, column=4, padx=5, pady=5, sticky=W)

label_exp_size = Label(f1, text='Experiment Volume:')
label_exp_size.grid(row=8, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
exp_ent = Entry(f1, width = 6)
exp_ent.grid(row=8, column=3, padx=5, pady=5, sticky=NW)
label_exp_size_units = Label(f1, text='mL')
label_exp_size_units.grid(row=8, column=4, padx=5, pady=5, sticky=W)


L_s1_concentration = Label(f1, text='Solution 1:')
L_s1_concentration.grid(row=9, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
E_s1_concentration = Entry(f1, width = 6)
E_s1_concentration.grid(row=9, column=3, padx=5, pady=5, sticky=NW)
L_s1_concentration_units = Label(f1, text='mM')
L_s1_concentration_units.grid(row=9, column=4, padx=5, pady=5, sticky=W)

L_s2_concentration = Label(f1, text='Solution 2:')
L_s2_concentration.grid(row=10, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
E_s2_concentration = Entry(f1, width = 6)
E_s2_concentration.grid(row=10, column=3, padx=5, pady=5, sticky=NW)
L_s2_concentration_units = Label(f1, text='mM')
L_s2_concentration_units.grid(row=10, column=4, padx=5, pady=5, sticky=W)

L_Setup_indicator = Label(f1)
L_Setup_indicator.grid(row=2, column=6, rowspan=10, columnspan=2)


###### Debug Tab
label_ep = Label(f2,text='Querys')
label_ep.grid(row=2, column=2, padx=5, pady=5, sticky=NW)
E_RS232_Query = Entry(f2, width = 80)
E_RS232_Query.grid(row=3, column=2, columnspan=10, padx=5, pady=5, sticky=SW)

version_check = Button(f2, text = "Version", command = ver_check)
version_check.grid(row=5, column = 2, padx=5, pady=5, sticky=W)

make_query = Button(f2, text = "Query", command = query_pump)
make_query.grid(row=4, column = 2, padx=5, pady=5, sticky=W)

just_write = Button(f2, text = "Write", command = write_from)
just_write.grid(row=4, column = 3, padx=5, pady=5, sticky=W)

just_read = Button(f2, text = "Read", command = read_from)
just_read.grid(row=4, column = 4, padx=5, pady=5, sticky=W)

Close_trig = Button(f2, text="Close", command=close_them)
Close_trig.grid(row=4, column=5, padx=5, pady=5, sticky=W)

Sc_RS232 = Scrollbar(f2, orient=VERTICAL)
Sc_RS232.grid(row=6, column=10, rowspan=10, sticky=NS) 

T_RS232_Respond = Text(f2, yscrollcommand= Sc_RS232.set, width = 60, height = 10)
T_RS232_Respond.grid(row=6, column = 1, padx=5, pady=5, rowspan=10, columnspan=9, sticky=N)


###### Test Tab

############ Test type
nbTest_type = ttkb.Notebook(f3)
nbTest_type.grid(row=2, column = 2, rowspan=5, columnspan= 5, padx=5, pady=5, sticky=N)
t1 = ttkb.Frame(nbtest)
nbTest_type.add(t1, text= "Preset Frequencys")
t2 = ttkb.Frame(nbtest)
nbTest_type.add(t2, text= "Custom")


L_freq_opt = Label(t1, text="Frequency:")
L_freq_opt.grid(row=1, column=1, padx=5, pady=5,sticky=NE)
freq_options = [
    "0.09 h", "0.18 h", "0.27 h", "0.36 h", "0.45 h", "0.54 h",
    "0.63 h", "0.72 h", "0.81 h", "0.90 h"
    ]
V_freq_type = StringVar()
V_freq_type.set(freq_options[0])
d_freq_type = OptionMenu(t1, V_freq_type, *freq_options)
d_freq_type.grid(row=1, column=2, padx=5, pady=5,sticky=NW)


L_Duration_opt = Label(t1, text="Estimated Duration:")
L_Duration_opt.grid(row=2, column=1, padx=5, pady=5,sticky=NE)
Duration_options = ["24 Hours", "48 hour",  "36 Hours", "72 Hours"]
V_Duration_type = StringVar()
V_Duration_type.set(Duration_options[1])
d_Duration_type = OptionMenu(t1, V_Duration_type, *Duration_options)
d_Duration_type.grid(row=2, column=2, padx=5, pady=5,sticky=NW)


L_Begin_opt = Label(t1, text="Start Options:")
L_Begin_opt.grid(row=3, column=1, padx=5, pady=5,sticky=NE)
Begin_options = ["1 hour Solution 1", "None"]
V_Begin_type = StringVar()
V_Begin_type.set(Begin_options[0])
d_Begin_type = OptionMenu(t1, V_Begin_type, *Begin_options)
d_Begin_type.grid(row=3, column=2, padx=5, pady=5,sticky=NW)


L_End_opt = Label(t1, text="Start Options:")
L_End_opt.grid(row=4, column=1, padx=5, pady=5,sticky=NE)
End_options = ["1 hour Solution 2", "None"]
V_End_type = StringVar()
V_End_type.set(End_options[0])
d_End_type = OptionMenu(t1, V_End_type, *End_options)
d_End_type.grid(row=4, column=2, padx=5, pady=5,sticky=NW)


L_Rate_opt = Label(t1, text="In Rate:")
L_Rate_opt.grid(row=5, column=1, padx=5, pady=5,sticky=NE)
Rate_options = ["0.010 mL Per Sec", "0.100 mL Per Sec", "1.000 mL Per Sec"]
V_Rate_type = StringVar()
V_Rate_type.set(Rate_options[0])
d_Rate_type = OptionMenu(t1, V_Rate_type, *Rate_options)
d_Rate_type.grid(row=5, column=2, padx=5, pady=5,sticky=NW)

L_Comp_opt = Label(t1, text="Drift Compensation\non Min and Max:")
L_Comp_opt.grid(row=6, column=1, padx=5, pady=5,sticky=NE)
Comp_options = ["None", "+0.5mL ","+1.0mL", "+1.5mL","+2.0mL","+2.5mL"]
V_Comp_type = StringVar()
V_Comp_type.set(Comp_options[0])
d_Comp_type = OptionMenu(t1, V_Comp_type, *Comp_options)
d_Comp_type.grid(row=6, column=2, padx=5, pady=5,sticky=NW)



B_Preset_verify = Button(t1, text = "Verify Test Conditions", font= 20, command = Load_presets)
B_Preset_verify.grid(row=2, column = 6, padx=5, pady=5, sticky=E)

B_Preset_test = Button(t1, text = "Run Experiment", font= 40, command = Run_Preset_test_thread)
B_Preset_test.grid(row=3, column = 6, rowspan=2, columnspan= 2, padx=20, pady=20, sticky=E)


######### Custom
wave_ent_label = Label(t2,text='Wavelength')
wave_ent_label.grid(row=2, column = 1, padx=5, pady=5, sticky=W)
wave_ent = Entry(t2, width=10)
wave_ent.grid(row=2, column = 2, padx=5, pady=5, sticky=W)
wave_ent_units = Label(t2,text='Seconds')
wave_ent_units.grid(row=2, column = 3, padx=5, pady=5, sticky=W)

dur_ent_label = Label(t2,text='Duration')
dur_ent_label.grid(row=3, column = 1, padx=5, pady=5, sticky=W)
dur_ent = Entry(t2, width=10)
dur_ent.grid(row=3, column = 2, padx=5, pady=5, sticky=W)
dur_ent_units = Label(t2,text='Cycles')
dur_ent_units.grid(row=3, column = 3, padx=5, pady=5, sticky=W)

label_Rate_set = Label(t2, text='Pump Rate:')
label_Rate_set.grid(row=5, column=1, columnspan=2, padx=5, pady=5, sticky=NW)
Rate_set_ent = Entry(t2, width = 6)
Rate_set_ent.grid(row=5, column=2, padx=5, pady=5, sticky=NW)
label_Rate_set_units = Label(t2, text='mL/S')
label_Rate_set_units.grid(row=5, column=3, padx=5, pady=5, sticky=W)

start_box = IntVar()
CB_Start_Blank = Checkbutton(t2, variable=start_box, text="1 hour start Solution 1")
CB_Start_Blank.grid(row=6, column=1,columnspan=3, padx=5, pady=5, sticky=NW)
CB_Start_Blank.select()

end_box = IntVar()
CB_End_full = Checkbutton(t2, variable=end_box ,text="1 hour end Solution 2")
CB_End_full.grid(row=7, column=1, columnspan=3, padx=5, pady=5, sticky=NW)
CB_End_full.select()

B_Custom_verify = Button(t2, text = "Verify Test Conditions", font= 20, command = test_con_ver_trig)
B_Custom_verify.grid(row=2, column = 6, padx=5, pady=5, sticky=E)

B_Custom_Test = Button(t2, text = "Run Experiment", font= 40, command = full_test_trig)
B_Custom_Test.grid(row=3, column = 6, rowspan=2, columnspan= 2, padx=20, pady=20, sticky=E)

L_Custom_comp_vol = Label(t2, text ="Compensate Volume:")
L_Custom_comp_vol.grid(row=8, column=1, columnspan=3, padx=5, pady=5, sticky=NW)
E_Custom_comp_vol = Entry(t2, width=6)
E_Custom_comp_vol.grid(row=9, column=1, padx=5, pady=5, sticky=NE)
L_Custom_comp_units = Label(t2, text ="mL at Peak and Valley")
L_Custom_comp_units.grid(row=9, column=2, columnspan=2, padx=5, pady=5, sticky=NW)


####Progress Boxes
LF_Progress_posts = LabelFrame(main, text="Progress", height=10, width=10)
LF_Progress_posts.grid(row=12, column=4, columnspan=7, padx=5, pady=5, sticky=SW)

L_Test_settings = Label(LF_Progress_posts, width=30, height=6)
L_Test_settings.grid(row=1, column=1, columnspan=2, padx=5, pady=5, sticky=NW)

L_Test_status = Label(LF_Progress_posts, width=30, height=6)
L_Test_status.grid(row=1, column=3, columnspan=2, padx=5, pady=5, sticky=NW)

L_Data_status = Label(LF_Progress_posts, width=30, height=6)
L_Data_status.grid(row=1, column=5, columnspan=2, padx=5, pady=5, sticky=NW)


########################################################################################################################################################################

#defaults

E_RS232_Query.insert(0,"")
#max volume of microscope size petri dish
mix_min_ent.insert(0,"5")
#max volume of what ever thing the solutions are being mixed in
mix_max_ent.insert(0,"100")
#volume that will be used in experiment chamber
exp_ent.insert(0,"3")
#default wavelegth in seconds
wave_ent.insert(0, "10000")
#default Number of cyles to run
dur_ent.insert(0, "16")
#default Pump Rate
Rate_set_ent.insert(0, "0.010")
E_Custom_comp_vol.insert(0, "0.500")
#Solution Concentrations
E_s1_concentration.insert(0, "0.0")
E_s2_concentration.insert(0, "50.0")
#Lock test Buttons behind verification
B_Preset_test["state"] = "disabled"
B_Custom_Test["state"] = "disabled"



#Status threads
d_update = threading.Thread(target=updater_data_status, args=(data_q,), daemon=True)
test_update = threading.Thread(target=updater_test_set_stat, args=(settings_q, status_q,), daemon=True)

d_update.start()
test_update.start()

mainloop()

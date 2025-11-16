from tkinter import*
import os
import math
import time
import openpyxl
import re
import json
import pyodbc
import datetime
import serial




def post(box, t):	#post a string, dictionary, or list at set box and inputs a \n
	if type(t) == str:
		box.insert(INSERT, str(t))
		box.insert(INSERT, '\n')
	elif type(t) == list:
		for n in range(len(t)):
			box.insert(INSERT,str(t[n]))
			box.insert(INSERT, '\n')
	elif type(t) == dict:
		for d, value in t.items():
			box.insert(INSERT, str(d) + ': ')
			box.insert(INSERT, '\n')
			box.insert(INSERT, str(t[d]))
			box.insert(INSERT, '\n')
	elif type(t) == int:
		box.insert(INSERT, 'Integer: ')
		box.insert(INSERT, str(t))
		box.insert(INSERT, '\n')
	elif type(t) == float:
		box.insert(INSERT, 'Float: ')
		box.insert(INSERT, str(t))
		box.insert(INSERT, '\n')
	elif type(t) == bytes:
		box.insert(INSERT, 'Bytes:')
		box.insert(INSERT, t)
		box.insert(INSERT, '\n')
		box.insert(INSERT, 'Decode:')
		box.insert(INSERT, t.decode('ascii'))
		box.insert(INSERT, '\n')
	elif type(t) == bytearray:
		box.insert(INSERT, 'Byte array:')
		box.insert(INSERT, t)
		box.insert(INSERT, '\n')
		box.insert(INSERT, 'Decode:')
		box.insert(INSERT, t.decode('ascii'))
		box.insert(INSERT, '\n')
	elif type(t) == tuple:
		box.insert(INSERT, 'Tuple:')
		box.insert(INSERT, t)
		box.insert(INSERT, '\n')
	elif type(t) == BooleanVar:
		box.insert(INSERT, "Boolean: ")
		if t == FALSE:
			box.insert(INSERT, "False")
		elif t == TRUE:
			box.insert(INSERT, "True")
		else:
			box.insert(INSERT, "Empty")
		box.insert(INSERT, '\n')
	else:
		box.insert(INSERT, str(type(t)))
		box.insert(INSERT, " : ")
		box.insert(INSERT, str(t))
		box.insert(INSERT, '\n')
	
def post_time(box): #function to post times to a text box to keep track ofhow long something is taking
	datetime_object = datetime.datetime.fromtimestamp(time.time())
	dt = datetime_object.strftime("%H:%M:%S")
	post(box,dt)

def Right_now(): #return a string of the curent time
	datetime_object = datetime.datetime.fromtimestamp(time.time())
	dt = datetime_object.strftime("%H:%M:%S")
	return dt

def get_sel(Lbox):	#gets the entry from a listbox set to single and returns it
	sel = Lbox.curselection()
	ents = Lbox.get(sel)
	return ents

def get_d_from_sel(Lbox, dic):	#returns the dictionary element by key taken from a listbox entry
	sel = Lbox.curselection()
	dic_key = Lbox.get(sel)
	dic_ent = dic[dic_key]
	return dic_ent


def input_list(ent, box):	#resets listbox and adds new list
	box.delete(0, END)
	for f in range(len(ent)):	
		box.insert(END, ent[f])
		

def input_list_dic(dic, box):	#resets listbox and loads the keys of a dictionary
	box.delete(0, END)
	for f, value in dic.items() :	
		box.insert(END, f)
				

def clear(box):
	box.delete(0.0, END)

def File_name_ext(fpath):		#returns name of file from a filepath
	sec = fpath.split('\\')
	fn = sec[-1]
	name, ext = os.path.splitext(fn)
	return name, ext


def build_filepath(p,nf):	#add to a filepath
	if type(nf) == str:
		np = p + '\\' + str(nf)
		if type(nf) == list:
			np = p
			for i in range(len(nf)):
				np += '\\' + str(nf[i])
		else:
			np = p
	return np
	
def pathcheck(thepath):
    if os.path.isdir(thepath):
        pass
    else:
        os.makedirs(thepath)

def get_sel(Lbox):	#gets the entry from a listbox set to single and returns them
	ents = ''
	sel = Lbox.curselection()
	if sel == ():
		ents = ''
	else:
		ents = Lbox.get(sel)
	return ents

def get_sel_all(Lbox):	#gets the entry from a listbox set to single and returns them
	ents = Lbox.get(0,END)
	return ents

def readablepaths(fpath):	#returns a list of filepaths from a folder given in an input path
	paths = []
	filelist = os.listdir(fpath)
	for f in range(len(filelist)):
		paths.append(fpath + "\\" + filelist[f])

	return paths

def nameit(path):	#return filename from path
	m = '?'
	if os.path.isdir(path) == True:
		n = os.path.basename(path)
		m = str(n)

	if os.path.isfile(path) ==  True:
		n = os.path.basename(path)
		a = n.split('.')
		if len(a) == 2 :
			m = a[0]
		else:
			for i in range(len(a)):	#return whole name except ext 
				if i is not (len(a) - 1):
					if i == 0:
						m = a[i]
					else:	#restores periods
						m += '.' + a[i]

	return m

def ext_str(path): #returns the extention (filetype as in file.png, this returns 'png') of the file from an absolute filepath
	a = path.split('\\')[-1]
	b = a.split('.')[-1]
	return b
	
def path_file_list(fpath): #A function to create a list of filenames and a 
	names = []				#dictionary to to link them to their files via filepath
	paths = {}
	try:
		filelist = os.listdir(fpath)
		for f in range(len(filelist)):
			now = fpath + "\\" + filelist[f]
			if os.path.isdir(now) == True:
				pass
			elif os.path.exists(now) == True:
				names.append(os.path.splitext(filelist[f])[0])
				paths.update({os.path.splitext(filelist[f])[0]:{'Filetype':os.path.splitext(filelist[f])[1], 'Path': now}})
			else:
				pass
	except NotADirectoryError:	#error handle for non-data folder objects
		pass
	return names, paths


def t_pathlist(fpath):	#return path list of all txt files in a walked path
	paths = []
	try:
		filelist = os.listdir(fpath)
		for f in range(len(filelist)):
			now = fpath + "\\" + filelist[f]
			if os.path.isdir(now) == True:
				rec = t_pathlist(now)
				for r in range(len(rec)):
					paths.append(rec[r])
			elif os.path.splitext(now)[1] == '.txt':
				paths.append(now)
			else:
				pass
	except NotADirectoryError:	#error handle for non-data folder objects
		pass
	return paths



def get_today():	#get current date and format it to needed string
	now = datetime.datetime.now()
	y = now.year
	m = str(now.month)
	d = str(now.day)
	
	final = str(y)
	
	if len(m) == 1:
		final += "-0"
		final += m
	else:
		final += "-"
		final += m
		
	if len(d) == 1:
		final += "-0"
		final += d
	else:
		final += "-"
		final += d
		
	return final


def get_PB_range(elements):	#function to return a noramlizes range for progressbar
	tot = elements
	n = tot/100
	return n


def Status_update(Label, status): #update the text of a tkInter Label
	if type(status) == str: #have the last few status shown w/ timestamps
		current_text = str(Label.cget('text'))
		new_text = ''
		if current_text == '':
			new_text = status + ' ' + Right_now() + '\n'
		else:
			update_lis = current_text.split('\n')
			if len(update_lis) < 5:	#five lines max
				new_text = current_text + status + ' ' + Right_now() + '\n'
			else:
				new_text = (
					update_lis[-4]  + '\n' +
					update_lis[-3]  + '\n' +
					update_lis[-2]  + '\n' +
					update_lis[-1]  + '\n' +
					status + ' ' + Right_now() 
				)

		Label.config(text = new_text)

	elif type(status) == list:
		upd = ""
		for i in range(len(status)):
			upd = upd + str(status[i]) + "\n"
		Label.config(text = upd)
	elif str(status).isdigit() == True:
		Label.config(text = str(status))
	elif type(status) == dict:
		upd = ""
		for d, value in status.items():
			upd = upd + str(d) + ': \n' + str(status[d]) + '\n'
		Label.config(text = upd)
	else:
		Label.config(text = str(type(status)))

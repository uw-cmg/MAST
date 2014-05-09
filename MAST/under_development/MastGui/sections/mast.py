##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk


class MastSection:
    def __init__(self, parent_window, grid_row):
        self.parent_window     = parent_window
        self.mast_window       = None
        self.window_properties = {\
                                     'system_name' : [None, None],\
                                 }
        self.summary_frame     = None
        self.has_summary       = False
        self.grid_row          = grid_row

    def content(self):
        content_start = "$mast\n"
        if self.window_properties['system_name'][0]:
            content_lines = "system_name %s\n" % self.window_properties['system_name'][0]
        else:
            content_lines = "system_name Unknown\n"
        content_end   = "$end\n\n"
 
        return content_start + content_lines + content_end
         
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="Mast")
        self.summary_frame.grid_rowconfigure(0, weight=1)
        self.summary_frame.grid_rowconfigure(1, weight=1)
        self.summary_frame.grid_columnconfigure(0, weight=1)
        self.summary_frame.grid_columnconfigure(1, weight=1)

        #display contents if present
        if self.has_summary:
            for p_name, p_value in self.window_properties.iteritems():
                label_key = Label(self.summary_frame, text=p_name.upper() + " : ")
                label_key.grid(row=0, column=0, sticky=E)
                label_value = Label(self.summary_frame, text=p_value[0])
                label_value.grid(row=0, column=1, sticky=W)

        #the frame contents
        button_text = "Add Mast properties"
        if self.has_summary:
            button_text = "Edit Mast properties"
        add_button = Button(self.summary_frame, text=button_text, command=self.create_window)
        add_button.grid(row=1, column=0, columnspan=2, padx=20, pady=10)
        
        self.summary_frame.grid(row=self.grid_row, column=0, sticky=W+E, padx=20, pady=5)

    def create_window(self):
        #property of new windows
        self.mast_window = ttk.Toplevel(self.parent_window)
        self.mast_window.title("Mast Properties")
        self.mast_window.minsize(200, 200)
        self.mast_window.grid_rowconfigure(0, weight=1)
        self.mast_window.grid_columnconfigure(0, weight=1)

        mast_frame = LabelFrame(self.mast_window, text="Mast Properties")
        mast_frame.grid(row=0, column=0, padx=20)
        mast_frame.grid_rowconfigure(0, weight=1)
        mast_frame.grid_rowconfigure(1, weight=1)
        mast_frame.grid_columnconfigure(0, weight=1)
        mast_frame.grid_rowconfigure(1, weight=1)

        for p_name, p_value in self.window_properties.iteritems():
            label = Label(mast_frame, text=p_name.upper())
            label.grid(row=0, column=0, sticky=E, pady=10)
            text_val = StringVar()
            if p_value[0] is not None:
                text_val.set(p_value[0])
            text  = Entry(mast_frame, textvariable=text_val)
            text.grid(row=0, column=1, sticky=W, pady=10)
            p_value[1] = text
        
        #save button
        save_button = Button(mast_frame, text="Save", command=self.print_value)
        save_button.grid(row=1, column=1, sticky=E, padx=20, pady=10)

    def print_value(self):
        for p_name, p_value in self.window_properties.iteritems():
            self.window_properties[p_name][0] = p_value[1].get()
        self.mast_window.destroy()
        self.mast_window = None
        self.has_summary = True
        self.summary_frame.grid_remove()
        self.summary()

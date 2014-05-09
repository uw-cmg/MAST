##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk


class PhononSection:
    def __init__(self, parent_window, grid_row):
        self.parent_window = parent_window
        self.grid_row = grid_row
        self.summary_frame = None
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="Phonon")

        #the frame contents
        add_button = Button(self.summary_frame, text="Add Phonon properties")
        add_button.pack(pady=20)

        self.summary_frame.grid(row=self.grid_row, column=0, sticky=W+E, padx=20, pady=5)

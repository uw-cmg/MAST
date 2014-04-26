##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk

from verticalScrolledFrame import VerticalScrolledFrame 
from sections.mast import MastSection 
from sections.structure import StructureSection
from sections.defects import DefectsSection
from sections.neb import NEBSection
from sections.ingredients import IngredientsSection
from sections.phonon import PhononSection

class MastGui:
    def __init__(self):
        self.root_window = Tk()
        self.vframe = VerticalScrolledFrame(self.root_window)
        self.vframe.pack(fill=BOTH, expand=1)
        self.window = self.vframe.interior
        self.sections = []
 
        self.initialize_window()

    def initialize_window(self):
        #window properties
        self.root_window.title("MAST")
        screensizeX = self.root_window.winfo_screenwidth()
        screensizeY = self.root_window.winfo_screenheight()
        self.root_window.minsize(screensizeX-100, screensizeY - 100)

        self.window.grid_rowconfigure(0, weight=1)
        self.window.grid_rowconfigure(1, weight=1)
        self.window.grid_rowconfigure(2, weight=1)
        self.window.grid_rowconfigure(3, weight=1)
        self.window.grid_rowconfigure(4, weight=1)
        self.window.grid_rowconfigure(5, weight=1)
        self.window.grid_rowconfigure(6, weight=1)
        self.window.grid_columnconfigure(0, weight=1)

        #add sections
        self.mast_section = MastSection(self.window, 0)
        self.mast_section.summary()
        self.sections.append(self.mast_section)

        self.structure_section = StructureSection(self.window, 1)
        self.structure_section.summary()
        self.sections.append(self.structure_section)

        self.defects_section = DefectsSection(self.window, 2, self.structure_section)
        self.defects_section.summary()
        self.sections.append(self.defects_section)

        self.neb_section = NEBSection(self.window, 3, self.structure_section, self.defects_section)
        self.neb_section.summary()
        self.sections.append(self.neb_section)

        self.ingredients_section = IngredientsSection(self.window, 4)
        self.ingredients_section.summary()
        self.sections.append(self.ingredients_section)

        self.phonon_section = PhononSection(self.window, 5)
        self.phonon_section.summary()

        self.create_btn = Button(self.window, text="Create Input File", command=self.create_input)
        self.create_btn.grid(row=6, column=0, pady=10, padx=10)


    def create_input(self):
        input_file = open("mast_input.inp", "w")

        for section_obj in self.sections:
            output_content =  section_obj.content()
            input_file.write(output_content)

        input_file.close()
        self.root_window.destroy()
       
        

    def run(self):
        self.root_window.mainloop()


if __name__ == "__main__":
    application_main = MastGui()
    application_main.run()

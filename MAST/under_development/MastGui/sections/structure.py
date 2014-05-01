##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk
from verticalScrolledFrame import VerticalScrolledFrame 


class StructureSection:
    def __init__(self, parent_window, grid_row):
        self.parent_window = parent_window
        self.grid_row = grid_row
        self.summary_frame = None
        self.structure_info = {\
                                  'coord_type' : [None, None],\
                                  'posfile'    : [None, None],\
                                  'elementmap' : [{}, {}],\
                                  'lattice'    : [None, None],\
                                  'coordinates': [{}, {}],\
                              }
        self.num_elt_fields = 1
        self.has_summary = False

    def get_elements_map(self):
        return self.structure_info['elementmap'][0]

    def content(self):
        content_start = "$structure\n"
        content_end   = "$end\n\n"

        content_lines = ""

        if self.structure_info['coord_type'][0]:
            content_lines += "coord_type %s\n\n" % self.structure_info['coord_type'][0]

        if self.structure_info['posfile'][0]:
            content_lines += "posfile %s\n\n" % self.structure_info['posfile'][0]

        if self.structure_info['elementmap'][0]:
            content_lines += "begin elementmap\n"
            for elt_key, elt_name in self.structure_info['elementmap'][0].iteritems():
                content_lines += "%s %s\n" % (elt_key, elt_name)
            content_lines += "end\n\n"

        if self.structure_info['coordinates'][0]:
            content_lines += "begin coordinates\n"
            for elt_key, elt_coord in self.structure_info['coordinates'][0].iteritems():
                coord_lines = elt_coord.split('\n')
                for line in coord_lines:
                    if line:
                        content_lines += "%s %s\n" % (elt_key, line)
            content_lines += "end\n\n"

        if self.structure_info['lattice'][0]:
            content_lines += "begin lattice\n"
            content_lines += "%s\n" % self.structure_info['lattice'][0] 
            content_lines += "end\n\n" 

        return content_start + content_lines + content_end
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="Structure")
        self.summary_frame.grid_rowconfigure(0, weight=1)
        self.summary_frame.grid_rowconfigure(1, weight=1)
        self.summary_frame.grid_columnconfigure(0, weight=1)
        self.summary_frame.grid_columnconfigure(1, weight=1)
        self.summary_frame.grid_columnconfigure(2, weight=1)

        #if it has summary show it
        if self.has_summary:
           #coord_type
           coord_type = Label(self.summary_frame, text="COORD_TYPE :")   
           coord_type.grid(row=0, column=1, sticky=E, pady=5)
           coord_type_value = Label(self.summary_frame, text=self.structure_info['coord_type'][0])
           coord_type_value.grid(row=0, column=2, sticky=W, pady=5)

           #posfile
           if self.structure_info['posfile'][0].strip() != "":
               posfile = Label(self.summary_frame, text="POSFILE :")   
               posfile.grid(row=1, column=1, sticky=E, pady=5)
               posfile_value = Label(self.summary_frame, text=self.structure_info['posfile'][0])
               posfile_value.grid(row=1, column=2, sticky=W, pady=5)

           #lattice structure
           lattice = Label(self.summary_frame, text="LATTICE : ")
           lattice.grid(row=2, column=1, sticky=E, pady=5)
           lattice_value = Label(self.summary_frame, text=self.structure_info['lattice'][0])
           lattice_value.grid(row=2, column=2, sticky=W, pady=5)

           for row_num in xrange(self.num_elt_fields + 1):
               if row_num == 0:
                   continue
               if self.structure_info['elementmap'][0]['X' + str(row_num)].strip() == "":
                   continue
               id_label = Label(self.summary_frame, text='X'+ str(row_num))
               id_label.grid(row=2 + row_num, column=0, padx=10, pady=2, sticky=E)
               element_label = Label(self.summary_frame, text=self.structure_info['elementmap'][0]['X'+str(row_num)])
               element_label.grid(row=2 + row_num, column=1, padx=10, pady=2, sticky=W)
               coordinates_label = Label(self.summary_frame, text=self.structure_info['coordinates'][0]['X'+str(row_num)])
               coordinates_label.grid(row=2 + row_num, column=2, padx=10, pady=2, sticky=W)
           


        #the frame contents
        add_button = Button(self.summary_frame, text="Add Structure properties", command=self.create_window)
        add_button.grid(row=3 + self.num_elt_fields, column=0, columnspan=3, pady=20, padx=20)

        self.summary_frame.grid(row=self.grid_row, column=0, padx=20, pady=5, sticky=W+E)

    def create_window(self):
        self.window = ttk.Toplevel(self.parent_window)
        self.window.title("Structure Information")
        self.window.minsize(500, 500)
        self.vframe = VerticalScrolledFrame(self.window)
        self.vframe.pack(fill=BOTH, expand=1)
        self.structure_window = self.vframe.interior
        self.structure_window.grid_rowconfigure(0, weight=1)
        self.structure_window.grid_columnconfigure(0, weight=1)
        self.structure_window.grid_rowconfigure(0, weight=1)
        self.structure_window.grid_columnconfigure(0, weight=1)

        structure_frame = LabelFrame(self.structure_window, text="Structure")
        structure_frame.grid(row=0, column=0, padx=20, pady=20)
        structure_frame.grid_rowconfigure(0, weight=1)
        structure_frame.grid_rowconfigure(1, weight=1)
        structure_frame.grid_columnconfigure(0, weight=1)
        structure_frame.grid_columnconfigure(1, weight=1)

        #add coordtype 
        label = Label(structure_frame, text="COORD_TYPE")
        label.grid(row=0, column=0, sticky=E, pady=10, padx=10)
        coord_type = StringVar(structure_frame)
        coord_type.set("fractional")
        if self.structure_info['coord_type'][0] is not None:
            coord_type.set(self.structure_info['coord_type'][0])
        coord_type_option = OptionMenu(structure_frame, coord_type, "fractional", "cartesian")  
        coord_type_option.grid(row=0, column=1, sticky=W, pady=10, padx=10)
        self.structure_info["coord_type"] = [coord_type.get(), coord_type]

        #pos file
        label = Label(structure_frame, text="POSFILE")
        label.grid(row=1, column=0, sticky=E, pady=10, padx=10)
        text_val = StringVar()
        if self.structure_info['posfile'][0] is not None:
            text_val.set(self.structure_info['posfile'][0])
        text  = Entry(structure_frame, textvariable=text_val)
        text.grid(row=1, column=1, sticky=W, pady=10, padx=10)
        self.structure_info['posfile'] = [text_val.get(), text_val]
        
      
        #pos file
        label = Label(structure_frame, text="LATTICE")
        label.grid(row=2, column=0, sticky=E, pady=10, padx=10)
        text  = Text(structure_frame, height=5, width=30)
        if self.structure_info['lattice'][0] is not None:
            text.insert('1.0', self.structure_info['lattice'][0])
        text.grid(row=2, column=1, sticky=W, pady=10, padx=10)
        self.structure_info['lattice'][1] = text

        #label frame to add elements info
        self.elements_frame =  LabelFrame(structure_frame, text="ELEMENTS INFORMATION")
        self.elements_frame.grid(row=3, column=0, columnspan=2, pady=10, padx=10)
        #create field headers 
        id_label = Label(self.elements_frame, text="ID")
        id_label.grid(row=0, column=0, padx=10, pady=10, sticky=W)
        element_label = Label(self.elements_frame, text="ELEMENT")
        element_label.grid(row=0, column=1, padx=10, pady=10, sticky=W)
        coordinates_label = Label(self.elements_frame, text="COORDINATES")
        coordinates_label.grid(row=0, column=2, padx=10, pady=10, sticky=W)

        for row_num in xrange(self.num_elt_fields + 1):
            if row_num == 0:
                continue
            self.create_element_fields(self.elements_frame, row_num)
        
        self.add_button = Button(self.elements_frame, text="Add more elements", command=self.add_elements)
        self.add_button.grid(row=self.num_elt_fields + 1, column=0, columnspan=3, pady=10, padx=10)

        #save button
        save_button = Button(structure_frame, text="Save", command=self.save)
        save_button.grid(row=4, column=0, columnspan=2, padx=20, pady=10)

    def save(self):
        self.structure_info['coord_type'][0] = self.structure_info['coord_type'][1].get()
        self.structure_info['posfile'][0] = self.structure_info['posfile'][1].get()
        self.structure_info['lattice'][0] = self.structure_info['lattice'][1].get("1.0", END)
        for row_num in xrange(self.num_elt_fields + 1):
            if row_num == 0:
                continue
            self.structure_info['elementmap'][0]['X'+str(row_num)] = self.structure_info['elementmap'][1]['X'+str(row_num)].get()
            self.structure_info['coordinates'][0]['X'+str(row_num)] = self.structure_info['coordinates'][1]['X'+str(row_num)].get("1.0", END)

        self.window.destroy()
        self.window = None
        self.has_summary = True
        self.summary_frame.grid_remove()
        self.summary()

    def add_elements(self):
        self.num_elt_fields += 1
        self.add_button.grid_remove()
        self.create_element_fields(self.elements_frame, self.num_elt_fields)
        self.add_button = Button(self.elements_frame, text="Add more elements", command=self.add_elements)
        self.add_button.grid(row=self.num_elt_fields + 1, column=0, columnspan=3, pady=10, padx=10)


    def create_element_fields(self, elements_frame, row_num):
        id_label = Label(elements_frame, text="X" + str(row_num))
        id_label.grid(row=row_num, column=0, padx=10, pady=10, sticky=W)
        text_val = StringVar()
        if self.structure_info['elementmap'][0] and ('X' + str(row_num)) in self.structure_info['elementmap'][0]:
            text_val.set(self.structure_info['elementmap'][0]['X'+str(row_num)])
        element_name  = Entry(elements_frame, textvariable=text_val)
        element_name.grid(row=row_num, column=1, pady=10, padx=10, sticky=W)
        self.structure_info['elementmap'][1]['X'+ str(row_num)] = text_val
        coordinates  = Text(elements_frame, height=5, width=30)
        if self.structure_info['coordinates'][0] and ('X' + str(row_num)) in self.structure_info['coordinates'][0]:
            coordinates.insert('1.0', self.structure_info['coordinates'][0]['X' + str(row_num)])
        coordinates.grid(row=row_num, column=2, pady=10, padx=10, sticky=W)
        self.structure_info['coordinates'][1]['X'+ str(row_num)] = coordinates
        
        
        

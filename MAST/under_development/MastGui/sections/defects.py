##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk
from verticalScrolledFrame import VerticalScrolledFrame 


class DefectsSection:
    def __init__(self, parent_window, grid_row, structure_section):
        self.parent_window = parent_window
        self.grid_row = grid_row
        self.summary_frame = None
        self.defects_info = {\
                                'coord_type'   : [None, None],\
                                'posfile'      : [None, None],\
                                'defecttype'   : [{}, {}],\
                                'namemap'      : [{}, {}],\
                                'elementmap'   : [{}, {}],\
                                'infomap'      : [{}, {}],\
                                'chargemap'    : [{}, {}],\
                            } 
        self.num_defect_fields = 1
        self.structure_section = structure_section
        self.has_summary=False
        self.defect_types = ['vacancy', 'interstitial', 'substitution', 'antisite']
        self.fields_hash  = {}

    def content(self):
        content_start = "$defects\n"
        content_end = "$end\n\n"

        content_lines = ""

        if self.defects_info['coord_type'][0]:
            content_lines += "coord_type %s\n\n" % self.defects_info['coord_type'][0]

        if self.defects_info['posfile'][0]:
            content_lines += "posfile %s\n\n" % self.defects_info['posfile'][0]

        if self.defects_info['namemap'][0]:
            for d_id, d_name in self.defects_info['namemap'][0].iteritems():
                content_lines += "begin %s\n" % d_name
                if self.defects_info['elementmap'][0][d_id]:
                    for idx, elt in enumerate(self.defects_info['elementmap'][0][d_id]):
                        elt = elt.split('(', 1)[1].replace(')', '').strip()
                        content_lines += "%s %s %s" % (self.defects_info['defecttype'][0][d_id][idx], self.defects_info['infomap'][0][d_id][idx], elt)
                        if self.defects_info['chargemap'][0][d_id][idx]:
                            content_lines += " charge=%s" % self.defects_info['chargemap'][0][d_id][idx]
                        content_lines += "\n"
                content_lines += "end\n\n"

        return content_start + content_lines + content_end

    def get_defects_map(self):
        return self.defects_info['namemap'][0]
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="Defects")
        self.summary_frame.grid_rowconfigure(0, weight=1)
        self.summary_frame.grid_rowconfigure(1, weight=1)
        self.summary_frame.grid_columnconfigure(0, weight=1)
        self.summary_frame.grid_columnconfigure(1, weight=1)
        self.summary_frame.grid_columnconfigure(2, weight=1)
        self.summary_frame.grid_columnconfigure(3, weight=1)
        self.summary_frame.grid_columnconfigure(4, weight=1)
        self.summary_frame.grid_columnconfigure(5, weight=1)

        #if it has summary show it
        if self.has_summary:
           #coord_type
           coord_type = Label(self.summary_frame, text="COORD_TYPE :")   
           coord_type.grid(row=0, column=2, sticky=E, pady=5)
           coord_type_value = Label(self.summary_frame, text=self.defects_info['coord_type'][0])
           coord_type_value.grid(row=0, column=3, sticky=W, pady=5)
           if self.defects_info['posfile'][0] is not None:
               pos_label = Label(self.summary_frame, text="POSFILE :")   
               pos_label.grid(row=1, column=2, sticky=E, pady=5)
               pos_value = Label(self.summary_frame, text=self.defects_info['posfile'][0])
               pos_value.grid(row=1, column=3, sticky=W, pady=5)

           for rownum in xrange(self.num_defect_fields):
               if rownum == 0:
                   continue
               if self.defects_info['elementmap'][0]['D' + str(rownum)][0].strip() == "":
                   continue
               id_label = Label(self.summary_frame, text='D'+ str(rownum))
               id_label.grid(row=2 + rownum, column=0, padx=10, pady=2, sticky=E)
               name_label = Label(self.summary_frame, text=self.defects_info['namemap'][0]['D'+str(rownum)])
               name_label.grid(row=2 + rownum, column=1, padx=10, pady=2, sticky=E)
               element_label = Label(self.summary_frame, text=self.defects_info['elementmap'][0]['D'+str(rownum)][0])
               element_label.grid(row=2 + rownum, column=2, padx=10, pady=2, sticky=W)
               type_label = Label(self.summary_frame, text=self.defects_info['defecttype'][0]['D'+str(rownum)][0])
               type_label.grid(row=2 + rownum, column=3, padx=10, pady=2, sticky=W)
               defect_info_label = Label(self.summary_frame, text=self.defects_info['infomap'][0]['D'+str(rownum)][0])
               defect_info_label.grid(row=2 + rownum, column=4, padx=10, pady=2, sticky=W)
               charge_label = Label(self.summary_frame, text=self.defects_info['chargemap'][0]['D'+str(rownum)][0])
               charge_label.grid(row=2 + rownum, column=5, padx=10, pady=2, sticky=W)

        #the frame contents
        add_button = Button(self.summary_frame, text="Add Defects", command=self.create_window)
        add_button.grid(row=3 + self.num_defect_fields, column=0, columnspan=6, pady=20, padx=20)

        self.summary_frame.grid(row=self.grid_row, column=0, padx=20, pady=5, sticky=W+E)


    def create_window(self):
        self.window = ttk.Toplevel(self.parent_window)
        self.window.title("Defects Information")
        self.window.minsize(500, 500)
        self.vframe = VerticalScrolledFrame(self.window)
        self.vframe.pack(fill=BOTH, expand=1)
        self.defects_window = self.vframe.interior
        self.defects_window.grid_rowconfigure(0, weight=1)
        self.defects_window.grid_columnconfigure(0, weight=1)


        defects_frame = LabelFrame(self.defects_window, text="Defects")
        defects_frame.grid(row=0, column=0, padx=20, pady=20)
        defects_frame.grid_rowconfigure(0, weight=1)
        defects_frame.grid_rowconfigure(1, weight=1)
        defects_frame.grid_columnconfigure(0, weight=1)
        defects_frame.grid_columnconfigure(1, weight=1)

        #add coordtype 
        label = Label(defects_frame, text="COORD_TYPE")
        label.grid(row=0, column=0, sticky=E, pady=10, padx=10)
        coord_type = StringVar(defects_frame)
        coord_type.set("fractional")
        if self.defects_info['coord_type'][0] is not None:
            coord_type.set(self.defects_info['coord_type'][0])
        coord_type_option = OptionMenu(defects_frame, coord_type, "fractional", "cartesian")  
        coord_type_option.grid(row=0, column=1, sticky=W, pady=10, padx=10)
        self.defects_info["coord_type"] = [coord_type.get(), coord_type]
        #add posfile
        label = Label(defects_frame, text="POS FILE")
        label.grid(row=1, column=0, sticky=E, pady=10, padx=10)
        pos_value = StringVar(defects_frame)
        pos_value.set("")
        if self.defects_info['posfile'][0] is not None:
            pos_value.set(self.defects_info['posfile'][0])
        pos_entry  = Entry(defects_frame, textvariable=pos_value)
        pos_entry.grid(row=1, column=1, sticky=W, pady=10, padx=10)
        self.defects_info["posfile"] = [pos_value.get(), pos_value]

        #a label frame to add defects related information
        self.defects_entities_frame = LabelFrame(defects_frame, text="DEFECT ENTITIES")  
        self.defects_entities_frame.grid(row=2, column=0, columnspan=2, padx=10, pady=10) 
        #create field headers 

        for rownum in xrange(self.num_defect_fields):
            self.create_defect_fields(self.defects_entities_frame, rownum)  

        self.add_button = Button(self.defects_entities_frame, text="Add more defects", command=self.add_elements)
        self.add_button.grid(row=self.num_defect_fields + 1, column=0, columnspan=2, pady=10, padx=10)

        #save button
        save_button = Button(defects_frame, text="Save", command=self.save)
        save_button.grid(row=3, column=0, columnspan=2, padx=20, pady=10)

    def save(self):
        self.defects_info['coord_type'][0] = self.defects_info['coord_type'][1].get()
        self.defects_info['posfile'][0]    = self.defects_info['posfile'][1].get() 
        for rownum in xrange(self.num_defect_fields):
            self.defects_info['namemap'][0]['D'+str(rownum)]    = self.defects_info['namemap'][1]['D'+str(rownum)].get()
            num_fields = self.fields_hash.get('D' + str(rownum), 0)
            for i in xrange(num_fields):
                if 'D' + str(rownum) in self.defects_info['elementmap'][0]:
                    if i < len(self.defects_info['elementmap'][0]['D' + str(rownum)]):
                        self.defects_info['elementmap'][0]['D' + str(rownum)][i] = self.defects_info['elementmap'][1]['D'+str(rownum)][i].get()
                    else:
                        self.defects_info['elementmap'][0]['D' + str(rownum)].append(self.defects_info['elementmap'][1]['D'+str(rownum)][i].get())
                else:
                    self.defects_info['elementmap'][0]['D'+str(rownum)] = [self.defects_info['elementmap'][1]['D'+str(rownum)][i].get()]
                if 'D' + str(rownum) in self.defects_info['infomap'][0]:
                    if i < len(self.defects_info['infomap'][0]['D' + str(rownum)]):
                        self.defects_info['infomap'][0]['D' + str(rownum)][i] = self.defects_info['infomap'][1]['D'+str(rownum)][i].get()
                    else:
                        self.defects_info['infomap'][0]['D' + str(rownum)].append(self.defects_info['infomap'][1]['D'+str(rownum)][i].get())
                else:
                    self.defects_info['infomap'][0]['D'+str(rownum)] = [self.defects_info['infomap'][1]['D'+str(rownum)][i].get()]
                if 'D' + str(rownum) in self.defects_info['chargemap'][0]:
                    if i < len(self.defects_info['chargemap'][0]['D' + str(rownum)]):
                        self.defects_info['chargemap'][0]['D' + str(rownum)][i] = self.defects_info['chargemap'][1]['D'+str(rownum)][i].get()
                    else:
                        self.defects_info['chargemap'][0]['D' + str(rownum)].append(self.defects_info['chargemap'][1]['D'+str(rownum)][i].get())
                else:
                    self.defects_info['chargemap'][0]['D'+str(rownum)] = [self.defects_info['chargemap'][1]['D'+str(rownum)][i].get()]
                if 'D' + str(rownum) in self.defects_info['defecttype'][0]:
                    if i < len(self.defects_info['defecttype'][0]['D' + str(rownum)]):
                        self.defects_info['defecttype'][0]['D' + str(rownum)][i] = self.defects_info['defecttype'][1]['D'+str(rownum)][i].get()
                    else:
                        self.defects_info['defecttype'][0]['D' + str(rownum)].append(self.defects_info['defecttype'][1]['D'+str(rownum)][i].get())
                else:
                    self.defects_info['defecttype'][0]['D'+str(rownum)] = [self.defects_info['defecttype'][1]['D'+str(rownum)][i].get()]

        self.window.destroy()
        self.window = None
        self.has_summary = True
        self.summary_frame.grid_remove()
        self.summary()
         

    def add_elements(self):
        self.num_defect_fields += 1
        self.add_button.grid_remove()
        self.create_defect_fields(self.defects_entities_frame, self.num_defect_fields - 1)
        self.add_button = Button(self.defects_entities_frame, text="Add more defects", command=self.add_elements)
        self.add_button.grid(row=self.num_defect_fields + 1, column=0, columnspan=6, pady=10, padx=10)
      

    def create_defect_fields(self, frame, rownum):
        elements_map = self.structure_section.get_elements_map()
        elements_options = ["%s (%s)" % (value, key) for key, value in elements_map.iteritems()]

        num_fields = self.fields_hash.get('D' + str(rownum), 1)
        self.fields_hash['D' + str(rownum)] = num_fields

        d_frame = LabelFrame(frame, text="D" + str(rownum))
        d_frame.grid(row=rownum, column=0, padx=10, pady=10, sticky=N+S+W+E)
        
        id_label = Label(d_frame, text="D" + str(rownum))
        id_label.grid(row=0, column=0, padx=10, pady=10, sticky=W)
        name_label = Label(d_frame, text="NAME")
        name_label.grid(row=0, column=0, padx=10, pady=10, sticky=W)
        element_label = Label(d_frame, text="ELEMENT")
        element_label.grid(row=0, column=1, padx=10, pady=10, sticky=W)
        vacancy_label = Label(d_frame, text="DEFECT_TYPE")
        vacancy_label.grid(row=0, column=2, padx=10, pady=10, sticky=W)
        vacancy_label = Label(d_frame, text="DEFECT_INFO")
        vacancy_label.grid(row=0, column=3, padx=10, pady=10, sticky=W)
        charge_label = Label(d_frame, text="CHARGE")
        charge_label.grid(row=0, column=4, padx=10, pady=10, sticky=W)

        name_text = StringVar()
        if self.defects_info['namemap'][0] and ('D' + str(rownum)) in self.defects_info['namemap'][0]:
            name_text.set(self.defects_info['namemap'][0]['D'+str(rownum)])
        name_val  = Entry(d_frame, textvariable=name_text)
        name_val.grid(row=1, column=0, pady=10, padx=10, sticky=W)
        self.defects_info['namemap'][1]['D' + str(rownum)] = name_text

        for i in xrange(num_fields):
            self.create_individual_fields(d_frame, elements_options, rownum, i + 1)

    def create_individual_fields(self, d_frame, elements_options, defect_id, rownum):
        element = StringVar(d_frame)
        element_val = None
        if self.defects_info['elementmap'][0] and ('D' + str(defect_id)) in self.defects_info['elementmap'][0] and \
                rownum <= len(self.defects_info['elementmap'][0]['D' + str(defect_id)]):
            element_val = self.defects_info['elementmap'][0]['D' + str(defect_id)][rownum - 1]
        option_menu = OptionMenu(d_frame, element, ())
        option_menu.grid(row=rownum, column=1, sticky=W, pady=10, padx=10)
        self.add_options(option_menu, elements_options, element, element_val)
        if 'D' + str(defect_id) in self.defects_info['elementmap'][1]:
            if rownum <= len(self.defects_info['elementmap'][1]['D' + str(defect_id)]):
                self.defects_info['elementmap'][1]['D' + str(defect_id)][rownum - 1] = element
            else:
                self.defects_info['elementmap'][1]['D' + str(defect_id)].append(element)
        else:
                self.defects_info['elementmap'][1]['D' + str(defect_id)] = [element]

        defect_type_text = StringVar(d_frame)
        defect_type_val  = None
        if self.defects_info['defecttype'][0] and ('D' + str(defect_id)) in self.defects_info['defecttype'][0] and \
                rownum <= len(self.defects_info['defecttype'][0]['D' + str(defect_id)]):
            defect_type_val = self.defects_info['defecttype'][0]['D' + str(defect_id)][rownum - 1] 
        defect_type_menu = OptionMenu(d_frame, defect_type_text, ())
        defect_type_menu.grid(row=rownum, column=2, sticky=W, pady=10, padx=10)
        self.add_options(defect_type_menu, self.defect_types, defect_type_text, defect_type_val)
        if 'D' + str(defect_id) in self.defects_info['defecttype'][1]:
            if rownum <= len(self.defects_info['defecttype'][1]['D' + str(defect_id)]):
                self.defects_info['defecttype'][1]['D' + str(defect_id)][rownum - 1] = defect_type_text
            else:
                self.defects_info['defecttype'][1]['D' + str(defect_id)].append(defect_type_text)
        else:
                self.defects_info['defecttype'][1]['D' + str(defect_id)] = [defect_type_text]

        vacancy_text = StringVar(d_frame)
        if self.defects_info['infomap'][0] and ('D' + str(defect_id)) in self.defects_info['infomap'][0] and \
                rownum <= len(self.defects_info['infomap'][0]['D' + str(defect_id)]):
            vacancy_text.set(self.defects_info['infomap'][0]['D'+str(defect_id)][rownum - 1])
        vacancy_val  = Entry(d_frame, textvariable=vacancy_text)
        vacancy_val.grid(row=rownum, column=3, pady=10, padx=10, sticky=W)
        if 'D' + str(defect_id) in self.defects_info['infomap'][1]:
            if rownum <= len(self.defects_info['infomap'][1]['D' + str(defect_id)]):
                self.defects_info['infomap'][1]['D' + str(defect_id)][rownum - 1] = vacancy_text
            else:
                self.defects_info['infomap'][1]['D' + str(defect_id)].append(vacancy_text)
        else:
                self.defects_info['infomap'][1]['D' + str(defect_id)] = [vacancy_text]

        charge_text = StringVar(d_frame)
        if self.defects_info['chargemap'][0] and ('D' + str(defect_id)) in self.defects_info['chargemap'][0] and \
                rownum <= len(self.defects_info['chargemap'][0]['D' + str(defect_id)]):
            charge_text.set(self.defects_info['chargemap'][0]['D'+str(defect_id)][rownum - 1])
        charge_val  = Entry(d_frame, textvariable=charge_text)
        charge_val.grid(row=rownum, column=4, pady=10, padx=10, sticky=W)
        if 'D' + str(defect_id) in self.defects_info['chargemap'][1]:
            if rownum <= len(self.defects_info['chargemap'][1]['D' + str(defect_id)]):
                self.defects_info['chargemap'][1]['D' + str(defect_id)][rownum - 1] = charge_text
            else:
                self.defects_info['chargemap'][1]['D' + str(defect_id)].append(charge_text)
        else:
                self.defects_info['chargemap'][1]['D' + str(defect_id)] = [charge_text]
        
        add_button = Button(d_frame, text="+")
        add_button.grid(row=rownum, column=5, pady=10, padx=10)
        add_button.configure(command=lambda:self.add_defect_field(add_button, d_frame, rownum, elements_options, defect_id))

    def add_defect_field(self, add_button, d_frame, rownum, elements_options, defect_id):
        add_button.grid_remove()
        num_fields = self.fields_hash.get('D' + str(defect_id), 1)
        self.fields_hash['D' + str(defect_id)] = num_fields + 1
        self.create_individual_fields(d_frame, elements_options, defect_id, rownum + 1)


    def add_options(self, option_menu, options, element_var, index_val):
        menu = option_menu["menu"]
        menu.delete(0, "end")
        for string in options:
            menu.add_command(label=string, command=lambda value=string:element_var.set(value))
        if index_val is not None:
            element_var.set(index_val)
        


        
 
        

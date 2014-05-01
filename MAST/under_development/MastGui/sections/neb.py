##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
import Tkinter as ttk
from verticalScrolledFrame import VerticalScrolledFrame 


class NEBSection:
    def __init__(self, parent_window, grid_row, structure_section, defects_section):
        self.parent_window = parent_window
        self.grid_row = grid_row
        self.structure_section = structure_section
        self.defects_section = defects_section
        self.has_summary = False
        self.neb_info = {\
                            'images'       : [None, None],\
                            'posfiles'     : [{}, {}],\
                            'startdefect'  : [{}, {}],\
                            'enddefect'    : [{}, {}],\
                            'elementmap'   : [{}, {}],\
                            'startpos'     : [{}, {}],\
                            'endpos'       : [{}, {}],\
                        }
        self.num_neb_fields = 1
        self.fields_hash    = {}

    def content(self):
        content_start = "$neb\n"
        content_end   = "$end\n\n"
        content_lines = ""

        if self.neb_info['startdefect'][0]:
            for n_id, s_defect in self.neb_info['startdefect'][0].iteritems():
                if n_id not in self.neb_info['enddefect'][0]:
                    continue
                content_lines += "begin %s-%s\n" % (s_defect.split('(', 1)[0].strip(), self.neb_info['enddefect'][0][n_id].split('(', 1)[0].strip())
                if n_id in self.neb_info['elementmap'][0]:
                    for idx, elt in enumerate(self.neb_info['elementmap'][0][n_id]):
                        elt = elt.split('(', 1)[1].replace(')', '').strip()
                        content_lines += "%s, %s, %s\n" % (elt, self.neb_info['startpos'][0][n_id][idx], self.neb_info['endpos'][0][n_id][idx])
                if n_id in self.neb_info['posfiles'][0]:
                    content_lines += "posfile %s\n" % ",".join(self.neb_info['posfiles'][0][n_id])
                    
                content_lines += "end\n\n"

        if self.neb_info['images'][0]:
            content_lines += "images %s\n" % self.neb_info['images'][0]

        return content_start + content_lines + content_end
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="NEB")
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
           images = Label(self.summary_frame, text="IMAGES :")   
           images.grid(row=0, column=2, sticky=E, pady=5)
           images_value = Label(self.summary_frame, text=self.neb_info['images'][0])
           images_value.grid(row=0, column=3, sticky=W, pady=5)

           for rownum in xrange(self.num_neb_fields):
               if self.neb_info['startdefect'][0]['N' + str(rownum)].strip() == "":
                   continue
               id_label = Label(self.summary_frame, text='N'+ str(rownum))
               id_label.grid(row=2 + rownum, column=0, padx=10, pady=2, sticky=E)
               s_defect_label = Label(self.summary_frame, text=self.neb_info['startdefect'][0]['N'+str(rownum)])
               s_defect_label.grid(row=2 + rownum, column=1, padx=10, pady=2, sticky=E)
               e_defect_label = Label(self.summary_frame, text=self.neb_info['enddefect'][0]['N'+str(rownum)])
               e_defect_label.grid(row=2 + rownum, column=2, padx=10, pady=2, sticky=E)
               element_label = Label(self.summary_frame, text=self.neb_info['elementmap'][0]['N'+str(rownum)][0])
               element_label.grid(row=2 + rownum, column=3, padx=10, pady=2, sticky=E)
               s_pos_label = Label(self.summary_frame, text=self.neb_info['startpos'][0]['N'+str(rownum)][0])
               s_pos_label.grid(row=2 + rownum, column=4, padx=10, pady=2, sticky=E)
               e_pos_label = Label(self.summary_frame, text=self.neb_info['endpos'][0]['N'+str(rownum)][0])
               e_pos_label.grid(row=2 + rownum, column=5, padx=10, pady=2, sticky=E)

        #the frame contents
        add_button = Button(self.summary_frame, text="Add NEB properties", command=self.create_window)
        add_button.grid(row=3 + self.num_neb_fields, column=0, columnspan=7, pady=20, padx=20)

        self.summary_frame.grid(row=self.grid_row, column=0, sticky=W+E, padx=20, pady=5)

    def create_window_with_posfiles(self):
        self.neb_info['images'][0] = self.neb_info['images'][1].get()
        self.window.destroy()
        self.window = None
        self.create_window()

    def create_window(self):
        self.window = ttk.Toplevel(self.parent_window)
        self.window.title("NEB Information")
        self.window.minsize(500, 500)
        self.vframe = VerticalScrolledFrame(self.window)
        self.vframe.pack(fill=BOTH, expand=1)
        self.neb_window = self.vframe.interior
        self.neb_window.grid_rowconfigure(0, weight=1)
        self.neb_window.grid_columnconfigure(0, weight=1)


        neb_frame = LabelFrame(self.neb_window, text="NEB")
        neb_frame.grid(row=0, column=0, padx=20, pady=20)
        neb_frame.grid_rowconfigure(0, weight=1)
        neb_frame.grid_columnconfigure(0, weight=1)

        #frame for images and button
        frame1 = Frame(neb_frame)
        frame1.grid(row=0, column=0)
        frame1.grid_rowconfigure(0, weight=1)
        frame1.grid_rowconfigure(1, weight=1)
        frame1.grid_columnconfigure(0, weight=1)
        frame1.grid_columnconfigure(1, weight=1)
        frame1.grid_columnconfigure(2, weight=1)

        #add images 
        label = Label(frame1, text="IMAGES")
        label.grid(row=0, column=0, sticky=E, pady=10, padx=10)
        images_val = StringVar(frame1)
        if self.neb_info['images'][0] is not None:
            images_val.set(self.neb_info['images'][0])
        images_menu  = Entry(frame1, textvariable=images_val)
        images_menu.grid(row=0, column=1, sticky=W, pady=10, padx=10)
        self.neb_info["images"] = [images_val.get(), images_val]
        images_btn = Button(frame1, text="Add Posfiles", command=self.create_window_with_posfiles)
        images_btn.grid(row=0, column=2, sticky=W, pady=10, padx=10)

        #frame for entities
        frame2 = Frame(neb_frame)
        frame2.grid(row=1, column=0)



        #a label frame to add neb related information
        self.neb_entities_frame = LabelFrame(frame2, text="NEB ENTITIES")  
        self.neb_entities_frame.grid(row=0, column=0, columnspan=3, padx=10, pady=10) 
        
         
        for rownum in xrange(self.num_neb_fields):
            self.create_neb_fields(self.neb_entities_frame, rownum)  

        self.add_button = Button(self.neb_entities_frame, text="Add more neb", command=self.add_elements)
        self.add_button.grid(row=self.num_neb_fields + 1, column=0, columnspan=6, pady=10, padx=10)

        #save button
        save_button = Button(neb_frame, text="Save", command=self.save)
        save_button.grid(row=3, column=0, columnspan=2, padx=20, pady=10)

    def add_elements(self):
        self.num_neb_fields += 1
        self.add_button.grid_remove()
        self.create_neb_fields(self.neb_entities_frame, self.num_neb_fields - 1)
        self.add_button = Button(self.neb_entities_frame, text="Add more neb", command=self.add_elements)
        self.add_button.grid(row=self.num_neb_fields + 1, column=0, pady=10, padx=10)

    def save(self):
        self.neb_info['images'][0] = self.neb_info['images'][1].get()
        for rownum in xrange(self.num_neb_fields):
            self.neb_info['startdefect'][0]['N'+str(rownum)]    = self.neb_info['startdefect'][1]['N'+str(rownum)].get()
            self.neb_info['enddefect'][0]['N'+str(rownum)]      = self.neb_info['enddefect'][1]['N'+str(rownum)].get()
            num_fields = self.fields_hash.get('N' + str(rownum), 0)
            for i in xrange(num_fields):
                if 'N' + str(rownum) in self.neb_info['elementmap'][0]:
                    if i < len(self.neb_info['elementmap'][0]['N' + str(rownum)]):
                        self.neb_info['elementmap'][0]['N' + str(rownum)][i] = self.neb_info['elementmap'][1]['N'+str(rownum)][i].get()
                    else:
                        self.neb_info['elementmap'][0]['N' + str(rownum)].append(self.neb_info['elementmap'][1]['N'+str(rownum)][i].get())
                else:
                    self.neb_info['elementmap'][0]['N'+str(rownum)] = [self.neb_info['elementmap'][1]['N'+str(rownum)][i].get()]
                if 'N' + str(rownum) in self.neb_info['startpos'][0]:
                    if i < len(self.neb_info['startpos'][0]['N' + str(rownum)]):
                        self.neb_info['startpos'][0]['N' + str(rownum)][i] = self.neb_info['startpos'][1]['N'+str(rownum)][i].get()
                    else:
                        self.neb_info['startpos'][0]['N' + str(rownum)].append(self.neb_info['startpos'][1]['N'+str(rownum)][i].get())
                else:
                    self.neb_info['startpos'][0]['N'+str(rownum)] = [self.neb_info['startpos'][1]['N'+str(rownum)][i].get()]
                if 'N' + str(rownum) in self.neb_info['endpos'][0]:
                    if i < len(self.neb_info['endpos'][0]['N' + str(rownum)]):
                        self.neb_info['endpos'][0]['N' + str(rownum)][i] = self.neb_info['endpos'][1]['N'+str(rownum)][i].get()
                    else:
                        self.neb_info['endpos'][0]['N' + str(rownum)].append(self.neb_info['endpos'][1]['N'+str(rownum)][i].get())
                else:
                    self.neb_info['endpos'][0]['N'+str(rownum)] = [self.neb_info['endpos'][1]['N'+str(rownum)][i].get()]

            #save the posfile values
            if not self.neb_info['posfiles'][1]:
                continue
            for entry in self.neb_info['posfiles'][1]['N' + str(rownum)]:
                val_list = self.neb_info['posfiles'][0].setdefault('N' + str(rownum), [])
                val_list.append(entry.get())

        self.window.destroy()
        self.window = None
        self.has_summary = True
        self.summary_frame.grid_remove()
        self.summary()

        


    def create_neb_fields(self, frame, rownum):
        elements_map = self.structure_section.get_elements_map()
        elements_options = ["%s (%s)" % (value, key) for key, value in elements_map.iteritems()]
        defects_map = self.defects_section.get_defects_map()
        defects_options = ["%s (%s)" % (value, key) for key, value in defects_map.iteritems()]

        num_fields = self.fields_hash.get('N' + str(rownum), 1)
        self.fields_hash['N' + str(rownum)] = num_fields


        #label frame for this entity
        m_frame = LabelFrame(frame, text="N" + str(rownum))
        m_frame.grid(row=rownum, column=0, padx=10, pady=10)

        e_frame = Frame(m_frame)
        e_frame.grid(row=0, column=0)
        p_frame = Frame(m_frame)
        p_frame.grid(row=1, column=0)

        #create field headers 
        start_label = Label(e_frame, text="START DEFECT")
        start_label.grid(row=0, column=0, padx=10, pady=10, sticky=W)
        end_label = Label(e_frame, text="END DEFECT")
        end_label.grid(row=0, column=1, padx=10, pady=10, sticky=W)
        element_label = Label(e_frame, text="ELEMENT")
        element_label.grid(row=0, column=2, padx=10, pady=10, sticky=W)
        start_pos = Label(e_frame, text="START POSITION")
        start_pos.grid(row=0, column=3, padx=10, pady=10, sticky=W)
        end_pos = Label(e_frame, text="END POSITION")
        end_pos.grid(row=0, column=4, padx=10, pady=10, sticky=W)



        start_defect = StringVar(e_frame)
        start_defect_val = None
        if self.neb_info['startdefect'][0] and ('N' + str(rownum)) in self.neb_info['startdefect'][0]:
            start_defect_val = self.neb_info['startdefect'][0]['N' + str(rownum)]
        s_option_menu = OptionMenu(e_frame, start_defect, ())
        s_option_menu.grid(row=1, column=0, sticky=W, pady=10, padx=10)
        self.add_options(s_option_menu, defects_options, start_defect, start_defect_val)
        self.neb_info['startdefect'][1]['N' + str(rownum)] = start_defect

        end_defect = StringVar(e_frame)
        end_defect_val = None
        if self.neb_info['enddefect'][0] and ('N' + str(rownum)) in self.neb_info['enddefect'][0]:
            end_defect_val = self.neb_info['enddefect'][0]['N' + str(rownum)]
        e_option_menu = OptionMenu(e_frame, end_defect, ())
        e_option_menu.grid(row=1, column=1, sticky=W, pady=10, padx=10)
        self.add_options(e_option_menu, defects_options, end_defect, end_defect_val)
        self.neb_info['enddefect'][1]['N' + str(rownum)] = end_defect

        for i in xrange(num_fields):
            self.create_individual_fields(e_frame, elements_options, rownum, i + 1)

        #put the pos file input boxes
        #add posfiles if needed
        if self.neb_info['images'][0] is not None and self.neb_info['images'][0] != "" and int(self.neb_info['images'][0]) > 0:
            label = Label(p_frame, text="POSFILES")
            label.grid(row=0, column=0, sticky=E, pady=10, padx=10)
            pos_files_frame = Frame(p_frame)
            pos_files_frame.grid(row=0, column=1, sticky=W, pady=10, padx=10)
            #create boxes
            for i in xrange(int(self.neb_info['images'][0])):
                entry_val = StringVar(pos_files_frame)
                if 'N' + str(rownum) in self.neb_info['posfiles'][0] and len(self.neb_info['posfiles'][0]['N' + str(rownum)]) >= i + 1:
                    entry_val.set(self.neb_info['posfiles'][0]['N' + str(rownum)][i])
                    self.neb_info['posfiles'][1]['N' + str(rownum)][i] = entry_val
                else:
                    if not self.neb_info['posfiles'][1] or 'N' + str(rownum) not in self.neb_info['posfiles'][1]:
                        self.neb_info['posfiles'][1]['N' + str(rownum)] = [entry_val]
                    else:
                        self.neb_info['posfiles'][1]['N' + str(rownum)].append(entry_val)
                pos_entry = Entry(pos_files_frame, textvariable=entry_val)
                pos_entry.grid(row=i, column=0, sticky=W)


    def create_individual_fields(self, e_frame, elements_options, neb_id, rownum):
        element = StringVar(e_frame)
        element_val = None
        if self.neb_info['elementmap'][0] and ('N' + str(neb_id)) in self.neb_info['elementmap'][0] and \
                rownum <= len(self.neb_info['elementmap'][0]['N' + str(neb_id)]):
            element_val = self.neb_info['elementmap'][0]['N' + str(neb_id)][rownum - 1]
        option_menu = OptionMenu(e_frame, element, ())
        option_menu.grid(row=rownum, column=2, sticky=W, pady=10, padx=10)
        self.add_options(option_menu, elements_options, element, element_val)
        if 'N' + str(neb_id) in self.neb_info['elementmap'][1]:
            if rownum <= len(self.neb_info['elementmap'][1]['N' + str(neb_id)]):
                self.neb_info['elementmap'][1]['N' + str(neb_id)][rownum - 1] = element
            else:
                self.neb_info['elementmap'][1]['N' + str(neb_id)].append(element)
        else:
                self.neb_info['elementmap'][1]['N' + str(neb_id)] = [element]
         
        startpos_text = StringVar()
        if self.neb_info['startpos'][0] and ('N' + str(neb_id)) in self.neb_info['startpos'][0] and \
                rownum <= len(self.neb_info['startpos'][0]['N' + str(neb_id)]):
            startpos_text.set(self.neb_info['startpos'][0]['N'+str(neb_id)][rownum - 1])
        startpos_val  = Entry(e_frame, textvariable=startpos_text)
        startpos_val.grid(row=rownum, column=3, pady=10, padx=10, sticky=W)
        if 'N' + str(neb_id) in self.neb_info['startpos'][1]:
            if rownum <= len(self.neb_info['startpos'][1]['N' + str(neb_id)]):
                self.neb_info['startpos'][1]['N' + str(neb_id)][rownum - 1] = startpos_text
            else:
                self.neb_info['startpos'][1]['N' + str(neb_id)].append(startpos_text)
        else:
                self.neb_info['startpos'][1]['N' + str(neb_id)] = [startpos_text]

        endpos_text = StringVar()
        if self.neb_info['endpos'][0] and ('N' + str(rownum)) in self.neb_info['endpos'][0] and \
                rownum <= len(self.neb_info['endpos'][0]['N' + str(neb_id)]):
            endpos_text.set(self.neb_info['endpos'][0]['N'+str(neb_id)][rownum - 1])
        endpos_val  = Entry(e_frame, textvariable=endpos_text)
        endpos_val.grid(row=rownum, column=4, pady=10, padx=10, sticky=W)
        if 'N' + str(neb_id) in self.neb_info['endpos'][1]:
            if rownum <= len(self.neb_info['endpos'][1]['N' + str(neb_id)]):
                self.neb_info['endpos'][1]['N' + str(neb_id)][rownum - 1] = endpos_text
            else:
                self.neb_info['endpos'][1]['N' + str(neb_id)].append(endpos_text)
        else:
                self.neb_info['endpos'][1]['N' + str(neb_id)] = [endpos_text]

        add_button = Button(e_frame, text="+")
        add_button.grid(row=rownum, column=5, pady=10, padx=10)
        add_button.configure(command=lambda:self.add_neb_field(add_button, e_frame, rownum, elements_options, neb_id))

    def add_neb_field(self, add_button, e_frame, rownum, elements_options, neb_id):
        add_button.grid_remove()
        num_fields = self.fields_hash.get('N' + str(neb_id), 1)
        self.fields_hash['N' + str(neb_id)] = num_fields + 1
        self.create_individual_fields(e_frame, elements_options, neb_id, rownum + 1)

    def add_options(self, option_menu, options, element_var, index_val):
        menu = option_menu["menu"]
        menu.delete(0, "end")
        for string in options:
            menu.add_command(label=string, command=lambda value=string:element_var.set(value))
        if index_val is not None:
            element_var.set(index_val)

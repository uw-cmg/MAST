##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Kumaresh Visakan Murugan
# Last updated: 2014-01-01
##############################################################
from Tkinter import *
from ttk import Combobox
from verticalScrolledFrame import VerticalScrolledFrame 


class IngredientsSection:
    def __init__(self, parent_window, grid_row):
        self.parent_window = parent_window
        self.grid_row = grid_row
        self.summary_frame = None
        self.keywords = [\
                           ('mast_program', 'vasp'), ('mast_kpoints', '2x2x2 M'), ('mast_pp_setup', 'La=La Mn=Mn_pv O=O_s'), ('mast_xc', 'PW91'), ('mast_multiplyencut', '1.5'),\
                           ('mast_setmagmom', '1 -1 1 -1 1 -1 1 -1'), ('mast_charge', '1'),('mast_coordinates', 'im1poscar,im2poscar,im3poscar'), ('mast_strain', '1.01 1.03 0.98'),\
                           ('mast_write_method', 'no_setup'), ('mast_ready_method', 'ready_singlerun'), ('mast_run_method', 'run_singlerun'), ('mast_complete_method', 'complete_singlerun'),\
                           ('mast_update_children_method', 'give_structure'), ('mast_nodes', '3'), ('mast_ppn', '8'), ('mast_queue', 'morgan1'),\
                           ('mast_exec', '//opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_CNEB'), ('mast_walltime', '12'), ('mast_memory', '2'),\
                           ('ibrion', '-1'), ('lcharge', 'False'), ('lwave', 'False'),\
                        ]
        self.ingredients_info = {\
                                    'global'       : [None, None],\
                                    'ingredients'  : [{}, {}],\
                                }
        self.key_boxes = {}
        self.num_ing_fields = 1

    def content(self):
        content_start = "$ingredients\n"
        content_end   = "$end\n\n"
        content_lines = ""

        if self.ingredients_info['global'][0]: 
            content_lines += "begin ingredients_global\n"
            content_lines += "%s\n" % self.ingredients_info['global'][0]
            content_lines += "end\n\n"

        if self.ingredients_info['ingredients'][0]:
            for ing_key, ing_data in self.ingredients_info['ingredients'][0].iteritems():
                content_lines += "begin %s\n" % ing_data[0]
                content_lines += "%s\n" % ing_data[1]
                content_lines += "end\n\n"

        return content_start + content_lines + content_end
    
    def summary(self):
        #the frame which encloses all the contents
        self.summary_frame = LabelFrame(self.parent_window, text="Ingredients")

        #the frame contents
        add_button = Button(self.summary_frame, text="Add Ingredients", command=self.create_window)
        add_button.pack(pady=20)

        self.summary_frame.grid(row=self.grid_row, column=0, sticky=W+E, padx=20, pady=5)

    def create_window(self):
        self.window = Toplevel(self.parent_window)
        self.window.title("Ingredients Information")
        self.window.minsize(500, 600)
        self.vframe = VerticalScrolledFrame(self.window)
        self.vframe.pack(fill=BOTH, expand=1)
        self.ingredients_window = self.vframe.interior
        self.ingredients_window.grid_rowconfigure(0, weight=1)
        self.ingredients_window.grid_rowconfigure(1, weight=1)
        self.ingredients_window.grid_rowconfigure(2, weight=1)
        self.ingredients_window.grid_columnconfigure(0, weight=1)


        #global param frame
        global_frame = LabelFrame(self.ingredients_window, text="Global Parameters")
        global_frame.grid(row=0, column=0, padx=20, pady=20, sticky=N+S+W+E)

        global_label = Label(global_frame, text="Global Parameters")
        global_label.grid(row=0, column=0, padx=20, sticky=W)
        global_box = Text(global_frame, height=5, width=80)
        if self.ingredients_info['global'][0] is not None:
            global_box.insert('1.0', self.ingredients_info['global'][0])
        else:
            global_box.insert('1.0', '\n'.join(["%s    %s" % (key, value) for key, value in self.keywords]))
        global_box.grid(row=1, column=0, padx=20, pady=5)
        self.ingredients_info['global'][1] = global_box 

        self.ingredients_frame = LabelFrame(self.ingredients_window, text="Ingredients information")
        self.ingredients_frame.grid(row=1, column=0, padx=20, pady=5, sticky=N+S+W+E)

        self.ingredients_frame.grid_rowconfigure(0, weight=1)
        self.ingredients_frame.grid_columnconfigure(0, weight=1)

        #ingredients frames
        for rownum in xrange(self.num_ing_fields):
            self.create_ingredient_fields(self.ingredients_frame, rownum) 

        self.add_button = Button(self.ingredients_frame, text="Add more Ingredients", command=self.add_ingredients)
        self.add_button.grid(row=self.num_ing_fields + 1, column=0, pady=10, padx=10)

        #save button
        save_button = Button(self.ingredients_window, text="Save", command=self.save)
        save_button.grid(row=2, column=0, padx=20, pady=10, sticky=N+S+W+E)
            

        #individual ingredients
        #param_box = self.create_param_box(global_frame, "global")
        #ind_frame = LabelFrame(self.ingredients_window, text="Ingredients")
        #ind_frame.grid(row=1, column=0, padx=20, pady=20)
        #for num in xrange(self.num_ing_fields):
        #    param_box = self.create_param_box(ind_frame, "individual", ingredient_index=num)
        #    param_box.grid(row=num, column=0, padx=20, pady=20)

    def save(self):
        self.ingredients_info['global'][0] = self.ingredients_info['global'][1].get("1.0", END)
        for rownum in xrange(self.num_ing_fields):
            self.ingredients_info['ingredients'][0][rownum] = [self.ingredients_info['ingredients'][1][rownum][0].get(), self.ingredients_info['ingredients'][1][rownum][1].get("1.0", END)]

        self.window.destroy()
        self.window = None
        self.has_summary = True
        self.summary_frame.grid_remove()
        self.summary()
       
    def add_ingredients(self):
        self.num_ing_fields += 1
        self.add_button.grid_remove()
        self.create_ingredient_fields(self.ingredients_frame, self.num_ing_fields-1)
        self.add_button = Button(self.ingredients_frame, text="Add more Ingredients", command=self.add_ingredients)
        self.add_button.grid(row=self.num_ing_fields + 1, column=0, pady=10, padx=10)

    def create_ingredient_fields(self, frame, rownum):
        if rownum not in self.ingredients_info['ingredients'][1]:
            self.ingredients_info['ingredients'][1][rownum] = [None, None]

        ing_frame = LabelFrame(frame, text="Ingredient %s" % (rownum + 1))
        ing_frame.grid(row=rownum, column=0, padx=10, pady=10, sticky=N+S+W+E)
        ing_frame.grid_rowconfigure(0, weight=1)
        ing_frame.grid_columnconfigure(0, weight=1)
        ing_frame.grid_columnconfigure(1, weight=10)

        space_label = Label(ing_frame, text="")
        space_label.grid(row=0, column=0, pady=5, sticky=W)

        name_label = Label(ing_frame, text="Ingredient Name")
        name_label.grid(row=1, column=0, columnspan=2, padx=10, pady=2, sticky=W)
        name_value = StringVar()
        if rownum in self.ingredients_info['ingredients'][0]:
            name_value.set(self.ingredients_info['ingredients'][0][rownum][0])
        name_box = Entry(ing_frame, textvariable=name_value)
        name_box.grid(row=2, column=0, columnspan=2, padx=10, pady=2, sticky=W)

        param_label = Label(ing_frame, text="Ingredient Parameters")
        param_label.grid(row=3, column=0, columnspan=2, padx=10, pady=2, sticky=W)

        key_option_val  = StringVar(ing_frame)
        key_option_val.set(self.keywords[0][0])
        key_option_menu = Combobox(ing_frame, textvariable=key_option_val, state='readonly')
        key_option_menu['values'] = [key for key, value in self.keywords]
        key_option_menu.grid(row=4, column=0, sticky=W, pady=10, padx=10)
        add_key_btn = Button(ing_frame, text="+")
        add_key_btn.grid(row=4, column=1, padx=10, pady=10, sticky=W)
        self.key_boxes[rownum] = key_option_val

        param_box = Text(ing_frame, height=5, width=80)
        if rownum in self.ingredients_info['ingredients'][0]:
            param_box.insert('1.0', self.ingredients_info['ingredients'][0][rownum][1])
        param_box.grid(row=5, column=0, columnspan=2, padx=10, pady=2, sticky=W)


        space_label = Label(ing_frame, text="")
        space_label.grid(row=5, column=0, pady=5, sticky=W)
        add_key_btn.configure(command=lambda:self.add_param_text(param_box, key_option_val))
        self.ingredients_info['ingredients'][1][rownum] = [name_value, param_box]


    def add_param_text(self, param_box, option_val):
        text = param_box.get("1.0", END)
        text += option_val.get()
        param_box.delete("1.0", END)
        param_box.insert("1.0", text)


    def create_param_box(self, parent, box_type, ingredient_index=-1): 
        values_map = {}
        object_map = {}

        #global param
        if box_type == 'global':
           values_map = self.ingredients_info['global'][0]
           object_map = self.ingredients_info['global'][1]
        #individual ingredients
        elif box_type == 'individual':
           if ingredient_index >= 0 and len(self.ingredients_info['ingredients'][0]) >= ingredient_index + 1:
               values_map = self.ingredients_info['ingredients'][0][ingredient_index]
               object_map = self.ingredients_info['ingredients'][1][ingredient_index]

        frame = Frame(parent)

        self.add_param_fields(frame, 0, object_map)
        add_button = Button(frame, text="+")
        add_button.configure(command=lambda:self.add_param(add_button, parent, box_type))
        add_button.grid(row=0, column=2, pady=20, padx=20)

        return frame

    def add_param(self, add_button, parent, box_type):
        print "success"
        add_button.grid_remove()
        self.create_param_box(parent, box_type)
       
       
        
    def add_param_fields(self, frame, rownum, object_map):
        #create the fields
        key_option_val  = StringVar(frame)
        key_option_val.set(self.keywords[0])
        key_option_menu = Combobox(frame, textvariable=key_option_val, state='readonly')
        key_option_menu['values'] = self.keywords
        key_option_menu.grid(row=0, column=0, sticky=W, pady=10, padx=10)

        val_text   = StringVar(frame)
        val_entry  = Entry(frame, textvariable=val_text)
        val_entry.grid(row=0, column=1, pady=10, padx=10, sticky=W)
        
        object_map.append((key_option_val, val_text))


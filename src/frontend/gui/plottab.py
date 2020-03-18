import tkinter as tk
import tkinter.ttk as ttk
import PIL.Image
import PIL.ImageTk

from frontend.gui.eventhandler import EventHandler
from frontend.gui.widgets import *  # TODO Remove wildcard import
from frontend.mod.constants import *
from frontend.mod.filehandler import FileHandler
from frontend.mod.plotter import plot_graph
from frontend.mod.session import Session, Property, PlotProperty


class PlotTab(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.parent = parent
        self._init_widgets()

    def _init_widgets(self):

        Session.write("Initializing the Plot Tab")
        self.display = PlotDisplay(self)
        self.display.pack(expand=True, fill=tk.BOTH, side=tk.LEFT, anchor=tk.E)
        self.sidebar = PlotSidebar(self)
        self.sidebar.pack(expand=True, fill=tk.Y, side=tk.TOP, anchor=tk.E)

    def create_tab(self, property_name):

        use_cur_prop = False
        
        # Check current property for fielname, plotname, and source
        if (self.sidebar.output_path.get().strip() != "") or \
            (self.sidebar.plot_name.get().strip() != "") or \
            (self.sidebar.input_path.get().strip() != ""):

            self.sidebar.create_new()

        self.sidebar.output_path.set(property_name)
        self.sidebar.plot_name.set(property_name + " Plot")
        self.sidebar.input_path.set(property_name)
        self.sidebar.load_file('../work-dir/output/' + property_name)
 

class PlotDisplay(tk.Canvas):

    def __init__(self, parent, **options):
        tk.Canvas.__init__(self, parent, borderwidth=0, background="#ffffff",
            **options)

        self.frame = tk.Frame(self, **options)
        self.frame.config(bg='white')
        self.scrollbar = tk.Scrollbar(parent, orient=tk.VERTICAL, 
            command=self.yview)
        self.configure(yscrollcommand=self.scrollbar.set)

        self.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scrollbar.pack(side=tk.LEFT, fill=tk.Y)

        self.create_window((0, 0), window=self.frame, anchor=tk.NW, 
            tags='self.frame')

        self.frame.bind('<Configure>', self._on_frame_configure)

        self.parent = parent
        self.result_list = {}
        self.image_list = {}
        self.large_image_list = {}
        self.image_counter = 0

    def add_image(self, filename, identifier):
        self._add_edit_image(filename, identifier, True)

    def edit_image(self, filename, identifier):
        self._add_edit_image(filename, identifier, False)
    
    def upload_new_image(self, filename, identifier):

        if not identifier in self.image_list:
            self.add_image(filename, identifier)
        else:
            self.edit_image(filename, identifier)

    def destory_image(self, filename, identifier):

        for widget in self.result_list[identifier].winfo_children():
            widget.destory()
        self.result_list[identifier].destroy()
        del self.result_list[identifier]
        del self.image_list[identifier]
        self.image_counter -= 1
        self._re_display_image()

    def destroy_all(self):
        """ Destroy all images """

        for key, value in self.result_list.items():
            value.grid_forget()
        
        self.result_list = {}
        self.image_list = {}
        self.large_image_list = {}
        self.image_counter = 0
        self._re_display_image()

    def _on_frame_configure(self, event):
        self.configure(scrollregion=self.bbox('all'))

    def _callback_btn_press_double(self, filename):

        import webbrowser
        webbrowser.open('../work-dir/plotresult/'+filename+'.html')

    def _add_edit_image(self, filename, identifier, action):
        
        path = '../work-dir/plotresult/' + filename + '.png'

        if platform == 'linux':
            small_image = PIL.Image.open(path)
        else:
            small_image = PIL.Image.open(path).convert('RGB')

        small_image = small_image.resize((200, 160), PIL.Image.ANTIALIAS)
        small_image = PIL.ImageTk.PhotoImage(small_image)
        self.image_list[identifier] = small_image

        # Same as previous similar conditional - unsure, but leaving in for now
        if platform == 'linux':
            display_image = PIL.Image.open(path)
        else:
            display_image = PIL.Image.open(path).convert('RGB')

        display_image = PIL.ImageTk.PhotoImage(display_image)
        self.large_image_list[filename] = display_image

        if action:  # Add Image
            self.result_list[identifier] = tk.LabelFrame(self.frame, 
                text=filename)
            self.result_list[identifier].config(bg='white', bd=5)
            self.result_list[identifier].grid(row=int(self.image_counter/3), 
                column=self.image_counter%3)
        else:  # Edit Image
            for widget in self.result_list[identifier].winfo_children():
                widget.destroy()
            self.result_list[identifier]['text'] = filename

        image = tk.Label(self.result_list[identifier], image=small_image)
        image.pack()
        image.config(bg='white')
        image.bind('<Double-Button-1>', 
            lambda event, x=filename: self._callback_btn_press_double(x))
        
        self.image_counter += 1

    def _re_display_image(self):

        counter = 0
        for key in sorted(self.result_list):
            self.result_list[key].grid(row=int(counter/3), column=counter%3)
            counter += 1


class PlotSidebar(tk.Frame):

    def __init__(self, parent, **options):
        tk.Frame.__init__(self, parent, **options)

        self.parent = parent

        # Selection Variables
        self.cur_prop = None  # PlotProperty()
        self.next_identifier = 0
        self.plot_prop_list = []  #  [PlotProperty()]
        self.sel_iid = None  # Int
        self.iid_list = []  # Int

        # Primary Widgets
        self.plot_view = None
        self.prop_view = None
        self.plot_list = None

        # Plot Variables
        self.input_path = tk.StringVar()
        self.output_path = tk.StringVar()
        self.plot_name = tk.StringVar()
        self.plot_status = tk.StringVar()
        self.variable_list = []  # tk.StringVar()
        self.mode_list = []  # tk.StringVar()
        self.horizontal_index = tk.IntVar()
        self.vertical_select = []  # tk.IntVar()
        self.vertical_indices = []  # Int
        
        # Plot button enable/disable conditions
        self.box_disable = True  
        self.plot_name_disable = True
        self.output_path_disable = True
        
        # Property Variables
        self.prop_name = tk.StringVar()
        self.prop_status = tk.StringVar()
        self.time_step = tk.DoubleVar()
        self.time_horizon = tk.DoubleVar()
        self.k_value = tk.DoubleVar()
        self.simulator = tk.StringVar()
        self.refine_strat = tk.StringVar()
        self.initial_set = tk.StringVar()
        self.unsafe_set = tk.StringVar()

        self._init_widgets()
        self._callback_new()

        self.OPEN_OPT = {
            'multiple': 0,
            'initialdir': '../work-dir/output/',
            'parent': self.parent,
            'title': 'Open File',
        }

    def create_new(self):
        """ Keep _callback_new private. """

        self._callback_new()

    def load_file(self, input_path_full):

        self.cur_prop.input_path_full = input_path_full
        self.cur_prop.input_path = input_path_full.split("/").pop()
        self.input_path.set(self.cur_prop.input_path)
        data_dict = FileHandler.parse_data_file(input_path_full)
               
        # Load Remaining property variables and set the appriopriate tk vars
        self.cur_prop.prop_name = data_dict[PROPNAME]
        self.cur_prop.prop_status = data_dict[SIMVER]
        self.cur_prop.time_step = data_dict[TIMESTEP]
        self.cur_prop.time_horizon = data_dict[TIMEHORIZON]
        self.cur_prop.k_value = data_dict[KVALUE]
        self.cur_prop.simulator = data_dict[SIMULATOR]
        self.cur_prop.refine_strat = data_dict[REFINEMENTSTRAT]
        self.cur_prop.initial_set = data_dict[INITIALSET]
        self.cur_prop.unsafe_set = data_dict[UNSAFESET]

        self.cur_prop.variable_list = data_dict[VARIABLENAMES]
        self.cur_prop.mode_list = data_dict[MODENAMES]

        self._load_cur_prop()
        self._create_boxes()

    def _init_widgets(self):

        self._init_plot_view()
        self.plot_view.pack(fill=tk.X)
        self.plot_view.columnconfigure(1, weight=1)

        self._init_prop_view()
        self.prop_view.pack(fill=tk.X)
        self.prop_view.columnconfigure(1, weight=1)

        self._init_plot_list()
        self.plot_list.pack(expand=True, fill=tk.BOTH)

    def _init_plot_view(self):

        # Plot View Frame
        self.plot_view = tk.LabelFrame(self, text="Plot")

        # Output File
        tk.Label(self.plot_view, text="Filename:")\
            .grid(row=0, column=0, sticky=tk.W)
        self.output_path.trace_variable('w', self._callback_output_path)
        tk.Entry(self.plot_view, textvariable=self.output_path)\
            .grid(row=0, column=1, sticky=tk.EW)
        self.output_path_vl = ValidLabel(self.plot_view)
        self.output_path_vl.grid(row=0, column=2, sticky=tk.E)

        # Plot Name
        tk.Label(self.plot_view, text="Plot Name:")\
            .grid(row=1, column=0, sticky=tk.W)
        self.plot_name.trace_variable('w', self._callback_plot_name)
        tk.Entry(self.plot_view, textvariable=self.plot_name)\
            .grid(row=1, column=1, sticky=tk.EW)
        self.plot_name_vl = ValidLabel(self.plot_view)
        self.plot_name_vl.grid(row=1, column=2, sticky=tk.E)

        # Input File
        tk.Label(self.plot_view, text="Source:")\
            .grid(row=2, column=0, sticky=tk.W)
        self.input_path.trace_variable('w', self._callback_input_path)
        tk.Entry(self.plot_view, textvariable=self.input_path, 
            state=tk.DISABLED)\
            .grid(row=2, column=1, sticky=tk.EW)
        tk.Button(self.plot_view, text="...", command=self._callback_open)\
            .grid(row=2, column=2, sticky=tk.EW)

        # Horizontal Axis
        self.horizontal_axis_frame = tk.LabelFrame(self.plot_view, 
            text="Horizontal Axis",height=100)
        self.horizontal_axis_frame.grid(row=3, column=0, columnspan=3,
             sticky=tk.EW)    
        self.horizontal_axis_frame.grid_propagate(0)

        self.h_canvas = tk.Canvas(self.horizontal_axis_frame, borderwidth=0,height = 100)
        # self.h_canvas.grid_propagate(0)
        self.h_frame = tk.Frame(self.h_canvas)
        # self.h_frame.grid_propagate(0)
        self.h_vsb = tk.Scrollbar(self.horizontal_axis_frame, orient="vertical", command=self.h_canvas.yview)
        self.h_canvas.configure(yscrollcommand=self.h_vsb.set)

        self.h_vsb.pack(side="right", fill="y")
        self.h_canvas.pack(side="left", fill="both", expand=True)
        self.h_canvas.create_window((4,4), window=self.h_frame, anchor="nw", 
                                  tags="self.h_frame")
        # self.h_frame.bind("<Configure>", self.onFrameConfigure)
        
        # Vertical Axis
        self.vertical_axis_frame = tk.LabelFrame(self.plot_view,
            text="Vertical Axis",height=100)
        self.vertical_axis_frame.grid(row=4, column=0, columnspan=3, 
            sticky=tk.EW)
        self.vertical_axis_frame.grid_propagate(0)

        self.v_canvas = tk.Canvas(self.vertical_axis_frame, borderwidth=0,height = 100)
        # self.v_canvas.grid_propagate(0)
        self.v_frame = tk.Frame(self.v_canvas)
        # self.v_frame.grid_propagate(0)
        self.v_vsb = tk.Scrollbar(self.vertical_axis_frame, orient="vertical", command=self.v_canvas.yview)
        self.v_canvas.configure(yscrollcommand=self.v_vsb.set)

        self.v_frame.bind("<Button-4>", self.v_mouseup)
        self.v_vsb.pack(side="right", fill="y")
        self.v_canvas.pack(side="left", fill="both", expand=True)
        self.v_canvas.create_window((4,4), window=self.v_frame, anchor="nw", 
                                  tags="self.v_frame")
        # self.v_frame.bind("<Configure>", self.onFrameConfigure)

        # self.vertical_axis_canvas = Canvas(self.vertical_axis_frame)
        # self.vertical_canvas_frame = Frame(self.vertical_axis_canvas)
        # self.vertical_axis_scrollbar=Scrollbar(self.vertical_axis_frame,orient="vertical",command=self.vertical_axis_canvas.yview)
        # self.vertical_axis_canvas.create_window((0,0),window=self.vertical_canvas_frame,anchor='nw')
        # self.vertical_axis_canvas.configure(scrollregion=self.vertical_axis_canvas.bbox('all'), yscrollcommand=self.vertical_axis_scrollbar.set)
        # self.vertical_axis_canvas.pack(fill='both', expand=True, side='left')
        # self.vertical_axis_scrollbar.pack(side="right",fill="y")
                
    # def onFrameConfigure(self, event):
    #     '''Reset the scroll region to encompass the inner frame'''
    #     self.canvas.configure(scrollregion=self.canvas.bbox("all"))


    def v_mouseup(self, event):
        self.v_frame.yview_scroll(-1*(event.delta/120), "units")

    def _init_prop_view(self):
        
        # Property View Frame
        self.prop_view = SimpleToggleFrame(self, text="Requirement")

        # Name
        tk.Label(self.prop_view.body_frame, text="Name:")\
            .grid(row=0, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.prop_name, 
            state=tk.DISABLED)\
            .grid(row=0, column=1, sticky=tk.EW)
        
        # Status
        tk.Label(self.prop_view.body_frame, text="Status:")\
            .grid(row=1, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.prop_status, 
            state=tk.DISABLED)\
            .grid(row=1, column=1, sticky=tk.EW)

        # Time Step
        tk.Label(self.prop_view.body_frame, text="Time Step:")\
            .grid(row=2, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.time_step,
            state=tk.DISABLED)\
            .grid(row=2, column=1, sticky=tk.EW)

        # Time Horizon
        tk.Label(self.prop_view.body_frame, text="Time Horizon:")\
            .grid(row=3, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.time_horizon,
            state=tk.DISABLED)\
            .grid(row=3, column=1, sticky=tk.EW)

        # K Value
        tk.Label(self.prop_view.body_frame, text="K Value:")\
            .grid(row=4, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.k_value,
            state=tk.DISABLED)\
            .grid(row=4, column=1, sticky=tk.EW)

        # Simulator
        tk.Label(self.prop_view.body_frame, text="Simulator:")\
            .grid(row=5, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.simulator,
            state=tk.DISABLED)\
            .grid(row=5, column=1, sticky=tk.EW)

        # Refinenment Strategy
        tk.Label(self.prop_view.body_frame, text="Refinement:")\
            .grid(row=6, column=0, sticky=tk.W)
        tk.Entry(self.prop_view.body_frame, textvariable=self.refine_strat,
            state=tk.DISABLED)\
            .grid(row=6, column=1, sticky=tk.EW)

        # Initial Set
        tk.Label(self.prop_view.body_frame, text="Initial Set:")\
            .grid(row=7, column=0, sticky=tk.W)
        self.initial_set_widget = Text(self.prop_view.body_frame, height=4, 
            width=42, wrap=tk.WORD, state=tk.DISABLED, 
            foreground='gray50', background='gray85')
        self.initial_set_widget.grid(row=8, column=0, columnspan=3)

        # Unsafe Set
        tk.Label(self.prop_view.body_frame, text="Unsafe Set:")\
            .grid(row=9, column=0, sticky=tk.W)
        self.unsafe_set_widget = Text(self.prop_view.body_frame, height=4, width=42,
            wrap=tk.WORD, state=tk.DISABLED,
            foreground='gray50', background='gray85')
        self.unsafe_set_widget.grid(row=10, column=0, columnspan=3)

    def _init_plot_list(self):
        
        # Plot List Frame
        self.plot_list = tk.LabelFrame(self, text='Plot List')

        # Display different plot names with a tree
        self.plot_tree = ttk.Treeview(self.plot_list)
        self.plot_tree.pack(fill=tk.BOTH, expand=True)
        self.plot_tree.bind('<Button-1>', self._callback_btn_press)

        self.plot_tree['show'] = 'headings'
        self.plot_tree['columns'] = ('name', 'status')
        self.plot_tree.column('name', width=150)
        self.plot_tree.column('status', width=150)
        self.plot_tree.heading('name', text="Name")
        self.plot_tree.heading('status', text="Status")

        # New, Copy, and Remove buttons
        btn_row = tk.Frame(self.plot_list)
        btn_row.pack(fill=tk.X)

        self.new_btn = tk.Button(btn_row, text="New", 
            command=self._callback_new)
        self.new_btn.pack(expand=True, fill=tk.X, side=tk.LEFT)

        self.copy_btn = tk.Button(btn_row, text="Copy",
            command=self._callback_copy)
        self.copy_btn.pack(expand=True, fill=tk.X, side=tk.LEFT)

        # Plot Buttons
        btn_row = tk.Frame(self.plot_list)  # Scrap variable reassigned
        btn_row.pack(fill=tk.X)

        self.plot_btn = tk.Button(btn_row, text="Plot", 
            command=self._callback_plot, state=tk.DISABLED)
        self.plot_btn.pack(expand=True, fill=tk.X, side=tk.LEFT)

        self.clear_btn = tk.Button(btn_row, text="Clear",
            command=self._callback_clear)
        self.clear_btn.pack(expand=True, fill=tk.X, side=tk.LEFT)

    def _callback_open(self):

        input_path_full = tk.filedialog.askopenfilename(**self.OPEN_OPT)
        if not input_path_full:
            return
        
        self.load_file(input_path_full)

    def _callback_output_path(self, *args):

        # Update Current Property
        self.cur_prop.output_path = self.output_path.get()

        # Validate Output Path
        valid = (self.cur_prop.output_path.strip() != "")
        # Update ValidLabel and Plot button
        self.output_path_vl.set_state(valid)
        self.output_path_disable = (not valid)
        self._enable_disable_plot_button()

    def _callback_plot_name(self, *args):

        # Update Current Property and Plot List display
        self.cur_prop.plot_name = self.plot_name.get()
        self.plot_tree.set(self.sel_iid, 0, self.cur_prop.plot_name)

        # Validate Plot Name
        valid = (self.cur_prop.plot_name.strip() != "")
        
        # Update ValidLabel and Plot button
        self.plot_name_vl.set_state(valid)
        self.plot_name_disable = (not valid)
        self._enable_disable_plot_button()

    def _callback_input_path(self, *args):

        self.cur_prop.input_path = self.input_path.get()

    def _callback_btn_press(self, event):

        iid = self.plot_tree.identify_row(event.y)
        if iid:
            self.sel_iid = iid
            index = self.plot_tree.index(self.sel_iid)
            self.cur_prop = self.plot_prop_list[index]
            self._load_cur_prop()
            self._update_boxes()
            self._enable_disable_plot_button()

    def _callback_new(self):

        self.cur_prop = PlotProperty(self.next_identifier)
        self.next_identifier += 1
        self.plot_prop_list.append(self.cur_prop)

        self.sel_iid = self._add_property(self.cur_prop)
        self.plot_tree.selection_set(self.sel_iid)
        self.iid_list.append(self.sel_iid)

        self._clear_all()

        self.box_disable = True  # New plots will never have boxes selected
        self._enable_disable_plot_button()

    def _callback_copy(self):

        input_path_full = self.cur_prop.input_path_full
        output_path_full = self.cur_prop.output_path_full

        # Copy current plot input
        plot_name = self.cur_prop.plot_name
        plot_status = NOT_PLOTTED
        horizontal_index = self.cur_prop.horizontal_index
        vertical_select = self.cur_prop.vertical_select
        prop_name = self.cur_prop.prop_name
        prop_status = self.cur_prop.prop_status
        time_step = self.cur_prop.time_step
        time_horizon = self.cur_prop.time_horizon
        k_value = self.cur_prop.k_value
        simulator = self.cur_prop.simulator
        refine_strat = self.cur_prop.refine_strat
        initial_set = self.cur_prop.initial_set
        unsafe_set = self.cur_prop.unsafe_set

        variable_list = self.cur_prop.variable_list
        mode_list = self.cur_prop.mode_list
        
        # Create new PlotProperty (self.cur_prop), select it, and load it
        self._callback_new()
        self.cur_prop.input_path_full = input_path_full
        self.cur_prop.output_path_full = output_path_full + "_Copy"
        self.cur_prop.plot_name = plot_name + "_Copy"
        self.cur_prop.horizontal_index = horizontal_index
        self.cur_prop.vertical_select = vertical_select
        self.cur_prop.plot_status = plot_status
        self.cur_prop.prop_name = prop_name
        self.cur_prop.prop_status = prop_status
        self.cur_prop.time_step = time_step
        self.cur_prop.time_horizon = time_horizon
        self.cur_prop.k_value = k_value
        self.cur_prop.simulator = simulator
        self.cur_prop.refine_strat = refine_strat
        self.cur_prop.initial_set = initial_set
        self.cur_prop.unsafe_set = unsafe_set
        self.cur_prop.variable_list = variable_list
        self.cur_prop.mode_list = mode_list
        self._load_cur_prop()
        self._update_boxes()

    def _callback_remove(self):

        print("Under construction")

    def _add_property(self, plot_prop):
        
        iid = self.plot_tree.insert('', 'end', 
            values=(plot_prop.plot_name, plot_prop.plot_status))
        return iid

    def _clear_all(self):

        self._clear_boxes()
        self._clear_entries()

    def _clear_boxes(self):
        """ Clear radio button and check box widgets """

        for widget in self.h_frame.winfo_children():
            widget.destroy()
        for widget in self.v_frame.winfo_children():
            widget.destroy()

    def _clear_entries(self):
        """ Clear entry widget values """

        self.input_path.set('')
        self.output_path.set('')
        self.plot_name.set('')
        self.plot_status.set('')
        self.variable_list = []
        self.mode_list = []
        self.horizontal_index.set(0)
        self.vertical_select = []
        self.vertical_indices = []

        self.prop_name.set('')
        self.prop_status.set('')
        self.time_step.set('')
        self.time_horizon.set('')
        self.k_value.set('')
        self.simulator.set('')
        self.refine_strat.set('')
        self.initial_set.set('')
        self.unsafe_set.set('')

        self.initial_set_widget.config(state=tk.NORMAL)
        self.initial_set_widget.delete('1.0', tk.END)
        self.initial_set_widget.config(state=tk.DISABLED)

        self.unsafe_set_widget.config(state=tk.NORMAL)
        self.unsafe_set_widget.delete('1.0', tk.END)
        self.unsafe_set_widget.config(state=tk.DISABLED)

    def _create_boxes(self):

        self._clear_boxes()

        self.horizontal_index = tk.IntVar()
        self.horizontal_index.set(0)
        self.vertical_select = []
        for i, var in enumerate(self.variable_list):
            tk.Radiobutton(self.h_frame, text=var, 
                variable=self.horizontal_index, value=i, 
                command=self._callback_box_select)\
                .grid(row=int(i/4), column=i%4, sticky=tk.W)
            self.vertical_select.append(tk.IntVar())
            self.vertical_select[i].set(0)
            
            tk.Checkbutton(self.v_frame, text=var, 
                variable=self.vertical_select[i], 
                command=self._callback_box_select)\
                .grid(row=int(i/4), column=i%4, sticky=tk.W)

        vertical_select = []
        for select in self.vertical_select:
            if select.get():
                vertical_select.append(1)
            else:
                vertical_select.append(0)
        self.cur_prop.vertical_select = vertical_select
        self.cur_prop.horizontal_index = self.horizontal_index.get()

    def _update_boxes(self):

        self._clear_boxes()

        self.horizontal_index = tk.IntVar()
        self.vertical_select = []
        for i, select in enumerate(self.cur_prop.vertical_select):
            self.vertical_select.append(tk.IntVar())
            self.vertical_select[i].set(select)

        for i, var in enumerate(self.variable_list):
            tk.Radiobutton(self.h_frame, text=var, 
                variable=self.horizontal_index, value=i, 
                command=self._callback_box_select)\
                .grid(row=int(i/4), column=i%4, sticky=tk.W)
            self.vertical_select[i].set(self.cur_prop.vertical_select[i])
            tk.Checkbutton(self.v_frame, text=var, 
                variable=self.vertical_select[i], 
                command=self._callback_box_select)\
                .grid(row=int(i/4), column=i%4, sticky=tk.W)
        self.horizontal_index.set(self.cur_prop.horizontal_index)
        # Enable/Disable plot button appropriately
        self._callback_box_select()  
    
    def _load_cur_prop(self):

        self.input_path.set(self.cur_prop.input_path)
        self.output_path.set(self.cur_prop.output_path)
        self.plot_name.set(self.cur_prop.plot_name)
        self.plot_status.set(self.cur_prop.plot_status)
        self.prop_name.set(self.cur_prop.prop_name)
        self.prop_status.set(self.cur_prop.prop_status)
        self.time_step.set(self.cur_prop.prop_status)
        self.time_horizon.set(self.cur_prop.time_horizon)
        self.k_value.set(self.cur_prop.k_value)
        self.simulator.set(self.cur_prop.simulator)
        self.refine_strat.set(self.cur_prop.refine_strat)
        self.initial_set.set(self.cur_prop.initial_set)
        self.unsafe_set.set(self.cur_prop.unsafe_set)

        self.variable_list = self.cur_prop.variable_list
        self.mode_list = self.cur_prop.mode_list

        self.initial_set_widget.config(state=tk.NORMAL)
        self.initial_set_widget.delete('1.0', tk.END)
        self.initial_set_widget.insert(tk.END, self.cur_prop.initial_set)
        self.initial_set_widget.config(state=tk.DISABLED)

        self.unsafe_set_widget.config(state=tk.NORMAL)
        self.unsafe_set_widget.delete('1.0', tk.END)
        self.unsafe_set_widget.insert(tk.END, self.cur_prop.unsafe_set)
        self.unsafe_set_widget.config(state=tk.DISABLED)

    def _callback_box_select(self):

        vertical_select = []
        self.box_disable = True  # Disable plot button unless one or more vertical index selected      
        for select in self.vertical_select:
            if select.get():
                vertical_select.append(1)
                self.box_disable = False
            else:
                vertical_select.append(0)
        
        self.cur_prop.vertical_select = vertical_select
        self.cur_prop.horizontal_index = self.horizontal_index.get()

        self._enable_disable_plot_button()

    def _callback_plot(self):

        self.vertical_indices = []
        for i, select in enumerate(self.vertical_select):
            if select.get():
                self.vertical_indices.append(i)

        # Make sure to use .get() on tk vars - send strings/ints to plot_graph
        plot_graph(self.cur_prop.input_path_full, self.horizontal_index.get(), 
            self.vertical_indices, self.variable_list, self.mode_list, 
            self.prop_status.get(), self.plot_name.get(), 
            self.cur_prop.output_path_full)

        self._update_property_status()
        self.parent.display.upload_new_image(self.cur_prop.output_path,
            self.cur_prop.identifier)

    def _callback_clear(self):

        self.parent.display.destroy_all()

        for iid in self.iid_list:
            plot_prop = self.plot_prop_list[self.plot_tree.index(iid)]
            plot_prop.plot_status = NOT_PLOTTED
            self.plot_tree.item(iid,
                values=(plot_prop.plot_name, plot_prop.plot_status))
            
    def _update_property_status(self, plotted=True):

        if plotted:
            self.cur_prop.plot_status = "Plotted"
        else:
            self.cur_prop.plot_status = "Not Plotted"

        self.plot_tree.item(self.sel_iid,
            values=(self.cur_prop.plot_name, self.cur_prop.plot_status))

    def _enable_disable_plot_button(self):

        if self.plot_name_disable or self.output_path_disable or \
            self.box_disable:
            self.plot_btn.config(state=tk.DISABLED)
        else:
            self.plot_btn.config(state=tk.NORMAL)
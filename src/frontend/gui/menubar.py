from tkinter import *
from tkinter import filedialog
from tkinter.ttk import *
import xml.dom.minidom

from frontend.gui.eventhandler import EventHandler
from frontend.gui.modeltab import *
from frontend.mod.filehandler import FileHandler
from frontend.mod.constants import *
from frontend.mod.session import *


class MenuBar(Menu):

    def __init__(self, parent, notebook):
        Menu.__init__(self, parent)
    
        self.parent = parent
        self.notebook = notebook
        self.tree = notebook.model_tab.tree
        self.model_sidebar = notebook.model_tab.sidebar
        self.plot_sidebar = notebook.plot_tab.sidebar

        # Open file constants
        self.OPEN_OPT = {
            'defaultextension': ['.hyxml','.mdl'],
            'filetypes': [('HyXML files', '.hyxml'), 
                ('Simulink files', '.mdl')],
            'multiple': 0,
            'parent': self.parent,
            'title': 'Open File',
        }

        """
        self.SAVE_OPT = {
            'defaultextension': '.hyxml',
            'filetypes': [('HyXML files', '.hyxml')],
            'parent': self.parent,
            'title': 'Save File',
        }
        """

        self._init_widgets()

    def disable_sim_ver(self, disable):

        if disable:
            state_ = DISABLED
        else:
            state_ = NORMAL
        # Simulate and Verify are the 3rd and 4th entires, respectively
        self.build_menu.entryconfig(2, state=state_)
        self.build_menu.entryconfig(3, state=state_)

    def _init_widgets(self):
        """ Initialize menus """

        # File Menu
        file_menu = Menu(self, tearoff=0)
        file_menu.add_command(label='Open', accelerator='Crtl+O', underline=0,
            command=self.open_callback)
        file_menu.add_command(label='New', accelerator='Ctrl+N', underline=0, 
            command=self.new_callback)
        file_menu.add_command(label='Save', accelerator='Ctrl+S', underline=0, 
            command=self.save_callback)
        file_menu.add_command(label='Save As', accelerator='Ctrl+Shift+S', 
            underline=5, command=self.save_as_callback)
        file_menu.add_command(label='Close', accelerator='Ctrl+L', underline=1, 
            command=self.close_callback)
        file_menu.add_command(label='Quit', accelerator='Ctrl+Q', underline=0, 
            command=self.quit_callback)
        self.add_cascade(label='File', menu=file_menu)

        # Build Menu
        self.build_menu = Menu(self, tearoff=0, postcommand=self._update_build_menu)
        self.build_menu.add_command(label="Parse",
            command=self.model_sidebar._callback_parse)
        self.build_menu.add_command(label="Compose",
            command=self.model_sidebar._callback_compose)
        self.build_menu.add_command(label="Verify",
            command=self.model_sidebar._callback_ver)
        self.build_menu.add_command(label="Simulate",
            command=self.model_sidebar._callback_sim)
        self.add_cascade(label="Build", menu=self.build_menu)

        # Model Menu
        model_menu = Menu(self, tearoff=0, postcommand=self._update_model_menu)
        
        #     Automata Submenu
        self.automata_submenu = Menu(self, tearoff=0)
        self.automata_submenu.add_command(label="Add Automaton",
            command=lambda: self.tree.launch_entry_popup(AUTOMATON, ADD))
        self.automata_submenu.add_command(label="Edit Automaton",
            command=lambda: self.tree.launch_entry_popup(AUTOMATON, EDIT))
        self.automata_submenu.add_command(label="Delete Automaton",
            command=lambda: self.tree.launch_entry_popup(AUTOMATON, DELETE))
                
        #     Variables Submenu
        self.variables_submenu = Menu(self, tearoff=0)
        self.variables_submenu.add_command(label="Edit Variables",
            command=lambda: self.tree.launch_entry_popup(VARIABLES, EDIT))
        
        #     Modes Submenu
        self.modes_submenu = Menu(self, tearoff=0)
        self.modes_submenu.add_command(label="Add Mode",
            command=lambda: self.tree.launch_entry_popup(MODES, ADD))
        self.modes_submenu.add_command(label="Edit Mode",
            command=lambda: self.tree.launch_entry_popup(MODES, EDIT))
        self.modes_submenu.add_command(label="Delete Mode",
            command=lambda: self.tree.launch_entry_popup(MODES, DELETE))

        #     Transitions Submenu
        self.transitions_submenu = Menu(self, tearoff=0)
        self.transitions_submenu.add_command(label="Add Transition",
            command=lambda: self.tree.launch_entry_popup(TRANSITIONS, ADD))
        self.transitions_submenu.add_command(label="Edit Transtion",
            command=lambda: self.tree.launch_entry_popup(TRANSITIONS, EDIT))
        self.transitions_submenu.add_command(label="Delete Transition",
            command=lambda: self.tree.launch_entry_popup(TRANSITIONS, DELETE))

        #     Add Cascade Menus
        model_menu.add_cascade(label="Automata", menu=self.automata_submenu)
        model_menu.add_cascade(label="Variables", menu=self.variables_submenu)
        model_menu.add_cascade(label="Modes", menu=self.modes_submenu)
        model_menu.add_cascade(label="Transitions", 
            menu=self.transitions_submenu)
        self.add_cascade(label="Model", menu=model_menu)
        
        # Requirements Menu
        self.requirements_menu = Menu(self, tearoff=0, 
            postcommand=self._update_requirements_menu)
        self.requirements_menu.add_command(label="New", 
            command=self.model_sidebar._callback_new)
        self.requirements_menu.add_command(label="Copy",
            command=self.model_sidebar._callback_cpy)
        self.requirements_menu.add_command(label="Remove",
            command=self.model_sidebar._callback_rmv)
        self.add_cascade(label="Requirements", menu=self.requirements_menu)

        # Plotter Menu
        self.plotter_menu = Menu(self, tearoff=0, 
            postcommand=self._update_plotter_menu)
        self.plotter_menu.add_command(label="New",
            command=self.plot_sidebar._callback_new)
        self.plotter_menu.add_command(label="Copy",
            command=self.plot_sidebar._callback_copy)
        self.plotter_menu.add_command(label="Remove",
            command=self.plot_sidebar._callback_remove)
        self.plotter_menu.add(SEPARATOR)
        self.plotter_menu.add_command(label="Plot",
            command=self.plot_sidebar._callback_plot)
        self.add_cascade(label="Plotter", menu=self.plotter_menu)

        # Help menu
        help_menu = Menu(self, tearoff=0)
        help_menu.add_command(label='About', command=self._callback_about)
        self.add_cascade(label='Help', menu=help_menu)
        self.parent.config(menu=self)

        # Bind accelerators
        self.parent.bind_all('<Control-o>', lambda event: self.open_callback())
        self.parent.bind_all('<Control-n>', lambda event: self.new_callback())
        self.parent.bind_all('<Control-s>', lambda event: self.save_callback())
        self.parent.bind_all('<Control-Shift-s>', 
            lambda event: self.save_as_callback())
        self.parent.bind_all('<Control-l>', 
            lambda event: self.close_callback())
        self.parent.bind_all('<Control-q>', lambda event: self.quit_callback())

    def _update_build_menu(self):
        """ Disable build menu options based on selected tab """

        # Enable/Disable based on selected tab
        if (not Session.file_opened) or (self.notebook.current_tab != MODEL):
            state_ = DISABLED
        else:
            state_ = NORMAL
        for i in range(4):
            self.build_menu.entryconfig(i, state=state_)

        # Enable/Disable based on model sidebar
        if self.model_sidebar.check_sim_ver_disable():
            self.build_menu.entryconfig(2, state=DISABLED)  # Simualte
            self.build_menu.entryconfig(3, state=DISABLED)  # Verify

    def _update_model_menu(self):
        """ Disable model menu options based on item selected in TreeView """

        self.variables_submenu.entryconfig(0, state=DISABLED)
        for i in range(3):
            self.automata_submenu.entryconfig(i, state=DISABLED)
            self.modes_submenu.entryconfig(i, state=DISABLED)
            self.transitions_submenu.entryconfig(i, state=DISABLED)

        if (not Session.file_opened) or (self.notebook.current_tab != MODEL):
            return
        
        self.automata_submenu.entryconfig(0, state=NORMAL) 
        if self.tree.slct_automaton is None:
            return
     
        self.automata_submenu.entryconfig(0, state=NORMAL)
        self.automata_submenu.entryconfig(1, state=NORMAL)
        self.automata_submenu.entryconfig(2, state=NORMAL)
        self.variables_submenu.entryconfig(0, state=NORMAL)
        self.modes_submenu.entryconfig(0, state=NORMAL)
        self.transitions_submenu.entryconfig(0, state=NORMAL)

        if self.tree.slct_mode is not None:       
            self.modes_submenu.entryconfig(1, state=NORMAL)
            self.modes_submenu.entryconfig(2, state=NORMAL)
        elif self.tree.slct_transition is not None:
            self.transitions_submenu.entryconfig(1, state=NORMAL)
            self.transitions_submenu.entryconfig(2, state=NORMAL)

    def _update_requirements_menu(self):
        """ Disable requirements menu options based on selected tab """

        if (not Session.file_opened) or (self.notebook.current_tab != MODEL):
            state_ = DISABLED
        else:
            state_ = NORMAL

        for i in range(3):
            self.requirements_menu.entryconfig(i, state=state_)

    def _update_plotter_menu(self):
        """ Disable requirements menu options based on selected tab """

        if (not Session.file_opened) or (self.notebook.current_tab != PLOT):
            state_ = DISABLED
        else:
            state_ = NORMAL

        for i in range(3): 
            self.plotter_menu.entryconfig(i, state=state_)
        # Entry 3 is the separator
        self.plotter_menu.entryconfig(4, state=state_)  

    def open_callback(self):
        """ Select and open file """

        # Forget welcome screen widgets
        self.parent.open_label.pack_forget() 
        self.parent.manual_label.pack_forget()
        self.parent.email_label.pack_forget()

        # If a file is already open, close it.
        if Session.file_opened:
            self.close_callback()
       
        file_path = filedialog.askopenfilename(**self.OPEN_OPT)
        if file_path:            
            Session.file_path = file_path
            status = FileHandler.open_file(file_path)
            if status:
                Session.file_opened = True
                Session.file_saved = True
                EventHandler.event_generate(OPEN_EVENT)

    def new_callback(self):
        """ Open template file """

        # Forget welcome screen widgets
        self.parent.open_label.pack_forget()
        self.parent.manual_label.pack_forget()
        self.parent.email_label.pack_forget()

        # If a file is already open, close it.
        if Session.file_opened:
            self.close_callback()

        Session.hybrid = HyIR.create_template()
        Session.file_path = None
        Session.file_opened = True
        Session.file_saved = False
        EventHandler.event_generate(OPEN_EVENT)

        return

    def save_callback(self):

        if self.notebook.current_tab == EDITOR:
            FileHandler.save(self.notebook.editor_tab.editor.get())
            # Load XML into HyIR object
            FileHandler.open_file(Session.file_path)
            # Refresh Model Tab
            self.tree._clear_model()
            self.tree._display_model()
            # Refresh Property Sidebar
            self.notebook.model_tab.sidebar._display_property(Session.cur_prop)
        else:
            FileHandler.save()
            # Refresh Editor Tab
            self.editor_tab.open_xml()
        
        return

    def save_as_callback(self, event=None):

        if self.notebook.current_tab == EDITOR:
            FileHandler.save_as(
                self.notebook.editor_tab.editor.get('1.0', 'end-1c'))
        else:
            FileHandler.save_as()
        
        return 

    def close_callback(self, event=None):

        if not Session.file_saved:
            if messagebox.askyesno('Save', "Save file before closing?"):
                self.save_callback()

        Session.lib_compiled = False
        Session.simulator = ODEINT_FIX
        Session.refine_strat = DEF_STRAT
        Session.file_saved = False
        Session.file_type = ''
        Session.file_path = ''

        if Session.file_opened:
            EventHandler.event_generate(CLOSE_EVENT)
        Session.file_opened = False

    def quit_callback(self, event=None):

        self.parent.destroy()

    def _callback_about(self):

        about = tk.Toplevel(self)
        about.title('About')
        about.config(bg='white')

        img = PIL.ImageTk.PhotoImage(PIL.Image.open('./frontend/res/logosmall.png'))
        img_label = tk.Label(about, image=img)
        img_label.image = img
        img_label.grid(row=0, column=0, columnspan=2)
        
        tk.Label(about, text="Version: ", bg='white').grid(row=1, column=0, sticky=tk.W)
        tk.Label(about, text="2.0.0", bg='white').grid(row=1, column=1, sticky=tk.E)
        tk.Label(about, text="Date: ", bg='white').grid(row=2, column=0, sticky=tk.W)
        tk.Label(about, text="August 31, 2018", bg='white').grid(row=2, column=1, sticky=tk.E)
        tk.Label(about, text="Support: ", bg='white').grid(row=3, column=0, sticky=tk.W) 
        tk.Label(about, text="c2e2help@gmail.com", bg='white').grid(row=3, column=1, sticky=tk.E)

        x = Session.window.winfo_pointerx()
        y = Session.window.winfo_pointery()
        about.geometry("+%d+%d" % (x, y))
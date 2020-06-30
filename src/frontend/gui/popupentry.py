from tkinter import *
from tkinter.ttk import *
from frontend.gui.widgets import ToggleFrame
from frontend.mod.automaton import *
from frontend.mod.constants import *
from frontend.mod.hyir import *
from frontend.mod.session import Session


class PopupEntry(Toplevel):

    def __init__(self, parent):
        Toplevel.__init__(self, parent)

        self.parent = parent
        self.resizable(width=False, height=False)

        self.title_label = Label(self, text='C2E2')
        self.title_label.grid(row=0, column=0, columnspan=2)
        
        self.TEXTBOX_HEIGHT = 10
        self.TEXTBOX_WIDTH = 30

        # Window appears by cursor position
        self.geometry("+%d+%d" % (Session.window.winfo_pointerx(), 
            Session.window.winfo_pointery()))

        # Prevent interaction with main window until Popup is Confirmed/Canceled
        self.wait_visibility()
        self.focus_set()
        self.grab_set()


class AutomatonEntry(PopupEntry):
    """
    Popup window for adding/deleting Automata from the hybrid system.

    Args:
        parent (obj): Popup's parent object
        hybrid (obj): Hybrid object - should always be Session.hybrid
        action (str): Action to be performed (ADD or DELETE)
    """

    def __init__(self, parent, hybrid, action, automaton=None):
        PopupEntry.__init__(self, parent)
        self.title_label.config(text="Automaton")

        if hybrid is not Session.hybrid:
            Session.write("ERROR: Attempting to edit non-Session hybrid.\n")
            self._cancel()

        self.parent = parent
        self.hybrid = hybrid
        self.automaton = automaton
        self.action = action
        self.changed = False

        self._init_widgets()

        if action == EDIT:
            self._load_session()

        if action == DELETE:
            self._disable_fields()
  
    def _init_widgets(self):
        """ Initialize GUI elements """

        # Name
        Label(self, text="Name:").grid(row=1, column=0, sticky=W)
        self.name = StringVar()
        self.name_entry = Entry(self, textvariable=self.name)
        self.name_entry.grid(row=1, column=1, sticky=E)

        # Buttons
        self.btn_frame = Frame(self)
        
        self.cancel_btn = Button(self.btn_frame, text="Cancel", 
                                 command=self._cancel)
        self.confirm_btn = Button(self.btn_frame, text="Confirm", 
                                  command=self._confirm)

        self.cancel_btn.grid(row=0, column=0)
        self.confirm_btn.grid(row=0, column=1)

        self.btn_frame.grid(row=2, column=0, columnspan=2)

        return

    def _load_session(self):

        # Name
        self.name.set(self.automaton.name)

        return

    def _disable_fields(self):

        # Name
        self.name_entry.config(state=DISABLED)

        self.confirm_btn.config(text="DELETE", command=self._delete)

        return

    def _confirm(self):

        if(self.action == ADD):
            self._confirm_add()
        else:
            self._confirm_edit()

        return

    def _confirm_add(self):

        self.hybrid.add_automaton(Automaton(self.name.get()))

        Session.write("Automaton Entry Confirmed.\n")
        self.changed = True
        self.destroy()

        return

    def _confirm_edit(self):

        self.automaton.name = self.name.get()
        
        Session.write("Automaton Entry Confirmed.\n")
        self.changed = True
        self.destroy()

        return

    def _delete(self):

        if messagebox.askyesno("Delete Automaton", 
            "Delete " + self.automaton.name + "?"):
            
            self.hybrid.remove_automaton(self.automaton)
            Session.write("Automaton Deleted.\n")
            self.changed = True
        else:
            Session.write("Automaton Deletion Canceled.\n")
            self.chagned = False
        self.destroy()

        return

    def _cancel(self):
        """ Cancels changes made in popup """

        Session.write("Automaton Entry Canceled.\n")
        self.changed = False
        self.destroy()

        return


class VariableEntry(PopupEntry):
    """ 
    Popup window for Variable editing.    

    The VariableEntry class is designed to be the popup displayed to users when 
    editing their model's variables. It controls the GUI elements of the popup, 
    and interacts with the Session variables to commit changes to the currently 
    active model
    
    Args:
        parent (obj): Popup's parent object
    """

    def __init__(self, parent, automaton):
        PopupEntry.__init__(self, parent)
        self.title_label.config(text="Variables")

        self.automaton = automaton
        self.changed = False

        # For readability, options differ from those stored in the var object
        self.scope_options = ('Local', 'Input', 'Output')

        self._init_widgets()
        self._load_session()


    def _init_widgets(self):
        """ Initialize GUI elements """

        self.title_label.grid(row=0, column=0, columnspan=4)

        Label(self, text="Name").grid(row=1, column=0)
        Label(self, text="Thin").grid(row=1, column=1)
        Label(self, text="Type").grid(row=1, column=2)
        Label(self, text="Scope").grid(row=1, column=3)

        # Variable lists for uknown number of inputs
        self.names = []  # StringVar()
        self.thins = []  # BoolVar()
        self.types = []  # StringVar()
        self.scopes = [] # StringVar()
        self.var_index = 0

        # Buttons
        self.btn_frame = Frame(self)

        self.cancel_btn = Button(self.btn_frame, text="Cancel", 
                                 command=self._cancel)
        self.add_btn = Button(self.btn_frame, text="Add", 
                              command=self._add_row)
        self.confirm_btn = Button(self.btn_frame, text="Confirm", 
                                  command=self._confirm)

        self.cancel_btn.grid(row=0, column=0)
        self.add_btn.grid(row=0, column=1)
        self.confirm_btn.grid(row=0, column=2)

        return
       
    def _load_session(self):
        """ Load current model's values. """

        scope_dict = {  # Convert Variable scopes to options displayed to user
            LOCAL: 'Local',  # LOCAL = 'LOCAL_DATA'
            INPUT: 'Input',  # INPUT = 'INPUT_DATA'
            OUTPUT: 'Output' # OUTPUT = 'OUTPUT_DATA'
        }

        # Add a blank row if there are no variables (happens with new automata)
        if len(self.automaton.vars) == 0 and len(self.automaton.thinvars) == 0:
            self._add_row()
            return

        for var in self.automaton.vars:
            self._add_row()
            self.names[self.var_index-1].set(var.name)
            self.thins[self.var_index-1].set(False)
            self.types[self.var_index-1].set(var.type)
            self.scopes[self.var_index-1].set(scope_dict[var.scope])

        for var in self.automaton.thinvars:
            self._add_row()
            self.names[self.var_index-1].set(var.name)
            self.thins[self.var_index-1].set(True)
            self.types[self.var_index-1].set(var.type)
            self.scopes[self.var_index-1].set(scope_dict[var.scope])

        return

    def _add_row(self):
        """ 
        Add a new variable row to VariableEntry popup. 
        Grid new entry widgets and regrid button frame.
        """

        self.names.append(StringVar())
        self.thins.append(BooleanVar())
        self.types.append(StringVar())
        self.scopes.append(StringVar())

        # Name
        Entry(self, textvariable=self.names[self.var_index])\
            .grid(row=self.var_index+2, column=0)

        # Thin
        Checkbutton(self, var=self.thins[self.var_index])\
            .grid(row=self.var_index+2, column=1)

        # Type
        self.types[self.var_index].set(REAL)
        OptionMenu(self, self.types[self.var_index], 
                   self.types[self.var_index].get(), 
                   *VARIABLE_TYPES)\
                   .grid(row=self.var_index+2, column=2)

        # Scope
        self.scopes[self.var_index].set('Local')
        OptionMenu(self, self.scopes[self.var_index], 
                   self.scopes[self.var_index].get(), 
                   *self.scope_options)\
                   .grid(row=self.var_index+2, column=3)

        self.btn_frame.grid(row=self.var_index+3, columnspan=4)

        self.var_index += 1

        return

    def _confirm(self):
        """ Commit changes to Session. Does NOT save these changes. """

        self.automaton.reset_vars()
        self.automaton.reset_thinvars()

        scope_dict = {  # Convert displayed scopes to values stored
            'Local': LOCAL,  # LOCAL = 'LOCAL_DATA'
            'Input': INPUT,  # INPUT = 'INPUT_DATA'
            'Output': OUTPUT # OUTPUT = 'OUTPUT_DATA'
        }

        for i in range(0, self.var_index):

            name = (self.names[i].get()).strip()
            thin = self.thins[i].get()
            type_ = self.types[i].get()  # Reserved word
            scope = scope_dict[self.scopes[i].get()]

            if not name:  # Delete variables by erasing their name 
                continue

            if thin:
                self.automaton.add_thinvar(
                    Variable(name=name, type=type_, scope=scope))
            else:
                self.automaton.add_var(
                    Variable(name=name, type=type_, scope=scope))
            
        Session.write("Variable Entry Confirmed.\n")
        self.changed = True
        self.destroy()

        return
 
    def _cancel(self):
        """ Cancels changes made in popup """

        Session.write("Variable Entry Canceled.")
        self.changed = False
        self.destroy()

        return


class ModeEntry(PopupEntry):
    """ 
    Popup window for Mode adding, editing, and deleting.

    The ModelEntry class is designed to be the popup displayed to users when
    editing their model's Modes, or adding/deleting Modes. It controls the GUI
    elements of the popup, and interacts with the Session variables to commit 
    changes to the currently active models.
    
    Args:
        parent (obj): Popup's parent object
        action (str): Action to be performed (constants ADD, EDIT, or DELETE)
        mode (Mode obj): Mode to be edited or deleted, not required for ADD 
    """

    def __init__(self, parent, automaton, action=ADD, mode=None):
        PopupEntry.__init__(self, parent)
        self.title_label.config(text='Mode')

        self.automaton = automaton
        self.mode = mode
        self.action = action
        self.mode_dict = automaton.mode_dict  # mode_dict[mode.id] = mode.name
        self.changed = False

        self._init_widgets()

        if(action == ADD): 
            self._load_new()
        else:
            self._load_session()
            if(action == DELETE):
                self._disable_fields()


    def _init_widgets(self):
        """ Initialize GUI elements """ 

        # Name
        Label(self, text='Name:').grid(row=1, column=0, sticky=W)
        self.name = StringVar()
        self.name_entry = Entry(self, textvariable=self.name)
        self.name_entry.grid(row=1, column=1, sticky=E)

        # ID
        Label(self, text='ID:').grid(row=2, column=0, sticky=W)
        self.mode_id = IntVar()
        self.id_entry = Entry(self, textvariable=self.mode_id, state=DISABLED)
        self.id_entry.grid(row=2, column=1, sticky=E)

        # Initial
        Label(self, text='Initial:').grid(row=3, column=0, sticky=W)
        self.initial = BooleanVar()
        self.initial_checkbutton = Checkbutton(self, var=self.initial)
        self.initial_checkbutton.grid(row=3, column=1)

        # Flows
        self.flow_toggle = ToggleFrame(self, text='Flows:')
        self.flow_toggle.grid(row=4, column=0, columnspan=2, sticky=E+W)

        # Invariants
        self.invariant_toggle = ToggleFrame(self, text='Invariants:')
        self.invariant_toggle.grid(row=5, column=0, columnspan=2, sticky=E+W)

        # Buttons
        
        self.btn_frame = Frame(self)

        self.cancel_btn = Button(self.btn_frame, 
                                 text='Cancel', 
                                 command=self._cancel)
        self.confirm_btn = Button(self.btn_frame, 
                                  text='Confirm', 
                                  command=self._confirm)
        
        self.cancel_btn.grid(row=0, column=0)
        self.confirm_btn.grid(row=0, column=1)

        self.btn_frame.grid(row=8, column=0, columnspan=2)

        return

    def _load_session(self):
        """ Load selected mode's Session values """

        # Name
        self.name.set(self.mode.name)

        # ID
        self.mode_id.set(self.mode.id)

        # Initial
        self.initial.set(self.mode.initial)

        # Flows
        if(len(self.mode.dais) < 1):
            self.flow_toggle.add_row()
        else:
            for dai in self.mode.dais:
                self.flow_toggle.add_row(text=dai.raw)
        self.flow_toggle.toggle()
        
        # Invariants
        if(len(self.mode.invariants) < 1):
            self.invariant_toggle.add_row()
        else:
            for invariant in self.mode.invariants:
                self.invariant_toggle.add_row(text=invariant.raw)
        self.invariant_toggle.toggle()

        return

    def _load_new(self):
        """ Load blank row and show toggle fields"""

        self.flow_toggle.add_row()
        self.flow_toggle.toggle()

        self.invariant_toggle.add_row()
        self.invariant_toggle.toggle()

        self.mode_id.set(self.automaton.next_mode_id)

        return

    def _disable_fields(self):
        """ Disable fields and reconfigure confirm button for deletion """

        self.name_entry.config(state=DISABLED)
        self.id_entry.config(state=DISABLED)
        self.initial_checkbutton.config(state=DISABLED)
        self.flow_toggle.disable_fields()
        self.invariant_toggle.disable_fields()

        self.confirm_btn.config(text='DELETE', command=self._delete)

        return

    def _confirm(self):
        """ Confirm button callback - call confirm method based on action """

        if(self.action == ADD):
            self._confirm_add()
        else:
            self._confirm_edit()

        return

    def _confirm_add(self):
        """ Confirm new mode addition """
        
        self.mode = Mode()
        self._confirm_edit()
        self.automaton.add_mode(self.mode)

        return

    def _confirm_edit(self):
        """ Commit changes to Session. Does NOT save changes """

        # Name
        self.mode.name = self.name.get()
        
        # ID
        self.mode.id = self.mode_id.get()
        
        # Initial
        self.mode.initial = self.initial.get()

        # Flows
        self.mode.clear_dais()
        for raw_text in self.flow_toggle.get_rows():
            if((raw_text.get()).strip()): 
                self.mode.add_dai(DAI(raw_text.get()))

        # Invariants
        self.mode.clear_invariants()
        for raw_text in self.invariant_toggle.get_rows():
            if((raw_text.get()).strip()): 
                self.mode.add_invariant(Invariant(raw_text.get()))
        
        Session.write("Mode Entry Confirmed.\n")
        self.changed = True
        self.destroy()

        return

    def _delete(self):
        """ Delete active Mode """

        # Build list of transitions that would be deleted
        del_trans = []
        for tran in self.automaton.transitions:
            if((tran.source == self.mode.id) or \
               (tran.destination == self.mode.id)):
                del_trans.append(tran)

        # Messagebox warning user of transitions that also will be deleted
        msg = "Delete " + self.mode.name + "(" + str(self.mode.id) + ") ?\n"
        msg += "WARNING: The following transitions will also be deleted:\n"
        for tran in del_trans:
            msg += tran.name + '\n'
                
        if(messagebox.askyesno('Delete Mode', msg)):
            self.automaton.remove_mode(self.mode)
            for tran in del_trans:
                self.automaton.remove_transition(tran)
        
            Session.write("Mode Deleted.\n")
            self.changed = True
        else:
            Session.write("Mode Deletion Canceled.\n")
            self.changed = False
        self.destroy()

        return

    def _cancel(self):
        """ Cancels changes made in popup """

        Session.write("Mode Entry Canceled.\n")
        self.changed = False
        self.destroy()

        return


class TransitionEntry(PopupEntry):
    """ 
    Popup window for Transition adding, editing, and deleting.

    The TransitionEntry class is designed to be the popup displayed to users 
    when editing their model's Modes, or adding/deleting Modes. It controls the
    GUI elements of the popup, and interacts with the Session variables to 
    commit changes to the currently active models.
    
    Args:
        parent (obj): Popup's parent object
        action (str): Action to be performed (constants ADD, EDIT, or DELETE)
        mode_dict (dictionary: int keys, str values): Dictionary connect mode 
                                                      IDs to mode names
        trans (Transition obj): Transition to be edited or deleted, not 
                                required for ADD action
    """    

    def __init__(self, parent, automaton, action=ADD, transition=None):
        PopupEntry.__init__(self, parent)
        self.title_label.config(text='Transition')

        self.automaton = automaton
        self.transition = transition
        self.mode_dict = automaton.mode_dict  # mode_dict[mode.id] = mode.name
        self.action = action
        self.changed = False

        # Load Mode list for Source/Destination Option Menus
        self.mode_list = []
        for mode_id in self.mode_dict:
            self.mode_list.append(self.mode_dict[mode_id])

        self._init_widgets()

        if(action == ADD): 
            self._load_new()
        else:
            self._load_session()
            if(action == DELETE):
                self._disable_fields()

    def _init_widgets(self):
        """ Initialize GUI elements """

        # Transition Label
        self.transition_str = StringVar()
        Label(self, textvariable=self.transition_str).grid(row=1, column=0, columnspan=2)

        # ID
        Label(self, text='ID:').grid(row=2, column=0, sticky=W)
        self.transition_id = IntVar()
        self.id_entry = Entry(self, textvariable=self.transition_id,
            state=DISABLED)
        self.id_entry.grid(row=2, column=1, sticky=E)

        # Source and Destination
        
        Label(self, text='Source:').grid(row=3, column=0, sticky=W)
        Label(self, text='Destination:').grid(row=4, column=0, sticky=W)

        self.source_str = StringVar()
        self.destination_str = StringVar()

        self.source_str.trace_variable('w', self._callback_mode_select)
        self.destination_str.trace_variable('w', self._callback_mode_select)

        # Arbitrarily set default source/destination. 
        # These are overwritten to be correct in _load_session when appropriate
        self.source_option_menu = OptionMenu(self, 
                                             self.source_str, 
                                             self.mode_list[0], 
                                             *self.mode_list)
        self.source_option_menu.grid(row=3, column=1, sticky=W+E)        
        self.destination_option_menu = OptionMenu(self, 
                                                  self.destination_str, 
                                                  self.mode_list[0], 
                                                  *self.mode_list)
        self.destination_option_menu.grid(row=4, column=1, sticky=W+E)

        # Guards
        Label(self, text='Guards:').grid(row=5, column=0, sticky=W)
        self.guard_str = StringVar()
        self.guard_entry = Entry(self, textvariable=self.guard_str)
        self.guard_entry.grid(row=5, column=1, sticky=E)

        # Actions
        self.action_toggle = ToggleFrame(self, text='Actions:')
        self.action_toggle.grid(row=6, column=0, columnspan=2, sticky=E+W)

        # Buttons
         
        self.btn_frame = Frame(self)

        self.cancel_btn = Button(self.btn_frame, 
                                 text='Cancel', 
                                 command=self._cancel)
        self.confirm_btn = Button(self.btn_frame, text='Confirm', command=self._confirm)

        self.cancel_btn.grid(row=0, column=0)
        self.confirm_btn.grid(row=0, column=1)

        self.btn_frame.grid(row=7, column=0, columnspan=2)

        return

    def _load_session(self):
        """ Load selected transition's Session values """

        # ID
        self.transition_id.set(self.transition.id)

        # Source and Destination
        self.source_str.set(self.mode_dict[self.transition.source])
        self.destination_str.set(self.mode_dict[self.transition.destination])

        # Guard
        self.guard_str.set(self.transition.guard.raw)

        # Actions
        if len(self.transition.actions) == 0:
            self.action_toggle.add_row()
        else:
            for action in self.transition.actions:
                self.action_toggle.add_row(text=action.raw)
        self.action_toggle.toggle()

        return
    
    def _load_new(self):
        """ Load blank rows and show toggle fields """

        self.action_toggle.add_row()
        self.action_toggle.toggle()

        self.transition_id.set(len(self.automaton.transitions))

        return

    def _disable_fields(self):
        """ Disable fields and reconfigure confirm button for deletion """

        self.id_entry.config(state=DISABLED)
        self.source_option_menu.config(state=DISABLED)
        self.destination_option_menu.config(state=DISABLED)
        self.guard_entry.config(state=DISABLED)
        self.action_toggle.disable_fields()

        self.confirm_btn.config(text='DELETE', command=self._delete)

        return

    def _callback_mode_select(self, *args):
        """ OptionMenu callback, updates transition label at top of window """

        self.transition_str.set(self.source_str.get() + " -> " + self.destination_str.get())

        return

    def _confirm(self):
        """ Confirm button callback - call confirm method based on action """
        
        if(self.action == ADD):
            self._confirm_add()
        else:
            self._confirm_edit()

        return

    def _confirm_add(self):
        """ Confirm new mode addition """

        # ID
        trans_id = self.transition_id.get()

        # Source and Destination
        for mode_id in self.mode_dict:
            if(self.mode_dict[mode_id] == self.source_str.get()):
                src = mode_id
            elif(self.mode_dict[mode_id] == self.destination_str.get()):
                dest = mode_id

        # Guard
        guard = Guard(self.guard_str.get())

        # Actions
        actions = []
        for action in self.action_toggle.get_rows():
            if((action.get()).strip()):
                actions.append(Action(action.get()))

        transition = Transition(guard, actions, trans_id, src, dest)
        self.automaton.add_transition(transition)

        Session.write("Transition Entry Confirmed.\n")
        self.changed = True
        self.destroy()     
        
        return

    def _confirm_edit(self):
        """" Commits changes to Session. Does NOT save changes """

        # ID
        self.transition.id = self.transition_id.get()
        
        # Source and Destination
        for mode_id in self.mode_dict:
            if(self.mode_dict[mode_id] == self.source_str.get()):
                self.transition.source = mode_id
            elif(self.mode_dict[mode_id] == self.destination_str.get()):
                self.transition.destination = mode_id

        # Guard
        self.transition.guard = Guard(self.guard_str.get())

        # Actions
        self.transition.clear_actions()
        for action in self.action_toggle.rows:
            if((action.get()).strip()):
                self.transition.add_action(Action(action.get()))
                
        Session.write("Transition Entry Confirmed.\n")
        self.changed = True
        self.destroy()

        return

    def _delete(self):
        """ Delete active Transiiton """
        
        if messagebox.askyesno('Delete Transition', 'Delete ' + \
           self.transition_str.get() + '?'):
            self.automaton.remove_transition(self.transition)
        
            Session.write("Transition Deleted.\n")
            self.changed = True
        else:
            Session.write("Transition Deletion Canceled.\n")
            self.changed = False
        self.destroy()

        return

    def _cancel(self):
        """ Cancels changes made in popup """

        Session.write("Transition Entry Canceled.\n")
        self.changed = False
        self.destroy()

        return
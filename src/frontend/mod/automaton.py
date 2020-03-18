import sympy
from scipy import optimize as opt

from frontend.mod.constants import *
from frontend.mod.session import Session, SymEq


class Automaton:

    def __init__(self, name="default_automaton"):

        self.name = name

        self.variables = Variables(self)
        self.thinvariables = ThinVariables(self)

        self.modes = []  # Mode objects
        self.transitions = []  # Transition objects

        self.prop_list = []
                
        self.next_mode_id = 0
        self.next_transition_id = 0

        self.mode_dict = {}

        self.parent = None

    # Variable Properties

    @property
    def vars(self):
        return self.variables.all

    @property
    def local_vars(self):
        return self.variables.local

    @property
    def output_vars(self):
        return self.variables.output

    @property
    def input_vars(self):
        return self.variables.input

    @property
    def var_names(self):
        return self.variables.names

    @property
    def local_var_names(self):
        return self.variables.local_names

    @property
    def output_var_names(self):
        return self.variables.output_names

    @property
    def input_var_names(self):
        return self.variables.input_names

    # Thinvariable Properties

    @property
    def thinvars(self):
        return self.thinvariables.all

    @property
    def local_thinvars(self):
        return self.thinvariables.local

    @property
    def output_thinvars(self):
        return self.thinvariables.output

    @property
    def input_thinvars(self):
        return self.thinvariables.input

    @property
    def thinvar_names(self):
        return self.thinvariables.names

    @property
    def local_thinvar_names(self):
        return self.thinvariables.local_names

    @property
    def output_thinvar_names(self):
        return self.thinvariables.output_names

    @property
    def input_thinvar_names(self):
        return self.thinvariables.input_names

    # Mode/Transition Properties
        
    @property
    def mode_names(self):
        names = []
        for mode in self.modes:
            names.append(mode.name)
        return names

    @property
    def trans(self):
        print("******************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - AUTOMATON.TRANS")
        print("******************************************************")
        return self.transitions

    def add_var(self, var):
        var.parent = self
        self.variables.add_var(var)
        return

    def remove_var(self, var):
        var.parent = None
        self.variables.remove_var(var)
        return
    
    def reset_vars(self):
        self.variables = Variables(self)
        return

    def add_thinvar(self, thinvar):
        thinvar.parent = self
        self.thinvariables.add_thinvar(thinvar)
        return

    def remove_thinvar(self, thinvar):
        thinvar.parent = None
        self.thinvariables.remove_thinvar(thinvar)
        return

    def reset_thinvars(self):
        self.thinvariables = ThinVariables(self)
        return

    def add_mode(self, mode, id_=None):

        if id_ is not None: 
            mode.id = id_
            if id_ > self.next_mode_id:
                self.next_mode_id = id_
        else:
            mode.id = self.next_mode_id

        self.next_mode_id += 1
        mode.parent = self
        self.mode_dict[mode.id] = mode.name
        self.modes.append(mode)
        return
    
    def remove_mode(self, mode):
        mode.parent = None
        self.mode_dict[mode.id] = "MODE NOT FOUND"
        self.modes.remove(mode)
        return

    def add_transition(self, tran):
        tran.parent = self
        self.transitions.append(tran)
        return

    def remove_transition(self, tran):
        tran.parent = None
        self.transitions.remove(tran)
        return

    def parse(self):
        
        Session.write("  Parsing Automaton " + self.name + "...\n")
        errors = []

        Session.write("    Parsing Modes...")
        if not self.verify_mode_names():
            errors.append(('Mode', self, "Mode names not unique", None))
        if not self.verify_mode_ids():
            errors.append(('Mode', self, "Mode IDs not unique", None))

        for mode in self.modes:
            errors += mode.parse()
        Session.write(" Modes Parsed.\n")

        Session.write("    Parsing Transitions...")
        if not self.verify_transition_src_dest():
            errors.append(('Transition', self, 
                           "Transition src/dest invalid", None))
        
        for transition in self.transitions:
            errors += transition.parse()
        Session.write(" Transitions Parsed.\n")

        Session.write("  Automaton Parsed.\n")
        return errors

    # Verify functions

    def verify_local_variables(self):
        """ Verify local variable names are unique within an Automaton """
        
        if len(self.local_var_names) > len(set(self.local_var_names)):
            return False
        
        if len(self.local_thinvar_names) > len(set(self.local_thinvar_names)):
            return False
        
        return True

    def verify_output_variables(self):
        """ Verify output variable names are unique within an Automaton """
        
        if len(self.output_var_names) > len(set(self.output_var_names)):
            return False
        
        if len(self.output_thinvar_names) >len(set(self.output_thinvar_names)):
            return False

        return True

    def verify_input_variables(self):
        """ Verify input variable names are unique within an Automaton """
        
        if len(self.input_var_names) > len(set(self.input_var_names)):
            return False
        
        if len(self.input_thinvar_names) > len(set(self.input_thinvar_names)):
            return False
        
        return True

    def verify_mode_ids(self):
        """ Verify mode ids are unique within an Automaton """
        id_set = set()
        for mode in self.modes:
            if mode.id in id_set:
                return False
            id_set.add(mode.id)
        return True

    def verify_mode_names(self):
        """ Verify mode names are unqiue with an Automaton """
        name_set = set()
        for mode in self.modes:
            name = mode.name.strip()
            if name in name_set:
                return False
            name_set.add(name)
        return True

    def verify_transition_src_dest(self):
        """ Verify transition src and destination IDs """

        id_set = set()
        for mode in self.modes:
            id_set.add(mode.id)

        for tran in self.transitions:
            if ((tran.source not in id_set) or 
                (tran.destination not in id_set)):
                return False
        return True

    # DEPRECATED
        
    def new_mode_id(self):
        m_id = self.next_mode_id
        self.next_mode_id += 1
        return m_id
    
    def new_transition_id(self):
        t_id = self.next_transition_id
        self.next_transition_id += 1
        return t_id

    # Prints
    
    def print_trans(self):
        print("Transition printing under construction")
        # print("--- Transitions ---")
        # for i in self.trans:
        #     print("%s: %s %s %s Actions: %s" % \
        #     (i.id, i.src, i.dest, i.guard.raw, 
        #     ", ".join(action.raw for action in i.actions)))
        return
    
    def print_modes(self):
        print("Mode printing under construction")
        # print("--- Modes ---")
        # for i in self.modes:
        #     print('DAIs:')
        #     print(str(i.id)+": "+i.name+" Initial: "+str(i.initial)+" Linear: "+str(i.linear))
        #     print('\n'.join(dai.raw for dai in i.dais))
        #     print('Invariants:')
        #     print('\n'.join(inv.raw for inv in i.invariants))
        return
    
    def print_all(self):
        print("%s:" % self.name)
        self.print_modes()
        self.print_trans()
        return


class Variables:
    
    def __init__(self, parent):
    
        self.local = []
        self.input = []
        self.output = []
        self.parent = parent

    @property
    def names(self):
        names = []
        for var in self.all:
            names.append(var.name)
        return names

    @property
    def local_names(self):
        names = []
        for var in self.local:
            names.append(var.name)
        return names

    @property
    def output_names(self):
        names = []
        for var in self.output:
            names.append(var.name)
        return names

    @property
    def input_names(self):
        names = []
        for var in self.input:
            names.append(var.name)
        return names

    @property
    def all(self):
        return (self.input + self.local + self.output)

    def add_var(self, v):
        """ Add Variable to appropriate list based on scope """
   
        if v.scope=='LOCAL_DATA':
            self.local.append(v)
        elif v.scope=='OUTPUT_DATA':
            self.output.append(v)
        elif v.scope=='INPUT_DATA':
            self.input.append(v)
    
        return

    def update_parents(self, parent):
        """ Update Variable parents. Needed for composition """

        for var in (self.local + self.input + self.output):
            var.parent = parent
        self.parent = parent

        return


class Variable:
    
    def __init__(self, name="default_variable", update_type="", 
                 type=REAL, scope='LOCAL_DATA'):
        self.name = name
        self.type = type
        self.scope = scope
        self.update_type = type
        self.parent = None

    def __eq__(self, other):
        return (self.name==other.name\
                and self.update_type==other.update_type\
                and self.type==other.type)


class ThinVariables:

    def __init__(self, parent):
        self.local = []
        self.input = []
        self.output = []
        self.parent = parent

    @property
    def names(self):
        names = []
        for var in self.all:
            names.append(var.name)
        return names
    
    @property
    def local_names(self):
        names = []
        for var in self.local:
            names.append(var.name)
        return names

    @property
    def output_names(self):
        names = []
        for var in self.output:
            names.append(var.name)
        return names

    @property
    def input_names(self):
        names = []
        for var in self.input:
            names.append(var.name)
        return names

    @property
    def all(self):
        return (self.input + self.local + self.output)

    def add_thinvar(self, v):

        if v.scope=='LOCAL_DATA':
            self.local.append(v)
        elif v.scope=='OUTPUT_DATA':
            self.output.append(v)
        elif v.scope=='INPUT_DATA':
            self.input.append(v)

        return


class ThinVariable:

    def __init__(self, name="default_thinvariable", update_type="", 
                 type=REAL, scope='LOCAL_DATA'):
        self.name = name
        self.type = type
        self.scope = scope
        self.update_type = type

    def __eq__(self, other):
        return (self.name==other.name\
                and self.update_type==other.update_type\
                and self.type==other.type)


class Mode:
    '''name - name of the mode (can be different from the one specified in 
              simulink/stateflow
    id - uniquely identifies mode. if mode created through HyIR, then id 
         corresponds to its index in the HyIR's modes list
    initial - True if this is the starting mode of the hybrid automata, 
              false otherwise
    invs - list of Invariant objects representing the mode's invariants
    dais - list of DAI objects representing the mode's governing differential 
           equations
    '''
    def __init__(self, name="default_mode", id=-1, initial=False):

        # Init the private "hidden" variables
        self._invariants = None
        self._dais = None
        # Call setters to ensure transition/dai parents are set
        self.invariants = []
        self.dais = []
        # Other variables, getters/setters not used
        self.name = name
        self.id = id
        self.initial = initial
        self.initialConditions = []
        self.linear = True
        self.parent = None

    @property
    def invariants(self):
        return self._invariants

    @invariants.setter
    def invariants(self, invs):
        for inv in invs:
            inv.parent = self
        self._invariants = invs
        return

    @property
    def dais(self):
        return self._dais

    @dais.setter
    def dais(self, ds):
        for d in ds:
            d.parent = self
        self._dais = ds
        return

    @property
    def invs(self):
        print("************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - MODE.INV")
        print("************************************************")
        return self.invariants

    @invs.setter
    def invs(self, invariants):
        print("************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - MODE.INV")
        print("************************************************")
        self.invariants = invariants
        return

    def add_invariant(self, inv):
        inv.parent = self
        self._invariants.append(inv)
        return

    def remove_invariant(self, inv):
        inv.parent = None
        self._invariants.remove(inv)
        return

    def clear_invariants(self):
        for inv in self._invariants:
            inv.parent = None
        self._invariants = []
        return

    def add_dai(self, dai):
        dai.parent = self
        self._dais.append(dai)
        return

    def remove_dai(self, dai):
        dai.parent = None
        self._dais.remove(dai)
        return

    def clear_dais(self):
        for dai in self._dais:
            dai.parent = None
        self._dais = []
        return

        
    def parse(self):
        """ Parse DAI equation and Invariant Equations """
        
        errors = []

        # Parse DAIs
        self.linear = True        
        for dai in self.dais:
            errors += dai.parse()
            if (self.linear and (dai.expr is not None)):
                self.linear = SymEq.is_linear(dai.expr.rhs)

        # Parse DAI variables
        errors += self.parse_dai_vars()

        # Parse Invariants
        for inv in self._invariants:
            p = inv.parse()
            errors += inv.parse()
            #if not inv.expr:
                #self.remove_invariants(inv)

        return errors

    def parse_dai_vars(self):
        """ 
        Parse variables used in DAI equations
        - Output variables on left-hand side must not end with _dot
        - Local variables on left-hand side must end with _dot
        - Input variables must not be used on the left-hand side
        
        Assumptions
        - <DAI>.expr contains only one symbol on the left-hand side
            - This check is done during <DAI>.parse()
        """

        errors = []
        local_vars = set(self.parent.local_var_names)
        output_vars = set(self.parent.output_var_names)
        input_vars = set(self.parent.input_var_names)

        for dai in self.dais:

            if dai.expr is None:
                continue

            lhs = str(dai.expr.lhs)
            var_ = lhs.replace('_dot', '')
            
            if var_ in local_vars:
                if not lhs.endswith('_dot'):
                    errors.append(('Flow', dai, "Incorrect Variable Usage",
                                   "Local variable equations must end with "
                                   + "_dot"))
                local_vars.remove(var_)
            elif var_ in output_vars:
                if lhs.endswith('_dot'):
                    errors.append(('Flow', dai, "Incorrect Variable Usage", 
                                   "Output variable equations must not end "
                                   + "with _dot"))
                output_vars.remove(var_)
            elif var_ in input_vars:
                errors.append(('Flow', dai, "Incorrect Variable Usage",
                               "Input variables should not appear on the "
                               "left-hand side of an equation"))

        for var_ in local_vars.union(output_vars):
            errors.append(('WARNING', self, "Unused variable: " + var_, None))
        
        return errors
        

class Transition:
    '''guard - node representing the guard for the transition
    actions - list of nodes representing the resets of the transition
    id - unique identifier for the transition
    src - the source of the transition
    dest - the destination of the transition
    ''' 
    def __init__(self, guard, actions, id=-1, source=-1, destination=-1):

        # Init the private "hidden" variables
        self._guard = None
        self._actions = None
        # Call setters to ensure guard/action parents are set
        self.guard = guard
        self.actions = actions
        # Other variables, getters/setters not used
        self.id = id
        self.source = source
        self.destination = destination
        self.parent = None

    @property
    def guard(self):
        return self._guard

    @guard.setter
    def guard(self, g):
        g.parent = self
        self._guard = g
        return

    @property
    def actions(self):
        return self._actions

    @actions.setter
    def actions(self, acts):
        for act in acts:
            act.parent = self
        self._actions = acts

    @property
    def name(self):
        return (self.parent.mode_dict[self.source] + " -> " 
                + self.parent.mode_dict[self.destination])

    def add_action(self, action):
        action.parent = self
        self._actions.append(action)
        return

    def clear_actions(self):
        for act in self._actions:
            act.parent = None
        self._actions = []
        return

    @property
    def src(self):
        print("********************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - TRANSITION.SOURCE")
        print("********************************************************")
        return self.source

    @src.setter
    def src(self, source):
        print("********************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - TRANSITION.SOURCE")
        print("********************************************************")
        self.source = source
        return

    @property
    def dest(self):
        print("*************************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - TRANSITION.DESTINATION")
        print("*************************************************************")
        return self.destination

    @dest.setter
    def dest(self, destination):
        print("*************************************************************")
        print(" WARNING: USING DEPRECATED PROPERTY - TRANSITION.DESTINATION")
        print("*************************************************************")
        self.destination = destination
        return

    def parse(self):

        errors = []
        errors += self.guard.parse()

        if not self.guard.expr:
            errors.append(('WARNING', self, "Transition will be deleted when "
                            + "composing. No guard expression.", None))
            
        for action in self.actions:
            errors += action.parse()

        return errors


class DAI:
    """Deterministic algebraic inequalities"""

    def __init__(self, raw=None):
        
        self.raw = raw
        self.expr = None
        self.parent = None

    @property
    def name(self):
        return self.raw
    
    def parse(self):

        errors = []
        if self.raw is None:
            errors.append(('Flow', self, "No Expression", None))
            self.expr = None
            return errors
        
        # Attempt to constructed a sympy equation
        constructed = SymEq.construct_eqn(self.raw, True, False)
        if constructed is None:  # construct_eqn returns None if it fails
            errors.append(('Flow', self, "Invalid Expression", None))
            self.expr = None
            return errors

        # Verify constructed equation has a single symbol on the left-hand side
        if type(constructed.lhs) is not sympy.Symbol:
            errors.append(('Flow', self, "Invalid Expession", "Left hand side "
                           + "of equation must be a single symbol"))
            self.expr = None
            return errors

        self.expr = constructed
        return errors


class Invariant:

    def __init__(self, raw=None):

        self.raw = raw
        self.expr = None
        self.parent = None

    @property
    def name(self):
        return self.raw

    def parse(self):
        
        errors = []
        if self.raw is None:
            errors.append(('Invariant', self, "No Expression", None))
            return errors

        eqns = self.raw.split('||')
        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, False, True)
            if constructed is None:
                errors.append(('Invariant', self, "Invalid Expression", eqn))
                return errors
            else:
                expr.append(constructed)
        self.expr = expr

        # Filter out equations that evaluate to False
        self.expr = filter(lambda eqn: eqn is not False, self.expr)
        self.expr = list(self.expr)
        if True in self.expr: 
            errors.append(('Invariant', self, "Redundant Invariant", self.raw))
            self.expr = []

        return errors


class Guard:

    def __init__(self, raw=None):

        self.raw = raw
        self.expr = None
        self.parent = None

    @property
    def name(self):
        return self.raw

    def parse(self):

        errors = []
        if self.raw is None:
            errors.append(('Guard', self, "No Expression", None)) 
            return errors

        eqns = self.raw.split('&&')
        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, False, True)
            if constructed is None:
                errors.append(('Guard', self, "Invalid Expression", eqn))
            else:
                expr.append(constructed)
        self.expr = expr

        # Filter out equations that evaluate to True
        self.expr = filter(lambda eqn: eqn is not True, self.expr)
        self.expr = list(self.expr)
        if False in self.expr: 
            errors.append(('Guard', self, "Redundant Guard", self.raw))
            self.expr = []
        
        return errors

        
class Action:

    def __init__(self, raw=None):

        self.raw = raw
        self.expr = None
        self.parent = None

    @property
    def name(self):
        return self.raw

    def parse(self):

        errors = []

        if self.raw is None:
            errors.append(('Action', self, "No Expression", None))
            return errors

        eqns = self.raw.split('&&')

        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, True, True)
            if constructed is None:
                errors.append(('Action', self, "Invalid Expression", eqn))
            else:
                expr.append(constructed)
        self.expr = expr

        return errors
import re
import sympy
import tkinter as tk
from scipy import optimize as opt
from tkinter import *
from tkinter.scrolledtext import *

from backend.lib.libc2e2 import Model
from frontend.mod.constants import *


class Property():

    def __init__(self):
        # Property values
        self.name = ''
        self.type = SAFETY
        self.time_step = 0.0
        self.time_horizon = 0.0
        self.k_value = 0.0
        self.initial_set_str = ''
        self.initial_set_obj = None
        self.unsafe_set_str = ''
        self.unsafe_set_obj = None
        self.simulator = ''

        # Property values valid
        self.name_valid = False
        self.type_valid = True
        self.time_step_valid = False
        self.time_horizon_valid = False
        self.k_value_valid = False
        self.initial_set_valid = False
        self.unsafe_set_valid = False

        # Misc
        self.status = "Not verified"
        self.result = ''
        self.is_visible = False

    def is_valid(self):
        
        return (self.name_valid and self.type_valid and self.time_step_valid \
                and self.time_horizon_valid and self.k_value_valid and \
                self.initial_set_valid and self.unsafe_set_valid)

    @staticmethod
    def validate_cur_prop():

        prop = Session.cur_prop

        Property.validate_name(prop.name)
        Property.validate_type(prop.type)
        Property.validate_time_step(prop.time_step)
        Property.validate_time_horizon(prop.time_horizon)
        Property.validate_k_value(prop.k_value)
        (is_status, is_error)= Property.validate_initial_set(prop.initial_set_str)
        Property.validate_unsafe_set(prop.unsafe_set_str)

        return

    @staticmethod
    def validate_name(name):
        
        valid = True
        if(name == ''):
            valid = False
        else:
            for prop in Session.hybrid.properties:
                if(prop == Session.cur_prop):
                    continue
                elif(name == prop.name):
                    valid = False
                    break

        Session.cur_prop.name_valid = valid
        return valid

    @staticmethod
    def validate_type(type):

        # TODO: LMB Not sure where this is ever invalidated (after being initialized to True). This function is a placeholder until I learn more.
        Session.cur_prop.type_valid = True
        return True
        
    @staticmethod
    def validate_time_step(time_step):
        
        #valid = (is_number(time_step) and (time_step >= 0))
        valid = (time_step >= 0)

        Session.cur_prop.time_step_valid = valid
        return valid

    @staticmethod
    def validate_time_horizon(horizon):
        
        #valid = (is_number(horizon) and (horizon >= 0))
        valid = (horizon >= 0)

        Session.cur_prop.time_horizon_valid = valid
        return valid

    @staticmethod
    def validate_k_value(k_value):

        #valid = (is_number(k_value) and (k_value >= 0))
        valid = (k_value >= 0)

        Session.cur_prop.k_value_valid = valid
        return valid

    @staticmethod
    def validate_initial_set(initial_set):
        
        # Regex building blocks
        flt = '(-?(\d*\.?\d+))'
        int = '(-?(\d+))'
        term = '(' + flt + '|(' + int + '\s*/\s*' + int + '))'
        var = '([a-zA-Z]\w*)'
        mode = '(' + var + ':)'
        eql = '((<=?)|(>=?)|(==))'
        expr = '(' + term + '?\s*' + var + '(\s*[-\+]?\s*((' + term + '?\s*' + var + ')|' + term + '))*)'
        eqn = '(' + expr + '\s*' + eql + '\s*' + term + ')'
        eqns = '(' + eqn + '(\s*&&\s*' + eqn + ')*)'
        init_set = '(\s*' + mode + '\s*' + eqns + '\s*){1}$'
        unsafe_set = '(\s*' + eqns + '\s*)+$'

        # Regex strings
        re_mode = mode
        re_var = var
        re_is = init_set
        re_us = unsafe_set

        match = re.match(re_is, initial_set)
        Session.cur_prop.initial_set_str = initial_set

        if(match == None):
            Session.cur_prop.initial_set_obj = None
            Session.cur_prop.initial_set_valid = False
            return (False, "Incorrect Syntax")
        else:
            is_sep = initial_set.split(':')

            # Validate Mode
            mode = re.search(re_var, is_sep[0]).group(0)
            mode_list = Session.hybrid.mode_names

            if(mode not in mode_list):
                Session.cur_prop.initial_set_obj = None
                Session.cur_prop.initial_set_valid = False
                return (False, "No matching modes")

            # Validate Vars
            vars_ = re.findall(re_var, is_sep[1])
            var_list = Session.hybrid.local_var_names
            var_union= set(vars_) | set(var_list)

            if(len(var_union) > len(var_list)):
                Session.cur_prop.initial_set_obj = None
                Session.cur_prop.initial_set_valid = False
                return (False, "Variable mismatch")

            # Parse equations
            a_m, b_m, eq_m = SymEq.get_eqn_matrix(is_sep[1], var_list)
            bounded = SymEq.check_boundedness(a_m, b_m, eq_m, var_list)

            if(is_sep[1].count('>')!= is_sep[1].count('<')):
                bounded = False
            
            if bounded:
                Session.cur_prop.initial_set_obj = [is_sep[0], a_m, b_m, eq_m]
                Session.cur_prop.initial_set_valid = True
                return (True, "")
            else:
                Session.cur_prop.initial_set_obj = None
                Session.cur_prop.initial_set_valid = False
                return (False, "Set unbounded")


    def validate_unsafe_set(unsafe_set):

        # Regex building blocks
        flt = '(-?(\d*\.?\d+))'
        int = '(-?(\d+))'
        term = '(' + flt + '|(' + int + '\s*/\s*' + int + '))'
        var = '([a-zA-Z]\w*)'
        mode = '(' + var + ':)'
        eql = '((<=?)|(>=?)|(==))'
        expr = '(' + term + '?\s*' + var + '(\s*[-\+]?\s*((' + term + '?\s*' + var + ')|' + term + '))*)'
        eqn = '(' + expr + '\s*' + eql + '\s*' + term + ')'
        eqns = '(' + eqn + '(\s*&&\s*' + eqn + ')*)'
        init_set = '(\s*' + mode + '\s*' + eqns + '\s*){1}$'
        unsafe_set_reg = '(\s*' + eqns + '\s*)+$'

        # Regex strings
        re_mode = mode
        re_var = var
        re_is = init_set
        re_us = unsafe_set_reg

        match = re.match(re_us, unsafe_set)
        Session.cur_prop.unsafe_set_str = unsafe_set

        # Check if input is valid
        if(match == None):
            Session.cur_prop.unsafe_set_obj = None
            Session.cur_prop.unsafe_set_valid = False
            return (False, "Incorrect Syntax")
        
        # Validate vars
        else:
            vars_ = re.findall(re_var, unsafe_set)
            var_list = Session.hybrid.local_var_names
            var_union = set(vars_) | set(var_list)

            if(len(var_union) > len(var_list)):
                Session.cur_prop.unsafe_set_obj = None
                Session.cur_prop.unsafe_set_valid = False
                return (False, "Variable mismatch")

            Session.cur_prop.unsafe_set_obj = SymEq.get_eqn_matrix(unsafe_set, var_list)
            Session.cur_prop.unsafe_set_valid = True
            return (True, "")


class PlotProperty():

    def __init__(self, identifier):

        self.identifier = identifier

        self._input_path = None
        self._input_path_full = None
        self._output_path = None
        self._output_path_full = None

        self.plot_name = ''
        self.plot_status = NOT_PLOTTED

        self.variable_list = []  # [String]
        self.mode_list = []  # [String]
        self.horizontal_index = ''  # Int
        self.vertical_indices = []  # Int
        self.vertical_select = []  # Int

        self.prop_name = ''
        self.prop_status = ''
        self.time_step = ''
        self.time_horizon = ''
        self.k_value = ''
        self.simulator = ''
        self.refine_strat = ''
        self.initial_set = ''
        self.unsafe_set = ''

    def load(data_dict):

        self.prop_name = data_dict[PROPNAME]
        self.prop_status = data_dict[SIMVER]
        self.mode_list = data_dict[MODENAMES]
        self.var_list = data_dict[VARIABLENAMES]    
        self.time_step = data_dict[TIMESTEP]
        self.time_horizon = data_dict[TIMEHORIZON]
        self.k_value = data_dict[KVALUE]
        self.simulator = data_dict[SIMULATOR]
        self.refine_strat = data_dict[REFINEMENTSTRAT]
        self.initial_set = data_dict[INITIALSET]
        self.unsafe_set = data_dict[UNSAFESET]
        
    @property
    def input_path_full(self):
        return self._input_path_full

    @input_path_full.setter
    def input_path_full(self, val):
        self._input_path_full = val
        self._input_path = self._input_path_full.split("/").pop()

    @property
    def input_path(self):
        return self._input_path

    @input_path.setter
    def input_path(self, val):
        self._input_path = val
        self._input_path_full = '../work-dir/output/' + self._input_path
        
    @property
    def output_path_full(self):
        return self._output_path_full

    @output_path_full.setter
    def output_path_full(self, val):
        self._output_path_full = val
        self._output_path = self._output_path_full.split("/").pop()

    @property
    def output_path(self):
        return self._output_path

    @output_path.setter
    def output_path(self, val):
        self._output_path = val
        self._output_path_full = '../work-dir/plotresult/' + self._output_path


class SymEq:
    """ Symbolic Equation Library """

    @staticmethod
    def construct_eqn(eqn, is_eq, rationalize):
        """ 
        Construct equation with sympy. Rationalize it if flag is set to True. is_eq is True for DAIs and Actions 
        
        NOTE: sympy.sympify(eqn) converts eqn into a type that can be used inside SymPy
        """

        try:
            if is_eq:
                eqn_split = eqn.split('=')
                lhs, rhs = eqn_split[0], eqn_split[1]
                #print (lhs, rhs)
                eqn = sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))
            else:
                eqn = sympy.sympify(eqn)
                if type(eqn) is bool:
                    return eqn
        except:
            Session.write("ERROR: Invalid expression.\n")
            return None

        if rationalize:
            eqn = SymEq.rationalize(eqn)
        return eqn

    def to_str(self):
        return str(self.eqn)

    # Given an equation convert all decimals to integers by multiplying by LCM
    @staticmethod
    def rationalize(expr):
        if expr.is_Relational:
            mult = max(SymEq.get_factor(expr.lhs), SymEq.get_factor(expr.rhs))
            return expr.func(mult*expr.lhs, mult*expr.rhs)

    @staticmethod
    def get_factor(expr):
        exp = 0
        if expr.is_Add:
            terms = expr.args
        else:
            terms = [expr]

        for term in terms:
            if term.is_Number:
                exp = max(exp, SymEq.get_term_exp(term))
            for unit in term.args:
                if unit.is_Number:
                    exp = max(exp, SymEq.get_term_exp(unit))

        return 10**exp

    @staticmethod
    def get_term_exp(unit):
        # print("****test****")
        # print(unit)
        text=str(float(unit))
        try:
            idx=text.index('e')
        except ValueError:
            val = len(str(float(unit))) - str((float(unit))).index('.') - 1
            # print(val)
            return val
        exp1 = abs(int(text[idx + 1:len(text)]))
        dec = text[0:idx]
        try:
            exp2 = len(dec)-dec.index('.')-1
        except ValueError:
            # print(exp1)
            return exp1
        val = exp1+exp2
        # print(val)
        return val
        
    @staticmethod
    def get_eqn_matrix(expressions, varList):
        """ Return A and B matrices for expression with x as varList """

        exprs = []
        for expr in expressions.split('&&'):
            if '==' in expr:
                eq = expr.split('==')
                exprs.append(sympy.sympify(eq[0]+'<='+eq[1]))
                exprs.append(sympy.sympify(eq[0]+'>='+eq[1]))
            else:
                sym_expr = sympy.sympify(expr)
                if sym_expr.func==sympy.StrictLessThan:
                    sym_expr = sympy.LessThan(sym_expr.lhs, sym_expr.rhs)
                elif sym_expr.func==sympy.StrictGreaterThan:
                    sym_expr = sympy.GreaterThan(sym_expr.lhs, sym_expr.rhs)
                exprs.append(sym_expr)
        aMatrix=[[0 for v in varList] for expr in exprs]
        bMatrix=[0 for expr in exprs]
        eqMatrix=[]
        for i, expr in enumerate(exprs):
            eqMatrix.append([expr.rel_op])
            expr = expr.lhs-expr.rhs
            expr *= SymEq.get_factor(expr)
            for v in expr.free_symbols:
                aMatrix[i][varList.index(str(v))] = expr.as_coefficients_dict()[v]
                bMatrix[i] = [-expr.as_coefficients_dict()[sympy.numbers.One()]]
        return [aMatrix, bMatrix, eqMatrix]

    # Return whether the initial set is bounded or not
    @staticmethod
    def check_boundedness(a_m, b_m, eq_m, var_list):

        A, B = [], []
        for r, a in enumerate(a_m):
            if eq_m[r][0] == '<=':
                A.append(a)
                B.append(b_m[r][0])
            else:
                A.append([-1 * x for x in a])
                B.append(-1 * b_m[r][0])

        # Solve minimum
        C = [0.0 for v in var_list]
        bounds = [(None, None) for v in var_list]
        for i, v in enumerate(var_list):
            C[i] = 1.0
            min_ = opt.linprog(C, A_ub=A, b_ub=B, bounds=bounds, 
                options=dict(tol=1e-11))  # Added tolerance to eliminate incorrect failures
            if not min_.success:
                return False

            C[i] = 0.0

        return True

    # Represent u**n as u*u*...u for the simulators.
    @staticmethod
    def convert_pow(expr):
        #SymEq.pow_to_mul(expr)
        return str(SymEq.convert_pow_helper(expr))

    @staticmethod
    def convert_pow_helper(expr):
        if not expr.args:
            return expr
        conv_args = [SymEq.convert_pow_helper(arg) for arg in expr.args]
        if expr.is_Pow and expr.exp>0:
            print("this expr is pow:", expr)
            print(expr.base, expr.exp)
            return sympy.Symbol('*'.join(['('+str(conv_args[0])+')']*int(conv_args[1])))
        return expr.func(*conv_args)

    @staticmethod
    def pow_to_mul(expr):
        """
        Convert integer powers in an expression to Muls, like a**2 => a*a.
        """
        pows = list(expr.atoms(sympy.Pow))
        if any(not e.is_Integer for b, e in (i.as_base_exp() for i in pows)):
            raise ValueError("A power contains a non-integer exponent")
        print(pows)
        for b,e in (i.as_base_exp() for i in pows):
            print(b,e)
        repl = zip(pows, (sympy.Mul(*[b]*e,evaluate=False) for b,e in (i.as_base_exp() for i in pows)))
        #print repl
        return expr.subs(repl)

    # Negate guard to construct invariant
    @staticmethod
    def construct_invariant(guard):
        inv_eqn = []
        for eqn in guard.expr:
            print('Eqn is: ' + str(eqn))
            if eqn.func==sympy.LessThan:
                inv = sympy.GreaterThan(*eqn.args)
            elif eqn.func==sympy.GreaterThan:
                inv = sympy.LessThan(*eqn.args)    
            elif eqn.func==sympy.StrictLessThan:
                inv = sympy.StrictGreaterThan(*eqn.args)    
            elif eqn.func==sympy.StrictGreaterThan:
                inv = sympy.StrictLessThan(*eqn.args)
            inv_eqn.append(inv)
        inv_eqn = [str(inv) for inv in inv_eqn]
        return Invariant('||'.join(inv_eqn))

    # Return all free vars in exprs
    @staticmethod
    def vars_used(exprs):
        var_set = set()
        for eqn in exprs:
            var_set = var_set.union(eqn.free_symbols)
        return [str(v) for v in var_set]

    # Return if expr is linear
    @staticmethod
    def is_linear(expr):
        syms = expr.free_symbols
        for x in syms:
            for y in syms:
                try:
                    if not sympy.Eq(sympy.diff(expr, x, y), 0):
                        return False
                except TypeError:
                    return False
        return True


class Feedback(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent, height=128)

        self.pack_propagate(0)

        # #EBEBEB = (235, 235, 235)
        # #DCDCDC = (220, 220, 220) = Disabled Entry widget bg color
        self.display = ScrolledText(self,state=DISABLED,bg='#EBEBEB',wrap=WORD)
        self.display.pack(fill=BOTH, side=TOP, expand=TRUE)


class StdRedirector(object):

    def __init__(self):
        self.widgets = []

    def add_feedback_frame(self, frame):
        self.widgets.append(frame)

    def write(self, string):
        # The Widgets are Feedback objects, <Feedback>.display = Text()
        for widget in self.widgets:
            widget.display.config(state=NORMAL)
            widget.display.insert(END, string)
            widget.display.see(END)
            widget.display.config(state=DISABLED)
        return

    def wait(self):
        for widget in self.widgets:
            widget.update()
    
    def flush(self):
        pass


class Session():

    file_opened = False
    file_modified = False
    file_saved = True
    file_path = ''
    file_type = None

    lib_compiled = False

    # Working Model
    hybrid = None
    cur_prop = None
    
    # Pre-composition (if hyrid was composed in this session)
    hybrid_precomp = None

    # Simulation parameters
    simulator = ODEINT_FIX
    refine_strat = DEF_STRAT
    linear_strat = GSTAR

    cpp_model = None

    # User Feedback
    writer = StdRedirector()
    window = None

    @classmethod
    def new_cpp_model(cls):
        cls.cpp_model = Model()
        return cls.cpp_model

    @classmethod
    def write(cls, string):
        cls.writer.write(string)
        print(string, end='')

    @classmethod
    def flush(cls):
        cls.writer.flush()

    @classmethod
    def add_feedback_frame(cls, parent):
        frame = Feedback(parent)
        cls.writer.add_feedback_frame(frame)
        return frame

    @classmethod
    def update(cls):
        if cls.window is not None:
            cls.window.update()
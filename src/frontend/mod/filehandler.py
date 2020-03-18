import os
import xml.dom.minidom
import xml.etree.ElementTree as ET
from tkinter import *
from tkinter.ttk import *

from frontend.mod.hyir import *
from frontend.mod.session import Session, Property
from frontend.mod.constants import * 
from frontend.mod.automaton import *
from frontend.gui.popupentry import PopupEntry


class SaveDialog(PopupEntry):

    def __init__(self, parent, hyxml_text=None):
        PopupEntry.__init__(self, parent)

        self.parent = parent
        self.hyxml_text = hyxml_text

        self._init_widgets()
        self.protocol('WM_DELETE_WINDOW', self._on_close)

    def _init_widgets(self):

        Label(self, text="Save changes before continuing?")\
            .grid(row=0, column=0)
        
        btn_row = Frame(self)
        btn_row.grid(row=1, column=0)
        Button(btn_row, text='Save', command=self._save_callback)\
            .pack(expand=TRUE, fill=X, side=LEFT)
        Button(btn_row, text='Save As', command=self._save_as_callback)\
            .pack(expand=TRUE, fill=X, side=LEFT)
        Button(btn_row, text='Discard', command=self._discard_callback)\
            .pack(expand=TRUE, fill=X, side=LEFT)

    def _save_callback(self):
        """ Save changes to model """
        FileHandler.save(self.hyxml_text)
        self.destroy()
    
    def _save_as_callback(self):
        """ Save changes to file dialog """
        FileHandler.save_as(self.hyxml_text)
        self.destroy()

    def _discard_callback(self):
        """ Reloaded file and return it to a saved state """
        Session.file_saved = True
        self.destroy()

    def _on_close(self):
        """ Discard changes if user closes the window """
        self._discard_callback()


class FileHandler:

    @staticmethod
    def save_as(hyxml_text=None):

        file_path = filedialog.asksaveasfile(**SAVE_OPT)
        if file_path:
            Session.file_path = file_path.name
            FileHandler.save(hyxml_text)

    @staticmethod
    def save(hyxml_text=None):

        if Session.file_path is None:
            Session.write("ERROR: SAVING WITH NO FILEPATH\n")

        if hyxml_text is None:
            FileHandler.save_tree()
        else:
            FileHandler.save_hyxml(hyxml_text)
            FileHandler.open_file(Session.file_path)

        Session.file_saved = True

        return

    @staticmethod
    def save_tree():

        if not Session.file_opened:
            return

        if Session.file_path is None:
            FileHandler.save_as()
            return

        if Session.file_type == MDL_FILE and not Session.file_saved:
            FileHandler.save_as()
            return

        FileHandler.save_model(
            Session.hybrid,
            Session.hybrid.properties,
            Session.file_path)
        
        Session.file_saved = True

        return

    @staticmethod
    def save_hyxml(hyxml_text):

        if not Session.file_path:
            Session.file_path = filedialog.asksaveasfile(**SAVE_OPT)

        if Session.file_path:
            with open(Session.file_path, 'w') as f:
                f.write(hyxml_text)
        
        Session.file_saved = True

        return

    @staticmethod
    def save_model(hybrid, property_list, file_path):

        Session.write("Saving filepath: " + str(file_path) + "\n")
        Session.write("Saving...")
        hyxml = ET.Element("hyxml", {"type":"Model"})

        for automaton in hybrid.automata:

            auto = ET.SubElement(hyxml, 'automaton', {"name": automaton.name})
            
            # Variables
            for var in automaton.vars:
                ET.SubElement(auto,"variable",{"name":var.name,"scope":var.scope,"type":var.type})

            # Thin Variables
            for thinvar in automaton.thinvars:
                ET.SubElement(auto, "thin_variable", { "name":thinvar.name, "scope":thinvar.scope, "type":thinvar.type })

            # Modes
            for mode in automaton.modes:
                m=ET.SubElement(auto,"mode",{"id":str(mode.id),"initial":str(mode.initial),"name":mode.name})
                for dai in mode.dais:
                    #equ = dai.raw[3:-1]
                    ET.SubElement(m, "dai", {"equation":dai.raw})
                for inv in mode.invariants:
                    ET.SubElement(m,"invariant",{"equation":inv.raw})
            
            # Transitions
            for tran in automaton.transitions:
                t = ET.SubElement(auto,"transition",{"id":str(tran.id),"destination":str(tran.destination),"source":str(tran.source)})
                ET.SubElement(t,"guard",{"equation":tran.guard.raw})
                for act in tran.actions:
                    ET.SubElement(t,"action",{"equation":act.raw})

        # Composition    
        #ET.SubElement(hyxml,"composition", {"automata":hybrid.automata[0].name})
            
        # Properties
        for prop in property_list:
            if not prop.initial_set_str:
                continue
            pt1 = ET.SubElement(hyxml,"property",{"name":prop.name,"type":str(prop.type),"initialSet":prop.initial_set_str,"unsafeSet":prop.unsafe_set_str})
            ET.SubElement(pt1,"parameters",{"timehorizon":str(prop.time_horizon),"timestep":str(prop.time_step),"kvalue":str(prop.k_value)})

        tree = ET.ElementTree(hyxml)
        
        def indent(elem,level=0):
            i = "\n" + level*"  "
            if len(elem):
                if not elem.text or not elem.text.strip():
                    elem.text = i + "  "
                if not elem.tail or not elem.tail.strip():
                    elem.tail = i
                for elem in elem:
                    indent(elem, level+1)
                if not elem.tail or not elem.tail.strip():
                    elem.tail = i
            else:
                if level and (not elem.tail or not elem.tail.strip()):
                    elem.tail = i

        indent(hyxml)
        tree.write(file_path)

        Session.write("  Done.\n")

        return    

    @staticmethod
    def open_file(file_path):
        """
        Loads a file and stores the model in Session.hybrid

        Args:
            file_path: The file path to load

        Returns:
            True/False status (True if file opened successfully)
        """

        Session.write("Opening File...\n")

        base_name = os.path.basename(file_path) 
        raw_name, ext = os.path.splitext(base_name)

        # Handle HyXML file
        if(ext == '.hyxml'):
            
            Session.file_type = HYXML_FILE

            # Get HyXML type and call corresponding function
            hyxml_tree = ET.parse(file_path)
            if(hyxml_tree == None):
                return False

            hyxml_root = hyxml_tree.getroot()
            hyxml_type = hyxml_root.get('type')

            thinvarprop = ""
            thinvarlist = ""

            if(hyxml_type == 'Model'):

                automata = FileHandler.open_hyxml_model(hyxml_root)
                properties = FileHandler.open_hyxml_properties(hyxml_root)

            else:
                return False

        # Handle MDL file
        elif(ext == '.mdl'):

            Session.file_type = MLD_FILE

            hybrid = FileHandler.open_mdl_model(file_path, raw_name)
            prop_list = []

        # Handle all other extensions
        else:
            return False
        
        hybrid = HyIR(file_name=file_path)

        hybrid.automata = automata
        hybrid.properties = properties

        Session.hybrid = hybrid
        Session.cur_prop = properties[0]

        Session.write("File Opened.\n")
        return True

    # Open HyXML Model file
    @staticmethod
    def open_mdl_model(file_path,file_name):
        print ("start to parse mdl file and construct hyir")
        model=open(file_path,"r")
        rawModel = model.read()
        x = re.search(r"Stateflow {",rawModel)
        sf_data=rawModel[x.start():]
        sf_tree= extract_sf(sf_data)
        if IsHierarchical(sf_tree):
            sf_tree=RemoveHierarchy(sf_tree)
        hybrid = HyIR(file_name = file_name)
        automaton = Automaton()
        hybrid.automata = [automaton]

        mapIDs = {}
        f_tran_list = []
        for i in sf_tree.children[1].children:
            if i.type == "SFState":
                m = Mode()
                acceptableMode = True
                f_tran_id = 0
                for j in i.children:
                    if j.type == "type" and j.value == "AND_STATE":
                        acceptableMode = False
                    if j.type == "id":
                        m.id = automaton.new_mode_id()
                        mapIDs[j.value] = m.id
                    if j.type == 'labelString':
                        m.name = j.children[0].value
                        for k in j.children[2].children:
                            rate = str(k.children[1].value)
                            if not re.search(r'^begin', rate) is None and \
                            not re.search(r'^end$',rate) is None:
                                rate = rate[5:-3]
                                rate = re.sub('_',',',rate)
                                rate = re.sub('neg', lambda x: '-', rate)
                                k.children[1].value = "[%s]" % rate
                                k.value = " in "
                            m.add_dai(DAI(collapse(k)))
                    if j.type == "firstTransition":
                        f_tran_id = j.value
                    if j.type == "position":
                        tmp3 = j.value
                        xp = int((tmp3[0]+tmp3[2]))
                        yp = int((tmp3[1]+tmp3[3]))
                        m.xpos = xp
                        m.ypos = yp
                if acceptableMode:
                    hybrid.automata[0].add_mode(m)
                else:
                    hybrid.automata.append(Automaton(m.name))
                    automaton.next_mode_id -= 1
                    f_tran_list.append(f_tran_id)
            
            if i.type == "SFAnnot":
                for j in i.children:
                    if j.type == "labelString":
                        bufferString = j.value.replace("@","")
                        hybrid.annotationsRaw += [bufferString]
                        hybrid.annotations += bufferString + "\n"

            if i.type == "SFData":
                v = Variable()
                for j in i.children:
                    if j.type == "name":
                        v.name = j.value
                    if j.type == "scope":
                        v.scope = j.value
                    if j.type == "props": 
                        if j.children[0].value == "SF_CONTINUOUS_TIME_DATA":
                            v.update_type = "Continuous"
                            v.type = "Real"
                        else:
                            v.update_type = "Discrete"
                            v.type = ""
                if v.scope != "INPUT_DATA":
                    hybrid.add_var(v)

            if i.type == "SFTransition":
                actions = []
                initialFlag = False
                for j in i.children:
                    if j.type == "id":
                        tid = automaton.new_transition_id()
                        sl_id = j.value #used for determining initial modes
                    elif j.type == "Source Block":
                        if j.children[0].value != 'None' and j.children[0].value !='':
                            srid = mapIDs[j.children[0].value]
                        else:
                            initialFlag = True
                    elif j.type == "Dest Block":
                            dsid = mapIDs[j.children[0].value]                           
                    elif j.type == "labelString":
                        for k in j.children:
                            if k.type == 'Logical' or k.type == 'Relational':
                                guard = Guard(collapse2(k))
                            elif k.type == 'Assignment':
                                actions.append(Action(collapseAction(k)))
                if not initialFlag: 
                    hybrid.automata[0].modes[srid].add_inv(SymEq.construct_invariant(guard))
                    t = Transition(guard,actions,tid,srid,dsid)
                    hybrid.automata[0].add_trans(t)
                else:
                    hybrid.automata[0].modes[dsid].initial = True
                    for action in actions:
                        hybrid.automata[0].modes[dsid].initialConditions.append(re.sub(r"'", lambda x: '', action.raw)) 
                    initialFlag = False

                    for x in range(0, len(f_tran_list)):
                        if sl_id == f_tran_list[x]:
                            hybrid.automata[x+1].initial_mode_id = dsid

        separateAutomata(hybrid)
        hybrid.populateInvGuards()

        return hybrid
        
    @staticmethod
    def open_hyxml_model(root):
        """
        Loads automata from input xml root

        Args:
            root: XML root

        Returns:
            List of Automaton() objects
        """

        Session.write("  Loading hyxml model...")

        automata = []

        for auto in root.iterfind("automaton"):

            name = auto.get("name")
            automaton = Automaton(name)
            
            for var in auto.iterfind("variable"):
                
                # Load variables
                v_name = var.get("name")
                v_scope = var.get("scope")
                v_type = var.get("type")
               
                v = Variable(name=v_name, type=v_type, scope=v_scope)
                automaton.add_var(v)

            for thinvar in auto.iterfind("thin_variable"):
                
                # Load thin variables
                v_name = thinvar.get("name")
                v_scope = thinvar.get("scope")
                v_type = thinvar.get("type")
                
                v = ThinVariable(name=v_name, type=v_type, scope=v_scope)
                automaton.add_thinvar(v)

            for mode in auto.iterfind("mode"):
                
                # Load modes
                mode_name = mode.get("name")
                mode_id = int(mode.get("id"))
                mode_init = (mode.get("initial") == "True")
                
                # if automaton.next_mode_id <= mode_id:
                #     automaton.next_mode_id = mode_id + 1
                    
                mode_obj = Mode(name=mode_name, id=mode_id, initial=mode_init)
        
                for dai in mode.iterfind("dai"):

                    # Load Flows 
                    raw_eq = dai.get("equation")                    
                    mode_obj.add_dai(DAI(raw_eq))

                for inv in mode.iterfind("invariant"):
                    
                    # Load Invariants
                    raw_eq = inv.get("equation")
                    # Equation 'cleaning' is needed for inequalities
                    clean_eq = FileHandler.clean_eq(raw_eq)
                    mode_obj.add_invariant(Invariant(clean_eq))
                
                automaton.add_mode(mode_obj, mode_id)
            
            if not automaton.verify_mode_ids():
                Session.write("  FILE READ ERROR: MODE IDS NOT UNIQUE\n")
                Session.write("  Automaton: " + automaton.name + "\n")
            if not automaton.verify_mode_names():
                Session.write("  FILE READ ERROR: MODE NAMES NOT UNIQUE\n")
                Session.write("  Automaton: " + automaton.name + "\n")

            for tran in auto.iterfind("transition"):

                # Load transitions
                g = tran.find("guard")
                guard = Guard(FileHandler.clean_eq(g.get("equation")))

                tran_id = int(tran.get("id"))
                tran_src = int(tran.get("source"))
                tran_dest = int(tran.get("destination"))
                
                # Actions
                actions = []
                for act in tran.iterfind("action"):

                    raw_eq = act.get("equation")
                    clean_eq = FileHandler.clean_eq(raw_eq)
                    actions.append(Action(clean_eq))

                transition = Transition(guard, actions, tran_id, 
                                        tran_src, tran_dest)
                automaton.add_transition(transition)

            if not automaton.verify_transition_src_dest():
                Session.write("  FILE READ ERROR: TRANSITION SOURCE/DESTINATION IDS " +
                      "NOT VALID MODE IDS\n")
                Session.write("  Automaton: " + automaton.name +"\n")

            automata.append(automaton)

        Session.write("  Done.\n")
        return automata

    @staticmethod
    def open_hyxml_properties(root):
        """ 
        Load properties from hyxml
        
        Args:
            root: XML root

        Returns:
            List of Proerty() objects
        """
        
        Session.write("  Loading hyxml properties...")

        prop_list = []
        for prop in root.iterfind('property'):
            
            p = Property()
            p.name = prop.get('name')

            p.type = SAFETY 
            p.initial_set_str = FileHandler.clean_eq(prop.get('initialSet'))
            p.unsafe_set_str = FileHandler.clean_eq(prop.get('unsafeSet'))

            # Handle properties parameters
            param = prop.find('parameters')
            if param is not None:
                time_step = param.get('timestep')
                if time_step is None:
                    p.time_step = 0.0
                else:
                    p.time_step = float(time_step)

                time_horizon = param.get('timehorizon')
                if time_horizon is None:
                    p.time_horizon = 0.0
                else:
                    p.time_horizon = float(time_horizon)

                k_value = param.get('kvalue')
                if k_value is None:
                    p.k_value = 0.0
                else:
                    p.k_value = float(k_value)

            prop_list.append(p)

        Session.write("  Done.\n")
        return prop_list

    @staticmethod
    def prepend(input, output, lines):
        with open(input, 'r+') as f:
            content = f.read()
        with open(output, 'w+') as f:
            f.seek(0, 0)
            f.write(lines.rstrip('\r\n') + '\n' + content)

    @staticmethod
    def prepend_cur_prop(simver):
        """ Prepend current property values to output file """
        
        # Property Name
        lines = Session.cur_prop.name + '\n'
        # Simulated or Verified
        lines += simver + '\n'
        # Modes Names
        for mode in Session.hybrid.mode_names:
            lines += mode + ' '
        lines += '\n'
        # Variable Names
        lines += 'time '  # Time is always first variable
        for var_ in Session.hybrid.local_var_names:
            lines += var_ + ' '
        lines += '\n'
        # Time Step
        lines += str(Session.cur_prop.time_step) + '\n'
        # Time Horizon
        lines += str(Session.cur_prop.time_horizon) + '\n'
        # K Value
        lines += str(Session.cur_prop.k_value) + '\n'
        # Simulator
        lines += Session.simulator + '\n'
        # Refinement Strategy
        lines += Session.refine_strat + '\n'
        # Initial Set
        lines += Session.cur_prop.initial_set_str + '\n'
        # Unsafe Set
        lines += Session.cur_prop.unsafe_set_str + '\n'

        FileHandler.prepend('../work-dir/' + Session.cur_prop.name, 
            '../work-dir/output/' + Session.cur_prop.name, lines)

    @staticmethod
    def parse_data_file(file_path):
        """ Parse data file, must be in C2E2 format """

        file_object = open(file_path, 'r')
        lines = file_object.readlines()
        data = {}

        # Property Name
        data[PROPNAME] = lines[0].strip()
        # Simulated or Verified
        data[SIMVER] = lines[1].strip()
        # Mode Names
        data[MODENAMES] = lines[2].split()
        # Variable Names
        data[VARIABLENAMES] = lines[3].split()
        # Time Step
        data[TIMESTEP] = float(lines[4].strip())
        # Time Horizon
        data[TIMEHORIZON] = float(lines[5].strip())
        # K Value
        data[KVALUE] = float(lines[6].strip())
        # Simulator
        data[SIMULATOR] = lines[7].strip()
        # Refinement Strategy
        data[REFINEMENTSTRAT] = lines[8].strip()
        # Initial Set
        data[INITIALSET] = lines[9].strip()
        # Unsafe Set
        data[UNSAFESET] = lines[10].strip()

        return data
   
    @staticmethod   
    def clean_eq(eq):
        r_dict = {'&lt;':'<', '&gt;':'>', '&amp;':'&', ' and ':'&&', ' or ':'||'}
        for term in r_dict:
            eq = eq.replace(term, r_dict[term])
        return eq

from tkinter import *
from tkinter.ttk import *

from frontend.gui.eventhandler import EventHandler
from frontend.gui.modeltab import ModelTab
from frontend.mod.filehandler import *
from frontend.mod.constants import *
from frontend.mod.session import Session
from frontend.gui.editortab import EditorTab 
from frontend.gui.plottab import PlotTab


class ModelNotebook(Notebook):

    def __init__(self, parent, **options):
        Notebook.__init__(self, parent, **options)

        self.parent = parent
        self.previous_tab = 'None'
        self.current_tab =  'None'
        self.plot_tab_dic = {}

        self._init_widgets()
        self._bind_events()

    def display_widgets(self, event=None):

        self.pack(fill=BOTH, expand=True)

        return

    def hide_widgets(self,event=None):

        keys = list(self.plot_tab_dic.keys())
        for key in keys:
            self.forget(self.plot_tab_dic[key])
            self.plot_tab_dic.pop(key,None)
        self.pack_forget()

        return

    def _bind_events(self):

        self.bind(OPEN_EVENT, self.display_widgets)
        EventHandler.add_event_listeners(self, OPEN_EVENT)
        self.bind(CLOSE_EVENT,self.hide_widgets)
        EventHandler.add_event_listeners(self, CLOSE_EVENT)
        self.bind('<<NotebookTabChanged>>', self._tab_changed)

        return

    def _init_widgets(self):

        self.model_tab = ModelTab(self)
        self.add(self.model_tab, text=MODEL)
        self.editor_tab = EditorTab(self)
        self.add(self.editor_tab, text=EDITOR)
        self.plot_tab = PlotTab(self)
        self.add(self.plot_tab, text=PLOT)

        return

    def _init_plot_widgets(self, *args):

        name = args[-1] + " plot"
        if name in self.plot_tab_dic:
            self.forget(self.plot_tab_dic[name])

        plot_model_tab = PlotterModelTab(self, *args)
        self.plot_tab_dic[name] = plot_model_tab
        self.add(plot_model_tab, text =name)

        return

    def _close_plot_tab(self, tabname):
        
        self.forget(self.plot_tab_dic[tabname])
        self.plot_tab_dic.pop(tabname,None)

        return

    def _tab_changed(self, event=None):
        
        self.previous_tab = self.current_tab
        self.current_tab = self.tab(self.select(), 'text')

        if not Session.file_saved:
            if self.previous_tab == EDITOR:
                hyxml_text = self.editor_tab.editor.get('1.0', 'end-1c')
            else:
                hyxml_text = None
            save_dialog = SaveDialog(self, hyxml_text)
            self.wait_window(save_dialog)

            FileHandler.open_file(Session.hybrid.file_name)

            # Refresh Edit Tab
            self.editor_tab.open_xml()
            # Refresh Model Tab
            self.model_tab.tree._clear_model()
            self.model_tab.tree._display_model()

        return
import os
import platform
import time
import tkinter as tk
import tkinter.ttk as ttk
import tkinter.scrolledtext as tkst
import PIL.Image
import PIL.ImageTk
import webbrowser

from frontend.gui.eventhandler import EventHandler
from frontend.gui.menubar import MenuBar
from frontend.gui.modelnotebook import ModelNotebook
from frontend.mod.constants import *
from frontend.mod.session import Session


class C2E2(tk.Tk):

    def __init__(self, parent):
        tk.Tk.__init__(self, parent)

        self.title("C2E2 v" + VERSION)
        
        width = self.winfo_screenwidth()
        height = self.winfo_screenheight()
        appwidth = (width * 2) // 3
        appheight = (height * 2) // 3
        appheight = 650 if (appheight < 650) else appheight
        offset_x = (width // 2) - (appwidth // 2)
        offset_y = (height // 2) - (appheight // 2)

        self.geometry('{}x{}+{}+{}'.format(appwidth, appheight, offset_x, offset_y))
        self.resizable(width=True, height=False)

        self.open_label = tk.Label(self,text="Open a new model (File -> Open)")
        self.open_label.pack(fill=tk.Y)
        self.open_label.config(font=('Courier', 20))

        self.manual_label = tk.Label(self, text="Read C2E2 manual (Click Me!)")
        self.manual_label.pack(fill=tk.Y)
        self.manual_label.config(font=('Courier', 20))
        self.manual_label.bind('<Button-1>', self._show_manual)
        
        self.email_label = tk.Label(self, text="Email c2e2help@gmail.com for support")
        self.email_label.pack(fill=tk.Y)
        self.email_label.config(font=('Courier', 20))

        self._init_widgets()
        self._bind_events()
    
    def _show_manual(self,event):
        webbrowser.open_new("http://publish.illinois.edu/c2e2-tool/files/2019/05/C2E2Manual_ver0.pdf")

    def _init_widgets(self):

        Session.window = self
        self.notebook = ModelNotebook(self)
        self.menu = MenuBar(self, self.notebook)

    def _bind_events(self):

        self.bind(CLOSE_EVENT, self._reset_name)
        EventHandler.add_event_listeners(self, CLOSE_EVENT)

        self.bind(OPEN_EVENT, self._set_name)
        EventHandler.add_event_listeners(self, OPEN_EVENT)

    def _set_name(self, event):

        self.title(Session.file_path.split("/").pop() + " - C2E2 v" + VERSION)

    def _reset_name(self, event):
        
        self.title("C2E2 v" + VERSION)
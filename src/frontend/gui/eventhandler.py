from tkinter import *
from tkinter.ttk import *

from frontend.mod.constants import *

class EventHandler():
    event_dict = {}

    @staticmethod
    def add_event_listeners(listeners, event):
        if not isinstance(listeners, list):
            listeners = [listeners]
        
        if event in EventHandler.event_dict.keys():
            EventHandler.event_dict[event] += listeners
        else:
            EventHandler.event_dict[event] = listeners

    @staticmethod
    def event_generate(event):
        #print(EventHandler.event_dict[event])
        for listener in EventHandler.event_dict[event]:
            listener.event_generate(event, when='tail')
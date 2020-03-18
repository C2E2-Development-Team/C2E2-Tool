from tkinter import *
from tkinter.ttk import *
import re
from sys import platform
from frontend.mod.constants import *
from frontend.mod.session import Session, Property
import PIL.Image
import PIL.ImageTk


class ValidLabel(Label):

    def __init__(self, parent, ses_var=None, **options):
        Label.__init__(self, parent, **options)
        
        self.parent = parent
        if platform == 'linux':
            self.inv_img = PhotoImage(file=CROSS_IMG)
            self.val_img = PhotoImage(file=TICK_IMG)
        else:
            self.inv_img = PIL.Image.open(CROSS_IMG).convert('RGB')
            self.val_img = PIL.Image.open(TICK_IMG).convert('RGB')

            self.inv_img = PIL.ImageTk.PhotoImage(self.inv_img)
            self.val_img = PIL.ImageTk.PhotoImage(self.val_img)

        self.config(image=self.inv_img)

    def set_state(self, state):

        if state == True:
            self.config(image=self.val_img)
        else:
            self.config(image=self.inv_img)

class SetText(Frame):
    
    def __init__(self, parent, width=0, height=0, callback=None, **options):
        Frame.__init__(self, parent, width=width, height=height)
        self.parent = parent
        self.callback = callback

        self.text = Text(self, **options)
        self.text.bind('<KeyRelease>', self._callback)
        self.text.pack(expand=TRUE, fill=Y)

    def _callback(self, event=None):
        txt = self.get()
        if self.callback:
            self.callback(txt)

    def pack(self, *args, **kwargs):
        Frame.pack(self, *args, **kwargs)
        self.pack_propagate(False)

    def grid(self, *args, **kwargs):
        Frame.grid(self, *args, **kwargs)
        self.pack_propagate(False)

    def get(self, end=END):
        txt = self.text.get('1.0', end).strip()
        return txt

    def set(self, str_):
        self.text.delete('1.0', END)
        self.text.insert(END, str_)
        self._callback()

    def delete(self):
        self.text.delete('1.0', END)

    def insert(self, str_):
        self.text.insert(INSERT, str_)

class SimpleToggleFrame(Frame):

    def __init__(self, parent, text):
        Frame.__init__(self, parent)

        # Title and Collapse Button

        self.header_frame = Frame(self)
        Label(self.header_frame, text=text).pack(side=LEFT, fill=X,expand=True)

        self.visible = BooleanVar()
        self.visible.set(False)
        self.toggle_btn = Checkbutton(self.header_frame, text="+", width=2,
            variable=self.visible, command=self._btn_toggle,style='Toolbutton')
        self.toggle_btn.pack(side=LEFT)

        self.header_frame.pack(fill=X, expand=True)

        self.body_frame = Frame(self)

    def _btn_toggle(self):

        if self.visible.get():
            self.body_frame.pack(fill=X, expand=True)
            self.toggle_btn.config(text="-")
        else:
            self.body_frame.forget()
            self.toggle_btn.config(text="+")


class ToggleFrame(Frame):

    def __init__(self, parent, text):
        Frame.__init__(self, parent)
              
        # Title / Collapse button
        
        self.header_frame = Frame(self)
        
        Label(self.header_frame, text=text).pack(side=LEFT, fill=X, expand=TRUE)
        
        self.visible = BooleanVar()
        self.visible.set(False)
        self.toggle_btn = Checkbutton(self.header_frame, text='+', width=2, 
            variable=self.visible, command=self._btn_toggle,style='Toolbutton')
        self.toggle_btn.pack(side=LEFT)

        self.header_frame.pack(fill=X, expand=TRUE)
    
        # Content (Entry Fields)

        self.body_frame = Frame(self)

        self.rows = []  # StringVar()
        self.entry_fields = []  # Entry fields 
        self.row_index = 0
               
        # Add Button

        self.btn_frame = Frame(self.body_frame)
        Button(self.btn_frame, text='Add Row', command=self.add_row).pack()
        
    def _btn_toggle(self):
        """ Toggle visibility when user clicks button """

        if(self.visible.get()):
            self.body_frame.pack(fill=X, expand=TRUE)
            self.toggle_btn.config(text='-')
        else:
            self.body_frame.forget()
            self.toggle_btn.config(text='+')


    def toggle(self):
        """ Toggle visibility (simulates a button click) """

        if(self.visible.get()):
            self.visible.set(False)
        else:
            self.visible.set(True)

        self._btn_toggle()       


    def add_row(self, text=''):
        """ Add entry row, optionally filled in with text """
        
        self.rows.append(StringVar())
        self.btn_frame.forget()

        if(text):
            self.rows[self.row_index].set(text)

        self.entry_fields.append(Entry(self.body_frame, textvariable=self.rows[self.row_index]))
        self.entry_fields[self.row_index].pack(fill=X, expand=TRUE)

        self.btn_frame.pack()

        self.row_index += 1

    def get_rows(self):
        return self.rows
    
    def disable_fields(self):
        """ Add readonly display row, must have text """

        self.btn_frame.forget()
        for entry in self.entry_fields:
            entry.config(state=DISABLED)
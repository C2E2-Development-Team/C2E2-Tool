import filecmp
import os
import platform
import subprocess
import sys
import time
import tkinter as tk
import tkinter.ttk as ttk
import tkinter.messagebox as tkmb
import tkinter.scrolledtext as tkst
import PIL.Image
import PIL.ImageTk

from frontend.mod.constants import *


class TestC2E2(tk.Tk):
    
    def __init__(self, parent, test=None):
        tk.Tk.__init__(self, parent)

        self.test = test
        self._init_widgets()

        if test is None:
            self._display_button_frame()
        else:
            switch_dict = {SHORT_TEST : self._short_test,
                           LONG_TEST : self._long_test,
                           DEVELOPER_TEST : self._developer_test,
                           STARTUP_TEST : self._startup_test} 
            test_result = switch_dict[test]()
            
    def _init_widgets(self):

        self.title('Test - C2E2 v' + VERSION)

        # Logo
        img = PIL.ImageTk.PhotoImage(PIL.Image.open('./frontend/res/logo.png'))
        img_width = 543
        img_height = 248
        self.logo = tk.Label(self, image=img)
        self.logo.image = img
        self.logo.grid(row=0, column=0)

        # Progress Bar
        self.progressbar = ttk.Progressbar(self, orient=tk.HORIZONTAL,
            mode='determinate')
        self.progressbar.grid(row=1, column=0, sticky=tk.NSEW)
        
        # Feedback Frame
        self.feedback_frame = tk.Frame(self)
        self.feedback_frame.grid(row=2, column=0, sticky=tk.NSEW)

        # Feedback Window
        self.feedback = tkst.ScrolledText(self.feedback_frame,
            state=tk.DISABLED, bg='#EBEBEB', wrap=tk.WORD, width=100)
        self.feedback.pack()

        # Butttons
        self.buttons = tk.Frame(self.feedback_frame)
        tk.Button(self.buttons, text="  Short  ", command=self._short_test)\
            .pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
        tk.Button(self.buttons, text="  Long   ", command=self._long_test)\
            .pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
        tk.Button(self.buttons, text="Developer", command=self._developer_test)\
            .pack(side=tk.LEFT, expand=True, fill=tk.BOTH)    

        # Geometry
        width = self.winfo_screenwidth()
        height = self.winfo_screenheight()
        img_width = 543
        img_height = 248
        img_x = (width//2) - (img_width//2)
        img_y = (height//2) - (img_height//2)

        self.logo.grid(row=0, column=0)
        self.progressbar.grid(row=1, column=0, sticky=tk.NSEW)
        self.feedback_frame.grid(row=2, column=0, sticky=tk.NSEW)
        self.feedback_frame.pack_propagate(False)
        self.feedback_frame.config(height=66, width=img_width)

        self.geometry('{}x{}+{}+{}'.format(
            img_width, img_height+88, img_x, img_y))
        self.resizable(width=False, height=False)

        self.update()

    def _display_button_frame(self):

        self.feedback.pack_forget()
        self.buttons.pack()

    def _display_feedback_frame(self):

        self.buttons.pack_forget()
        self.feedback.pack()       

    def _startup_test(self):
        """ Run short tests, then automatically close. """
        
        self._short_test()
        time.sleep(2)
        self.destroy()

    def _short_test(self):
        """" Startup tests. Check OS, libraries, files, and folders. """

        self._display_feedback_frame()

        required_dirs = ["../work-dir/",
                         "../work-dir/output/",
                         "../work-dir/plotresult/"]
        required_files = ["backend/lib/libc2e2.so"]

        result = 1

        system = platform.system()
        version = platform.version()
        # system = "Windows"
        # version = "Debian 9.5"

        # Soft stops for system (Linux) and version (Ubuntu 16.04)
        self.write("Checking system... ")
        if system != 'Linux':
            message = "C2E2 currently only supports Ubuntu 16.04. " + system +\
            " was detected as teh current operating system. You may " +\
            "experience unintended behavior."
            tkmb.showwarning("Unsupported System", message)
            self.write("\nWARNING: Unsupported System\n")
            self.write(message + "\n")
        else:
            self.write("Done!\n")
        
        self.write("Checking system version... ")
        if ('16.04' not in version) or ('Ubuntu' not in version):
            message = "An unsupported version of Linux was detected. C2E2 " +\
            "was built and tested using Ubuntu 16.04. You may experience " +\
            "unintended behavior."
            tkmb.showwarning("Unsupported System Version", message)
            self.write("\nWARNING: Unsupported System Version\n")
            self.write(message)
        else:
            self.write("Done!\n")

        # Hard stops for missing files and directories
        self.write("Checking required files and directories...")
        for dir_ in required_dirs:
            if not os.path.isdir(dir_):
                message = "Required path does not exist."
                message += "Missing path: " + dir_
                tkmb.showerror("Fatal Error", message)
                self.write("\nERROR: " + message + "\n")
                result = 0

        for file_ in required_files:
            if not os.path.exists(file_):
                message = "Required file does not exist."
                message += "Missing files: " + file_
                tkmb.showerror("Fatal Error", message)
                self.write("\nERROR: " + message + "\n")
                result = 0
            
        self.write("Test Complete.\n")
        self.write("-----------------\n")
        if result:
            self.write("All tests passed!\n")
        else:
            self.write("One or more failures occured. See testlog.txt")

        with open('../test-dir/TestResult.dat', 'w+') as f:
                f.write(str(result))

        return 

    def _long_test(self):
        """ Install tests. Run a few examples and compare to expected results. """

        self._display_feedback_frame ()

        test_directory = '../test-dir/'
        filename_list = ['TotalMotion40s', 'cardiacCell', 'linThermo', 'TotalMotionV2', 'bruss']
        sim_result_list = ['1', '1', '-1', '1', '1']
        ver_result_list = ['1', '1', '-1', '1', '1']

        return self._run_examples(test_directory, filename_list, 
            sim_result_list, ver_result_list)

    def _developer_test(self):
        """ Developer tests. Tests and direction compare output files against expected. """

        self._display_feedback_frame()

        filename_list = ['bruss', 'cardiacCell', 'TotalMotionV2', 
            '7d_model_flowstar', 'hybrid_inverter_ramp']    


        result = 1
        for i, filename in enumerate(filename_list):
            for action in ['simulate', 'verify']:
                self.write(filename + " " + action + " test starting... ")
                sts = subprocess.Popen('python3 main.py -f ' + \
                    ('../test-dir/' + filename + '.hyxml') + ' -a ' + action, 
                    shell=True).wait()
                test_file = '../work-dir/' + filename + '_test'
                control_file = '../test-dir/' + filename + '_' + action +\
                    '_control'
                if filecmp.cmp(test_file, control_file):
                    self.write("Passed!\n")
                else:
                    self.write("FAILED\n")
                    self.write("Verification results do not match with expected results.\n")
                    result = 0

        with open('../test-dir/TestResult.dat', 'w+') as f:
            f.write(str(result))

        self.write("Tests Complete.")
        if result:
            self.write(" All tests passed!\n")
            #time.sleep(2)
        else:
            self.write("\nOne or more failures occured. See testlog.txt\n")
            #time.sleep(5)

    def _run_examples(self, directory, filenames, sim_results, ver_results):

        progressbar_step = 100 / (len(sim_results) + len(ver_results))
        result = 1

        for i, filename in enumerate(filenames):

            self.write(filename + " simulation test starting... ")
            sts = subprocess.Popen("cd ../src/; python3 main.py -f " + \
                (directory + filename + '.hyxml ') + '-a simulate', shell=True)\
                .wait()
            with open('../work-dir/Result.dat', 'r') as f:
                result = f.readline().strip()
                if result == sim_results[i]:
                    self.write("Passed!\n")
                else:
                    self.write("FAILED\n")
                    self.write("Expected: " + sim_result_list[i] + "\n")
                    self.write("result: " + result + "\n")
                    result = 0
            
            self.progressbar.step(progressbar_step)
            self.update()

            self.write(filename + " verification test starting...")
            sts = subprocess.Popen("cd ../src/; python3 main.py -f " + \
                (directory + filename + '.hyxml ') + '-a simulate', shell=True)\
                .wait()
            with open('../work-dir/Result.dat', 'r') as f:
                result = f.readline().strip()
                if result == ver_results[i]:
                    self.write("Passed!\n")
                else:
                    self.write("FAILED\n")
                    self.write("Expected: " + ver_results[i] + "\n")
                    self.write("Result: " + result + "\n")
                    result = 0
            self.progressbar.step(progressbar_step)
            self.update()

        self.write("Tests Complete.")
        self.write("-----------------\n")
        if result:
            self.write("All tests passed!")
        else:
            self.write("One or more failures occured. See testlog.txt")

        with open('../test-dir/TestResult.dat', 'w+') as f:
            f.write(str(result))

        return 

    def write(self, string):

        print(string, end='')
        with open('../test-dir/testlog.txt', 'w+') as f:
            f.write(string)
        self.feedback.config(state=tk.NORMAL)
        self.feedback.insert(tk.END, string)
        self.feedback.see(tk.END)
        self.feedback.config(state=tk.DISABLED)
        self.update()


def getopts(argv):

    if (len(argv) == 1):
        return None
    elif (argv[1] == '-s'):
        return SHORT_TEST
    elif (argv[1] == '-l'):
        return LONG_TEST
    elif (argv[1] == '-d'):
        return DEVELOPER_TEST
    elif (argv[1] == '-startup'):
        return STARTUP_TEST

    return None


if __name__ == '__main__':

    arg = getopts(sys.argv)
    if arg is not None:
        app = TestC2E2(None, arg)
    else:
        app = TestC2E2(None)
    app.mainloop()
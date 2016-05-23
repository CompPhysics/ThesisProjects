#pythonw wlagui.py

from Tkinter import *
import os, time, Pmw, tkFileDialog, tkMessageBox, tkSimpleDialog, master
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class Gui:
    def __init__(self, parent):
        self.master = parent

        #Set title of program
        self.program_name = "Wla 0.0 "
        self.master.title(self.program_name)

        ######
        #     Need some "global" variables
        ######

        self.cwd = os.getcwd()
        
        self.signal = master.Data(0)
        self.wadata = master.Data(1)
        self.mradata = master.Data(2)
        self.ahurstdata = master.Data(3)
        self.ihurstdata = master.Data(4)
        self.genset = master.Settings(0)
        self.waset = master.Settings(1)
        self.mraset = master.Settings(2)
        self.ahurstset = master.Settings(3)
        self.ihurstset = master.Settings(4)
        self.openfile = False
        self.activeepoch = IntVar()
        #Stupid IntVar() shit!!!
        self.filetype = IntVar()
        self.channel = IntVar()
        self.waactive = IntVar()
        self.mraactive = IntVar()
        self.ahurstactive = IntVar()
        self.ihurstactive = IntVar()
        self.watranstype = IntVar()
        self.mratranstype = IntVar()
        self.ahursttranstype = IntVar()
        self.ihursttranstype = IntVar()
        self.walevel = IntVar() #Bruker j_0 og level om hverandre
        self.mralevel = IntVar()
        self.ahurstlevel = IntVar()
        self.ihurstlevel = IntVar()
        self.waboundary = IntVar()
        self.mraboundary = IntVar()
        self.ahurstboundary = IntVar()
        self.ihurstboundary = IntVar()
        self.wawavelet = IntVar()
        self.mrawavelet = IntVar()
        self.ahurstwavelet = IntVar()
        self.ihurstwavelet = IntVar()
        self.walength = IntVar()
        self.mralength = IntVar()
        self.ahurstlength = IntVar()
        self.ihurstlength = IntVar()
        self.ahurstjmin = IntVar()
        self.ihurstjmin = IntVar()
        self.ahurstjmax = IntVar()
        self.ihurstjmax = IntVar()
        self.ahurstbias = IntVar()
        self.ahurstp = IntVar()
        self.ihurstsmooth = IntVar()

        self.resolution = 1
        self.l_min = 2
        self.l_max = 30
        self.j_max = 10

        #Geometrical settings
        self.slider_length = 150
        
        #Create right hand side frame for "Compute buttons" and for displaying relevant numerics
        self.right_frame = Frame(self.master, width=300, height=640, bg="red")
        self.right_frame.pack(side='right', expand=True, fill='both')
        self.right_frame.propagate(0)
        
        self.button_frame = Frame(self.right_frame, width=300, height=30)
        self.button_frame.pack(side=BOTTOM, anchor=S, fill=BOTH)

        #Create a frame for tools like new, open, options etc 
        self.tool_frame = Frame(self.master, width = 900, height = 60, bg="black")
        self.tool_frame.pack(side='top', expand=False, fill=None)
        self.tool_frame.propagate(0)
        
        #Create main frame for displaying graphs
        self.main_frame = Frame(self.master, width=900, height=600, bg="blue")
        self.main_frame.pack(side='top', expand=True, fill='both')
        self.main_frame.propagate(0)

        #Run button
        rbutton = Button(self.button_frame, text="Run", command=self.calculate, width=16)
        rbutton.pack(side=LEFT, anchor=S)
        #rbutton.config(bg="red", bd=8)

        #Quit button
        qbutton = Button(self.button_frame, text="Quit", command=self.quit, width=15)
        qbutton.pack(side=RIGHT, anchor=S)

        #New button
        new_icon = PhotoImage(file=self.cwd+'/icons/new32.gif')
        nbutton = Button(self.tool_frame, image=new_icon, command=self.new)
        nbutton.pack(side='left', anchor='w')
        nbutton.image = new_icon       #Protect image from the garbage collector

        #Open button
        open_icon = PhotoImage(file=self.cwd+'/icons/open32.gif')
        obutton = Button(self.tool_frame, image=open_icon, command=self.open)
        obutton.pack(side='left', anchor='w')
        obutton.image = open_icon

        #Save button
        save_icon = PhotoImage(file=self.cwd+'/icons/save32.gif')
        sbutton = Button(self.tool_frame, image=save_icon, command=self.save)
        sbutton.pack(side='left', anchor='w')
        sbutton.image = save_icon

        #Zoom button
        zoom_icon = PhotoImage(file=self.cwd+'/icons/zoom32.gif')
        zbutton = Button(self.tool_frame, image=zoom_icon, command=self.zoom)
        zbutton.pack(side='left', anchor='w')
        zbutton.image = zoom_icon

        #Preferences button
        pref_icon = PhotoImage(file=self.cwd+'/icons/pref32.gif')
        pbutton = Button(self.tool_frame, image=pref_icon, command=self.preferences)
        pbutton.pack(side='left', anchor='w')
        pbutton.image = pref_icon

    def epoch_buttons(self):
        self.activeepoch.set(1)
        self.p1 = Radiobutton(self.tool_frame, text='1', variable=self.activeepoch, value=1, 
                              command=self.regraph, bg="black", fg="red")
        self.p1.pack(side=LEFT, anchor=SW, pady=7)
        self.p2 = Radiobutton(self.tool_frame, text='2', variable=self.activeepoch, value=2, 
                              command=self.regraph, bg="black", fg="red")
        self.p2.pack(side=LEFT, anchor=SW, pady=7)
        self.p3 = Radiobutton(self.tool_frame, text='3', variable=self.activeepoch, value=3, 
                              command=self.regraph, bg="black", fg="red")
        self.p3.pack(side=LEFT, anchor=SW, pady=7)
        if self.genset.filetype is 0:
            self.p4 = Radiobutton(self.tool_frame, text='4', variable=self.activeepoch, value=4, 
                                  command=self.regraph, bg="black", fg="red")
            self.p4.pack(side=LEFT, anchor=SW, pady=7)

    def new(self):
        self.signal = master.Data(0)
        self.wadata = master.Data(1)
        self.mradata = master.Data(2)
        self.ahurstdata = master.Data(3)
        self.ihurstdata = master.Data(4)
        try:
            self.waoutgeneral.destroy()
        except:
            pass
        try:
            self.mraoutgeneral.destroy()
        except:
            pass
        try:
            self.ahurstoutgeneral.destroy()
            self.showahurstout1.destroy()
            self.showahurstout2.destroy()
            self.showahurstout3.destroy()
        except:
            pass
        try:
            self.showahurstout4.destroy()
        except:
            pass
        try:
            self.ihurstoutgeneral.destroy()
            self.showihurstout1.destroy()
            self.showihurstout2.destroy()
            self.showihurstout3.destroy()
        except:
            pass
        try:
            self.showihurstout4.destroy()
        except:
            pass
        try:
            self.graphs.destroy()
        except:
            pass
        try:
            self.p1.destroy()
            self.p2.destroy()
            self.p3.destroy()
        except:
            pass
        try:
            self.p1.destroy()
            self.p2.destroy()
            self.p3.destroy()
            self.p4.destroy()
        except:
            pass
        self.openfile = False
        self.master.title(self.program_name)
    
    def open2(self):
        self.ifilename = '/Users/vegard/Master/Hurst05.dat'
        for value in open(self.ifilename):
            self.signal.signal.append(float(value.split()[1]))
            #self.signal.signal.append(float(value))
        self.signal.e1 = self.signal.signal
        self.genset.f = 1
        self.signal.set_time(self.genset.f)
        self.master.title(self.program_name + ':' + self.ifilename)
        self.graphs = Pmw.NoteBook(self.main_frame) 
        self.graphs.pack(expand=True, fill=BOTH)
        self.signalplot = self.graphs.add("Signal")
        self.show_graph(self.signalplot, self.signal.time, self.signal.signal)
        self.openfile = True
        self.activeepoch.set(1)

    def open(self): #Make sure that if a file is opened, it is impossible to open another file
        def validtime(input):
            try:
                time.strptime(input, "%H:%M:%S")
                return True
            except:
                tkMessageBox.showerror("IO Error", "Time format not valid")
                return False

        if not self.openfile:
            self.ifilename = tkFileDialog.askopenfilename(initialdir=self.cwd+'/data')
            if(self.genset.filetype == 0):
                try:
                    self.signal.signal = master.read_oldfile(self.ifilename, self.genset.channel)
                    self.genset.f = 1
                    valid = False
                    while not valid:
                        ep1 = tkSimpleDialog.askstring("Epoch 1", "Enter start time epoch 1 (hh::mm::ss):")
                        valid = validtime(ep1)
                    valid = False
                    while not valid:
                        ep2 = tkSimpleDialog.askstring("Epoch 2", "Enter start time epoch 2 (hh::mm::ss):")
                        valid = validtime(ep2)
                    valid = False
                    while not valid:
                        ep3 = tkSimpleDialog.askstring("Epoch 3", "Enter start time epoch 3 (hh::mm::ss):")
                        valid = validtime(ep3)
                    valid = False
                    while not valid:
                        ep4 = tkSimpleDialog.askstring("Epoch 4", "Enter start time epoch 4 (hh::mm::ss):")
                        valid = validtime(ep4)
                    master.get_epochtimes(self.ifilename, self.signal, ['1','2','3','4'], [ep1,ep2,ep3,ep4])
                    master.get_epochs(self.genset, self.signal)
                    self.signal.set_time(self.genset.f)
                    self.master.title(self.program_name + ':' + self.ifilename)
                    self.graphs = Pmw.NoteBook(self.main_frame) #This fuck up the blue color
                    self.graphs.pack(expand=True, fill=BOTH)
                    self.signalplot = self.graphs.add("Signal")
                    self.show_graph(self.signalplot, self.signal.time, self.signal.signal)
                    self.epochplot = self.graphs.add("Epoch 1")
                    self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
                    self.epoch_buttons()
                    self.openfile = True
                except:
                    tkMessageBox.showerror("IO Error", "Error in file")
            elif(self.genset.filetype == 1):
                try:
                    self.signal.signal, self.signal.epoch_times = master.read_newfile(self.ifilename)
                    master.get_epochs(self.genset, self.signal)
                    self.signal.set_time(self.genset.f)
                    self.master.title(self.program_name + ':' + self.ifilename)
                    self.graphs = Pmw.NoteBook(self.main_frame) #This fuck up the blue color
                    self.graphs.pack(expand=True, fill=BOTH)
                    self.signalplot = self.graphs.add("Signal")
                    self.show_signal(self.signalplot, self.signal.time, self.signal.signal, self.signal.epoch_times)
                    self.epochplot = self.graphs.add("Epoch 1")
                    self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
                    self.epoch_buttons()
                    self.openfile = True
                except:
                    tkMessageBox.showerror("IO Error", "Error in file")
        else:
            tkMessageBox.showerror("IO Error", "A file is already open")

    def save(self):
        ofilename = tkFileDialog.asksaveasfilename(initialdir=self.cwd+"/saves")
        try:
            master.save(ofilename, self.program_name+self.ifilename, self.genset, 
                        self.ahurstset, self.ihurstset, self.ahurstdata, self.ihurstdata)
        except:
            tkMessageBox.showerror("IO Error", "Nothing to save...")

    def zoom(self): 
        active = self.graphs.getcurselection()
        f = mpl.figure(figsize=(10,5),dpi=100)
        f.canvas.set_window_title(active)
        if(active == "Signal"):
            self.zoom_graph(f, self.signal.time, self.signal.signal)
        elif(active == "Epoch 1"):
            self.zoom_graph(f, self.signal.etime, self.signal.e1)
        elif(active == "Epoch 2"):
            self.zoom_graph(f, self.signal.etime, self.signal.e2)
        elif(active == "Epoch 3"):
            self.zoom_graph(f, self.signal.etime, self.signal.e3)
        elif(active == "Epoch 4"):
            self.zoom_graph(f, self.signal.etime, self.signal.e4)
        elif active is "WA": 
            if self.activeepoch.get() is 1:
                self.zoom_multiple(f, self.signal.etime, self.wadata.wv1, self.wadata.boundaries)
            elif self.activeepoch.get() is 2:
                self.zoom_multiple(f, self.signal.etime, self.wadata.wv2, self.wadata.boundaries)
            elif self.activeepoch.get() is 3:
               self.zoom_multiple(f, self.signal.etime, self.wadata.wv3, self.wadata.boundaries)
            elif self.activeepoch.get() is 4:
               self.zoom_multiple(f, self.signal.etime, self.wadata.wv4, self.wadata.boundaries)
        elif active is "MRA":
            if self.activeepoch.get() is 1:
               self.zoom_multiple(f, self.signal.etime, self.mradata.ds1, self.mradata.boundaries)
            elif self.activeepoch.get() is 2:
               self.zoom_multiple(f, self.signal.etime, self.mradata.ds2, self.mradata.boundaries)
            elif self.activeepoch.get() is 3:
              self.zoom_multiple(f, self.signal.etime, self.mradata.ds3, self.mradata.boundaries)
            elif self.activeepoch.get() is 4:
              self.zoom_multiple(f, self.signal.etime, self.mradata.ds4, self.wadata.boundaries)
        elif active is "aHurst":
            if self.activeepoch.get() is 1:
                self.zoom_logplot(f, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max, 
                                  self.ahurstdata.y_hat1, self.ahurstdata.y1, self.ahurstdata.ci1)
            elif self.activeepoch.get() is 2:
                self.zoom_logplot(f, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max, 
                                  self.ahurstdata.y_hat2, self.ahurstdata.y2, self.ahurstdata.ci2)
            elif self.activeepoch.get() is 3:
                self.zoom_logplot(f, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max, 
                                  self.ahurstdata.y_hat3, self.ahurstdata.y3, self.ahurstdata.ci3)
            elif self.activeepoch.get() is 4:
                self.zoom_logplot(f, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max, 
                                  self.ahurstdata.y_hat4, self.ahurstdata.y4, self.ahurstdata.ci4)
        elif active is "iHurst":
            if self.activeepoch.get() is 1:
                self.zoom_ihurst(f, self.signal.etime, self.ihurstdata.hurst1, self.ihurstdata.out1)
            elif self.activeepoch.get() is 2:
                self.zoom_ihurst(f, self.signal.etime, self.ihurstdata.hurst2, self.ihurstdata.out2)
            elif self.activeepoch.get() is 3:
                self.zoom_ihurst(f, self.signal.etime, self.ihurstdata.hurst3, self.ihurstdata.out3)
            elif self.activeepoch.get() is 4:
                self.zoom_ihurst(f, self.signal.etime, self.ihurstdata.hurst4, self.ihurstdata.out4)

    def preferences(self): #Consider putting this in a class on its own
        def build_generalpage(rebuild):
            self.generalfgroup = Pmw.Group(self.generalpage, tag_text="Filetype")
            self.generalfgroup.pack(fill=BOTH, expand = 1)
            if not rebuild:
                self.filetype.set(self.genset.filetype)
            rb1 = Radiobutton(self.generalfgroup.interior(), text="Old file type", variable=self.filetype, 
                              value=0, command=rebuild_generalpage)
            rb1.pack(side=TOP, anchor=NW)
            rb2 = Radiobutton(self.generalfgroup.interior(), text="New file type", variable=self.filetype, 
                              value=1, command=rebuild_generalpage)
            rb2.pack(side=TOP, anchor=NW)
            
            self.state = NORMAL
            if(self.filetype.get() == 1):
                self.state= DISABLED
            self.channel.set(self.genset.channel)
            self.generalcgroup = Pmw.Group(self.generalpage, tag_text="Channel")
            self.generalcgroup.pack(fill=BOTH, expand=1)
            rb3 = Radiobutton(self.generalcgroup.interior(), text="Channel 1", variable=self.channel, value=1, 
                              state=self.state)
            rb3.pack(side=TOP, anchor=NW)
            rb4 = Radiobutton(self.generalcgroup.interior(), text="Channel 2", variable=self.channel, value=2, 
                              state=self.state)
            rb4.pack(side=TOP, anchor=NW)
            rb5 = Radiobutton(self.generalcgroup.interior(), text="Channel 3", variable=self.channel, value=3, 
                              state=self.state)
            rb5.pack(side=TOP, anchor=NW)
            rb6 = Radiobutton(self.generalcgroup.interior(), text="Channel 4", variable=self.channel, value=4, 
                              state=self.state)
            rb6.pack(side=TOP, anchor=NW)

            self.generalagroup = Pmw.Group(self.generalpage, tag_text="Analysis")
            self.generalagroup.pack(fill=BOTH, expand = 1)
            self.waactive.set(self.waset.active)
            cb1 = Checkbutton(self.generalagroup.interior(), text="Wavelet analysis", 
                              variable=self.waactive, command=rebuild_wapage)
            cb1.pack(anchor=NW)
            self.mraactive.set(self.mraset.active)
            cb2 = Checkbutton(self.generalagroup.interior(), text="Multiresolution analysis", 
                              variable = self.mraactive, command=rebuild_mrapage)
            cb2.pack(anchor=NW)
            self.ahurstactive.set(self.ahurstset.active)
            cb3 = Checkbutton(self.generalagroup.interior(), text="Averaged Hurst estimation", 
                              variable = self.ahurstactive, command=rebuild_ahurstpage)
            cb3.pack(anchor=NW)
            self.ihurstactive.set(self.ihurstset.active)
            cb4 = Checkbutton(self.generalagroup.interior(), text="Instant Hurst estimation", 
                              variable = self.ihurstactive, command=rebuild_ihurstpage)
            cb4.pack(anchor=NW)
        
        def build_wapage():
            self.state = NORMAL
            if(self.waactive.get() == 0):
                self.state = DISABLED

            self.watgroup = Pmw.Group(self.wapage, tag_text="Transform")
            self.watgroup.pack(fill=BOTH, expand=1)
            self.watranstype.set(self.waset.transtype)
            rb1 = Radiobutton(self.watgroup.interior(), text="DWT", variable=self.watranstype, value=0, 
                              state=self.state, command=update_waslider1b)
            rb1.pack(side =TOP, anchor=NW)
            rb2 = Radiobutton(self.watgroup.interior(), text="MODWT", variable=self.watranstype, value=1, 
                              state=self.state, command=update_waslider1b)
            rb2.pack(side =TOP, anchor=NW)
            self.walevel.set(self.waset.j_0)
            self.resolution = 1
            if not self.openfile:
                self.state = DISABLED
            self.waslider1 = Scale(self.watgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Transformation level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.walevel, state=self.state)
            if(self.waactive.get() == 1):
                self.state = NORMAL
            self.waslider1.pack(anchor=W, pady=20)
            text = Label(self.watgroup.interior(), text="Boundary conditions:")
            text.pack(side=TOP, anchor=NW)
            self.waboundary.set(self.waset.refl)
            rb3 = Radiobutton(self.watgroup.interior(), text="Circular", variable=self.waboundary, 
                              value=0, state=self.state)
            rb3.pack(side=LEFT, anchor=NW)
            rb4 = Radiobutton(self.watgroup.interior(), text="Reflection", variable=self.waboundary, 
                              value=1, state=self.state)
            rb4.pack(side=LEFT, anchor=NW)
        
            self.wawgroup = Pmw.Group(self.wapage, tag_text="Wavelet")
            self.wawgroup.pack(fill=BOTH, expand=1)
            self.wawavelet.set(self.waset.type)
            rb1 = Radiobutton(self.wawgroup.interior(), text="Daubechies", variable=self.wawavelet, value=1, 
                              command=set_walength, state=self.state)
            rb1.pack(anchor=W)
            rb2 = Radiobutton(self.wawgroup.interior(), text="Least Asymmetrical", variable=self.wawavelet, value=2, 
                              command=set_walength, state=self.state)
            rb2.pack(anchor=W)
            rb3 = Radiobutton(self.wawgroup.interior(), text="Best Localized", variable=self.wawavelet, value=3, 
                              command=set_walength, state=self.state)
            rb3.pack(anchor=W)
            rb4 = Radiobutton(self.wawgroup.interior(), text="Coiflet", variable=self.wawavelet, value=4, 
                              command=set_walength, state=self.state)
            rb4.pack(anchor=W)
            self.walength.set(self.waset.l)
            if(self.wawavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.wawavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.wawavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.wawavelet.get() == 4):
                self.resolution = 6
                self.l_min = 2
                self.l_max = 30
            self.waslider2 = Scale(self.wawgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                                   label="Wavelet filter length:", from_=self.l_min,to=self.l_max, 
                                   resolution=self.resolution, variable=self.walength, state=self.state, 
                                   command=update_waslider1)
            self.waslider2.pack(anchor=W, pady=20)

        def build_mrapage():
            self.state = NORMAL
            if(self.mraactive.get() == 0):
                self.state = DISABLED

            self.mratgroup = Pmw.Group(self.mrapage, tag_text="Transform")
            self.mratgroup.pack(fill=BOTH, expand=1)
            self.mratranstype.set(self.mraset.transtype)
            rb1 = Radiobutton(self.mratgroup.interior(), text="DWT", variable=self.mratranstype, value=0, 
                              state=self.state, command=update_mraslider1b)
            rb1.pack(side =TOP, anchor=NW)
            rb2 = Radiobutton(self.mratgroup.interior(), text="MODWT", variable=self.mratranstype, value=1, 
                              state=self.state, command=update_mraslider1b)
            rb2.pack(side =TOP, anchor=NW)
            self.mralevel.set(self.mraset.j_0)
            self.resolution = 1
            if not self.openfile:
                self.state = DISABLED
            self.mraslider1 = Scale(self.mratgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Transformation level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.mralevel, state=self.state)
            if(self.mraactive.get() == 1):
                self.state = NORMAL 
            self.mraslider1.pack(anchor=W, pady=20)
            text = Label(self.mratgroup.interior(), text="Boundary conditions:")
            text.pack(side=TOP, anchor=NW)
            self.mraboundary.set(self.mraset.refl)
            rb3 = Radiobutton(self.mratgroup.interior(), text="Circular", variable=self.mraboundary, value=0, 
                              state=self.state)
            rb3.pack(side=LEFT, anchor=NW)
            rb4 = Radiobutton(self.mratgroup.interior(), text="Reflection", variable=self.mraboundary, value=1, 
                              state=self.state)
            rb4.pack(side=LEFT, anchor=NW)
        
            self.mrawgroup = Pmw.Group(self.mrapage, tag_text="Wavelet")
            self.mrawgroup.pack(fill=BOTH, expand=1)
            self.mrawavelet.set(self.mraset.type)
            rb1 = Radiobutton(self.mrawgroup.interior(), text="Daubechies", variable=self.mrawavelet, value=1, 
                              command=set_mralength, state=self.state)
            rb1.pack(anchor=W)
            rb2 = Radiobutton(self.mrawgroup.interior(), text="Least Asymmetrical", variable=self.mrawavelet, value=2, 
                              command=set_mralength, state=self.state)
            rb2.pack(anchor=W)
            rb3 = Radiobutton(self.mrawgroup.interior(), text="Best Localized", variable=self.mrawavelet, value=3, 
                              command=set_mralength, state=self.state)
            rb3.pack(anchor=W)
            rb4 = Radiobutton(self.mrawgroup.interior(), text="Coiflet", variable=self.mrawavelet, value=4, 
                              command=set_mralength, state=self.state)
            rb4.pack(anchor=W)
            self.mralength.set(self.mraset.l)
            if(self.mrawavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.mrawavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.mrawavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.mrawavelet.get() == 4):
                self.resolution = 6
                self.l_min = 2
                self.l_max = 30
            self.mraslider2 = Scale(self.mrawgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                                    label="Wavelet filter length:", from_=self.l_min,to=self.l_max, 
                                    resolution=self.resolution, variable=self.mralength, state=self.state,
                                    command=update_mraslider1)
            self.mraslider2.pack(anchor=W, pady=20)

        def build_ahurstpage():
            self.state = NORMAL
            if(self.ahurstactive.get() == 0):
                self.state = DISABLED
                
            self.ahursttgroup = Pmw.Group(self.ahurstpage, tag_text="Transform")
            self.ahursttgroup.pack(fill=BOTH, expand=1)
            self.ahursttranstype.set(self.ahurstset.transtype)
            rb1 = Radiobutton(self.ahursttgroup.interior(), text="DWT", variable=self.ahursttranstype, value=0, 
                              state=self.state, command=update_ahurstslider1b)
            rb1.pack(side =TOP, anchor=NW)
            rb2 = Radiobutton(self.ahursttgroup.interior(), text="MODWT", variable=self.ahursttranstype, value=1, 
                              state=self.state, command=update_ahurstslider1b)
            rb2.pack(side =TOP, anchor=NW)
            self.ahurstlevel.set(self.ahurstset.j_0)
            #Set ahurstlevel if file is opened
            self.resolution = 1
            if not self.openfile:
                self.state = DISABLED
            self.ahurstslider1 = Scale(self.ahursttgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Transformation level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ahurstlevel, state=self.state, command=update_ahurstslider4)
            if(self.ahurstactive.get() == 1):
                self.state = NORMAL
            self.ahurstslider1.pack(anchor=W, pady=20)
            text = Label(self.ahursttgroup.interior(), text="Boundary conditions:")
            text.pack(side=TOP, anchor=NW)
            self.ahurstboundary.set(self.ahurstset.refl)
            rb3 = Radiobutton(self.ahursttgroup.interior(), text="Circular", variable=self.ahurstboundary, 
                              value=0, state=self.state, command=update_ahurstslider1b)
            rb3.pack(side=LEFT, anchor=NW)
            rb4 = Radiobutton(self.ahursttgroup.interior(), text="Reflection", variable=self.ahurstboundary, 
                              value=1, state=self.state, command=update_ahurstslider1b)
            rb4.pack(side=LEFT, anchor=NW)
            
            self.ahurstwgroup = Pmw.Group(self.ahurstpage, tag_text="Wavelet")
            self.ahurstwgroup.pack(fill=BOTH, expand=1)
            self.ahurstwavelet.set(self.ahurstset.type)
            rb1 = Radiobutton(self.ahurstwgroup.interior(), text="Daubechies", variable=self.ahurstwavelet, 
                              value=1, command=set_ahurstlength, state=self.state)
            rb1.pack(anchor=W)
            rb2 = Radiobutton(self.ahurstwgroup.interior(), text="Least Asymmetrical", variable=self.ahurstwavelet, 
                              value=2, command=set_ahurstlength, state=self.state)
            rb2.pack(anchor=W)
            rb3 = Radiobutton(self.ahurstwgroup.interior(), text="Best Localized", variable=self.ahurstwavelet, 
                              value=3, command=set_ahurstlength, state=self.state)
            rb3.pack(anchor=W)
            rb4 = Radiobutton(self.ahurstwgroup.interior(), text="Coiflet", variable=self.ahurstwavelet, value=4, 
                              command=set_ahurstlength, state=self.state)
            rb4.pack(anchor=W)
            self.ahurstlength.set(self.ahurstset.l)
            if(self.ahurstwavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 4):
                self.resolution = 6
                self.l_min = 2
                self.l_max = 30
            self.ahurstslider2 = Scale(self.ahurstwgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                                       label="Wavelet filter length:", from_=self.l_min,to=self.l_max, 
                                       resolution=self.resolution, variable=self.ahurstlength, state=self.state,
                                       command=update_ahurstslider1)
            self.ahurstslider2.pack(anchor=W, pady=20)

            self.ahurstwavgroup = Pmw.Group(self.ahurstpage, tag_text="Wavelet Variance")
            self.ahurstwavgroup.pack(fill=BOTH, expand=1)
            self.ahurstbias.set(self.ahurstset.bias)
            rb1_frame = Frame(self.ahurstwavgroup.interior())
            rb1_frame.pack(anchor=NW)
            rb1 = Radiobutton(rb1_frame, text="Unbiased", variable=self.ahurstbias, value=0, state=self.state, 
                              command=update_ahurstslider1b)
            rb1.pack(side=LEFT, anchor=NW)
            rb2 = Radiobutton(rb1_frame, text="Biased", variable=self.ahurstbias, value=1, state=self.state, 
                              command=update_ahurstslider1b)
            rb2.pack(side=LEFT, anchor=NW)
            empty_frame = Frame(self.ahurstwavgroup.interior())
            empty_frame.pack(anchor=NW, pady=10)
            text_frame = Frame(self.ahurstwavgroup.interior())
            text_frame.pack(anchor=NW)
            rb2_frame = Frame(self.ahurstwavgroup.interior())
            rb2_frame.pack(anchor=NW)
            text=Label(text_frame, text="Confidence interval:")
            text.pack(side=LEFT, anchor=NW)
            self.ahurstp.set(self.ahurstset.p)
            rb3 = Radiobutton(rb2_frame, text="90", variable=self.ahurstp, value=90, state=self.state)
            rb3.pack(side=LEFT, anchor=NW)
            rb4 = Radiobutton(rb2_frame, text="95", variable=self.ahurstp, value=95, state=self.state)
            rb4.pack(side=LEFT, anchor=NW)
            rb5 = Radiobutton(rb2_frame, text="99", variable=self.ahurstp, value=99, state=self.state)
            rb5.pack(side=LEFT, anchor=NW)

            self.ahurstwlsegroup = Pmw.Group(self.ahurstpage, tag_text="WLSE")
            self.ahurstwlsegroup.pack(fill=BOTH, expand=1)
            self.ahurstjmin.set(self.ahurstset.j_min)
            self.resolution = 1
            self.slider_length = 160
            if not self.openfile:
                self.state = DISABLED
            ahurstslider3 = Scale(self.ahurstwlsegroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Lowest regression level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ahurstjmin, state=self.state)
            ahurstslider3.pack(side=LEFT, anchor=NW)
            self.ahurstjmax.set(self.ahurstset.j_max)
            self.ahurstslider4 = Scale(self.ahurstwlsegroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Highest regression level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ahurstjmax, state=self.state)
            self.ahurstslider4.pack(side=LEFT, anchor=NW, padx=20)
            self.state = NORMAL
            self.slider_length = 150

        def build_ihurstpage():
            self.state = NORMAL
            if(self.ihurstactive.get() == 0):
                self.state = DISABLED

            self.ihursttgroup = Pmw.Group(self.ihurstpage, tag_text="Transform")
            self.ihursttgroup.pack(fill=BOTH, expand=1)
            self.ihursttranstype.set(self.ihurstset.transtype)
            rb1 = Radiobutton(self.ihursttgroup.interior(), text="DWT", variable=self.ihursttranstype, 
                              value=0, state=DISABLED)
            rb1.pack(side =TOP, anchor=NW)
            rb2 = Radiobutton(self.ihursttgroup.interior(), text="MODWT", variable=self.ihursttranstype, 
                              value=1, state=DISABLED)
            rb2.pack(side =TOP, anchor=NW)
            self.ihurstlevel.set(self.ihurstset.j_0)
            self.resolution = 1
            if not self.openfile:
                self.state = DISABLED
            self.ihurstslider1 = Scale(self.ihursttgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Transformation level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ihurstlevel, state=self.state, command=update_ihurstslider45)
            if(self.ihurstactive.get() == 1):
                self.state = NORMAL
            self.ihurstslider1.pack(anchor=W, pady=20)
            text = Label(self.ihursttgroup.interior(), text="Boundary conditions:")
            text.pack(side=TOP, anchor=NW)
            self.ihurstboundary.set(self.ihurstset.refl)
            rb3 = Radiobutton(self.ihursttgroup.interior(), text="Circular", variable=self.ihurstboundary, 
                              value=0, state=self.state)
            rb3.pack(side=LEFT, anchor=NW)
            rb4 = Radiobutton(self.ihursttgroup.interior(), text="Reflection", variable=self.ihurstboundary, 
                              value=1, state=self.state)
            rb4.pack(side=LEFT, anchor=NW)
            
            self.ihurstwgroup = Pmw.Group(self.ihurstpage, tag_text="Wavelet")
            self.ihurstwgroup.pack(fill=BOTH, expand=1)
            self.ihurstwavelet.set(self.ihurstset.type)
            rb1 = Radiobutton(self.ihurstwgroup.interior(), text="Daubechies", variable=self.ihurstwavelet, 
                              value=1, command=set_ihurstlength, state=DISABLED)
            rb1.pack(anchor=W)
            rb2 = Radiobutton(self.ihurstwgroup.interior(), text="Least Asymmetrical", variable=self.ihurstwavelet, 
                              value=2, command=set_ihurstlength, state=self.state)
            rb2.pack(anchor=W)
            rb3 = Radiobutton(self.ihurstwgroup.interior(), text="Best Localized", variable=self.ihurstwavelet, 
                              value=3, command=set_ihurstlength, state=self.state)
            rb3.pack(anchor=W)
            rb4 = Radiobutton(self.ihurstwgroup.interior(), text="Coiflet", variable=self.ihurstwavelet, value=4, 
                              command=set_ihurstlength, state=self.state)
            rb4.pack(anchor=W)
            self.ihurstlength.set(self.ihurstset.l)
            if(self.ihurstwavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 4):
                self.resolution = 6
                self.l_min = 2
                self.l_max = 30
            self.ihurstslider2 = Scale(self.ihurstwgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                                       label="Wavelet filter length:", from_=self.l_min,to=self.l_max, 
                                       resolution=self.resolution, variable=self.ihurstlength, state=self.state,
                                       command=update_ihurstslider1)
            self.ihurstslider2.pack(anchor=W, pady=20)

            self.ihurstiwlsegroup = Pmw.Group(self.ihurstpage, tag_text="iWLSE")
            self.ihurstiwlsegroup.pack(fill=BOTH, expand=1)
            self.ihurstjmin.set(self.ihurstset.j_min)
            self.resolution = 1
            if not self.openfile:
                self.state = DISABLED
            ihurstslider3 = Scale(self.ihurstiwlsegroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Lowest regression level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ihurstjmin, state=self.state)
            ihurstslider3.pack(side=LEFT, anchor=NW)
            self.ihurstjmax.set(self.ihurstset.j_max)
            self.slider_length = 160
            self.ihurstslider4 = Scale(self.ihurstiwlsegroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                            label="Highest regression level", from_=1,to=self.j_max, resolution=self.resolution, 
                            variable=self.ihurstjmax, state=self.state)
            self.ihurstslider4.pack(side=LEFT, anchor=NW, padx=20)
            self.state = NORMAL
            self.slider_length = 150
            self.ihurstsmoothgroup = Pmw.Group(self.ihurstpage, tag_text="MRA smooth")
            self.ihurstsmoothgroup.pack(fill=BOTH, expand=1)
            self.resolution = 1
            self.ihurstsmooth.set(self.ihurstset.smoothlevel)
            self.ihurstslider5 = Scale(self.ihurstsmoothgroup.interior(), orient=HORIZONTAL, length=self.slider_length, 
                                  label="Smooth level:", from_=0,to=self.j_max, resolution=self.resolution, 
                                  variable=self.ihurstsmooth, state=self.state)
            self.ihurstslider5.pack(side=LEFT, anchor=NW)

        def rebuild_generalpage():
            self.generalfgroup.destroy()
            self.generalcgroup.destroy()
            self.generalagroup.destroy()
            build_generalpage(True)

        def rebuild_wapage():
            self.watgroup.destroy()
            self.wawgroup.destroy()
            build_wapage()
        
        def rebuild_mrapage():
            self.mratgroup.destroy()
            self.mrawgroup.destroy()
            build_mrapage()
        
        def rebuild_ahurstpage():
            self.ahursttgroup.destroy()
            self.ahurstwgroup.destroy()
            self.ahurstwavgroup.destroy()
            self.ahurstwlsegroup.destroy()
            build_ahurstpage()
        
        def rebuild_ihurstpage():
            self.ihursttgroup.destroy()
            self.ihurstwgroup.destroy()
            self.ihurstiwlsegroup.destroy()
            build_ihurstpage()

        def update_waslider1(length):
            self.waslider1.configure(to=master.get_maxlevel(self.watranstype.get(), 
                                                            len(self.signal.etime), int(length)))
        
        def update_waslider1b():
            self.waslider1.configure(to=master.get_maxlevel(self.watranstype.get(), 
                                                            len(self.signal.etime), self.walength.get()))

        def update_mraslider1(length):
            self.mraslider1.configure(to=master.get_maxlevel(self.mratranstype.get(), 
                                                            len(self.signal.etime), int(length)))
        
        def update_mraslider1b():
            self.mraslider1.configure(to=master.get_maxlevel(self.watranstype.get(), 
                                                            len(self.signal.etime), self.mralength.get()))

        def update_ahurstslider1(length):
            if self.ahurstbias.get() is 0:
                if self.ahurstboundary.get() is 0:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              len(self.signal.etime), 
                                                                              int(length),
                                                                              self.ahurstbias.get()))
                elif self.ahurstboundary.get() is 1:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              2*len(self.signal.etime), 
                                                                              int(length),
                                                                              self.ahurstbias.get()))
            elif self.ahurstbias.get() is 1:
                if self.ahurstboundary.get() is 0:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              len(self.signal.etime), 
                                                                              int(length),
                                                                              self.ahurstbias.get()))
                elif self.ahurstboundary.get() is 1:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              2*len(self.signal.etime), 
                                                                              int(length),
                                                                              self.ahurstbias.get()))
                                             
        def update_ahurstslider1b():
            if self.ahurstbias.get() is 0:
                if self.ahurstboundary.get() is 0:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              len(self.signal.etime), 
                                                                              self.ahurstlength.get(),
                                                                              self.ahurstbias.get()))
                elif self.ahurstboundary.get() is 1:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              2*len(self.signal.etime), 
                                                                              self.ahurstlength.get(),
                                                                              self.ahurstbias.get()))
            elif self.ahurstbias.get() is 1:
                if self.ahurstboundary.get() is 0:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              len(self.signal.etime), 
                                                                              self.ahurstlength.get(),
                                                                              self.ahurstbias.get()))
                elif self.ahurstboundary.get() is 1:
                    self.ahurstslider1.configure(to=master.get_ahurstmaxlevel(self.ahursttranstype.get(), 
                                                                              2*len(self.signal.etime), 
                                                                              self.ahurstlength.get(),
                                                                              self.ahurstbias.get()))

        def update_ihurstslider1(length):
            if self.ihurstboundary.get() is 0:
                print "Circular"
                self.ihurstslider1.configure(to=master.get_ahurstmaxlevel(self.ihursttranstype.get(), 
                                                                          len(self.signal.etime), 
                                                                          self.ihurstlength.get(),
                                                                          1))
   #10 is max... (Need to be able to shift...)
           #     self.ihurstslider1.configure(to=master.get_maxlevel(self.ihursttranstype.get(), 
            #                                                        len(self.signal.etime), int(length)))
           
        def update_ahurstslider4(j_0):
            if self.ahurstbias.get() is 0:
                self.ahurstslider4.configure(to=j_0)
            elif self.ahurstbias.get() is 1:
                print "Biased" #N-L_j >= 0 Check which transformation levels that are allowed now
                self.ahurstslider4.configure(to=j_0)

        def update_ihurstslider45(j_0):
            self.ihurstslider4.configure(to=j_0)
            self.ihurstslider5.configure(to=j_0)

        def set_walength():
            self.waslider2.destroy()
            if(self.wawavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.wawavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.wawavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.wawavelet.get() == 4):
                self.resolution = 6
                self.l_min = 6
                self.l_max = 30
            self.waslider2 = Scale(self.wawgroup.interior(), orient=HORIZONTAL, length=self.slider_length ,
                           label="Wavelet filter length:", from_=self.l_min, to=self.l_max, resolution=self.resolution, 
                           variable=self.walength, command=update_waslider1)
            self.waslider2.pack(anchor=W, pady=20)

        def set_mralength():
            self.mraslider2.destroy()
            if(self.mrawavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.mrawavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.mrawavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.mrawavelet.get() == 4):
                self.resolution = 6
                self.l_min = 6
                self.l_max = 30
            self.mraslider2 = Scale(self.mrawgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                           label="Wavelet filter length:", from_=self.l_min, to=self.l_max, resolution=self.resolution, 
                           variable=self.mralength, command=update_mraslider1)
            self.mraslider2.pack(anchor=W, pady=20)

        def set_ahurstlength():
            self.ahurstslider2.destroy()
            if(self.ahurstwavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ahurstwavelet.get() == 4):
                self.resolution = 6
                self.l_min = 6
                self.l_max = 30
            self.ahurstslider2 = Scale(self.ahurstwgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                           label="Wavelet filter length:", from_=self.l_min, to=self.l_max, resolution=self.resolution, 
                           variable=self.ahurstlength, command=update_ahurstslider1)
            self.ahurstslider2.pack(anchor=W, pady=20)

        def set_ihurstlength():
            self.ihurstslider2.destroy()
            if(self.ihurstwavelet.get() == 1):
                self.resolution = 2
                self.l_min = 2
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 2):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 3):
                self.resolution = 2
                self.l_min = 8
                self.l_max = 20
            elif(self.ihurstwavelet.get() == 4):
                self.resolution = 6
                self.l_min = 6
                self.l_max = 30
            self.ihurstslider2 = Scale(self.ihurstwgroup.interior(), orient=HORIZONTAL, length=self.slider_length , 
                           label="Wavelet filter length:", from_=self.l_min, to=self.l_max, resolution=self.resolution, 
                           variable=self.ihurstlength, command=update_ihurstslider1)
            self.ihurstslider2.pack(anchor=W, pady=20)

        def ok():
            self.genset.filetype = self.filetype.get()
            self.genset.channel = self.channel.get()
            self.waset.active = self.waactive.get()
            self.mraset.active = self.mraactive.get()
            self.ahurstset.active = self.ahurstactive.get()
            self.ihurstset.active = self.ihurstactive.get()
            if(self.waactive.get() == 1):
                self.waset.transtype = self.watranstype.get()
                self.waset.j_0 = self.walevel.get()
                self.waset.refl = self.waboundary.get()
                self.waset.type = self.wawavelet.get()
                self.waset.l = self.walength.get()
            if(self.mraactive.get() == 1):
                self.mraset.transtype = self.mratranstype.get()
                self.mraset.j_0 = self.mralevel.get()
                self.mraset.refl = self.mraboundary.get()
                self.mraset.type = self.mrawavelet.get()
                self.mraset.l = self.mralength.get()
            if(self.ahurstactive.get() == 1):
                self.ahurstset.transtype = self.ahursttranstype.get()
                self.ahurstset.j_0 = self.ahurstlevel.get()
                self.ahurstset.refl = self.ahurstboundary.get()
                self.ahurstset.type = self.ahurstwavelet.get()
                self.ahurstset.l = self.ahurstlength.get()
                self.ahurstset.bias = self.ahurstbias.get()
                self.ahurstset.p = self.ahurstp.get()
                self.ahurstset.j_min = self.ahurstjmin.get()
                self.ahurstset.j_max = self.ahurstjmax.get()
            if(self.ihurstactive.get() == 1):
                print self.ihurstlevel.get()
                print self.ihurstjmax.get()
                self.ihurstset.j_0 = self.ihurstlevel.get()
                self.ihurstset.refl = self.ihurstboundary.get()
                self.ihurstset.type = self.ihurstwavelet.get()
                self.ihurstset.l = self.ihurstlength.get()
                self.ihurstset.j_min = self.ihurstjmin.get()
                self.ihurstset.j_max = self.ihurstjmax.get()
                self.ihurstset.smoothlevel = self.ihurstsmooth.get()
            self.preferences.destroy()

        def cancel():
            self.preferences.destroy()
        
        def set_title(name):
            self.preferences.title("Preferences: " + name)
        
        self.preferences = Toplevel()
        self.preferences.title("Preferences: General")
        self.preferences.geometry("440x712")
        self.nb = Pmw.NoteBook(self.preferences, raisecommand=set_title)
        self.nb.pack(expand=True, fill=BOTH)
        self.generalpage = self.nb.add("General")
        build_generalpage(False)
        self.wapage = self.nb.add("WA")
        build_wapage()
        self.mrapage = self.nb.add("MRA")
        build_mrapage()
        self.ahurstpage = self.nb.add("aHurst")
        build_ahurstpage()
        self.ihurstpage = self.nb.add("iHurst")
        build_ihurstpage()
        cancelbutton = Button(self.preferences, text="Cancel", command=cancel)
        cancelbutton.pack(side=RIGHT)
        okbutton = Button(self.preferences, text="Ok", command=ok)
        okbutton.pack(side=RIGHT)
  
    def zoom_graph(self, f, x, y):
        a = f.add_subplot(111)
        a.plot(x, y)
        a.set_xlabel("Time (min)")
        a.set_ylabel("Conductance (G)")
        f.show()
        
    def zoom_multiple(self, f, x, y, boundaries):
        i = len(y)
        k = 1
        low = 0
        up = len(boundaries)/2
        for list in y:
            a = f.add_subplot(i,1,k)
            a.set_axis_off()
            a.plot(x, list)
            lowvalue = boundaries[low]/float(self.genset.f*60)
            upvalue = boundaries[up]/float(self.genset.f*60)
            a.axvline(x=lowvalue,ymin=0, ymax=2, color='r')
            a.axvline(x=upvalue,ymin=0, ymax=2, color='r')
            low += 1
            up += 1
            k += 1
        f.show()

    def zoom_logplot(self, f, j_0, j_min, j_max, y, wvar, ci):
        tau1 = []
        tau2 = []
        low = []
        up = []
        for i in range(0,j_0):
            tau1.append(i+1)
            low.append(wvar[i]-ci[i])
            up.append(wvar[i]-ci[i])
        for j in range(j_min-1, j_max):
            tau2.append(j+1)
        a = f.add_subplot(111)    
        a.plot(tau2, y)
        a.errorbar(tau1, wvar, yerr=[low, up], fmt='o')
        a.set_xlim(0, j_0+1)
        a.set_xlabel("Scale")
        a.set_ylabel("Wavelet variance")
        f.show()

    def zoom_ihurst(self, f, x, y, data):
        avg = data[0]
        min = avg-data[2]
        max = avg+data[2]
        a = f.add_subplot(111)
        a.plot(x, y)
        a.axhline(y=avg, color='r')
        a.axhline(y=min, color='g')
        a.axhline(y=max, color='g')
        a.set_xlabel("Time (min)")
        a.set_ylabel("Instant Hurst estimate")
        f.show()

    def show_signal(self, parent, x, y, et):
        ep1_start = et['1'][1]/float(self.genset.f*60)
        ep1_end = ep1_start+10
        ep2_start = et['2'][1]/float(self.genset.f*60)
        ep2_end = ep2_start+10
        ep3_start = et['3'][1]/float(self.genset.f*60)
        ep3_end = ep3_start+10
        f = mpl.figure(figsize=(10,5),dpi=100) 
        a = f.add_subplot(111)
        a.plot(x, y)
        a.axvline(x=ep1_start, color='black')
        a.axvline(x=ep1_end, color='r')
        a.axvline(x=ep2_start, color='black')
        a.axvline(x=ep2_end, color='r')
        a.axvline(x=ep3_start, color='black')
        a.axvline(x=ep3_end, color='r')
        if self.genset.filetype is 0:
            ep4_start = et['4'][1]/float(self.genset.f*60)
            ep4_end = ep4_start+10
            a.axvline(x=ep4_start, color='black')
            a.axvline(x=ep4_end, color='r')
        a.set_xlabel("Time (min)")
        a.set_ylabel("Conductance (G)")
        canvas = FigureCanvasTkAgg(f, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()

    def show_graph(self, parent, x, y):
        f = mpl.figure(figsize=(10,5),dpi=100) 
        a = f.add_subplot(111)
        a.plot(x, y)
        a.set_xlabel("Time (min)")
        a.set_ylabel("Conductance (G)")
        canvas = FigureCanvasTkAgg(f, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()

    def show_igraph(self, parent, x, y, data):
        avg = data[0]
        min = avg-data[2]
        max = avg+data[2]
        f = mpl.figure(figsize=(10,5),dpi=100)
        a = f.add_subplot(111)
        a.plot(x, y)
        a.axhline(y=avg, color='r')
        a.axhline(y=min, color='g')
        a.axhline(y=max, color='g')
        a.set_xlabel("Time (min)")
        a.set_ylabel("Instant Hurst estimate")
        canvas = FigureCanvasTkAgg(f, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()

    def multiple_graphs(self, parent, x, y, boundaries):
        f = mpl.figure(figsize=(10,5),dpi=100)
        i = len(y)
        k = 1
        low = 0
        up = len(boundaries)/2
        for list in y:
            a = f.add_subplot(i,1,k)
            a.set_axis_off()
            a.plot(x, list)
            lowvalue = boundaries[low]/float(self.genset.f*60)
            upvalue = boundaries[up]/float(self.genset.f*60)
            a.axvline(x=lowvalue,ymin=0, ymax=2, color='r')
            a.axvline(x=upvalue,ymin=0, ymax=2, color='r')
            low += 1
            up += 1
            k += 1
        canvas = FigureCanvasTkAgg(f, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()

    def logplot(self, parent, j_0, j_min, j_max, y, wvar, ci):
        f = mpl.figure(figsize=(10,5),dpi=100)
        a = f.add_subplot(111)
        tau1 = []
        tau2 = []
        low = []
        up = []
        for i in range(0,j_0):
            tau1.append(i+1)
            low.append(wvar[i]-ci[i])
            up.append(wvar[i]-ci[i])
        for j in range(j_min-1, j_max):
            tau2.append(j+1)
        print "---"
        print len(tau2), len(y)
        print len(tau1), len(wvar)
        a.plot(tau2, y)
        a.errorbar(tau1, wvar, yerr=[low, up], fmt='o')
        a.set_xlim(0, j_0+1)
        a.set_xlabel("Scale")
        a.set_ylabel("Wavelet variance")
        canvas = FigureCanvasTkAgg(f, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()

    def regraph(self):
        activepage = self.graphs.index(self.graphs.getcurselection())
        if self.activeepoch.get() is 1:
            if "Epoch 1" in self.graphs.pagenames():
                self.graphs.delete("Epoch 1")
                self.epochplot = self.graphs.add("Epoch 1")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
            elif "Epoch 2" in self.graphs.pagenames():
                self.graphs.delete("Epoch 2")
                self.epochplot = self.graphs.add("Epoch 1")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
            elif "Epoch 3" in self.graphs.pagenames():
                self.graphs.delete("Epoch 3")
                self.epochplot = self.graphs.add("Epoch 1")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
            elif "Epoch 4" in self.graphs.pagenames():
                self.graphs.delete("Epoch 4")
                self.epochplot = self.graphs.add("Epoch 1")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e1)
            if "WA" in self.graphs.pagenames():
                self.graphs.delete("WA")
                self.waplot = self.graphs.add("WA")
                self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv1, self.wadata.boundaries)
            if "MRA" in self.graphs.pagenames():
                self.graphs.delete("MRA")
                self.mraplot = self.graphs.add("MRA")
                self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds1, self.mradata.boundaries)
            if "aHurst" in self.graphs.pagenames():
                self.graphs.delete("aHurst")
                self.ahurstplot = self.graphs.add("aHurst")
                self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                         self.ahurstdata.y_hat1, self.ahurstdata.y1, self.ahurstdata.ci1)
            if "iHurst" in self.graphs.pagenames():
                self.graphs.delete("iHurst")
                self.ihurstplot = self.graphs.add("iHurst")
                self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst1, self.ihurstdata.out1)
        elif self.activeepoch.get() is 2:
            if "Epoch 1" in self.graphs.pagenames():
                self.graphs.delete("Epoch 1")
                self.epochplot = self.graphs.add("Epoch 2")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e2)
            elif "Epoch 2" in self.graphs.pagenames():
                self.graphs.delete("Epoch 2")
                self.epochplot = self.graphs.add("Epoch 2")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e2)
            elif "Epoch 3" in self.graphs.pagenames():
                self.graphs.delete("Epoch 3")
                self.epochplot = self.graphs.add("Epoch 2")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e2)
            elif "Epoch 4" in self.graphs.pagenames():
                self.graphs.delete("Epoch 4")
                self.epochplot = self.graphs.add("Epoch 2")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e2)
            if "WA" in self.graphs.pagenames():
                self.graphs.delete("WA")
                self.waplot = self.graphs.add("WA")
                self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv2, self.wadata.boundaries)
            if "MRA" in self.graphs.pagenames():
                self.graphs.delete("MRA")
                self.mraplot = self.graphs.add("MRA")
                self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds2, self.mradata.boundaries)
            if "aHurst" in self.graphs.pagenames():
                self.graphs.delete("aHurst")
                self.ahurstplot = self.graphs.add("aHurst")
                self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                             self.ahurstdata.y_hat2, self.ahurstdata.y2, self.ahurstdata.ci2)
            if "iHurst" in self.graphs.pagenames():
                self.graphs.delete("iHurst")
                self.ihurstplot = self.graphs.add("iHurst")
                self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst2, self.ihurstdata.out2)
        elif self.activeepoch.get() is 3:
            if "Epoch 1" in self.graphs.pagenames():
                self.graphs.delete("Epoch 1")
                self.epochplot = self.graphs.add("Epoch 3")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e3)
            elif "Epoch 2" in self.graphs.pagenames():
                self.graphs.delete("Epoch 2")
                self.epochplot = self.graphs.add("Epoch 3")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e3)
            elif "Epoch 3" in self.graphs.pagenames():
                self.graphs.delete("Epoch 3")
                self.epochplot = self.graphs.add("Epoch 3")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e3)
            elif "Epoch 4" in self.graphs.pagenames():
                self.graphs.delete("Epoch 4")
                self.epochplot = self.graphs.add("Epoch 3")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e3)
            if "WA" in self.graphs.pagenames():
                self.graphs.delete("WA")
                self.waplot = self.graphs.add("WA")
                self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv3, self.wadata.boundaries)
            if "MRA" in self.graphs.pagenames():
                self.graphs.delete("MRA")
                self.mraplot = self.graphs.add("MRA")
                self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds3, self.mradata.boundaries)
            if "aHurst" in self.graphs.pagenames():
                self.graphs.delete("aHurst")
                self.ahurstplot = self.graphs.add("aHurst")
                self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                             self.ahurstdata.y_hat3, self.ahurstdata.y3, self.ahurstdata.ci3)
            if "iHurst" in self.graphs.pagenames():
                self.graphs.delete("iHurst")
                self.ihurstplot = self.graphs.add("iHurst")
                self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst3, self.ihurstdata.out3)
        elif self.activeepoch.get() is 4:
            if "Epoch 1" in self.graphs.pagenames():
                self.graphs.delete("Epoch 1")
                self.epochplot = self.graphs.add("Epoch 4")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e4)
            elif "Epoch 2" in self.graphs.pagenames():
                self.graphs.delete("Epoch 2")
                self.epochplot = self.graphs.add("Epoch 4")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e4)
            elif "Epoch 3" in self.graphs.pagenames():
                self.graphs.delete("Epoch 3")
                self.epochplot = self.graphs.add("Epoch 4")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e4)
            elif "Epoch 4" in self.graphs.pagenames():
                self.graphs.delete("Epoch 4")
                self.epochplot = self.graphs.add("Epoch 4")
                self.show_graph(self.epochplot, self.signal.etime, self.signal.e4)
            if "WA" in self.graphs.pagenames():
                self.graphs.delete("WA")
                self.waplot = self.graphs.add("WA")
                self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv4, self.wadata.boundaries)
            if "MRA" in self.graphs.pagenames():
                self.graphs.delete("MRA")
                self.mraplot = self.graphs.add("MRA")
                self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds4, self.mradata.boundaries)
            if "aHurst" in self.graphs.pagenames():
                self.graphs.delete("aHurst")
                self.ahurstplot = self.graphs.add("aHurst")
                self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                             self.ahurstdata.y_hat4, self.ahurstdata.y4, self.ahurstdata.ci4)
            if "iHurst" in self.graphs.pagenames():
                self.graphs.delete("iHurst")
                self.ihurstplot = self.graphs.add("iHurst")
                self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst4, self.ihurstdata.out4)
        self.graphs.selectpage(activepage)

    def calculate(self):
        #master.shurst(self.signal)
        #If file is opened first destroy existing stuff...(and set Data to zero)
        try:
            if self.waset.active is 1:
                master.wa(self.waset, self.signal, self.wadata)
                self.waplot = self.graphs.add("WA")
                if self.activeepoch.get() is 1:
                    self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv1, self.wadata.boundaries)
                elif self.activeepoch.get() is 2:
                    self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv2, self.wadata.boundaries)
                elif self.activeepoch.get() is 3:
                    self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv3, self.wadata.boundaries)
                elif self.activeepoch.get() is 4:
                    self.multiple_graphs(self.waplot, self.signal.etime, self.wadata.wv4, self.wadata.boundaries)
            if self.mraset.active is 1:
                master.mra(self.mraset, self.signal, self.mradata)
                self.mraplot = self.graphs.add("MRA")
                if self.activeepoch.get() is 1:
                    self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds1, self.mradata.boundaries)
                elif self.activeepoch.get() is 2:
                    self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds2, self.mradata.boundaries)
                elif self.activeepoch.get() is 3:
                    self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds3, self.mradata.boundaries)
                elif self.activeepoch.get() is 5:
                    self.multiple_graphs(self.mraplot, self.signal.etime, self.mradata.ds4, self.mradata.boundaries)  
            if self.ahurstset.active is 1:
                master.ahurst(self.ahurstset, self.signal, self.ahurstdata)
                self.ahurstplot = self.graphs.add("aHurst")
                if self.activeepoch.get() is 1:
                    self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                                 self.ahurstdata.y_hat1, self.ahurstdata.y1, self.ahurstdata.ci1)
                elif self.activeepoch.get() is 2:
                    self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                                 self.ahurstdata.y_hat2, self.ahurstdata.y2, self.ahurstdata.ci2)
                elif self.activeepoch.get() is 3:
                    self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                                 self.ahurstdata.y_hat3, self.ahurstdata.y3, self.ahurstdata.ci3)
                elif self.activeepoch.get() is 4:
                    self.logplot(self.ahurstplot, self.ahurstset.j_0, self.ahurstset.j_min, self.ahurstset.j_max,
                                 self.ahurstdata.y_hat4, self.ahurstdata.y4, self.ahurstdata.ci4)
                
            if self.ihurstset.active is 1:
                master.ihurst(self.ihurstset, self.signal, self.ihurstdata)
                self.ihurstplot = self.graphs.add("iHurst")
                if self.activeepoch.get() is 1:
                    self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst1, self.ihurstdata.out1)
                elif self.activeepoch.get() is 2:
                    self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst2, self.ihurstdata.out2)
                elif self.activeepoch.get() is 3:
                    self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst3, self.ihurstdata.out3)
                elif self.activeepoch.get() is 4:
                    self.show_igraph(self.ihurstplot, self.signal.etime, self.ihurstdata.hurst4, self.ihurstdata.out4)
        except:
            tkMessageBox.showerror("Run Error", "Cannot calculate!")
        self.show_results()

    def quit(self):
        self.master.quit()
        
    def show_results(self):
        if self.waset.active is 1:
            if self.waset.transtype is 0:
                algo = "DWT"
            elif self.waset.transtype is 1:
                algo = "MODWT"
            if self.waset.type is 1:
                wavelet = "Daubechies"
            elif self.waset.type is 2:
                wavelet = "Least Asymmetrical"
            elif self.waset.type is 3:
                wavelet = "Best Localized"
            elif self.waset.type is 4:
                wavelet = "Coiflet"
            if self.waset.refl is 0:
                boundary = "Circular"
            elif self.waset.refl is 1:
                boundary = "Reflection"
            text="WA settings\nAlgorithm: %s\nWavelet: %s %d\nBoundary: %s\nLevel: %d\n"%(algo, wavelet, self.waset.l, boundary, self.waset.j_0)
            self.waoutgeneral = Label(self.right_frame, text=text, justify=LEFT, bg="red")
            self.waoutgeneral.pack(side=TOP, anchor=NW, padx=15)
        if self.mraset.active is 1:
            if self.mraset.transtype is 0:
                algo = "DWT"
            elif self.mraset.transtype is 1:
                algo = "MODWT"
            if self.mraset.type is 1:
                wavelet = "Daubechies"
            elif self.mraset.type is 2:
                wavelet = "Least Asymmetrical"
            elif self.mraset.type is 3:
                wavelet = "Best Localized"
            elif self.mraset.type is 4:
                wavelet = "Coiflet"
            if self.mraset.refl is 0:
                boundary = "Circular"
            elif self.mraset.refl is 1:
                boundary = "Reflection"
            text="MRA settings\nAlgorithm: %s\nWavelet: %s %d\nBoundary: %s\nLevel: %d\n"%(algo, wavelet, self.mraset.l, boundary, self.mraset.j_0)
            self.mraoutgeneral = Label(self.right_frame, text=text, justify=LEFT, bg="red")
            self.mraoutgeneral.pack(side=TOP, anchor=NW, padx=15)
        if self.ahurstset.active is 1:
            if self.ahurstset.transtype is 0:
                algo = "DWT"
            elif self.ahurstset.transtype is 1:
                algo = "MODWT"
            if self.ahurstset.type is 1:
                wavelet = "Daubechies"
            elif self.ahurstset.type is 2:
                wavelet = "Least Asymmetrical"
            elif self.ahurstset.type is 3:
                wavelet = "Best Localized"
            elif self.ahurstset.type is 4:
                wavelet = "Coiflet"
            if self.ahurstset.refl is 0:
                boundary = "Circular"
            elif self.ahurstset.refl is 1:
                boundary = "Reflection"
            if self.ahurstset.bias is 0:
                bias = "Unbiased"
            elif self.ahurstset.bias is 1:
                bias = "Biased"
            text="Averaged estimates\nAlgorithm: %s\nWavelet: %s %d\nBoundary: %s\nBias: %s\nWLSE: %d,%d\n"%(algo, wavelet, self.ahurstset.l, boundary, bias, self.ahurstset.j_min, self.ahurstset.j_max)
            self.ahurstoutgeneral = Label(self.right_frame, text=text, justify=LEFT, bg="red")
            self.ahurstoutgeneral.pack(side=TOP, anchor=NW, padx=15)
            ahurstvalue1 = self.ahurstdata.out1[0]
            ahurstvar1 = self.ahurstdata.out1[1]
            ahurstsd1 = self.ahurstdata.out1[2]
            text = "1:\t %.5f \t %.5f \t %.5f"%(ahurstvalue1,ahurstvar1,ahurstsd1)
            self.showahurstout1 = Label(self.right_frame, text=text, bg="red")
            self.showahurstout1.pack(side=TOP, anchor=NW, padx=15)
            ahurstvalue2 = self.ahurstdata.out2[0]
            ahurstvar2 = self.ahurstdata.out2[1]
            ahurstsd2 = self.ahurstdata.out2[2]
            text = "2:\t %.5f \t %.5f \t %.5f"%(ahurstvalue2,ahurstvar2,ahurstsd2)
            self.showahurstout2 = Label(self.right_frame, text=text, bg="red")
            self.showahurstout2.pack(side=TOP, anchor=NW, padx=15)
            ahurstvalue3 = self.ahurstdata.out3[0]
            ahurstvar3 = self.ahurstdata.out3[1]
            ahurstsd3 = self.ahurstdata.out3[2]
            if self.genset.filetype is 0:
                text = "3:\t %.5f \t %.5f \t %.5f"%(ahurstvalue3,ahurstvar3,ahurstsd3)
            elif self.genset.filetype is 1:
                text = "3:\t %.5f \t %.5f \t %.5f\n"%(ahurstvalue3,ahurstvar3,ahurstsd3)
            self.showahurstout3 = Label(self.right_frame, text=text, bg="red")
            self.showahurstout3.pack(side=TOP, anchor=NW, padx=15)
            if self.genset.filetype is 0:
                ahurstvalue4 = self.ahurstdata.out4[0]
                ahurstvar4 = self.ahurstdata.out4[1]
                ahurstsd4 = self.ahurstdata.out4[2]
                text = "4:\t %.5f \t %.5f \t %.5f\n"%(ahurstvalue4,ahurstvar4,ahurstsd4)
                self.showahurstout4 = Label(self.right_frame, text=text, bg="red")
                self.showahurstout4.pack(side=TOP, anchor=NW, padx=15)

        if self.ihurstset.active is 1:
            if self.ihurstset.type is 1:
                wavelet = "Daubechies"
            elif self.ihurstset.type is 2:
                wavelet = "Least Asymmetrical"
            elif self.ihurstset.type is 3:
                wavelet = "Best Localized"
            elif self.ihurstset.type is 4:
                wavelet = "Coiflet"
            if self.ihurstset.refl is 0:
                boundary = "Circular"
            elif self.ihurstset.refl is 1:
                boundary = "Reflection"
            if self.ihurstset.bias is 0:
                bias = "Unbiased"
            elif self.ihurstset.bias is 1:
                bias = "Biased"
            text="Instant estimates\nAlgorithm: MODWT\nWavelet: %s %d\nBoundary: %s\nBias: %s\niWLSE: %d,%d\n"%(wavelet, self.ihurstset.l, boundary, bias, self.ihurstset.j_min, self.ihurstset.j_max)
            self.ihurstoutgeneral = Label(self.right_frame, text=text, justify=LEFT, bg="red")
            self.ihurstoutgeneral.pack(side=TOP, anchor=NW, padx=15)
            ihurstvalue1 = self.ihurstdata.out1[0]
            ihurstvar1 = self.ihurstdata.out1[1]
            ihurstsd1 = self.ihurstdata.out1[2]
            text = "1:\t %.5f \t %.5f \t %.5f"%(ihurstvalue1,ihurstvar1,ihurstsd1)
            self.showihurstout1 = Label(self.right_frame, text=text, bg="red")
            self.showihurstout1.pack(side=TOP, anchor=NW, padx=15)
            ihurstvalue2 = self.ihurstdata.out2[0]
            ihurstvar2 = self.ihurstdata.out2[1]
            ihurstsd2 = self.ihurstdata.out2[2]
            text = "2:\t %.5f \t %.5f \t %.5f"%(ihurstvalue2,ihurstvar2,ihurstsd2)
            self.showihurstout2 = Label(self.right_frame, text=text, bg="red")
            self.showihurstout2.pack(side=TOP, anchor=NW, padx=15)
            ihurstvalue3 = self.ihurstdata.out3[0]
            ihurstvar3 = self.ihurstdata.out3[1]
            ihurstsd3 = self.ihurstdata.out3[2]
            text = "3:\t %.5f \t %.5f \t %.5f"%(ihurstvalue3,ihurstvar3,ihurstsd3)
            self.showihurstout3 = Label(self.right_frame, text=text, bg="red")
            self.showihurstout3.pack(side=TOP, anchor=NW, padx=15)
            if self.genset.filetype is 0:
                ihurstvalue4 = self.ihurstdata.out4[0]
                ihurstvar4 = self.ihurstdata.out4[1]
                ihurstsd4 = self.ihurstdata.out4[2]
                text = "4:\t %.5f \t %.5f \t %.5f"%(ihurstvalue4,ihurstvar4,ihurstsd4)
                self.showihurstout4 = Label(self.right_frame, text=text, bg="red")
                self.showihurstout4.pack(side=TOP, anchor=NW, padx=15)

if __name__ == "__main__": 
    root = Tk()
    Pmw.initialise(root)
    Gui(root)
    root.mainloop()

#Update epoch buttons...
#Dimensionerror aHurst regraph

#Finish off ihurst function in c++
#Check aHurst DWT results if I initialise to 0 first

from tkinter import *
import numpy as np
from os import system
from tkinter.filedialog import askopenfilename

#from tkFileDialog import askopenfilename

def readinput(input_obj,main_frm):
    path = askopenfilename()
    if path == '':
        return
    input_obj.destroy()
    input_obj = input_frames(main_frm)

    input_obj.pump_amp,input_obj.inj_time,input_obj.seed_amp,\
    input_obj.seed_time,input_obj.seed_delay,input_obj.seed_freq,input_obj.seed_profile,\
        input_obj.dens,input_obj.gamma,input_obj.beam_tempX,input_obj.beam_tempY,input_obj.beam_top,\
        input_obj.beam_left,input_obj.beam_right,input_obj.beam_bottom,input_obj.num_part,\
        input_obj.simul_time,input_obj.steps,input_obj.box_top,input_obj.box_left,input_obj.box_right,\
        input_obj.box_bottom,input_obj.snap_time,input_obj.diag_phs,input_obj.diag_fld,\
    input_obj.diag_dens,input_obj.num_node,input_obj.seed_width,input_obj.estat,\
    input_obj.neutr,input_obj.out_folder = np.loadtxt(path,dtype='str')

    input_obj.simulation_opts()
    input_obj.field_frames()
    input_obj.particles_frames()


class input_frames:
    def __init__(self,parent):
        self.main_frm = parent
        self.save_buttom = Button(parent,text='save input',command=self.save_params)
        self.save_buttom.grid(row=1,column=1)
        self.main_buttom = Button(parent,text='launch ewok',command=self.get_params)
        self.main_buttom.grid(row=1,column=2)
        return

    def save_params(self):
        self.pump_amp =  self.pump_amp_entr.get()
        self.inj_time = self.inj_time_entr.get()
        self.seed_amp =  self.seed_amp_entr.get()
        self.seed_time = self.seed_time_entr.get()
        self.seed_delay = self.seed_delay_entr.get()
        self.seed_freq = self.seed_freq_entr.get()
        self.seed_profile = self.seed_profile_lstbx.curselection()[0]
        self.dens = self.dens_entr.get()
        self.gamma =  self.gamma_entr.get()
        self.beam_tempX = self.beam_tempX_entr.get()
        self.beam_tempY = self.beam_tempY_entr.get()
        self.beam_top = self.beam_top_entr.get()
        self.beam_left = self.beam_left_entr.get()
        self.beam_right = self.beam_right_entr.get()
        self.beam_bottom = self.beam_bottom_entr.get()
        self.num_part = self.num_part_entr.get()
        self.simul_time = self.simul_time_entr.get()
        self.steps = self.steps_entr.get()
        self.box_top = self.box_top_entr.get()
        self.box_left = self.box_left_entr.get()
        self.box_right = self.box_right_entr.get()
        self.box_bottom = self.box_bottom_entr.get()
        self.snap_time = self.snap_time_entr.get()
        self.diag_phs = self.diag_phs_var.get()
        self.diag_fld = self.diag_fld_var.get()
        self.diag_dens = self.diag_dens_var.get()
        self.num_node = self.num_node_entr.get()
        self.seed_width = self.seed_width_entr.get()
        self.out_folder = self.out_folder_entr.get()
        self.estat = self.estat_var.get()
        self.neutr = self.neutr_var.get()

        input = np.array([self.pump_amp,self.inj_time,self.seed_amp,\
        self.seed_time,self.seed_delay,self.seed_freq,self.seed_profile,\
        self.dens,self.gamma,self.beam_tempX,self.beam_tempY,self.beam_top,\
        self.beam_left,self.beam_right,self.beam_bottom,self.num_part,\
        self.simul_time,self.steps,self.box_top,self.box_left,self.box_right,\
        self.box_bottom,self.snap_time,self.diag_phs,self.diag_fld,\
        self.diag_dens,self.num_node,self.seed_width,self.estat,self.neutr,self.out_folder],dtype='str')
        np.savetxt('input.txt',input,fmt="%s")

    def get_params(self):
        self.save_params()
        system('rm -f nohup.out')
        nodes2launch =  str(int(float(self.num_node)))
        system('nohup mpirun -np '+nodes2launch+'python ./ewok_launch.py < /dev/null &')

    def default_input(self):
        self.pump_amp,self.inj_time,self.seed_amp,\
          self.seed_time,self.seed_delay,self.seed_freq,self.seed_profile,\
          self.dens,self.gamma,self.beam_tempX,self.beam_tempY,self.beam_top,\
          self.beam_left,self.beam_right,self.beam_bottom,self.num_part,\
          self.simul_time,self.steps,self.box_top,self.box_left,self.box_right,\
          self.box_bottom,self.snap_time,self.diag_phs,self.diag_fld,self.diag_dens,\
          self.num_node,self.seed_width,self.estat,self.neutr,self.out_folder = np.loadtxt('input.txt',dtype='str')

        self.pump_amp,self.inj_time,self.seed_amp,\
          self.seed_time,self.seed_delay,self.seed_freq,self.seed_profile,\
          self.dens,self.gamma,self.beam_tempX,self.beam_tempY,self.beam_top,\
          self.beam_left,self.beam_right,self.beam_bottom,self.num_part,\
          self.simul_time,self.steps,self.box_top,self.box_left,self.box_right,\
          self.box_bottom,self.snap_time,self.diag_phs,self.diag_fld,self.diag_dens,\
          self.num_node, self.seed_width, self.estat, self.neutr = np.array([\
          self.pump_amp,self.inj_time,self.seed_amp,\
          self.seed_time,self.seed_delay,self.seed_freq,self.seed_profile,\
          self.dens,self.gamma,self.beam_tempX,self.beam_tempY,self.beam_top,\
          self.beam_left,self.beam_right,self.beam_bottom,self.num_part,\
          self.simul_time,self.steps,self.box_top,self.box_left,self.box_right,\
          self.box_bottom,self.snap_time,self.diag_phs,self.diag_fld,self.diag_dens,\
          self.num_node,self.seed_width,self.estat,self.neutr] ,dtype='float')


    def simulation_opts(self):
        self.sim_n_diag_frm = Frame(self.main_frm)

        simul_frm = Frame(self.sim_n_diag_frm,relief=SUNKEN,borderwidth=1)

        num_node_lab = Label(simul_frm,text='number of nodes')
        num_node_lab.grid(row=0,column=0,sticky=N)
        self.num_node_entr = Entry(simul_frm,width=6,justify=CENTER)
        self.num_node_entr.insert(0,self.num_node)
        self.num_node_entr.grid(row=1,column=0,sticky=N)

        num_part_lab = Label(simul_frm,text='particles per cell') 
        num_part_lab.grid(row=2,column=0,sticky=N)
        self.num_part_entr = Entry(simul_frm,width=6,justify=CENTER) 
        self.num_part_entr.insert(0,self.num_part) 
        self.num_part_entr.grid(row=3,column=0,sticky=N)

        simul_time_lab = Label(simul_frm,text='simulation time') 
        simul_time_lab.grid(row=4,column=0,sticky=N)
        self.simul_time_entr = Entry(simul_frm,width=6,justify=CENTER) 
        self.simul_time_entr.insert(0,self.simul_time) 
        self.simul_time_entr.grid(row=5,column=0,sticky=N)

        steps_lab = Label(simul_frm,text='steps per laser period') 
        steps_lab.grid(row=6,column=0,sticky=N)
        self.steps_entr = Entry(simul_frm,width=6,justify=CENTER) 
        self.steps_entr.insert(0,self.steps) 
        self.steps_entr.grid(row=7,column=0,sticky=N)

        simul_frm.grid(row=0,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        box_geom_frm = Frame(self.sim_n_diag_frm,relief=SUNKEN,borderwidth=1)
        box_geom_lab = Label(box_geom_frm,text='box geometry') 
        box_geom_lab.grid(row=0,column=0,sticky=N, columnspan=2)

        box_top_lab = Label(box_geom_frm,text='top') 
        box_top_lab.grid(row=1,column=0,sticky=N, columnspan=2)
        self.box_top_entr = Entry(box_geom_frm,width=6,justify=CENTER) 
        self.box_top_entr.insert(0,self.box_top) 
        self.box_top_entr.grid(row=2,column=0,sticky=N, columnspan=2)

        box_left_lab = Label(box_geom_frm,text='left') 
        box_left_lab.grid(row=3,column=0,sticky=N)
        self.box_left_entr = Entry(box_geom_frm,width=6,justify=CENTER) 
        self.box_left_entr.insert(0,self.box_left) 
        self.box_left_entr.grid(row=4,column=0,sticky=N)

        box_right_lab = Label(box_geom_frm,text='right') 
        box_right_lab.grid(row=3,column=1,sticky=N)
        self.box_right_entr = Entry(box_geom_frm,width=6,justify=CENTER) 
        self.box_right_entr.insert(0,self.box_right) 
        self.box_right_entr.grid(row=4,column=1,sticky=N)

        box_bottom_lab = Label(box_geom_frm,text='bottom') 
        box_bottom_lab.grid(row=5,column=0,sticky=N, columnspan=2)
        self.box_bottom_entr = Entry(box_geom_frm,width=6,justify=CENTER) 
        self.box_bottom_entr.insert(0,self.box_bottom) 
        self.box_bottom_entr.grid(row=6,column=0,sticky=N, columnspan=2)

        box_geom_frm.grid(row=1,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        self.sim_n_diag_frm.grid(row=0,column=2,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)


    def field_frames(self):

        self.pump_diag_frm = Frame(self.main_frm)
        
        pump_frm = Frame(self.pump_diag_frm,relief=SUNKEN,borderwidth=1)
        
        pump_amp_lab = Label(pump_frm,text='pump amplitude')
        pump_amp_lab.grid(row=0,column=0,sticky=N)
        self.pump_amp_entr = Entry(pump_frm,width=12,justify=CENTER)
        self.pump_amp_entr.insert(0,self.pump_amp)
        self.pump_amp_entr.grid(row=1,column=0,sticky=N)

        inj_time_lab = Label(pump_frm,text='injection time')
        inj_time_lab.grid(row=2,column=0,sticky=N)
        self.inj_time_entr = Entry(pump_frm,width=12,justify=CENTER) 
        self.inj_time_entr.insert(0,self.inj_time) 
        self.inj_time_entr.grid(row=3,column=0,sticky=N)

        gamma_lab = Label(pump_frm,text='lattice angle')
        gamma_lab.grid(row=4,column=0,sticky=N)
        self.gamma_entr = Entry(pump_frm,width=6,justify=CENTER)
        self.gamma_entr.insert(0,self.gamma)
        self.gamma_entr.grid(row=5,column=0,sticky=N)

        pump_frm.grid(row=0,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        self.diag_frm = Frame(self.pump_diag_frm,relief=SUNKEN,borderwidth=1)
        snap_time_lab = Label(self.diag_frm,text='snaps interval') 
        snap_time_lab.grid(row=0,column=0,sticky=N)

        self.snap_time_entr = Entry(self.diag_frm,width=6,justify=CENTER) 
        self.snap_time_entr.insert(0,self.snap_time) 
        self.snap_time_entr.grid(row=1,column=0,sticky=N)

        diag_lab = Label(self.diag_frm,text='diagnostics') 
        diag_lab.grid(row=2,column=0,sticky=N)

        self.diag_phs_var = IntVar() 
        self.diag_phs_var.set(int(self.diag_phs))
        diag_phs_chkbut= Checkbutton(self.diag_frm,text='phase',variable=self.diag_phs_var) 
        diag_phs_chkbut.grid(row=3,column=0,sticky=NW)

        self.diag_fld_var = IntVar()
        self.diag_fld_var.set(int(self.diag_fld))
        diag_fld_chkbut= Checkbutton(self.diag_frm,text='field',variable=self.diag_fld_var) 
        diag_fld_chkbut.grid(row=4,column=0,sticky=NW)

        self.diag_dens_var = IntVar() 
        self.diag_dens_var.set(int(self.diag_dens)) 
        diag_dens_chkbut= Checkbutton(self.diag_frm,text='density',variable=self.diag_dens_var) 
        diag_dens_chkbut.grid(row=5,column=0,sticky=NW)

        self.diag_frm.grid(row=1,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)


        self.other_opt_frm = Frame(self.pump_diag_frm,relief=SUNKEN,borderwidth=1)

        self.estat_var = IntVar()
        self.estat_var.set(int(self.estat))
        estat_chkbut= Checkbutton(self.other_opt_frm,text='electrostatics',variable=self.estat_var)
        estat_chkbut.grid(row=0,column=0,sticky=NW)

        self.neutr_var = IntVar()
        self.neutr_var.set(int(self.neutr))
        neutr_chkbut= Checkbutton(self.other_opt_frm,text='neutralize',variable=self.neutr_var)
        neutr_chkbut.grid(row=1,column=0,sticky=NW)

#                out_folder_lab = Label(self.other_opt_frm,text='out folder')
#                out_folder_lab.grid(row=2,column=0,sticky=N)
#                self.out_folder_entr = Entry(self.other_opt_frm,width=15,justify=CENTER)
#                self.out_folder_entr.insert(0,self.out_folder)
#                self.out_folder_entr.grid(row=3,column=0,sticky=N)

#                out_folder_lab = Label(self.main_frm,text='out folder')
#                out_folder_lab.grid(row=1,column=3,sticky=N,columnspan=2)

        self.out_folder_entr = Entry(self.main_frm,justify=CENTER)
        self.out_folder_entr.insert(0,self.out_folder)
        self.out_folder_entr.grid(row=1,column=3,sticky=EW,columnspan=3)

        self.other_opt_frm.grid(row=2,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        self.pump_diag_frm.grid(row=0,column=3, ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)


        self.seed_frm = Frame(self.main_frm,relief=SUNKEN,borderwidth=1)

        seed_amp_lab = Label(self.seed_frm,text='signal amplitude') 
        seed_amp_lab.grid(row=0,column=0,sticky=N)
        self.seed_amp_entr = Entry(self.seed_frm,width=12,justify=CENTER) 
        self.seed_amp_entr.insert(0,self.seed_amp) 
        self.seed_amp_entr.grid(row=1,column=0,sticky=N)

        seed_time_lab = Label(self.seed_frm,text='signal duration') 
        seed_time_lab.grid(row=2,column=0,sticky=N)
        self.seed_time_entr = Entry(self.seed_frm,width=12,justify=CENTER) 
        self.seed_time_entr.insert(0,self.seed_time) 
        self.seed_time_entr.grid(row=3,column=0,sticky=N)

        seed_delay_lab = Label(self.seed_frm,text='signal delay') 
        seed_delay_lab.grid(row=4,column=0,sticky=N)
        self.seed_delay_entr = Entry(self.seed_frm,width=12,justify=CENTER)
        self.seed_delay_entr.insert(0,self.seed_delay)
        self.seed_delay_entr.grid(row=5,column=0,sticky=N)

        seed_freq_lab = Label(self.seed_frm,text='signal frequency')
        seed_freq_lab.grid(row=6,column=0,sticky=N)
        self.seed_freq_entr = Entry(self.seed_frm,width=12,justify=CENTER) 
        self.seed_freq_entr.insert(0,self.seed_freq) 
        self.seed_freq_entr.grid(row=7,column=0,sticky=N)

        seed_profile_lab = Label(self.seed_frm,text='signal profile') 
        seed_profile_lab.grid(row=8,column=0,sticky=N)
        self.seed_profile_lstbx = Listbox(self.seed_frm,exportselection=0,height=3,width=12) 
        self.seed_profile_lstbx.insert(END,'cosinus') 
        self.seed_profile_lstbx.insert(END,'flat') 
        self.seed_profile_lstbx.insert(END,'gauss') 
        self.seed_profile_lstbx.selection_set(int(self.seed_profile)) 
        self.seed_profile_lstbx.grid(row=9,column=0,sticky=N)

        seed_width_lab = Label(self.seed_frm,text='profile width')
        seed_width_lab.grid(row=10,column=0,sticky=N)
        self.seed_width_entr = Entry(self.seed_frm,width=12,justify=CENTER)
        self.seed_width_entr.insert(0,self.seed_width)
        self.seed_width_entr.grid(row=11,column=0,sticky=N)

        self.seed_frm.grid(row=0,column=4,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

    def particles_frames(self):

        self.eons_frm = Frame(self.main_frm)

        eons_gen_frm = Frame(self.eons_frm,relief=SUNKEN,borderwidth=1)

        dens_lab = Label(eons_gen_frm,text='electron density') 
        dens_lab.grid(row=2,column=0,sticky=N)
        self.dens_entr = Entry(eons_gen_frm,width=6,justify=CENTER) 
        self.dens_entr.insert(0,self.dens) 
        self.dens_entr.grid(row=3,column=0,sticky=N)

        eons_gen_frm.grid(row=0,column=0, ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        beam_temps_frm = Frame(self.eons_frm,relief=SUNKEN,borderwidth=1)

        beam_temps_lab = Label(beam_temps_frm,text='thermal velocities') 
        beam_temps_lab.grid(row=0,column=0,sticky=N, columnspan=2)
    
        beam_tempX_lab = Label(beam_temps_frm,text='Vx') 
        beam_tempX_lab.grid(row=1,column=0,sticky=N)
        self.beam_tempX_entr = Entry(beam_temps_frm,width=6,justify=CENTER) 
        self.beam_tempX_entr.insert(0,self.beam_tempX) 
        self.beam_tempX_entr.grid(row=2,column=0,sticky=N)

        beam_tempY_lab = Label(beam_temps_frm,text='Vy') 
        beam_tempY_lab.grid(row=1,column=1,sticky=N)
        self.beam_tempY_entr = Entry(beam_temps_frm,width=6,justify=CENTER) 
        self.beam_tempY_entr.insert(0,self.beam_tempY) 
        self.beam_tempY_entr.grid(row=2,column=1,sticky=N)

        beam_temps_frm.grid(row=1,column=0, ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        beam_geom_frm = Frame(self.eons_frm,relief=SUNKEN,borderwidth=1)

        beam_geom_lab = Label(beam_geom_frm,text='beam geometry') 
        beam_geom_lab.grid(row=0,column=0,sticky=N, columnspan=2)

        beam_top_lab = Label(beam_geom_frm,text='top') 
        beam_top_lab.grid(row=1,column=0,sticky=N, columnspan=2)
        self.beam_top_entr = Entry(beam_geom_frm,width=6,justify=CENTER) 
        self.beam_top_entr.insert(0,self.beam_top) 
        self.beam_top_entr.grid(row=2,column=0,sticky=N, columnspan=2)

        beam_left_lab = Label(beam_geom_frm,text='left') 
        beam_left_lab.grid(row=3,column=0,sticky=N)
        self.beam_left_entr = Entry(beam_geom_frm,width=6,justify=CENTER) 
        self.beam_left_entr.insert(0,self.beam_left) 
        self.beam_left_entr.grid(row=4,column=0,sticky=N)

        beam_right_lab = Label(beam_geom_frm,text='right') 
        beam_right_lab.grid(row=3,column=1,sticky=N)
        self.beam_right_entr = Entry(beam_geom_frm,width=6,justify=CENTER) 
        self.beam_right_entr.insert(0,self.beam_right) 
        self.beam_right_entr.grid(row=4,column=1,sticky=N)

        beam_bottom_lab = Label(beam_geom_frm,text='bottom') 
        beam_bottom_lab.grid(row=5,column=0,sticky=N, columnspan=2)
        self.beam_bottom_entr = Entry(beam_geom_frm,width=6,justify=CENTER) 
        self.beam_bottom_entr.insert(0,self.beam_bottom) 
        self.beam_bottom_entr.grid(row=6,column=0,sticky=N, columnspan=2)

        beam_geom_frm.grid(row=2,column=0,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

        self.eons_frm.grid(row=0,column=5,ipadx=3,ipady=3,padx=3,pady=3,sticky=NW)

    def destroy(self):
        self.diag_frm.destroy()
        self.sim_n_diag_frm.destroy()
        self.eons_frm.destroy()
        self.pump_diag_frm.destroy()
        self.seed_frm.destroy()

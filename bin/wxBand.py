#!/usr/bin/env python
import os
import pprint
import random
import wx
from numpy import *
from wx import *
#import lattice as Lat
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib
matplotlib.use('WXAgg')
#from matplotlib.figure import Figure
from matplotlib.pyplot import figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar


class Band:
	def __init__(self,filename,shift=0):
		f = open(filename,"r")

		self.kx = []
		self.ky = []
		self.kz = []
		self.Elist= []

		for line in f.readlines():
			sp = line.split()
			self.kx.append(float(sp[1]))
			self.ky.append(float(sp[2]))
			self.ky.append(float(sp[3]))
			En = []
			for ss in sp[6:-1] :
				En.append(float(ss)-shift)

			self.Elist.append(En)

		self.Elist = mat(self.Elist)

class BarsFrame(wx.Frame):
	title = 'Demo: wxPython with matplotlib'
	def __init__(self):
		wx.Frame.__init__(self, None, -1, self.title)
		
		self.create_menu()
		self.create_status_bar()
		self.create_main_panel()
		self.textbox.SetValue('0 10 -10 0 200 0')
		self.draw_flag, self.Y_high, self.Y_low, self.mu, self.X_high, self.X_low= 0, 10, -10, 0.0, 200, 0
		self.draw_flag=0
		
	def create_menu(self):
		self.menubar = wx.MenuBar()
		
		menu_file = wx.Menu()
		menu_file.AppendSeparator()
		m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
		self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
		
		menu_help = wx.Menu()
		m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
		self.Bind(wx.EVT_MENU, self.on_about, m_about)
		
		self.menubar.Append(menu_file, "&File")
		self.menubar.Append(menu_help, "&Help")
		self.SetMenuBar(self.menubar)
	def create_main_panel(self):
		""" Creates the main panel with all the controls on it:
		     * mpl canvas 
		     * mpl navigation toolbar
		     * Control panel for interaction
		"""
		self.panel = wx.Panel(self)

		self.dirCtrl= wx.GenericDirCtrl( self.panel, size=(500,500),  dir=os.getcwd())
		self.dirCtrl.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSel)
		
		# Create the mpl Figure and FigCanvas objects. 
		# 5x4 inches, 100 dots-per-inch
		#
		self.dpi = 200
		self.fig = figure(figsize=(2.5, 3.0), dpi=self.dpi)
		self.canvas = FigCanvas(self.panel, -1, self.fig)
		

		# Since we have only one plot, we can use add_axes 
		# instead of add_subplot, but then the subplot
		# configuration tool in the navigation toolbar wouldn't
		# work.
		#
		self.fig.subplots_adjust(wspace=0.6)
		from mpl_toolkits.axes_grid import make_axes_locatable
		self.axes11 = self.fig.add_subplot(111)
		
		self.canvas.mpl_connect('pick_event', self.on_pick)

		# Create the navigation toolbar, tied to the canvas
		#
		self.toolbar = NavigationToolbar(self.canvas)
		
		self.textbox = wx.TextCtrl( self.panel, size=(300,-1), style=wx.TE_PROCESS_ENTER)
		self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
		# Layout with box sizers
		#
		flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
		self.vbox = wx.BoxSizer(wx.VERTICAL)
		
		self.vbox = wx.BoxSizer(wx.VERTICAL)
		self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
		self.vbox.Add(self.toolbar, 0, wx.EXPAND)
		self.vbox.Add(self.textbox, 0, wx.EXPAND)
		self.vbox.AddSpacer(10)

		self.hbox = wx.BoxSizer(wx.HORIZONTAL)
		self.vbox.Add(self.hbox, 0, wx.EXPAND)
		self.hbox0 = wx.BoxSizer(wx.HORIZONTAL)


		self.hbox0.Add(self.dirCtrl, 0, wx.EXPAND )
		self.hbox0.Add(self.vbox, 1, wx.RIGHT| wx.TOP | wx.GROW)
		

		#self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
		
		self.panel.SetSizer(self.hbox0)
		self.hbox0.Fit(self)
	def OnSel(self,event):
		path = self.dirCtrl.GetFilePath()

		self.axes11.clear()

		mu = self.mu
		uL = Band(path,mu)
		uL.kx = array(uL.kx)

		
		uEk = uL.Elist.T.tolist()
		_l= len(uEk)/2
		Xpoint = arange(len(uEk[0]));

		self.axes11.axhline(y=0,linestyle='--',color='black')
		#self.axes11.axvline(x=10,linestyle='--',color='gray')
		#self.axes11.axvline(x=20,linestyle='--',color='gray')


		lw=1
		if self.draw_flag == 0:

			for i in xrange(len(uEk)):
				ek = array(uEk[i])
				self.axes11.plot(Xpoint,ek,"-",color=(0,0.1,0.5),linewidth=lw)

		#self.axes11.axis([self.X_low,self.X_high,self.Y_low,self.Y_high])
		#self.axes11.set_xticklabels([],fontsize=30)
		#self.axes11.set_yticklabels([],fontsize=30)

		self.canvas.draw()

	def create_status_bar(self):
		self.statusbar = self.CreateStatusBar()
	def on_cb_grid(self, event):
		pass
	def on_slider_height(self, event):
		pass
	def on_pick(self, event):
		pass
	def on_text_enter(self, event):
		#self.draw_flag, self.Y_high, self.Y_low, self.mu, self.fold = map(float,self.textbox.GetValue().split())
		#self.draw_flag, self.Y_high, self.Y_low, self.mu, self.X_high, self.X_low= map(float,self.textbox.GetValue().split())
		self.OnSel(self)

	def on_exit(self, event):
		self.Destroy()
		#self.FS_frame.Destroy()
	def on_about(self, event):
		msg = """ A demo using wxPython with matplotlib:
		 * Use the matplotlib navigation bar
		 * Add values to the text box and press Enter (or click "Draw!")
		 * Show or hide the grid
		 * Drag the slider to modify the width of the bars
		 * Save the plot to a file using the File menu
		 * Click on a bar to receive an informative message
		"""
		dlg = wx.MessageDialog(self, msg, "About", wx.OK)
		dlg.ShowModal()
		dlg.Destroy()
	def flash_status_message(self, msg, flash_len_ms=1500):
		self.statusbar.SetStatusText(msg)
		self.timeroff = wx.Timer(self)
		self.Bind(
		    wx.EVT_TIMER, 
		    self.on_flash_status_off, 
		    self.timeroff)
		self.timeroff.Start(flash_len_ms, oneShot=True)
	def on_flash_status_off(self, event):
		self.statusbar.SetStatusText('')

if __name__ == '__main__':
	app = wx.PySimpleApp()
	app.frame = BarsFrame()
	app.frame.Show()
	app.MainLoop()

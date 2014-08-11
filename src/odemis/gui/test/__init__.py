# -*- coding: utf-8 -*-
"""
Created on 1 Jul 2013

@author: Éric Piel

Copyright © 2013 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License version 2 as published by the Free Software
Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

"""

# Common configuration and code for the GUI test cases
import random
import os.path
import unittest
import numpy

import wx

import odemis.gui.test.test_gui
import odemis.gui.model as gmodel
import odemis.model as omodel
from odemis.gui.xmlh import odemis_get_test_resources


MANUAL = False
INSPECT = False

SLEEP_TIME = 50 # ms: time to sleep between actions (to slow down the tests)

def goto_manual():
    """ Call this function as soon as possible, to go to manual mode, where
    the test GUI will stay open after finishing the test case. """
    global MANUAL
    MANUAL = True

def goto_inspect():
    global INSPECT
    INSPECT = True

def gui_loop(slp=None):
    """
    Execute the main loop for the GUI until all the current events are processed
    """
    app = wx.GetApp()
    if app is None:
        return

    while True:
        wx.CallAfter(app.ExitMainLoop)
        app.MainLoop()
        if not app.Pending():
            break

    wx.MilliSleep(slp or SLEEP_TIME)

def sleep(ms=None):
    wx.MilliSleep(ms or SLEEP_TIME)

def set_sleep_time(slp_tm):
    global SLEEP_TIME
    SLEEP_TIME = slp_tm

# Default wxPython App that can be used as a basis for testing
class GuiTestApp(wx.App):

    test_frame = None

    def __init__(self, frame):
        odemis.gui.test.test_gui.get_resources = odemis_get_test_resources
        self.test_frame = frame
        self.module_name = ""

        # gen_test_data()
        wx.App.__init__(self, redirect=False)

    def OnInit(self):
        self.test_frame = self.test_frame(None) # odemis.gui.test.test_gui.xrccanvas_frame(None)
        self.test_frame.SetSize((400, 400))
        self.test_frame.Center()
        self.test_frame.Layout()

        import __main__
        self.module_name = os.path.basename(__main__.__file__)
        self.test_frame.SetTitle(self.module_name)

        self.test_frame.Show()

        return True

    def panel_finder(self, win=None):
        """ Find the first child panel of win """

        win = win or self.test_frame

        for c in win.GetChildren():
            if isinstance(c, wx.Panel):
                return c
            else:
                return self.panel_finder(c)
        return None

# TestCase base class, with GuiTestApp support
class GuiTestCase(unittest.TestCase):

    frame_class = None
    app_class = None

    @classmethod
    def setUpClass(cls):
        if not cls.frame_class:
            raise ValueError("No frame_class set!")
        cls.app_class = cls.app_class or GuiTestApp
        cls.app = cls.app_class(cls.frame_class)
        cls.frame = cls.app.test_frame
        cls.panel = cls.app.panel_finder(cls.app.test_frame)
        cls.sizer = cls.panel.GetSizer()

        # NOTE!: Call Layout on the panel here, because otherwise the
        # controls layed out using XRC will not have the right sizes!
        gui_loop()

    @classmethod
    def tearDownClass(cls):
        if not MANUAL:
            cls.app.test_frame.Destroy()
            wx.CallAfter(cls.app.Exit)
        elif INSPECT:
            from wx.lib import inspection
            inspection.InspectionTool().Show()
        cls.app.MainLoop()

    def setUp(self):
        self.app.test_frame.SetTitle(
            "%s > %s" % (self.app.module_name, self._testMethodName))

    @classmethod
    def add_control(cls, ctrl, flags=0, border=10, proportion=0, clear=False):
        if clear:
            cls.remove_all()

        cls.sizer.Add(ctrl, flag=flags, border=border, proportion=proportion)
        cls.sizer.Layout()
        return ctrl

    @classmethod
    def remove_all(cls):
        for child in cls.sizer.GetChildren():
            cls.sizer.Remove(child.Window)
            child.Window.Destroy()
        cls.sizer.Layout()

    def create_simple_tab_model(self):
        main = gmodel.MainGUIData(None) # no microscope backend
        tab = gmodel.MicroscopyGUIData(main)

        # Add one view
        fview = gmodel.MicroscopeView("fakeview")
        tab.views.value.append(fview)
        tab.focussedView.value = fview

        return tab


# # Dummy clases for testing purposes
# class Object(object):
#     pass
#
# class MainGUIData(object):
#     """
#     Imitates a MainGUIData wrt stream entry: it just needs a focussedView
#     """
#     def __init__(self):
#         fview = gmodel.MicroscopeView("fakeview")
#         self.focussedView = omodel.VigilantAttribute(fview)
#
#         self.main = Object() #pylint: disable=E0602
#         self.main.light = None
#         self.main.ebeam = None
#         self.main.debug = omodel.VigilantAttribute(fview)
#         self.focussedView = omodel.VigilantAttribute(fview)
#
#         self.light = None
#         self.light_filter = None
#         self.ccd = None
#         self.sed = None
#         self.ebeam = None
#         self.tool = None
#         self.subscribe = None

# Utility functions

def set_img_meta(img, pixel_size, pos):
    img.metadata[omodel.MD_PIXEL_SIZE] = pixel_size
    img.metadata[omodel.MD_POS] = pos

def generate_img_data(width, height, depth, alpha=255):
    """ Create an image of the given dimensions """

    shape = (height, width, depth)
    rgb = numpy.empty(shape, dtype=numpy.uint8)

    if width > 100 or height > 100:
        tl = random_color(alpha=alpha)
        tr = random_color(alpha=alpha)
        bl = random_color(alpha=alpha)
        br = random_color(alpha=alpha)

        rgb = numpy.zeros(shape, dtype=numpy.uint8)

        rgb[..., -1, 0] = numpy.linspace(tr[0], br[0], height)
        rgb[..., -1, 1] = numpy.linspace(tr[1], br[1], height)
        rgb[..., -1, 2] = numpy.linspace(tr[2], br[2], height)

        rgb[..., 0, 0] = numpy.linspace(tl[0], bl[0], height)
        rgb[..., 0, 1] = numpy.linspace(tl[1], bl[1], height)
        rgb[..., 0, 2] = numpy.linspace(tl[2], bl[2], height)

        for i in xrange(height):
            sr, sg, sb = rgb[i, 0, :3]
            er, eg, eb = rgb[i, -1, :3]

            rgb[i, :, 0] = numpy.linspace(int(sr), int(er), width)
            rgb[i, :, 1] = numpy.linspace(int(sg), int(eg), width)
            rgb[i, :, 2] = numpy.linspace(int(sb), int(eb), width)

        if depth == 4:
            rgb[..., 3] = min(255, max(alpha, 0))

    else:
        for w in xrange(width):
            for h in xrange(height):
                rgb[h, w] = random_color((230, 230, 255), alpha)

    return omodel.DataArray(rgb)


def random_color(mix_color=None, alpha=255):
    """ Generate a random color, possibly tinted using mix_color """
    red = random.randint(0, 255)
    green = random.randint(0, 255)
    blue = random.randint(0, 255)

    if mix_color:
        red = (red - mix_color[0]) / 2
        green = (green - mix_color[1]) / 2
        blue = (blue - mix_color[2]) / 2

    a = alpha / 255.0

    return red * a, green * a, blue * a, alpha

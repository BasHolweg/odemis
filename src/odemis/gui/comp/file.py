# -*- coding: utf-8 -*-

"""

@author: Rinze de Laat

Copyright © 2014Rinze de Laat, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License version 2 as published by the Free
Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

Content:

    This module contains controls for file selection.

"""

import os
import logging

import wx
import wx.lib.newevent

import odemis.gui
from .buttons import ImageTextButton, ImageButton
from odemis.gui.img import data

FileSelectEvent, EVT_FILE_SELECT = wx.lib.newevent.NewEvent()

class FileBrowser(wx.Panel):
    """
    Widget that displays a file name and allows to change it by selecting a 
    different file.
    It will generate a EVT_FILE_SELECT when the file changes.
    Note that like most of the wx widgets, SetValue does not generate an event.
    """

    def __init__(self, parent, id=wx.ID_ANY,
                  pos=wx.DefaultPosition,
                  size=wx.DefaultSize,
                  style=wx.TAB_TRAVERSAL,
                  tool_tip=None,
                  clear_btn=False,
                  label="",
                  dialog_title="Browse for file",
                  wildcard="*.*",
                  name='fileBrowser',
                  default_dir=None
        ):

        self.file_path = None
        self.default_dir = default_dir or os.path.curdir

        self.dialog_title = dialog_title
        self.wildcard = wildcard
        self.clear_btn = clear_btn # Add clear buttons
        self.label = label # Text to show when the control is cleared

        self.text_ctrl = None
        self.btn_ctrl = None
        self._btn_clear = None

        self.create_dialog(parent, id, pos, size, style, name)

    def create_dialog(self, parent, id, pos, size, style, name):
        """Setup the graphic representation of the dialog"""
        wx.Panel.__init__(self, parent, id, pos, size, style, name)
        self.SetBackgroundColour(parent.GetBackgroundColour())

        box = wx.BoxSizer(wx.HORIZONTAL)

        self.text_ctrl = wx.TextCtrl(self,
                            style=wx.BORDER_NONE|wx.TE_READONLY)
        self.text_ctrl.SetForegroundColour(odemis.gui.FOREGROUND_COLOUR_EDIT)
        self.text_ctrl.SetBackgroundColour(odemis.gui.BACKGROUND_COLOUR)
        self.text_ctrl.Bind(wx.EVT_TEXT, self.on_changed)

        box.Add(self.text_ctrl, 1)

        if self.clear_btn:
            self._btn_clear = ImageButton(self,
                                          wx.ID_ANY,
                                          data.getico_clearBitmap(),
                                          (10, 8),
                                          (18, 18),
                                          background_parent=parent)
            self._btn_clear.SetBitmaps(data.getico_clear_hBitmap())
            self._btn_clear.SetToolTipString("Clear calibration") # FIXME: do not hard code
            self._btn_clear.Hide()
            self._btn_clear.Bind(wx.EVT_BUTTON, self._on_clear)
            box.Add(self._btn_clear, 0, wx.LEFT, 10)

        self.btn_ctrl = ImageTextButton(self, -1, data.getbtn_64x16Bitmap(),
                                        label_delta=1,
                                        style=wx.ALIGN_CENTER)
        self.btn_ctrl.SetBitmaps(data.getbtn_64x16_hBitmap(),
                                 data.getbtn_64x16_aBitmap())
        self.btn_ctrl.SetForegroundColour("#000000")
        self.btn_ctrl.SetLabel("change...")
        self.btn_ctrl.Bind(wx.EVT_BUTTON, self._on_browse)

        box.Add(self.btn_ctrl, 0, wx.LEFT, 5)

        self.SetAutoLayout(True)
        self.SetSizer(box)
        self.Layout()
        if isinstance(size, tuple):
            size = wx.Size(size)
        self.SetDimensions(-1, -1, size.width, size.height, wx.SIZE_USE_EXISTING)

    def on_changed(self, evt):
        evt.SetEventObject(self)
        evt.Skip()

    def _SetValue(self, file_path, raise_event):

        if file_path:
            logging.debug("Setting file control to %s", file_path)

            self.file_path = file_path

            if not os.path.exists(self.file_path):
                self.text_ctrl.SetForegroundColour(odemis.gui.ALERT_COLOUR)
            else:
                self.text_ctrl.SetForegroundColour(
                                            odemis.gui.FOREGROUND_COLOUR_EDIT)

            self.text_ctrl.SetValue(self.file_path)

            self.text_ctrl.SetToolTipString(self.file_path)
            self.text_ctrl.SetInsertionPointEnd()
            self._btn_clear.Show()
        else:
            logging.debug("Clearing file control")

            self.file_path = None
            self.text_ctrl.SetForegroundColour(odemis.gui.FOREGROUND_COLOUR_DIS)

            self.text_ctrl.SetValue(self.label)

            self.text_ctrl.SetToolTipString("")
            self._btn_clear.Hide()

        self.Layout()

        if raise_event:
            wx.PostEvent(self, FileSelectEvent(selected_file=self.file_path))

    def SetValue(self, file_path):
        logging.debug("File set to '%s' by Odemis", file_path)
        self._SetValue(file_path, raise_event=False)

    def GetValue(self):
        return self.file_path

    @property
    def basename(self):
        """
        the base name of the file
        """
        return os.path.basename(self.file_path or "")

    @property
    def path(self):
        """
        the name of the directory containing the file
        """
        return os.path.dirname(self.file_path or "")

    def SetWildcard(self, wildcard):
        self.wildcard = wildcard

    def _on_clear(self, evt):
        self._SetValue(None, raise_event=True)

    def clear(self):
        self.SetValue(None)

    def _on_browse(self, evt):
        current = self.GetValue() or ""
        directory = os.path.split(current)

        if os.path.isdir(current):
            directory = current
            current = ""
        elif directory and os.path.isdir(directory[0]):
            current = directory[1]
            directory = directory[0]
        else:
            directory = self.default_dir
            current = ""

        dlg = wx.FileDialog(self, self.dialog_title, directory, current,
                            wildcard=self.wildcard,
                            style=wx.FD_OPEN)


        if dlg.ShowModal() == wx.ID_OK:
            self._SetValue(dlg.GetPath(), raise_event=True)
        dlg.Destroy()

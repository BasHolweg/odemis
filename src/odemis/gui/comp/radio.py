#-*- coding: utf-8 -*-

"""
@author: Rinze de Laat

Copyright © 2012 Rinze de Laat, Delmic

Custom (graphical) radio button control.

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 2 of the License, or (at your option) any later
version.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

"""

from odemis.gui.comp.buttons import GraphicRadioButton
import logging
import odemis.gui.img.data as img
import wx


class GraphicalRadioButtonControl(wx.Panel):

    def __init__(self, *args, **kwargs)    :

        #self.bnt_width = kwargs.pop("bnt_width", 32)

        self.choices = kwargs.pop("choices", [])
        self.buttons = []
        self.labels = kwargs.pop("labels", [])
        self.units = kwargs.pop("units", None)



        wx.Panel.__init__(self, *args, **kwargs)

        self.SetBackgroundColour(self.Parent.GetBackgroundColour())

        sizer = wx.BoxSizer(wx.HORIZONTAL)

        for choice, label in zip(self.choices, self.labels):
            btn = GraphicRadioButton(self,
                                     -1,
                                     img.getbtn_32x16Bitmap(),
                                     value=choice,
                                     #size=(self.bnt_width, 16),
                                     style=wx.ALIGN_CENTER,
                                     label=label,
                                     label_delta=1)

            btn.SetForegroundColour("#111111")

            btn.SetBitmaps(img.getbtn_32x16_hBitmap(),
                           img.getbtn_32x16_aBitmap(),
                           img.getbtn_32x16_aBitmap())

            self.buttons.append(btn)

            sizer.Add(btn, flag=wx.RIGHT, border=5)
            btn.Bind(wx.EVT_BUTTON, self.OnClick)
            btn.Bind(wx.EVT_KEY_UP, self.OnKeyUp)

        if self.units:
            lbl = wx.StaticText(self, -1, self.units)
            lbl.SetForegroundColour("#DDDDDD")
            sizer.Add(lbl, flag=wx.RIGHT, border=5)

        self.SetSizer(sizer)

    def _reset_buttons(self, btn=None):
        for button in [b for b in self.buttons if b != btn]:
            button.SetToggle(False)

    def SetValue(self, value):
        logging.debug("Set radio button control to %s", value)
        self._reset_buttons()

        for i, btn in enumerate(self.buttons):
            if btn.value == value:
                self.buttons[i].SetToggle(True)

    def GetValue(self):
        for btn in self.buttons:
            if btn.GetToggle():
                return btn.value

    def OnKeyUp(self, evt):
        btn = evt.GetEventObject()
        if btn.hasFocus and evt.GetKeyCode() == ord(" "):
            self._reset_buttons(btn)
            btn.up = False
            btn.Notify()
            btn.Refresh()

    def OnClick(self, evt):
        btn = evt.GetEventObject()
        self._reset_buttons(btn)
        #if not btn.GetToggle():
        evt.Skip()



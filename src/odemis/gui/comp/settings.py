#-*- coding: utf-8 -*-
"""
.. codeauthor:: Rinze de Laat <laat@delmic.com>

Copyright © 2014 Rinze de Laat, Delmic

This file is part of Odemis.

.. license::
    Odemis is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License version 2 as published by the Free
    Software Foundation.

    Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with
    Odemis. If not, see http://www.gnu.org/licenses/.

"""

from decorator import decorator

import wx

import odemis.gui as gui
from odemis.gui.comp.combo import ComboBox
from odemis.gui.comp.foldpanelbar import FoldPanelItem
from odemis.gui.comp.radio import GraphicalRadioButtonControl
from odemis.gui.comp.slider import UnitIntegerSlider, UnitFloatSlider
from odemis.gui.comp.text import UnitIntegerCtrl, UnitFloatCtrl


@decorator
def control_bookkeeper(f, self, *args, **kwargs):
    """ Clear the default message, if needed, and advance the row count """
    self.clear_default_message()
    result = f(self, *args, **kwargs)
    # Redo FoldPanelBar layout
    self.Parent.Parent.Layout()
    self.num_rows += 1
    return result


class SettingsPanel(wx.Panel):

    def __init__(self, *args, **kwargs):
        default_msg = kwargs.pop('default_msg', "")
        super(SettingsPanel, self).__init__(*args, **kwargs)

        assert isinstance(self.Parent, FoldPanelItem)

        self.SetBackgroundColour(gui.BG_COLOUR_MAIN)
        self.SetForegroundColour(gui.FG_COLOUR_MAIN)

        self.num_rows = 0  # The number of rows in the GridBagSizer

        # The main sizer is used to create a margin on the inside of the panel
        self._main_sizer = wx.BoxSizer(wx.VERTICAL)
        # The GridBagSizer is use to create a 2-column lay-out for the settings controls
        self._gb_sizer = wx.GridBagSizer()
        self._gb_sizer.SetEmptyCellSize((0, 0))

        # The default message text is added here, because at least control needs to be present
        # before the growable column can be added.
        self.message_text = None
        # A default message need to be inserted here, because otherwise AddGrowableCol will cause
        # an exception, since there wouldn't be any columns yet.
        self.set_default_message(default_msg)

        # Make the 2nd column expand
        self._gb_sizer.AddGrowableCol(1)

        self.SetSizer(self._main_sizer)
        self._main_sizer.Add(self._gb_sizer, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)

    def set_default_message(self, msg):
        """ Set the default message in the settings panel """

        if not self.message_text:
            self.message_text = wx.StaticText(self, -1, msg)
            self._gb_sizer.Add(self.message_text, (0, 0), (1, 2))
            self.num_rows = 1
        else:
            self.message_text.SetLabel(msg)
        self.message_text.Show()

    def clear_default_message(self):
        """ Remove the default message if it exists """
        if self.message_text:
            self.message_text.Hide()

    def clear_all(self):
        """ Remove all children """
        for c in self.GetChildren():
            c.Destroy()
        self.num_rows = 0
        self.Parent.Parent.Layout()

    # Control methods

    def _add_side_label(self, label_text):
        """ Add a static text label to the left column at the current row

        This method should only be called from another control adding method!

        """

        self.clear_default_message()

        # Create label
        lbl_ctrl = wx.StaticText(self, -1, unicode(label_text))
        self._gb_sizer.Add(lbl_ctrl, (self.num_rows, 0),
                           flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)

        return lbl_ctrl

    @control_bookkeeper
    def add_divider(self):
        """ Add a horizontal divider to the panel """
        line_ctrl = wx.StaticLine(self, size=(-1, 1))
        self._gb_sizer.Add(line_ctrl, (self.num_rows, 0), span=(1, 2),
                           flag=wx.ALL | wx.EXPAND, border=5)

    @control_bookkeeper
    def add_readonly_field(self, label_text, value=None, selectable=True):
        """ Adds a value to the control panel that cannot directly be changed by the user

        :param label_text: (str) Label text to display
        :param value: (None or object) Value to display next to the label
        :param selectable: (boolean) whether the value can be selected for copying by the user

        :return: (Ctrl, Ctrl or None) Label and value control

        """

        lbl_ctrl = self._add_side_label(label_text)

        if value:
            if selectable:
                value_ctrl = wx.TextCtrl(self, value=unicode(value),
                                         style=wx.BORDER_NONE | wx.TE_READONLY)
                value_ctrl.SetForegroundColour(gui.FG_COLOUR_DIS)
                value_ctrl.SetBackgroundColour(gui.BG_COLOUR_MAIN)
                self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                                   flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)
            else:
                value_ctrl = wx.StaticText(self, label=unicode(value))
                value_ctrl.SetForegroundColour(gui.FG_COLOUR_DIS)
                self._gb_sizer.Add(value_ctrl, (self.num_rows, 1), flag=wx.ALL, border=5)
        else:
            value_ctrl = None

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_text_field(self, label_text, value=None, readonly=False):
        """ Add a label and text control to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or str) Value to display
        :param readonly: (boolean) Whether the value can be changed by the user

        :return: (Ctrl, Ctrl) Label and text control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = wx.TextCtrl(self, value=unicode(value or ""),
                                 style=wx.BORDER_NONE | (wx.TE_READONLY if readonly else 0))
        if readonly:
            value_ctrl.SetForegroundColour(gui.FG_COLOUR_DIS)
        else:
            value_ctrl.SetForegroundColour(gui.FG_COLOUR_EDIT)
        value_ctrl.SetBackgroundColour(gui.BG_COLOUR_MAIN)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_integer_slider(self, label_text, value=None, conf=None):
        """ Add an integer value slider to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or int) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = UnitIntegerSlider(self, value=value, **conf)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_float_slider(self, label_text, value=None, conf=None):
        """ Add a float value slider to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or float) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = UnitFloatSlider(self, value=value, **conf)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_float_slider(self, label_text, value=None, conf=None):
        """ Add a float value slider to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or float) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = UnitFloatSlider(self, value=value, **conf)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_int_field(self, label_text, value=None, conf=None):
        """ Add an integer value field to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or int) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        return self._add_num_field(UnitIntegerCtrl, label_text, value, conf)

    @control_bookkeeper
    def add_float_field(self, label_text, value=None, conf=None):
        """ Add a float value field to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or float) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        return self._add_num_field(UnitFloatCtrl, label_text, value, conf)

    def _add_num_field(self, klass, label_text, value, conf):

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = klass(self, value=value, style=wx.NO_BORDER, **conf)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)
        value_ctrl.SetForegroundColour(gui.FG_COLOUR_EDIT)
        value_ctrl.SetBackgroundColour(gui.BG_COLOUR_MAIN)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_radio_control(self, label_text, value=None, conf=None):
        """ Add a series of radio buttons to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or float) Value to display
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = GraphicalRadioButtonControl(self, -1, style=wx.NO_BORDER, **conf)
        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)
        if value:
            value_ctrl.SetValue(value)

        return lbl_ctrl, value_ctrl

    @control_bookkeeper
    def add_combobox_control(self, label_text, value=None, conf=None):
        """ Add a combo box to the settings panel

        :param label_text: (str) Label text to display
        :param value: (None or float) Value to display *NOT USED ATM*
        :param conf: (None or dict) Dictionary containing parameters for the control

        """

        lbl_ctrl = self._add_side_label(label_text)
        value_ctrl = ComboBox(self, -1, value='', pos=(0, 0), size=(-1, 12),
                              style=wx.NO_BORDER | wx.TE_PROCESS_ENTER, **conf)

        self._gb_sizer.Add(value_ctrl, (self.num_rows, 1),
                           flag=wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, border=5)

        # TODO: move this to ComboBox?
        def _eat_event(evt):
            """ Quick and dirty empty function used to 'eat'
            mouse wheel events
            """
            # TODO: This solution only makes sure that the control's value
            # doesn't accidentally get altered when it gets hit by a mouse
            # wheel event. However, it also stop the event from propagating
            # so the containing scrolled window will not scroll either.
            # (If the event is skipped, the control will change value again)
            pass

        value_ctrl.Bind(wx.EVT_MOUSEWHEEL, _eat_event)

        if value:
            value_ctrl.SetValue(unicode(value))

        return lbl_ctrl, value_ctrl


# END Control methods

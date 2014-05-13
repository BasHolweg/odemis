# -*- coding: utf-8 -*-
'''
Created on 6 Sep 2013

@author: Éric Piel

Copyright © 2013 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.
'''
from __future__ import division
from odemis.gui.util.widgets import VigilantAttributeConnector
from odemis.gui.comp.combo import ComboBox
import logging
import wx



# Known good key bindings
# WXK -> (args for tab_data_model.step())
KB_SECOM = {
        wx.WXK_LEFT: ("stage", "x", -1),
        wx.WXK_RIGHT: ("stage", "x", 1),
        wx.WXK_DOWN: ("stage", "y", -1),
        wx.WXK_UP: ("stage", "y", 1),
        wx.WXK_PAGEDOWN: ("focus", "z", -1),
        wx.WXK_PAGEUP: ("focus", "z", 1),
        wx.WXK_NUMPAD_HOME: ("aligner", "b", -1),
        wx.WXK_NUMPAD_PAGEDOWN: ("aligner", "b", 1),
        wx.WXK_NUMPAD_END: ("aligner", "a", -1),
        wx.WXK_NUMPAD_PAGEUP: ("aligner", "a", 1),
        # same but with NumLock
        wx.WXK_NUMPAD7: ("aligner", "b", -1),
        wx.WXK_NUMPAD3: ("aligner", "b", 1),
        wx.WXK_NUMPAD1: ("aligner", "a", -1),
        wx.WXK_NUMPAD9: ("aligner", "a", 1),
}

KB_SPARC = {
        wx.WXK_LEFT: ("mirror", "x", -1), # so that image goes in same direction
        wx.WXK_RIGHT: ("mirror", "x", 1),
        wx.WXK_DOWN: ("mirror", "y", 1),
        wx.WXK_UP: ("mirror", "y", -1),
        wx.WXK_PAGEDOWN: ("focus_sem", "z", -1),
        wx.WXK_PAGEUP: ("focus_sem", "z", 1),
        wx.WXK_NUMPAD_LEFT: ("mirror", "rz", -1),
        wx.WXK_NUMPAD_RIGHT: ("mirror", "rz", 1),
        wx.WXK_NUMPAD_DOWN: ("mirror", "ry", 1),
        wx.WXK_NUMPAD_UP: ("mirror", "ry", -1),
        # same but with NumLock
        wx.WXK_NUMPAD4: ("mirror", "rz", -1),
        wx.WXK_NUMPAD6: ("mirror", "rz", 1),
        wx.WXK_NUMPAD2: ("mirror", "ry", 1),
        wx.WXK_NUMPAD8: ("mirror", "ry", -1),
}

KEY_BINDINGS = {
            "sparc": KB_SPARC,
            "secom": KB_SECOM,
}

class ActuatorController(object):
    """ This controller manages the buttons to manually move the actuators.
    """

    def __init__(self, tab_data, main_frame, tab_prefix):
        """ Binds the step and axis buttons to their appropriate
        Vigilant Attributes in the model.ActuatorGUIData. It only connects the
        buttons which exists with the actuators which exists.

        tab_data (ActuatorGUIData): the data model of the tab
        main_frame: (wx.Frame): the main frame of the GUI
        tab_prefix (string): common prefix of the names of the buttons
        """
        self._tab_data_model = tab_data
        self._main_frame = main_frame
        # Check which widgets and VAs exist. Bind the matching ones.

        # Bind size steps (= sliders)
        self._va_connectors = []
        for an, ss in tab_data.stepsizes.items():
            slider_name = tab_prefix + "slider_" + an
            try:
                slider = getattr(main_frame, slider_name)
            except AttributeError:
                continue

            slider.SetRange(*ss.range)

            vac = VigilantAttributeConnector(ss, slider, events=wx.EVT_SLIDER)
            self._va_connectors.append(vac)
        if not self._va_connectors:
            logging.warning("No slider found for tab %s", tab_prefix)

        # Bind buttons
        for actuator, axis in tab_data.axes:
            for suffix, factor in [("m", -1), ("p", 1)]:
                # something like "lens_align_btn_p_mirror_rz"
                btn_name = "%sbtn_%s_%s_%s" % (tab_prefix, suffix, actuator, axis)
                try:
                    btn = getattr(main_frame, btn_name)
                except AttributeError:
                    logging.debug("No button in GUI found for axis %s", axis)
                    continue

                def btn_action(evt, tab_data=tab_data, actuator=actuator, axis=axis, factor=factor):
                    # Button events don't contain key state, so check ourselves
                    if wx.GetKeyState(wx.WXK_SHIFT):
                        factor /= 10
                    tab_data.step(actuator, axis, factor)

                btn.Bind(wx.EVT_BUTTON, btn_action)


    def bind_keyboard(self, tab_frame):
        """
        Bind the keyboard keys to the actuator axes
        """
        role = self._tab_data_model.main.role
        try:
            self.key_bindings = KEY_BINDINGS[role]
        except KeyError:
            logging.warning("No known key binding for microscope %s", role)
            return
        # Remove keys for axes not available
        for key, (actuator, axis, _) in self.key_bindings.items():
            if not (actuator, axis) in self._tab_data_model.axes:
                del self.key_bindings[key]

        # Keybinding is difficult:
        # evt_key_* and evt_char are not passed to their parents, even if
        # skipped. Only evt_char_hook is propagated, the problem is that it's
        # not what the children bind to, so we always get it, even if the child
        # handles the key events.
        # http://article.gmane.org/gmane.comp.python.wxpython/50485
        # http://wxpython.org/Phoenix/docs/html/KeyEvent.html
        tab_frame.Bind(wx.EVT_CHAR_HOOK, self._on_key)


    def _on_key(self, event):
        key = event.GetKeyCode()
        if key in self.key_bindings:
            # check the focus is not on some children that'll handle the key
            focusedWin = wx.Window.FindFocus()
            # TODO: need to check for more widget types?
            if not isinstance(focusedWin, (wx.TextCtrl, ComboBox)):
                actuator, axis, size = self.key_bindings[key]
                if event.ShiftDown():
                    size /= 10
                self._tab_data_model.step(actuator, axis, size)
                return # keep it for ourselves

        # everything else we don't process
        event.Skip()

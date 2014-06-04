#-*- coding: utf-8 -*-
"""
@author: Rinze de Laat

Copyright © 2012 Rinze de Laat, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License version 2 as published by the Free Software
Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/. """

import logging
import os.path
from logging.handlers import RotatingFileHandler

import wx

LOG_FILE = "odemis-gui.log"

LOG_LINES = 500 # maximum lines in the GUI logger
log = logging.getLogger() # for compatibility only


def init_logger(level=logging.DEBUG):
    """
    Initializes the logger to some nice defaults
    To be called only once, at the initialisation
    """
    logging.basicConfig(format=" - %(levelname)s \t%(message)s")
    l = logging.getLogger()
    l.setLevel(level)
    frm = "%(asctime)s  %(levelname)-7s %(module)-15s: %(message)s"
    l.handlers[0].setFormatter(logging.Formatter(frm))

def create_gui_logger(log_field, error_va=None):
    """
    log_field (wx text field)
    error_va (Boolean VigilantAttribute)
    """
    # Create file handler

    # Path to the log file
    logfile_path = os.path.join(os.path.expanduser("~"), LOG_FILE)
    # Formatting string for logging messages to file
    frm = "%(asctime)s %(levelname)-7s %(module)s:%(lineno)d: %(message)s"
    file_format = logging.Formatter(frm)

    # Max 5 log files of 10Mb
    file_handler = RotatingFileHandler(logfile_path,
                                       maxBytes=10 * (2 ** 20),
                                       backupCount=5)

    file_handler.setFormatter(file_format)

    # Create gui handler
    frm = "%(asctime)s %(levelname)-7s %(module)-15s: %(message)s"
    gui_format = logging.Formatter(frm, '%H:%M:%S')
    text_field_handler = TextFieldHandler()
    text_field_handler.setTextField(log_field)
    if error_va is not None:
        text_field_handler.setErrorVA(error_va)
    text_field_handler.setFormatter(gui_format)
    logging.debug("Switching to GUI logger")

    # remove whatever handler was already there
    for handler in log.handlers:
        log.removeHandler(handler)

    try:
        log.addHandler(text_field_handler)
        log.addHandler(file_handler)
    except:
        # Use print here because log probably doesn't work
        print("Failed to set-up logging handlers")
        logging.exception("Failed to set-up logging handlers")
        raise

def stop_gui_logger():
    """
    Stop the logger from displaying logs to the GUI.
    Use just before ending the GUI.
    """

    # remove whatever handler was already there
    for handler in log.handlers:
        if isinstance(handler, TextFieldHandler):
            log.removeHandler(handler)

class TextFieldHandler(logging.Handler):
    """ Custom log handler, used to output log entries to a text field. """
    def __init__(self):
        """ Call the parent constructor and initialize the handler """
        logging.Handler.__init__(self)
        self.textfield = None
        self.error_va = None

    def setTextField(self, textfield):
        self.textfield = textfield
        self.textfield.Clear()

    def setErrorVA(self, error_va):
        self.error_va = error_va

    def emit(self, record):
        """ Write a record, in colour, to a text field. """
        if self.textfield is not None:
            if record.levelno >= logging.ERROR:
                colour = "#B00B2C"
            elif record.levelno == logging.WARNING:
                colour = "#C87000"
            elif record.levelno == logging.INFO:
                colour = "#555555"
            else:
                colour = "#777777"

            # FIXME: still seems to be possible to completely hog the GUI by
            # logging too much. A way to fix it would be to run the textfield
            # update at a maximum frequency (10Hz), and queue the logs in
            # between.

            # Do the actual writing in a CallAfter, so logging won't interfere
            # with the GUI drawing process.
            wx.CallAfter(self.write_to_field, record, colour)

        # Will typically ensure that the text field is displayed
        if self.error_va is not None and record.levelno >= logging.ERROR:
            self.error_va.value = True

    def write_to_field(self, record, colour):
        nb_lines = self.textfield.GetNumberOfLines()
        nb_old = nb_lines - LOG_LINES
        if nb_old > 0:
            # Removes the characters from position 0 up to and including the
            # Nth line break
            first_new = 0
            txt = self.textfield.Value
            for i in range(nb_old):
                first_new = txt.find('\n', first_new) + 1

            self.textfield.Remove(0, first_new)

        self.textfield.SetDefaultStyle(wx.TextAttr(colour, None))
        self.textfield.AppendText("\n" + self.format(record))


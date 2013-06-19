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


def init_logger():
    """
    Initializes the logger to some nice defaults
    To be called only once, at the initialisation
    """
    logging.basicConfig(format=" - %(levelname)s \t%(message)s")
    l = logging.getLogger()
    l.setLevel(logging.DEBUG)
    frm = "%(asctime)s  %(levelname)-7s %(module)-15s: %(message)s"
    l.handlers[0].setFormatter(logging.Formatter(frm))

def create_gui_logger(log_field):
    # Create file handler

    # Path to the log file
    logfile_path = os.path.join(os.path.expanduser("~"), LOG_FILE)
    # Maximum size of the log file before it's rotated
    max_logfile_size = 1024**2
    # Maximum number of (rotated) log files
    max_logfile_count = 5
    # Formatting string for logging messages to file
    frm = "%(asctime)s %(levelname)-7s %(module)s:%(lineno)d: %(message)s"
    file_format = logging.Formatter(frm)

    file_handler = RotatingFileHandler(logfile_path, 'w',
                                       max_logfile_size,
                                       max_logfile_count)

    file_handler.setFormatter(file_format)

    # Create gui handler
    frm = "%(asctime)s %(levelname)-7s %(module)-15s: %(message)s"
    gui_format = logging.Formatter(frm, '%H:%M:%S')
    text_field_handler = TextFieldHandler()
    text_field_handler.setTextField(log_field)
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

    def setTextField(self, textfield):
        self.textfield = textfield
        self.textfield.Clear()

    def emit(self, record):
        """ Write a record, in color, to a text field. """
        if self.textfield is not None:
            color = "#777777"

            if record.levelno > logging.WARNING:
                color = "#B00B2C"
            elif record.levelno > logging.INFO:
                color = "#C87000"
            elif record.levelno > logging.DEBUG:
                color = "#555555"
            else:
                color = "#777777"

            # Do the actual writing in a CallAfter, so logging won't interfere
            # with the GUI drawing process.
            wx.CallAfter(self.write_to_field, record, color)

    def write_to_field(self, record, color):

        while self.textfield.GetNumberOfLines() > LOG_LINES:
            # Removes the characters from position 0 up to and including the
            # first line break
            self.textfield.Remove(0, self.textfield.GetValue().find('\n') + 1)

        self.textfield.SetDefaultStyle(wx.TextAttr(color, None))
        self.textfield.AppendText(''.join([self.format(record), '\n']))
        self.textfield.LineUp()



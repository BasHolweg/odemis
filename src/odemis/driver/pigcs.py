# -*- coding: utf-8 -*-
'''
Created on 7 Aug 2012

@author: Éric Piel

Copyright © 2012 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms 
of the GNU General Public License version 2 as published by the Free Software 
Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
Odemis. If not, see http://www.gnu.org/licenses/.
'''
from __future__ import division

from _curses import baudrate
import collections
from concurrent import futures
import glob
import logging
from odemis import model
import odemis
from odemis.model import isasync
from odemis.util import driver
import os
import re
import serial
import socket
import threading
import time


# Driver to handle PI's piezo motor controllers that follow the 'GCS' (General
# Command Set). In particular it handle the PI E-861 controller. Information can
# be found the manual E-861_User_PZ205E121.pdf (p.107). See PIRedStone for the PI C-170.
#
# In a daisy-chain, connected via USB or via RS-232, there must be one
# controller with address 1 (=DIP 1111). There is also a broadcast address: 255.
#
# The controller contains many parameters in flash memory. These parameters must
# have been previously written correctly to fit the stage it is driving. This
# driver uses and expects the parameters to be correct.
# The configuration database is available in a file called pistages2.dat. The
# PIMikroMove Windows program allows to load it, but by default doesn't copy it to
# the non-volatile memory, so you need to explicitly save them in persistent
# memory. (Can also be done with the WPA command and password "100".)
# In particular, for closed-loop stages, ensure that the settling windows and time
# are correctly set so that a move is considered on target quickly (< 1 s).
#
# The controller support closed-loop mode (i.e., absolute positioning) but only
# if it is associated to a sensor (not software detectable). It can also work in
# open-loop mode but to avoid damaging the hardware (which is moved by this
# actuator):
# * Do not switch servo on (SVO command)
# * Do not send commands for closed-loop motion, like MOV or MVR
# * Do not send the open-loop commands OMA and OMR, since they
#    use a sensor, too
#
# The controller accepts several baud rates. We choose 38400 (DIP=01) as it's fast
# and it seems accepted by every version. Other settings are 8 data, 1 stop,
# no parity.
#
#
# In open-loop, the controller has 2 ways to move the actuators:
#  * Nanostepping: high-speed, and long distance
#       1 step ~ 10 μm without load (less with load)
#  * Analog: very precise, but moves maximum ~5μm
#      "40 volts corresponds to a motion of approx. 3.3μm"
#      "20 volts corresponds to a motion of approx. 1μm"
#
# As an exception, the C-867 only supports officially closed-loop. However, there
# is a "testing" command, SMO, that allows to move in open-loop by simulating the input
# to the PID controller. PI assured us that as long as the stage used can reach
# the limit without getting damaged, it is safe. It's pretty straightforward to
# use the command. The voltage defines the speed (and direction) of the move. The
# voltage should be set to 0 again when the position desired is reached. 3V is
# approximately the minimum to move, and 10V is the maximum. Voltage is more or
# less linear between -32766 and 32766 -> -10 and 10V. So the distance moved
# depends on the time the SMO is set, which is obviously very imprecise.
#
# In closed-loop, it's almost all automagical.
# There are two modes in closed-loop: before and after referencing. Referencing
# consists in going to at least one "reference" point so that the actual position
# is known.
#  * Non referenced: that's the only one possible just after boot. It's only
#    possible to do relative moves.
#  * Referenced: both absolute and relative moves are possible. It's the default
#    mode.
# The problem with referencing is that for some cases, it might be dangerous to
# move the actuator, so a user feedback is needed. This means an explicit request
# via the API must be done before this is going on, and stopping must be possible.
# In addition, in many cases, relative move is sufficient.
# Note that for the closed-loop to work, in addition to the actuator there must
# be:
#   * a sensor (which indicates a distance)
#   * a reference switch (which indicates a point, usually in the middle)
#   * 2 limit switches (which indicate the borders)
# Most of the time, all of these are present, but the controller can do with just
# either the ref switch or the limit switches.
#
# The recommended maximum step frequency is 800 Hz.
#
# The architecture of the driver relies on four main classes:
#  * Controller: represent one controller with one or several axes (E-861 has only
#     one). Each subclass defines a different type of control. The subclass is
#     picked dynamically according to the physical controller.
#  * Bus: represent the whole group of controllers daisy-chained from the same
#     serial port. It's also the Actuator interface for the rest of Odemis.
#  * SerialBusAccesser: wrapper around the serial port which handles the low-level GCS
#     protocol.
#  * ActionManager: handles all the actions (move/stop) sent to the controller so
#     that the asynchronous ones are ordered.
#
# In the typical usage, Odemis ask to moveRel() an axis to the Bus. The Bus converts
# it into an action, returns a Future and queue the action on the ActionManager.
# When the Controller is free, the ActionManager pick the next action and convert
# it into a command for the Controller, which sends it to the actual PI controller
# and waits for it to finish.
#
# Note: in some rare cases, the controller might not answer to commands correctly,
# reporting error 555. In that case, it's possible to do a factory reset with the
# hidden command (which must be followed by the reconfiguration of the parameters):
# zzz 100 parameter
class PIGCSError(Exception):

    def __init__(self, errno):
        self.errno = errno
        self.strerror = self._errordict.get(errno, "Unknown error")

    def __str__(self):
        return "PIGCS error %d: %s" % (self.errno, self.strerror)

    _errordict = {
        0: "No error",
        1: "Parameter syntax error",
        2: "Unknown command",
        3: "Command length out of limits or command buffer overrun",
        4: "Error while scanning",
        5: "Unallowable move attempted on unreferenced axis, or move attempted with servo off",
        6: "Parameter for SGA not valid",
        7: "Position out of limits",
        8: "Velocity out of limits",
        9: "Attempt to set pivot point while U,V and W not all 0",
        10: "Controller was stopped by command",
        11: "Parameter for SST or for one of the embedded scan algorithms out of range",
        12: "Invalid axis combination for fast scan",
        13: "Parameter for NAV out of range",
        14: "Invalid analog channel",
        15: "Invalid axis identifier",
        16: "Unknown stage name",
        17: "Parameter out of range",
        18: "Invalid macro name",
        19: "Error while recording macro",
        20: "Macro not found",
        21: "Axis has no brake",
        22: "Axis identifier specified more than once",
        23: "Illegal axis",
        24: "Incorrect number of parameters",
        25: "Invalid floating point number",
        26: "Parameter missing",
        27: "Soft limit out of range",
        28: "No manual pad found",
        29: "No more step-response values",
        30: "No step-response values recorded",
        31: "Axis has no reference sensor",
        32: "Axis has no limit switch",
        33: "No relay card installed",
        34: "Command not allowed for selected stage(s)",
        35: "No digital input installed",
        36: "No digital output configured",
        37: "No more MCM responses",
        38: "No MCM values recorded",
        39: "Controller number invalid",
        40: "No joystick configured",
        41: "Invalid axis for electronic gearing, axis cannot be slave",
        42: "Position of slave axis is out of range",
        43: "Slave axis cannot be commanded directly when electronic gearing is enabled",
        44: "Calibration of joystick failed",
        45: "Referencing failed",
        46: "OPM (Optical Power Meter) missing",
        47: "OPM (Optical Power Meter) not initialized or cannot be initialized",
        48: "OPM (Optical Power Meter) Communication Error",
        49: "Move to limit switch failed",
        50: "Attempt to reference axis with referencing disabled",
        51: "Selected axis is controlled by joystick",
        52: "Controller detected communication error",
        53: "MOV! motion still in progress",
        54: "Unknown parameter",
        55: "No commands were recorded with REP",
        56: "Password invalid",
        57: "Data Record Table does not exist",
        58: "Source does not exist; number too low or too high",
        59: "Source Record Table number too low or too high",
        60: "Protected Param: current Command Level (CCL) too low",
        61: "Command execution not possible while Autozero is running",
        62: "Autozero requires at least one linear axis",
        63: "Initialization still in progress",
        64: "Parameter is read-only",
        65: "Parameter not found in non- volatile memory",
        66: "Voltage out of limits",
        67: "Not enough memory available for requested wave curve",
        68: "Not enough memory available for DDL table; DDL cannot be started",
        69: "Time delay larger than DDL table; DDL cannot be started",
        70: "The requested arrays have different lengths; query them separately",
        71: "Attempt to restart the generator while it is running in single step mode",
        72: "Motion commands and wave generator activation are not allowed when analog target is active",
        73: "Motion commands are not allowed when wave generator is active",
        74: "No sensor channel or no piezo channel connected to selected axis (sensor and piezo matrix)",
        75: "Generator started (WGO) without having selected a wave table (WSL).",
        76: "Interface buffer did overrun and command couldn't be received correctly",
        77: "Data Record Table does not hold enough recorded data",
        78: "Data Record Table is not configured for recording",
        79: "Open-loop commands (SVA, SVR) are not allowed when servo is on",
        80: "Hardware error affecting RAM",
        81: "Not macro command",
        82: "Macro counter out of range",
        83: "Joystick is active",
        84: "Motor is off",
        85: "Macro-only command",
        86: "Invalid joystick axis",
        87: "Joystick unknown",
        88: "Move without referenced stage",
        89: "Command not allowed in current motion mode",
        90: "No tracing possible while digital IOs are used on this HW revision. Reconnect to switch operation mode.",
        91: "Move not possible, would cause collision",
        92: "Stage is not capable of following the master. Check the gear ratio.",
        100: "PI LabVIEW driver reports error. See source control for details.",
        200: "No stage connected to axis",
        201: "File with axis parameters not found",
        202: "Invalid axis parameter file",
        203: "Backup file with axis parameters not found",
        204: "PI internal error code 204",
        205: "SMO with servo on",
        206: "uudecode: incomplete header",
        207: "uudecode: nothing to decode",
        208: "uudecode: illegal UUE format",
        209: "CRC32 error",
        210: "Illegal file name (must be 8-0 format)",
        211: "File not found on controller",
        212: "Error writing file on controller",
        213: "VEL command not allowed in DTR Command Mode",
        214: "Position calculations failed",
        215: "The connection between controller and stage may be broken",
        216: "The connected stage has driven into a limit switch, some controllers need CLR to resume operation",
        217: "Strut test command failed because of an unexpected strut stop",
        218: "While MOV! is running position can only be estimated!",
        219: "Position was calculated during MOV motion",
        230: "Invalid handle",
        231: "No bios found",
        232: "Save system configuration failed",
        233: "Load system configuration failed",
        301: "Send buffer overflow",
        302: "Voltage out of limits",
        303: "Open-loop motion attempted when servo ON",
        304: "Received command is too long",
        305: "Error while reading/writing EEPROM",
        306: "Error on I2C bus",
        307: "Timeout while receiving command",
        308: "A lengthy operation has not finished in the expected time",
        309: "Insufficient space to store macro",
        310: "Configuration data has old version number",
        311: "Invalid configuration data",
        333: "Internal hardware error",
        400: "Wave generator index error",
        401: "Wave table not defined",
        402: "Wave type not supported",
        403: "Wave length exceeds limit",
        404: "Wave parameter number error",
        405: "Wave parameter out of range",
        406: "WGO command bit not supported",
        502: "Position consistency check failed",
        503: "Hardware collision sensor(s) are activated",
        504: "Strut following error occurred, e.g. caused by overload or encoder failure",
        555: "BasMac: unknown controller error",
        601: "not enough memory",
        602: "hardware voltage error",
        603: "hardware temperature out of range",
        1000: "Too many nested macros",
        1001: "Macro already defined",
        1002: "Macro recording not activated",
        1003: "Invalid parameter for MAC",
        1004: "PI internal error code 1004",
        1005: "Controller is busy with some lengthy operation (e.g. reference move, fast scan algorithm)",
        2000: "Controller already has a serial number",
        4000: "Sector erase failed",
        4001: "Flash program failed",
        4002: "Flash read failed",
        4003: "HW match code missing/invalid",
        4004: "FW match code missing/invalid",
        4005: "HW version missing/invalid",
        4006: "FW version missing/invalid",
        4007: "FW update failed",
        5000: "PicoCompensation scan data is not valid",
        5001: "PicoCompensation is running, some actions cannot be executed during scanning/recording",
        5002: "Given axis cannot be defined as PPC axis",
        5003: "Defined scan area is larger than the travel range",
        5004: "Given PicoCompensation type is not defined",
        5005: "PicoCompensation parameter error",
        5006: "PicoCompensation table is larger than maximum table length",
        5100: "Common error in NEXLINE® firmware module",
        5101: "Output channel for NEXLINE® cannot be redefined for other usage",
        5102: "Memory for NEXLINE® signals is too small",
        5103: "RNP cannot be executed if axis is in closed loop",
        5104: "Relax procedure (RNP) needed",
        5200: "Axis must be configured for this action",
        - 1024: "Motion error: position error too large, servo is switched off automatically",
    }

# constants for model number
MODEL_C867 = 867
MODEL_E861 = 861
MODEL_UNKNOWN = 0

# TODO: Have a separate Controller class for each type of move mode. For now, we
# are overriding the method dynamically depending on what the controller can do,
# but that's too hard to read.
class Controller(object):
    def __new__(cls, busacc, address=None, axes=None, *args, **kwargs):
        """
        Takes care of selecting the right class of controller depending on the
        hardware.
        For the arguments, see __init__()
        """
        busacc.flushInput()
        # Three types of controllers: Closed-loop (detected just from the axes
        # arguments), normal open-loop, and open-loop via SMO test command.
        # Difference between the 2 open-loop is hard-coded on the model as it's
        # faster than checking for the list of commands available.
        if address is None:
            subcls = Controller # just for tests/scan
        elif any(axes.values()):
            if not all(axes.values()):
                raise ValueError("Controller %d, mix of closed-loop and "
                                 "open-loop axes is not supported", address)
            # TODO: don't ask user for CL/OL and check if has limit and sensor?
            subcls = CLController
        else:
            # Check controller model by asking it, but cannot rely on the
            # normal commands as nothing is ready, so do all "manually"
            # Note: IDN works even if error is set
            idn = busacc.sendQueryCommand(address, "*IDN?\n")
            if re.search(cls.idn_matches[MODEL_C867], idn):
                subcls = SMOController
            else:
                subcls = OLController

        return super(Controller, cls).__new__(subcls, busacc, address, axes,
                                              *args, **kwargs)
    
    
    def __init__(self, busacc, address=None, axes=None,
                 dist_to_steps=None, min_dist=None):
        """
        busacc: a SerialBusAccesser (with at least a serial port opened)
        address 1<int<16: address as configured on the controller
        If not address is given, it just allows to do some raw commands
        axes (dict int -> boolean): determine which axis will be used and whether
          it will be used closed-loop (True) or open-loop (False).
        Next 2 parameters are calibration values for E-861
        dist_to_steps (0 < float): allows to calibrate how many steps correspond
          to a given distance (in step/m). Default is 1e5, a value that could 
          make sense.
        min_dist (0 <= float < 1): minimum distance required for the axis to 
          even move (in m). Below this distance, a command will be sent, but it
          is expected that the actuator doesn't move at all. Default is 0.01 
          step (= 0.01 / dist_to_steps).
        Next parameter is calibration value for C-867
        vpms (0 < float): calibration value voltage -> speed, in V/(m/s), 
          default is a not too bad value of 87 V/(m/s). Note: it's not linear
          at all actually, but we tend to try to always go at lowest speed (near 3V)
        """
        # TODO: calibration values should be per axis (but for now we only have controllers with 1 axis)
        self.busacc = busacc
        self.address = address
        self._try_recover = False # for now, fully raw access
        # did the user asked for a raw access only?
        if address is None:
            self._channels = set(range(1, 17)) # allow commands to work on any axis
            return
        if axes is None:
            raise LookupError("Need to have at least one axis configured")

        self.GetErrorNum() # make it happy again (in case it wasn't)
        # We don't reboot by default because:
        # * in almost any case we can assume it's in a good state (if not, it's
        #   up to the controller to force it.
        # * if it's in a really bad state, GetErrorNum will fail and cause a
        #   reboot anyway
        # * it sometimes actually cause more harm by putting the controller in
        #   a zombie state
        # * it's slow (especially if you have 5 controllers)

        version = self.GetSyntaxVersion()
        if version != "2.0":
            logging.warning("Controller %d announces untested GCS %s", address, version)

        self._model = self.getModel()
        if self._model == MODEL_UNKNOWN:
            logging.warning("Controller %d is an unsupported version (%s)", self.address, self.GetIdentification())

        self._channels = self.GetAxes() # available channels (=axes)
        # dict axis -> boolean
        self._hasLimitSwitches = dict([(a, self.HasLimitSwitches(a)) for a in self._channels])
        # dict axis -> boolean
        self._hasRefSwitch = dict([(a, self.HasRefSwitch(a)) for a in self._channels])
        self._position = {} # m (dict axis-> position)

        # only for interpolated position (on open-loop)
        self._target = {} # m (dict axis-> future position when a move is over)
        self._end_move = {} # m (dict axis -> time the move will finish)
        self._start_move = {} # m (dict axis -> time the move started)

        # If the controller is mis-configured for the actuator, things can go quite
        # wrong, so make it clear
        for c in self._channels:
            logging.info("Controller %d is configured for actuator %s", address, self.GetStageName(c))
            logging.info("Axis %s has %slimit switches and has %sreference switch",
                         c,
                         "" if self._hasLimitSwitches[c] else "no ",
                         "a " if self._hasRefSwitch[c] else "no ")
        self._avail_params = self.GetAvailableParameters()

        self._try_recover = True # full feature only after init

    def terminate(self):
        pass

    def _sendOrderCommand(self, com):
        """
        Send a command which does not expect any report back
        com (string): command to send (including the \n if necessary)
        """
        self.busacc.sendOrderCommand(self.address, com)

    def _sendQueryCommand(self, com):
        """
        Send a command and return its report
        com (string): the command to send (without address prefix but with \n)
        return (string or list of strings): the report without prefix 
           (e.g.,"0 1") nor newline. If answer is multiline: returns a list of each line 
        """
        try:
            lines = self.busacc.sendQueryCommand(self.address, com)
        except IOError:
            if not self._try_recover:
                raise

            success = self.recoverTimeout()
            if success:
                logging.warning("Controller %d timeout after '%s', but recovered.",
                                self.address, com.encode('string_escape'))
                # try one more time
                lines = self.busacc.sendQueryCommand(self.address, com)
            else:
                raise IOError("Controller %d timeout after '%s', not recovered." %
                              (self.address, com.encode('string_escape')))

        return lines

    err_ans_re = r"(-?\d+)$" # ex: ("0 1 ")[-54](\n)
    def recoverTimeout(self):
        """
        Try to recover from error in the controller state
        return (boolean): True if it recovered
        raise PIGCSError: if the timeout was due to a controller error (in which
            case the controller will be set back to working state if possible)
        """
        self.busacc.flushInput()

        # It makes the controller more comfortable...
        try:
            resp = self.busacc.sendQueryCommand(self.address, "ERR?\n")
            m = re.match(self.err_ans_re, resp)
            if m: # looks like an answer to err?
                err = int(m.group(1))
                if err == 0:
                    return True
                else:
                    raise PIGCSError(err)
        except IOError:
            pass

        # We timed out again, try harder: reboot
        self.Reboot()
        self.busacc.sendOrderCommand("ERR?\n")
        try:
            resp = self.busacc.sendQueryCommand(self.address, "ERR?\n")
            if re.match(self.err_ans_re, resp): # looks like an answer to err?
                # TODO Check if error == 307 or 308?
                return True
        except IOError:
            pass

        # that's getting pretty hopeless
        return False

    # The following are function directly mapping to the controller commands.
    # In general it should not be need to use them directly from outside this class
    def GetIdentification(self):
        #*IDN? (Get Device Identification):
        #ex: 0 2 (c)2010 Physik Instrumente(PI) Karlsruhe,E-861 Version 7.2.0
        version = self._sendQueryCommand("*IDN?\n")
        return version

    def GetSyntaxVersion(self):
        #CSV? (Get Current Syntax Version)
        #GCS version, can be 1.0 (for GCS 1.0) or 2.0 (for GCS 2.0)
        return self._sendQueryCommand("CSV?\n")

    def GetStageName(self, axis):
        """
        return (str) the name of the stage for which the controller is configured.
        Note that the actual stage might be different.
        """
        # CST? does this as well
        #parameter 0x3c
        return self.GetParameter(axis, 0x3C)

    def GetAxes(self):
        """
        returns (set of int): all the available axes
        """
        #SAI? (Get List Of Current Axis Identifiers)
        #SAI? ALL: list all axes (included disabled ones)
        answer = self._sendQueryCommand("SAI? ALL\n")
        # TODO check it works with multiple axes
        axes = set([int(a) for a in answer.split(" ")])
        return axes

    def GetAvailableCommands(self):
        #HLP? (Get List Of Available Commands)
        # first line starts with \x00
        lines = self._sendQueryCommand("HLP?\n")
        lines[0] = lines[0].lstrip("\x00")
        return lines

    def GetAvailableParameters(self):
        """
        Returns the available parameters
        return (dict param -> list of strings): parameter number and strings 
         used to describe it (typically: 0, 1, FLOAT, description)
        """
        #HPA? (Get List Of Available Parameters)
        lines = self._sendQueryCommand("HPA?\n")
        lines[0] = lines[0].lstrip("\x00")
        params = {}
        # first and last lines are typically just user-friendly text
        # look for something like '0x412=\t0\t1\tINT\tmotorcontroller\tI term 1'
        # (and old firmwares report like: '0x412 XXX')
        for l in lines:
            m = re.match(r"0x(?P<param>[0-9A-Fa-f]+)[= ]\w*(?P<desc>(\t?\S+)+)", l)
            if not m:
                logging.debug("Line doesn't seem to be a parameter: '%s'", l)
                continue
            param, desc = int(m.group("param"), 16), m.group("desc")
            params[param] = tuple(filter(bool, desc.split("\t")))
        return params

    def GetParameter(self, axis, param):
        """
        axis (1<int<16): axis number
        param (0<int): parameter id (cf p.35)
        returns (string): the string representing this parameter 
        """
        # SPA? (Get Volatile Memory Parameters)
        assert((1 <= axis) and (axis <= 16))
        assert(0 <= param)
        if hasattr(self, "_avail_params") and  not param in self._avail_params:
            raise ValueError("Parameter %d %d not available" % (axis, param))

        answer = self._sendQueryCommand("SPA? %d %d\n" % (axis, param))
        try:
            value = answer.split("=")[1]
        except IndexError:
            # no "=" => means the parameter is unknown
            raise ValueError("Parameter %d %d unknown" % (axis, param))
        return value

    def SetParameter(self, axis, param, val):
        """
        axis (1<int<16): axis number
        param (0<int): parameter id (cf p.35)
        val (str): value to set (if not a string, it will be converted)
        Raises ValueError if hardware complains
        """
        # SPA (Set Volatile Memory Parameters)
        assert((1 <= axis) and (axis <= 16))
        assert(0 <= param)
        self._sendOrderCommand("SPA %d %d %s\n" % (axis, param, val))
        err = self.GetErrorNum()
        if err:
            raise ValueError("Error %d: setting param 0x%X with val %s failed." %
                             (err, param, val), err)

    def _readAxisValue(self, com, axis):
        """
        Returns the value for a command with axis.
        Ex: POS? 1 -> 1=25.3
        com (str): the 4 letter command (including the ?)
        axis (1<int<16): axis number
        returns (int or float or str): value returned depending on the type detected
        """
        assert(axis in self._channels)
        assert(2 < len(com) < 8)
        resp = self._sendQueryCommand("%s %d\n" % (com, axis))
        try:
            value_str = resp.split("=")[1]
        except IndexError:
            raise ValueError("Failed to parse answer from %s %d: '%s'" %
                             (com, axis, resp))
        try:
            value = int(value_str)
        except ValueError:
            try:
                value = float(value_str)
            except ValueError:
                value = value_str
        
        return value

    def HasLimitSwitches(self, axis):
        """
        Report whether the given axis has limit switches (is able to detect 
         the ends of the axis).
        Note: It's just read from a configuration value in flash 
        memory. Can be configured easily with PIMikroMove 
        axis (1<int<16): axis number
        returns (bool)
        """
        #LIM? (Indicate Limit Switches)
        # 1 => True, 0 => False
        return self._readAxisValue("LIM?", axis) == 1

    def HasRefSwitch(self, axis):
        """
        Report whether the given axis has a reference switch (is able to detect
         the "middle" of the axis).
        Note: apparently it's just read from a configuration value in flash 
        memory. Can be configured easily with PIMikroMove
        axis (1<int<16): axis number
        returns (bool)
        """
        # TODO: Rename to has RefSwitch?
        # TRS? (Indicate Reference Switch)
        # 1 => True, 0 => False
        return self._readAxisValue("TRS?", axis) == 1

    def GetMotionStatus(self):
        """
        returns (set of int): the set of moving axes
        Note: it seems the controller doesn't report moves when using OL via PID
        """
        # "\x05" (Request Motion Status)
        # hexadecimal number bitmap of which axis is moving => 0 if everything is stopped
        # Ex: 4 => 3rd axis moving
        answer = self._sendQueryCommand("\x05")
        bitmap = int(answer, 16)
        # convert to a set
        i = 1
        mv_axes = set()
        while bitmap > 0:
            if bitmap & 1:
                mv_axes.add(i)
            i += 1
            bitmap >>= 1
        return mv_axes

    def GetStatus(self):
        #SRG? = "\x04" (Query Status Register Value)
        #SRG? 1 1
        #Check status
        # hexadecimal number bitmap of which axis is moving => 0 if everything is stopped
        # Ex: 0x9004
        bitmap = self._sendQueryCommand("\x04")
        assert(bitmap.startswith("0x"))
        value = int(bitmap[2:], 16)
        # TODO change to constants
        return value

    def IsReady(self):
        """
        return (boolean): True if ready for new command
        """
        # "\x07" (Request Controller Ready Status)
        # returns 177 if ready, 178 if not
        ans = self._sendQueryCommand("\x07")
        if ans == "\xb1":
            return True
        elif ans == "\xb2":
            return False

        logging.warning("Controller %d replied unknown ready status '%s'", self.address, ans)

    def IsReferenced(self, axis):
        """
        Report whether the given axis has been referenced
        Note: setting position with RON disabled will also put it in this mode
        axis (1<int<16): axis number
        returns (bool)
        """
        #FRF? (Get Referencing Result)
        # 1 => True, 0 => False
        return self._readAxisValue("FRF?", axis) == 1

    def IsOnTarget(self, axis):
        """
        Report whether the given axis is considered on target (for closed-loop 
          moves only)
        axis (1<int<16): axis number
        returns (bool)
        """
        #ONT? (Get On Target State)
        # 1 => True, 0 => False
        # cf parameters 0x3F (settle time), and 0x4D (algo), 0x406 (window size)
        # 0x407 (window off size)
        return self._readAxisValue("ONT?", axis) == 1

    def GetErrorNum(self):
        """
        return (int): the error number (can be negative) of last error
        See p.192 of manual for the error codes
        """
        #ERR? (Get Error Number): get error code of last error
        answer = self._sendQueryCommand("ERR?\n")
        error = int(answer)
        return error

    def Reboot(self):
        self._sendOrderCommand("RBT\n")

        # empty the serial buffer
        self.busacc.flushInput()

        # Sending commands before it's fully rebooted can seriously mess it up.
        # It might end up in a state where only power cycle can reset it.
        # Give it some time to reboot before it's accessible again.
        time.sleep(2)

        self.busacc.flushInput()

    # TODO: use it when terminating?
    def RelaxPiezos(self, axis):
        """
        Call relaxing procedure. Reduce voltage, to increase lifetime and needed
          to change between modes
        axis (1<int<16): axis number
        """
        #RNP (Relax PiezoWalk Piezos): reduce voltage when stopped to increase lifetime
        #Also needed to change between nanostepping and analog
        assert(axis in self._channels)
        self._sendOrderCommand("RNP %d 0\n" % axis)

    def Halt(self, axis=None):
        """
        Stop motion with deceleration
        Note: see Stop
        axis (1<int<16): axis number, 
        """
        #HLT (Stop All Axes): immediate stop (high deceleration != HLT)
        # set error code to 10
        if axis is None:
            self._sendOrderCommand("HLT\n")
        else:
            assert(axis in self._channels)
            self._sendOrderCommand("HLT %d\n" % axis)
#        time.sleep(1) # give it some time to stop before it's accessible again

        # need to recover from the "error", otherwise nothing works
        error = self.GetErrorNum()
        if error != 10: #PI_CNTR_STOP
            logging.warning("Stopped controller %d, but error code is %d instead of 10", self.address, error)

    def Stop(self):
        """
        Stop immediately motion on all axes
        Note: it's not efficient enough with SMO commands
        """
        #STP = "\x18" (Stop All Axes): immediate stop (high deceleration != HLT)
        # set error code to 10
        self._sendOrderCommand("\x18")

        # need to recover from the "error", otherwise nothing works
        error = self.GetErrorNum()
        if error != 10: #PI_CNTR_STOP
            logging.warning("Stopped controller %d, but error code is %d instead of 10", self.address, error)

    def SetServo(self, axis, activated):
        """
        Activate or de-activate the servo. 
        Note: only activate it if there is a sensor (cf .HasRefSwitch and ._hasRefSwitch)
        axis (1<int<16): axis number
        activated (boolean): True if the servo should be activated (closed-loop)
        """
        #SVO (Set Servo State)
        assert(axis in self._channels)

        if activated:
            assert(self._hasRefSwitch[axis])
            state = 1
        else:
            state = 0
        # FIXME: on E861 it seems recommended to first relax piezo.
        # On C867, it's RNP doesn't exists
        self._sendOrderCommand("SVO %d %d\n" % (axis, state))

    def SetReferenceMode(self, axis, absolute):
        """
        Select the reference mode. 
        Note: only useful for closed-loop moves 
        axis (1<int<16): axis number
        absolute (bool): If True, absolute moves can be used, but needs to have
          been referenced.
          If False only relative moves can be used, but only needs a sensor to
          be used.
        """
        #RON (Set Reference Mode)
        assert(axis in self._channels)

        if absolute:
            assert(self._hasLimitSwitches[axis] or self._hasRefSwitch[axis])
            state = 1
        else:
            state = 0
        self._sendOrderCommand("RON %d %d\n" % (axis, state))

    # Functions for relative move in open-loop (no sensor)
    def OLMoveStep(self, axis, steps):
        """
        Moves an axis for a number of steps. Can be done only with servo off.
        axis (1<int<16): axis number
        steps (float): number of steps to do (can be a float). If negative, goes
          the opposite direction. 1 step is about 10µm.
        """
        #OSM (Open-Loop Step Moving): move using nanostepping
        assert(axis in self._channels)
        if steps == 0:
            return
        self._sendOrderCommand("OSM %d %.5g\n" % (axis, steps))

    def SetStepAmplitude(self, axis, amplitude):
        """
        Set the amplitude of one step (in nanostep mode). It affects the velocity
        of OLMoveStep.
        Note: probably it's best to set it to 55 and use OVL to change speed.
        axis (1<int<16): axis number
        amplitude (0<=float<=55): voltage applied (the more the further)
        """
        #SSA (Set Step Amplitude) : for nanostepping
        assert(axis in self._channels)
        assert((0 <= amplitude) and (amplitude <= 55))
        self._sendOrderCommand("SSA %d %.5g\n" % (axis, amplitude))

    def GetStepAmplitude(self, axis):
        """
        Get the amplitude of one step (in nanostep mode).
        Note: mostly just for self-test
        axis (1<int<16): axis number
        returns (0<=float<=55): voltage applied
        """
        #SSA? (Get Step Amplitude), returns something like:
        # 1=10.0000
        assert(axis in self._channels)
        answer = self._sendQueryCommand("SSA? %d\n" % axis)
        amp = float(answer.split("=")[1])
        return amp

    def OLAnalogDriving(self, axis, amplitude):
        """
        Use analog mode to move the axis by a given amplitude.
        axis (1<int<16): axis number
        amplitude (-55<=float<=55): Amplitude of the move. It's only a small move.
          55 is approximately 5 um.
        """
        #OAD (Open-Loop Analog Driving): move using analog
        assert(axis in self._channels)
        assert((-55 <= amplitude) and (amplitude <= 55))
        self._sendOrderCommand("OAD %d %.5g\n" % (axis, amplitude))

    def SetOLVelocity(self, axis, velocity):
        """
        Set velocity for open-loop nanostepping motion.
        axis (1<int<16): axis number
        velocity (0<float): velocity in step-cycles/s. Default is 200 (~ 0.002 m/s)
        """
        #OVL (Set Open-Loop Velocity)
        assert(axis in self._channels)
        assert(velocity > 0)
        self._sendOrderCommand("OVL %d %.5g\n" % (axis, velocity))

    def SetOLAcceleration(self, axis, value):
        """
        Set open-loop acceleration of given axis.
        axis (1<int<16): axis number
        value (0<float): acceleration in step-cycles/s². Default is 2000 
        """
        #OAC (Set Open-Loop Acceleration)
        assert(axis in self._channels)
        assert(value > 0)
        self._sendOrderCommand("OAC %d %.5g\n" % (axis, value))

    def SetOLDeceleration(self, axis, value):
        """
        Set the open-loop deceleration.
        axis (1<int<16): axis number
        value (0<float): deceleration in step-cycles/s². Default is 2000 
        """
        #ODC (Set Open-Loop Deceleration)
        assert(axis in self._channels)
        assert(value > 0)
        self._sendOrderCommand("ODC %d %.5g\n" % (axis, value))

    # Methods for closed-loop functionality. For all of them, servo must be on
    def MoveAbs(self, axis, pos):
        """
        Start an absolute move of an axis to specific position.
         Can only be done with servo on and referenced.
        axis (1<int<16): axis number
        pos (float): position in "user" unit (~m)
        """
        #MOV (Set Target Position)
        assert(axis in self._channels)
        self._sendOrderCommand("MOV %d %.5g\n" % (axis, pos))

    def MoveRel(self, axis, shift):
        """
        Start an relative move of an axis to specific position.
         Can only be done with servo on and referenced.
        axis (1<int<16): axis number
        shift (float): change of position in "user" unit (~m)
        """
        #MVR (Set Target Relative To Current Position)
        assert(axis in self._channels)
        self._sendOrderCommand("MVR %d %.5g\n" % (axis, shift))

    def ReferenceToLimit(self, axis, lim=1):
        """
        Start to move the axis to the switch position (typically, the center)
        Note: Servo and referencing must be on
        See IsReferenced()
        axis (1<int<16): axis number
        lim (-1 or 1): -1 for negative limit and 1 for positive limit
        """
        #FNL (Fast Reference Move To Negative Limit)
        #FPL (Fast Reference Move To Positive Limit)
        assert(axis in self._channels)
        assert(lim in [-1, 1])
        if lim == 1:
            self._sendOrderCommand("FPL %d\n" % axis)
        else:
            self._sendOrderCommand("FNL %d\n" % axis)

    def ReferenceToSwitch(self, axis):
        """
        Start to move the axis to the switch position (typically, the center)
        Note: Servo and referencing must be on
        See IsReferenced()
        axis (1<int<16): axis number
        """
        #FRF (Fast Reference Move To Reference Switch)
        assert(axis in self._channels)
        self._sendOrderCommand("FRF %d\n" % axis)

    def GetPosition(self, axis):
        """
        Get the position (in "user" units)
        axis (1<int<16): axis number
        return (float): pos can be negative
        Note: after referencing, a constant is added by the controller
        """
        #POS? (GetRealPosition)
        return self._readAxisValue("POS?", axis)

    def GetMinPosition(self, axis):
        """
        Get the minimum reachable position (in "user" units)
        axis (1<int<16): axis number
        return (float): pos can be negative
        """
        #TMN? (Get Minimum Commandable Position)
        return self._readAxisValue("TMN?", axis)

    def GetMaxPosition(self, axis):
        """
        Get the maximum reachable position (in "user" units)
        axis (1<int<16): axis number
        return (float): pos can be negative
        """
        #TMX? (Get Maximum Commandable Position)
        assert(axis in self._channels)
        return self._readAxisValue("TMX?", axis)

    def SetCLVelocity(self, axis, velocity):
        """
        Set velocity for closed-loop montion.
        axis (1<int<16): axis number
        velocity (0<float): velocity in units/s
        """
        #VEL (Set Closed-Loop Velocity)
        assert(axis in self._channels)
        assert(velocity > 0)
        self._sendOrderCommand("VEL %d %.5g\n" % (axis, velocity))

    def SetCLAcceleration(self, axis, value):
        """
        Set closed-loop acceleration of given axis.
        axis (1<int<16): axis number
        value (0<float): acceleration in units/s² 
        """
        #ACC (Set Closed-Loop Acceleration)
        assert(axis in self._channels)
        assert(value > 0)
        self._sendOrderCommand("ACC %d %.5g\n" % (axis, value))

    def SetCLDeceleration(self, axis, value):
        """
        Set the closed-loop deceleration.
        axis (1<int<16): axis number
        value (0<float): deceleration in units/s² 
        """
        #DEC (Set Closed-Loop Deceleration)
        assert(axis in self._channels)
        assert(value > 0)
        self._sendOrderCommand("DEC %d %.5g\n" % (axis, value))


# Different from OSM because they use the sensor and are defined in physical unit.
# Servo must be off! => Probably useless... compared to MOV/MVR
#OMR (Relative Open-Loop Motion)
#OMA (Absolute Open-Loop Motion)
#

    # Below are methods for manipulating the controller
    idn_matches = {
               MODEL_C867: "Physik Instrumente.*,.*C-867",
               MODEL_E861: "Physik Instrumente.*,.*E-861", 
               }
    def getModel(self):
        """
        returns a model constant
        """
        idn = self.GetIdentification()
        for c, m in self.idn_matches.items():
            if re.search(m, idn):
                return c
        return MODEL_UNKNOWN

    def _storeMove(self, axis, shift, duration):
        """
        Save move information for interpolating the position
        To be called when a new move is started
        axis (int): the channel
        shift (float): relative change in position (im m)
        duration (0<float): time it will take (in s)
        """
        now = time.time()
        cur_pos = self._interpolatePosition(axis)
        self._position[axis] = cur_pos
        self._start_move[axis] = now
        self._end_move[axis] = now + duration
        self._target[axis] = cur_pos + shift

    def _storeStop(self, axis):
        """
        Save the fact that a move was stop immediately (maybe not achieved)
        """
        self._position[axis] = self._interpolatePosition(axis)
        self._target[axis] = self._position[axis]
        self._end_move[axis] = 0

    def _storeMoveComplete(self, axis):
        """
        Save the fact that the current move is complete (even if the end time
        is not yet achieved)
        """
        self._position[axis] = self._target.get(axis, self._position[axis])
        self._end_move[axis] = 0

    def _interpolatePosition(self, axis):
        """
        return (float): interpolated position at the current time
        """
        now = time.time()
        end_move = self._end_move.get(axis, 0)
        if now > end_move:
            target = self._target.get(axis, self._position[axis])
            logging.debug("Interpolating move by reporting target position: %g",
                           target)
            return target
        else:
            start = self._start_move[axis]
            completion = (now - start) / (end_move - start)
            pos = self._position[axis]
            cur_pos = pos + (self._target[axis] - pos) * completion
            logging.debug("Interpolating move to %g %% of complete move: %g",
                          completion * 100, cur_pos)
            return cur_pos

    def getPosition(self, axis):
        """
        Note: in open-loop mode it's very approximate (and interpolated)
        return (float): the current position of the given axis
        """
        # This is using interpolation, closed-loop must override this method
        assert(axis in self._channels)

        # make sure that if a move finished early, we report the final position
        if not self.isMoving(set([axis])):
            self._storeMoveComplete(axis)

        return self._interpolatePosition(axis)

    def setSpeed(self, axis, speed):
        """
        Changes the move speed of the motor (for the next move).
        Note: in open-loop mode, it's very approximate.
        speed (0<float<10): speed in m/s.
        axis (1<=int<=16): the axis
        """
        assert (self.min_speed <= speed <= self.max_speed)
        assert (axis in self._channels)
        self._speed[axis] = speed

    def getSpeed(self, axis):
        return self._speed[axis]

    def setAccel(self, axis, accel):
        """
        Changes the move acceleration (and deceleration) of the motor (for the next move).
        Note: in open-loop mode, it's very approximate.
        accel (0<float<100): acceleration in m/s².
        axis (1<=int<=16): the axis
        """
        assert (0 < accel <= self.max_accel)
        assert (axis in self._channels)
        self._accel[axis] = accel

    def moveRel(self, axis, distance):
        """
        Move on a given axis for a given distance.
        It's asynchronous: the method might return before the move is complete.
        axis (1<=int<=16): the axis
        distance (float): the distance of move in m (can be negative)
        returns (float): approximate distance actually moved
        """
        raise NotImplementedError("This method must be overridden by a subclass")

    def moveAbs(self, axis, position):
        """
        Move on a given axis to a given position.
        It's asynchronous: the method might return before the move is complete.
        axis (1<=int<=16): the axis
        position (float): the target position in m (can be negative)
        returns (float): approximate distance actually moved
        """
        # This is a default implementation relying on the moveRel
        # It should be overriden, if it can be done directly
        old_pos = self.getPosition(axis)
        shift = position - old_pos
        return self.moveRel(axis, shift)

    def isMoving(self, axes=None):
        """
        Indicate whether the motors are moving. 
        axes (None or set of int): axes to check whether for move, or all if None
        return (boolean): True if at least one of the axes is moving, False otherwise
        """
        if axes is None:
            axes = self._channels
        else:
            assert axes.issubset(self._channels)

        return not axes.isdisjoint(self.GetMotionStatus())

    def stopMotion(self):
        """
        Stop the motion on all axes immediately
        """
        self.Stop()

    def waitEndMotion(self, axes=None):
        """
        Wait until the motion of all the given axis is finished.
        Note: there is a 5 s timeout
        axes (None or set of int): axes to check whether for move, or all if None
        """
        #TODO use the time, distance, and speed of last move to evaluate the timeout
        # approximately the time for the longest move
        timeout = 5 #s
        end = time.time() + timeout
        while self.isMoving(axes):
            if time.time() >= end:
                raise IOError("Timeout while waiting for end of motion")
            time.sleep(0.005)

    def selfTest(self):
        """
        check as much as possible that it works without actually moving the motor
        return (boolean): False if it detects any problem
        """
        try:
            error = self.GetErrorNum()
            if error:
                logging.warning("Controller %d had error status %d", self.address, error)

            version = self.GetSyntaxVersion()
            logging.info("GCS version: '%s'", version)
            ver_num = float(version)
            if ver_num < 1 or ver_num > 2:
                logging.error("Controller %d has unexpected GCS version %s", self.address, version)
                return False

            axes = self.GetAxes()
            if len(axes) == 0 or len(axes) > 16:
                logging.error("Controller %d report axes %s", self.address, str(axes))
                return False

            if self._model in (MODEL_E861,): # support open-loop mode
                for a in self._channels:
                    self.SetStepAmplitude(a, 10)
                    amp = self.GetStepAmplitude(a)
                    if amp != 10:
                        logging.error("Failed to modify amplitude of controller %d (%f instead of 10)", self.address, amp)
                        return False

            if self._model in (MODEL_C867,): # support temperature reading
                # No support for direct open-loop mode
                # TODO put the temperature as a RO VA?
                current_temp = float(self.GetParameter(1, 0x57))
                max_temp = float(self.GetParameter(1, 0x58))
                if current_temp >= max_temp:
                    logging.error("Motor of controller %d too hot (%f C)", self.address, current_temp)
                    return False
        except Exception:
            return False

        return True

    @staticmethod
    def scan(busacc, max_add=16):
        """
        Scan the serial network for all the PI GCS compatible devices available.
        Note this is the low-level part, you probably want to use Controller.scan()
         for scanning devices on a computer.
        bus: the bus
        max_add (1<=int<=16): maximum address to scan
        return (dict int -> tuple): addresses of available controllers associated
            to number of axes, and presence of limit switches/sensor
        """
        ctrl = Controller(busacc)

        present = {}
        for i in range(1, max_add + 1):
            # ask for controller #i
            logging.debug("Querying address %d", i)

            # is it answering?
            try:
                ctrl.address = i
                axes = {}
                for a in ctrl.GetAxes():
                    axes = {a: ctrl.HasRefSwitch(a)}
                if not axes:
                    logging.info("Found controller %d with no axis", i)
                else:
                    version = ctrl.GetIdentification()
                    logging.info("Found controller %d with ID '%s'.", i, version)
                    present[i] = axes
            except IOError:
                pass

        ctrl.address = None
        return present

class CLController(Controller):
    """
    Controller managed via closed-loop commands (ex: C-867 with reference)
    TODO: For now, only relative moves are supported.
          For supporting absolute moves, we need to add querying and requesting
          "homing" procedure. Then the position would reset to 0 (and that's it
          from the user's point of view).
    """
    def __init__(self, busacc, address=None, axes=None):
        super(CLController, self).__init__(busacc, address, axes)

        self._speed = {} # m/s dict axis -> speed
        self._accel = {} # m/s² dict axis -> acceleration/deceleration
        self.pos_rng = {} # m, dict axis -> min,max position

        for a, cl in axes.items():
            if not a in self._channels:
                raise LookupError("Axis %d is not supported by controller %d" % (a, address))

            if not cl: # want open-loop?
                raise ValueError("Initialising CLController with request for open-loop")
            if not self._hasRefSwitch[a]:
                raise ValueError("Closed-loop control requested but controller "
                                 "%d reports no sensor for axis %d" % (address, a))


            # TODO:
            # * if not referenced => disable reference mode to be able to
            # move relatively (until referencing happens). Use the position as
            # is, so that if no referencing ever happens, at least the position
            # is correctly as long as the controller is powered.
            # * if referenced => beleive it and stay in this mode.

            # Start with servo, non-referenced mode, assume the position is
            # correct, and set units to meters
            self.SetServo(a, True)
            self.SetReferenceMode(a, False)
            # It seems the internal device unit is in 100 nm, but "user" unit
            # is usually set to report in mm (but can be changed).
            # => change to be always in meters
            self.SetParameter(a, 0xE, 10000000) # numerator
            self.SetParameter(a, 0xF, 1) # denumerator

            # Movement range before refenrencing is max range in both directions
            pos = self.GetPosition(a)
            width = self.GetMaxPosition(a) - self.GetMinPosition(a)
            self.pos_rng[a] = (pos - width, pos + width)

            # Read speed/accel ranges
            self._speed[a] = self._readAxisValue("VEL?", a) # m/s
            self._accel[a] = self._readAxisValue("ACC?", a) # m/s²

            # TODO: also use per-axis info
            try:
                self.max_speed = float(self.GetParameter(1, 0xA)) # m/s
                self.max_accel = float(self.GetParameter(1, 0x4A)) # m/s²
            except (IOError, ValueError):
                self.max_speed = self._speed[a]
                self.max_accel = self._accel[a]

        self.min_speed = 10e-6 # m/s (default low value)
        self._prev_speed_accel = ({}, {})

        # Note: setting position only works if ron is disabled. It's possible
        # also indirectly set it after referencing, but then it will conflict
        # with TMN/TMX and some correct moves will fail.
        # So referencing could look like:
        # ron 1 1
        # frf -> go home and now know the position officially
        # orig_pos = pos?
        #

        # Use FRF? to know about the referencing status
        # TMN? TMX?  give range (careful, it's typically 0 -> BIG Number)
        # Note: the actuator actively control the motor to stay at the last 
        # position requested as long as the servo is on. So need to disable
        # servo when stopping the program to allow the user to move the stage
        # manually (or maybe even after any move?)

    def terminate(self):
        super(CLController, self).terminate()
        # Disable servo, to allow the user to move the axis manually
        for a in self._channels:
            self.SetServo(a, False)

    def _updateSpeedAccel(self, axis):
        """
        Update the speed and acceleration values for the given axis. 
        It's only done if necessary, and only for the current closed- or open-
        loop mode.
        axis (1<=int<=16): the axis
        """
        prev_speed = self._prev_speed_accel[0].get(axis, None)
        new_speed = self._speed[axis]
        if prev_speed != new_speed:
            # TODO: check it's within range
            self.SetCLVelocity(axis, new_speed)
            self._prev_speed_accel[0][axis] = new_speed

        prev_accel = self._prev_speed_accel[1].get(axis, None)
        new_accel = self._accel[axis]
        if prev_accel != new_accel:
            # TODO: check it's within range
            self.SetCLAcceleration(axis, new_accel)
            self.SetCLDeceleration(axis, new_accel)
            self._prev_speed_accel[1][axis] = new_accel

    def moveRel(self, axis, distance):
        """
        See Controller.moveRel
        """
        assert(axis in self._channels)
        self._updateSpeedAccel(axis)
        # We trust the caller that it knows it's in range
        # (worst case the hardware will not go further)
        self.MoveRel(axis, distance)

        return distance

    def moveAbs(self, axis, position):
        """
        See Controller.moveAbs
        """
        assert(axis in self._channels)
        self._updateSpeedAccel(axis)
        # We trust the caller that it knows it's in range
        # (worst case the hardware will not go further)
        old_pos = self.getPosition(axis)
        distance = position - old_pos

        # Absolute move is only legal if already referenced.
        if not self.IsReferenced(axis): # TODO: cache, or just always do relative?
            self.MoveRel(axis, distance)
        else:
            self.MoveAbs(axis, position)

        return distance

    def getPosition(self, axis):
        """
        Find current position as reported by the sensor
        return (float): the current position of the given axis
        """
        return self.GetPosition(axis)

    # Warning: if the settling window is too small or settling time too big,
    # it might take several seconds to reach target (or even never reach it)
    def isMoving(self, axes=None):
        """
        Indicate whether the motors are moving (ie, last requested move is over)
        axes (None or set of int): axes to check whether for move, or all if None
        return (boolean): True if at least one of the axes is moving, False otherwise
        """
        if axes is None:
            axes = self._channels
        else:
            assert axes.issubset(self._channels)

        # With servo on, it might constantly be _slightly_ moving (around the
        # target), so it's much better to use IsOnTarget info. The controller
        # needs to be correctly configured with the right window size.
        for c in axes:
            if not self.IsOnTarget(c):
                return True
        return False

    def startReferencing(self, axis):
        """
        Start a referencing move. Use isMoving() or isReferenced to know if
        the move is over. Position will change, as well as absolute positions.
        axis (1<=int<=16)
        """
        if self._hasRefSwitch[axis]:
            self.ReferenceToSwitch(axis)
        else:
            raise NotImplementedError("Don't know how to reference to limit yet")
            self.ReferenceToLimit(axis)
            # TODO: need to do that after the move is complete
            self.waitEndMotion(set(axis))
            # Go to 0 (="home")
            self.MoveAbs(axis, 0)

    def isReferenced(self, axis):
        """
        returns (bool): True if the axis is referenced
        """
        return self.IsReferenced(axis)


class OLController(Controller):
    """
    Controller managed via open-loop commands (ex: E-861)
    """
    def __init__(self, busacc, address=None, axes=None,
                 dist_to_steps=None, min_dist=None):
        if dist_to_steps and not (0 < dist_to_steps):
            raise ValueError("dist_to_steps (%s) must be > 0" % dist_to_steps)
        if min_dist and not (0 <= min_dist < 1):
            raise ValueError("min_dist (%s) must be between 0 and 1 m" % min_dist)

        super(OLController, self).__init__(busacc, address, axes,
                                           dist_to_steps=dist_to_steps,
                                           min_dist=min_dist)
        for a, cl in axes.items():
            if not a in self._channels:
                raise LookupError("Axis %d is not supported by controller %d" % (a, address))

            if cl: # want closed-loop?
                raise ValueError("Initialising OLController with request for closed-loop")
            # that should be the default, but for safety we force it
            self.SetServo(a, False)
            self.SetStepAmplitude(a, 55) # maximum is best
            self._position[a] = 0


        # TODO: allow to pass a polynomial
        self._dist_to_steps = dist_to_steps or 1e5 # step/m
        if min_dist is None:
            self.min_stepsize = 0.01 # step, under this, no move at all
        else:
            self.min_stepsize = min_dist * self._dist_to_steps

        self.min_speed = 10e-6 # m/s (default low value)

        # FIXME 0x7000204 seems specific to E-861. => use CL info if not available?
        self.max_speed = 0.5 # m/s
        self.max_accel = 0.01 # m/s²
        try:
            # (max m/s) = (max step/s) * (step/m)
            self.max_speed = float(self.GetParameter(1, 0x7000204)) / self._dist_to_steps # m/s
            # Note: the E-861 claims max 0.015 m/s but actually never goes above 0.004 m/s
            # (max m/s²) = (max step/s²) * (step/m)
            self.max_accel = float(self.GetParameter(1, 0x7000205)) / self._dist_to_steps # m/s²
        except (IOError, ValueError) as err:
            # TODO detect better that it's just a problem of sending unsupported command/value
            # Put default (large values)
            self.GetErrorNum() # reset error (just in case)
            logging.debug("Using default speed and acceleration value after error '%s'", err)

        self._speed = dict([(a, (self.min_speed + self.max_speed) / 2) for a in axes]) # m/s
        self._accel = dict([(a, self.max_accel) for a in axes]) # m/s² (both acceleration and deceleration)
        self._prev_speed_accel = ({}, {})

    def _convertDistanceToDevice(self, distance):
        """
        converts meters to the unit for this device (steps) in open-loop.
        distance (float): meters (can be negative)
        return (float): number of steps, <0 if going opposite direction
            0 if too small to move.
        """
        steps = distance * self._dist_to_steps
        if abs(steps) < self.min_stepsize:
            return 0

        return steps

    def _convertSpeedToDevice(self, speed):
        """
        converts meters/s to the unit for this device (steps/s) in open-loop.
        distance (float): meters/s (can be negative)
        return (float): number of steps/s, <0 if going opposite direction
        """
        steps_ps = speed * self._dist_to_steps
        return max(1, steps_ps) # don't go at 0 m/s!

    # in linear approximation, it's the same
    _convertAccelToDevice = _convertSpeedToDevice

    def _updateSpeedAccel(self, axis):
        """
        Update the speed and acceleration values for the given axis. 
        It's only done if necessary, and only for the current closed- or open-
        loop mode.
        axis (1<=int<=16): the axis
        """
        prev_speed = self._prev_speed_accel[0].get(axis, None)
        new_speed = self._speed[axis]
        if prev_speed != new_speed:
            steps_ps = self._convertSpeedToDevice(new_speed)
            self.SetOLVelocity(axis, steps_ps)
            self._prev_speed_accel[0][axis] = new_speed

        prev_accel = self._prev_speed_accel[1].get(axis, None)
        new_accel = self._accel[axis]
        if prev_accel != new_accel:
            steps_pss = self._convertAccelToDevice(new_accel)
            self.SetOLAcceleration(axis, steps_pss)
            self.SetOLDeceleration(axis, steps_pss)
            self._prev_speed_accel[1][axis] = new_accel

    def moveRel(self, axis, distance):
        """
        See Controller.moveRel
        """
        # TODO: also report expected time for the move?
        assert(axis in self._channels)

        self._updateSpeedAccel(axis)
        steps = self._convertDistanceToDevice(distance)
        if steps == 0: # if distance is too small, report it
            return 0
            # TODO: try to move anyway, just in case it works

        self.OLMoveStep(axis, steps)
        # TODO use OLAnalogDriving for very small moves (< 5µm)?

        duration = abs(distance) / self._speed[axis]
        self._storeMove(axis, distance, duration)
        return distance

    def stopMotion(self):
        # Controller.stopMotion(self)
        super(OLController, self).stopMotion()
        for c in self._channels:
            self._storeStop(c)

class SMOController(Controller):
    """
    Controller managed via the test open-loop command "SMO" (ex: C-867)
    """
    def __init__(self, busacc, address=None, axes=None, vmin=2., speed_base=0.03):
        """
        vmin (0.5 < float < 10): lowest voltage at which the actuator moves 
          reliably in V. This is the voltage used when performing small moves
          (~< 50 µm).
        speed_base (0<float<10): speed in m/s at the base voltage (3.5V). The
          base voltage is used for long moves (~> 50 µm). 
        """
        # TODO: need 4 settings:
        # vmin/speed_min: voltage/speed for the smallest moves (will be used < 50 µm)
        # vmax/speed_max: voltage/speed for the big moves (will be used above)

        if not (0.5 < vmin < 10):
            raise ValueError("vmin (%s) must be between 0.5 and 10 V" % vmin)
        if not (0 < speed_base < 10):
            raise ValueError("speed_base (%s) must be between 0 and 10 m/s" % speed_base)

        super(SMOController, self).__init__(busacc, address, axes)
        for a, cl in axes.items():
            if not a in self._channels:
                raise LookupError("Axis %d is not supported by controller %d" % (a, address))

            if cl: # want closed-loop?
                raise ValueError("Initialising OLController with request for closed-loop")
            # that should be the default, but for safety we force it
            self.SetServo(a, False)
            self._position[a] = 0

        # Get maximum motor output parameter (0x9) allowed
        # Because some type of stages cannot bear as much as the full maximum
        # The maximum output voltage is calculated following this formula:
        # 200 Vpp*Maximum motor output/32767
        self._max_motor_out = int(self.GetParameter(1, 0x9))
        # official approx. min is 3V, but from test, it can sometimes go down to 1.5V.
        self._min_motor_out = int((vmin / 10) * 32767) # encoded as a ratio of 10 V * 32767
        if self._max_motor_out < self._min_motor_out:
            raise ValueError("Controller report max voltage lower than vmin=%g V" % vmin)

        # FIXME: the manual somehow implies that writing macros writes on the
        # flash, which is only possible a number of times. So only write macro
        # if it's different?

        # Set up a macro that will do the job
        # To be called like "MAC START OS 16000 500"
        # First param is voltage between -32766 and 32766
        # Second param is delay in ms between 1 and 9999
        # Note: it seems it doesn't work to have a third param is the axis
        # WARNING: old firmware (<1.2) don't support macro arguments, but there
        # is not a clear way to detect them (but checking firmware version)
        # Note: macro name is short to be sure to have a short command line
        mac = "MAC BEG OS\n" \
              "%(n)d SMO 1 $1\n" \
              "%(n)d DEL $2\n"   \
              "%(n)d SMO 1 0\n"  \
              "%(n)d MAC END\n" % {"n": self.address}
        self._sendOrderCommand(mac)

#        # For short moves
#        mac = "MAC BEG OS0\n" \
#              "%(n)d SMO 1 $1\n" \
#              "%(n)d SAI? ALL\n" \
#              "%(n)d SMO 1 0\n"  \
#              "%(n)d MAC END\n" % {"n": self.address}
#        self._sendOrderCommand(mac)

        # Don't authorize different speeds or accels
        self._speed_base = speed_base
        self.min_speed = speed_base # m/s
        self.max_speed = speed_base # m/s
        self.max_accel = 0.01 # m/s² (actually I've got no idea)

        self._speed = dict([(a, self.min_speed) for a in axes]) # m/s
        self._accel = dict([(a, self.max_accel) for a in axes]) # m/s² (both acceleration and deceleration)

    def StopOLViaPID(self, axis):
        """
        Stop the fake PID driving when doing open-loop
        """
        self._sendOrderCommand("SMO %d 0\n" % axis)

    def OLMovePID(self, axis, voltage, t):
        """
        Moves an axis for a number of steps. Can be done only with servo off.
        axis (1<int<16): axis number
        voltage (-32766<=int<=32766): voltage for the PID control. <0 to go towards
          the negative direction. 32766 is 10V
        t (0<int <= 9999): time in ms.
        """
        # Uses MAC OS, based on SMO
        assert(axis == 1) # seems not possible to have 3 parameters?!
        assert(-32768 <= voltage <= 32767)
        assert(0 < t <= 9999)

        # From experiment: a delay of 0 means actually 2**16, and >= 10000 it's 0
        self._sendOrderCommand("MAC START OS %d %d\n" % (voltage, t))

#    def OLMovePID0(self, axis, voltage):
#        """
#        Moves an axis a very little bit. Can be done only with servo off.
#        Warning: it's completely hacky, there is no idea if it even moves
#        axis (1<int<16): axis number
#        voltage (-32766<=int<=32766): voltage for the PID control. <0 to go towards
#          the negative direction. 32766 is 10V
#        """
#        # Uses MAC OS0, based on SMO
#        assert(axis == 1)
#        assert(-32768 <= voltage <= 32767)
#
#        self._sendOrderCommand("MAC START OS0 %d\n" % (voltage,))


    def _isAxisMovingOLViaPID(self, axis):
        """
        axis (1<int<16): axis number
        returns (boolean): True moving axes for the axes controlled via PID
        """
#        # FIXME: SMO doesn't seem to work any more with the fine grain firmware
#        now = time.time()
#        end_move = self._end_move.get(axis, 0)
#        if now > end_move:
#            return False
#        else:
#            return True

        # "SMO?" (Get Control Value)
        # Reports the speed set. If it's 0, it's not moving, otherwise, it is.
        answer = self._sendQueryCommand("SMO? %d\n" % axis)
        value = answer.split("=")[1]
        if value == "0":
            return False
        else:
            return True

    cycles_per_s = 20000 # 20 kHz for new experimental firmware
    #cycles_per_s = 1000 # 1 kHz for normal firmware

    _base_motor_out = int((3.5 / 10) * 32767) # High enough to ensures it's never stuck

    def _convertDistanceSpeedToPIDControl(self, distance, speed):
        """
        converts meters and speed to the units for this device (~V, ms) in
        open-loop via PID control.
        distance (float): meters (can be negative)
        speed (0<float): meters/s UNUSED
        return (tuple: int, 0<number, float): PID control (in device unit),
          duration (in device cycles), expected actual distance (in m)
        """
        # TODO: smooth transition of voltage between 50µm and 20µm?
        if abs(distance) > 50e-6: # big move
            # Use the fastest speed
            voltage_u = max(self._min_motor_out, self._base_motor_out)
            speed = self._speed_base # m/s
        else: # small move
            # Consider the speed linear of the voltage with:
            # * vmin - 0.5 -> 0 m/s
            # * base voltage -> speed_base
            # * vmin -> 0.5 * speed_base / (bv - vmin +0.5) 
            # It's totally wrong, but approximately correct
            voltage_u = self._min_motor_out
            bv = self._base_motor_out * (10 / 32767)
            vmin = self._min_motor_out * (10 / 32767)
            speed = 0.5 * self._speed_base / (bv - vmin + 0.5)

        mv_time = abs(distance) / speed # s
        mv_time_cy = int(round(mv_time * self.cycles_per_s)) # cycles
        if mv_time_cy < 1:
            # really no hope
            return 0, 0, 0
        elif mv_time_cy < 5:
            # On fine grain firmware, below 6 cyles gives always the same time
            mv_time_cy = 5
        elif mv_time_cy >= 10000:
            logging.debug("Too big distance of %f m, shortening it", distance)
            mv_time_cy = 9999

        act_dist = (mv_time_cy / self.cycles_per_s) * speed # m (very approximate)

        if distance < 0:
            voltage_u = -voltage_u
            act_dist = -act_dist

        return voltage_u, mv_time_cy, act_dist

#    def _OLDconvertDistanceSpeedToPIDControl(self, distance, speed):
#        """
#        converts meters and speed to the units for this device (~V, ms) in
#        open-loop via PID control.
#        distance (float): meters (can be negative)
#        speed (0<float): meters/s (can be negative)
#        return (tuple: int, 0<number, float): PID control (in device unit),
#          duration (in device cycles), distance (in m)
#        """
#        voltage_u = round(int(speed * self._vpms)) # uV
#        # clamp it to the possible values
#        voltage_u = min(max(self._min_motor_out, voltage_u), self._max_motor_out)
#        act_speed = voltage_u / self._vpms # m/s
#
#        mv_time = abs(distance) / act_speed # s
#        mv_time_cy = int(round(mv_time * self.cycles_per_s)) # cycles
#        if mv_time < 1e-3 and act_speed > self.min_speed:
#            # small distance => try with the minimum speed to have a better precision
#            return self._convertDistanceSpeedToPIDControl(distance, self.min_speed)
#        # FIXME: unused with new firmware
#        elif 0.1e-3 < mv_time and mv_time_cy < 1:
#            logging.error("Should not be called")
#            # try our special super small step trick if at least 0.1 ms
#            mv_time_cy = 0.5e-3 * self.cycles_per_s # 0.5 ms
#            voltage_u = self._min_motor_out # TODO: change according to requested distance?
#        # TODO: on new experimental firmware, below 6 cyles might not be possible
#        # according to PI
#        elif mv_time_cy < 1:
#            # really no hope
#            return 0, 0, 0
#        elif mv_time_cy >= 10000:
#            logging.debug("Too big distance of %f m, shortening it", distance)
#            mv_time_cy = 9999
#
#        if distance < 0:
#            voltage_u = -voltage_u
#
#        act_dist = (mv_time_cy / self.cycles_per_s) * voltage_u / self._vpms # m (very approximate)
#        return voltage_u, mv_time_cy, act_dist

    def moveRel(self, axis, distance):
        """
        See Controller.moveRel
        """
        assert(axis in self._channels)
        
        speed = self._speed[axis]
        v, t, ad = self._convertDistanceSpeedToPIDControl(distance, speed)
        if t == 0: # if distance is too small, report it
            logging.debug("Move of %g µm too small, not moving", distance * 1e-6)
            return 0
#        elif t < 1: # special small move command
#            logging.error("Should never be called")
#            self.OLMovePID0(axis, v)
#            duration = 0.5e-3 # = 1/2 ms
        else:
            self.OLMovePID(axis, v, t)
            duration = t / self.cycles_per_s
        logging.debug("Moving axis at %f V, for %f ms", v * (10 / 32767), duration)

        self._storeMove(axis, ad, duration)
        return ad

    def isMoving(self, axes=None):
        """
        See Controller.isMoving
        """
        if axes is None:
            axes = self._channels
        else:
            assert axes.issubset(self._channels)

        for c in self._channels:
            if self._isAxisMovingOLViaPID(c):
                return True
        return False

    def stopMotion(self):
        """
        Stop the motion on all axes immediately
        Implementation for open-loop PID control 
        """
        self.Stop() # doesn't seem to be effective with SMO macro
        for c in self._channels:
            self.StopOLViaPID(c)
            self._storeStop(c)

class Bus(model.Actuator):
    """
    Represent a chain of PIGCS controllers over a serial port
    """
    def __init__(self, name, role, port, axes, baudrate=38400,
                 dist_to_steps=None, min_dist=None,
                 vmin=None, speed_base=None,
                 _addresses=None, **kwargs):
        """
        port (string): name of the serial port to connect to the controllers
         (starting with /dev on Linux or COM on windows) or "autoip" for
         automatically finding an ip address or a "host[:port]".
        axes (dict string -> 3-tuple(1<=int<=16, 1<=int, boolean): the configuration
         of the network. For each axis name associates the controller address,
         channel, and whether it's closed-loop (absolute positioning) or not.
         Note that even if it's made of several controllers, each controller is 
         _not_ seen as a child from the odemis model point of view.
        baudrate (int): baudrate of the serial port (default is the recommended 
          38400). Use .scan() to detect it.
        Next 3 parameters are for calibration, see Controller for definition
        dist_to_steps (dict string -> (0 < float)): axis name -> value
        min_dist (dict string -> (0 <= float < 1)): axis name -> value
        vpms (dict string -> (0 < float)): axis name -> value
        """
        self.accesser = self._openPort(port, baudrate, _addresses)

        dist_to_steps = dist_to_steps or {}
        min_dist = min_dist or {}
        vmin = vmin or {}
        speed_base = speed_base or {}

        # Prepare initialisation by grouping axes from the same controller
        ac_to_axis = {} # address, channel -> axis name
        controllers = {} # address -> kwargs (axes, dist_to_steps, min_dist, vpms...)
        for axis, (add, channel, isCL) in axes.items():
            if not add in controllers:
                controllers[add] = {"axes":{}}
            elif channel in controllers[add]:
                raise ValueError("Cannot associate multiple axes to controller %d:%d" % (add, channel))
            ac_to_axis[(add, channel)] = axis
            kwc = controllers[add]
            kwc["axes"].update({channel: isCL})
            # FIXME: for now we rely on the fact 1 axis = 1 controller for the calibration values
            if axis in dist_to_steps:
                kwc["dist_to_steps"] = dist_to_steps[axis]
            if axis in min_dist:
                kwc["min_dist"] = min_dist[axis]
            if axis in vmin:
                kwc["vmin"] = vmin[axis]
            if axis in speed_base:
                kwc["speed_base"] = speed_base[axis]

        # Init each controller
        self._axis_to_cc = {} # axis name => (Controller, channel)
        axes_def = {} # axis name => axis definition
        # TODO also a rangesRel : min and max of a step
        position = {}
        speed = {}
        referenced = {}
        for address, kwc in controllers.items():
            try:
                controller = Controller(self.accesser, address, **kwc)
            except IOError:
                logging.exception("Failed to find a controller with address %d on %s", address, port)
                raise
            except LookupError:
                logging.exception("Failed to initialise controller %d on %s", address, port)
                raise
            channels = kwc["axes"]
            for c, isCL in channels.items():
                axis = ac_to_axis[(address, c)]
                self._axis_to_cc[axis] = (controller, c)

                position[axis] = controller.getPosition(c)
                # TODO if closed-loop, the ranges should be updated after homing
                try:
                    rng = controller.pos_rng[c]
                except (IndexError, AttributeError):
                    # Unknown? Give room
                    rng = (-1, 1) # m
                # Just to make sure it doesn't go too fast
                speed[axis] = 0.001 # m/s
                speed_rng = (controller.min_speed, controller.max_speed)
                ad = model.Axis(unit="m", range=rng, speed=speed_rng,
                                canAbs=isCL)
                axes_def[axis] = ad
                if ad.canAbs:
                    referenced[axis] = controller.isReferenced(c)

        # this set ._axes
        model.Actuator.__init__(self, name, role, axes=axes_def, **kwargs)

        # TODO: allow to override the unit (per axis)
        # RO, as to modify it the client must use .moveRel() or .moveAbs()
        self.position = model.VigilantAttribute(self._applyInversionAbs(position),
                                                unit="m", readonly=True)

        # RO VA dict axis -> bool: True if the axis has been referenced
        # Only axes which can be referenced are listed
        self.referenced = model.VigilantAttribute(referenced, readonly=True)

        # min speed = don't be crazy slow. max speed from hardware spec
        self.speed = model.MultiSpeedVA(speed, range=[0, 10.],
                                        unit="m/s", setter=self._setSpeed)
        self._setSpeed(speed)

        # set HW and SW version
        self._swVersion = self.accesser.driverInfo
        hwversions = []
        for axis, (ctrl, channel) in self._axis_to_cc.items():
            hwversions.append("'%s': %s (GCS %s) for %s" %
                              (axis, ctrl.GetIdentification(),
                               ctrl.GetSyntaxVersion(), ctrl.GetStageName(channel))
                             )
        self._hwVersion = ", ".join(hwversions)

        self._action_mgr = ActionManager(self)
        self._action_mgr.start()
    
    def _updatePosition(self):
        """
        update the position VA
        Note: it should not be called while holding the lock to the serial port
        """
        position = {}
        # TODO: improve efficiency (don't ask all controllers every times)
        # request position from each controller
        for axis, (controller, channel) in self._axis_to_cc.items():
            position[axis] = controller.getPosition(channel)

        pos = self._applyInversionAbs(position)
        logging.debug("Reporting new position at %s", pos)

        # it's read-only, so we change it via _value
        self.position._value = pos
        self.position.notify(self.position.value)

    def _setSpeed(self, value):
        """
        value (dict string-> float): speed for each axis
        returns (dict string-> float): the new value
        """
        for axis, v in value.items():
            rng = self._axes[axis].speed
            if not rng[0] <= v <= rng[1]:
                raise ValueError("Requested speed of %f for axis %s not within %f->%f" %
                                 (v, axis, rng[0], rng[1]))
            controller, channel = self._axis_to_cc[axis]
            controller.setSpeed(channel, v)
        return value

    @isasync
    def moveRel(self, shift):
        """
        Move the stage the defined values in m for each axis given.
        shift dict(string-> float): name of the axis and shift in m
        returns (Future): future that control the asynchronous move
        """
        self._checkMoveRel(shift)
        shift = self._applyInversionRel(shift)
        # converts the request into one action (= a dict controller -> channels + distance)
        action_axes = collections.defaultdict(list)
        for axis, distance in shift.items():
            controller, channel = self._axis_to_cc[axis]
            action_axes[controller].append((channel, distance))

        action = ActionFuture(MOVE_REL, action_axes)
        self._action_mgr.append_action(action)
        return action

    # TODO implement moveAbs fallback if not canAbs
    @isasync
    def moveAbs(self, pos):
        """
        Move the stage to the defined position in m for each axis given. This is an
        asynchronous method.
        pos dict(string-> float): name of the axis and new position in m
        returns (Future): object to control the move request
        """
        self._checkMoveAbs(pos)
        pos = self._applyInversionAbs(pos)
        # converts the request into one action (= a dict controller -> channels + distance)
        action_axes = collections.defaultdict(list)
        for axis, p in pos.items():
            controller, channel = self._axis_to_cc[axis]
            action_axes[controller].append((channel, p))

        action = ActionFuture(MOVE_ABS, action_axes)
        self._action_mgr.append_action(action)
        return action


    # TODO reference(self, axes)
    # start a move to reference. axes: set of str
    # returns a Future, done once the axis is referenced

    def stop(self, axes=None):
        """
        stops the motion on all axes
        Warning: this might stop the motion even of axes not managed (it stops
          all the axes of the related controllers).
        axes (set of str)
        """
        if self._action_mgr:
            self._action_mgr.cancel_all()

        axes = axes or self._axes
        # TODO: use the broadcast address to request a stop to all
        # controllers on the bus at the same time?

        # send stop to all controllers (including the ones not in action)
        controllers = set()
        for axis, (controller, channel) in self._axis_to_cc.items():
            if not axis in axes:
                continue
            if controller not in controllers:
                controller.stopMotion()
                controllers.add(controller)

        # wait all controllers are done moving
        for controller in controllers:
            controller.waitEndMotion()

    def terminate(self):
        if not hasattr(self, "_action_mgr"):
            # not even fully initialised
            return

        self.stop()

        if self._action_mgr:
            self._action_mgr.terminate()
            self._action_mgr = None

        controllers = set()
        for axis, (controller, channel) in self._axis_to_cc.items():
            if controller not in controllers:
                controller.terminate()

        self.accesser.terminate()

    def selfTest(self):
        """
        No move should be going one while doing a self-test
        """
        passed = True
        controllers = set([c for c, a in self._axis_to_cc.values()])
        for controller in controllers:
            logging.info("Testing controller %d", controller.address)
            passed &= controller.selfTest()

        return passed

    def _openPort(self, port, baudrate, _addresses):
        if port.startswith("/dev/") or port.startswith("COM"):
            ser = self._openSerialPort(port, baudrate, _addresses)
            return SerialBusAccesser(ser)
        else: # ip address
            if port == "autoip": # Search for IP (and hope there is only one result)
                ipmasters = self._scanIPMasters()
                if not ipmasters:
                    raise IOError("Failed to find any PI network master controller")
                host, ipport = ipmasters[0]
                logging.info("Will connect to %s:%d", host, ipport)
            else:
                # split the (IP) port, separated by a :
                if ":" in port:
                    host, ipport_str = port.split(":")
                    ipport = int(ipport_str)
                else:
                    host = port
                    ipport = 50000 # default
            
            sock = self._openIPSocket(host, ipport)
            return IPBusAccesser(sock)

    @classmethod
    def scan(cls, port=None):
        """
        port (string): name of the serial port. If None, all the serial ports are tried
        returns (list of 2-tuple): name, kwargs (port, axes(channel -> CL?)
        Note: it's obviously not advised to call this function if moves on the motors are ongoing
        """
        if port:
            ports = [port]
        else:
            # TODO: use serial.tools.list_ports.comports() (but only availabe in pySerial 2.6)
            if os.name == "nt":
                ports = ["COM" + str(n) for n in range (0, 8)]
            else:
                ports = glob.glob('/dev/ttyS?*') + glob.glob('/dev/ttyUSB?*')

        logging.info("Serial network scanning for PI-GCS controllers in progress...")
        found = []  # (list of 2-tuple): name, args (port, axes(channel -> CL?)
        axes_names = "xyzabcdefghijklmnopqrstuvw"
        for p in ports:
            try:
                # check all possible baud rates, in the most likely order
                for br in [38400, 9600, 19200, 115200]:
                    logging.debug("Trying port %s at baud rate %d", p, br)
                    ser = cls._openSerialPort(p, br)
                    controllers = Controller.scan(SerialBusAccesser(ser))
                    if controllers:
                        axis_num = 0
                        arg = {}
                        for add, axes in controllers.items():
                            for a, cl in axes.items():
                                arg[axes_names[axis_num]] = (add, a, cl)
                                axis_num += 1
                        found.append(("Actuator " + os.path.basename(p),
                                     {"port": p, "baudrate": br, "axes": arg}))
                        # it implies the baud rate was correct, and as it's impossible
                        # to have devices on different baud rate, so we are done
                        break
            except serial.SerialException:
                # not possible to use this port? next one!
                pass

        # Scan for controllers via each IP master controller
        ipmasters = cls._scanIPMasters()
        for ipadd in ipmasters:
            try:
                logging.debug("Scanning controllers on master %s:%d", ipadd[0], ipadd[1])
                sock = cls._openIPSocket(*ipadd)
                controllers = Controller.scan(IPBusAccesser(sock))
                if controllers:
                    axis_num = 0
                    arg = {}
                    for add, axes in controllers.items():
                        for a, cl in axes.items():
                            arg[axes_names[axis_num]] = (add, a, cl)
                            axis_num += 1
                    found.append(("Actuator IP",
                                 {"port": "%s:%d" % ipadd, "axes": arg}))
            except IOError:
                logging.info("Failed to scan on master %s:%d", ipadd[0], ipadd[1])
        
        return found
    
    @classmethod
    def _scanIPMasters(cls):
        """
        Scans the IP network for master controllers
        return (list of tuple of str, int): list of ip add and port of the master
          controllers found.
        """
        logging.info("Ethernet network scanning for PI-GCS controllers in progress...")
        found = set()  # (set of 2-tuple): ip address, ip port

        # Find all the broadcast addresses possible (one or more per network interfaces)
        # In the ideal world, we could just use '<broadcast>', but apprently if
        # there is not gateway to WAN, it will not work.
        bdc = []
        try:
            import netifaces
            for itf in netifaces.interfaces():
                try:
                    for addrinfo in netifaces.ifaddresses(itf)[socket.AF_INET]:
                        bdc.append(addrinfo["broadcast"])
                except KeyError:
                    pass # no INET or no "broadcast"
        except ImportError:
            bdc = ['<broadcast>']

        for bdcaddr in bdc:
            for port in [50000]: # TODO: the PI program tries on more ports
                # Special protocol by PI (reversed-engineered):
                # * Broadcast "PI" on a (known) port
                # * Listen for an answer
                try:
                    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
                    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                    s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
                    s.bind(('', 0))
                    logging.debug("Broadcasting on %s:%d", bdcaddr, port)
                    s.sendto('PI', (bdcaddr, port))
                    s.settimeout(1.0)  # It should take less than 1 s to answer

                    while True:
                        data, fulladdr = s.recvfrom(1024)
                        if not data:
                            break
                        # data should contain something like "PI C-863K016 SN 0 -- listening on port 50000 --"
                        if data.startswith("PI"):
                            found.add(fulladdr)
                        else:
                            logging.info("Received %s from %s", data.encode('string_escape'), fulladdr)
                except socket.timeout:
                    pass
                except socket.error:
                    logging.info("Couldn't broadcast on %s:%d", bdcaddr, port)
                except Exception:
                    logging.exception("Failed to broadcast on %s:%d", bdcaddr, port)

        return list(found)

    @staticmethod
    def _openSerialPort(port, baudrate=38400, _addresses=None):
        """
        Opens the given serial port the right way for the PI controllers.
        port (string): the name of the serial port (e.g., /dev/ttyUSB0)
        baudrate (int): baudrate to use, default is the recommended 38400
        _addresses (unused): only for testing (cf FakeBus)
        return (serial): the opened serial port
        """
        ser = serial.Serial(
            port=port,
            baudrate=baudrate,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=0.5 # s
        )

        return ser

    @staticmethod
    def _openIPSocket(host, port=50000):
        """
        Opens a socket connection to an PI master controller over IP.
        host (string): the IP address or host name of the master controller
        port (int): the (IP) port number
        return (socket): the opened socket connection
        """
        sock = socket.create_connection((host, port), timeout=5)
        sock.settimeout(0.5) # s
        return sock

class SerialBusAccesser(object):
    """
    Manages connections to the low-level bus
    """
    def __init__(self, serial):
        self.serial = serial
        # to acquire before sending anything on the serial port
        self.ser_access = threading.Lock()
        self.driverInfo = "serial driver: %s" % (driver.getSerialDriver(serial.port),)

    def terminate(self):
        self.serial.close()

    def sendOrderCommand(self, addr, com):
        """
        Send a command which does not expect any report back
        addr (None or 1<=int<=16): address of the controller. If None, no address
        is used (and it's typically controller 1 answering)
        com (string): command to send (including the \n if necessary)
        """
        assert(len(com) <= 100) # commands can be quite long (with floats)
        assert(1 <= addr <= 16 or addr == 254 or addr == 255)
        if addr is None:
            full_com = com
        else:
            full_com = "%d %s" % (addr, com)
        with self.ser_access:
            logging.debug("Sending: '%s'", full_com.encode('string_escape'))
            self.serial.write(full_com)
            # We don't flush, as it will be done anyway if an answer is needed

    def sendQueryCommand(self, addr, com):
        """
        Send a command and return its report (raw)
        addr (None or 1<=int<=16): address of the controller
        com (string): the command to send (without address prefix but with \n)
        return (string or list of strings): the report without prefix 
           (e.g.,"0 1") nor newline. 
           If answer is multiline: returns a list of each line
        Note: multiline answers seem to always begin with a \x00 character, but
         it's left as is.
        """
        assert(len(com) <= 100) # commands can be quite long (with floats)
        assert(1 <= addr <= 16 or addr == 254)
        if addr is None:
            full_com = com
        else:
            full_com = "%d %s" % (addr, com)
        with self.ser_access:
            logging.debug("Sending: '%s'", full_com.encode('string_escape'))
            self.serial.write(full_com)

            # ensure everything is received, before expecting an answer
            self.serial.flush()

            char = self.serial.read() # empty if timeout
            line = ""
            lines = []
            while char:
                if char == "\n":
                    if (len(line) > 0 and line[-1] == " " and  # multiline: "... \n"
                        not re.match(r"0 \d+ $", line)):  # excepted empty line "0 1 \n"
                        lines.append(line[:-1]) # don't include the space
                        line = ""
                    else:
                        # full end
                        lines.append(line)
                        break
                else:
                    # normal char
                    line += char
                char = self.serial.read()

        if not char:
            raise IOError("Controller %d timeout." % addr)

        assert len(lines) > 0

        logging.debug("Received: '%s'", "\n".join(lines).encode('string_escape'))
        if addr is None:
            prefix = ""
        else:
            prefix = "0 %d " % addr
        if not lines[0].startswith(prefix):
            raise IOError("Report prefix unexpected after '%s': '%s'." % (com, lines[0]))
        lines[0] = lines[0][len(prefix):]

        if len(lines) == 1:
            return lines[0]
        else:
            return lines

    def flushInput(self):
        """
        Ensure there is no more data queued to be read on the bus (=serial port)
        """
        with self.ser_access:
            # Flush buffer + give it some time to recover from whatever
            self.serial.flush()
            self.serial.flushInput()
            while self.serial.read():
                pass


class IPBusAccesser(object):
    """
    Manages connections to the low-level bus
    """
    def __init__(self, socket):
        self.socket = socket
        # to acquire before sending anything on the socket
        self.ser_access = threading.Lock()

        # recover the main controller from previous errors (just in case)
        err = self.sendQueryCommand(254, "ERR?\n")

        # Get the master controller version
        version = self.sendQueryCommand(254, "*IDN?\n")
        self.driverInfo = "%s" % (version.encode('string_escape'),)

    def terminate(self):
        self.socket.close()

    def sendOrderCommand(self, addr, com):
        """
        Send a command which does not expect any report back
        addr (None or 1<=int<=16): address of the controller. If None, no address
        is used (and it's typically controller 1 answering)
        com (string): command to send (including the \n if necessary)
        """
        assert(len(com) <= 100) # commands can be quite long (with floats)
        assert(1 <= addr <= 16 or addr == 254 or addr == 255)
        if addr is None:
            full_com = com
        else:
            full_com = "%d %s" % (addr, com)
        with self.ser_access:
            logging.debug("Sending: '%s'", full_com.encode('string_escape'))
            self.socket.sendall(full_com)

    def sendQueryCommand(self, addr, com):
        """
        Send a command and return its report (raw)
        addr (None or 1<=int<=16): address of the controller
        com (string): the command to send (without address prefix but with \n)
        return (string or list of strings): the report without prefix 
           (e.g.,"0 1") nor newline. 
           If answer is multiline: returns a list of each line
        """
        assert(len(com) <= 100) # commands can be quite long (with floats)
        assert(1 <= addr <= 16 or addr == 254)
        if addr is None:
            full_com = com
        else:
            full_com = "%d %s" % (addr, com)

        with self.ser_access:
            logging.debug("Sending: '%s'", full_com.encode('string_escape'))
            self.socket.sendall(full_com)

            # read the answer
            end_time = time.time() + 0.5
            ans = ""
            while True:
                try:
                    data = self.socket.recv(4096)
                except socket.timeout:
                    raise IOError("Controller %d timeout." % addr)
                # If the master is already accessed from somewhere else it will just
                # immediately answer an empty message
                if not data:
                    if time.time() > end_time:
                        raise IOError("Master controller not answering. It might be already connected with another client.")
                    time.sleep(0.01)
                    continue

                ans += data
                # does it look like we received the end of an answer?
                # To be really sure we'd need to wait until timeout, but that
                # would slow down a lot. Normally, if we've received one full
                # answer, there's 99% chance we've received everything.
                # An answer ends with \n (and not " \n", which indicates multi-
                # line).
                if (ans[-1] == "\n" and (
                     ans[-2:-1] != " " or  # multiline: "... \n"
                     re.match(r"0 \d+ $", ans))):  # excepted empty line "0 1 \n"
                    break

        logging.debug("Received: '%s'", ans.encode('string_escape'))

        # remove the prefix and last newline
        if addr is None:
            prefix = ""
        else:
            prefix = "0 %d " % addr
        if not ans.startswith(prefix):
            raise IOError("Report prefix unexpected after '%s': '%s'." % (com, ans))
        ans = ans[len(prefix):-1]

        # Interpret the answer
        lines = []
        for i, l in enumerate(ans.split("\n")):
            if l[-1] == " ": # remove the spaces indicating multi-line
                l = l[:-1]
            elif i != len(lines):
                logging.warning("Skipping previous answer from hardware %s",
                                "\n".join(lines + [l]).encode('string_escape'))
                lines = []
                continue
            lines.append(l)

        if len(lines) == 1:
            return lines[0]
        else:
            return lines

    def flushInput(self):
        """
        Ensure there is no more data queued to be read on the bus
        """
        with self.ser_access:
            try:
                while True:
                    data = self.socket.recv(4096)
            except socket.timeout:
                pass
            except Exception:
                logging.exception("Failed to flush correctly the socket")


# TODO: simplify, by using CancellableFuture, as in tmcm
class ActionManager(threading.Thread):
    """
    Thread running the requested actions (=moves)
    Provides a queue (deque) of actions (action_queue)
    For each action in the queue: performs and wait until the action is finished
    At the end of the action, call all the callbacks
    """
    def __init__(self, bus, name="PIGCS action manager"):
        threading.Thread.__init__(self, name=name)
        self.daemon = True # If the backend is gone, just die

        self.action_queue_cv = threading.Condition()
        self.action_queue = collections.deque()
        self.current_action = None
        self._bus = bus

    def run(self):
        while True:
            # Pick the next action
            with self.action_queue_cv:
                while not self.action_queue:
                    self.action_queue_cv.wait()
                self.current_action = self.action_queue.popleft()

            # Special action "None" == stop
            if self.current_action is None:
                return

            try:
                self.current_action._start_action()
            except futures.CancelledError:
                # cancelled in the mean time: skip the action
                continue

            while not self.current_action._wait_action(0.2):
                # regularly update position (5 Hz)
                self._bus._updatePosition()

            # Update position one last time
            # FIXME: should update position just _before_ calling the callbacks
            self._bus._updatePosition()

    def cancel_all(self):
        must_terminate = False
        with self.action_queue_cv:
            # cancel current action
            if self.current_action:
                self.current_action.cancel()

            # cancel every action in the queue
            while self.action_queue:
                action = self.action_queue.popleft()
                if action is None: # asking to terminate the thread
                    must_terminate = True
                    continue
                action.cancel()

        if must_terminate:
            self.append_action(None)

    def append_action(self, action):
        """
        appends an action in the doer's queue
        action (Action)
        """
        with self.action_queue_cv:
            self.action_queue.append(action)
            self.action_queue_cv.notify()

    def terminate(self):
        """
        Ask the action manager to terminate (once all the queued actions are done)
        """
        self.append_action(None)


MOVE_REL = "moveRel"
MOVE_ABS = "moveAbs"
MOVE_REF = "moveRef"

PENDING = 'PENDING'
RUNNING = 'RUNNING'
CANCELLED = 'CANCELLED'
FINISHED = 'FINISHED'

class ActionFuture(object):
    """
    Provides the interface for the clients to manipulate an (asynchronous) action 
    they requested.
    It follows http://docs.python.org/dev/library/concurrent.futures.html
    The result is always None, or raises an Exception.
    Internally, it has a reference to the action manager thread.
    """
    possible_types = [MOVE_REL, MOVE_ABS]

    # TODO handle exception in action
    def __init__(self, action_type, args):
        """
        type (str): name of the action (only supported so far is "moveRel"
        args (tuple): arguments to pass to the action
        """
        assert(action_type in self.possible_types)

        logging.debug("New action of type %s with arguments %s", action_type, args)
        self._type = action_type
        self._args = args
        self._expected_end = None # when it expects to finish (only during RUNNING)
        self._timeout = None # really too late to be running normally

        # acquire to modify the state, wait to wait for it to be done
        self._condition = threading.Condition()
        self._state = PENDING
        self._callbacks = []

    def _invoke_callbacks(self):
        # do not call with _condition! And ensure it's called only once
        for callback in self._callbacks:
            try:
                callback(self)
            except Exception:
                logging.exception('exception calling callback for %r', self)

    def cancel(self):
        with self._condition:
            if self._state == CANCELLED:
                return True
            elif self._state == FINISHED:
                return False
            elif self._state == RUNNING:
                self._stop_action()
                # go through, like for state == PENDING

            self._state = CANCELLED
            self._condition.notify_all()

        self._invoke_callbacks()
        return True

    def cancelled(self):
        with self._condition:
            return self._state == CANCELLED

    def running(self):
        with self._condition:
            return self._state == RUNNING

    def done(self):
        with self._condition:
            return self._state in [CANCELLED, FINISHED]

    def result(self, timeout=None):
        with self._condition:
            if self._state == CANCELLED:
                raise futures.CancelledError()
            elif self._state == FINISHED:
                return None

            self._condition.wait(timeout)

            if self._state == CANCELLED:
                raise futures.CancelledError()
            elif self._state == FINISHED:
                return None
            else:
                raise futures.TimeoutError()

    def exception(self, timeout=None):
        """
        return None or return what result raises
        """
        try:
            return self.result(timeout)
        except (futures.TimeoutError, futures.CancelledError) as exp:
            raise exp
        except Exception as exp:
            return exp

    def add_done_callback(self, fn):
        with self._condition:
            if self._state not in [CANCELLED, FINISHED]:
                self._callbacks.append(fn)
                return
        fn(self)

    def _start_action(self):
        """
        Start the physical action, and immediately return. It also set the 
        state to RUNNING.
        Note: to be called without the lock (._condition) acquired.
        """
        with self._condition:
            if self._state == CANCELLED:
                raise futures.CancelledError()

            # Do the action
            if self._type == MOVE_REL:
                duration = self._moveRel(self._args)
            elif self._type == MOVE_ABS:
                duration = self._moveAbs(self._args)
            else:
                raise Exception("Unknown action %s" % self._type)

            self._state = RUNNING
            duration = min(duration, 60) # => wait maximum 2 min
            self._expected_end = time.time() + duration
            self._timeout = self._expected_end + duration + 1 # 2 *duration + 1s

    def _wait_action(self, timeout=None):
        """
        Wait for the action to finish normally. If the action finishes normally
        it's also in charge of calling all the callbacks.
        timeout (float or None): maximum time to wait
        return True if the action is finished/cancel, and False if it's still
         running (so more waiting is needed) 
        Note: to be called without the lock (._condition) acquired.
        """
        # create a dict of controllers => channels
        controllers = {}
        for controller, moves in self._args.items():
            channels = [c for c, d in moves]
            controllers[controller] = channels

        with self._condition:
            assert(self._expected_end is not None)
            # if it has been cancelled in the mean time
            if self._state != RUNNING:
                return True

            if timeout is None:
                end_wait = self._timeout + 1 # => will never be triggered
            else:
                end_wait = time.time() + timeout

            # it's over when either all axes are finished moving, it's too late,
            # or the action was cancelled
            logging.debug("Waiting %f s for the move to finish", self._expected_end - time.time())
            while self._state == RUNNING:
                now = time.time()
                if now > self._timeout:
                    logging.debug("Giving up waiting for move %g s late",
                                  now - self._expected_end)
                    break # move took too much time
                if now > end_wait:
                    return False # could wait more but the caller is not interested
                duration = (self._expected_end - now) / 2
                duration = max(0.01, duration)
                self._condition.wait(duration)
                if not self._isMoving(controllers): # we are done!
                    self._checkError(controllers)
                    break

            # if cancelled, we don't update state
            if self._state != RUNNING:
                return True

            self._state = FINISHED
            self._condition.notify_all()

        self._invoke_callbacks()
        return True

    def _stop_action(self):
        """
        Stop the action. Do not call directly, call cancel()
        Note: to be called with the lock (._condition) acquired.
        """
        # The only two possible actions are stopped the same way

        # create a dict of controllers => channels
        controllers = {}
        for controller, moves in self._args.items():
            channels = [c for c, d in moves]
            controllers[controller] = channels

        self._stopMotion(controllers)

    def _isMoving(self, axes):
        """
        axes (dict: Controller -> list (int)): controller to channel which must be check for move
        """
        moving = False
        for controller, channels in axes.items():
            if len(channels) == 0:
                logging.warning("Asked to check move on a controller without any axis")
            else:
                moving |= controller.isMoving(set(channels))
        return moving

    def _checkError(self, axes):
        """
        axes (dict: Controller -> list (int)): controller to channel which must be check for error
        return nothing
        raise PIGCSError if an error on a controller happened
        """
        for controller in axes.keys():
            err = controller.GetErrorNum()
            if err:
                raise PIGCSError(err)

    def _stopMotion(self, axes):
        """
        axes (dict: Controller -> list (int)): controller to channel which must be stopped
        """
        for controller in axes:
            # it can only stop all axes (that's the point anyway)
            controller.stopMotion()

    def _moveRel(self, axes):
        """
        axes (dict: Controller -> list (tuple(int, double)): 
            controller to list of channel/distance to move (m)
        returns (float): approximate time in s it will take (optimistic)
        """
        max_duration = 0 #s
        for controller, channels in axes.items():
            for channel, distance in channels:
                actual_dist = controller.moveRel(channel, distance)
                duration = abs(actual_dist) / controller.getSpeed(channel)
                max_duration = max(max_duration, duration)

        return max_duration

    def _moveAbs(self, axes):
        """
        axes (dict: Controller -> list (tuple(int, double)): 
            controller to list of channel/position to move (m)
        returns (float): approximate time in s it will take (optimistic)
        """
        max_duration = 0 #s
        for controller, channels in axes.items():
            for channel, pos in channels:
                actual_dist = controller.moveAbs(channel, pos)
                duration = abs(actual_dist) / controller.getSpeed(channel)
                max_duration = max(max_duration, duration)

        return max_duration


# All the classes below are for the simulation of the hardware

class SimulatedError(Exception):
    """
    Special exception class to simulate error in the controller
    """
    pass

class E861Simulator(object):
    """
    Simulates a GCS controller (+ serial port at 38400). Only used for testing.
    1 axis, open-loop only, very limited behaviour
    Same interface as the serial port
    """
    _idn = "(c)2013 Delmic Fake Physik Instrumente(PI) Karlsruhe, E-861 Version 7.2.0"
    _csv = "2.0"
    def __init__(self, port, baudrate=9600, timeout=0, address=1, 
                 closedloop=False, *args, **kwargs):
        """
        parameters are the same as a serial port
        address (1<=int<=16): the address of the controller
        closedloop (bool): whether it simulates a closed-loop actuator or not
        """
        self.port = port
        self._address = address
        self._has_encoder = closedloop
        # we don't care about the actual parameters but timeout
        self.timeout = timeout

        self._init_mem()
        
        self._end_move = 0 # time the last requested move is over

        # only used in closed-loop
        # If move is over:
        #   position == target
        # else:
        #   position = original position
        #   target = requested position
        #   current position = weigthed average (according to time)
        self._position = 0.012 # m
        self._target = self._position # m 
        self._start_move = 0

        self._output_buf = "" # what the commands sends back to the "host computer"
        self._input_buf = "" # what we receive from the "host computer"

        # special trick to only answer if baudrate is correct
        if baudrate != 38400:
            logging.debug("Baudrate incompatible: %d", baudrate)
            self.write = (lambda s: "")

    def _init_mem(self):
        # internal values to simulate the device
        # Note: the type is used to know how it should be decoded, so it's
        # important to differentiate between float and int.
        # Parameter table: address -> value
        self._parameters = {0x14: 1 if self._has_encoder else 0, # 0 = no ref switch, 1 = ref switch
                            0x32: 0 if self._has_encoder else 1, # 0 = limit switches, 1 = no limit switches
                            0x3c: "DEFAULT-FAKE", # stage name
                            0x15: 0.025, # TMX (in m)
                            0x30: 0.0, # TMN (in m)
                            0x16: 0.012, # value at ref pos
                            0x49: 10.0, # VEL
                            0x0B: 3.2, # ACC
                            0x0C: 0.9, # DEC
                            0x0A: 50.0, # max vel
                            0x4A: 5.0, # max acc
                            0x4B: 5.0, # max dec
                            0x0E: 10000000, # unit num (note: normal default is 10000)
                            0x0F: 1,       # unit denum
                            0x7000003: 10.0, # SSA
                            0x7000201: 3.2, # OVL
                            0x7000202: 0.9, # OAC
                            0x7000204: 15.3, # max step/s
                            0x7000205: 1.2, # max step/s²
                            0x7000206: 0.9, # ODC
                            }
        self._servo = 0 # servo state
        self._ready = True # is ready?
        self._referenced = 0
        self._ref_mode = 1
        self._errno = 0 # last error set

    _re_command = ".*?[\n\x04\x05\x07\x08\x18\x24]"
    def write(self, data):
        self._input_buf += data
        # process each commands separated by a "\n" or is short command
        while len(self._input_buf) > 0:
            m = re.match(self._re_command, self._input_buf)
            if not m:
                return # no more full command available
            c = m.group(0)
            self._processCommand(c)
            self._input_buf = self._input_buf[m.end(0):] # all the left over

    def read(self, size=1):
        # simulate timeout
        end_time = time.time() + self.timeout

        # FIXME: to be correct, we'd need to take a lock
        ret = self._output_buf[:size]
        self._output_buf = self._output_buf[len(ret):]

        while len(ret) < size:
            time.sleep(0.01)
            left = size - len(ret)
            ret += self._output_buf[:left]
            self._output_buf = self._output_buf[len(ret):]
            if self.timeout and time.time() > end_time:
                break

        return ret

    def close(self):
        # using read or write will fail after that
        del self._output_buf
        del self._input_buf

    def _get_cur_pos_cl(self):
        """
        Computes the current position, in closed loop mode
        """
        if self._servo != 1:
            raise SimulatedError(2)
        now = time.time()
        if now > self._end_move:
            self._position = self._target
            return self._position
        else:
            completion = (now - self._start_move) / (self._end_move - self._start_move)
            cur_pos = self._position + (self._target - self._position) * completion
            return cur_pos
    
    # TODO: some commands are read-only
    # Command name -> parameter number
    _com_to_param = {# "LIM": 0x32, # LIM actually report the opposite of 0x32
                     "TRS": 0x14,
                     "CTS": 0x3c,
                     "TMN": 0x30,
                     "TMX": 0x15,
                     "VEL": 0x49,
                     "ACC": 0x0B,
                     "DEC": 0x0C,
                     "OVL": 0x7000201,
                     "OAC": 0x7000202,
                     "ODC": 0x7000206,
                     "SSA": 0x7000003,
                     }
    _re_addr_com = r"((?P<addr>\d+) (0 )?)?(?P<com>.*)"
    def _processCommand(self, com):
        """
        process the command, and put the result in the output buffer
        com (str): command
        """
        logging.debug("Fake controller %d processing command '%s'",
                       self._address, com.encode('string_escape'))
        out = None # None means error while decoding command

        # command can start with a prefix like "5 0 " or "5 "
        m = re.match(self._re_addr_com, com)
        assert m # anything left over should be in com
        if m.group("addr"):
            addr = int(m.group("addr"))
            prefix = "0 %d " % addr
        else:
            addr = 1 # default is address == 1
            prefix = ""

        if addr != self._address and addr != 255: # message is for us?
            logging.debug("Controller %d skipping message for %d",
                          self._address, addr)
            return

        com = m.group("com") # also removes the \n at the end if it's there
        # split into arguments separated by spaces (not including empty strings)
        args = filter(bool, com.split(" "))
        logging.debug("Command decoded: %s", args)

        if self._errno:
            # if errno is not null, most commands don't work any more
            if com not in ["*IDN?", "RBT", "ERR?", "CSV?"]:
                logging.debug("received command %s while errno = %d",
                              com.encode('string_escape'), self._errno)
                return

        # TODO: to support more commands, we should have a table, with name of
        # the command + type of arguments (+ number of optional args)
        try:
            if com == "*IDN?": # identification
                out = self._idn
            elif com == "CSV?": # command set version
                out = self._csv
            elif com == "ERR?": # last error number
                out = "%d" % self._errno
                self._errno = 0 # reset error number
            elif com == "RBT": # reboot
                self._init_mem()
                time.sleep(0.1)
            elif com == "\x04": # Query Status Register Value
                # return hexadecimal bitmap of moving axes
                # TODO: to check, much more info returned
                val = 0
                if time.time() < self._end_move:
                    val |= 0x400  #  first axis moving
                out = "0x%x" % val
            elif com == "\x05": # Request Motion Status
                # return hexadecimal bitmap of moving axes
                if time.time() > self._end_move:
                    val = 0
                else:
                    val = 1 # first axis moving
                out = "%x" % val
            elif com == "\x07": # Request Controller Ready Status
                if self._ready: # TODO: when is it not ready??
                    out = "\xb1"
                else:
                    out = "\xb2"
            elif com == "\x18" or com == "STP": # Stop immediately
                self._end_move = 0
                self._errno = 10 # PI_CNTR_STOP
            elif args[0].startswith("HLT"): # halt motion with deceleration: axis (optional)
                self._end_move = 0
            elif args[0][:3] in self._com_to_param:
                param = self._com_to_param[args[0][:3]]
                logging.debug("Converting command %s to param %d", args[0], param)
                axis = int(args[1])
                if axis != 1:
                    raise SimulatedError(15)
                if args[0][3:4] == "?" and len(args) == 2: # query
                    out = "%s=%s" % (args[1], self._parameters[param])
                elif len(args[0]) == 3 and len(args) == 3: # set
                    # convert according to the current type of the parameter
                    typeval = type(self._parameters[param])
                    self._parameters[param] = typeval(args[2])
                else:
                    raise SimulatedError(15)
            elif args[0] == "SPA?" and len(args) == 3: # GetParameter: axis, address
                # TODO: when no arguments -> list all parameters
                axis, addr = int(args[1]), int(args[2])
                if axis != 1:
                    raise SimulatedError(15)
                try:
                    out = "%d=%s" % (addr, self._parameters[addr])
                except KeyError:
                    logging.debug("Unknown parameter %d", addr)
                    raise SimulatedError(56)
            elif args[0] == "SPA" and len(args) == 4: # SetParameter: axis, address, value
                axis, addr = int(args[1]), int(args[2])
                if axis != 1:
                    raise SimulatedError(15)
                if addr in [0x0E, 0x0F] and self._parameters[addr] != int(args[3]):
                    # TODO: have a list of parameters to update
                    raise NotImplementedError("Simulator cannot change unit")
                try:
                    typeval = type(self._parameters[addr])
                    self._parameters[addr] = typeval(args[3])
                except KeyError:
                    logging.debug("Unknown parameter %d", addr)
                    raise SimulatedError(56)
            elif args[0] == "LIM?" and len(args) == 2: # Get Limit Switches
                axis = int(args[1])
                if axis == 1:
                    # opposite of param 0x32
                    out = "%s=%s" % (args[1], 1 - self._parameters[0x32])
                else:
                    self._errno = 15
            elif args[0] == "SVO" and len(args) == 3: # Set Servo State
                axis, state = int(args[1]), int(args[2])
                if axis == 1:
                    self._servo = state
                else:
                    self._errno = 15
            elif args[0] == "RON" and len(args) == 3: # Set Reference mode
                axis, state = int(args[1]), int(args[2])
                if axis == 1:
                    self._ref_mode = state
                else:
                    self._errno = 15
            elif args[0] == "OSM" and len(args) == 3: #Open-Loop Step Moving
                axis, steps = int(args[1]), float(args[2])
                if axis != 1:
                    raise SimulatedError(15)
                speed = self._parameters[self._com_to_param["OVL"]]
                duration = abs(steps) / speed
                logging.debug("Simulating a move of %f s", duration)
                self._end_move = time.time() + duration # current move stopped
            elif args[0] == "MOV" and len(args) == 3: # Closed-Loop absolute move
                axis, pos = int(args[1]), float(args[2])
                if axis != 1:
                    raise SimulatedError(15)
                if self._ref_mode and not self._referenced:
                    raise SimulatedError(8)
                speed = self._parameters[self._com_to_param["VEL"]]
                cur_pos = self._get_cur_pos_cl()
                distance = cur_pos - pos
                duration = abs(distance) / speed
                logging.debug("Simulating a move of %f s", duration)
                self._start_move = time.time()
                self._end_move = self._start_move + duration
                self._position = cur_pos
                self._target = pos
            elif args[0] == "MVR" and len(args) == 3: # Closed-Loop relative move
                axis, distance = int(args[1]), float(args[2])
                if axis != 1:
                    raise SimulatedError(15)
                if self._ref_mode and not self._referenced:
                    raise SimulatedError(8)
                speed = self._parameters[self._com_to_param["VEL"]]
                duration = abs(distance) / speed
                logging.debug("Simulating a move of %f s", duration)
                cur_pos = self._get_cur_pos_cl()
                self._start_move = time.time()
                self._end_move = self._start_move + duration
                self._position = cur_pos
                self._target = cur_pos + distance
            elif args[0] == "POS?" and len(args) == 2: # Closed-Loop position query
                axis = int(args[1])
                if axis != 1:
                    raise SimulatedError(15)
                out = "%s=%s" % (args[1], self._get_cur_pos_cl())
            elif args[0] == "ONT?" and len(args) == 2: # on target
                axis = int(args[1])
                if axis != 1:
                    raise SimulatedError(15)
                ont = time.time() > self._end_move
                out = "%s=%d" % (args[1], 1 if ont else 0)
            elif args[0] == "FRF?" and len(args) == 2: # is referenced?
                axis = int(args[1])
                if axis != 1:
                    raise SimulatedError(15)
                out = "%s=%d" % (args[1], self._referenced)
            elif args[0] == "FRF" and len(args) == 2: # reference to ref switch
                axis = int(args[1])
                if axis != 1:
                    raise SimulatedError(15)
                self._referenced = 1
                self._end_move = 0
                self._position = self._parameters[0x16] # value at reference
            elif args[0] == "SAI?" and len(args) <= 2: # List Of Current Axis Identifiers
                # Can be followed by "ALL", but for us, it's the same
                out = "1"
            elif com == "HLP?":
                # The important part is " \n" at the end of each line
                out = ("\x00The following commands are available: \n" +
                       "#4 request status register \n" +
                       "HLP list the available commands \n" +
                       "ERR? get error number \n" +
                       "VEL {<AxisId> <Velocity>} set closed-loop velocity \n" +
                       "end of help"
                       )
            elif com == "HPA?":
                out = ("\x00The following parameters are valid: \n" +
                       "0x1=\t0\t1\tINT\tmotorcontroller\tP term 1 \n" +
                       "0x32=\t0\t1\tINT\tmotorcontroller\thas limit\t(0=limitswitchs 1=no limitswitchs) \n" +
                       "0x3C=\t0\t1\tCHAR\tmotorcontroller\tStagename \n" +
                       "0x7000000=\t0\t1\tFLOAT\tmotorcontroller\ttravel range minimum \n" +
                       "end of help"
                       )
            else:
                logging.debug("Unknown command '%s'", com)
                self._errno = 1
        except SimulatedError as ex:
            logging.debug("Error detected while processing command '%s'", com)
            self._errno = ex.args[0]
        except Exception:
            logging.debug("Failed to process command '%s'", com)
            self._errno = 1

        # add the response header
        if out is None:
            #logging.debug("Fake controller %d doesn't respond", self._address)
            pass
        else:
            out = "%s%s\n" % (prefix, out)
            logging.debug("Fake controller %d responding '%s'", self._address,
                          out.encode('string_escape'))
            self._output_buf += out

class DaisyChainSimulator(object):
    """
    Simulated serial port that can simulate daisy chain on the controllers
    Same interface as the serial port + list of (fake) serial ports to connect 
    """
    def __init__(self, port, timeout=0, *args, **kwargs):
        """
        subports (list of open ports): the ports to receive the data
        """
        self.port = port
        self.timeout = timeout
        self._subports = kwargs["subports"]
        self._output_buf = "" # TODO: probably cleaner to user lock to access it

        # TODO: for each port, put a thread listening on the read and push to output
        self._is_terminated = False
        for p in self._subports:
            t = threading.Thread(target=self._thread_read_serial, args=(p,))
            t.daemon = True
            t.start()

    def flush(self):
        return

    def flushInput(self):
        return

    def write(self, data):
        # just duplicate
        for p in self._subports:
            p.write(data)

    def read(self, size=1):
        # simulate timeout
        end_time = time.time() + self.timeout

        ret = self._output_buf[:size]
        self._output_buf = self._output_buf[len(ret):]

        while len(ret) < size:
            time.sleep(0.01)
            left = size - len(ret)
            ret += self._output_buf[:left]
            self._output_buf = self._output_buf[len(ret):]
            if self.timeout and time.time() > end_time:
                break

        return ret

    def _thread_read_serial(self, ser):
        """
        Push the output of the given serial port into our output
        """
        try:
            while not self._is_terminated:
                c = ser.read(1)
                if len(c) == 0:
                    time.sleep(0.01)
                else:
                    self._output_buf += c
        except Exception:
            logging.exception("Fake daisy chain thread received an exception")

    def close(self):
        self._is_terminated = True
        # using read or write will fail after that
        del self._output_buf
        del self._subports

class FakeBus(Bus):
    """
    Same as the normal Bus, but connects to simulated controllers
    """
    def __init__(self, name, role, port, axes, baudrate=38400, **kwargs):
        # compute the addresses from the axes declared
        addresses = dict([(d[0], d[2]) for d in axes.values()])
        Bus.__init__(self, name, role, port, axes, baudrate=baudrate,
                     _addresses=addresses, **kwargs)

    @classmethod
    def scan(cls, port=None):
        # force only one port
        return super(FakeBus, cls).scan(port="/dev/fake")

    @staticmethod
    def _openSerialPort(port, baudrate=38400, _addresses=None):
        """
        Opens a fake serial port
        port (string): the name of the serial port (e.g., /dev/ttyUSB0)
        _addresses (dict of int -> bool, or None): list of each address that should have 
         a simulated controller created and wheter it is closed-loop or not.
         Default to {1: False, 2:False} (used for scan).
        return (serial): the opened serial port
        """
        _addresses = _addresses or {1: False, 2: False}
        simulators = []
        for addr, cl in _addresses.items():
            sim = E861Simulator(
                    port=port,
                    baudrate=baudrate,
                    bytesize=serial.EIGHTBITS,
                    parity=serial.PARITY_NONE,
                    stopbits=serial.STOPBITS_ONE,
                    timeout=0.5, #s
                    address=addr,
                    closedloop=cl,
                   )
            simulators.append(sim)

        # link all of them in daisy chain
        ser = DaisyChainSimulator(
                port=port,
                baudrate=baudrate,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE,
                timeout=0.5, #s
                subports=simulators,
              )

        return ser

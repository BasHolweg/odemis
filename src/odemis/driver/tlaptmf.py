# -*- coding: utf-8 -*-
'''
Created on 25 Mar 2014

@author: Éric Piel

Copyright © 2014 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.
'''

# Driver for the Thorlabs "MFF10X" motorised filter flipper mounts. It uses the APT
# protocol (over serial/USB).
# Most of the protocol is documented in APT_Communications_Protocol_Rev_9.pdf
# http://www.thorlabs.co.uk/software/apt/APT_Communications_Protocol_Rev_9.pdf
# (provided by Thorlabs on request). This protocol allows to manage a very wide
# variety of devices.

# For now, we have a simple implementation of APT directly here, but if more
# devices are to be supported, it should be move to a APT library layer.
# cf pyAPT: https://github.com/freespace/pyAPT
# The typical way distinguish Thorlabs devices is to indicate the serial number
# of the device (which is clearly physically written on it too). This can be
# then easily compared with the USB attribute cf /sys/bus/usb/devices/*/serial

from __future__ import division

from Pyro4.core import isasync
import glob
import logging
import math
from odemis import model
import odemis
from odemis.model._futures import CancellableThreadPoolExecutor
from odemis.util import driver
import os
import serial
import struct
import sys
import threading
import time


# Classes for defining the messages
class APTMessage(object):
    # TODO: also indicates whether the command expect p1, p2, or the length of the data
    def __init__(self, mid):
        """
        mid (int): Message ID
        """
        assert 1 <= mid <= 0xffff
        self.id = mid

class APTSet(APTMessage):
    """
    Represent a command message which does not expect a response
    """
    pass

class APTReq(APTMessage):
    """
    Represent a request message, which expects a response
    """
    def __init__(self, mid, rid):
        """
        mid (int): Message ID
        rid (int): Message ID of the response
        """
        assert 1 <= rid <= 0xffff
        APTMessage.__init__(self, mid)
        self.rid = rid

# Messages
MOD_IDENTIFY = APTSet(0x0223)
HW_REQ_INFO = APTReq(0x0005, 0x0006)
HW_START_UPDATEMSGS = APTSet(0x0011)
HW_STOP_UPDATEMSGS = APTSet(0x0012)
HW_NO_FLASH_PROGRAMMING = APTSet(0x0018) # Or 0x0017???
MOT_MOVE_JOG = APTSet(0x046a)
MOT_MOVE_STOP = APTSet(0x0465)
MOT_SUSPEND_ENDOFMOVEMSGS = APTSet(0x046b)
MOT_RESUME_ENDOFMOVEMSGS = APTSet(0x046c)
MOT_REQ_STATUSUPDATE = APTReq(0x0480, 0x0481)
MOT_REQ_DCSTATUSUPDATE = APTReq(0x0490, 0x0491)
MOT_ACK_DCSTATUSUPDATE = APTSet(0x0492)
MOT_SET_AVMODES = APTSet(0x04b3)
MOT_REQ_POWERPARAMS = APTReq(0x0427, 0x0428)
MOT_REQ_JOGPARAMS = APTReq(0x0417, 0x0418)
# FIXME: these ones are event messages from the device
MOT_MOVE_COMPLETED = APTSet(0x0464)
MOT_MOVE_STOPPED = APTSet(0x0466)
# TODO: Use this to change the "transit time" (= speed)
MOT_SET_MFF_OPERPARAMS = APTSet(0x0510)
MOT_REQ_MFF_OPERPARAMS = APTReq(0x0511, 0x0512)

# Status flags (for MOT_REQ_*STATUSUPDATE)
# There are more, but we don't use them for now (cf p.90)
STA_FWD_HLS = 0x0001
STA_RVS_HLS = 0x0002
STA_FWD_MOT = 0x0010
STA_RVS_MOT = 0x0020
STA_FWD_JOG = 0x0040
STA_RVS_JOG = 0x0080
STA_CHA_ENB = 0x80000000

STA_IN_MOTION = (STA_FWD_MOT | STA_RVS_MOT | STA_FWD_JOG | STA_RVS_JOG)

# All MFFxxx have serial number starting with 37
SN_PREFIX_MFF = "37"

POS_UP = 0
POS_DOWN = math.radians(90)

class MFF(model.Actuator):
    """
    Represents one Thorlabs Motorized Filter Flipper (ie: MFF101 or MFF102)
    """
    def __init__(self, name, role, sn=None, port=None, axis="rz", inverted=None, **kwargs):
        """
        sn (str): serial number (recommended)
        port (str): port name (only if sn is not specified)
        axis (str): name of the axis
        inverted (set of str): names of the axes which are inverted (IOW, either
         empty or the name of the axis) 
        """
        if (sn is None and port is None) or (sn is not None and port is not None):
            raise ValueError("sn or port argument must be specified (but not both)")
        if sn is not None:
            if not sn.startswith(SN_PREFIX_MFF) or len(sn) != 8:
                logging.warning("Serial number '%s' is unexpected for a MFF "
                                "device (should be 8 digits starting with %s).",
                                sn, SN_PREFIX_MFF)
            self._port = self._getSerialPort(sn)
        else:
            self._port = port

        self._serial = self._openSerialPort(self._port)
        self._ser_access = threading.Lock()

        # Ensure we don't receive anything
        self.SendMessage(HW_STOP_UPDATEMSGS)
        self._serial.flushInput()

        # Documentation says it should be done first, though it doesn't seem
        # required
        self.SendMessage(HW_NO_FLASH_PROGRAMMING)

        # will take care of executing axis move asynchronously
        self._executor = CancellableThreadPoolExecutor(max_workers=1) # one task at a time

        # TODO: have the standard inverted Actuator functions work on enumerated
        # use a different format than the standard Actuator
        if inverted and axis in inverted:
            self._pos_to_jog = {POS_UP: 2,
                                POS_DOWN: 1}
            self._status_to_pos = {STA_RVS_HLS: POS_UP,
                                   STA_FWD_HLS: POS_DOWN}
        else:
            self._pos_to_jog = {POS_UP: 1,
                                POS_DOWN: 2}
            self._status_to_pos = {STA_FWD_HLS: POS_UP,
                                   STA_RVS_HLS: POS_DOWN}

        # TODO: add support for speed
        axes = {axis: model.Axis(unit="rad",
                                 choices=set(self._pos_to_jog.keys()))
                }
        model.Actuator.__init__(self, name, role, axes=axes, **kwargs)

        driver_name = driver.getSerialDriver(self._port)
        self._swVersion = "%s (serial driver: %s)" % (odemis.__version__, driver_name)
        snd, modl, typ, fmv, notes, hwv, state, nc = self.GetInfo()
        self._hwVersion = "%s v%d (firmware %s)" % (modl, hwv, fmv)

        self.position = model.VigilantAttribute({}, readonly=True)
        self._updatePosition()

        # It'd be nice to know when a move is over, but it seems the MFF10x
        # never report ends of move.
        # self.SendMessage(MOT_RESUME_ENDOFMOVEMSGS)

        # If we need constant status updates, then, we'll need to answer them
        # with MOT_ACK_DCSTATUSUPDATE at least once per second.
        # For now we don't track the current device status, so it's easy.
        # When requesting update messages, messages are sent at ~10Hz, even if
        # no change has happened.
        # self.SendMessage(HW_START_UPDATEMSGS) # Causes a lot of messages

        # We should make sure that the led is always off, but apparently, it's
        # always off without doing anything (cf MOT_SET_AVMODES)

    def terminate(self):
        if self._executor:
            self.stop()
            self._executor.shutdown()
            self._executor = None

        with self._ser_access:
            if self._serial:
                self._serial.close()
                self._serial = None

    def SendMessage(self, msg, dest=0x50, src=1, p1=None, p2=None, data=None):
        """
        Send a message to a device and possibility wait for its response
        msg (APTSet or APTReq): the message definition
        dest (0<int): the destination ID (always 0x50 if directly over USB)
        p1 (None or 0<=int<=255): param1 (passed as byte2)
        p2 (None or 0<=int<=255): param2 (passed as byte3)
        data (None or bytes): data to be send further. Cannot be mixed with p1
          and p2
        return (None or bytes): the content of the response or None if it was
          an APTSet message
        raise:
           IOError: if failed to send or receive message
        """
        assert 0 <= dest < 0x80

        # create the message
        if data is None: # short message
            p1 = p1 or 0
            p2 = p2 or 0
            com = struct.pack("<HBBBB", msg.id, p1, p2, dest, src)
        else: # long message
            com = struct.pack("<HHBB", msg.id, len(data), dest | 0x80, src) + data

        logging.debug("Sending: '%s'", ", ".join("%02X" % ord(c) for c in com))
        with self._ser_access:
            self._serial.write(com)

            if isinstance(msg, APTReq):  # read the response
                # ensure everything is sent, before expecting an answer
                self._serial.flush()

                # Read until end of answer
                while True:
                    rid, res = self._ReadMessage()
                    if rid == msg.rid:
                        return res
                    logging.debug("Skipping unexpected message %X", rid)

    # Note: unused
    def WaitMessage(self, msg, timeout=None):
        """
        Wait until a specified message is received
        msg (APTMessage)
        timeout (float or None): maximum amount of time to wait
        return (bytes): the 2 params or the data contained in the message
        raise:
            IOError: if timeout happened
        """
        start = time.time()
        # Read until end of answer
        with self._ser_access:
            while True:
                if timeout is not None:
                    left = time.time() - start + timeout
                    if left <= 0:
                        raise IOError("No message %d received in time" % msg.id)
                else:
                    left = None

                mid, res = self._ReadMessage(timeout=left)
                if mid == msg.id:
                    return res
                # TODO: instead of discarding the message, it could go into a
                # queue, to be handled later
                logging.debug("Skipping unexpected message %X", mid)

    def _ReadMessage(self, timeout=None):
        """
        Reads the next message
        timeout (0 < float): maximum time to wait for the message
        return:
             mid (int): message ID
             data (bytes): bytes 3&4 or the data of the message
        raise:
           IOError: if failed to send or receive message
        """
        old_timeout = self._serial.timeout
        if timeout is not None:
            # Should be only for the first byte, but doing it for the first 6
            # should rarely matter
            self._serial.timeout = timeout
        try:
            # read the first (required) 6 bytes
            msg = b""
            for i in range(6):
                char = self._serial.read() # empty if timeout
                if not char:
                    raise IOError("Controller timeout, after receiving %s" % msg)

                msg += char
        finally:
            self._serial.timeout = old_timeout

        mid = struct.unpack("<H", msg[0:2])[0]
        if not (ord(msg[4]) & 0x80): # short message
            logging.debug("Received: '%s'", ", ".join("%02X" % ord(c) for c in msg))
            return mid, msg[2:4]

        # long message
        length = struct.unpack("<H", msg[2:4])[0]
        for i in range(length):
            char = self._serial.read() # empty if timeout
            if not char:
                raise IOError("Controller timeout, after receiving %s" % msg)

            msg += char

        logging.debug("Received: '%s'", ", ".join("%02X" % ord(c) for c in msg))
        return mid, msg[6:]

    # Low level functions
    def GetInfo(self):
        """
        returns:
            serial number (int)
            model number (str)
            type (int)
            firmware version (str)
            notes (str)
            hardware version (int)
            hardware state (int)
            number of channels (int)
        """
        res = self.SendMessage(HW_REQ_INFO)
        # Expects 0x54 bytes
        values = struct.unpack('<I8sHI48s12xHHH', res)
        sn, modl, typ, fmv, notes, hwv, state, nc = values

        # remove trailing 0's
        modl = modl.rstrip("\x00")
        notes = notes.rstrip("\x00")

        # Convert firmware version to a string
        fmvs = "%d.%d.%d" % ((fmv & 0xff0000) >> 16,
                             (fmv & 0xff00) >> 8,
                             fmv & 0xff)

        return sn, modl, typ, fmvs, notes, hwv, state, nc

    def MoveJog(self, pos):
        """
        Move the position. Note: this is asynchronous.
        pos (int): 1 or 2
        """
        assert pos in [1, 2]
        # p1 is chan ident, always 1
        self.SendMessage(MOT_MOVE_JOG, p1=1, p2=pos)

    def GetStatus(self):
        """
        return:
            pos (int): position count
            status (int): status, as a flag of STA_*
        """
        res = self.SendMessage(MOT_REQ_STATUSUPDATE)
        # expect 14 bytes
        c, pos, enccount, status = struct.unpack('<HiiI', res)

        return pos, status
    

    # high-level methods (interface)
    def _updatePosition(self):
        """
        update the position VA
        """
        _, status = self.GetStatus()
        pos = {}
        for axis in self.axes: # axes contains precisely one axis
            # status' flags should never be present simultaneously
            for f, p in self._status_to_pos.items():
                if f & status:
                    pos[axis] = p
                    break
            else:
                logging.warning("Status %X doesn't contain position information", status)
                return # don't change position

        # it's read-only, so we change it via _value
        self.position._value = pos
        self.position.notify(self.position.value)

    def _waitNoMotion(self, timeout=None):
        """
        Block as long as the controller reports motion
        timeout (0 < float): maximum time to wait for the end of the motion
        """
        start = time.time()

        # Read until end of motion
        while True:
            _, status = self.GetStatus()
            if not (status & STA_IN_MOTION):
                return

            if timeout is not None and (time.time() > start + timeout):
                raise IOError("Device still in motion after %g s" % (timeout,))

            # Give it a small break
            time.sleep(0.05) # 20Hz

    @isasync
    def moveRel(self, shift):
        if not shift:
            return model.InstantaneousFuture()
        self._checkMoveRel(shift)
        # TODO move to the +N next position? (and modulo number of axes)
        raise NotImplementedError("Relative move on enumerated axis not supported")

    @isasync
    def moveAbs(self, pos):
        if not pos:
            return model.InstantaneousFuture()
        self._checkMoveAbs(pos)

        return self._executor.submit(self._doMovePos, pos.values()[0])

    def stop(self, axes=None):
        self._executor.cancel()

    def _doMovePos(self, pos):
        jogp = self._pos_to_jog[pos]
        self.MoveJog(jogp)
        self._waitNoMotion(10) # by default, a move lasts ~0.5 s
        self._updatePosition()

    @staticmethod
    def _openSerialPort(port):
        """
        Opens the given serial port the right way for a Thorlabs APT device.
        port (string): the name of the serial port (e.g., /dev/ttyUSB0)
        return (serial): the opened serial port
        """
        # For debugging purpose
        if port == "/dev/fake":
            return MFF102Simulator(timeout=1)

        ser = serial.Serial(
            port=port,
            baudrate=115200,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            rtscts=True,
            timeout=1 #s
        )

        # Purge (as recommended in the documentation)
        time.sleep(0.05) # 50 ms
        ser.flush()
        ser.flushInput()
        time.sleep(0.05) # 50 ms

        # Prepare the port
        ser.setRTS()

        return ser

    def _getSerialPort(self, sn):
        """
        sn (str): serial number of the device
        return (str): serial port name (eg: "/dev/ttyUSB0" on Linux)
        """
        if sys.platform.startswith('linux'):
            # Look for each USB device, if the serial number is good
            sn_paths = glob.glob('/sys/bus/usb/devices/*/serial')
            for p in sn_paths:
                try:
                    f = open(p)
                    snp = f.read().strip()
                except IOError:
                    logging.debug("Failed to read %s, skipping device", p)
                if snp == sn:
                    break
            else:
                raise ValueError("No USB device with S/N %s" % sn)

            # Deduce the tty:
            # .../3-1.2/serial => .../3-1.2/3-1.2:1.0/ttyUSB1
            sys_path = os.path.dirname(p)
            usb_num = os.path.basename(sys_path)
            tty_paths = glob.glob("%s/%s/ttyUSB?*" % (sys_path, usb_num + ":1.0"))
            if not tty_paths:
                raise ValueError("Failed to find tty for device with S/N %s" % sn)
            tty = os.path.basename(tty_paths[0])

            # Convert to /dev
            # Note: that works because udev rules create a dev with the same name
            # otherwise, we would need to check the char numbers
            return "/dev/%s" % (tty,)
        else:
            # TODO: Windows version
            raise NotImplementedError("OS not yet supported")

    @classmethod
    def scan(cls):
        """
        returns (list of 2-tuple): name, args (sn)
        Note: it's obviously not advised to call this function if a device is already under use
        """
        logging.info("Serial ports scanning for Thorlabs MFFxxx in progress...")
        found = []  # (list of 2-tuple): name, kwargs

        if sys.platform.startswith('linux'):
            # Look for each USB device, if the serial number is potentially good
            sn_paths = glob.glob('/sys/bus/usb/devices/*/serial')
            for p in sn_paths:
                try:
                    f = open(p)
                    snp = f.read().strip()
                except IOError:
                    logging.debug("Failed to read %s, skipping device", p)
                if not (snp.startswith(SN_PREFIX_MFF) and len(snp) == 8):
                    continue

                # Deduce the tty:
                # .../3-1.2/serial => .../3-1.2/3-1.2:1.0/ttyUSB1
                sys_path = os.path.dirname(p)
                usb_num = os.path.basename(sys_path)
                logging.info("Looking at device %s with S/N=%s", usb_num, snp)
                tty_paths = glob.glob("%s/%s/ttyUSB?*" % (sys_path, usb_num + ":1.0"))
                if not tty_paths: # 0 or 1 paths
                    continue
                tty = os.path.basename(tty_paths[0])

                # Convert to /dev
                # Note: that works because udev rules create a dev with the same name
                # otherwise, we would need to check the char numbers
                port = "/dev/%s" % (tty,)

                # open and try to communicate
                try:
                    dev = cls(name="test", role="test", port=port)
                    _, modl, typ, fmv, notes, hwv, state, nc = dev.GetInfo()
                    found.append((modl, {"sn": snp, "axis": "rz"}))
                except Exception:
                    pass
        else:
            # TODO: Windows version
            raise NotImplementedError("OS not yet supported")

        return found


class MFF102Simulator(object):
    """
    Simulates a MFF102 (+ serial port). Only used for testing.
    Same interface as the serial port
    """
    def __init__(self, timeout=0, *args, **kwargs):
        # we don't care about the actual parameters but timeout
        self.timeout = timeout
        self._output_buf = "" # what the commands sends back to the "host computer"
        self._input_buf = "" # what we receive from the "host computer"

        # internal values
        self._state = {"jog": 1, # 1 or 2
                       }
        self._end_motion = 0 # time at which the current motion end(ed)
        self._add = 0x50 # the address of this device
        self._sn = 37000001
        self._model = "MPP002"
        self._fmv = 0x020304
        self._hwv = 2
        self._nchans = 1

    def write(self, data):
        self._input_buf += data

        self._parseMessages() # will update _input_buf

    def read(self, size=1):
        ret = self._output_buf[:size]
        self._output_buf = self._output_buf[len(ret):]

        if len(ret) < size:
            # simulate timeout
            time.sleep(self.timeout)
        return ret

    def flush(self):
        pass

    def flushInput(self):
        self._output_buf = ""

    def close(self):
        # using read or write will fail after that
        del self._output_buf
        del self._input_buf

    def _parseMessages(self):
        """
        Parse as many messages available in the buffer
        """
        while len(self._input_buf) >= 6:
            # Similar to MFF._ReadMessage()
            # read the first (required) 6 bytes
            msg = self._input_buf[0:7]

            if ord(msg[4]) & 0x80: # long message
                length = struct.unpack("<H", msg[2:4])[0]
                if len(self._input_buf) < 6 + length:
                    return # not yet all the message received
                msg += self._input_buf[6:6 + length + 1]

            # remove the bytes we've just read
            self._input_buf = self._input_buf[len(msg):]
            
            self._processMessage(msg)

    def _createMessage(self, mid, dest=0x1, src=0x50, p1=None, p2=None, data=None):
        """
        msg (APTSet or APTReq): the message definition
        dest (0<int): the destination ID (always 0x50 if directly over USB)
        p1 (None or 0<=int<=255): param1 (passed as byte2)
        p2 (None or 0<=int<=255): param2 (passed as byte3)
        data (None or bytes): data to be send further. Cannot be mixed with p1
          and p2
        return (bytes): full message
        """
        # create the message
        if data is None: # short message
            p1 = p1 or 0
            p2 = p2 or 0
            msg = struct.pack("<HBBBB", mid, p1, p2, dest, src)
        else: # long message
            msg = struct.pack("<HHBB", mid, len(data), dest | 0x80, src) + data

        return msg

    def _processMessage(self, msg):
        """
        process the msg, and put the result in the output buffer
        msg (str): raw message (including header)
        """
        logging.debug("Simulator received: '%s'", ", ".join("%02X" % ord(c) for c in msg))

        mid = struct.unpack("<H", msg[0:2])[0]
        dest = ord(msg[4]) & 0x7f
        if dest != self._add:
            logging.debug("Simulator (add = %X) skipping message for %X",
                          self._add, dest)
            return
        src = ord(msg[5]) & 0x7f
        
        res = None
        try:
            if mid == HW_REQ_INFO.id:
                data = struct.pack('<I8sHI48s12xHHH', self._sn, self._model, 2,
                                   self._fmv, "APT Fake Filter Flipper",
                                   self._hwv, 0, self._nchans)
                res = self._createMessage(HW_REQ_INFO.rid, src, self._add, data=data)
            elif mid == HW_STOP_UPDATEMSGS.id:
                # good, because we don't support update messages ;-)
                pass
            elif mid == HW_START_UPDATEMSGS.id:
                logging.warning("Simulator doesn't support updates messages")
                pass
            elif mid == HW_NO_FLASH_PROGRAMMING.id:
                # nothing to do
                pass
            elif mid == MOT_REQ_STATUSUPDATE.id:
                # compute status from the state
                if self._end_motion < time.time(): # stopped
                    jog_to_sta = {1: STA_FWD_HLS,
                                  2: STA_RVS_HLS}
                    status = jog_to_sta[self._state["jog"]]
                else: # moving
                    # that's what the hardware reports when jog moves!
                    status = STA_FWD_MOT
                status |= STA_CHA_ENB
                data = struct.pack('<HiiI', 1, 0, 0, status)
                res = self._createMessage(MOT_REQ_STATUSUPDATE.rid,
                                          src, self._add, data=data)
            elif mid == MOT_MOVE_JOG.id:
                chan, jog = struct.unpack('BB', msg[2:4])
                if chan != 1:
                    raise ValueError("Channel = %d" % chan)
                if not jog in [1, 2]:
                    raise ValueError("jog = %d" % jog)
                # simulate a move
                self._state["jog"] = jog
                self._end_motion = time.time() + 1 # 1s move
                # no output
            else:
                logging.warning("Message '%X' unknown", mid)
        except Exception:
            logging.exception("Simulator failed on message %X", mid)

        # add the response end
        if res is not None:
            self._output_buf += res

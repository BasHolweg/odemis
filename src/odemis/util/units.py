# -*- coding: utf-8 -*-
"""
Created on 20 Feb 2012

@author: Éric Piel

Various utility functions for displaying numbers (with and without units).

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

"""
from __future__ import division
import collections
import logging
import math

SI_PREFIXES = {9: u"G",
               6: u"M",
               3: u"k",
               0: u"",
               -3: u"m",
               -6: u"µ",
               -9: u"n",
               -12: u"p"}

# The following units should ignore SI formatting
IGNORE_UNITS = (None, "", "px", "C", u"°C", "rad", "%", "nm")

def round_significant(x, n):
    """
    Round a number to n significant figures
    """
    if x == 0:
        return 0

    return round(x, int(n - math.ceil(math.log10(abs(x)))))

def round_down_significant(x, n):
    """
    Round a number to n significant figures making sure it's smaller
    """
    if x == 0:
        return 0

    exp = n - math.ceil(math.log10(abs(x)))
    if x > 0:
        ret = math.floor(x * 10 ** exp) / (10 ** exp)
    else:
        ret = math.ceil(x * 10 ** exp) / (10 ** exp)
    # assert(abs(ret) <= abs(x))
    return ret

def get_si_scale(x):
    """ This function returns the best fitting SI scale for the given numerical
    value x.
    Returns a (float, string) tuple: (divisor , SI prefix)
    """
    if x == 0:
        return (1, u"")

    most_significant = math.floor(math.log10(abs(x)))
    prefix_order = (most_significant // 3) * 3 # rounding to multiple of 3
    prefix_order = max(-12, min(prefix_order, 9)) # clamping
    return (10 ** prefix_order), SI_PREFIXES[int(prefix_order)]

def to_si_scale(x):
    """ Scale the given value x to the best fitting metric prefix.
    Return a tuple: (scaled value of x, prefix)
    """
    divisor, prefix = get_si_scale(x)
    return x / divisor, prefix

def si_scale_list(values):
    """ Scales a list of numerical values using the same metrix scale """
    if values:
        marker = max(values)
        divisor, prefix = get_si_scale(marker)
        return [v / divisor for v in values], prefix
    return None, u""

def to_string_si_prefix(x, sig=None):
    """
    Convert a number to a string with the most appropriate SI prefix appended
    ex: 0.0012 -> "1.2 m"
    x (float): number
    return (string)
    """
    value, prefix = to_si_scale(x)
    return u"%s %s" % (to_string_pretty(value, sig), prefix)

def to_string_pretty(x, sig=None, unit=None):
    """ Convert a number to a string as int or float as most appropriate

    :param sig: (int) The number of significant decimals

    """

    if x == 0:
        # don't consider this a float
        return u"0"

    if sig is not None:
        x = round_significant(x, sig)

    # so close from an int that it's very likely one?
    if abs(x - round(x)) < 1e-5 and abs(x) >= 1:
        x = int(round(x)) # avoid the .0

    if isinstance(x, float):

        str_val = "%r" % x

        if unit in IGNORE_UNITS:
            return str_val
        else:
            # Get the scale that a readable (formatted) string would use
            eo, _ = get_si_scale(x)
            scale = int(round(math.log(eo, 10)))

            fn, _, ep = str_val.partition('e')
            ep = int(ep or 0)

            dot_move = ep - scale

            if dot_move and '.' in fn:
                dot_pos = fn.index('.')
                new_dot_pos = dot_pos + dot_move
                fn = fn.replace(".", "")

                if new_dot_pos > len(fn):
                    fn = fn.ljust(new_dot_pos, '0')

                fn = ".".join([fn[:new_dot_pos], fn[new_dot_pos:]])
                return u"%se%d" % (fn.strip('0').strip('.'), scale)
            else:
                return str_val


    return u"%s" % x

def readable_str(value, unit=None, sig=None):
    """
    Convert a value with a unit into a displayable string for the user

    :param value: (number or [number...]): value(s) to display
    :param unit: (None or string): unit of the values. If necessary a SI prefix
        will be used to make the value more readable, unless None is given.
    :param sig: (int or None) The number of significant numbers

    return (string)
    """
    # TODO: convert % to ‰ when small value?
    # check against our black list of units which don't support SI prefix

    if value is None:
        return ""

    if unit in IGNORE_UNITS:
        # don't put SI scaling prefix
        if unit in (None, ""):
            sunit = u""
        else:
            sunit = u" %s" % unit
        if isinstance(value, collections.Iterable):
            # Could use "×" , but less readable than "x"
            return u"%s%s" % (u" x ".join([to_string_pretty(v, sig) for v in value]), sunit)
        else:
            return u"%s%s" % (to_string_pretty(value, sig), sunit)

    # TODO: special case for s: only if < 10

    if isinstance(value, collections.Iterable):
        values, prefix = si_scale_list(value)
        return u"%s %s%s" % (u" x ".join([to_string_pretty(v, sig) for v in values]), prefix, unit)
    else:
        return u"%s%s" % (to_string_si_prefix(value, sig), unit)


def readable_time(seconds):
    """This function translates intervals given in seconds into human readable
    strings.
    seconds (float)
    """
    # TODO: a way to indicate some kind of significant number? (If it's going to
    # last 5 days, the number of seconds is generally pointless)
    result = []

    sign = 1
    if seconds < 0:
        # it's just plain weird, but let's do as well as we can
        logging.warning("Asked to display negative time %f", seconds)
        sign = -1
        seconds = -seconds

    if seconds > 60 * 60 * 24 * 30:
        # just for us to remember to extend the function
        logging.debug("Converting time longer than a month.")

    second, subsec = divmod(seconds, 1)
    msec = round(subsec * 1e3)
    if msec == 1000:
        msec = 0
        second += 1
    if second == 0 and msec == 0:
        # exactly 0 => special case
        return u"0 second"

    minute, second = divmod(second, 60)
    hour, minute = divmod(minute, 60)
    day, hour = divmod(hour, 24)

    if day:
        result.append(u"%d day%s" % (day, u"" if day == 1 else u"s"))

    if hour:
        result.append(u"%d hour%s" % (hour, u"" if hour == 1 else u"s"))

    if minute:
        result.append(u"%d minute%s" % (minute, u"" if minute == 1 else u"s"))

    if second:
        result.append(u"%d second%s" % (second, u"" if second == 1 else u"s"))

    if msec:
        result.append(u"%d ms" % msec)

    if len(result) == 1:
        # simple case
        ret = result[0]
    else:
        # make them "x, x, x and x"
        ret = u"{} and {}".format(u", ".join(result[:-1]), result[-1])

    if sign == -1:
        ret = u"minus " + ret

    return ret

# vim:tabstop=4:shiftwidth=4:expandtab:spelllang=en_gb:spell:

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 02 2022

@author: Kornee Kleijwegt

Copyright © 2022 Kornee Kleijwegt, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.

----------------------------------------------------------------------------------------------------

This file allows to filter a pytest test report such that only the errors, failures etc. are outputted.
This means that the part of the report regarding the passed tests and on the current progress are filtered out.

Run the script by calling it with as first (and only) input argument the location of a text based file where the pytest test report can be found.
    > ./pytest_log_filter.py ~/location/test_file.txt
"""
import re
import sys


def filter_test_log(log_txt):
    """
    Filters a log file to only include the parts with a failure warning etc.
    Text regarding progress and passed tests are excluded.
    If the output log is too long the middle part of the logged is not returned.

    :param log (str): String with the log of a single test file.
    :return (str): filtered log

    Example of filtering:
    This function expects a pytest test report of single test file as input, a string such as the example below:

        ============================= test session starts ==============================
        platform linux -- Python 3.6.9, pytest-3.3.2, py-1.5.2, pluggy-0.6.0 -- /usr/bin/python3
        cachedir: ../../../../.cache
        rootdir: /home/kleijwegt/development/odemis, inifile:
        collecting ... collected 9 items

        comp_canvas_test.py::TestDblMicroscopeCanvas::test_basic_display PASSED  [ 11%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_basic_move PASSED     [ 22%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_conversion_functions PASSED [ 33%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_conversion_methods PASSED [ 44%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_crosshair PASSED      [ 55%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_pyramidal_3x2 FAILED  [ 66%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_pyramidal_one_tile FAILED [ 77%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_pyramidal_zoom FAILED [ 88%]
        comp_canvas_test.py::TestDblMicroscopeCanvas::test_zoom_move PASSED      [100%]


        =================================== FAILURES ===================================
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:742: AssertionError: Tuples differ: (0, 128, 127) != (0, 76, 179)
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:501: AssertionError: Tuples differ: (0, 128, 127) != (0, 76, 179)
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:604: AssertionError: Tuples differ: (0, 127, 128) != (0, 179, 76)
        ====================== 3 failed, 6 passed in 7.43 seconds ======================

    For a full test report only failures, errors etc. are of interest, therefore only this part is outputted, the rest is
    filtered. An example of such an output string is the following example:

        =================================== FAILURES ===================================
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:742: AssertionError: Tuples differ: (0, 128, 127) != (0, 76, 179)
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:501: AssertionError: Tuples differ: (0, 128, 127) != (0, 76, 179)
        /home/kleijwegt/development/odemis/src/odemis/gui/test/comp_canvas_test.py:604: AssertionError: Tuples differ: (0, 127, 128) != (0, 179, 76)
        ====================== 3 failed, 6 passed in 7.43 seconds ======================

    """
    test_results = re.search("\n\n===(.+?)((.|\n)*)==", log_txt).group()  # Filter to only get the usefull parts of the test report
    test_results = test_results.lstrip("\n")  # Remove preceding empty lines
    # Only display when there is more to tell than just passed test cases (meaning failures, warning etc.)
    if "\n" in test_results:  # Only a message with multiple lines contains interesting information.
        if len(test_results) > 2000:
            test_results = test_results[0:1000] + \
                           "\n\n \t--The middle part of the logging message is not displayed because it is too long for the default summary report. -- \n\n" + \
                            test_results[-2000:] + \
                           "\n \t--Logging message is too long for the default summary report, the middle part is not displayed -- \n\n"
        return test_results


if __name__ == "__main__":
    # Run with as first argument the path to a log txt file of a single test file.
    with open(sys.argv[1]) as f:
        log_txt = f.read()
    output = filter_test_log(log_txt)
    if output:
        print(output)

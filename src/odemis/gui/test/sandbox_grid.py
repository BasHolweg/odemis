#-*- coding: utf-8 -*-

"""
.. codeauthor:: Rinze de Laat <delaat@delmic.com>

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

import unittest

import odemis.gui.test as test
from odemis.gui.test import gui_loop


test.goto_manual()


class GridPanelTestCase(test.GuiTestCase):

    frame_class = test.test_gui.xrcgrid_frame

    def test(self):

        sizer = self.frame.grid_panel.GetSizer()
        sizer.SetEmptyCellSize((0, 0))

        self.frame.blue.Hide()
        self.frame.purple.Hide()
        self.frame.red.Hide()
        self.frame.brown.SetSize(sizer.GetSize())
        self.frame.brown.SetPosition((0, 0))

        # sizer.RemoveGrowableRow(0)
        # sizer.RemoveGrowableCol(0)
        sizer.Layout()

        gui_loop()

        # self.frame.brown.Show()
        # self.frame.purple.Show()


if __name__ == "__main__":
    unittest.main()

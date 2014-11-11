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

import logging
import unittest

import wx
from odemis.gui.cont.views import ViewPortController

import odemis.gui.test as test
from odemis.gui.test import gui_loop


test.goto_manual()


class GridPanelTestCase(test.GuiTestCase):

    frame_class = test.test_gui.xrcgrid_frame
    # test.set_log_level(logging.DEBUG)

    def test_grid_view(self):

        test.set_sleep_time(200)

        gp = self.frame.grid_panel

        gui_loop()
        csize = self.frame.ClientSize

        # Hide 1 windows

        # Hide top left
        gp.hide_viewport(self.frame.red)
        self.assertEqual(self.frame.blue.Size, (csize.x, gp.grid_layout.tr.size.y))
        self.assertEqual(self.frame.purple.Size, gp.grid_layout.bl.size)
        self.assertEqual(self.frame.brown.Size, gp.grid_layout.br.size)
        gui_loop()
        gp.show_viewport(self.frame.red)

        # Hide top right
        gp.hide_viewport(self.frame.blue)
        self.assertEqual(self.frame.red.Size, (csize.x, gp.grid_layout.tl.size.y))
        self.assertEqual(self.frame.purple.Size, gp.grid_layout.bl.size)
        self.assertEqual(self.frame.brown.Size, gp.grid_layout.br.size)
        gui_loop()
        gp.show_viewport(self.frame.blue)

        # Hide bottom left
        gp.hide_viewport(self.frame.purple)
        self.assertEqual(self.frame.red.Size, gp.grid_layout.tl.size)
        self.assertEqual(self.frame.blue.Size, gp.grid_layout.tr.size)
        self.assertEqual(self.frame.brown.Size, (csize.x, gp.grid_layout.br.size.y))
        gui_loop()
        gp.show_viewport(self.frame.purple)

        # Hide bottom right
        gp.hide_viewport(self.frame.brown)
        self.assertEqual(self.frame.red.Size, gp.grid_layout.tl.size)
        self.assertEqual(self.frame.blue.Size, gp.grid_layout.tr.size)
        self.assertEqual(self.frame.purple.Size, (csize.x, gp.grid_layout.bl.size.y))
        gui_loop()
        gp.show_viewport(self.frame.brown)

        # Hide 2 windows

        # Hide top
        gp.hide_viewport(self.frame.red)
        gp.hide_viewport(self.frame.blue)
        self.assertEqual(self.frame.purple.Size, (gp.grid_layout.bl.size.x, csize.y))
        self.assertEqual(self.frame.brown.Size, (gp.grid_layout.br.size.x, csize.y))
        gui_loop()
        gp.show_viewport(self.frame.red)
        gp.show_viewport(self.frame.blue)

        # Hide right
        gp.hide_viewport(self.frame.blue)
        gp.hide_viewport(self.frame.brown)
        self.assertEqual(self.frame.red.Size, (csize.x, gp.grid_layout.tl.size.y))
        self.assertEqual(self.frame.purple.Size, (csize.x, gp.grid_layout.bl.size.y))
        gui_loop()
        gp.show_viewport(self.frame.brown)
        gp.show_viewport(self.frame.blue)

        # Hide bottom
        gp.hide_viewport(self.frame.purple)
        gp.hide_viewport(self.frame.brown)
        self.assertEqual(self.frame.red.Size, (gp.grid_layout.tl.size.x, csize.y))
        self.assertEqual(self.frame.blue.Size, (gp.grid_layout.tr.size.x, csize.y))
        gui_loop()
        gp.show_viewport(self.frame.brown)
        gp.show_viewport(self.frame.purple)

        # Hide left
        gp.hide_viewport(self.frame.red)
        gp.hide_viewport(self.frame.purple)
        self.assertEqual(self.frame.blue.Size, (csize.x, gp.grid_layout.tr.size.y))
        self.assertEqual(self.frame.brown.Size, (csize.x, gp.grid_layout.br.size.y))
        gui_loop()
        gp.show_viewport(self.frame.purple)
        gp.show_viewport(self.frame.red)

        # Hide 3 windows

        gp.set_shown_viewports(self.frame.red)
        self.assertEqual(self.frame.red.Size, csize)
        gui_loop()

        gp.set_shown_viewports(self.frame.blue)
        self.assertEqual(self.frame.blue.Size, csize)
        gui_loop()

        gp.set_shown_viewports(self.frame.purple)
        self.assertEqual(self.frame.purple.Size, csize)
        gui_loop()

        gp.set_shown_viewports(self.frame.brown)
        self.assertEqual(self.frame.brown.Size, csize)
        gui_loop()

        gp.set_shown_viewports(self.frame.yellow)
        self.assertEqual(self.frame.yellow.Size, csize)
        gui_loop()

        gp.show_grid_viewports()

    def test_grid_edit(self):

        test.set_sleep_time(200)

        gp = self.frame.grid_panel
        gui_loop()

        self.assertRaises(ValueError, gp.swap_viewports, self.frame.red, self.frame.yellow)
        gp.swap_viewports(self.frame.red, self.frame.yellow)

        self.frame.red.Hide()
        self.frame.yellow.Show()
        gp.swap_viewports(self.frame.red, self.frame.yellow)
        gui_loop()
        self.assertEqual(self.frame.yellow.Position, (0, 0))

        self.frame.red.Show()
        self.frame.yellow.Hide()
        gp.swap_viewports(self.frame.red, self.frame.yellow)
        gui_loop()
        self.assertEqual(self.frame.red.Position, (0, 0))

    def test_grid_resize(self):

        test.set_sleep_time(200)

        gp = self.frame.grid_panel
        gui_loop()

        self.frame.SetSize((600, 600))
        self.frame.Center()
        gui_loop()



if __name__ == "__main__":
    unittest.main()

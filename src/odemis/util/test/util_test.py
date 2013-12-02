#-*- coding: utf-8 -*-

"""
@author: Rinze de Laat

Copyright © 2013 Rinze de Laat, Delmic

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

from odemis import util
import os
import unittest


# Bounding box clipping test data generation

def tp(trans, ps):
    """ Translate points ps using trans """
    r = []
    i = 0
    for p in ps:
        r.append(p + trans[i])
        i = (i + 1) % len(trans)
    return tuple(r)

# First we define a bounding boxes, at different locations
bounding_boxes = [(-2, -2, 0, 0),
                  (-1, -1, 1, 1),
                  (0, 0, 2, 2),
                  (2, 2, 4, 4)]

# From this, we generate boxes that are situated all around these
# bounding boxes, but that do not touch or overlap them.

def relative_boxes(bb):

    t_left = [(-3, i) for i in range(-3, 4)]
    to_the_left = [tp(t, bb) for t in t_left]

    t_top = [(i, -3) for i in range(-3, 4)]
    to_the_top = [tp(t, bb) for t in t_top]

    t_right = [(3, i) for i in range(-3, 4)]
    to_the_right = [tp(t, bb) for t in t_right]

    t_bottom = [(i, 3) for i in range(-3, 4)]
    to_the_bottom = [tp(t, bb) for t in t_bottom]

    outside_boxes = to_the_left + to_the_top + to_the_right + to_the_bottom

    # Selection boxes that touch the outside of the bounding box
    touch_left = [tp((1, 0), b) for b in to_the_left[1:-1]]
    touch_top = [tp((0, 1), b) for b in to_the_top[1:-1]]
    touch_right = [tp((-1, 0), b) for b in to_the_right[1:-1]]
    touch_bottom = [tp((0, -1), b) for b in to_the_bottom[1:-1]]

    touching_boxes = touch_left + touch_top + touch_right + touch_bottom

    # Partial overlapping boxes
    overlap_left = [tp((1, 0), b) for b in touch_left[1:-1]]
    overlap_top = [tp((0, 1), b) for b in touch_top[1:-1]]
    overlap_right = [tp((-1, 0), b) for b in touch_right[1:-1]]
    overlap_bottom = [tp((0, -1), b) for b in touch_bottom[1:-1]]

    overlap_boxes = overlap_left + overlap_top + overlap_right + overlap_bottom

    return outside_boxes, touching_boxes, overlap_boxes

class AlmostEqualTestCase(unittest.TestCase):

    def test_simple(self):
        in_exp = {(0., 0): True,
                  (-5, -5.): True,
                  (1., 1. - 1e-9): True,
                  (1., 1. - 1e-3): False,
                  (1., 1. + 1e-3): False,
                  (-5e-8, -5e-8 + 1e-19): True,
                  (5e18, 5e18 + 1): True,
                  }
        for i, eo in in_exp.items():
            o = util.almost_equal(*i)
            self.assertEqual(o, eo, "Failed to get correct output for %s" % (i,))


class CanvasTestCase(unittest.TestCase):

    def test_clipping(self):

        tmp = "{}: {} - {} -> {}"

        for bb in bounding_boxes:
            outside, touching, overlap = relative_boxes(bb)

            for b in outside:
                r = util.rect_intersect(b, bb)
                msg = tmp.format("outside", b, bb, r)
                self.assertIsNone(r, msg)

            for b in touching:
                r = util.rect_intersect(b, bb)
                msg = tmp.format("touching", b, bb, r)
                self.assertIsNone(r, msg)

            for b in overlap:
                r = util.rect_intersect(b, bb)
                msg = tmp.format("overlap", b, bb, r)
                self.assertIsNotNone(r, msg)

                # 'Manual' checks
                if bb == (-1, -1, 1, 1):
                    if b[:2] == (-2, -2):
                        self.assertEqual(r, (-1, -1, 0, 0), msg)
                    elif b[:2] == (0, -1):
                        self.assertEqual(r, (0, -1, 1, 1), msg)
                    elif b[:2] == (0, 0):
                        self.assertEqual(r, (0, 0, 1, 1), msg)

            # full and exact overlap
            b = bb
            r = util.rect_intersect(b, bb)
            self.assertEqual(r, bb)

            # inner overlap
            b = (bb[0] + 1, bb[1] + 1, bb[2], bb[3])
            r = util.rect_intersect(b, bb)
            self.assertEqual(r, b)

            # overflowing overlap
            b = (bb[0] - 1, bb[1] - 1, bb[2] + 1, bb[2] + 1)
            r = util.rect_intersect(b, bb)
            self.assertEqual(r, bb)

if __name__ == "__main__":
    unittest.main()

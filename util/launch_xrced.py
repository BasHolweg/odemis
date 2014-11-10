#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This custom XRCED launcher allows a small wx function to be wrapped
so it provides a little extra needed functionality.

XRC sometimes need to check if a node contains a filename. It does so by
checking node types. This works fine, until we start working with custom
controls, of which XRC knows nothing by default.

The little wrapper added to the pywxrc.XmlResourceCompiler.NodeContainsFilename
method,  will return true if it contains a value ending with '.png', indicating
the content is an PNG image.

"""

import os
import sys

#pylint: disable=W0105
from wx.tools.XRCed.globals import set_debug


if __name__ == '__main__':
    try:
        from wx.tools.XRCed.xrced import main
        set_debug(True)
    except ImportError:
        try:
            from XRCed.xrced import main
        except ImportError:
            print >> sys.stderr, 'XRCed parent directory must be in PYTHONPATH'
            raise
        sys.modules['wx.tools.XRCed'] = sys.modules['XRCed']

    from wx.tools import pywxrc

    # The XRCEDPATH environment variable is used to define additional plugin directories
    xrced_path = os.getenv('XRCEDPATH')
    this_path = os.path.dirname(__file__)
    os.environ['XRCEDPATH'] = xrced_path or os.path.join(this_path, "../src/odemis/gui/xmlh")
    print "'XRCEDPATH' is set to %s" % os.getenv('XRCEDPATH')

    # Move this to a separate launcher so it can be spread with
    # odemis

    def ncf_decorator(ncf):
        def wrapper(self, node):
            if (
                node.firstChild and
                node.firstChild.nodeType == 3 and
                node.firstChild.nodeValue.lower().endswith(".png")
            ):
                # print node.firstChild.nodeValue
                return True
            return ncf(self, node)
        return wrapper

    pywxrc.XmlResourceCompiler.NodeContainsFilename = ncf_decorator(pywxrc.XmlResourceCompiler.NodeContainsFilename)

    main()

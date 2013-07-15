# -*- coding: utf-8 -*-

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
Odemis. If not, see http://www.gnu.org/licenses/.

"""

# This module is used to construct Delmic's custom FoldPanelBar according to
# the definition in the XRC configuration file.
#
# This module is used both by Odemis' GUI and XRCED.

import ast
import logging
import odemis.gui.comp.buttons as btns
import odemis.gui.comp.foldpanelbar as fpb
import odemis.gui.comp.viewport as vport
import odemis.gui.comp.slider as slide
import odemis.gui.comp.stream as strm
import odemis.gui.comp.text as txt
import odemis.gui.img.data as img
import wx
import wx.combo
import wx.lib.buttons
import wx.xrc as xrc


HANDLER_CLASS_LIST = []

##################################
# Fold Panel Bar related Handlers
##################################

class StreamPanelXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddStyle("wxTAB_TRAVERSAL", wx.TAB_TRAVERSAL)
        #self.AddStyle("CP_GTK_EXPANDER", CP_GTK_EXPANDER)
        #self.AddStyle("CP_DEFAULT_STYLE", CP_DEFAULT_STYLE)
        #self.AddStyle("CP_NO_TLW_RESIZE", CP_NO_TLW_RESIZE)
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "StreamPanel")

        return capable

    def DoCreateResource(self):
        assert self.GetInstance() is None

        parent_window = self.GetParentAsWindow()
        # Now create the object
        panel = strm.StreamPanel(parent_window,
                                         self.GetID(),
                                         self.GetText('label'),
                                         self.GetPosition(),
                                         self.GetSize(),
                                         self.GetStyle("style", wx.TAB_TRAVERSAL),
                                         #self.GetStyle('exstyle'),
                                         name=self.GetName(),
                                         collapsed=self.GetBool('collapsed')
                                         )

        # These two things should be done in either case:
        # Set standard window attributes
        self.SetupWindow(panel)
        panel.finalize()

        parent_window.add_stream(panel)
        # Create any child windows of this node
        # deprecated: all children are hard-coded
        #self.CreateChildren(panel.get_panel())

        return panel
HANDLER_CLASS_LIST.append(StreamPanelXmlHandler)


class FoldPanelBarXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        #self._isInside = False

    def CanHandle(self, node):
        # return not self._isInside and self.IsOfClass(node, 'wx.lib.foldpanelbar.FoldPanelBar') or \
        #        self._isInside and self.IsOfClass(node, 'foldpanel')
        return self.IsOfClass(node, 'FoldPanelBar')


    # Process XML parameters and create the object
    def DoCreateResource(self):

        if self.GetClass() == 'FoldPanelBar':
            #print "Creating FoldpanelBar"
            w = fpb.FoldPanelBar(self.GetParentAsWindow(),
                                 self.GetID(),
                                 self.GetPosition(),
                                 self.GetSize(),
                                 self.GetStyle())

            self.SetupWindow(w)

            parent = w.GetParent()
            if parent.__class__ == wx.ScrolledWindow:
                parent.EnableScrolling(False, True)
                parent.SetScrollbars(-1, 10, 1, 1)

            self.CreateChildren(w, False)

            return w
HANDLER_CLASS_LIST.append(FoldPanelBarXmlHandler)


class CaptionBarXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'CaptionBar')


    # Process XML parameters and create the object
    def DoCreateResource(self):

        if self.GetClass() == 'CaptionBar':
            #print "Creating CaptionBar"
            parent = self.GetParentAsWindow()
            w = fpb.CaptionBar(
                            parent,
                            self.GetText('label'),
                            self.GetBool('collapsed')
            )
            self.SetupWindow(w)
            return w
HANDLER_CLASS_LIST.append(CaptionBarXmlHandler)

class FoldPanelItemXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'FoldPanelItem')


    # Process XML parameters and create the object
    def DoCreateResource(self):

        if self.GetClass() == 'FoldPanelItem':
            #print "Creating FoldpanelItem"
            parent = self.GetParentAsWindow()
            w = fpb.FoldPanelItem(parent,
                                  self.GetID(),
                                  self.GetPosition(),
                                  self.GetSize(),
                                  self.GetStyle(),
                                  self.GetText('label'),
                                  self.GetBool('collapsed'),
                                  nocaption=self.GetBool('nocaption'))
            self.SetupWindow(w)

            self.CreateChildren(w, False)

            # Move all the FoldPanelItem children to the main sizer
            w.children_to_sizer()
            parent.add_item(w)
            return w
HANDLER_CLASS_LIST.append(FoldPanelItemXmlHandler)


class StreamBarXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'StreamBar')


    # Process XML parameters and create the object
    def DoCreateResource(self):

        if self.GetClass() == 'StreamBar':
            #print "Creating FoldpanelItem"
            parent = self.GetParentAsWindow()
            w = strm.StreamBar(parent,
                                 self.GetID(),
                                 self.GetPosition(),
                                 self.GetSize(),
                                 self.GetStyle(),
                                 add_button=self.GetBool('add_button'))
            self.SetupWindow(w)
            # 'Dirty' fix for the hard coded 'add stream' child button
            if self.GetBool('add_button'):
                w.btn_add_stream.SetBackgroundColour(w.GetBackgroundColour())
            parent.add_item(w)
            return w
HANDLER_CLASS_LIST.append(StreamBarXmlHandler)


################################
# ImageButton sub class handlers
################################

class GenBitmapButtonHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'wx.lib.buttons.GenBitmapButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = wx.lib.buttons.GenBitmapButton(self.GetParentAsWindow(),
                                    self.GetID(),
                                    bmp,
                                    self.GetPosition(),
                                    self.GetSize(),
                                    self.GetStyle())

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)

        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(GenBitmapButtonHandler)


class ImageButtonHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'ImageButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = btns.ImageButton(self.GetParentAsWindow(),
                            self.GetID(),
                            bmp,
                            pos=self.GetPosition(),
                            size=self.GetSize(),
                            style=self.GetStyle(),
                            label_delta=self.GetLong('delta'))

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("hover"):
            bmp = self.GetBitmap("hover")
            w.SetBitmapHover(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)


        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(ImageButtonHandler)


class ImageToggleButtonHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'ImageToggleButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = btns.ImageToggleButton(self.GetParentAsWindow(),
                                   self.GetID(),
                                   bmp,
                                   pos=self.GetPosition(),
                                   size=self.GetSize(),
                                   style=self.GetStyle(),
                                   label_delta=self.GetLong('delta'))

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("hover"):
            bmp = self.GetBitmap("hover")
            w.SetBitmapHover(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)

        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(ImageToggleButtonHandler)


class ImageTextButtonHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles
        self.AddStyle('wxALIGN_LEFT', wx.ALIGN_LEFT)
        self.AddStyle('wxALIGN_RIGHT', wx.ALIGN_RIGHT)
        self.AddStyle('wxALIGN_CENTRE', wx.ALIGN_CENTRE)

    def CanHandle(self, node):
        return self.IsOfClass(node, 'ImageTextButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = btns.ImageTextButton(self.GetParentAsWindow(),
                                 self.GetID(),
                                 bmp,
                                 pos=self.GetPosition(),
                                 size=self.GetSize(),
                                 style=self.GetStyle(),
                                 label=self.GetText('label'),
                                 label_delta=self.GetLong('delta'))

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("hover"):
            bmp = self.GetBitmap("hover")
            w.SetBitmapHover(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)

        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(ImageTextButtonHandler)


class ImageTextToggleButtonHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles
        self.AddStyle('wxALIGN_LEFT', wx.ALIGN_LEFT)
        self.AddStyle('wxALIGN_RIGHT', wx.ALIGN_RIGHT)
        self.AddStyle('wxALIGN_CENTRE', wx.ALIGN_CENTRE)

        self.klass = btns.ImageTextToggleButton

    def CanHandle(self, node):
        return self.IsOfClass(
            node, 'ImageTextToggleButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = self.klass(self.GetParentAsWindow(),
                       self.GetID(),
                       bmp,
                       pos=self.GetPosition(),
                       size=self.GetSize(),
                       style=self.GetStyle(),
                       label=self.GetText('label'),
                       label_delta=self.GetLong('delta'))

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("hover"):
            bmp = self.GetBitmap("hover")
            w.SetBitmapHover(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)

        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(ImageTextToggleButtonHandler)


class TabButtonHandler(ImageTextToggleButtonHandler):

    def __init__(self):
        ImageTextToggleButtonHandler.__init__(self)
        self.klass = btns.TabButton

    def CanHandle(self, node):
        return self.IsOfClass(node, 'TabButton')
HANDLER_CLASS_LIST.append(TabButtonHandler)


class ViewButtonHandler(ImageTextToggleButtonHandler):

    def __init__(self):
        ImageTextToggleButtonHandler.__init__(self)
        self.klass = btns.ViewButton

    def CanHandle(self, node):
        return self.IsOfClass(node, 'ViewButton')
HANDLER_CLASS_LIST.append(ViewButtonHandler)


class PopupImageButtonHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'PopupImageButton')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        bmp = wx.NullBitmap
        if self.GetParamNode("bitmap"):
            bmp = self.GetBitmap("bitmap")

        w = btns.PopupImageButton(self.GetParentAsWindow(),
                                  self.GetID(),
                                  bmp,
                                  pos=self.GetPosition(),
                                  size=self.GetSize(),
                                  style=self.GetStyle())

        if self.GetParamNode("selected"):
            bmp = self.GetBitmap("selected")
            w.SetBitmapSelected(bmp)

        if self.GetParamNode("hover"):
            bmp = self.GetBitmap("hover")
            w.SetBitmapHover(bmp)

        if self.GetParamNode("focus"):
            bmp = self.GetBitmap("focus")
            w.SetBitmapFocus(bmp)

        if self.GetParamNode("disabled"):
            bmp = self.GetBitmap("disabled")
            w.SetBitmapDisabled(bmp)

        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(PopupImageButtonHandler)


class SuggestTextCtrlHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'SuggestTextCtrl')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        w = txt.SuggestTextCtrl(self.GetParentAsWindow(),
                                id=self.GetID(),
                                value=self.GetText('value'),
                                pos=self.GetPosition(),
                                size=self.GetSize(),
                                style=self.GetStyle(),
                                choices=[str(i) for i in range(2)])
        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(SuggestTextCtrlHandler)


class UnitIntegerCtrlHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'UnitIntegerCtrl')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        val = int(self.GetText('value'))

        w = txt.UnitIntegerCtrl(self.GetParentAsWindow(),
                                id=self.GetID(),
                                value=val,
                                pos=self.GetPosition(),
                                size=self.GetSize(),
                                style=self.GetStyle(),
                                unit=self.GetText('unit'),
                                min_val=self.GetLong('min'),
                                max_val=self.GetLong('max'))
        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(UnitIntegerCtrlHandler)


class UnitFloatCtrlHandler(xrc.XmlResourceHandler):

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        return self.IsOfClass(node, 'UnitFloatCtrl')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        val = float(self.GetText('value'))

        w = txt.UnitFloatCtrl(self.GetParentAsWindow(),
                              id=self.GetID(),
                              value=val,
                              pos=self.GetPosition(),
                              size=self.GetSize(),
                              style=self.GetStyle(),
                              unit=self.GetText('unit'),
                              min_val=self.GetLong('min'),
                              max_val=self.GetLong('max'))
        self.SetupWindow(w)
        return w
HANDLER_CLASS_LIST.append(UnitFloatCtrlHandler)


##################################
# Canvas Handlers
##################################

class MicroscopeViewportXmlHandler(xrc.XmlResourceHandler):

    klass = vport.MicroscopeViewport

    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddStyle("wxTAB_TRAVERSAL", wx.TAB_TRAVERSAL)
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        return self.IsOfClass(node, "MicroscopeViewport")

    def DoCreateResource(self):
        assert self.GetInstance() is None

        # Now create the object
        panel = self.klass(self.GetParentAsWindow(),
                           id=self.GetID(),
                           pos=self.GetPosition(),
                           size=self.GetSize(),
                           style=self.GetStyle())
        self.SetupWindow(panel)
        return panel
HANDLER_CLASS_LIST.append(MicroscopeViewportXmlHandler)


class SecomViewportXmlHandler(MicroscopeViewportXmlHandler):

    klass = vport.SecomViewport

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        return self.IsOfClass(node, "SecomViewport")
HANDLER_CLASS_LIST.append(SecomViewportXmlHandler)


class SparcAcquisitionViewportXmlHandler(MicroscopeViewportXmlHandler):

    klass = vport.SparcAcquisitionViewport

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        return self.IsOfClass(node, "SparcAcquisitionViewport")
HANDLER_CLASS_LIST.append(SparcAcquisitionViewportXmlHandler)


class SparcAnalysisViewportXmlHandler(MicroscopeViewportXmlHandler):

    klass = vport.SparcAnalysisViewport

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        return self.IsOfClass(node, "SparcAnalysisViewport")
HANDLER_CLASS_LIST.append(SparcAnalysisViewportXmlHandler)

class SparcAlignViewportXmlHandler(MicroscopeViewportXmlHandler):

    klass = vport.SparcAlignViewport

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        return self.IsOfClass(node, "SparcAlignViewport")
HANDLER_CLASS_LIST.append(SparcAlignViewportXmlHandler)


##################################
# Sliders
##################################

class UnitIntegerSliderHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "UnitIntegerSlider")
        return capable

    def DoCreateResource(self):
        assert self.GetInstance() is None

        # parent, id=wx.ID_ANY, value=0.0, val_range=(0.0, 1.0),
        #          size=(-1, -1), pos=wx.DefaultPosition, style=wx.NO_BORDER,
        #          name="Slider", scale=None, t_class=NumberTextCtrl,
        #          t_size=(50, -1), unit="", accuracy=0):

        val = int(self.GetText('value') or 0)
        rng = (self.GetLong('min'), self.GetLong('max'))
        text_size = ast.literal_eval(self.GetText('text_size') or "50, -1")

        if rng[0] == rng[1]:
            rng = (rng[0], rng[1] + 1)

        # Now create the object
        slider = slide.UnitIntegerSlider(self.GetParentAsWindow(),
                                        id=self.GetID(),
                                        pos=self.GetPosition(),
                                        size=self.GetSize(),
                                        style=self.GetStyle(),
                                        value=val,
                                        unit=self.GetText('unit'),
                                        val_range=rng,
                                        scale=self.GetText('scale'),
                                        t_size=text_size)
        self.SetupWindow(slider)
        return slider
HANDLER_CLASS_LIST.append(UnitIntegerSliderHandler)


class UnitFloatSliderHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "UnitFloatSlider")
        return capable

    def GetFloat(self, param, defaultv=0):
        # there is a bug in wxWidgets, which doesn't export GetFloat
        # => recreate in Python
        # self, String param, long defaultv=0

        string = self.GetParamValue(param)

        try:
            value = float(string)
        except ValueError:
            logging.error("Float param incorrect %s", string)
        return value

    def DoCreateResource(self):
        assert self.GetInstance() is None

        val = self.GetFloat('value')
        rng = (self.GetFloat('min'), self.GetFloat('max'))
        text_size = ast.literal_eval(self.GetText('text_size') or "50, -1")

        if rng[0] == rng[1]:
            logging.warning("Incorrect range between %r and %r", rng[0], rng[1])
            rng = (rng[0], rng[1] + 1.0)

        accuracy = self.GetLong('accuracy', -1)
        if accuracy == -1:
            accuracy = None

        # Now create the object
        slider = slide.UnitFloatSlider(self.GetParentAsWindow(),
                                       id=self.GetID(),
                                       pos=self.GetPosition(),
                                       size=self.GetSize(),
                                       style=self.GetStyle(),
                                       value=val,
                                       unit=self.GetText('unit'),
                                       val_range=rng,
                                       scale=self.GetText('scale'),
                                       accuracy=accuracy,
                                       t_size=text_size)

        self.SetupWindow(slider)
        return slider
HANDLER_CLASS_LIST.append(UnitFloatSliderHandler)

class VisualRangeSliderHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "VisualRangeSlider")
        return capable

    def DoCreateResource(self):
        assert self.GetInstance() is None
        # Now create the object
        slider = slide.VisualRangeSlider(self.GetParentAsWindow(),
                                         id=self.GetID(),
                                         pos=self.GetPosition(),
                                         size=self.GetSize(),
                                         style=self.GetStyle())
        self.SetupWindow(slider)
        slider.SetForegroundColour(slider.GetForegroundColour())
        return slider
HANDLER_CLASS_LIST.append(VisualRangeSliderHandler)

class BandwidthSliderHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "BandwidthSlider")
        return capable

    def DoCreateResource(self):
        assert self.GetInstance() is None
        # Now create the object
        slider = slide.BandwidthSlider(self.GetParentAsWindow(),
                                       id=self.GetID(),
                                       pos=self.GetPosition(),
                                       size=self.GetSize(),
                                       style=self.GetStyle())
        self.SetupWindow(slider)
        # FIXME: this shouldn't be needed, but without it, the content colour is not set
        slider.SetForegroundColour(slider.GetForegroundColour())
        return slider
HANDLER_CLASS_LIST.append(BandwidthSliderHandler)

####################################################################
# OwnerDrawnComboBox Handlers
#
# This handler was needed because the one included with wxPython
# and XRCed did not allow for the alteration of the down button.
####################################################################

class OwnerDrawnComboBoxHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Specify the styles recognized by objects of this type
        self.AddStyle("wxCB_SIMPLE", wx.CB_SIMPLE)
        self.AddStyle("wxCB_DROPDOWN", wx.CB_DROPDOWN)
        self.AddStyle("wxCB_READONLY", wx.CB_READONLY)
        self.AddStyle("wxCB_SORT", wx.CB_SORT)
        self.AddStyle("wxODCB_STD_CONTROL_PAINT", wx.combo.ODCB_STD_CONTROL_PAINT)
        self.AddStyle("wxODCB_DCLICK_CYCLES", wx.combo.ODCB_DCLICK_CYCLES)
        self.AddStyle("wxTE_PROCESS_ENTER", wx.TE_PROCESS_ENTER)

        self.AddWindowStyles()

    # This method and the next one are required for XmlResourceHandlers
    def CanHandle(self, node):
        capable = self.IsOfClass(node, "OwnerDrawnComboBox")
        return capable

    def DoCreateResource(self):
        assert self.GetInstance() is None
        print self.GetParamValue('content')

        # Now create the object
        new_ctrl = wx.combo.OwnerDrawnComboBox(self.GetParentAsWindow(),
                                            id=self.GetID(),
                                            pos=self.GetPosition(),
                                            size=self.GetSize(),
                                            style=self.GetStyle())
        new_ctrl.SetButtonBitmaps(img.getbtn_downBitmap(), pushButtonBg=False)
        self.SetupWindow(new_ctrl)
        return new_ctrl
HANDLER_CLASS_LIST.append(OwnerDrawnComboBoxHandler)

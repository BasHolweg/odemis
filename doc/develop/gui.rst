************************************
Graphical User Interface (GUI) Layer
************************************

Introduction
============

The Odemis graphical user interface (GUI) is implemented based on the
`Hierarchical Model–view–controller architectural pattern
<https://en.wikipedia.org/wiki/Model%E2%80%93view%E2%80%93controller>`_.
It consists of several main controllers and subcontrollers,
each related to a specific aspect of Odemis functionality.
This data is stored as vigilant attributes, so that the state can be subscribed
to by components to notify them when changes occur.

wxPython is used as an underlying framework for the GUI. Various resources are
defined via XRC. Many of the stock wxWidget components have been extended in
Odemis to better suit Odemis-specific use cases.

.. figure:: stream_view_diag.*
    :width: 100 %
    :align: center

    Simplified UML class diagram of the main GUI classes.

Models
======

When a session of the GUI loads, configuration and components from the backend
are preloaded into a data structure representing the model, or shared state, of the GUI.
This object, inherited from ``MainGUIData``, is defined in ``gui/model/__init__.py``,
and contains attributes to all components loaded in the current microscope configuration.
This object is accessible from many GUI components as the argument ``main_data``.

*MicroscopyGUIData* corresponds to a specific GUI tab, and contains data related
to the shared state of that tab. The following classes refine the model for each
type of tab:

    - **LiveViewGUIData:**

      Represents an interface used to only show the current data from the microscope.

    - **SparcAcquisitionGUIData:**

      Represents an interface used to select a precise area to scan and
      acquire signal. It allows fine control of the shape and density of the scan.
      It is specifically made for the SPARC system.

    - **ChamberGUIData:**

        Represents an interface used by the chamber view tab. 

    - **AnalysisGUIData:**

        Represents an interface used to show the recorded microscope data. Typically
        it represents all the data present in a specific file. All the streams should be StaticStreams

    - **ActuatorGUIData:**

        Represents an interface used to move the actuators of a microscope. It might
        also display one or more views, but it's not required.
        Typically used for the SECOM and SPARC(v2) alignment tabs.

    - **SecomAlignGUIData**, **SparcAlignGUIData**, and **Sparc2AlignGUIData**:

        Represents an interface used for alignment tabs. 

Controllers
===========

.. figure:: gui_cont.*
    :width: 100 %
    :align: center

    UML class diagram of the main controller classes.

The Odemis GUI has several main controllers that allow a user to directly control functionality. 

The top level controllers include:
    #. **TabBarController** (``gui/cont/tabs.py``)
        Controller to handle display of the tab bar at the top of the window. 
    #. **MenuController** (``gui/cont/menu.py``)
        Controller for the Odemis application menu bar. 
    #. **SnapshotController:** (``gui/cont/acquisition.py``)
        Controller to handle snapshot acquisition in a 'global' context.

TabBarController
----------------

    Odemis contains many different tabs in its interface for different modes of operation. Each tab is controlled by a corresponding tab controller, wherein each inheriting from the *Tab* base class. These include:

    - **SecomStreamTab:**

        Provides views for streams and corresponding controls for the SECOM and DELPHI platform.
        To controller the main "Optical" and "SEM" buttons as well as the sample loading,
        it relies on one of these two controllers:

        - **SecomStateController** (``gui/cont/microscope.py``)
        - **DelphiStateController** (``gui/cont/microscope.py``)

    - **SparcAcquisitionTab:**

        Provides views for streams and acquisition controls for the SPARC and SPARCv2 platform. 

    - **ChamberTab:**
    
        Handles control of the measurement chamber state. 
    
    - **AnalysisTab:** 
    
        Handles the loading and displaying of acquisition files.
        
    - **SecomAlignTab:**
    
        Tab for the lens alignment on the SECOM platform.
        
    - **SparcAlignTab:**
    
        Tab for the mirror/fiber alignment on the SPARCv1.
        
    - **Sparc2AlignTab:**
    
        Tab for the mirror/fiber alignment on the SPARCv2. Note that the basic idea is similar to the SPARCv1, but the actual procedure is entirely different.
   
The following controllers are subcontrollers of a tab controller. 

    #. **StreamController** (``gui/cont/streams.py``)
    #. **StreamBarController** (``gui/cont/streams.py``)
    #. **ViewPortController** (``gui/cont/views.py``)
    #. **ViewButtonController** (``gui/cont/views.py``)
    #. Acquisition Controllers (``gui/cont/acquisition.py``)



StreamController
----------------

A controller that is created for each stream. It controls the playing and pausing of a stream, and the display of the stream in the visible views. It also determines whether or not the stream is visible in the stream bar that are displayed, and generates widget controls that control stream VA's. The widgets which are created are generated based on the stream type, and are determined in *conf.data.STREAM_SETTINGS_CONFIG* and *conf.data.HW_SETTINGS_CONFIG*. 


StreamBarController
-------------------

Manages the stream bar, which is a side view tab which shows current streams.
This allows a user to add and remove streams to and from the workspace.
Variants of this controller exist for SECOM and SPARC configurations. 

   - The controller has member functions which correspond to the stream types which can be added, based on the configuration. 
   - ``add_stream()``: This function adds the stream by creating it and its representative *StreamCont*. 
   - Functionality such as the repetition overlay and ROI selection are handled by this controller since they apply to all current streams. 

ViewPortController
------------------

The viewport controller creates a view layout based on the list of available views in the configuration. It is created by the Tab controller, which also has definitions for which views are created based on the main data model configuration. 
        

ViewButtonController
--------------------

The view button controller manages the view button thumbnails on the Odemis left side panel. 


Acquisition Controllers
-----------------------

The acquisition controllers found in ``gui.cont.acquisition`` handle the acquisition process for the SPARC systems.
SECOM acquisition is handled by a separate window in ``gui.win.acquisition``. 


Components
==========
Components in Odemis are extended versions of wxPython components. 

.. TODO * Add the rest of the components 

- **Viewport(wx.Panel):**

- **miccanvas.DblMicroscopyCanvas:**

        Provides a space where microscopy images can be displayed.
        It is draggable, can provide display of various overlays (such as ROA and FOV),
        and handles the display of user tools and custom cursors. 

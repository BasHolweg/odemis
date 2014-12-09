# -*- coding: utf-8 -*-
'''
Created on 18 Jun 2012

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
import Pyro4
from Pyro4.core import oneway
import collections
import logging
import multiprocessing
import os
import threading
import urllib
import weakref


# Pyro4.config.COMMTIMEOUT = 30.0 # a bit of timeout
# There is a problem with threadpool: threads have a timeout on waiting for a
# request. That obviously doesn't make much sense, but also means it's not
# possible to put a global timeout with the current version and threadpool.
# One possibility is to change ._pyroTimeout on each proxy.
# thread is restricted: it can handle at the same time only
# MAXTHREADS concurrent connections.
# After that it simply blocks. As there is one connection per object, it goes fast.
# Multiplex can handle a much larger number of connections, but will always
# execute the requests one at a time, which can cause deadlock when handling
# callbacks.
#Pyro4.config.SERVERTYPE = "multiplex"
Pyro4.config.THREADPOOL_MINTHREADS = 48 # big, because it can block when increasing the pool
Pyro4.config.THREADPOOL_MAXTHREADS = 128
# TODO make sure Pyro can now grow the pool: it used to allocate a huge static
# number of threads. It seems also that when growing the pool it sometimes blocks

# TODO needs a different value on Windows
# TODO try a user temp directory if /var/run/odemisd doesn't exist (and cannot be created)
BASE_DIRECTORY="/var/run/odemisd"
BASE_GROUP="odemis" # user group that is allowed to access the backend


BACKEND_FILE = BASE_DIRECTORY + "/backend.ipc" # the official ipc file for backend (just to detect status)
BACKEND_NAME = "backend" # the official name for the backend container

_microscope = None

def getMicroscope():
    """
    return the microscope component managed by the backend
    Note: if a connection has already been set up, it will reuse it, unless
    you reset _microscope to None
    """
    global _microscope # cached at the module level
    if _microscope is None:
        backend = getContainer(BACKEND_NAME, validate=False)
        _microscope = backend.getRoot()
    return _microscope

def getComponent(name=None, role=None):
    """
    Find a component, according to its name or role.
    At least a name or a role should be provided
    name (str): name of the component to look for
    role (str): role of the component to look for
    return (Component): the component with the given name
    raise LookupError: if no component with such a name is given
    """
    # Note: we could have a "light-weight" version which directly connects to
    # the right container (by-passing the backend), but it's probably going to
    # save time only if just one component is ever used (and immediately found)

    if name is None and role is None:
        raise ValueError("Need to specify at least a name or a role")

    for c in getComponents():
        if name is not None and c.name != name:
            continue
        if role is not None and c.role != role:
            continue
        return c
    else:
        errors = []
        if name is not None:
            errors.append("name %s" % name)
        if role is not None:
            errors.append("role %s" % role)
        raise LookupError("No component with the %s" % (" and ".join(errors),))

def getComponents():
    """
    return (set of Component): all the HwComponents (alive) managed by the backend
    """
    microscope = getMicroscope()
    return microscope.alive.value | {microscope}
    # return _getChildren(microscope)

def _getChildren(root):
    """
    Return the set of components which are referenced from the given component
     (via children)
    root (HwComponent): the component to start from
    returns (set of HwComponents)
    """
    ret = set([root])
    for child in root.children.value:
        ret |= _getChildren(child)

    return ret


#TODO special attributes, which are just properties that are explicitly duplicated
# on the proxy. Getting/setting them always access the actual object remotely.
# declarator is like a property. Two possible implementations:
# * special message types (get/set) instead of method call
# * create special methods on the object, to handle these attributes (when the parent object is registered or shared)

# The special read-only attribute which are duplicated on proxy objects
class roattribute(property):
    """
    A member of an object which will be cached in the proxy when remotely shared.
    It can be modified only before the object is ever shared. (Technically, it
    can still be written afterwards but the values will not be synchronised
    between the containers).
    """
    # the implementation is just a (python) property with only a different name
    # TODO force to not have setter, but I have no idea how to, override __init__?
    pass

def get_roattributes(self):
    """
    list all roattributes of an instance
    Note: this only works on an original class, not on a proxy
    """
#    members = inspect.getmembers(self.__class__)
#    return [name for name, obj in members if isinstance(obj, roattribute)]
    klass = self.__class__
    roattributes = []
    for key in dir(klass):
        try:
            if isinstance(getattr(klass, key), roattribute):
                roattributes.append(key)
        except AttributeError:
            continue

    return roattributes

def dump_roattributes(self):
    """
    list all the roattributes and their value
    """
    # if it is a proxy, use _odemis_roattributes
    roattr = getattr(self, "_odemis_roattributes", [])
    roattr += get_roattributes(self)

    return dict([[name, getattr(self, name)] for name in roattr])

def load_roattributes(self, roattributes):
    """
    duplicate the given roattributes into the instance.
    useful only for a proxy class
    """
    for a, value in roattributes.items():
        setattr(self, a, value)

    # save the list in case we need to pickle the object again
    self._odemis_roattributes = roattributes.keys()


# Container management functions and class

class ContainerObject(Pyro4.core.DaemonObject):
    """Object which represent the daemon for remote access"""

    # it'll never be able to answer back if everything goes fine
    @oneway
    def terminate(self):
        """
        stops the server
        """
        self.daemon.terminate()

    def instantiate(self, klass, kwargs):
        """
        instantiate a component and publish it
        klass (class): component class
        kwargs (dict (str -> value)): arguments for the __init__() of the component
        returns the new component instantiated
        """
        return self.daemon.instantiate(klass, kwargs)

    def getRoot(self):
        """
        returns the root object, if it has been defined in the container
        """
        return self.getObject(self.daemon.rootId)

# Basically a wrapper around the Pyro Daemon
class Container(Pyro4.core.Daemon):
    def __init__(self, name):
        """
        name: name of the container (must be unique)
        """
        assert not "/" in name
        self._name = name
        # all the sockets are in the same directory so it's independent from the PWD
        self.ipc_name = BASE_DIRECTORY + "/" + urllib.quote(name) + ".ipc"

        if not os.path.isdir(BASE_DIRECTORY + "/."): # + "/." to check it's readable
            logging.error("Directory " + BASE_DIRECTORY + " is not accessible, "
                          "which is needed for creating the container %s", name)
        elif os.path.exists(self.ipc_name):
            try:
                os.remove(self.ipc_name)
                logging.warning("The file '%s' was deleted to create container '%s'.", self.ipc_name, name)
            except OSError:
                logging.error("Impossible to delete file '%s', needed to create container '%s'.", self.ipc_name, name)

        Pyro4.Daemon.__init__(self, unixsocket=self.ipc_name, interface=ContainerObject)

        # To be set by the user of the container
        self.rootId = None # objectId of a "Root" component

    def run(self):
        """
        runs and serve the objects registered in the container.
        returns only when .terminate() is called
        """
        # wrapper to requestLoop() just because the name is strange
        self.requestLoop()

    def terminate(self):
        """
        stops the server
        Can be called remotely or locally
        """
        # wrapper to shutdown(), in order to be more consistent with the vocabulary
        self.shutdown()
        # All the cleaning is done in the original thread, after the run()

    def close(self):
        """
        Cleans up everything behind, once the container is already done with running
        Has to be called locally, at the end.
        """
        # unregister every object still around, to be sure everything gets
        # deallocated from the memory (but normally, it's up to the client to
        # terminate() every component before)
        for obj in self.objectsById.values():
            if hasattr(obj, "_unregister"):
                try:
                    obj._unregister()
                except Exception:
                    logging.exception("Failed to unregister object %s when terminating container", str(obj))
            else:
                self.unregister(obj)

        Pyro4.Daemon.close(self)

    def instantiate(self, klass, kwargs):
        """
        instantiate a Component and publish it
        klass (class): component class
        kwargs (dict (str -> value)): arguments for the __init__() of the component
        returns the new component instantiated
        """
        kwargs["daemon"] = self # the component will auto-register
        comp = klass(**kwargs)
        return comp

    def setRoot(self, component):
        """
        sets the root object. It has to be one of the component handled by the
         container.
        component (Component)
        """
        self.rootId = component._pyroId

# helper functions
def getContainer(name, validate=True):
    """
    returns (a proxy to) the container with the given name
    validate (boolean): if the connection should be validated
    raises an exception if no such container exist
    """
    # detect when the base directory doesn't even exists and is readable
    if not os.path.isdir(BASE_DIRECTORY + "/."): # + "/." to check it's readable
        raise IOError("Directory " + BASE_DIRECTORY + " is not accessible.")

    # the container is the default pyro daemon at the address named by the container
    container = Pyro4.Proxy("PYRO:Pyro.Daemon@./u:"+BASE_DIRECTORY+"/"+urllib.quote(name)+".ipc")
    container._pyroTimeout = 120  # s
    container._pyroOneway.add("terminate")

    # A proxy doesn't connect until the first remote call, check the connection
    if validate:
        container.ping() # raise an exception if connection fails
    return container

def getObject(container_name, object_name):
    """
    returns (a proxy to) the object with the given name in the given container
    raises an exception if no such object or container exist
    """
    container = getContainer(container_name, validate=False)
    return container.getObject(urllib.quote(object_name))

def createNewContainer(name, validate=True, in_own_process=True):
    """
    creates a new container in an independent and isolated process
    validate (bool): if the connection should be validated
    in_own_process (bool): if True, creates the container in a separate process
     (so can run fully asynchronously). Otherwise, it is run in a thread. 
    returns the (proxy to the) new container
    """
    # create a container separately
    if in_own_process:
        isready = multiprocessing.Event()
        p = multiprocessing.Process(name="Container " + name, target=_manageContainer,
                                    args=(name, isready))
    else:
        isready = threading.Event()
        p = threading.Thread(name="Container " + name, target=_manageContainer,
                             args=(name, isready))
    p.start()
    if not isready.wait(5): # wait maximum 5s
        logging.error("Container %s is taking too long to get ready", name)
        raise IOError("Container creation timeout")

    # connect to the new container
    return getContainer(name, validate)

def createInNewContainer(container_name, klass, kwargs):
    """
    creates a new component in a new container
    container_name (string)
    klass (class): component class
    kwargs (dict (str -> value)): arguments for the __init__() of the component
    returns:
        (Container) the new container
        (Component) the (proxy to the) new component
    """
    container = createNewContainer(container_name, validate=False)
    try:
        comp = container.instantiate(klass, kwargs)
    except Exception as exp:
        try:
            container.terminate()
        except Exception:
            logging.exception("Failed to stop the container %s after component failure",
                              container_name)
        raise exp
    return container, comp

def _manageContainer(name, isready=None):
    """
    manages the whole life of a container, from birth till death
    name (string)
    isready (Event): set when the container is (almost) ready to publish objects
    """
    container = Container(name)
    # TODO: also change the process name/arguments to easily known which process
    # is what? cf py-setproctitle
    logging.debug("Container %s runs in PID %d", name, os.getpid())
    if isready is not None:
        isready.set()
    container.run()
    container.close()

# Special functions and class to manage method/function with weakref
# wxpython.pubsub has something similar

class WeakRefLostError(Exception):
    pass

class WeakMethodBound(object):
    def __init__(self, f):
        self.f = f.__func__
        self.c = weakref.ref(f.__self__)
        # cache the hash so that it's the same after deref'd
        self.hash = hash(f.__func__) + hash(f.__self__)

    def __call__(self, *arg, **kwargs):
        ins = self.c()
        if ins == None:
            raise WeakRefLostError('Method called on dead object')
        return self.f(ins, *arg, **kwargs)

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        try:
            return (type(self) is type(other) and self.f == other.f
                    and self.c() == other.c())
        except:
            return False

    # def __ne__(self, other):
    #     return not self == other

class WeakMethodFree(object):
    def __init__(self, f):
        self.f = weakref.ref(f)
        # cache the hash so that it's the same after deref'd
        self.hash = hash(f)

    def __call__(self, *arg, **kwargs):
        fun = self.f()
        if fun == None:
            raise WeakRefLostError('Function no longer exist')
        return fun(*arg, **kwargs)

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        try:
            return type(self) is type(other) and self.f() == other.f()
        except:
            return False

    # def __ne__(self, other):
    #    return not self == other

def WeakMethod(f):
    try:
        # Check if the paramater has a function object, which is the case
        # if it's a bound function (ie.e a method)
        f.__func__
    except AttributeError:
        return WeakMethodFree(f)
    return WeakMethodBound(f)

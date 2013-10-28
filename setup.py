#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob
import os
import subprocess
import sys

# To be updated to the current version
VERSION = "1.3"
# We cannot use the git version because it's not (always) available when building
# the debian package

# Trick from openshot
# Boolean: running as root?
ROOT = os.geteuid() == 0
# For Debian packaging it could be a fakeroot so reset flag to prevent execution of
# system update services for Mime and Desktop registrations.
# The debian/odemis.postinst script must do those.
if not os.getenv("FAKEROOTKEY") == None:
    print "NOTICE: Detected execution in a FakeRoot so disabling calls to system update services."
    ROOT = False

# almost copy from odemis.__init__.py, but we cannot load it as it's not installed yet
def _get_version_git():
    """
    Get the version via git
    raises LookupError if no version info found
    """
    # change directory to root
    rootdir = os.path.dirname(__file__) # .

#    if not os.path.isdir(rootdir) or not os.path.isdir(os.path.join(rootdir, ".git")):
#        raise LookupError("Not in a git directory")

    try:
        out = subprocess.check_output(args=["git", "describe", "--tags", "--dirty", "--always"],
                                      cwd=rootdir)

        return out.strip()
    except EnvironmentError:
        raise LookupError("Unable to run git")

# Check version
try:
    gver = _get_version_git()
    if "-" in gver:
            sys.stderr.write("Warning: packaging a non-tagged version: %s\n" % gver)
    if VERSION != gver:
        sys.stderr.write("Warning: package version and git version don't match:"
                         " %s <> %s\n" % (VERSION, gver))
except LookupError:
    pass


if sys.platform.startswith('linux'):
    data_files = [('/etc/', ['install/linux/etc/odemis.conf']),
                  # Not copying sudoers file, as we are not sure there is a sudoers.d directory
                  # TODO udev rules might actually be better off in /lib/udev/rules.d/
                  ('/etc/udev/rules.d', glob.glob('install/linux/etc/udev/rules.d/*.rules')),
                  ('share/odemis/', glob.glob('install/linux/usr/share/odemis/*.odm.yaml')),
                  # TODO: need to run desktop-file-install in addition to update-desktop-database?
                  ('share/applications/', ['install/linux/usr/share/applications/odemis.desktop']),
                  ('share/icons/hicolor/32x32/apps/', ['install/linux/usr/share/icons/hicolor/32x32/apps/odemis.png']),
                  ('share/icons/hicolor/64x64/apps/', ['install/linux/usr/share/icons/hicolor/64x64/apps/odemis.png']),
                  ('share/icons/hicolor/128x128/apps/', ['install/linux/usr/share/icons/hicolor/128x128/apps/odemis.png']),
                  ('share/doc/odemis/', glob.glob('doc/*.txt')),
                  ('share/doc/odemis/scripts/', glob.glob('scripts/*.py')),
                  ]
    # TODO: see if we could use entry_points instead
    scripts = ['install/linux/usr/local/bin/odemisd',
               'install/linux/usr/local/bin/odemis-cli',
               'install/linux/usr/local/bin/odemis-convert',
               'install/linux/usr/local/bin/odemis-gui',
               'install/linux/usr/local/bin/odemis-start',
               'install/linux/usr/local/bin/odemis-stop'
               ]
else:
    data_files = []
    scripts = []
    sys.stderr.write("Warning: Platform %s not supported" % sys.platform)

dist = setup(name='Odemis',
             version=VERSION,
             description='Open Delmic Microscope Software',
             author='Éric Piel, Rinze de Laat',
             author_email='piel@delmic.com, laat@delmic.com',
             url='https://github.com/delmic/odemis',
             classifiers=["License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                          "Operating System :: POSIX :: Linux",
                          "Programming Language :: Python",
                          "Intended Audience :: Science/Research",
                          "Topic :: Scientific/Engineering",
                          "Environment :: Console",
                          "Environment :: X11 Applications :: GTK",
                         ],
             package_dir={'': 'src'},
             packages=find_packages('src', exclude=["*.test"]),
             package_data={'odemis.gui.img': ["example/*.png", "example/*.mat",
                                              "calibration/*.png"]
                          },
             scripts=scripts,
             data_files=data_files # not officially in setuptools, but works as for distutils
            )

if ROOT and dist != None:
    # for mime file association, see openshot's setup.py
    # update the XDG .desktop file database
    try:
        print "Updating the .desktop file database."
        subprocess.check_output(["update-desktop-database"])
    except Exception:
        sys.stderr.write("Failed to update.\n")

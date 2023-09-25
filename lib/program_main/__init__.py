# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os
import sys
import shutil
import re
import json
import datetime

# PyMOL.
from pymol import cmd

from pymol.Qt import QtWidgets, QtCore, QtGui

import os
import shutil

from .program_gui.main_window import PCN_Miner_main_window_qt

def program_launcher(app, plugin_name, program_version, program_revision):
    pcn_miner_plugin = PCN_Miner(app, plugin_name, program_version, program_revision)


class PCN_Miner():

    def __init__(self, app, plugin_name, program_version, program_revision):

        self.plugin_name = plugin_name
        self.program_version = program_version
        self.program_revision = program_revision

        # Set to 'True' when developing, useful for debugging.
        self.DEVELOP = True
        # Set to 'True' to perform some tests on sequences/structures from the GUI.
        self.TEST = False

        self.app = app

         # If set to 'True' the most time consuming protocols will be run in a thread so
        # that the GUI is not freezed. When developing the code, it is better to set it to 'False',
        # in order to better track exceptions.
        if self.DEVELOP:
            self.use_protocol_threads = False
        else:
            self.use_protocol_threads = True


      #### MAIN WINDOW OF THE PLUGIN ####

        self.main_window = PCN_Miner_main_window_qt(self)

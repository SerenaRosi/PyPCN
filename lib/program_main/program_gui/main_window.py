# Copyright 2023 by Serena Rosignoli. All rights reserved.
# This code is part of the PyPCN package and governed by its license. Please
# see the LICENSE file.


import os
import sys
import shutil
import re
import json
import datetime
import io as io
# PyMOL.
import pymol
from pymol import cmd

from pymol.Qt import QtWidgets, QtCore, QtGui

# Import tabs
from .tabs import CHILD_TAB, INPUTStab_widgets, CENTRALITYtab_widgets, SPECTRALtab_widgets, EMBEDDINGStab_widgets, COMMUNITYtab_widgets, RESULTStab_widgets, OTHERtab_widgets

import warnings
import math


class PCN_Miner_main_window_main_menu:

    is_program_main_window = True

    def make_main_menu(self):

        """
        A method to create the Main Window Main Menu.

        note: it is currently not used.
        """

        self.menubar = self.menuBar()
        self.menubar.setNativeMenuBar(False)

        #---------------
        # "File" menu. -
        #---------------

        self.file_menu = self.menubar.addMenu('File')

        # Workspace submenu.
        self.sessions_submenu = QtWidgets.QMenu('Sessions', self)
        self.file_menu.addMenu(self.sessions_submenu)

        self.file_menu.addSeparator()
        self.exit_submenu = QtWidgets.QMenu('Exit', self)
        self.file_menu.addMenu(self.exit_submenu)



class PCN_Miner_main_window_qt(QtWidgets.QMainWindow, PCN_Miner_main_window_main_menu):

    is_program_main_window = True

    def __init__(self, program, parent=None):
        super(PCN_Miner_main_window_qt, self).__init__(parent)

        """

        MAIN WINDOW

        """

        # Initial settings.
        self.program = program
        self.title = self.program.plugin_name
        self.statusBar_message = "Welcome to PCN Miner"

        self.left = 550
        self.top = 50
        self.width = 550
        self.height = 400

        self.vertical_spacing = 1

        # Creating a menu bar.
        #self.make_main_menu()

        ####
        # Main Window Widget
        ####
        self.main_window_widget = Main_Window_Widget(self)

        ####
        # Central Widget
        ####
        self.central_widget = Centralwid(self)
        self.setCentralWidget(self.central_widget)

        # Add the Main Window Widget to the Central Widget and set the layout
        self.central_widget.central_layout.addWidget(self.main_window_widget.widget)
        self.main_window_widget.widget.setLayout(self.main_window_widget.main_widget_layout)

        # Set the layout of the Central Widget.
        self.central_widget.setLayout(self.central_widget.central_layout)
        self.central_widget.central_layout.setFormAlignment(QtCore.Qt.AlignLeft)
        self.central_widget.central_layout.setVerticalSpacing(self.vertical_spacing)

        # Creating status bar.
        self.statusBar().showMessage(self.statusBar_message, 3000)
        self.statusBar().setSizeGripEnabled(1)

        # Initialize User Interface.
        self.initUI()

        self.check_dependencies()


    def check_dependencies(self):

        # TODO: FUNCTION TO CHECK IF PIP INSTALL IS FUNCTIONAL

        from contextlib import redirect_stderr

        # If a simple command like 'pip config list' raise Exception, the PyMOL version in use is not compatible with in-PyMOL installation,
        # (Likely is an Open-source PyMOl)
        # with io.StringIO() as stderr, redirect_stderr(stderr):
        #     r = cmd.do("pip config list")
        #     self.s = stderr.getvalue()
        #
        # if re.search("SyntaxError", self.s):
        #     notfunctional = True

        list_of_dependencies = ["configparser", "scipy", "fcmeans", "gem", "cdlib", "pytz", "node2vec", "leidenalg", "sklearn", "matplotlib", "infomap", "wurlitzer",
        "karateclub", "ASLPAw", "Graph", "pandas"]

        dict_dep = {}
        dict_dep["configparser"] = "pip install configparser"
        dict_dep["scipy"] = "pip install scipy"
        dict_dep["regex"] = "pip install regex"
        dict_dep["fcmeans"] = "pip install fuzzy-c-means"

        dict_dep["gem"] = "pip install git+https://github.com/palash1992/GEM.git"
        dict_dep["cdlib"] = "pip install cdlib"
        dict_dep["pytz"] = "pip install pytz"

        dict_dep["node2vec"] = "pip install node2vec"
        dict_dep["leidenalg"] = "pip install leidenalg"
        dict_dep["sklearn"] = "pip install sklearn"

        dict_dep["matplotlib"] = "pip install matplotlib"
        dict_dep["infomap"] = "pip install infomap"
        dict_dep["wurlitzer"] = "pip install wurlitzer"
        dict_dep["karateclub"] = "pip install karateclub"
        dict_dep["ASLPAw"] = "pip install ASLPAw"
        dict_dep["Graph"] = "pip install graph-tools"
        dict_dep["pandas"] = "pip install pandas"

        list_of_missing_dep = []

        for dep in list_of_dependencies:
            if dep == "ASLPAw":
                try:
                    from ASLPAw_package import ASLPAw
                except:
                    print(dep + " not found")
                    list_of_missing_dep.append([dep, dict_dep[dep]])
            elif dep == "Graph":
                try:
                    from graph_tools import Graph
                except:
                    print(dep + " not found")
                    list_of_missing_dep.append([dep, dict_dep[dep]])
            else:
                try:
                    __import__(dep)
                except:
                    #print(dep + " not found")
                    list_of_missing_dep.append([dep, dict_dep[dep]])

        if list_of_missing_dep:
            string = "\n"
            for dep in list_of_missing_dep:
                string = string + ': '.join(dep) + "\n"

            message = "To use all features of PyPCN the next missing dependencies must be installed.\n\nHere are some suggestions with 'pip install':\n" + string

            button_1_text = "Install Automatically"
            button_2_text = "Install Manually"

            choice = self.installwindow("Missing dependencies",
                                 message = message,
                                 parent = self,
                                 buttons_text = [button_1_text, button_2_text])

            if choice:
                q = QtWidgets.QMessageBox.question(self, "Automatic Installation Warning", "The automatic installation will proceed with the installation of all the missing dependencies in the PyMOL environment.\nThis may take some time. Please be patient.\n\nDo you really want to continue?")
                if q == QtWidgets.QMessageBox.Yes:
                    for dep in list_of_missing_dep:
                        cmd.do(dep[1])

                    QtWidgets.QMessageBox.warning(self, "Automatic Installation Completed", "Please restart PyMOL and check if all dependecies have been installed.\n\nThe automatic installation may fail because:\n- your PyMOL setup does not support it (e.g. Open-Source or Incentive versions < 2.5)\n- only some dependencies are failing to be installed due to not widely known reasons (check output messages on PyMOL prompt).\n\nFor further information, please check requirements in the PyPCN User's Guide.")

            else:
                QtWidgets.QMessageBox.warning(self, "Manual Installation Warning", "Please read the instructions in the PyPCN User's Guide at: https://github.com/pcnproject/PyPCN/releases/download/utilities/Supplementary-PyPCN_User_Guide.pdf")

# fcmeans: pip install fuzzy-c-means
# gem: pip install git+https://github.com/palash1992/GEM.git
# cdlib: pip install cdlib
# node2vec: pip install node2vec
# leidenalg: pip install leidenalg
# infomap: pip install infomap
# wurlitzer: pip install wurlitzer
# karateclub: pip install karateclub
# ASLPAw: pip install ASLPAw
# Graph: pip install graph-tools

    def installwindow(self, title, message, parent=None, buttons_text=None):

        """
        Wrapper to a Yes/no dialog in PyQt. If 'buttons_text' is 'None', the default
        "Yes" and "No" buttons will be used. If 'buttons_text' is a list with two
        strings, the first string will be the text of the "Yes" button and the second
        one will be the text of the "No" button.
        """

        dialog = QtWidgets.QMessageBox(parent)
        dialog.setWindowTitle(title)
        dialog.setText(message)

        count = dialog.layout().count()

        # for i in range(count):
        #     item = dialog.layout().itemAt(i)
        #     print(item)

        if len(buttons_text) == 1:
            start_button = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.AcceptRole)
            dialog.setStyleSheet('QLabel { font-size: 13pt; padding: 15px}')
            answer = dialog.exec_()

        else:

            yesbutton = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.YesRole)
            nobutton = dialog.addButton(buttons_text[1], QtWidgets.QMessageBox.NoRole)
            dialog.setStyleSheet('QLabel { font-size: 13pt; padding: 15px}')
            yesbutton.setStyleSheet('QPushButton {font-weight: bold; border-color: green}')

            answer = dialog.exec_()
            return dialog.clickedButton() is yesbutton



    def set_pymol_visualization_options(self):

        pymol_version = float(".".join(cmd.get_version()[0].split(".")[0:2]))

        if pymol_version >= 2.5:
            try:
                cmd.undo_disable()
            except:
                pass

        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_highlight_color", "grey50")
        cmd.set("sphere_scale", 0.15)
        cmd.set("orthoscopic", "on")
        cmd.set("cartoon_dumbbell_length", 1.5)
        cmd.set("cartoon_dumbbell_width", 0.4)
        cmd.set("cartoon_dumbbell_radius", 0.3)
        cmd.set("cartoon_ring_mode", 1)
        cmd.set("cartoon_ring_transparency", 0.5)

    def initUI(self):

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # Loads and sets the Qt stylesheet.
        module_path = sys.modules[__name__].__file__
        self.qss = QSSHelper.open_qss(os.path.join(os.path.dirname(module_path), 'aqua', 'aqua.qss'))
        self.setStyleSheet(self.qss)
        self.show()



class QSSHelper:
    def __init__(self):
        pass

    @staticmethod
    def open_qss(path):
        """
        opens a Qt stylesheet with a path relative to the project

        Note: it changes the urls in the Qt stylesheet (in memory), and makes these urls relative to the project
        Warning: the urls in the Qt stylesheet should have the forward slash ('/') as the pathname separator
        """
        with open(path) as f:
            qss = f.read()
            pattern = r'url\((.*?)\);'
            for url in sorted(set(re.findall(pattern, qss)), key=len, reverse=True):
                directory, basename = os.path.split(path)
                new_url = os.path.join(directory, *url.split('/'))
                new_url = os.path.normpath(new_url)
                new_url = new_url.replace(os.path.sep, '/')
                qss = qss.replace(url, new_url)
            return qss



class Centralwid(QtWidgets.QWidget):

    """
    A class to reppresent the Central Widget of the Main Window
    """

    def __init__(self, main_window):
        super(Centralwid, self).__init__(main_window)
        self.style = "background-color: rgb(0, 0, 0); color: rgb(255, 255, 255); font-weight: bold"
        self.main_window = main_window
        self.initUI()

    def initUI(self):

        self.central_layout = QtWidgets.QFormLayout()


class Main_Window_Widget(QtWidgets.QWidget):

    """
    Main Window Widget of the plugin.
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window

        self.create_main_window_widget()

    def create_main_window_widget(self):

        # Customize PyMOl visualization.
        self.main_window.set_pymol_visualization_options()

        # if os.system == "darwin":
        #     self.USE_THREADS = False
        # else:
        #     self.USE_THREADS = True

        # Path where the plugin is located
        self.current_path = (os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
        self.pcn_main_path = os.path.join(self.current_path, "program_scripts", "pcn", "pcn_main.py")
        self.tools_path = os.path.join(self.current_path, "program_scripts", "pcn", "tools")
        self.config_path = os.path.join(self.tools_path, "config.ini")
        self.working_dir_path = ""

        self.dict_of_algorithms = {}
        self.dict_of_algorithms["centrality"] = ["closeness", "betweenness", "eigenvector_c", "degree_c"]
        self.dict_of_algorithms["spectral"] = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
        self.dict_of_algorithms["embedd"] = ["fuzzycmeans_hope", "kmeans_hope", "fuzzycmeans_laplacianeigenmaps", "kmeans_laplacianeigenmaps", "fuzzycmeans_node2vec", "kmeans_node2vec"]
        self.dict_of_algorithms["community"] = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]

        # Dict of results for PyMOL sessions
        self.algorithms_results_dict = {}

        for key in self.dict_of_algorithms:
            for alg in self.dict_of_algorithms[key]:
                self.algorithms_results_dict[alg] = {}
                self.algorithms_results_dict[alg]["results"] = []
                self.algorithms_results_dict[alg]["location"] = ""

        '''
        # Dict of results for txt files
        {'closeness': {'results': ['5xh6_closeness.txt', '7o4g_closeness.txt', '3uoh_closeness.txt'],
                      'location': '/home/ale/Scrivania/pcn miner/outputCentralities/closeness/Txt'},
        'betweenness': {'results': [], 'location': ''}, 'eigenvector_c': {'results': [],
                      'location': ''}, 'degree_c': {'results': [], 'location': ''},
        'unnorm_ssc': {'results': [], 'location': ''},
        'norm_ssc': {'results': [], 'location': ''},
        'unnorm_hsc': {'results': [], 'location': ''},
        '''
        self.algorithms_results_dict_txt = {}

        for key in self.dict_of_algorithms:
            for alg in self.dict_of_algorithms[key]:
                self.algorithms_results_dict_txt[alg] = {}
                self.algorithms_results_dict_txt[alg]["results"] = []
                self.algorithms_results_dict_txt[alg]["location"] = ""

        '''
        # Dict of results for participation coeff files

        '''
        self.part_coeff_results_dict = {}

        for key in self.dict_of_algorithms:
            for alg in self.dict_of_algorithms[key]:
                self.part_coeff_results_dict[alg] = {}
                self.part_coeff_results_dict[alg]["results"] = []
                self.part_coeff_results_dict[alg]["location"] = ""

        self.algorithms_results_dict_clusters = {}

        for key in self.dict_of_algorithms:
            for alg in self.dict_of_algorithms[key]:
                self.algorithms_results_dict_clusters[alg] = {}
                self.algorithms_results_dict_clusters[alg]["results"] = []
                self.algorithms_results_dict_clusters[alg]["location"] = ""

        #
        self.algorithms_results_dict_summary = {}

        for key in self.dict_of_algorithms:
            for alg in self.dict_of_algorithms[key]:
                self.algorithms_results_dict_summary[alg] = {}
                self.algorithms_results_dict_summary[alg]["results"] = []
                self.algorithms_results_dict_summary[alg]["location"] = ""

        # Separators
        if sys.platform == "win32":
            self.path_sep = "\\"

        elif sys.platform == "linux":
            self.path_sep = "/"

        elif sys.platform == "darwin":
            self.path_sep = "/"

        # Main Window widget
        self.widget = QtWidgets.QWidget()

        # Main Window layout
        self.main_widget_layout = QtWidgets.QGridLayout()

        # Tab Widget
        self.TABS = QtWidgets.QTabWidget()

        # INPUTS tab and layout
        self.INPUTS = QtWidgets.QWidget()
        self.INPUTS_layout = CHILD_TAB(self).tab_layout

        # SPECTRAL tab and layout
        self.SPECTRAL = QtWidgets.QWidget()
        self.SPECTRAL_layout = CHILD_TAB(self).tab_layout

        # EMEDDINGS tab and layout
        self.EMBEDDINGS = QtWidgets.QWidget()
        self.EMBEDDINGS_layout = CHILD_TAB(self).tab_layout

        # COMMUNITY tab and layout
        self.COMMUNITY = QtWidgets.QWidget()
        self.COMMUNITY_layout = CHILD_TAB(self).tab_layout

        # CENTRALITY tab and Layout
        self.CENTRALITY = QtWidgets.QWidget()
        self.CENTRALITY_layout = CHILD_TAB(self).tab_layout

        # CENTRALITY tab and Layout
        self.RESULTS = QtWidgets.QWidget()
        self.RESULTS_layout = CHILD_TAB(self).tab_layout

        # CENTRALITY tab and Layout
        self.OTHER = QtWidgets.QWidget()
        self.OTHER_layout = CHILD_TAB(self).tab_layout

        # Add child tabs to TABS and set Layouts
        self.TABS.addTab(self.INPUTS, "INPUTS")
        self.TABS.addTab(self.CENTRALITY, "CENTRALITY ANALYSIS")
        self.TABS.addTab(self.SPECTRAL, "SPECTRAL CLUSTERING")
        self.TABS.addTab(self.EMBEDDINGS, "EMBEDDINGS + CLUSTERING")
        self.TABS.addTab(self.COMMUNITY, "COMMUNITY EXTRACTION")
        self.TABS.addTab(self.RESULTS, "RESULTS")
        self.TABS.addTab(self.OTHER, "Other")

        self.INPUTS.setLayout(self.INPUTS_layout)
        self.CENTRALITY.setLayout(self.CENTRALITY_layout)
        self.SPECTRAL.setLayout(self.SPECTRAL_layout)
        self.EMBEDDINGS.setLayout(self.EMBEDDINGS_layout)
        self.COMMUNITY.setLayout(self.COMMUNITY_layout)
        self.RESULTS.setLayout(self.RESULTS_layout)
        self.OTHER.setLayout(self.OTHER_layout)

        self.INPUTS_widgets = INPUTStab_widgets(self, self.INPUTS_layout)
        self.CENTRALITY_widgets = CENTRALITYtab_widgets(self, self.CENTRALITY_layout)
        self.SPECTRAL_widgets = SPECTRALtab_widgets(self, self.SPECTRAL_layout)
        self.EMBEDDINGS_widgets = EMBEDDINGStab_widgets(self, self.EMBEDDINGS_layout)
        self.COMMUNITY_widgets = COMMUNITYtab_widgets(self, self.COMMUNITY_layout)
        self.RESULTS_widgets = RESULTStab_widgets(self, self.RESULTS_layout)
        self.OTHER_widgets = OTHERtab_widgets(self, self.OTHER_layout)

        self.main_widget_layout.addWidget(self.TABS)

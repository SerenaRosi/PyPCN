import os
import sys
import subprocess
import shutil
from pathlib import Path
import numpy as np
import warnings
import json
import urllib.request
import json


# PyMOL.
import pymol
from pymol import cmd

from pymol.Qt import QtWidgets, QtCore, QtGui

import lib.program_main.program_scripts.pcn.pcn_main as pcn_module
from lib.program_main.program_scripts.threads import Protocol_exec_dialog
from lib.program_main.program_gui.frames import Frame
from lib.program_main.program_gui.cgo_arrow import *
from lib.program_main.program_gui.contact_map_visualization import *
from lib.program_main.program_gui.plots import *

try:
    from networkx import from_numpy_matrix
except:
    pass

def get_adj_pdb_choice(cb):

    if cb.isChecked():
        return "adj"
    else:
        return "pdb"


def unselect_all(list, all_box):

    for boxes in list:
        if not boxes.isChecked():
            all_box.setChecked(False)

def click_all(list, all_box):

    for boxes in list:
        if all_box.isChecked():
            boxes.setChecked(True)
        else:
            boxes.setChecked(False)


def check_inputs(main, pdb_line_edit, working_dir_line_edit):

    input_correct = True

    if main.INPUTS_widgets.use_precomputed_matrix_cb.isChecked():
        if not main.INPUTS_widgets.precomputed_matrix_linedit.text():
            QtWidgets.QMessageBox.warning(main, "Warning", "'Use pre-computed adjacency matrix' option is selected\nPlease select an adjacency matrix file or uncheck the option to continue")
            input_correct = False

    if main.INPUTS_widgets.use_pymol_protein_cb.isChecked():
        pymol_objs_list = cmd.get_object_list()
        if not pdb_line_edit.text() in pymol_objs_list:
            QtWidgets.QMessageBox.warning(main, "Warning", pdb_line_edit.text() + " not found in PyMOL workspace\nPlease insert a valid PyMOL object name")
            input_correct = False

    # Get the text from the configuration_line_edit
    # Control if it's a path leading to a folder, otherwise give a warning message
    input_line = working_dir_line_edit.text()
    if not os.path.isdir(input_line):
        QtWidgets.QMessageBox.warning(main, "Warning", "Please choose the working directory")
        input_correct = False

    # Get the text from the pdb_line_edit
    # If the pdb code is not present, give a warning message
    if not pdb_line_edit.text():
        QtWidgets.QMessageBox.warning(main, "Warning", "There aren't enough parameters\n Please enter the PDB ID")
        input_correct = False

    return input_correct


class CHILD_TAB(QtWidgets.QWidget):

    """
    Layout for INPUTS tab.
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window

        self.initUI()

    def initUI(self):

        # Layout
        self.tab_layout = QtWidgets.QGridLayout()

######################## MODIFIED ##########################
class INPUTStab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        # Create two Group Boxes
        self.config_box = QtWidgets.QGroupBox("Configuration")
        self.input_parameters_box = QtWidgets.QGroupBox("Input parameters")

        # Add the Group Boxes to the INPUTS tab Layout
        self.layout.addWidget(self.config_box)
        self.layout.addWidget(self.input_parameters_box)

        # Create and set Group Boxes layouts
        self.config_box_layout = QtWidgets.QGridLayout()
        self.input_parameters_box_layout = QtWidgets.QGridLayout()
        self.config_box.setLayout(self.config_box_layout)
        self.input_parameters_box.setLayout(self.input_parameters_box_layout)

        # Create and add input, output and adj to the configuration box
        self.input_label = QtWidgets.QLabel("Working Directory:")
        self.config_box_layout.addWidget(self.input_label, 0, 0)

        # self.output_label = QtWidgets.QLabel("OUTPUT:")
        # self.config_box_layout.addWidget(self.output_label, 1, 0)
        #
        # self.adj_label = QtWidgets.QLabel("ADJ:")
        # self.config_box_layout.addWidget(self.adj_label, 2, 0)

        # Create the configuration QLineEdit Widget and disable it
        self.working_dir_line_edit = QtWidgets.QLineEdit("")
        # self.configuration_line1_edit = QtWidgets.QLineEdit("")
        # self.configuration_line2_edit = QtWidgets.QLineEdit("")
        # self.configuration_line3_edit = QtWidgets.QLineEdit("")

        self.working_dir_line_edit.setEnabled(False)
        # self.configuration_line1_edit.setEnabled(False)
        # self.configuration_line2_edit.setEnabled(False)
        # self.configuration_line3_edit.setEnabled(False)

        # A function to control if the path is present in the configuration_line_edit
        self.read_config(self.working_dir_line_edit)
        # self.read_config(self.configuration_line1_edit, "input")
        # self.read_config(self.configuration_line2_edit, "output")
        # self.read_config(self.configuration_line3_edit, "adj")

        # Add the configuration QLineEdit Widget to the configuration box
        self.config_box_layout.addWidget(self.working_dir_line_edit, 0, 1)
        #self.config_box_layout.addWidget(self.configuration_line1_edit, 0, 1)
        # self.config_box_layout.addWidget(self.configuration_line2_edit, 1, 1)
        # self.config_box_layout.addWidget(self.configuration_line3_edit, 2, 1)

        # Create and add open .txt file button
        self.working_dir_button = QtWidgets.QPushButton("Choose directory ...")
        #self.clear_dir_button = QtWidgets.QPushButton("Clear Directory")
        #self.clear_dir_button.setStyleSheet("border-color: red")
        # self.open_file1_button = QtWidgets.QPushButton("Choose directory ...")
        # self.open_file2_button = QtWidgets.QPushButton("Choose directory ...")
        # self.open_file3_button = QtWidgets.QPushButton("Choose directory ...")

        self.config_box_layout.addWidget(self.working_dir_button, 0, 2)
        #self.config_box_layout.addWidget(self.clear_dir_button, 0, 3)
        # self.config_box_layout.addWidget(self.open_file1_button, 0, 2)
        # self.config_box_layout.addWidget(self.open_file2_button, 1, 2)
        # self.config_box_layout.addWidget(self.open_file3_button, 2, 2)

        self.working_dir_button.clicked.connect(lambda: self.open_file_button_func(self.working_dir_line_edit))
        # self.open_file2_button.clicked.connect(lambda: self.open_file_button_func(self.configuration_line2_edit, "output"))
        # self.open_file3_button.clicked.connect(lambda: self.open_file_button_func(self.configuration_line3_edit, "adj"))

        # Create and add pdb id to the input parameters box
        self.use_precomputed_matrix_cb = QtWidgets.QCheckBox("Use pre-computed adjacency matrix")
        self.input_parameters_box_layout.addWidget(self.use_precomputed_matrix_cb, 0, 0)
        self.use_precomputed_matrix_cb.clicked.connect(self.show_precomputed_options)

        self.precomputed_matrix_linedit = QtWidgets.QLineEdit()
        self.input_parameters_box_layout.addWidget(self.precomputed_matrix_linedit, 0, 1)
        self.precomputed_matrix_linedit.setEnabled(False)

        self.precomputed_matrix_button = QtWidgets.QPushButton("Choose file ...")
        self.input_parameters_box_layout.addWidget(self.precomputed_matrix_button, 0, 2)
        self.precomputed_matrix_button.setEnabled(False)
        self.precomputed_matrix_button.clicked.connect(lambda: self.open_adj_file_func(self.precomputed_matrix_linedit))

        self.use_pymol_protein_cb = QtWidgets.QCheckBox("Use protein loaded in PyMOL")
        self.input_parameters_box_layout.addWidget(self.use_pymol_protein_cb, 1, 0)

        self.pdb_id_label = QtWidgets.QLabel("PyMOL object name or PDB-ID:")
        self.input_parameters_box_layout.addWidget(self.pdb_id_label, 2, 0)

        self.pdb_line_edit = QtWidgets.QLineEdit()
        self.pdb_line_edit.setPlaceholderText("e.g. 1bj4")
        self.input_parameters_box_layout.addWidget(self.pdb_line_edit, 2, 1)

        # Create and add the covalent bonds label to the parameters box
        self.non_covalent_label = QtWidgets.QLabel("Min threshold: ")
        self.input_parameters_box_layout.addWidget(self.non_covalent_label, 3, 0)

        # Create and add a non covalent QSpinBox to the parameters box to select a non covalent bonds threshold distance
        self.non_covalent_box = QtWidgets.QSpinBox()
        self.non_covalent_box.setRange(0, 1000)
        self.non_covalent_box.setSingleStep(1)
        self.non_covalent_box.setValue(4)
        self.input_parameters_box_layout.addWidget(self.non_covalent_box, 3, 1)

        # Create and add the only significant bonds label to the parameters box
        self.only_significant_label = QtWidgets.QLabel("Max threshold: ")
        self.input_parameters_box_layout.addWidget(self.only_significant_label, 4, 0)

        # Create and add a significant bonds QSpinBox to the parameters box to select an only significant bonds threshold distance
        self.only_significant_box = QtWidgets.QSpinBox()
        self.only_significant_box.setRange(0, 1000)
        self.only_significant_box.setSingleStep(1)
        self.only_significant_box.setValue(8)
        self.input_parameters_box_layout.addWidget(self.only_significant_box, 4, 1)

        # Create ComboBox to select the type of adjacency matrix to compute
        self.adj_mat_type_cb_label = QtWidgets.QLabel("Adjacency Matrix Type: ")
        self.input_parameters_box_layout.addWidget(self.adj_mat_type_cb_label, 5, 0)

        self.adj_mat_type_cb = QtWidgets.QComboBox()
        self.adj_mat_type_cb.addItems(["alpha-Carbons", "beta-Carbons", "Centroids"])
        self.input_parameters_box_layout.addWidget(self.adj_mat_type_cb, 5, 1)

        self.USE_THREADS = QtWidgets.QCheckBox("Use Threads")
        self.input_parameters_box_layout.addWidget(self.USE_THREADS, 6, 0)
        if sys.platform == "darwin":
            self.USE_THREADS.setChecked(False)
        else:
            self.USE_THREADS.setChecked(True)

        self.USE_THREADS.setToolTip('Able/Disable the use of threads.')

        # Create and add the coef plot label to the parameters box
        self.coef_plot_label = QtWidgets.QLabel("Compute the participation coefficient plot:")
        #self.input_parameters_box_layout.addWidget(self.coef_plot_label, 4, 0)

        # Create and add a coef plot QCheckBox to the parameters box to compute or not a participation coef plot
        self.coef_plot_box = QtWidgets.QCheckBox()
        self.coef_plot_box.setChecked(True)
        #self.input_parameters_box_layout.addWidget(self.coef_plot_box, 4, 1)

        # Create a Next Page button, add it to the layout and connect it with the "next_page_func"
        self.next_page_button = QtWidgets.QPushButton("Next")
        self.layout.addWidget(self.next_page_button)
        self.next_page_button.clicked.connect(self.next_page_func)

    def show_precomputed_options(self):

        if self.use_precomputed_matrix_cb.isChecked():
            self.precomputed_matrix_button.setEnabled(True)
        else:
            self.precomputed_matrix_button.setEnabled(False)


    def read_config(self, line_edit):

        # Find the path to the config.ini file
        path_to_config = os.path.join(self.main_window.tools_path, "config.ini")

        # Open the config.ini file
        config_ini = open(path_to_config, "rt")

        # # Assign a variable to the path name
        # if label == "input":
        #     user_path = "proteins_path ="
        # elif label == "output":
        #     user_path = "output_path ="
        # elif label == "adj":
        #     user_path = "adj_filespath ="
        # elif label == "working_dir":
        #     user_path = "working_dir_path ="

        # Read the config.ini file line by line
        for line in config_ini:
            # Control if the line starts with the path name
            if line.startswith("working_dir_path ="):
                # Convert the line into a string
                line_string = str(line)
                # Split the string using "=" as a divider, take what is after the "=" and remove "\n"
                string_to_line_edit = (line_string.split("=")[1]).replace("\n", "")
                # If the config.ini file is empty
                if not string_to_line_edit:
                    # Replace "" with "Not Found" in the configuration_line_edit
                    line_edit.setText("Not Found")
                # If the config.ini file is filled in
                else:
                    # Check if the adj, input and output directories exist
                    if os.path.isdir(os.path.join(string_to_line_edit, "adj")) and os.path.isdir(os.path.join(string_to_line_edit, "input")) and os.path.isdir(os.path.join(string_to_line_edit, "output")):
                        # Take the path present after "=" and use it to replace "Not Found" in the configuration_line_edit
                        line_edit.setText(string_to_line_edit)
                        self.main_window.working_dir_path = string_to_line_edit
                    else:
                        line_edit.setText("Not Found")

    def open_adj_file_func(self, line_edit):

        self.file_path = QtWidgets.QFileDialog.getOpenFileName(self, "Open Text File", "","text files (*.txt)")
        if self.file_path:
            line_edit.setText(self.file_path[0])
            QtWidgets.QMessageBox.warning(self.main_window, "Attention!", "Please select the Min and Max thresholds according to those used to compute the selected adjacent matrix.")
        else:
            pass

    def open_file_button_func(self, line_edit):

        ## File dialog to select a directory
        self.file_path = QtWidgets.QFileDialog.getExistingDirectory(self.main_window, 'Select Directory')

        ## Check and continue only if the directory is empty
        if os.path.isdir(self.file_path) and not os.listdir(self.file_path):

            # set path to the line edit
            line_edit.setText(self.file_path)
            self.main_window.working_dir_path = self.file_path

            ## Make directories: input, output, adj
            adj_filespath = os.path.join(self.file_path, "adj")
            os.mkdir(adj_filespath)

            proteins_path = os.path.join(self.file_path, "input")
            os.mkdir(proteins_path)

            output_path = os.path.join(self.file_path, "output")
            os.mkdir(os.path.join(self.file_path, "output"))

            ## Update config.ini file
            os.rename(os.path.join(self.main_window.tools_path, "config.ini"), os.path.join(self.main_window.tools_path, "config_temp.ini"))

            file_input = open(os.path.join(self.main_window.tools_path, "config_temp.ini"), "rt")
            file_output = open(os.path.join(self.main_window.tools_path, "config.ini"), "wt")

            for line in file_input:
                if line.startswith("working_dir_path ="):
                    file_output.write("working_dir_path =" + self.file_path + "\n")

                elif line.startswith("proteins_path ="):
                    file_output.write("proteins_path =" + proteins_path + "\n")

                elif line.startswith("output_path ="):
                    file_output.write("output_path =" + output_path + "\n")

                elif line.startswith("adj_filespath ="):
                    file_output.write("adj_filespath =" + adj_filespath + "\n")
                else:
                    file_output.write(line)

            file_input.close()
            file_output.close()

            os.remove(os.path.join(self.main_window.tools_path, "config_temp.ini"))

        # Error message if the directory is not empty
        elif os.path.isdir(self.file_path) and os.listdir(self.file_path):
            QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please chose an empty directory")

        # Pass if the directory is not selected
        elif not self.file_path:
            pass


    def next_page_func(self):
        self.main_window.TABS.setCurrentIndex(1)


class ScrollArea:

    def __init__(self):

        # Widget that contains the scroll Area
        self.widget = QtWidgets.QWidget()

        # Create the Scroll Area and add it to the widget
        self.scroll_area = QtWidgets.QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setWidget(self.widget)

        # Set the layout of the Scroll Area
        self.scroll_layout = QtWidgets.QGridLayout()
        self.widget.setLayout(self.scroll_layout)


class CENTRALITYtab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        # Additional text for centrality analysis
        add_text_centrality = ("Centrality Analysis implements main centrality algorithms. Supported algorithms: Closeness Centrality provides information about how close a node is to all other nodes; Betweenness Centrality provides a measure of how information can flow between nodes in a network; Eigenvector Centrality takes into account both the number of connections of a given node and its relevance in terms of information flow; Degree Centrality gives a measure of the relative connectivity of a node in the network.")

        self.centrality_info_text_area = QtWidgets.QPlainTextEdit()
        self.centrality_info_text_area.setReadOnly(True)

        self.centrality_info_text_area.setPlainText(add_text_centrality)

        self.layout.addWidget(self.centrality_info_text_area, 0, 0)


        # Create a Group Box for the algorithms to use in the centrality analysis
        self.node_centrality_box = QtWidgets.QGroupBox("Algorithm/s")
        self.layout.addWidget(self.node_centrality_box)

        # Add the Group Box to the CENTRALITY tab layout
        self.node_centrality_box_layout = QtWidgets.QVBoxLayout()
        self.node_centrality_box.setLayout(self.node_centrality_box_layout)

        # Add the algorithms to the node_centrality_box
        self.list_centrality_algorithms = []

        self.centrality_all = QtWidgets.QCheckBox("all")
        self.node_centrality_box_layout.addWidget(self.centrality_all)
        self.centrality_all.setChecked(True)

        self.centrality_clos = QtWidgets.QCheckBox("closeness")
        self.node_centrality_box_layout.addWidget(self.centrality_clos)

        self.centrality_betw = QtWidgets.QCheckBox("betweenness")
        self.node_centrality_box_layout.addWidget(self.centrality_betw)

        self.centrality_eigen = QtWidgets.QCheckBox("eigenvector_c")
        self.node_centrality_box_layout.addWidget(self.centrality_eigen)

        self.centrality_degree = QtWidgets.QCheckBox("degree_c")
        self.node_centrality_box_layout.addWidget(self.centrality_degree)

        # Create a list of centrality algorithms
        self.list_centrality_algorithms.extend([self.centrality_all,
        self.centrality_clos,
        self.centrality_betw,
        self.centrality_eigen,
        self.centrality_degree])

        for boxes in self.list_centrality_algorithms:
            boxes.setChecked(True)

        for boxes in self.list_centrality_algorithms[1:]:
            boxes.clicked.connect(lambda a: unselect_all(self.list_centrality_algorithms, self.centrality_all))

        self.centrality_all.clicked.connect(lambda a: click_all(self.list_centrality_algorithms, self.centrality_all))

        # Create a Run Analysis button, add it to the layout and connect it with the "run_analysis_func"
        self.run_analysis_button = QtWidgets.QPushButton("Run Analysis")
        self.layout.addWidget(self.run_analysis_button)
        self.run_analysis_button.clicked.connect(self.run_analysis_func)

    def run_analysis_func(self):

        # A function to check if the Inputs are present
        input_correct = check_inputs(self.main_window, self.main_window.INPUTS_widgets.pdb_line_edit, self.main_window.INPUTS_widgets.working_dir_line_edit)

        if not input_correct:
            pass
        else:
            if self.main_window.INPUTS_widgets.use_precomputed_matrix_cb.isChecked():
                adj_file_path = self.main_window.INPUTS_widgets.precomputed_matrix_linedit.text()
                adj_file_name = os.path.basename(adj_file_path)
                adj_mat_type_selection = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText()

                if adj_mat_type_selection == "alpha-Carbons":
                    adj_mat_type = "CA"
                elif adj_mat_type_selection == "beta-Carbons":
                    adj_mat_type = "CB"
                elif adj_mat_type_selection == "Centroids":
                    adj_mat_type = "centroid"

                new_adj_file_name = self.main_window.INPUTS_widgets.pdb_line_edit.text() + "_adj_{}_{}_{}.txt".format(adj_mat_type, self.main_window.INPUTS_widgets.non_covalent_box.value(), self.main_window.INPUTS_widgets.only_significant_box.value())
                new_adj_file_path = os.path.join(self.main_window.working_dir_path, "adj", new_adj_file_name)

                if os.path.isfile(new_adj_file_path):
                    os.remove(new_adj_file_path)

                if os.path.isfile(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name)):
                    os.remove(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name))

                shutil.copy(adj_file_path, os.path.join(self.main_window.working_dir_path, "adj"))
                os.rename(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name), new_adj_file_path)

                if not os.path.isfile(os.path.join(self.main_window.working_dir_path, "outputAdj", new_adj_file_name)):
                    shutil.copy(new_adj_file_path, os.path.join(self.main_window.working_dir_path, "outputAdj"))

            if self.main_window.INPUTS_widgets.use_pymol_protein_cb.isChecked():
                cmd.save(os.path.join(self.main_window.working_dir_path, "input", self.main_window.INPUTS_widgets.pdb_line_edit.text() + ".pdb"), self.main_window.INPUTS_widgets.pdb_line_edit.text())

            if self.main_window.INPUTS_widgets.USE_THREADS.isChecked():
                p_dialog = Protocol_exec_dialog(app=self, program=self,
                                                function=self.run,
                                                args=(),
                                                wait_start=0.4, wait_end=0.4,
                                                lock=True,
                                                stdout_silence=True,
                                                title="Running",
                                                label_text="PyPCN is running. Please wait")
                p_dialog.exec_()
            else:
                print("THREADS DISABLED - PLEASE WAIT")
                self.run()

        self.main_window.RESULTS_widgets.check_centrality()
        self.main_window.RESULTS_widgets.outputAdj_list = self.main_window.RESULTS_widgets.check_adj_directory()
        self.main_window.RESULTS_widgets.adj_frame.update_adj_frame(self.main_window.RESULTS_widgets.outputAdj_list)

        if input_correct:
            QtWidgets.QMessageBox.information(self.main_window, "Process Completed", "PyMOL Loaded only the last analysis.\nPlease explore the RESULTS tab to see all computations")



    def run(self):

        choice = get_adj_pdb_choice(self.main_window.INPUTS_widgets.use_precomputed_matrix_cb)
        # A function to create the string for the algorithms to be used
        self.get_algorithms()

        pcn_main = pcn_module.PCN_MAIN(self,
        pdb_input = self.main_window.INPUTS_widgets.pdb_line_edit,
        covalent_bonds_threshold = self.main_window.INPUTS_widgets.non_covalent_box,
        significant_bonds_threshold = self.main_window.INPUTS_widgets.only_significant_box,
        adj_mat_type = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText(),
        type_of_analysis = "centrality",
        algorithms_to_use = self.algorithms_to_use,
        initial_choice = choice,
        participation_coef_plot = self.main_window.INPUTS_widgets.coef_plot_box)

    def get_algorithms(self):

        list_of_algorithms = ["all", "closeness", "betweenness", "eigenvector_c", "degree_c"]

        self.algorithms_to_use = ""

        for box in self.list_centrality_algorithms:
            if box.text() == "all" and box.isChecked():
                for a in list_of_algorithms:
                    self.algorithms_to_use = self.algorithms_to_use + "," + a

            else:
                if box.isChecked():
                    self.algorithms_to_use = self.algorithms_to_use + "," + box.text()

        if self.algorithms_to_use[0] == ",":
            self.algorithms_to_use = self.algorithms_to_use[1:]


class SPECTRALtab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        # Additional text for Spectral Clustering
        add_text_spectral = ("Spectral Clustering extracts clusters from a graph with a clustering approach based on the Laplacian matrix eigenvectors. Both Hard (K-Means) and Soft (Fuzzy C-Means) clustering approaches are used on the eigenvectors of the Laplacian matrix (both normalized or unnormalized form). Shi Malik spectral clustering approach is also supported, to resolve the generalized eigenvalues problem. The parameter k”, which indicates the number of clusters, is required for all algorithms.")

        self.spectral_info_text_area = QtWidgets.QPlainTextEdit()
        self.spectral_info_text_area.setReadOnly(True)

        self.spectral_info_text_area.setPlainText(add_text_spectral)

        self.layout.addWidget(self.spectral_info_text_area, 0, 0, 1, 2)

        # Create a Group Box
        self.algorithm_box = QtWidgets.QGroupBox("Algorithm/s")

        # Add the Group Box to the SPECTRAL tab Layout
        self.layout.addWidget(self.algorithm_box, 1, 0, 1, 2)

        # Create and set the Group Box layout
        self.algorithm_box_layout = QtWidgets.QVBoxLayout()
        self.algorithm_box.setLayout(self.algorithm_box_layout)

        # Add the algorithms to the algorithm_box
        self.list_of_algorithms_boxes = []

        self.spectral_1_box = QtWidgets.QCheckBox("all")
        self.algorithm_box_layout.addWidget(self.spectral_1_box)
        self.spectral_1_box.setChecked(True)

        self.spectral_2_box = QtWidgets.QCheckBox("unnorm_ssc")
        self.algorithm_box_layout.addWidget(self.spectral_2_box)

        self.spectral_3_box = QtWidgets.QCheckBox("norm_ssc")
        self.algorithm_box_layout.addWidget(self.spectral_3_box)

        self.spectral_4_box = QtWidgets.QCheckBox("unnorm_hsc")
        self.algorithm_box_layout.addWidget(self.spectral_4_box)

        self.spectral_5_box = QtWidgets.QCheckBox("norm_hsc")
        self.algorithm_box_layout.addWidget(self.spectral_5_box)

        self.spectral_6_box = QtWidgets.QCheckBox("hsc_shimalik")
        self.algorithm_box_layout.addWidget(self.spectral_6_box)

        self.spectral_7_box = QtWidgets.QCheckBox("ssc_shimalik")
        self.algorithm_box_layout.addWidget(self.spectral_7_box)

        self.spectral_8_box = QtWidgets.QCheckBox("skl_spectral_clustering")
        self.algorithm_box_layout.addWidget(self.spectral_8_box)

        # Create a list of Spectral algorithms
        self.list_of_algorithms_boxes.extend([self.spectral_1_box,
        self.spectral_2_box,
        self.spectral_3_box,
        self.spectral_4_box,
        self.spectral_5_box,
        self.spectral_6_box,
        self.spectral_7_box,
        self.spectral_8_box])

        for boxes in self.list_of_algorithms_boxes:
            boxes.setChecked(True)

        for boxes in self.list_of_algorithms_boxes[1:]:
            boxes.clicked.connect(lambda a: unselect_all(self.list_of_algorithms_boxes, self.spectral_1_box))

        self.spectral_1_box.clicked.connect(lambda a: click_all(self.list_of_algorithms_boxes, self.spectral_1_box))

        self.cluster_box = ClusterK(self, self.layout)

        # Create a Run Analysis button, add it to the layout and connect it with the "run_analysis_func"
        self.run_analysis_button = QtWidgets.QPushButton("Run Analysis")
        self.layout.addWidget(self.run_analysis_button, 3, 0, 1, 2)
        self.run_analysis_button.clicked.connect(self.run_analysis_func)

    def run_analysis_func(self):

        # A function to check if the Inputs are present
        input_correct = check_inputs(self.main_window, self.main_window.INPUTS_widgets.pdb_line_edit, self.main_window.INPUTS_widgets.working_dir_line_edit)

        if not input_correct:
            pass
        else:
            if self.main_window.INPUTS_widgets.use_precomputed_matrix_cb.isChecked():
                adj_file_path = self.main_window.INPUTS_widgets.precomputed_matrix_linedit.text()
                adj_file_name = os.path.basename(adj_file_path)
                adj_mat_type_selection = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText()

                if adj_mat_type_selection == "alpha-Carbons":
                    adj_mat_type = "CA"
                elif adj_mat_type_selection == "beta-Carbons":
                    adj_mat_type = "CB"
                elif adj_mat_type_selection == "Centroids":
                    adj_mat_type = "centroid"

                new_adj_file_name = self.main_window.INPUTS_widgets.pdb_line_edit.text() + "_adj_{}_{}_{}.txt".format(adj_mat_type, self.main_window.INPUTS_widgets.non_covalent_box.value(), self.main_window.INPUTS_widgets.only_significant_box.value())
                new_adj_file_path = os.path.join(self.main_window.working_dir_path, "adj", new_adj_file_name)

                if os.path.isfile(new_adj_file_path):
                    os.remove(new_adj_file_path)

                if os.path.isfile(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name)):
                    os.remove(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name))

                shutil.copy(adj_file_path, os.path.join(self.main_window.working_dir_path, "adj"))
                os.rename(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name), new_adj_file_path)

                if not os.path.isfile(os.path.join(self.main_window.working_dir_path, "outputAdj", new_adj_file_name)):
                    shutil.copy(new_adj_file_path, os.path.join(self.main_window.working_dir_path, "outputAdj"))

            if self.main_window.INPUTS_widgets.use_pymol_protein_cb.isChecked():
                cmd.save(os.path.join(self.main_window.working_dir_path, "input", self.main_window.INPUTS_widgets.pdb_line_edit.text() + ".pdb"), self.main_window.INPUTS_widgets.pdb_line_edit.text(), format = "pdb")

            if self.main_window.INPUTS_widgets.USE_THREADS.isChecked():
                p_dialog = Protocol_exec_dialog(app=self, program=self,
                                                function=self.run,
                                                args=(),
                                                wait_start=0.4, wait_end=0.4,
                                                lock=True,
                                                stdout_silence=True,
                                                title="Running",
                                                label_text="PyPCN is running. Please wait")
                p_dialog.exec_()
            else:
                print("THREADS DISABLED - PLEASE WAIT")
                self.run()

        self.main_window.RESULTS_widgets.check_results(self.main_window.RESULTS_widgets.spectral_frame, "spectral")
        self.main_window.RESULTS_widgets.outputAdj_list = self.main_window.RESULTS_widgets.check_adj_directory()
        self.main_window.RESULTS_widgets.adj_frame.update_adj_frame(self.main_window.RESULTS_widgets.outputAdj_list)

        self.main_window.RESULTS_widgets.check_part_coeff_results(self.main_window.RESULTS_widgets.spectral_frame, "spectral")

        if input_correct:
            QtWidgets.QMessageBox.information(self.main_window, "Process Completed", "PyMOL Loaded only the last analysis.\nPlease explore the RESULTS tab to see all computations")


    def run(self):

        choice = get_adj_pdb_choice(self.main_window.INPUTS_widgets.use_precomputed_matrix_cb)

        # A function to create the string for the algorithms to be used
        self.get_algorithms()
        self.get_k()

        pcn_main = pcn_module.PCN_MAIN(self,
        pdb_input = self.main_window.INPUTS_widgets.pdb_line_edit,
        covalent_bonds_threshold = self.main_window.INPUTS_widgets.non_covalent_box,
        significant_bonds_threshold = self.main_window.INPUTS_widgets.only_significant_box,
        adj_mat_type = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText(),
        type_of_analysis = "spectral",
        initial_choice = choice,
        algorithms_to_use = self.algorithms_to_use,
        participation_coef_plot = self.main_window.INPUTS_widgets.coef_plot_box,
        k_value = self.k_value)

    # def check_inputs(self):
    #     self.input_correct = True
    #
    #     # Get the text from the pdb_line_edit
    #     # If the pdb code is not present, give a warning message
    #     if not self.main_window.INPUTS_widgets.pdb_line_edit.text():
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "There aren't enough parameters\n Please enter the PDB ID")
    #         self.input_correct = False
    #
    #     # Get the text from the configuration_line_edit
    #     # Control if it's a path leading to a folder, otherwise give a warning message
    #     input_line = self.main_window.INPUTS_widgets.configuration_line1_edit.text()
    #     if not os.path.isdir(input_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the input directory")
    #         self.input_correct = False
    #
    #     output_line = self.main_window.INPUTS_widgets.configuration_line2_edit.text()
    #     if not os.path.isdir(output_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the output directory")
    #         self.input_correct = False
    #
    #     adj_line = self.main_window.INPUTS_widgets.configuration_line3_edit.text()
    #     if not os.path.isdir(adj_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the adj directory")
    #         self.input_correct = False

    def get_algorithms(self):

        dict = {}
        dict["all"] = "0"
        dict["unnorm_ssc"] = "1"
        dict["norm_ssc"] = "2"
        dict["unnorm_hsc"] = "3"
        dict["norm_hsc"] = "4"
        dict["hsc_shimalik"] = "5"
        dict["ssc_shimalik"] = "6"
        dict["skl_spectral_clustering"] = "7"

        self.algorithms_to_use = ""

        for box in self.list_of_algorithms_boxes:
            if box.text() == "all" and box.isChecked():
                self.algorithms_to_use = "0"

            else:
                if box.isChecked():
                    number = dict[box.text()]
                    self.algorithms_to_use = self.algorithms_to_use + "," + number

        if self.algorithms_to_use[0] == ",":
            self.algorithms_to_use = self.algorithms_to_use[1:]


    def get_k(self):

        if self.cluster_box.best_k_radiobutton.isChecked():
            self.k_value = "best_k"

        if self.cluster_box.choose_k_radiobutton.isChecked():
            self.k_value = self.cluster_box.choose_k_box.value()



class ClusterK(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        # Create a Group Box
        self.cluster_k_box = QtWidgets.QGroupBox("Clustering parameters")

        # Add the Group Box to the SPECTRAL or EMBEDDINGS tab Layout
        self.layout.addWidget(self.cluster_k_box, 2, 0, 1, 2)

        # Create and set the Group Box layout
        self.cluster_k_box_layout = QtWidgets.QGridLayout()
        self.cluster_k_box.setLayout(self.cluster_k_box_layout)

        # Add k to the cluster_k_box
        self.best_k_radiobutton = QtWidgets.QRadioButton("Best k")
        self.n_best_ks_label = QtWidgets.QLabel("N° of best ks to try:")
        self.choose_k_radiobutton = QtWidgets.QRadioButton("Choose a k")

        self.best_k_radiobutton.setChecked(True)
        self.best_k_radiobutton.clicked.connect(self.change_state)
        self.choose_k_radiobutton.clicked.connect(self.change_state)

        self.cluster_k_box_layout.addWidget(self.best_k_radiobutton, 0, 0)
        self.cluster_k_box_layout.addWidget(self.n_best_ks_label, 0, 1)
        self.cluster_k_box_layout.addWidget(self.choose_k_radiobutton, 1, 0)

        self.n_best_ks_box = QtWidgets.QSpinBox()
        self.n_best_ks_box.setRange(1, 50)
        self.n_best_ks_box.setSingleStep(1)
        self.n_best_ks_box.setEnabled(True)
        self.cluster_k_box_layout.addWidget(self.n_best_ks_box, 0, 2)

        self.choose_k_box = QtWidgets.QSpinBox()
        self.choose_k_box.setRange(1, 50)
        self.choose_k_box.setSingleStep(1)
        self.choose_k_box.setEnabled(False)
        self.cluster_k_box_layout.addWidget(self.choose_k_box, 1, 2)


    def change_state(self):

        if self.best_k_radiobutton.isChecked():
            self.choose_k_box.setEnabled(False)
            self.n_best_ks_box.setEnabled(True)

        if self.choose_k_radiobutton.isChecked():
            self.choose_k_box.setEnabled(True)
            self.n_best_ks_box.setEnabled(False)


class EMBEDDINGStab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        # Additional text for embeddings+clustering
        add_text_embeddings = ("Embeddings + Clustering approach uses one of the embedding algorithms in the GEM library and then applies clustering. Supported algorithms: Node2vec (additional parameters random walk length and number of random walks for each node are required), HOPE (the additional parameter beta, which indicates the decay factor, is required) and Laplacianeigenmaps embedding, followed by a supported spectral clustering algorithm.")

        self.embeddings_info_text_area = QtWidgets.QPlainTextEdit()
        self.embeddings_info_text_area.setReadOnly(True)

        self.embeddings_info_text_area.setPlainText(add_text_embeddings)

        self.layout.addWidget(self.embeddings_info_text_area, 0, 0, 1, 2)

        # Create a Group Box
        self.algorithm_box = QtWidgets.QGroupBox("Algorithm/s")

        # Add the Group Box to the EMBEDDINGS tab Layout
        self.layout.addWidget(self.algorithm_box, 1, 0)

        # Add the algorithms to the algorithm_box
        self.list_of_algorithms_boxes = []

        # Create and set the Group Box layout
        self.algorithm_box_layout = QtWidgets.QVBoxLayout()
        self.algorithm_box.setLayout(self.algorithm_box_layout)

        # Add the algorithms to the algorithm_box
        self.embeddings_1_box = QtWidgets.QCheckBox("all")
        self.algorithm_box_layout.addWidget(self.embeddings_1_box)
        self.embeddings_1_box.setChecked(True)

        self.embeddings_2_box = QtWidgets.QCheckBox("fuzzycmeans_hope")
        self.algorithm_box_layout.addWidget(self.embeddings_2_box)

        self.embeddings_3_box = QtWidgets.QCheckBox("kmeans_hope")
        self.algorithm_box_layout.addWidget(self.embeddings_3_box)

        self.embeddings_4_box = QtWidgets.QCheckBox("fuzzycmeans_laplacianeigenmaps")
        self.algorithm_box_layout.addWidget(self.embeddings_4_box)

        self.embeddings_5_box = QtWidgets.QCheckBox("kmeans_laplacianeigenmaps")
        self.algorithm_box_layout.addWidget(self.embeddings_5_box)

        self.embeddings_6_box = QtWidgets.QCheckBox("fuzzycmeans_node2vec")
        self.algorithm_box_layout.addWidget(self.embeddings_6_box)

        self.embeddings_7_box = QtWidgets.QCheckBox("kmeans_node2vec")
        self.algorithm_box_layout.addWidget(self.embeddings_7_box)

        # Create a list of Embeddings algorithms
        self.list_of_algorithms_boxes.extend([self.embeddings_1_box,
        self.embeddings_2_box,
        self.embeddings_3_box,
        self.embeddings_4_box,
        self.embeddings_5_box,
        self.embeddings_6_box,
        self.embeddings_7_box])

        for boxes in self.list_of_algorithms_boxes:
            boxes.setChecked(True)

        for boxes in self.list_of_algorithms_boxes[1:]:
            boxes.clicked.connect(lambda a: unselect_all(self.list_of_algorithms_boxes, self.embeddings_1_box))

        self.embeddings_1_box.clicked.connect(lambda a: click_all(self.list_of_algorithms_boxes, self.embeddings_1_box))

        self.cluster_box = ClusterK(self, self.layout)

        # Create an additional parameters box
        self.add_parameters_box = QtWidgets.QGroupBox("Additional parameters")

        # Add the Group Box to the EMBEDDINGS tab Layout
        self.layout.addWidget(self.add_parameters_box, 1, 1)

        # Create and set the Group Box layout
        self.add_parameters_box_layout = QtWidgets.QGridLayout()
        self.add_parameters_box.setLayout(self.add_parameters_box_layout)

        # Create and add the widgets to the additional parameters box
        self.beta_label = QtWidgets.QLabel("Decay factor beta:")
        self.add_parameters_box_layout.addWidget(self.beta_label, 0, 0)
        self.beta_box = QtWidgets.QDoubleSpinBox()
        self.beta_box.setRange(0, 50)
        self.beta_box.setSingleStep(0.01)
        self.beta_box.setValue(0.01)
        self.add_parameters_box_layout.addWidget(self.beta_box, 0, 1)

        self.walk_len_label = QtWidgets.QLabel("Random walk lenght:")
        self.add_parameters_box_layout.addWidget(self.walk_len_label, 1, 0)
        self.walk_len_box = QtWidgets.QSpinBox()
        self.walk_len_box.setRange(0, 1000)
        self.walk_len_box.setSingleStep(1)
        self.walk_len_box.setValue(100)
        self.add_parameters_box_layout.addWidget(self.walk_len_box, 1, 1)

        self.num_walks_label = QtWidgets.QLabel("Number of random walks:")
        self.add_parameters_box_layout.addWidget(self.num_walks_label, 2, 0)
        self.num_walks_box = QtWidgets.QSpinBox()
        self.num_walks_box.setRange(0, 1000)
        self.num_walks_box.setSingleStep(1)
        self.num_walks_box.setValue(100)
        self.add_parameters_box_layout.addWidget(self.num_walks_box, 2, 1)

        # Create a d-parameter Box
        self.d_parameter_box = QtWidgets.QGroupBox("d parameter")

        # Add the d-parameter Box to the EMEDDINGS tab Layout
        self.layout.addWidget(self.d_parameter_box, 3, 0, 1, 2)

        # Create and set the d-parameter Box layout
        self.d_parameter_box_layout = QtWidgets.QGridLayout()
        self.d_parameter_box.setLayout(self.d_parameter_box_layout)

        # Create and Add the widgets to the d-parameter Box
        self.enter_d_label = QtWidgets.QLabel("Enter d parameter for d-dimensional embedding:")
        self.d_parameter_box_layout.addWidget(self.enter_d_label, 1, 0)

        self.enter_d_box = QtWidgets.QSpinBox()
        self.enter_d_box.setRange(2, 1000)
        self.enter_d_box.setSingleStep(1)
        self.enter_d_box.setValue(0)
        self.d_parameter_box_layout.addWidget(self.enter_d_box, 1, 1)

        # Create a Run Analysis button, add it to the layout and connect it with the "run_analysis_func"
        self.run_analysis_button = QtWidgets.QPushButton("Run Analysis")
        self.layout.addWidget(self.run_analysis_button, 4, 0, 1, 2)
        self.run_analysis_button.clicked.connect(self.run_analysis_func)


    def run_analysis_func(self):

        # A function to check if the Inputs are present
        input_correct = check_inputs(self.main_window, self.main_window.INPUTS_widgets.pdb_line_edit, self.main_window.INPUTS_widgets.working_dir_line_edit)

        if not input_correct:
            pass
        else:
            if self.main_window.INPUTS_widgets.use_precomputed_matrix_cb.isChecked():
                adj_file_path = self.main_window.INPUTS_widgets.precomputed_matrix_linedit.text()
                adj_file_name = os.path.basename(adj_file_path)
                adj_mat_type_selection = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText()

                if adj_mat_type_selection == "alpha-Carbons":
                    adj_mat_type = "CA"
                elif adj_mat_type_selection == "beta-Carbons":
                    adj_mat_type = "CB"
                elif adj_mat_type_selection == "Centroids":
                    adj_mat_type = "centroid"

                new_adj_file_name = self.main_window.INPUTS_widgets.pdb_line_edit.text() + "_adj_{}_{}_{}.txt".format(adj_mat_type, self.main_window.INPUTS_widgets.non_covalent_box.value(), self.main_window.INPUTS_widgets.only_significant_box.value())
                new_adj_file_path = os.path.join(self.main_window.working_dir_path, "adj", new_adj_file_name)

                if os.path.isfile(new_adj_file_path):
                    os.remove(new_adj_file_path)

                if os.path.isfile(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name)):
                    os.remove(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name))

                shutil.copy(adj_file_path, os.path.join(self.main_window.working_dir_path, "adj"))
                os.rename(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name), new_adj_file_path)

                if not os.path.isfile(os.path.join(self.main_window.working_dir_path, "outputAdj", new_adj_file_name)):
                    shutil.copy(new_adj_file_path, os.path.join(self.main_window.working_dir_path, "outputAdj"))

            if self.main_window.INPUTS_widgets.use_pymol_protein_cb.isChecked():
                cmd.save(os.path.join(self.main_window.working_dir_path, "input", self.main_window.INPUTS_widgets.pdb_line_edit.text() + ".pdb"), self.main_window.INPUTS_widgets.pdb_line_edit.text())

            if self.main_window.INPUTS_widgets.USE_THREADS.isChecked():
                p_dialog = Protocol_exec_dialog(app=self, program=self,
                                                function=self.run,
                                                args=(),
                                                wait_start=0.4, wait_end=0.4,
                                                lock=True,
                                                stdout_silence=True,
                                                title="Running",
                                                label_text="PyPCN is running. Please wait")
                p_dialog.exec_()
            else:
                print("THREADS DISABLED - PLEASE WAIT")
                self.run()

        self.main_window.RESULTS_widgets.check_results(self.main_window.RESULTS_widgets.embedd_clust_frame, "embedd")
        self.main_window.RESULTS_widgets.outputAdj_list = self.main_window.RESULTS_widgets.check_adj_directory()
        self.main_window.RESULTS_widgets.adj_frame.update_adj_frame(self.main_window.RESULTS_widgets.outputAdj_list)
        self.main_window.RESULTS_widgets.check_part_coeff_results(self.main_window.RESULTS_widgets.embedd_clust_frame, "embedd")

        if input_correct:
            QtWidgets.QMessageBox.information(self.main_window, "Process Completed", "PyMOL loaded only the last analysis.\nPlease explore the RESULTS tab to see all computations")


    def run(self):

        choice = get_adj_pdb_choice(self.main_window.INPUTS_widgets.use_precomputed_matrix_cb)

        # A function to create the string for the algorithms to be used
        self.get_algorithms()
        self.get_k()

#         beta_box
# walk_len_box
# num_walks_box

        pcn_main = pcn_module.PCN_MAIN(self,
        pdb_input = self.main_window.INPUTS_widgets.pdb_line_edit,
        covalent_bonds_threshold = self.main_window.INPUTS_widgets.non_covalent_box,
        significant_bonds_threshold = self.main_window.INPUTS_widgets.only_significant_box,
        adj_mat_type = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText(),
        type_of_analysis = "embeddings",
        algorithms_to_use = self.algorithms_to_use,
        initial_choice = choice,
        participation_coef_plot = self.main_window.INPUTS_widgets.coef_plot_box,
        k_value = self.k_value,
        d_value = self.enter_d_box,
        beta = self.beta_box.value(),
        walk_len = self.walk_len_box.value(),
        num_walks = self.num_walks_box.value())

    # def check_inputs(self):
    #     self.input_correct = True
    #
    #     # Get the text from the pdb_line_edit
    #     # If the pdb code is not present, give a warning message
    #     if not self.main_window.INPUTS_widgets.pdb_line_edit.text():
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "There aren't enough parameters\n Please enter the PDB ID")
    #         self.input_correct = False
    #
    #     # Get the text from the configuration_line_edit
    #     # Control if it's a path leading to a folder, otherwise give a warning message
    #     input_line = self.main_window.INPUTS_widgets.configuration_line1_edit.text()
    #     if not os.path.isdir(input_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the input directory")
    #         self.input_correct = False
    #
    #     output_line = self.main_window.INPUTS_widgets.configuration_line2_edit.text()
    #     if not os.path.isdir(output_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the output directory")
    #         self.input_correct = False
    #
    #     adj_line = self.main_window.INPUTS_widgets.configuration_line3_edit.text()
    #     if not os.path.isdir(adj_line):
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Please choose the adj directory")
    #         self.input_correct = False


    def get_algorithms(self):

        dict = {}
        dict["all"] = "0"
        dict["fuzzycmeans_hope"] = "1"
        dict["kmeans_hope"] = "2"
        dict["fuzzycmeans_laplacianeigenmaps"] = "3"
        dict["kmeans_laplacianeigenmaps"] = "4"
        dict["fuzzycmeans_node2vec"] = "5"
        dict["kmeans_node2vec"] = "6"

        self.algorithms_to_use = ""

        for box in self.list_of_algorithms_boxes:
            if box.text() == "all" and box.isChecked():
                self.algorithms_to_use = "0"

            else:
                if box.isChecked():
                    number = dict[box.text()]
                    self.algorithms_to_use = self.algorithms_to_use + "," + number

        # To delete commas at the beginning of the string
        if self.algorithms_to_use[0] == ",":
            self.algorithms_to_use = self.algorithms_to_use[1:]


    def get_k(self):

        if self.cluster_box.best_k_radiobutton.isChecked():
            self.k_value = "best_k"

        if self.cluster_box.choose_k_radiobutton.isChecked():
            self.k_value = self.cluster_box.choose_k_box.value()


class COMMUNITYtab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):


        # Additional text for community extraction
        add_text_community = ("Community Detection uses one of the community detection algorithms in the cdlib library. Supported algorithms: Louvain, Leiden, Walktrap, Infomap, Asyn FluidC (requires a fixed number of communities “k” for each protein), Greedy Modularity, Spinglass.")

        self.community_info_text_area = QtWidgets.QPlainTextEdit()
        self.community_info_text_area.setReadOnly(True)

        self.community_info_text_area.setPlainText(add_text_community)

        self.layout.addWidget(self.community_info_text_area, 0, 0, 1, 2)

        # Create a Group Box
        self.algorithm_box = QtWidgets.QGroupBox("Algorithm/s")

        # Add the Group Box to the COMMUNITY tab Layout
        self.layout.addWidget(self.algorithm_box, 1, 0)

        # Add the algorithms to the algorithm_box
        self.list_of_algorithms_boxes = []

        # Create and set the Group Box layout
        self.algorithm_box_layout = QtWidgets.QVBoxLayout()
        self.algorithm_box.setLayout(self.algorithm_box_layout)

        # Add the algorithms to the algorithm_box
        self.community_1_box = QtWidgets.QCheckBox("all")
        self.algorithm_box_layout.addWidget(self.community_1_box)
        self.community_1_box.setChecked(True)

        self.community_2_box = QtWidgets.QCheckBox("louvain")
        self.algorithm_box_layout.addWidget(self.community_2_box)

        self.community_3_box = QtWidgets.QCheckBox("leiden")
        self.algorithm_box_layout.addWidget(self.community_3_box)

        self.community_4_box = QtWidgets.QCheckBox("walktrap")
        self.algorithm_box_layout.addWidget(self.community_4_box)

        self.community_5_box = QtWidgets.QCheckBox("asyn_fluidc")
        self.algorithm_box_layout.addWidget(self.community_5_box)

        self.community_6_box = QtWidgets.QCheckBox("greedy_modularity")
        self.algorithm_box_layout.addWidget(self.community_6_box)

        self.community_7_box = QtWidgets.QCheckBox("infomap")
        self.algorithm_box_layout.addWidget(self.community_7_box)

        if sys.platform == "win32":
            self.community_7_box.setEnabled(False)
            self.community_7_box.setChecked(False)

        self.community_8_box = QtWidgets.QCheckBox("spinglass")
        self.algorithm_box_layout.addWidget(self.community_8_box)

        # Create a list of Community algorithms
        self.list_of_algorithms_boxes.extend([self.community_1_box,
        self.community_2_box,
        self.community_3_box,
        self.community_4_box,
        self.community_5_box,
        self.community_6_box,
        self.community_7_box,
        self.community_8_box])

        for boxes in self.list_of_algorithms_boxes:
            boxes.setChecked(True)

        for boxes in self.list_of_algorithms_boxes[1:]:
            boxes.clicked.connect(lambda a: unselect_all(self.list_of_algorithms_boxes, self.community_1_box))

        self.community_1_box.clicked.connect(lambda a: click_all(self.list_of_algorithms_boxes, self.community_1_box))

        # Create an additional parameters box
        self.add_parameters_box = QtWidgets.QGroupBox("Additional parameter for asyn_fluidc")

        # Add the Group Box to the EMBEDDINGS tab Layout
        self.layout.addWidget(self.add_parameters_box, 1, 1)

        # Create and set the Group Box layout
        self.add_parameters_box_layout = QtWidgets.QGridLayout()
        self.add_parameters_box.setLayout(self.add_parameters_box_layout)

        # Create and add the widgets to the add_parameters_box
        self.n_of_communities_label = QtWidgets.QLabel("Number of communities:")
        self.add_parameters_box_layout.addWidget(self.n_of_communities_label, 1, 0)
        self.n_of_communities_box = QtWidgets.QSpinBox()
        self.n_of_communities_box.setRange(0, 50)
        self.n_of_communities_box.setSingleStep(1)
        self.n_of_communities_box.setValue(1)
        self.add_parameters_box_layout.addWidget(self.n_of_communities_box, 1, 1)


        # Create a Run Analysis button, add it to the layout and connect it with the "run_analysis_func"
        self.run_analysis_button = QtWidgets.QPushButton("Run Analysis")
        self.layout.addWidget(self.run_analysis_button, 2, 0, 1, 2)
        self.run_analysis_button.clicked.connect(self.run_analysis_func)

    def run_analysis_func(self):

        # A function to check if the Inputs are present
        input_correct = check_inputs(self.main_window, self.main_window.INPUTS_widgets.pdb_line_edit, self.main_window.INPUTS_widgets.working_dir_line_edit)

        if not input_correct:
            pass
        else:
            if self.main_window.INPUTS_widgets.use_precomputed_matrix_cb.isChecked():
                adj_file_path = self.main_window.INPUTS_widgets.precomputed_matrix_linedit.text()
                adj_file_name = os.path.basename(adj_file_path)
                adj_mat_type_selection = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText()

                if adj_mat_type_selection == "alpha-Carbons":
                    adj_mat_type = "CA"
                elif adj_mat_type_selection == "beta-Carbons":
                    adj_mat_type = "CB"
                elif adj_mat_type_selection == "Centroids":
                    adj_mat_type = "centroid"

                new_adj_file_name = self.main_window.INPUTS_widgets.pdb_line_edit.text() + "_adj_{}_{}_{}.txt".format(adj_mat_type, self.main_window.INPUTS_widgets.non_covalent_box.value(), self.main_window.INPUTS_widgets.only_significant_box.value())
                new_adj_file_path = os.path.join(self.main_window.working_dir_path, "adj", new_adj_file_name)

                if os.path.isfile(new_adj_file_path):
                    os.remove(new_adj_file_path)

                if os.path.isfile(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name)):
                    os.remove(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name))

                shutil.copy(adj_file_path, os.path.join(self.main_window.working_dir_path, "adj"))
                os.rename(os.path.join(self.main_window.working_dir_path, "adj", adj_file_name), new_adj_file_path)

                if not os.path.isfile(os.path.join(self.main_window.working_dir_path, "outputAdj", new_adj_file_name)):
                    shutil.copy(new_adj_file_path, os.path.join(self.main_window.working_dir_path, "outputAdj"))

            if self.main_window.INPUTS_widgets.use_pymol_protein_cb.isChecked():
                cmd.save(os.path.join(self.main_window.working_dir_path, "input", self.main_window.INPUTS_widgets.pdb_line_edit.text() + ".pdb"), self.main_window.INPUTS_widgets.pdb_line_edit.text())

            if self.main_window.INPUTS_widgets.USE_THREADS.isChecked():
                p_dialog = Protocol_exec_dialog(app=self, program=self,
                                                function=self.run,
                                                args=(),
                                                wait_start=0.4, wait_end=0.4,
                                                lock=True,
                                                stdout_silence=True,
                                                title="Running",
                                                label_text="PyPCN is running. Please wait")
                p_dialog.exec_()
            else:
                print("THREADS DISABLED - PLEASE WAIT")
                self.run()

        self.main_window.RESULTS_widgets.check_results(self.main_window.RESULTS_widgets.community_frame, "community")
        self.main_window.RESULTS_widgets.outputAdj_list = self.main_window.RESULTS_widgets.check_adj_directory()
        self.main_window.RESULTS_widgets.adj_frame.update_adj_frame(self.main_window.RESULTS_widgets.outputAdj_list)

        self.main_window.RESULTS_widgets.check_part_coeff_results(self.main_window.RESULTS_widgets.community_frame, "community")

        if input_correct:
            QtWidgets.QMessageBox.information(self.main_window, "Process Completed", "PyMOL Loaded only the last analysis.\nPlease explore the RESULTS tab to see all computations")


    def run(self):

        choice = get_adj_pdb_choice(self.main_window.INPUTS_widgets.use_precomputed_matrix_cb)

        # A function to create the string for the algorithms to be used
        self.get_algorithms()
        n_comm = self.n_of_communities_box.value()

        pcn_main = pcn_module.PCN_MAIN(self,
        pdb_input = self.main_window.INPUTS_widgets.pdb_line_edit,
        covalent_bonds_threshold = self.main_window.INPUTS_widgets.non_covalent_box,
        significant_bonds_threshold = self.main_window.INPUTS_widgets.only_significant_box,
        adj_mat_type = self.main_window.INPUTS_widgets.adj_mat_type_cb.currentText(),
        type_of_analysis = "community",
        algorithms_to_use = self.algorithms_to_use,
        initial_choice = choice,
        participation_coef_plot = self.main_window.INPUTS_widgets.coef_plot_box,
        k_value = n_comm)

    # def check_inputs(self):
    #     self.input_correct = True
    #
    #     # Get the text from the pdb_line_edit
    #     # If the pdb code is not present, give a warning message
    #     if not self.main_window.INPUTS_widgets.pdb_line_edit.text():
    #         QtWidgets.QMessageBox.warning(self.main_window, "Warning", "There aren't enough parameters\n Please enter the PDB ID")
    #

    def get_algorithms(self):

        dict = {}
        dict["all"] = "0"
        dict["louvain"] = "1"
        dict["leiden"] = "2"
        dict["walktrap"] = "3"
        dict["asyn_fluidc"] = "4"
        dict["greedy_modularity"] = "5"

        if sys.platform == "win32":
            dict["spinglass"] = "6"
        else:
            dict["infomap"] = "6"
            dict["spinglass"] = "7"

        self.algorithms_to_use = ""

        for box in self.list_of_algorithms_boxes:
            if box.text() == "all" and box.isChecked():
                self.algorithms_to_use = "0"

            else:
                if box.isChecked():
                    number = dict[box.text()]
                    self.algorithms_to_use = self.algorithms_to_use + "," + number

        if self.algorithms_to_use[0] == ",":
            self.algorithms_to_use = self.algorithms_to_use[1:]



class RESULTStab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()


    def add_widgets(self):

        # Create a Group Box for information about Spectral Clustering Analysis
        self.group_boxes_dict = {}
        self.contact_map_gb = QtWidgets.QGroupBox("Contact Maps")
        self.centrality_analysis_gb = QtWidgets.QGroupBox("CENTRALITY ANALYSIS")
        self.spectral_clust_gb = QtWidgets.QGroupBox("SPECTRAL CLUSTERING")
        self.embedd_blust_gb = QtWidgets.QGroupBox("EMBEDDINGS + CLUSTERING")
        self.comm_extr_gb = QtWidgets.QGroupBox("COMMUNITY EXTRACTION")

        # Set layouts
        self.contact_map_gb.setLayout(QtWidgets.QGridLayout())
        self.centrality_analysis_gb.setLayout(QtWidgets.QGridLayout())
        self.spectral_clust_gb.setLayout(QtWidgets.QGridLayout())
        self.embedd_blust_gb.setLayout(QtWidgets.QGridLayout())
        self.comm_extr_gb.setLayout(QtWidgets.QGridLayout())

        # Add the Group Box to the SPECTRAL tab layout
        self.layout.addWidget(self.contact_map_gb, 0, 0)
        self.layout.addWidget(self.centrality_analysis_gb, 1, 0)
        self.layout.addWidget(self.spectral_clust_gb, 2, 0)
        self.layout.addWidget(self.embedd_blust_gb, 3, 0)
        self.layout.addWidget(self.comm_extr_gb, 4, 0)

        # Create and set the Group Box layout
        self.info_spectral_box_layout = QtWidgets.QVBoxLayout()

        outputAdj_list = self.check_adj_directory()
        self.adj_frame = Frame(parent=None, main_window=self.main_window, list_of_algorithms = outputAdj_list)
        self.adj_frame.build_adj_matrices_frames()
        self.contact_map_gb.layout().addWidget(self.adj_frame)

        # Add wigets to CENTRALITY ANALYSIS
        self.centrality_frame = Frame(parent=None, main_window=self.main_window,
        list_of_algorithms = self.main_window.dict_of_algorithms["centrality"])
        self.centrality_frame.build_analysis_frame(alg_type = "centrality")
        self.centrality_analysis_gb.layout().addWidget(self.centrality_frame)

        self.check_centrality()

        # Add wigets to SPECTRAL CLUSTERING
        self.spectral_frame = Frame(parent=None, main_window=self.main_window,
        list_of_algorithms = self.main_window.dict_of_algorithms["spectral"])
        self.spectral_frame.build_analysis_frame(alg_type = "spectral")
        self.spectral_clust_gb.layout().addWidget(self.spectral_frame)

        self.check_results(self.spectral_frame, "spectral")
        self.check_part_coeff_results(self.spectral_frame, "spectral")


        # Add wigets to EMBEDDINGS + CLUSTERING
        self.embedd_clust_frame = Frame(parent=None, main_window=self.main_window,
        list_of_algorithms = self.main_window.dict_of_algorithms["embedd"])
        self.embedd_clust_frame.build_analysis_frame(alg_type = "embedd")
        self.embedd_blust_gb.layout().addWidget(self.embedd_clust_frame)

        self.check_results(self.embedd_clust_frame, "embedd")
        self.check_part_coeff_results(self.embedd_clust_frame, "embedd")

        # Add wigets to COMMUNITY EXTRACTION
        self.community_frame = Frame(parent=None, main_window=self.main_window,
        list_of_algorithms = self.main_window.dict_of_algorithms["community"])
        self.community_frame.build_analysis_frame(alg_type = "community")
        self.comm_extr_gb.layout().addWidget(self.community_frame)

        self.check_results(self.community_frame, "community")
        self.check_part_coeff_results(self.community_frame, "community")


    # def add_widgets(self):
    #
    #     # Create a Group Box for information about Spectral Clustering Analysis
    #     self.group_boxes_dict = {}
    #     self.contact_map_gb = QtWidgets.QGroupBox("Contact Maps")
    #     self.centrality_analysis_gb = QtWidgets.QGroupBox("CENTRALITY ANALYSIS")
    #     self.spectral_clust_gb = QtWidgets.QGroupBox("SPECTRAL CLUSTERING")
    #     self.embedd_blust_gb = QtWidgets.QGroupBox("EMBEDDINGS + CLUSTERING")
    #     self.comm_extr_gb = QtWidgets.QGroupBox("COMMUNITY EXTRACTION")
    #
    #     # Set layouts
    #     self.contact_map_gb.setLayout(QtWidgets.QGridLayout())
    #     self.centrality_analysis_gb.setLayout(QtWidgets.QGridLayout())
    #     self.spectral_clust_gb.setLayout(QtWidgets.QGridLayout())
    #     self.embedd_blust_gb.setLayout(QtWidgets.QGridLayout())
    #     self.comm_extr_gb.setLayout(QtWidgets.QGridLayout())
    #
    #     self.group_boxes_dict["Contact Maps"] = self.contact_map_gb
    #     self.group_boxes_dict["CENTRALITY ANALYSIS"] = self.centrality_analysis_gb
    #     self.group_boxes_dict["SPECTRAL CLUSTERING"] = self.spectral_clust_gb
    #     self.group_boxes_dict["EMBEDDINGS + CLUSTERING"] = self.embedd_blust_gb
    #     self.group_boxes_dict["COMMUNITY EXTRACTION"] = self.comm_extr_gb
    #
    #     # Add the Group Box to the SPECTRAL tab layout
    #     self.layout.addWidget(self.contact_map_gb, 0, 0)
    #
    #     self.analysis_type_cb = QtWidgets.QComboBox()
    #     self.analysis_type_cb.addItems(["CENTRALITY ANALYSIS", "SPECTRAL CLUSTERING", "EMBEDDINGS + CLUSTERING", "COMMUNITY EXTRACTION"])
    #     self.analysis_type_cb.currentTextChanged.connect(self.on_show_groupbox)
    #     # self.analysis_type_cb.view().pressed.connect(self.on_show_groupbox)
    #
    #     self.layout.addWidget(self.analysis_type_cb)
    #     self.layout.addWidget(self.centrality_analysis_gb)
    #
    #     # Create and set the Group Box layout
    #     self.info_spectral_box_layout = QtWidgets.QVBoxLayout()
    #     self.contact_map_gb.setLayout(self.info_spectral_box_layout)
    #
    #     self.check_adj_directory()
    #
    #     # Add wigets to CENTRALITY ANALYSIS
    #     pdb_id = "prova"
    #     self.betweenness_frame = Frame(parent=None, main_window=self, pdb_id = pdb_id)
    #     self.betweenness_frame.build_analysis_frame()
    #     self.centrality_analysis_gb.layout().addWidget(self.betweenness_frame)


    # def on_show_groupbox(self):
    #
    #     current_gb = self.analysis_type_cb.currentText()
    #
    #     for i in reversed(range((self.group_boxes_dict[current_gb]).layout().count())):
    #         self.group_boxes_dict[current_gb].layout().itemAt(i).widget().deleteLater()
    #
    #     # self.group_boxes_dict[current_gb].layout().itemAt(i).widget().deleteLater()
    #     self.group_boxes_dict[current_gb].layout().deleteLater()

    def check_part_coeff_results(self, frame, analysis):

        path = os.path.join(self.main_window.working_dir_path, "output")

        # If the "outputCentralities" directory exists
        if os.path.isdir(path):
            for alg in self.main_window.dict_of_algorithms[analysis]:
                # If the algorithm-specific sub-directory exists
                tmp_path = os.path.join(path, alg, "Part_coefs_Sessions")
                if os.path.isdir(tmp_path):
                    self.main_window.part_coeff_results_dict[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.part_coeff_results_dict[alg]["location"] = tmp_path

        if self.main_window.part_coeff_results_dict:
            frame.update_frame(analysis)


    def check_results(self, frame, analysis):

        for alg in self.main_window.dict_of_algorithms[analysis]:
            path = os.path.join(self.main_window.working_dir_path, "output" + alg)
            # If the "output<alg>" directory exists
            if os.path.isdir(path):
                # If the algorithm-specific sub-directory exists
                tmp_path = os.path.join(path, "Sessions")
                if os.path.isdir(tmp_path):
                    self.main_window.algorithms_results_dict[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.algorithms_results_dict[alg]["location"] = tmp_path

                tmp_path = os.path.join(path, "Clusters")
                if os.path.isdir(tmp_path):
                    self.main_window.algorithms_results_dict_clusters[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.algorithms_results_dict_clusters[alg]["location"] = tmp_path

                tmp_path = os.path.join(path, "Summary")
                if os.path.isdir(tmp_path):
                    self.main_window.algorithms_results_dict_summary[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.algorithms_results_dict_summary[alg]["location"] = tmp_path

        if self.main_window.algorithms_results_dict and self.main_window.algorithms_results_dict_clusters and self.main_window.algorithms_results_dict_summary:
            frame.update_frame(analysis)


    def check_centrality(self):

        path = os.path.join(self.main_window.working_dir_path, "outputCentralities")

        # If the "outputCentralities" directory exists
        if os.path.isdir(path):
            for alg in self.main_window.dict_of_algorithms["centrality"]:
                # If the algorithm-specific sub-directory exists
                tmp_path = os.path.join(path, alg, "Sessions")
                if os.path.isdir(tmp_path):
                    self.main_window.algorithms_results_dict[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.algorithms_results_dict[alg]["location"] = tmp_path
                tmp_path = os.path.join(path, alg, "Txt")
                if os.path.isdir(tmp_path):
                    self.main_window.algorithms_results_dict_txt[alg]["results"] = os.listdir(tmp_path)
                    self.main_window.algorithms_results_dict_txt[alg]["location"] = tmp_path

        if self.main_window.algorithms_results_dict and self.main_window.algorithms_results_dict_txt:
            self.centrality_frame.update_frame("centrality")


    def check_adj_directory(self):
        ## This function checks wether there are some files already stored in the working directory

        # os.path.basename(os.path.normpath(self.main_window.working_dir_path))
        outputAdj_path = os.path.join(self.main_window.working_dir_path, "outputAdj")

        if os.path.isdir(outputAdj_path):
            outputAdj_list = [file for file in os.listdir(outputAdj_path) if file != "adj_matrix_dict.json"]
        else:
            outputAdj_list = []

        return outputAdj_list


    # def add_adj_matrices_frames(self, outputAdj_list):
    #
    #     pdb_frame = Frame(parent=self, main_window=self.main_window, list_of_algorithms = outputAdj_list)
    #     pdb_frame.build_adj_matrices_frames()
    #     self.contact_map_gb.layout().addWidget(pdb_frame)

class OTHERtab_widgets(QtWidgets.QWidget):

    def __init__(self, parent, layout):
        super().__init__(parent)

        self.main_window = parent
        self.layout = layout

        self.add_widgets()

    def add_widgets(self):

        self.computed_consensus = {}

        # Additional text for community extraction
        add_text_community = ("")

        self.other_info_text_area = QtWidgets.QPlainTextEdit()
        self.other_info_text_area.setReadOnly(True)

        self.other_info_text_area.setPlainText(add_text_community)

        self.layout.addWidget(self.other_info_text_area, 0, 0, 1, 1)

        ############################ Advanced Options
        # Create advanced options group box
        self.advanced_options_group = QtWidgets.QGroupBox("Advanced Options")
        self.layout.addWidget(self.advanced_options_group, 2, 0)
        # Create and set the Group Box layout
        self.advanced_group_layout = QtWidgets.QGridLayout()
        self.advanced_options_group.setLayout(self.advanced_group_layout)

        # Add widgets probability
        self.prob_label = QtWidgets.QLabel("Probability cutoff: ")
        self.advanced_group_layout.addWidget(self.prob_label, 1, 0)
        self.prob_label_box = QtWidgets.QDoubleSpinBox()
        self.prob_label_box.setRange(0, 1)
        self.prob_label_box.setSingleStep(1)
        self.prob_label_box.setValue(0.9)
        self.advanced_group_layout.addWidget(self.prob_label_box, 1, 1)

        self.frames_radio_automatic = QtWidgets.QRadioButton("Derive from data")
        self.frames_radio_manual = QtWidgets.QRadioButton("Set manually")
        self.frames_radio_automatic.clicked.connect(self.get_radiostate)
        self.frames_radio_manual.clicked.connect(self.get_radiostate)
        self.advanced_group_layout.addWidget(self.frames_radio_automatic, 3, 0)
        self.advanced_group_layout.addWidget(self.frames_radio_manual, 4, 0)
        self.frames_label = QtWidgets.QLabel("MD frames: ")
        self.advanced_group_layout.addWidget(self.frames_label, 2, 0)
        self.frames_box = QtWidgets.QSpinBox()
        self.frames_box.setRange(0, 10000000)
        self.frames_box.setSingleStep(1)
        self.frames_box.setValue(2000000)
        self.advanced_group_layout.addWidget(self.frames_box, 4, 1)
        self.frames_box.setEnabled(False)
        self.frames_radio_automatic.setChecked(True)


        ############################ Map Contacts
        # Create a Group Box
        self.map_contacts = QtWidgets.QGroupBox("Map Contacts")

        # Add the Group Box to the COMMUNITY tab Layout
        self.layout.addWidget(self.map_contacts, 1, 0)

        # Create and set the Group Box layout
        self.map_contacts_layout = QtWidgets.QGridLayout()
        self.map_contacts.setLayout(self.map_contacts_layout)

        # Widgets for simulation 1
        self.map_contacts_data_label = QtWidgets.QLabel("Simulation 1: ")
        self.map_contacts_layout.addWidget(self.map_contacts_data_label, 0, 0)

        self.simulation_1_line = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.simulation_1_line, 0, 1)

        self.map_contacts_data_file = QtWidgets.QPushButton("Choose file ...")
        self.map_contacts_layout.addWidget(self.map_contacts_data_file, 0, 2)
        self.map_contacts_data_file.clicked.connect(lambda: self.open_external_data(self.simulation_1_line, "data", 1))

        # Widgets for simulation 2
        self.map_contacts_data_label = QtWidgets.QLabel("Simulation 2: ")
        self.map_contacts_layout.addWidget(self.map_contacts_data_label, 1, 0)

        self.simulation_2_line = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.simulation_2_line, 1, 1)

        self.map_contacts_data_file = QtWidgets.QPushButton("Choose file ...")
        self.map_contacts_layout.addWidget(self.map_contacts_data_file, 1, 2)
        self.map_contacts_data_file.clicked.connect(lambda: self.open_external_data(self.simulation_2_line, "data", 2))

        # Widgets for pdb file
        self.map_contacts_prot_label = QtWidgets.QLabel("Open PDB file: ")
        self.map_contacts_layout.addWidget(self.map_contacts_prot_label, 2, 0)
        # self.use_precomputed_matrix_cb.clicked.connect(self.show_precomputed_options)

        self.map_contacts_prot_edit = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.map_contacts_prot_edit, 2, 1)

        self.map_contacts_prot_file = QtWidgets.QPushButton("Choose file ...")
        self.map_contacts_layout.addWidget(self.map_contacts_prot_file, 2, 2)
        self.map_contacts_prot_file.clicked.connect(lambda: self.open_external_data(self.map_contacts_prot_edit, "prot"))

        # Widgets for computing consensus
        self.compute_consensus_button = QtWidgets.QPushButton("Compute Consensus")
        self.consensus_path_line = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.compute_consensus_button, 3, 0)
        self.map_contacts_layout.addWidget(self.consensus_path_line, 3, 1)
        self.compute_consensus_button.clicked.connect(self.compute_prob_consensus)
        self.show_consensus_map_button = QtWidgets.QPushButton("Show Contact Map")
        self.map_contacts_layout.addWidget(self.show_consensus_map_button, 3, 2)
        self.show_consensus_map_button.clicked.connect(self.show_consensus_map_func)

        # Widgets for computing difference
        self.compute_difference_button = QtWidgets.QPushButton("Compute Difference")
        self.difference_path_line = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.compute_difference_button, 4, 0)
        self.map_contacts_layout.addWidget(self.difference_path_line, 4, 1)
        self.compute_difference_button.clicked.connect(self.compute_difference)
        self.show_difference_map_button = QtWidgets.QPushButton("Show Contact Map")
        self.map_contacts_layout.addWidget(self.show_difference_map_button, 4, 2)
        self.show_difference_map_button.clicked.connect(self.show_difference_map_func)

        self.difference_path_line_map = QtWidgets.QLineEdit()
        self.map_contacts_layout.addWidget(self.difference_path_line_map, 5, 1)
        self.show_map_pymol_button = QtWidgets.QPushButton("Map in PyMOL")
        self.map_contacts_layout.addWidget(self.show_map_pymol_button, 5, 2)
        self.show_map_pymol_button.clicked.connect(self.map_external_data)
        # Create a Run Analysis button, add it to the layout and connect it with the "run_analysis_func"
        # self.run_analysis_button = QtWidgets.QPushButton("Go")
        # self.layout.addWidget(self.run_analysis_button, 2, 0, 1, 2)
        # self.run_analysis_button.clicked.connect(self.map_external_data)

    def get_radiostate(self):
        if not self.frames_radio_automatic.isChecked():
            self.frames_box.setEnabled(True)
        else:
            self.frames_box.setEnabled(False)

    def open_external_data(self, line_edit, data, num = int()):

        if data == "data":
            self.file_path = QtWidgets.QFileDialog.getOpenFileName(self, "Open Text File", "","text files (*.txt)")
            if self.file_path:

                QtWidgets.QMessageBox.warning(self.main_window, "Attention!", "Please load the protein for which data have been computed")
                if num == 1:
                    valid, data  = self.parse_data(self.file_path[0])
                    self.simulation_1 = data
                elif num == 2:
                    valid, data = self.parse_data(self.file_path[0])
                    self.simulation_2 = data
                if valid:
                    line_edit.setText(self.file_path[0])
            else:
                pass

        if data == "prot":
            self.file_path = QtWidgets.QFileDialog.getOpenFileName(self, "Open Text File", "","text files (*.pdb)")
            if self.file_path:
                line_edit.setText(self.file_path[0])
            else:
                pass


    def parse_data(self, file_path):

        # Define a list of possible separators to try
        possible_separators = [',', ';', '\t', ' ', '|', ':', '~']

        # Initialize variables to store the best result
        best_separator = None
        best_data = None

        for separator in possible_separators:
            try:
                with open(file_path, 'r') as file:
                    lines = file.readlines()

                # Split the first line using the current separator to check the column count
                sample_line = lines[0].strip()
                columns = sample_line.split(separator)

                # Ensure there are at least three columns
                if len(columns) == 3:
                    # If we reach this point, the separator is valid and the file has at least three columns
                    best_separator = separator
                    best_data = [line.strip().split(separator) for line in lines]
                    break  # Stop searching for separators
                else:
                    raise ValueError("File does not have three columns.")

            except Exception as e:
                pass
                # Handle any exceptions (e.g., file not found, permission denied, etc.)
                # print(f"Error with separator '{separator}': {str(e)}")

        if best_separator is not None:
            print(f"Best separator: '{best_separator}'")
            valid = True
            # print("Parsed data:")
            # for row in best_data:
            #     print(row)
        else:
            print("No valid separator found or file doesn't have three columns.")
            valid = False

        return valid, best_data


    def map_external_data(self):

        map_contacts_prot_path = self.map_contacts_prot_edit.text()
        data_path = self.difference_path_line_map.text()

        if os.path.isfile(data_path):

            with open(data_path, 'r') as f:
                data = [[float(num) for num in line.split()] for line in f]
                data_array = np.array(data)

        cmd.reinitialize()
        cmd.load(map_contacts_prot_path)

        print(data_array)
        # Extract the third column
        third_column = data_array[:, 3]

        # Calculate the minimum and maximum values of the third column
        min_value = third_column.min()
        max_value = third_column.max()

        # Normalize the third column between 0 and 1
        normalized_column = (third_column - min_value) / (max_value - min_value)

        # Create a copy of the original array and replace the third column with the normalized values
        normalized_array = data_array.copy()
        normalized_array[:, 3] = normalized_column


        for idx, i in enumerate(data_array):
            res1 = int(i[0])
            res2 = int(i[1])
            res3 = float(i[3])

            if res3 != 0.0:

                cmd.select("res1", "resi {} and n. CA".format(res1))
                cmd.select("res2", "resi {} and n. CA".format(res2))

                if i[3] > 0:

                    cgo_arrow("(res1)", "(res2)", color = "blue", radius = res3)

                elif i[3] < 0:
                    cgo_arrow("(res1)", "(res2)", color = "red", radius = res3)

                cmd.delete("res1")
                cmd.delete("res2")


    def compute_prob_consensus(self):

        # Get input files
        simulation_1 = self.simulation_1_line.text()
        simulation_2 = self.simulation_2_line.text()
        map_contacts_prot_path = self.map_contacts_prot_edit.text()
        directory = os.path.dirname(simulation_1)
        simulation_1_name = os.path.basename(os.path.normpath(simulation_1)).replace(".txt", "")
        simulation_2_name = os.path.basename(os.path.normpath(simulation_2)).replace(".txt", "")

        if simulation_1 and simulation_2 and map_contacts_prot_path:

            # Compute probabilities
            self.array_simulation_1 = self.compute_probabilities(self.simulation_1, frames = self.frames_radio_manual.isChecked())
            self.array_simulation_2 = self.compute_probabilities(self.simulation_2, frames = self.frames_radio_manual.isChecked())
            # Save to file
            prob_sim_1 = os.path.join(directory, "{}_prob.txt".format(simulation_1_name))
            prob_sim_2 = os.path.join(directory, "{}_prob.txt".format(simulation_2_name))
            np.savetxt(prob_sim_1, self.array_simulation_1)
            np.savetxt(prob_sim_2, self.array_simulation_2)

            # Compute consensus
            self.consensus_simulations = self.compute_consensus(self.array_simulation_1, self.array_simulation_2, self.prob_label_box.value())
            consensus_simulations_path = os.path.join(directory, "{}_{}_consensus.txt".format(simulation_1_name, simulation_2_name))
            np.savetxt(consensus_simulations_path, self.consensus_simulations)

            # Convert consensus to adjacency matrix format
            self.consensus_binary = self.convert_to_binary_adjacency_matrix(self.consensus_simulations, col = 2, binary = True)
            consensus_binary_path = os.path.join(directory, "{}_{}_adj_mat_cons.txt".format(simulation_1_name, simulation_2_name))
            self.consensus_path_line.setText(consensus_binary_path)
            np.savetxt(consensus_binary_path, self.consensus_binary)

            # Show interactive consensus contact map
            self.show_consensus_map_func()

        else:
            QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Missing Inputs")


    def compute_difference(self):

        # Get input files
        simulation_1 = self.simulation_1_line.text()
        simulation_2 = self.simulation_2_line.text()
        map_contacts_prot_path = self.map_contacts_prot_edit.text()
        directory = os.path.dirname(simulation_1)
        simulation_1_name = os.path.basename(os.path.normpath(simulation_1)).replace(".txt", "")
        simulation_2_name = os.path.basename(os.path.normpath(simulation_2)).replace(".txt", "")

        # self.array_simulation_1 = self.convert_to_array(self.simulation_1)
        # self.array_simulation_2 = self.convert_to_array(self.simulation_2)
        if simulation_1 and simulation_2 and map_contacts_prot_path:

            # Compute probabilities
            self.array_simulation_1 = self.compute_probabilities(self.simulation_1, frames = self.frames_radio_manual.isChecked())
            self.array_simulation_2 = self.compute_probabilities(self.simulation_2, frames = self.frames_radio_manual.isChecked())
            # Save to file
            prob_sim_1 = os.path.join(directory, "{}_prob.txt".format(simulation_1_name))
            prob_sim_2 = os.path.join(directory, "{}_prob.txt".format(simulation_2_name))
            np.savetxt(prob_sim_1, self.array_simulation_1)
            np.savetxt(prob_sim_2, self.array_simulation_2)

            # Normalize the fourth column of f1 and f2
            self.array_simulation_1[:, 3] /= np.max(self.array_simulation_1[:, 3])
            self.array_simulation_2[:, 3] /= np.max(self.array_simulation_2[:, 3])

            # Create a DataFrame 'df' based on f1
            df_data = self.array_simulation_1.copy()

            # Calculate the difference between the fourth column of f2 and f1
            df_data[:, 3] = self.array_simulation_1[:, 3] - self.array_simulation_2[:, 3]

            # Save to file
            difference_path = os.path.join(directory, "{}_{}_difference_contacts.txt".format(simulation_1_name, simulation_2_name))
            np.savetxt(difference_path, df_data)
            self.difference_path_line_map.setText(difference_path)

            # Convert consensus to adjacency matrix format
            self.difference_binary = self.convert_to_binary_adjacency_matrix(df_data, col = 3, binary = False)
            difference_binary_path = os.path.join(directory, "{}_{}_adj_mat_diff.txt".format(simulation_1_name, simulation_2_name))
            self.difference_path_line.setText(difference_binary_path)
            np.savetxt(difference_binary_path, self.difference_binary)

            self.map_external_data()
            self.show_difference_map_func()

        else:
            QtWidgets.QMessageBox.warning(self.main_window, "Warning", "Missing Inputs")


# from_numpy_matrix

    def convert_to_array(self, data):
        # Convert the list of lists to a NumPy array
        numpy_array = np.array(data, dtype=np.float)

        # Reshape the NumPy array to have 3 columns
        numpy_array = numpy_array.reshape(-1, 3)

        return numpy_array


    def compute_probabilities(self, data, frames = int()):

        # Convert the list of lists to a NumPy array
        numpy_array = np.array(data, dtype=np.int)

        # Reshape the NumPy array to have 3 columns
        numpy_array = numpy_array.reshape(-1, 3)

        # Define the cutoff value
        if frames:
            cutoff_value = self.frames_box.value()
        else:
            cutoff_value = self.get_frames(numpy_array)
        #cutoff_value = 2000000


        # # Calculate the new column by dividing the third column by the cutoff value as a float
        new_column = numpy_array[:, 2] / cutoff_value
        result_array = np.column_stack((numpy_array, new_column.astype(np.float)))

        return result_array


    def get_frames(self, numpy_array):

        cutoff = numpy_array[:, 2].max()
        return cutoff


    def compute_consensus(self, array1, array2, cutoff):

        # Check if both input arrays have values in the fourth column higher than the cutoff
        mask1 = array1[:, 3] > cutoff
        mask2 = array2[:, 3] > cutoff

        # Create the consensus array
        consensus_array = np.column_stack((array1[:, :2], (mask1 & mask2).astype(int)))

        return consensus_array

    def convert_to_binary_adjacency_matrix(self, data, col, binary):

        tot_lista = []
        tmp_dict = {}
        prev = "dummy"
        tmp_list = []

        for row in data:
            res1 = int(row[0])
            res2 = int(row[1])
            if binary:
                bin = int(row[col])
            else:
                bin = float(row[col])

            if res1 == 1:
                tmp_dict[res2] = bin

            if res1 != prev and res1 != 1:
                tot_lista.append(tmp_list)
                tmp_list = []
                for i in range(res2):
                    if i > 0:
                        tmp_list.append(tmp_dict[i])

            prev = res1
            tmp_list.append(bin)

        tot_lista.append(tmp_list)
        arr_adj = np.array(tot_lista)

        return arr_adj

    def show_consensus_map_func(self):

        adj = self.consensus_path_line.text()
        prot = self.map_contacts_prot_edit.text()

        if adj and prot:

            self.contact_map_vis = Contactmap(self, self.main_window)
            self.contact_map_vis.create_contact_map(adj, prot, "CA", threshold = 20, feature_type = "contact")

        else:
            QtWidgets.QMessageBox.warning(self.main_window, "Warning", "To continue, compute a consensus adjacency matrix first.")


    def show_difference_map_func(self):

        adj = self.difference_path_line.text()
        prot = self.map_contacts_prot_edit.text()

        if adj and prot:

            self.contact_map_vis = Contactmap(self, self.main_window)
            self.contact_map_vis.create_contact_map(adj, prot, "CA", threshold = 20, feature_type = "probdiff")

        else:
            QtWidgets.QMessageBox.warning(self.main_window, "Warning", "To continue, compute a consensus adjacency matrix first.")

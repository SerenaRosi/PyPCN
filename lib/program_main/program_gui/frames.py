# PyMOL.
import pymol
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol.Qt import QtWidgets, QtCore, QtGui

from pymol import Qt
from pymol import stored
from pymol import viewing
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QProcess

import json
import ast
import re
try:
    import pandas as pd
except:
    pass

from lib.program_main.program_gui.contact_map_visualization import *
from lib.program_main.program_gui.plots import *


class NewWindow(QtWidgets.QMainWindow):

    middle_layout_type = "qform"

    def __init__(self, parent,
                 title="New Window",
                 upper_frame_title="New Window Sub-title",
                 submit_command=None, submit_button_text="Submit",
                 with_scroll=True,
                 # geometry=None
                 ):

        super().__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)

        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()

        self.middle_scroll = QtWidgets.QScrollArea()
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        self.main_vbox.addWidget(self.middle_scroll)

        self.middle_layout_type = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_layout_type)

        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())
            self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)


class Frame(QtWidgets.QFrame):

    def __init__(self, parent, main_window,
                 pdb_id = "",
                 list_of_algorithms = [],
                 dict_of_results = {},
                 *args, **configs):

        super(Frame, self).__init__(main_window, *args, **configs)

        self.main_window = main_window

        self.pdb_id = pdb_id
        self.list_of_algorithms = list_of_algorithms

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 0

        # Builds a frame for each template structure and all its options.
        #self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.structure_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


    def build_adj_matrices_frames(self):

        self.adj_file_cb = QtWidgets.QComboBox()
        self.adj_file_cb.addItems(self.list_of_algorithms)

        self.show_contact_map = QtWidgets.QPushButton("Show Contact Map")
        self.show_contact_map.clicked.connect(self.show_contact_map_interactive_plot)

        self.structure_frame_layout.addWidget(self.adj_file_cb, 0, 0)
        self.structure_frame_layout.addWidget(self.show_contact_map, 0, 1)


    def build_analysis_frame(self, alg_type = ""):

        self.alg_type = alg_type

        self.algorithms_cb = QtWidgets.QComboBox()
        self.algorithms_cb.addItems(self.list_of_algorithms)
        self.algorithms_cb.currentTextChanged.connect(self.update_frame)
        self.algorithms_cb.view().pressed.connect(self.update_frame)
        self.results_cb = QtWidgets.QComboBox()

        self.view_in_pymol_pb = QtWidgets.QPushButton("View in PyMOL")
        self.view_in_pymol_pb.clicked.connect(self.view_in_pymol_func)
        self.save_to_file_pb = QtWidgets.QPushButton("Save to file")
        # self.save_to_file_pb.clicked.connect(self.save_to_file_pb)

        if self.alg_type == "centrality":
            self.show_plot_cb = QtWidgets.QPushButton("Show Plot")
            self.show_plot_cb.clicked.connect(self.show_plot_centrality)
            self.structure_frame_layout.addWidget(self.show_plot_cb, 0, 2)

        else:
            self.part_coeff_plot_cb = QtWidgets.QPushButton("Participation Coefficient Plot")
            self.part_coeff_plot_cb.clicked.connect(self.show_particpation_coeff_plot_options)
            self.structure_frame_layout.addWidget(self.part_coeff_plot_cb, 0, 2)

            # self.clusters_plot_cb = QtWidgets.QPushButton("Clusters Plot")
            # self.clusters_plot_cb.clicked.connect(self.show_clusters_plot)
            # self.structure_frame_layout.addWidget(self.clusters_plot_cb, 1, 2)


        self.structure_frame_layout.addWidget(self.algorithms_cb, 0, 0)
        self.structure_frame_layout.addWidget(self.results_cb, 1, 0)
        self.structure_frame_layout.addWidget(self.view_in_pymol_pb, 0, 1)
        #self.structure_frame_layout.addWidget(self.save_to_file_pb, 1, 1)


    def show_particpation_coeff_plot_options(self):

        self.part_coeff_options_window = NewWindow(parent = self.main_window,
        title = "Participation Coefficient Plot", upper_frame_title = "",
        submit_command = self.plot_p, submit_button_text="Plot", with_scroll = True)

        # GroupBox
        self.plotting_gb = QtWidgets.QGroupBox("Plotting Options")
        self.part_coeff_options_window.middle_layout_type.addWidget(self.plotting_gb)
        self.plotting_gb.setLayout(QtWidgets.QGridLayout())

        # Add single-multiple radiobuttons to layout
        self.number_of_plots = QtWidgets.QButtonGroup(self)

        self.single_mode = QtWidgets.QRadioButton("Single")
        self.single_mode.setChecked(True)
        self.number_of_plots.addButton(self.single_mode)

        self.multiple_mode = QtWidgets.QRadioButton("Multiple")
        self.number_of_plots.addButton(self.multiple_mode)

        self.single_mode.clicked.connect(self.get_single_multiple_state)
        self.multiple_mode.clicked.connect(self.get_single_multiple_state)

        self.plotting_gb.layout().addWidget(self.single_mode)
        self.plotting_gb.layout().addWidget(self.multiple_mode)

        self.single_mode_cb = QtWidgets.QComboBox()
        self.plotting_gb.layout().addWidget(self.single_mode_cb)

        self.multiple_mode_cb = QtWidgets.QListWidget()
        self.plotting_gb.layout().addWidget(self.multiple_mode_cb)

        self.multiple_mode_cb.setSelectionMode(3)

        self.single_mode_cb.setEnabled(True)
        self.multiple_mode_cb.setEnabled(False)

        self.plot_type = QtWidgets.QButtonGroup(self)

        self.vsNode_gb = QtWidgets.QRadioButton("vs_Node")
        self.vsZ_gb = QtWidgets.QRadioButton("vs_Z-intraconnectivity")

        self.plot_type.addButton(self.vsNode_gb)
        self.plot_type.addButton(self.vsZ_gb)

        self.plotting_gb.layout().addWidget(self.vsNode_gb)
        self.plotting_gb.layout().addWidget(self.vsZ_gb)

        self.vsNode_gb.setChecked(True)

        list_of_coeff = []

        for alg in self.main_window.part_coeff_results_dict:
            if alg in self.main_window.dict_of_algorithms[self.alg_type]:
                tmp = [a for a in self.main_window.part_coeff_results_dict[alg]["results"] if a.endswith('.txt') and not re.search("z_intraconn", a)]
                list_of_coeff.extend(tmp)

        self.single_mode_cb.addItems(list_of_coeff)
        self.multiple_mode_cb.addItems(list_of_coeff)

        self.part_coeff_options_window.show()


    def get_single_multiple_state(self):

        if self.single_mode.isChecked():
            self.single_mode_cb.setEnabled(True)
            self.multiple_mode_cb.setEnabled(False)

        elif self.multiple_mode.isChecked():
            self.multiple_mode_cb.setEnabled(True)
            self.single_mode_cb.setEnabled(False)


    def show_clusters_plot(self):

        # Initialize Plot Window
        plot_centrality = PlotCentrality(self, self.main_window)
        plot_centrality.plot_window()

        # Get algorithm name
        algo = self.algorithms_cb.currentText()
        # Get current selection
        current_selection = self.results_cb.currentText()
        # Get protein name
        protein_name = current_selection.split("_")[0]
        # Get location of the results
        location = self.main_window.algorithms_results_dict_clusters[algo]['location']

        # Modify session file name to clusters file name
        tmp = current_selection.replace("part_coefs", "Clusters")
        file_name = re.sub(r'_session.*?\.pse', '.txt', tmp)
        #file_name = tmp.replace("_session.pse", ".txt")
        #file_name = tmp.replace("_session.pse", ".txt")

        # Make path
        path_to_cluster_file = os.path.join(location, file_name)

        # Make dataframe
        file = open(path_to_cluster_file, "r")
        test_string = file.read()
        file.close()

        dictio = ast.literal_eval(test_string)

        df = pd.DataFrame.from_dict(dictio.items())
        df.columns = ['Residue', 'Cluster']
        x_data = df['Residue'].to_list()
        y_data = df['Cluster'].to_list()

        if os.path.isfile(path_to_cluster_file):
            plot_centrality.add_plot("", "", "", type = "scatter_plot",
                                     scatter_plot_data = [x_data, y_data],
                                     original_data = [],
                                     algorithm = algo,
                                     pdb_name = protein_name)

        plot_centrality.show_plot()


    def plot_p(self):

        plot_centrality = PlotCentrality(self, self.main_window)
        plot_centrality.plot_window()

        part_coeff_name_list = self.get_parameters()

        cmd.reinitialize()

        if self.vsNode_gb.isChecked():
            #if self.single_mode.isChecked():
            for part_coeff_name in part_coeff_name_list:
                data, algo, pdb = self.read_part_coeff_data(part_coeff_name)
                list_of_data = self.init_values(data, algo)
                plot_centrality.add_plot("Nodes (residues)", "Participation Coefficient", "", type = "line_plot",
                                         line_plot_data = list_of_data,
                                         original_data = data,
                                         algorithm = algo,
                                         pdb_name = pdb)

        if self.vsZ_gb.isChecked():
            for part_coeff_name in part_coeff_name_list:
                data, algo, pdb = self.read_part_coeff_data(part_coeff_name)
                z_data, algo = self.read_z_data(part_coeff_name)
                list_of_data = self.init_values(data, algo)
                list_of_z_data = self.init_values(z_data, algo)
                plot_centrality.add_plot("Participation Coefficient", "Z-intraconnectivity", "", type = "scatter_plot",
                                         scatter_plot_data = [list_of_data, list_of_z_data],
                                         original_data = data,
                                         algorithm = algo,
                                         pdb_name = pdb)

        plot_centrality.show_plot()


    def init_values(self, data_dict, algorithm):

        values_list = []

        for keys in data_dict:
            tmp_list = []
            resnum = keys.split()[0]
            res = resnum[0:3]
            num = resnum[3:]
            chain = keys.split()[1]
            value = data_dict[keys]

            values_list.append(value)

        return values_list


    def get_parameters(self):

        if self.single_mode.isChecked():
            part_coeff_name = [self.single_mode_cb.currentText()]

        if self.multiple_mode.isChecked():
            part_coeff_name = [a.text() for a in self.multiple_mode_cb.selectedItems()]

        return part_coeff_name


    def read_part_coeff_data(self, part_coeff_name):

        for alg in self.main_window.dict_of_algorithms[self.alg_type]:
            if re.search('_' + alg, part_coeff_name):
                algo = alg
                file_path = os.path.join(self.main_window.part_coeff_results_dict[alg]["location"], part_coeff_name)
                part_coeff_session_name = part_coeff_name.replace(".txt", "_session.pse")
                pse_file_path = os.path.join(self.main_window.part_coeff_results_dict[alg]["location"], part_coeff_session_name)

        cmd.load(pse_file_path, partial=1)

        pymol_obj = cmd.get_names('objects')
        for i in pymol_obj:
            if re.search(i, pse_file_path):
                pdb = i

        cmd.zoom(pdb)

        file = open(file_path, "r")
        test_string = file.read()
        file.close()

        dictio = ast.literal_eval(test_string)

        return dictio, algo, pdb


    def read_z_data(self, part_coeff_name):

        for alg in self.main_window.dict_of_algorithms[self.alg_type]:
            if re.search('_' + alg, part_coeff_name):
                algo = alg
                file_path = os.path.join(self.main_window.part_coeff_results_dict[alg]["location"], part_coeff_name.replace("part_coefs", "z_intraconn"))

        file = open(file_path, "r")
        test_string = file.read()
        file.close()

        dictio = ast.literal_eval(test_string)

        return dictio, algo

    def show_plot_centrality(self):

        # Get Parameters
        alg = self.algorithms_cb.currentText()
        result = self.results_cb.currentText()
        pdb_id = result.split("_")[0]
        adj_mat_type = (result.split("_")[2].replace("session", "")).replace(".pse", "")
        path = self.main_window.algorithms_results_dict[alg]["location"]
        path_to_selected = os.path.join(path, result)

        # Load in PyMOL
        if os.path.isfile(path_to_selected):
            cmd.reinitialize()
            cmd.load(path_to_selected)

        file_name = pdb_id + "_" + alg + adj_mat_type + ".txt"

        path_to_file = os.path.join(self.main_window.algorithms_results_dict_txt[alg]["location"], file_name)

        # Read data from file to dict
        with open(path_to_file) as f:
            data = f.read()

        data_dict = ast.literal_eval(data)

        plot_centrality = PlotCentrality(self, self.main_window)
        plot_centrality.plot_window()

        list_of_data = self.init_values(data_dict, alg)
        plot_centrality.add_plot("Participation Coefficient", "Nodes (residues)", "", type = "line_plot",
                                 line_plot_data = list_of_data,
                                 original_data = data_dict,
                                 algorithm = alg,
                                 pdb_name = pdb_id)

        plot_centrality.show_plot()


    def view_in_pymol_func(self):

        selected = self.results_cb.currentText()
        path = self.main_window.algorithms_results_dict[self.algorithms_cb.currentText()]["location"]

        path_to_selected = os.path.join(path, selected)

        if os.path.isfile(path_to_selected):
            cmd.reinitialize()
            cmd.load(path_to_selected)

        path = self.main_window.part_coeff_results_dict[self.algorithms_cb.currentText()]["location"]

        path_to_selected = os.path.join(path, selected)

        if os.path.isfile(path_to_selected):
            cmd.reinitialize()
            cmd.load(path_to_selected)

    def update_frame(self, algorithm):

        list_to_add = []
        for alg in self.main_window.part_coeff_results_dict:
            if self.algorithms_cb.currentText() == alg:
                list_to_add.extend(a for a in self.main_window.part_coeff_results_dict[alg]["results"] if a.endswith('.pse'))

        for alg in self.main_window.algorithms_results_dict:
            if self.algorithms_cb.currentText() == alg:
                list_to_add.extend(self.main_window.algorithms_results_dict[alg]["results"])

        self.results_cb.clear()
        self.results_cb.addItems(list_to_add)

    def update_adj_frame(self, outputAdj_list):

        self.adj_file_cb.clear()
        self.adj_file_cb.addItems(outputAdj_list)

    def show_contact_map_interactive_plot(self):

        if self.adj_file_cb.currentText():

            adj_file_name = self.adj_file_cb.currentText()
            pdb_name = self.adj_file_cb.currentText().split("_")[0]
            adj_mat_type = self.adj_file_cb.currentText().split("_")[2]

            self.contact_map_vis = Contactmap(self, self.main_window)
            self.contact_map_vis.create_contact_map(adj_file_name, pdb_name, adj_mat_type)

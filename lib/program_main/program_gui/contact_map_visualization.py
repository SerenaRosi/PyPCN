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

import numpy as np
import warnings
import sys
from io import StringIO
from .center_of_mass import *

# dict for 20 standard amino acids: 3 -> 1
code_standard = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
    'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
    'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G'}

# dict for 20 standard amino acids: 1 -> 3
code_standard_rev = {v: k for k, v in code_standard.items()}

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

#---------------------------------------------
# Viridis colormap, adapted from matplotlib. -
#---------------------------------------------

viridis_colors = [
    [0.267004, 0.004874, 0.329415],
    [0.272594, 0.025563, 0.353093],
    [0.277018, 0.050344, 0.375715],
    [0.280267, 0.073417, 0.397163],
    [0.282327, 0.094955, 0.417331],
    [0.283197, 0.11568, 0.436115],
    [0.282884, 0.13592, 0.453427],
    [0.281412, 0.155834, 0.469201],
    [0.278826, 0.17549, 0.483397],
    [0.275191, 0.194905, 0.496005],
    [0.270595, 0.214069, 0.507052],
    [0.265145, 0.232956, 0.516599],
    [0.258965, 0.251537, 0.524736],
    [0.252194, 0.269783, 0.531579],
    [0.244972, 0.287675, 0.53726],
    [0.237441, 0.305202, 0.541921],
    [0.229739, 0.322361, 0.545706],
    [0.221989, 0.339161, 0.548752],
    [0.214298, 0.355619, 0.551184],
    [0.206756, 0.371758, 0.553117],
    [0.19943, 0.387607, 0.554642],
    [0.192357, 0.403199, 0.555836],
    [0.185556, 0.41857, 0.556753],
    [0.179019, 0.433756, 0.55743],
    [0.172719, 0.448791, 0.557885],
    [0.166617, 0.463708, 0.558119],
    [0.160665, 0.47854, 0.558115],
    [0.154815, 0.493313, 0.55784],
    [0.149039, 0.508051, 0.55725],
    [0.143343, 0.522773, 0.556295],
    [0.13777, 0.537492, 0.554906],
    [0.132444, 0.552216, 0.553018],
    [0.127568, 0.566949, 0.550556],
    [0.123463, 0.581687, 0.547445],
    [0.120565, 0.596422, 0.543611],
    [0.119423, 0.611141, 0.538982],
    [0.120638, 0.625828, 0.533488],
    [0.12478, 0.640461, 0.527068],
    [0.132268, 0.655014, 0.519661],
    [0.143303, 0.669459, 0.511215],
    [0.157851, 0.683765, 0.501686],
    [0.175707, 0.6979, 0.491033],
    [0.196571, 0.711827, 0.479221],
    [0.220124, 0.725509, 0.466226],
    [0.24607, 0.73891, 0.452024],
    [0.274149, 0.751988, 0.436601],
    [0.304148, 0.764704, 0.419943],
    [0.335885, 0.777018, 0.402049],
    [0.369214, 0.788888, 0.382914],
    [0.404001, 0.800275, 0.362552],
    [0.440137, 0.811138, 0.340967],
    [0.477504, 0.821444, 0.318195],
    [0.515992, 0.831158, 0.294279],
    [0.555484, 0.840254, 0.269281],
    [0.595839, 0.848717, 0.243329],
    [0.636902, 0.856542, 0.21662],
    [0.678489, 0.863742, 0.189503],
    [0.720391, 0.87035, 0.162603],
    [0.762373, 0.876424, 0.137064],
    [0.804182, 0.882046, 0.114965],
    [0.845561, 0.887322, 0.099702],
    [0.886271, 0.892374, 0.095374],
    [0.926106, 0.89733, 0.104071],
    [0.964894, 0.902323, 0.123941]
]

#viridis_colors_hex = [convert_rgb_to_hex(c) for c in viridis_colors]

viridis_colors_rev = list(reversed(viridis_colors))


class Res():

    def __init__(self, main, res):
        # e.g. res = ('SER', ' 123', 'A')
        if res[0] in code_standard:
            resid = code_standard[res[0]]
        else:
            resid = res[0]
        num = res[1]
        self.res_name = resid + "__" + num
        self.db_index = int(num)
        self.three_letter_code = res[0]

        self.chain = res[2]


class Contactmap():

    def __init__(self, main, main_window):
        self.main = main
        self.main_window = main_window

    def create_contact_map(self, adj_matrix_file, pdb_name, adj_mat_type, threshold = "", feature_type = "contact"):

        outputAdj_path = os.path.join(self.main_window.working_dir_path, "outputAdj")
        adj_path = os.path.join(outputAdj_path, adj_matrix_file)

        if os.path.isfile(adj_path):

            with open(adj_path, 'r') as f:
                data_array = [[float(num) for num in line.split()] for line in f]

            # Reference residues and selectors lists. Used to interact with the 3D structures loaded in
            # PyMOL.
            self.target_map = []

            title = adj_matrix_file + " Interactive Contact Map"
            self.pixel_size = 5
            self.feature_type = feature_type
            if threshold:
                self.dist_threshold = threshold
            else:
                self.dist_threshold = float(adj_matrix_file.split("_")[4].replace(".txt", ""))
            self.interaction_center = adj_mat_type

            self.residues_list = []

            if os.path.isfile(pdb_name):
                path_to_pdb = pdb_name
                pdb_name = os.path.basename(os.path.normpath(path_to_pdb)).replace(".pdb", "")
                self.objsel = pdb_name
            else:
                path_to_pdb = os.path.join(self.main_window.working_dir_path, "input", pdb_name + ".pdb")
                self.objsel = pdb_name

            if os.path.isfile(path_to_pdb):

                res_list = []
                atoms = []
                num = 0

                # check if it has a chain
                # with open(path_to_pdb) as pdbfile:
                #     for line in pdbfile:
                #         if line[:4] == 'ATOM':
                #             if line[21].strip():
                #                 has_chain = True
                #                 break
                #             else:
                #                 has_chain = False

                # if not has_chain:
                #     cmd.reinitialize()
                #     cmd.load(path_to_pdb)
                #     cmd.do("alter (all),chain='/A/A/'")
                #     cmd.save("tmp.pdb", '(all)', 0, "pdb")
                #     cmd.reinitialize()
                #
                #     os.remove(path_to_pdb)
                #     os.rename("tmp.pdb", path_to_pdb)
                #     os.remove("tmp.pdb")

                with open(path_to_pdb) as pdbfile:
                    for line in pdbfile:
                        if line[:4] == 'ATOM':
                            if line[22:26] != num:
                                res_list.append((str(line[17:20]), line[22:26], line[21]))
                            if line[21].strip():
                                has_chain = True
                            else:
                                has_chain = False
                            num = line[22:26]
                            #split the line
                            splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                            atoms.append(splitted_line)


            for res in res_list:
                res_obj = Res(self, res)
                self.residues_list.append(res_obj)

            # Load the corresponding PDB file in PyMOL - TODO: if not present
            cmd.reinitialize()
            cmd.load(path_to_pdb)

            residues, coords_array, selectors = self.get_coords_array(self.objsel, self.residues_list, self.interaction_center)

            self.ref_residues = residues
            self.ref_selectors = selectors

            #----------------------------------------------------------------------------------
            # Builds the array with the contact/distance maps to show for a single structure. -
            #----------------------------------------------------------------------------------

            # if self.feature_type in ("contact", "distance"):
            #
            #     # Get the coordinates of the residues.
            #     # residues, coords_array, selectors = self.get_coords_array(self.target_sequence, self.interaction_center)
            #
            #     self.target_map = np.zeros((len(residues), len(residues)))
            #     for res_idx_i, res_i in enumerate(residues):
            #         for res_idx_j, res_j in enumerate(residues):
            #             if res_idx_i >= res_idx_j:
            #                 continue
            #             dist = get_distance(coords_array[res_idx_i], coords_array[res_idx_j])
            #             self.target_map[res_idx_i][res_idx_j] = self.get_map_pixel(dist)
            #             self.target_map[res_idx_j][res_idx_i] = self.target_map[res_idx_i][res_idx_j]
            #
            #     self.ref_residues = residues
            #     self.ref_selectors = selectors


            # Initializes the contact/distance map window.
            cp = Contact_map_analysis_window_qt(self.main_window)
            cp.initialize_map(prog=self.main_window,
                              data_array=data_array,
                              ref_residues=self.ref_residues, ##
                              ref_selectors=self.ref_selectors, ##
                              title=title,
                              pixel_size=self.pixel_size,
                              feature_type=self.feature_type,
                              threshold=self.dist_threshold,
                              interaction_center=self.interaction_center,
                              centroid = adj_mat_type)
            cp.show()

            if not has_chain:
                QtWidgets.QMessageBox.warning(self.main_window, "Warning with PDB file", "The provided PDB does not have a valid assigned chain.\nSelection from the interactive contacts map is not possible")




    def get_coords_array(self, objsel, residues_list, interaction_center, get_selectors=True):

        objsel = objsel
        residues = residues_list

        if interaction_center == "CA":
            atm_types = ("CA", )
        elif interaction_center == "CB":
            atm_types = ("CB", )
        elif interaction_center == "centroid":
            atm_types = ("CA", )
        else:
            raise KeyError("Invalid 'interaction_center': %s" % interaction_center)

        # Gets the coordinates from PyMOL.
        from pymol import stored

        # Builds a dictionary in which the keys are tuples with (residue_id, atom_name) and the values
        # are coordinates.
        stored.element_coords = {}

        # The 'cmd.iterate_state' function is particularly fast in PyMOL, much faster than iterating to
        # every residue and calling 'cmd.get_coords' for them.
        cmd.iterate_state(1, objsel, "stored.element_coords[(int(resv), name)] = [x, y, z]")

        # Store the coordinates of residues having an interaction center atom (CA or CB).
        coords = [] # Coordinates list.
        _residues = [] # New residues list, which will substitute the 'residues' one.
        selectors = []
        for res in residues:
            # Attempts to get various atoms.
            for atm_type in atm_types:
                if (res.db_index, atm_type) in stored.element_coords:
                    coords_i = stored.element_coords[(res.db_index, atm_type)]
                    coords.append(coords_i)
                    _residues.append(res)
                    if get_selectors:
                        sel = "object %s and n. %s and i. %s and c. %s" % (objsel, atm_type, res.db_index, res.chain)
                        selectors.append(sel)
                    break

        coords = np.array(coords)

        if get_selectors:
            return _residues, coords, selectors
        else:
            return _residues, coords




class Contact_map_analysis_window_qt(QtWidgets.QMainWindow):
    """
    Class for a window containing a PyQt widget in which a contact/distance map
    will be drawn. Minimal PyQt implementation of a graphical contact map. NumPy
    is required.
    Note: for large protein, building every rectangle in the canvas takes a lot
    of time and there is room for optimization.
    """

    distance_count = 0

    def initialize_map(self, prog, data_array,
                       ref_residues,
                       ref_selectors,
                       title=None,
                       pixel_size=5,
                       feature_type="contact",
                       threshold=8.0,
                       interaction_center="ca",
                       centroid = False):

        # Sets the attributes.
        self.data_array = data_array
        self.pixel_size = pixel_size
        self.feature_type = feature_type
        self.threshold = threshold
        self.interaction_center = interaction_center
        self.ref_residues = ref_residues
        self.ref_selectors = ref_selectors
        if centroid == "centroid":
            self.centroid = True
        else:
            self.centroid = False


        # Assign the methods to get the labels.
        if self.feature_type == "contact":
            self.get_value_label = self._get_value_label_contact
        elif self.feature_type == "probdiff":
            self.get_value_label = self._get_value_label_probdiff
        elif self.feature_type == "distance":
            self.get_value_label = self._get_value_label_distance
        elif self.feature_type == "distances_difference":
            self.get_value_label = self._get_value_label_distance_diff
        elif self.feature_type == "distances_mean":
            self.get_value_label = self._get_value_label_distance_mean
        elif self.feature_type == "distances_std":
            self.get_value_label = self._get_value_label_distance_std
        else:
            raise KeyError(self.feature_type)

        # Set the canvas size.
        min_size = 150
        h = self.pixel_size*len(self.data_array)
        win_size = min((910, h))
        win_size = max((min_size, win_size))

        if title:
            self.setWindowTitle(title)


        # Set some appearance parameters.
        self.controls_padding = 4
        if self.feature_type in ("contact", "distance", "probdiff"):
            self.controls_font = "helvetica 11 bold"
        else:
            self.controls_font = "helvetica 10 bold"
        self.controls_config = {"fg": "black", "font": self.controls_font,
                                "padx": self.controls_padding,
                                "pady": self.controls_padding}
        self.labels_pack_config = {"side": "left", "pady": (0, 5), "padx": (5, 0)}
        self.buttons_pack_config = {"side": "left", "pady": (0, 5), "padx": (1, 0)}


        # Frame of the window containing a row for some control buttons, a row for
        # the plot and a row for a messagebar.
        self.plot_frame = QtWidgets.QWidget()
        self.plot_frame_layout = QtWidgets.QGridLayout()
        self.plot_frame.setLayout(self.plot_frame_layout)
        self.setCentralWidget(self.plot_frame)


        # Control frame.
        self.controls_frame = QtWidgets.QWidget()
        self.controls_frame_layout = QtWidgets.QGridLayout()
        self.controls_frame.setLayout(self.controls_frame_layout)
        self.plot_frame_layout.addWidget(self.controls_frame)

        self.delete_distances_button = QtWidgets.QPushButton("Delete all distances in PyMOL")
        self.delete_distances_button.setEnabled(False)
        self.delete_distances_button.clicked.connect(lambda a=None: self.clear_plot())
        self.controls_frame_layout.addWidget(self.delete_distances_button, 0, 0)

        self.scale_factor = 0
        self.scale_down_button = QtWidgets.QPushButton("Zoom out")
        try:
            self.scale_down_button.setIcon(QtGui.QIcon.fromTheme("go-down"))
        except:
            pass
        self.scale_down_button.clicked.connect(lambda a=None: self.scale_plot_down())
        self.controls_frame_layout.addWidget(self.scale_down_button, 0, 1)

        self.scale_up_button = QtWidgets.QPushButton("Zoom in")
        try:
            self.scale_up_button.setIcon(QtGui.QIcon.fromTheme("go-up"))
        except:
            pass
        self.scale_up_button.clicked.connect(lambda a=None: self.scale_plot_up())
        self.controls_frame_layout.addWidget(self.scale_up_button, 0, 2)

        self.controls_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


        # Frame containing the plot (with a scrollbar).
        self.canvas_plot_frame = QtWidgets.QWidget()
        self.canvas_plot_frame.setStyleSheet("background-color: white")
        self.canvas_plot_frame_layout = QtWidgets.QGridLayout()
        self.canvas_plot_frame.setLayout(self.canvas_plot_frame_layout)

        self.canvas_plot_scrollarea = QtWidgets.QScrollArea()
        self.canvas_plot_scrollarea.setWidgetResizable(True)
        self.canvas_plot_scrollarea.setWidget(self.canvas_plot_frame)
        self.plot_frame_layout.addWidget(self.canvas_plot_scrollarea)


        # Builds the scene where to draw the contact map.
        self.canvas_plot_scene = QtWidgets.QGraphicsScene()
        # Builds the graphics view containing the scene above.
        self.canvas_plot_view = Contact_map_graphics_view(self.canvas_plot_scene)
        self.canvas_plot_frame_layout.addWidget(self.canvas_plot_view)


        # A bottom frame fo the window, containing some buttons to interact with the graph.
        self.message_frame = QtWidgets.QFrame()
        self.message_frame_layout = QtWidgets.QHBoxLayout()
        self.message_frame.setLayout(self.message_frame_layout)
        self.plot_frame_layout.addWidget(self.message_frame)

        # Label to show which residue/position pair is currently being hovered by the mouse pointer.
        if self.feature_type in ("contact", "distance", "probdiff"):
            view_label_text = "Couple:"
        else:
            view_label_text = "Alignment positions:"
        self.view_label = QtWidgets.QLabel(view_label_text)
        # self.view_label.setStyleSheet(self.controls_config)
        self.message_frame_layout.addWidget(self.view_label)


        # Actually draws the contact map.
        self.draw_map()
        # self._test_plot()


    def draw_map(self):
        """
        Methods that actually draw the contact/distance map on a canvas widget.
        """
        w = self.pixel_size

        # Prepare the brushes.
        self.viridis_brushes = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for c in viridis_colors_rev:
                brush = QtGui.QBrush(QtGui.QColor(c[0]*255, c[1]*255, c[2]*255))
                self.viridis_brushes.append(brush)
        self.default_brush = QtGui.QBrush(QtGui.QColor(242, 242, 242))
        self.highlight_brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))

        # Function to get the color of pixels.
        if self.feature_type == "contact":
            self.get_color = self._get_color_contact

        elif self.feature_type == "probdiff":
            self.get_color = self._get_color_probdiff

        elif self.feature_type in ("distance", "distances_mean"):
            self._dist_bins = np.linspace(2.5, self.threshold, 63)
            self.get_color = self._get_color_distance

        elif self.feature_type in ("distances_difference", "distances_std"):
            self._dist_bins = np.linspace(0.0, self.threshold, 63)
            self.get_color = self._get_color_distance

        else:
            raise KeyError(self.feature_type)


        #-------------------------------
        # Draws the map on the canvas. -
        #-------------------------------

        # Use a transparent pen for the border pf the pixels.
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 0), 0)

        for i in range(0, len(self.data_array)):

            for j in range(0, len(self.data_array)):

                if i <= j:

                    # Gets the color brush.
                    color_brush = self.get_color(self.data_array[i][j])

                    # Builds rectangles for both the upper and lower part of the
                    # matrix and adds them to the graphics scene.
                    pid = Contact_map_pixel(j*w, i*w, w, w, i=i, j=j,
                                            contact_map_window=self,
                                            pen=pen, brush=color_brush)
                    self.canvas_plot_scene.addItem(pid)

                    pid = Contact_map_pixel(i*w, j*w, w, w, i=j, j=i,
                                            contact_map_window=self,
                                            pen=pen, brush=color_brush)
                    self.canvas_plot_scene.addItem(pid)


        # Add events to the canvas.
        if self.feature_type in ("contact", "distance", "probdiff"):
            # Draws the map of a single structure on the canvas.
            self.canvas_plot_move_event = self.move_on_plot
            self.canvas_plot_left_click_event = self.click_on_plot

        else:
            # Draws the map of multiple structures on the canvas.
            self.canvas_plot_move_event = self.move_on_plot_ali
            self.canvas_plot_left_click_event = self.click_on_plot_ali


    # Click on the plot.
    def click_on_plot(self, i, j):
        """
        When clicking on a square, draw a distance in the PyMOL viewer between the corresponding
        residues. Used when analyzing a single structure.
        """
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]

        res_name_1 = res_1.three_letter_code + res_1.res_name.split()[1]
        res_name_2 = res_2.three_letter_code + res_2.res_name.split()[1]

        if res_1 is res_2:
            return None

        sel_1 = self.ref_selectors[i]
        sel_2 = self.ref_selectors[j]

        if self.centroid:
            sel_1_res = sel_1.replace(" and n. CA ", "")
            cmd.select("tmp1", sel_1_res)
            com("(tmp1)")
            sel_2_res = sel_2.replace(" and n. CA ", "")
            cmd.select("tmp2", sel_2_res)
            com("(tmp2)")
            cmd.distance("{}__dist_{}_{}".format(self.distance_count, res_name_1, res_name_2), "tmp1_COM", "tmp2_COM")
            cmd.center("{}__dist_{}_{}".format(self.distance_count, res_name_1, res_name_2))
            cmd.delete("tmp1")
            cmd.delete("tmp1_COM")
            cmd.delete("tmp2")
            cmd.delete("tmp2_COM")

        else:
            cmd.distance("{}__dist_{}_{}".format(self.distance_count, res_name_1, res_name_2), sel_1, sel_2)
            cmd.center("{}__dist_{}_{}".format(self.distance_count, res_name_1, res_name_2))

        self.distance_count += 1
        self.delete_distances_button.setEnabled(True)


    def click_on_plot_ali(self, i, j):
        """
        Used when analyzing multiple structures. The residues of the reference structures will be
        used in PyMOL.
        """
        if i == j:
            return None

        sel_1 = self.ref_selectors[i]
        sel_2 = self.ref_selectors[j]
        if sel_1 is None or sel_2 is None:
            return None
        cmd.distance("{}__dist_{}_{}".format(self.distance_count, res_name_1, res_name_2), sel_1, sel_2)
        self.distance_count += 1
        self.delete_distances_button.setEnabled(True)


    # Move the mouse on the plot.
    def move_on_plot(self, i, j):
        """
        Used when showing the contact/distance map of a single structure.
        """

        val = self.data_array[i][j]
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]
        self.view_label.setText("Couple: %s %s - %s %s (%s)." % (res_1.three_letter_code, res_1.db_index,
                                                                 res_2.three_letter_code, res_2.db_index,
                                                                 self.get_value_label(val)))

    def move_on_plot_ali(self, i, j):
        """
        Used when showing the distance map of multiple structures.
        """
        val = self.data_array[i][j]
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]

        label = "Alignment positions: %s - %s (%s)." % (i+1, j+1, self.get_value_label(val))
        if res_1 is not None and res_2 is not None:
            label += " Reference: %s %s - %s %s" % (res_1.three_letter_code, res_1.db_index,
                                                    res_2.three_letter_code, res_2.db_index)
        self.view_label.setText(label)


    # Get the labels to show on the bottom of the window.
    def _get_value_label_contact(self, value):
        if value == 1:
            return "contact"
        else:
            return "non contact"

    def _get_value_label_probdiff(self, value):
        if value == 0:
            return ""
        else:
            return ""

    def _get_value_label_distance(self, value):
        return str(round(value, 2)) + " \u212B"

    def _get_value_label_distance_diff(self, value):
        return self._get_value_label_distance_ali(value, "diff")

    def _get_value_label_distance_mean(self, value):
        return self._get_value_label_distance_ali(value, "mean")

    def _get_value_label_distance_std(self, value):
        return self._get_value_label_distance_ali(value, "std")

    def _get_value_label_distance_ali(self, value, label):
        if value == -1:
            return "not aligned"
        else:
            return "%s = %s %s" % (label, round(value, 2), " \u212B")


    # Get the colors for pixels in the map.
    def _get_color_contact(self, v):
        if v == 0:
            return self.viridis_brushes[-1] # "white"
        else:
            return self.viridis_brushes[0] # "gray"

    def _get_color_probdiff(self, v):
        if v == 0.0:
            return self.viridis_brushes[1] # "white"
        if v > 0:
            return self.viridis_brushes[40] # "gray"
        if v < 0:
            return self.viridis_brushes[55]

    def _get_color_distance(self, v):
        """
        Assignes the bin for the distance value and returns the corresponding color.
        """
        if v == -1:
            return self.default_brush
        return self.viridis_brushes[np.digitize(v, self._dist_bins)]


    # Events influecing the whole plot.
    def clear_plot(self):
        """
        Remove all distances which have been drawn.
        """
        cmd.delete("*__dist_*")
        self.delete_distances_button.setEnabled(False)

    def scale_plot_up(self):
        if self.scale_factor > 10:
            return None
        self.canvas_plot_view.scale(1.25, 1.25)
        self.scale_factor += 1
        if self.scale_factor > 10:
            self.scale_up_button.setEnabled(False)
        self.scale_down_button.setEnabled(True)

    def scale_plot_down(self):
        if self.scale_factor < -10:
            return None
        self.canvas_plot_view.scale(0.8, 0.8)
        self.scale_factor -=1
        if self.scale_factor < -10:
            self.scale_down_button.setEnabled(False)
        self.scale_up_button.setEnabled(True)


###############################################################################
# Results window.                                                             #
###############################################################################

class Contact_map_graphics_view(QtWidgets.QGraphicsView):
    pass


class Contact_map_pixel(QtWidgets.QGraphicsRectItem):
    """
    Custom class for drawing on a graphic scene of PyQt rectangles corresponding
    to pixels of a contact/distance map.
    """

    def __init__(self, x, y, w, h, contact_map_window, pen, brush, i, j):

        super(Contact_map_pixel, self).__init__(x, y, w, h)
        self.contact_map_window = contact_map_window
        self.setPen(pen)
        self.original_brush = brush
        self.setBrush(brush)
        self.setAcceptHoverEvents(True)
        self.i_idx = i
        self.j_idx = j


    def mousePressEvent(self, e):
        self.contact_map_window.canvas_plot_left_click_event(self.i_idx, self.j_idx)

    def hoverEnterEvent(self, e):
        self.setBrush(self.contact_map_window.highlight_brush)
        self.contact_map_window.canvas_plot_move_event(self.i_idx, self.j_idx)

    def hoverLeaveEvent(self, e):
        self.setBrush(self.original_brush)

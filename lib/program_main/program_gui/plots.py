import sys

from PyQt5.QtWidgets import QDialog, QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    import matplotlib.pyplot as plt
except:
    pass


import random
from pymol import cmd


class PlotCentrality():

    def __init__(self, tab, main):

        self.tab = tab
        self.main = main

    def plot_window(self):

        self.window = Window(self.main)

    def add_plot(self, x_label, y_label, title,
                original_data,
                type = "",
                line_plot_data = [],
                scatter_plot_data = [],
                algorithm = "",
                pdb_name = "",
                cmap_scatter = []):

        self.alg = algorithm

        if type == "line_plot":

            self.window.line_plot(original_data=original_data,
                             x_list = line_plot_data,
                             y_list = [],
                             alg = self.alg,
                             x_label = x_label,
                             y_label = y_label,
                             title = title,
                             pdb_name = pdb_name)

        if type == "scatter_plot":
            self.window.scatter_plot(original_data = original_data,
                             x_list = scatter_plot_data[0],
                             y_list = scatter_plot_data[1],
                             alg = self.alg,
                             x_label = x_label,
                             y_label = y_label,
                             title = title,
                             pdb_name = pdb_name,
                             cmap_scatter = cmap_scatter)

    def show_plot(self):
        self.window.show()


class Window(QDialog):

    # constructor
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that
        # displays the 'figure'it takes the
        # 'figure' instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('button_press_event', self.onclick_lineplot)
        self.canvas.mpl_connect('motion_notify_event', self.hover_on_plot)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # # adding action to the button
        # self.button.clicked.connect(self.plot)

        # creating a Vertical Box layout
        layout = QVBoxLayout()

        # adding tool bar to the layout
        layout.addWidget(self.toolbar)

        # adding canvas to the layout
        layout.addWidget(self.canvas)

        # setting layout to the main window
        self.setLayout(layout)

        # clearing old figure
        self.figure.clear()

        # create an axis
        self.ax = self.figure.add_subplot(111)

        self.all_subplots = []

        self.original_data = {}

        plt.rcParams['font.size'] = 10
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        annotation = True

        if annotation:
            self.annotation = self.ax.annotate(text = '',
                             xy = (0,0),
                             xytext=(15,15),
                             textcoords='offset points',
                             bbox = {'boxstyle': 'round', 'fc': 'w'},
                             arrowprops = {'arrowstyle': '->'}
                             )

            self.annotation.set_visible(False)
        else:
            self.annotation = False


    def onclick_lineplot(self, event):

        if self.annotation:

            if event.inaxes == self.ax:
                for subplot in self.all_subplots:
                    if type(subplot) is list: # Line Plot
                        plot_instance = subplot[0] # e.g. plot_instance = Line2D(1OL5:unnorm_ssc)
                        plot_instance_name = str(plot_instance).split("(")[1].replace(")", "")
                        pdb_name = (str(plot_instance).split("(")[1].replace(")", "")).split(":")[0]

                        is_contained, annotation_index = plot_instance.contains(event)

                        if is_contained:
                            posx, posy = [plot_instance.get_xdata()[annotation_index['ind'][0]], plot_instance.get_ydata()[annotation_index['ind'][0]]]
                            text_label = str(list(self.original_data[plot_instance_name].keys())[int(posx)])
                            res = text_label.split()[0][:3]
                            res_num = text_label.split()[0][3:]
                            chain = text_label.split()[1]
                            sele_name = "{}{}{}".format(res, res_num, chain)
                            self.pymol_sele_list.append(sele_name)
                            cmd.select(pdb_name + "_sele", pdb_name)
                            cmd.select(sele_name, "{}/{}`{}/ & {}".format(chain, res, res_num, pdb_name + "_sele"))
                            cmd.show("dots", sele_name)

                # #self.show_in_pymol()
                # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                #       ('double' if event.dblclick else 'single', event.button,
                #        event.x, event.y, event.xdata, event.ydata))


    def hover_on_plot(self, event):

        if self.annotation:

            annotation_visible = self.annotation.get_visible()

            if event.inaxes == self.ax:
                for subplot in self.all_subplots:
                    if type(subplot) is list: # Line Plot
                        plot_instance = subplot[0] # e.g. plot_instance = Line2D(1OL5:unnorm_ssc)
                        plot_instance_name = str(plot_instance).split("(")[1].replace(")", "")

                        is_contained, annotation_index = plot_instance.contains(event)

                        if is_contained:
                            posx, posy = [plot_instance.get_xdata()[annotation_index['ind'][0]], plot_instance.get_ydata()[annotation_index['ind'][0]]]
                            self.annotation.xy = (posx, posy)
                            text_label = str(list(self.original_data[plot_instance_name].keys())[int(posx)])
                            self.annotation.set_text(text_label)
                            self.annotation.set_visible(True)
                            self.figure.canvas.draw_idle()
                        else:
                            if annotation_visible:
                                self.annotation.set_visible(False)
                                self.figure.canvas.draw_idle()

                    ## TODO -  also with scatter plot
                    # else:
                    #     # print(self.plot_instance[0])
                    #     # plot_instance_name = str(self.plot_instance).split("(")[1].replace(")", "")
                    #     is_contained, annotation_index = subplot.contains(event)
                    #
                    #     if is_contained: # Scatter Plot
                    #         data_point_location = subplot.get_offsets()[annotation_index['ind'][0]]
                    #         self.annotation.xy = data_point_location
                    #         #text_label = str(list(self.original_data[plot_instance_name].keys())[int(data_point_location[0])])
                    #         text_label = ""
                    #         self.annotation.set_text(text_label)
                    #         self.annotation.set_visible(True)
                    #         self.figure.canvas.draw_idle()
                    #
                    #     else:
                    #         if annotation_visible:
                    #             self.annotation.set_visible(False)
                    #             self.figure.canvas.draw_idle()


    def scatter_plot(self, x_list, y_list, alg, original_data, x_label = "", y_label = "", title = "", pdb_name = "", cmap_scatter = []):

        self.original_data["{}:{}".format(pdb_name.upper(), alg)] = original_data

        # plot data
        if cmap_scatter:
            self.plot_instance = self.ax.scatter(x_list, y_list, c = cmap_scatter, alpha = 0.6, s=10)
        else:
            self.plot_instance = self.ax.scatter(x_list, y_list, label = "{}:{}".format(pdb_name.upper(), alg))

        self.all_subplots.append(self.plot_instance)
        self.ax.set_ylabel(y_label, fontsize=12)
        self.ax.set_xlabel(x_label, fontsize=12)
        self.ax.legend(fontsize=12)
        self.ax.set_title(title, fontsize=12)
        # refresh canvas
        self.canvas.draw()


    def line_plot(self, x_list, y_list, alg, original_data, x_label = "", y_label = "", title = "", pdb_name = ""):

        # Create list of PyMOL selection
        self.pymol_sele_list = []

        self.original_data["{}:{}".format(pdb_name.upper(), alg)] = original_data

        # plot data
        self.plot_instance = self.ax.plot(x_list, label = "{}:{}".format(pdb_name.upper(), alg), marker='o')
        self.all_subplots.append(self.plot_instance)
        self.ax.set_ylabel(y_label, fontsize=12)
        self.ax.set_xlabel(x_label, fontsize=12)
        self.ax.legend(fontsize=12)
        self.ax.set_title(title, fontsize=12)

        # refresh canvas
        self.canvas.draw()

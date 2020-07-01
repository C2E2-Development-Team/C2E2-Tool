from bokeh.io import save, export_png
from bokeh.plotting import figure

from frontend.mod.constants import *
from frontend.mod.session import Session


def plot_graph(data_filepath, horizontal_index, vertical_indices, 
    variable_list, mode_list, type_, title, filename):
    """
    File Format for input files:
        Line 0 - 10: Header information
        Line 11+: Simulation/Verification Data

    data_filepath:      File path containing information to plot
    horizontal_index:   Index of horizontal axis variable in variable_list
    vertical_indices:   List of indices of variables in variable_list to be 
                            plotted on the vertical axis
    type_:              VERIFIED or SIMULATED
                            VERIFIED = Reach Tube plotting
                            SIMULATED = Line plotting
    variable_list:      List of variables
    mode_list:          List of modes
    title:              Plot title
    filename:           Output file name. 
                            filename.png and filename.html will be output
    """

    # Load data from file into 2D list
    file_object = open(data_filepath, 'r')
    lines = file_object.readlines()
    data_points = []
    for i in range(11, len(lines)):  # Output starts on line 11
        points = lines[i].split()
        if points[0] == '%':
            continue
        else:
            data_points.append(points)

    if type_ == SIMULATED:
        plot_line(data_points, horizontal_index, vertical_indices, 
            variable_list, mode_list, title, filename)
    else:
        plot_quad(data_points, horizontal_index, vertical_indices, 
            variable_list, mode_list, title, filename)

def plot_line(data_points, horizontal_index, vertical_indices, variable_list, 
    mode_list, title, filename):

    Session.write("Generating simulation plot... ")

    bokeh_plot = figure(title=title)

    x_axis = []
    y_axes = []
    for i in range(len(vertical_indices)):
        y_axes.append([])
    
    toggle = True
    plot_index = -1
    for line in data_points:

        if float(line[0]) == 0:
            plot_index += 1
            x_axis.append([])
            for y_index in range(len(vertical_indices)):
                y_axes[y_index].append([])
            toggle = True

        if toggle:
            x_axis[plot_index].append(float(line[horizontal_index]))
            for y_index, vert_index in enumerate(vertical_indices):
                y_axes[y_index][plot_index].append(float(line[vert_index]))
            toggle = False
        else:
            toggle = True

    for i in range(plot_index+1):
        for j, y_axis in enumerate(y_axes):
            bokeh_plot.line(x_axis[i], y_axis[i], line_width=2, 
                color=PLOT_COLORS[j-1], 
                legend_label=variable_list[vertical_indices[j]])

    # X - axis
    bokeh_plot.xaxis.axis_label = variable_list[horizontal_index]
    # Y - axis
    y_axis_label = variable_list[vertical_indices[0]]
    for i in range(1, len(vertical_indices)):
        y_axis_label += ", " + variable_list[vertical_indices[i]]
    bokeh_plot.yaxis.axis_label = y_axis_label

    save(bokeh_plot, filename=filename+'.html', title=title)
    export_png(bokeh_plot, filename=filename+'.png' )

    Session.write("Done.\n")

def plot_quad(data_points, horizontal_index, vertical_indices, variable_list, 
    mode_list, title, filename):

    Session.write("Generating verification plot...")

    bokeh_plot = figure(title=title)

    for i, vertical_index in enumerate(vertical_indices):
        top = []
        bottom = []
        left = []
        right = []

        for j in range(1, len(data_points), 2):
            # Left and Right side defined by horizontal axis variable
            left.append(float(data_points[j-1][horizontal_index]))
            right.append(float(data_points[j][horizontal_index]))
            # Top and Bottom are defined by vertical axis variable
            bottom.append(min(float(data_points[j-1][vertical_index]),
                float(data_points[j][vertical_index])))
            top.append(max(float(data_points[j-1][vertical_index]),
                float(data_points[j][vertical_index])))

        bokeh_plot.quad(left=left, right=right, bottom=bottom, top=top, 
            color=PLOT_COLORS[i], line_width=0, line_alpha=0.35, fill_alpha=0.35,
            legend_label=variable_list[vertical_indices[i]])

    # X - axis
    bokeh_plot.xaxis.axis_label = variable_list[horizontal_index]
    # Y - axis
    y_axis_label = variable_list[vertical_indices[0]]
    for i in range(1, len(vertical_indices)):
        y_axis_label += ", " + variable_list[vertical_indices[i]]
    bokeh_plot.yaxis.axis_label = y_axis_label

    save(bokeh_plot, filename=filename + '.html', title=title)
    export_png(bokeh_plot, filename=filename + '.png')

    Session.write("Done.\n")

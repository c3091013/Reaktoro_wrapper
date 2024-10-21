# Reaktoro_wrapper
wrappers around the reaktoro code that make generating figures easier.

wrapper(lst)
wrapper to allow parralelisation of functions used in this script. input is a list with the name of the function followed by a list of the inputs
lst[0] = name of function to be used
lst[1,:] = list of typical input for that function
returns the output of the function. if no known function is used it just returns the string "missed"

get_1D_plot(ax,data,xlim,imask,title,xtit,ytit,colors)
generates a 1D plot of data and casts it to an axis object from matplotlib
ax = an axis object the data will be attached too
data = the data to be plotted
xlim = limits of the x axis
imask = list of species to be plotted, used hide species not being targeteed that are included in the data.
title = axis title
xtit = x axis title
ytit = y axis title
colors = dictionary of rgb colors associated with each species the imask.
returns the axis object after modification.

get_2D_plot(ax,data,xlim,ylim,imask,title,xtit,ytit,colors)
generates a 2D patches plotfor a given set of predominance data and casts it to an axis object.
ax = an axis object the data will be attached too
data = the data to be plotted
xlim = limits of the x axis
ylim = limits of the y axis
imask = list of species to be plotted, used hide species not being targeteed that are included in the data.
title = axis title
xtit = x axis title
ytit = y axis title
colors = dictionary of rgb colors associated with each species the imask.
returns the axis object after modification.

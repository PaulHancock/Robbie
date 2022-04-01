#! /usr/bin/env python3

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh import palettes
from bokeh.models import (ColumnDataSource, Circle,
                          DataTable, TableColumn, NumberFormatter)
from bokeh.layouts import gridplot
import pandas as pd
import numpy as np
from astropy.io.votable import parse


def load_mean_image(filename):
    hdu = fits.open(filename)#'results/mean_image.fits')
    return hdu

def mean_image_data(hdu):
    wcs = WCS(hdu[0].header)
    x,y = np.indices(hdu[0].data.shape)
    sky_map = wcs.pixel_to_world(x.ravel(),y.ravel())
    ra = sky_map.ra.degree.reshape(x.shape)
    dec = sky_map.dec.degree.reshape(x.shape)
    return hdu[0].data, ra, dec

def get_imdata(data, ra, dec):
    imdata ={'image':[data],
        'ra':[ra],
        'dec':[dec],
        'x':[ra[0,0]],
        'y':[dec[0,0]],
        'dw':[abs(ra[-1,-1]-ra[0,0])],
        'dh':[abs(dec[-1,-1]-dec[0,0])]
        }
    return imdata


def get_data(fname):
    tab = parse(fname).get_first_table().to_table()
    df = tab.to_pandas()
    return df


def get_joined_table_source():
    # most likely there is a solution that doesn't involve so many intermediaries!
    flux_table = get_data("/home/paulhancock/alpha/ADACS/MAP22A-TGalvin/Robbie-ADACS/results/flux_table.vot")
    stats_table = get_data("/home/paulhancock/alpha/ADACS/MAP22A-TGalvin/Robbie-ADACS/results/stats_table.vot")
    joined = flux_table.join(stats_table.set_index('uuid'), on='uuid')
    df = joined
    source = ColumnDataSource(data=dict( (i,df[i]) for i in df.columns))
    return source

def get_scatter_plots():
    selected_circle = Circle(fill_alpha=1, fill_color="firebrick", line_color=None)
    nonselected_circle = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)


    source = get_joined_table_source()
    TOOLS = "box_select,lasso_select,help,pan,tap,wheel_zoom"

    # create a new plot and add a renderer
    left = figure(tools=TOOLS, width=300, height=300, title=None)
    lr = left.circle('ref_ra', 'ref_dec', source=source)

    # create another new plot and add a renderer
    right = figure(tools=TOOLS, width=300, height=300, title=None)
    rr = right.circle('md', 'pval_peak_flux_ks', source=source)

    for r in [lr,rr]:
        r.selection_glyph = selected_circle
        r.nonselection_glyph = nonselected_circle

    #p = gridplot([[left, right]])
    columns = [
        TableColumn(field="uuid", title="UUID"),
        TableColumn(field="ref_ra", title="RA", formatter=NumberFormatter(format="0.0000")),
        TableColumn(field="ref_dec", title="DEC", formatter=NumberFormatter(format="0.0000")),
        TableColumn(field="md", title="Debiased modulation index"),
        TableColumn(field="pval_peak_flux_ks", title="PVal")
    ]
    data_table = DataTable(source=source, columns=columns, scroll_to_selection=True, sortable=True, background='#111111')
    return left,right,data_table


def get_mean_image_plot():
    hdu = load_mean_image('results/mean_image.fits')
    data, ra, dec = mean_image_data(hdu)
    imdata = get_imdata(data,ra,dec)

    p = figure(tooltips=[
                            ("value", "@image Jy/beam"),
                            ("RA", "@ra{0.00}°"),
                            ("DEC", "@dec{0.00}°")
                        ],
            x_axis_label='RA',
            y_axis_label='DEC')
    p.x_range.range_padding = p.y_range.range_padding = 0

    # must give a vector of image data for image parameter
    p.image(source=imdata, 
            image='image',
            palette=palettes.mpl['Cividis'][256])#, level="image")
    #p.grid.grid_line_width = 0.5
    return p

def main():
    output_file(filename='viewer/RV.html', title="Robbie Viewer")
    mean_image = get_mean_image_plot()
    scatter_plots = get_scatter_plots()
    p = gridplot([ scatter_plots,[mean_image]])
    show(p)


if __name__ == "__main__":
    main()
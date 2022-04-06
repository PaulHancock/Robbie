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

TOOLS = "box_select,lasso_select,help,pan,tap,wheel_zoom,reset"

def load_mean_image(filename):
    hdu = fits.open(filename)#'results/mean_image.fits')
    return hdu

def mean_image_data(hdu):
    wcs = WCS(hdu[0].header)
    data = np.fliplr(hdu[0].data)
    x,y = np.indices(data.shape[::-1])
    sky_map = wcs.pixel_to_world(x.reshape(-1),y.reshape(-1))
    ra = np.fliplr(sky_map.ra.degree.reshape(x.shape).T)
    dec = sky_map.dec.degree.reshape(y.shape).T
    return data, ra, dec

def get_imdata(data, ra, dec):
    imdata ={'image':[data],
        'ra':[ra],
        'dec':[dec],
        'x':[ra[0,0]],
        'y':[dec[0,0]],
        'dw':[abs(ra[0,0]-ra[0,-1])],
        'dh':[abs(dec[0,0]-dec[-1,0])]
        }
    return imdata


def get_tabdata(fname):
    tab = parse(fname).get_first_table().to_table()
    df = tab.to_pandas()
    return df


def get_joined_table_source():
    # most likely there is a solution that doesn't involve so many intermediaries!
    flux_table = get_tabdata("/home/paulhancock/alpha/ADACS/MAP22A-TGalvin/Robbie-ADACS/results/flux_table.vot")
    stats_table = get_tabdata("/home/paulhancock/alpha/ADACS/MAP22A-TGalvin/Robbie-ADACS/results/stats_table.vot")
    joined = flux_table.join(stats_table.set_index('uuid'), on='uuid')
    df = joined
    source = ColumnDataSource(data=dict( (i,df[i]) for i in df.columns))
    return source

def get_scatter_plots(source):
    selected_circle = Circle(fill_alpha=1, fill_color="firebrick", line_color=None)
    nonselected_circle = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)

    

    # create a new plot and add a renderer
    left = figure(tools=TOOLS, 
                  #width=300, height=300,
                  title='Sky Plot',
                  x_axis_label='RA (deg)', y_axis_label='DEC (deg)')
    lr = left.circle('ref_ra', 'ref_dec', source=source)

    # create another new plot and add a renderer
    right = figure(tools=TOOLS, 
                   #width=300, height=300,
                   title='Variables Plot',
                   x_axis_label='md', y_axis_label='pval_peak_flux_ks')
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


def get_mean_image_plot(source):
    hdu = load_mean_image('results/mean_image.fits')
    data, ra, dec = mean_image_data(hdu)
    imdata = get_imdata(data,ra,dec)

    p = figure(tools=TOOLS,
               tooltips=[
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
    p.circle(source=source, 
             x='ref_ra', y='ref_dec',
             radius=0.03, fill_color=None,
             line_width=1.5, line_color='yellow')
    #p.grid.grid_line_width = 0.5
    return p

def main():
    output_file(filename='viewer/RV.html', title="Robbie Viewer")
    source = get_joined_table_source()
    mean_image = get_mean_image_plot(source)
    scatter_plots = get_scatter_plots(source)
    p = gridplot([ scatter_plots,[mean_image]])
    show(p)


if __name__ == "__main__":
    main()
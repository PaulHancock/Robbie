from functools import lru_cache
from os import stat
from astropy.io import fits
from astropy.io.votable import parse
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

import astropy.units as u

import numpy as np
import datetime

from bokeh.plotting import figure, show, output_file
from bokeh import palettes
from bokeh.layouts import gridplot, layout
from bokeh.io import curdoc
from bokeh.models import (ColumnDataSource, Circle, Whisker,
                          DataTable, TableColumn, NumberFormatter,
                          CustomJS, Slider, DatetimeTickFormatter,
                          HoverTool)
TOOLS = "box_select,lasso_select,help,pan,tap,wheel_zoom,reset"
selected_circle = Circle(fill_alpha=1, fill_color="firebrick", line_color=None)
nonselected_circle = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)
    
# Testing
# ref_ra_min, ref_ra_max = 335, 340
# ref_dec_min, ref_dec_max = -35, -30

cutout_position = SkyCoord(335*u.deg, -15*u.deg, frame='icrs')

# x_lim, y_lim = ((ref_ra_min, ref_ra_max), (ref_dec_min, ref_dec_max))


@lru_cache
def load_mean_image(filename):
    hdu = fits.open(filename)#'results/mean_image.fits')
    hdu[0].data = hdu[0].data
    return hdu

def mean_image_data(hdu):
    wcs = WCS(hdu[0].header)
    cutout = Cutout2D(hdu[0].data, cutout_position, 20*u.deg, wcs=wcs)
    wcs = cutout.wcs
    data = np.fliplr(cutout.data)   
    x,y = np.indices(data.shape[::-1])
    sky_map = wcs.pixel_to_world(x.reshape(-1),y.reshape(-1))
    ra = np.fliplr(sky_map.ra.degree.reshape(x.shape).T)
    dec = sky_map.dec.degree.reshape(y.shape).T

    # Decimation
    # data, ra, dec = data[::5, ::5], ra[::5, ::5], dec[::5, ::5]

    # Create mask array based on provided percentage of data to cut off
    # ra_mask = (ra < x_lim[-1]) & (ra > x_lim[0])
    # dec_mask = (dec < y_lim[-1]) & (dec > y_lim[0])
    # mask = ra_mask * dec_mask
    # sqrt_size = int(np.sqrt(mask.flatten().sum()))
    

    # # Try to reshape based on the sqrt to ensure a square array
    # try:        
    #     data, ra, dec  = data[mask].reshape((sqrt_size, sqrt_size)), ra[mask].reshape((sqrt_size, sqrt_size)), dec[mask].reshape((sqrt_size, sqrt_size))
    # except:
    #     # Number of extra elements needed to mask to ensure a square matrix
    #     mask_diff =  mask.flatten().sum() - sqrt_size**2
    #     # Make True elements False based on the number of mask diff
    #     new_mask = mask.copy().flatten()
    #     true_idxs = np.where(new_mask == True)[0][-mask_diff:]
    #     new_mask[true_idxs] = False
    #     # Reshape data back to square 2D array
    #     new_mask = new_mask.reshape((mask.shape[0], mask.shape[-1]))
    #     data, ra, dec  = data[new_mask].reshape((sqrt_size, sqrt_size)), ra[new_mask].reshape((sqrt_size, sqrt_size)), dec[new_mask].reshape((sqrt_size, sqrt_size))
    #     print(f"Actual ra range: {ra.min()} to {ra.max()} degrees")
    #     print(f"Actual dec range: {dec.min()} to {dec.max()} degrees")

    return data, ra, dec

def get_imdata(data, ra, dec, ei=None):
    if ra[0][0] > ra[0][-1]:
        # needs inverting
        ra_new = []
        for ra_row in ra:
            ra_row = list(ra_row)
            ra_row.reverse()
            ra_new.append(ra_row)
        ra = np.array(ra_new)
        data_new = []
        for data_row in data:
            data_row = list(data_row)
            data_row.reverse()
            data_new.append(data_row)
        data = np.array(data_new)
    if ei is None:
        imdata ={
            'image':[data],
            'ra':[ra],
            'dec':[dec],
            'x':[ra[0,0]],
            'y':[dec[0,0]],
            'dw':[abs(ra[0,0]-ra[0,-1])],
            'dh':[abs(dec[0,0]-dec[-1,0])]
        }
    else:
        imdata ={
            f'image_{ei}':[data],
            f'ra_{ei}':[ra],
            f'dec_{ei}':[dec],
            f'x_{ei}':[ra[0,0]],
            f'y_{ei}':[dec[0,0]],
            f'dw_{ei}':[abs(ra[0,0]-ra[0,-1])],
            f'dh_{ei}':[abs(dec[0,0]-dec[-1,0])]
        }
    return imdata

def get_tabdata(fname):
    tab = parse(fname).get_first_table().to_table()
    df = tab.to_pandas()
    return df

@lru_cache
def get_joined_table_source(result_dir):
    # most likely there is a solution that doesn't involve so many intermediaries!
    flux_table = get_tabdata(f"{result_dir}/flux_table.vot")
    stats_table = get_tabdata(f"{result_dir}/stats_table.vot")

    # flux_table = flux_table[::10]
    # stats_table = stats_table[::10]

    # # Mask data based on provided RA and DEC ranges
    # flux_table = flux_table.drop(flux_table[(flux_table['ref_ra'] > x_lim[-1]) + (flux_table['ref_ra'] < x_lim[0])].index)
    # flux_table = flux_table.drop(flux_table[(flux_table['ref_dec'] > y_lim[-1]) + (flux_table['ref_dec'] < y_lim[0])].index)
   
    # # Drop NaN
    flux_table = flux_table.dropna()
    stats_table = stats_table.dropna()

    joined = flux_table.join(stats_table.set_index('uuid'), on='uuid')    
    df = joined
    source = ColumnDataSource(data=dict( (i,df[i]) for i in df.columns))
    
    flux_cols = list(sorted([i for i in df.columns if i.startswith('peak_flux_')]))
    err_flux_cols = [f'err_{i}' for i in flux_cols]
    fluxes = df.set_index('uuid')[flux_cols]
    err_fluxes = df.set_index('uuid')[err_flux_cols]
    fluxes = fluxes.T
    err_fluxes = err_fluxes.fillna(0).T

    lc = dict( (i,fluxes[i].values) for i in fluxes.columns)
    err_lc = dict( (f'err_{i}',err_fluxes[i].values) for i in err_fluxes.columns)
    lc.update(err_lc)

    epochs = [int(e.split('_')[-1]) for e in flux_cols]
    dates = [ df[e].unique()[0] for e in df.columns if e.startswith('epoch_')]
    datetimes = [datetime.datetime.strptime(a, "%Y-%m-%dT%H:%M:%S") for a in dates]
    lc.update({
        'epoch':epochs,
        'date':dates,
        'datetimes':datetimes,
        'current':fluxes[df['uuid'].iloc[0]],
        'current_upper':fluxes[df['uuid'].iloc[0]].values+err_fluxes[df['uuid'].iloc[0]].values,
        'current_lower':fluxes[df['uuid'].iloc[0]].values-err_fluxes[df['uuid'].iloc[0]].values,
    })

    lc_source = ColumnDataSource(data=lc)
    return source, lc_source

def get_scatter_plots(source):

    # create a new plot and add a renderer
    left = figure(
        tools=TOOLS,
        title='Sky Plot',
        x_axis_label='RA (deg)', y_axis_label='DEC (deg)',
        x_range=(0,1),
        y_range=(0,1),
    )
    lr = left.circle('ref_ra', 'ref_dec', source=source)

    # create another new plot and add a renderer
    right = figure(
        tools=TOOLS,
        title='Variables Plot',
        x_axis_label='md', y_axis_label='pval_peak_flux_ks',
    )
    rr = right.circle('md', 'pval_peak_flux_ks', source=source)

    for r in [lr,rr]:
        r.selection_glyph = selected_circle
        r.nonselection_glyph = nonselected_circle

    columns = [
        TableColumn(field="uuid", title="UUID"),
        TableColumn(field="ref_ra", title="RA", formatter=NumberFormatter(format="0.0000")),
        TableColumn(field="ref_dec", title="DEC", formatter=NumberFormatter(format="0.0000")),
        TableColumn(field="mean_peak_flux", title="μ(flux)"),
        TableColumn(field="std_peak_flux", title="σ(flux)"),
        TableColumn(field="md", title="Debiased modulation index"),
        TableColumn(field="pval_peak_flux_ks", title="PVal")
    ]
    data_table = DataTable(source=source, columns=columns, scroll_to_selection=True, sortable=True, background='#111111')
    return left,right,data_table

def get_mean_image_plot(source, result_dir):
    hdu = load_mean_image(f"{result_dir}/mean_image_reprojection_test.fits")
    data, ra, dec = mean_image_data(hdu)
    imdata = get_imdata(data,ra,dec)

    p = figure(
        tools=TOOLS,
        title='Mean Image',
        tooltips=[
            ("value", "@image Jy/beam"),
            ("RA", "@ra{0.00}°"),
            ("DEC", "@dec{0.00}°")
        ],
        x_axis_label='RA',
        y_axis_label='DEC',
        x_range=(min(source.data['ref_ra']),  max(source.data['ref_ra'])),
        y_range=(min(source.data['ref_dec']), max(source.data['ref_dec'])),
    )

    # Adjust hover tool to only work on image data
    hover_tool = p.select(type=HoverTool)
    hover_tool.names = ["mean_image"]

    # must give a vector of image data for image parameter
    p.image(source=imdata,
            image='image',
            palette=palettes.mpl['Cividis'][256], name="mean_image")#, level="image")
    p.circle(source=source,
             x='ref_ra', y='ref_dec',
             radius=0.03, fill_color=None,
             line_width=1.5, line_color='yellow')
    #p.grid.grid_line_width = 0.5
    return p

def get_epoch_image_plots(source, epoch_source, num_epochs):
    # Set up the figure
    epochs = figure(
        tools=TOOLS,
        title="Epoch_0",
        tooltips=[
            ("value", "@image Jy/beam"),
            ("RA", "@ra{0.00}°"),
            ("DEC", "@dec{0.00}°")
        ],
        x_axis_label='RA',
        y_axis_label='DEC',
        x_range=(min(source.data['ref_ra']),  max(source.data['ref_ra'])),
        y_range=(min(source.data['ref_dec']), max(source.data['ref_dec'])),
    )



    # Adjust hover tool to only work on image data
    hover_tool = epochs.select(type=HoverTool)
    hover_tool.names = ["epoch_image"]
    epoch_slider = Slider(start=0, end=num_epochs-1, value=0, step=1, title="Epoch")

    # must give a vector of image data for image parameter
    epochs.image(source=epoch_source,
            palette=palettes.mpl['Cividis'][256], name="epoch_image")


    # Add the source circles
    epochs.circle(source=source,
                x='ref_ra', y='ref_dec',
                radius=0.03, fill_color=None,
                line_width=1.5, line_color='yellow')

    return epochs, epoch_slider

def get_light_curve_plot(source):
    tooltips = [("epoch","@date"),]
    lc_plot = figure(
        tools=TOOLS,
        title='Light curve:',
        tooltips=tooltips,
        x_axis_type="datetime",
        x_axis_label='Epoch',
        y_axis_label='Peak flux (Jy/beam)',
        y_range=(0,1),
    )
    lc_plot.xaxis[0].formatter = DatetimeTickFormatter(
        years = ["%Yy"],
        months = ["%b"],
        days = ["%dd"],
        hours=["%Hh"],
        hourmin=['%Hh %Mmin'],
        minutes=["%Mmin"],
        minsec = ['%Mmin %Ss'],
        seconds=["%Ss"],
    )
    lines = lc_plot.line(source=source, y='current', x='datetimes')
    circles = lc_plot.circle(source=source, y='current', x='datetimes', size=6)
    whiskers = lc_plot.add_layout(Whisker(source=source, base='datetimes', upper='current_upper', lower='current_lower', line_color='gray'))
    for r in [circles,]:
        r.selection_glyph = selected_circle
        r.nonselection_glyph = nonselected_circle

    return lc_plot
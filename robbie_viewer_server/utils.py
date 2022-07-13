from copy import copy
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

@lru_cache
def load_mean_image(filename):
    hdu = fits.open(filename)#'results/mean_image.fits')    
    return hdu

def mean_image_data(hdu, ra_ref=None, dec_ref=None, degrees_around_ref_coords=None):
    wcs = WCS(hdu[0].header)  
    data = np.fliplr(hdu[0].data)
    if ra_ref is not None:
        cutout_position = SkyCoord(ra_ref*u.deg, dec_ref*u.deg, frame='icrs')
        cutout = Cutout2D(hdu[0].data, cutout_position, degrees_around_ref_coords*u.deg, wcs=wcs, mode='partial')
        wcs = cutout.wcs
        data = np.fliplr(cutout.data)
        
    x,y = np.indices(data.shape[::-1])
    sky_map = wcs.pixel_to_world(x.reshape(-1),y.reshape(-1))
    ra = np.fliplr(sky_map.ra.degree.reshape(x.shape).T)
    dec = sky_map.dec.degree.reshape(y.shape).T
    # As we have reprojected we only need ra and dec extents
    ra_dec = np.array([[ra[0,0], ra[0,-1]], [dec[0][0], dec[-1][0]]])
    return data, ra_dec

def get_imdata(data, ra_dec, ei=None):
    if ra_dec[0][0] > ra_dec[0][-1]:
        # needs inverting
        ra_new = []
        ra_dec_row = list(ra_dec[0])
        ra_dec_row.reverse()
        ra_new.append(ra_dec_row)
        ra_dec[0] = np.array(ra_new)
        data_new = []
        for data_row in data:
            data_row = list(data_row)
            data_row.reverse()
            data_new.append(data_row)
        data = np.array(data_new)
    if ei is None:
        imdata ={
            'image':[data],
            # 'ra':[ra],
            # 'dec':[dec],
            'x':[ra_dec[0,0]],
            'y':[ra_dec[-1,0]],
            'dw':[abs(ra_dec[0,0]-ra_dec[0,-1])],
            'dh':[abs(ra_dec[-1,0]-ra_dec[-1,-1])]
        }
    else:
        imdata ={
            f'image_{ei}':[data],
            # f'ra_{ei}':[ra],
            # f'dec_{ei}':[dec],
            f'x_{ei}':[ra_dec[0,0]],
            f'y_{ei}':[ra_dec[-1,0]],
            f'dw_{ei}':[abs(ra_dec[0,0]-ra_dec[0,-1])],
            f'dh_{ei}':[abs(ra_dec[-1,0]-ra_dec[-1,-1])]
        }
    return imdata

def get_tabdata(fname):
    tab = parse(fname).get_first_table().to_table()
    df = tab.to_pandas()
    return df

@lru_cache
def get_joined_table_source(result_dir, ra_ref=None, dec_ref=None, degrees_around_ref_coords=None):
    # most likely there is a solution that doesn't involve so many intermediaries!
    flux_table = get_tabdata(f"{result_dir}/flux_table.vot")
    stats_table = get_tabdata(f"{result_dir}/stats_table.vot")

    # Mask data based on provided RA and DEC ranges
    if ra_ref is not None:
        dec_min, dec_max = dec_ref - degrees_around_ref_coords/2, dec_ref + degrees_around_ref_coords/2
        ra_min, ra_max = ra_ref - degrees_around_ref_coords/2, ra_ref + degrees_around_ref_coords/2
        flux_table = flux_table.drop(flux_table[(flux_table['ref_ra'] > ra_max) + (flux_table['ref_ra'] < ra_min)].index)
        flux_table = flux_table.drop(flux_table[(flux_table['ref_dec'] > dec_max) + (flux_table['ref_dec'] < dec_min)].index)
    
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
        x_axis_label='RA (hour)', y_axis_label='DEC (deg)',
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

def get_mean_image_plot(source, result_dir, ra_ref=None, dec_ref=None, degrees_around_ref_coords=None):
    hdu = load_mean_image(f"{result_dir}/mean_image_reprojected.fits")
    if ra_ref:
        data, ra_dec = mean_image_data(hdu, ra_ref, dec_ref, degrees_around_ref_coords)
    else:
        data, ra_dec = mean_image_data(hdu)

    imdata = get_imdata(data,ra_dec)

    p = figure(
        tools=TOOLS,
        title='Mean Image',
        tooltips=[
            ("value", "@image Jy/beam"),
            ("RA", "$x{0.00}°"),
            ("DEC", "$y{0.00}°")
        ],
        x_axis_label='RA (hour)', y_axis_label='DEC (deg)',
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


def get_epoch_image_plots(epoch_files, mean_source, ra_ref=None, dec_ref=None, degrees_around_ref_coords=None):
    # Set up the figure
    p = figure(
        tools=TOOLS,
        title="Epoch_0",
        tooltips=[
            ("value", "@image Jy/beam"),
            ("RA", "$x{0.00}°"),
            ("DEC", "$y{0.00}°")
        ],
        x_axis_label='RA (hour)', y_axis_label='DEC (deg)',
        x_range=(min(mean_source.data['ref_ra']),  max(mean_source.data['ref_ra'])),
        y_range=(min(mean_source.data['ref_dec']), max(mean_source.data['ref_dec'])),
    )

    # Adjust hover tool to only work on image data
    hover_tool = p.select(type=HoverTool)
    hover_tool.names = ["epoch_image"]

    # Make a data dictionary of each epoch with the _(int) format keys
    data_dict = {}
    for ei, epoch_file in enumerate(epoch_files):
        hdu = load_mean_image(epoch_file)
        if ra_ref:
            data, ra_dec = mean_image_data(hdu, ra_ref, dec_ref, degrees_around_ref_coords)
        else:
            data, ra_dec = mean_image_data(hdu)
        imdata = get_imdata(data, ra_dec, ei=ei)
        data_dict.update(imdata)

    # Point to first image first
    data_dict.update({
            'image':data_dict["image_0"],
            # 'ra':data_dict["ra_0"],
            # 'dec':data_dict["dec_0"],
            'x':data_dict["x_0"],
            'y':data_dict["y_0"],
            'dw':data_dict["dw_0"],
            'dh':data_dict["dh_0"]})
    source = ColumnDataSource(data=data_dict)

    # must give a vector of image data for image parameter
    p.image(source=source,
            palette=palettes.mpl['Cividis'][256], name="epoch_image")
    slider = Slider(start=0, end=len(epoch_files)-1, value=0, step=1, title="Epoch")
    # Add the source circles
    p.circle(source=mean_source,
             x='ref_ra', y='ref_dec',
             radius=0.03, fill_color=None,
             line_width=1.5, line_color='yellow')

    callback = CustomJS(args=dict(source=source, slider=slider, p=p),
                        code="""
        const data = source.data;
        const i = slider.value;
        p.title.text = 'Epoch_'+i.toString(10);
        const image = data['image_'+i.toString(10)];
        data['image'] = image;
        source.change.emit();
    """)
    slider.js_on_change("value", callback)
    return p, slider
# def get_epoch_image_plots(source, epoch_source, num_epochs):
#     # Set up the figure
#     epochs = figure(
#         tools=TOOLS,
#         title="Epoch_0",
#         tooltips=[
#             ("value", "@image Jy/beam"),
#             ("RA", "$x{0.00}°"),
#             ("DEC", "$y{0.00}°")
#         ],
#         x_axis_label='RA',
#         y_axis_label='DEC',
#         x_range=(min(source.data['ref_ra']),  max(source.data['ref_ra'])),
#         y_range=(min(source.data['ref_dec']), max(source.data['ref_dec'])),
#     )



#     # Adjust hover tool to only work on image data
#     hover_tool = epochs.select(type=HoverTool)
#     hover_tool.names = ["epoch_image"]
#     epoch_slider = Slider(start=0, end=num_epochs-1, value=0, step=1, title="Epoch")

#     # must give a vector of image data for image parameter
#     epochs.image(source=epoch_source,
#             palette=palettes.mpl['Cividis'][256], name="epoch_image")


#     # Add the source circles
#     epochs.circle(source=source,
#                 x='ref_ra', y='ref_dec',
#                 radius=0.03, fill_color=None,
#                 line_width=1.5, line_color='yellow')

#     return epochs, epoch_slider

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
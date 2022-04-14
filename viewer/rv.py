#! /usr/bin/env python3

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import glob
import pandas as pd
from astropy.io.votable import parse
from functools import lru_cache
import datetime

from bokeh.plotting import figure, show, output_file
from bokeh import palettes
from bokeh.layouts import gridplot, layout
from bokeh.io import curdoc
from bokeh.models import (ColumnDataSource, Circle, Whisker,
                          DataTable, TableColumn, NumberFormatter,
                          CustomJS, Slider, DatetimeTickFormatter)


TOOLS = "box_select,lasso_select,help,pan,tap,wheel_zoom,reset"
selected_circle = Circle(fill_alpha=1, fill_color="firebrick", line_color=None)
nonselected_circle = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)

@lru_cache
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
def get_joined_table_source():
    # most likely there is a solution that doesn't involve so many intermediaries!
    flux_table = get_tabdata("results/flux_table.vot")
    stats_table = get_tabdata("results/stats_table.vot")
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
    print(f"epochs:{epochs}")
    print(f"dates:{dates}")
    lc.update({
        'epoch':epochs,
        'date':dates,
        'datetimes':datetimes,
        'current':fluxes[df['uuid'][0]],
        'current_upper':fluxes[df['uuid'][0]].values+err_fluxes[df['uuid'][0]].values,
        'current_lower':fluxes[df['uuid'][0]].values-err_fluxes[df['uuid'][0]].values
    })

    lc_source = ColumnDataSource(data=lc)
    return source, lc_source


def get_scatter_plots(source):

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


def get_mean_image_plot(source):
    hdu = load_mean_image('results/mean_image.fits')
    data, ra, dec = mean_image_data(hdu)
    imdata = get_imdata(data,ra,dec)

    p = figure(tools=TOOLS,
               title='Mean Image',
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

def get_light_curve_plot(source):
    tooltips = [("epoch","@date"),]
    lc_plot = figure(tools=TOOLS, width=300, height=300, title='Light curve:',
                    tooltips=tooltips,
                    x_axis_type="datetime",
                    x_axis_label='Epoch',
                    y_axis_label='Peak flux (Jy/beam)',
                    y_range=(0,1))
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

def get_epoch_image_plots(epoch_files, mean_source):
    # Set up the figure
    p = figure(
        tools=TOOLS,
        title="Epoch_0",
        tooltips=[
            ("value", "@image Jy/beam"),
            ("RA", "@ra{0.00}°"),
            ("DEC", "@dec{0.00}°")
        ],
        x_axis_label='RA',
        y_axis_label='DEC')
    p.x_range.range_padding = p.y_range.range_padding = 0

    # Make a data dictionary of each epoch with the _(int) format keys
    data_dict = {}
    for ei, epoch_file in enumerate(epoch_files):
        hdu = load_mean_image(epoch_file)
        data, ra, dec = mean_image_data(hdu)
        imdata = get_imdata(data, ra, dec, ei=ei)
        data_dict.update(imdata)
    # Point to first image first
    data_dict.update({
            'image':data_dict["image_0"],
            'ra':data_dict["ra_0"],
            'dec':data_dict["dec_0"],
            'x':data_dict["x_0"],
            'y':data_dict["y_0"],
            'dw':data_dict["dw_0"],
            'dh':data_dict["dh_0"]})
    source = ColumnDataSource(data=data_dict)

    # must give a vector of image data for image parameter
    p.image(source=source,
            palette=palettes.mpl['Cividis'][256])
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


def main():
    output_file(filename='viewer/RV.html', title="Robbie Viewer")
    #curdoc().theme = 'night_sky'
    source, lc_source = get_joined_table_source()
    mean_image = get_mean_image_plot(source)
    sky, variable, table = get_scatter_plots(source)
    lc = get_light_curve_plot(lc_source)
    epochs, epoch_slider = get_epoch_image_plots(glob.glob("results/E*.fits"), source)

    # Make the sky plot and mean image zoom/pan together
    sky.x_range = mean_image.x_range = epochs.x_range
    sky.y_range = mean_image.y_range = epochs.y_range

    # When a tranisent is selected
    source.selected.js_on_change('indices',
        CustomJS(args=dict(lc_source=lc_source, lc=lc, source=source),
            code="""
            if (cb_obj.indices.length == 0)
                return;
            var idx = cb_obj.indices[0];
            var uuid = source.data['uuid'][idx];
            var data = lc_source.data;
            var plot = lc;

            var lower = [];
            var upper = [];
            var min = 999;
            var max = -999;
            for (let i = 0; i<data[uuid].length; i++){
                lower.push(data[uuid][i]-data['err_'+uuid][i]);
                upper.push(data[uuid][i]+data['err_'+uuid][i]);
                if (min>lower[i]){min=lower[i]};
                if (max<upper[i]){max=upper[i]};
            }

            data['current'] = data[uuid];
            data['current_lower'] = lower;
            data['current_upper'] = upper;

            // update some of the plot visuals
            plot.title.text = ''+uuid;
            plot.y_range.start = min- Math.abs(min)*0.8;
            plot.y_range.end = max/0.8;

            // emit changes so that we can update all the other
            lc_source.change.emit();
            source.change.emit();
            """
        )
    )
    # When an epoch is selected
    lc_source.selected.js_on_change('indices',
        CustomJS(args=dict(lc_source=lc_source, epoch_slider=epoch_slider, source=source),
            code="""
            if (cb_obj.indices.length == 0)
                return;
            var idx = cb_obj.indices[0];
            console.log(idx);
            var slider = epoch_slider;

            slider.value = Number(idx);

            // emit changes so that we can update all the other
            epoch_slider.change.emit();
            """
        )
    )
    p = layout([ [sky, variable,mean_image,[epochs, epoch_slider], lc], [table]],
                 sizing_mode='scale_width')
    show(p)


if __name__ == "__main__":
    main()
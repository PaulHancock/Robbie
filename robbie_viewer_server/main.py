#! /usr/bin/env python

from astropy.io import fits
from astropy.io.votable import parse
from astropy.wcs import WCS
from astropy.coordinates import Angle
import astropy.units as u
import numpy as np
import glob
import pandas as pd
import datetime
import argparse
import sys
import os
import numpy as np

from bokeh.layouts import layout
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, Circle, CustomJS
from bokeh.models.formatters import FuncTickFormatter
from utils import get_joined_table_source, get_scatter_plots, get_mean_image_plot,\
    get_light_curve_plot, get_light_curve_plot, get_epoch_image_plots, load_mean_image, mean_image_data, get_imdata


# Initialise tools and Glyphs
TOOLS = "box_select,lasso_select,help,pan,tap,wheel_zoom,reset"
selected_circle = Circle(fill_alpha=1, fill_color="firebrick", line_color=None)
nonselected_circle = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)
    
# Get results from the command line args



if len(sys.argv) > 2:
    result_dir = sys.argv[1]
    ra_loc_centre = float(sys.argv[2])
    dec_loc_centre = float(sys.argv[3])
    degrees_around_centre = float(sys.argv[-1])
else:
    result_dir = sys.argv[-1]

# Find input fits files
reproject_img_dir = f"{result_dir}/reprojected_images"
all_fits = glob.glob(f"{reproject_img_dir}/*.fits")
epoch_fits = []
for fits in all_fits:
    if not ("mean" in fits or "transients" in fits):
        epoch_fits.append(fits)

# Get data
if len(sys.argv) > 2:
    source, lc_source = get_joined_table_source(result_dir, ra_loc_centre, dec_loc_centre, degrees_around_centre)
    mean_image = get_mean_image_plot(source, reproject_img_dir, ra_loc_centre, dec_loc_centre, degrees_around_centre)
else:
    source, lc_source = get_joined_table_source(result_dir)
    mean_image = get_mean_image_plot(source, reproject_img_dir)

sky, variable, table = get_scatter_plots(source)
lc = get_light_curve_plot(lc_source)

# Initialise empty CDS for Epochs
epochs_cds = ColumnDataSource()

# Get Epoch and slider
#epochs, epoch_slider = get_epoch_image_plots(source, epochs_cds, len(epoch_fits)) # New
if len(sys.argv) > 2:
    epochs, epoch_slider = get_epoch_image_plots(epoch_fits, source, ra_loc_centre, dec_loc_centre, degrees_around_centre)
else:
    epochs, epoch_slider = get_epoch_image_plots(epoch_fits, source)

# Make the sky plot and mean image zoom/pan together
sky.x_range = mean_image.x_range = epochs.x_range
sky.y_range = mean_image.y_range = epochs.y_range

# make the plot coords be in hms/dms format
dms = """
var s = Math.sign(tick);
var t = Math.abs(tick);
var deg = Math.floor(t);
t = (t -deg)*60;
var min = Math.floor(t);
t = (t-min)*60;
var sec = Math.floor(t);
var tk = ""
if (s < 0) {
  tk = String(s)[0]
}
tk = tk + String(deg).padStart(2,'0') + ":"+String(min).padStart(2,'0')+":"+String(sec).padStart(2,'0');
return tk;
"""
hms = """
var s = Math.sign(tick);
var t = Math.abs(tick/15);
var deg = Math.floor(t);
t = (t -deg)*60;
var min = Math.floor(t);
t = (t-min)*60;
var sec = Math.floor(t);
var tk = ""
if (s < 0) {
  tk = String(s)[0]
}
tk = tk + String(deg).padStart(2,'0') + ":"+String(min).padStart(2,'0')+":"+String(sec).padStart(2,'0');
return tk;
"""

sky.yaxis.formatter = mean_image.yaxis.formatter = epochs.yaxis.formatter = FuncTickFormatter(code=dms)
sky.xaxis.formatter = mean_image.xaxis.formatter = epochs.xaxis.formatter = FuncTickFormatter(code=hms)


# # Callback to update epochs
# def update_epochs():

#     # Make a data dictionary of each epoch with the _(int) format keys
#     data_dict = {}
#     for ei, epoch_file in enumerate([epoch_fits[epoch_slider.value]]):
#         print(ei, epoch_file)
#         hdu = load_mean_image(epoch_file)
#         if 'RA' in os.environ:
#             data, ra_dec = mean_image_data(hdu, ra_loc_centre, dec_loc_centre, degrees_around_centre)
#         else:
#             data, ra_dec = mean_image_data(hdu)

#         imdata = get_imdata(data, ra_dec, ei=ei)

#         data_dict.update(imdata)

#     # Point to first image first
#     data_dict.update({
#             'image':data_dict["image_0"],
#             # 'ra':data_dict["ra_0"],
#             # 'dec':data_dict["dec_0"],
#             'x':data_dict["x_0"],
#             'y':data_dict["y_0"],
#             'dw':data_dict["dw_0"],
#             'dh':data_dict["dh_0"]
#             }
#             )

#     epochs_cds.data = data_dict
  
  
# # Add on change to epoch slider
# epoch_slider.on_change('value', lambda attr, old, new: update_epochs())

# Epoch slider callback
callback = CustomJS(args=dict(source=source, slider=epoch_slider, p=epochs),
                    code="""
    const i = slider.value;
    p.title.text = 'Epoch_'+i.toString(10);
 
""")
epoch_slider.js_on_change("value", callback)

# When a tranisent is selected
source.selected.js_on_change('indices',
    CustomJS(args=dict(lc_source=lc_source, lc=lc, source=source, mean_image=mean_image),
        code="""
        if (cb_obj.indices.length == 0)
            return;
        var idx = cb_obj.indices[0];
        var uuid = source.data['uuid'][idx];
        var data = lc_source.data;
        var plot = lc;
        var image_plot = mean_image;

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
# When a tranisent is selected
source.selected.js_on_change('indices',
    CustomJS(args=dict(source=source, xr=mean_image.x_range, yr=mean_image.y_range, mean_image=mean_image),
        code="""
        if (cb_obj.indices.length == 0)
            return;
        var idx = cb_obj.indices[0];

        // zoom into the selected tranisent
        xr.start = source.data['ref_ra'][idx] - 0.2;
        xr.end = source.data['ref_ra'][idx] + 0.2;
        yr.start = source.data['ref_dec'][idx] - 0.2;
        yr.end = source.data['ref_dec'][idx] + 0.2;
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

# Initial data load
#update_epochs() 

# Add plots to the current document root
p = layout([ [sky, variable,mean_image,[epochs, epoch_slider], lc], [table]],
                 sizing_mode='scale_width')
curdoc().add_root(p)
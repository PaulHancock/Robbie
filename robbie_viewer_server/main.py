#! /usr/bin/env python3

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import glob
import pandas as pd
import datetime
import argparse
import sys

from bokeh.layouts import layout
from bokeh.io import curdoc
from bokeh.models import CustomJS
from utils import get_joined_table_source, get_scatter_plots, get_mean_image_plot,\
    get_light_curve_plot, get_light_curve_plot, get_epoch_image_plots


def main(result_dir):
 
    # Find input fits files
    reproject_img_dir = f"{result_dir}/reprojected_images"
    all_fits = glob.glob(f"{reproject_img_dir}/*.fits")
    epoch_fits = []
    for fits in all_fits:
        if not ("mean" in fits or "transients" in fits):
            epoch_fits.append(fits)

    source, lc_source = get_joined_table_source(result_dir)
    mean_image = get_mean_image_plot(source, reproject_img_dir)
    sky, variable, table = get_scatter_plots(source)
    lc = get_light_curve_plot(lc_source)
    epochs, epoch_slider = get_epoch_image_plots(epoch_fits, source)

    # Make the sky plot and mean image zoom/pan together
    sky.x_range = mean_image.x_range = epochs.x_range
    sky.y_range = mean_image.y_range = epochs.y_range

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
    p = layout([ [sky, variable,mean_image,[epochs, epoch_slider], lc], [table]],
                 sizing_mode='scale_width')

    return p


# Get results from the command line args
result_dir = sys.argv[-1]

# Get the outputted plots
p = main(result_dir)

# Add plots to the current document root
curdoc().add_root(p)
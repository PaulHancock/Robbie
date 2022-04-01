#! /usr/bin/env python3

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh import palettes


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

def main():
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

    output_file(filename='viewer/RV.html', title="Robbie Viewer")
    show(p)

if __name__ == "__main__":
    main()
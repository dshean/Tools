#! /usr/bin/env python
"""
Precompute valid tiles for dem_mosaic
"""

import argparse
import math
from collections import OrderedDict

from osgeo import gdal, ogr

from pygeotools.lib import geolib
from pygeotools.lib import warplib

def getparser():
    parser = argparse.ArgumentParser(description='Utility to compute valid dem_mosaic tiles for input res/extent/proj')
    parser.add_argument('-tr', default='min', help='Output resolution (default: %(default)s)')
    parser.add_argument('-te', default='union', help='Output extent (default: %(default)s)')
    parser.add_argument('-t_srs', default='first', help='Output projection (default: %(default)s)')
    parser.add_argument('-tile_width', default=100000., type=float, help='Output tile width (meters)')
    parser.add_argument('src_fn_list', nargs='+', help='Input filenames (img1.tif img2.tif ...)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    #Input filelist
    fn_list = args.src_fn_list
    #Might hit OS open file limit here
    print("Loading all datasets")
    ds_list = [gdal.Open(fn) for fn in fn_list]

    #Mosaic t_srs
    print("Parsing t_srs")
    t_srs = warplib.parse_t_srs(args.t_srs, ds_list)
    print(t_srs)

    #Mosaic res
    print("Parsing tr")
    tr = warplib.parse_res(args.tr, ds_list, t_srs=t_srs) 
    print(tr)

    #Mosaic extent 
    #xmin, ymin, xmax, ymax
    print("Parsing te")
    te = warplib.parse_extent(args.te, ds_list, t_srs=t_srs) 
    print(te)
    mos_xmin, mos_ymin, mos_xmax, mos_ymax = te

    #Tile dimensions in output projected units (meters)
    tile_width = args.tile_width
    tile_height = tile_width

    #Compute extent geom for all input datsets
    print("Computing extent geom for all input datasets")
    input_geom_dict = OrderedDict()
    for ds in ds_list:
        ds_geom = geolib.ds_geom(ds, t_srs)
        #Could use filename as key here
        input_geom_dict[ds] = ds_geom

    #Mosaic tile size
    #Should have float extent and tile dim here
    ntiles_w = int(math.ceil((mos_xmax - mos_xmin)/tile_width))
    ntiles_h = int(math.ceil((mos_ymax - mos_ymin)/tile_height))
    ntiles = ntiles_w * ntiles_h
    print("%i (%i cols x %i rows) tiles required for full mosaic" % (ntiles, ntiles_w, ntiles_h))
    #Use this for zero-padding of tile number
    ntiles_digits = len(str(ntiles))

    print("Computing extent geom for all output tiles")
    tile_dict = OrderedDict()
    for i in range(ntiles_w):
        for j in range(ntiles_h):
            tilenum = j*ntiles_w + i
            tile_xmin = mos_xmin + i*tile_width
            tile_xmax = mos_xmin + (i+1)*tile_width
            tile_ymax = mos_ymax - j*tile_height
            tile_ymin = mos_ymax - (j+1)*tile_height
            #Corner coord needed for geom
            x = [tile_xmin, tile_xmax, tile_xmax, tile_xmin, tile_xmin]
            y = [tile_ymax, tile_ymax, tile_ymin, tile_ymin, tile_ymax]
            tile_geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(x,y)]))
            tile_geom = ogr.CreateGeometryFromWkt(tile_geom_wkt)
            tile_geom.AssignSpatialReference(t_srs)
            tile_dict[tilenum] = tile_geom

    out_tile_list = []
    print("Computing valid intersections between input dataset geom and tile geom")
    for tilenum, tile_geom in tile_dict.iteritems():
        for ds, ds_geom in input_geom_dict.iteritems():
            if tile_geom.Intersects(ds_geom):
                out_tile_list.append(tilenum)
                #Write out shp for debugging
                #geolib.geom2shp(tile_geom, 'tile_%03i.shp' % tilenum)
                break 

    print("Valid output tiles")
    out_tile_list.sort()
    out_tile_list_str = ' '.join(map(str, out_tile_list))
    print(out_tile_list_str)

    out_fn = 'mos_tilenum_list.txt'
    with open(out_fn, 'w') as f:
        f.write(out_tile_list_str)
        
if __name__ == "__main__":
    main()

#! /usr/bin/env python
"""
Precompute valid tiles for dem_mosaic
"""
#res=32
#res=8
#mkdir conus_${res}m_tile
#lfs setstripe conus_${res}m_tile --count 64
#~/src/Tools/dem_mosaic_validtiles.py --tr $res --t_projwin 'union' --t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o conus_${res}m_tile/conus_${res}m *00/*/*DEM_${res}m.tif

#res=8
#mkdir hma_${res}m_tile
#lfs setstripe hma_${res}m_tile --count 64
#~/src/Tools/dem_mosaic_validtiles.py --tr $res --t_projwin 'union' --t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o hma_${res}m_tile/hma_${res}m */*00/*/*DEM_${res}m.tif

import os
import argparse
import math
from collections import OrderedDict
import time
import subprocess

from osgeo import gdal, ogr

from pygeotools.lib import geolib, warplib, iolib

def get_dem_mosaic_cmd(fn_list, o, tr, t_srs, t_projwin, georef_tile_size, threads, tile):
    """
    Call dem_mosaic command with valid tile number 
    Eventually, could call with multiple tiles
    """
    cmd = ['dem_mosaic', '--threads', threads, '--tr', tr, '--t_srs', t_srs.ExportToProj4(), \
           '--georef-tile-size', georef_tile_size, '-o', o, '--t_projwin']
    cmd.extend(t_projwin)
    #Not yet implemented
    #cmd.extend(tile_list)
    cmd.append('--tile-index')
    cmd.append(tile)
    cmd.extend(fn_list)
    cmd = [str(i) for i in cmd]
    #print(cmd)
    #return subprocess.call(cmd)
    return cmd

def getparser():
    parser = argparse.ArgumentParser(description='Wrapper for dem_mosaic that will only write valid tiles')
    parser.add_argument('--tr', default='min', help='Output resolution (default: %(default)s)')
    parser.add_argument('--t_projwin', default='union', help='Output extent (default: %(default)s)')
    parser.add_argument('--t_srs', default='first', help='Output projection (default: %(default)s)')
    parser.add_argument('--georef_tile_size', type=float, default=100000., help='Output tile width (meters)')
    parser.add_argument('--threads', type=int, default=iolib.cpu_count(), help='Number of threads')
    parser.add_argument('-o', type=str, default=None, help='Output mosaic prefix')
    parser.add_argument('src_fn_list', type=str, nargs='+', help='Input filenames (img1.tif img2.tif ...)')
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
    t_srs = warplib.parse_srs(args.t_srs, ds_list)
    print(t_srs)

    #Mosaic res
    print("Parsing tr")
    tr = warplib.parse_res(args.tr, ds_list, t_srs=t_srs) 
    print(tr)

    #Mosaic extent 
    #xmin, ymin, xmax, ymax
    print("Parsing t_projwin")
    t_projwin = warplib.parse_extent(args.t_projwin, ds_list, t_srs=t_srs) 
    print(t_projwin)
    #This could trim off some fraction of a pixel around margins
    t_projwin = geolib.extent_round(t_projwin, tr)
    mos_xmin, mos_ymin, mos_xmax, mos_ymax = t_projwin

    #Tile dimensions in output projected units (meters)
    #Assume square
    tile_width = args.georef_tile_size
    tile_height = tile_width

    o = args.o
    if o is None:
        o = 'mos_%im/mos' % tr
    threads = args.threads

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
    
    #Could also preserve list of input files that intersect tile
    #Then only process those files for given tile bounds 
    #Avoid loading all files for each dem_mosaic call

    print("%i valid output tiles" % len(out_tile_list))
    out_tile_list.sort()
    out_tile_list_str = ' '.join(map(str, out_tile_list))
    print(out_tile_list_str)

    out_fn = o+'_tilenum_list.txt'
    with open(out_fn, 'w') as f:
        f.write(out_tile_list_str)

    print("Running dem_mosaic in parallel")
    dem_mosaic_args = (fn_list, o, tr, t_srs, t_projwin, tile_width, 1)
    processes = []
    log = False
    delay = 1
    outf = open(os.devnull, 'w') 
    for tile in out_tile_list:
        cmd = get_dem_mosaic_cmd(*dem_mosaic_args, tile=tile)
        if log:
            outf = open('%s-log-dem_mosaic-tile-%i.log' % (o, tile), 'w')
            #can add close_fds=True to Popen to close logfile
        processes.append(subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT))
        time.sleep(delay)
    for p in processes: p.wait()
    outf = None

if __name__ == "__main__":
    main()

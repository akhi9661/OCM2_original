import rasterio, os, re, shutil
import numpy as np
from osgeo import gdal, osr, gdalconst

def ExportSubdatasets(path, hdf_file):
    
    '''
    This function takes the folder path and the HDF file as input and exports individual layers to TIFF (named GeoTIFF)
    
    '''
    opf_tif = os.path.join(path, 'GeoTiff')
    if os.path.exists(opf_tif):
        shutil.rmtree(opf_tif)
    os.makedirs(opf_tif)
    
    inp_hdf = os.path.join(path, hdf_file)
    hdf_ds = gdal.Open(inp_hdf, gdal.GA_ReadOnly)
    subdatasets = hdf_ds.GetSubDatasets()
    
    for i in range(0, len(subdatasets)):
        subdataset_name = subdatasets[i][0]
        band_ds = gdal.Open(subdataset_name, gdal.GA_ReadOnly)
        band_path = os.path.join(opf_tif, 'band{0}.TIF'.format(i))
        if band_ds.RasterCount > 1:
            for j in range(1,band_ds.RasterCount + 1):
                band = band_ds.GetRasterBand(j)
                band_array = band.ReadAsArray()
        else:
            band_array = band_ds.ReadAsArray()
        
        out_ds = gdal.GetDriverByName('GTiff').Create(band_path,
                                                      band_ds.RasterXSize,
                                                      band_ds.RasterYSize,
                                                      1,
                                                      gdal.GDT_Float64,
                                                      ['COMPRESS=LZW', 'TILED=YES'])
        
        
        out_ds.SetGeoTransform(band_ds.GetGeoTransform())
        out_ds.SetProjection(band_ds.GetProjection())
        out_ds.GetRasterBand(1).WriteArray(band_array)
        out_ds.GetRasterBand(1).SetNoDataValue(-32768)
        
    out_ds = None
        
    return opf_tif
    
def metaInfo(path, hdf_file):
    
    inp = gdal.Open(os.path.join(path, hdf_file))
    meta = inp.GetMetadata()
    
    ulx, uly = float(meta['Upper Left Longitude']), float(meta['Upper Left Latitude'])
    urx, ury = float(meta['Upper Right Longitude']), float(meta['Upper Right Latitude'])
    blx, bly = float(meta['Lower Left Longitude']), float(meta['Lower Left Latitude'])
    brx, bry = float(meta['Lower Right Longitude']), float(meta['Lower Right Latitude'])
    
    sun_elev = float(meta['Sun Elevation Angle'])
    
    return (ulx,  uly), (urx, ury), (brx, bry), (blx, bly), (sun_elev)
    
def GetExtent(ds):
    
    ''' Return list of corner coordinates from a gdal Dataset '''
    
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)


def Georeference(inpf, gtif, meta, opf_ref):
    
    inp_file = os.path.join(inpf, gtif)
    band_tif = gdal.Open(inp_file)
    ext = GetExtent(band_tif)
    
    ext_pos = [(abs(val[0]), abs(val[1])) for val in ext]

    out_file = os.path.join(opf_ref, os.path.basename(gtif).split('.')[0] + '_georef.TIF')
    shutil.copy(inp_file, out_file)
    ds = gdal.Open(out_file, gdal.GA_Update)
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326) 
    
    
    '''
    Enter the GCPs
    Format: [map x-coordinate(longitude)], [map y-coordinate (latitude)], [elevation],
    [image column index(x)], [image row index (y)]
    If map pixel is negative, multiply it with -1 to make positive since GDAL can't handle negative pixel position that well.
    
    '''
    
    gcps = [gdal.GCP(meta[0][0], meta[0][1], 0, ext_pos[0][0], ext_pos[0][1]), 
            gdal.GCP(meta[1][0], meta[1][1], 0, ext_pos[1][0], ext_pos[1][1]),
            gdal.GCP(meta[2][0], meta[2][1], 0, ext_pos[2][0], ext_pos[2][1]), 
            gdal.GCP(meta[3][0], meta[3][1], 0, ext_pos[3][0], ext_pos[3][1])]
    
    ds.SetGCPs(gcps, sr.ExportToWkt())
    ds = None
    
    return 'Done'


def calc_toa(rad, sun_elev, band_no):
    esol = [1.72815, 1.85211, 1.9721, 1.86697, 1.82781, 1.65765, 1.2897, 0.952073]
    toa_reflectance = (np.pi * 1 * rad * 10) / (esol[band_no] * 1000 * math.sin(math.radians(sun_elev)))
    return toa_reflectance

def toa_convert(inpf, inp_name, opf, sun_elev):
    
    band_no = int(''.join(list(filter(str.isdigit, inp_name.split('.')[0].split('_')[0]))))
    with rasterio.open(os.path.join(inpf, inp_name)) as (r):
        rad = r.read(1).astype('float32')
        profile = r.profile
    
    toa = calc_toa(rad, sun_elev, band_no)
    toa[toa < 0] = 0.0
    toa[toa > 2] = 0.0   
    op_name = os.path.basename(inp_name).split('.')[0] + '.TIF'
    with (rasterio.open)((os.path.join(opf, op_name)), 'w', **profile) as (dataset):
        dataset.write(toa, 1)
    dataset.close()
        
    return 'done'

def list_files(inpf, inp_name, files):
    
    files.append(os.path.join(inpf, inp_name))
    return files
        
    
def sum_toa(filelist):
    
    with rasterio.open(filelist[0]) as r:
        arr = r.read()
        profile = r.profile
        
    for f in filelist[1:]:
        with rasterio.open(f) as r:
            assert profile == r.profile, 'stopping, file {} and  {} do not have matching profiles'.format(filelist[0], f)
            arr = arr + r.read()

    return (arr)

def toa_other(filelist):
     
    one, two, seven = [filelist[check] for check in [0,1,6]]
    
    with rasterio.open(one) as r:
        band1 = r.read()
        profile = r.profile
    shape = band1.shape
    with rasterio.open(two) as r:
        band2 = r.read()
    with rasterio.open(seven) as r:
        band7 = r.read()
        
    
    band2[band2 == 0.0] = np.nan
    band7[band7 == 0.0] = np.nan

    toa_diff = band2 - band1
    toa_ratio = band2/band7 

    return (toa_diff, toa_ratio, shape, profile)

def cloudmask_ocm(inpf, filelist):

'''
Cloud mask based on Mishra et. al., 2018 (https://doi.org/10.1007/s12524-017-0715-5)
'''
    
    toa_sum = sum_toa(filelist)
    toa_diff, toa_ratio, shape, profile = toa_other(filelist) 
    
    cldmsk = np.zeros(shape, dtype = 'float32')
    cldmsk = np.where(((toa_sum > 2.7) & (toa_ratio > 1.5) & (toa_diff < 0)), 1, 0)
    
    with (rasterio.open)((os.path.join(inpf, 'cloud_mask.TIF')), 'w', **profile) as (dst):
        dst.write(cldmsk)
    dst.close()
    
    return cldmsk


meta = metaInfo(path, hdf_file)
opf_tif = ExportSubdatasets(path, hdf_file)
print('Done: Layers converted to GeoTIFF. Wait.')

opf_ref = os.path.join(path, 'Reflectance')
if os.path.exists(opf_ref):
    shutil.rmtree(opf_ref)
os.makedirs(opf_ref)

opf_georef = os.path.join(path, 'Georeferenced')
if os.path.exists(opf_georef):
    shutil.rmtree(opf_georef)
os.makedirs(opf_georef)

def do_ref(opf_tif, meta, opf_ref):
    
    original = os.listdir(opf_tif)
    gtif = list(filter(lambda x: x.endswith('TIF'), original))
    for band_name in gtif:
        if (int(''.join(list(filter(str.isdigit, band_name.split('.')[0].split('_')[0]))))) <= 7:
            toa_convert(opf_tif, band_name, opf_ref, meta[4])
        else:
            shutil.copy(os.path.join(opf_tif, band_name), os.path.join(opf_ref, band_name))
            
    return None
    
def do_georef(geo_ref, meta, opf_georef):
    
    original = os.listdir(opf_ref)
    gtif = list(filter(lambda x: x.endswith('TIF'), original))
    for band_name in gtif:
        Georeference(opf_ref, band_name, meta, opf_georef)
        
    return None
    
def do_cldmsk(opf_ref):
    
    files = []

    original = os.listdir(opf_ref)
    gtif = list(filter(lambda x: x.endswith('TIF'), original))
    for band_name in gtif:
        if (int(''.join(list(filter(str.isdigit, band_name.split('.')[0].split('_')[0]))))) <= 7:
            filelist = list_files(opf_ref, band_name, files)

    cldmsk = cloudmask_ocm(opf_ref, filelist)
    
    return None


do_ref(opf_tif, meta, opf_ref)
print('Done: Reflectance conversion. Wait.')
do_cldmsk(opf_ref)
print('Done: Cloudmasking. Wait.')
do_georef(opf_ref, meta, opf_georef)
print('Done: Georeferncing')

print('\nRemoving temporary folders. Wait.')
if os.path.exists(opf_tif):
    shutil.rmtree(opf_tif)
    
if os.path.exists(opf_ref):
    shutil.rmtree(opf_ref)
    
print('\nDone! Open Georeferenced folder for results')

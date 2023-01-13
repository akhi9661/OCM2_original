This python code module takes OCM2 HDF file as input and generates georeferenced reflectance, angle files, and cloud mask layer

Input folder: path
Input file: hdf_file

Outputs:
Final output folder: Georeferenced

Two temporary output folders are also generated which are automatically deleted after the completion.
GeoTiff: Contains geotiff files of each layer in HDF file
Reflectance: Contains reflectance (of seven bands) and cloudmask layer. 

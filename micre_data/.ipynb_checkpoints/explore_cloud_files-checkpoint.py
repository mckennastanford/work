import xarray
import glob

path = '/Volumes/STANFORD_6/cloud_precip_properties/'

files = glob.glob(path+'*.nc')
files = sorted(files)
len(files)

key_list = []

for file in files:
    ncfile = xarray.open_dataset(file)
    for key in ncfile.keys():
        key_list.append(key)
    print(aaaa)
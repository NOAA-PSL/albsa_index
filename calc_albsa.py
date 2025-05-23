#!/usr/bin/env python3
# -*- coding: utf-8 -*-  

# ################################################################################################################
#
# AUTHORS:
#
#   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
#
# ################################################################################################################
#
# PURPOSE:
#
# Create a NetCDF4 file containing daily calculations of the ALBSA index.
# ALBSA: "Aleutian Low - Beaufort Sea Anticyclone" (Cox et al., 2019)
# It is a 4-pt index on the geopotential height field a 850 hPa where the points are 
#   North (“N” at 75°N/170°W), 
#   South (“S” at 50°N/170°W), 
#   East (“E” at 55°N/150°W), and 
#   West (“W” at 55°N/160°E)
#
#   ALBSA = ( E-W ) - ( N - S ) and is expressed in meters of GPH
#
# The index has been most commonly defined based on the NCEP/NCAR Reanalysis 1, but in depreciation 
#   this routine calculates ALBSA using the ECMWF ERA5 reanalysis.
# Testing in May 2025 shows differences between the reanalyses of mean X m, std X m, averaged for all
#   days overlapping 1948-2025.
# This code handles the download from ECMWF, the calculation, and the writes the new NetCDF.
#
# Cox, C. J., R. S. Stone, D. C. Douglas, D. M. Stanitski, and D. C. Douglas (2019), The Aleutian Low - 
#   Beaufort Sea Anticyclone: A climate index correlated with the timing of springtime melt in the
#   Pacific Arctic cryosphere. Geophysical Research Letters, 46(13), 7464-7473, 
#   https://doi.org/10.1029/2019GL083306
#
# ################################################################################################################
#
# DATA and ACKNOWLEDGEMENT:
#
# The product is generated using Copernicus Climate Change Service information, 2025. 
# ECMWF Reanalysis v5 (ERA5) (Hersbach, et al. 2020)
#
# Hersbach, H., B. Bell, P. Berrisford, S. Hirahara, A. Horányi, J. Muñoz-Sabater, J. Nicolas, C. Peubey, 
#   R.Radu, D. Schepers, A. Simmons, C. Soci, S. Abdalla, X. Abellan, G. Balsamo, P. Bechtold, G. Biavati, 
#   J. Bidlot, M. Bonavita, G. De Chiara, P. Dahlgren, D. Dee, M. Diamantakis, R. Dragani, J. Flemming, 
#   R. Forbes, M. Fuentes, A. Geer, L. Haimberger, S. Healy, R. J. Hogan, E. Hólm, M. Janisková, S. Keeley, 
#   P. Laloyaux, P. Lopez, C. Lupu, G. Radnoti, P. de Rosnay, I. Rozum, F. Vamborg, S. Villaume, and J.-N. Thépaut 
#   (2020) The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), 1999-2049,
#   https://doi.org/10.1002/qj.3803
#
# ################################################################################################################
#
# HOW TO:
#
#   python3 ./calc_albsa.py -p /path/to/working/directory/ -s 1940 -e 2025
#
#   -p: path to working directory (REQUIRED)
#   -s: starting year for the desired file contents (OPTIONAL; default to 1940)
#   -e: ending year for the desired file contents (OPTIONAL; default to current year)
#
# DEPENDENCIES:
#
#   python  ≥ 3.13.2
#   netCDF4 ≥ 1.7.2
#   xarray ≥ 2024.11.0
#   cdsapi ≥ 0.7.6
#
# ################################################################################################################

# import modules
import os, argparse, cdsapi, datetime, xarray as xr

# you are required to provide a directory path
# you are permitted to specify beginning and end dates (year granularity)
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', metavar='str', help='path to working directory, include trailing /')
parser.add_argument('-s', '--start_date', metavar='str', help='beginning of processing period, yyyy')
parser.add_argument('-e', '--end_date', metavar='str', help='ending of processing period, yyyy')
args = parser.parse_args()

try: working_dir = args.path
except: print("Please specify a working directory path using -e")
if working_dir[-1] != "/": working_dir = working_dir+"/"

# create a list of years to access

if args.start_date: yrbeg = int(args.start_date)
else: yrbeg = 1940 # earliest date in the data set

if args.end_date: yrend = int(args.end_date)
else: yrend = datetime.datetime.now().year # to the current year

if yrbeg > yrend: raise ValueError("beginning year must >= ending year")

yrlist = []
for i in range(yrbeg,yrend+1):
    yrlist.append(str(i))

# uses the Climate Data Store (CDS) API to download data from ECMWF
def era5_downloader(yrlist):

    # create a request for ERA5 presure level 850 hPa, 4x daily, all years for the AK domain, 2.5d resolution
    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": [
            "reanalysis",
        ],
        "variable": ["geopotential"],
        "year": yrlist,
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "day": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"
        ],
        "time": [
            "00:00", "06:00", "12:00",
            "18:00"
        ],
        "pressure_level": ["850"],
        "grid": ["2.5","2.5"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [85, 150, 40, -120]
    }

    client = cdsapi.Client()
    client.retrieve(dataset, request, working_dir+'tmp.nc') # tmp.nc is a temporary file that will be deleted at the end of the code


def main():

    # download data
    era5_downloader(yrlist)

    # open the temporary ERA5 file using xarray
    file = xr.open_dataset(working_dir+'tmp.nc')

    # resample from 4x/day to daily
    file = file.resample(valid_time="D").mean() 

    file['longitude'] = file['longitude'] + 180

    # calculate albsa

    # coordinates.
    laS = 50
    loS = 360-170
    laN = 75
    loN = 360-170
    laE = 55
    loE = 360-150
    laW = 55
    loW = 160

    # 850 hPa daily mean GPH at 4 ALBSA coordinates.
    # The method is NN, but the data has been downloaded on pre-processed 2.5d grid, 
    #   which is a coarser, but divisible resampling of the native 0.25d grid.  
    Scoord = file.sel(latitude=laS, longitude=loS, method="nearest")/10
    Ncoord = file.sel(latitude=laN, longitude=loN, method="nearest")/10
    Ecoord = file.sel(latitude=laE, longitude=loE, method="nearest")/10
    Wcoord = file.sel(latitude=laW, longitude=loW, method="nearest")/10

    # calculate ALBSA index
    #   Eq. (3) from Cox et al. (2019)
    albi = ( Ecoord.z - Wcoord.z ) - ( Ncoord.z - Scoord.z )
    file = file.assign(index=albi.transpose())

    print(file)
    
    # remove some variables
    albi['number'].drop_attrs()
    albi = albi.drop_vars('number') 

    # write back to new file
    file.to_netcdf(working_dir+'albsa_index.nc')

    # delete the intermediary file
    #os.remove(working_dir+'tmp.nc')


# executes main():
if __name__ == "__main__":
    main()
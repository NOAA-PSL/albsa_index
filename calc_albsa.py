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
# Testing in May 2025 shows differences between the reanalyses of mean -0.24 m, std 19.5 m, averaged for all 
#   days overlapping 1948-2025 (n = 28259). The mean absolute value of the error variance is ~0.035%.
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
# The product is generated using Copernicus Climate Change Service information (CCCS/CDS, 2023). 
# ECMWF Reanalysis v5 (ERA5) (Hersbach, et al. 2020; Hersbach et al., 2023)
#
# Copernicus Climate Change Service, Climate Data Store, (2023): ERA5 hourly data on pressure levels from 1940 to present. 
#   Copernicus Climate Change Service (C3S) Climate Data Store (CDS), https://doi.org/10.24381/cds.bd0915c6
#
# Hersbach, H., B. Bell, P. Berrisford, G. Biavati, A. Horányi, J. Muñoz Sabater, J. Nicolas, C. Peubey, C., 
#   R. Radu, I. Rozum, D. Schepers, A. Simmons, C. Soci, D. Dee, and J.-N. Thépaut (2023): ERA5 hourly data on 
#   pressure levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), 
#   https://doi.org/10.24381/cds.bd0915c6
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
import os, shutil, argparse, cdsapi, datetime, glob, xarray as xr

# filenames
fname = "albsa_index.nc" # the new file name

# you are required to provide a directory path
# you are permitted to specify beginning and end dates (year granularity)
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", metavar="str", help="path to working directory, include trailing /")
parser.add_argument("-s", "--start_date", metavar="str", help="beginning of processing period, yyyy")
parser.add_argument("-e", "--end_date", metavar="str", help="ending of processing period, yyyy")
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

# for the NetCDF history
def code_version():
    cv = '1.0, 5/28/2025'
    return cv

# uses the Climate Data Store (CDS) API to download data from ECMWF
def era5_downloader(yrlist,tname):

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
        "grid": ["2.5/2.5"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [85, 150, 40, -120]
    }

    client = cdsapi.Client()
    client.retrieve(dataset, request, working_dir+tname) # the tmp file is a temporary file that will be deleted at the end of the code

# encode desired global attributes (in order of appearance) here
def define_glob_atts():

    # key-value pairs are attributes names / attributes
    global_atts = {
            "title"            :"Aleutian Low-Beaufort Sea Anticyclone (ALBSA) index",
            "institution"      :"National Oceanic and Atmospheric Administration (NOAA) Physical Sciences Laboratory (PSL)",
            "history"          :"created "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" with calc_albsa.py version "+code_version(),
            "file_creator"     :"Christopher J. Cox",
            "creator_email"    :"christopher.j.cox@noaa.gov", 
            "funding"          :"NOAA Physical Sciences Laboratory (PSL)",
            "source"           :"ECMWF Reanalysis v5 (ERA5)",
            "license"          :"Creative Commons Attribution 4.0 License, CC 4.0",  
            "methods"          :"Refer to Cox et al. (2019)",
            "conventions"      :"CF-1.12",  
            "commment"         :"",
            "keywords"         :"Arctic, climate index, snow, Alaska, atmospheric advection, sea ice",
            "acknowledgements" :"This product was generated using the Copernicus Climate Change Service information (CCCS/CDS, 2023) using data ECMWF Reanalysis v5 (ERA5) (Hersbach, et al. 2020; Hersbach et al., 2023).",
            "references"       :"""References:
        
                                    Copernicus Climate Change Service, Climate Data Store, (2023): ERA5 hourly data on pressure levels from 1940 to present. 
                                        Copernicus Climate Change Service (C3S) Climate Data Store (CDS), https://doi.org/10.24381/cds.bd0915c6

                                    Cox, C. J., R. S. Stone, D. C. Douglas, D. M. Stanitski, and D. C. Douglas (2019), The Aleutian Low - 
                                        Beaufort Sea Anticyclone: A climate index correlated with the timing of springtime melt in the
                                        Pacific Arctic cryosphere. Geophysical Research Letters, 46(13), 7464-7473, 
                                        https://doi.org/10.1029/2019GL083306

                                    Hersbach, H., B. Bell, P. Berrisford, G. Biavati, A. Horányi, J. Muñoz Sabater, J. Nicolas, C. Peubey, C., 
                                        R. Radu, I. Rozum, D. Schepers, A. Simmons, C. Soci, D. Dee, and J.-N. Thépaut (2023): ERA5 hourly data on 
                                        pressure levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), 
                                        https://doi.org/10.24381/cds.bd0915c6

                                    Hersbach, H., B. Bell, P. Berrisford, S. Hirahara, A. Horányi, J. Muñoz-Sabater, J. Nicolas, C. Peubey, 
                                        R.Radu, D. Schepers, A. Simmons, C. Soci, S. Abdalla, X. Abellan, G. Balsamo, P. Bechtold, G. Biavati, 
                                        J. Bidlot, M. Bonavita, G. De Chiara, P. Dahlgren, D. Dee, M. Diamantakis, R. Dragani, J. Flemming, 
                                        R. Forbes, M. Fuentes, A. Geer, L. Haimberger, S. Healy, R. J. Hogan, E. Hólm, M. Janisková, S. Keeley, 
                                        P. Laloyaux, P. Lopez, C. Lupu, G. Radnoti, P. de Rosnay, I. Rozum, F. Vamborg, S. Villaume, and J.-N. Thépaut 
                                        (2020) The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), 1999-2049,
                                        https://doi.org/10.1002/qj.3803"""
    }

    return global_atts

# encode additional variable attributes here
def define_var_atts():
       
    # keys are place holders      {"variable","att_name","att_value"}
    var_atts = {

        "index_att1"            :("index", "long_name", "ALBSA index"),
        "index_att2"            :("index", "units", "meters"),
        "index_att3"            :("index", "standard_name", ""),

        "scoord_att1"           :("s_coord", "long_name", "geopotential height at S coordinate"),
        "scoord_att2"           :("s_coord", "units", "meters"),
        "scoord_att3"           :("s_coord", "standard_name", "geopotential_height"),

        "ncoord_att1"           :("n_coord", "long_name", "geopotential height at N coordinate"),
        "ncoord_att2"           :("n_coord", "units", "meters"),
        "ncoord_att3"           :("n_coord", "standard_name", "geopotential_height"),

        "ecoord_att1"           :("e_coord", "long_name", "geopotential height at E coordinate"),
        "ecoord_att2"           :("e_coord", "units", "meters"),
        "ecoord_att3"           :("e_coord", "standard_name", "geopotential_height"),

        "wcoord_att1"           :("w_coord", "long_name", "geopotential height at W coordinate"),
        "wcoord_att2"           :("w_coord", "units", "meters"),
        "wcoord_att3"           :("w_coord", "standard_name", "geopotential_height"),
    }

    return var_atts

def main():

    # download data
    #   this needs to be batched. limits for CDS are 60K "items", which is 1 field x 1 var x 1 level x n time steps
    #   for this application at 4x/day, 20 years is about 30K. So lets do 20 year batches.
    batch_size = 5
    for i in range(0, len(yrlist), batch_size):
        era5_downloader(yrlist[i:i+batch_size],"tmp_"+yrlist[i]+".nc")

    # open the temporary ERA5 files using xarray and merge them into one xr.Dataset()
    file = xr.Dataset()
    for filename in glob.glob(working_dir+"tmp*.nc"):
        filetmp = xr.open_dataset(filename)
        file = xr.merge([file,filetmp])

    # resample from 4x/day to daily
    file = file.resample(valid_time="D").mean() 

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

    # convert the geopotential field [m^2/s^2] to geopotential height [m]
    z_attrs = file['z'].attrs.copy() # save the attributes
    grav = 9.80665 # gravity constant used by ECMWF
    file['z'] = file['z']/grav
    file['z'].attrs = z_attrs # reassign the attributes

    # 850 hPa daily mean GPH at 4 ALBSA coordinates.
    # The method is NN, but the data has been downloaded on pre-processed 2.5d grid, 
    #   which is a coarser, but divisible resampling of the native 0.25d grid.  
    Scoord = file.sel(latitude=laS, longitude=loS, method="nearest")
    Ncoord = file.sel(latitude=laN, longitude=loN, method="nearest")
    Ecoord = file.sel(latitude=laE, longitude=loE, method="nearest")
    Wcoord = file.sel(latitude=laW, longitude=loW, method="nearest")

    # calculate ALBSA index
    #   Eq. (3) from Cox et al. (2019)
    albi = ( Ecoord.z - Wcoord.z ) - ( Ncoord.z - Scoord.z )
    
    # create vars for index variables
    file = file.assign(index=albi.transpose())
    file = file.assign(s_coord=Scoord.z.transpose())
    file = file.assign(n_coord=Ncoord.z.transpose())
    file = file.assign(e_coord=Ecoord.z.transpose())
    file = file.assign(w_coord=Wcoord.z.transpose())

    # housekeeping
        # placeholder from CDS for ensembles. deceptive for reanalysis, in this case = 1
    file["number"].drop_attrs() 
    file = file.drop_vars("number") 
        # 4d to 3d bv squeezing the single-level pressure field
    file["z"] = file["z"].squeeze()
        # use single 
    file["latitude"] = file["latitude"].astype("float32")
    file["longitude"] = file["longitude"].astype("float32")

    # if an existing data set is found, politely move it
    if os.path.isfile(working_dir+fname): shutil.move(working_dir+fname, working_dir+fname+".arch")

    # we will carry most of the attributes, but "index" is a new var, so we created some in define_var_atts()
    var_atts = define_var_atts()
    for value in var_atts.values(): file[value[0]].attrs = {} # clear current atts
    for value in var_atts.values(): file[value[0]].attrs[value[1]] = value[2] # replace with desired atts

    # global attributes
    global_atts = define_glob_atts()
    for key, value in global_atts.items(): file.attrs[key] = value

    # write back to new file
    file.to_netcdf(working_dir+fname)
    file.close()

    # delete the intermediary files
    #for filename in glob.glob(working_dir+"tmp*.nc"):
    #    os.remove(filename)


# executes main():
if __name__ == "__main__":
    main()
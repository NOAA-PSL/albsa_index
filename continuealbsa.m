
cd /Users/ccox/Documents/Projects/2024/albsa/

laS = 50;
loS = 360-170;
laN = 75;
loN = 360-170;
laE = 55;
loE = 360-150;
laW = 55;
loW = 160;

d = dir('*.nc');
file   = rd_netcdf(d(1).name);
indlon = find(file.lon >= 150 & file.lon <= 240);
indlat = find(file.lat >= 40 & file.lat <= 85);

albsa.i   = NaN(length(d).*366*4,1);
albsa.ew  = NaN(length(d).*366*4,1);
albsa.ns  = NaN(length(d).*366*4,1);
albsa.N   = NaN(length(d).*366*4,1);
albsa.S   = NaN(length(d).*366*4,1);
albsa.E   = NaN(length(d).*366*4,1);
albsa.W   = NaN(length(d).*366*4,1);
albsa.h   = NaN(length(d).*366*4,1);
albsa.d   = NaN(length(d).*366*4,1);
albsa.yd  = NaN(length(d).*366*4,1);
albsa.m   = NaN(length(d).*366*4,1);
albsa.y   = NaN(length(d).*366*4,1);
albsa.map = NaN(length(indlon),length(indlat),length(d).*366);
albsa.lat = file.lat(indlat);
albsa.lon = file.lon(indlon);

count = 1;

for k = 1:length(d)
    
    disp(d(k).name);
    
    file  = rd_netcdf(d(k).name);
    a     = attribute_list(d(k).name);
    
    years = str2num(datestr(double(file.time./24 + datenum([1800 01 01 00 00 00])),'yyyy'));
    mos   = str2num(datestr(double(file.time./24 + datenum([1800 01 01 00 00 00])),'mm'));
    days  = str2num(datestr(double(file.time./24 + datenum([1800 01 01 00 00 00])),'dd'));
    hours  = str2num(datestr(double(file.time./24 + datenum([1800 01 01 00 00 00])),'hh'));
    doys  = double(file.time./24 + datenum([1800 01 01 00 00 00])) - datenum([years years.*0+1 years.*0 years.*0 years.*0 years.*0]);
    len   = length(mos);

    Scoord = [];
    for la = laS
        for lo = loS
            Scoord = [Scoord; [findnearest(file.lat,la,1) findnearest(file.lon,lo,1)]];
        end
    end
    Ncoord = [];
    for la = laN
        for lo = loN
            Ncoord = [Ncoord; [findnearest(file.lat,la,1) findnearest(file.lon,lo,1)]];
        end
    end
    Ecoord = [];
    for la = laE
        for lo = loE
            Ecoord = [Ecoord; [findnearest(file.lat,la,1) findnearest(file.lon,lo,1)]];
        end
    end
    Wcoord = [];
    for la = laW
        for lo = loW
            Wcoord = [Wcoord; [findnearest(file.lat,la,1) findnearest(file.lon,lo,1)]];
        end
    end

    Nmean = squeeze(file.hgt(Ncoord(2),Ncoord(1),3,:));
    Smean = squeeze(file.hgt(Scoord(2),Scoord(1),3,:));
    Emean = squeeze(file.hgt(Ecoord(2),Ecoord(1),3,:));
    Wmean = squeeze(file.hgt(Wcoord(2),Wcoord(1),3,:));
    albsa.i(count:count+len-1)       = (Emean-Wmean) - (Nmean-Smean);
    albsa.ew(count:count+len-1)      = (Emean-Wmean);
    albsa.ns(count:count+len-1)      = (Nmean-Smean);
    albsa.N(count:count+len-1)       = Nmean;
    albsa.S(count:count+len-1)       = Smean;
    albsa.E(count:count+len-1)       = Emean;
    albsa.W(count:count+len-1)       = Wmean;
    albsa.d(count:count+len-1)       = days;
    albsa.h(count:count+len-1)       = hours;
    albsa.yd(count:count+len-1)      = doys;
    albsa.m(count:count+len-1)       = mos;
    albsa.y(count:count+len-1)       = years;
    albsa.map(:,:,count:count+len-1) = squeeze(file.hgt(indlon,indlat,3,:));
    clear Emean Wmean Nmean Smean EmeanV WmeanV NmeanV SmeanV
    count = count + len;

end
albsa.i(count:end)       = [];
albsa.ew(count:end)      = [];
albsa.ns(count:end)      = [];
albsa.N(count:end)       = [];
albsa.S(count:end)       = [];
albsa.E(count:end)       = [];
albsa.W(count:end)       = [];
albsa.h(count:end)       = [];
albsa.d(count:end)       = [];
albsa.yd(count:end)      = [];
albsa.m(count:end)       = [];
albsa.y(count:end)       = [];
albsa.map(:,:,count:end) = [];
albsa.dn = datenum([albsa.y albsa.m albsa.d albsa.h albsa.h.*0 albsa.h.*0]);

albsanew = albsa;

load('/Users/ccox/Documents/Projects/2024/albsa/albsadaily.mat')

albsa.dn = datenum([albsa.y albsa.m albsa.d albsa.d.*0 albsa.d.*0 albsa.d.*0]);

fnames = fieldnames(albsa);

ind = find(albsanew.dn > albsa.dn(end));

for k = 1:length(fnames)
    
    if strcmp(char(fnames(k)),'map')
        albsa.(char(fnames(k))) = cat(3,albsa.(char(fnames(k))),albsanew.(char(fnames(k)))(:,:,ind));
    elseif strcmp(char(fnames(k)),'lat') | strcmp(char(fnames(k)),'lon')
        % do nothing
    else
        albsa.(char(fnames(k))) = [albsa.(char(fnames(k))); albsanew.(char(fnames(k)))(ind)];
    end

end

clearallbut albsa






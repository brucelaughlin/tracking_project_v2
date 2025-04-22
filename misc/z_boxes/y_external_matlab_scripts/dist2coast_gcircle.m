function dist2cst=dist2coast_gcircle(lonr,latr,maskr,lonc,latc);
% dist2coast_gcircle: finds the shortest distance to coast 
%                     (defined by a coastline) using ROMS gcircle.m routine.
%
% SYNTAX
% ------
% dist2cst=dist2coast_gcircle(lonr,latr,maskr,lonc,latc);
%
% where lonr,latr are the lon/lat points of the grid 
%       maskr is the mask of the grid (used to only find dist for ocean points
%       lonc,latc define the coastline (can and should probably have NaNs in it)
%  and  dist is what's returned (in meters);
%
% Chris Edwards		                          June 2018

% get cells that have coastline in them
[ni,nj]=size(lonr);

dist2cst=NaN*zeros([ni nj]);

% find minimum distance to that set of coastlines
[iocean,jocean]=find(maskr==1);
nocean=length(jocean);
for kocean=1:nocean
  j0=jocean(kocean);
  i0=iocean(kocean);
  lon0=lonr(i0,j0);
  lat0=latr(i0,j0);
  distt = gcircle(lon0,lat0,lonc,latc);
  dist2cst(i0,j0)=min(distt);
end

return

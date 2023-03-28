import numpy as np

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE

def ll2xy(lon, lat, lon0, lat0):
    """
    This converts lon, lat into meters relative to lon0, lat0.
    It should work for lon, lat scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y

def lon2x(lon, lon0, lat0):
    """
    This converts lon into meters relative to lon0, along lat0.
    It should work for lon scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    return x

def lat2y(lat, lat0):
    """
    This converts lat into meters relative to lat0.
    It should work for lon, lat scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    y = R * np.pi * (lat - lat0) / 180
    return y

def x2lon(x, lon0, lat0):
    """
    This converts a distance in meters to a longitude relative to lon0, along lat0.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    lon = (180 * x) / (R * clat * np.pi) + lon0
    return lon

def y2lat(y, lat0):
    """
    This converts a distance in meters to a latitude relative to lat0.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    lat = (180 * y) / (R * np.pi) + lat0
    return lat

def xy2ll(x, y, lon0, lat0):
    """
    This converts a position x,y in meters into lon, lat relative to lon0, lat0.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    lon = (180 * x) / (R * clat * np.pi) + lon0
    lat = (180 * y) / (R * np.pi) + lat0
    return lon, lat
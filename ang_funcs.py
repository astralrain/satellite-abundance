from libs import *
from consts import *

def angular_separation(ra1: float, dec1: float, ra2: float, 
                       dec2: float, degrees=True) -> float:
    """
    Calculate angular separation of coordinates (ra1, dec1) and (ra2, dec2) using the
        haversine formula.
    
    Parameters:
        ra1 (float): Right ascension of the first location.
        dec1 (float): Declination of the first location.
        ra2 (float): Right ascension of the second location.
        dec2 (float): Declination of second location.
        degrees (bool): Whether the values were given in degrees, default True.
    
    Returns:
        c (float): Angular separation in degrees.
    """
    if degrees:
        ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
        
    dra = ra2 - ra1
    # wrap ra 
    dra = (dra + np.pi) % (2 * np.pi) - np.pi  
    # eg. 1 degree and 359 degrees are only 2 degrees apart
    # and not 358 degrees apart
    
    ddec = dec2 - dec1
        
    a = np.sin(ddec / 2.0)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra / 2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.degrees(c) if degrees else c

def angular_diameter_distance(z: float) -> float:
    """
    Calculates the angular diameter distance for a given redshift z
        in Plank cosmology.

    Parameters:
        z (float): Redshift of the object.

    Returns:
        float: Angular diameter distance.
    """
    return cosmo.angular_diameter_distance(z).to(u.Mpc).value

def projected_distance(sep_deg: float, ang_dia_dist: float) -> float:
    """
    Convert angular separation, sep_deg to physical projected distance
        using angular diameter distance, ang_dia_dist.
    
    Parameters:
        separation_deg (float): Angular separation in degrees.
        ang_dia_dist (float): Angular diameter distance in Mpc.
    
    Returns:
        float: Projected distance in Mpc.
    """
    return ang_dia_dist * np.radians(sep_deg)
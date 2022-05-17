
# -*- coding: utf-8 -*-
"""
This scripts contains all the needed T functions

@author: zhaoliang
"""

import numpy as np

def arrangeT(waypts,T):
    """
    This function is to calculate the initial 
        time series for each waypoint section.

    Parameters
    ----------
    waypts : 2d np.array
        2d way points.
    T : float
        Total execution time.

    Returns
    -------
    ts : np.array[float]
        Time to reach each waypoint.

    """
    x = waypts[:,1:] - waypts[:,0:-1] 
    
    dist = sum(x**2)**0.5
    k = T/sum(dist)
    new = np.cumsum(dist*k)
    ts = np.insert(new,0,0)
    
    return ts

def init_T(waypts,ave_speed):
    """
    This functino is to initialize total estimated T
    The idea in this function is to assign T by 
        using total path length / average speed 

    Parameters
    ----------
    waypts :  2 * n np.array
        A series of 2d waypoints.
    ave_speed : float
        A estimated average speed.

    Returns
    -------
    T : float 
        Initial total estimated time.

    Example
    -------
    waypts = np.array([[1,2,3],[1,2,3]])
    ave_speed = 2 # m/s
    
    """
    x = waypts[:,1:] - waypts[:,0:-1] 
    dist = sum(x**2)**0.5
    total_dist = sum(dist)
    T = total_dist/ave_speed
    return T
# wave-interp

A MATLAB data structure for pre-computed free fluid surface profiles, _η(x,y,t)_, and the underlying fluid particle kinematics, _**u**(x,y,z,t)_, complete with functions to interpolate between the pre-computed data points.

Released under the GNU General Public License v3.0. In addition, the author reserves the rights to later publications of the project under a different licence. All usage must incorporate appropriate citation as listed below.

Disclaimer: the codes are provided 'as-is' without warranty and is used solely at the risk of the end-user.

Copyright (c) 2013-2020 lm808. All rights reserved.

## Cite as

L. Ma, wave-interp (2020), open-source repository, https://github.com/lm808/wave-interp

## Input data structure

The input data file is a \*.mat file for MATLAB R2014b or newer. The Cartesian coordinate convention is such that:
* _x_ is the mean direction of wave propagation.
* _y_ is the transverse direction along the wave crest
* _z_ is the vertical direction, which is positive upwards and originates at the still water level (SWL).

**1. Shared data**

Variable | Size | Type | Description
------------- | ---- | ---- | -----------
`d` | 1 * 1 | double | water depth with respect to the SWL
`nt` | 1 * 1 |  numeric | number of time-steps
`t` | `nt` * 1 |  double | time vector

**2. Surface elevation data**

Variable | Size | Type | Description
------------- | ---- | ---- | -----------
`nx` | 1 * 1 | numeric | number of data points in the *x*-direction
`ny` | 1 * 1 | numeric | number of data points in the *y*-direction
`X` | `ny` * `nx` | double | array of *x*-coordinates (`meshgrid` layout)
`Y` | `ny` * `nx` | double | array of *y*-coordinates (`meshgrid` layout)
`ETA` | `ny` * `nx` * `nt` | double | surface elevation data (`meshgrid` layout)
`maxEta` | 1 * 1 | double | maximum wave crest elevation
`minEta` | 1 * 1 | double | minimum free surface elevation

**3. Underlying fluid particle kinematics**

Variable | Size | Type | Description
------------- | ---- | ---- | -----------
`x` | 1 * `nt` | cell array | *x*-location vectors for kinematics
`y` | 1 * `nt` | cell array | *y*-locations vectors for kinematics
`z` | 1 * `nt` | cell array | *z*-locations vectors for kinematics
`ux` | 1 * `nt` | cell array | vectors containing velocity in the *x*-direction
`uy` | 1 * `nt` | cell array | vectors containing velocity in the *y*-direction
`uz` | 1 * `nt` | cell array | vectors containing velocity in the *z*-direction

For the `x`, `y`, `z`, `ux`, `uy` and `uz` cell arrays, each cell element corresponds to a time-step. For each time-step the cell should contain a simple 1-D vector. Taken together, they describe the (*ux*, *uy*, *uz*) velocities for a given set of (*x*, *y*, *z*) points. As such, for the same time-step (or cell address) the vector sizes should be identical across all 6 cell-arrays.

**Additional notes**

* In general, it is a good convention to shift the highest wave crest to occur at `x=0`,`y=0`, and `t=0`. In this way it is very easy to locate the maximum.
* To improve speed, do not include data ranges beyond what is absolutely required.
* To ensure correct pressure calculations based on the interpolated velocity for multi-valued free surfaces (e.g. an over-turning wave):
  * Provide the highest elevation where water exists in `ETA` (i.e. an envelope).
  * Set velocities to zero for air-filled cells in `ux`, `uy` and `uz`.

## Limitation

Interpolation cannot and should not be performed in between pre-computed time-steps to avoid unnecessary errors or the propagation of NaN's. This is enforced by an _1e-10_ tolerance check on the requested time against the pre-computed time values. As such, it is advisable to extract the time vector directly from the data file and index into it to retrieve required time values.

The interpolation, however, can be performed at any spatial location within a given time-step; provided that the requested point _(xq,yq,zq)_ is within the convex hull defined by the `x`, `y` and `z` vectors for each time-step.

## List of functions

`fInterpEta.m`

Returns the free surface elevation η for given _(x,y)_ locations, the time-step, _t_, and the input data file.

`fInterpVel.m`

Returns the three-component velocities _**u**=[ux,uy,uz]_ for given _(x,y,z)_ locations, the time-step, _t_, and the input data file.

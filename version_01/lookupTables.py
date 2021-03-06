"""
Copyright 2017, 2017 nowaythatsok

This file is part of GNU_lapser 0.1.

    GNU_lapser is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GNU_lapser is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GNU_lapser.  If not, see <http://www.gnu.org/licenses/>.
"""

#this file contains the extracted look up tables

#-------------------------------Official INCLUDES-------------------------------
import numpy as np

#-------------------------------Function DEFINITIONS----------------------------



def percentile_61_3ColorLookUp(pxl):
#D3s
	x = np.array( [0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.200, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.050, 1.100, 1.200, 1.300, 1.350, 1.400, 1.500, 1.550, 1.600, 1.650, 1.700, 1.750, 1.800, 1.850, 1.900, 1.950, 2.000, 2.050, 2.100, 2.133, 2.167, 2.200, 2.250, 2.300, 2.333, 2.367, 2.400, 2.433, 2.467, 2.500, 2.550, 2.600, 2.633, 2.667, 2.700, 2.733, 2.767, 2.800, 2.833, 2.867, 2.900, 2.933, 2.967, 3.000, 3.025, 3.050, 3.075, 3.100, 3.133, 3.167, 3.200, 3.225, 3.250, 3.275, 3.300, 3.333, 3.367, 3.400, 3.425, 3.450, 3.475, 3.500, 3.525, 3.550, 3.575, 3.600, 3.620, 3.640, 3.660, 3.680, 3.700, 3.725, 3.750, 3.775, 3.800, 3.820, 3.840, 3.860, 3.880, 3.900, 3.920, 3.940, 3.960, 3.980, 4.000, 4.020, 4.040, 4.060, 4.080, 4.100, 4.120, 4.140, 4.160, 4.180, 4.200, 4.217, 4.233, 4.250, 4.267, 4.283, 4.300, 4.317, 4.333, 4.350, 4.367, 4.383, 4.400, 4.420, 4.440, 4.460, 4.480, 4.500, 4.517, 4.533, 4.550, 4.567, 4.583, 4.600, 4.617, 4.633, 4.650, 4.667, 4.683, 4.700, 4.717, 4.733, 4.750, 4.767, 4.783, 4.800, 4.820, 4.840, 4.860, 4.880, 4.900, 4.917, 4.933, 4.950, 4.967, 4.983, 5.000, 5.020, 5.040, 5.060, 5.080, 5.100, 5.117, 5.133, 5.150, 5.167, 5.183, 5.200, 5.220, 5.240, 5.260, 5.280, 5.300, 5.320, 5.340, 5.360, 5.380, 5.400, 5.417, 5.433, 5.450, 5.467, 5.483, 5.500, 5.520, 5.540, 5.560, 5.580, 5.600, 5.620, 5.640, 5.660, 5.680, 5.700, 5.725, 5.750, 5.775, 5.800, 5.825, 5.850, 5.875, 5.900, 5.925, 5.950, 5.975, 6.000, 6.025, 6.050, 6.075, 6.100, 6.133, 6.167, 6.200, 6.233, 6.267, 6.300, 6.333, 6.367, 6.400, 6.433, 6.467, 6.500, 6.533, 6.567, 6.600, 6.650, 6.700, 6.750, 6.800, 6.833, 6.867, 6.900, 7.000, 7.050, 7.100, 7.150, 7.200, 7.233, 7.267, 7.300, 7.350, 7.400, 7.500, 7.600, 7.700, 7.800, 8.000, 8.400, 8.600, 9.100, 9.100] )

	return x[pxl]


def percentile_61_3ColorLookUp_n(pxl):
#D3100
	x = np.array( [0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.200, 0.250, 0.300, 0.400, 0.500, 0.621, 0.705, 0.800, 0.867, 0.920, 0.961, 1.003, 1.068, 1.134, 1.203, 1.286, 1.354, 1.411, 1.447, 1.483, 1.523, 1.567, 1.611, 1.657, 1.702, 1.751, 1.800, 1.837, 1.874, 1.912, 1.952, 1.992, 2.030, 2.068, 2.108, 2.168, 2.217, 2.251, 2.286, 2.319, 2.352, 2.385, 2.414, 2.440, 2.466, 2.492, 2.523, 2.556, 2.589, 2.624, 2.661, 2.698, 2.732, 2.765, 2.799, 2.826, 2.853, 2.880, 2.907, 2.934, 2.962, 2.989, 3.017, 3.045, 3.073, 3.101, 3.126, 3.150, 3.174, 3.199, 3.226, 3.253, 3.280, 3.307, 3.337, 3.366, 3.395, 3.418, 3.439, 3.460, 3.481, 3.503, 3.524, 3.545, 3.566, 3.587, 3.610, 3.635, 3.659, 3.684, 3.709, 3.732, 3.756, 3.780, 3.803, 3.822, 3.841, 3.860, 3.879, 3.898, 3.918, 3.938, 3.958, 3.978, 3.998, 4.017, 4.035, 4.053, 4.071, 4.089, 4.108, 4.129, 4.149, 4.169, 4.189, 4.210, 4.232, 4.253, 4.275, 4.297, 4.315, 4.332, 4.350, 4.367, 4.385, 4.403, 4.423, 4.443, 4.463, 4.482, 4.502, 4.521, 4.539, 4.558, 4.577, 4.596, 4.615, 4.634, 4.653, 4.672, 4.692, 4.711, 4.730, 4.750, 4.769, 4.788, 4.807, 4.826, 4.846, 4.865, 4.884, 4.903, 4.923, 4.942, 4.962, 4.981, 5.001, 5.022, 5.042, 5.063, 5.083, 5.104, 5.125, 5.147, 5.168, 5.189, 5.209, 5.227, 5.245, 5.264, 5.282, 5.300, 5.319, 5.338, 5.357, 5.376, 5.395, 5.415, 5.436, 5.457, 5.478, 5.499, 5.521, 5.542, 5.563, 5.585, 5.606, 5.628, 5.649, 5.671, 5.692, 5.714, 5.736, 5.759, 5.781, 5.804, 5.833, 5.861, 5.889, 5.915, 5.939, 5.964, 5.988, 6.014, 6.041, 6.069, 6.096, 6.128, 6.160, 6.192, 6.230, 6.271, 6.308, 6.336, 6.364, 6.393, 6.437, 6.487, 6.530, 6.571, 6.622, 6.702, 6.753, 6.804, 6.854, 6.912, 7.027, 7.103, 7.169, 7.234, 7.297, 7.400, 7.595, 7.700, 7.750, 7.800, 7.883, 8.108, 8.300, 8.600, 8.600] )

	return x[pxl]

"""Simple Equal Earth projection realisation in Python"""

import numpy as np
import scipy.optimize as spo
import matplotlib.pylab as plt

__author__ = "Ilya Oshchepkov"
__copyright__ = "Copyright 2018, Ilya Oshchepkov"
__credits__ = ["Ilya Oshchepkov"]
__license__ = "GPL3"
__version__ = "1.0"
__maintainer__ = "Ilya Oshchepkov"
__status__ = "Prototype"


class EqualEarth:
    """
    Class for working with the Equal Earth projection.

    Parameters
    ----------
    auth_radius : float
        Radius of the authlatic sphere, in metres.

    Methods
    -------
    latlon_to_xy(lot, lan, degrees=True)
        Convert geocentric coordinates to Equal Earth coordinates.
    xy_to_latlon(x, y, degrees=True, tol=1e-9)
        Convert Equal Earth coordinates to geocentric coordinates.
    plot_grid(ax=None, meridians=None, parallels=None)
        Init matplotlib axis for plotting in the Equal Earth projection.

    References
    ----------
    .. [1] Bojan Šavrič, Tom Patterson & Bernhard Jenny (2018)
        The Equal Earth map projection, International
        Journal of Geographical Information Science,
        DOI: 10.1080/13658816.2018.1504949
    """
    a1 = 1.340264
    a2 = -0.081106
    a3 = 0.000893
    a4 = 0.003796
    m = np.sqrt(3) / 2
    tol = 1e-9

    def __init__(self, auth_radius=1.0):
        self.auth_radius = auth_radius

    def _x(self, lat):
        """Auxiliary function for x.
        """
        x1 = 2 / np.sqrt(3) * np.cos(lat)
        x2 = self.a1 + 3*self.a2*lat**2 + \
            lat**6 * (7*self.a3 + 9*self.a4*lat**2)
        return x1 / x2

    def _y(self, lat):
        """Auxiliary function for y.

        """
        return lat * (self.a1 + self.a2*lat**2 + lat**6 * (self.a3 + self.a4*lat**2))

    def latlon_to_xy(self, lat, lon, degrees=True):
        """Convert geocentric coordinates to Equal Earth coordinates.

        Parameters
        ----------
        lat, lon : float or array_like of floats
            Geocentric latitude and longitudes, in degrees or radians.
        degrees : bool, optional
            If True, the input geocentric ccordinates are given
            in degrees, otherwise radians. Default is True.

        Returns
        -------
        x, y : float or array_like of floats
            Planar coordinates in the Equal Earth projection, in metres.
        """
        if degrees:
            lat = np.radians(lat)
            lon = np.radians(lon)

        lat = np.arcsin(self.m * np.sin(lat))

        x = self.auth_radius*lon*self._x(lat)
        y = self.auth_radius*self._y(lat)

        return x, y

    def xy_to_latlon(self, x, y, degrees=True):
        """Convert Equal Earth coordinates to geocentric coordinates.

        Parameters
        ----------
        x, y : float or array_like of floats
            Planar coordinates in the Equal Earth projection, in metres.
        degrees : bool, optional
            If True, the output geocentric ccordinates will be given
            in degrees, otherwise radians. Default is True.

        Returns
        -------
        lat, lon : float or array_like of floats
            Geocentric latitude and longitudes
        """
        def func(lat): return self._y(lat) - y / self.auth_radius

        lat = spo.root(func, y / self.auth_radius, tol=self.tol).x
        lon = x*self._x(lat)**-1

        lat = np.arcsin(self.m**-1 * np.sin(lat))

        if degrees:
            lat = np.degrees(lat)
            lon = np.degrees(lon)

        return lat, lon

    def plot_grid(self, ax=None, meridians=None, parallels=None):
        """Init matplotlib axis for plotting in the Equal Earth projection.

        """
        if ax is None:
            fig, ax = plt.subplots()

        ax.axis('off')
        ax.set_aspect('equal')

        if parallels is None:
            parallels = np.linspace(-90, 90, 7, endpoint=True)
        for lat in parallels:
            x, y = prj.latlon_to_xy([lat, lat], [-180, +180])
            ax.plot(x, y, 'k', linewidth=0.5)

        if meridians is None:
            meridians = np.linspace(-180, 180, 13, endpoint=True)

        for lon in meridians:
            p = np.linspace(-90, 90, 100)
            x, y = prj.latlon_to_xy(p, lon)
            ax.plot(x, y, 'k', linewidth=0.5)

        x_w, x_e = self.latlon_to_xy(0.0, [-180, 180])[0]
        y_n, y_s = self.latlon_to_xy([-90, 90], 0.0)[1]

        ax.set_xlim(x_w, x_e)
        ax.set_ylim(y_n, y_s)

        return ax

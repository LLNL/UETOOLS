class Utilities:
    """Class providing useful utilities for UETOOLS

    These methods are based on functions provided by UETOOLS,
    mainly Case.get and Case.setue. To utilize these outside of
    UETOOLS a class accessing UEDGE data by a get and setue method
    has to be duck-typed.

    Methods
    -------
    calc_lcon()
        stores the connection length in Case.lcon
    potent_1dsol()
        sets the theoretical 1D potential for the UEDGE case
    squareinterp(data, r=None, z=None, method='linear',
            resolution=(500j, 800j) mask=False, fill=float('NaN')
        interpolates data from UEDGE grid to square grid
    psinormc(simagx=None, sibdrys=None)
        returns the normalized-psi values at the OMP cell centers
    psinormf(simagx=None, sibdrys=None)
        returns the normalized-psi values at the OMP east cell faces
    """

    def __init__(self, case):
        """Connects UeUtilitie.Utilities to Case"""
        self.get = case.get
        self.setue = case.setue

    def calc_lcon(self):
        """Stores connection length in Case.lcon"""
        from numpy import cumsum

        geom = self.get("geometry")[0].decode("UTF-8").strip()
        if geom in ["snull", "uppersn"]:
            # TODO: How to extend to dnulls????
            # Calc and store connection lengths
            lcon = cumsum(1 / (self.get("rr") * self.get("gx") + 1e-20), axis=0)
            # No connlen for core cells
            lcon[
                self.get("ixpt1")[0] + 1 : self.get("ixpt2")[0] + 1,
                : self.get("iysptrx") + 1,
            ] = 0
            lcon[self.get("ixpt2")[0] + 1 :, : self.get("iysptrx") + 1] = cumsum(
                1
                / (self.get("rr") * self.get("gx") + 1e-20)[
                    self.get("ixpt2")[0] + 1 :, : self.get("iysptrx") + 1
                ],
                axis=0,
            )
            return lcon
        else:
            raise NotImplementedError(
                "Connection lengt calculation for "
                + "geometry {} not yet implemented".format(geom)
            )

    def potent_1dsol(self):
        """Sets the theoretical 1D potential for the UEDGE case"""
        # TODO: How to extend to dnulls????
        from numpy import zeros_like

        phi = self.get("phi")
        phis = zeros_like(phi)
        gx = self.get("gx")
        ixp1 = self.get("ixp1")
        ex = self.get("ex")

        phis[1, :] = self.get("kappal")[:, 0] * self.get("te")[1] / self.get("qe")
        if self.get("isudsym") == 0:
            for iy in range(1, self.get("ny") + 1):
                for ix in range(1, self.get("nx") + 1):
                    ix1 = ixp1[ix, iy]
                    dxf = 0.5 * (gx[ix, iy] + gx[ix1, iy]) / (gx[ix, iy] * gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf

            for ix in range(self.get("ixpt1")[0] + 1, self.get("ixpt2")[0] + 1):
                for iy in range(self.get("iysptrx") + 1):
                    phis[ix, iy] = phis[ix, self.get("iysptrx") + 1]
                phis[ix, self.get("ny") + 1] = phis[ix, self.get("ny")]
        else:
            nxc = self.get("nxc")
            for iy in range(1, self.get("ny") + 1):
                for ix in range(0, nxc):
                    ix1 = ixp1[ix, iy]
                    dxf = 0.5 * (gx[ix, iy] + gx[ix1, iy]) / (gx[ix, iy] * gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf
                # NOTE: no overlap for the two subdomains in poloidal direction
                for ix in range(self.get("nx"), nxc, -1):
                    ix1 = ixm1[ix, iy]
                    ix2 = ixp1[ix, iy]
                    dxf = 0.5 * (gx[ix, iy] + gx[ix1, iy]) / (gx[ix, iy] * gx[ix1, iy])
                    phis[ix, iy] = phis[ix2, iy] - ex[ix, iy] * dxf

        self.setue("phis", phis)
        self.setue("phi", phis)

    def squareinterp(
        self,
        data,
        r=None,
        z=None,
        method="linear",
        resolution=(500j, 800j),
        mask=False,
        fill=float("NaN"),
        zshift = 0,
    ):
        """Interpolates date from UEDGE grid to uniform square

        Arguments
        ---------
        data - 2D array of values to interpolate on UEDGE grid

        Keyword arguments
        -----------------
        r : len-2 tuple of floats (default = None)
            tuple of (r_begin, r_end) R-coordinates of new grid. If
            None, bases range on rm
        z : len-2 tuple of floats (default = None)
            tuple of (z_begin, z_end) Z-coordinates of new grid. If
            None, bases range on zm
        method : str (default = 'linear')
            intepolation method, passed to scipy.interpolate.griddata
        resolution : len-2 tuple of im ints (default = (500j, 800j))
            r and z resoltuion of interpoalted grid
        mask : bool (default = False)
            if True, masks out points external to the UEDGE grid
        fill : float (default = float('NaN'))
            What to fill masked/extrapolated cells with

        Returns
        -------
        (gx, gy, interp)
        gx - uniform mgrid data of x-coordinates for resolution
        gy - uniform mgrid data of y-coordinates for resolution
        interp - interpolated data on new grid dedined by resolution
        """
        from scipy.interpolate import griddata
        from numpy import floor, ceil, mgrid, array, nan, concatenate
        from shapely.geometry import Point
        from shapely.geometry.polygon import Polygon
        from copy import deepcopy

        # Get R,Z coords
        rm = self.get("rm")
        zm = self.get("zm")
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "     uppersn":
            zm = self.disp - zm
        zm = deepcopy(zm)
        zm += zshift
        if mask:
            # Mask out external points
            x = []
            y = []
            for j in range(rm.shape[1]):
                x.append(rm[1, j, 1])
                y.append(zm[1, j, 1])
            for i in range(rm.shape[0]):
                x.append(rm[i, -1, 1])
                y.append(zm[i, -1, 1])
            for j in range(rm.shape[1] - 1, -1, -1):
                x.append(rm[-1, j, 1])
                y.append(zm[-1, j, 1])
            for i in range(rm.shape[0] - 1, self.get("ixpt2")[0], -1):
                x.append(rm[i, 0, 4])
                y.append(zm[i, 0, 4])
            for i in range(self.get("ixpt1")[0], -1, -1):
                x.append(rm[i, 0, 4])
                y.append(zm[i, 0, 4])
            if self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
                y = [self.disp - a for a in y]
            outline = Polygon(concatenate([[array(x), array(y)]], axis=0).T)

        # Automatically set sqare grid boundaries
        if r is None:
            r = (0.01 * floor(rm.min() * 100), 0.01 * ceil(rm.max() * 100))
        if z is None:
            z = (0.01 * floor(zm.min() * 100), 0.01 * ceil(zm.max() * 100))
        # Get cell centers and unravel data
        gx, gy = mgrid[r[0] : r[1] : resolution[0], z[0] : z[1] : resolution[1]]

        interp = griddata(
            (rm[:, :, 0].ravel(), zm[:, :, 0].ravel()),
            data.ravel(),
            (gx, gy),
            method=method,
            fill_value=fill,
        )
        if mask:
            for i in range(gx.shape[0]):
                for j in range(gx.shape[1]):
                    if not outline.contains(Point(gx[i, j], gy[i, j])):
                        interp[i, j] = nan

        # Perform interpolation
        return gx, gy, interp

    def psinormc(self, simagxs=None, sibdrys=None):
        """Returns normalized-psi values at the OMP cell centers

        Keyword arguments
        -----------------
        simagxs : float (default = None)
            the psi-value at the magnetic axis. If None, accesses
            simagxs value in UEDGE memory
        sibdrys : float (default = None)
            the psi-value at the inner separatrix. If None, accesses
            simagxs value in UEDGE memory

        Returns
        -------
        1D ndarray of normalized cell-center psi values at OMP
        """
        if simagxs is None:
            simagxs = self.get("simagxs")
        if sibdrys is None:
            sibdrys = self.get("sibdrys")
        return (self.get("psi")[self.get("ixmp"), :, 0] - simagxs) / (sibdrys - simagxs)

    def psinormf(self, simagxs=None, sibdrys=None):
        """Returns normalized-psi values at the OMP east cell faces

        Keyword arguments
        -----------------
        simagxs : float (default = None)
            the psi-value at the magnetic axis. If None, accesses
            simagxs value in UEDGE memory
        sibdrys : float (default = None)
            the psi-value at the inner separatrix. If None, accesses
            simagxs value in UEDGE memory

        Returns
        -------
        1D ndarray of normalized east cell face psi values at OMP
        """
        if simagxs is None:
            simagxs = self.get("simagxs")
        if sibdrys is None:
            sibdrys = self.get("sibdrys")
        psi = self.get("psi")[self.get("ixmp")]
        return ((0.5 * (psi[:, 3] + psi[:, 4])) - simagxs) / (sibdrys - simagxs)

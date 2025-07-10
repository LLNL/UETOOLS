class VacuumTests:
 
    def twoSurfacePlot(self):
        S1 = Surface((2, 5), (4, 2), 0)
        S2 = Surface((2, 1), (3, 2), 1)
        S1.showTwoSurfacePlot(S2, r_offset=1)
    
    def outerCirclePlot(self):
        S1 = Surface((2, 5), (4, 1,), 0)
        S1.showOuterCirclePlot(r_offset=1)

    def analyticPlot(self):
        S1 = Surface((2, 5), (4, 1), 0)
        S1.showAnalyticPlot(True, r_offset=1)
    
    def trianglePlot(self):
        from shapely import Point, LineString
        import math
        import numpy as np

        S1 = Surface((1, 3), (2, 6), 0)
        height = math.sqrt(3) * S1.surfaceLength / 2 # height of the et
        vertex = S1.normal.interpolate(height)
        test = VacuumRegion( [
            (self.start.x, self.start.y), 
            (self.end.x, self.end.y), 
            (vertex.x, vertex.y)
        ])

    def squarePlot(self):
        from shapely import Point, LineString
        S1 = Surface((1, 2), (3, 6), 0)
        # # # Make the sides of the square perpendicular to self # # #
        side1Start = (S1.end.x, S1.end.y)
        side1End = S1.normalHelper(abs(S1.dx), abs(S1.dy), side1Start[0], side1Start[1], False)

        side3End = (S1.start.x, S1.start.y)
        side3Start = self.normalHelper(abs(S1.dx), abs(S1.dy), side3End[0], side3End[1], False)

        test = VacuumRegion([side1Start, side1End, side3Start, side3End])

    def shadedSquarePlot(self):
        from shapely import Point, LineString

        """To be used alongside the source surface S1-- Surface((2, 1), (1, 1))-- which is defined in the lineOfSightPlot function
        in the test functions."""

        geometryVertices = [(3, 1), (1, 1), (1, 5), (2, 6), (1, 7), (7, 7), (7, 1), (5, 1), (5, 3), (3, 3)]

        test = VacuumRegion(geometryVertices, P=1)
        for i in test.errors:
            f = test.plotGeometry(labels=True, testsurf=i)
            test.surfaces[i].plotSelf(ax=f)
            test.surfaces[i].printReport()
        f = test.plotGeometry(labels=True, testsurf=6)
        return test


    def tokamakPlot(self, savefile):
        from uetools import Case
        from numpy import zeros
        c = Case(savefile, inplace=True)
        (main, pf) = c.coupling.get_snull_vacuum_regions(maxlength = 1)
        nobug = zeros((main[0].shape[0]-1, main[0].shape[1]))
        nobug[:66] = main[0][:66]
        nobug[66:] = main[0][67:]
        test = VacuumRegion(nobug, P=main[1])
        for i in test.errors:
            if (i > 50) and (i<90):
                f = test.plotGeometry(labels=False, testsurf=i, markers='.')
                f.get_axes()[0].set_title(f"Surface {i}")
        f = test.plotGeometry(labels=False, testsurf=60, showCircle=True)
        f = test.plotGeometry(labels=False, testsurf=4)
        return test
       


class VacuumRegion:
    def __init__(self, nodeList, P=0):
        from shapely import Point, Polygon
        from tqdm import tqdm
        # Generate all surfaces and create dictionary
        # Create neigbors dictionaries by checking LOS for each surface
        # Dictionary containing all surfaces making up the geometry
        self.surfaces = {}
        # # # Set up surfaces of geometry and the polygon object # # #
        for i in range(len(nodeList)):
            startNode = Point(nodeList[i])
            if i == len(nodeList) - 1:
                endNode = Point(nodeList[0])
            else:
                endNode = Point(nodeList[i + 1])
            self.surfaces[i] = Surface((startNode.x, startNode.y), (endNode.x, endNode.y), i)

            # later if n < i < m: s.reflecting = 0, something like this
 
        # Create Polygon of Vacuum region for intersect checks
        self.geometry = Polygon(nodeList) 
        self.P = P # Number of plasma surfaces
    
        # Iterate surfaces to identify surface neigbors
        for _, surface in tqdm(self.surfaces.items()):
            surface.getNeighbors(self.surfaces, self.geometry)

        if not self.checkContinuity(False):
            print("Warning! Continuity violated for surfaces:", self.errors)
        


    def checkContinuity(self, verbose=True):
        self.errors = []
        for surfid, surface in self.surfaces.items():
            if abs(surface.totflux - 1) > 1e-6:
                if verbose:
                    surface.printReport()
                self.errors.append(surfid)
        return len(self.errors)==0


    def plotGeometry(self, ax=None, labels=False, testsurf=[], 
        showCircle=False, markers=None, connectionLineWidth=0.5,**kwargs):
        from matplotlib.pyplot import subplots, Figure, Axes, ioff
        import matplotlib.pyplot as plt

        if isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif ax is None:
            f, ax = subplots(figsize=(5,10))
        elif not isinstance(ax, Axes):
            raise Exception('ax not a valid Figure or Axes object')
        
        if isinstance(testsurf, int):
            testsurf = [testsurf]
        
        ioff()
            
        for _, surface in self.surfaces.items():
            color = 'k'
            if (surface.ID <  self.P):
                color = 'r'
            surface.plotSelf(color=color, ax=ax, label=labels, showCircle=showCircle)

        for itest in testsurf:
            self.surfaces[itest].plotConnections(
                            ax=ax, 
                            linewidth=connectionLineWidth,
                            **kwargs
            )
    
        for line in ax.lines:
            line.set_marker(markers)

        plt.show(block=False)

        ax.set_aspect('equal')
        ax.grid(False)
        return ax.get_figure()
   
class Surface:

    def __init__(self, start, end, ID, material=1, emitting=0, absorbing=0, r_offset=1): # creates the surface
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.start (Point), self.end (Point), self.ID, self.segment (LineString), self.surfaceLength, self.midpoint (Point), self.normalStart (Point), 
        self.normalEnd (Point), self.normal (LineString), self.material, self.emitting, self.absorbing
        """

        """Start and end passed into the constructor are tuples (x, y)"""

        # # # Start and end points of the surface and a segment representation of the surface # # # 
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])
        # LineString representing the surface
        self.segment = LineString([start, end]) 
        self.ID = ID 

        self.surfaceLength = math.sqrt(
                                (self.end.x - self.start.x)**2 \
                                + (self.end.y - self.start.y)**2
                            )

        # # # Midpoint of surface # # #
        self.midpoint = self.segment.centroid

        # # # Construct normal vector # # #
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        
        # ending x coordinate for normal (don't use for reference, only to 
        # determine which direction the normal vector should point in the 
        # normal helper method)
        self.normalEndX = self.midpoint.x 
        # ending y coordinate for normal
        self.normalEndY = self.midpoint.y 

        # # # In order to get the correct slope for the normal vector # # #
        self.dx = self.end.x - self.start.x
        self.dy = self.end.y - self.start.y

        self.normalHelper(
                abs(self.dx), 
                abs(self.dy),
                self.normalEndX, 
                self.normalEndY, 
                True
        )

        # # # Start and end points, and the normal line itself (use these for reference) # # #
        self.normalStart = Point(self.midpoint.x, self.midpoint.y)
        self.normalEnd = Point(self.normalEndX, self.normalEndY)
        self.normal = LineString([self.normalStart, self.normalEnd])

        # # # Properties of the surface # # #
        self.circle = None 
        self.material = material
        self.emitting = emitting
        self.absorbing = absorbing

        """ Coupling to surfaces with LOS """
        self.neighbors = {}
        self.totflux = 0
        self.issues = {}

        self.distributionCircle(r_offset)
        return

    def distributionCircle(self, r_offset): # creates the distribution circle
        from shapely import Point, LineString, plotting, LinearRing, Polygon
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.dCircleCenter (Point), self.circle (Polygon), self.r_offset
        """

        # # # Finding radius and center of circle # # #
        self.r_offset = r_offset 
        radius = self.surfaceLength / 8
        #if self.r_offset == 1:
            #r_offset = radius
        self.dCircleCenter = self.normal.interpolate(self.r_offset * radius) # the center of the distribution circle along the normal line

        # # # for labeling the plots # # #
        if r_offset == 0:
            self.distType = "Uniform Distribution"
        elif r_offset == 1:
            self.distType = "Cosine Distribution"
        else:
            self.distType = "Distribution"
        

        numPoints = 500
        circlePoints = []

        # # # Getting the Points of the circle # # #
        ''' at pi, none of the x values for the circle will have been repeated,
            so we can add the point and then sort based on x values, descending
            order for the surfaces that have the normal going from the bottom 
            to top: if endX > startX, angle is negative '''
        # builds circle clockwise rather than counter clockwise
        if self.end.x < self.start.x: 
            for i in range(numPoints, -1, -1):
                # calculates every angle from 0-2pi, placing 500 points to 
                # create the circle
                angle = (2 * math.pi) * (i / numPoints) 
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))
        else: 
            for i in range(numPoints + 1):
                angle = (2 * math.pi) * (i / numPoints)
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))

        # # # Object representation of the distribution circle # # #
        self.circle = Polygon(circlePoints) 

        return
    
    def printReport(self):
        print("Surface {}: {}".format(self.ID, self.totflux))
        for neighid, neighbor in self.neighbors.items():
            print("    -> {}: {}".format(f"{neighid}".rjust(4), neighbor['flux']))

    def intersectionArea(self, s2): # finds the overlapping area (flux) between two surfaces, given that self has a distribution circle generated
        from shapely import Point, LineString, plotting, Polygon, is_closed
        import math

        """Reference Variables/Important:
        self.triangle, self.leg1, self.leg2, self.overlapShape, overlapArea, fractionalArea
        """

        if self.circle == None:
            print("Call distributionCircle on self before finding intersection area!")
            return

        # # # Creates/draws the relevant shapes for finding the flux (fractional area) and other reference # # #
        triangle = Polygon([s2.start, s2.end, self.midpoint, s2.start]) # triangle from midpoint of self to the endpoints of the other surface, s2
        self.leg1 = LineString([self.midpoint, s2.start]) # legs of the triangle
        self.leg2 = LineString([self.midpoint, s2.end])
        self.overlapShape = triangle.intersection(self.circle)

        self.vLeg1 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2.start.x, s2.start.y)) #vector representation

        self.vLeg2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2.end.x, s2.end.y)) #vector representation
        
        
        # # # Getting the correct area of the distribution circle on one side of the normal line # # #
        overlapArea = self.overlapShape.area
        if (0 <= self.r_offset and self.r_offset < 1) and self.circle.intersects(self.segment): # if we're dealing with a circle shifted between uniform and cosine
            buffLine = self.segment.buffer(0.000000000001)
            splitdCircle = self.circle.difference(buffLine) # split the distribution circle and the surface line
            if splitdCircle.geoms[0].intersects(self.normal): 
                self.totalAreaCircle = splitdCircle.geoms[0]
            else:
                self.totalAreaCircle = splitdCircle.geoms[1]
            circleArea = self.totalAreaCircle.area 
        else:
            circleArea = self.circle.area

        # # # Calculate the flux (fractional area) # # #
        fractionalArea = overlapArea / circleArea

        return fractionalArea, triangle 


    def getSmallestIntersectAngle(self, neighbor, geometry, polygon):
        from shapely import intersects, crosses, contains, buffer, intersection
        # # # If both legs of the triangle are blocked by the same outside region # # #
        if crosses(geometry, self.leg1) and crosses(geometry, self.leg2): 
            self.issues[neighbor.ID] = "Intersects both legs of triangle (polygon)"
            return False

        # # # If Leg 1 (from S1 midpoint to S2 start) intersects the outside region # # #
        if (    intersects(polygon, buffer(self.leg1, 0.000001)) \
                and (contains(geometry, self.leg1) == False)
        ):
            coordinates1 = list(polygon.exterior.coords)
            leg1Start = buffer(self.midpoint, 0.00001)
            leg1End = buffer(neighbor.start, 0.00001)

            if (    (contains(leg1Start, intersection(polygon, buffer(self.leg1, 0.000001))) == False) \
                    and (contains(leg1End, intersection(polygon, buffer(self.leg1, 0.000001))) == False)
            ):
                # print(f"Surface {surface1.ID} to Surface {surface2.ID} leg 1 interrupted by polygon.")
                self.issues[neighbor.ID] = "Leg 1 interrupted by polygon"
                for pair in coordinates1:
                    self.vNewLeg1 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y),
                                pair
                    )
                    self.vNewLeg2 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y), 
                                (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1])
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg1SmallestPoint = pair
            
        # # # If Leg 2 (S1 midpoint to S2 end) intersects the outside region # # #
        if (    intersects(polygon, buffer(self.leg2, 0.000001)) \
                and (contains(geometry, self.leg2) == False)
        ):
            coordinates2 = list(polygon.exterior.coords)
            leg2Start = buffer(self.midpoint, 0.00001)
            leg2End = buffer(neighbor.end, 0.00001)

            if (    (contains(leg2Start, intersection(polygon, buffer(self.leg2, 0.000001))) == False) \
                    and (contains(leg2End, intersection(polygon, buffer(self.leg2, 0.000001))) == False)
            ):
                # print(f"Surface {surface1.ID} to Surface {surface2.ID} leg 2 interrupted by polygon.")
                self.issues[neighbor.ID] = "Leg 2 interrupted by polygon"
                for pair in coordinates2:
                    self.vNewLeg1 = self.vectorHelper(
                            (self.midpoint.x, self.midpoint.y), 
                            (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1])
                    )
                    self.vNewLeg2 = self.vectorHelper(
                            (self.midpoint.x, self.midpoint.y), 
                            pair
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)
                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg2SmallestPoint = pair

        return True

    def getSmallestIntersectAngleMulti(self, neighbor, geometry, polygon):
        from shapely import intersects, crosses, contains, buffer, intersection

        # # # Both legs of the triangle are blocked by the same outer region # # #
        if (    crosses(polygon, buffer(self.leg1, 0.000001)) \
                and crosses(polygon, buffer(self.leg2, 0.000001))
        ):
            # print(f"Both legs blocked by multipolygon, Surface {surface1.ID} to Surface {surface2.ID}.")
            self.issues[ineighbor.ID] = "Both legs blocked by multipolygon"
            return False

        # # # Leg 1 (S1 midpoint to S2 start) is blocked by this polygon in the outer region # # #
        if (    intersects(polygon, buffer(self.leg1, 0.000001)) \
                and (contains(geometry, self.leg1) == False)
        ): # and crosses(geometry, surface1.leg1) 
            coordinates1 = list(polygon.exterior.coords)
            leg1Start = buffer(self.midpoint, 0.00001)
            leg1End = buffer(neighbor.start, 0.00001)
            
            if (    (contains(leg1Start, intersection(polygon, buffer(self.leg1, 0.000001))) == False) \
                    and (contains(leg1End, intersection(polygon, buffer(self.leg1, 0.000001))) == False)
            ):
                # print(f"Surface {surface1.ID} to Surface {surface2.ID} leg 1 interrupted by multipolygon.")
                self.issues[neighbor.ID] = "Leg 1 blocked by multipolygon"
                for pair in coordinates1:
                    self.vNewLeg1 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y), 
                        pair
                    )
                    self.vNewLeg2 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y), 
                        (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1])
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)
                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg1SmallestPoint = pair

        # # # Leg 2 (S2 midpoint to S2 end) is blocked by this polygon in the outer region # # #
        if (    intersects(polygon, buffer(self.leg2, 0.000001)) \
                and (contains(geometry, self.leg2) == False)
        ):  #and crosses(geometry, surface1.leg2) 
            coordinates2 = list(polygon.exterior.coords)
            self.leg2Start = buffer(self.midpoint, 0.00001)
            self.leg2End = buffer(neighbor.end, 0.00001)
            
            if (    (contains(self.leg2Start, intersection(polygon, buffer(self.leg2, 0.000001))) == False) \
                    and (contains(self.leg2End, intersection(polygon, buffer(self.leg2, 0.000001))) == False)
            ):
                # print(f"Surface {surface1.ID} to Surface {surface2.ID} leg 2 interrupted by multipolygon.")
                self.issues[neighbor.ID] = "Leg 2 blocked by multipolygon"
                for pair in coordinates2:
                    self.vNewLeg1 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y), 
                                (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1])
                    )
                    self.vNewLeg2 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y), 
                                pair
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)
                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg2SmallestPoint = pair

        return True


    def getNeighbors(self, surfaces, geometry):
        from shapely import intersects, difference
 
        # # # Fractional Area Calculations # # # 
        for neighid, neighbor in surfaces.items():
            if self.ID != neighid:
                flux, triangle = self.intersectionArea(neighbor)

                # Flux would go to the wrong (back) side of the surface
                try:
                    if not intersects(triangle, neighbor.normal):
                        self.issues[neighid] = "Wrong side of normal"
                        continue
                except:
                    continue
               
                if not intersects(triangle, self.normal):
                    self.issues[neighid] = "Wrong side of self surface"
                    continue

                # polygon/multipolygon of the intersections of the outside area of the geometry w/ the legs
                AllIntersections = difference(triangle, geometry) 

                self.leg1SmallestPoint = (neighbor.start.x, neighbor.start.y)
                self.leg2SmallestPoint = (neighbor.end.x, neighbor.end.y)

                self.vNewLeg1 = self.vLeg1
                self.vNewLeg2 = self.vLeg2
                self.smallestAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                # # # Only one area of intersection with the triangle # # #
                if AllIntersections.geom_type == 'Polygon':
                    if not self.getSmallestIntersectAngle(neighbor, geometry, AllIntersections):
                        continue
                # # # Multiple intersections with the outside region and the triangle/triangle legs # # #
                elif AllIntersections.geom_type == 'MultiPolygon':
                    # # # Check all possible polygons # # #
                    for polygon in AllIntersections.geoms:
                        if not self.getSmallestIntersectAngleMulti(neighbor, geometry, polygon):
                            continue

                # # # Set new S2 Surface and calculate the new flux. Update the total flux out of Surface 1. # # #
                newS2 = Surface(
                                (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1]), 
                                (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1]), 
                                neighbor.ID
                )
                flux, triangle = self.intersectionArea(newS2)
                self.totflux += flux
                if flux > 0:
                    if neighid not in self.neighbors:
                        self.neighbors[neighid] = {}
                    self.neighbors[neighid]['flux'] = flux
                    self.neighbors[neighid]['los'] = triangle # Get the new triangle shape and add it here



    def drawOuterCircle(self): # creates the distribution circle that will be compared to the analytic cosine (does not show plot)
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        import math
        
        """Reference Variables/Important:
            self.outerCircle (Polygon)
        """

        # # # Creating the outer circle of S2 surfaces # # #
        outerRadius = self.surfaceLength / 2 # outer circle has diameter equal to the length of self surface (S1)
        outerCenter = self.midpoint 
        outerCirclePoints = []
        numPoints = 500

        for i in range(numPoints):
                angle = (2 * math.pi) * (i / numPoints)
                x = outerCenter.x + outerRadius * math.cos(angle) # uses same center as the distribution circle
                y = outerCenter.y + outerRadius * math.sin(angle)
                outerCirclePoints.append((x, y))
    
        self.fullCircle = Polygon(outerCirclePoints) # polygon object of the full outer circle, before splitting to correct side of normal

        buffLine = self.segment.buffer(0.000000000001)
        splitCircles = self.fullCircle.difference(buffLine)
        if splitCircles.geoms[0].intersects(self.normal):
            self.outerCircle = splitCircles.geoms[0] # desired circle
        else:
            self.outerCircle = splitCircles.geoms[1] # desired circle

        return


    # # # # # # # # # # # 
    # PLOTTING FUNCTIONS #
    # # # # # # # # # # # 

    def plotSelf(self, ax=None, color='k', label=False, showCircle=True): # plots just the self surface
        from matplotlib.pyplot import subplots, ioff, Figure, Axes
        from shapely import plotting

        ioff()
        if ax is None:
            fig, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif not isinstance(ax, Axes):
            raise Exception('ax not a valid Figure or Axes object')
        ax.set_aspect('equal')

        plotting.plot_line(self.segment, ax, color=color, linewidth=2) # plots surface
        if label:
            ax.text(self.start.x, self.start.y, f"Surface {self.ID}", color='black')

        if showCircle:
            plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
        return ax

    def plotConnections(self,  ax=None, linewidth=2, 
            colorseed=1, **kwargs):
        from matplotlib.pyplot import subplots, ioff, Figure, Axes, get_cmap
        from shapely import plotting
        import random
        from numpy import linspace

        ioff()
        if ax is None:
            fig, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif not isinstance(ax, Axes):
            raise Exception('ax not a valid Figure or Axes object')
        ax.set_aspect('equal')

        cmap=get_cmap('jet')
        cols = linspace(0,1,len(self.neighbors))
        random.Random(colorseed).shuffle(cols)
        colors = iter(cmap(cols))
        for neighid, neighbor in self.neighbors.items():
            plotting.plot_polygon(
                    neighbor['los'], 
                    add_points=False,
                    color=next(colors), 
                    linewidth=linewidth,
                    **kwargs
            )
    

    def showAnalyticPlot(self, comparison, r_offset=1): # Plots the curve from the outer circle and compares it to the analytic equation plot
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        """The 'comparison' variable determines if we are comparing the plotted distribution to the known cosine distribution or not. 
        True means yes, do the comparison and plot both the generated and plot for comparison.
        False means only plot the points being generated by my code, and not the cosine equation.
        Only set 'comparison' equal to true if you are modeling a COSINE distribution (offset = 1) 
        and you want to compare it to the standard."""

        # # # Plotting setup # # #
        ioff()
        fig, ax = subplots()
        ax.set_aspect('auto')
        plt.xlabel("Angle (Radians)")
        plt.ylabel("Fractional Area")

        self.distributionCircle(r_offset)
        self.drawOuterCircle()
        
        # # # CREATING THE PLOT OF ANGLE VS AREA # # #
        s2Points = get_coordinates(self.outerCircle) # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0] # The starting point of the first S2 surface
        plotPoints = [] # Points to plot for the outer circle plot
        pdfArea = 0 # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):
            # # # End point of S2 # # # 
            s2End = s2Points[i]

            # # # S2 # # #
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]), i)

            # # # Removes some outlier points # # #
            if (s2Surface.segment.length >= (self.surfaceLength - 1)) and (s2Surface.segment.length <= (self.surfaceLength + 1)): 
                s2Start = s2End
                continue
            
            # # # Edge case outliers (mainly for uniform distribution) # # #
            startTuple = (s2Start[0], s2Start[1])
            endTuple = (s2End[0], s2End[1])
            if startTuple not in list(self.fullCircle.exterior.coords) or endTuple not in list(self.fullCircle.exterior.coords):
                #print("Removal", startTuple, endTuple)
                s2Start = s2End
                continue

            # # # Vector representations of the normal and the line from the midpoint of S2 to the midpoint of S1 # # #
            vNormal = self.vectorHelper((self.normalStart.x, self.normalStart.y), (self.normalEnd.x, self.normalEnd.y))

            vS2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2Surface.midpoint.x, s2Surface.midpoint.y))

            # # # Call helper function that uses dot product to calculate angle between vectors # # #
            angle = self.dotProductAngle(vNormal, vS2) # Plot on x-axis

            areaValue, _ = self.intersectionArea(s2Surface) # Plot on y-axis

            # # # If making a comparison to the formula plot, need to normalize the areaValue values with pdfArea and dTheta # # #
            # # # This section updates the pdfArea, not the areaValue yet # # #
            # # # Vector representations of the borders of the segments being swept out by S2 # # #
            #self.vLeg1 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2Start[0], s2Start[1]))

            #self.vLeg2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2End[0], s2End[1]))

            dTheta = self.dotProductAngle(self.vLeg1, self.vLeg2)
            pdfArea += areaValue * dTheta

            #if areaValue < 0.003: # single out the points for the outliers
                #print(s2Start, s2End)

            # # # Add the point and continue to the next iteration of the loop # # #
            plotPoints.append(Point(angle, areaValue))
            s2Start = s2End
            # # # End of for loop # # #

        # # # Just plotting the outer circle plot points-- no need to normalize # # #
        if comparison == False:
            for point in plotPoints:
                plot_points(point, ax, color='black')
             
        # # # Plotting the analytic cosine distribution from the equation y = (1/2pi) * (1 + cos(x)), scaled to be from -pi/2 to pi/2 # # #
        # # # Normalizing the outer circle area value points as well using pdfArea # # #
        else: # if comparison == True
            cosPoints = []
            adjustedPlotPoints = []
            for plotPoint in plotPoints:
                # # # Red points (not normalized) # # #
                cosXval = plotPoint.x
                cosYval = (1/math.pi)*(1 + math.cos(cosXval*2)) #(normalized from solving integral to be = to 1 -- see notebook)

                #cosYval = (1/(2*math.pi))*(1 + math.cos(cosXval*2)) (integrated this from -pi/2 to pi/2, C = 1/2)
                cosPoints.append(Point(cosXval, cosYval))

                # # # Normalized outer circle points # # #
                adjustedArea = plotPoint.y / abs(pdfArea)
                adjustedPlotPoints.append(Point(plotPoint.x, adjustedArea))
            
            for adjPoint in adjustedPlotPoints:
                plot_points(adjPoint, ax, color='black')

            for cosPoint in cosPoints:
                plot_points(cosPoint, ax, color='red', marker='1')
        

        # # # Sanity check of total fractional Area and the PDF Area # # #
        total = 0
        for point in plotPoints:
            total += point.y 

        print("Total Fractional Area: ", total)
        print("PDF Area: ", pdfArea)

        # # # Generate the plot # # #
        if self.r_offset == 0:
            yUpperLim = plotPoints[-1].y + plotPoints[-1].y * 0.1
            ax.set_ylim(bottom=0, top=yUpperLim)
        plt.show(block=False)

        return

    def showTwoSurfacePlot(self, s2, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')

        self.distributionCircle(r_offset)
        self.intersectionArea(s2)

        plotting.plot_line(self.segment, ax, color='red', linewidth=2) # plots surface
        ax.text(self.start.x, self.start.y, "Surface 1", color='red')

        plotting.plot_line(self.normal, ax, color='blue', linewidth=2) # plots normal line to surface
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color='blue')

        plotting.plot_line(self.dCircleCenter, ax, color='green') # plots point at center of distribution circle
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
        ax.text(self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color='green')

        plotting.plot_polygon(self.triangle, ax, color='orange', linewidth=2) # plots the triangle connecting to another surface
        ax.text(self.triangle.centroid.x, self.triangle.centroid.y, "Surface 2", color='orange')

        plotting.plot_polygon(self.overlapShape, ax, add_points=False, color='black', linewidth=2) # displays the overlapping area
        ax.text(self.overlapShape.centroid.x, self.overlapShape.centroid.y, "Overlap Area", color='black')

        plt.show(block=False)

        return

    def showOuterCirclePlot(self, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')

        self.distributionCircle(r_offset)
        self.drawOuterCircle()

        plotting.plot_line(self.segment, ax, color='red', linewidth=2) # plots surface
        ax.text(self.start.x, self.start.y, "Surface 1", color='red')

        plotting.plot_line(self.normal, ax, color='blue', linewidth=2) # plots normal line to surface
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color='blue')

        plotting.plot_line(self.dCircleCenter, ax, color='green') # plots point at center of distribution circle
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
        ax.text(self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color='green')

        plotting.plot_polygon(self.outerCircle, ax, add_points=True, color='gray') # displays the "outerCircle"

        plt.show(block=False)

        return



    # # # # # # # # # # # #
    # # HELPER FUNCTIONS # # 
    # # # # # # # # # # # #
    def normalHelper(self, dx, dy, endX, endY, init): # helper function to find normal direction
        if self.end.y > self.start.y:
            endX += dy
            if self.end.x > self.start.x:
                endY -= dx # (+x, -y)
            else:
                endY += dx # (+x, +y)
        else: # -x
            endX -= dy
            if self.end.x < self.start.x:
                endY += dx # (-x, +y)
            else:
                endY -= dx # (-x, -y)
        if init:
            self.normalEndX = endX
            self.normalEndY = endY
        else:
            return (endX, endY)

    def vectorHelper(self, start, end): # Finds the vector representation of a segment (surface object)
        from shapely import Point, LineString
        import numpy as np

        """Takes in tuples (x, y) that represent the start and end points of a surface."""

        iSurface = end[0] - start[0]
        jSurface = end[1] - start[1]
        vSurface = np.array([iSurface, jSurface])

        return vSurface
    
    def dotProductAngle(self, v1, v2): # Calculates the angle between two vectors using the dot product
        import numpy as np

        """v1 and v2 must be np.array objects that represent the surfaces: [iSurface, jSurface]. Use vectorHelper on the start and end points
        of a surface before passing anything into dotProductAngle."""

        dotProduct = np.dot(v1, v2)
        magnitude1 = np.linalg.norm(v1)
        magnitude2 = np.linalg.norm(v2)

        cosineAngle = dotProduct / (magnitude1 * magnitude2 + 1e-10)

        angle = np.arccos(cosineAngle) 

        crossProduct = v1[0]*v2[1] - v1[1]*v2[0]
        if crossProduct < 0:
            angle = -angle

        return angle
 


        








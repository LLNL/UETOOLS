class VacuumRegion:
    def __init__(self, case):
        return 

    def testFunction(self):
        s = Surface((5, 1), (5, 8))
        s.distributionCircle(1)
        s.distributionPlot()
        #s.fullPlot()

class Surface:

    def __init__(self, start, end): # creates the surface
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.start, self.end, self.segment, self.midpoint, self.normalStart, self.normalEnd
        """

        # Start and end should be shapely point objects?
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])
        self.segment = LineString([start, end]) # segment representing the surface

        self.surfaceLength = math.sqrt((self.end.x - self.start.x)**2 + (self.end.y - self.start.y)**2)

        # Midpoint of surface
        self.midpoint = self.segment.centroid

        # Construct normal vector
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        
        self.normalEndX = self.midpoint.x # ending x coordinate for normal (don't use for reference, only to determine which direction the normal vector should point in the normal helper method)
        self.normalEndY = self.midpoint.y # ending y coordinate for normal

        # in order to get the correct slope for the normal vector
        dx = abs(self.end.x - self.start.x)
        dy = abs(self.end.y - self.start.y)

        self.normalHelper(dx, dy)

        self.normalStart = Point(self.midpoint.x, self.midpoint.y) # Start point of normal # USE THIS FOR REFERENCE
        self.normalEnd = Point(self.normalEndX, self.normalEndY) # End point of normal USE THIS FOR REFERENCE
        self.normal = LineString([self.normalStart, self.normalEnd])

        return

    def normalHelper(self, dx, dy): # helper function to find normal direction

        if self.end.y > self.start.y:
            self.normalEndX += dy
            if self.end.x > self.start.x:
                self.normalEndY -= dx # (+x, -y)
            else:
                self.normalEndY += dx # (+x, +y)
        else: # -x
            self.normalEndX -= dy
            if self.end.x < self.start.x:
                self.normalEndY += dx # (-x, +y)
            else:
                self.normalEndY -= dx # (-x, -y)


    def distributionCircle(self, r_offset): # draws the distribution circle
        from shapely import Point, LineString, plotting, LinearRing, Polygon
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.dCircleCenter, self.circle
        """

        radius = 1.0
        self.dCircleCenter = self.normal.interpolate(r_offset) # the center of the circle along the normal line

        # for labeling the plots
        if r_offset == 0:
            self.distType = "Uniform Distribution"
        if r_offset == 1:
            self.distType = "Cosine Distribution"

        numPoints = 500
        circlePoints = []
        # at pi, none of the x values for the circle will have been repeated, so we can add the point and then sort based on x values, descending order
        # for the surfaces that have the normal going from the bottom to top:
            # if endX > startX, angle is negative

        if self.end.x < self.start.x: # builds circle clockwise rather than counter clockwise
            for i in range(numPoints, -1, -1):
                angle = (2 * math.pi) * (i / numPoints) # calculates every angle from 0-2pi, placing 500 points to create the circle
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and r_offset != 0:
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))
        else: 
            for i in range(numPoints + 1):
                angle = (2 * math.pi) * (i / numPoints)
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and r_offset != 0:
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))

        
        self.circle = Polygon(circlePoints) # object representation of the circle

        return

    def intersectionArea(self, s2): # finds the overlapping area (flux) between two surfaces
        from shapely import Point, LineString, plotting, Polygon
        import math
        """ For now: pass in a surface object, s2, so that we can access the end points of the surface and make a triangle from the midpoint of self to the endpoints of s2.
            The parameters s2Start and s2End would represent the start and end points of the second surface (the one we are measuring flux into from self).
            This parameter might be changed instead to an array that represents the start and end points of all relevant surfaces to self, rather than just 2 input points.
            """

        """Reference Variables/Important:
        self.triangle, self.overlapShape, overlapArea, fractionalArea
        """
        # find area that overlaps between circle and triangle 
        # areaOverlap / circleArea = flux on surface 2 from surface 1

        # need 3 segements (LineStrings) to make the triangle
        # s2 start to end, s2 end to s1 midpoint, s1 midpoint to s2 end

        self.triangle = Polygon([s2.start, s2.end, self.midpoint, s2.start]) # triangle from midpoint of self to the endpoints of the other surface, s2
        self.overlapShape = self.triangle.intersection(self.circle)
        
        overlapArea = self.overlapShape.area
        circleArea = self.circle.area

        fractionalArea = overlapArea / circleArea

        return fractionalArea

    def distributionPlot(self): # creates a plot of the distribution circle that will be compared to the analytic cosine
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        
        """Reference Variables/Important:
            self.outerCircle

        """

        # creating the outer circle of S2 surfaces
        outerRadius = self.surfaceLength / 2
        outerCenter = self.midpoint
        outerCirclePoints = []
        numPoints = 500

        for i in range(numPoints + 1):
                angle = (2 * math.pi) * (i / numPoints)
                x = outerCenter.x + outerRadius * math.cos(angle) # uses same center as the distribution circle
                y = outerCenter.y + outerRadius * math.sin(angle)
                outerCirclePoints.append((x, y))
    
        fullCircle = Polygon(outerCirclePoints) # polygon object of the full outer circle

        buffLine = self.segment.buffer(0.000000000001)
        splitCircles = fullCircle.difference(buffLine) # make this not self after done testing things out
        if splitCircles.geoms[0].intersects(self.normal):
            self.outerCircle = splitCircles.geoms[0]
        else:
            self.outerCircle = splitCircles.geoms[1]



        # # CREATING THE PLOT OF ANGLE VS AREA # #
        s2Points = get_coordinates(self.outerCircle) 
        s2Start = s2Points[0]
        plotPoints = []

        for i in range(1, len(s2Points), 1):
            s2End = s2Points[i]
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1])) # the surface to pass in as s2 into intersectionArea
            s2MidpointLine = LineString([self.midpoint, s2Surface.midpoint])
            # take dot product to find the angle from the normal (?)
            iNormal = self.normalEnd.x - self.normalStart.x 
            jNormal = self.normalEnd.y - self.normalStart.y
            magnitudeNormal = math.sqrt(iNormal**2 + jNormal**2)

            iS2 = s2Surface.midpoint.x - self.midpoint.x # line from midpoint of self to s2 to calculate the angle between the surface and the normal
            jS2 = s2Surface.midpoint.y - self.midpoint.y 
            magnitudeS2 = math.sqrt(iS2**2 + jS2**2)

            #angle = math.acos((iNormal * iS2) + (jNormal * jS2) / (magnitudeNormal * magnitudeS2)) # store as x in tuple
            angle = math.atan2(jS2, iS2) - math.atan2(jNormal, iNormal)
            areaValue = self.intersectionArea(s2Surface) # store as y in tuple
            plotPoints.append(Point(angle, areaValue))

            s2Start = s2End

        """plotting the actual distribution"""
        ioff()
        fig, ax = subplots()
        #ax.set_aspect('equal')         
        for point in plotPoints:
           plot_points(point, ax, color='black')
        plt.xlabel("Angle (Radians)")
        plt.ylabel("Fractional Area")
        plt.show()

        return

    def fullPlot(self): # plots all relevant objects (surface, normal, distribution circle, etc.) (excluding the plot of the distribution itself (angle vs area plot))
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')
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

        plotting.plot_polygon(self.outerCircle, ax, add_points=True, color='gray') # displays the "outerCircle" which is needed to create the analytic distributions (angle vs fractionalArea plot)

        plt.show()

        return

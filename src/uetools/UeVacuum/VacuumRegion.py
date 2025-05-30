class VacuumRegion:
    def __init__(self, case):
        return 

    def twoSurfacePlot(self):
        S1 = Surface((1, 2), (4, 5))
        S1.distributionCircle(1)
        S2 = Surface((7, 4), (6, 2))
        S1.intersectionArea(S2)
        S1.showTwoSurfacePlot()
    
    def outerCirclePlot(self):
        S1 = Surface((1, 2), (4, 5))
        S1.distributionCircle(1)
        S1.drawOuterCircle()
        S1.showOuterCirclePlot()

    def analyticPlot(self):
        S1 = Surface((1, 2), (4, 5))
        S1.distributionCircle(1)
        S1.drawOuterCircle()
        S1.showAnalyticPlot()
    
    def fullPlotTest(self):
        S1 = Surface((1, 2), (4, 5))
        S2 = Surface((7, 4), (6, 2))
        S1.distributionCircle(1)
        S1.intersectionArea(S2)
        S1.drawOuterCircle()
        S1.fullPlot()

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
        self.dCircleCenter, self.circle, self.r_offset
        """

        self.r_offset = r_offset
        radius = 1.0
        self.dCircleCenter = self.normal.interpolate(r_offset) # the center of the distribution circle along the normal line

        # for labeling the plots
        if r_offset == 0:
            self.distType = "Uniform Distribution"
        elif r_offset == 1:
            self.distType = "Cosine Distribution"
        else:
            self.distType = "Distribution"
        

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

        
        self.circle = Polygon(circlePoints) # object representation of the circle

        return

    def intersectionArea(self, s2): # finds the overlapping area (flux) between two surfaces
        # also draws the triangle between s1 and s2
        from shapely import Point, LineString, plotting, Polygon
        import math

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
        if (0 < self.r_offset and self.r_offset < 1): # if we're dealing with a circle shifted between uniform and cosine
            buffLine = self.segment.buffer(0.000000000001)
            splitdCircle = self.circle.difference(buffLine) # split the distribution circle and the surface line
            if splitdCircle.geoms[0].intersects(self.normal): 
                self.totalAreaCircle = splitdCircle.geoms[0]
            else:
                self.totalAreaCircle = splitdCircle.geoms[1]
            circleArea = self.totalAreaCircle.area 
        else:
            circleArea = self.circle.area

        fractionalArea = overlapArea / circleArea

        return fractionalArea 

    def drawOuterCircle(self): # creates a plot of the distribution circle that will be compared to the analytic cosine
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
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
        
        return


    def showAnalyticPlot(self): # plots the curve
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math

        # # CREATING THE PLOT OF ANGLE VS AREA # #
        s2Points = get_coordinates(self.outerCircle) # points defining the outer circle
        s2Start = s2Points[0]
        plotPoints = [] # points to plot for the curve

        for i in range(1, len(s2Points), 1):
            s2End = s2Points[i] # end point of the surface

            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1])) # the surface to pass in as s2 into intersectionArea

            # take dot product to find the angle from the normal
            iNormal = self.normalEnd.x - self.normalStart.x 
            jNormal = self.normalEnd.y - self.normalStart.y

            iS2 = s2Surface.midpoint.x - self.midpoint.x # line from midpoint of self to s2 to calculate the angle between the surface and the normal
            jS2 = s2Surface.midpoint.y - self.midpoint.y 

            angle = math.atan2(jS2, iS2) - math.atan2(jNormal, iNormal)
            areaValue = self.intersectionArea(s2Surface) # store as y in tuple
            plotPoints.append(Point(angle, areaValue))
    
            s2Start = s2End

        """plotting the actual distribution"""
        ioff()
        fig, ax = subplots()

        for point in plotPoints: # uncomment
           plot_points(point, ax, color='black')
        plt.xlabel("Angle (Radians)")
        plt.ylabel("Fractional Area")
        
        total = 0
        for point in plotPoints:
            total += point.y 

        print("Total Area: ", total)

        plt.show()

        return

    def showTwoSurfacePlot(self):
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

        plt.show()

        return

    def showOuterCirclePlot(self):
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

        plotting.plot_polygon(self.outerCircle, ax, add_points=True, color='gray') # displays the "outerCircle"

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

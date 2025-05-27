class VacuumRegion:
    def __init__(self, case):
        return 

    def testFunction(self):
        s = Surface((2, 1), (1, 2))
        s.distribution(1.0)
        s.fullPlot()

class Surface:

    def __init__(self, start, end):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        # Start and end should be shapely point objects?
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])
        self.segment = LineString([start, end])

        # Midpoint of surface
        self.midpoint = self.segment.centroid

        # Construct normal vector
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        
        self.normalEndX = self.midpoint.x # ending x coordinate for normal
        self.normalEndY = self.midpoint.y # ending y coordinate for normal

        # in order to get the correct slope for the normal vector
        dx = abs(self.end.x - self.start.x)
        dy = abs(self.end.y - self.start.y)

        self.normalHelper(dx, dy)

        self.normalStart = Point(self.midpoint.x, self.midpoint.y) # Start point of normal
        self.normalEnd = Point(self.normalEndX, self.normalEndY) # End point of normal
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


    def distribution(self, r_offset): # for now r_offset can only equal 0 (uniform) or 1 (cosine)
        from shapely import Point, LineString, plotting, LinearRing, Polygon
        from matplotlib.pyplot import subplots
        import math

        radius = 1.0
        self.center = self.normal.interpolate(r_offset) # the center of the circle along the normal line

        numPoints = 500
        circlePoints = []
        # at pi, none of the x values for the circle will have been repeated, so we can add the point and then sort based on x values, descending order
        # for the surfaces that have the normal going from the bottom to top:
            # if endX > startX, angle is negative

        if self.end.x < self.start.x: # builds circle clockwise rather than counter clockwise
            for i in range(500, -1, -1):
                angle = (2 * math.pi) * (i / numPoints) # calculates every angle from 0-2pi, placing 500 points to create the circle
                x = self.center.x + radius * math.cos(angle)
                y = self.center.y + radius * math.sin(angle)
                if angle == math.pi:
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))
        else: 
            for i in range(500 + 1):
                angle = (2 * math.pi) * (i / numPoints) # calculates every angle from 0-2pi, placing 500 points to create the circle
                x = self.center.x + radius * math.cos(angle)
                y = self.center.y + radius * math.sin(angle)
                if angle == math.pi:
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))

        
        self.circle = Polygon(circlePoints) # a representation of the circle
        
        #self.circle = self.center.buffer(radius, 5000)
        centerDist = math.sqrt((self.center.x - self.midpoint.x)**2 + (self.center.y - self.midpoint.y)**2) 
        print("Distance from center to midpoint: ", centerDist) # checks that the r_offset is correct for the specified distribution
        print("Center of circle: (", self.center.x, ", ", self.center.y, ")")
        return

    
    def fullPlot(self): # plots all relevant objects (surface, normal, distribution circle)
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots

        fig, ax = subplots()
        ax.set_aspect('equal')
        plotting.plot_line(self.segment, ax, color='red', linewidth=2)
        plotting.plot_line(self.normal, ax, color='blue', linewidth=2)
        plotting.plot_line(self.center, ax, color='green')
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1)
        return

    def intersectionArea(self, angles):
        from shapely import Point, LineString, plotting
        import math
        """ The parameter angles is a 1-D array of all the angles from the normal of self to the endpoints of the other 
                surfaces. 
        # find area that overlaps between circle and triangle 
        # areaOverlap / circleArea = flux on surface 2 from surface 1
        # triangle is based on the angles we are covering from the midpoint (0 to pi)

        circleArea = self.circle.area
        
        return




    
class VacuumRegion:
    def __init__(self, case):
        return 

    def testFunction(self):
        s = Surface((3, 4), (6, 5))
        s.distribution(1)
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
        self.normalStartX = self.midpoint.x # starting x coordinate for normal
        self.normalStartY = self.midpoint.y # starting y coordinate for normal
        self.normalEndX = self.midpoint.x # ending x coordinate for normal
        self.normalEndY = self.midpoint.y # ending y coordinate for normal

        # in order to get the correct slope for the normal vector
        dx = abs(self.end.x - self.start.x)
        dy = abs(self.end.y - self.start.y)

        self.normalHelper(dx, dy)

        self.normalStart = Point(self.normalStartX, self.normalStartY) # Start point of normal
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
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        radius = 1.0
        self.center = self.normal.interpolate(r_offset)
        self.circle = self.center.buffer(radius, 5000)
        centerDist = math.sqrt((self.center.x - self.midpoint.x)**2 + (self.center.y - self.midpoint.y)**2) # check that the offset is correct
        print("Distance from center to midpoint: ", centerDist)
        return

        # 1 unitize normal vector
        # 2 multiply normal points by r_offset
        # 3 add midpoint + multiplied points for location of center of circle

        ''' #1
        normFactor = math.sqrt((self.normalEndX-self.normalStartX)**2+(self.normalEndY-self.normalStartY)**2)
        unitNormX = (self.normalEndX-self.normalStartX) / normFactor
        unitNormY = (self.normalEndY-self.normalStartY) / normFactor
        #2
        unitNormX = unitNormX * r_offset
        unitNormY = unitNormY * r_offset
        #3
        centerX = self.midpoint.x + unitNormX
        centerY = self.midpoint.y + unitNormY

        self.center = Point(centerX, centerY)
        self.circle = self.center.buffer(radius, 5000)

        centerDist = math.sqrt((centerX - self.midpoint.x)**2 + (centerY - self.midpoint.y)**2) # check that the offset is correct
        print("Distance from center to midpoint: ", centerDist) '''
    
    def fullPlot(self):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots

        fig, ax = subplots()
        ax.set_aspect('equal')
        plotting.plot_line(self.segment, ax, color='red', linewidth=2)
        plotting.plot_line(self.normal, ax, color='blue', linewidth=2)
        plotting.plot_line(self.center, ax, color='green')
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', facecolor=None, linewidth=1)
        return

    def intersectionArea(self):
        from shapely import Point, LineString, plotting
        import math
        # find area that overlaps between circle and triangle 
        # areaOverlap / circleArea = flux on surface

        circleArea = self.circle.area
        minAngle = 0
        maxAngle = math.pi

        return




    
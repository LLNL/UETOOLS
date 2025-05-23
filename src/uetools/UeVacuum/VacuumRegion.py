class VacuumRegion:
    def __init__(self, case):
        return 

    def testFunction(self):
        s = Surface((3, 2), (4, 5))

class Surface:
    def __init__(self, start, end):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math
        # Start and end should be shapely point objects?
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])
        segment = LineString([start, end])

        # Midpoint of surface
        self.midpoint = segment.centroid

        # Construct normal vector
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        normalStartX = self.midpoint.x # starting x coordinate for normal
        normalStartY = self.midpoint.y # starting y coordinate for normal
        normalEndX = self.midpoint.x # ending x coordinate for normal
        normalEndY = self.midpoint.y # ending y coordinate for normal

        # in order to get the correct slope for the normal vector
        dx = abs(self.end.x - self.start.x)
        dy = abs(self.end.y - self.start.y)

        if self.end.y > self.start.y:
            normalEndX += dy
            if self.end.x > self.start.x:
                normalEndY -= dx # (+x, -y)
            else:
                normalEndY += dx # (+x, +y)
        else: # -x
            normalEndX -= dy
            if self.end.x < self.start.x:
                normalEndY += dx # (-x, +y)
            else:
                normalEndY -= dx # (-x, -y)

    
        normalStart = Point(normalStartX, normalStartY) # Start point of normal
        normalEnd = Point(normalEndX, normalEndY) # End point of normal
        normal = LineString([normalStart, normalEnd])
        
        # Plotting the segment and the normal
        fig, ax = subplots()
        plotting.plot_line(segment, ax, color='red', linewidth=2)
        plotting.plot_line(normal, ax, color='blue', linewidth=2)
    
        print(dy/dx)
        print( (normalEndX - normalStartX)/(normalEndY - normalStartY))
        return

    def distribution(r_offset):
        radius = 1
        centerX = self.midpoint.x + r_offset
        centerY = self.midpoint.x + r_offset

        center = Point(centerX, centerY)
        circle = center.buffer()

        return


    
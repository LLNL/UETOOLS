class VacuumRegion:
    def __init__(self, case):
        return 

    def testFunction(self):
        s = Surface((3, 5), (4, 6))

class Surface:
    def __init__(self, start, end):
        # Start and end should be shapely point objects?
        self.start = start
        self.end = end
        segment = LineString([start, end])

        # Midpoint of surface
        midpoint = segment.centroid

        # Construct normal vector
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        normalStartX = midpoint.x # starting x coordinate for normal
        normalStartY = midpoint.y # starting y coordinate for normal
        normalEndX = normalStartX + end.x - start.x # ending x coordinate for normal
        normalEndY = normalStartY + end.y - start.y # ending y coordinate for normal

        if end.y > start.y: # Normal has +x
            normalEndX = abs(normalEndX)
        else: # Normal has -x
            normalEndX = -abs(normalEndX)

        if end.x < start.x: # Normal has +y
            normalEndY = abs(normalEndY)
        else: # Normal has -y
            normalEndY = -abs(normalEndY)


        normalStart = Point(normalStartX, normalStartY) # Start point of normal
        normalEnd = Point(normalEndX, normalEndY) # End point of normal
        normal = LineString([normalStart, normalEnd])
        

        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        # Plotting the segment and the normal
        fig, ax = subplots()
        plotting.plot_line(segment, ax, color='red', linewidth=2)
        plotting.plot_line(normal, ax, color='blue', linewidth=2)
        return

    def distribution(r_offset):
        return


    
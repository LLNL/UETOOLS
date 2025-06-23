class VacuumRegion:
    def __init__(self, case):
        return 

    def twoSurfacePlot(self):
        S1 = Surface((2, 5), (4, 2))
        S1.distributionCircle(1)
        S2 = Surface((2, 1), (3, 1))
        S1.intersectionArea(S2)
        S1.showTwoSurfacePlot()
    
    def outerCirclePlot(self):
        S1 = Surface((2, 5), (4, 1))
        S1.distributionCircle(1)
        S1.drawOuterCircle()
        S1.showOuterCirclePlot()

    def analyticPlot(self):
        S1 = Surface((2, 5), (4, 1))
        S1.distributionCircle(1)
        S1.drawOuterCircle()
        S1.showAnalyticPlot(False)
    
    def fullPlotTest(self):
        S1 = Surface((5, 2), (4, 1))
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
        self.start, self.end, self.segment, self.surfaceLength, self.midpoint, self.normalStart, self.normalEnd
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
        radius = self.surfaceLength / 8
        #if self.r_offset == 1:
            #r_offset = radius
        self.dCircleCenter = self.normal.interpolate(self.r_offset * radius) # the center of the distribution circle along the normal line

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
        
        # Getting the correct area of the circle on one side of the normal line
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

        fractionalArea = overlapArea / circleArea

        return fractionalArea 

    def drawOuterCircle(self): # creates a plot of the distribution circle that will be compared to the analytic cosine
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        import math
        
        """Reference Variables/Important:
            self.outerCircle
        """

        # creating the outer circle of S2 surfaces
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
            self.outerCircle = splitCircles.geoms[0]
        else:
            self.outerCircle = splitCircles.geoms[1]
        
        return

    def showAnalyticPlot(self, comparison): # Plots the curve from the outer circle and compares it to the analytic equation plot
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
        

        # # # CREATING THE PLOT OF ANGLE VS AREA # # #
        s2Points = get_coordinates(self.outerCircle) # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0] # The starting point of the first S2 surface
        plotPoints = [] # Points to plot for the outer circle plot
        pdfArea = 0 # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):
            # # # End point of S2 # # # 
            s2End = s2Points[i]

            # # # S2 # # #
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]))

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
            iNormal = self.normalEnd.x - self.normalStart.x
            jNormal = self.normalEnd.y - self.normalStart.y
            vNormal = np.array([iNormal, jNormal])

            iS2 = s2Surface.midpoint.x - self.midpoint.x
            jS2 = s2Surface.midpoint.y - self.midpoint.y 
            vS2 = np.array([iS2, jS2])

            # # # Call helper function that uses dot product to calculate angle between vectors # # #
            angle = self.dotProductAngle(vNormal, vS2) # Plot on x-axis

            areaValue = self.intersectionArea(s2Surface) # Plot on y-axis

            # # # If making a comparison to the formula plot, need to normalize the areaValue values with pdfArea and dTheta # # #
            # # # This section updates the pdfArea, not the areaValue yet # # #
            # # # Vector representations of the borders of the segments being swept out by S2 # # #
            iLeg1 = s2Start[0] - self.midpoint.x 
            jLeg1 = s2Start[1] - self.midpoint.y
            vLeg1 = np.array([iLeg1, jLeg1])

            iLeg2 = s2End[0] - self.midpoint.x 
            jLeg2 = s2End[1] - self.midpoint.y 
            vLeg2 = np.array([iLeg2, jLeg2])
            
            dTheta = self.dotProductAngle(vLeg1, vLeg2)
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
        plt.show(block=False)

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

        plt.show(block=False)

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

        plt.show(block=False)

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

        plt.show(block=False)

        return

    def dotProductAngle(self, v1, v2):
        import numpy as np

        dotProduct = np.dot(v1, v2)
        magnitude1 = np.linalg.norm(v1)
        magnitude2 = np.linalg.norm(v2)

        cosineAngle = dotProduct / (magnitude1 * magnitude2)

        angle = np.arccos(cosineAngle) 

        crossProduct = v1[0]*v2[1] - v1[1]*v2[0]
        if crossProduct < 0:
            angle = -angle

        return angle


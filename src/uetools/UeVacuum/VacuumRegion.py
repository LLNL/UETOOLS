class VacuumRegion:
    def __init__(self, case):
        return 

    def twoSurfacePlot(self):
        S1 = Surface((2, 5), (4, 2))
        S1.distributionCircle(1)
        S2 = Surface((2, 1), (3, 2))
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
        S1.showAnalyticPlot(True)
    
    def trianglePlot(self):
        S1 = Surface((1, 3), (2, 6))
        S1.geometries(S1.equilateral(), 1)

    def squarePlot(self):
        S1 = Surface((1, 2), (3, 6))
        S1.geometries(S1.square(), 1)
    
    def lineOfSightPlot(self):
        source = Surface((3, 1), (1, 1))
        source.geometries(source.lineOfSightVertices(), 1)

class Surface:

    def __init__(self, start, end, material=1, emitting=0, absorbing=0): # creates the surface
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.start (Point), self.end (Point), self.segment (LineString), self.surfaceLength, self.midpoint (Point), self.normalStart (Point), 
        self.normalEnd (Point), self.normal (LineString), self.material, self.emitting, self.absorbing
        """

        """Start and end passed into the constructor are tuples (x, y)"""

        # # # Start and end points of the surface and a segment representation of the surface # # # 
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])
        self.segment = LineString([start, end]) # LineString representing the surface

        self.surfaceLength = math.sqrt((self.end.x - self.start.x)**2 + (self.end.y - self.start.y)**2)

        # # # Midpoint of surface # # #
        self.midpoint = self.segment.centroid

        # # # Construct normal vector # # #
            # slope of normal vector = -slope of line segment
            # add x and y components of the line segment to the midpoint to get the correct
                # x and y for the normal
            # need to make sure that we set the normal in the right direction:
                # if end.Y > start.Y --> positive x direction
                # if end.X < start.X --> positive y direction
        
        self.normalEndX = self.midpoint.x # ending x coordinate for normal (don't use for reference, only to determine which direction the normal vector should point in the normal helper method)
        self.normalEndY = self.midpoint.y # ending y coordinate for normal

        # # # In order to get the correct slope for the normal vector # # #
        self.dx = self.end.x - self.start.x
        self.dy = self.end.y - self.start.y

        self.normalHelper(abs(self.dx), abs(self.dy), self.normalEndX, self.normalEndY, True)

        # # # Start and end points, and the normal line itself (use these for reference) # # #
        self.normalStart = Point(self.midpoint.x, self.midpoint.y)
        self.normalEnd = Point(self.normalEndX, self.normalEndY)
        self.normal = LineString([self.normalStart, self.normalEnd])

        # # # Properties of the surface # # #
        self.circle = None 
        self.material = material
        self.emitting = emitting
        self.absorbing = absorbing

        return

    def distributionCircle(self, r_offset): # draws the distribution circle
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
        '''at pi, none of the x values for the circle will have been repeated, so we can add the point and then sort based on x values, descending order
        # for the surfaces that have the normal going from the bottom to top:
            # if endX > startX, angle is negative '''
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

        # # # Object representation of the distribution circle # # #
        self.circle = Polygon(circlePoints) 

        return

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
        self.triangle = Polygon([s2.start, s2.end, self.midpoint, s2.start]) # triangle from midpoint of self to the endpoints of the other surface, s2
        self.leg1 = LineString([self.midpoint, s2.start]) # legs of the triangle
        self.leg2 = LineString([self.midpoint, s2.end])
        self.overlapShape = self.triangle.intersection(self.circle)

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

        return fractionalArea 

    def drawOuterCircle(self): # creates the plot of the distribution circle that will be compared to the analytic cosine (does not show plot)
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
            vNormal = self.vectorHelper((self.normalStart.x, self.normalStart.y), (self.normalEnd.x, self.normalEnd.y))

            vS2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2Surface.midpoint.x, s2Surface.midpoint.y))

            # # # Call helper function that uses dot product to calculate angle between vectors # # #
            angle = self.dotProductAngle(vNormal, vS2) # Plot on x-axis

            areaValue = self.intersectionArea(s2Surface) # Plot on y-axis

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

    def geometries(self, nodeList, distributionType): # Does flux calculations for a geometry, such as a triangle or a square or something more complex
        from shapely import Point, LineString, plotting, Polygon, intersects, contains, intersection, is_closed, LinearRing, buffer
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        """The parameter nodeList WILL include the start and end points of S1 (or self) to define the endpoints of the first and last surface.
        Ex: For a triangle with base (0, 0) to (1, 1), nodeList would be [(0, 0), (2, 2), (1, 1)]. 
        DistributionType determines if we are using cosine or uniform or other type of distribution circle with a number like
        0 (uniform) or 1 (cosine)."""

        """nodeList should be an array of tuples like this: [(x1, y1), (x2, y2), ..., (xn, yn)]"""

        # # # Plot set-up # # #
        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')

        if len(nodeList) == 3:
            geometryType = "Triangle"
        elif len(nodeList) == 4:
            geometryType = "Square"
        else:
            geometryType = "Other Geometry"
    
        # # # Set up surfaces of geometry and the polygon object # # #
        surfaces = [self]
        startNode = Point(nodeList[0])
        for i in range(1, len(nodeList), 1):
            endNode = Point(nodeList[i])
            s = Surface((startNode.x, startNode.y), (endNode.x, endNode.y))
            s.distributionCircle(distributionType)
            surfaces.append(s)
            startNode = endNode
            # later if n < i < m: s.reflecting = 0, something like this
        self.distributionCircle(distributionType)

        nodeList.append(nodeList[0]) # to "close" the geometry
        geometry = Polygon(nodeList) # boundary of the figure

        # # # Plot geometry # # #
        for surface in surfaces:
            if surface == self:
                plotting.plot_line(surface.segment, ax, color='blue', linewidth=2)
            else:
                plotting.plot_line(surface.segment, ax, color='black', linewidth=2) # plots surfaces
            plotting.plot_polygon(surface.circle, ax, add_points=False, color='green', linewidth=1) # plots distribution circles
            
        plotting.plot_line(self.normal, ax, color='gray', linewidth=2) # plots normal line
        ax.text(self.midpoint.x, self.midpoint.y, "Surface 1 (Source)", color='blue') # labels the self/S1 surface

        # # # Fractional Area Calculations # # # --> add self to surfaces and have double for loops running
        i = 1 # labeling 
        s1test = 9
        s2test = 4
        for surface1 in surfaces:
            j = 1 # labeling

            for surface2 in surfaces:
                if surface1 == surface2: # don't want self to self flux
                    j += 1
                    continue

                flux = surface1.intersectionArea(surface2)

                if intersects(surface1.triangle, surface2.normal) == False: # flux would go to the wrong (back) side of the surface
                    flux = 0
                    j += 1
                    continue

                fullIntersection = surface1.triangle.difference(geometry) # polygon/multipolygon of the intersections w/ the legs
                leg1SmallestPoint = (surface2.start.x, surface2.start.y)
                leg2SmallestPoint = (surface2.end.x, surface2.end.y)

                vNewLeg1 = surface1.vLeg1
                vNewLeg2 = surface1.vLeg2
                smallestAngle = self.dotProductAngle(vNewLeg1, vNewLeg2)
                
                if i == s1test and j == s2test: # testing
                    print(f"polygon coords from {i} to {j}: {fullIntersection}")

                if fullIntersection.geom_type == 'Polygon': # only one polygon intersects the legs of the polygon
                    if i == s1test and j == s2test: # testing
                        print("polygon type, first if statement")
            
                    if intersects(fullIntersection, buffer(surface1.leg1, 0.01)) and intersects(fullIntersection, buffer(surface1.leg2, 0.01)): # intersects both legs of triangle
                        flux = 0
                        j += 1
                        continue

                    elif intersects(fullIntersection, buffer(surface1.leg1, 0.01)): # single polygon intersects leg1 of triangle
                        coordinates1 = list(fullIntersection.exterior.coords)
                        if i == s1test and j == s2test:
                            print("Coordinates 1: ", coordinates1)

                        for pair in coordinates1:
                            vNewLeg1 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), pair)
                            vNewLeg2 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), (leg2SmallestPoint[0], leg2SmallestPoint[1]))
                            newAngle = self.dotProductAngle(vNewLeg1, vNewLeg2)
                            if abs(newAngle) < abs(smallestAngle):
                                smallestAngle = newAngle
                                leg1SmallestPoint = pair
                            if i == s1test and j == s2test:
                                print("points:", leg1SmallestPoint, "and", leg2SmallestPoint)
                                print("current point:", pair)
                                print("current angle:", newAngle)
                                print("smallest angle polygon leg1:", smallestAngle)
                        
                    elif intersects(fullIntersection, buffer(surface1.leg2, 0.01)): # single polygon intersects leg2 of triangle
                        coordinates2 = list(fullIntersection.exterior.coords)

                        for pair in coordinates2:
                            vNewLeg1 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), (leg1SmallestPoint[0], leg1SmallestPoint[1]))
                            vNewLeg2 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), pair)
                            newAngle = self.dotProductAngle(vNewLeg1, vNewLeg2)
                            if abs(newAngle) < abs(smallestAngle):
                                smallestAngle = newAngle
                                leg2SmallestPoint = pair
                            if i == s1test and j == s2test:
                                print("points:", leg1SmallestPoint, "and", leg2SmallestPoint)
                                print("current point:", pair)
                                print("current angle:", newAngle)
                                print("smallest angle polygon leg2:", smallestAngle)

                elif fullIntersection.geom_type == 'MultiPolygon': # multiple polygons that intersect the legs
                    if i == s1test and j == s2test: # testing
                        print("multipolygon type, second if statement")
                        k = 1 # remove after testing
                    for polygon in fullIntersection.geoms:
                        if i == s1test and j == s2test:
                            print("multipolygon number:", k)
                        #same logic as above (repetitive for now)
                        if intersects(polygon, buffer(surface1.leg1, 0.01)) and intersects(polygon, buffer(surface1.leg2, 0.01)): # intersects both legs of triangle
                            flux = 0
                            j += 1
                            continue

                        elif intersects(polygon, buffer(surface1.leg1, 0.01)): # single polygon intersects leg1 of triangle
                            coordinates1 = list(polygon.exterior.coords)

                            for pair in coordinates1:
                                vNewLeg1 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), pair)
                                vNewLeg2 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), (leg2SmallestPoint[0], leg2SmallestPoint[1]))
                                newAngle = self.dotProductAngle(vNewLeg1, vNewLeg2)
                                if abs(newAngle) < abs(smallestAngle):
                                    smallestAngle = newAngle
                                    leg1SmallestPoint = pair
                                if i == s1test and j == s2test: # testing
                                    print("points:", leg1SmallestPoint, "and", leg2SmallestPoint)
                                    print("current point:", pair)
                                    print("current angle:", newAngle)
                                    print("smallest angle mp leg1:", smallestAngle)
                            if i == s1test and j == s2test: # testing
                                k += 1

                        elif intersects(polygon, buffer(surface1.leg2, 0.01)): # single polygon intersects leg2 of triangle
                            coordinates2 = list(polygon.exterior.coords)
                            for pair in coordinates2:
                                vNewLeg1 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), (leg1SmallestPoint[0], leg1SmallestPoint[1]))
                                vNewLeg2 = self.vectorHelper((surface1.midpoint.x, surface1.midpoint.y), pair)
                                newAngle = self.dotProductAngle(vNewLeg1, vNewLeg2)
                                if abs(newAngle) < abs(smallestAngle):
                                    smallestAngle = newAngle
                                    leg2SmallestPoint = pair
                                if i == s1test and j == s2test:
                                    print("points:", leg1SmallestPoint, "and", leg2SmallestPoint)
                                    print("current point:", pair)
                                    print("current angle:", newAngle)
                                    print("smallest angle mp leg2:", smallestAngle)
                            if i == s1test and j == s2test:
                                k += 1

                if i == s1test and j == s2test: # find intersection points
                    #plotting.plot_line(newS2.segment, ax, add_points=True, color='orange', linewidth=1)
                    #plotting.plot_line(surface1.leg1, ax, add_points=True, color='orange', linewidth=1)
                    #plotting.plot_line(surface1.leg2, ax, add_points=True, color='orange', linewidth=1)

                    newS2 = Surface((leg1SmallestPoint[0], leg1SmallestPoint[1]), (leg2SmallestPoint[0], leg2SmallestPoint[1]))
                    flux = surface1.intersectionArea(newS2)
                    plotting.plot_polygon(surface1.triangle, ax, add_points=True, color='orange', linewidth=1)
                    intersectionPoints = intersection(geometry, surface1.triangle)
                    plotting.plot_points(intersectionPoints, ax, color='red', marker='1')
                    print(f"({geometryType}) Surface {i} onto Surface {j}: {flux}")
                    print(newS2.start, newS2.end)

                j += 1

            ax.text(surface1.midpoint.x, surface1.midpoint.y, f"Surface {i}", color='black')
            i += 1

        # # # Generates the Plot # # #
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

        cosineAngle = dotProduct / (magnitude1 * magnitude2)

        angle = np.arccos(cosineAngle) 

        crossProduct = v1[0]*v2[1] - v1[1]*v2[0]
        if crossProduct < 0:
            angle = -angle

        return angle
    
    def equilateral(self): # Generates points that form an equilateral triangle with self
        from shapely import Point, LineString
        import math
        import numpy as np

        height = math.sqrt(3) * self.surfaceLength / 2 # height of the et
        vertex = self.normal.interpolate(height)
        ETvertices = [(self.end.x, self.end.y), (vertex.x, vertex.y), (self.start.x, self.start.y)]

        return ETvertices

    def square(self): # Generates points that form a square geometry with self
        from shapely import Point, LineString
        
        # # # Make the sides of the square perpendicular to self # # #
        side1Start = (self.end.x, self.end.y)
        side1End = self.normalHelper(abs(self.dx), abs(self.dy), side1Start[0], side1Start[1], False)

        side3End = (self.start.x, self.start.y)
        side3Start = self.normalHelper(abs(self.dx), abs(self.dy), side3End[0], side3End[1], False)

        squareVertices = [side1Start, side1End, side3Start, side3End]

        return squareVertices

    def lineOfSightVertices(self): # Store the test vertices for the line of sight shape
        from shapely import Point, LineString

        """To be used alongside the source surface S1-- Surface((2, 1), (1, 1))-- which is defined in the lineOfSightPlot function
        in the test functions."""

        #geometryVertices = [(1, 1), (1, 7), (7, 7), (7, 1), (5, 1), (5, 3), (3, 3), (3, 1)]
        geometryVertices = [(1, 1), (1, 5), (2, 6), (1, 7), (7, 7), (7, 1), (5, 1), (5, 3), (3, 3), (3, 1)]


        return geometryVertices



        








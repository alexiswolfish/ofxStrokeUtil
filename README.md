ofxStrokeUtil
=============

an openFrameworks addon for comparing/analyzing 2D gestures and contours using ofPath

created in collaboration with Golan Levin @ the Frank-Ratchye Studio for Creative Inquiry

-------------------------------
<img src="https://raw.github.com/alexiswolfish/ofxStrokeUtil/master/example-GraffitiAnalyzer/exampleScreen2.png" width="600">
<br>
##Description
***
ofxStrokeUtil provides a number of ways to compare and quantify lines, drawings, or shapes represented by an ofPath or ofPolyline. Some useful features include…<br>
#####Convex Hull
the smallest set of points that completely envelop your gesture. 
#####Orientation
the axis the gesture is most strongly aligned to, and a metric to determine the strength of that alignment
#####BoundingBox + Oriented Bounding Box
The smallest box that the gesture can be contained in, aligned to either the XY axis, or the axis of orientation.  
<br>
…intersections, corners, various centroid functions, moments, and more. Check out the .cpp file for more detailed descriptions of the various functions and what they do. 
<br><br>
<img src="https://raw.github.com/alexiswolfish/ofxStrokeUtil/master/example-GraffitiAnalyzer/exampleScreen1.png" width ="600">



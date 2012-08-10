
#include "ofxStrokeUtil.h"
#include <math.h>



float ofxStrokeUtil::getStdDevDistanceFromCentroid(ofPath p){
    float stdDevDistance = 0;
    float meanDistance = getMeanDistanceFromCentroid(p);
    ofPoint centroid = getCentroid(p);
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            stdDevDistance += getStdDevDistanceFromPoint(polyline, centroid, meanDistance);
        }
    }
    return stdDevDistance/p.getOutline().size();
}
float ofxStrokeUtil::getStdDevDistanceFromCentroid(ofPolyline p){
    float stdDevDistance = 0;
    float meanDistance = getMeanDistanceFromCentroid(p);
    ofPoint centroid = p.getCentroid2D();
    if(p.getVertices().size() > 2){
        for(ofPoint vertex : p.getVertices()){
            ofPoint dist = vertex - centroid;
            float difFromMean = dist.length() - meanDistance;
            stdDevDistance += difFromMean*difFromMean;
        }
        
        stdDevDistance /= p.getVertices().size();
        stdDevDistance = sqrt(stdDevDistance);
    }
    return stdDevDistance;
}
float ofxStrokeUtil::getStdDevDistanceFromPoint(ofPolyline p, ofPoint point, float meanDistance){
    float stdDevDistance = 0;
    if(p.getVertices().size() > 2){
        for(ofPoint vertex : p.getVertices()){
            ofPoint dist = vertex - point;
            float difFromMean = dist.length() - meanDistance;
            stdDevDistance += difFromMean*difFromMean;
        }
        
        stdDevDistance /= p.getVertices().size();
        stdDevDistance = sqrt(stdDevDistance);
    }
    return stdDevDistance;

}
/*---------------------------------------------------
 *---------------------------------------------------*/
float ofxStrokeUtil::getMeanDistanceFromCentroid(ofPath p){
    ofPoint centroid = getCentroid(p);
    float meanDist = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            meanDist += getMeanDistanceFromPoint(polyline, centroid);
        }
    }
    return meanDist/p.getOutline().size(); //DEBUG::is this correct?
}
float ofxStrokeUtil::getMeanDistanceFromCentroid(ofPolyline p){
    ofPoint centroid = p.getCentroid2D();
    float meanDist =0;
    if(p.getVertices().size() > 2){
        for(ofPoint vertex : p.getVertices()){
            ofPoint dist = vertex - centroid;
            meanDist += dist.length();
        }
        meanDist/= p.getVertices().size();
    }
    return meanDist;
}
float ofxStrokeUtil::getMeanDistanceFromPoint(ofPolyline p, ofPoint point){
    float meanDist =0;
    if(p.getVertices().size() > 2){
        for(ofPoint vertex : p.getVertices()){
            ofPoint dist = vertex - point;
            meanDist += dist.length();
        }
        meanDist/= p.getVertices().size();
    }
    return meanDist;
}

/*---------------------------------------------------*
 *---------------------------------------------------*/
float ofxStrokeUtil::getStdDevSpeed(ofPath p){
    float speed = 0;
    if(p.getOutline().size() > 0){
        int meanSpeed = getMeanSpeed(p);
        for(ofPolyline polyline: p.getOutline()){
            speed += getStdDevSpeed(polyline, meanSpeed);
        }
    }
    return speed/p.getOutline().size();
}

//for a single polyline, input getMeanSpeed(p);
float ofxStrokeUtil::getStdDevSpeed(ofPolyline p, float meanSpeed){
    float speed = 0;
    float count = 0;
    
    if(p.getVertices().size() > 0){
        for(int i=1; i<p.getVertices().size(); i++){
            ofPoint a = p.getVertices()[i-1];
            ofPoint b = p.getVertices()[i];
            
            ofPoint velocity = b-a;
            float differenceFromMean = velocity.length() - meanSpeed;
            speed += differenceFromMean*differenceFromMean;
            count ++;
        }
    }
    if(count>0){
        speed/=count;
        speed = sqrt(speed);
    }
    return speed;
}

/*---------------------------------------------------
 *---------------------------------------------------*/
float ofxStrokeUtil::getMeanSpeed(ofPath p){
    float meanSpeed = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            meanSpeed += getMeanSpeed(polyline);
        }
    }
    return meanSpeed; //DEBUG :: divide out by #polylines?
}

float ofxStrokeUtil::getMeanSpeed(ofPolyline p){ //DEBUG :: get rid of count variable
    float meanSpeed = 0;
    if(p.getVertices().size() > 1){ //DEBUG :: is this necessary?
        float count = 0;
        for(int i=1; i<p.getVertices().size(); i++){
            ofPoint a = p.getVertices()[i-1];
            ofPoint b = p.getVertices()[i];
            
            ofPoint velocity = b-a;
            meanSpeed += velocity.length();
            count++;
        }
        if(count>0){
            return meanSpeed/= float(count);
        }
    }
    return meanSpeed;
}

/*---------------------------------------------------
 *---------------------------------------------------*/
ofPoint ofxStrokeUtil::getMeanVelocity(ofPath p){
    ofPoint meanVel;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            meanVel += getMeanVelocity(polyline);
        }
    }
    return meanVel/p.getOutline().size(); 
}

ofPoint ofxStrokeUtil::getMeanVelocity(ofPolyline p){
    
    ofPoint meanVel;
    for(int i=1; i<p.getVertices().size(); i++){
        ofPoint a = p.getVertices()[i-1];
        ofPoint b = p.getVertices()[i];
        
        ofPoint velocity = b-a;
        meanVel += velocity;
    }
    if(p.getVertices().size() > 0){
        meanVel /= p.getVertices().size();
    }
    return meanVel;
}

/*---------------------------------------------------
 Returns the orientation of the shape, as an angle from
 the horizontal axis, and the "orientedness", how strongly
 the shape follows that axis, as a float. 
 
 **FLAG || return as struct
 *---------------------------------------------------*/
ofVec2f ofxStrokeUtil::getOrientation(ofPath p){
    ofVec2f eigenvector;
    if(p.getOutline().size() > 0){
        ofPoint centroid = getCentroid(p); //DEBUG:: adjust if function moved to ofPath
        
        float XXsum = 0;
        float YYsum = 0;
        float XYsum = 0;
        
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                ofPoint d = vertex - centroid;
                XXsum += d.x*d.x;
                YYsum += d.y*d.y;
                XYsum += d.x*d.y;
            }
        }
       ofVec2f orientation = calculateMajorAxis(YYsum, -XYsum, -XYsum, XXsum);
        return orientation;
    }
    return ofVec2f(0,0);
}


ofVec2f ofxStrokeUtil::getVelocityOrientation(ofPath p){
    
    float XXsum = 0; 
    float YYsum = 0; 
    float XYsum = 0;
    
    if(p.getOutline().size() > 0){
        ofPoint meanVel = getMeanVelocity(p);
        for(ofPolyline polyline: p.getOutline()){
            for(int i=1; i<polyline.getVertices().size(); i++){
                ofPoint a = polyline.getVertices()[i-1];
                ofPoint b = polyline.getVertices()[i];
                ofPoint velocity = b-a;
                ofPoint difVel = velocity - meanVel;
                
                XXsum += difVel.x*difVel.x;
                YYsum += difVel.y*difVel.y;
                XYsum += difVel.x*difVel.y;
            }
        }
        
        ofVec2f orientation = calculateMajorAxis(YYsum, -XYsum, -XYsum, XXsum);
        return orientation;
    }
    return ofVec2f(0,0);
}


/*---------------------------------------------------
 PRIVATE
 MAJOR AXIS
 
 Finds the major axis that descripes the shape, and how
 strongly oriented to the axis that shape is
 
 returns: ofVec2f containing the orientation axis's
          angle from the horizontal axis and then
          the strength of the orientation.
 *---------------------------------------------------*/
ofVec2f ofxStrokeUtil::calculateMajorAxis(float A, float B, float C, float D){
    
    //solve for eigenvalues using the Quadratic formula
    //eigenvalues are the roots of the equation det(lambda*I-T) = 0;
    float eigenvalue1, eigenvalue2;
    float a = 1.0;
    float b = (0.0-A)-D;
    float c = (A*D)-(B*C);
    float Q = (b*b)-(4.0*a*c);
    
    if(Q>=0){
        float magnitude;
        eigenvalue1 = ((0.0-b)+sqrt(Q))/(2.0*a);
        eigenvalue2 = ((0.0-b)-sqrt(Q))/(2.0*a);
        
        eigenvalue1 = (min(eigenvalue1, eigenvalue2) - A)/B;
        magnitude = sqrt(1.0 + eigenvalue1*eigenvalue1);
        
        if((magnitude ==0) || (isnan(magnitude))){
            return ofVec2f(0,0);
        }
        else{
            float orientation = atan2(float(1.0/magnitude), float(eigenvalue1/magnitude)); //in angle from horizontal axis
            float anisitrophy = log(1.0+eigenvalue2); //how strongly oriented it is along the angle
            return(ofVec2f(orientation,anisitrophy));
        }
    }
    else {
        return ofVec2f(0,0);
    }
}

/*---------------------------------------------------*
 *---------------------------------------------------*/
float ofxStrokeUtil::getArcLength(ofPath p){
    
    float arcLength;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline : p.getOutline()){
            arcLength += getArcLength(polyline);
        }  
    }
    return arcLength;
}
float ofxStrokeUtil::getArcLength(ofPolyline p){
    
    float arcLength;
    if(p.getVertices().size()>1){
        for(int i=1; i<p.getVertices().size(); i++){
            ofPoint a = p.getVertices()[i-1];
            ofPoint b = p.getVertices()[i];
            ofPoint d = b-a;
            arcLength += d.length();
        }
    }
    return arcLength;
}

/*---------------------------------------------------*
 *---------------------------------------------------*/
float ofxStrokeUtil::getArea(ofPath p){
    
    float area = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            area += polyline.getArea();
        }
    }
    return area;
}

/*---------------------------------------------------*
 returns the true "average" center point of the given path
 *---------------------------------------------------*/
ofPoint ofxStrokeUtil::getCentroid(ofPath p){
    ofPoint centroid;
    float numPts = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                centroid += vertex;
                numPts++;
            }
        }
        if(numPts > 0)
            centroid/=numPts;
    }
    return centroid;
}

/*---------------------------------------------------*
 *---------------------------------------------------*/
ofPoint ofxStrokeUtil::getCenterOfMass(ofPath p){
    
}
    
//this might make more sense to put in ofPath
ofRectangle ofxStrokeUtil::getBoundingBox(ofPath p){
    //based off of the corresponding ofPolyline function
    ofRectangle box;
    if(p.getOutline().size() > 0){
        box.set(MAXFLOAT,MAXFLOAT,0,0);
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                if(vertex.x < box.x){
                    box.x = vertex.x;
                }
                if(vertex.x > box.width){
                    box.width = vertex.x;
                }
                if(vertex.y < box.y){
                    box.y = vertex.y;
                }
                if(vertex.y > box.height){
                    box.height = vertex.y;
                }
            }
        }
        box.width -= box.x;
        box.height -= box.y;
    }
    return box;
}

/*---------------------------------------------------*
 *---------------------------------------------------*/
float ofxStrokeUtil::getAspectRatio(ofPath p){
    //DEBUG :: adjust accordingly if getBoundingBox is moved into ofPath
    ofRectangle bounds = getBoundingBox(p);
    return (bounds.width)/(bounds.height);
}

float ofxStrokeUtil::getAspectRatio(ofPolyline p){
    ofRectangle bounds = p.getBoundingBox();
    return (bounds.width)/(bounds.height);
}


/*---------------------------------------------------*
 Returns the cumulative angle through which the mark moves
 *---------------------------------------------------*/
float ofxStrokeUtil::getTotalAngle(ofPath p){
    
    float totalAngle = 0;
    if(p.getOutline().size() > 0){
        float meanAbsAngle = 0;
        for(ofPolyline polyLine: p.getOutline()){
            if(polyLine.getVertices().size() > 2){
                for(int i=0; i<polyLine.getVertices().size(); i++){
                    ofPoint a, b, c;
                    if(i == 0){
                        a = ofVec2f(0,0);
                        b = ofVec2f(0,0);
                    }
                    else if(i == 1){
                        a = ofVec2f(0,0);
                        b = polyLine.getVertices()[i-1];
                    }
                    else{
                        a = polyLine.getVertices()[i-2];
                        b = polyLine.getVertices()[i-1];
                    }
                    c = polyLine.getVertices()[i];
                    
                    float theta = getJointAngle(a,b,c);
                    totalAngle += theta;
                }
            }
        }
    }
    return totalAngle;
}

float ofxStrokeUtil::getTotalAbsoluteAngle(ofPath p){
    float totalAngle = 0;
    float nCorners = 0;
    float cornerThreshold = (30*PI)/180;
    
    if(p.getOutline().size() > 0){
        float meanAbsAngle = 0;
        for(ofPolyline polyLine: p.getOutline()){
            if(polyLine.getVertices().size() > 2){
                for(int i=0; i<polyLine.getVertices().size(); i++){
                    ofPoint a, b, c;
                    if(i == 0){
                        a = ofVec2f(0,0);
                        b = ofVec2f(0,0);
                    }
                    else if(i == 1){
                        a = ofVec2f(0,0);
                        b = polyLine.getVertices()[i-1];
                    }
                    else{
                        a = polyLine.getVertices()[i-2];
                        b = polyLine.getVertices()[i-1];
                    }
                    c = polyLine.getVertices()[i];
                    
                    float theta = abs(getJointAngle(a,b,c));
                    
                    totalAngle += theta;
                    
                    if((theta > cornerThreshold) && (i>1) && (i<(polyLine.getVertices().size()-1)))
                        nCorners++;
                }
            }
        }
    }
    corners = nCorners;
    return totalAngle;
}
float ofxStrokeUtil::getMeanAngle(ofPath p){
    float angle = getTotalAngle(p);
    float count = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            count += polyline.getVertices().size();
        }
    }
    if(count > 2)
        return (angle/count)*100; //DEBUG :: ask golan for logic of /100
    else
        return 0;
}
float ofxStrokeUtil::getMeanAbsoluteAngle(ofPath p){
    float angle = getTotalAbsoluteAngle(p);
    float count = 0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            count += polyline.getVertices().size();
        }
    }
    if(count > 2)
        return (angle/count)*10; //DEBUG :: ask golan for logic of /10
    else
        return 0;

}
float ofxStrokeUtil::getStdDevAbsoluteAngle(ofPath p){
    float stdDevAngle = 0;
    float count = 0;
    float meanAbsAngle = getMeanAbsoluteAngle(p);
    
    if(p.getOutline().size() > 0){
        float meanAbsAngle = 0;
        for(ofPolyline polyLine: p.getOutline()){
            if(polyLine.getVertices().size() > 2){
                for(int i=0; i<polyLine.getVertices().size(); i++){
                    ofPoint a, b, c;
                    if(i == 0){
                        a = ofVec2f(0,0);
                        b = ofVec2f(0,0);
                    }
                    else if(i == 1){
                        a = ofVec2f(0,0);
                        b = polyLine.getVertices()[i-1];
                    }
                    else{
                        a = polyLine.getVertices()[i-2];
                        b = polyLine.getVertices()[i-1];
                    }
                    c = polyLine.getVertices()[i];
                    
                    float theta = abs(getJointAngle(a,b,c));
                    float difFromMean = theta - meanAbsAngle;
                    stdDevAngle += (difFromMean*difFromMean);
                    count ++;
                }
            }
        }
        if(count > 2)
            return (sqrt(stdDevAngle/count));
        else
            return 0;
    }
    return stdDevAngle;
    
}
float ofxStrokeUtil::getNumberofCorners(ofPath p){
    getTotalAbsoluteAngle(p);
    return corners;
}
       
//Return angle between two segments
float ofxStrokeUtil::getJointAngle(ofVec2f a, ofVec2f b, ofVec2f c){
    float angle = 0;
    ofVec2f dBA = b-a;
    ofVec2f dCB = c-b;
    dBA.x += 0.00001; //DEBUG :: why do we do this?
    dCB.x += 0.00001;
    
    if((dBA.length() > 0) && (dCB.length() > 0)){
        float slopeBA = dBA.y/dBA.x;
        float slopeCB = dCB.y/dCB.x;
        float angleBA = atan(slopeBA);
        float angleCB = atan(slopeCB);
        
        if(dBA.x >0){
            if(dBA.y >0){
                if(dCB.x > 0){
                    if(dCB.y > 0)
                        angle = angleBA + (PI - angleCB);
                    else
                        angle = PI + angleBA - angleCB;
                }
                else if(dCB.x < 0){
                    if(dCB.y > 0)
                        angle = angleBA - angleCB;
                    else if(dCB.y <= 0){
                        if((angleBA - angleCB) > 0)
                            angle = angleBA - angleCB;
                        else
                            angle = TWO_PI + angleBA - angleCB;
                    }
                }
            }
            else if(dBA.y <= 0){
                if(dCB.x > 0){
                    if(dCB.y > 0)
                        angle = PI + angleBA - angleCB;
                    else
                        angle = PI + angleBA - angleCB;
                }
                else if(dCB.x <0){
                    if(dCB.y >0){
                        float dif = (0-angleBA)+angleCB;
                        if(dif > 0)
                            angle = TWO_PI - dif;
                        else
                            angle = 0 - dif;
                    }
                    else{
                        angle = TWO_PI + angleBA - angleCB;
                    }
                }
            }
        }
        else if(dBA.x <0){
            if(dBA.y <= 0){
                if(dCB.x < 0){
                    if(dCB.y <= 0)
                        angle = angleBA + PI - angleCB;
                    else 
                        angle = PI + angleBA - angleCB;
                }
                else if(dCB.x > 0){
                    if(dCB.y <=0)
                        angle = angleBA - angleCB;
                    else{
                        float dif = angleBA - angleCB;
                        if(dif > 0)
                            angle = dif;
                        else
                            angle = TWO_PI + dif;
                    }
                }
            }
            else{
                if(dCB.x < 0){
                    if(dCB.y > 0)
                        angle = PI + angleBA - angleCB;
                    else 
                        angle = PI + angleBA - angleCB;
                }
                else if(dCB.x > 0){
                    if(dCB.y > 0)
                        angle = TWO_PI - angleCB + angleBA;
                    else if(dCB.y <=0){
                        float dif = (0-angleBA) + angleCB;
                        if(dif>0)
                            angle = TWO_PI - dif;
                        else
                            angle = 0 - dif;
                    }
                }
            }
        }
    }
    angle = PI - angle;
    return angle;
}
       


/*---------------------------------------------------*
 Returns the number of times an ofPath intersects with itself
 *---------------------------------------------------*/
vector<ofPoint> ofxStrokeUtil::getSelfIntersections(ofPath tag){
    
    vector<ofPolyline> strokes = tag.getOutline();
    vector<ofPoint> intersections;
    
    for(ofPolyline p1 : strokes){
        for(ofPolyline p2 : strokes){
            for(ofPoint i : getIntersections(p1, p2)){
                if(!contains(i,intersections))
                    intersections.push_back(i);}
        }
    }
    return intersections;
}

bool ofxStrokeUtil::contains(ofPoint a, vector<ofPoint> points){
    for(ofPoint p: points){
        if(a==p)
            return true;
    }
    return false;
}

/*---------------------------------------------------*
 Returns the number of times a single Polyline intersects
 with itself.
 *---------------------------------------------------*/
vector<ofPoint> ofxStrokeUtil::getSelfIntersections(ofPolyline tag){
    vector<ofPoint> intersections;
    vector<ofPoint> aPoints = tag.getVertices();
    if(aPoints.size() > 3){
        for(int i=1; i<aPoints.size(); i++){
            ofPoint a1 = aPoints[i-1];
            ofPoint a2 = aPoints[i];
            for(int j=i+1; j<aPoints.size(); j++){
                ofPoint b1 = aPoints[j-1];
                ofPoint b2 = aPoints[j];
                if(a1!=b1 && a2!=b2){
                    ofPoint test;
                    if(testIntersect(a1, a2, b1, b2, test))
                        intersections.push_back(test);
                }
            }
        }
    }
    return intersections;
}

/*---------------------------------------------------*
 Tests the line a-b against c-d for intersection.
 *---------------------------------------------------*/
bool ofxStrokeUtil::testIntersect(ofPoint a1, ofPoint a2, ofPoint b1, ofPoint b2, ofPoint &i){
   
    ofVec2f a = a2-a1;
    ofVec2f b = b2-b1;
    ofVec2f c = b2-a2;
    
    double f = perpDot(a,b);
    if(!f)
        return false; //lines are parallel

    double aa = perpDot(a,c);
    double bb = perpDot(b,c);
    
    if(f<0){
        if( aa>0 || bb>0 || aa<f || bb<f)
            return false;
    }
    else{
        if(aa<0 || bb<0 || aa>f || bb>f)
            return false;
    }
    double mu = 1.0-(aa/f);
    ofVec2f intersection = ((b2-b1)*mu) + b1;
    i.set(intersection.x, intersection.y);
    if(i!= a1 && i!= a2 && i != b1 && i !=b2)
        return true;
    else {
        return false;
    }
}

double ofxStrokeUtil::perpDot(ofVec2f a , ofVec2f b){
    return (a.y*b.x)-(a.x*b.y);
}
vector<ofPoint> ofxStrokeUtil::getIntersections(ofPolyline a, ofPolyline b){
    vector<ofPoint> intersections;
    vector<ofPoint> aPoints = a.getVertices();
    vector<ofPoint> bPoints = b.getVertices();
    
    for(int j=1; j<aPoints.size(); j++){
        ofPoint a1 = aPoints[j-1];
        ofPoint a2 = aPoints[j];
        for(int k=1; k<bPoints.size(); k++){
            ofPoint b1 = bPoints[k-1];
            ofPoint b2 = bPoints[k];
            ofPoint test;
                if(testIntersect(a1, a2, b1, b2, test)){
                    intersections.push_back(test);
                }
        }
    }
    return intersections;
}

/*---------------------------------------------------*
 Moments
 Adapted from from http://kenai.com/projects/audiveris/sources/hg/content/src/main/omr/math/Moments.java?rev=3011
 
 * Compute the moments for a set of points whose x and y coordinates are
 * provided, all values being normed by the provided unit value.
 * @param x    the array of abscissa values
 * @param y    the array of ordinate values
 * @param dim  the number of points
 * @param unit the length (number of pixels, for example 20) of norming unit
 
[0] =  "weight"
[1] =  "width"
[2] =  "height"
[3] =  "n20"    //x absolute eccentricity
[4] =  "n11"    //xy covariance
[5] =  "n02"    //y absolute eccentricity
[6] =  "n30"    //x signed eccentricity
[7] =  "n21"
[8] =  "n12"
[9] =  "n03"    //y signed eccentricity
[10]=  "h1"     //Assign non-orthogonal centralized moments
[11]=  "h2"     //(invariant to translation/scaling/rotation)
[12]=  "h3"
[13]=  "h4"
[14]=  "h5"
[15]=  "h6"
[16]=  "h7"
[17]=  "xBar"
[18]=  "yBar"
 *---------------------------------------------------*/
vector<float> ofxStrokeUtil::getMoments(ofPath p){
    vector<float> moments;
    vector<ofPoint> points;
    for(int i=0; i<19; i++)
        moments.push_back(0);
    
    int numPoints =0;
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                points.push_back(vertex);
            }
        }
    }
    
    if(points.size() > 3){
        // Normalized Moments
        float unit = 1.0; 
        float n00 = (float) points.size() / (float) (unit * unit);
        float n01 = 0.0;
        float n02 = 0.0;
        float n03 = 0.0;
        float n10 = 0.0;
        float n11 = 0.0;
        float n12 = 0.0;
        float n20 = 0.0;
        float n21 = 0.0;
        float n30 = 0.0;
        
        // Total weight
        float w = points.size(); // For p+q = 0
        float w2 = w * w;  // For p+q = 2
        float w3 = sqrt(w * w * w * w * w); // For p+q = 3
        
        for(ofPoint vertex: points){
            n10 += vertex.x;
            n01 += vertex.y;
        }
        n10/=points.size();
        n01/=points.size();
        
        // width & height
        float  xMin = FLT_MAX;
        float  xMax = FLT_MIN;
        float  yMin = FLT_MAX;
        float  yMax = FLT_MIN;
        for (int i = points.size()-1; i >= 0; i--) {
            float xx = points[i].x;
            if (xx < xMin) {
                xMin = xx;
            }
            if (xx > xMax) {
                xMax = xx;
            }
            
            float yy = points[i].y;
            if (yy < yMin) {
                yMin = yy;
            }
            if (yy > yMax) {
                yMax = yy;
            }
        }
        
        float dx, dy;
        for (int i = points.size()-1; i >= 0; i--) {
            dx = points[i].x - n10;
            dy = points[i].y - n01;
            
            n11 += (dx * dy);
            n12 += (dx * dy * dy);
            n21 += (dx * dx * dy);
            n20 += (dx * dx);
            n02 += (dy * dy);
            n30 += (dx * dx * dx);
            n03 += (dy * dy * dy);
        }
        
        // Normalize
        // p + q = 2
        n11 /= w2;
        n20 /= w2;
        n02 /= w2;
        
        // p + q = 3
        n12 /= w3;
        n21 /= w3;
        n30 /= w3;
        n03 /= w3;
        
        // Assign non-orthogonal centralized moments
        // (invariant to translation & scaling)
        moments[0] = (float) n00; // Weight
        moments[1] = (float) ((xMax - xMin) / unit); // Width
        moments[2] = (float) ((yMax - yMin) / unit); // Height
        moments[3] = (float) n20; // X  absolute eccentricity
        moments[4] = (float) n11; // XY covariance
        moments[5] = (float) n02; // Y  absolute eccentricity
        moments[6] = (float) n30; // X  signed eccentricity
        moments[7] = (float) n21; //
        moments[8] = (float) n12; //
        moments[9] = (float) n03; // Y signed eccentricity
        
        
        // Assign orthogonals moments (Hu set)
        // (Invariant to translation / scaling / rotation)
        moments[10] = (float) (n20 + n02);
        moments[11] = (float) (((n20 - n02) * (n20 - n02)) + (4 * n11 * n11));
        moments[12] = (float) (((n30 - (3 * n12)) * (n30 - (3 * n12))) + ((n03 - (3 * n21)) * (n03 - (3 * n21))));
        moments[13] = (float) (((n30 + n12) * (n30 + n12)) + ((n03 + n21) * (n03 + n21)));
        moments[14] = (float) (((n30 - (3 * n12)) * (n30 + n12) * (((n30 + n12) * (n30 + n12)) - (3 * (n21 + n03) * (n21 + n03)))) + 
                               ((n03 - (3 * n21)) * (n03 + n21) * (((n03 + n21) * (n03 + n21)) - (3 * (n12 + n30) * (n12 + n30)))));
        moments[15] = (float) (((n20 - n02) * (((n30 + n12) * (n30 + n12)) - ((n03 + n21) * (n03 + n21)))) + (4 * n11 * (n30 + n12) * (n03 + n21)));
        moments[16] = (float) ((((3 * n21) - n03) * (n30 + n12) * (((n30 + n12) * (n30 + n12)) - (3 * (n21 + n03) * (n21 + n03)))) - 
                               (((3 * n12) - n30) * (n03 + n21) * (((n03 + n21) * (n03 + n21)) - (3 * (n12 + n30) * (n12 + n30)))));
        
        // Mass center placed here
        moments[17] = (float) (n10); // xBar
        moments[18] = (float) (n01); // yBar
    }
    return moments;
    
}

float ofxStrokeUtil::getHullPointPercentage(ofPath p){
    if(p.getOutline().size() > 0){
        float numPts = 0;
        ofHull *h = new ofHull(p);
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                numPts++;
            }
        }
        if((h->hull.size() > 3) && (numPts > 3)){
            float dMin = 0;
            float dMax = 1;
            float percentage = (float)(h->hull.size())/numPts;
            return min(dMax, max(dMin, percentage));
        }
    }
    return 0;
}

float ofxStrokeUtil::getCompactness(ofPath p){
    float compactness = 0;
    if(p.getOutline().size() > 0){
        ofHull *h = new ofHull(p);
        float arcLength = getArcLength(p);
        float area = h->getArea();
        if((area > 0) && (arcLength > 0)){
            compactness = (arcLength*arcLength)/area;
            compactness = sqrt(compactness);
        }
    }
    return compactness;
}

float ofxStrokeUtil::getPointDensity(ofPath p){
    if(p.getOutline().size() > 0){
        float numPts = 0;
        ofHull *h = new ofHull(p);
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                numPts++;
            }
        }
        if(numPts > 3){
            float area = h->getArea();
            if(area > 0)
                return float(numPts/area);
        }
    }
    return 0;
}

/*-----------------------------------------------*/
ofHull::ofHull(ofPath p){
    
    depth = 0;
    vector<ofPoint> points, above, below;
    
    if(p.getOutline().size() > 0){
        for(ofPolyline polyline: p.getOutline()){
            for(ofPoint vertex: polyline.getVertices()){
                points.push_back(vertex);
            }
        }
    }
    if(points.size() > 3){
        ofPoint left, right;
        //just sort the vector via the x coord rather than
        //looking individually at each point to find the
        //min and max
        std::sort(points.begin(), points.end(), ofHull::xCoordComparator);
        left = points.front();
        right = points.back();
        
        //test each point to see if it is over or under the
        //line formed by min and max. 
        for(ofPoint point: points){
            if( (point != left) && (point != right)){
                //if(onLeft(left, right, point)){
                if(direction(left, right, point) < 0){
                    above.push_back(point);
                }
                else {
                    below.push_back(point);
                }
            }
        }
        
        //add left and right to each set
        above.push_back(left);
        above.push_back(right);
        below.push_back(left);
        below.push_back(right);
        
        //calculate the upper hull
        vector<ofPoint> tempHull;
        golanQuickHull(above, hull, left, right, 0, 0);
        hull.push_back(left);
        hull.push_back(right);
        std::sort(hull.begin(), hull.end(), ofHull::xCoordComparator);
        golanQuickHull(below, tempHull, left, right, 1, 0);
        std::sort(tempHull.begin(), tempHull.end(), ofHull::xCoordComparator);
        for(int i=tempHull.size()-1; i>=0; i--){
            hull.push_back(tempHull[i]);
        }
    }
}

ofHull::ofHull(ofPolyline p){
    
}

void ofHull::renderHull(){
    if(hull.size() > 0){
        ofBeginShape();
        for(ofPoint p : hull)
            ofVertex(p.x, p.y);
        ofEndShape();
    }
}
void ofHull::renderHullPoints(){
    if(hull.size() > 0){
        for(ofPoint p: hull)
            ofCircle(p.x, p.y, 2);
    }
}

float ofHull::getArea(){
    float area = 0;
    if(hull.size() > 2){
        for(int i=0; i < hull.size()-1; i++){
            ofPoint p0 = hull[i];
            ofPoint p1 = hull[i+1];
            area += (p0.x*p1.y)-(p0.y*p1.x);
        }
        area = abs(0.5*area);
    }
    return area;
}

//Alex's slightly cleaner implementation of the onLeft implementation
//DEBUG:: needs testing

//translate the rightmost point to the origin, and then take the cross
//product of left and p. The function will return positive if p is
//on the right, negative if its on the left, or 0 if it is on the line

int ofHull::direction(ofPoint left, ofPoint right, ofPoint p){
    return((left.x-right.x)*(p.y-right.y)) - ((p.x-right.x)*(left.y-right.y));
}

//determine if p is lying on the left or right side of the
//first point lying on the line formed by a and b
bool ofHull::onLeft(ofPoint a, ofPoint b, ofPoint p){
    if(a.x == b.x){
        if(p.x < a.x)
            return true;
        else{
            if(p.x == a.x){
                if (((p.y > a.y) && (p.y < b.y)) || ((p.y > b.y) && (p.y < a.y)))
                    return true;
                else
                    return false;
            }
            else {
                return false;
            }
        }
    }
    else{
        //DEBUG :: super sketch
        float slope = (b.y - a.y) / (b.x - a.x);
        float x3    = (((p.x + slope * (slope * a.x - a.y + p.y)) /  (1.0 + slope * slope)) * 10000.0 );
        float y3    = ((slope * (x3 / 10000 - a.x) + a.y) * 10000.0 );
        
        if (slope == 0.0) {
            if ((p.y * 10000.0 ) > y3) 
                return true; 
            else 
                return false;
        }
        
        else { 
            if (slope > 0.0 ) {
                if (x3 > (p.x * 10000.0 )) 
                    return true; 
                else 
                    return false;
            }
            
            else {
                if ((p.x * 10000.0 ) > x3) 
                    return true; 
                else 
                    return false;
            }
        }
    }
}

//recursive hull finding method
//let's call the line segment formed by left and right, line LR.
//from a given set of points, all on the same side of LR, find
//the farthest, F. Partition the remaining points into 3 subsets
//depending on whether they are inside the triangle formed by
//LRF, on the right side of the line formed by L and F, or on
//on the right of F to R

void ofHull::golanQuickHull(vector<ofPoint> p, vector<ofPoint> &output, ofPoint left, ofPoint right, int dir, int depth){
    int MAX_RECURSION_DEPTH = 1000; //meh
    if((p.size() == 2) || (depth > MAX_RECURSION_DEPTH)){
        ofPoint point = p[0];
        output.push_back(point);
        p.clear();
        return;
    }
    else{
        depth++;
        
        vector<ofPoint> p1; 
        vector<ofPoint> p2;
        ofPoint curPoint;
        int farthestIndex = splitAt(p, left, right);
        ofPoint farthest = p[farthestIndex];
        
        for(int i=0; i<p.size()-2; i++){
            if(i != farthestIndex){
                curPoint = p[i];
                if(dir == 0){ //finding top Hull
                    if(onLeft(left, farthest, curPoint)){
                        p1.push_back(curPoint);
                    }
                    if(onLeft(farthest, right, curPoint)){
                        p2.push_back(curPoint);
                    }
                }
                else{ //finding bottom hull
                    if(!onLeft(left, farthest, curPoint)){
                        p1.push_back(curPoint);
                    }
                    if(!onLeft(farthest, right, curPoint)){
                        p2.push_back(curPoint);
                    }
                }
            }
        }
        
        p1.push_back(left);
        p1.push_back(farthest); //DEBUG:: looks like to me, P1 will only
        p2.push_back(farthest); //ever add left to the output array, and
        p2.push_back(right); //not the farthest it finds
        
        if(dir == 0){
            golanQuickHull(p1, output, left, farthest, 0, depth);
            golanQuickHull(p2, output, farthest, right, 0, depth);
        }
        else{
            golanQuickHull(p1, output, left, farthest, 1, depth);
            golanQuickHull(p2, output, farthest, right, 1, depth);
        }
        
    }
}
void ofHull::quickHull(vector<ofPoint> p, vector<ofPoint> &output, ofPoint left, ofPoint right, int dir, int depth){
    int MAX_RECURSION_DEPTH = 1000; //meh
    if((p.size() == 0) || (depth > MAX_RECURSION_DEPTH)){
        return;
    }
    else{
        depth++;
        
        vector<ofPoint> p1; 
        vector<ofPoint> p2;
        ofPoint curPoint;
        int farthestIndex = splitAt(p, left, right);
        ofPoint farthest = p[farthestIndex];
        output.push_back(farthest);
        
        for(int i=0; i<p.size()-1; i++){
            if(i != farthestIndex){
                curPoint = p[i];
                if(dir == 0){ //finding top Hull
                    if(onLeft(left, farthest, curPoint)){
                        p1.push_back(curPoint);
                    }
                    if(onLeft(farthest, right, curPoint)){
                        p2.push_back(curPoint);
                    }
                }
                else{ //finding bottom hull
                    if(!onLeft(left, farthest, curPoint)){
                        p1.push_back(curPoint);
                    }
                    if(!onLeft(farthest, right, curPoint)){
                        p2.push_back(curPoint);
                    }
                }
            }
        }
        if(dir == 0){
            quickHull(p1, output, left, farthest, 0, depth);
            quickHull(p2, output, farthest, right, 0, depth);
        }
        else{
            quickHull(p1, output, left, farthest, 1, depth);
            quickHull(p2, output, farthest, right, 1, depth);
        }
    }
}


//Find  a point which is certain to be in the Hull, from among a group of points
//All the given points are on the same side of the line formed by l and r,
//so the point with the longest distance perpendicular to this line is 
//the point we are looking for.  Return the index of this point in the Vector.
int ofHull::splitAt(vector<ofPoint> p, ofPoint l, ofPoint r){
    float denom = r.x-l.x;
    if (r.x == l.x) { 
        denom = 0.00001; //DEBUG:: is this neccessary?
    }
    
    float slope = (r.y-l.y)/denom;
    float farthestIndex = 0;
    float maxDist = 0;
    ofPoint perpPoint;
    
    for(int i=0; i<p.size()-1; i++){
        ofPoint testPoint = p[i];
        //find the point on the line between l and r that
        //is perpendicular to the one we're testing
        if(r.y == l.y){ //horizontal line
            perpPoint.x = testPoint.x;
            perpPoint.y = l.y;
        }
        else{
            perpPoint.x = (testPoint.x + slope*(slope*l.x-l.y+testPoint.y))/(1+slope*slope);
            perpPoint.y = slope*(perpPoint.x - l.x) + l.y;
        }
        
        //check if it is farther than our current farthestPoint
        float distance = testPoint.distance(perpPoint);
        if(distance > maxDist){
            maxDist = distance;
            farthestIndex = i;
        }
    }
    
    return farthestIndex;
}








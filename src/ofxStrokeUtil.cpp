
#include "ofxStrokeUtil.h"
#include <math.h>

float ofxStrokeUtil::getStdDevDistanceFromCentroid(ofPath p){
    float stdDevDistance = 0;
    float meanDistance = getMeanDistanceFromCentroid(p);
    ofPoint centroid = getCentroid(p);
    if(p.hasOutline()){
        for(ofPolyline polyline: p.getOutline()){
            stdDevDistance += getStdDevDistanceFromPoint(polyline, centroid, meanDistance);
        }
    }
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
    if(p.hasOutline()){
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
float ofxStrokeUtil::getStandardDeviatedVelocity(ofPath p){
    float velocity = 0;
    if(p.hasOutline()){
        for(ofPolyline polyline: p.getOutline()){
            velocity += getStandardDeviatedVelocity(polyline);
        }
    }
    return velocity;
}

float ofxStrokeUtil::getStandardDeviatedVelocity(ofPolyline p){
    float velocity = 0;
    float count = 0;
    float meanVelocity = getMeanVelocity(p);
    
    if(p.getVertices().size() > 0){
        for(int i=0; i<p.getVertices().size(); i++){
            ofPoint a = p.getVertices()[i-1];
            ofPoint b = p.getVertices()[i];
            
            ofPoint c = b-a;
            float differenceFromMean = c.length() - meanVelocity;
            velocity += differenceFromMean*differenceFromMean;
            count ++;
        }
    }
    if(count>0){
        velocity/=count;
        velocity = sqrt(velocity);
    }
    return velocity;
}

/*---------------------------------------------------
 *---------------------------------------------------*/
float ofxStrokeUtil::getMeanVelocity(ofPath p){
    float meanVel = 0;
    if(p.hasOutline()){
        for(ofPolyline polyline: p.getOutline()){
            meanVel += getMeanVelocity(polyline);
        }
    }
    return meanVel;
}

float ofxStrokeUtil::getMeanVelocity(ofPolyline p){
    float meanVel = 0;
    if(p.getVertices().size() > 1){
        float count = 0;
        for(int i=1; i<p.getVertices().size(); i++){
            ofPoint a = p.getVertices()[i-1];
            ofPoint b = p.getVertices()[i];
            
            ofPoint c = b-a;
            meanVel += c.length();
            count++;
        }
        if(count>0){
            return meanVel/= float(count);
        }
    }
    return meanVel;
}

/*---------------------------------------------------
 Returns: ofVec2f containing the eigenvector of the ofPath's
 tensor matrix****** 
 
 DEBUG :: put this in english. why is this useful? ask golan about COM
 *---------------------------------------------------*/
ofVec2f ofxStrokeUtil::getOrientation(ofPath p){
    ofVec2f eigenvector;
    if(p.hasOutline()){
        ofPoint centroid = getCentroid(p); //DEBUG:: adjust if function moved to ofPath
        
        //DEBUG:: this seems suspect, talk to golan
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
    }
}

/*---------------------------------------------------
 PRIVATE
 MAJOR AXIS
 
 orientation/orientedness is how strongly distrubed along that
 axis the shape is
 
 find a way to seperate these
 *---------------------------------------------------*/
ofVec2f ofxStrokeUtil::calculateMajorAxis(float A, float B, float C, float D){
    
    //solve for eigenvalues using the Quadratic formula
    //eigenvalues are the roots of the equation det(lambda*I-T) = 0;
    float eigenvalue1, eigenvalue2, eigenvector1, eigenvector2;
    float a = 1.0;
    float b = (0.0-A)-D;
    float c = (A*D)-(B*C);
    float Q = (b*b)-(4.0*a*c);
    
    if(Q>=0){
        float factor2, magnitude; //this is a terrible name, what is this variable? eigenvector 1
        eigenvalue1 = ((0.0-b)+sqrt(Q))/(2.0*a);
        eigenvalue2 = ((0.0-b)-sqrt(Q))/(2.0*a);
        
        factor2 = (min(eigenvalue1, eigenvalue2) - A)/B;
        magnitude = sqrt(1.0 + factor2*factor2);
        
        if((magnitude ==0) || (isnan(magnitude))){
            return ofVec2f(0,0);
        }
        else{
            eigenvector1 = atan2((1.0/magnitude), (factor2/magnitude));
            eigenvector2 = log(1.0+eigenvalue2);
            return(ofVec2f(eigenvector1,eigenvector2));
        }
    }
    else {
        return ofVec2f(0,0);
    }
}

/*---------------------------------------------------
 *---------------------------------------------------*/
float ofxStrokeUtil::getArcLength(ofPath p){
    
    float arcLength;
    if(p.hasOutline()){
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

/*---------------------------------------------------
 *---------------------------------------------------*/
float ofxStrokeUtil::getArea(ofPath p){
    
    float area = 0;
    if(p.hasOutline()){
        for(ofPolyline polyline: p.getOutline()){
            area += polyline.getArea();
        }
    }
    return area;
}
    
ofPoint ofxStrokeUtil::getCentroid(ofPath p){
    ofPoint centroid; //DEBUG :: check if this autoinitalizes to 0
    if(p.hasOutline()){
        for(ofPolyline polyline: p.getOutline()){
            centroid += polyline.getCentroid2D();
        }
        centroid/=p.getOutline().size();
    }
    return centroid;
}

//this might make more sense to put in ofPath
ofRectangle ofxStrokeUtil::getBoundingBox(ofPath p){
    //based off of the corresponding ofPolyline function
    ofRectangle box;
    if(p.hasOutline()){
        box.set(0,0,0,0);
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

/*---------------------------------------------------
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

/*---------------------------------------------------
 *---------------------------------------------------*/
int ofxStrokeUtil::getIntersections(ofPolyline a, ofPolyline b){
    int intersections = 0;
    vector<ofPoint> aPoints = a.getVertices();
    vector<ofPoint> bPoints = b.getVertices();
    
    for(int j=1; j<aPoints.size(); j++){
        ofPoint a1 = aPoints[j-1];
        ofPoint a2 = aPoints[j];
        for(int k=1; j<bPoints.size(); k++){
            ofPoint b1 = bPoints[k-1];
            ofPoint b2 = bPoints[k];
            
            if(testIntersect(a1, a2, b1, b2))
                intersections++;
        }
    }
    return intersections;
}

/*---------------------------------------------------
 Returns the number of times an ofPath intersects with itself
 *---------------------------------------------------*/
int ofxStrokeUtil::getSelfIntersections(ofPath tag){
    
    //Paths verses polylines.
    //Either use a template with generics or ALWAYS convert to polyline?
    
    //this needs to be tested to MAKE SURE it doesn't break,
    //regardless of the MODE. I think we're good though
    vector<ofPolyline> strokes = tag.getOutline(); //converts to polyline regardless of MODE
    int intersections = 0;
    
    for(ofPolyline p1 : strokes){
        for(ofPolyline p2 : strokes){
            intersections += getIntersections(p1, p2);
        }
    }
    return intersections;
}

/*---------------------------------------------------
 Returns the number of times a single Polyline intersects
 with itself.
 *---------------------------------------------------*/
int ofxStrokeUtil::getSelfIntersections(ofPolyline tag){
    return(getIntersections(tag, tag));
}

/*---------------------------------------------------
 Tests the line a-b against c-d for intersection.
 *---------------------------------------------------*/
bool ofxStrokeUtil::testIntersect(ofPoint a, ofPoint b, ofPoint c, ofPoint d){
    
    float denominator, numerator, alpha, beta;
    
    float Ax = b.x - a.x;
    float Ay = b.y - a.y;
    float Bx = c.x - d.x;
    float By = c.y - d.y;
    float Cx = a.x - c.y;
    float Cy = a.y - c.y;
    
    denominator = (Ay*Bx) - (Ax*By);
    numerator = (Ax*Cy) - (Ay*Cx);
    
    if(denominator != 0){
        alpha = numerator/denominator;
        if((alpha > 0.0) && (alpha < 1.0)){
            numerator = (By*Cx) - (Bx*Cy);
            beta = numerator/denominator;
            if((beta > 0.0) && (beta < 1.0)){
                return true;
            }
        }
    }
    
    return false;
    
}
/*-----------------------------------------------------------*
 ofStrokeUtil.h
 
 algorithm : Golan Levin 
 OF port   : Alex Wolfe
 @ the studio for creative inquiry
 *-----------------------------------------------------------*/

#pragma once

#include "ofMain.h"

class ofxStrokeUtil{
    
    public:
        ofxStrokeUtil();
        ~ofxStrokeUtil();
    
        //should go into ofPath
        ofRectangle getBoundingBox(ofPath p);
        ofPoint getCentroid(ofPath p);
        float getArea(ofPath p);
        
        //finished
        float getMeanDistanceFromCentroid(ofPath p);
        float getMeanDistanceFromCentroid(ofPolyline p);
        float getMeanDistanceFromPoint(ofPolyline p, ofPoint centroid);
    
        float getStdDevDistanceFromCentroid(ofPath p);
        float getStdDevDistanceFromCentroid(ofPolyline p);
        float getStdDevDistanceFromPoint(ofPolyline p, ofPoint point, float meanDistance);
    
        float getStandardDeviatedVelocity(ofPath p);
        float getStandardDeviatedVelocity(ofPolyline p);
        float getMeanVelocity(ofPath p);
        float getMeanVelocity(ofPolyline p);
        float getAspectRatio(ofPath p);
        float getAspectRatio(ofPolyline p);
        float getArcLength(ofPath p);
        float getArcLength(ofPolyline p);
        int getSelfIntersections(ofPolyline tag);
        int getSelfIntersections(ofPath tag);
    
        //suspect
        ofVec2f getOrientation(ofPath p);
        ofVec4f getVelocityOrientation(ofPath p); //return type?? ofVec4f or vector<float>
    
        float getHullPointPercentage(ofPath p); // if the hull is just an outline, usefullness?
        float getPointDensity(ofPath p);  //number of points/hullArea, usefullness/accuracy?
        float getCompactness(ofPath p); //relies on hull, useful?
        vector<float> getTotalAngle(ofPath p); //split into multiple function? Terrible style to return huge vector of info
    
        //toDo
        float getStdDevAbsoluteAngle(ofPath p);
        float getJointAngle(ofVec2f a, ofVec2f b, ofVec2f c); //changed input from seperate floats
        float getMoments(ofPath p);
    
    private:
        bool testIntersect(ofPoint a, ofPoint b, ofPoint c, ofPoint d);  
        int getIntersections(ofPolyline a, ofPolyline b);
        ofVec2f calculateMajorAxis(float a, float b, float c, float d); 
        
};
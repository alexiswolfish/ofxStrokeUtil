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
        //ofxStrokeUtil();
        //~ofxStrokeUtil();
    
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
    
        float getStdDevSpeed(ofPath p);
        float getStdDevSpeed(ofPolyline p);
        float getMeanSpeed(ofPath p);
        float getMeanSpeed(ofPolyline p);
    
        ofPoint getMeanVelocity(ofPath p);
        ofPoint getMeanVelocity(ofPolyline p);
    
        float getAspectRatio(ofPath p);
        float getAspectRatio(ofPolyline p);
        float getArcLength(ofPath p);
        float getArcLength(ofPolyline p);
        vector<ofPoint> getSelfIntersections(ofPolyline tag);
        vector<ofPoint> getSelfIntersections(ofPath tag);
    
        ofVec2f getOrientation(ofPath p);
        ofVec2f getVelocityOrientation(ofPath p);
    
        float getTotalAngle(ofPath p); 
        float getTotalAbsoluteAngle(ofPath p);
        float getMeanAngle(ofPath p);
        float getMeanAbsoluteAngle(ofPath p);
        float getStdDevAbsoluteAngle(ofPath p);
        float getNumberofCorners(ofPath p);
    
        vector<float> getMoments(ofPath p);
    
        //hull functions
        float getHullPointPercentage(ofPath p); 
        float getPointDensity(ofPath p); 
        float getCompactness(ofPath p); 
    
    
    private:
        float corners;
        float stdDevAbsAngle;
    
        ofPoint* testIntersect(ofPoint a, ofPoint b, ofPoint c, ofPoint d);  
        vector<ofPoint> getIntersections(ofPolyline a, ofPolyline b);
        ofVec2f calculateMajorAxis(float a, float b, float c, float d); 
        float getJointAngle(ofVec2f a, ofVec2f b, ofVec2f c);
        
};

class ofHull{
    
    public:
        ofHull(ofPath p);
        ofHull(ofPolyline p);
    
        void renderHull();
        void renderHullPoints();
        void buildHalfHull(vector<ofPoint> sortedPoints, ofPoint left, ofPoint right,vector<ofPoint> &output, int dir);
        void golanQuickHull(vector<ofPoint> input, vector<ofPoint> &output, ofPoint left, ofPoint right, int dir, int depth);
        void quickHull(vector<ofPoint> input, vector<ofPoint> &output, ofPoint left, ofPoint right, int dir, int depth);
        int direction(ofPoint left, ofPoint right, ofPoint p);
        bool onLeft(ofPoint left, ofPoint right, ofPoint p);
        int splitAt(vector<ofPoint> p, ofPoint l, ofPoint r);
        float getArea();
        vector<ofPoint> hull;
        int depth;
    
    private:
    static bool xCoordComparator(ofPoint a, ofPoint b){return (a.x < b.x);}

};
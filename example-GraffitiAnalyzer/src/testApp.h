#pragma once

#include "ofMain.h"
#include "ofxXmlSettings.h"
#include "ofxUI.h"
#include "ofxStrokeUtil.h"


class testApp : public ofBaseApp{
    
public:
    
    void setup();
    void update();
    void draw();
    
    void guiSetup();
    
    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    void guiEvent(ofxUIEventArgs &e);
    
    /*-------Drawing--------*/
    float strokeWidth;
    ofPath curPath;
    ofRectangle canvas;
    ofColor curColor;
    
    /*---------Gui----------*/
    ofxUICanvas *gui;
    ofImage samplerImg;
    float guiWidth, canvasHeight;
    
    /*---------Utils--------*/
    
    bool useofPath;
    bool newPath;
    bool drawOutline;
    
    bool testHull;
    
    /*----READYSTEADYGO!----*/
    
    ofxStrokeUtil strokeUtil;
    
    
};


#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
    
    guiWidth = 200;
    guiSpace = 8;
    
    //OF colorTheme
    main.set(185,187,189);
    comp1.set(58,56,57);
    comp2.set(238,57,135);
    comp3.set(253,164,5);
    comp4.set(208,83,153);
    comp5.set(253,187,28);
    comp6.set(123,223,232);
    curColor.set(comp4);
    
    //GUI setup
    panelMain.setup("", "main.xml", guiSpace, guiSpace*4);
    panelMain.add(analyzeToggle.setup("<---ANALYZE", false, guiWidth));
    panelMain.add(newButton.setup("NEW",guiWidth));
    panelMain.add(fillToggle.setup("FILL", false, guiWidth));
    panelMain.add(pathToggle.setup("CLOSE PATHS", false, guiWidth));
    panelMain.add(strokeWidthSlider.setup("STROKE WIDTH", 2, 0.001, 5));
    
    newButton.addListener(this,&testApp::clearDisplay);
    
    panelDisplay.setup("display options", "display.xml", guiSpace, panelMain.getHeight()+guiSpace*5);
    panelDisplay.add(hullToggle.setup("CONVEX HULL", false, guiWidth));
    panelDisplay.add(intersectToggle.setup("INTERSECTIONS", false, guiWidth));
    panelDisplay.add(orientToggle.setup("ORIENTATION", false, guiWidth));
    panelDisplay.add(bBoxToggle.setup("BOUNDING BOX", false, guiWidth));
    panelDisplay.add(centroidToggle.setup("CENTROID", false, guiWidth));
    
    panelColor.setup("","color.xml", guiSpace, panelMain.getHeight()+panelDisplay.getHeight()+guiSpace*6);
    panelColor.add(rSlider.setup("R",curColor.r,0,255));
    panelColor.add(gSlider.setup("G",curColor.g,0,255));
    panelColor.add(bSlider.setup("B",curColor.b,0,255));
    
    canvas = *new ofRectangle(guiWidth+guiSpace*3, guiSpace*4, ofGetWindowWidth() - guiWidth-32, ofGetWindowHeight()-16-guiSpace*3);
    swatch = *new ofRectangle(guiSpace, panelColor.getPosition().y + panelColor.getHeight() + guiSpace, guiWidth,guiWidth/2);
    float displayHeight = swatch.y + swatch.height + guiSpace;
    resultsRec = *new ofRectangle(guiSpace, displayHeight, guiWidth, guiSpace*2); 
    
    curPath.setFilled(false);
    cleared = false;
    
    ofEnableSmoothing();


}

//--------------------------------------------------------------
void testApp::update(){
    
    curColor.set(rSlider.value, gSlider.value, bSlider.value);
    curPath.setStrokeWidth(strokeWidthSlider.value);
    curPath.setStrokeColor(curColor);
    curPath.setFillColor(curColor);
    
    curPath.setFilled(fillToggle.value);
    

}

//--------------------------------------------------------------
void testApp::draw(){
    
    ofBackground(main);
    ofSetColor(comp1);
    ofRect(0, 0, ofGetWindowWidth(), guiSpace*3);
    ofSetColor(comp2);
    ofDrawBitmapString("ofxStrokeUtil GRAFFITI ANALYZER EXAMPLE", guiSpace, guiSpace*2);
    ofSetColor(255, 255, 255);
    ofRect(canvas);
    ofSetColor(curColor);
    ofRect(swatch);
    
    panelMain.draw();
    panelDisplay.draw();
    panelColor.draw();
    
    curPath.draw();
    
    if(hullToggle.value){
        ofPushStyle();
        ofFill();
        ofSetColor(comp5, 180);
        ofHull *h = new ofHull(curPath);
        h->renderHull();
        ofSetColor(0,255,255);
        h->renderHullPoints();
        ofPopStyle();
    }
    if(bBoxToggle.value){
        ofPushStyle();
        ofNoFill();
        ofSetColor(comp6);
        ofSetLineWidth(4);
        ofRect(strokeUtil.getBoundingBox(curPath));
        ofPopStyle();
    }
    if(analyzeToggle.value){
        ofPushStyle();
        ofSetColor(comp1);
        ofRect(resultsRec);
        ofSetColor(comp2);
        ofDrawBitmapString("ANALYSIS", resultsRec.x, resultsRec.y+12);
        ofVec2f pos = *new ofVec2f(resultsRec.x, resultsRec.y + resultsRec.height + guiSpace*2);
        for(int i=0; i<results.size(); i++){
            ofDrawBitmapString(results[i], pos.x, pos.y + i*12);
        }
        ofPopStyle();
    }
    if(centroidToggle.value){
        ofPushStyle();
        ofSetColor(242,31,13);
        ofFill();
        ofPoint centroid = strokeUtil.getCentroid(curPath);
        ofCircle(centroid.x, centroid.y, 10);
        
        ofSetColor(comp6);
        for(ofPolyline p : curPath.getOutline()){
            //ofPolyline's normal centroid function fails on self intersecting
            //non closed shapes (even though they force a "close" by connecting
            //the beginning and the end
            ofPoint c = p.getCentroid2D();
            ofCircle(c.x, c.y, 5);
        }
        ofPopStyle();
    }
    if(intersectToggle.value){
        ofPushStyle();
        ofSetColor(35,0,97);
        ofFill();
        vector<ofPoint> intersections = strokeUtil.getSelfIntersections((curPath));
        for(ofPoint p : intersections){
            ofCircle(p.x,p.y,3);
        }
        ofPopStyle();
    }
    
    if(cleared)
        clear();

}

void testApp::clearDisplay(bool &pressed){
    cleared = true;
}
void testApp::clear(){
    ofPushStyle();
    ofSetColor(255, 255, 255);
    ofFill();
    ofRect(canvas);
    curPath.clear();
    ofPopStyle();
}

void testApp::analyze(){
    
    ofPoint meanVel = strokeUtil.getMeanVelocity(curPath);
    results.push_back("area: " + ofToString(strokeUtil.getArea(curPath)));
    results.push_back("corners: " + ofToString(strokeUtil.getNumberofCorners(curPath),2));
    results.push_back("intersections: " + ofToString(strokeUtil.getSelfIntersections(curPath).size()));
    results.push_back("aspectRatio: " + ofToString(strokeUtil.getAspectRatio(curPath),2));
    results.push_back("arcLength: " + ofToString(strokeUtil.getArcLength(curPath),2));
    results.push_back("meanVelocity: (" + ofToString(meanVel.x,2) + ", " + ofToString(meanVel.y,2)+")");
    results.push_back("meanSpeed: " + ofToString(strokeUtil.getMeanSpeed(curPath),2));
    results.push_back("stdDevSpeed: " + ofToString(strokeUtil.getStdDevSpeed(curPath),2));
    results.push_back("totalAngle: " + ofToString(strokeUtil.getTotalAngle(curPath),2));
    results.push_back("totalAbsAngle: " + ofToString(strokeUtil.getTotalAbsoluteAngle(curPath),2));
    results.push_back("meanAngle: " + ofToString(strokeUtil.getMeanAngle(curPath),2));
    results.push_back("meanAbsAngle: " + ofToString(strokeUtil.getMeanAbsoluteAngle(curPath),2));
    results.push_back("stdDevAbsAngle: " + ofToString(strokeUtil.getStdDevAbsoluteAngle(curPath),2));
    results.push_back("meanDistFromCenter: " + ofToString(strokeUtil.getMeanDistanceFromCentroid(curPath),2));
    results.push_back("stdDevDistFromCenter: " + ofToString(strokeUtil.getStdDevDistanceFromCentroid(curPath)));
    results.push_back("pointDensity: " + ofToString(strokeUtil.getPointDensity(curPath)));
    results.push_back("compactness: " + ofToString(strokeUtil.getCompactness(curPath),2));
    results.push_back("hullPointPercentage: " + ofToString(strokeUtil.getHullPointPercentage(curPath),2));
    
}
//--------------------------------------------------------------
void testApp::keyPressed(int key){
    cleared = true;
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
    if(canvas.inside(mouseX, mouseY))
        curPath.lineTo(mouseX, mouseY);
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    if(canvas.inside(mouseX, mouseY)){
        cleared = false;
        curPath.moveTo(mouseX, mouseY);
    }
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
    if(canvas.inside(mouseX, mouseY) && pathToggle.value){
        curPath.close();
    }
    
    if(analyzeToggle.value){
        results.clear();
        analyze();
    }
}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}
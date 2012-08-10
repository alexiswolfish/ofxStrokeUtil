#include "testApp.h"


//--------------------------------------------------------------
void testApp::setup(){
	ofBackground(127);

    useofPath = true;
    newPath = false;
    drawOutline = false;
    testHull = false;
    
    guiWidth = 200;
    strokeWidth = 2;
    
    samplerImg.loadImage("landscape.jpg");
    float resize = guiWidth/samplerImg.width;
    samplerImg.resize(guiWidth, samplerImg.height*resize);
    
    guiSetup();
    
    canvas = *new ofRectangle(guiWidth+24, 8, ofGetWindowWidth() - guiWidth-32, ofGetWindowHeight()-16);
    
    curPath.setFilled(false);
    curColor.set(0,0,0);
    
	}

//--------------------------------------------------------------
void testApp::update(){
    
  //  curPath.getSubPaths()[0];
    curPath.setStrokeWidth(strokeWidth);
    curPath.setStrokeColor(curColor);
    curPath.setFillColor(curColor);

}

//--------------------------------------------------------------
void testApp::draw(){
    ofBackground(201,255,255);
    ofSetColor(255, 255, 255);
    ofRect(canvas);
    
    if(!drawOutline){
        curPath.draw();
    }
    else{
        if(curPath.hasOutline()){
            ofSetLineWidth( strokeWidth );
            ofSetColor(curColor);
			for(ofPolyline p: curPath.getOutline()){
                p.draw();
            }
            cout << "end" << endl;
        }
    }
    if(testHull){
        ofPushStyle();
        ofFill();
        ofSetColor(255, 255, 0, 150);
        ofHull *h = new ofHull(curPath);
        h->renderHull();
        ofSetColor(0,255,255);
        h->renderHullPoints();
        ofPopStyle();
    }

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
        curPath.moveTo(mouseX, mouseY);
    }
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
    if(canvas.inside(mouseX, mouseY)){
      //  curPath.close();
    }
}

//--------------------------------------------------------------

void testApp::guiSetup(){
    float dim = 16;
    vector<string> toggleNames;
    
    toggleNames.push_back("ofPath");
    toggleNames.push_back("ofPolyLine");
    
    gui = new ofxUICanvas(0,0,guiWidth+OFX_UI_GLOBAL_WIDGET_SPACING*2, ofGetHeight());
    gui->addWidgetDown(new ofxUILabel("SIMPLE EXAMPLE", OFX_UI_FONT_LARGE));
    gui->addWidgetDown(new ofxUISpacer(guiWidth, 2)); 
    
    gui->addWidgetDown(new ofxUIRadio( dim, dim, "ofxStrokeUtil datatype",toggleNames,OFX_UI_ORIENTATION_HORIZONTAL));
    gui->addWidgetDown(new ofxUISpacer(guiWidth, 2));
    gui->addWidgetDown(new ofxUIButton(dim, dim, false, "NEW", OFX_UI_FONT_SMALL));
    
    toggleNames.clear();
    toggleNames.push_back("FILL");
    toggleNames.push_back("OUTLINE");
    toggleNames.push_back("STRAIGHT UP");
    gui->addWidgetDown(new ofxUIRadio(dim, dim, "drawStyle", toggleNames, OFX_UI_ORIENTATION_VERTICAL));
    gui->addWidgetDown(new ofxUISpacer(guiWidth, 2)); 
    gui->addWidgetDown(new ofxUISlider(guiWidth, dim, 0, 30, 
                                        strokeWidth, "stroke width")); 
    gui->addWidgetDown(new ofxUIImageSampler(samplerImg.getWidth(), samplerImg.getHeight(),
                                             &samplerImg, "sampler"));
    gui->addWidgetDown(new ofxUISpacer(guiWidth, 2)); 
    gui->addWidgetDown(new ofxUIToggle(dim, dim, testHull, "DRAW HULL"));
    
    ofAddListener(gui->newGUIEvent,this,&testApp::guiEvent);	
    
    
    
}
void testApp::guiEvent(ofxUIEventArgs &e){
    string name = e.widget->getName(); 
	int kind = e.widget->getKind(); 
    
    if(name == "ofPath")
    {
        curPath.setMode(curPath.PATHS);
    }
    else if(name == "ofPolyLine")
    {
        curPath.setMode(curPath.POLYLINES);
    }
    else if(name == "NEW")
	{
        curPath.clear();
        cout << "cleared " << curPath.hasOutline() << endl;
	}
    else if(name == "FILL"){
        //ofxUIToggle *toggle = (ofxUIToggle *)e.widget;
        //curPath.setFilled(toggle->getValue());
        //ofxUIRadio * rad = (ofxUIRadio *)e.widget;
        curPath.setFilled(true);
    }
    else if(name == "OUTLINE"){
        drawOutline = true;
        curPath.setFilled(false);
    }
    else if(name == "STRAIGHT UP"){
        drawOutline = false;
        curPath.setFilled(false);
    }
    else if(name == "stroke width")
	{
		ofxUISlider *slider = (ofxUISlider *) e.widget; 
        strokeWidth = (slider->getScaledValue());
    }
    else if(name == "sampler")
	{
		ofxUIImageSampler *sampler = (ofxUIImageSampler *) e.widget; 
        curColor = sampler->getColor();
	} 
    else if(name == "DRAW HULL"){
        ofxUIToggle *toggle = (ofxUIToggle *) e.widget;
        testHull = toggle->getValue();
    }
}
//--------------------------------------------------------------
void testApp::keyPressed  (int key){

}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

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





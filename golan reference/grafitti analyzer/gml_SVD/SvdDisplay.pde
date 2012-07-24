
//=====================================================================================
void computeProjectedSvdBounds () {

  projectedSvdDataBounds = new Bounds[nProjectionDimensions];

  // compute the min/max bounds of the projected data ...for display purposes.
  int nSvdDatums = SVDData.size();
  for (int n=0; n<nProjectionDimensions; n++) {
    projectedSvdDataBounds[n] = new Bounds(); 
    float minVal = 99999;
    float maxVal =-99999;

    for (int d=0; d<nSvdDatums; d++) {
      float[] p = ((SVDDatum3D)SVDData.get(d)).projectedPoint;
      float dataVal = p[n]; 
      if (dataVal < minVal) {
        minVal = dataVal;
      }
      if (dataVal > maxVal) {
        maxVal = dataVal;
      }
    }
    projectedSvdDataBounds[n].id   = n;
    projectedSvdDataBounds[n].low  = minVal; 
    projectedSvdDataBounds[n].high = maxVal;
  }

  boolean bVerbose = false;
  if (bVerbose) {
    for (int n=0; n<nProjectionDimensions; n++) {
      projectedSvdDataBounds[n].print();
    }
  }
}



//=====================================================================================
void computeSvdDisplayCoordinates() {


  float minX = projectedSvdDataBounds[0].low;
  float maxX = projectedSvdDataBounds[0].high; 

  float minY = projectedSvdDataBounds[1].low; 
  float maxY = projectedSvdDataBounds[1].high; 

  float minZ = 0; 
  float maxZ = 0; 
  if (nProjectionDimensions == 3) {
    minZ = projectedSvdDataBounds[2].low; 
    maxZ = projectedSvdDataBounds[2].high;
  }


  int nSvdDatums = SVDData.size();
  for (int d=0; d<nSvdDatums; d++) {
    SVDDatum3D aSvdDatum = (SVDDatum3D)SVDData.get(d);

    float[] projectedPoint = aSvdDatum.projectedPoint;
    float dataX = projectedPoint[0];
    float dataY = projectedPoint[1]; 

    displayPoint[0] = map(dataX, minX, maxX, svdDrawRectR, svdDrawRectX);
    displayPoint[1] = map(dataY, minY, maxY, svdDrawRectY, svdDrawRectB);

    float dataZ = 0.0; 
    if (nProjectionDimensions == 3) {
      dataZ  = projectedPoint[2]; 
      displayPoint[2] = map(dataZ, minZ, maxZ, svdDrawRectZ0, svdDrawRectZ1);
    }

    aSvdDatum.setDisplayPoint(displayPoint, nProjectionDimensions);
  }
}





//=====================================================================================
void displayProjectedSvdPoints() {

  smooth();
  noFill();
  stroke(0, 0, 0, 64);
  strokeWeight(4); 

  int nSvdDatums = SVDData.size();

  if (bDrawSvdIn2D || (nProjectionDimensions == 2)) {
    beginShape(POINTS);
    for (int d=0; d<nSvdDatums; d++) {
      SVDDatum3D aSvdDatum = (SVDDatum3D) SVDData.get(d);

      float[] displayPoint = aSvdDatum.displayPoint;  
      float plotX = displayPoint[0];
      float plotY = displayPoint[1];
      point(plotX, plotY);
    }
    endShape();
  } 
  
  
  else if ((!bDrawSvdIn2D) && (nProjectionDimensions == 3)) {
    
    pushMatrix();
    float x2 = svdDrawRectX + svdDrawRectW/2;
    float y2 = svdDrawRectY + svdDrawRectH/2;
    translate(-x2,-y2,0);
    
    float sx,sy;
    float sM = 2;
    float sL = svdDrawRectX + sM;
    float sR = svdDrawRectR - sM;
    float sT = svdDrawRectY + sM;
    float sB = svdDrawRectB - sM;
    float plotX, plotY, plotZ;
    float[] displayPoint;
    SVDDatum3D aSvdDatum;
    
    beginShape(POINTS);
    for (int d=0; d<nSvdDatums; d++) {
      
      aSvdDatum = (SVDDatum3D) SVDData.get(d);
      displayPoint = aSvdDatum.displayPoint;  
      
      plotX = displayPoint[0];
      plotY = displayPoint[1];
      plotZ = displayPoint[2];
      
      sx = screenX(plotX, plotY, plotZ);
      sy = screenY(plotX, plotY, plotZ);
      aSvdDatum.setScreenPoint(sx,sy);
      
      if ((sx > sL) && 
          (sx < sR) && 
          (sy > sT) && 
          (sy < sB)){
        vertex(plotX, plotY, plotZ);
      }
    }
    endShape();
    popMatrix();
  }

  
  strokeWeight(1);
}


//=====================================================================================
void displayProjectedSvdAxes() {
  // draw axes projected into SVD space.
  stroke(0, 0, 0); 
  strokeWeight(1); 

  float minX = projectedSvdDataBounds[0].low;
  float maxX = projectedSvdDataBounds[0].high; 

  float minY = projectedSvdDataBounds[1].low; 
  float maxY = projectedSvdDataBounds[1].high; 

  float xAxis = map(0, minX, maxX, svdDrawRectR, svdDrawRectX);
  float yAxis = map(0, minY, maxY, svdDrawRectY, svdDrawRectB);

  if (nProjectionDimensions == 2) {
    line (xAxis, svdDrawRectY, xAxis, svdDrawRectB);
    line (svdDrawRectX, yAxis, svdDrawRectR, yAxis);
  } 
  else if (nProjectionDimensions == 3) {
    
    float minZ = projectedSvdDataBounds[2].low; 
    float maxZ = projectedSvdDataBounds[2].high; 
    
    ///line(minX, minY, minZ,   maxX, maxY, maxZ); 
    
    /*
    
    
    float zAxis = map(0, minZ, maxZ, svdDrawRectZ0, svdDrawRectZ1);
    line (xAxis, svdDrawRectY, svdDrawRectZ0, xAxis, svdDrawRectB, svdDrawRectZ0);
    line (svdDrawRectX, yAxis, svdDrawRectZ0, svdDrawRectR, yAxis, svdDrawRectZ0);
    // remaining axis?
    */
    
  }

  stroke(0);
}

//=====================================================================================
void getMouseFilename() {
  mouseFileName = ""; 

  float minDist = 99999; 
  int nSvdDatums = SVDData.size();
  
  
  //if (nProjectionDimensions == 2) {
  
  
  for (int d=0; d<nSvdDatums; d++) {
    SVDDatum3D aSvdDatum = (SVDDatum3D) SVDData.get(d);
    //float posx = aSvdDatum.displayPoint[0];
    //float posy = aSvdDatum.displayPoint[1];
    float posx = aSvdDatum.screenPoint[0];
    float posy = aSvdDatum.screenPoint[1];

    float dx = mouseX - posx;
    float dy = mouseY - posy;
    float dh2 = dx*dx + dy*dy; 
    if (dh2 < 16) {
      if (dh2 < minDist) {
        minDist = dh2;
        mouseFileName = aSvdDatum.name + ".gml";
        theCurrentDatumId = d;
      }
    }
  }
  boolean bVerbose = false; 
  if (bVerbose) {
    println("mouseFileName = " + mouseFileName);
  }
}


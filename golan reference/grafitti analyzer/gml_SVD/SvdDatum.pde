class SVDDatum3D {
  // Class to contain SVD data loaded from the server. 
  // Each object contains the coordinates computed by the SVD for each tag. 
  
  String   name; 
  float[]  projectedPoint;
  float[]  displayPoint; 
  float[]  screenPoint;
  
  SVDDatum3D (String gmlFilename, float svdx, float svdy, float svdz){
    name = gmlFilename;
    projectedPoint = new float[3]; 
    displayPoint   = new float[3];
    screenPoint    = new float[2];
    
    projectedPoint[0] = svdx;
    projectedPoint[1] = svdy;
    projectedPoint[2] = svdz;
  }
  
  void print(){
    float svdx = projectedPoint[0];
    float svdy = projectedPoint[1];
    float svdz = projectedPoint[2];
    println(name + "\t" + svdx + "\t" + svdy + "\t" + svdz); 
  }
  
 
  void setDisplayPoint(float[] p, int nDimensions) {
    for (int i=0; i<nDimensions; i++) {
      displayPoint[i] = (float) p[i];
    }
  }
  
  void setScreenPoint(float sx, float sy){
    screenPoint[0] = sx;
    screenPoint[1] = sy;
  }

}


// http://www.igvita.com/2007/01/15/svd-recommendation-system-in-ruby/
// http://math.nist.gov/javanumerics/jama/



//==================================================
class Datum {
  
  int nValues;
  String name; 
  
  double[] values;
  float[] projectedPoint; // the coordinate of the Datum in the new projective space. 
  float[] displayPoint;   // the location of the Datum on-screen. 
  
  
  //-------------------------------------------
  Datum (int n, int ndp) {
    nValues = n;
    values = new double[nValues];
    projectedPoint = new float[ndp];
    displayPoint   = new float[ndp];
  }

  //-------------------------------------------
  Datum (String n, int ndp, float[] inVals) {
    name = n;
    projectedPoint = new float[ndp];
    displayPoint   = new float[ndp];

    nValues = inVals.length;
    values = new double[nValues];
    for (int i=0; i<nValues; i++){
      values[i] = (double) inVals[i];
    }
  }

  //==========================================
  void setProjectedPoint(double[] p, int nDimensions) {
    for (int i=0; i<nDimensions; i++) {
      projectedPoint[i] = (float)p[i];
    }
  }
  
  //-------------------------------------------
  float[] getProjectedPoint() {
    return projectedPoint;
  }

  //==========================================
  void setDisplayPoint(float[] p, int nDimensions) {
    for (int i=0; i<nDimensions; i++) {
      displayPoint[i] = (float) p[i];
    }
  }
  
  void setDisplayPoint2d (float x, float y) {
      displayPoint[0] = x;
      displayPoint[1] = y;
  }
  
  void setDisplayPoint3d (float x, float y, float z){
      displayPoint[0] = x;
      displayPoint[1] = y;
      displayPoint[2] = z;
  }
  

  float[] getDisplayPoint() {
    return displayPoint;
  }

  //==========================================
  double getValue (int which) {
    double out = 0; 
    if ((which < nValues) && (which >= 0)) {
      out = values[which];
    }
    return out;
  }

  void setValue (int which, float val) {
    if ((which < nValues) && (which >= 0)) {
      values[which] = (double) val;
    }
  }

  int getNValues() {
    return nValues;
  }

  double[] getValues() {
    return values;
  }
}





/*
  double[][] bobValues = {
 {
 5,5,0,0,0,5
 }
 };
 Matrix bob = new Matrix(bobValues);
 Matrix bobEmbed = bob.times(u2).times(eig2.inverse());
 
 
 
 int printSpacing   = 5; 
 int printPrecision = 4;
 bobEmbed.print(printSpacing, printPrecision);
 double d00 = bobEmbed.get(0,0);
 double d01 = bobEmbed.get(0,1);
 println("00 = " + d00); 
 println("01 = " + d01); 
 */

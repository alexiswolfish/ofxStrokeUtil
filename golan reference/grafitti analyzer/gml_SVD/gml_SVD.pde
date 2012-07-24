

// http://www.igvita.com/2007/01/15/svd-recommendation-system-in-ruby/
// http://math.nist.gov/javanumerics/jama/
// Uses Jama-1.0.2.jar in "code" folder

import processing.opengl.*;
Arcball arcball;


import Jama.*; 
SingularValueDecomposition svd;
Matrix m;

String gmlAnalysisFilename = "_ANALYSIS.tsv"; // or 10000 500

ArrayList<Datum> DatumList;
ArrayList<Datum> DatumsSimilarToCurrentTag;
ArrayList<GMLTag> TagsSimilarToCurrentTag;
int nSimilarTagsToShow = 30; 

double[][] matrixRawData;
float[] projectedDataMinBounds;
float[] projectedDataMaxBounds;
float[] displayPoint;
int  nProjectionDimensions; 
boolean bDrawSvdIn2D;

ArrayList<SVDDatum3D> SVDData;
Bounds projectedSvdDataBounds[];

float   svdDrawMargin;
float   svdDrawRectX;
float   svdDrawRectY;
float   svdDrawRectW;
float   svdDrawRectH;
float   svdDrawRectR;
float   svdDrawRectB;
float   svdDrawRectZ0;
float   svdDrawRectZ1;

boolean bViewComparison = false;
String  mouseFileName = ""; 
String  prevMouseFileName="";

GMLTag theCurrentTag = null; 
int    theCurrentDatumId = -1;

String fieldNames[] = {
  "Filename     ", //      0
  "totalLength  ", //      1
  "aspectRatio  ", //      2 *
  "nStrokes     ", //      3
  "duration     ", //      4
  "meanDistCent ", //      5
  "stdvDistCent ", //      6
  "totalAbsAngle", //      7
  "totalAngle   ", //      8 *
  "meanAngle    ", //      9 *
  "meanAbsAngle ", //     10 *
  "stdvAbsAng   ", //     11 *
  "orientation  ", //     12
  "orientedness ", //     13 *
  "orientDot    ", //     14 *
  "vorientation ", //     15 
  "vorientedness", //     16
  "meanVelocity ", //     17
  "stdvVelocity ", //     18 
  "nIntersects  ", //     19
  "nCorners     ", // 20
  "compactness  ", // 21
  "hullPointPct ", // 22
  "moment03     ", // 23
  "moment04     ", // 24
  "moment05     ", // 25
  "moment06     ", // 26
  "moment07     ", 
  "moment08     ", 
  "moment09     ", 
  "moment10     ", 
  "moment11     ", 
  "moment12     ", 
  "moment13     ", 
  "moment14     ", 
  "moment15     ", 
  "moment16     "
};



// not 2,8,9,10,11,13,14
int fieldSwitches[] = {
  1, 3, 4, 5, 6, 7, 12, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36
};

int recenterMeanOnZero[] = {
  0, // "Filename", //      0
  0, // "totalLength", //   1
  0, // "aspectRatio", //   2
  0, // "nStrokes", //      3
  0, // "duration", //      4
  0, // "meanDistCent", //  5
  0, // "stdvDistCent", //  6
  0, // "totalAbsAngle", // 7
  0, // "totalAngle", //    8
  0, // "meanAngle", //     9
  0, // "meanAbsAngle", // 10
  0, // "stdvAbsAng", //   11
  0, // "orientation", //  12
  0, // "orientedness", // 13
  0, // "orientDot", //    14
  0, // "vorientation", // 15
  0, // "vorientedness", //16
  0, // "meanVelocity", // 17
  0, // "stdvVelocity", // 18 
  0, // "nIntersects", //  19
  0, // "nCorners", //     20
  0, // "compactness", //  21
  0, // "hullPointPct", // 22
  0, // "moment03", //     23
  0, // "moment04", //     24
  0, // "moment05", //     25
  1, // "moment06", //     26
  1, // "moment07", //     27
  1, // "moment08", //     28
  1, // "moment09", //     29
  0, // "moment10", //     30
  0, // "moment11", //     31
  0, // "moment12", //     32
  0, // "moment13", //     33
  0, // "moment14", //     34
  0, // "moment15", //     35
  0  // "moment16"  //     36
};

int zeroesAllowed[] = {
  0, // "Filename", //      0
  0, // "totalLength", //   1
  0, // "aspectRatio", //   2
  0, // "nStrokes", //      3
  0, // "duration", //      4
  0, // "meanDistCent", //  5
  0, // "stdvDistCent", //  6
  1, // "totalAbsAngle", // 7
  1, // "totalAngle", //    8
  1, // "meanAngle", //     9
  1, // "meanAbsAngle", // 10
  1, // "stdvAbsAng", //   11
  1, // "orientation", //  12
  1, // "orientedness", // 13
  1, // "orientDot", //    14
  1, // "vorientation", // 15
  0, // "vorientedness", //16
  1, // "meanVelocity", // 17
  1, // "stdvVelocity", // 18 
  1, // "nIntersects", //  19
  1, // "nCorners", 
  1, // "compactness", 
  1, // "hullPointPct", 
  1, // "moment03", 
  1, // "moment04", 
  1, // "moment05", 
  1, // "moment06", 
  1, // "moment07", 
  1, // "moment08", 
  1, // "moment09", 
  1, // "moment10", 
  1, // "moment11", 
  1, // "moment12", 
  1, // "moment13", 
  1, // "moment14", 
  1, // "moment15", 
  1  // "moment16"
};


int    nTotalFields; 
Bounds fieldBoundsArray[]; 

int borderBoxGray = 180;
int xFieldIndex = 0;
int yFieldIndex = 1;

// Arcball in JS http://www.math.tamu.edu/~romwell/arcball_js/index.html



//=====================================================================================
void setup() {

  bDrawSvdIn2D = false;
  nProjectionDimensions  = 3;


  if ((bDrawSvdIn2D == true) || (nProjectionDimensions == 2)) {
    size(1024, 768);
  } 
  else if ((!bDrawSvdIn2D) && (nProjectionDimensions == 3)) {
    size(1024, 768, OPENGL);
    hint(DISABLE_OPENGL_2X_SMOOTH);
  }


  svdDrawMargin = 34;
  svdDrawRectW = 380;
  svdDrawRectH = 380;
  svdDrawRectX = svdDrawMargin;
  svdDrawRectY = svdDrawMargin;
  svdDrawRectR = svdDrawRectX + svdDrawRectW;
  svdDrawRectB = svdDrawRectY + svdDrawRectH;
  svdDrawRectZ0 =  200;
  svdDrawRectZ1 = -200;

  float arcballCenterX = svdDrawRectX + svdDrawRectW/2.0;
  float arcballCenterY = svdDrawRectY + svdDrawRectH/2.0;
  arcball = new Arcball(arcballCenterX, arcballCenterY, svdDrawRectW/2.0); 


  textMode(SCREEN);
  displayPoint = new float[nProjectionDimensions]; 


  nTotalFields = fieldSwitches.length;
  DatumList = new ArrayList<Datum>(); 


  DatumsSimilarToCurrentTag = new ArrayList<Datum>(); 
  TagsSimilarToCurrentTag   = new ArrayList<GMLTag>();

  // Pretending to be the server ------------------------------------------------------------
  loadGMLAnalysisFile();               // Obtain the Tag analyses, the 36 metrics.
  preprocessDatumsForSvd();            // Compute the metrics' stdev's and constrain outliers
  performSvdAndSaveSvdAnalysisFile();  // Perform SVD and save the "_SVDAnalysis.tsv";

  // This is the client ---------------------------------------------------------------------
  loadSvdAnalysisFile();               // Load "_SVDAnalysis.tsv"; populate SVDData
  computeProjectedSvdBounds();         // Compute the bounds, for display purposes. 
  computeSvdDisplayCoordinates();      // Within those bounds, map the data for displauy.



  // This is only if we intend to plot the data according to individual metrics
  boolean bDiagnostic = true;
  if (bDiagnostic) {
    computeBoundsAndStatsForAllDatumFields();
  }
}




//=====================================================================================
void loadGMLAnalysisFile() {
  // this loads the "ANALYSIS" file which synopsizes the ~36 different statistical metrics of the marks. 

  float inVals[] = new float[nTotalFields]; 

  String analysisLines[] = loadStrings(gmlAnalysisFilename);
  int nLines = analysisLines.length;
  int nFaulty = 0; 

  for (int i=0; i<nLines; i++) {
    boolean bGMLAnalysisIsOK = true;
    String aLine = analysisLines[i]; 
    String aLineValStrs[] = split(aLine, '\t'); 


    String gmlFilename = aLineValStrs[0];
    if (gmlFilename.endsWith(".gml")) {
      gmlFilename = gmlFilename.substring(0, gmlFilename.indexOf('.'));
    }

    for (int j=0; j<nTotalFields; j++) {
      int fieldId = fieldSwitches[j];

      float aVal = 0.0;
      try {
        aVal = (float)(Float.parseFloat (aLineValStrs[fieldId]));
      } 
      catch (NumberFormatException e) {
        bGMLAnalysisIsOK = false;
      }

      aVal = getFieldSpecificPreprocessedValue (aVal, fieldId);
      inVals[j] = aVal;

      if (zeroesAllowed[fieldId] == 0) {
        if (aVal == 0.0) {
          bGMLAnalysisIsOK = false;
        }
      }
    }


    if (bGMLAnalysisIsOK == false) {
      nFaulty++;
    }

    if (bGMLAnalysisIsOK) {
      Datum D = new Datum (gmlFilename, nProjectionDimensions, inVals);
      DatumList.add (D);
    }
  }
  println("nLines = " + nLines);
  println("nFaulty = " + nFaulty);
}



//=====================================================================================
void preprocessDatumsForSvd() {
  computeBoundsAndStatsForAllDatumFields();
  constrainAnalysisOutliers();
}

//=====================================================================================
void computeBoundsAndStatsForAllDatumFields() {
  // create an Array of Bounds, representing the min/max of all Datum fields. 

  fieldBoundsArray = new Bounds[nTotalFields];
  for (int i=0; i<nTotalFields; i++) {
    fieldBoundsArray[i] = new Bounds();
  }

  int nDatums = DatumList.size();
  for (int i=0; i< nTotalFields; i++) {
    float fieldMax = -99999;
    float fieldMin =  99999;

    float M = 0.0; 
    float S = 0.0; 
    int k = 1;

    for (int d=0; d<nDatums; d++) {
      Datum aDatum = (Datum)DatumList.get(d);
      float val = (float) aDatum.getValue(i);
      if (val > fieldMax) {
        fieldMax = val;
      }
      if (val < fieldMin) {
        fieldMin = val;
      }

      // running standard deviation, from 
      // http://stackoverflow.com/questions/895929/how-do-i-determine-the-standard-deviation-stddev-of-a-set-of-values
      float tmpM = M;
      M += (val - tmpM) / (float) k;
      S += (val - tmpM) * (val - M);
      k++;
    }

    float fieldMean = M;
    float fieldStdv = (k > 1) ? (sqrt(S / (k-1))) : 0.0; 

    fieldBoundsArray[i].name = fieldNames[fieldSwitches[i]];
    fieldBoundsArray[i].set(fieldSwitches[i], fieldMin, fieldMax, fieldMean, fieldStdv); 
    // fieldBoundsArray[i].print();
  }
}


//=====================================================================================
void constrainAnalysisOutliers() {

  float maxStdvs = 6.0;

  int nDatums = DatumList.size();
  for (int d=0; d<nDatums; d++) {
    Datum aDatum = (Datum)DatumList.get(d);

    for (int j=0; j<nTotalFields; j++) {

      float aVal = (float) aDatum.getValue(j);
      float mean = fieldBoundsArray[j].mean;
      float stdv = fieldBoundsArray[j].stdv;
      float hiOutlier = mean + maxStdvs * stdv;
      float loOutlier = mean - maxStdvs * stdv;

      aVal = constrain(aVal, loOutlier, hiOutlier); 
      if (recenterMeanOnZero[fieldSwitches[j]] > 0) {
        aVal = map(aVal, loOutlier, hiOutlier, 0, 1); 
        aVal = doubleExponentialSigmoid (aVal, 0.98);
      }

      aDatum.setValue(j, aVal);
    }
  }
}


//------------------------------------------------
float doubleExponentialSigmoid (float x, float a) {

  float epsilon = 0.00001;
  float min_param_a = 0.0 + epsilon;
  float max_param_a = 1.0 - epsilon;
  a = min(max_param_a, max(min_param_a, a));
  a = 1.0-a; // for sensible results

  float y = 0;
  if (x<=0.5) {
    y = (pow(2.0*x, 1.0/a))/2.0;
  } 
  else {
    y = 1.0 - (pow(2.0*(1.0-x), 1.0/a))/2.0;
  }
  return y;
}






//=====================================================================================
float getFieldSpecificPreprocessedValue (float aVal, int fieldId) {

  if (false) {//(fieldId >= 23) && (fieldId <= 36)) {
  } 
  else {
    switch (fieldId) {
    case 15: // vorientation
      if (aVal < HALF_PI) {
        aVal += PI;
      }
      break;

    case 2:  // aspect ratio
      aVal = atan2(aVal, 1.0); 
      break;

    case 3:  // nStrokes
    case 4:  // duration
    case 7:  // totalAbsAngle
    case 17: // meanVelocity
    case 18: // stdvVelocity
    case 19: // nIntersects
      aVal = max(0, aVal);   // sanity
      aVal = log(1.0+ aVal); 
      break;

    case 1:  // totalLength
      aVal = pow(aVal, 0.500); 
      break;
    case 16: //vorientedness
      aVal = pow(aVal, 0.250); 
      break;
    case 20: // nCorners
    case 21: // compactness
      aVal = pow(aVal, 0.300); 
      break;


    case 23:
    case 25:
    case 30: 
    case 31:
    case 32:
    case 33:

      if (aVal <= 0) {
        //println(negCount++);
        aVal = 16;
      }
      else if (aVal == 0) {
        aVal = 16;//log(MAX_FLOAT);
      } 
      else {
        aVal = log(1.0 + (1.0 / aVal));
      }

      break;
    }
  }
  return aVal;
}

int negCount =0; 

//=====================================================================================
void drawFieldFieldComparison() {

  int nDatums = DatumList.size();
  float xMin = fieldBoundsArray[xFieldIndex].low;
  float xMax = fieldBoundsArray[xFieldIndex].high;
  float yMin = fieldBoundsArray[yFieldIndex].low;
  float yMax = fieldBoundsArray[yFieldIndex].high;

  noStroke();
  noSmooth();

  smooth();
  noFill();
  stroke(0, 0, 0, 128);
  strokeWeight(2); 

  //fill(0, 0, 126, 120);
  beginShape(POINTS); 
  for (int d=0; d<nDatums; d++) {
    Datum aDatum = (Datum)DatumList.get(d);
    float valx = (float) aDatum.getValue(xFieldIndex);
    float valy = (float) aDatum.getValue(yFieldIndex);
    float posx = map(valx, xMin, xMax, svdDrawRectX, svdDrawRectR ); 
    float posy = map(valy, yMin, yMax, svdDrawRectB, svdDrawRectY ); 
    aDatum.setDisplayPoint2d (posx, posy); 

    vertex(posx, posy);// 3, 3);
    //text(aDatum.name, posx, posy-10);
  }
  endShape();

  fill(0, 0, 0); 
  int xField = fieldSwitches[xFieldIndex];
  int yField = fieldSwitches[yFieldIndex];
  String heading = ("x = [" + xField + "] " +  fieldNames[xField] + "\ny = [" + yField + "] " +  fieldNames[yField]);
  text (heading, 20, 20) ;
}




void searchForSimilarDatums () {

  // cosine similarity. 
  if (theCurrentDatumId >= 0) {

    HashMap<Datum, Float> similarities = new HashMap<Datum, Float>();
    similarities.clear();


    float Ax, Ay, Az, Amag;
    float Bx, By, Bz, Bmag;
    Ax = Ay = Az = Amag = 0; 
    int nDatums = DatumList.size();

    Datum queryDatum = (Datum)DatumList.get(theCurrentDatumId);
    Ax   = queryDatum.projectedPoint[0];
    Ay   = queryDatum.projectedPoint[1];

    if (nProjectionDimensions == 2) {
      Amag = sqrt(Ax*Ax + Ay*Ay);
      if (Amag > 0) {
        for (int d=0; d<nDatums; d++) {
          Datum aDatum = (Datum)DatumList.get(d);
          Bx = aDatum.projectedPoint[0];
          By = aDatum.projectedPoint[1];
          Bmag = sqrt(Bx*Bx + By*By); 

          if (Bmag > 0) {
            float dotProduct = Ax*Bx + Ay*By;
            float similarity = dotProduct / (Amag * Bmag); 
            similarities.put(aDatum, similarity);
          }
        }
      }
    } 
    else if (nProjectionDimensions == 3) {
      Az = queryDatum.projectedPoint[2];
      Amag = sqrt(Ax*Ax + Ay*Ay + Az*Az);
      if (Amag > 0) {
        for (int d=0; d<nDatums; d++) {
          Datum aDatum = (Datum)DatumList.get(d);
          Bx = aDatum.projectedPoint[0];
          By = aDatum.projectedPoint[1];
          Bz = aDatum.projectedPoint[2];
          Bmag = sqrt(Bx*Bx + By*By + Bz*Bz); 

          if (Bmag > 0) {
            float dotProduct = Ax*Bx + Ay*By + Az*Bz;
            float similarity = dotProduct / (Amag * Bmag); 
            similarities.put(aDatum, similarity);
          }
        }
      }
    }




    List<Map.Entry> entrySet = new ArrayList<Map.Entry>(similarities.entrySet());

    Collections.sort(entrySet, new Comparator<Map.Entry>() {
      public int compare(Map.Entry e0, Map.Entry e1) {
        return ((Float)(e1.getValue())).compareTo((Float)(e0.getValue()));
      }
    }
    );


    DatumsSimilarToCurrentTag.clear();
    TagsSimilarToCurrentTag.clear();

    boolean bAllowSelfMatching = true;
    if (bAllowSelfMatching) {
      int nEntries = entrySet.size();
      int nSearch = min(nSimilarTagsToShow, nEntries);  
      for (int i=0; i<nSearch; i++) {
        Datum aDatum = (Datum) entrySet.get(i).getKey();
        if ((aDatum != null)) {
          float howSimilar = ((Float) entrySet.get(i).getValue());
          DatumsSimilarToCurrentTag.add(aDatum);
        }
      }
    } 
    else {
      int nEntries = entrySet.size();
      int nSearch = min(nSimilarTagsToShow+1, nEntries);  
      for (int i=0; i<nSearch; i++) {
        Datum aDatum = (Datum) entrySet.get(i).getKey();
        if ((aDatum != null) && (aDatum.equals(queryDatum) == false)) {
          float howSimilar = ((Float) entrySet.get(i).getValue());
          DatumsSimilarToCurrentTag.add(aDatum);
        }
      }
    }

    if (DatumsSimilarToCurrentTag.size() > 0) {
      nLoadedSimilarTags = 0; 
      bLoadingSimilarTags = true;
    }

    //-----------------
    // don't use this; currently replaced by the progressive loader, loadNextSimilarTag();
    // loadAllSimilarTags();
  }
}


//=====================================================================================


int     nLoadedSimilarTags = 0; 
boolean bLoadingSimilarTags = false;

void loadNextSimilarTagIfNecessary () {

  if (bLoadingSimilarTags) {
    if ((nLoadedSimilarTags < nSimilarTagsToShow) && (nLoadedSimilarTags < DatumsSimilarToCurrentTag.size())) {
      Datum similarDatum = (Datum)DatumsSimilarToCurrentTag.get(nLoadedSimilarTags);
      if (similarDatum != null) {
        String similarGmlFileName = similarDatum.name + ".gml";
        GMLTag aSimilarTag = loadGMLTag (similarGmlFileName); 
        if (aSimilarTag != null) {
          aSimilarTag.normalizeTag();
          TagsSimilarToCurrentTag.add(aSimilarTag); 
          nLoadedSimilarTags++;
        }
      }
    } 
    else {
      bLoadingSimilarTags = false;
    }
  }
}

//--------------------------------------------------------------
void loadAllSimilarTags() {
  // load the Tags corresponding to the Datums.  
  // currently, not used; replaced by loadNextSimilarTag(), a progressive loader. 

  int similarTagCount = 0; 
  while ( (similarTagCount < nSimilarTagsToShow) && (similarTagCount < DatumsSimilarToCurrentTag.size())) {
    Datum similarDatum = (Datum)DatumsSimilarToCurrentTag.get(similarTagCount);
    if (similarDatum != null) {
      String similarGmlFileName = similarDatum.name + ".gml";
      GMLTag aSimilarTag = loadGMLTag (similarGmlFileName); 
      if (aSimilarTag != null) {
        aSimilarTag.normalizeTag();
        TagsSimilarToCurrentTag.add(aSimilarTag); 
        similarTagCount++;
      }
    }
  }
  nLoadedSimilarTags = similarTagCount;
  bLoadingSimilarTags = false;
}

//=====================================================================================
void drawCurrentTag() {

  noFill();  
  pushMatrix();
  translate(svdDrawRectR + svdDrawMargin, svdDrawMargin); 

  pushMatrix();
  translate(10,10); 
  strokeWeight(1.0); 
  smooth();
  stroke(0, 0, 0); 
  theCurrentTag.renderWithScale(svdDrawRectW-20, svdDrawRectH-20);
  popMatrix();
  
  noSmooth();
  stroke(borderBoxGray);
  rect(0,0, svdDrawRectW, svdDrawRectH); 
  popMatrix();
}

//=====================================================================================
void drawTagsSimilarToCurrentTag() {

  if (theCurrentDatumId >= 0) {

    loadNextSimilarTagIfNecessary();
    int nSimilarTags = TagsSimilarToCurrentTag.size();
    if (nSimilarTags > 0) {

      float omx = svdDrawMargin;
      float omy = svdDrawMargin + svdDrawRectY + svdDrawRectH; 
      int nRow = 10; 
      float rm = 0; // margin
      float rw = (width - 2*omx - (nRow-1)*rm)/(float)nRow;
      float rh = rw; 
      float bm = 10;


      for (int i=0; i<nSimilarTags; i++) {
        GMLTag aSimilarTag = (GMLTag)TagsSimilarToCurrentTag.get(i);
        if (aSimilarTag != null) {

          // compute position of the tag in this row. 
          float rx = (i%nRow)*(rw+rm);
          float ry = omy + (rh+rm)*(i/nRow); 

          pushMatrix();
          translate (omx, 0); 
          translate (rx, ry);

          boolean bDrawRect = false;
          if (bDrawRect || (i==0)) {
            noSmooth();
            stroke(borderBoxGray);
            rect(0, 0, rw, rh);
          }

          noFill();
          smooth(); 
          strokeWeight(0.5);
          stroke(0, 0, 0); 
          pushMatrix();
          translate (bm, bm);
          aSimilarTag.renderWithScale(rw-bm*2, rh-bm*2);
          popMatrix();
          popMatrix();
        }
      }
    }
  }
}

boolean bMouseDragged = false;
int mouseClickX;
int mouseClickY; 

void mousePressed() {
  bMouseDragged = false;
  arcball.mousePressed();
}

void mouseDragged() {
  bMouseDragged = true;
  arcball.mouseDragged();
}



void mouseReleased() {

  if (bMouseDragged == false) {
    if (mouseFileName != "") {
      println("Loading " + mouseFileName); 
      theCurrentTag = loadGMLTag (mouseFileName); 
      if (theCurrentTag != null) {
        theCurrentTag.normalizeTag();

        nLoadedSimilarTags = 0; 
        long then = millis();
        searchForSimilarDatums();
        //long elapsed = millis() - then;
        //println("elapsed = " + elapsed);
      }
    }
  }
}


//=====================================================================================
void draw() {
  background(255); 


  if (bViewComparison) {
    drawFieldFieldComparison();
    getMouseFilename();
  } 
  else {
    // The usual case. 
    hint(DISABLE_DEPTH_TEST);
    pushMatrix();
    float arcballCenterX = svdDrawRectX + svdDrawRectW/2.0;
    float arcballCenterY = svdDrawRectY + svdDrawRectH/2.0;
    translate(arcballCenterX, arcballCenterY); 
    arcball.run();

    displayProjectedSvdPoints ();
    displayProjectedSvdAxes ();
    getMouseFilename ();
    popMatrix();
  }

  if (!mousePressed) {
    if ((mouseFileName.equals("") == false)  && (mouseFileName.equals(prevMouseFileName) == false)) {
      prevMouseFileName = mouseFileName;
      GMLTag aMouseTag = loadGMLTag (mouseFileName); 
      if (aMouseTag != null) {
        theCurrentTag = aMouseTag;
        theCurrentTag.normalizeTag();
      }
    }
  }


  if (theCurrentTag != null) {
    drawCurrentTag();
    drawTagsSimilarToCurrentTag();
  }

  noFill();
  stroke(borderBoxGray);
  noSmooth();
  rect(svdDrawRectX, svdDrawRectY, svdDrawRectW, svdDrawRectH);
}






//=====================================================================================
void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) {
      yFieldIndex = min(yFieldIndex+1, nTotalFields-1);
    } 
    else if (keyCode == DOWN) {
      yFieldIndex = max(0, yFieldIndex-1);
    }

    if (keyCode == LEFT) {
      xFieldIndex = max(0, xFieldIndex-1);
    } 
    else if (keyCode == RIGHT) {
      xFieldIndex = min(xFieldIndex+1, nTotalFields-1);
    }

    //println("Fields: " + xFieldIndex + " " + yFieldIndex) ;
    //println("x = " + fieldNames[fieldSwitches[xFieldIndex]] + "\ty = " + fieldNames[fieldSwitches[yFieldIndex]]);
  } 
  else {
    if (key == 'v') {
      bViewComparison = !bViewComparison;
    }
    if(key == ' ')
      exit();
  }
}


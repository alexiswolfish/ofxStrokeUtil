// Javascript Array: http://www.w3schools.com/jsref/jsref_obj_array.asp
// http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/MORSE/region-props-and-moments.pdf
// http://kenai.com/projects/audiveris/sources/hg/content/src/main/omr/math/Moments.java?rev=3011

// PRESS 's' TO ANALYZE ALL DATA, AND SAVE ANALYSIS FILE. 
// See: saveCSVDataFromAllTags();

Vector gmlFilenameVector;
int currentGMLFileIndex;
Vector data;
GMLTag aTag;
GMLAnalyzer analyzer;
ConvexHull  huller;

boolean bRotateTag90 = true;
boolean bIsTag3D = false;

HashMap<Integer, GMLTag> allTagMap; //this should really be a hashset

PrintWriter fileOutput;
float calculationTimeTotal;
int nCalculatedFiles;

//=====================================================================
void setup() {

  size(600, 600, P3D);

  allTagMap = new HashMap<Integer, GMLTag>();

  findGMLFiles();
  analyzer = new GMLAnalyzer();
  huller   = new ConvexHull();

  String gmlFilename = "";
  currentGMLFileIndex = 0;
  if (gmlFilenameVector.size() > 0) {
    gmlFilename = (String) gmlFilenameVector.elementAt(currentGMLFileIndex);
  }

  loadData (gmlFilename);
  aTag.normalizeTag();
  aTag.makeFlatCopyStroke(); 
  processFile (gmlFilename);


  fileOutput = createWriter("_ANALYSIS.tsv"); 
  calculationTimeTotal = 0.0;
  nCalculatedFiles = 0;


  GMLTag aBogusTag = makeBogusTag();
  aTag = aBogusTag;
  processCurrentTag();
  //aTag.

}






//=====================================================
GMLTag makeBogusTag() {

  GMLTag outTag = new GMLTag();

  GMLStroke s0 = new GMLStroke();
  float tim = 0; 

  float ra = 1.123;
  float rb = 2.317;	
  for (int i=0; i<100; i++) {
    float t = (float)i/100.0 * TWO_PI;
    float r = 50; 
    float x = 100 + r*1.2 * cos(ra*t);
    float y = 100 + r     * sin(rb*t);
    tim += 25.0; // "millis"
    s0.addXYTPoint(x, y, tim);
  }
  outTag.addStroke(s0);

  ra = 1.23;
  rb = 0.71;
  GMLStroke s1 = new GMLStroke();
  for (int i=0; i<75; i++) {
    float t = (float)i/75.0 * TWO_PI; 
    float x = 100 + 80 * cos(ra*t);
    float y = 100 + 40 * sin(rb*t);
    tim += 25.0; // "millis"
    s1.addXYTPoint(x, y, tim);
  }
  outTag.addStroke(s1);

  return outTag;
}


//=====================================================================
// Identifies all GML files in the 'data' folder.
void findGMLFiles() {
  gmlFilenameVector = new Vector();
  String path = sketchPath + "/data";
  String[] filenames = listFileNames(path);
  for (int f=0; f<filenames.length; f++) {
    if (filenames[f].endsWith(".gml")) {
      // println("Found GML file: " + filenames[f]);
      gmlFilenameVector.addElement(filenames[f]);
    }
  }
  println("Found # GML files: " + gmlFilenameVector.size());
}

//=====================================================================
// This function returns all the files in a directory as an array of Strings  
String[] listFileNames(String dir) {
  File file = new File(dir);
  if (file.isDirectory()) {
    String names[] = file.list();
    return names;
  } 
  else {
    // If it's not a directory
    return null;
  }
}


//=====================================================================
// Loads the data from the GML file into a Vector of GMLPoint arrays. 
// Each GMLPoint array is used to store the information of one stroke.
void loadData(String gmlFilename) {

  XMLElement gmlXML = null;
  try {
    gmlXML = new XMLElement(this, gmlFilename);
  } 
  catch (Exception e) {
    println ("Problem loading GML file " + gmlFilename + "!");
  }

  boolean bRotated90 = true; // the classic GMLs from 000000book are rotated!
  boolean bHasBounds = false;
  float   screenBoundsX = 0; 
  float   screenBoundsY = 0; 
  float   aspectRatio = 1.0;

  if (gmlXML != null) {
    // data = new Vector(10);
    aTag = new GMLTag();
    aTag.fileName = gmlFilename;

    int numCh0 = gmlXML.getChildCount();
    if (numCh0 > 0) {
      println("Loading GML data from " + gmlFilename);

      for (int i=0; i<numCh0; i++) {
        XMLElement ch0 = gmlXML.getChild(i);
        String ch0name = ch0.getName();
        if (ch0name.equals("tag")) {

          //------------------
          // Extract header information: upDirection, screenBounds
          XMLElement gmlHeader = ch0.getChild("header");
          if (gmlHeader != null) {
            XMLElement gmlEnvironment = gmlHeader.getChild("environment");
            if (gmlEnvironment != null) {

              // simple detection of rotation request.
              XMLElement gmlUpDirection = gmlEnvironment.getChild("up");
              if (gmlUpDirection != null) {
                XMLElement upX = gmlUpDirection.getChild("x");
                XMLElement upY = gmlUpDirection.getChild("y");
                XMLElement upZ = gmlUpDirection.getChild("z");
                if ((upX != null) && (upY != null) && (upZ != null)) {
                  float upXf = Float.valueOf((upX.getContent()).trim()).floatValue();
                  float upYf = Float.valueOf((upY.getContent()).trim()).floatValue();
                  float upZf = Float.valueOf((upZ.getContent()).trim()).floatValue();
                  if (upYf != 0.0) {
                    bRotated90 = false;
                  }
                }
              }

              XMLElement gmlScreenBounds = gmlEnvironment.getChild("screenBounds");
              if (gmlScreenBounds != null) {
                XMLElement boundXelt = gmlScreenBounds.getChild("x");
                XMLElement boundYelt = gmlScreenBounds.getChild("y");
                if ((boundXelt != null) && (boundYelt != null)) {
                  bHasBounds = true;
                  screenBoundsX = Float.valueOf((boundXelt.getContent()).trim()).floatValue();
                  screenBoundsY = Float.valueOf((boundYelt.getContent()).trim()).floatValue();
                  if ((screenBoundsX > 0) && (screenBoundsY > 0)) {
                    aspectRatio = screenBoundsY / screenBoundsX;
                  }
                }
              }
            }
          }

          //------------------
          // Extract the drawing
          XMLElement drawing = ch0.getChild("drawing");
          int numCh1 = drawing.getChildCount();
          int nStrokes = 0;

          for (int j=0; j<numCh1; j++) {
            XMLElement ch1 = drawing.getChild(j);
            String ch1name = ch1.getName();
            if (ch1name.equals("stroke")) {
              int nPts = ch1.getChildCount();
              // println("Loading data: stroke #" + j + " has " + nPts + " points.");

              GMLStroke aStroke = new GMLStroke(); // Now we read in each stroke.
              for (int p=0; p<nPts; p++) {         // For each point in the stroke
                XMLElement pt = ch1.getChild(p);   // Fetch the p'th point of the stroke

                // Search for the field in that point which has the right name. 
                // Usually, these will be in the order "x", "y", "z", "time", 
                // But we need to do this search in case they're out of order for whatever reason.
                int nPointFields = pt.getChildCount();

                float x, y, z, t; 
                x=y=z=t= 0.0;
                for (int n=0; n<nPointFields; n++) {
                  XMLElement pointField = pt.getChild(n);
                  String pointFieldName = pointField.getName();
                  pointFieldName = (pointFieldName.toLowerCase()).trim();
                  float value = Float.valueOf((pointField.getContent()).trim()).floatValue();

                  if      (pointFieldName.equals("x")) {
                    x = value;
                  } 
                  else if (pointFieldName.equals("y")) {
                    y = value;
                  } 
                  else if (pointFieldName.equals("z")) {
                    z = value;
                  } 
                  else if (pointFieldName.equals("t")) {
                    t = value;
                  } 
                  else if (pointFieldName.equals("time")) {
                    t = value;
                  }
                }
                if (bRotated90) {
                  float w = y;
                  y = 1.0-x;
                  x = w;
                } 
                y *= aspectRatio;
                XYTPoint aPoint = new XYTPoint(x, y, t);
                aStroke.addXYTPoint(aPoint);
              }
              aTag.addStroke(aStroke);
            }
          }
        }
      }
    }

    gmlXML = null;
    System.gc();
  }
}





//=====================================================================
void draw() {
  background(0);
  noFill();
  smooth();

  pushMatrix();
  scale(600, 600);
  stroke(255, 255, 0);
  aTag.render();
  //stroke(255,0,0);
  //aTag.renderBlurredCopy();
  // aTag.renderVelocities();

 // stroke(0, 200, 0); 
 // ArrayList hull = huller.quickHull (aTag);
 // huller.renderHull(hull);


  XYTPoint centroid = analyzer.getCentroid (aTag);
  fill(255, 0, 0); 
  ellipse(centroid.x, centroid.y, 0.03, 0.03); 

  float lr, lx, ly;

  float orientation = (aTag.descriptors.get("orientation"));
  lr = 0.1;
  lx = centroid.x + lr*sin( (orientation));
  ly = centroid.y + lr*cos( (orientation));
  stroke(0, 255, 0);
  line(centroid.x, centroid.y, lx, ly);

  float vorientation = (aTag.descriptors.get("vorientation"));
  lr = 0.1;
  lx = centroid.x + lr*sin( (vorientation));
  ly = centroid.y + lr*cos( (vorientation));
  stroke(255, 0, 255);
  line(centroid.x, centroid.y, lx, ly);


  popMatrix();


  fill(255);
  //String instructions = "Current GML file: " + gmlFilename + "\n";
  //text(instructions, width/2, height-80);
}


//=====================================================================
void mousePressed() {
  processNextFile();
}

void processNextFile() {
  if (gmlFilenameVector.size() > 0) {
    currentGMLFileIndex = (currentGMLFileIndex+1)%(gmlFilenameVector.size());
    String gmlFilename = (String) gmlFilenameVector.elementAt(currentGMLFileIndex);
    processFile(gmlFilename);
  }
}


//=====================================================================
void processCurrentTag () {


  if (true) {
  
    aTag.normalizeTag();
    aTag.makeFlatCopyStroke();
    aTag.makeBlurredCopy();


    float totalLength  = analyzer.getArcLength(aTag); 
    float aspectRatio  = analyzer.getAspectRatio (aTag);
    float nStrokes     = analyzer.getNStrokes(aTag); 
    float duration     = analyzer.getDuration (aTag); 

    XYTPoint C = analyzer.getCentroid (aTag);
    float meanDistCent = analyzer.getMeanDistanceFromCentroid (aTag, C);
    float stdvDistCent = analyzer.getStDevDistanceFromCentroid (aTag, C, meanDistCent);

    float totAngle[]    = analyzer.getTotalAngle (aTag);
    float totalAngle    = totAngle[0];
    float totalAbsAngle = totAngle[1];
    float meanAngle     = totAngle[2];
    float meanAbsAngle  = totAngle[3];
    float nCorners      = totAngle[4];
    float stdvAbsAng    = analyzer.getStDevAbsAngle ( aTag, meanAbsAngle);

    float orient[]      = analyzer.getOrientation (aTag, C); 
    float orientation   = orient[0];
    float orientedness  = orient[1];

    float vorient[]     = analyzer.getVelocityOrientation (aTag);
    float vorientation  = vorient[0];
    float vorientedness = vorient[1];
    float meanVelocity  = vorient[2];
    float stdvVelocity  = vorient[3];

    float orientDot     = analyzer.getOrientationDot (orientation, vorientation); 
    float nIntersects   = (float) getSelfIntersectionCount (aTag);
    float moments[]    = analyzer.getMoments(aTag); 

    ArrayList hull = huller.quickHull (aTag);
    float compactness   = analyzer.getCompactness ( totalLength, hull);
    float hullPointPct  = analyzer.getHullPointPercentage (aTag, hull); 



    aTag.descriptors.put("totalLength", totalLength);
    aTag.descriptors.put("aspectRatio", aspectRatio);
    aTag.descriptors.put("nStrokes", nStrokes);
    aTag.descriptors.put("duration", duration);
    aTag.descriptors.put("meanDistCent", meanDistCent);
    aTag.descriptors.put("stdvDistCent", stdvDistCent);
    aTag.descriptors.put("totalAbsAngle", totalAbsAngle);
    aTag.descriptors.put("totalAngle", totalAngle);
    aTag.descriptors.put("stdvAbsAng", stdvAbsAng); 
    aTag.descriptors.put("meanAngle", meanAngle); 
    aTag.descriptors.put("meanAbsAngle", meanAbsAngle); 
    aTag.descriptors.put("orientation", orientation);  
    aTag.descriptors.put("orientedness", orientedness); 
    aTag.descriptors.put("vorientation", vorientation);   
    aTag.descriptors.put("vorientedness", vorientedness);
    aTag.descriptors.put("orientDot", orientDot);
    aTag.descriptors.put("meanVelocity", meanVelocity);
    aTag.descriptors.put("stdvVelocity", stdvVelocity);
    aTag.descriptors.put("nIntersects", nIntersects);
    aTag.descriptors.put("nCorners", nCorners); 
    aTag.descriptors.put("compactness", compactness); 
    aTag.descriptors.put("hullPointPct", hullPointPct); 

    for (int m=3; m<=16; m++) {
      String momentName = "moment" + nf(m, 2);
      aTag.descriptors.put(momentName, moments[m]);
    }

    if ((allTagMap.containsKey(aTag.hashCode())) == false) {
      allTagMap.put(aTag.hashCode(), aTag);
    }

    long processFileEndTime = millis();
    
    //println("PRINTING ATAG"); 
    //aTag.printDescriptors();

  }
}



//=====================================================================
void processFile (String gmlFilename) {

  long processFileStartTime = millis();
  if (gmlFilename != null) {
    loadData (gmlFilename);
    long calcStartTime = millis();

    aTag.normalizeTag();
    aTag.makeFlatCopyStroke();
    aTag.makeBlurredCopy();


    float totalLength  = analyzer.getArcLength(aTag); 
    float aspectRatio  = analyzer.getAspectRatio (aTag);
    float nStrokes     = analyzer.getNStrokes(aTag); 
    float duration     = analyzer.getDuration (aTag); 

    XYTPoint C = analyzer.getCentroid (aTag);
    float meanDistCent = analyzer.getMeanDistanceFromCentroid (aTag, C);
    float stdvDistCent = analyzer.getStDevDistanceFromCentroid (aTag, C, meanDistCent);

    float totAngle[]    = analyzer.getTotalAngle (aTag);
    float totalAngle    = totAngle[0];
    float totalAbsAngle = totAngle[1];
    float meanAngle     = totAngle[2];
    float meanAbsAngle  = totAngle[3];
    float nCorners      = totAngle[4];
    float stdvAbsAng    = analyzer.getStDevAbsAngle ( aTag, meanAbsAngle);

    float orient[]      = analyzer.getOrientation (aTag, C); 
    float orientation   = orient[0];
    float orientedness  = orient[1];

    float vorient[]     = analyzer.getVelocityOrientation (aTag);
    float vorientation  = vorient[0];
    float vorientedness = vorient[1];
    float meanVelocity  = vorient[2];
    float stdvVelocity  = vorient[3];

    float orientDot     = analyzer.getOrientationDot (orientation, vorientation); 
    float nIntersects   = (float) getSelfIntersectionCount (aTag);
    float moments[]    = analyzer.getMoments(aTag); 

    ArrayList hull = huller.quickHull (aTag);
    float compactness   = analyzer.getCompactness ( totalLength, hull);
    float hullPointPct  = analyzer.getHullPointPercentage (aTag, hull); 

    long calcEndTime = millis();
    long calculationTime = calcEndTime - calcStartTime;

    calculationTimeTotal += calculationTime;
    nCalculatedFiles++; 



    aTag.descriptors.put("totalLength", totalLength);
    aTag.descriptors.put("aspectRatio", aspectRatio);
    aTag.descriptors.put("nStrokes", nStrokes);
    aTag.descriptors.put("duration", duration);
    aTag.descriptors.put("meanDistCent", meanDistCent);
    aTag.descriptors.put("stdvDistCent", stdvDistCent);
    aTag.descriptors.put("totalAbsAngle", totalAbsAngle);
    aTag.descriptors.put("totalAngle", totalAngle);
    aTag.descriptors.put("stdvAbsAng", stdvAbsAng); 
    aTag.descriptors.put("meanAngle", meanAngle); 
    aTag.descriptors.put("meanAbsAngle", meanAbsAngle); 
    aTag.descriptors.put("orientation", orientation);  
    aTag.descriptors.put("orientedness", orientedness); 
    aTag.descriptors.put("vorientation", vorientation);   
    aTag.descriptors.put("vorientedness", vorientedness);
    aTag.descriptors.put("orientDot", orientDot);
    aTag.descriptors.put("meanVelocity", meanVelocity);
    aTag.descriptors.put("stdvVelocity", stdvVelocity);
    aTag.descriptors.put("nIntersects", nIntersects);
    aTag.descriptors.put("nCorners", nCorners); 
    aTag.descriptors.put("compactness", compactness); 
    aTag.descriptors.put("hullPointPct", hullPointPct); 

    for (int m=3; m<=16; m++) {
      String momentName = "moment" + nf(m, 2);
      aTag.descriptors.put(momentName, moments[m]);
    }

    if ((allTagMap.containsKey(aTag.hashCode())) == false) {
      allTagMap.put(aTag.hashCode(), aTag);
    }

    long processFileEndTime = millis();

    aTag.printDescriptors();

    long processFileTime = processFileEndTime - processFileStartTime; 
    float averageCalculationTime = calculationTimeTotal / (float)nCalculatedFiles;
    println("calculationTime = " + calculationTime + "\t" + processFileTime + "\t" + averageCalculationTime + "\tnFiles = " + nCalculatedFiles);
  }
}

//=====================================================================
void saveCSVDataFromAllTags() {

  int nTags = gmlFilenameVector.size();
  for (int i=0; i<nTags; i++) {
    processNextFile();
  }

  int count = 0; 
  for (GMLTag tag: allTagMap.values()) {
    fileOutput.println (tag.toCSV());
    if (count%50 == 0) {
      println("Handling " + count);
    }
    count++;
  }
  println("Finished fileOutput"); 

  fileOutput.flush(); // Writes the remaining data to the file
  fileOutput.close(); // Finishes the file
}

//=====================================================================
// Based on keystrokes, issue commands to export DXF, etc.
void keyPressed() {

  if (key == 's') {
    saveCSVDataFromAllTags();
  }
  if(key == ' '){
    exit();
  }

  //println(aTag.toCSV());
}



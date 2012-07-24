

//=====================================================================
// Loads the data from the GML file into a Vector of GMLPoint arrays. 
// Each GMLPoint array is used to store the information of one stroke.
GMLTag loadGMLTag (String gmlFilename) {

  GMLTag loadedTag = new GMLTag();
  XMLElement gmlXML = null;
  try {
    gmlXML = new XMLElement(this, gmlFilename);
  } 
  catch (Exception e) {
    println ("Problem loading GML file " + gmlFilename + "!");
    return null;
  }

  boolean bRotated90 = true; // the classic GMLs from 000000book are rotated!
  boolean bHasBounds = false;
  float   screenBoundsX = 0; 
  float   screenBoundsY = 0; 
  float   aspectRatio = 1.0;

  if (gmlXML != null) {
    // data = new Vector(10);
    loadedTag.fileName = gmlFilename;

    int numCh0 = gmlXML.getChildCount();
    if (numCh0 > 0) {
      // println("Loading GML data from " + gmlFilename);

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
              loadedTag.addStroke(aStroke);
            }
          }
        }
      }
    }

    gmlXML = null;
    System.gc();
    return loadedTag;
  }
  return null;
}


//=====================================================================
// A data structure for storing GML point data.
class XYTPoint {
  float x;
  float y;
  float t;

  XYTPoint () {
    x = y = t = 0;
  }

  XYTPoint (float _x, float _y, float _t) {
    x = _x;
    y = _y;
    t = _t;
  }

  XYTPoint (XYTPoint src) {
    x = src.x;
    y = src.y;
    t = src.t;
  }

  void set (float _x, float _y, float _t) {
    x = _x;
    y = _y;
    t = _t;
  }

  void set (float _x, float _y) {
    x = _x;
    y = _y;
  }
}





//=====================================================================
class GMLStroke {

  ArrayList pts;
  //---------------------------------
  GMLStroke () {
    pts = new ArrayList();
  }

  //---------------------------------
  void addXYTPoint (XYTPoint p) {
    pts.add(p);
  }

  //---------------------------------
  void render () {
    noFill();
    int nPoints = pts.size();

    beginShape();
    for (int i=0; i<nPoints; i++) {
      XYTPoint p = (XYTPoint) pts.get(i);
      float x = p.x;
      float y = p.y;
      vertex(x, y);
    }
    endShape();
  }

  //---------------------------------
  void renderWithScale (float sx, float sy) {
    noFill();
    int nPoints = pts.size();

    beginShape();
    for (int i=0; i<nPoints; i++) {
      XYTPoint p = (XYTPoint) pts.get(i);
      float x = sx * p.x;
      float y = sy * p.y;
      vertex(x, y);
    }
    endShape();
  }
}

//=====================================================================
class GMLTag {
  String fileName;
  ArrayList strokes;
  
  // For animated drawing
  float durationMillis;
  long  birthdayMillis;

  //---------------------------------
  GMLTag () {
    strokes  = new ArrayList();
  }



  //---------------------------------
  void addStroke (GMLStroke s) {
    strokes.add(s);
  }


  //---------------------------------
  int getNStrokes() {
    return strokes.size();
  }
  //---------------------------------
  void render () {
    stroke(0, 0, 0); 
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke)strokes.get(i);
      s.render();
    }
  }

  //---------------------------------
  void renderWithScale (float sx, float sy) {
    stroke(0, 0, 0); 
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke)strokes.get(i);
      s.renderWithScale (sx, sy);
    }
  }

  //---------------------------------
  void animateTagWithScale (float sx, float sy){
    
    // For animated drawing
    // float durationMillis;
   // long  birthdayMillis;
  }

  //---------------------------------
  void normalizeTag() {

    // Initialize range bounds
    float maxY = -99999;
    float maxX = -99999;
    float maxT = -99999;
    float minY =  99999;
    float minX =  99999;
    float minT =  99999;

    // Compute the bounding ranges of the tag.
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke) strokes.get(i);
      ArrayList pts = s.pts;
      int nPoints = pts.size();
      for (int j=0; j<nPoints; j++) {
        XYTPoint pt = (XYTPoint) pts.get(j);
        float x = pt.x;
        float y = pt.y;
        float t = pt.t;

        if (y > maxY) {
          maxY = y;
        }
        if (x > maxX) {
          maxX = x;
        }
        if (t > maxT) {
          maxT = t;
        }


        if (y < minY) {
          minY = y;
        }
        if (x < minX) {
          minX = x;
        }
        if (t < minT) {
          minT = t;
        }
      }
    }

    // Normalize the marks to the range 0..1, 
    // without changing the aspect ratio of the tag.
    float rangeX = maxX - minX;
    float rangeY = maxY - minY;
    float rangeT = maxT - minT;
    float normalizingRange = (rangeX > rangeY) ? rangeX : rangeY;
    float lesserRange      = (rangeX > rangeY) ? rangeY : rangeX;

    // heuristic for guessing whether the time units are seconds or millis:
    boolean bTimeUnitsAreMillis = (rangeT > 60); 

    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke) strokes.get(i);
      ArrayList pts = s.pts;
      int nPoints = pts.size();
      for (int j=0; j<nPoints; j++) {
        XYTPoint pt = (XYTPoint) pts.get(j);
        pt.x = (pt.x - minX)/normalizingRange;
        pt.y = (pt.y - minY)/normalizingRange;

        if (bTimeUnitsAreMillis) {
          pt.t = pt.t / 1000.0;
        }
      }
    }
  }
}




/*
    int nEntries = entrySet.size();
 for (int i=0; i<nEntries; i++){
 Datum aDatum =     (Datum) entrySet.get(i).getKey();
 float howSimilar = ((Float) entrySet.get(i).getValue());
 println(aDatum.name + ".gml  " + howSimilar); 
 }
 */

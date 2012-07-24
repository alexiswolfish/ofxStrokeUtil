




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
  void addXYTPoint (float _x, float _y, float _t) {
    XYTPoint p = new XYTPoint(_x, _y, _t);
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

    /*
    for (int i=0; i<nPoints; i++) {
     XYTPoint p = (XYTPoint) pts.get(i);
     float x = p.x;
     float y = p.y;
     ellipse(x, y, 0.01, 0.01);
     }
     */
  }
}

//=====================================================================
class GMLTag {
  String fileName;
  ArrayList strokes;
  ArrayList strokesSmoothed;
  GMLStroke flatCopyStroke;
  GMLStroke resampledTag;
  HashMap<String, Float> descriptors;

  //---------------------------------
  GMLTag () {
    strokes         = new ArrayList();
    resampledTag    = new GMLStroke();
    flatCopyStroke  = new GMLStroke();
    descriptors     = new HashMap<String, Float>();
  }

  int hashCode() {
    return descriptors.hashCode();
  }

  boolean equals (GMLTag T) {
    println(this + " Called equals on " + T); 
    return (descriptors.equals(T.descriptors));
  }


  String toString() {
    return "<Tag " + hashCode() + " >\n";
  }


  String toCSV() {
    String out = fileName;

    out += "\t" + nf(descriptors.get("totalLength"), 1, 5); 
    out += "\t" + nf(descriptors.get("aspectRatio"), 1, 5); 
    out += "\t" +    descriptors.get("nStrokes"   ); 
    out += "\t" + nf(descriptors.get("duration"   ), 1, 5); 
    out += "\t" + nf(descriptors.get("meanDistCent"), 1, 5);
    out += "\t" + nf(descriptors.get("stdvDistCent"), 1, 5);
    out += "\t" + nf(descriptors.get("totalAbsAngle"), 1, 5);
    out += "\t" + nf(descriptors.get("totalAngle"), 1, 5);
    out += "\t" + nf(descriptors.get("meanAngle"), 1, 5);
    out += "\t" + nf(descriptors.get("meanAbsAngle"), 1, 5);
    out += "\t" + nf(descriptors.get("stdvAbsAng"), 1, 5);
    out += "\t" + nf(descriptors.get("orientation"), 1, 5);
    out += "\t" + nf(descriptors.get("orientedness"), 1, 5);
    out += "\t" + nf(descriptors.get("orientDot"), 1, 5); 
    out += "\t" + nf(descriptors.get("vorientation"), 1, 5);
    out += "\t" + nf(descriptors.get("vorientedness"), 1, 5);
    out += "\t" + nf(descriptors.get("meanVelocity"), 1, 5);
    out += "\t" + nf(descriptors.get("stdvVelocity"), 1, 5);
    out += "\t"    + descriptors.get("nIntersects");
    out += "\t"    + descriptors.get("nCorners");
    out += "\t" + nf(descriptors.get("compactness"), 1, 5);
    out += "\t" + nf(descriptors.get("hullPointPct"), 1, 5);

    for (int m=3; m<=16; m++) {
      String momentName = "moment" + nf(m, 2);
      out += "\t" + descriptors.get(momentName);
    }

    return out;
  }

  void printDescriptors() {

    println("totalLength   = " + descriptors.get("totalLength")); 
    println("aspectRatio   = " + descriptors.get("aspectRatio")); 
    println("nStrokes      = " + descriptors.get("nStrokes")); 
    println("duration      = " + descriptors.get("duration")); 
    println("meanDistCent  = " + descriptors.get("meanDistCent"));
    println("stdvDistCent  = " + descriptors.get("stdvDistCent"));
    println("totalAbsAngle = " + descriptors.get("totalAbsAngle"));
    println("totalAngle    = " + descriptors.get("totalAngle"));
    println("meanAngle     = " + descriptors.get("meanAngle"));
    println("meanAbsAngle  = " + descriptors.get("meanAbsAngle"));
    println("stdvAbsAng    = " + descriptors.get("stdvAbsAng"));
    println("orientation   = " + descriptors.get("orientation"));
    println("orientedness  = " + descriptors.get("orientedness"));
    println("orientDot     = " + descriptors.get("orientDot")); 
    println("vorientation  = " + descriptors.get("vorientation"));
    println("vorientedness = " + descriptors.get("vorientedness"));
    println("meanVelocity  = " + descriptors.get("meanVelocity"));
    println("stdvVelocity  = " + descriptors.get("stdvVelocity"));
    println("nIntersects   = " + descriptors.get("nIntersects"));
    println("nCorners      = " + descriptors.get("nCorners"));
    println("compactness   = " + descriptors.get("compactness"));
    println("hullPointPct  = " + descriptors.get("hullPointPct"));

    for (int m=3; m<=16; m++) {
      String momentName = "moment" + nf(m, 2);
      println(momentName + "      = " + descriptors.get(momentName));
    }

    println("--------------");
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
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke)strokes.get(i);
      s.render();
    }
  }

  //---------------------------------
  void makeFlatCopyStroke () {
    // Make a single Stroke that contains all of the points. 
    // (Just make a simple copy if there was only one stroke to begin with.)

    flatCopyStroke.pts.clear();
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke srcStroke = (GMLStroke)strokes.get(i);
      ArrayList srcPts = srcStroke.pts;

      int nPoints = srcPts.size();
      for (int j=0; j<nPoints; j++) {
        XYTPoint srcPt = (XYTPoint) srcPts.get(j);
        flatCopyStroke.addXYTPoint(srcPt);
      }
    }
  }

  //---------------------------------
  void makeBlurredCopy() {

    int nStrokes = strokes.size();
    strokesSmoothed = new ArrayList(nStrokes);
    for (int i=0; i<nStrokes; i++) {
      GMLStroke srcStroke = (GMLStroke)strokes.get(i);
      GMLStroke dstStroke = new GMLStroke();

      ArrayList srcPts = srcStroke.pts;
      int nPoints = srcPts.size();
      for (int j=0; j<nPoints; j++) {
        XYTPoint srcPt = (XYTPoint) srcPts.get(j);
        XYTPoint dstPt = new XYTPoint(srcPt);
        dstStroke.addXYTPoint(dstPt);
      }
      strokesSmoothed.add(dstStroke);
    }

    int nRuns = 10;
    for (int n=0; n<nRuns; n++) {
      float A = 0.30;
      float B = 1.0 - 2.0*A;
      for (int i=0; i<nStrokes; i++) {
        GMLStroke aStroke = (GMLStroke)strokesSmoothed.get(i);

        ArrayList aStrokePts = aStroke.pts;
        int nPoints = aStrokePts.size();
        if (nPoints >= 7) {

          for (int j=1; j<(nPoints-1); j++) {
            XYTPoint pt0 = (XYTPoint) aStrokePts.get(j-1);
            XYTPoint pt1 = (XYTPoint) aStrokePts.get(j  );
            XYTPoint pt2 = (XYTPoint) aStrokePts.get(j+1);

            float x = (pt0.x + pt1.x + pt2.x)/3.0;
            float y = (pt0.y + pt1.y + pt2.y)/3.0;
            ((XYTPoint) aStrokePts.get(j)).set(x, y);
          }
        }
      }
    }
  }

  //---------------------------------
  void renderBlurredCopy () {
    int nStrokes = strokesSmoothed.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke)strokesSmoothed.get(i);
      s.render();
    }
  }


  //---------------------------------
  void renderVelocities () {
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke)strokes.get(i);

      ArrayList pts = s.pts;
      int nPoints = pts.size();

      for (int j=1; j<nPoints; j++) {
        XYTPoint pt0 = (XYTPoint) pts.get(j-1);
        XYTPoint pt1 = (XYTPoint) pts.get(j  );
        float vx = 10*(pt1.x - pt0.x);
        float vy = 10*(pt1.y - pt0.y);
        ellipse(0.5+vx, 0.5+vy, 0.01, 0.01);
      }
    }
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
    
    println("normalizingRange= " + normalizingRange + " lesserRange" + lesserRange);

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




  void  resampleTag ( Vector path, XYTPoint resampledPath[], int nResamples, float totalPathLength) {

    // takes a vector of Vec3f's, 
    // an array for nResampledPoints Vec3f's,
    // and the int nResampledPoints, which is the number of resamples

    XYTPoint	lower, upper;
    int		nResampledPoints= nResamples;
    int 	nPathPoints = path.size();
    float	RSL = totalPathLength / (float) nResampledPoints; // "resampledSegmentLength"
    float	Dx, Dy;
    float	dx, dy;
    float       RSLdx, RSLdy; 
    float	px, py;
    float	segLength, ASL; // available segment length
    float	remainder;

    float	prevRemainder=RSL;
    int 	nsegs;
    int 	p=0;
    int 	i;

    if (nPathPoints <= 1) { // special case for one-point path
      for (p=0; p<nResampledPoints; p++) {
        lower = (XYTPoint) path.elementAt(0);
        px = lower.x+ (float)p*0.0001f;
        py = lower.y+ (float)p*0.0001f;
        resampledPath[p].set(px, py, 0);
      }
    } 
    else {

      float neededSpace;
      for (i=0; i< nPathPoints-1; i++) {
        lower = (XYTPoint) path.elementAt(i);
        upper = (XYTPoint) path.elementAt(i+1);
        Dx = upper.x - lower.x;
        Dy = upper.y - lower.y;

        segLength = sqrt (Dx*Dx + Dy*Dy);
        ASL = segLength; // available segment length
        dx = (Dx/segLength); // unit vector components
        dy = (Dy/segLength);
        RSLdx = dx * RSL; // resampled segment vector components
        RSLdy = dy * RSL;

        neededSpace = RSL-prevRemainder;
        if (ASL >= neededSpace) { 
          // if there is enough room to place the first point
          // then place the first resample point in the latest segment
          remainder = ASL;
          px = lower.x + (neededSpace*dx);
          py = lower.y + (neededSpace*dy);

          if (p < nResampledPoints) { 
            resampledPath[p].set(px, py, 0);
            remainder -= neededSpace;
            p++;
          }

          int nPtsToDo = (int)(remainder/RSL);
          for (int d=0; d<nPtsToDo; d++) {

            px += RSLdx; 
            py += RSLdy;
            if (p < nResampledPoints) { 
              resampledPath[p].set(px, py, 0);
              remainder -= RSL;
              p++;
            }
          }

          prevRemainder = remainder;
        } 
        else { 
          // if there is not enough room to place the first point 
          prevRemainder += ASL;
        }
      }
    }
  }
}




//=====================================================================
class GMLAnalyzer {

  XYTPoint centroid;
  XYTPoint tmpPoint;
  float matrix2x2[][];
  float multiPartData[];
  float moments[];

  //---------------------------------
  GMLAnalyzer () {
    centroid      = new XYTPoint();
    tmpPoint      = new XYTPoint();
    matrix2x2     = new float[2][2];
    multiPartData = new float[8];


    moments       = new float[19];
  }

  //---------------------------------
  float getNothing (GMLTag tag) {
    if (tag != null) {
      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();

      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          float x = pt.x;
          float y = pt.y;
        }
      }
    }
    return 0.0;
  }

  //-----------------------------------------------------------------------------------------------------------------
  float getNothingFromFlatCopy (GMLTag tag) {

    if (tag != null) {
      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList pts         = aStroke.pts;
      int nPoints           = pts.size();

      for (int j=0; j<nPoints; j++) {
        XYTPoint pt = (XYTPoint) pts.get(j);
        float x = pt.x;
        float y = pt.y;
      }
    }
    return 0.0;
  }


  //---------------------------------
  float getNStrokes (GMLTag tag) {
    float out = 0; 
    if (tag != null) {
      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();
      out = nStrokes;
    }
    return out;
  }

  //---------------------------------
  float getAspectRatio (GMLTag tag) {
    float out = 0.0;
    if (tag != null) {
      GMLStroke aStroke = tag.flatCopyStroke;
      ArrayList pts = aStroke.pts;
      int nPoints = pts.size();
      if (nPoints > 2) {

        float minX =  99999;
        float maxX = -99999;
        float minY =  99999;
        float maxY = -99999;

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          float x = pt.x;
          float y = pt.y;

          if (x > maxX) {
            maxX = x;
          }
          if (x < minX) {
            minX = x;
          }
          if (y > maxY) {
            maxY = y;
          }
          if (y < minY) {
            minY = y;
          }
        }

        float dx = maxX - minX;
        float dy = maxY - minY;
        if (dx > 0) {
          out = dy / dx;
        }
      }
    }
    return out;
  }

  //---------------------------------
  float getOrientationDot (float orientation, float vorientation) {
    return (cos(orientation - vorientation));
  }


  //---------------------------------
  float[] getOrientation (GMLTag tag,  XYTPoint COM) {

    float orientation  = 0.0; 
    float orientedness = 0.0;

    if (tag != null) {

      GMLStroke aStroke = tag.flatCopyStroke;
      ArrayList pts = aStroke.pts;
      int nPoints = pts.size();
      if (nPoints > 2) {

        //arguments: an array of pixels, the array's width & height, and the location of the center of mass (com).
        //this function calculates the elements of a point set's tensor matrix,
        //calls the function calcEigenvector() to get the best eigenvector of this matrix
        //and returns this eigenVector as a pair of doubles

        //first we look at all the pixels, determine which ones contribute mass (the black ones),
        // and accumulate the sums for the tensor matrix
        float dX, dY; 
        float XXsum, YYsum, XYsum;

        XXsum = 0; 
        YYsum = 0; 
        XYsum = 0; 

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          dX = pt.x - COM.x;
          dY = pt.y - COM.y;
          XXsum += dX * dX;
          YYsum += dY * dY;
          XYsum += dX * dY;
        }

        // here's the tensor matrix
        matrix2x2[0][0] =  YYsum;
        matrix2x2[0][1] = -XYsum;
        matrix2x2[1][0] = -XYsum;
        matrix2x2[1][1] =  XXsum;

        // get the orientation of the bounding box
        float[] response = calcEigenvector ( matrix2x2 );
        orientation  = response[0];
        orientedness = response[1];
      }
    }

    multiPartData[0] = (orientation);
    multiPartData[1] = (orientedness);
    return multiPartData;
  }

  //---------------------------------
  float[] calcEigenvector ( float[][] matrix ) {
    //this function takes a 2x2 matrix, and returns a pair of angles which are the eigenvectors
    float A = matrix[0][0]; 
    float B = matrix[0][1];
    float C = matrix[1][0];
    float D = matrix[1][1];

    //because we assume a 2x2 matrix,
    //we can solve explicitly for the eigenValues using the Quadratic formula.
    //the eigenvalues are the roots of the equation  det( lambda * I  - T) = 0
    float a, b, c, root1, root2;
    a = 1.0;
    b = (0.0 - A) - D;
    c = (A * D) - (B * C);
    float Q = (b * b) - (4.0 * a * c);
    if (Q >= 0) {
      root1 = ((0.0 - b) + sqrt ( Q)) / (2.0 * a);
      root2 = ((0.0 - b) - sqrt ( Q)) / (2.0 * a);

      //assume x1 and x2 are the elements of the eigenvector.  Then, because Ax1 + Bx2 = lambda * x1, 
      //we know that x2 = x1 * (lambda - A) / B.
      float factor2 = ( min (root1, root2) - A) / B;

      //we arbitrarily set x1 = 1.0 and compute the magnitude of the eigenVector with respect to this assumption
      float magnitude2 = sqrt (1.0 + factor2*factor2);

      //we now find the exact components of the eigenVector by scaling by 1/magnitude
      if ((magnitude2 == 0) || (Float.isNaN(magnitude2))) {
        multiPartData[0] = 0;
        multiPartData[1] = 0;
      } 
      else {
        float orientedBoxOrientation = atan2 ( (1.0 / magnitude2), (factor2 / magnitude2));
        float orientedBoxEigenvalue  = log (1.0+root2); // orientedness
        multiPartData[0] = orientedBoxOrientation;
        multiPartData[1] = orientedBoxEigenvalue;
      }
    } 
    else {
      multiPartData[0] = 0;
      multiPartData[1] = 0;
    }

    return multiPartData;
  }


  //---------------------------------
  float[] getVelocityOrientation (GMLTag tag) {

    float velOrientation  = 0.0; 
    float velOrientedness = 0.0;
    float meanVelocity    = 0.0;
    float stDevVelocity   = 0.0; 

    if (tag != null) {

      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();

      XYTPoint velMean = tmpPoint;
      velMean.x = 0; 
      velMean.y = 0;
      float vx, vy, vh;
      int count = 0;
      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();
        for (int j=1; j<nPoints; j++) {
          XYTPoint pt0 = (XYTPoint) pts.get(j-1);
          XYTPoint pt1 = (XYTPoint) pts.get(j  );
          vx = pt1.x - pt0.x;
          vy = pt1.y - pt0.y;
          velMean.x += vx;
          velMean.y += vy;
          meanVelocity += sqrt(vx*vx + vy*vy);
          count++;
        }
      }
      if (count > 0) {
        velMean.x    /= (float) count;
        velMean.y    /= (float) count;
        meanVelocity /= (float) count;
      }


      float dvX, dvY; 
      float XXsum = 0; 
      float YYsum = 0; 
      float XYsum = 0;
      count = 0;  
      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        for (int j=1; j<nPoints; j++) {
          XYTPoint pt0 = (XYTPoint) pts.get(j-1);
          XYTPoint pt1 = (XYTPoint) pts.get(j  );

          vx = pt1.x - pt0.x;
          vy = pt1.y - pt0.y;
          dvX = vx - velMean.x;
          dvY = vy - velMean.y;

          vh = sqrt(vx*vx + vy*vy);
          float difFromMean = vh - meanVelocity;
          stDevVelocity += sq(difFromMean); 
          count++;

          XXsum += dvX * dvX;
          YYsum += dvY * dvY;
          XYsum += dvX * dvY;
        }
      }

      if (count > 0) { 
        stDevVelocity /= (float)count;
        stDevVelocity = sqrt(stDevVelocity);
      }

      // here's the tensor matrix
      matrix2x2[0][0] =  YYsum;
      matrix2x2[0][1] = -XYsum;
      matrix2x2[1][0] = -XYsum;
      matrix2x2[1][1] =  XXsum;

      // get the orientation of the bounding box
      float[] response = calcEigenvector ( matrix2x2 );
      velOrientation  = response[0];
      velOrientedness = response[1];
    }

    multiPartData[0] =         velOrientation;
    multiPartData[1] = 10.0  * velOrientedness;
    multiPartData[2] = 100.0 * meanVelocity;
    multiPartData[3] = 100.0 * stDevVelocity;

    return multiPartData;
  }

  //---------------------------------
  float getDuration (GMLTag tag) {

    float duration = 0; 

    if (tag != null) {
      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();
      float minT =  99999;
      float maxT = -99999;

      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          float t = pt.t;
          if (t > maxT) {
            maxT = t;
          }
          if (t < minT) {
            minT = t;
          }
        }
      }
      duration = maxT - minT;
    }
    return duration;
  }





  //--------------------------------- 
  float getMeanFVelocity (GMLTag tag) {
    float meanVel = 0; 

    if (tag != null) {
      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();
      int count = 0; 

      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        if (nPoints > 1) {
          for (int j=1; j<nPoints; j++) {
            XYTPoint pt0 = (XYTPoint) pts.get(j-1);
            XYTPoint pt1 = (XYTPoint) pts.get(j  );

            float dx = pt1.x - pt0.x;
            float dy = pt1.y - pt0.y;

            float dh = sqrt(dx*dx + dy*dy);
            meanVel += dh; 
            count++;
          }
        }
      }

      if (count > 0) { 
        meanVel /= (float)count;
      }
    }
    return meanVel;
  }


  //--------------------------------- 
  float getStDevFVelocity (GMLTag tag, float meanVel) {
    float stDevVel = 0; 

    if (tag != null) {
      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();
      int count = 0; 

      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        if (nPoints > 1) {
          for (int j=1; j<nPoints; j++) {
            XYTPoint pt0 = (XYTPoint) pts.get(j-1);
            XYTPoint pt1 = (XYTPoint) pts.get(j  );

            float dx = pt1.x - pt0.x;
            float dy = pt1.y - pt0.y;

            float dh = sqrt(dx*dx + dy*dy);
            float difFromMean = dh - meanVel;
            stDevVel += sq(difFromMean); 
            count++;
          }
        }
      }

      if (count > 0) { 
        stDevVel /= (float)count;
        stDevVel = sqrt(stDevVel);
      }
    }
    return stDevVel;
  }





  //---------------------------------
  XYTPoint getCentroid (GMLTag tag) {
    float cx = 0; 
    float cy = 0; 
    float ct = 0; 
    int nCentroidPts = 0; 

    if (tag != null) {

      ArrayList strokes = tag.strokes;
      int nStrokes = strokes.size();

      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPoints = pts.size();

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          cx += pt.x;
          cy += pt.y;
          ct += pt.t;
          nCentroidPts++;
        }
      }

      if (nCentroidPts > 0) {
        cx /= (float)nCentroidPts;
        cy /= (float)nCentroidPts;
        ct /= (float)nCentroidPts;
      }
    }

    centroid.set(cx, cy, ct); 
    return centroid;
  }



  //-----------------------------------------------------------------------------------------------------------------
  float getHullPointPercentage (GMLTag tag, ArrayList hull) {

    float out = 0;
    if ((tag != null) && (hull != null)) {
      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList strokePts   = aStroke.pts;
      int nStrokePoints     = strokePts.size();
      int nHullPoints       = hull.size() - 1;

      if ((nStrokePoints > 3) && (nHullPoints > 3)) {
        out = (float) nHullPoints / (float) nStrokePoints;
        out = min(1, max(0, out)); 
      }
    }
    return out;
  }

  //---------------------------------
  float getPointDensity (GMLTag tag, ArrayList hull) {
    float out = 0; 
    if ((tag != null) && (hull != null)) {
      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList strokePts   = aStroke.pts;
      int nStrokePoints     = strokePts.size();
      if (nStrokePoints > 3) {
        float hullArea = getHullArea( hull); 
        if (hullArea > 0) {
          out = (float) nStrokePoints / hullArea;
        }
      }
    }
    return out;
  }

  //---------------------------------
  float getCompactness (float arcLength, ArrayList hull) {
    float out = 0; 
    if (hull != null) {

      float hullArea = getHullArea( hull); 
      if ((hullArea > 0) && (arcLength > 0)) {
        out = (arcLength*arcLength) /hullArea;
        out = sqrt(out);
      }
    }
    return out;
  }

  //---------------------------------
  // return area of polygon
  float getHullArea( ArrayList hull) {
    float out = 0; 
    if (hull != null) {
      float sum = 0.0;
      int nHullPoints = hull.size();
      for (int i = 0; i < (nHullPoints-1); i++) {
        XYTPoint hullpt0 = (XYTPoint) hull.get(i); 
        XYTPoint hullpt1 = (XYTPoint) hull.get(i+1); 
        sum += (hullpt0.x * hullpt1.y) - (hullpt0.y * hullpt1.x);
      }
      out = abs(0.5 * sum);
    }
    return out;
  }

  //---------------------------------
  float getMeanDistanceFromCentroid  (GMLTag tag, XYTPoint cent) {
    float meanDistance = 0.0;

    if (tag != null) {

      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList pts         = aStroke.pts;
      int nPoints           = pts.size();

      if (nPoints > 2) {
        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          float dx = pt.x - cent.x;
          float dy = pt.y - cent.y;
          float dh = sqrt(dx*dx + dy*dy); 
          meanDistance += dh;
        }
        meanDistance /= (float)nPoints;
      }
    }
    return meanDistance;
  }


  //---------------------------------
  float getStDevDistanceFromCentroid  (GMLTag tag, XYTPoint cent, float meanDist) {
    float stDevDistance = 0;

    if (tag != null) {

      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList pts         = aStroke.pts;
      int nPoints           = pts.size();
      if (nPoints > 2) {

        for (int j=0; j<nPoints; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          float dx = pt.x - cent.x;
          float dy = pt.y - cent.y;
          float dh = sqrt(dx*dx + dy*dy); 

          float difFromMean = dh - meanDist;
          stDevDistance += sq(difFromMean);
        }

        stDevDistance /= (float)nPoints;
        stDevDistance  = sqrt(stDevDistance);
      }
    } 
    return stDevDistance;
  }


  //---------------------------------
  float getArcLength (GMLTag tag) {
    float tagLength = 0.0; 

    if (tag != null) {

      GMLStroke aStroke     = tag.flatCopyStroke;
      ArrayList pts         = aStroke.pts;
      int nPoints           = pts.size();

      if (nPoints > 1) {
        for (int j=1; j<nPoints; j++) {
          XYTPoint pt0 = (XYTPoint) pts.get(j-1);
          XYTPoint pt1 = (XYTPoint) pts.get(j  );
          float dx = pt1.x - pt0.x;
          float dy = pt1.y - pt0.y;
          float dh = sqrt(dx*dx + dy*dy);
          tagLength += dh;
        }
      }
    }
    return tagLength;
  }



  //---------------------------------
  float[] getTotalAngle (GMLTag tag) {

    float totalAngle         = 0.0; // the cumulative angle through which the mark moves.
    float totalAbsoluteAngle = 0.0;
    float meanAngle          = 0.0;
    float meanAbsAngle       = 0.0;
    float nCorners           = 0.0;

    if (tag != null) {
      int count = 0;
      float cornerThreshold = radians(30.0); 

      ArrayList strokes = tag.strokesSmoothed;
      int nStrokes = strokes.size();
      if (nStrokes > 0) {

        float x0 = 0;
        float x1 = 0;
        float x2 = 0;

        float y0 = 0;
        float y1 = 0;
        float y2 = 0;

        float x1b,y1b;

        for (int i=0; i<nStrokes; i++) {
          GMLStroke s = (GMLStroke) strokes.get(i);
          ArrayList pts = s.pts;
          int nPoints = pts.size();

          for (int j=0; j<nPoints; j++) {
            XYTPoint pt = (XYTPoint) pts.get(j);
            x0 = x1;
            x1 = x2;
            x2 = pt.x;
            y0 = y1;
            y1 = y2;
            y2 = pt.y;

            float theta = getJointAngle (x0, y0, x1, y1, x2, y2);
            float absTheta     = abs(theta);
            totalAngle         +=    theta;
            totalAbsoluteAngle += absTheta;

            if ((absTheta > cornerThreshold) && (j>1) && (j<(nPoints-1))) {
              nCorners++;
            }
            count++;
          }
        }
      }
      if (count > 2) {      
        meanAngle    = totalAngle / (float) count;
        meanAbsAngle = totalAbsoluteAngle / (float) count;
      } 
      else {
        totalAngle = totalAbsoluteAngle = 0;
      }
    }

    multiPartData[0] =  (totalAngle);
    multiPartData[1] =  (totalAbsoluteAngle);
    multiPartData[2] =  100.0 * meanAngle;
    multiPartData[3] =  10.0  * meanAbsAngle;
    multiPartData[4] =   nCorners;
    return multiPartData;
  }


  //---------------------------------
  float getStDevAbsAngle (GMLTag tag, float meanAbsAngle) {
    float stDevAbsAngle = 0;

    if (tag != null) {

      ArrayList strokes = tag.strokesSmoothed;
      int nStrokes = strokes.size();
      int count = 0; 
      if (nStrokes > 0) {

        float x0 = 0;
        float x1 = 0;
        float x2 = 0;
        float y0 = 0;
        float y1 = 0;
        float y2 = 0;

        float x1b,y1b;

        for (int i=0; i<nStrokes; i++) {
          GMLStroke s = (GMLStroke) strokes.get(i);
          ArrayList pts = s.pts;
          int nPoints = pts.size();

          for (int j=0; j<nPoints; j++) {
            XYTPoint pt = (XYTPoint) pts.get(j);
            x0 = x1;
            x1 = x2;
            x2 = pt.x;
            y0 = y1;
            y1 = y2;
            y2 = pt.y;

            float absTheta = abs(getJointAngle (x0, y0, x1, y1, x2, y2));
            float difFromMean = absTheta - meanAbsAngle;
            stDevAbsAngle += sq(difFromMean); 
            count++;
          }
        }

        if (count > 2) { 
          stDevAbsAngle /= (float)count;
          stDevAbsAngle = sqrt(stDevAbsAngle);
        } 
        else {
          stDevAbsAngle = 0;
        }
      }
    }
    return stDevAbsAngle;
  }


  //---------------------------------
  float getJointAngle (float x0, float y0, float x1, float y1, float x2, float y2) {
    // get the angle between two segments. 

    float anglei = 0;
    float dxBA = x1 - x0 + 0.00001;
    float dyBA = y1 - y0;
    float dhBA = (dxBA*dxBA + dyBA*dyBA);
    float dxCB = x2 - x1 + 0.00001;
    float dyCB = y2 - y1;
    float dhCB = (dxCB*dxCB + dyCB*dyCB);

    if ((dhBA > 0) && (dhCB > 0)) {
      float slopeCB = dyCB /  dxCB;
      float slopeBA = dyBA /  dxBA;
      float angleBA = atan(slopeBA);
      float angleCB = atan(slopeCB);


      if (dxBA > 0) {
        if (dyBA > 0) {
          if (dxCB > 0) {
            if (dyCB > 0) {
              anglei = angleBA + (PI - angleCB);
            }  
            else if (dyCB <= 0) { 
              anglei = PI + angleBA - angleCB   ;
            }
          }
          else if (dxCB < 0) {
            if (dyCB > 0) {
              anglei = angleBA - angleCB;
            }  
            else if (dyCB <= 0) {
              float dif = 	 angleBA - angleCB;
              if (dif > 0) { 
                anglei = dif;
              }
              else {
                anglei = TWO_PI + dif;
              }
            }
          }
        }
        else if (dyBA <= 0) {
          if (dxCB > 0) {
            if (dyCB > 0) {
              anglei =  PI + angleBA - angleCB;
            } 				 
            else if (dyCB <= 0) {
              anglei = PI + angleBA - angleCB;
            }
          }
          else if (dxCB < 0) {
            if (dyCB > 0) {
              float dif = (0 - angleBA) + angleCB;
              if (dif > 0) {
                anglei = TWO_PI - dif;
              }
              else { 
                anglei = 0 - dif;
              }
            }
            else if (dyCB <= 0) {
              anglei = TWO_PI + angleBA  - angleCB;
            }
          }
        }
      }
      else if (dxBA < 0) {
        if (dyBA <= 0) {
          if (dxCB < 0) {
            if (dyCB <= 0) { 
              anglei = angleBA + PI - angleCB;
            }
            else if (dyCB > 0) { 
              anglei = PI + angleBA - angleCB ;
            }
          }
          else if (dxCB > 0) {
            if (dyCB <= 0) { 
              anglei = angleBA - angleCB;
            }
            else if (dyCB > 0) { 
              float dif = angleBA - angleCB;
              if (dif > 0) {
                anglei = dif;
              }
              else { 
                anglei = TWO_PI + dif;
              }
            }
          }
        }
        else if (dyBA > 0) {
          if (dxCB < 0) {
            if (dyCB > 0) { 
              anglei = PI + angleBA - angleCB;
            }
            else if (dyCB <= 0) { 
              anglei = PI + angleBA - angleCB;
            }
          }
          else if (dxCB > 0) {
            if (dyCB > 0) { 
              anglei = TWO_PI - angleCB + angleBA;
            }
            else if (dyCB <= 0) {
              float dif = 	(0 - angleBA) + angleCB;
              if ( dif > 0) {
                anglei = TWO_PI - dif;
              }
              else {
                anglei = 0 - dif;
              }
            }
          }
        }
      }
    }

    anglei = (PI - anglei);
    return anglei;
  }









  //---------//
  // Moments //
  // Adapted from http://kenai.com/projects/audiveris/sources/hg/content/src/main/omr/math/Moments.java?rev=3011
  //---------//

  /**
   
   * Compute the moments for a set of points whose x and y coordinates are
   * provided, all values being normed by the provided unit value.
   * @param x    the array of abscissa values
   * @param y    the array of ordinate values
   * @param dim  the number of points
   * @param unit the length (number of pixels, for example 20) of norming unit
   */

  /*
 "weight", // 0
   "width", // 1 
   "height", //  2
   "n20", // 3
   "n11", // 4
   "n02", // 5
   "n30", // 6
   "n21", // 7
   "n12", // 8
   "n03", // 9
   "h1", // 10
   "h2", // 11
   "h3", // 12
   "h4", // 13
   "h5", // 14
   "h6", // 15
   "h7", // 16
   "xBar", // 17
   "yBar" // 18
   */


  float[] getMoments (GMLTag tag) { 

    for (int i=0; i<19; i++) {
      moments[i] = 0;
    }

    // how many points altogether?
    int nPoints = 0; 
    ArrayList strokes = tag.strokes;
    int nStrokes = strokes.size();
    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke) strokes.get(i);
      ArrayList pts = s.pts;
      int nPointsInStroke = pts.size();
      nPoints += nPointsInStroke;
    }

    // create arrays to store points. 
    if (nPoints > 3) {
      double x[] = new double[nPoints];
      double y[] = new double[nPoints];

      int count = 0; 
      for (int i=0; i<nStrokes; i++) {
        GMLStroke s = (GMLStroke) strokes.get(i);
        ArrayList pts = s.pts;
        int nPointsInStroke = pts.size();

        for (int j=0; j<nPointsInStroke; j++) {
          XYTPoint pt = (XYTPoint) pts.get(j);
          x[count] = (double) pt.x;
          y[count] = (double) pt.y;
          count++;
        }
      }


      double dx;
      double dy;

      // Normalized Moments
      double unit = 1.0; 
      double n00 = (double) nPoints / (double) (unit * unit);
      double n01 = 0.0;
      double n02 = 0.0;
      double n03 = 0.0;
      double n10 = 0.0;
      double n11 = 0.0;
      double n12 = 0.0;
      double n20 = 0.0;
      double n21 = 0.0;
      double n30 = 0.0;

      // Total weight
      double w = nPoints; // For p+q = 0
      double w2 = w * w;  // For p+q = 2
      double w3 = Math.sqrt(w * w * w * w * w); // For p+q = 3

      // Mean x & y
      for (int i = nPoints-1; i >= 0; i--) {
        n10 += x[i];
        n01 += y[i];
      }
      n10 /= (double) nPoints;
      n01 /= (double) nPoints;


      // width & height
      double  xMin = Double.MAX_VALUE;
      double  xMax = Double.MIN_VALUE;
      double  yMin = Double.MAX_VALUE;
      double  yMax = Double.MIN_VALUE;
      for (int i = nPoints-1; i >= 0; i--) {
        double xx = x[i];
        if (xx < xMin) {
          xMin = xx;
        }
        if (xx > xMax) {
          xMax = xx;
        }

        double yy = y[i];
        if (yy < yMin) {
          yMin = yy;
        }
        if (yy > yMax) {
          yMax = yy;
        }
      }

      for (int i = nPoints-1; i >= 0; i--) {
        dx = x[i] - n10;
        dy = y[i] - n01;

        n11 += (dx * dy);
        n12 += (dx * dy * dy);
        n21 += (dx * dx * dy);
        n20 += (dx * dx);
        n02 += (dy * dy);
        n30 += (dx * dx * dx);
        n03 += (dy * dy * dy);
        
        //println(i + "dx3 = " + (dx * dx * dx));
      }

      // Normalize
      // p + q = 2
      n11 /= w2;
      n20 /= w2;
      n02 /= w2;
      
      // p + q = 3
      n12 /= w3;
      n21 /= w3;
      n30 /= w3;
      n03 /= w3;

      // Assign non-orthogonal centralized moments
      // (invariant to translation & scaling)
      moments[0] = (float) n00; // Weight
      moments[1] = (float) ((xMax - xMin) / unit); // Width
      moments[2] = (float) ((yMax - yMin) / unit); // Height
      moments[3] = (float) n20; // X  absolute eccentricity
      moments[4] = (float) n11; // XY covariance
      moments[5] = (float) n02; // Y  absolute eccentricity
      moments[6] = (float) n30; // X  signed eccentricity
      moments[7] = (float) n21; //
      moments[8] = (float) n12; //
      moments[9] = (float) n03; // Y signed eccentricity


      // Assign orthogonals moments (Hu set)
      // (Invariant to translation / scaling / rotation)
      moments[10] = (float) (n20 + n02);
      moments[11] = (float) (((n20 - n02) * (n20 - n02)) + (4 * n11 * n11));
      moments[12] = (float) (((n30 - (3 * n12)) * (n30 - (3 * n12))) + ((n03 - (3 * n21)) * (n03 - (3 * n21))));
      moments[13] = (float) (((n30 + n12) * (n30 + n12)) + ((n03 + n21) * (n03 + n21)));
      moments[14] = (float) (((n30 - (3 * n12)) * (n30 + n12) * (((n30 + n12) * (n30 + n12)) - (3 * (n21 + n03) * (n21 + n03)))) + 
        ((n03 - (3 * n21)) * (n03 + n21) * (((n03 + n21) * (n03 + n21)) - (3 * (n12 + n30) * (n12 + n30)))));
      moments[15] = (float) (((n20 - n02) * (((n30 + n12) * (n30 + n12)) - ((n03 + n21) * (n03 + n21)))) + (4 * n11 * (n30 + n12) * (n03 + n21)));
      moments[16] = (float) ((((3 * n21) - n03) * (n30 + n12) * (((n30 + n12) * (n30 + n12)) - (3 * (n21 + n03) * (n21 + n03)))) - 
        (((3 * n12) - n30) * (n03 + n21) * (((n03 + n21) * (n03 + n21)) - (3 * (n12 + n30) * (n12 + n30)))));

      // Mass center placed here
      moments[17] = (float) (n10); // xBar
      moments[18] = (float) (n01); // yBar
      
    }

    return moments;
  }
}


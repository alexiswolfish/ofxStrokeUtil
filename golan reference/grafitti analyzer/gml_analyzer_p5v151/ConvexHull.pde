

class ConvexHull {

  ArrayList tempHull;
  ArrayList hull;
  int depth; 

  public ConvexHull() {
    tempHull    = new ArrayList();
    hull        = new ArrayList();
  }





  //-----------------------------------------------------------------------------------------------------------------
  // QuickHull Algorithm implementation.
  public ArrayList quickHull (GMLTag tag) {

    GMLStroke aStroke   = tag.flatCopyStroke;
    ArrayList strokePts = aStroke.pts;
    int numberOfPoints  = strokePts.size();
    depth = 0; 


    if (numberOfPoints > 3) {

      hull.clear();
      tempHull.clear();

      ArrayList P1 = new ArrayList();
      ArrayList P2 = new ArrayList();

      XYTPoint l = (XYTPoint) strokePts.get(0);
      XYTPoint r = (XYTPoint) strokePts.get(0);

      float minX = l.x;
      float maxX = l.x;
      int minAt = 0;
      int maxAt = 0;	


      // find the max and min x-coord point
      XYTPoint currPt;
      for (int i = 1; i < numberOfPoints; i++) {
        currPt = (XYTPoint) strokePts.get(i);
        if ( currPt.x > maxX) {
          r = currPt;
          maxX = currPt.x;
          maxAt = i;
        };
        if ( currPt.x < minX) {
          l = currPt;
          minX = currPt.x;
          minAt = i;
        };
      }


      //find out whether each point is over or under the line formed by the two points 
      // with min and max x-coord, and put them in 2 group according to whether they are above or under 
      for (int i = 0; i < numberOfPoints; i++) {
        if ((i != maxAt) && (i != minAt)) {
          currPt = (XYTPoint) strokePts.get(i);
          if (onLeft( l, r, currPt)) {
            P1.add( currPt );
          } 
          else {
            P2.add( currPt );
          }
        }
      };

      //put the max and min x-cord points in each group 
      P1.add( l);
      P1.add( r);
      P2.add( l );
      P2.add( r );

      //calculate the upper hull; put the upper hull result in final result  
      quick (P1, l, r, 0);
      for (int k=0; k<tempHull.size(); k++) {
        hull.add( tempHull.get(k));
      }	

      // important for proper shape
      hull.add(r);

      //critical! or else you get duplication
      tempHull.clear();

      //calculate the lower hull; append the result from lower hull to final result 
      quick(P2, l, r, 1);
      for (int k=(tempHull.size()-1); k>=0; k--) {
        hull.add( tempHull.get(k));
      }
    }
    System.gc();
    return hull;
  }


  //-----------------------------------------------------------------------------------------------------------------
  void renderHull (ArrayList hull) {

    if (hull != null) {
      int nHullPoints = hull.size();
      fill(0, 100, 0, 100); 

      beginShape();
      for (int j=0; j<nHullPoints; j++) {
        XYTPoint hullpt = (XYTPoint) hull.get(j); 
        vertex(hullpt.x, hullpt.y);
      }
      endShape();
    }
  }





  //-----------------------------------------------------------------------------------------------------------------
  boolean onLeft (XYTPoint aLinep0, XYTPoint aLinep1, XYTPoint chkpt) {
    //Given a Check point, determine if this check point is lying on the
    //left side or right side of the first point of the line.

    if (aLinep0.x == aLinep1.x ) {
      if (chkpt.x < aLinep0.x) return true;

      else {
        if (chkpt.x == aLinep0.x) {
          if (((chkpt.y > aLinep0.y) && (chkpt.y < aLinep1.y)) || ((chkpt.y > aLinep1.y) && (chkpt.y < aLinep0.y))) {
            return true;
          } 
          else {
            return false;
          }
        }
        else return false;
      }
    }

    else {            

      //stupid shit left over from dealing with integer points & rounding error, feh
      float slope = (aLinep1.y - aLinep0.y) / (aLinep1.x - aLinep0.x);
      float x3    = (((chkpt.x + slope * (slope * aLinep0.x - aLinep0.y + chkpt.y)) /  (1.0 + slope * slope)) * 10000.0 );
      float y3    = ((slope * (x3 / 10000 - aLinep0.x) + aLinep0.y) * 10000.0 );

      if (slope == 0.0) {
        if ((chkpt.y * 10000.0 ) > y3) return true; 
        else return false;
      }

      else { 
        if (slope > 0.0 ) {
          if (x3 > (chkpt.x * 10000.0 )) return true; 
          else return false;
        }

        else {
          if ((chkpt.x * 10000.0 ) > x3) return true; 
          else return false;
        }
      }
    }
  }




  //-----------------------------------------------------------------------------------------------------------------
  void quick (ArrayList P, XYTPoint l, XYTPoint r, int faceDir) {

    // Recursive method to find the Hull.
    // faceDir is 0 if we are calculating the upper hull.
    // faceDir is 1 if we are calculating the lower hull.

    int MAX_RECURSION_DEPTH = 1000;
    if ((P.size() == 2) || (depth > MAX_RECURSION_DEPTH)) {
      XYTPoint p0 = (XYTPoint) P.get(0);
      tempHull.add(p0); 
      P.clear();
      return;
    } 

    else {
      
      depth++;
      int hAt = splitAt(P, l, r);
      XYTPoint phAt = (XYTPoint) P.get(hAt);
      XYTPoint currPt;

      ArrayList P1 = new ArrayList();
      ArrayList P2 = new ArrayList();

      for (int i = 0; i < (P.size() - 2); i++) {
        if (i != hAt) {
          currPt = (XYTPoint) P.get(i);
          if (faceDir == 0) {
            if (onLeft(l, phAt, currPt)) {  
              P1.add(   currPt     );
            }
            if (onLeft(phAt, r, currPt)) {  
              P2.add(   currPt     );
            }
          } 

          else {
            if (!(onLeft(l, phAt, currPt))) {  
              P1.add(   currPt     );
            };
            if (!(onLeft(phAt, r, currPt))) { 
              P2.add(   currPt     );
            };
          };
        }
      }

      P1.add( l );  
      P1.add( phAt );
      P2.add( phAt );
      P2.add( r );

      if (faceDir == 0) {
        quick(P1, l, phAt, 0);
        quick(P2, phAt, r, 0);
      } 
      else {
        quick(P1, l, phAt, 1);
        quick(P2, phAt, r, 1);
      }

      return;
    }
  }



  //Find  a point which is certain to be in the Hull, from among a group of points
  //All the given points are on the same side of the line formed by l and r,
  //so the point with the longest distance perpendicular to this line is 
  //the point we are looking for.  Return the index of this point in the Vector.

  int splitAt (ArrayList P, XYTPoint l, XYTPoint r) {

    float denom;  
    denom = r.x - l.x;
    if (r.x == l.x) { 
      denom = 0.00001;
    }

    float newLnSlope = (r.y - l.y) / denom;
    float x1 = 0; 
    float y1 = 0;
    float x3 = 0; 
    float y3 = 0;

    float maxDist = 0; 
    float distance = 0;
    int farPtIndex = 0;

    for (int i = 0; i < (P.size() -2); i++) {

      XYTPoint aPoint = (XYTPoint) P.get(i);
      x1 = aPoint.x;
      y1 = aPoint.y;

      if (r.y == l.y) {
        x3 = x1;
        y3 = l.y;
      } 
      else {
        x3 =   (x1 + newLnSlope * (newLnSlope * l.x - l.y +  y1)) / (1.0 + newLnSlope * newLnSlope);
        y3 =   newLnSlope * (x3 - l.x) + l.y ;
      }

      distance = sqrt(pow((y1-y3), 2) + pow((x1-x3), 2));
      if (distance > maxDist) {
        maxDist = distance;
        farPtIndex = i;
      }
    }

    return farPtIndex;
  }
}


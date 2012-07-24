

//---------------------------------
int getSelfIntersectionCount (GMLTag tag) {
  
  int nSelfIntersections = 0; 
  if (tag != null) {
    ArrayList strokes = tag.strokes;
    int nStrokes = strokes.size();

    for (int i=0; i<nStrokes; i++) {
      GMLStroke s = (GMLStroke) strokes.get(i);
      ArrayList spts = s.pts;
      int snPoints = spts.size();
      for (int j=1; j<snPoints; j++) {
        XYTPoint spt0 = (XYTPoint) spts.get(j-1);
        XYTPoint spt1 = (XYTPoint) spts.get(j  );

        for (int k=0; k<nStrokes; k++) {
          GMLStroke t = (GMLStroke) strokes.get(k);
          ArrayList tpts = t.pts;
          int tnPoints = tpts.size();
          for (int l=1; l<tnPoints; l++) {
            XYTPoint tpt0 = (XYTPoint) tpts.get(l-1);
            XYTPoint tpt1 = (XYTPoint) tpts.get(l  );

            boolean bSelfIntersects = testTwoLinesIntersection (spt0, spt1, tpt0, tpt1);
            if (bSelfIntersects) { 
              nSelfIntersections++;
            }
          }
        }
      }
    }
  }
  nSelfIntersections /= 2;
  return nSelfIntersections;
}


//-----------------------------------------------------------------------------------------------------------------
boolean testTwoLinesIntersection( XYTPoint p1, XYTPoint p2, XYTPoint p3, XYTPoint p4) {

  //determine whether two line segments intersect
  boolean intersectsP = false;
  float denominator, numerator1, numerator2, alph, beta;
  	
  float Ax = p2.x - p1.x; 
  float Ay = p2.y - p1.y; 		 

  float Bx = p3.x - p4.x; 
  float By = p3.y - p4.y;

  float Cx = p1.x - p3.x; 
  float Cy = p1.y - p3.y;

  denominator = (Ay*Bx) - (Ax*By);
  numerator1  = (Ax*Cy) - (Ay*Cx);

  if (denominator !=0) {
    alph = numerator1 / denominator;
    if ((alph > 0.0) & (alph < 1.0)) {

      numerator2 = (By * Cx) - (Bx * Cy);
      beta = numerator2 / denominator;
      if ((beta > 0.0) && (beta < 1.0)) {
        intersectsP = true;
      }
    }
  }

  return intersectsP;
}


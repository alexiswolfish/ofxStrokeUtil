// CREATES and LOADS the SVD analysis file, in the format: 
//    #####.gml  coord0  coord1  coord2
//    #####.gml  coord0  coord1  coord2
// The creator will be replaced by the server-side python script. 
// The loader  will be replaced by the Processing.js client-side script. 





String SVDAnalysisFilename3D = "_SVDAnalysis3D.js";
String SVDAnalysisFilename2D = "_SVDAnalysis2D.js";
//=============================================================================
void performSvdAndSaveSvdAnalysisFile () {

  // Using the DatumList, 
  // Populate a raw data array (matrix), for processing by the SVD algorithm. 
  int nDatums = DatumList.size(); // count of the above
  int nValues = nTotalFields; // values per Datum
  matrixRawData = new double[nValues][nDatums];
  for (int d=0; d<nDatums; d++) {
    Datum aDatum = (Datum)DatumList.get(d);
    for (int v=0; v<nValues; v++) {
      matrixRawData[v][d] = (double) aDatum.getValue(v); // transposing happens here.
    }
  }


  // Compute the actual SVD.
  m = new Matrix (matrixRawData);
  long svdStartTime = millis();
  svd = new SingularValueDecomposition(m); // m.transpose()
  long svdEndTime = millis();
  println("The SVD took " + (svdEndTime - svdStartTime) + " millis.");


  // Extract the projected point data
  Matrix v                 = svd.getV();
  double[][] projectedData = v.getArray(); 
  for (int d=0; d<nDatums; d++) {
    double[] p = projectedData[d];
    ((Datum)DatumList.get(d)).setProjectedPoint(p, nProjectionDimensions);
  }

  // Write the projected point data to a JSON file. 
  PrintWriter output;
  String SVDAnalysisFilename = (nProjectionDimensions == 2) ? SVDAnalysisFilename2D : SVDAnalysisFilename3D;
  output = createWriter("data/" + SVDAnalysisFilename);
  output.println("[");
  
  for (int d=0; d<nDatums; d++) {
    String aLine = ""; 
    String filename =  ((Datum)DatumList.get(d)).name;
    aLine += "[" + filename; 
    float projectedPoint[] = ((Datum)DatumList.get(d)).getProjectedPoint();
    for (int i=0; i<nProjectionDimensions; i++) {
      aLine += "," + (projectedPoint[i]); // eliminate minus signs to save 3% of file size!
    }
    aLine += "]";
    if (d < (nDatums-1)) {
      aLine += ",";
    }
    output.println(aLine);
  }
  output.println("]");
  output.flush(); // Writes the remaining data to the file
  output.close(); // Finishes the file
  println("Wrote " + SVDAnalysisFilename);
}

//=============================================================================
void loadSvdAnalysisFile() {
  // Loads "SVDAnalysis.tsv" (eventually, from the server) 

  boolean bVerbose = false;
  String SVDAnalysisFilename = (nProjectionDimensions == 2) ? SVDAnalysisFilename2D : SVDAnalysisFilename3D;
  String svdAnalysisStrings[] = loadStrings(SVDAnalysisFilename); 
  int nSvdAnalysisStrings = svdAnalysisStrings.length;

  float tempFloats3[] = {
    0.0, 0.0, 0.0
  };

  // Reading in this file line-by-line,
  // Put all SVD data into the ArrayList, SVDData
  if (nSvdAnalysisStrings > 1) {
    SVDData = new ArrayList<SVDDatum3D>();
    for (int i=1; i<nSvdAnalysisStrings; i++) {

      String aLine = svdAnalysisStrings[i];
      String aLinePieces[] = split(aLine, ","); 
      int nLinePieces = aLinePieces.length;
      if (nLinePieces > 1) {

        boolean bSvdLineOK = true;

        String gmlFilename = ""; 
        gmlFilename = aLinePieces[0]; 
        int nGmlFilenameChars = gmlFilename.length(); 
        gmlFilename = gmlFilename.substring(1, nGmlFilenameChars);

        for (int j=1; j<(nLinePieces-1); j++) {

          float aVal = 0; 
          try {
            String aPiece = aLinePieces[j];
            if (aPiece.endsWith("]")) {
              aPiece = aPiece.substring(0, aPiece.length()-1);
            }
            //println(j + " " + aPiece); 
            aVal = (float)(Float.parseFloat (aPiece));
            tempFloats3[j-1] = aVal;
          } 
          catch (NumberFormatException e) {
            bSvdLineOK = false;
          }
        }
        float svdx = tempFloats3[0];
        float svdy = tempFloats3[1]; 
        float svdz = (nProjectionDimensions == 2) ? 0 : tempFloats3[2];  

        if (bSvdLineOK) {
          SVDDatum3D aSvdDatum = new SVDDatum3D(gmlFilename, svdx, svdy, svdz);
          SVDData.add(aSvdDatum); 
          if (bVerbose) {
            aSvdDatum.print();
          }
        }
      }
    }

    if (bVerbose) {
      println("SVDData has " + (SVDData.size()));
    }
  } 
  else {
    println(SVDAnalysisFilename + " is invalid");
  }
}


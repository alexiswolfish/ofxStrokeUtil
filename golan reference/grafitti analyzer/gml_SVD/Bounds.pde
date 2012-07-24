
// A class which stores the min and max for some range of data. 
class Bounds {
  String name = ""; 
  int    id;
  float  low ;
  float  high;
  float  mean;
  float  stdv;
  
  Bounds (){
    name = ""; 
    id = 0; 
    low = high = 0.0;
    mean = stdv = 0.0; 
  }
  
  void set(int i, float l, float h){ 
    id = i;
    low = l;
    high = h;
  }
  
  void set(int i, float l, float h, float m, float s){ 
    id = i;
    low = l;
    high = h;
    mean = m; 
    stdv = s; 
  }
  
  void print(){
    println("Bound "+ id + ":\t" + name +  "\t" + low + "\t" + high + "\t" + mean + "\t" + stdv); 
  }
}



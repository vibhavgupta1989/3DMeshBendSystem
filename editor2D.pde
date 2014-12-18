// //<>//
//**************************** 2D Graph Library Definitions ***********************
//

import java.io.*;

/*enum StrokeMode {
  Normal, 
  Radial, 
  Ball,
};*/

Point2D calcStrokeOffset(int mode, Point2D p, Point2D pd, float d, float dd, int dir) {
  //--
  Point2D O;  
  
  //--
  Vector2D T = new Vector2D(pd.x, pd.y);  
  Vector2D N = rotate90(T);
  Vector2D n = normalize(N);
  float k = normSq(T);

  switch (mode) {
  default:
  case 0:
    O = translate(p, n, dir * d); 
    break;
  case 1:
    O = translate(p, add(scale(T, -dd), scale(N, dir * sqrt(k - sq(dd)))), d / k);    
    break;
  }
  return O;
}

Point2D[] _generateStroke(Point2D[] vertices, float[] radii, int strokeMode, float scale, int dir) {
  //--
  int n = vertices.length;
  float delta = 1.0f / 50;
  Point2D[] stroke = new Point2D[n];
      
  //--
  for (int i = 0; i < n; ++i) {
    //--
    Point2D p1 = vertices[i];
    float r1 = radii[i];
    
    //--
    Point2D p0, p2;
    float r0, r2;
    if (i == 0) {
      p0 = p1;
      r0  = r1;
    } else {
      p0 = vertices[i-1];
      r0 = radii[i-1];
    }
    
    //--
    if (i == n-1) {
      p2 = p1;
      r2  = r1; 
    } else {
      p2 = vertices[i+1];
      r2 = radii[i+1];
    }
    
    //--
    Point2D pd = calcNevilleDerivative(p0, p1, p2, delta);
    float rd   = calcNevilleDerivative(r0, r1, r2, delta);
    Point2D O  = calcStrokeOffset(strokeMode, p1, pd, r1 * scale, rd * scale, dir);
    //println("i="+i+",r0"+r0+",r1="+r1+",r2="+r2+",rd="+rd);
    stroke[i] = O;
  }
  
  return stroke;
}

Point2D[] _generateStroke(Point2D[] vertices, float[] radii, int strokeMode) {
  //--
  Point2D[] l, r;
  if (strokeMode == 2) {
    l = _generateStroke(vertices, radii, 0, 0.5, 1);
    l = _generateStroke(l, radii, 1, 0.5, 1);
    r = _generateStroke(vertices, radii, 0, 0.5, -1);
    r = _generateStroke(r, radii, 1, 0.5, -1);
  } else {
    l = _generateStroke(vertices, radii, strokeMode, 1.0, 1);
    r = _generateStroke(vertices, radii, strokeMode, 1.0, -1);
  }
  
  //--
  int n = vertices.length;
  Point2D[] ouput = new Point2D[2 * n];
  for (int i = 0; i < n; ++i) {
    ouput[i] = r[i];
    ouput[2 * n - 1 - i] = l[i];
  }
  
  //--
  return ouput;
}

class Curve {
  //--
  Point2D[] vertices;
  Point2D[] curvePts;
  Point2D[] spinePts;
  float[] radii;
  int[] radiiVertices;
  
  //--
  ArrayList<Point2D> controlPts;
  FloatList controlRadii;

  //--
  int selectedControlPt = -1;
  int selectedControlRadius = -1;
  float selectedControlRadiusR;
  float selectedControlRadiusD;

  Curve() {
    //--   
    this.controlPts   = new ArrayList<Point2D>();
    this.controlRadii = new FloatList();   
  }

  boolean isValid() {
    return this.controlPts.size() >= bspline_degree;
  }

  boolean hasControlPtSelected() {
    return (this.selectedControlPt != -1);
  }

  boolean hasControlRadiusSelected() {
    return (this.selectedControlRadius != -1);
  }

  void clear() {
    //--
    this.controlPts.clear();
    this.controlRadii.clear();
    this.selectedControlPt = -1;
  }

  void addControlPoint(float x, float y, float r) {
    this.controlPts.add(new Point2D(x, y));
    this.controlRadii.append(r);
  } 

  void generateCurve() {
    //--    
    int n = this.controlPts.size();
    this.radii = new float[num_curve_points];
    this.curvePts = new Point2D[num_curve_points];    
    this.radiiVertices = new int[n];
    
    //--
    BSpline bspline = new BSpline(bspline_degree, n, true);  

    //--
    float step = 1.0f / (num_curve_points-1);    
    float dist = 0.0f;
    int k = 0;
    int j = 0;
    float t = 0.0;

    //--
    for (; j < num_curve_points-1; ++j, t += step) {
      //--
      Point2D p = new Point2D(0, 0);
      Point2D pd = new Point2D(0, 0);
      float r = 0.0f;
      float rd = 0.0f;

      //--
      for (int i = 0; i < n; ++i) {
        //--
        float w = bspline.calcWeight(i, t);
        p = add(p, scale(this.controlPts.get(i), w));
        r += w * this.controlRadii.get(i);
        //println("j="+j+",i="+i+",t="+t+",w="+w+",p="+p.toString()+",r="+r);
      }           
      //println("t="+t+"r="+r+",p="+p.toString());

      //--
      float d = distance(p, this.controlPts.get(k));
      //println("j="+j+",d="+d+",dist="+dist);
      if (d < dist) {
        dist = d;
      } else {
        int s = max(0, j - 1);
        this.radiiVertices[k] = s;
        //println("k="+k+",val="+s+",j="+j+",t="+t);
        ++k;
        dist = distance(p, this.controlPts.get(k));
      }

      //--
      this.curvePts[j] = p;      
      this.radii[j] = r;
    }

    // generate last sample pointl
    Point2D p = this.controlPts.get(n-1);
    float r = this.controlRadii.get(n-1);
    //println("t="+t+",r="+r+",p="+p.toString());

    //--
    this.curvePts[j] = new Point2D(p.x, p.y);    
    this.radii[j] = r;   

    if (k < n-1) {
      this.radiiVertices[k] = j-1;
      //println("k="+k+",val="+(j-1));
      ++k;
    }
    this.radiiVertices[k] = j;
    //println("k="+k+",val="+(num_curve_points-1)+"#c="+n);

    //--
    int selectedControlPt = -1;
    int selectedControlRadius = -1;
  }
  
  void generateStroke(int strokeMode) {
    //--
    this.vertices = _generateStroke(this.spinePts, this.radii, strokeMode);
  } //<>//
  
  Vector2D[] generateTangents() {
    //--
    int n = this.curvePts.length;
    Vector2D[] T = new Vector2D[n];    

    //--   
    float delta = 1.0f / 50;
    for(int i=0; i < n; ++i) {
      //--
      Point2D p1 = this.curvePts[i];
      Point2D p0 = (i > 0) ? this.curvePts[i-1] : p1;
      Point2D p2 = (i < n-1) ? this.curvePts[i+1] : p1;
      Point2D pd = calcNevilleDerivative(p0, p1, p2, delta);
      T[i] = normalize(new Vector2D(pd.x, pd.y));
    }
    
    /*for(int i=0; i < n; ++i){
      if (i==0) {
        T[i] = normalize(new Vector2D(this.curvePts[i], this.curvePts[i+1]));
      } else if (i == n-1) {
        T[i] = normalize(new Vector2D(this.curvePts[i-1], this.curvePts[i]));
      } else {
        T[i] = normalize(new Vector2D(this.curvePts[i-1], this.curvePts[i+1]));
      }
    }*/    
    
    //--
    return T;
  }

  void renderCurve(int _color, boolean showRadii) {
    //--
    stroke(black);
    fill(_color);
    beginShape();    
    for (Point2D p : this.vertices) {
      vertex(p);
    }
    vertex(this.vertices[0]);
    endShape();

    noFill();
    stroke(black);
    for (int i = 0; i < this.curvePts.length-1; ++i) {
      Point2D p0 = this.curvePts[i];
      Point2D p1 = this.curvePts[i+1];
      line(p0, p1);
    }

    if (showRadii) {
      noFill();      
      for (int i = 0, n = this.controlPts.size (); i < n; ++i) {
        stroke((i == this.selectedControlRadius) ? red : blue);        
        float r = this.controlRadii.get(i);        
        int j = this.radiiVertices[i];
        Point2D p = this.curvePts[j];         
        circle(p.x, p.y, r);
      }
    }
  }

  void renderControlPts() {
    //--
    noStroke();    
    for (int i = 0, n = this.controlPts.size (); i < n; ++i) {
      fill((i == this.selectedControlPt) ? red : black);
      Point2D p = this.controlPts.get(i);
      circle(p.x, p.y, vertex_radius);
    }
  }

  void selectControlPt(boolean showCurves) {
    //--    
    this.selectedControlPt = -1;
    this.selectedControlRadius = -1;
    Point2D pm = new Point2D(mouseX, mouseY);
    for (int i = 0, n = this.controlPts.size (); i < n; ++i) {
      {
        Point2D p = this.controlPts.get(i);      
        float dist = distance(p, pm);      
        if (dist < vertex_radius) {           
          this.selectedControlPt = i;
          break;
        }
      }

      if (showCurves) {
        float r = this.controlRadii.get(i);
        int j = radiiVertices[i];
        //println("xxx i="+i+",j="+j);
        
        // check the size because the curve length 
        // could shrink after the correspondance mapping
        if (j < this.curvePts.length) {
          Point2D p = this.curvePts[j];
          float dist = distance(p, pm);
          if (dist < r) {
            this.selectedControlRadius = i;
            this.selectedControlRadiusR = r;
            this.selectedControlRadiusD = dist;
            break;
          }
        }
      }
    }
  }

  void moveVertex() {
    int i = this.selectedControlPt;
    if (i != -1) {
      Point2D p = this.controlPts.get(i);
      p.x = mouseX;
      p.y = mouseY;
    }
  }  

  void changeRadius() {
    int i = this.selectedControlRadius;
    if (i != -1) {
      float r = this.selectedControlRadiusR;
      Point2D pm = new Point2D(mouseX, mouseY);
      int j = this.radiiVertices[i];
      Point2D p = this.curvePts[j];      
      float dist = distance(p, pm);
      float diff = dist - this.selectedControlRadiusD;      
      float newRadius = max(min_stroke_radius, min(r + diff, max_stroke_radius));
      this.controlRadii.set(i, newRadius);
    }
  }  

  void save(DataOutputStream dos) {
    try {
      dos.writeInt(this.controlPts.size());
      for (int i = 0, n = this.controlPts.size (); i < n; ++i) {
        Point2D p = this.controlPts.get(i);
        float r = this.controlRadii.get(i);
        dos.writeFloat(p.x);
        dos.writeFloat(p.y);
        dos.writeFloat(r);
      }
    } 
    catch (IOException e) {
      println("failed to save the curve");
    }
  }

  void load(DataInputStream dis) {
    try {
      int size = dis.readInt();    
      while (size-- > 0) {
        float x = dis.readFloat();
        float y = dis.readFloat();
        float r = dis.readFloat();
        this.addControlPoint(x, y, r);
      }
    } 
    catch (IOException e) {
      println("failed to load the curve");
    }
  }
}

class BallMorphInfo {
  int correspondence;
  Point2D center;
  float radius;
};

class Face {
  Point2D[] spine;
  Point2D[] stroke;
};

class Editor2D implements IControlView {

  Curve c0, c1;  // curves

  Point2D[] spine;
  Face morphfaces[];
  ArrayList<BallMorphInfo> ballMorphTable;

  int strokeMode = 2;
  boolean showCurves = false; 
  boolean showRadii = false;
  boolean showBallMorph = false;

  Editor2D() {
    this.c0 = new Curve();
    this.c1 = new Curve();
    this.ballMorphTable = new ArrayList<BallMorphInfo>();
    this.morphfaces = new Face[num_morph_faces-2];
  }

  void init() {
    //--
  }

  boolean isReady() {
    if (!c0.isValid() || !c1.isValid()) {
      println("*** curves insufficient control points!");
      return false;
    }    
    return true;
  }

  void clear() {
    this.c0.clear();
    this.c1.clear();
    this.ballMorphTable.clear();
    this.showCurves = false; 
    this.showRadii = false;
    this.showBallMorph = false;
  }

  void render() { 
    //--
    hint(DISABLE_DEPTH_TEST);
    camera();
    noLights();
  
    // Render the split bar
    stroke(black);
    noFill();
    line(width/2, 0, width/2, height);

    // Render curves
    if (showCurves) {
      c0.renderCurve(green, this.showRadii);
      c1.renderCurve(red, this.showRadii);      
      // render spine
      stroke(red);
      noFill();
      line(spine[0], spine[1]);
      noStroke();
      fill(black); 
      circle(spine[0].x, spine[0].y, vertex_radius);
      circle(spine[1].x, spine[1].y, vertex_radius);
    }
    c0.renderControlPts();
    c1.renderControlPts();
    
    if (this.showBallMorph) { //<>//
      this.drawBallMorph();
    }
  }  
  
  void drawBallMorph(){
    //--
    stroke(black, 32);    
    
    //--
    for(Face face: morphfaces) {
      // render stroke
      fill(blue, 32);
      beginShape();      
      for (Point2D p : face.stroke) {
        vertex(p);
      }
      vertex(face.stroke[0]);
      endShape();
      noFill();
      
      // render curve
      for (int i=0; i < face.spine.length-1; ++i) {
        Point2D p = face.spine[i];
        Point2D q = face.spine[i+1];
        line(p.x, p.y, q.x, q.y);
      }
    }
  }

  void computeSpine() {
    //--
    float minY = MAX_FLOAT;
    float maxY = MIN_FLOAT;

    for (Point2D p : c0.vertices) {
      minY = min(minY, p.y);
      maxY = max(maxY, p.y);
    }

    for (Point2D p : c1.vertices) {
      minY = min(minY, p.y);
      maxY = max(maxY, p.y);
    }

    spine = new Point2D[] {
      new Point2D(width/2, minY), 
      new Point2D(width/2, maxY)
    };
  }
  
  void buildCorrespondance() {
    //--
    ballMorphTable.clear();
    
    //--
    Vector2D[] T0 = c0.generateTangents();
    Vector2D[] T1 = c1.generateTangents();
    
    //--
    int firstCorrespondence = -1;
    int lastCorrespondence = -1;    
    int n = c0.curvePts.length;
    int i = 0;
    
    //--
    for (;i < n; ++i) {
      //println("i="+i);
      Point2D centerBall = new Point2D(c0.curvePts[i]);
      float rBall = Float.MAX_VALUE;
      Vector2D N = normalize(rotate(T0[i], -PI/2));
      //println("N = " + N.x+","+N.y);
      
      int correspondence = -1;
      for (int j = lastCorrespondence+1; j < n; ++j) {
        //--
        Vector2D RL = new Vector2D(c1.curvePts[j], c0.curvePts[i]); 
        float rBallCurrent = -1.0 * dot(RL, RL)/(2.0 * dot(RL, N));
        if (rBallCurrent < rBall) {
          rBall = rBallCurrent;
          centerBall = translate(c0.curvePts[i], N, rBall);         
          correspondence = j;
        }
      }
      
      if (correspondence != -1) {      
        //--
        if (firstCorrespondence == -1) {
          firstCorrespondence = correspondence;
        }
        lastCorrespondence = correspondence;
        
        //--
        Vector2D PQ = new Vector2D(c0.curvePts[i], c1.curvePts[correspondence]);
        float radius = dot(PQ, sub(T0[i],T1[correspondence]))/dot(sub(T0[i],T1[correspondence]),sub(T0[i],T1[correspondence]));
        Point2D R = translate(c0.curvePts[i], T0[i], radius);      
        //println("Correspondence ="+correspondence);
        //println("R = "+R.x+","+R.y);
        //println("radius = " + radius);       
        
        //--
        BallMorphInfo info = new BallMorphInfo();
        info.correspondence = correspondence - firstCorrespondence;
        info.center = R;
        info.radius = radius;      
        ballMorphTable.add(info);
      } else {
        println("*** correspondance cutoff index="+i+",c1-start="+firstCorrespondence+",c1-end="+lastCorrespondence);
        break;
      }
    }
    
    if (i != n) {
      // trim excess vertices on c0
      c0.spinePts = Arrays.copyOf(c0.curvePts, i);
      c0.radii = Arrays.copyOf(c0.radii, i);
      
      // trim excess vertices on c1
      c1.spinePts = new Point2D[(lastCorrespondence - firstCorrespondence) + 1];
      arrayCopy(c1.curvePts, firstCorrespondence, c1.spinePts, 0, c1.spinePts.length);
      float[] radii = new float[(lastCorrespondence - firstCorrespondence) + 1];
      arrayCopy(c1.radii, firstCorrespondence, radii, 0, radii.length);
      c1.radii = radii;
    } else {
      c0.spinePts = c0.curvePts;
      c1.spinePts = c1.curvePts;
    }
  }
  
  Face createMorphFace(float t) {
    //--
    int n = c0.spinePts.length;    
    Point2D[] curve = new Point2D[n];
    float[] radii = new float[n];
    Point2D origin = this.spine[0];
    
    for (int i = 0; i < n; ++i) {
      //--
      BallMorphInfo info = this.ballMorphTable.get(i);
      int correspondence = info.correspondence;
      //--      
      Point2D p = c0.spinePts[i];
      Point2D q = c1.spinePts[correspondence];
      Point2D o = new Point2D(origin.x, q.y);
      Point2D qr= translate(o, new Vector2D(q, o));
      
      curve[i] = lerp(c0.spinePts[i], qr, t);
      radii[i] = lerp(c0.radii[i], c1.radii[correspondence], t);
    }      
    //--
    Face face = new Face();
    face.spine = curve;
    face.stroke = _generateStroke(curve, radii, editor.strokeMode);
    return face;
  }
  
  Face createBallMorphFace(float t) {
    //--
    int n = c0.spinePts.length;    
    Point2D[] curve = new Point2D[n];
    float[] radii = new float[n];
    
    for (int i = 0; i < n; ++i) {
      //--
      BallMorphInfo info = this.ballMorphTable.get(i);
      Point2D R = info.center;
      float radius = info.radius;
      int correspondence = info.correspondence;
      //--
      Vector2D RP = new Vector2D(R, c0.spinePts[i]);
      Vector2D RQ = new Vector2D(R, c1.spinePts[correspondence]);    
      float angle = angle(RP, RQ);
      if (angle < 0){
        angle += PI;
      }        
      //println("angle = "+angle);
      
      //--
      curve[i] = translate(R, lerp(RP, RQ, t));
      radii[i] = lerp(c0.radii[i], c1.radii[correspondence], t);
    }      
    //--
    Face face = new Face();
    face.spine = curve;
    face.stroke = _generateStroke(curve, radii, editor.strokeMode);
    return face;
  } 
  
  void createBallMorphAnimation() {
    //--
    float delta = 1.0 / (num_morph_faces-1);    
    for (int i = 1; i < num_morph_faces-1; ++i) {
      this.morphfaces[i-1] = this.createBallMorphFace(i * delta);
    }
  }

  void processCurves() {
    //--
    if (!this.isReady()) {
      return;
    }
    
    //--
    this.c0.generateCurve();
    this.c1.generateCurve();
    this.buildCorrespondance();
    this.c0.generateStroke(this.strokeMode);
    this.c1.generateStroke(this.strokeMode);
    this.computeSpine();
    
    //--
    this.showBallMorph = false;
  }

  void onKeyPressed() {
    //--
    if (mousePressed) {
      return;
    }

    switch (key) {
    case 'c':
      this.clear();
      break;
    case ' ':
      if (this.isReady()) {
        if (!this.showCurves) {
          this.processCurves();
          this.showCurves= true;
        } else {
          this.showCurves = false;
        }
      }
      break;
    case 'r':
      this.showRadii = !this.showRadii;
      break;
    case 'k':
      if (this.strokeMode == 0) {
        this.strokeMode = 1;
        displayText("Radial Offset");
      } else 
        if (this.strokeMode == 1) {
        this.strokeMode = 2;
        displayText("Ball Offset");
      } else {
        this.strokeMode = 0;
        displayText("Normal Offset");
      }  
      this.processCurves();
      break;
    case 'a':
      if (this.showCurves && !this.showBallMorph) {
        this.createBallMorphAnimation();
        this.showBallMorph= true;
      } else {
        this.showBallMorph = false;
      }
      break;
    case 's':
      this.save("data/save.dat");
      break; 
    case 'l':            
      this.load(models[model_index++ % models.length]);
      break;
    }
  }
  
  void onKeyReleased() {
    //--
  }

  void onMouseMoved() {
    //--
    c0.selectControlPt(this.showCurves);
    c1.selectControlPt(this.showCurves);
  }

  void onMouseDragged() {
    //--
    if (c0.hasControlRadiusSelected()) {
      c0.changeRadius();
    } else
      if (c1.hasControlRadiusSelected()) {
      c1.changeRadius();
    } else
      if (c0.hasControlPtSelected()) {
      c0.moveVertex();
    } else
      if (c1.hasControlPtSelected()) {
      c1.moveVertex();
    }
  }

  void onMousePressed() {
    //--
    if (!this.showCurves && !this.showRadii) {
      if (mouseX < width/2) {
        if (!c0.hasControlPtSelected()) {
          c0.addControlPoint(mouseX, mouseY, default_stroke_radius);
        }
      } else {
        if (!c1.hasControlPtSelected()) {
          c1.addControlPoint(mouseX, mouseY, default_stroke_radius);
        }
      }
    }
  }

  void onMouseReleased() {
    //--
    if (this.showCurves) {
      boolean applyUpdate = false;
      if (c0.hasControlPtSelected() || c0.hasControlRadiusSelected()) {
        applyUpdate = true;
      } else
        if (c1.hasControlPtSelected() || c1.hasControlRadiusSelected()) {
        applyUpdate = true;
      }
      if (applyUpdate) {
        this.processCurves();
      } 
    }
  }

  void save(String filename) {
    if (!this.isReady()) {
      return;
    }
    try {
      DataOutputStream dos = new DataOutputStream(createOutput(filename));
      c0.save(dos);
      c1.save(dos);
      dos.close();
    } 
    catch (IOException e) {
      println("failed to save the model");
    }
  }

  void load(String filename) {  
    this.clear();
    try {    
      DataInputStream dis = new DataInputStream(createInput(filename));    
      c0.load(dis);
      c1.load(dis);
      dis.close();
    } 
    catch (IOException e) {
      println("failed to load the model");
    }
  }
};

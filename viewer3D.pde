//
//**************************** 3D Graph Library Definitions ***********************
//  

class Face3D {
  IntList vIndices;
  Vector3D normal;
  Point2D[] UVs;
  boolean isCap;
  
  Face3D() {
    this.vIndices = new IntList();
    this.isCap = false;
  }
  
  void clear() {
    this.vIndices.clear();
  }
  
  void render(ArrayList<Point3D> vertices, Vector3D[] normals, PImage texture, Point3D eye, int _color, boolean wireframe) {
    //--
    if (wireframe) {
      int alpha = this.isBackfacing(eye, vertices) ? 32 : 255;      
      stroke(_color, alpha); 
      noFill(); 
    } else {
      fill(_color);
      noStroke();
    }
    
    //--
    beginShape();
    
    if (!wireframe && texture != null) {
      textureMode(IMAGE);
      textureWrap(REPEAT);
      texture(texture);
    }
    
    for (int i = 0, n = this.vIndices.size(); i < n; ++i) {
      //--
      int vIndex = this.vIndices.get(i);
      
      //--
      if (!wireframe) {
        normal((normals != null) ? normals[vIndex] : this.normal);
      }
      
      //--
      if (!wireframe && texture != null) {
        Point2D uv = this.UVs[i];
        vertex(vertices.get(vIndex), uv.x, uv.y);
      } else {
        vertex(vertices.get(vIndex));
      }
    }    
    
    //--
    if (wireframe) {
      vertex(vertices.get(this.vIndices.get(0)));
    }
    endShape();
    
    /*// debug normals
    if (normals  != null) {
      noFill();
      stroke(black);
      for (int i = 0, n = this.vIndices.size(); i < n; ++i) {
        int vIndex = this.vIndices.get(i);
        Point3D P1 = vertices.get(vIndex);
        Vector3D N = normals[vIndex];
        Point3D P2 = translate(P1, N, 10);
        line(P1.x, P1.y, P1.z, P2.x, P2.y, P2.z);
      }
    }*/   
  }
  
  float computeVolume(ArrayList<Point3D> vertices) {
    // reference: http://research.microsoft.com/en-us/um/people/chazhang/publications/icip01_ChaZhang.pdf
    float volume = 0;
    Point3D O = new Point3D(0, 0, 0);
    Vector3D v0 = new Vector3D(O, vertices.get(this.vIndices.get(0))); 
    for (int i = 1, n = this.vIndices.size(); i+1 < n; ++i) {
      Vector3D v1 = new Vector3D(O, vertices.get(this.vIndices.get(i)));
      Vector3D v2 = new Vector3D(O, vertices.get(this.vIndices.get(i+1)));
      volume += dot(v0, cross(v1, v2));
    }
    return volume;
  }
  
  void generateNormal(ArrayList<Point3D> vertices, Vector3D[] normals, int[] nFaces) {
    // using Newell algorithm
    // reference: https://www.opengl.org/wiki/Calculating_a_Surface_Normal
    Vector3D N = new Vector3D(0, 0, 0);
    for (int i = 0, n = this.vIndices.size(); i < n; ++i) {
      Point3D v  = vertices.get(this.vIndices.get(i));
      Point3D vn = vertices.get(this.vIndices.get((i+1) % n));
      N.x += (v.y - vn.y) * (v.z + vn.z);
      N.y += (v.z - vn.z) * (v.x + vn.x);
      N.z += (v.x - vn.x) * (v.y + vn.y);      
    }
    this.normal = normalize(N);
    
    //--
    for (int vIndex: this.vIndices) {
      int nFace = nFaces[vIndex];
      Vector3D normal = normals[vIndex];
      if (nFace == 0) {
        normal.x = N.x;
        normal.y = N.y;
        normal.z = N.z;
      } else {
        normal.x += N.x;
        normal.y += N.y;
        normal.z += N.z;
      }
      nFaces[vIndex] = nFace + 1;
    }
  }
  
  void generateUVs(Point2D[] vertices) {
    //-
    this.UVs = new Point2D[vertices.length];
    float[] bb = computeBoundingBall(vertices);
    float minx = bb[0];
    float miny = bb[2];
    for (int i = 0; i < vertices.length; ++i) {
      Point2D p = vertices[i];
      this.UVs[i] = new Point2D(p.x - minx, p.y - miny);
    }
  }
  
  boolean isBackfacing(Point3D eye, ArrayList<Point3D> vertices) {
    //--
    int i0 = this.vIndices.get(0);
    Vector3D vView = new Vector3D(eye, vertices.get(i0));
    //println("vView="+vView.toString()+",normal="+this.normal);
    return dot(vView, this.normal) >= 0; 
  }
};

class Model3D {
  //--
  ArrayList<Point3D> vertices;
  Vector3D[] normals;
  Vector2D[] uvs;
  ArrayList<Point3D> bendVertices;
  ArrayList<Point3D> _vertices;
  ArrayList<Face3D> faces;
  Point3D[] spine;
  int[] controlPtsIndices;
  int num_spine_points;
  float spineDistance;
  PImage texture;
  
  Model3D() {
    //--
    this.vertices = new ArrayList<Point3D>();
    this.faces = new ArrayList<Face3D>();
    this._vertices = vertices;
  }
  
  void init() {
    texture = loadImage("data/texture.jpg");
  }
  
  void clear() {
    //--
    this.vertices.clear();
    this.normals = null;
    this.faces.clear();
    this._vertices = vertices;
  }
  
  void render(Point3D eye, boolean wireframe, boolean textureMode) {
    //--
    int n = this.faces.size();
    for (int i = 1; i < n-1; ++i) {
      Face3D face = this.faces.get(i); 
      Vector3D[] _normals = (textureMode && !face.isCap) ? this.normals : null;  
      face.render(this._vertices, _normals, null, eye, face.isCap ? grey : yellow, wireframe);
    }
    //--
    PImage _texture = textureMode ? texture : null;
    this.faces.get(0).render(this._vertices, null, _texture, eye, green, wireframe);
    this.faces.get(n-1).render(this._vertices, null, _texture, eye, red, wireframe);
  }
  
  void restoreVertices() {
    //--
    this._vertices = this.vertices;
    this.generateNormals();
  }
  
  void applyBending(Point3D A, Point3D C, float theta) {
    int j;
    Vector3D I, J, K;
    Vector3D Is, Js, Ks;
    Vector3D T;
    Vector3D AP;
    Vector3D AC = new Vector3D(A,C);
    I = normalize(AC);
    K = new Vector3D(0,0,1);
    J = normalize(cross(K,I));
    //println("I = " + I.x+","+I.y+","+I.z);
    //println("J = " + J.x+","+J.y+","+J.z);
    //println("K = " + K.x+","+K.y+","+K.z);
    Point3D pP;
    float spineTravelled;
    
    float toTravelZ, toTravelX, toTravelY;
    
    for (int i = 0, n = this.vertices.size(); i < n; ++i){
      //println("vertex = " + this.vertices.get(i).x+","+this.vertices.get(i).y+","+this.vertices.get(i).z);
      AP = new Vector3D(this.spine[0], this.vertices.get(i));
      toTravelZ = abs(dot(AP, K));
      toTravelX = dot(AP, I);
      toTravelY = dot(AP, J);
      //println("toTravelX = " + toTravelX +",toTravelY="+toTravelY+",toTravelZ="+toTravelZ);
      spineTravelled = 0.0;
      pP = new Point3D(this.vertices.get(i));
      for(j = 0; j <num_spine_points; j++){
        
        if(j < num_spine_points - 1){
          spineTravelled += distance(this.spine[j], this.spine[j+1]);
        }
        if(spineTravelled > toTravelZ){
          break;
        }
      }
      
      if(j==num_spine_points) j--; // ugly hack
      //println("stop at spineIndex = " + j);
      if(j == 0){
        Ks = rotate(normalize(new Vector3D(this.spine[j], this.spine[j+1])), theta, I, J);
      }
      else if(j == num_spine_points - 1){
        Ks = rotate(normalize(new Vector3D(this.spine[j-1], this.spine[j])), theta, I, J);
      }
      else{
        Ks = rotate(normalize(new Vector3D(this.spine[j-1], this.spine[j+1])), theta, I, J);
      }
      //println("Ks = " + Ks.x+","+Ks.y+","+Ks.z);
      Js = neg(normalize(cross(normalize(AC), Ks)));
      //println("Js = " + Js.x+","+Js.y+","+Js.z);
      Is = normalize(cross(Js, Ks));
      //println("Is = " + Is.x+","+Is.y+","+Is.z);
      pP = translate(translate(this.spine[j], Is, toTravelX), Js, toTravelY);
       
      this.bendVertices.set(i,pP);
       
    }
    
    //--
    this._vertices = this.bendVertices;
    this.generateNormals();
  }
  void refine(){
    Point3D[] spineNew = new Point3D[2 * num_spine_points - 1];
    int j = 0;
    for (int i = 0; i < num_spine_points - 1; i++){
      spineNew[j++] = new Point3D(spine[i]);
      spineNew[j++] = new Point3D((spine[i].x + spine[i+1].x)/2,(spine[i].y + spine[i+1].y)/2,(spine[i].z + spine[i+1].z)/2);
    }
    spineNew[j++] = new Point3D(spine[num_spine_points - 1]);
    num_spine_points = 2*num_spine_points - 1;
    spine = spineNew;
    
    for(int i=0; i<num_spine_control_points; i++){
      controlPtsIndices[i] = controlPtsIndices[i] * 2;
    }
  }
  
  void dual(){
    Point3D[] spineNew = new Point3D[num_spine_points + 1];
    spineNew[0] = new Point3D(spine[0]);
    spineNew[num_spine_points] = new Point3D(spine[num_spine_points - 1]);
    
    for (int i = 0; i < num_spine_points - 1; i++){
      spineNew[i+1] = new Point3D((spine[i].x + spine[i+1].x)/2,(spine[i].y + spine[i+1].y)/2,(spine[i].z + spine[i+1].z)/2);
    }
    num_spine_points = num_spine_points + 1;
    spine = spineNew; 
    
    controlPtsIndices[num_spine_control_points - 1] = num_spine_points - 1;
  }
  void calcStraightSpine() {
    //--
    float delta = this.spineDistance / (num_spine_control_points - 1);
    for(int i = 1; i < num_spine_control_points; ++i) {
      this.spine[i] = new Point3D(0, 0, this.spine[0].z + i * delta);
    }
    for(int i=0; i < num_spine_control_points; i++){
      controlPtsIndices[i] = i;
    }
    num_spine_points = num_spine_control_points;
    println("num_spine_points = " + num_spine_points);
  }
  void calcBendedSpine(boolean isBended) {
    //--
    if(!isBended){
    float delta = this.spineDistance / (num_spine_control_points - 1);
    int random;
    for(int i = 1; i < num_spine_control_points; ++i) {
      random = (int)(Math.random() * 40 + (-20));
      this.spine[i] = new Point3D(this.spine[0].x + random, 0, this.spine[0].z + i * delta);
      //this.spine[i] = new Point3D(0, 0, this.spine[0].z + i * delta);
    }
    for(int i=0; i < num_spine_control_points; i++){
      controlPtsIndices[i] = i;
    }
    num_spine_points = num_spine_control_points;
    refine();dual();dual();dual();dual();refine();dual();dual();dual();dual();
    refine();dual();dual();dual();dual();
    refine();dual();dual();dual();dual();
    refine();dual();dual();dual();dual();
    }
  }
  
  void generateNormals() {
    //--
    if (this.normals == null) {
      this.normals = new Vector3D[this._vertices.size()];
      for (int i = 0; i < this.normals.length; ++i) {
        this.normals[i] = new Vector3D(0, 0, 0);
      } 
    }
    
    //--
    int[] nfaces = new int[this.normals.length];
    for (int i = 0; i < nfaces.length; ++i) {
      nfaces[i] = 0;
    }
    
    //--
    for (Face3D face : this.faces) {
      face.generateNormal(this._vertices, this.normals, nfaces);
    }
    
    //--
    for (int i = 0; i < this.normals.length; ++i) {      
      this.normals[i] = normalize(scale(this.normals[i], 1.0f / nfaces[i]));
      //println("i="+i+",#n="+nfaces[i]+",n="+this.normals[i].toString()+",ns="+norm(this.normals[i]));
    }
  }
  void build(Editor2D editor, float sweep_angle) {
    //--
    this.clear();
    
    // generate all 3D faces
    Face3D fp = new Face3D();
    float delta = 1.0 / (num_sweep_faces-1);
    for (int i = 0; i < num_sweep_faces-1; ++i) {
      float t = i * delta;
      this.createFaces(editor, fp, t, t * sweep_angle);
    }
    this.createFaces(editor, fp, 1.0, sweep_angle);
    
    // center mesh coordinates
    float[] bb = computeBoundingBall(this.vertices);
    float centerx = (bb[0] + bb[1]) / 2;
    float centerz = (bb[4] + bb[5]) / 2;
    for (Point3D p: vertices) {
      p.x -= centerx;    
      p.z -= centerz;
    }
    
    // Allocate the bend vertices
    this.bendVertices = new ArrayList<Point3D>(this.vertices.size());
    for (Point3D p : this.vertices) {
      this.bendVertices.add(p);
    }
    
    float offsetX = (editor.spine[0].x + editor.spine[1].x) / 2;
    float offsetY = (editor.spine[0].y + editor.spine[1].y) / 2;
    this.spine = new Point3D[num_spine_control_points];
    this.controlPtsIndices = new int[num_spine_control_points];
    //this.spineBended = new Point3D[num_spine_points];    
    this.spine[0] = new Point3D(editor.spine[0].x - offsetX, 0, editor.spine[0].y - offsetY);
    this.spineDistance = distance(editor.spine[0], editor.spine[1]);
    this.calcStraightSpine();
    
    
    // generate normals
    this.generateNormals();
  }
  
  void createFaces(Editor2D editor, Face3D fp, float t, float angle) {  
    //--
    Point2D origin = editor.spine[0];    
    Face3D endFace = null;
    
    //--
    Point2D[] vertices;    
    if (t == 0.0) {
      vertices = editor.c0.vertices;
      endFace = new Face3D();
    } else if (t == 1.0) {
      Face face2D = editor.createMorphFace(1.0);
      vertices = face2D.stroke;
      endFace = new Face3D();
    } else {
      Face face2D = editor.createMorphFace(t);
      vertices = face2D.stroke;
    }
    
    //--    
    int n = vertices.length;
    int vfpp = (t > 0.0) ? fp.vIndices.get(n-1) : 0;
    int vn = this.createVertex(this.createPoint(origin, vertices[n-1], angle));
    int vp = vn;
        
    //--
    for (int i = 0; i < n; ++i) {
      //--
      int vfp = (t > 0.0) ? fp.vIndices.get(i) : 0;
      int v = (i < n-1) ? this.createVertex(this.createPoint(origin, vertices[i], angle)) : vn;
            
      //--
      if (endFace != null) {
        endFace.vIndices.append(v);
      }
      
      //--
      if (t > 0.0) {
        Face3D face3D = new Face3D();
        face3D.vIndices.append(vfp);
        face3D.vIndices.append(vfpp);
        face3D.vIndices.append(vp);
        face3D.vIndices.append(v);
        face3D.isCap = (i == 0) || (i == n/2);
        this.faces.add(face3D);
      }     
      
      //--
      vp = v;
      vfpp = vfp;
      if (t > 0.0) {
        fp.vIndices.set(i, v);
      } else {
        fp.vIndices.append(v);
      }
    }
    
    if (endFace != null) {      
      if (t == 1.0) {
        // change the indices orientation for backfacing test
        endFace.vIndices.reverse();
      }
      this.faces.add(endFace);
      
      // Generate texture coordinates
      endFace.generateUVs(vertices);
    }
  }
  
  Point3D createPoint(Point2D origin, Point2D p, float angle) {
    //--
    Point2D oX = new Point2D(origin.x, 0);
    Point2D pX = new Point2D(p.x, 0);
    Vector2D v = rotate(new Vector2D(oX, pX), angle);
    Point2D pn = translate(oX, v);    
    //--
    return new Point3D(pn.x, pn.y, height - p.y);
  } 
  
  int createVertex(Point3D p) {
    int index = this.vertices.size();
    this.vertices.add(p);
    return index;
  } 
  
  float computeVolume() {
    float volume = 0;
    for (Face3D face: this.faces) {
      volume += face.computeVolume(this._vertices);
    }
    return abs(volume);
  }
};

class Viewer3D implements IControlView {
  //--
  Model3D model;
  //Point3D[] spine;
  //float spineDistance;
  float radius = default_bend_radius;
  float theta = 0;  
  int theta_x = 0; 
  float sweep = default_sweep_angle;  
  int sweep_x = (int)((sweep * 180) / PI); 
  
  //--
  float dz = 0; // distance to camera. Manipulated with wheel or when 
  float rx =-0.06 * TWO_PI;
  float ry =-0.04 * TWO_PI; // view angles manipulated when space pressed but not mouse
  
  //--
  boolean wireframe = false;
  boolean textureMode = false;
  boolean showFloor = true;
  boolean isBended = false;
  boolean changeSpinePlane = false;
  boolean changeSweepAngle = false;
  boolean showSpineControlPtsBalls = false;
  
  //--
  Vector3D I = new Vector3D(1,0,0);
  Vector3D J = new Vector3D(0,1,0);
  Vector3D K = new Vector3D(0,0,1);
  
  Viewer3D() {
    //--
    model = new Model3D();
  }
  
  void init() {
    //--
    model.init();
  }
  
  void clear() {
    //--
    model.clear();
  }
  
  void render() {    
    //--
    hint(ENABLE_DEPTH_TEST);
    pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
    camera();       // sets a standard perspective
    translate(width/2,height/2,dz); // puts origin of model at screen center and moves forward/away by dz
    lights();  // turns on view-dependent lighting
    rotateX(rx); rotateY(ry); // rotates the model around the new origin (center of screen)
    rotateX(PI/2); // rotates frame around X to make X and Y basis vectors parallel to the floor
    
    if (showFloor) {
      noStroke();
      // X-red, Y-green, Z-blue reference frame
      pushMatrix(); translate(0,0,+this.model.spine[0].z); drawFrame(50); popMatrix();
     
      if (changeSpinePlane) {
        fill(black); 
        pushMatrix(); 
        translate(0,0,+this.model.spine[0].z);
        rotateZ(this.theta);
        rotateY(PI/2); 
        drawArrow(50,50/10); 
        popMatrix();
      }
      
      if (changeSweepAngle) {
        fill(black); 
        pushMatrix(); 
        translate(0,0,+this.model.spine[0].z);
        rotateZ(this.sweep - PI);
        rotateY(PI/2); 
        drawArrow(50,50/10); 
        popMatrix();
      }
  
      fill(0xfff5f5f5, 100); pushMatrix(); translate(0,0,+this.model.spine[0].z); box(100,100,0.5); popMatrix(); // draws floor as thin plate
    }
    
    //--
    //directionalLight(0.8, 0.8, 0.8, 0, 0, -1);
    
    //--
    Point3D eye = getCameraPosition();
    //println("eye="+eye);
    
    if(fDragging != true){
      g_center = pick( mouseX, mouseY );
    }
    
    // render the spine
    pushMatrix();
    rotateZ(theta);
    stroke(grey);noFill();
    for (int i = 0, n = this.model.num_spine_points-1; i < n; ++i) {
      line(this.model.spine[i], this.model.spine[i+1]);
    }
    popMatrix();
    
    pushMatrix();
    rotateZ(theta);
    if(this.showSpineControlPtsBalls){
      for(int i=0; i <num_spine_control_points; i++){
        Point3D p = this.model.spine[this.model.controlPtsIndices[i]];
        pushMatrix(); translate(p.x, p.y, p.z); fill(red);sphere(5); popMatrix();
      }
    }
    popMatrix();   
    
    // render the model 
    model.render(eye, this.wireframe, this.textureMode);

    // render the 3D model here! 
    popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas   
  }
  
  void build() {
    //--
    this.clear();
    
    //--
    this.model.build(editor, this.sweep);
    
    //--
    /*float offsetX = (editor.spine[0].x + editor.spine[1].x) / 2;
    float offsetY = (editor.spine[0].y + editor.spine[1].y) / 2;
    this.spine = new Point3D[num_spine_points];    
    this.spine[0] = new Point3D(editor.spine[0].x - offsetX, 0, editor.spine[0].y - offsetY);
    this.spineDistance = distance(editor.spine[0], editor.spine[1]);*/
    
    //--
    //this.generateStraightSpine();
  } 
  
  /*void generateStraightSpine() {
    //--
    float delta = this.spineDistance / (num_spine_points - 1);
    for(int i = 1; i < num_spine_points; ++i) {
      this.spine[i] = new Point3D(this.spine[0].x, 0, this.spine[0].z + i * delta);
    }
  }*/
  
  /*void generateBendedSpine() {
    //--
    float delta = this.spineDistance / (this.radius * (num_spine_points - 1));
    for(int i = 1; i < num_spine_points; ++i) {
      float alpha = i * delta;
      Point3D p = new Point3D(this.spine[0].x + this.radius - this.radius * cos(alpha), this.spine[0].y, this.spine[0].z + this.radius * sin(alpha));
      this.spine[i] = translate(this.spine[0], rotate(new Vector3D(this.spine[0], p), theta, I, J));
    }
  }*/
  
  void restoreModel() {
    //--
    this.model.calcStraightSpine();   
    this.model.restoreVertices();
  }
  
  void applyBending() {  
    //--
    this.model.calcBendedSpine(this.isBended);
    //---
    Point3D A = this.model.spine[0];
    Point3D C = translate(A, rotate(new Vector3D(A, translate(A, I, radius)), theta, I, J));
    model.applyBending(A, C, theta);
  }
  
  void onKeyPressed() {
    //--
    if (mousePressed) {
      return;
    }
    
    switch (key) {
    case '_':
      this.showFloor=!this.showFloor;
      break;
    case 'w':
      this.wireframe=!this.wireframe;
      break;
    case 'b':
      if (this.isBended) {
        this.restoreModel();
        this.isBended = false;        
      } else {
        this.applyBending();
        this.isBended = true;
      }
      break;
    case 'c':
      if(this.showSpineControlPtsBalls){
        this.showSpineControlPtsBalls = false;
      }
      else{
        this.showSpineControlPtsBalls = true;
      }
      break;
    case 'n':
      this.changeSpinePlane = true;
      break;
    case 'm':
      this.changeSweepAngle = true;
      break;
    case 't':
      this.textureMode=!this.textureMode;
      break;
    case 'v':
      float volume = this.model.computeVolume();
      displayText("volume="+volume);
      break;
    }
  }
  
  void onKeyReleased() {
    //--
    if (changeSweepAngle) {
      this.model.build(editor, this.sweep);
      if (this.isBended) {
        this.applyBending();
      }
      this.changeSweepAngle = false;
    }    
  }
  
  void onMouseMoved() {
    //--
    if (keyPressed && key==' ') {
      rx -= PI * (mouseY-pmouseY) / height; 
      ry += PI * (mouseX-pmouseX) / width;
    }    
    if (keyPressed && key=='z') {
      dz += (float)(mouseY-pmouseY); // approach view (same as wheel)
    }
    //--
    if (changeSweepAngle) {
      this.sweep_x = min(max(this.sweep_x + pmouseX - mouseX, 0), 360);
      this.sweep = this.sweep_x * (PI / 180);
    }
  }

  void onMouseDragged() {    
    //--
    fDragging = true;
    
    if (this.isBended) {
      if (mouseButton == LEFT) {
        this.radius = min(max(this.radius + pmouseX - mouseX, 200), 900);        
      } else {
        this.changeSpinePlane = true;
        this.theta_x = min(max(this.theta_x + pmouseX - mouseX, 0), 180);
        this.theta = this.theta_x * (PI / 180);
      }
      this.applyBending();
    }
  }
  
  void onMousePressed() {
    //--
  }
  
  void onMouseReleased() {
    //--
    fDragging = false;
    this.changeSpinePlane = false;
  }  
};

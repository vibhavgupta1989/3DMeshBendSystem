//
//**************************** 3D Library Definitions ***********************
//

import java.nio.*;

class Point3D {
  float x;
  float y;
  float z;
 
 Point3D(Point3D p) {
    this.x = p.x;
    this.y = p.y;
    this.z = p.z;
  } 
  
  Point3D(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  } 
  
  boolean equal(Point3D p) {
    return this.x == p.x && this.y == p.y && this.z == p.z;
  }
  
  String toString() {
    return "("+x+","+y+","+z+")";
  }
  
  Point2D toPoint2D() {
    return new Point2D(this.x, this.y);
  }
}

class Vector3D {
  float x;
  float y;
  float z;
 
 Vector3D(Vector3D v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
  } 
  
  Vector3D(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
 
  Vector3D(Point3D v0, Point3D v1) {
    this.x = v1.x - v0.x;
    this.y = v1.y - v0.y;
    this.z = v1.z - v0.z;
  }
  
  String toString() {
    return "("+x+","+y+","+z+")";
  }
}

Vector3D neg(Vector3D v) {
  return new Vector3D(-v.x, -v.y, -v.z);
}

Vector3D cross(Vector3D u, Vector3D v) {
  // UxV cross product (normal to both)
  return new Vector3D(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x); 
};

Vector3D axisNormal(Vector3D v) {
  if (abs(v.z) <= min(abs(v.x),abs(v.y))) {
    return new Vector3D(-v.y,v.x,0);
  }
  if (abs(v.x) <= min(abs(v.z),abs(v.y))) {
    return new Vector3D(0,-v.z,v.y);
  }
  return new Vector3D(v.z,0,-v.x);
}

Vector3D normal(Vector3D u, Vector3D v) {
  return cross(u, v); 
};  

Vector3D normal(Point3D A, Point3D B, Point3D C) {
  // normal to triangle (A,B,C), not normalized (proportional to area)
  return cross(new Vector3D(A,B), new Vector3D(A,C)); 
};                                                   

Vector3D add(Vector3D u, Vector3D v) {
  return new Vector3D(u.x + v.x, u.y + v.y, u.z + v.z);
}

Vector3D sub(Vector3D u, Vector3D v) {
  return new Vector3D(u.x - v.x, u.y - v.y, u.z - v.z);
}

Point3D translate(Point3D p, Vector3D dir, float dist) {
  return new Point3D(p.x + dist * dir.x, p.y + dist * dir.y, p.z + dist * dir.z); 
}

Point3D translate(Point3D p, Vector3D dir) {
  return translate(p, dir, 1);
}

float distanceSq(Point3D u, Point3D v) {
  return sq(v.x - u.x) + sq(v.y - u.y) + sq(v.z - u.z);
}

float distance(Point3D u, Point3D v) {
  return sqrt(distanceSq(u, v));
}

float dot(Vector3D u, Vector3D v) {
  //u*v dot product
  return u.x * v.x + u.y * v.y + u.z * v.z; 
}

float det2(Vector3D u, Vector3D v) {
  return u.x * v.y - u.y * v.x; 
};

float det3(Vector3D u, Vector3D v) {
  // u|v det product
  return sqrt(dot(u,u) * dot(v,v) - sq(dot(u,v))); 
}

float mixed(Vector3D u, Vector3D v, Vector3D W) {
  // (UxV)*W  mixed product, determinant  
  return dot(u, cross(v,W)); 
}

float normSq(Vector3D v) { 
  return sq(v.x) + sq(v.y) + sq(v.z);
}

float norm(Vector3D v) { 
  return sqrt(normSq(v));
}

Vector3D normalize(Vector3D v) {
  float inv_size = 1.0 / norm(v);
  return new Vector3D(v.x * inv_size, v.y * inv_size, v.z * inv_size);
}

Vector3D scale(Vector3D v, float s) {
  return new Vector3D(v.x * s, v.y * s, v.z * s);
}

Vector3D rotate(Vector3D v, float angle, Vector3D I, Vector3D J) {
  float x = dot(v, I), y = dot(v, J); 
  float c = cos(angle), s = sin(angle); 
  return add(v, add(scale(I, x * c - x - y * s), scale(J, x * s + y * c - y)));
}

// rendering ============================================================================

void normal(Vector3D v) {
  // emit normal for smooth shading
  normal(v.x,v.y,v.z);
}

void vertex(Point3D p) {
  // emit vertex for shading or drawing
  vertex(p.x,p.y,p.z);
}

void vertex(Point3D p, float u, float v) {
  // emit vertex with texture coordinates
  vertex(p.x,p.y,p.z,u,v);
}                          

void line(Point3D p, Point3D q) {
  // draws edge (p,q)
  line(p.x,p.y,p.z,q.x,q.y,q.z); 
}

void line(Point3D p, Vector3D v) {
  // shows edge from p to p+v
  line(p.x,p.y,p.z,p.x+v.x,p.y+v.y,p.z+v.z); 
};

void line(Point3D p, Vector3D v, float d) {
  // shows edge from p to p+dV
  line(p.x,p.y,p.z,p.x+d*v.x,p.y+d*v.y,p.z+d*v.z); 
}; 

void drawTriangle(Point3D A, Point3D B, Point3D C) {
  // draw triangle
  beginShape(); 
  vertex(A);
  vertex(B); 
  vertex(C); 
  endShape(CLOSE);
}

void drawquad(Point3D A, Point3D B, Point3D C, Point3D D) {
  // draw quad
  beginShape(); 
  vertex(A); 
  vertex(B); 
  vertex(C); 
  vertex(D); 
  endShape(CLOSE);
} 

void drawSphere(Point3D p, float r) {
  // render sphere of radius r and center p
  pushMatrix(); 
  translate(p.x,p.y,p.z); 
  sphere(r); 
  popMatrix();
}

void drawShadow(Point3D p, float r) {
  pushMatrix(); 
  translate(p.x,p.y,0); 
  scale(1,1,0.01); 
  sphere(r); 
  popMatrix();
}

void drawDisk(Point3D p, Vector3D v, float r) {  
  Vector3D I = normalize(axisNormal(v));
  Vector3D J = normalize(normal(I,v));
  drawDisk(p,I,J,r);
}

void drawDisk(Point3D p, Vector3D I, Vector3D J, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
  vertex(p);
  for(float a=0; a<=TWO_PI+da; a+=da) {
    vertex(translate(p, add(scale(I, r*cos(a)), scale(J, r*sin(a)))));
  }
  endShape();
}
  
void drawCollar(Point3D p, Vector3D v, float r, float rd) {
  Vector3D I = normalize(axisNormal(v));
  Vector3D J = normalize(normal(I,v));
  drawCollar(p,v,I,J,r,rd);
}
 
void drawCollar(Point3D p, Vector3D v, Vector3D I, Vector3D J, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
  for(float a=0; a<=TWO_PI+da; a+=da) {
    vertex(translate(p, add(scale(I, r*cos(a)), scale(J, r*sin(a))))); 
    vertex(translate(p, add(add(scale(I, rd*cos(a)), scale(J, rd*sin(a))), v)));
  }
  endShape();
}

void drawStub(Point3D p, Vector3D v, float r, float rd) {
  drawCollar(p,v,r,rd); 
  drawDisk(p,v,r); 
  drawDisk(translate(p,v),v,rd); 
}

Point3D toScreen(Point3D p) {
  return new Point3D(screenX(p.x,p.y,p.z), screenY(p.x,p.y,p.z), screenZ(p.x,p.y,p.z));
}

Point3D toModel(Point3D p) {
  return new Point3D(modelX(p.x,p.y,p.z), modelY(p.x,p.y,p.z), modelZ(p.x,p.y,p.z));
}

void drawText(Point3D p, String s) {
  // prints string s in 3D at p
  text(s, p.x, p.y, p.z); 
}

void drawText(Point3D p, String s, Vector3D D) {
  // prints string s in 3D at p+D
  text(s, p.x+D.x, p.y+D.y, p.z+D.z);  
}

Point3D pick(int mX, int mY)
{
  PGL pgl = beginPGL();
  FloatBuffer depthBuffer = ByteBuffer.allocateDirect(1 << 2).order(ByteOrder. nativeOrder()).asFloatBuffer();
  pgl.readPixels(mX, height - mY - 1, 1, 1, PGL.DEPTH_COMPONENT, PGL.FLOAT, depthBuffer);
  float depthValue = depthBuffer.get(0);
  depthBuffer.clear();
  endPGL();

  //get 3d matrices
  PGraphics3D p3d = (PGraphics3D)g;
  PMatrix3D proj = p3d.projection.get();
  PMatrix3D modelView = p3d.modelview.get();
  PMatrix3D modelViewProjInv = proj; modelViewProjInv.apply( modelView ); modelViewProjInv.invert();
  
  float[] viewport = {0, 0, p3d.width, p3d.height};
  
  float[] normalized = new float[4];
  normalized[0] = ((mX - viewport[0]) / viewport[2]) * 2.0f - 1.0f;
  normalized[1] = ((height - mY - viewport[1]) / viewport[3]) * 2.0f - 1.0f;
  normalized[2] = depthValue * 2.0f - 1.0f;
  normalized[3] = 1.0f;
  
  float[] unprojected = new float[4];
  
  modelViewProjInv.mult( normalized, unprojected );
  return new Point3D(unprojected[0]/unprojected[3], unprojected[1]/unprojected[3], unprojected[2]/unprojected[3]);
}

Point3D getCameraPosition() {
  //--
  PGraphics3D p3D = (PGraphics3D)g;
  PVector p = p3D.modelviewInv.mult(new PVector(1.0, 1.0, 1.0), null);
  return new Point3D(p.x, p.y, p.z);
}

float[] computeBoundingBall(ArrayList<Point3D> vertices) {
  //--
  float minx = MAX_FLOAT;
  float miny = MAX_FLOAT;
  float minz = MAX_FLOAT;
  float maxx = MIN_FLOAT;
  float maxy = MIN_FLOAT;
  float maxz = MIN_FLOAT;
  //--
  for (Point3D p: vertices) {
    if (p.x < minx) minx = p.x;
    if (p.y < miny) miny = p.y;
    if (p.z < minz) minz = p.z;
    if (p.x > maxx) maxx = p.x;
    if (p.y > maxy) maxy = p.y;
    if (p.z > maxz) maxz = p.z;      
  }
  return new float[] {minx, maxx, miny, maxy, minz, maxz};
}

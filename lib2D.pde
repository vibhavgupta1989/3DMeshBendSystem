//
//**************************** 2D Library Definitions ***********************
//

class Point2D {
  float x;
  float y;
 
 Point2D(Point2D p) {
    this.x = p.x;
    this.y = p.y;
  } 
  
  Point2D(float x, float y) {
    this.x = x;
    this.y = y;
  } 
  
  boolean equal(Point2D p) {
    return this.x == p.x && this.y == p.y;
  }
  
  String toString() {
    return "("+this.x+","+this.y+")";
  }
}

class Vector2D {
  float x;
  float y;
 
 Vector2D(Vector2D v) {
    this.x = v.x;
    this.y = v.y;
  } 
  
  Vector2D(float x, float y) {
    this.x = x;
    this.y = y;
  }
 
  Vector2D(Point2D v0, Point2D v1) {
    this.x = v1.x - v0.x;
    this.y = v1.y - v0.y;
  }
  
  String toString() {
    return "("+this.x+","+this.y+")";
  }
}

boolean isNormalized(float v) {
  return (v >= 0.0 && v <= 1.0);
}

float distanceSq(Point2D a, Point2D v) {
  return sq(v.x - a.x) + sq(v.y - a.y);
}

float distance(Point2D u, Point2D v) {
  return sqrt(distanceSq(u, v));
}

float clamp(float val, float min, float max) {
  if (val < min) return min; 
  if (val > max) return max; 
  return val;
}

Vector2D rotate90(Vector2D v) {
  // rotate 90 degree ccw on cartesian coordinates system
  return new Vector2D(-v.y, v.x);
}

Vector2D normalize(Vector2D v) {
  float inv_size = 1.0 / norm(v);
  return new Vector2D(v.x * inv_size, v.y * inv_size);
}

Vector2D rotate(Vector2D v, float angle) {
  float c = cos(angle);
  float s = sin(angle);
  return new Vector2D(v.x * c - v.y * s,  v.x * s + v.y * c);
}

Vector2D add(Vector2D u, Vector2D v) {
  return new Vector2D(u.x + v.x, u.y + v.y); 
}

Vector2D sub(Vector2D u, Vector2D v) {
  return new Vector2D(u.x - v.x, u.y - v.y);
}

Vector2D neg(Vector2D v) {
  return new Vector2D(-v.x, -v.y);
}

Vector2D lerp(Vector2D u, Vector2D v, float t) {
  return new Vector2D(u.x + t * (v.x - u.x), u.y + t * (v.y - u.y)); 
}

Point2D translate(Point2D p, Vector2D dir, float dist) {
  return new Point2D(p.x + dist * dir.x, p.y + dist * dir.y); 
}

Point2D translate(Point2D p, Vector2D dir) {
  return translate(p, dir, 1);
}

Point2D add(Point2D p, Point2D q) {
  return new Point2D(p.x + q.x, p.y + q.y); 
}

Point2D sub(Point2D p, Point2D q) {
  return new Point2D(p.x - q.x, p.y - q.y); 
}

Point2D scale(Point2D p, float s) {
  return new Point2D(p.x * s, p.y * s);
}

Point2D lerp(Point2D p, Point2D q, float t) {
  return new Point2D(p.x + t * (q.x - p.x), p.y + t * (q.y - p.y)); 
}

float dot(Vector2D u, Vector2D v) {
  //
  // this operator expresses combined magniture of two vectors factoring their same orientation 
  // u . v = norm(u) * norm(v) * cos(angle)u,v))
  // if u . v == 0 then u and v are orthogonal
  // if norm(u) = 1, then u.v = norm(v) * cos(angle)u,v)) which is the scalar projection of v onto u.
  //
  return u.x * v.x + u.y * v.y; 
}

float det(Vector2D u, Vector2D v) {
  //
  // this operato expresses the signed magnitude of the bivector (u,v), the parallelogram defined by vectors u and v. 
  // u x u = det(u,u) == 0
  //
  return u.x * v.y - v.x * u.y; 
}

float angle(Vector2D u, Vector2D v) {
  return atan2(det(u,v), dot(u,v));
}

float normSq(Vector2D v) { 
  return sq(v.x) + sq(v.y);
}

float norm(Vector2D v) { 
  return sqrt(normSq(v));
}

Vector2D scale(Vector2D v, float s) {
  return new Vector2D(v.x * s, v.y * s);
}

Vector2D spiralLerp(float t, Vector2D u, Vector2D v) {
  float m = norm(v) / norm(u);
  float a = angle(u,v);
  return scale(rotate(u, t * a), pow(m, t)); 
}

Point2D polygonCentroid(ArrayList<Point2D> points) {
  Point2D c = new Point2D(0, 0);
  int size = points.size();  
  for (int i = 0; i < size; ++i) {    
    Point2D P = points.get(i);
    c.x += P.x;
    c.y += P.y;   
  }     
  float is = 1.0 / size;
  c.x *= is;
  c.y *= is;  
  return c;
}

float polygonArea(ArrayList<Point2D> points) {  
  float area = 0; 
        
  int n = points.size();
  if (n > 1) {
    Point2D O = new Point2D(0, 0);
    Vector2D vp = new Vector2D(O, points.get(n-1));
    for (int i = 0; i < n; ++i) {
      Vector2D v = new Vector2D(O, points.get(i));
      area += det(vp, v); 
      vp = v;
    }
  }
  
  return abs(area / 2);
}

float triangleArea(Point2D p1, Point2D p2, Point2D p3) {
  Point2D O = new Point2D(0, 0);
  Vector2D vOP1 = new Vector2D(O, p1);
  Vector2D vOP2 = new Vector2D(O, p2);
  Vector2D vOP3 = new Vector2D(O, p3);
  float area = det(vOP1, vOP2) + det(vOP2, vOP3) + det(vOP3, vOP1);
  return abs(area / 2);
}

Point2D midPoint(Point2D pA, Point2D pB) {
  return new Point2D((pA.x + pB.x) / 2, (pA.y + pB.y) / 2);
}

// Return point M on line (A, B) given distance (A, M).
Point2D pointOnSegmentAB(Point2D pA, Point2D pB, float distAM) {
  // AM = t * AB and dist(AM) = t * dist(AB)
  float t = distAM / distance(pA, pB);
  Vector2D vAB = new Vector2D(pA, pB);
  return translate(pA, vAB, t);
}

// Return 'true' if line segment (p1,p2) contains p 
boolean segmentIntersectPoint(Point2D p1, Point2D p2, Point2D p) {
  float ds1 = distanceSq(p, p1);
  float ds2 = distanceSq(p, p2);
  float ds3 = distanceSq(p1, p2);
  return abs((ds1 + ds2) - ds3) < EPSILON;
}

// Return the lerp coefficients (t1,t2) for segments s1(A, B) and s2(C,D) where s1 and s2 intersects each other
float[] segmentIntersectionPoint(Point2D pA, Point2D pB, Point2D pC, Point2D pD) {
  //
  // AM = t1 * AB and CM = t2 * CD
  // given the basis origin O, AM = AO + OM and CM  = CO + OM
  // then OM = OA + t1 * AB = OC + t2 * CD
  // given u x u = det(u,u) = 0 
  // (OA + t1 * AB) x CD = (OC + t2 * CD) x CD
  // then OA x CD + t1 * (AB x CD) = OC x CD
  // t1 = (AC x CD) / (AB x CD)
  // (OA + t1 * AB) x AB = (OC + t2 * CD) x AB
  // then OA x AB = OC x AB + t2 * CD x AB
  // t2 = (CA x AB) / (CD x AB)
  // CD x AB = -(AB x CD)
  // t2 = (AC x AB) / (AB x CD)
  //
  
  Vector2D vAB = new Vector2D(pA, pB);
  Vector2D vCD = new Vector2D(pC, pD);
    
  //--
  float _det = det(vAB, vCD);
  if (abs(_det) < EPSILON) {
    return null; // parallel segments
  }
  
  Vector2D vAC = new Vector2D(pA, pC);
  float inv_det = 1.0 / _det;
  float t1 = det(vAC, vCD) * inv_det;
  float t2 = det(vAC, vAB) * inv_det;
  
  return new float[] {t1, t2};
}

//--
float NevilleLerp(float a, float b, float c, float ta, float tb, float tc, float t) {
  float x = a + (b - a) * t / (tb - ta);
  float y = b + (c - b) * t / (tc - tb);
  return x + (y - x) * t / (tc - ta);
}

//--
float NevilleLerp(float a, float b, float c, float t) {
  return NevilleLerp(a, b, c, 0.0f, 0.5f, 1.0f, t);
}

//--
Point2D NevilleLerp(Point2D pA, Point2D pB, Point2D pC, float t) {
  //--
  return new Point2D(
    NevilleLerp(pA.x, pB.x, pC.x, t), 
    NevilleLerp(pA.y, pB.y, pC.y, t)
    );
}

//--
Point2D calcNevilleDerivative(Point2D pA, Point2D pB, Point2D pC, float delta) {
  //--
  Point2D p = NevilleLerp(pA, pB, pC, 0.5f - delta/2);
  Point2D q = NevilleLerp(pA, pB, pC, 0.5f + delta/2);
  return scale(sub(q, p), 1 / delta);
}

//--
float calcNevilleDerivative(float a, float b, float c, float delta) {
  //--
  float p = NevilleLerp(a, b, c, 0.5f - delta/2);
  float q = NevilleLerp(a, b, c, 0.5f + delta/2);
  return (q - p) / delta;
}

// Return point M where lines l1(A, B) and l2(C,D) intersects each other
Point2D lineIntersectionPoint(Point2D pA, Point2D pB, Point2D pC, Point2D pD) {  
  float[] lerps = segmentIntersectionPoint(pA, pB, pC, pD);
  if (lerps == null) {
    return null;
  }
  return translate(pA, new Vector2D(pA, pB), lerps[0]);
}

// Return a subset of points in a list, given the start and end indices, as well as the advance step. 
<T> ArrayList<T> subset(ArrayList<T> list, int start, int end, int step) {
  ArrayList<T> result = new ArrayList<T>();
  for (int i = start, n = list.size(); i != end; i = pmod(i + step, n)) {
    result.add(list.get(i));   
  }  
  return result;
}

boolean isConvexPolygon(ArrayList<Point2D> polygon) {
  //
  // Traverse pairs of consective edges on the polygon and ensure that the angle between them is always less than 180 degree
  // that is, the angle direction is always clockwise or always counter-clockwise
  // that is, the sign of det(v1,v2) is the same
  //
  
  //--
  int n = polygon.size();
  if (n < 3) {
    return false;
  }

  int flag = 0;

  for (int i = 0; i < n; ++i) {   
    //--
    int j = (i + 1) % n;
    int k = (i + 2) % n;
    
    //--
    Point2D p0 = polygon.get(i);
    Point2D p1 = polygon.get(j);
    Point2D p2 = polygon.get(k);
     
    //--
    Vector2D vA = new Vector2D(p0, p1); 
    Vector2D vB = new Vector2D(p1, p2);
    float z  = det(vA, vB);
    if (z < 0) {
      flag |= 1;
    } else 
    if (z > 0) {
      flag |= 2;
    }
    if (flag == 3) {
      return false;
    }
  }
  
  return true;
}

// Return the two points that intersect a line and a circle
Point2D[] segmentCircleIntersection(Point2D p1, Point2D p2, Point2D C, float r) {
  //
  // AM = t x AB => AC + CM = t x AB => CM = t x AB + CA
  // dist(C,M) = r => CM.CM = r^2
  // then (t x AB + CA) . (t x AB + CA) = r^2
  //      (AB.AB)xt^2 + 2x(AB.CA)xt + CA.CA - r^2 = 0
  // Solve quadratic equation a*t^2+b*t+c = 0, then t = (-b+sqrt(delta)/2a, -b-sqrt(delta)/2a) where delta = b^2-4ac;
  //
  Vector2D vd = new Vector2D(p1, p2);
  Vector2D vf = new Vector2D(C, p1);
  
  float a = dot(vd, vd);
  float b = 2 * dot(vf, vd);
  float c = dot(vf, vf) - r * r;
  
  float delta = b*b - 4*a*c;
  if (delta < 0) {
    // there is no real solution, so no intersection
    return null;
  }
  
  float deltaSqrt = sqrt(delta);
  float inv_2a = 0.5 / a;
  float t1 = (-b - deltaSqrt) * inv_2a;
  float t2 = (-b + deltaSqrt) * inv_2a;
  
  if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
    return new Point2D[]{translate(p1, vd, t1), translate(p1, vd, t2)}; 
  }
  
  if (t1 >= 0 && t1 <= 1) {
    return new Point2D[]{translate(p1, vd, t1)}; 
  }
  
  if (t2 >= 0 && t2 <= 1) {
    return new Point2D[]{translate(p1, vd, t2)}; 
  }
  
  return null;
}

int convertHSVToRGB(int h, int s, int v) {
  float _h = clamp(h, 0, 360) / 60;
  float _s = clamp(s, 0, 100) / 100;
  float _v = clamp(v, 0, 100) / 100;    
  
  if (_s == 0) {
      int x = round(_v * 255);
      return (x << 16) | (x << 8) | x;
  }    
  
  int e = floor(_h);
  float l = _h - e;
  float d = _v * (1 - _s);
  float c = _v * (1 - _s * l);
  float n = _v * (1 - _s * (1 - l));
  
  float a, k, m;
  switch (e) {
  case 0: a = _v; k = n; m = d; break;
  case 1: a = c; k = _v; m = d; break;
  case 2: a = d; k = _v; m = n; break;
  case 3: a = d; k = c; m = _v; break;
  case 4: a = n; k = d; m = _v; break;
  default:
      a = _v; k = d; m = c;
  }
  
  int r = round(a * 255);
  int g = round(k * 255);
  int b = round(m * 255);
  
  return (r << 16) | (g << 8) | b;
}

Point2D mirrorPoint(Point2D p, Point2D midpoint) {
  return new Point2D(2 * midpoint.x - p.x, 2 * midpoint.y - p.y);
}

float[] computeBoundingBall(Point2D[] vertices) {
  //--
  float minx = MAX_FLOAT;
  float miny = MAX_FLOAT;
  float maxx = MIN_FLOAT;
  float maxy = MIN_FLOAT;
  //--
  for (Point2D p: vertices) {
    if (p.x < minx) minx = p.x;
    if (p.y < miny) miny = p.y;
    if (p.x > maxx) maxx = p.x;
    if (p.y > maxy) maxy = p.y;      
  }
  return new float[] {minx, maxx, miny, maxy};
}

// rendering ============================================================================

void vertex(Point2D p) {
  // emit vertex for shading or drawing
  vertex(p.x,p.y);
}

void line(Point2D p, Point2D q) {
  // draws edge (p,q)
  line(p.x,p.y,q.x,q.y); 
}

// LecturesInGraphics: utilities
// Colors, pictures, text formatting
// Author: Jarek ROSSIGNAC, last edited on September 10, 2012

// ************************************************************************ COLORS 
color black=#000000, white=#FFFFFF, // set more colors using Menu >  Tools > Color Selector
   red=#FF0000, green=#00FF01, blue=#0300FF, yellow=#FEFF00, cyan=#00FDFF, magenta=#FF00FB,
   grey=#818181, orange=#FFA600, brown=#B46005, metal=#B5CCDE, dgreen=#157901;
 
// ************************************************************************ GRAPHICS 
void pen(color c, float w) {stroke(c); strokeWeight(w);}
void circle(float x, float y, float r) {ellipse(x,y,r*2,r*2);}

// ************************************************************************ IMAGES & VIDEO 
int pictureCounter=0;
PImage myFace1, myFace2; // picture of author's face, should be: data/pic.jpg in sketch folder
int face_width = 80, face_height = 80;
void snapPicture() {saveFrame("PICTURES/P"+nf(pictureCounter++,3)+".jpg"); }

// ************************************************************************ TEXT 
Boolean scribeText=true; // toggle for displaying of help text
void scribe(String S, float x, float y) {fill(0); text(S,x,y); noFill();} // writes on screen at (x,y) with current fill color
void scribeHeader(String S, int i) {fill(0); text(S,10,20+i*20); noFill();} // writes black at line i
void scribeHeaderRight(String S) {fill(0); text(S,width-textWidth(S)-4,20); noFill();} // writes black on screen top, right-aligned
void scribeFooter(String S, int i) {fill(0); text(S,10,height-10-i*20); noFill();} // writes black on screen at line i from bottom
void scribeAtMouse(String S) {fill(0); text(S,mouseX,mouseY); noFill();} // writes on screen near mouse
void scribeMouseCoordinates() {fill(black); text("("+mouseX+","+mouseY+")",mouseX+7,mouseY+25); noFill();}
void displayHeader() { // Displays title and authors face on screen
    scribeHeader(title,0); 
    scribeHeader(instructions,1);
    scribeHeaderRight(name);     
    image(myFace1, width-2*(face_width+4)-20,25,face_width,face_height);
    image(myFace2, width-(face_width+4),25,face_width,face_height);
    }
void displayFooter() { // Displays help text at the bottom
    scribeFooter(guide,2); 
    scribeFooter(menu,1); 
    }
void scribeText(String S, float x, float y) { int l = S.length(); text(S,x-l*4,y+4); }

// ************************************************************************ 3D 

void drawFrame(float d) { 
  noStroke(); 
  fill(metal); 
  sphere(d/10);
  fill(blue);  
  drawArrow(d,d/10);
  fill(red); 
  pushMatrix(); 
  rotateY(PI/2); 
  drawArrow(d,d/10); 
  popMatrix();
  fill(green); 
  pushMatrix(); 
  rotateX(-PI/2); 
  drawArrow(d,d/10); 
  popMatrix();
}

void drawFan(float d, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
  vertex(0,0,d);
  for(float a=0; a<=TWO_PI+da; a+=da) {
    vertex(r*cos(a),r*sin(a),0);
  }
  endShape();
}

void drawCollar(float d, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
  for(float a=0; a<=TWO_PI+da; a+=da) {
    vertex(r*cos(a),r*sin(a),0); 
    vertex(rd*cos(a),rd*sin(a),d);
  }
  endShape();
}

void drawCone(float d, float r) {
  drawFan(d,r);  
  drawFan(0,r);
}

void drawStub(float d, float r, float rd) {
  drawCollar(d,r,rd); 
  drawFan(0,r);  
  pushMatrix(); 
  translate(0,0,d); 
  drawFan(0,rd); 
  popMatrix();
}
  
void drawArrow(float d, float r) {
  float dd=d/5;
  drawStub(d-dd,r*2/3,r/3); 
  pushMatrix(); 
  translate(0,0,d-dd); 
  drawCone(dd,r); 
  popMatrix();
}  
  
void drawBlock(float w, float d, float h, float x, float y, float z, float a) {
  pushMatrix(); 
  translate(x,y,h/2); 
  rotateZ(TWO_PI*a); 
  box(w,d,h); 
  popMatrix(); 
}

//************************ capturing frames for a movie ************************
boolean filming=false;  // when true frames are captured in FRAMES for a movie
int frameCounter=0;     // count of frames captured (used for naming the image files)
boolean change=false;   // true when the user has presed a key or moved the mouse
boolean animating=false; // must be set by application during animations to force frame capture

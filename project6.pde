/**************************** HEADER ****************************
 Class: CS6491 Fall 2014
 Students: Blaise Tine & Gupta Vibhav
 Project number: Number 6
 Project title: 3D Mesh Bend System
 Date of submission: Nov 11, 2014
*****************************************************************/

import java.util.*;

interface IControlView {
  
  void render();
  
  void onKeyPressed();
  
  void onKeyReleased();
  
  void onMouseMoved();

  void onMouseDragged();
  
  void onMousePressed();
  
  void onMouseReleased();
};  

//**************************** global variables ****************************

int fps = 30;
int fontsize = 12;
String[] models = {"data/model0.dat", "data/model1.dat"};
int model_index = 0;

Point3D g_center = null;
boolean fDragging = false;
int default_stroke_radius = 20;
int max_stroke_radius = 50;
int min_stroke_radius = 8;
int bspline_degree = 5;
int num_curve_points = 50;
//int num_spine_points = 30;
int num_spine_control_points = 7;
float vertex_radius = 6;
int num_sweep_faces = 12;
int num_morph_faces = 6;
float default_sweep_angle = PI;
float default_bend_radius = 400;

String text_msg = "";
int text_timer = 0;

Editor2D editor = new Editor2D();
Viewer3D viewer = new Viewer3D();
IControlView controlView = editor;

boolean displayLogging = true;

void log(String s) {
  if (displayLogging) {
    print(s);
  }
}

void logln(String s) {
  if (displayLogging) {
    println(s);
  }
}

void displayText(String msg) {
  text_msg = msg;
  text_timer = millis();
}

void renderText() {
  int fw = 16;
  int pad = 4;
  textSize(fw);
  float tw = textWidth(text_msg);
  float x0 = (width-tw)/2;
  float y0 = (height-fw)/2;
  noStroke();
  fill(blue);
  rect(x0-pad, y0-fw, tw+2*pad, fw+2*pad, pad);
  fill(white);
  text(text_msg, x0, y0, 0);
  int elapsed = millis() - text_timer;
  if (elapsed > 1000) {
    text_timer = 0;
  }
  textSize(fontsize);
}

//**************************** initialization ****************************
void setup() {           // executed once at the begining
  textFont(createFont("Courier", fontsize));
  size(600, 600, P3D);
  frameRate(fps);        // render 30 frames per second
  smooth();              // turn on antialiasing
  cursor(CROSS);
  
  myFace1 = loadImage("data/blaise.jpg");
  myFace2 = loadImage("data/vibhav.jpg");
  
  //--
  editor.init();
  viewer.init();
}

//**************************** display current frame ****************************
void draw() {        // executed at each frame
  //--
  background(white); // clear screen and paints white background
  pen(black,2);      // sets stroke color (to blackk) and width (to 3 pixels)
  
  // Render the UI
  displayHeader();
  if (scribeText && !filming) displayFooter(); // shows title, menu, and my face & name 
  
  // render model graph
  controlView.render();
  
  //--
  if (text_timer != 0) {
    renderText();
  }
  
  //--
  displayLogging = false;
    
  if (filming && (animating || change)) saveFrame("FRAMES/"+nf(frameCounter++,4)+".png");  
  change = false; // to avoid capturing frames when nothing happens
  // make sure that animating is set to true at the beginning of an animation and to false at the end
}  // end of draw()
  
//************************* mouse and key actions ****************************
void keyPressed() { // executed each time a key is pressed: the "key" variable contains the correspoinding char,
  //--
  int ascii = key - 0;
  //println("key="+key+" #"+ascii+")");
  
  // Keys mapping
  switch (key) {
  case '?':
    scribeText=!scribeText; // toggle display of help text and authors picture
    break;
  case '!':
    snapPicture(); // make a picture of the canvas
    break;
  case '~':
    filming=!filming; // filming on/off capture frames into folder FRAMES  
    break;
  case 'Q':
    exit();  // quit application
    break;
  case '3':
    if (controlView == editor) {
      if (editor.isReady()) {
        editor.processCurves();
        viewer.build();
        controlView = viewer;
      }
    } else {
      controlView = editor;
    }
    break; 
  }
  
  //--
  controlView.onKeyPressed();
  
  //--
  displayLogging = true;
  
  //--
  change=true;
}

void keyReleased() { // executed each time a key is released
  //--
  controlView.onKeyReleased(); 
  
  //--  
  change=true;
}

void mouseWheel(MouseEvent event) {
  //--
  viewer.dz -= event.getAmount(); 
  
  //--
  change=true;
}

void mouseMoved() { // executed when mouse is pressed and moved
  //--
  controlView.onMouseMoved(); 
  
  //--
  change=true;
}

void mouseDragged() { // when mouse is moved
  //--
  controlView.onMouseDragged();
    
  //--
  change=true;
}
  
void mousePressed() { // when mouse key is pressed
  //--
  controlView.onMousePressed();
  
  //--
  change=true;
}
  
void mouseReleased() { // when mouse key is released
  //--
  controlView.onMouseReleased();
  
  //--
  displayLogging = true;

  //--
  change=true;
}

//*************** text drawn on the canvas for name, title and help  *******************
String title = "CS6491, Fall 2014, Project 6: '3D Mesh Bend System";
String name = "Blaise Tine & Gupta Vibhav";
String instructions = "Draw two curves to describe a 3D model for bending\n";
String menu = "tab: show curves; c: reset; m: move vertex; s: save; l: load; w: wireframe\n?:(show/hide) help; !:snap picture; ~:toggle recording; Q:quit";
String guide = "Use the mouse to draw the curves control points"; // help info

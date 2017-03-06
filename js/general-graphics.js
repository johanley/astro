/* Tools of general utility for HTML5 canvas graphics.*/
var GRAPH = (function(){ 
  
  /**
   Set a standard font for a canvas.
   The intent is that this is called for all canvas objects.
   (Setting the canvas font can't be done in CSS.)
   If the browser has a hard minimum size of font, then it will take 
   precedence over any size specified here.
  */
  function myFont(ctx){
    //12 is chosen here since it's the default min for Chrome.
    ctx.font = '12px Courier New';
  }
  function myBigFont(ctx){
    ctx.font = '16px Courier New';
  }
  
  function text(ctx,string,x,y){
    ctx.fillText(string,x,y); 
  }
  
  function text(ctx,string,x,y,color){
    ctx.save();
    ctx.fillStyle = color; 
    ctx.fillText(string,x,y); 
    ctx.restore();
  }
  
  /** Single line segment. */
  function line(ctx, fromx, fromy, tox, toy){
    ctx.beginPath();
    ctx.moveTo(fromx,fromy);
    ctx.lineTo(tox,toy);
    ctx.stroke();
    ctx.closePath();
  }
  
  /** A polyline (N line segments.) */
  function lines(ctx, points){
    var i;
    ctx.beginPath();
    ctx.moveTo(points[0].x,points[0].y);
    for (i=1;i<points.length;i++){
      ctx.lineTo(points[i].x,points[i].y);
    }
    ctx.stroke();
    ctx.closePath();
  }
  
  /** Arrow-head up, with the point on the given xy. */
  function arrowUp(ctx, x, y){
    ctx.beginPath();
    ctx.moveTo(x-5,y+5); 
    ctx.lineTo(x,y);
    ctx.lineTo(x+5,y+5); 
    ctx.stroke();
    ctx.closePath();
  }
  
  /** Arrow-head down, with the point on the given xy. */
  function arrowDown(ctx, x, y){
    ctx.beginPath();
    ctx.moveTo(x-5,y-5); 
    ctx.lineTo(x,y);
    ctx.lineTo(x+5,y-5); 
    ctx.stroke();
    ctx.closePath();
  }
  
  /** Arrow-head to the right, with the point on the given xy. */
  function arrowRight(ctx, x, y){
    ctx.beginPath();
    ctx.moveTo(x-5,y-5); 
    ctx.lineTo(x,y);
    ctx.lineTo(x-5,y+5); 
    ctx.stroke();
    ctx.closePath();
  }
  
  /** Arrow-head to the left, with the point on the given xy. */
  function arrowLeft(ctx, x, y){
    ctx.beginPath();
    ctx.moveTo(x+5,y-5); 
    ctx.lineTo(x,y);
    ctx.lineTo(x+5,y+5); 
    ctx.stroke();
    ctx.closePath();
  }
  
  /** Filled triangle. */
  function triangle(ctx,ax,ay,bx,by,cx,cy){
    ctx.beginPath();
    ctx.moveTo(ax,ay);
    ctx.lineTo(bx,by);
    ctx.lineTo(cx,cy);
    ctx.lineTo(ax,ay);
    ctx.fill();
    ctx.stroke();
    ctx.closePath();
  }
  
  /** 
   Circle, filled in. Circles are less performant than rectangles.
   Firefox is currently the worst performer.
  */
  function spot(ctx,x,y,r){
    ctx.beginPath();
    ctx.arc(x,y,r, 0, Math.PI*2, false);
    ctx.stroke();
    ctx.fill();
    ctx.closePath();
  }
  
  /** A single pixel. Smaller than a spot.*/
  function point(ctx,x,y){
    ctx.beginPath();
    ctx.fillRect(x,y,1,1);
    ctx.closePath();
  }
  
  function square(ctx,x,y,size){
    ctx.beginPath();
    ctx.fillRect(x-size/2,y-size/2,size,size); //this avoids small shifts when a square morphs into a spot
    ctx.closePath();
  }
  
  /** Circle, unfilled. */
  function circle(ctx,x,y,r){
    ctx.beginPath();
    ctx.arc(x,y,r, 0, Math.PI*2, false);
    ctx.stroke();
    ctx.closePath();
  }
  
  function tickMarkVertical(ctx,x,y,size){
    line(ctx,x,y-size,x,y+size);
  }
  
  function tickMarkHorizontal(ctx,x,y,size){
    line(ctx,x-size,y,x+size,y);
  }
  
  /** Filled rectangle of the given color. */
  function rect(ctx,x,y,width,height,color){
    ctx.save();
    ctx.fillStyle=color;
    ctx.beginPath();
    ctx.fillRect(x,y,width,height);
    ctx.stroke();
    ctx.closePath();
    ctx.restore();
  }
  
  /** Unfilled rectangle of the given color. */
  function rectUnfilled(ctx,x,y,width,height,color){
    ctx.save();
    ctx.strokeStyle=color;
    ctx.beginPath();
    ctx.strokeRect(x,y,width,height);
    ctx.stroke();
    ctx.closePath();
    ctx.restore();
  }
  
  /**  
   Return the mouse position relative to the canvas. 
   Used in canvas event listeners.
  */
  function mousePosition(canvas,event) {
    var rect=canvas.getBoundingClientRect();
    return {
      x:event.clientX-rect.left,
      y:event.clientY-rect.top
    };
  }
  
  return {
   myFont:myFont,
   myBigFont:myBigFont,
   text:text,
   line:line,
   lines:lines,
   arrowUp:arrowUp,
   arrowDown:arrowDown,
   arrowRight:arrowRight,
   arrowLeft:arrowLeft,
   triangle:triangle,
   spot:spot,
   point:point,
   square:square,
   circle:circle,
   tickMarkVertical:tickMarkVertical,
   tickMarkHorizontal:tickMarkHorizontal,
   rect:rect,
   rectUnfilled:rectUnfilled,
   mousePosition:mousePosition
  };  

}()); // the top level function is invoked here; its return value is stored in UTIL, a global variable
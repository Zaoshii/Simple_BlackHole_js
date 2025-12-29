let particles = [];
let delp = false;
let moving = false;
let direction = "left";
let bhx;
let bhy;
let LmagMax = 0;
const bhmass = 8000;
const G = 2000;
const softening = 350;
let innerRadius = 400;
let outerRadius = 1050;



function createParticles() {
  let r = random(innerRadius,outerRadius);
  let theta = random(0, TWO_PI);
  let x = bhx + r*cos(theta);
  let y = bhy + r*sin(theta);
  let ux = -sin(theta);
  let uy = cos(theta);
  let v = sqrt(G*bhmass/r);
  let vx = ux*v;
  let vy = uy*v;
  return {
    x,y,
    vx,vy,
    ax:0,
    ay:0,
    mass: random(5,50),
  };
}

function getDistance(x1,y1) {
  return sqrt( (x1 - bhx) ** 2 + (y1 - bhy) ** 2 )
}

function getAngularMomentum(p) {
  let rx = p.x - bhx;
  let ry = p.y - bhy;
  let vx = p.vx;
  let vy = p.vy;

  let L = p.mass * (rx*vy - ry*vx);
  return L;
}

function applyForces(p, dt) {
  let distance = getDistance(p.x,p.y);

  distance = max(distance, 1);

  let dx = (bhx - p.x) / distance;
  let dy = (  bhy - p.y) / distance;

  let F = G * p.mass * bhmass / ( distance ** 2 + softening);

  p.ax = dx * F / p.mass;
  p.ay = dy * F / p.mass;

  p.vx += p.ax * dt;
  p.vy += p.ay * dt;

  p.x += p.vx * dt;
  p.y += p.vy * dt;

  p.ax = 0;
  p.ay = 0;

}

function accelerationAt(x,y,mass) {
  let dx = bhx - x ;
  let dy = bhy - y ;
  let r = sqrt( dx*dx + dy*dy );
  r = max(r,1);

  dx /= r;
  dy /= r;

  let F = ( G * mass * bhmass ) / ( r*r + softening );

  return {
    ax: dx * F / mass,
    ay: dy * F / mass,
  };


}

function RK4Integrate(p,dt) {
  
  let a1 = accelerationAt(p.x,p.y,p.mass);
  let k1 = {
    dx: p.vx,
    dy: p.vy,
    dvx: a1.ax,
    dvy: a1.ay
  };

  let a2 = accelerationAt(
    p.x + k1.dx * dt/2,
    p.y + k1.dy * dt/2,
    p.mass
  );
  
  let k2 = {
    dx: p.vx + k1.dvx * dt/2,
    dy: p.vy + k1.dvy * dt/2,
    dvx: a2.ax,
    dvy: a2.ay
  };

  let a3 = accelerationAt(
    p.x + k2.dx * dt/2,
    p.y + k2.dy * dt/2,
    p.mass
  );

  let k3 = {
    dx: p.vx + k2.dvx * dt/2,
    dy: p.vy + k2.dvy * dt/2,
    dvx: a3.ax,
    dvy: a3.ay
  };

  let a4 = accelerationAt(
    p.x + k3.dx * dt,
    p.y + k3.dy * dt,
    p.mass
  );

  let k4 = {
    dx: p.vx + k3.dvx * dt,
    dy: p.vy + k3.dvy * dt,
    dvx: a4.ax,
    dvy: a4.ay
  };

  p.x += dt/6 * ( k1.dx + 2*k2.dx + 2*k3.dx + k4.dx );
  p.y += dt/6 * ( k1.dy + 2*k2.dy + 2*k3.dy + k4.dy );

  p.vx += dt/6 * ( k1.dvx + 2*k2.dvx + 2*k3.dvx + k4.dvx );
  p.vy += dt/6 * ( k1.dvy + 2*k2.dvy + 2*k3.dvy + k4.dvy );
}

function reset() {
  innerRadius = 400;
  outerRadius = 1050;
  bhx = windowWidth/2;
  bhy = windowHeight/2;
  for (let i = 0; i < 1500; i++) {
    let randomParticle = createParticles();
    particles.push(randomParticle)
    console.log(particles[i]);
  } 
  for (let i = 0; i < 50; i ++) {
    innerRadius = 200;
    outerRadius = 250;
    let randomParticle = createParticles();
    particles.push(randomParticle)
    console.log(particles[i]);
  }
}

function setup() {
  createCanvas(windowWidth, windowHeight);
  colorMode(HSB,255);
  reset();
  for (let p of particles) {
    let L_ = abs(getAngularMomentum(p));
    if (L_ > LmagMax) LmagMax = L_;
  }
}

function draw() {
  background(0, 30);
  let dt = deltaTime/2000;
  fill(255);
  if (mouseIsPressed == true) {
    bhx = mouseX;
    bhy = mouseY;
  }
  else {
    bhx = windowWidth/2;
    bhy = windowHeight/2;
    circle(bhx,bhy,150);
  }

  for (let i = particles.length - 1; i >=0; i--) {
    let p = particles[i];
    RK4Integrate(p,dt);
    let L = getAngularMomentum(p);
    let Lmag = abs(L); 
    if (frameCount%120 == 0) {
      for (let p of particles) {
        let L_ = abs(getAngularMomentum(p));
        if (L_ > LmagMax) LmagMax = L_;
        }
      }
  console.log(LmagMax); 
    let Huevalue = map(Lmag,0,LmagMax,0,240);
    Huevalue = constrain(Huevalue,0,240);
    stroke(Huevalue, 200, 255);
    if (getDistance(p.x,p.y) < 10 && delp == true) {
      particles.splice(i,1);
      continue;
    }
    strokeWeight(p.mass * 0.1);
    point(p.x,p.y);
  }
  
  if (moving) {
    if (direction == "right") {
      bhx++;
    }
    if (direction == "left") {
      bhx--;
    }
    if (bhx > windowWidth) {
      direction = "left";
    }
    if (bhx < 0) {
      direction = "right";
    }
  }
}

function keyPressed() {
  if (key == "r") {
    particles = [];
    reset();
  }
}

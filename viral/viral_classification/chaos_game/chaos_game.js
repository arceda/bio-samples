let points = []
let current;
let previous;
let percent = 0.5;
let n_vertex = 5;

function setup(){

    createCanvas(windowWidth,windowHeight);
    points = [];    

    for (let i = 0; i < n_vertex; i++){
        //the vector point surrond a circle
        let angle   = i * TWO_PI / n_vertex;
        let v = p5.Vector.fromAngle(angle);
        v.mult(width/2);
        v.add(width/2, height/2)

        //let v = createVector(random(width), random(height))
    
        points.push(v)
    }


    reset();   

}

function reset(){   
    current = createVector(random(width), random(height))

    background(0);
    stroke(255);
    strokeWeight(8);

    for (let p of points){
        point(p.x, p.y)
    }
}

function draw(){
    if (frameCount % 200 == 0){
        reset();
    }

    for (let i = 0; i < 1000; i++){
        strokeWeight(2)
        stroke(255,0,255, 100);

        let next = random(points)

        
        // for pentagon
        
        if (n_vertex == 5 && next !== previous){
            current.x = lerp(current.x, next.x, percent);
            current.y = lerp(current.y, next.y, percent);
            point(current.x, current.y)
        }
        else if (n_vertex != 5){
            current.x = lerp(current.x, next.x, percent);
            current.y = lerp(current.y, next.y, percent);
            point(current.x, current.y)
        }

        

        previous  = next;

    }
    
}
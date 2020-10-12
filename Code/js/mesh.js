// In this file you will implement traversal and analysis for your assignment.
// Make sure to familiarize yourself with the utility functions in meshUtils.js
// they might be useful for the second part of your assignment!

////////////////////////////////////////////////////////////////////////////////
// Traversal
////////////////////////////////////////////////////////////////////////////////

// Return all vertices on face f
Mesh.prototype.verticesOnFace = function(face) {
  const vertices = [];
  let run = face.halfedge;
  const first = run;
  while (true) {
    vertices.push(run.vertex);
    run = run.next;
    if (run === first) {
      break;
    }
  }
  return vertices;
};

// Return all halfedges on face f
Mesh.prototype.edgesOnFace = function(face) {
  var halfedges = [];
  var run = face.halfedge;
  var first = run;
  while (true) {
      halfedges.push(run);
      run = run.next;
      if (run === first ) {
        break;
      }
  }
  return halfedges;
};

// Return all faces adjacent to input face, not
// including input face.
Mesh.prototype.facesOnFace = function(face) {

  var faces = [];
  var run = face.halfedge;
  var first = run;
  while ( true ) {
    faces.push( run.opposite.face );
    run = run.next;
    if (run === first) {
       break;
    }
  }

  return faces;
};

// Return one-ring neighbors of input vertex, not
// including the input vertex itself
Mesh.prototype.verticesOnVertex = function(v) {
  var vertices = [];
  var run = v.halfedge;
  var first = run;
  while (true) {
    vertices.push(run.vertex);
    run = run.opposite.next;
    if (run === first ) {
      break;
    }
  }

  return vertices;
};

// Return all halfedges that point away from v
Mesh.prototype.edgesOnVertex = function(v) {
  var halfedges = [];
  var run = v.halfedge;
  var first = run;
  while (true) {
    halfedges.push( run );
    run = run.opposite.next;
    if (run === first) {
      break;
    }
  }
  return halfedges;
};

// Return all faces that include v as a vertex.
Mesh.prototype.facesOnVertex = function(v) {
  var faces = [];

  var run = v.halfedge;
  var first = run;
  while ( true ) {
    faces.push( run.face );
    run = run.opposite.next;
    if ( run === first ) {
      break;
    }
  }
  return faces;
};

// Return the vertices that form the endpoints of a given edge
Mesh.prototype.verticesOnEdge = function(e) {
  const vertices = [];

  faces.push( e.vertex );
  faces.push( e.opposite.vertex );

  return vertices;
};

// Return the faces that include a given edge
Mesh.prototype.facesOnEdge = function(e) {
  const faces = [];
  faces.push( e.face );
  faces.push( e.opposite.face );
  return faces;
};

// Return the edge pointing from v1 to v2
Mesh.prototype.edgeBetweenVertices = function(v1, v2) {
  let outEdge = undefined;
  edges = this.edgesOnVertex(v1)
  var index = 0;
  var run = edges[index]
  var first = run;
  while (true) {
    if (run.vertex === v2) { outEdge = run }
    index = index + 1
    run = edges[index]
    if (index > edges.length - 1) {
      break;
    }
  }
  return outEdge;
};

////////////////////////////////////////////////////////////////////////////////
// Analysis
////////////////////////////////////////////////////////////////////////////////

// Return the surface area of a provided face f.
Mesh.prototype.calculateFaceArea = function(f) {
  let area = 0.0;
  var he = f.halfedge;
  var v1 = he.vertex;
  var v2 = he.next.vertex;
  var v3 = he.next.next.vertex;
  var first = he.vertex;
  var second = he.next.vertex;
  while (true) {
    var a = this.dist(v1.position, v2.position);
    var b = this.dist(v2.position, v3.position);
    var c = this.dist(v1.position, v3.position);
    var s = 0.5 * (a + b + c);
    var newArea = Math.sqrt(s * (s - a) * (s - b) * (s - c));
		area = area + newArea;
    he = he.next;
    v1 = he.vertex;
    v2 = he.next.vertex;
    v3 = he.next.next.vertex;
    if ( v1 === first || v1 === second) {
      break;
    }
  }

  f.area = area;
  return area;
};

// Update the area attributes of all faces in the mesh
Mesh.prototype.calculateFacesArea = function() {
  for (let i = 0; i < this.faces.length; ++i) {
    this.calculateFaceArea(this.faces[i]);
  }
};

// Calculate the vertex normal at a given vertex,
// using the face normals of bordering faces, weighted by face area
Mesh.prototype.calculateVertexNormal = function(v) {
  var v_faces = Mesh.prototype.facesOnVertex(v);
	var v_normal = new THREE.Vector3( 0, 0, 0 );
	for (var i = 0; i < v_faces.length; ++i) {
		var temp = new THREE.Vector3(v_faces[i].normal.x, v_faces[i].normal.y, v_faces[i].normal.z);
		temp.multiplyScalar(Mesh.prototype.calculateFaceArea(v_faces[i]));
        v_normal.add(temp);;
  }
	v_normal = v_normal.normalize();
	
  return v_normal;
};

// update the vertex normals of every vertex in the mesh
Mesh.prototype.updateVertexNormals = function() {
  for (let i = 0; i < this.vertices.length; ++i) {
    this.calculateVertexNormal(this.vertices[i]);
  }
};

// compute the average length of edges touching v
Mesh.prototype.averageEdgeLength = function(v) {
  let avg = 0.0;

  var edges = this.edgesOnVertex(v);
	var lenSum = 0;

	for (var i = 0; i < edges.length; ++i) {
		var edgeLen = this.dist(edges[i].vertex.position, edges[i].opposite.vertex.position);
		lenSum = lenSum + edgeLen;
	}

	if (edges.length != 0) {
		avg = lenSum/edges.length;
	}

  return avg;
};

////////////////////////////////////////////////////////////////////////////////
// Topology
////////////////////////////////////////////////////////////////////////////////

// Given a face in the shape of an arbitrary polygon,
// split that face so it consists only of several triangular faces. 
Mesh.prototype.triangulateFace = function(f) {
  var he = f.halfedge;
  var first = he.vertex;
  var second = he.next.vertex;
  var v1 = he.next.vertex;
  var v2 = he.next.next.vertex;
  var v3 = he.next.next.next.vertex;
  while (true) {
    if (v3 === first || v2 === first || v1 === first) {
      break;
    }
    this.splitFaceMakeEdge (f, first, v2, second, true)
    he = he.next;
    var v1 = he.vertex;
    var v2 = he.next.vertex;
    var v3 = he.next.next.vertex;
  }
};

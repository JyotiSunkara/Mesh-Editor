// In this file you will implement traversal and analysis for your assignment.
// Make sure to familiarize yourself with the utility functions in meshUtils.js
// they might be useful for the second part of your assignment!

////////////////////////////////////////////////////////////////////////////////
// Traversal
////////////////////////////////////////////////////////////////////////////////

// Return all vertices on face f
Mesh.prototype.verticesOnFace = function(f) {
  const vertices = [];
  let runEdges = f.halfedge;
  const first = runEdges;
  while (true) {
    vertices.push(runEdges.vertex);
    runEdges = runEdges.next;
    if (runEdges === first) {
      break;
    }
  }
  return vertices;
};

// Return all halfedges on face f
Mesh.prototype.edgesOnFace = function(f) {
  var halfedges = [];
  var runEdges = f.halfedge;
  var first = runEdges;
  while (true) {
      halfedges.push(runEdges);
      runEdges = runEdges.next;
      if (runEdges === first ) {
        break;
      }
  }

  return halfedges;
};

// Return all faces adjacent to input face, not
// including input face.
Mesh.prototype.facesOnFace = function(f) {

  var faces = [];
  var runEdges = f.halfedge;
  var first = runEdges;
  while ( true ) {
    faces.push(runEdges.opposite.face);
    runEdges = runEdges.next;
    if (runEdges === first) {
       break;
    }
  }

  return faces;
};

// Return one-ring neighbors of input vertex, not
// including the input vertex itself
Mesh.prototype.verticesOnVertex = function(v) {
  var vertices = [];
  var runEdges = v.halfedge;
  var first = runEdges;
  while (true) {
    vertices.push(runEdges.vertex);
    runEdges = runEdges.opposite.next;
    if (runEdges === first ) {
      break;
    }
  }

  return vertices;
};

// Return all halfedges that point away from v
Mesh.prototype.edgesOnVertex = function(v) {
  var halfedges = [];
  var runEdges = v.halfedge;
  var first = runEdges;
  while (true) {
    halfedges.push(runEdges);
    runEdges = runEdges.opposite.next;
    if (runEdges === first) {
      break;
    }
  }
   return halfedges;
};

// Return all faces that include v as a vertex.
Mesh.prototype.facesOnVertex = function(v) {
  var faces = [];

  var runEdges = v.halfedge;
  var first = runEdges;
  while (true) {
    faces.push(runEdges.face);
    runEdges = runEdges.opposite.next;
    if (runEdges === first) {
      break;
    }
  }
  return faces;
};

// Return the vertices that form the endpoints of a given edge
Mesh.prototype.verticesOnEdge = function(e) {
  const vertices = [];
  faces.push(e.vertex);
  faces.push(e.opposite.vertex);
  return vertices;
};

// Return the faces that include a given edge
Mesh.prototype.facesOnEdge = function(e) {
  const faces = [];
  faces.push(e.face);
  faces.push(e.opposite.face);
  return faces;
};

// Return the edge pointing from v1 to v2
Mesh.prototype.edgeBetweenVertices = function(v1, v2) {
  let outEdge = undefined;
  edges = this.edgesOnVertex(v1)
  var index = 0;
  var runEdges = edges[index]
  var first = runEdges;
  while (true) {
    if (runEdges.vertex === v2) { 
      outEdge = runEdges 
    }
    index = index + 1
    runEdges = edges[index]
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
    var area = 0.0;
    var runEdge, first, second;
    var v1, v2, v3;

    // Herons Formula for triangle
    var a, b, c, s;

    
    runEdge = f.halfedge;
    first = runEdge.vertex;
    second = runEdge.next.vertex;

    v1 = runEdge.vertex;
    v2 = runEdge.next.vertex;
    v3 = runEdge.next.next.vertex;


    while (true) {
        a = this.dist(v1.position, v2.position);
        b = this.dist(v2.position, v3.position);
        c = this.dist(v1.position, v3.position);
        s = 0.5 * (a + b + c);
        var newArea = Math.sqrt(s * (s - a) * (s - b) * (s - c));
		    
        area = area + newArea;
        runEdge = runEdge.next;
        v1 = runEdge.vertex;
        v2 = runEdge.next.vertex;
        v3 = runEdge.next.next.vertex;
        if (v1 === first || v1 === second) {
          break;
        }
    }
    return area;
};

// Update the area attributes of all faces in the mesh
Mesh.prototype.calculateFacesArea = function() {
    for ( var i = 0; i < this.faces.length; ++i ) {
      this.faces[i].area = this.calculateFaceArea(this.faces[i]);
    }
};

// Calculate the vertex normal at a given vertex,
// using the face normals of bordering faces, weighted by face area
Mesh.prototype.calculateVertexNormal = function(v) {
  var vFaces = Mesh.prototype.facesOnVertex(v);
	var vNormal = new THREE.Vector3(0, 0, 0);
	for (var i = 0; i < vFaces.length; ++i) {
		var normal = new THREE.Vector3(vFaces[i].normal.x, 
                                   vFaces[i].normal.y, 
                                   vFaces[i].normal.z);

		normal.multiplyScalar(Mesh.prototype.calculateFaceArea(vFaces[i]));
    vNormal.add(normal);
  }

	vNormal = vNormal.normalize();
  return vNormal;
};

// update the vertex normals of every vertex in the mesh
Mesh.prototype.updateVertexNormals = function() {
   for ( var i = 0; i < this.vertices.length; ++i ) {
      this.vertices[i].normal = this.calculateVertexNormal(this.vertices[i]);
   }
};



// compute the average length of edges touching v
Mesh.prototype.averageEdgeLength = function ( v ) {
  var avg = 0.0;

	var edges = this.edgesOnVertex(v);
	var lenSum = 0;
	for (let i = 0; i < edges.length; ++i) {
		var edgeLen = this.dist(edges[i].vertex.position, 
                            edges[i].opposite.vertex.position);
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
  var runEdge = f.halfedge;
  var first = runEdge.vertex;
  var second = runEdge.next.vertex;
  var v1, v2, v3;
  v1 = runEdge.next.vertex;
  v2 = runEdge.next.next.vertex;
  v3 = runEdge.next.next.next.vertex;
  while (true) {
    if (v1 === first || v2 === first || v3 === first) {
      break;
    }
    this.splitFaceMakeEdge(f, first, v2, second, true)
    runEdge = runEdge.next;
    v1 = runEdge.vertex;
    v2 = runEdge.next.vertex;
    v3 = runEdge.next.next.vertex;
  }
};



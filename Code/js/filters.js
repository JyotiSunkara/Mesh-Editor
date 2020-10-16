var Filters = Filters || {};

// Space for your helper functions
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 105 lines of code.
// ----------- STUDENT CODE END ------------

// Translate all selected vertices in the mesh by the given x,y,z offsets.
Filters.translation = function(mesh, x, y, z) {
  const T = new THREE.Vector3(x, y, z);

  const verts = mesh.getModifiableVertices();
  const N = verts.length;
  for (var i = 0; i < N; ++i) {
    verts[i].position.add(T);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Given x,y,z, the desired rotation around each axis, in radians,
// apply this rotation to all selected vertices in the mesh.
Filters.rotation = function(mesh, x, y, z) {
  const verts = mesh.getModifiableVertices();
  const N = verts.length;
  
  for (var i = 0; i < N; i++){
    if (x !== 0){
      var axis = new THREE.Vector3(1, 0, 0);
      verts[i].position.applyAxisAngle(axis, x);
    }

    if (y !== 0) {
      var axis = new THREE.Vector3(0, 1, 0);
      verts[i].position.applyAxisAngle(axis, y);
    }

    if (z !== 0) {
      var axis = new THREE.Vector3(0, 0, 1);
      verts[i].position.applyAxisAngle(axis, z);
    }
  }

  mesh.updateNormals();
  mesh.calculateFacesArea();

};

// Uniformly scale the position of all selected vertices in the mesh
// by the provided scale factor s
Filters.scale = function(mesh, s) {
  const verts = mesh.getModifiableVertices();

  const N = verts.length;
  for (var i = 0; i < N; ++i) {
    verts[i].position.multiplyScalar(s);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// estimate the per-vertex gaussian vurvature of the mesh at each vertex.
// set that vertex's color to some value based on its curvature value.
// (the precise mapping of curvature to color is left to you)
Filters.curvature = function(mesh) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 102 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Curvature is not implemented yet");
};

// Apply a random offset to each selected vertex in the direction of its normal
// scale the random offset by the provided factor and by
// the average length of edges at that vertex
Filters.noise = function(mesh, factor) {
  var verts = mesh.getModifiableVertices();
  var N = verts.length;
  for ( var i = 0 ; i < N ; ++i ) {
    verts[i].position.addScalar(mesh.averageEdgeLength(verts[i]) * factor * (Math.random() * 2 - 1));
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Smooth the mesh using the specified weighting scheme.
// In the standard case, this is done using uniform Laplacian smoothing,
// by moving each vertex towards the average position of its neighbors.
//
// Arguments:
//  - mesh: the mesh to smooth
//  - iter: the number of iterations of smoothing to apply
//  - delta: a scaling factor for the amount of smoothing
//  - curvFlow: a bool. if true, use cotangent weights instead of uniform (requires triangular mesh)
//  - scaleDep: a bool. if true, scale offsets differently for each vertex (see spec.)
//  - implicit: a bool. if true, perform implicit smoothing (see spec.)
//
// Note that the reference solution calls a giant utility function so the line
// count is not terribly representative of the true solution
//
// For implicit, you will want to compute I - M*L*delta, where I is the identity
// matrix, M is a diagonal "mass" matrix, and L is a Laplacian matrix. Then
// you will want to call math.lup() on your result in order to decompose the
// matrix. Finally, call math.lusolve() to compute the X,Y, and Z positions of
// vertices. Note that the decomposition step allows for fast solving multiple
// times. It would be possible to replace a few of these steps with simple matrix
// inversion; however, matrix inversion is far slower and less numerically stable
//
Filters.smooth = function(mesh, iter, delta, curvFlow, scaleDep, implicit) {
  const verts = mesh.getModifiableVertices(); 

  if (curvFlow){
      this.triangulate(mesh);
  }   
  var neighbors = [];
  for(var i = 0; i < verts.length; i++) {
      neighbors.push(mesh.verticesOnVertex(verts[i]));
  }
  var he, v1, v2, edge2, new1;
  var cosAlpha, sinAlpha, cotAlpha, cosBeta, sinBeta, cotBeta, w;
  var vClone, nClone;
  var sumNeighbors = new THREE.Vector3(0, 0, 0)
  var weightsSum = 0;
  while(iter--) {
      var vertexArea = [], sumArea = 0;
      var newVerts = [];
      for(var i = 0; i < verts.length; i++) {
        sumNeighbors = new THREE.Vector3(0, 0, 0)
        weightsSum = 0;
        vClone = verts[i].position.clone();
        newVerts.push(vClone);
        for(j in neighbors[i]) {
            nClone = neighbors[i][j].position.clone();
            if(!curvFlow) { 
              sumNeighbors.add(nClone);
              weightsSum++;  
              
            } 
            
            if(curvFlow) {
              he = mesh.edgeBetweenVertices(verts[i], neighbors[i][j]);
              v1 = he.next.vertex.position.clone(); 
              v2 = he.opposite.next.vertex.position.clone();
              edge2 = nClone.sub(v1);
              new1 = vClone.sub(v1).clone();

              cosAlpha = new1.dot(edge2);
              sinAlpha = new1.cross(edge2).length();
              cotAlpha = cosAlpha / sinAlpha;

              cosBeta = new1.dot(edge2);
              sinBeta = new1.cross(edge2).length();
              cotBeta = cosBeta / sinBeta;

              w = 0.5 * (cotAlpha + cotBeta);
              sumNeighbors.add(nClone.multiplyScalar(w));
              weightsSum = weightsSum + w;
              
            }
        }

        faces = mesh.facesOnVertex(verts[i]);
        var associatedArea = 0
        for(var l = 0; l < faces.length; l++) {
          if(scaleDep)
            associatedArea += mesh.calculateFaceArea(faces[l]);
        }

        if(!scaleDep) {
          vertexArea.push(1);
          sumArea++;

        } 
        
        if(scaleDep) {
          vertexArea.push(associatedArea);
          sumArea = sumArea + associatedArea;
        }

        newVerts[i].add(sumNeighbors.add(verts[i].position.clone().negate().multiplyScalar(weightsSum)).multiplyScalar(delta / vertexArea[i]));
      }

      for(var i = 0; i < newVerts.length; i++) {
        var avgArea = sumArea / newVerts.length;
        newVerts[i].sub(verts[i].position).multiplyScalar(avgArea).add(verts[i].position);
        verts[i].position.set(newVerts[i].x, newVerts[i].y, newVerts[i].z);
      }
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Sharpen the mesh by moving selected vertices away from the average position
// of their neighbors (i.e. Laplacian smoothing in the negative direction)
Filters.sharpen = function(mesh, iter, delta) {

	Filters.smooth(mesh, iter, -delta);
	
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Move every selected vertex along its normal direction
// Scale the amount by the provided factor and average edge length at that vertex
Filters.inflate = function (mesh, factor) {

  var verts = mesh.getModifiableVertices();
	var N = verts.length;
  for ( var i = 0 ; i < N ; ++i ) {
		var avgLen = mesh.averageEdgeLength(verts[i]);
		var normal = new THREE.Vector3(verts[i].normal.x, verts[i].normal.y, verts[i].normal.z);    
    verts[i].position.add(normal.multiplyScalar(factor * avgLen));
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// rotate selected vertices around the Y axis by an amount
// proportional to its Y value times the scale factor.
Filters.twist = function(mesh, factor) {
  
  var verts = mesh.getModifiableVertices();
  var N = verts.length;

  var vx, vy, vz, rot;
	for ( var i = 0 ; i < N ; ++i ) { 
		vx = verts[i].position.getComponent(0);
		vy = verts[i].position.getComponent(1);
		vz = verts[i].position.getComponent(2);
		rot = vy * factor;
    verts[i].position = new THREE.Vector3(Math.cos(rot) * vx + Math.sin(rot) * vz, 
                                          vy,
                                          -Math.sin(rot) * vx + Math.cos(rot) * vz);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// warp a mesh using a nonlinear mapping of your choice
Filters.wacky = function(mesh, factor) {
	var verts = mesh.getModifiableVertices();
	var N = verts.length;

  for ( var i = 0 ; i < N ; ++i ) {
		  var shear = new THREE.Vector3(Math.sin(((verts[i].position.y) * 22)/(7 * factor)), 
                                    Math.cos(((verts[i].position.y) * 22)/(7 * factor)), 
                                    0)
      verts[i].position.add(shear);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Convert the selected faces from arbitrary polygons into all triangles
Filters.triangulate = function(mesh) {
  const faces = mesh.getModifiableFaces();
  var N = faces.length;
  for ( var i = 0 ; i < N ; i++ ) { 
        mesh.triangulateFace(faces[i]);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitEdgeMakeVert in mesh.js
Filters.splitEdge = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 2) {
    mesh.splitEdgeMakeVert(verts[0], verts[1], 1/2);
  } else {
    console.log("ERROR: to use split edge, select exactly 2 adjacent vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinEdgeKillVert in mesh.js
Filters.joinEdges = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 3) {
    var v0 = verts[0],
      v1 = verts[1],
      v2 = verts[2];

    const he01 = mesh.edgeBetweenVertices(v0, v1);
    const he12 = mesh.edgeBetweenVertices(v1, v2);

    if (he01) {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[1], verts[2]);
      } else {
        mesh.joinEdgeKillVert(verts[1], verts[0], verts[2]);
      }
    } else {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[2], verts[1]);
      } else {
        console.log(
          "ERROR: to use join edge, select exactly 3 vertices such that one only has edges to the other two"
        );
      }
    }
  } else {
    console.log("ERROR: to use join edge, select exactly 3 vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitFaceMakeEdge in mesh.js
Filters.splitFace = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 1) {
    mesh.splitFaceMakeEdge(faces[0], verts[0], verts[1]);
  } else {
    console.log("ERROR: to use split face, select exactly 1 face and 2 nonadjacent vertices on it");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinFaceKillEdge in mesh.js
Filters.joinFaces = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 2) {
    mesh.joinFaceKillEdge(faces[0], faces[1], verts[0], verts[1]);
  } else {
    console.log(
      "ERROR: to use split face, select exactly 2 adjacent faces the 2 vertices between them"
    );
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// extrude the selected faces from the mesh in the direction of their normal
// vector, scaled by the provided factor.
// See the spec for more detail.
Filters.extrude = function(mesh, factor) {
  const faces = mesh.getModifiableFaces();
  var N = faces.length;

  for(var i = 0; i < N ; i++) {
    var newVerts = [];
		var oldVerts = mesh.verticesOnFace(faces[i]);
    var oldEdges = mesh.edgesOnFace(faces[i]);
    
    for(var j = 0; j < oldEdges.length; j++) {
      var o1 = oldEdges[j].vertex;
      var o2 = oldEdges[j].opposite.vertex;
      var newVert = mesh.splitEdgeMakeVert(o1, o2, 0);

      facesOne = mesh.facesOnVertex(o1);
      facesTwo = mesh.facesOnVertex(o2);
      var commonFace = undefined;
      for (var fi = 0; fi < facesOne.length; fi++) {
        if (facesOne[fi] !== faces[i]) {
          for (var fj = 0; fj < facesTwo.length; fj++) {
                  if (facesOne[fi] === facesTwo[fj])
                      commonFace = facesOne[fi];
              }
        }
    }
      
      mesh.splitFaceMakeEdge(commonFace, newVert.halfedge.vertex, newVert.halfedge.opposite.next.vertex);
      newVerts.push(newVert);
    }
    
    for(var j = 0; j < newVerts.length; j++) {
      
      mesh.splitFaceMakeEdge(faces[i], newVerts[j], newVerts[(j + 1) % newVerts.length]);
      var f1 = newVerts[j].halfedge.opposite.next.next.opposite.face;
      var f2 = newVerts[j].halfedge.opposite.face;
      var v1 = newVerts[j].halfedge.opposite.next.vertex;
      var v2 = newVerts[j].halfedge.vertex;
      mesh.joinFaceKillEdge(f1, f2, v1, v2);
    }
      
    for(var j = 0; j < newVerts.length; j++)  {
        newVerts[j].position.add(faces[i].normal.normalize().multiplyScalar(factor));
    }
  }
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Truncate the selected vertices of the mesh by "snipping off" corners
// and replacing them with faces. factor specifies the size of the truncation.
// See the spec for more detail.
Filters.truncate = function (mesh, factor) {

	origMesh = new Mesh();
	origMesh.copy(mesh);
	var origVerts = origMesh.getModifiableVertices(), verts = mesh.getModifiableVertices(), otherVerts, newVerts;
	var N = verts.length, adjFace;
	for (var i = 0 ; i < N ; ++i ) {

		otherVerts = mesh.verticesOnVertex(verts[i]);
		newVerts = []
		for (var j = 0; j < otherVerts.length-1; ++j) {
			newVerts[j] = mesh.splitEdgeMakeVert(verts[i], otherVerts[j], factor);
		}
			
		for (var j = 0; j < newVerts.length - 1; j++) {
      var v1 = newVerts[j];
      var v2 = newVerts[(j + 1) % newVerts.length];
			A = mesh.facesOnVertex(v1);
			B = mesh.facesOnVertex(v2);
			for (var k = 0; k < A.length; k++) {
				for (var l = 0; l < B.length; l++) {
					if (A[k] === B[l])
						adjFace = A[k];
				}
			}
			mesh.splitFaceMakeEdge(adjFace, v1, v2);
			
			if (j > 0) { 
				var edge = mesh.edgeBetweenVertices(v1, verts[i]);
				mesh.joinFaceKillEdgeSimple(edge);
			}
		}
		

		origOtherVerts = origMesh.verticesOnVertex(origVerts[i]);
		for (var j = 0; j < newVerts.length; j++) {
			var newPos = new THREE.Vector3(0, 0, 0);
			var p1 = new THREE.Vector3(0, 0 ,0);
			var p2 = new THREE.Vector3(0, 0 ,0);
			p1.copy(origVerts[i].position);
			p2.copy(origOtherVerts[j].position);

			newPos.add(p1.multiplyScalar(1 - factor)).add(p2.multiplyScalar(factor));
			newVerts[j].position = newPos;
		}
		
		
		var newPos = new THREE.Vector3(0, 0, 0);
		var p1 = new THREE.Vector3(0, 0 ,0);
		var p2 = new THREE.Vector3(0, 0 ,0);
		p1.copy(origVerts[i].position);
		p2.copy(origOtherVerts[otherVerts.length-1].position);

		newPos.add(p1.multiplyScalar(1 - factor)).add(p2.multiplyScalar(factor));
		verts[i].position = newPos;
		
	}

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Apply the bevel operation to the mesh, scaling the degree of bevelling by factor
Filters.bevel = function ( mesh, factor ) {

    var verts = mesh.getModifiableVertices();
    this.truncate(mesh, factor);
    var origFaces = mesh.getModifiableFaces();
    var faces = [];
    var newFaces = [];
    var removedEdge = new Set();

    for(var i = 0; i < origFaces.length; i++) {
      if(mesh.verticesOnFace(origFaces[i]).length == 3) {
        faces.push(origFaces[i]);
      }
    }
    
    for(var i = 0; i < faces.length; i++) {
      var halfedges = mesh.edgesOnFace(faces[i]);
      for(var j = 0; j < halfedges.length; j++) {
        var v1 = halfedges[j].vertex;
        var v2 = halfedges[j].opposite.vertex;
        var f1 = halfedges[j].opposite.face
        newFaces.push([mesh.splitEdgeMakeVert(v1, v2, 1/2), f1]);
      }
    }

    
    for(var i = 0; i < newFaces.length; i++) {
      var f1 = newFaces[i][1]
      var f2 = newFaces[i][0]
      var f3 = f2.halfedge.opposite.next.next.next.vertex
      mesh.splitFaceMakeEdge(f1, f2, f3);
      var edge = mesh.edgeBetweenVertices(f2.halfedge.opposite.next.next.vertex,
                                          f2.halfedge.opposite.next.vertex)
      removedEdge.add(edge);
    }
    // Set iteration
    for(var edge of removedEdge) {
      mesh.joinFaceKillEdge(edge.face, edge.opposite.face, edge.vertex, edge.opposite.vertex);
      mesh.joinEdgeKillVert(edge.vertex.halfedge.vertex, edge.vertex, edge.vertex.halfedge.opposite.next.vertex)
    }

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Split the longest edges in the mesh into shorter edges.
// factor is a float in [0,1]. it tells the proportion
// of the total number of edges in the mesh that should be split.
Filters.splitLong = function ( mesh, factor  ) {
  var initFaces, intEdges = [], faces, edges, longEdge, longDist, vertices, count;

	initFaces = mesh.getModifiableFaces();
	for (var i = 0; i < initFaces.length; i++) {
		var fEdges = mesh.edgesOnFace(initFaces[i]);
		for (var j = 0; j < fEdges.length; j++) {
			if(intEdges.indexOf(fEdges[j]) == -1) 
				intEdges.push(fEdges[j]);
		}
	}
  var iterations = factor * intEdges.length;
	
	while (iterations--) {
		var faces = mesh.getModifiableFaces();
		var edges = [];
		for (var i = 0; i < faces.length; i++) {
			var fEdges = mesh.edgesOnFace(faces[i]);
			for (var j = 0; j < fEdges.length; j++) {
				if(edges.indexOf(fEdges[j]) == -1) 
					edges.push(fEdges[j]);
			}
		}

		
		longDist = Number.MIN_SAFE_INTEGER;
    var initFaces, intEdges = [], faces, edges, longEdge, longDist, vertices, count;
    
		for (var j = 0; j < edges.length; j++) {
			var d = mesh.dist(edges[j].vertex.position, edges[j].opposite.vertex.position);
			if (d > longDist) {
				longDist = d;
				longEdge = edges[j];
			}
		}
		
		vertices = mesh.verticesOnFace(longEdge.face);
		count = 0;
		var v = vertices[count], v1 = longEdge.vertex, v2 = longEdge.opposite.vertex;
		while (true) {
			if (v !== v1){
        if(v !== v2) {
				    break;
        }
      }
			v = vertices[++count];
		}
		
		var newVert = mesh.splitEdgeMakeVert(v1, v2, 0.5);
		mesh.splitFaceMakeEdge(longEdge.face, newVert, v);
	}

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Triangulate a mesh, and apply triangular subdivision to its faces.
// Repeat for the specified number of levels.
Filters.triSubdiv = function (mesh, levels) {
    Filters.triangulate( mesh );
    while(levels--) {
        var faces = mesh.getModifiableFaces();
        var N = faces.length
        var vertI = 0, oldI = 0, newI = 0;
        var v1, v2, v3, v4, v5, v6, he1, he2, he3;
        var vertsOld = [], newVerts = [], verts = [];

        for (var i = 0; i < N; i++) {
            var vertices = mesh.verticesOnFace( faces[i] )
            for ( var j = 0; j < 3; ++j, vertI++) { 
                vertsOld[vertI] = vertices[j]; 
            }
        }

        for ( var i = 0; i < N; i++, newI++) {
            verts = []
            for ( var j = 0; j < 3; j++, oldI++) { 
                verts.push(vertsOld[oldI]); 
            }

            v1 = verts[0]
            v2 = verts[1]
            v3 = verts[2]
            he1 = mesh.edgeBetweenVertices(v1, v2);
            he2 = mesh.edgeBetweenVertices(v1, v3);
            he3 = mesh.edgeBetweenVertices(v2, v3);

            if (he1) v4 = mesh.splitEdgeMakeVert(v1, v2, 1/2) 
            if (he2) v5 = mesh.splitEdgeMakeVert(v1, v3, 1/2) 
            if (he3) v6 = mesh.splitEdgeMakeVert(v2, v3, 1/2) 
            // Then
            if (!he1) v4 = mesh.vertBetweenVertices (v1,v2)  
            if (!he2) v5 = mesh.vertBetweenVertices (v1,v3)  
            if (!he3) v6 = mesh.vertBetweenVertices (v2,v3)  

            newVerts[newI] = v4;
            newI = newI + 1;
            newVerts[newI] = v5;
            newI = newI + 1;
            newVerts[newI] = v6;   
        }

        for (var i = 0, vi = 0; i < N; i ++) {
            verts = []
            for ( var j = 0; j < 3; j++, vi++) { 
                verts.push(newVerts[vi]); 
            }
            v4 = verts[0]
            v5 = verts[1]
            v6 = verts[2]
            f1 = mesh.splitFaceMakeEdge(faces[i], v4, v6, v5, true)
            f2 = mesh.splitFaceMakeEdge(f1, v4, v5, v6, true)
            mesh.splitFaceMakeEdge(f2, v5, v6, v4, true) 
        }
    }

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

//Triangulate the mesh and apply loop subdivision to the faces
//repeat for the specified number of levels.
Filters.loop = function (mesh, levels) {
    Filters.triangulate(mesh);
    betaOne = 3 / 8
    betaTwo = 1 /8

    while(levels--) {
        var faces = mesh.getModifiableFaces();
        var N = faces.length;    

        var oldVerts = [];
        var vertI = 0;

        for (var i = 0; i < N; i++) {
            var vertices = mesh.verticesOnFace(faces[i])
            for (var j = 0; j < 3; ++j, vertI++) { 
                oldVerts[vertI] = vertices[j]; 
            }
        }

        var iOld = 0, iNew = 0, he1, he2, he3;
        var v1, v2, v3, v4, v5, v6;
        var newVerts = [], verts = [];
        var oldX = [], oldY = [], oldZ = [], newX = [], newY = [], newZ = [];
        var selfX, selfY, selfZ, totalX, totalY, totalZ, otherX, otherY, otherZ;
        var neighbors, beta;

        for (var i = 0; i < N; i++, iNew++) {
            verts = []
            
            for (var j = 0; j < 3; j++, iOld++) { 
                verts[j] = oldVerts[iOld]; 
            }

            v1 = verts[0]
            v2 = verts[1]
            v3 = verts[2]

            he1 = mesh.edgeBetweenVertices(v1, v2);
            he2 = mesh.edgeBetweenVertices(v1, v3);
            he3 = mesh.edgeBetweenVertices(v2, v3);

            if (he1) v4 = mesh.splitEdgeMakeVert(v1, v2, 1/2) 
            if (he2) v5 = mesh.splitEdgeMakeVert(v1, v3, 1/2) 
            if (he3) v6 = mesh.splitEdgeMakeVert(v2, v3, 1/2) 

            // Then

            if (!he1) v4 = mesh.vertBetweenVertices (v1,v2)  
            if (!he2) v5 = mesh.vertBetweenVertices (v1,v3)  
            if (!he3) v6 = mesh.vertBetweenVertices (v2,v3)  

            newVerts[iNew] = v4;
            iNew = iNew + 1;
            newVerts[iNew] = v5;
            iNew = iNew + 1;
            newVerts[iNew] = v6;   
        }

        var opp1, opp2
        var twoFaces = mesh.facesOnVertices(v1, v2)
        var face1 = twoFaces[0]
        var face2 = twoFaces[1]
        var verts1 = mesh.verticesOnFace(face1)
        var verts2 = mesh.verticesOnFace(face2)
        var index1 = 0, index2 = 0;
        var ringVerts, locSum;



        for (var i = 0; i < 1; i ++) {
            totalX = 0;
            totalY = 0; 
            totalZ = 0;
            neighbors = mesh.verticesOnVertex(oldVerts[i])

            if (oldVerts.length > 3) { 
              beta = 3 / (8 * oldVerts.length) 
            } else { 
              beta = 3 / 16 
            }

            for (var j = 0; j < neighbors.length; j ++) {
                totalX = totalX + beta * neighbors[j].position.x
                totalY = totalY + beta * neighbors[j].position.y
                totalZ = totalZ + beta * neighbors[j].position.z
            }

            selfX = (1 - (neighbors.length * beta)) * oldVerts[i].position.x
            selfY = (1 - (neighbors.length * beta)) * oldVerts[i].position.y
            selfZ = (1 - (neighbors.length * beta)) * oldVerts[i].position.z

            oldX[i] = totalX + selfX
            oldY[i] = totalY + selfX
            oldZ[i] = totalZ + selfX
        }

        for (var i = 0; i < newVerts.length; i++) {
            totalX = 0; 
            totalY = 0; 
            totalZ = 0;
            neighbors = mesh.verticesOnVertex(newVerts[i])

            totalX = totalX + beta * neighbors[0].position.x + beta * neighbors[1].position.x;
            totalY = totalY + beta * neighbors[0].position.y + beta * neighbors[1].position.y;
            totalZ = totalZ + beta * neighbors[0].position.z + beta * neighbors[1].position.z;
            
            index1 = 0;
            index2 = 0;

            do {
              if (verts1[index1] !== neighbors[0]) {
                  if(verts1[index1] !== neighbors[1]) { 
                      opp1 = verts1[index1]; 
                  }
              }
              index1++;
              if (index1 > verts1.length - 1) {
                  break;
              }
            } while(true)


            do {
                if (verts2[index2] !== neighbors[0]){
                    if (verts2[index2] !== neighbors[1]) { 
                        opp2 = verts2[index2]; 
                    }
                }
                index2++;

                if (index2 > verts2.length - 1) {
                    break;
                }
            } while(true)

            otherX = betaTwo * opp1.position.x + betaTwo * opp2.position.x
            otherY = betaTwo * opp1.position.y + betaTwo * opp2.position.y
            otherZ = betaTwo * opp1.position.z + betaTwo * opp2.position.z
           
            newX[i] = totalX + otherX
            newY[i] = totalY + otherY
            newZ[i] = totalZ + otherZ
        }

        for (var i = 0, vi = 0; i < N; i ++) {
            var verts = []
            
            for (var j = 0; j < 3; j++, vi++) { 
                verts[j] = newVerts[vi]; 
            }

            v4 = verts[0]
            v5 = verts[1]
            v6 = verts[2]

            f1 = mesh.splitFaceMakeEdge(faces[i], v4, v6, v5, true)
            f2 = mesh.splitFaceMakeEdge(f1, v4, v5, v6, true)
            f3 = mesh.splitFaceMakeEdge(f2, v5, v6, v4, true) 

        }
        for (var i = 0; i < newVerts.length; i ++) {
            ringVerts = mesh.verticesOnVertex(newVerts[i]);
            locSum = new THREE.Vector3(0,0,0);
            for (var j = 0; j < ringVerts.length; j++) {
                locSum.add(ringVerts[j].position);
            }
            newVerts[i].position = locSum.divideScalar(ringVerts.length);
        }

        for (var i = 0; i < oldVerts.length; i ++) {
            ringVerts = mesh.verticesOnVertex(oldVerts[i]);
            locSum = new THREE.Vector3(0,0,0);
            for (var j = 0; j < ringVerts.length; j++) {
                locSum.add(ringVerts[j].position);
            }
            oldVerts[i].position = locSum.divideScalar(ringVerts.length);
        }

    }

    mesh.calculateFacesArea();
    mesh.updateNormals();
};


// Requires a quad mesh. Apply quad subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.quadSubdiv = function(mesh, levels) {
  for (var l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 55 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Quad subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply catmull clark subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.catmullClark = function(mesh, levels) {
  for (var l = 0; l < levels; l++) {
    const faces = mesh.faces;
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 102 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Catmull-Clark subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// ================= internal functions =======================

// internal function for selecting faces in the form of a loop
Filters.procFace = function(mesh, f) {
  const faceFlags = new Array(mesh.faces.length);
  for (var i = 0; i < mesh.faces.length; i++) {
    faceFlags[i] = 0;
  }
  var sum = f.area;
  const start_he = f.halfedge.opposite.next;
  var curr_he = start_he;
  do {
    if (faceFlags[curr_he.face.id] > 0) {
      break;
    }
    sum += curr_he.face.area;
    curr_he.face.selected = true;
    faceFlags[curr_he.face.id]++;
    const last_he = curr_he;
    curr_he = curr_he.opposite.next;
    if (curr_he.face == f) {
      curr_he = last_he.next.opposite.next;
    }
  } while (true);
};

Filters.parseSelected = function(sel) {
  if (sel === undefined || sel.replace === undefined) {
    return [];
  }
  if (typeof sel === "number") {
    return [sel];
  }
  // sel = sel.replace(/[\(\)]/g,'');
  sel = sel.split(",");
  const parsedSel = [];
  for (var i = 0; i < sel.length; i++) {
    const idx = parseInt(sel[i]);
    if (!isNaN(idx)) {
      parsedSel.push(idx);
    }
  }
  return parsedSel;
};

// internal filter for updating selection
Filters.selection = function(mesh, vertIdxs, faceIdxs) {
  mesh.setSelectedVertices(Filters.parseSelected(vertIdxs));
  mesh.setSelectedFaces(Filters.parseSelected(faceIdxs));
};

// internal filter for setting display settings
Filters.displaySettings = function(
  mesh,
  showLabels,
  showHalfedge,
  shading,
  showVN,
  showFN,
  showGrid,
  showVertDots,
  showAxes,
  showVC,
  meshColor
) {
  Main.displaySettings.showIdLabels = showLabels;
  Main.displaySettings.wireframe = showHalfedge;
  Main.displaySettings.shading = shading;
  Main.displaySettings.showVN = showVN;
  Main.displaySettings.showFN = showFN;
  Main.displaySettings.showGrid = showGrid;
  Main.displaySettings.showVertDots = showVertDots;

  Main.displaySettings.showAxes = showAxes;
  Main.displaySettings.showVC = showVC;
  // Main.displaySettings.meshColor = meshColor;

  // Main.refreshDisplaySettings();
};

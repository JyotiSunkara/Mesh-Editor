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
  for (let i = 0; i < N; ++i) {
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

    if (y !== 0){
      var axis = new THREE.Vector3(0, 1, 0);
      verts[i].position.applyAxisAngle(axis, y);
    }

    if (z !== 0){
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
  for (let i = 0; i < N; ++i) {
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
Filters.smooth = function(mesh, iter, delta, curvFlow = false, scaleDep = false, implicit = false) {
     for (var t = 0; t < iter; t++){
        var verts = mesh.getModifiableVertices();

    //Calculate weighted averages
    var averagePos = [];


    //Calculates weighted for each vertex
    for (var i = 0; i < verts.length; i++){
        var sigma = mesh.averageEdgeLength(verts[i]);
        var sumGaussian = 1;

        var curPos = new THREE.Vector3(0, 0, 0);
        curPos.copy(verts[i].position);
 
        var neighbors = Mesh.prototype.verticesOnVertex(verts[i]);

        for (var j = 0; j < neighbors.length; j++){
            var gaussian = Math.exp(-verts[i].position.distanceToSquared(neighbors[j].position) / (2 * sigma * sigma));
            sumGaussian += gaussian;

            var tmpNeighbor = new THREE.Vector3(0, 0, 0);
            tmpNeighbor.copy(neighbors[j].position);

            //console.log("Tmp " + tmpNeighbor.x + ", " + tmpNeighbor.y + ", " + tmpNeighbor.z);

            curPos.add(tmpNeighbor.multiplyScalar(gaussian));
            //console.log("New " + curPos.x + ", " + curPos.y + ", " + curPos.z);
        }

        curPos.multiplyScalar(1/sumGaussian);
        //console.log("PPos " + curPos.x + ", " + curPos.y + ", " + curPos.z);
        //console.log("Pos " + curPos.x + ", " + curPos.y + ", " + curPos.z);
        averagePos.push(curPos);
    }

    
    for (var i = 0; i < verts.length; i++){
        //console.log("Before " + verts[i].position.x + ", " + verts[i].position.y + ", " + verts[i].position.z);
        verts[i].position = averagePos[i];
        //console.log("After " + verts[i].position.x + ", " + verts[i].position.y + ", " + verts[i].position.z);
    }

    mesh.updateNormals();
    mesh.calculateFacesArea();
    }
};

// Sharpen the mesh by moving selected vertices away from the average position
// of their neighbors (i.e. Laplacian smoothing in the negative direction)
Filters.sharpen = function(mesh, iter, delta) {
  var verts = mesh.getModifiableVertices();

	var smoothMesh = new Mesh();
	smoothMesh.copy(mesh);
	Filters.smooth(smoothMesh, iter, delta);
	var smoothVerts = smoothMesh.getModifiableVertices();
	
	var N = verts.length;
	for ( var i = 0 ; i < N ; ++i ) {
    var pos = new THREE.Vector3(verts[i].position.x, verts[i].position.y, verts[i].position.z);
		var diff = pos.sub(smoothVerts[i].position);
		verts[i].position.add(diff);
	}
  
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Move every selected vertex along its normal direction
// Scale the amount by the provided factor and average edge length at that vertex
Filters.inflate = function (mesh, factor) {

  var verts = mesh.getModifiableVertices();

	var origMesh = new Mesh();
	origMesh.copy(mesh);
	var origVerts = origMesh.getModifiableVertices();

	var N = verts.length;
  for ( var i = 0 ; i < N ; ++i ) {
		var avgLen = origMesh.averageEdgeLength(origVerts[i]);
		var temp = new THREE.Vector3(origVerts[i].normal.x, origVerts[i].normal.y, origVerts[i].normal.z)
		temp.multiplyScalar(factor * avgLen);
    verts[i].position.add(temp);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// rotate selected vertices around the Y axis by an amount
// proportional to its Y value times the scale factor.
Filters.twist = function(mesh, factor) {
  var verts = mesh.getModifiableVertices();

  var verts = mesh.getModifiableVertices();

  var n_vertices = verts.length;
	for ( var i = 0 ; i < n_vertices ; ++i ) { 
		var vx = verts[i].position.getComponent(0);
		var vy = verts[i].position.getComponent(1);
		var vz = verts[i].position.getComponent(2);
		var rot = vy*factor;
		var t_y = new THREE.Vector3(Math.cos(rot)*vx+Math.sin(rot)*vz, vy, -Math.sin(rot)*vx+Math.cos(rot)*vz);
        verts[i].position = t_y;
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// warp a mesh using a nonlinear mapping of your choice
Filters.wacky = function(mesh, factor) {
	var verts = mesh.getModifiableVertices();
	var N = verts.length;
	
	var maxH, minH;
	var minH = maxH = verts[0].position.y;

	for ( var i = 0 ; i < N ; ++i ) {
		if (verts[i].position.y > maxH) {
			maxH = verts[i].position.y;
		}
		if (verts[i].position.y < minH) {
			minH = verts[i].position.y;
		}
	}

	var H = maxH - minH;

  for ( var i = 0 ; i < N ; ++i ) {
		  var temp = new THREE.Vector3(H/2*Math.sin(((verts[i].position.y/H) * 22)/(7 * factor)), 0, 0)
      verts[i].position.add(temp);
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
    mesh.splitEdgeMakeVert(verts[0], verts[1], 0.5);
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
    let v0 = verts[0],
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

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 32 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Extrude is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Truncate the selected vertices of the mesh by "snipping off" corners
// and replacing them with faces. factor specifies the size of the truncation.
// See the spec for more detail.
Filters.truncate = function ( mesh, factor ) {

  origMesh = new Mesh();
	origMesh.copy(mesh);
	var orig_verts = origMesh.getModifiableVertices();
	var verts = mesh.getModifiableVertices();
	var n_verts = verts.length;
	for ( var i = 0 ; i < n_verts ; ++i ) {

		//make new verts
		var otherVerts = mesh.verticesOnVertex(verts[i]);
		var newVerts = []
		for (var v_i = 0; v_i < otherVerts.length-1; ++v_i) {
			newVerts[v_i] = mesh.splitEdgeMakeVert(verts[i], otherVerts[v_i], factor);
		}
			
		//make new face
		for (var v_i = 0; v_i < newVerts.length-1; v_i++) {
			faces_a = mesh.facesOnVertex(newVerts[v_i]);
			faces_b = mesh.facesOnVertex(newVerts[(v_i+1)%newVerts.length]);
			var common_face = undefined;
			for (var f_i = 0; f_i < faces_a.length; f_i++) {
				for (var f_j = 0; f_j < faces_b.length; f_j++) {
					if (faces_a[f_i] === faces_b[f_j])
						common_face = faces_a[f_i];
				}
			}
			mesh.splitFaceMakeEdge(common_face, newVerts[v_i], newVerts[(v_i+1)%newVerts.length]);
			
			if (v_i > 0) { //only make one new face
				var edge = mesh.edgeBetweenVertices(newVerts[v_i], verts[i]);
				mesh.joinFaceKillEdgeSimple(edge);
			}
		}
		

		//move to correct spot base on origMesh
		origOtherVerts = origMesh.verticesOnVertex(orig_verts[i]);
		for (var v_i = 0; v_i < newVerts.length; v_i++) {
			var new_pos = new THREE.Vector3( 0, 0, 0 );
			var p1 = new THREE.Vector3( 0, 0 ,0 );
			p1.copy( orig_verts[i].position );
			var p2 = new THREE.Vector3( 0, 0 ,0 );
			p2.copy( origOtherVerts[v_i].position );
			new_pos.add( p1.multiplyScalar( 1 - factor ) );
			new_pos.add( p2.multiplyScalar( factor ) );
			newVerts[v_i].position = new_pos;
		}
		
		
		//move the last vertex
		var new_pos = new THREE.Vector3( 0, 0, 0 );
		var p1 = new THREE.Vector3( 0, 0 ,0 );
		p1.copy( orig_verts[i].position );
		var p2 = new THREE.Vector3( 0, 0 ,0 );
		p2.copy( origOtherVerts[otherVerts.length-1].position );
		new_pos.add( p1.multiplyScalar( 1 - factor ) );
		new_pos.add( p2.multiplyScalar( factor ) );
		verts[i].position = new_pos;
		
	}
    // ----------- Our reference solution uses 54 lines of code.
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('Truncate is not implemented yet');

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Apply the bevel operation to the mesh, scaling the degree of bevelling by factor
Filters.bevel = function ( mesh, factor ) {

    var verts = mesh.getModifiableVertices();

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 104 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('Bevel is not implemented yet');

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Split the longest edges in the mesh into shorter edges.
// factor is a float in [0,1]. it tells the proportion
// of the total number of edges in the mesh that should be split.
Filters.splitLong = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 35 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Split Long Edges is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate a mesh, and apply triangular subdivision to its faces.
// Repeat for the specified number of levels.
Filters.triSubdiv = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 43 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate the mesh and apply loop subdivision to the faces
// repeat for the specified number of levels.
Filters.loop = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 123 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Requires a quad mesh. Apply quad subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.quadSubdiv = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
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
  for (let l = 0; l < levels; l++) {
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
  for (let i = 0; i < mesh.faces.length; i++) {
    faceFlags[i] = 0;
  }
  let sum = f.area;
  const start_he = f.halfedge.opposite.next;
  let curr_he = start_he;
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
  for (let i = 0; i < sel.length; i++) {
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

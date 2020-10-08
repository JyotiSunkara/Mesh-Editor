// Empty mesh constructor
function Mesh() {
  this.clear();
}

// Vertex constructor, only stores position
// Should not be used directly, vertices should be added via addVertex method or Euler operators
function Vertex(x, y, z) {
  this.id = -1;
  this.position = new THREE.Vector3(x, y, z);
  this.normal = new THREE.Vector3();
  this.color = new THREE.Vector3(Math.random(), Math.random(), Math.random());
  this.curvature = 0.0;
  this.selected = false;

  this.halfedge = undefined;

  this.removed = undefined;
}

// Face Constructor
// Should not be used directly, faces should be added via addFace method or Euler operators
function Face() {
  this.id = -1;
  this.normal = new THREE.Vector3();
  this.area = 0.0;
  this.centroid = undefined;
  this.selected = false;

  this.halfedge = undefined;

  this.removed = undefined;
}

// HalfEdge Constructor
// Should not be used directly, faces should be added via addFace method or Euler operators
function HalfEdge() {
  this.id = -1;
  this.selected = false;
  this.midpoint = undefined;

  this.vertex = undefined;
  this.next = undefined;
  this.opposite = undefined;
  this.face = undefined;

  this.removed = undefined;
}

////////////////////////////////////////////////////////////////////////////////
// Data Structure Modification (lowest, lowest level; students need not modify)
////////////////////////////////////////////////////////////////////////////////

Mesh.prototype.clear = function() {
  this.vertices = [];
  this.halfedges = [];
  this.faces = [];
};

Mesh.prototype.copy = function(mesh) {
  this.vertices = [];
  this.faces = [];
  this.halfedges = [];
  let v_id, f_id, he_id;

  assert(mesh, "Mesh.copy must take a mesh to copy from");

  // create list of vertices
  for (v_id in mesh.vertices) {
    const v_org = mesh.vertices[v_id];
    const v_cpy = this.addVertex(v_org.position);
    v_cpy.normal.copy(v_org.normal);
    v_cpy.color.copy(v_org.color);
    v_cpy.curvature = v_org.curvature;
    v_cpy.selected = v_org.selected;
  }

  // create list of faces
  for (f_id in mesh.faces) {
    var f_org = mesh.faces[f_id];
    var f_cpy = this.addFace();
    if (f_org.normal !== undefined) {
      f_cpy.normal = new THREE.Vector3(f_org.normal.x, f_org.normal.y, f_org.normal.z);
    }
    f_cpy.selected = f_org.selected;
  }

  // create all half edges
  for (he_id in mesh.halfedges) {
    var he_org = mesh.halfedges[he_id];
    var f_org = he_org.face;
    var f_cpy = this.faces[f_org.id];

    var he_cpy = this.addHalfEdge(
      this.vertices[he_org.opposite.vertex.id],
      this.vertices[he_org.vertex.id],
      f_cpy
    );
    he_cpy.selected = he_org.selected;
  }

  // relink halfedges
  for (he_id in mesh.halfedges) {
    var he_org = mesh.halfedges[he_id];
    var he_cpy = this.halfedges[he_id];
    he_cpy.next = this.halfedges[he_org.next.id];
    he_cpy.opposite = this.halfedges[he_org.opposite.id];
  }

  // add face reference to halfedges
  for (f_id in mesh.faces) {
    var f_org = mesh.faces[f_id];
    var f_cpy = this.faces[f_id];
    f_cpy.halfedge = this.halfedges[f_org.halfedge.id];
  }
};

Mesh.prototype.addVertex = function(position) {
  const v = new Vertex(position.x, position.y, position.z);
  v.id = this.vertices.length;

  this.vertices.push(v);

  return v;
};

Mesh.prototype.addFace = function() {
  const f = new Face();
  f.id = this.faces.length;
  this.faces.push(f);
  f.halfedge = undefined;

  return f;
};

Mesh.prototype.addHalfEdge = function(origin, destination, face) {
  const he = new HalfEdge();
  he.id = this.halfedges.length;
  this.halfedges.push(he);

  he.vertex = destination;
  he.face = face;
  origin.halfedge = he;

  he.next = undefined;
  he.opposite = undefined;

  return he;
};

// TODO make removal less slow
Mesh.prototype.removeVertex = function(v) {
  const idx = this.vertices.indexOf(v);
  if (idx > -1) {
    this.vertices.splice(idx, 1);
  }
  v.removed = true;
  for (let i = 0; i < this.vertices.length; i++) {
    this.vertices[i].id = i;
  }
};

Mesh.prototype.removeFace = function(f) {
  const idx = this.faces.indexOf(f);
  if (idx > -1) {
    this.faces.splice(idx, 1);
  }
  f.removed = true;
  for (let i = 0; i < this.faces.length; i++) {
    this.faces[i].id = i;
  }
};

Mesh.prototype.removeHalfEdge = function(he) {
  const idx = this.halfedges.indexOf(he);
  if (idx > -1) {
    this.halfedges.splice(idx, 1);
  }
  he.removed = true;
  for (let i = 0; i < this.halfedges.length; i++) {
    this.halfedges[i].id = i;
  }
};

////////////////////////////////////////////////////////////////////////////////
// Euler Operators
////////////////////////////////////////////////////////////////////////////////
// these operations always return a valid mesh (when called on a valid mesh with valid parameters)
// for further descriptions of functions, see the section on Euler Operators at: http://wiki.blender.org/index.php/Dev:2.6/Source/Modeling/BMesh/Design

// splits the edge between v1 and v2 into two, creating a new vert
// v1 -- v2   =>   v1 -- newVert -- v2
// the previously existing halfedges will point to newVert and new ones will be added pointing out from newVert
// returns the new vert, newVert
Mesh.prototype.splitEdgeMakeVert = function(v1, v2, factor) {
  if (factor === undefined) {
    factor = 0.5;
  }
  // get all relevant info
  const he1 = this.edgeBetweenVertices(v1, v2);
  if (!he1) {
    return false;
  }
  const he2 = he1.opposite;
  const f1 = he1.face;
  const f2 = he2.face;
  const he1_next = he1.next;
  const he2_next = he2.next;

  // compute new vertex position
  const new_pos = new THREE.Vector3(0, 0, 0);
  const p1 = new THREE.Vector3(0, 0, 0);
  p1.copy(v1.position);
  const p2 = new THREE.Vector3(0, 0, 0);
  p2.copy(v2.position);

  new_pos.add(p1.multiplyScalar(1 - factor));
  new_pos.add(p2.multiplyScalar(factor));

  // create new vertex and halfedges
  const newVert = this.addVertex(new_pos);
  const he3 = this.addHalfEdge(newVert, v2, f1);
  const he4 = this.addHalfEdge(newVert, v1, f2);

  he1.vertex = newVert;
  he2.vertex = newVert;

  // relink everything
  he3.next = he1_next;
  he1.next = he3;
  he1.opposite = he4;
  he4.opposite = he1;

  he4.next = he2_next;
  he2.next = he4;
  he2.opposite = he3;
  he3.opposite = he2;

  f1.halfedge = he3;
  f2.halfedge = he4;

  return newVert;
};

// takes 3 verts, v2 adjacent only to v1 and v3
// v1 -- v2 -- v3   =>   v1 -- v3
// removes v2 and the edges pointing out from v2 (and relinks the others)
// returns true if successful, false otherwise
Mesh.prototype.joinEdgeKillVert = function(v1, v2, v3) {
  // get the halfedges to update
  let he12 = this.edgeBetweenVertices(v1, v2);
  let he32 = this.edgeBetweenVertices(v3, v2);
  // assert(he12, "v1 and v2 must share an edge");
  // assert(he32, "v2 and v3 must share an edge");
  if (he12 === undefined) {
    return false;
  }
  if (he32 === undefined) {
    return false;
  }

  // make sure they're pointing at v2
  // we'll be removing the halfedges pointing away from v2
  he12 = he12.vertex == v2 ? he12 : he12.opposite;
  he32 = he32.vertex == v2 ? he32 : he32.opposite;

  // update vertices to skip over v2
  he12.vertex = he32.opposite.vertex;
  he32.vertex = he12.opposite.vertex;

  // update next halfedges to skip over v2
  he12.next = he32.opposite.next;
  he32.next = he12.opposite.next;

  // make sure the faces don't point to the removed halfedges
  // don't necessarily need to update, but easier to do it than to check
  he12.face.halfedge = he12;
  he32.face.halfedge = he32;

  this.removeHalfEdge(he12.opposite);
  this.removeHalfEdge(he32.opposite);

  // make them opposite each other
  he12.opposite = he32;
  he32.opposite = he12;

  this.removeVertex(v2);

  return true;
};

// takes a face and two vertices on it
// creates a new edge between v1 and v2, splitting the face in two
// if optional arg vertOnF is supplied, ensures that vertOnF remains on face f after the split (otherwise, it is undefined which way the face is split)
// if switchFaces is true, faces will be assigned so that vertOnF lies on the new face
// returns the new face (which points to one of the new halfedges)
Mesh.prototype.splitFaceMakeEdge = function(f, v1, v2, vertOnF, switchFaces) {
  // Go around f and find halfedges adjacent to v1 and v2
  let heToV2, heToV1;
  let heFromV1, heFromV2;
  // also build a string to track relative locations of v1, v2, and vertOnF
  let vertOrder = "";
  var he = f.halfedge;
  var first = he;
  do {
    if (he.vertex === v1) {
      heToV1 = he;
      vertOrder += "1";
    }
    if (he.vertex === v2) {
      heToV2 = he;
      vertOrder += "2";
    }
    if (he.vertex === vertOnF) {
      vertOrder += "3";
    }
    if (he.opposite.vertex === v1) {
      heFromV1 = he;
    }
    if (he.opposite.vertex === v2) {
      heFromV2 = he;
    }
    he = he.next;
  } while (he !== first);

  // create new halfedges
  const he12 = this.addHalfEdge(v1, v2);
  const he21 = this.addHalfEdge(v2, v1);
  he21.opposite = he12;
  he12.opposite = he21;

  // link new halfedges with preexisting ones
  heToV1.next = he12;
  he12.next = heFromV2;
  heToV2.next = he21;
  he21.next = heFromV1;

  // create new face
  const newF = this.addFace();

  // assign new halfedges to faces (and vice versa)
  let heOnF;
  // for heOnF to be he21, vertOnF must occur after v1 and before v2
  if (vertOrder === "132" || vertOrder === "321" || vertOrder === "213") {
    heOnF = he21;
  } else {
    // the alternative (also the default if vertOnF isn't given )
    heOnF = he12;
  }
  if (switchFaces) {
    heOnF = heOnF == he12 ? he21 : he12;
  }
  const heOnNewF = heOnF == he12 ? he21 : he12;
  f.halfedge = heOnF;
  heOnF.face = f;
  newF.halfedge = heOnNewF;
  heOnNewF.face = newF;

  // go around each face and make sure halfedges point to it
  for (let i = 0; i < 2; i++) {
    const curFace = i == 0 ? f : newF;

    var he = curFace.halfedge;
    var first = he;
    do {
      he.face = curFace;
      he = he.next;
    } while (he !== first);
  }

  return newF;
};

// takes two faces and two vertices which denote an edge between them
// f1 and f2 must be separated only by the edge (v1, v2)
// removes the edge (v1, v2) and the face f2, then expands f1 to fill the space of both faces
// returns true if successful, false otherwise
Mesh.prototype.joinFaceKillEdge = function(f1, f2, v1, v2) {
  let he = this.edgeBetweenVertices(v1, v2);
  if (!he) {
    return false;
  }

  // make sure he points from v1 to v2
  if (he.vertex !== v2) {
    he = he.opposite;
  }

  // make other halfedges point past the edge to be removed
  for (var i = 0; i < 2; i++) {
    const startingHe = i === 0 ? he : he.opposite;

    let heToAdjust = startingHe.next;
    while (true) {
      if (heToAdjust.next === startingHe) {
        heToAdjust.next = startingHe.opposite.next;
        break;
      }
      heToAdjust = heToAdjust.next;
    }
  }

  // make sure faces don't point at the edge to be removed
  f1.halfedge = he.next;
  f2.halfedge = he.next;

  // make sure verts don't point at the edge to be removed
  v1.halfedge = he.opposite.next;
  v2.halfedge = he.next;

  // make remaining edges point to f1 rather than f2
  const edges = this.edgesOnFace(f2);
  for (var i = 0; i < edges.length; i++) {
    const e = edges[i];
    e.face = f1;
  }

  this.removeFace(f2);
  this.removeHalfEdge(he);
  this.removeHalfEdge(he.opposite);

  for (var i = 0; i < this.vertices.length; ++i) {
    var v = this.vertices[i];
    if (v.removed) {
      continue;
    }
    assert(!v.halfedge.removed);
  }
  for (var i = 0; i < this.faces.length; ++i) {
    var v = this.faces[i];
    if (v.removed) {
      continue;
    }
    assert(!v.halfedge.removed);
  }
  for (var i = 0; i < this.halfedges.length; ++i) {
    var v = this.halfedges[i];
    if (v.removed) {
      continue;
    }
    assert(!v.next.removed);
    assert(!v.face.removed);
    assert(!v.vertex.removed);
    assert(!v.opposite.removed);
  }

  return true;
};

// a simpler interface for joinFaceKillEdge
// doesn't give control over which face is removed
Mesh.prototype.joinFaceKillEdgeSimple = function(he) {
  return this.joinFaceKillEdge(he.face, he.opposite.face, he.opposite.vertex, he.vertex);
};

////////////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////////////

// return a copy of vec, leaving the original unchanged
function CopyVec(vec) {
  return new THREE.Vector3(0, 0, 0).copy(vec);
}

// takes two verts sharing an edge, and a face touching that edge
// returns the halfedge between the two verts that belongs to the specificed face
Mesh.prototype.edgeBetweenVerticesFacingFace = function(v1, v2, f) {
  const e = this.edgeBetweenVertices(v1, v2);
  if (e.face == f) {
    return e;
  } else if (e.opposite.face == f) {
    return e.opposite;
  } else {
    console.log("edge does not lie on face!", e, f);
  }
};

// takes two halfedges directed away from a vertex v
// return the angle between them, in radians
Mesh.prototype.angleBetweenEdges = function(v, he1, he2) {
  const p0 = v.position;
  const p1 = he1.vertex.position;
  const p2 = he2.vertex.position;
  const v1 = new THREE.Vector3();
  const v2 = new THREE.Vector3();
  v1.subVectors(p1, p0);
  v2.subVectors(p2, p0);

  // Return angle between vectors
  const d1 = v1.length();
  if (Math.abs(d1) < 0.000001) {
    return 0.0;
  }
  const d2 = v2.length();
  if (Math.abs(d2) < 0.000001) {
    return 0.0;
  }
  const cosine = v1.dot(v2) / (d1 * d2);
  if (cosine >= 1.0) {
    return 0.0;
  } else if (cosine <= -1.0) {
    return Math.PI;
  } else {
    return Math.acos(cosine);
  }
};

Mesh.prototype.updateEdgeMidpoints = function() {
  for (let i = 0; i < this.halfedges.length; i++) {
    const he = this.halfedges[i];
    he.midpoint = CopyVec(he.vertex.position)
      .add(he.opposite.vertex.position)
      .multiplyScalar(0.5);
  }
};

Mesh.prototype.updateFaceCentroids = function() {
  for (let i = 0; i < this.faces.length; i++) {
    const f = this.faces[i];
    // get the centroid (inlined for performance)
    const verts = this.verticesOnFace(f);
    const centroid = new THREE.Vector3(0, 0, 0);
    for (let j = 0; j < verts.length; ++j) {
      centroid.add(verts[j].position);
    }
    centroid.divideScalar(verts.length);
    f.centroid = centroid;
  }
};

Mesh.prototype.calculateFaceCentroid = function(f) {
  const verts = this.verticesOnFace(f);
  const centroid = new THREE.Vector3(0, 0, 0);
  for (let i = 0; i < verts.length; ++i) {
    centroid.add(verts[i].position);
  }
  centroid.divideScalar(verts.length);
  return centroid;
};

Mesh.prototype.calculateFaceNormal = function(f) {
  // Get vertices of queried face
  const vertices = this.verticesOnFace(f);

  // Since every face has at least three vertices, we can getpositions of first three.
  // Here we assume that face is planar, even if number of vertices is greater than 3.

  // We find two edges that are not colinear, and their cross product is the normal.
  const nverts = vertices.length;
  const normal = new THREE.Vector3();

  for (let start = 0; start < nverts; start++) {
    const i0 = (start + 0) % nverts;
    const i1 = (start + 1) % nverts;
    const i2 = (start + 2) % nverts;
    const p0 = vertices[i0].position;
    const p1 = vertices[i1].position;
    const p2 = vertices[i2].position;

    const vec1 = new THREE.Vector3();
    vec1.subVectors(p0, p1);

    const vec2 = new THREE.Vector3();
    vec2.subVectors(p0, p2);

    // compute average edge length
    const vec3 = new THREE.Vector3();
    vec3.subVectors(p1, p2);
    const avelen = (vec1.length() + vec2.length() + vec3.length()) / 3.0;

    normal.crossVectors(vec1, vec2);

    // area of square would be like avelen^2 so compare cross product to that
    const square = avelen * avelen;
    const area = normal.length();
    const ratio = area / square;
    if (!isNaN(ratio) && ratio > 1e-6) {
      normal.normalize();
      return normal;
    }
  }
  console.log("oops - tried to find normal of degenerate face");
  return normal;
};

Mesh.prototype.updateFaceNormals = function() {
  for (let i = 0; i < this.faces.length; ++i) {
    this.faces[i].normal = this.calculateFaceNormal(this.faces[i]);
  }
};

Mesh.prototype.updateNormals = function() {
  this.updateFaceNormals();
  this.updateVertexNormals();
};

// number of faces that are selected
Mesh.prototype.numSelectedFaces = function() {
  let count = 0;
  for (let i = 0; i < this.faces.length; ++i) {
    if (this.faces[i].selected) {
      count++;
    }
  }
  return count;
};

// list of faces that are selected
Mesh.prototype.getModifiableFaces = function() {
  const faces = [];
  for (let i = 0; i < this.faces.length; ++i) {
    if (this.faces[i].selected) {
      faces.push(this.faces[i]);
    }
  }

  if (faces.length === 0) {
    return this.faces;
  }

  return faces;
};

// list of selected vertices and vertices on faces that are selected
Mesh.prototype.getModifiableVertices = function() {
  const verts_movable = [];
  for (var i = 0; i < this.faces.length; ++i) {
    var verts = this.verticesOnFace(this.faces[i]);
    for (let j = 0; j < verts.length; ++j) {
      if (this.faces[i].selected || verts[j].selected) {
        verts_movable[verts[j].id] = true;
      }
    }
  }
  var verts = [];
  if (verts_movable.length === 0) {
    verts = this.vertices;
  } else {
    for (var i = 0; i < verts_movable.length; ++i) {
      if (verts_movable[i] === true) {
        verts.push(this.vertices[i]);
      }
    }
  }
  return verts;
};

// list of selected vertices
Mesh.prototype.getSelectedVertices = function() {
  const selectedVerts = [];
  for (let i = 0; i < this.vertices.length; ++i) {
    if (this.vertices[i].selected) {
      selectedVerts.push(this.vertices[i]);
    }
  }
  return selectedVerts;
};

// set list of faces that are selected
Mesh.prototype.setSelectedFaces = function(sel) {
  for (var i = 0; i < this.faces.length; ++i) {
    this.faces[i].selected = false;
  }
  if (sel === undefined) {
    return;
  }
  for (var i = 0; i < sel.length; ++i) {
    const id = sel[i];
    if (id >= 0) {
      // just skip negative ids, to allow easies toggling
      this.faces[id].selected = true;
    }
  }
};

// set list of verts that are selected
Mesh.prototype.setSelectedVertices = function(sel) {
  for (var i = 0; i < this.vertices.length; ++i) {
    this.vertices[i].selected = false;
  }
  if (sel === undefined) {
    return;
  }
  for (var i = 0; i < sel.length; ++i) {
    const id = sel[i];
    if (id >= 0) {
      // just skip negative ids, to allow easies toggling
      this.vertices[id].selected = true;
    }
  }
};

// returns list of selected face ids
Mesh.prototype.getSelectedFaceIds = function() {
  const sel = [];
  for (let i = 0; i < this.faces.length; i++) {
    if (this.faces[i].selected) {
      sel.push(i);
    }
  }
  return sel;
};

// returns list of selected vert ids
Mesh.prototype.getSelectedVertexIds = function() {
  const sel = [];
  for (let i = 0; i < this.vertices.length; i++) {
    if (this.vertices[i].selected) {
      sel.push(i);
    }
  }
  return sel;
};

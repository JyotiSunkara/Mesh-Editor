////////////////////////////////////////////////////////////////////////////////
// Conversions / Constructors
// These are functions used to create meshes - students should have no need to
// modify this code
////////////////////////////////////////////////////////////////////////////////

Mesh.prototype.buildFromVerticesAndFaces = function(vertices, faces) {
  const start = new Date().getTime();

  var n_vertices = vertices.length;
  const n_faces = faces.length;

  // lets have an edge map, mapping the i,j vertex to a half edge
  function emIJ2Key(i, j) {
    return i + "_" + j;
  }

  function emKey2IJ(key) {
    const parts = key.split("_");
    const ret = [parts[0], parts[1]];
    return ret;
  }

  const edgeMap = {};

  for (var i = 0; i < n_vertices; ++i) {
    this.addVertex(vertices[i]);
  }

  for (var i = 0; i < n_faces; ++i) {
    const cur_face_ind = faces[i];
    const cur_vertices = [];

    for (var j = 0; j < cur_face_ind.length; j++) {
      cur_vertices.push(this.vertices[cur_face_ind[j]]);
    }

    const f = this.addFace();

    // add halfedges between consecutive vertices
    var n_vertices = cur_vertices.length;
    for (var j = 0; j < n_vertices; j++) {
      var next_j = (j + 1) % n_vertices;
      const he = this.addHalfEdge(cur_vertices[j], cur_vertices[next_j], f);
      edgeMap[emIJ2Key(cur_vertices[j].id, cur_vertices[next_j].id)] = he;
      f.halfedge = he;
    }

    // relink edges around face
    for (var j = 0; j < n_vertices; j++) {
      var next_j = (j + 1) % n_vertices;
      cur_vertices[j].halfedge.next = cur_vertices[next_j].halfedge;
    }
  }

  for (const key in edgeMap) {
    const he1 = edgeMap[key];
    const ind = emKey2IJ(key);
    const key2 = emIJ2Key(ind[1], ind[0]);
    const he2 = edgeMap[key2];
    he1.opposite = he2;
  }

  this.calculateFacesArea();
  this.updateNormals();

  const end = new Date().getTime();
  const elapsed = end - start;

  // console.log( "Conversion took " + elapsed + " ms. Mesh contains " + this.vertices.length + " vertices" );
};

Mesh.prototype.fromOBJ = function(filename, meshLoadCallback) {
  this.filename = filename;

  filename = "obj/" + filename; // all obj files are in the obj folder

  const start = new Date().getTime();

  const manager = new THREE.LoadingManager();

  // load using the three js loading manager plus pass reference to current mesh
  const loader = new OBJLoader(manager);
  const mesh = this;

  loader.load(filename, function(vertices, faces) {
    mesh.buildFromVerticesAndFaces(vertices, faces);
    meshLoadCallback();
  });
};

// this code is ugly since translateions and rotations need to be treated separately
Mesh.prototype.applyFilters = function(values) {
  //console.log(values);
  if (values == undefined) {
    return;
  }

  // first parse translations & rotations
  const translation = new THREE.Vector3(values.translateX, values.translateY, values.translateZ);
  const rotation = new THREE.Vector3(values.rotateX, values.rotateY, values.rotateZ);
  if (translation.x !== 0 || translation.y !== 0 || translation.z !== 0) {
    Filters.translate(this, translation);
  }

  if (rotation.x !== 0 || rotation.y !== 0 || rotation.z !== 0) {
    Filters.rotate(this, rotation);
  }

  //then all other values
  for (const prop in values) {
    if (prop in Gui.defaults && values[prop] !== Gui.defaults[prop]) {
      const val = values[prop];
      if (
        prop === "translateX" ||
        prop === "translateY" ||
        prop === "translateZ" ||
        prop === "rotateX" ||
        prop === "rotateY" ||
        prop === "rotateZ"
      ) {
        continue;
      }
      Filters[prop](this, val);
    }
  }
};

// generates the THREE.js geometry that's actually rendered (base mesh, selections, halfedge vis, etc)
Mesh.prototype.toBufferGeometry = function() {
  // useful variables and default settings
  let v_pos, v_nor, v_col;
  var i = 0,
    j = 0,
    k = 0,
    sel_k = 0,
    l = 0,
    idx = 0;

  const arrowSize = 0.005;

  this.updateFaceCentroids();
  this.updateEdgeMidpoints();

  // precalculate stuff for halfedge vis
  for (var i = 0; i < this.halfedges.length; i++) {
    var he = this.halfedges[i];
    // he.avgNormal = CopyVec(he.face.normal).add(he.opposite.face.normal).normalize();
    const edgeDir = CopyVec(he.vertex.position)
      .sub(he.opposite.vertex.position)
      .normalize();
    he.shortTangent = CopyVec(he.face.normal)
      .cross(edgeDir)
      .normalize()
      .multiplyScalar(0.05);
    he.shortDir = edgeDir.multiplyScalar(arrowSize);
  }

  const hideSelected = Gui.hidden;

  const n_faces = this.faces.length;

  let n_triangles = 0;
  let n_selected_triangles = 0;
  for (i = 0; i < n_faces; ++i) {
    const vertices = this.verticesOnFace(this.faces[i]);
    // need to figure number of vertices for a face
    if (this.faces[i].selected && !hideSelected) {
      n_selected_triangles += vertices.length - 2;
    } else {
      n_triangles += vertices.length - 2;
    }
  }

  // geometries
  const meshGeometry = new THREE.BufferGeometry();
  const selectedFacesGeo = new THREE.BufferGeometry();
  const faceNormalsGeo = new THREE.Geometry();
  const vertexNormalsGeo = new THREE.Geometry();
  const wireframeGeo = new THREE.Geometry();
  const verticesGeo = new THREE.Geometry();
  const selectedVerticesGeo = new THREE.Geometry();
  const structureVisGeo = new THREE.Geometry();
  const edgeStructureVisGeo = new THREE.Geometry();

  const unselectedVertexColors = [];

  // buffers for bufferGeometries
  // each face - 3 vertices - 3 attributes
  const vertex_positions = new Float32Array(n_triangles * 3 * 3);
  const sel_vertex_positions = new Float32Array(n_selected_triangles * 3 * 3);
  const vertex_normals = new Float32Array(n_triangles * 3 * 3);
  const sel_vertex_normals = new Float32Array(n_selected_triangles * 3 * 3);
  const vertex_colors = new Float32Array(n_triangles * 3 * 3);
  const sel_vertex_colors = new Float32Array(n_selected_triangles * 3 * 3);

  // text billboards
  const textSprites = [];

  for (i = 0; i < n_faces; i++) {
    const f = this.faces[i];

    // Face normals and halfedge vis
    centroid = f.centroid;

    // arrow from face to edge
    if (Main.displaySettings.wireframe) {
      var he = f.halfedge;
      var v1 = he.vertex.position;
      var v2 = he.opposite.vertex.position;

      //var midpoint = he.midpoint.add(he.shortTangent);
      var midpoint = new THREE.Vector3();
      midpoint.addVectors(v1, v2);
      midpoint.divideScalar(2);
      const dir = new THREE.Vector3();
      dir.subVectors(midpoint, centroid);
      dir.multiplyScalar(0.9);
      midpoint.addVectors(centroid, dir);

      structureVisGeo.vertices.push(centroid);
      structureVisGeo.vertices.push(midpoint);
      const arrowDirInv = CopyVec(centroid)
        .sub(midpoint)
        .normalize()
        .multiplyScalar(arrowSize);
      structureVisGeo.vertices.push(midpoint);
      structureVisGeo.vertices.push(
        CopyVec(midpoint)
          .add(arrowDirInv)
          .add(he.shortDir)
      );
      structureVisGeo.vertices.push(midpoint);
      structureVisGeo.vertices.push(
        CopyVec(midpoint)
          .add(arrowDirInv)
          .sub(he.shortDir)
      );
    }

    // face idx billboard
    if (Main.displaySettings.showIdLabels) {
      var sprite = makeTextSprite(f.id, {
        borderColor: "black",
        backgroundColor: { r: 255, g: 0, b: 0, a: 1 },
      });
      var pos = CopyVec(f.centroid).add(
        CopyVec(f.normal)
          .normalize()
          .multiplyScalar(0.1)
      );
      sprite.position.set(pos.x, pos.y, pos.z);
      textSprites.push(sprite);
    }

    // face normals
    if (Main.displaySettings.showFN) {
      const fn_p2 = new THREE.Vector3();
      fn_p2.copy(f.normal);
      if (fn_p2.length() > 1e-6) {
        fn_p2.normalize();
      }
      fn_p2.normalize();
      fn_p2.multiplyScalar(0.2);
      fn_p2.add(centroid);

      faceNormalsGeo.vertices.push(centroid);
      faceNormalsGeo.vertices.push(fn_p2);
    }

    // triangulate faces
    const verts = [];
    const f_verts = this.verticesOnFace(f);

    verts[0] = f_verts[0];

    let positions_ptr = vertex_positions;
    let normals_ptr = vertex_normals;
    let colors_ptr = vertex_colors;
    idx = k;

    if (f.selected && !hideSelected) {
      positions_ptr = sel_vertex_positions;
      normals_ptr = sel_vertex_normals;
      colors_ptr = sel_vertex_colors;
      idx = sel_k;
    }

    for (j = 1; j < f_verts.length - 1; ++j) {
      const next_j = j + 1;
      verts[1] = f_verts[j];
      verts[2] = f_verts[next_j];

      for (l = 0; l < 3; ++l) {
        v_pos = verts[l].position;
        if (Main.displaySettings.shading === "smooth") {
          v_nor = verts[l].normal;
        } else {
          v_nor = f.normal;
        }
        v_col = verts[l].color;

        colors_ptr[idx] = v_col.x;
        colors_ptr[idx + 1] = v_col.y;
        colors_ptr[idx + 2] = v_col.z;

        normals_ptr[idx] = v_nor.x;
        normals_ptr[idx + 1] = v_nor.y;
        normals_ptr[idx + 2] = v_nor.z;

        positions_ptr[idx] = v_pos.x;
        positions_ptr[idx + 1] = v_pos.y;
        positions_ptr[idx + 2] = v_pos.z;

        // vertex normals
        if (Main.displaySettings.showVN) {
          const vn_p1 = new THREE.Vector3();
          const vn_p2 = new THREE.Vector3();

          vn_p1.copy(v_pos);
          vn_p2.copy(verts[l].normal);
          if (vn_p2.length() > 1) {
            vn_p2.normalize();
          }
          vn_p2.multiplyScalar(0.2);
          vn_p2.add(vn_p1);
          vertexNormalsGeo.vertices.push(vn_p1);
          vertexNormalsGeo.vertices.push(vn_p2);
        }
        if (f.selected && !hideSelected) {
          sel_k += 3;
          idx = sel_k;
        } else {
          k += 3;
          idx = k;
        }
      }
    }
  }

  for (j = 0; j < this.halfedges.length; ++j) {
    var he = this.halfedges[j];

    // edge arrows
    if (Main.displaySettings.wireframe) {
      wireframeGeo.vertices.push(he.vertex.position);
      wireframeGeo.vertices.push(he.opposite.vertex.position);

      const edgeVec = new THREE.Vector3();
      edgeVec.subVectors(he.vertex.position, he.opposite.vertex.position);
      edgeVec.normalize();
      //var edgeVec = CopyVec(he.shortDir).normalize();
      edgeVec.multiplyScalar(0.1);

      var centroid = he.face.centroid;
      var v1 = he.vertex.position;
      var v2 = he.opposite.vertex.position;

      //var midpoint = he.midpoint.add(he.shortTangent);
      var midpoint = new THREE.Vector3();
      midpoint.addVectors(v1, v2);
      midpoint.divideScalar(2);
      const fDir = new THREE.Vector3();
      fDir.subVectors(midpoint, centroid);
      fDir.multiplyScalar(0.9);
      midpoint.addVectors(centroid, fDir);

      const eDir = new THREE.Vector3();
      eDir.subVectors(v1, v2);
      const p2 = CopyVec(eDir)
        .multiplyScalar(-0.425)
        .add(midpoint);
      const p1 = CopyVec(eDir)
        .multiplyScalar(0.45)
        .add(midpoint);

      // edge points to p1
      //var p1 = CopyVec(he.vertex.position).add(he.shortTangent).sub(edgeVec);
      //var p2 = CopyVec(he.opposite.vertex.position).add(he.shortTangent).add(edgeVec.normalize().multiplyScalar(0.05));
      edgeStructureVisGeo.vertices.push(p1);
      edgeStructureVisGeo.vertices.push(p2);

      const arrowOffsetDir = CopyVec(edgeVec)
        .cross(he.face.normal)
        .normalize()
        .multiplyScalar(arrowSize);
      const p1ArrowBase = CopyVec(p1).sub(edgeVec.normalize().multiplyScalar(arrowSize));
      edgeStructureVisGeo.vertices.push(p1);
      edgeStructureVisGeo.vertices.push(CopyVec(p1ArrowBase).add(arrowOffsetDir));
      edgeStructureVisGeo.vertices.push(p1);
      edgeStructureVisGeo.vertices.push(CopyVec(p1ArrowBase).sub(arrowOffsetDir));
    }

    // edge billboard
    if (Main.displaySettings.showIdLabels) {
      var sprite = makeTextSprite(he.id, {
        borderColor: "black",
        backgroundColor: { r: 255, g: 255, b: 0, a: 1 },
      });
      var pos = CopyVec(he.vertex.position)
        .multiplyScalar(0.66)
        .add(CopyVec(he.opposite.vertex.position).multiplyScalar(0.33))
        .add(
          CopyVec(he.face.normal)
            .normalize()
            .multiplyScalar(0.1)
        );
      sprite.position.set(pos.x, pos.y, pos.z);
      textSprites.push(sprite);
    }
  }

  for (j = 0; j < this.vertices.length; ++j) {
    const vert = this.vertices[j];

    // vertex billboard
    if (Main.displaySettings.showIdLabels) {
      var sprite = makeTextSprite(vert.id, {
        borderColor: "black",
        backgroundColor: { r: 100, g: 100, b: 255, a: 1 },
      });
      // bias it back a bit so selected vertices show over
      var pos = CopyVec(vert.position).add(
        CopyVec(vert.normal)
          .normalize()
          .multiplyScalar(-0.001)
      );
      sprite.position.set(pos.x, pos.y, pos.z);
      textSprites.push(sprite);
    }

    // vertex dots
    if (vert.selected && !hideSelected) {
      selectedVerticesGeo.vertices.push(vert.position);
    } else {
      verticesGeo.vertices.push(vert.position);
      // var vc = vert.color;
      // unselectedVertexColors.push( new THREE.Color(vc.x, vc.y, vc.z) );
      unselectedVertexColors.push(new THREE.Color(255, 255, 255));
    }
  }

  verticesGeo.colors = unselectedVertexColors;

  // assign BufferGeometry buffers
  meshGeometry.addAttribute("position", new THREE.BufferAttribute(vertex_positions, 3));
  meshGeometry.addAttribute("color", new THREE.BufferAttribute(vertex_colors, 3));
  meshGeometry.addAttribute("normal", new THREE.BufferAttribute(vertex_normals, 3));

  selectedFacesGeo.addAttribute("position", new THREE.BufferAttribute(sel_vertex_positions, 3));
  selectedFacesGeo.addAttribute("color", new THREE.BufferAttribute(sel_vertex_colors, 3));
  selectedFacesGeo.addAttribute("normal", new THREE.BufferAttribute(sel_vertex_normals, 3));

  return {
    textSprites: textSprites,
    geometries: {
      main: meshGeometry,
      faceNormals: faceNormalsGeo,
      vertexNormals: vertexNormalsGeo,
      wireframe: wireframeGeo,
      selectedFaces: selectedFacesGeo,
      vertices: verticesGeo,
      selectedVertices: selectedVerticesGeo,
      structureVisualization: structureVisGeo,
      edgeStructureVisualization: edgeStructureVisGeo,
    },
  };
};

Mesh.prototype.toOBJ = function(mesh) {
  function destroyClickedElement(event) {
    document.body.removeChild(event.target);
  }

  let objContent = "# Princeton COS426 OBJ model\n\n";

  for (var i = 0; i < this.vertices.length; ++i) {
    const p = this.vertices[i].position;
    objContent += "v " + p.x + " " + p.y + " " + p.z + "\n";
  }

  for (var i = 0; i < this.faces.length; ++i) {
    objContent += "f";
    const face = this.faces[i];
    const verts = this.verticesOnFace(face);
    for (let j = 0; j < verts.length; ++j) {
      objContent += " " + (verts[j].id + 1);
    }
    objContent += "\n";
  }

  const textFileAsBlob = new Blob([objContent], { type: "text/obj" });
  const fileNameToSaveAs = "mesh.obj";

  const downloadLink = document.createElement("a");
  downloadLink.download = fileNameToSaveAs;
  downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
  downloadLink.onclick = destroyClickedElement;
  downloadLink.style.display = "none";
  document.body.appendChild(downloadLink);
  downloadLink.click();
};

// add event listener that will cause 'O' key to download mesh
window.addEventListener("keyup", function(event) {
  // only respond to 'O' key
  if (event.which == 79) {
    Main.mesh.toOBJ();
  }
});

////////////////////////////////////////////////////////////////////////////////
// Misc Utilities
////////////////////////////////////////////////////////////////////////////////

// from http://stackoverflow.com/questions/15313418/javascript-assert
function assert(condition, message) {
  if (!condition) {
    message = message || "Assertion failed";
    if (typeof Error !== "undefined") {
      throw new Error(message);
    }
    throw message; // Fallback
  }
  return condition;
}

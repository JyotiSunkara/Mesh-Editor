"use strict";

// NOTE: requires a version of dat.gui.js with the following additions:
//  - callbacks on folder open/close
//  - guiFolder.removeFolder(folderName)

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

assert(GuiConfig);

var Gui = Gui || {
  hidden: false,
  filterDefs: GuiConfig.filterDefs,
};

Gui.init = function(controlsChangeCallback) {
  const genFuncAddFilter = function(filterDef) {
    return function() {
      Gui.addHistoryEntry(filterDef);
    };
  };

  // save off callbacks for easy access
  this.controlsChangeCallback = controlsChangeCallback;

  if (this.hidden) {
    this.historyFilters = [];

    this.loadFromUrl();

    this.fullyInitialized = true;

    this.handleControlsChange();

    return;
  }

  // make and save off dat.gui gui objects
  this.historyDatGui = new dat.GUI();
  this.filterListDatGui = new dat.GUI();

  // make a button for each filterDef, in the correct folder
  this.filterFolders = {};
  this.applyFuncs = {};
  for (let filterIdx = 0; filterIdx < this.filterDefs.length; filterIdx++) {
    const filterDef = this.filterDefs[filterIdx];

    // get the gui folder where the button goes (make it if doesn't yet exist)
    if (!filterDef.hidden) {
      let filterFolder = undefined;
      if (filterDef.folderName) {
        filterFolder = this.filterFolders[filterDef.folderName];

        if (filterFolder === undefined) {
          filterFolder = this.filterListDatGui.addFolder(
            filterDef.folderName,
            function() {},
            function() {}
          );
          this.filterFolders[filterDef.folderName] = filterFolder;
        }
      } else {
        filterFolder = this.filterListDatGui;
      }

      // generate the button function and make the button
      filterDef.applyFunc = filterDef.applyFunc || genFuncAddFilter(filterDef);
      this.applyFuncs[filterDef.name] = filterDef.applyFunc;
      const filterButton = filterFolder.add(this.applyFuncs, filterDef.name);
    }
  }

  // CONSTRUCT THE HISTORY PANE
  this.historyFolder = this.historyDatGui.addFolder("History");
  this.historyFolder.open();
  this.historyFilters = [];

  this.fullyInitialized = true;

  this.loadFromUrl();

  if (GuiConfig.onInit) {
    GuiConfig.onInit(this);
  }
};

Gui.handleControlsChange = function() {
  if (!Gui.fullyInitialized || Gui.suspendDisplayUpdate) {
    return;
  }
  this.controlsChangeCallback();
  this.updateUrl();
};

Gui.getFilterHistoryData = function() {
  return this.historyFilters;
};

Gui.getFilterDef = function(filterName) {
  for (let i = 0; i < Gui.filterDefs.length; i++) {
    const def = Gui.filterDefs[i];
    if (def.name === filterName) {
      return def;
    }
  }
  assert(false, "failed to find filter: " + filterName);
};

Gui.addHistoryEntry = function(filterDef, argVals) {
  this.suspendDisplayUpdate = true;

  // make the folder for the history entry
  // HACK: number all history entries to keep folder names unique (required by dat.gui)
  if (!Gui.hidden) {
    const folderLabel = (this.historyFilters.length + 1).toString() + ": " + filterDef.name;
    var folder = this.historyFolder.addFolder(folderLabel);
    folder.open();
  }

  // make the filter instance
  const filterInst = {
    filterDef: filterDef,
    paramVals: {},
    argsList: [],
    guiControls: [],
    folder: folder,
  };

  let animatedValFound = false;

  // make param controls in gui folder
  for (let paramIdx = 0; paramIdx < filterDef.paramDefs.length; paramIdx++) {
    const paramDef = filterDef.paramDefs[paramIdx];

    var paramVal = paramDef.defaultVal;
    if (argVals !== undefined && argVals[paramIdx] !== undefined)
      paramVal = argVals[paramIdx];

    // booleans and undefined become empty string
    if (paramDef.isString) {
      if (typeof paramVal !== "string") {
        if (typeof paramVal === "number") {
          paramVal = paramVal.toString();
        } else {
          paramVal = "";
        }
      }
    }

    // HACK(drew): dat.gui needs to be initialized with a float for float sliders
    filterInst.paramVals[paramDef.name] = paramDef.isFloat ? 0.5 : paramVal;

    filterInst.argsList[paramIdx] = paramVal;

    // mess of nested functions to generate a valid closure
    const defaultOnChangeFunc = (function(paramIdx) {
      return function(val) {
        filterInst.argsList[paramIdx] = val;
        Gui.handleControlsChange();
      };
    })(paramIdx);

    if (Gui.hidden) {
      this.suspendDisplayUpdate = false;
      if (paramDef.isColor) {
        defaultOnChangeFunc(new Pixel(paramVal[0] / 255, paramVal[1] / 255, paramVal[2] / 255, 1));
      } else {
        if (paramVal.isAnimated) {
          assert(!animatedValFound, "only one animated parameter is allowed at a time");
          animatedValFound = true;

          Main.animatedValData = {
            current: parseFloat(paramVal.start),
            start: parseFloat(paramVal.start),
            end: parseFloat(paramVal.end),
            step: parseFloat(paramVal.step),
            changeFunc: defaultOnChangeFunc,
          };

          paramVal = parseFloat(paramVal.start);
        }

        defaultOnChangeFunc(paramVal);
      }
    } else {
      let paramControl = undefined;
      if (paramDef.sliderRange) {
        paramControl = folder.add(
          filterInst.paramVals,
          paramDef.name,
          paramDef.sliderRange[0],
          paramDef.sliderRange[1]
        );
        paramControl.step(
          paramDef.step ||
            (!paramDef.isFloat && 1) ||
            (paramDef.sliderRange[1] - paramDef.sliderRange[0]) / 50
        );

        paramControl.onChange(defaultOnChangeFunc);
      } else if (paramDef.dropdownOptions) {
        paramControl = folder.add(filterInst.paramVals, paramDef.name, paramDef.dropdownOptions);
        paramControl.onChange(defaultOnChangeFunc);
      } else if (paramDef.isColor) {
        paramControl = folder.addColor(filterInst.paramVals, paramDef.name);
        paramControl.onChange(
          (function(filterInst, paramIdx) {
            return function(val) {
              filterInst.argsList[paramIdx] = new Pixel(
                val[0] / 255,
                val[1] / 255,
                val[2] / 255,
                1
              );
              Gui.handleControlsChange();
            };
          })(filterInst, paramIdx)
        );
      } else if (paramDef.isString) {
        paramControl = folder.add(filterInst.paramVals, paramDef.name);
        paramControl.onChange(defaultOnChangeFunc);
      } else if (paramDef.isBoolean) {
        paramControl = folder.add(filterInst.paramVals, paramDef.name);
        paramControl.onChange(defaultOnChangeFunc);
      } else {
        assert(false, "unsupported control type in paramDef: " + paramDef);
      }

      assert(!paramVal.isAnimated, "animated parameters are only allowed in batch mode");

      // in case we initialized with a float dummy value, set default again
      paramControl.setValue(paramVal);

      filterInst.guiControls[paramIdx] = paramControl;
    }
  }

  if (!Gui.hidden) {
    if (filterDef.permanent) {
      // no delete button
    } else if (filterDef.hasDeleteBelowButton) {
      filterInst["Delete Below"] = function() {
        Gui.deleteHistoryBelow(filterInst);
        Gui.deleteHistoryEntry(filterInst);
        Gui.handleControlsChange();
      };
      folder.add(filterInst, "Delete Below");
    } else {
      filterInst.Delete = function() {
        Gui.deleteHistoryEntry(filterInst);
        Gui.handleControlsChange();
      };
      folder.add(filterInst, "Delete");
    }

    this.historyFolder.open();
  }

  this.historyFilters.push(filterInst);

  this.suspendDisplayUpdate = false;
  Gui.handleControlsChange();

  return filterInst;
};

Gui.deleteHistoryEntry = function(filterInst) {
  for (let i = 0; i < this.historyFilters.length; i++) {
    const maybeInputFilterInst = this.historyFilters[i];

    if (maybeInputFilterInst === filterInst) {
      this.historyFolder.removeFolder(filterInst.folder.name);
      this.historyFilters.splice(i, 1);
      this.updateUrl();
      return;
    }
  }
};

Gui.deleteHistoryBelow = function(filterInst) {
  for (let i = 0; i < this.historyFilters.length; i++) {
    const maybeInputFilterInst = this.historyFilters[i];

    if (maybeInputFilterInst === filterInst) {
      for (let j = i + 1; j < this.historyFilters.length; j++) {
        const filterInstToRemove = this.historyFilters[j];
        this.historyFolder.removeFolder(filterInstToRemove.folder.name);
      }
      this.historyFilters.splice(i + 1, this.historyFilters.length - (i + 1));
      this.updateUrl();
      return;
    }
  }
};

// gets rid of the ".0000000001" etc when stringifying floats
// from http://stackoverflow.com/questions/1458633/how-to-deal-with-floating-point-number-precision-in-javascript
function stripFloatError(number) {
  if (number && number.toPrecision) {
    return parseFloat(number.toPrecision(12));
  } else {
    return number;
  }
}

Gui.loadFromUrl = function() {
  for (let i = 0; i < Parser.commands.length; i++) {
    const cmd = Parser.commands[i];

    if (cmd.name == "imageFile") {
      const fileName = cmd.args[0];
    } else {
      let filterDef = undefined;
      for (let filterIdx = 0; filterIdx < this.filterDefs.length; filterIdx++) {
        const curFilterDef = this.filterDefs[filterIdx];
        if (curFilterDef.name === cmd.name) {
          filterDef = curFilterDef;
          break;
        }
      }

      if (filterDef !== undefined) {
        // is known command
        this.addHistoryEntry(filterDef, cmd.args);
      }
    }
  }
};

Gui.getUrl = function() {
  let url = "";
  for (let i = 0; i < this.historyFilters.length; i++) {
    const filterInst = this.historyFilters[i];

    // filter separator
    url += (i > 0 && "&") || "";

    // filter name
    url += filterInst.filterDef.name + "=";

    for (let j = 0; j < filterInst.argsList.length; j++) {
      // parameter separator
      if (j > 0) {
        url += ";";
      }
      const val = filterInst.argsList[j];

      if (isNaN(val) || val === "") {
        if (typeof Pixel != "undefined" && Pixel.prototype.isPrototypeOf(val)) {
          url +=
            "[" +
            stripFloatError(val.data[0]) +
            "," +
            stripFloatError(val.data[1]) +
            "," +
            stripFloatError(val.data[2]) +
            "," +
            stripFloatError(val.a) +
            "]";
        } else {
          url += val;
        }
      } else {
        url += stripFloatError(val);
      }
    }
  }
  url = url.replace(/ /g, "_");
  return url;
};

var stateObj = {};
Gui.updateUrl = function() {
  if (Gui.batchMode) {
    return;
  }

  let url = (Gui.batchMode && "batch.html?") || "index.html?";
  url += Gui.getUrl();
  history.pushState(stateObj, "", url);
};

Gui.alertOnce = function(msg) {
  const mainDiv = document.getElementById("main_div");
  mainDiv.style.opacity = "0.3";
  const alertDiv = document.getElementById("alert_div");
  alertDiv.innerHTML = "<p>" + msg + '</p><button id="ok" onclick="Gui.closeAlert()">ok</button>';
  alertDiv.style.display = "inline";
};

Gui.closeAlert = function() {
  const mainDiv = document.getElementById("main_div");
  mainDiv.style.opacity = "1";
  const alertDiv = document.getElementById("alert_div");
  alertDiv.style.display = "none";
};

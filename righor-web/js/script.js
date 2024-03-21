import init, { WrappedModel } from "../pkg/righor_web.js";

var column_index = {
  'sequence': 0, 'pgen': 5, 'ncdr3': 1, 'v': 2, 'j':3, 'acdr3': 4
}

var results = {};
var global = {not_a_click:false, table: undefined, currentSpeciesId: 'human'
              //, worker: undefined
             };




function alignReprString(str1, str2) {
    if (str1.length !== str2.length) {
        throw new Error("Strings must be of the same length.");
    }

    let result = "";
    for (let i = 0; i < str1.length; i++) {
        if (str1[i] === str2[i]) {
            result += '.';
        } else {
            result += str2[i];
        }
    }
    return result;
}

function switch_model(species, chain) {
  WrappedModel.load_model(species, chain)
}


function generateWrapper(n) {
  let promises = [];
  const checkbox = document.getElementById('checkBoxFunctional');
  for (let step = 0; step < n; step++){
    promises.push(WrappedModel.generate(checkbox.checked).then(
      (rjs) => {update_table(rjs, step, false)}));
  }
  return Promise.all(promises);
}


function generate(n) {
  document.getElementById('loadingWheel').style.visibility = 'visible';
  setTimeout(() => {
    generateWrapper(n).then(document.getElementById('loadingWheel').style.visibility = 'hidden')
  }, 100);

}


async function analyze_sequence(rowIndex) {
    var rowData = global.table.getRowData(rowIndex)
    // full sequence
    if (rowData[0].trim() != "") {
      var seq = rowData[0].trim()
      if (!(seq in results)) {
        try {
          console.log(seq)
          console.log(WrappedModel.evaluate(seq))
          return WrappedModel.evaluate(seq).then((res) => {  results[seq] = res; return update_table(res, rowIndex, true)}).catch((e) =>{ console.log(e);})

        } catch(e) { // in case something goes wrong (avoid a freeze)
          console.log(e)
        }
      }

    }
    // just cdr3 nt and V/J genes
  else if (rowData[1].trim() != "" && rowData[2].trim() != "" && rowData[3].trim() != "") {
    var seq = rowData[1].trim();
    var vname = rowData[2].trim();
    var jname = rowData[3].trim();
    try {
      return WrappedModel.evaluate_from_cdr3(seq, vname, jname).then((res) => { results[seq] = res; return update_table(res, rowIndex, true)}).catch((e) =>{ console.log(e);})

    } catch(e) { // in case something goes wrong (avoid a freeze)
      console.log(e)
    }
  }
  // just cdr3 aa and V/J genes
  else if (rowData[4].trim() != "" && rowData[2].trim() != "" && rowData[3].trim() != "") {
    console.log("Not implemented yet")
  }
}



function update_table(r, idx, with_pgen) {
  var currentRowCount = global.table.getData().length;
  if (idx >= currentRowCount) {
    // Add a new row. Adjust parameters according to your needs
    global.table.insertRow([]);
  }

  if(with_pgen){
    global.table.setValueFromCoords(column_index['pgen'], idx, r.pgen.toExponential(2));
    global.table.getCellFromCoords(column_index['pgen'], idx).setAttribute("style", "background-color:#e6fce7;")
  }
  global.table.setValueFromCoords(column_index['sequence'], idx, r.full_seq);
  global.table.setValueFromCoords(column_index['j'], idx, r.j_name.split('-')[0].split('*')[0]);
  global.table.setValueFromCoords(column_index['v'], idx, r.v_name.split('-')[0].split('*')[0]);
  global.table.setValueFromCoords(column_index['ncdr3'], idx, r.n_cdr3);

  if (r.aa_cdr3 === undefined || r.aa_cdr3.length == 0) {
    global.table.setValueFromCoords(column_index['acdr3'], idx, '<em> Out of Frame </em>');
  }
  else {
    global.table.setValueFromCoords(column_index['acdr3'], idx, r.aa_cdr3);
  }

  printSequence(idx);
  global.not_a_click = true;

  // remove gray color, in case the line contained the og example
  if (idx == 0) {
    global.table.rows[0].style.color = '';
  }
}


function load_menu() {

  const speciesButton = document.getElementById('modelChoiceLink');
  const speciesMenu = document.getElementById('speciesMenu');
  const chainsMenu = document.getElementById('chainsMenu');
  const selectionDisplay = document.getElementById('modelName'); // Get the display div

  let currentSpecies = ''; // Variable to keep track of the current species

  // Mapping of species to chains
  const speciesToChains = {
    'human': ['TRα', 'TRβ', 'IGH', 'IGλ', 'IGκ'],
    'mouse': ['TRα', 'TRβ'],
    // Add more mappings as needed
  };

  const chainToId = {
    'TRα': 't_alpha',
    'TRβ': 't_beta',
    'IGH': 'b_heavy',
    'IGλ': 'b_lambda',
    'IGκ': 'b_lambda'
  }

  // Toggle species menu display
  speciesButton.addEventListener('click', () => {
    speciesMenu.style.display = speciesMenu.style.display === 'none' ? 'block' : 'none';
    chainsMenu.style.display = 'none'; // Hide chains menu if species menu is toggled
  });

  // Show chains menu for selected species
  speciesMenu.addEventListener('click', e => {
    if (e.target.dataset.species) {
      currentSpecies = e.target.textContent; // Store the current species name
      global.currentSpeciesId = e.target.dataset.species;
      const chains = speciesToChains[e.target.dataset.species];
      chainsMenu.innerHTML = ''; // Clear previous chains
      chains.forEach(chain => {
        const div = document.createElement('div');
        div.textContent = chain;
        chainsMenu.appendChild(div);
      });

      // Calculate and set position
      const speciesRect = e.target.getBoundingClientRect();
      const menuRect = speciesMenu.getBoundingClientRect();

      chainsMenu.style.display = 'block'; // Ensure the menu is displayed before positioning
      chainsMenu.style.left = `${menuRect.right}px`; // Position to the right of the speciesMenu
      chainsMenu.style.top = `${speciesRect.top + window.scrollY}px`; // Align with the clicked species

    }
  });

  // Event listener for clicks on chains to update the display div
  chainsMenu.addEventListener('click', e => {
    if (e.target.textContent) { // Ensure the click is on a chain option
      const chainSelected = e.target.textContent;
      selectionDisplay.textContent = `${currentSpecies} / ${chainSelected}`; // Update the display div
      switch_model(global.currentSpeciesId, chainToId[e.target.textContent]);
      speciesMenu.style.display = 'none';
      chainsMenu.style.display = 'none';
    }
  });

  // Optional: Hide menus when clicking outside
  document.addEventListener('click', (e) => {
    if (!e.target.matches('#modelChoiceLink, #speciesMenu, #speciesMenu div, #chainsMenu, #chainsMenu div')) {
      speciesMenu.style.display = 'none';
    }
    if (!e.target.matches('#speciesMenu div, #chainsMenu, #chainsMenu div')) {
      chainsMenu.style.display = 'none';
    }
  }, true);

}

function choose_evaluate() {
  let promises = [];
    // for each occupied / selected row
    let selected = global.table.selectedCell;
    console.log(selected)

    if ( selected != null) {
      // Iterate over selected rows
      for(let rowIndex = selected[1]; rowIndex <= selected[3]; rowIndex++){
        promises.push(analyze_sequence(rowIndex));
      }
    } else {
      // No selection, just do the first row
      let rowCount = global.table.getData().length;
      for (let rowIndex = 0; rowIndex < 1; rowIndex++) {
        promises.push(analyze_sequence(rowIndex));
      }
    }
  return Promise.all(promises);
}


async function evaluate_click() {

  document.getElementById('loadingWheel').style.visibility = 'visible';
  setTimeout(() => {
    choose_evaluate().then(() => {
      document.getElementById('loadingWheel').style.visibility = 'hidden';
    }).catch(() => { document.getElementById('loadingWheel').style.visibility = 'hidden';});
  }, 100);
}

function onRowClick(instance, colNumber, rowNumber, cells, e) {
  if (global.not_a_click) {
    global.not_a_click = false;
    return;
  }
  // Remove highlighting from all rows
  Array.from(document.querySelectorAll('#spreadsheet .jexcel tbody tr')).forEach(function(row) {
      row.style.backgroundColor = ''; // Reset background color
  });
  // remove gray color, in case the line contained the og example
  if (rowNumber == 0) {
    instance.jexcel.rows[0].style.color = '';
  }


  // Highlight the clicked row
  var clickedRow = instance.jexcel.rows[rowNumber];
  if (clickedRow) {
    clickedRow.style.backgroundColor = '#D3E5FF'; // Set your highlight color
  }

  printSequence(rowNumber)
}

function eraseSequence() {
  var div = document.getElementById('alignment');
  div.innerHTML = "";
}

function printSequence(rowNumber) {
  eraseSequence();
  var data = global.table.getRowData(rowNumber)
  var div = document.getElementById('alignment');
  console.log(data)
  if (data[0].trim() in results) {
    var obj = results[data[0].trim()];
    div.innerHTML += "<span style='color:#0056b3;'>" + obj.cdr3_pos_string + "</span>\n"
    div.innerHTML += obj.seq + "\n"
    div.innerHTML += alignReprString(obj.seq, obj.aligned_v) + "\n"
    div.innerHTML += alignReprString(obj.seq, obj.aligned_j) + "\n"
    // scroll to the cdr3

    // console.log(div.scrollWidth);
    // console.log(obj.pos_start_cdr3);
    // console.log(obj.seq)
    div.scrollLeft = div.scrollWidth * (obj.pos_start_cdr3/obj.seq.length) - 30;
  }
}

document.addEventListener('DOMContentLoaded', () => {



  // global.worker = new Worker('js/wasm_worker.js', { type: 'module' });
  // global.worker.onerror = function(event) {
//   console.error("Worker error event:", event);
//   // Log specific properties that might contain useful information
//   console.error(`Error filename: ${event.filename}, lineno: ${event.lineno}, colno: ${event.colno}`);
// };

  init().then(() => {
    WrappedModel.initialize();


    var example_data = [
      ["TCTCAGACTATTCATCAATGGCCAGCGACCCTGGTGCAGCCTGTGGGCAGCCCGCTCTCTCTGGAGTGCACTGTGGAGGGAACATCAAACCCCAACCTATACTGGTACCGACAGGCTGCAGGCAGGGGCCTCCAGCTGCTCTTCTACTCCGTTGGTATTGGCCAGATCAGCTCTGAGGTGCCCCAGAATCTCTCAGCCTCCAGACCCCAGGACCGGCAGTTCATCCTGAGTTCTAAGAAGCTCCTTCTCAGTGACTCTGGCTTCTATCTCTGTGCCTGGAATGGGACTAGCGGGGGGAAGGTTTCTGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGG", "TGTGCCTGGAATGGGACTAGCGGGGGGAAGGTTTCTGAAAAACTGTTTTTT  ", "TRBV30", "TRBJ1", "CAWNGTSGGKVSEKLFF", ""],

        ["", "", "", "", "", ""]

      ]


      global.table = jspreadsheet(document.getElementById('spreadsheet'), {
        data: example_data,
        onselection: onRowClick,
        columns: [
          {
            type: 'text',
            title: 'Sequence',
            tooltip: 'Complete V(D)J sequence.',
	    width: '185px',
          },
          {
            type: 'text',
            title: 'nCDR3',
            tooltip: 'CDR3 sequence (in nucleotides).',
	    width: '150px',
          },
          {
            type: 'text',
            title: 'V family',
            tooltip: 'Most likely V (variable) gene family.',
	    width: '100px',
          },
          {
            type: 'text',
            title: 'J family',
            tooltip: 'Most likely J (joining) gene family.',
	    width: '100px',
          },
          {
            type: 'html',
            title: 'aaCDR3',
            tooltip: 'CDR3 sequence (in amino-acids).',
	    width: '150px',
          },
          {
            type: 'numeric',
            title: 'pgen',
            tooltip: 'Probability of being generated according to the model.',
	    width: '100px',
          }],
        minDimensions: [6, 4],
        tableWidth: "800px",
        tableHeight: "300px",
      });
      // easier to debug
      window.table = global.table

      // Apply gray background for the examples
      var numCols = global.table.getHeaders().length;
      global.table.rows[0].style.color = '#cccccc';

      document.getElementById('generateLink').addEventListener('click', function() {
        var inputValue = document.getElementById('nbGeneratedInput').value;
        generate(inputValue);
      });

      document.getElementById('nbGeneratedInput').addEventListener('click', function() {
        global.table.resetSelection();
      });

      document.getElementById('evaluateLink').addEventListener('click', function() {
        evaluate_click();
      });
    document.getElementById('loadingWheel').style.visibility = 'hidden';
    load_menu();

    document.getElementById('help-window').style.display = 'none'

    document.getElementById('help-button').addEventListener('click', function() {
      document.getElementById('help-window').style.display = 'block'
    });

    document.addEventListener('click', function(event) {
      var isClickInside = document.getElementById('help-window').contains(event.target);
      var isClickButton = document.getElementById('help-button').contains(event.target);

      if (!isClickInside && !isClickButton) {
        document.getElementById('help-window').style.display = 'none'
      }
    });

  });

});


// function sendTaskToWorker(task, data = {}) { // Default parameter in case no data is passed
//   return new Promise((resolve, reject) => {
//     global.worker.onmessage = (e) => {
//       const { status, result, message } = e.data;

//       if (status === 'success') {
//         console.log(result);
//         resolve(result);
//       } else {
//         reject(new Error(message));
//       }
//     };

//     global.worker.postMessage({ task, ...data });
//   });
// }

// // Specific functions for each task
// function loadModelInWorker(a, b) {
//   return sendTaskToWorker('loadModel', { a, b });
// }

// function evaluateInWorker(seq) {
//   return sendTaskToWorker('evaluate', { seq });
// }

// function generateInWorker(functional) {
//   return sendTaskToWorker('generate', { functional });
// }

// function initializeInWorker() {
//   return sendTaskToWorker('initialize');
// }

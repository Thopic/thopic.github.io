import init, { WrappedModel } from "../pkg/righor_web.js";

var column_index = {
  'sequence': 0, 'pgen': 5, 'ncdr3': 1, 'v': 2, 'j':3, 'acdr3': 4
}

var results = {};
var global = {not_a_click:false, table: undefined };



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

function generate(n) {
  for (let step = 0; step < n; step++){
    let rjs = WrappedModel.generate(step)
    update_table(rjs, step, false)
  }
}


function analyze_sequence(rowIndex) {
  var rowData = global.table.getRowData(rowIndex)
  // full sequence
  if (rowData[0].trim() != "") {
    var seq = rowData[0].trim()
    if (!(seq in results)) {
      try {
        results[seq] = WrappedModel.evaluate(seq);
      } catch(e) { // in case something goes wrong (avoid a freeze)
        console.log(e)
        return
      }
      update_table(results[seq], rowIndex, true)
    }
  }
  // just cdr3 nt and V/J genes
  else if (rowData[1].trim() != "" && rowData[2].trim() != "" && rowData[3].trim() != "") {
    console.log("Not implemented yet")

  }
  // just cdr3 aa and V/J genes
  else if (rowData[4].trim() != "" && rowData[2].trim() != "" && rowData[3].trim() != "") {
    console.log("Not implemented yet")
  }
}


function update_table(r, idx, with_pgen) {
  console.log(r)
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


function evaluate_click() {

  // for each occupied / selected row

  let selected = global.table.selectedCell;
  console.log(selected)

  if ( selected != null) {
    // Iterate over selected rows
    for(let rowIndex = selected[1]; rowIndex <= selected[3]; rowIndex++){
      analyze_sequence(rowIndex);
    }
  } else {
    // No selection, iterate over all rows
    let rowCount = global.table.getData().length;
    for (let rowIndex = 0; rowIndex < rowCount; rowIndex++) {
      analyze_sequence(rowIndex);
    }
  }
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
  if (data[0].trim() in results) {
    var obj = results[data[0].trim()];
    console.log(obj.aligned_j);
    div.innerHTML += "<span style='color:#0056b3;'>" + obj.cdr3_pos_string + "</span>\n"
    div.innerHTML += obj.seq + "\n"
    div.innerHTML += alignReprString(obj.seq, obj.aligned_v) + "\n"
    div.innerHTML += alignReprString(obj.seq, obj.aligned_j) + "\n"
    // scroll to the cdr3

    console.log(div.scrollWidth);
    console.log(obj.pos_start_cdr3);
    console.log(obj.seq)
    div.scrollLeft = div.scrollWidth * (obj.pos_start_cdr3/obj.seq.length) - 30;
  }
}

document.addEventListener('DOMContentLoaded', () => {
    init().then(() => {
      WrappedModel.initialize()

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

      document.getElementById('evaluateLink').addEventListener('click', function() {
        evaluate_click();
      });



    });

});

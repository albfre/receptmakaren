function solveQP(Q, c, Aeq, beq, Aineq, bineq, variables = []) {
  let solutionElement = document.getElementById("status");

  const start = performance.now();
  const {x, f, res, gap, iter} = interiorPointQP(Q, c, Aeq, beq, Aineq, bineq);
  const end = performance.now();

  let tableStr = '<table>';
  function addRow(str, val) {
    tableStr += `<tr><td>${str}</td><td>${val}</td></tr>`;
  }

  addRow('Objective value', f);
  addRow('Number of iterations', iter);
  addRow('Residual', res);
  addRow('Gap', gap);
  addRow('Elapsed time', `${end - start} milliseconds`);
  for (let i = 0; i < x.length; i++) {
    addRow(variables.length == x.length ? variables[i] : `x${i}`, x[i]);
  }
  addRow('Variable vector', x);
  tableStr += '</table>';

  solutionElement.innerHTML = tableStr;
  return x;
}

// Functions relating to buttons on the html page
const qp_consideredNaringsvarden = ["Energi (kJ)", "Fett", "Summa mättade fettsyror", "Kolhydrater", "Socker totalt", "Fibrer", "Protein", "Salt"];
const shortNaringsvarden = ["energi", "fett", "mattat-fett", "kolhydrat", "socker", "fibrer", "protein", "salt"];
//const qp_consideredNaringsvarden = ["Energi (kJ)"];
//const shortNaringsvarden = ["energi"];
const naringsvardenMap = {}
for (let i = 0; i < qp_consideredNaringsvarden.length; i++) {
  naringsvardenMap[qp_consideredNaringsvarden[i]] = shortNaringsvarden[i];
}

function parseSelectedFoodOption(option) {
  console.log(option.value); // or option.text to get the text content of the option
  const values = option.dataset.value.split(",");
  if (values.length % 2 != 0) {
    throw new Error('Value error: ' + values);
  }
  const varden = {};
  for (let i = 0; i < values.length; i += 2) {
    const name = values[i].trim();
    const value = parseFloat(values[i + 1]);
    if (name in naringsvardenMap) {
      varden[naringsvardenMap[name]] = value;
    }
    else {
      //throw new Error('Unknown name: ' + name);
    }
  }
  return varden;
}

function parseSelectedFoods() {
  if (qp_consideredNaringsvarden.length != shortNaringsvarden.length) {
    throw new Error('consideredNaringsvarden.length != shortNaringsvarden.length: ' + qp_consideredNaringsvarden + ", " + shortNaringsvarden);
  }

  const naringsvarden = {};
  const selectedFoodsList = document.querySelector("#selected-foods");
  for (const option of selectedFoodsList.options) {
    naringsvarden[option.value] = parseSelectedFoodOption(option);
  }
  return naringsvarden;
}

function parseInput() {
  const naringsvarden = {};
  for (const naringsvarde of shortNaringsvarden) {
    const targetValue = parseFloat(document.getElementById(naringsvarde).value);
    naringsvarden[naringsvarde] = targetValue
  }
  return naringsvarden;
}

function setOutput(naringsvarden, prefix) {
  for (const naringsvarde of shortNaringsvarden) {
    const value = naringsvarden[naringsvarde];
    const target = document.getElementById(`${prefix}-${naringsvarde}`);
    target.textContent = value;
  }
}

function solve() {
  const selectedFoodsList = document.querySelector("#selected-foods");
  const options = selectedFoodsList.options;
  const foods = zeroVector(options.length);
  for (let i = 0; i < options.length; i++) {
    foods[i] = options[i].value;
  }
  const n = options.length;
  let Q = zeroMatrix(n, n);
  const c = zeroVector(n);

  // e' x = 1
  const Aeq = filledMatrix(1, n, 1.0);
  const beq = filledVector(1, 1.0);

  // x >= 0
  const Aineq = diag(filledVector(n, 1.0));
  const bineq = zeroVector(n);

  // (t - k0 x0 - k1 x1 ...)^2 = t^2 - 2 t ki xi + ki^2 xi^2 + 2 ki kj xi xj
  const selectedFoodsNaringsvarden = parseSelectedFoods();
  console.log(selectedFoodsNaringsvarden);
  const targetNaringsvarden = parseInput();
  for (const naringsvarde of shortNaringsvarden) {
    const target = targetNaringsvarden[naringsvarde];
    for (let i = 0; i < foods.length; i++) {
      const ki = selectedFoodsNaringsvarden[foods[i]][naringsvarde];
      c[i] -= target * ki / target ** 2;
      Q[i][i] += ki ** 2 / target ** 2;
      for (let j = i + 1; j < foods.length; j++) {
        const kj = selectedFoodsNaringsvarden[foods[j]][naringsvarde];
        Q[i][j] += ki * kj / target ** 2;
        Q[j][i] += ki * kj / target ** 2;
      }
    }
  }
  console.log('Q')
  console.log(Q)
  //Q = diag(filledVector(n, 1));

  try {
    const x = solveQP(Q, c, Aeq, beq, Aineq, bineq);
    console.log(x)
    const naringsvarden = {};
    for (const naringsvarde of shortNaringsvarden) {
      naringsvarden[naringsvarde] = 0.0;
    }
    
    for (const naringsvarde of shortNaringsvarden) {
      for (let i = 0; i < foods.length; i++) {
        naringsvarden[naringsvarde] += selectedFoodsNaringsvarden[foods[i]][naringsvarde] * x[i];
      }
    }
    setOutput(naringsvarden, 'result');
  } catch (error) {
    let solutionElement = document.getElementById("status");
    solutionElement.innerHTML = `Error ${error.lineNumber}: ${error.message}`;
  }
}

function printResults() {
  const selectedFoodsList = document.querySelector("#selected-foods");
  for (const food of selectedFoods) {
    const listItem = document.createElement("li");
    listItem.textContent = `${food.name}: Energy=${food.energy}, Fat=${food.fat}, Sugar=${food.sugar}`;
    selectedFoodsList.appendChild(listItem);
  }
}

document.getElementById("optimize").addEventListener("click", solve);

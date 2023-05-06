function solveQP(Q, c, Aeq, beq, Aineq, bineq, variables = []) {
  let solutionElement = document.getElementById("status");

  const start = performance.now();
  const {x, f, res, gap, iter} = interiorPointQP(Q, c, Aeq, beq, Aineq, bineq);
  const end = performance.now();

  let tableStr = "<table>";
  function addRow(str, val) {
    tableStr += `<tr><td>${str}</td><td>${val}</td></tr>`;
  }

  addRow("Objective value", f);
  addRow("Number of iterations", iter);
  addRow("Residual", res);
  addRow("Gap", gap);
  addRow("Elapsed time", `${end - start} milliseconds`);
  for (let i = 0; i < x.length; i++) {
    addRow(variables.length == x.length ? variables[i] : `x${i}`, x[i]);
  }
  addRow("Variable vector", x);
  tableStr += "</table>";

  solutionElement.innerHTML = tableStr;
  return x;
}

// Functions relating to buttons on the html page
const naringsvardenMap = {}
for (let i = 0; i < consideredNaringsvarden.length; i++) {
  naringsvardenMap[consideredNaringsvarden[i]] = shortNaringsvarden[i];
}

function parseSelectedFoodOption(option) {
  //console.log(option);
  const values = option.dataset.value.split(",");
  if (values.length % 2 != 0) {
    throw new Error("Value error: " + values);
  }
  const varden = {};
  for (let i = 0; i < values.length; i += 2) {
    const name = values[i].trim();
    const value = parseFloat(values[i + 1]);
    if (name in naringsvardenMap) {
      varden[naringsvardenMap[name]] = value;
    }
    else {
      //throw new Error("Unknown name: " + name);
    }
  }
  return varden;
}

function parseSelectedFoods() {
  if (consideredNaringsvarden.length != shortNaringsvarden.length) {
    throw new Error("consideredNaringsvarden.length != shortNaringsvarden.length: " + consideredNaringsvarden + ", " + shortNaringsvarden);
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
    const str = document.getElementById("mimic-" + naringsvarde).value.trim();
    const targetValue = str.length === 0 ? 0.0 : parseFloat(str);
    naringsvarden[naringsvarde] = targetValue
  }
  return naringsvarden;
}

function setOutput(naringsvarden, prefix) {
  for (const naringsvarde of shortNaringsvarden) {
    const value = naringsvarden[naringsvarde];
    const formattedValue = value.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2, useGrouping: false });
    const target = document.getElementById(`${prefix}-${naringsvarde}`);
    target.textContent = formattedValue;
  }
}

function solve(decreasing = false) {
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
  let Aineq = diag(filledVector(n, 1.0));
  let bineq = zeroVector(n);
  if (decreasing) {
    // x[i] >= x[i + 1] => (1 -1 0 ... ) x >= 0
    for (let i = 0; i + 1 < options.length; i++) {
      const row = zeroVector(n);
      row[i] = 1;
      row[i + 1] = -1;
      Aineq.push(row);
      bineq.push(0.0);
    }
  }
  console.log(Aineq);

  // (t - k0 x0 - k1 x1 ...)^2 = t^2 - 2 t ki xi + ki^2 xi^2 + 2 ki kj xi xj
  const selectedFoodsNaringsvarden = parseSelectedFoods();
  const targetNaringsvarden = parseInput();
  //console.log("Näringsämnen");
  //console.log(selectedFoodsNaringsvarden);
  //console.log(targetNaringsvarden);
  for (const naringsvarde of shortNaringsvarden) {
    const target = targetNaringsvarden[naringsvarde];
    const denominator = Math.max(target, 0.1) ** 2;
    for (let i = 0; i < foods.length; i++) {
      const ki = selectedFoodsNaringsvarden[foods[i]][naringsvarde];
      c[i] -= target * ki / denominator;
      Q[i][i] += ki ** 2 / denominator;
      for (let j = i + 1; j < foods.length; j++) {
        const kj = selectedFoodsNaringsvarden[foods[j]][naringsvarde];
        Q[i][j] += ki * kj / denominator;
        Q[j][i] += ki * kj / denominator;
      }
    }
  }
  console.log("Q");
  console.log(Q);

  try {
    const x = solveQP(Q, c, Aeq, beq, Aineq, bineq);
    console.log(x);
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
    printResults(options, x);
  } catch (error) {
    let solutionElement = document.getElementById("status");
    solutionElement.innerHTML = `Error ${error.lineNumber}: ${error.message}`;
  }
}

function solveWithoutOrder() {
  console.log("apa");
  solve(false);
}

function solveDecreasing() {
  solve(true);
}

function printResults(options, weights) {
  const weightList = document.querySelector("#resulting-weights");
  while (weightList.firstChild) {
    weightList.removeChild(weightList.firstChild);
  }

  for (let i = 0; i < options.length; i++) {
    const listItem = document.createElement("li");
    const formattedNumber = weights[i].toLocaleString(undefined, { minimumFractionDigits: 3, maximumFractionDigits: 3, useGrouping: false });
    listItem.textContent = `${options[i].value}: ${formattedNumber}`;
    weightList.appendChild(listItem);
  }
}

document.getElementById("optimize").addEventListener("click", solveWithoutOrder);
document.getElementById("optimize-decreasing").addEventListener("click", solveDecreasing);

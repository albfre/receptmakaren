function interiorPointQP(H, c, Aeq, beq, Aineq, bineq, tol=1e-8, maxIter=500) {
  /* minimize 0.5 x' H x + c' x
   *   st    Aeq x = beq
   *         Aineq x >= bineq
   */

  // Matrix sizes
  const n = H.length;
  const mIneq = Aineq.length;
  const mEq = Aeq.length;

  // Preconditions
  if (H.some(row => row.length != n)) {
    throw new Error('H is not a square matrix');
  }
  if (Aineq.some(row => row.length != n)) {
    throw new Error('All rows of Aineq must have the same length as H');
  }
  if (Aeq.some(row => row.length != n)) {
    throw new Error('All rows of Aeq must have the same length as H');
  }
  if (bineq.length !== mIneq) {
    throw new Error('Aineq and bineq must have the same length. Aineq.length = ' + mIneq + ', bineq.length = ' + bineq.length);
  }
  if (beq.length !== mEq) {
    throw new Error('Aeq and beq must have the same length. Aeq.length = ' + mEq + ', beq.length = ' + beq.length);
  }

  const AineqT = transpose(Aineq);
  const AeqT = transpose(Aeq);

  // Define the function for evaluating the objective and constraints
  function evalFunc(x, s, y, z, mu) {
    const Hx = matrixTimesVector(H, x);
    const Aeqx = matrixTimesVector(Aeq, x);
    const Aineqx = matrixTimesVector(Aineq, x);
    const AeqTy = matrixTimesVector(AeqT, y);
    const AineqTz = matrixTimesVector(AineqT, z);

    // Objective
    const f = 0.5 * dot(x, Hx) + dot(c, x); // 0.5 x' H x + c' x

    // Residuals
    let rGrad = add(Hx, c); // Hx + c + Aeq' y - Aineq' z
    if (mEq > 0 ) {
      rGrad = add(rGrad, AeqTy);
    }
    if (mIneq > 0 ) {
      rGrad = subtract(rGrad, AineqTz);
    }
    const rEq = subtract(Aeqx, beq); // Aeq x - beq
    const rIneq = subtract(subtract(Aineqx, s), bineq); // Aineq x - s - bineq
    const rS = subtract(elementwiseVectorProduct(s, z), new Array(mIneq).fill(mu)); // SZe - mu e

    return { f, rGrad, rEq, rIneq, rS };
  }


  // Construct the augmented KKT system
  /*  [ H       Aeq'   Aineq' ]
   *  [ Aeq      0      0     ]
   *  [ Aineq    0   -Z^-1 S  ]
  */
  const m = n + mEq + mIneq;
  const KKT = zeroMatrix(m, m);
  setSubmatrix(KKT, H, 0, 0);
  setSubmatrix(KKT, AeqT, 0, n);
  setSubmatrix(KKT, AineqT, 0, n + mEq);
  setSubmatrix(KKT, Aeq, n, 0);
  setSubmatrix(KKT, Aineq, n + mEq, 0);
  function updateMatrix(s, z) {
    const minusZinvS = negate(elementwiseVectorDivision(s, z));
    setSubdiagonal(KKT, minusZinvS, n + mEq, n + mEq);
  }

  // Define the function for computing the search direction
  function computeSearchDirection(s, z, L, D, rGrad, rEq, rIneq, rS) {
    const rIneqMinusYinvrS = add(rIneq, elementwiseVectorDivision(rS, z)); // Aineq x - s - bineq + Z^-1 (SZe - mue)
    const rhs = negate(rGrad.concat(rEq).concat(rIneqMinusYinvrS));

    // Solve the KKT system
    const d = solveUsingFactorization(L, D, rhs);

    // Extract the search direction components
    const dx = d.slice(0, n);
    const dy = d.slice(n, n + mEq);
    const dz = negate(d.slice(n + mEq, n + mEq + mIneq));
    const ds = negate(elementwiseVectorDivision(add(rS, elementwiseVectorProduct(s, dz)), z)); // -Z^-1 (rS + S dz)

    return { dx, ds, dy, dz };
  }

  // Define the function for computing the step size
  function getMaxStep(v, dv) {
    return v.reduce((m, value, index) => dv[index] < 0 ? Math.min(-value/dv[index], m) : m, 1.0);
  }

  // Initialize primal and dual variables
  const x = new Array(n).fill(1.0);        // Primal variables
  const s = new Array(mIneq).fill(1.0);    // Slack variables for inequality constraints
  const y = new Array(mEq).fill(1.0);    // Multipliers for equality constraints
  const z = new Array(mIneq).fill(1.0); // Multipliers for inequality constraints
  
  function getMu(s, z) {
    return mIneq > 0 ? dot(s, z) / mIneq : 0;
  }

  function getResidualAndGap(s, z, rGrad, rEq, rIneq) {
    const res = norm(rGrad.concat(rEq).concat(rIneq));
    const gap = getMu(s, z);
    return { res, gap };
  }

  // Perform the interior point optimization
  let iter = 0;
  for (; iter < maxIter; iter++) {
    const { f, rGrad, rEq, rIneq, rS } = evalFunc(x, s, y, z, 0);

    // Check the convergence criterion
    const { res, gap } = getResidualAndGap(s, z, rGrad, rEq, rIneq);
    console.log('f: ' + f + ', res: ' + res + ', gap: ' + gap )
    if (res <= tol && gap <= tol) {
      break;
    }

    // Update and factorize KKT matrix
    updateMatrix(s, z);
    [L, D] = symmetricIndefiniteFactorization(KKT);

    // Use the predictor-corrector method

    // Compute affine scaling step
    const { dx : dxAff, ds : dsAff, dy : dyAff, dz : dzAff } = computeSearchDirection(s, z, L, D, rGrad, rEq, rIneq, rS);
    const alphaAffP = getMaxStep(s, dsAff);
    const alphaAffD = getMaxStep(z, dzAff);
    const zAff = Array.from(z);
    const sAff = Array.from(s);
    vectorPlusEqScalarTimesVector(sAff, alphaAffP, dsAff);
    vectorPlusEqScalarTimesVector(zAff, alphaAffD, dzAff);
    const muAff = getMu(zAff, sAff);

    // Compute aggregated centering-corrector direction
    const mu = getMu(s, z);
    const sigma = mu > 0 ? Math.pow(muAff / mu, 3.0) : 0;
    const { rS : rSCenter } = evalFunc(x, s, y, z, sigma * mu);
    const rSCenterCorr = add(elementwiseVectorProduct(dzAff, dsAff), rS);
    const { dx, ds, dy, dz } = computeSearchDirection(s, z, L, D, rGrad, rEq, rIneq, rSCenterCorr);
    const alphaP = getMaxStep(s, ds);
    const alphaD = getMaxStep(z, dz);

    // Update the variables
    const fractionToBoundary = 0.995;
    vectorPlusEqScalarTimesVector(x, fractionToBoundary * alphaP, dx);
    vectorPlusEqScalarTimesVector(s, fractionToBoundary * alphaP, ds);
    vectorPlusEqScalarTimesVector(y, fractionToBoundary * alphaD, dy);
    vectorPlusEqScalarTimesVector(z, fractionToBoundary * alphaD, dz);
  }

  // Return the solution and objective value
  const { f, rGrad, rEq, rIneq, rS } = evalFunc(x, s, y, z, 0);
  const { res, gap } = getResidualAndGap(s, z, rGrad, rEq, rIneq);
  return { x, f, res, gap, iter };
}


// Helper functions for linear algebra operations
function zeroVector(n) {
  return new Array(n).fill(0.0);
}

function zeroMatrix(m, n) {
  return new Array(m).fill().map(() => new Array(n).fill(0.0));
}

function setSubmatrix(M, X, startI, startJ) {
  const m = X.length;
  if (M.length < m + startI) {
    throw new Error('Invalid submatrix row');
  }
  for (let i = 0; i < m; i++) {
    const si = i + startI;
    const n = X[i].length;
    if (M[si].length < n + startJ) {
      throw new Error('Invalid submatrix column');
    }
    for (let j = 0; j < n; j++) {
      M[si][j + startJ] = X[i][j];
    }
  }
}

function setSubdiagonal(M, d, startI, startJ) {
  const m = d.length;
  if (M.length < m + startI) {
    throw new Error('Invalid submatrix row');
  }
  for (let i = 0; i < m; i++) {
    M[i + startI][i + startJ] = d[i];
  }
}

function isVector(x) {
  return Array.isArray(x) && x.every(xi => typeof xi == 'number');
}

function assertIsVector(x, name) {
  if (!isVector(x)) {
    throw new Error('Invalid input type: ' + name + ' must be an array. ' + name + ': ' + x);
  }
}

function assertIsMatrix(A) {
  if (!Array.isArray(A) || A.some(row => !isVector(row))) {
    throw new Error('Invalid input type: A must be a matrix. A: ' + A);
  }
}

function diag(x) {
  assertIsVector(x, 'x');
  m = x.length;
  X = zeroMatrix(m, m);
  for (let i = 0; i < m; i++) {
    X[i][i] = x[i];
  }
  return X;
}

function assertAreEqualLengthVectors(x, y) {
  assertIsVector(x, 'x');
  assertIsVector(y, 'y');

  if (x.length !== y.length) {
    throw new Error('Invalid input shape: x and y must have the same length. x.length = ' + x.length + ', y.length = ' + y.length);
  }
}

function transpose(A) {
  assertIsMatrix(A);
  const m = A.length;
  const n = m > 0 ? A[0].length : 0;
  const B = zeroMatrix(n, m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      B[j][i] = A[i][j];
    }
  }
  return B;
}

function negate(x) {
  assertIsVector(x, 'x');
  return x.map(value => -value);
}

function elementwiseVectorProduct(x, y) {
  assertAreEqualLengthVectors(x, y);
  return x.map((value, index) => value * y[index]);
}

function elementwiseVectorDivision(x, y) {
  assertAreEqualLengthVectors(x, y);
  return x.map((value, index) => value / y[index]);
}

function vectorPlusEqScalarTimesVector(x, s, y) {
  assertAreEqualLengthVectors(x, y);
  for (let i = 0; i < x.length; i++) {
    x[i] += s * y[i];
  }
}

function matrixTimesVector(A, x) {
  assertIsMatrix(A);
  A.every(row => assertAreEqualLengthVectors(row, x));
  return A.map(ai => dot(ai, x));
}

function add(x, y) {
  assertAreEqualLengthVectors(x, y);
  return x.map((value, index) => value + y[index]);
}

function addVectors(...vectors) {
  return vectors.reduce((acc, vec) => add(acc, vec));
}

function subtract(x, y) {
  assertAreEqualLengthVectors(x, y);
  return x.map((value, index) => value - y[index]);
}

function subtractVectors(...vectors) {
  return vectors.reduce((acc, vec) => subtract(acc, vec));
}

function norm(x) {
  return Math.sqrt(dot(x, x));
}

function dot(x, y) {
  assertAreEqualLengthVectors(x, y);
  return x.reduce((sum, value, index) => sum + value * y[index], 0);
}

function symmetricIndefiniteFactorization(A) {
  const n = A.length;
  const L = zeroMatrix(n, n);
  const D = zeroVector(n);

  for (let i = 0; i < n; i++) {
    // Compute the (i,i) entry of D
    let d_ii = A[i][i];
    for (let k = 0; k < i; k++) {
      d_ii -= L[i][k] ** 2 * D[k];
    }

    // Check for singularity
    if (d_ii === 0) {
      throw new Error('Matrix is singular');
    }

    D[i] = d_ii;

    // Compute the entries of L
    for (let j = i + 1; j < n; j++) {
      let l_ij = A[i][j];
      for (let k = 0; k < i; k++) {
        l_ij -= L[i][k] * D[k] * L[j][k];
      }
      L[j][i] = l_ij / d_ii;
    }
  }

  return [L, D];
}

function solveUsingFactorization(L, D, b) {
  const n = L.length;
  const x = new Array(n).fill(0);
  const y = new Array(n).fill(0);

  // Forward substitution: solve Ly = b
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < i; j++) {
      sum += L[i][j] * y[j];
    }
    y[i] = b[i] - sum;
  }

  // Backward substitution: solve L^Tx = y
  for (let i = n - 1; i >= 0; i--) {
    let sum = 0;
    for (let j = i + 1; j < n; j++) {
      sum += L[j][i] * x[j];
    }
    x[i] = (y[i] - sum) / D[i];
  }

  return x;
}

// Parsing of objectives and constraints

function solveQP(Q, c, Aeq, beq, Aineq, bineq, variables = []) {
  let solutionElement = document.getElementById("solution");
  
  try {
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

  } catch (error) {
    solutionElement.innerHTML = `Error ${error.lineNumber}: ${error.message}`;
  }
}

// Functions relating to buttons on the html page

function solve() {
  const objective = document.getElementById("objective").value;
  const table = document.getElementById("optimization-problem");
  let constraints = [];
  for (let i = 2; i < table.rows.length; i++) {
    const constraint = document.getElementById(`constraint-${i}`).textContent;
    constraints.push(constraint);
  }
  try {
    const { Q, c, variables } = parseObjective(objective);
    const { Aeq, beq, Aineq, bineq } = parseConstraints(variables, constraints)
    solveQP(Q, c, Aeq, beq, Aineq, bineq, variables);
  } catch (error) {
    let solutionElement = document.getElementById("solution");
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

document.getElementById("solve").addEventListener("click", solve);

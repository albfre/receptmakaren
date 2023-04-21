const consideredNaringsvarden = ["Energi (kJ)", "Fett", "Summa mättade fettsyror", "Kolhydrater", "Socker totalt", "Fibrer", "Protein", "Salt"];
const displayNaringsvarden = ["Energi (kJ)", "Fett (g)", "- varav mättat fett (g)" , "Kolhydrater (g)", "Socker (g)", "Fibrer (g)", "Protein (g)", "Salt (g)"]
const shortNaringsvarden = ["energi", "fett", "mattat-fett", "kolhydrat", "socker", "fibrer", "protein", "salt"];

async function parseLivsmedelFromXML(filename) {
  const parser = new DOMParser();
  const data = await readFile(filename);
  const xmlDoc = parser.parseFromString(data, "text/xml");

  const result = [];
  for (const entry of xmlDoc.getElementsByTagName("Livsmedel")) {
    const livsmedelsnamn = entry.getElementsByTagName("Namn")[0].textContent;
    const naringsvarden = {};
    for (const naringsvarde of entry.getElementsByTagName("Naringsvarde")) {
      const namn = naringsvarde.getElementsByTagName("Namn")[0].textContent;
      const varde = naringsvarde.getElementsByTagName("Varde")[0].textContent.replace(/\s/g, "").replace(",", ".");
      const enhet = naringsvarde.getElementsByTagName("Enhet")[0].textContent;
      if (consideredNaringsvarden.includes(namn)) {
        const val = parseFloat(varde);
        if (enhet !== (namn === "Energi (kJ)" ? "kJ" : "g")) {
          throw new Error(`Unknown unit: ${enhet}`);
        }
        naringsvarden[namn] = [parseFloat(varde), enhet];
      }
    }
    result.push({ livsmedelsnamn, naringsvarden });
  }
  result.sort(function(a, b) { return a.livsmedelsnamn.localeCompare(b.livsmedelsnamn, 'sv'); });

  //console.log(result);
  return result;
}

async function readFile(filename) {
  const response = await fetch(`./data/${filename}`, {cache: "no-store"});
  if (!response.ok) {
    throw new Error("Could not read file");
  }
  const data = await response.text();
  return data;
}

async function populateFoods() {
  const livsmedel = await parseLivsmedelFromXML("compressed");
  const foodsSelect = document.querySelector("#foods");
  foodsSelect.innerHTML = "";
  for (const { livsmedelsnamn, naringsvarden } of livsmedel) {
    const option = document.createElement("option");
    option.textContent = livsmedelsnamn;
    const mappedArray = Object.entries(naringsvarden).map(([key, value]) => `${key}, ${value[0]}`);
    option.setAttribute("data-value", mappedArray.join(", "));
    foodsSelect.appendChild(option);
  }

  const preselect = ["Havregryn fullkorn", "Socker", "Rapsolja", "Ljus sirap", "Veteflingor ångprep. fullkorn", "Vetekli", "Salt m. jod"];
  const selected = document.querySelector("#selected-foods");
  for (const pre of preselect) {
    for (const option of foodsSelect.options) {
      if (option.value === pre) {
        selected.appendChild(option);
      }
    }
  }
}

async function populateDivs(prefix, valType) {
  const div = document.querySelector(`#${prefix}-div`);
  for (let i = 0; i < shortNaringsvarden.length; i++) {
    // Create label
    const label = document.createElement("label");
    label.textContent = displayNaringsvarden[i] + ":";
    label.classList.add("label-width");
    div.appendChild(label);

    // Create label or input
    const val = document.createElement(valType);
    if (valType === "input") {
      val.type = "text";
      const valMap = {0: 1800, 1: 15, 2: 2, 3: 63, 4: 20, 5: 7.8, 6: 8.9, 7: 0.5 };
      if (i in valMap) {
        val.value = valMap[i];
      }
    }
    val.classList.add("input-width");
    val.setAttribute("id", prefix + "-" + shortNaringsvarden[i]);
    div.appendChild(val);

    // Create br
    const br = document.createElement("br");
    div.appendChild(br);
  }
}

async function populate() {
  populateFoods();
  populateDivs("mimic", "input");
  populateDivs("vald", "label");
  populateDivs("result", "label");
}

function addFoods() {
  const selectedFoods = [];
  const foodsSelect = document.querySelector("#foods");
  for (const option of foodsSelect.selectedOptions) {
    const foodName = option.textContent;
    selectedFoods.push([foodName, option.dataset.value]);
  }
  const selectedFoodsList = document.querySelector("#selected-foods");
  for (const food of selectedFoods) {
    const option = document.createElement("option");
    option.textContent = food[0];
    option.setAttribute("data-value", food[1]);
    selectedFoodsList.appendChild(option);
  }
}

function removeFoods() {
  const select = document.querySelector("#selected-foods");
  for (var i = select.options.length - 1; i >= 0; i--) {
    if (select.options[i].selected) {
      select.remove(i);
    }
  }
}

function searchFoods() {
  const input = document.getElementById("searchInput");
  const filter = input.value.toUpperCase();
  const foodsList = document.querySelector("#foods");
  const foods = Array.from(foodsList.options);

  // Sort options by whether they start with filter
  foods.sort((a, b) => {
    const atext = a.textContent.toUpperCase();
    const btext = b.textContent.toUpperCase();
    const a1 = atext.startsWith(filter);
    const b1 = btext.startsWith(filter);
    if (a1 && !b1) return -1;
    if (b1 && !a1) return 1;
    return atext.localeCompare(btext, 'sv');
  });

  // Remove all options from the select element
  while (foodsList.options.length > 0) {
    foodsList.remove(0);
  }

  // Add the sorted options back to the select element and hide the ones not containing filter
  for (const option of foods) {
    option.style.display = option.textContent.toUpperCase().indexOf(filter) > -1 ? '' : 'none';
    foodsList.add(option);
  }
}

function showFood() {
  const foodsSelect = document.getElementById("selected-foods");
  for (const option of foodsSelect.selectedOptions) {
    const name = option.value;
    const varde = parseSelectedFoodOption(option)
    setOutput(varde, 'vald');
    break;
  }
}

function up() {
  const select = document.querySelector("#selected-foods");
  for (var i = 0; i < select.options.length; i++) {
    if (select.options[i].selected && i > 0 && !select.options[i - 1].selected) {
      select.insertBefore(select.options[i], select.options[i - 1]);
    }
  }
}

function down() {
  const select = document.querySelector("#selected-foods");
  for (var i = select.options.length - 1; i >= 0; i--) {
    if (select.options[i].selected && i + 1 < select.options.length && !select.options[i + 1].selected) {
      select.insertBefore(select.options[i + 1], select.options[i]);
    }
  }
}

document.getElementById("add-food").addEventListener("click", addFoods);
document.getElementById("remove-food").addEventListener("click", removeFoods);
document.getElementById("up-btn").addEventListener("click", up);
document.getElementById("down-btn").addEventListener("click", down);
window.onload = populate;
